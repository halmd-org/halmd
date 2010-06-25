/* Lennard-Jones fluid simulation using CUDA
 *
 * Copyright © 2008-2010  Peter Colberg
 *                        Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_MDSIM_LJFLUID_GPU_SQUARE_HPP
#define HALMD_MDSIM_LJFLUID_GPU_SQUARE_HPP

#include <halmd/algorithm/reduce.hpp>
#include <halmd/mdsim/ljfluid_gpu_base.hpp>
#include <halmd/mdsim/gpu/lattice.hpp>

#define foreach BOOST_FOREACH

namespace halmd
{

template <typename ljfluid_impl, int dimension>
class ljfluid;

template <int dimension>
class ljfluid<ljfluid_impl_gpu_square, dimension>
    : public ljfluid_gpu_base<ljfluid_impl_gpu_square, dimension>
{
public:
    typedef ljfluid_gpu_base<ljfluid_impl_gpu_square, dimension> _Base;
    typedef gpu::ljfluid<ljfluid_impl_gpu_square, dimension> _gpu;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_vector_type gpu_vector_type;
    typedef typename _Base::host_sample_type host_sample_type;
    typedef typename _Base::gpu_sample_type gpu_sample_type;
    typedef gpu_sample_type trajectory_sample_type;
    typedef boost::variant<host_sample_type, gpu_sample_type> trajectory_sample_variant;
    typedef typename _Base::energy_sample_type energy_sample_type;
    typedef typename _Base::virial_tensor virial_tensor;

public:
    /** set number of particles in system */
    template <typename T>
    void particles(T const& value);
    /** set number of CUDA execution threads */
    void threads(unsigned int value);

    /** restore system state from phase space sample */
    void state(host_sample_type& sample, float_type box);
    /** rescale particle velocities */
    void rescale_velocities(double coeff);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(float_type temp);

    /** MD integration step on GPU */
    void mdstep();
    /** sample phase space on host */
    void sample(host_sample_type& sample) const;
    /** sample phase space on GPU */
    void sample(gpu_sample_type& sample) const;
    /** sample thermodynamic equilibrium properties */
    void sample(energy_sample_type& sample) const;

    /** returns number of particles */
    unsigned int particles() const { return npart; }
    /** get number of CUDA execution threads */
    unsigned int threads() const { return dim_.threads_per_block(); }

private:
    /** assign particle positions */
    void assign_positions();
    /** generate Maxwell-Boltzmann distributed velocities */
    void boltzmann(float temp);
    /** first leapfrog step of integration of differential equations of motion */
    void velocity_verlet();
    /** Lennard-Jones force calculation */
    void update_forces();

private:
    using _Base::box_;
    using _Base::density_;
    using _Base::dim_;
    using _Base::m_times;
    using _Base::mixture_;
    using _Base::mpart;
    using _Base::npart;
    using _Base::potential_;
    using _Base::r_cut;
    using _Base::timestep_;
    using _Base::thermostat_steps;
    using _Base::thermostat_count;
    using _Base::thermostat_temp;

    /** CUDA execution dimensions for phase space sampling */
    std::vector<cuda::config> dim_sample;

    using _Base::reduce_squared_velocity;
    using _Base::reduce_velocity;
    using _Base::reduce_en;
    using _Base::reduce_virial;
    using _Base::reduce_helfand;

    /** system state in page-locked host memory */
    struct {
        /** tagged periodically reduced particle positions */
        cuda::host::vector<float4> mutable r;
        /** periodic box traversal vectors */
        cuda::host::vector<gpu_vector_type> mutable R;
        /** particle velocities */
        cuda::host::vector<gpu_vector_type> mutable v;
    } h_part;

    /** system state in global device memory */
    struct {
        /** tagged periodically reduced particle positions */
        cuda::vector<float4> r;
        /** periodic box traversal vectors */
        cuda::vector<gpu_vector_type> R;
        /** particle velocities */
        cuda::vector<gpu_vector_type> v;
        /** particle forces */
        cuda::vector<gpu_vector_type> f;
        /** particle tags */
        cuda::vector<unsigned int> tag;
        /** potential energies per particle */
        cuda::vector<float> en;
        /** virial equation sums per particle */
        cuda::vector<gpu_vector_type> virial;
        /** time integral of virial stress tensor to calculate Helfand moment */
        cuda::vector<gpu_vector_type> helfand;
    } g_part;

};

template <int dimension>
template <typename T>
void ljfluid<ljfluid_impl_gpu_square, dimension>::particles(T const& value)
{
    _Base::particles(value);

    // allocate page-locked host memory for system state
    try {
        h_part.r.resize(npart);
        h_part.R.resize(npart);
        h_part.v.resize(npart);
        // particle forces reside only in GPU memory
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to allocate page-locked host memory for system state");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::threads(unsigned int value)
{
    _Base::threads(value);

    try {
#ifdef USE_VERLET_DSFUN
        LOG("using double-single arithmetic in Verlet integration");
        g_part.r.reserve(2 * dim_.threads());
#else
        g_part.r.reserve(dim_.threads());
#endif
        // allocate sufficient memory for binary mixture sampling
        for (size_t n = 0, i = 0; n < npart; n += mpart[i], ++i) {
            cuda::config dim((mpart[i] + threads() - 1) / threads(), threads());
            g_part.r.reserve(n + dim.threads());
            dim_sample.push_back(dim);
        }
        g_part.r.resize(npart);
        g_part.R.reserve(g_part.r.capacity());
        g_part.R.resize(npart);
        g_part.v.reserve(g_part.r.capacity());
        g_part.v.resize(npart);
        g_part.f.reserve(dim_.threads());
        g_part.f.resize(npart);
        g_part.tag.reserve(dim_.threads());
        g_part.tag.resize(npart);
        g_part.en.reserve(dim_.threads());
        g_part.en.resize(npart);
        g_part.virial.reserve(dim_.threads());
        g_part.virial.resize(npart);
        g_part.helfand.reserve(dim_.threads());
        g_part.helfand.resize(npart);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to allocate global device memory for system state");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::state(host_sample_type& sample, float_type box)
{
    _Base::state(sample, box, h_part.r, h_part.v);
#ifdef USE_VERLET_DSFUN
    cuda::memset(g_part.r, 0, g_part.r.capacity());
#endif
    cuda::copy(h_part.r, g_part.r);
    assign_positions();
#ifdef USE_VERLET_DSFUN
    cuda::memset(g_part.v, 0, g_part.v.capacity());
#endif
    cuda::copy(h_part.v, g_part.v);
    // init accumulator for Helfand moment
    cuda::memset(g_part.helfand, 0, g_part.helfand.capacity());
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::rescale_velocities(double coeff)
{
    LOG("rescaling velocities with coefficient: " << coeff);
    _Base::rescale_velocities(g_part.v, coeff, dim_);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::lattice()
{
#ifdef USE_VERLET_DSFUN
    cuda::memset(g_part.r, 0, g_part.r.capacity());
#endif
    // place particles on an fcc lattice
    _Base::lattice(g_part.r);
    // randomly permute particle coordinates for binary mixture
    _Base::random_permute(g_part.r);
    // init accumulator for Helfand moment
    cuda::memset(g_part.helfand, 0, g_part.helfand.capacity());

    assign_positions();
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::temperature(float_type temp)
{
    LOG("initialising velocities from Boltzmann distribution at temperature: " << temp);

    boost::array<high_resolution_timer, 2> timer;
    cuda::thread::synchronize();
    timer[0].record();
    try {
        boltzmann(temp);
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to compute Boltzmann distributed velocities on GPU");
    }
    timer[1].record();
    m_times["boltzmann"] += timer[1] - timer[0];
}

/**
 * MD integration step on GPU
 */
template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::mdstep()
{
    boost::array<high_resolution_timer, 5> timer;
    cuda::thread::synchronize();
    timer[1].record();

    // first leapfrog step of integration of differential equations of motion
    try {
        velocity_verlet();
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to stream first leapfrog step on GPU");
    }
    timer[2].record();

    // Lennard-Jones force calculation
    try {
        update_forces();
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to stream force calculation on GPU");
    }
    timer[3].record();

    // heat bath coupling
    if (thermostat_steps && ++thermostat_count > thermostat_steps) {
        try {
            boltzmann(thermostat_temp);
            cuda::thread::synchronize();
        }
        catch (cuda::error const& e) {
            LOG_ERROR("CUDA: " << e.what());
            throw exception("failed to compute Boltzmann distributed velocities on GPU");
        }
    }
    timer[4].record();

    // potential energy sum calculation
    try {
        reduce_en(g_part.en);
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to stream potential energy sum calculation on GPU");
    }
    timer[0].record();

    // CUDA time for MD integration step
    m_times["mdstep"] += timer[0] - timer[1];
    // GPU time for velocity-Verlet integration
    m_times["velocity_verlet"] += timer[2] - timer[1];
    // GPU time for Lennard-Jones force update
    m_times["update_forces"] += timer[3] - timer[2];
    // GPU time for potential energy sum calculation
    m_times["potential_energy"] += timer[0] - timer[4];

    if (thermostat_steps && thermostat_count > thermostat_steps) {
        // reset MD steps since last heatbath coupling
        thermostat_count = 0;
        // GPU time for Maxwell-Boltzmann distribution
        m_times["boltzmann"] += timer[4] - timer[3];
    }

    if (!std::isfinite(reduce_en.value())) {
        throw potential_energy_divergence();
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::sample(host_sample_type& sample) const
{
    typedef typename host_sample_type::value_type sample_type;
    typedef typename sample_type::position_sample_vector position_sample_vector;
    typedef typename sample_type::position_sample_ptr position_sample_ptr;
    typedef typename sample_type::velocity_sample_vector velocity_sample_vector;
    typedef typename sample_type::velocity_sample_ptr velocity_sample_ptr;

    boost::array<high_resolution_timer, 2> timer;
    cuda::thread::synchronize();
    timer[0].record();
    try {
        cuda::copy(g_part.r, h_part.r);
        cuda::copy(g_part.R, h_part.R);
        cuda::copy(g_part.v, h_part.v);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to copy MD integration step results from GPU to host");
    }
    timer[1].record();
    m_times["sample_memcpy"] += timer[1] - timer[0];

    for (size_t n = 0, i = 0; n < npart; ++i) {
        // allocate memory for trajectory sample
        position_sample_ptr r(new position_sample_vector);
        r->reserve(mpart[i]);
        velocity_sample_ptr v(new velocity_sample_vector);
        v->reserve(mpart[i]);
        sample.push_back(sample_type(r, v));
        // assign particle positions and velocities of homogenous type
        for (size_t j = 0; j < mpart[i]; ++j, ++n) {
            r->push_back(h_part.r[n] + box_ * static_cast<vector_type>(h_part.R[n]));
            v->push_back(h_part.v[n]);
        }
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::sample(gpu_sample_type& sample) const
{
    typedef typename gpu_sample_type::value_type sample_type;
    typedef typename sample_type::position_sample_vector position_sample_vector;
    typedef typename sample_type::position_sample_ptr position_sample_ptr;
    typedef typename sample_type::velocity_sample_vector velocity_sample_vector;
    typedef typename sample_type::velocity_sample_ptr velocity_sample_ptr;

    boost::array<high_resolution_timer, 2> timer;
    cuda::thread::synchronize();
    timer[0].record();

    for (size_t n = 0, i = 0; n < npart; n += mpart[i], ++i) {
        // allocate global device memory for phase space sample
        position_sample_ptr r(new position_sample_vector);
        r->reserve(dim_sample[i].threads());
        r->resize(mpart[i]);
        velocity_sample_ptr v(new velocity_sample_vector);
        v->reserve(dim_sample[i].threads());
        v->resize(mpart[i]);
        sample.push_back(sample_type(r, v));
        // sample trajectories
        cuda::configure(dim_sample[i].grid, dim_sample[i].block);
        _gpu::sample(g_part.r.data() + n, g_part.R.data() + n, g_part.v.data() + n, *r, *v);
    }

    cuda::thread::synchronize();
    timer[1].record();
    m_times["sample"] += timer[1] - timer[0];
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::sample(energy_sample_type& sample) const
{
    boost::array<high_resolution_timer, 2> timer;
    cuda::thread::synchronize();

    // mean potential energy per particle
    sample.en_pot = reduce_en.value() / npart;

    // virial tensor trace and off-diagonal elements for particle species
    try {
        timer[0].record();
        if (mixture_ == BINARY) {
            reduce_virial(g_part.virial, g_part.tag, mpart);
            reduce_helfand(g_part.helfand, g_part.tag, mpart);
        }
        else {
            reduce_virial(g_part.virial);
            reduce_helfand(g_part.helfand);
        }
        cuda::thread::synchronize();
        timer[1].record();
        m_times["virial_sum"] += timer[1] - timer[0];
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to calculate virial equation sum on GPU");
    }
    sample.virial = reduce_virial.value();
    for (size_t i = 0; i < sample.virial.size(); ++i) {
        sample.virial[i] /= mpart[i];
    }
    sample.helfand = reduce_helfand.value();
    for (size_t i = 0; i < sample.helfand.size(); ++i) {
        sample.helfand[i] /= mpart[i];
    }

    // mean squared velocity per particle
    try {
        timer[0].record();
        reduce_squared_velocity(g_part.v);
        cuda::thread::synchronize();
        timer[1].record();
        m_times["reduce_squared_velocity"] += timer[1] - timer[0];
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to calculate mean squared velocity on GPU");
    }
    sample.vv = reduce_squared_velocity.value() / npart;

    // mean velocity per particle
    try {
        timer[0].record();
        reduce_velocity(g_part.v);
        cuda::thread::synchronize();
        timer[1].record();
        m_times["reduce_velocity"] += timer[1] - timer[0];
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to calculate mean velocity on GPU");
    }
    sample.v_cm = reduce_velocity.value() / npart;
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::assign_positions()
{
    // assign ascending particle numbers
    _Base::init_tags(g_part.r, g_part.tag);

    try {
        // set periodic box traversal vectors to zero
        cuda::memset(g_part.R, 0);
        // calculate forces
        update_forces();
        // calculate potential energy
        reduce_en(g_part.en);

        // wait for CUDA operations to finish
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to assign particle positions on GPU");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::boltzmann(float temp)
{
#ifdef USE_VERLET_DSFUN
    cuda::memset(g_part.v, 0, g_part.v.capacity());
#endif
    _Base::boltzmann(g_part.v, temp);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::velocity_verlet()
{
    cuda::configure(dim_.grid, dim_.block);
    _gpu::inteq(g_part.r, g_part.R, g_part.v, g_part.f, g_part.virial);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square, dimension>::update_forces()
{
    cuda::configure(dim_.grid, dim_.block, dim_.threads_per_block() * (dimension + 1) * sizeof(int));
    _Base::update_forces(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial, g_part.helfand);
}

} // namespace halmd

#undef foreach

#endif /* ! HALMD_MDSIM_LJFLUID_GPU_SQUARE_HPP */
