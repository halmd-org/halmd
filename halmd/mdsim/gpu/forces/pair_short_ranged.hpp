/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_SHORT_RANGED_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_SHORT_RANGED_HPP

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/forces/pair_short_ranged_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{

/**
 * class template for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type>
class pair_short_ranged
  : public mdsim::gpu::force<dimension, float_type>
{
public:
    typedef mdsim::gpu::force<dimension, float_type> _Base;
    typedef typename _Base::matrix_type matrix_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_stress_tensor_type gpu_stress_tensor_type;
    typedef gpu::particle<dimension, float> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef utility::profiler profiler_type;

    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_short_ranged_wrapper<dimension, gpu_potential_type> gpu_wrapper;

    boost::shared_ptr<potential_type> potential;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    inline pair_short_ranged(
        boost::shared_ptr<potential_type> potential
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
    );
    inline virtual void compute();
    inline void register_runtimes(profiler_type& profiler);

    //! returns potential cutoff distance
    virtual matrix_type const& cutoff()
    {
        return potential->r_cut();
    }

    //! returns potential energies of particles
    virtual cuda::vector<float> const& potential_energy()
    {
        return g_en_pot_;
    }

    /** potential part of stress tensors of particles */
    virtual cuda::vector<gpu_stress_tensor_type> const& potential_stress()
    {
        return g_stress_pot_;
    }

    // module runtime accumulator descriptions
    HALMD_PROFILE_TAG( compute_, "force computation" );

private:
    /** potential energy for each particle */
    cuda::vector<float> g_en_pot_;
    /** potential part of stress tensor for each particle */
    cuda::vector<gpu_stress_tensor_type> g_stress_pot_;

    boost::fusion::map<
        boost::fusion::pair<compute_, accumulator<double> >
    > runtime_;
};

template <int dimension, typename float_type, typename potential_type>
pair_short_ranged<dimension, float_type, potential_type>::pair_short_ranged(
    boost::shared_ptr<potential_type> potential
  , boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<box_type> box
  // FIXME , boost::shared_ptr<smooth_type> smooth
)
  // dependency injection
  : potential(potential)
  , particle(particle)
  , box(box)
  // memory allocation
  , g_en_pot_(particle->dim.threads())
  , g_stress_pot_(particle->dim.threads())
{
    cuda::copy(static_cast<vector_type>(box->length()), gpu_wrapper::kernel.box_length);
}

/**
 * Compute pair forces, potential energy, and potential part of stress tensor
 */
template <int dimension, typename float_type, typename potential_type>
void pair_short_ranged<dimension, float_type, potential_type>::compute()
{
#ifdef USE_FORCE_DSFUN
#endif /* HALMD_VARIANT_FORCE_DSFUN */

    scoped_timer<timer> timer_(boost::fusion::at_key<compute_>(runtime_));

    cuda::copy(particle->neighbour_size, gpu_wrapper::kernel.neighbour_size);
    cuda::copy(particle->neighbour_stride, gpu_wrapper::kernel.neighbour_stride);
    gpu_wrapper::kernel.r.bind(particle->g_r);
    potential->get_kernel_param().bind(potential->g_param());

    cuda::configure(particle->dim.grid, particle->dim.block);
    gpu_wrapper::kernel.compute(particle->g_f, particle->g_neighbour, g_en_pot_, g_stress_pot_);
    cuda::thread::synchronize();
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type, typename potential_type>
void pair_short_ranged<dimension, float_type, potential_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

}}} // namespace mdsim::gpu::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_LJ_HPP */
