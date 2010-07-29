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

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <string>

#include <cuda_wrapper.hpp>
#include <halmd/deprecated/mdsim/backend/exception.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/forces/lj.hpp>
#include <halmd/mdsim/gpu/forces/lj_kernel.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace boost::assign;
using namespace boost::fusion;
using namespace boost::numeric::ublas;

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void lj<dimension, float_type>::options(po::options_description& desc)
{
    po::options_description group("Lennard-Jones potential");
    group.add_options()
        ("cutoff", po::value<boost::array<float, 3> >()->default_value(list_of(2.5f)(2.5f)(2.5f)),
         "truncate potential at cutoff radius")
        ("epsilon", po::value<boost::array<float, 3> >()->default_value(list_of(1.0f)(1.5f)(0.5f)),
         "potential well depths AA,AB,BB")
        ("sigma", po::value<boost::array<float, 3> >()->default_value(list_of(1.0f)(0.8f)(0.88f)),
         "collision diameters AA,AB,BB")
        ;
    desc.add(group);
}

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void lj<dimension, float_type>::select(po::options const& vm)
{
    if (vm["force"].as<std::string>() != "lj") {
        throw unsuitable_module("mismatching option force");
    }
}

/**
 * Initialize Lennard-Jones potential parameters
 */
template <int dimension, typename float_type>
lj<dimension, float_type>::lj(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  // allocate potential parameters
  , epsilon_(scalar_matrix<float_type>(particle->ntype, particle->ntype, 1))
  , sigma_(scalar_matrix<float_type>(particle->ntype, particle->ntype, 1))
  , r_cut_sigma_(particle->ntype, particle->ntype)
  , r_cut_(particle->ntype, particle->ntype)
  , rr_cut_(particle->ntype, particle->ntype)
  , sigma2_(particle->ntype, particle->ntype)
  , en_cut_(particle->ntype, particle->ntype)
  , g_ljparam_(epsilon_.data().size())
{
    // register module runtime accumulators
    profiler->register_map(runtime_);

    // parse deprecated options
    boost::array<float, 3> epsilon = vm["epsilon"].as<boost::array<float, 3> >();
    boost::array<float, 3> sigma = vm["sigma"].as<boost::array<float, 3> >();
    boost::array<float, 3> r_cut_sigma;
    try {
        r_cut_sigma = vm["cutoff"].as<boost::array<float, 3> >();
    }
    catch (boost::bad_any_cast const&) {
        // backwards compatibility
        std::fill(r_cut_sigma.begin(), r_cut_sigma.end(), vm["cutoff"].as<float>());
    }
    for (size_t i = 0; i < std::min(particle->ntype, 2U); ++i) {
        for (size_t j = i; j < std::min(particle->ntype, 2U); ++j) {
            epsilon_(i, j) = epsilon[i + j];
            sigma_(i, j) = sigma[i + j];
            r_cut_sigma_(i, j) = r_cut_sigma[i + j];
        }
    }

    // precalculate derived parameters
    for (size_t i = 0; i < particle->ntype; ++i) {
        for (size_t j = i; j < particle->ntype; ++j) {
            r_cut_(i, j) = r_cut_sigma_(i, j) * sigma_(i, j);
            rr_cut_(i, j) = std::pow(r_cut_(i, j), 2);
            sigma2_(i, j) = std::pow(sigma_(i, j), 2);
            // energy shift due to truncation at cutoff length
            float_type rri_cut = std::pow(r_cut_sigma_(i, j), -2);
            float_type r6i_cut = rri_cut * rri_cut * rri_cut;
            en_cut_(i, j) = 4 * epsilon_(i, j) * r6i_cut * (r6i_cut - 1);
        }
    }

    LOG("potential well depths: ε = " << epsilon_);
    LOG("potential pair separation: σ = " << sigma_);
    LOG("potential cutoff length: r = " << r_cut_sigma_);
    LOG("potential cutoff energy: U = " << en_cut_);

    cuda::host::vector<float4> ljparam(g_ljparam_.size());
    for (size_t i = 0; i < ljparam.size(); ++i) {
        fixed_vector<float, 4> p;
        p[lj_kernel::EPSILON] = epsilon_.data()[i];
        p[lj_kernel::RR_CUT] = rr_cut_.data()[i];
        p[lj_kernel::SIGMA2] = sigma2_.data()[i];
        p[lj_kernel::EN_CUT] = en_cut_.data()[i];
        ljparam[i] = p;
    }
    cuda::copy(ljparam, g_ljparam_);
    cuda::copy(static_cast<vector_type>(box->length()), get_lj_kernel<dimension>().box_length);
/* FIXME
    // initialise CUDA symbols
    typedef lj_wrapper<dimension> _gpu;

    get_lj_kernel<dimension>().r = // cuda::texture<float4>
    get_lj_kernel<dimension>().box_length =  // cuda::symbol<vector_type>
    get_lj_kernel<dimension>().neighbour_size = // cuda::symbol<unsigned int> ;
    get_lj_kernel<dimension>().neighbour_stride = // cuda::symbol<unsigned int> ;
    get_lj_kernel<dimension>().ljparam = // cuda::texture<float4> ;
*/
}

/**
 * Compute Lennard-Jones forces
 */
template <int dimension, typename float_type>
void lj<dimension, float_type>::compute()
{
#ifdef USE_FORCE_DSFUN
#endif /* HALMD_VARIANT_FORCE_DSFUN */

    scoped_timer<timer> timer_(at_key<compute_>(runtime_));
    cuda::copy(particle->neighbour_size, get_lj_kernel<dimension>().neighbour_size);
    cuda::copy(particle->neighbour_stride, get_lj_kernel<dimension>().neighbour_stride);
    get_lj_kernel<dimension>().r.bind(particle->g_r);
    get_lj_kernel<dimension>().ljparam.bind(g_ljparam_);
    cuda::configure(particle->dim.grid, particle->dim.block);
    get_lj_kernel<dimension>().compute(
        particle->g_f
      , particle->g_neighbour
      , thermodynamics->g_en_pot_
      , thermodynamics->g_virial_
    );
    cuda::thread::synchronize();
//                   thermodynamics->en_pot, thermodynamics->virial);

/*    // initialize particle forces to zero
    std::fill(particle->f.begin(), particle->f.end(), 0);

    // potential energy
    float_type& en_pot_ = thermodynamics->en_pot_;
    en_pot_ = 0;

    // virial equation sum
    typedef typename _Base::thermodynamics_type thermodynamics_type;
    typedef typename thermodynamics_type::virial_type virial_type;
    std::vector<virial_type>& virial_ = thermodynamics->virial_;
    std::fill(virial_.begin(), virial_.end(), 0);

    for (size_t i = 0; i < particle->nbox; ++i) {
        // calculate pairwise Lennard-Jones force with neighbour particles
        BOOST_FOREACH(size_t j, particle->neighbour[i]) {
            // particle distance vector
            vector_type r = particle->r[i] - particle->r[j];
            box->reduce_periodic(r);
            // particle types
            size_t a = particle->type[i];
            size_t b = particle->type[j];
            // squared particle distance
            float_type rr = inner_prod(r, r);

            // truncate potential at cutoff length
            if (rr >= rr_cut_(a, b))
                continue;

            // compute Lennard-Jones force in reduced units
            float_type sigma2 = sigma2_(a, b);
            float_type rri = sigma2 / rr;
            float_type r6i = rri * rri * rri;
            float_type epsilon = epsilon_(a, b);
            float_type fval = 48 * rri * r6i * (r6i - 0.5) * (epsilon / sigma2);
            float_type en_pot = 4 * epsilon * r6i * (r6i - 1) - en_cut_(a, b);

            // optionally smooth potential yielding continuous 2nd derivative
            // FIXME test performance of template versus runtime bool
            if (smooth) {
                smooth->compute(std::sqrt(rr), r_cut_(a, b), fval, en_pot);
            }

            // add force contribution to both particles
            particle->f[i] += r * fval;
            particle->f[j] -= r * fval;

            // add contribution to potential energy
            en_pot_ += en_pot;

            // add contribution to virial
            virial_type vir = 0.5f * fval * thermodynamics_type::virial_tensor(rr, r);
            virial_[a] += vir;
            virial_[b] += vir;
        }
    }

    en_pot_ /= particle->nbox;
    for (size_t i = 0; i < virial_.size(); ++i) {
        virial_[i] /= particle->ntypes[i];
    }

    // ensure that system is still in valid state
    if (std::isinf(en_pot_)) {
        throw potential_energy_divergence();
    }*/
}

// explicit instantiation
template class lj<3, float>;
template class lj<2, float>;

}}} // namespace mdsim::gpu::forces

template class module<mdsim::gpu::forces::lj<3, float> >;
template class module<mdsim::gpu::forces::lj<2, float> >;

} // namespace halmd
