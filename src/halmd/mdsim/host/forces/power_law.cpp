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

#include <halmd/io/logger.hpp>
#include <halmd/math/pow.hpp>
#include <halmd/mdsim/backend/exception.hpp>
#include <halmd/mdsim/host/forces/power_law.hpp>
#include <halmd/utility/module.hpp>

using namespace boost;
using namespace boost::assign;
using namespace boost::numeric::ublas;

using namespace halmd::math;

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**boost/algorithm/string/predicate.hpp>
 * Assemble module options
 */
template <int dimension, typename float_type>
void power_law<dimension, float_type>::options(po::options_description& desc)
{
    po::options_description group("Power law potential");
    group.add_options()
        ("index", po::value<int>()->default_value(12),
         "index of soft power-law potential")
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
void power_law<dimension, float_type>::resolve(po::options const& vm)
{
    if (vm["force"].as<std::string>() != "power-law") {
        throw unsuitable_module("mismatching option force");
    }
}

/**
 * Initialize Lennard-Jones potential parameters
 */
template <int dimension, typename float_type>
power_law<dimension, float_type>::power_law(po::options const& vm)
  : _Base(vm)
  // allocate potential parameters
  , index_(vm["index"].as<int>())
  , epsilon_(scalar_matrix<float_type>(particle->ntype, particle->ntype, 1))
  , sigma_(scalar_matrix<float_type>(particle->ntype, particle->ntype, 1))
  , r_cut_sigma_(particle->ntype, particle->ntype)
  , r_cut_(particle->ntype, particle->ntype)
  , rr_cut_(particle->ntype, particle->ntype)
  , sigma2_(particle->ntype, particle->ntype)
  , en_cut_(particle->ntype, particle->ntype)
{
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
            en_cut_(i, j) = epsilon_(i, j) * std::pow(1 / r_cut_sigma_(i, j), static_cast<int>(index_));
        }
    }

    LOG("potential power-law index: n = " << index_);
    LOG("potential strength: ε = " << epsilon_);
    LOG("potential pair separation: σ = " << sigma_);
    LOG("potential cutoff length: r = " << r_cut_sigma_);
    LOG("potential cutoff energy: U = " << en_cut_);
}

/**
 * Compute power-law forces
 *
 * call index-dependent template implementations
 * for efficiency of pow() function
 */
template <int dimension, typename float_type>
void power_law<dimension, float_type>::compute()
{
    switch (index_) {
        case 6:  compute_impl<6>();  break;
        case 12: compute_impl<12>(); break;
        case 24: compute_impl<24>(); break;
        case 48: compute_impl<48>(); break;
        default:
            LOG_WARNING_ONCE("Using non-optimised force routine for index " << index_);
            compute_impl<0>();
            break;
    }
}

template <int dimension, typename float_type>
template <int index>
void power_law<dimension, float_type>::compute_impl()
{
    // initialize particle forces to zero
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

            // compute power-law force in reduced units,
            // choose arbitrary index_ if template parameter index = 0
            float_type rni;
            if (index > 0) {
                rni = pow<index>(sigma_(a, b) / std::sqrt(rr));
            }
            else {
                rni = std::pow(sigma_(a, b) / std::sqrt(rr), index_);
            }
            float_type en_pot = epsilon_(a, b) * rni;      // U(r)
            float_type fval = (index > 0 ? index : index_) * en_pot / rr;
                                                           // F(r) / r
            en_pot -= en_cut_(a, b);                       // shift potential

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
    }
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class power_law<3, double>;
template class power_law<2, double>;
#else
template class power_law<3, float>;
template class power_law<2, float>;
#endif

}}} // namespace mdsim::host::forces

#ifndef USE_HOST_SINGLE_PRECISION
template class module<mdsim::host::forces::power_law<3, double> >;
template class module<mdsim::host::forces::power_law<2, double> >;
#else
template class module<mdsim::host::forces::power_law<3, float> >;
template class module<mdsim::host::forces::power_law<2, float> >;
#endif

} // namespace halmd
