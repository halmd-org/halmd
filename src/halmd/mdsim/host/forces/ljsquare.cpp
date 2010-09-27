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

#include <halmd/deprecated/mdsim/backend/exception.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/forces/ljsquare.hpp>
#include <halmd/utility/module.hpp>

using namespace boost;
using namespace boost::assign;
using namespace boost::numeric::ublas;

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void lj_square<dimension, float_type>::select(po::variables_map const& vm)
{
    if (vm["force"].as<std::string>() != "lj_square") {
        throw unsuitable_module("mismatching option force");
    }
}

/**
 * Initialize Lennard-Jones potential parameters
 */
template <int dimension, typename float_type>
lj_square<dimension, float_type>::lj_square(modules::factory& factory, po::variables_map const& vm)
  : _Base(factory, vm)
  // allocate potential parameters
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
            float_type rri_cut = std::pow(r_cut_sigma_(i, j), -2);
            float_type r6i_cut = rri_cut * rri_cut * rri_cut;
            en_cut_(i, j) = 4 * epsilon_(i, j) * r6i_cut * (r6i_cut - 1);
        }
    }

    LOG("potential well depths: ε = " << epsilon_);
    LOG("potential pair separation: σ = " << sigma_);
    LOG("potential cutoff length: r = " << r_cut_sigma_);
    LOG("potential cutoff energy: U = " << en_cut_);
}

/**
 * Compute Lennard-Jones forces
 */
template <int dimension, typename float_type>
void lj_square<dimension, float_type>::compute()
{
    // initialise particle forces to zero
    std::fill(particle->f.begin(), particle->f.end(), 0);

    // initialise potential energy and stress tensor
    en_pot_ = 0;
    stress_pot_ = 0;

    // calculate pairwise Lennard-Jones force with O(N²) algorithm
    for (size_t i = 0; i < particle->nbox - 1; ++i) {
        for (size_t j = i + 1; j < particle->nbox; ++j) {
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

            // ... and potential part of stress tensor
            stress_pot_ += fval * make_stress_tensor(rr, r);
        }
    }

    en_pot_ /= particle->nbox;
    stress_pot_ /= particle->nbox;

    // ensure that system is still in valid state
    if (std::isinf(en_pot_)) {
        throw potential_energy_divergence();
    }
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class lj_square<3, double>;
template class lj_square<2, double>;
#else
template class lj_square<3, float>;
template class lj_square<2, float>;
#endif

}}} // namespace mdsim::host::forces

#ifndef USE_HOST_SINGLE_PRECISION
template class module<mdsim::host::forces::lj_square<3, double> >;
template class module<mdsim::host::forces::lj_square<2, double> >;
#else
template class module<mdsim::host::forces::lj_square<3, float> >;
template class module<mdsim::host::forces::lj_square<2, float> >;
#endif

} // namespace halmd
