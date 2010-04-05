/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <halmd/mdsim/backend/exception.hpp>
#include <halmd/mdsim/host/forces/lj.hpp>
#include <halmd/util/log.hpp>

using namespace boost;
using namespace boost::numeric::ublas;

namespace halmd { namespace mdsim { namespace host { namespace forces
{

/**
 * Initialize Lennard-Jones potential parameters
 */
template <int dimension, typename float_type>
lj<dimension, float_type>::lj(particle_ptr const& particle, box_ptr const& box, options const& vm)
    // dependency injection
    : particle(static_pointer_cast<particle_type>(particle))
    , box(static_pointer_cast<box_type>(box))
    // allocate potential parameters
    , epsilon_(scalar_matrix<float_type>(particle->ntype, particle->ntype, 1))
    , sigma_(scalar_matrix<float_type>(particle->ntype, particle->ntype, 1))
    , r_cut_(particle->ntype, particle->ntype)
    , r_cut_sigma_(particle->ntype, particle->ntype)
    , r_cut_sigma2_(particle->ntype, particle->ntype)
    , sigma2_(particle->ntype, particle->ntype)
    , en_cut_(particle->ntype, particle->ntype)
{
    // parse options
    if (!vm["binary"].empty() && vm["particles"].defaulted()) {
        boost::array<float, 3> epsilon(vm["epsilon"].as<boost::array<float, 3> >());
        std::copy(epsilon.begin(), epsilon.end(), epsilon_.data().begin());
        boost::array<float, 3> sigma(vm["sigma"].as<boost::array<float, 3> >());
        std::copy(sigma.begin(), sigma.end(), sigma_.data().begin());
    }
    try {
        boost::array<float, 3> r_cut(vm["cutoff"].as<boost::array<float, 3> >());
        std::copy(r_cut.begin(), r_cut.end(), r_cut_.data().begin());
    }
    catch (boost::bad_any_cast const&) {
        // backwards compatibility
        std::fill(r_cut_.data().begin(), r_cut_.data().end(), vm["cutoff"].as<float>());
    }

    // precalculate derived parameters
    for (size_t i = 0; i < particle->ntype; ++i) {
        for (size_t j = i; j < particle->ntype; ++j) {
            r_cut_sigma_(i, j) = r_cut_(i, j) * sigma_(i, j);
            r_cut_sigma2_(i, j) = std::pow(r_cut_sigma_(i, j), 2);
            sigma2_(i, j) = std::pow(sigma_(i, j), 2);
            // energy shift due to truncation at cutoff length
            float_type rri_cut = std::pow(r_cut_(i, j), -2);
            float_type r6i_cut = rri_cut * rri_cut * rri_cut;
            en_cut_(i, j) = 4 * epsilon_(i, j) * r6i_cut * (r6i_cut - 1);
        }
    }

    LOG("potential well depths: ε = " << epsilon_);
    LOG("potential pair separation: σ = " << sigma_);
    LOG("potential cutoff length: r = " << r_cut_);
    LOG("potential cutoff energy: U = " << en_cut_);
}

/**
 * Compute Lennard-Jones forces
 */
template <int dimension, typename float_type>
void lj<dimension, float_type>::compute()
{
    // initialize particle forces to zero
    std::fill(particle->f.begin(), particle->f.end(), 0);

    // potential energy
    en_pot_ = 0;
    // virial equation sum
    std::fill(virial_.begin(), virial_.end(), 0);
    // half periodic box edge lengths for nearest mirror-image particle
    vector_type box_half = static_cast<float_type>(0.5) * box->length();

    for (size_t i = 0; i < particle->nbox; ++i) {
        // calculate pairwise Lennard-Jones force with neighbor particles
        BOOST_FOREACH (size_t j, particle->neighbor[i]) {
            // particle distance vector
            vector_type r = particle->r[i] - particle->r[j];
            // particle types
            size_t a = particle->type[i];
            size_t b = particle->type[j];
            // enforce periodic boundary conditions
            // FIXME use ghost particles to implement minimum image convention
            for (size_t k = 0; k < dimension; ++k) {
                if (r[k] > box_half[k]) {
                    r[k] -= box->length()[k];
                }
                else if (r[k] < -box_half[k]) {
                    r[k] += box->length()[k];
                }
            }
            // squared particle distance
            float_type rr = r * r;

            // truncate potential at cutoff length
            if (rr >= r_cut_sigma2_(a, b))
                continue;

            // compute Lennard-Jones force in reduced units
            float_type sigma2 = sigma2_(a, b);
            float_type rri = sigma2 / rr;
            float_type r6i = rri * rri * rri;
            float_type epsilon = epsilon_(a, b);
            float_type fval = 48 * rri * r6i * (r6i - 0.5) * (epsilon / sigma2);
            float_type en_pot = 4 * epsilon * r6i * (r6i - 1) - en_cut_(a, b);

            // add force contribution to both particles
            particle->f[i] += r * fval;
            particle->f[j] -= r * fval;

            // add contribution to potential energy
            en_pot_ += en_pot;

            // add contribution to virial
            float_type virial = 0.5 * rr * fval;
            virial_[a][0] += virial;
            virial_[b][0] += virial;

            // compute off-diagonal virial stress tensor elements
            if (dimension == 3) {
                virial = 0.5 * r[1] * r[2] * fval;
                virial_[a][1] += virial;
                virial_[b][1] += virial;

                virial = 0.5 * r[2] * r[0] * fval;
                virial_[a][2] += virial;
                virial_[b][2] += virial;

                virial = 0.5 * r[0] * r[1] * fval;
                virial_[a][3] += virial;
                virial_[b][3] += virial;
            }
            else {
                virial = 0.5 * r[0] * r[1] * fval;
                virial_[a][1] += virial;
                virial_[b][1] += virial;
            }
        }
    }

    en_pot_ /= particle->nbox;

    // ensure that system is still in valid state
    if (std::isinf(en_pot_)) {
        throw potential_energy_divergence();
    }
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class lj<3, double>;
template class lj<2, double>;
#else
template class lj<3, float>;
template class lj<2, float>;
#endif

}}}} // namespace halmd::mdsim::host::forces
