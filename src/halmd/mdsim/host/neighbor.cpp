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

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/range/iterator_range.hpp>
#include <exception>

#include <halmd/mdsim/host/neighbor.hpp>
#include <halmd/util/log.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim { namespace host
{

template <int dimension, typename float_type>
neighbor<dimension, float_type>::neighbor(particle_ptr particle, force_ptr force, box_ptr box, options const& vm)
    // dependency injection
    : particle(dynamic_pointer_cast<particle_type>(particle))
    , force(dynamic_pointer_cast<force_type>(force))
    , box(dynamic_pointer_cast<box_type>(box))
    // allocate parameters
    , rr_cut_skin_(particle->ntype, particle->ntype)
    , r0_(particle->nbox)
{
    // parse options
    r_skin_ = vm["skin"].as<float>();

    matrix_type r_cut = force->cutoff();
    matrix_type r_cut_skin(particle->ntype, particle->ntype);
    float_type r_cut_max = 0;
    for (size_t i = 0; i < particle->ntype; ++i) {
        for (size_t j = i; j < particle->ntype; ++j) {
            r_cut_skin(i, j) = r_cut(i, j) + r_skin_;
            rr_cut_skin_(i, j) = std::pow(r_cut_skin(i, j), 2);
            r_cut_max = max(r_cut_skin(i, j), r_cut_max);
        }
    }
    vector_type L = box->length();
    for (size_t i = 0; i < dimension; ++i) {
        ncell_[i] = static_cast<size_t>(L[i] / r_cut_max);
    }
    if (*min_element(ncell_.begin(), ncell_.end()) < 3) {
        throw logic_error("less than least 3 cells per dimension");
    }
    cell_.resize(ncell_);
    for (size_t i = 0; i < dimension; ++i) {
        cell_length_[i] = L[i] / ncell_[i];
    }
    r_skin_half_ = r_skin_ / 2;

    LOG("neighbor list skin: " << r_skin_);
    LOG("number of cells per dimension: " << ncell_);
    LOG("cell edge lengths: " << cell_length_);
}

/**
 * Update neighbor lists
 */
template <int dimension, typename float_type>
void neighbor<dimension, float_type>::update()
{
    cell_index i;
    for (i[0] = 0; i[0] < ncell_[0]; ++i[0]) {
        for (i[1] = 0; i[1] < ncell_[1]; ++i[1]) {
            if (dimension == 3) {
                for (i[2] = 0; i[2] < ncell_[2]; ++i[2]) {
                    update_cell_neighbors(i);
                }
            }
            else {
                update_cell_neighbors(i);
            }
        }
    }
    copy(particle->r.begin(), particle->r.end(), r0_.begin());
}

/**
 * Check if neighbour list update is needed
 */
template <int dimension, typename float_type>
bool neighbor<dimension, float_type>::check()
{
    float_type rr_max = 0;
    for (size_t i = 0; i < particle->nbox; ++i) {
        vector_type r = particle->r[i] - r0_[i];
        rr_max = max(rr_max, r * r);
    }
    return sqrt(rr_max) > r_skin_half_;
}

/**
 * Update cell lists
 */
template <int dimension, typename float_type>
void neighbor<dimension, float_type>::update_cells()
{
    // empty cell lists without memory reallocation
    BOOST_FOREACH(cell_list& c, make_iterator_range(cell_.data(), cell_.data() + cell_.num_elements())) {
        c.clear();
    }
    // add particles to cells
    for (size_t i = 0; i < particle->nbox; ++i) {
        compute_cell(particle->r[i]).push_back(i);
    }
}

/**
 * Returns cell list which a particle belongs to
 */
template <int dimension, typename float_type>
typename neighbor<dimension, float_type>::cell_list&
neighbor<dimension, float_type>::compute_cell(vector_type const& r)
{
    cell_index i;
    for (size_t j = 0; j < dimension; ++j) {
        i[j] = static_cast<size_t>(r[j] / cell_length_[j]) % ncell_[j];
    }
    return cell_(i);
}

/**
 * Update neighbor lists for a single cell
 */
template <int dimension, typename float_type>
void neighbor<dimension, float_type>::update_cell_neighbors(cell_index const& i)
{
    BOOST_FOREACH(size_t p, cell_(i)) {
        // empty neighbor list of particle
        particle->neighbor[p].clear();

        cell_index j;
        for (j[0] = -1; j[0] <= 1; ++j[0]) {
            for (j[1] = -1; j[1] <= 1; ++j[1]) {
                if (dimension == 3) {
                    for (j[2] = -1; j[2] <= 1; ++j[2]) {
                        // visit half of 26 neighbor cells due to pair potential
                        if (j[0] == 0 && j[1] == 0 && j[2] == 0) {
                            goto self;
                        }
                        // update neighbor list of particle
                        cell_index k;
                        for (int n = 0; n < dimension; ++n) {
                            k[n] = (i[n] + ncell_[n] + j[n]) % ncell_[n];
                        }
                        compute_cell_neighbors<false>(p, cell_(k));
                    }
                }
                else {
                    // visit half of 8 neighbor cells due to pair potential
                    if (j[0] == 0 && j[1] == 0) {
                        goto self;
                    }
                    // update neighbor list of particle
                    cell_index k;
                    for (int n = 0; n < dimension; ++n) {
                        k[n] = (i[n] + ncell_[n] + j[n]) % ncell_[n];
                    }
                    compute_cell_neighbors<false>(p, cell_(k));
                }
            }
        }
self:
        // visit this cell
        compute_cell_neighbors<true>(p, cell_(i));
    }
}

/**
 * Update neighbor list of particle
 */
template <int dimension, typename float_type>
template <bool same_cell>
void neighbor<dimension, float_type>::compute_cell_neighbors(size_t i, cell_list& c)
{
    vector_type L = box->length();
    vector_type L_half = static_cast<float_type>(0.5) * L;

    BOOST_FOREACH(size_t j, c) {
        // skip identical particle and particle pair permutations if same cell
        if (same_cell && particle->tag[j] <= particle->tag[i])
            continue;

        // particle distance vector
        vector_type r = particle->r[i] - particle->r[j];
        // particle types
        size_t a = particle->type[i];
        size_t b = particle->type[j];
        // enforce periodic boundary conditions
        for (size_t k = 0; k < dimension; ++k) {
            if (r[k] > L_half[k]) {
                r[k] -= L[k];
            }
            else if (r[k] < -L_half[k]) {
                r[k] += L[k];
            }
        }
        // squared particle distance
        float_type rr = r * r;

        // enforce cutoff radius with neighbor list skin
        if (rr >= rr_cut_skin_(a, b))
            continue;

        // add particle to neighbor list
        particle->neighbor[i].push_back(j);
    }
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class neighbor<3, double>;
template class neighbor<2, double>;
#else
template class neighbor<3, float>;
template class neighbor<2, float>;
#endif

}}} // namespace halmd::mdsim::host
