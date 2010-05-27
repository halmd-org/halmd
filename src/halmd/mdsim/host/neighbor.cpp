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
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/range/iterator_range.hpp>
#include <exception>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/neighbor.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host
{

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void neighbor<dimension, float_type>::options(po::options_description& desc)
{
    desc.add_options()
        ("skin", po::value<float>()->default_value(0.5),
         "neighbour list skin")
        ;
}

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void neighbor<dimension, float_type>::resolve(po::options const& vm)
{
    module<particle_type>::required(vm);
    module<box_type>::required(vm);
    module<force_type>::required(vm);
}

template <int dimension, typename float_type>
neighbor<dimension, float_type>::neighbor(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(module<particle_type>::fetch(vm))
  , force(module<force_type>::fetch(vm))
  , box(module<box_type>::fetch(vm))
  // allocate parameters
  , rr_cut_skin_(particle->ntype, particle->ntype)
  , r0_(particle->nbox)
{
    // parse options
    r_skin_ = vm["skin"].as<float>();

    matrix_type r_cut = force->cutoff();
    matrix_type r_cut_skin(particle->ntype, particle->ntype);
    typename matrix_type::value_type r_cut_max = 0;
    for (size_t i = 0; i < particle->ntype; ++i) {
        for (size_t j = i; j < particle->ntype; ++j) {
            r_cut_skin(i, j) = r_cut(i, j) + r_skin_;
            rr_cut_skin_(i, j) = std::pow(r_cut_skin(i, j), 2);
            r_cut_max = max(r_cut_skin(i, j), r_cut_max);
        }
    }
    vector_type L = box->length();
    ncell_ = static_cast<cell_size_type>(L / r_cut_max);
    if (*min_element(ncell_.begin(), ncell_.end()) < 3) {
        throw logic_error("less than least 3 cells per dimension");
    }
    cell_.resize(ncell_);
    cell_length_ = element_div(L, static_cast<vector_type>(ncell_));
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
    // rebuild cell lists
    update_cells();
    // rebuild neighbor lists
    cell_size_type i;
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
    // make snapshot of absolute particle displacements
    copy(particle->r.begin(), particle->r.end(), r0_.begin());
}

/**
 * Check if neighbor list update is needed
 */
template <int dimension, typename float_type>
bool neighbor<dimension, float_type>::check()
{
    float_type rr_max = 0;
    for (size_t i = 0; i < particle->nbox; ++i) {
        vector_type r = particle->r[i] - r0_[i];
        rr_max = max(rr_max, inner_prod(r, r));
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
    for_each(cell_.data(), cell_.data() + cell_.num_elements(), bind(&cell_list::clear, _1));
    // add particles to cells
    for (size_t i = 0; i < particle->nbox; ++i) {
        vector_type const& r = particle->r[i];
        cell_size_type index = element_mod(static_cast<cell_size_type>(element_div(r, cell_length_)), ncell_);
        cell_(index).push_back(i);
    }
}

/**
 * Update neighbor lists for a single cell
 */
template <int dimension, typename float_type>
void neighbor<dimension, float_type>::update_cell_neighbors(cell_size_type const& i)
{
    BOOST_FOREACH(size_t p, cell_(i)) {
        // empty neighbor list of particle
        particle->neighbor[p].clear();

        cell_diff_type j;
        for (j[0] = -1; j[0] <= 1; ++j[0]) {
            for (j[1] = -1; j[1] <= 1; ++j[1]) {
                if (dimension == 3) {
                    for (j[2] = -1; j[2] <= 1; ++j[2]) {
                        // visit half of 26 neighbor cells due to pair potential
                        if (j[0] == 0 && j[1] == 0 && j[2] == 0) {
                            goto self;
                        }
                        // update neighbor list of particle
                        cell_size_type k = element_mod(static_cast<cell_size_type>(static_cast<cell_diff_type>(i + ncell_) + j), ncell_);
                        compute_cell_neighbors<false>(p, cell_(k));
                    }
                }
                else {
                    // visit half of 8 neighbor cells due to pair potential
                    if (j[0] == 0 && j[1] == 0) {
                        goto self;
                    }
                    // update neighbor list of particle
                    cell_size_type k = element_mod(static_cast<cell_size_type>(static_cast<cell_diff_type>(i + ncell_) + j), ncell_);
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
    BOOST_FOREACH(size_t j, c) {
        // skip identical particle and particle pair permutations if same cell
        if (same_cell && particle->tag[j] <= particle->tag[i])
            continue;

        // particle distance vector
        vector_type r = particle->r[i] - particle->r[j];
        box->reduce_periodic(r);
        // particle types
        size_t a = particle->type[i];
        size_t b = particle->type[j];
        // squared particle distance
        float_type rr = inner_prod(r, r);

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

}} // namespace mdsim::host

#ifndef USE_HOST_SINGLE_PRECISION
template class module<mdsim::host::neighbor<3, double> >;
template class module<mdsim::host::neighbor<2, double> >;
#else
template class module<mdsim::host::neighbor<3, float> >;
template class module<mdsim::host::neighbor<2, float> >;
#endif

} // namespace halmd
