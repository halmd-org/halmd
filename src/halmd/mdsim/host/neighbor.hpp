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

#ifndef HALMD_MDSIM_HOST_NEIGHBOR_HPP
#define HALMD_MDSIM_HOST_NEIGHBOR_HPP

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/mdsim/neighbor.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim { namespace host
{

template <int dimension, typename float_type>
class neighbor : public mdsim::neighbor<dimension, float_type>
{
public:
    typedef mdsim::neighbor<dimension, float_type> _Base;
    typedef vector<float_type, dimension> vector_type;

    typedef host::particle<dimension, float_type> particle_type;
    typedef mdsim::force<dimension, float_type> force_type;
    typedef mdsim::box<dimension, float_type> box_type;

    typedef typename particle_type::neighbor_list cell_list;
    typedef boost::multi_array<cell_list, dimension> cell_lists;
    typedef vector<size_t, dimension> cell_index;
    typedef vector<ssize_t, dimension> cell_diff;

    typedef typename force_type::matrix_type matrix_type;

    neighbor(options const& vm);
    virtual ~neighbor() {}
    void update();
    bool check();

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<force_type> force;
    boost::shared_ptr<box_type> box;

protected:
    void update_cells();
    void update_cell_neighbors(cell_index const& i);
    template <bool same_cell>
    void compute_cell_neighbors(size_t i, cell_list& c);

    /** neighbor list skin in MD units */
    float_type r_skin_;
    /** (cutoff lengths + neighbor list skin)² */
    matrix_type rr_cut_skin_;
    /** cell lists */
    cell_lists cell_;
    /** number of cells per dimension */
    cell_index ncell_;
    /** cell edge lengths */
    vector_type cell_length_;
    /** half neighbor list skin */
    float_type r_skin_half_;
    /** particle positions at last neighbor list update */
    std::vector<vector_type> r0_;
};

}}} // namespace halmd::mdsim::host

#endif /* ! HALMD_MDSIM_HOST_NEIGHBOR_HPP */
