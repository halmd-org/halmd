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

#ifndef HALMD_MDSIM_HOST_NEIGHBOUR_HPP
#define HALMD_MDSIM_HOST_NEIGHBOUR_HPP

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/neighbour.hpp>

namespace halmd
{
namespace mdsim { namespace host
{

namespace sort
{
template <int dimension, typename float_type>
class hilbert;
}

template <int dimension, typename float_type>
class neighbour
  : public mdsim::neighbour<dimension>
{
public:
    typedef mdsim::neighbour<dimension> _Base;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef mdsim::box<dimension> box_type;

    typedef typename particle_type::neighbour_list cell_list;
    typedef boost::multi_array<cell_list, dimension> cell_lists;
    typedef fixed_vector<size_t, dimension> cell_size_type;
    typedef fixed_vector<ssize_t, dimension> cell_diff_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    static void luaopen(lua_State* L);

    neighbour(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , matrix_type const& r_cut
      , double skin
    );
    virtual void update();
    virtual bool check();

    //! returns neighbour list skin in MD units
    float_type r_skin() const
    {
        return r_skin_;
    }

protected:
    friend class sort::hilbert<dimension, float_type>; //< public interface

    void update_cells();
    void update_cell_neighbours(cell_size_type const& i);
    template <bool same_cell>
    void compute_cell_neighbours(size_t i, cell_list& c);

    /** neighbour list skin in MD units */
    float_type r_skin_;
    /** (cutoff lengths + neighbour list skin)² */
    matrix_type rr_cut_skin_;
    /** cell lists */
    cell_lists cell_;
    /** number of cells per dimension */
    cell_size_type ncell_;
    /** cell edge lengths */
    vector_type cell_length_;
    /** half neighbour list skin */
    float_type r_skin_half_;
    /** particle positions at last neighbour list update */
    std::vector<vector_type> r0_;
};

}} // namespace mdsim::host

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_NEIGHBOUR_HPP */
