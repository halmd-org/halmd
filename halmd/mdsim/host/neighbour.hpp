/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/binning.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/neighbour.hpp>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class neighbour
  : public mdsim::neighbour
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef mdsim::box<dimension> box_type;
    typedef host::binning<dimension, float_type> binning_type;
    typedef std::vector<unsigned int> neighbour_list;
    typedef typename neighbour::signal_type signal_type;
    typedef typename neighbour::slot_function_type slot_function_type;
    typedef typename neighbour::connection_type connection_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    neighbour(
        boost::shared_ptr<particle_type const> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<binning_type const> binning
      , matrix_type const& r_cut
      , double skin
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    virtual void update();

    virtual connection_type on_update(slot_function_type const& slot)
    {
        return on_update_.connect(slot);
    }

    //! returns neighbour list skin in MD units
    float_type r_skin() const
    {
        return r_skin_;
    }

    //! returns neighbour lists
    std::vector<neighbour_list> const& lists() const
    {
        return neighbour_;
    }

private:
    typedef typename binning_type::cell_size_type cell_size_type;
    typedef typename binning_type::cell_diff_type cell_diff_type;
    typedef typename binning_type::cell_list cell_list;
    typedef typename binning_type::cell_lists cell_lists;

    boost::shared_ptr<particle_type const> particle_;
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<binning_type const> binning_;
    boost::shared_ptr<logger_type> logger_;

    void update_cell_neighbours(cell_size_type const& i);
    template <bool same_cell>
    void compute_cell_neighbours(size_t i, cell_list const& c);

    /** neighbour lists */
    std::vector<neighbour_list> neighbour_;
    /** neighbour list skin in MD units */
    float_type r_skin_;
    /** (cutoff lengths + neighbour list skin)² */
    matrix_type rr_cut_skin_;
    /** signal emitted before neighbour list update */
    signal_type on_update_;
};

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_NEIGHBOUR_HPP */
