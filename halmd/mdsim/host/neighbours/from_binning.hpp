/*
 * Copyright © 2015       Nicolas Höft
 * Copyright © 2015       Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_NEIGHBOURS_FROM_BINNING_HPP
#define HALMD_MDSIM_HOST_NEIGHBOURS_FROM_BINNING_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/binning.hpp>
#include <halmd/mdsim/host/max_displacement.hpp>
#include <halmd/mdsim/host/neighbour.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/profiler.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <lua.hpp>

#include <memory>
#include <vector>

namespace halmd {
namespace mdsim {
namespace host {
namespace neighbours {

template <int dimension, typename float_type>
class from_binning
  : public mdsim::host::neighbour
{
private:
    typedef mdsim::host::neighbour _Base;

public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef mdsim::box<dimension> box_type;
    typedef host::binning<dimension, float_type> binning_type;
    typedef typename _Base::neighbour_list neighbour_list;
    typedef max_displacement<dimension, float_type> displacement_type;

    typedef _Base::array_type array_type;

    static void luaopen(lua_State* L);

    from_binning(
        std::pair<std::shared_ptr<particle_type const>, std::shared_ptr<particle_type const>> particle
      , std::pair<std::shared_ptr<binning_type>, std::shared_ptr<binning_type>> binning
      , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>> displacement
      , std::shared_ptr<box_type const> box
      , matrix_type const& r_cut
      , double skin
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    connection on_prepend_update(std::function<void ()> const& slot)
    {
        return on_prepend_update_.connect(slot);
    }

    connection on_append_update(std::function<void ()> const& slot)
    {
        return on_append_update_.connect(slot);
    }

    //! returns neighbour list skin in MD units
    float_type r_skin() const
    {
        return r_skin_;
    }

    //! returns neighbour lists
    virtual cache<array_type> const& lists();

    //! returns true if the binning modules are compatible with the neighbour list module
    static bool is_binning_compatible(
        std::shared_ptr<binning_type const> binning1
      , std::shared_ptr<binning_type const> binning2
    );

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::reverse_tag_array_type reverse_tag_array_type;
    typedef typename particle_type::species_array_type species_array_type;
    typedef typename particle_type::species_type species_type;
    typedef typename particle_type::size_type size_type;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type update;
    };

    typedef typename binning_type::cell_size_type cell_size_type;
    typedef typename binning_type::cell_diff_type cell_diff_type;
    typedef typename binning_type::cell_list cell_list;
    typedef typename binning_type::array_type cell_array_type;

    std::shared_ptr<particle_type const> particle1_;
    std::shared_ptr<particle_type const> particle2_;
    std::shared_ptr<binning_type> binning1_;
    std::shared_ptr<binning_type> binning2_;
    std::shared_ptr<displacement_type> displacement1_;
    std::shared_ptr<displacement_type> displacement2_;
    std::shared_ptr<box_type const> box_;
    std::shared_ptr<logger> logger_;

    void update();
    void update_cell_neighbours(cell_size_type const& i);
    template <bool same_cell>
    void compute_cell_neighbours(size_t i, cell_list const& c);

    /** neighbour lists */
    cache<array_type> neighbour_;
    /** cache observer for neighbour list update */
    std::tuple<cache<>, cache<>> neighbour_cache_;
    /** neighbour list skin in MD units */
    float_type r_skin_;
    /** (cutoff lengths + neighbour list skin)² */
    matrix_type rr_cut_skin_;
    /** signal emitted before neighbour list update */
    signal<void ()> on_prepend_update_;
    /** signal emitted after neighbour list update */
    signal<void ()> on_append_update_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace neighbours
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_NEIGHBOURS_FROM_BINNING_HPP */
