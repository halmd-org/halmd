/*
 * Copyright © 2008-2012  Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_BINNING_HPP
#define HALMD_MDSIM_HOST_BINNING_HPP

#include <halmd/algorithm/multi_range.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/profiler.hpp>

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <lua.hpp>

#include <algorithm>
#include <memory>
#include <vector>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class binning
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef mdsim::box<dimension> box_type;
    typedef logger logger_type;

    typedef std::vector<unsigned int> cell_list;
    typedef boost::multi_array<cell_list, dimension> array_type;
    typedef fixed_vector<size_t, dimension> cell_size_type;
    typedef fixed_vector<ssize_t, dimension> cell_diff_type;

    static void luaopen(lua_State* L);

    binning(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<box_type const> box
      , matrix_type const& r_cut
      , float_type skin
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    //! returns neighbour list skin in MD units
    float_type r_skin() const
    {
        return r_skin_;
    }

    //! cell edge length
    vector_type const& cell_length() const
    {
        return cell_length_;
    }

    //! number of cells per dimension
    cell_size_type ncell() const
    {
        return ncell_;
    }

    //! get cell lists
    cache<array_type> const& cell();

private:
    typedef typename particle_type::size_type size_type;
    typedef typename particle_type::position_array_type position_array_type;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type update;
    };

    void update();

    //! system state
    std::shared_ptr<particle_type const> particle_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;
    /** neighbour list skin in MD units */
    float_type r_skin_;
    /** cell lists */
    cache<array_type> cell_;
    /** cache observer for cell list update */
    cache<> cell_cache_;
    /** number of cells per dimension */
    cell_size_type ncell_;
    /** cell edge lengths */
    vector_type cell_length_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

/**
 * Copy cells to multi-range output iterator.
 *
 * @param output multi-range output iterator
 *
 * A multi-range iterator is a functor that accepts a multi-dimensional
 * index of array type, and returns an output iterator for the given
 * index. The particle indices in the cell of the given index are
 * then copied to the returned output iterator.
 */
template <typename binning_type, typename output_iterator>
inline void
get_cell(binning_type& binning, output_iterator output)
{
    typedef typename binning_type::array_type array_type;
    typedef typename binning_type::cell_size_type cell_size_type;
    array_type const& cell = read_cache(binning.cell());
    multi_range_for_each(
        cell_size_type(0)
      , binning.ncell()
      , [&](cell_size_type const& index) {
            std::copy(cell(index).begin(), cell(index).end(), output(index));
        }
    );
}

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_BINNING_HPP */
