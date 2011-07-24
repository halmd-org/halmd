/*
 * Copyright © 2011  Peter Colberg
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

#ifndef HALMD_IO_READERS_H5MD_APPEND_HPP
#define HALMD_IO_READERS_H5MD_APPEND_HPP

#include <boost/function.hpp>
#include <lua.hpp>

#include <h5xx/h5xx.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace io {
namespace readers {
namespace h5md {

/**
 * H5MD dataset reader (append mode)
 *
 * This module implements collective reading from one or multiple H5MD
 * datasets, where each dataset is a time series. Upon initialisation,
 * the reader is assigned a collective H5MD group. A dataset within this
 * group is created by connecting a data slot to the on_read signal.
 * All datasets share common step and time datasets, which are linked
 * into each dataset group upon connection.
 *
 * The reader provides a common write slot, which may be connected to
 * the sampler to write to the datasets at a fixed interval. Further
 * signals on_prepend_read and on_append_read are provided to call
 * arbitrary slots before and after reading.
 */
class append
{
private:
    typedef mdsim::clock clock_type;
    typedef signal<void ()> signal_type;

public:
    typedef clock_type::step_type step_type;
    typedef clock_type::time_type time_type;
    typedef signal_type::slot_function_type slot_function_type;
    typedef H5::Group subgroup_type;

    /** open reader group and create time and step datasets */
    append(
        H5::Group const& root
      , std::vector<std::string> const& location
    );
    /** connect data slot for reading dataset */
    template <typename T>
    void on_read(
        subgroup_type& group
      , boost::function<T ()> const& slot
      , std::vector<std::string> const& location
    );
    /** connect slot called before reading */
    void on_prepend_read(slot_function_type const& slot);
    /** connect slot called after reading */
    void on_append_read(slot_function_type const& slot);
    /** read at given step offset */
    void read_step(step_type step);
    /** read at given time offset */
    void read_time(time_type time);
    /** Lua bindings */
    static void luaopen(lua_State* L);

private:
    typedef boost::function<hsize_t (H5::Group const& group)> index_function_type;

    template <typename T>
    static void read_dataset(
        H5::DataSet dataset
      , boost::function<T ()> const& slot
      , index_function_type const& index
      , H5::Group const& group
    );
    static hsize_t read_step_index(
        step_type step
      , H5::Group const& group
    );
    static hsize_t read_time_index(
        time_type time
      , H5::Group const& group
    );

    /** reader group */
    H5::Group group_;
    /** signal emitted for reading datasets */
    signal<void (index_function_type const&)> on_read_;
    /** signal emitted before reading datasets */
    signal_type on_prepend_read_;
    /** signal emitted before after datasets */
    signal_type on_append_read_;
};

} // namespace h5md
} // namespace readers
} // namespace io
} // namespace halmd

#endif /* ! HALMD_IO_READERS_H5MD_APPEND_HPP */
