/*
 * Copyright © 2011  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_IO_READERS_H5MD_TRUNCATE_HPP
#define HALMD_IO_READERS_H5MD_TRUNCATE_HPP

#include <functional>
#include <lua.hpp>

#include <h5xx/h5xx.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace io {
namespace readers {
namespace h5md {

/**
 * H5MD dataset reader (truncate mode)
 *
 * This module implements collective reading from one or multiple H5MD
 * datasets. Upon initialisation, the reader is assigned a collective
 * H5MD group. A dataset within this group is created by connecting a
 * data slot to the on_read signal.
 *
 * The reader provides a common read slot, which may be connected to
 * another module to read the datasets. Further signals on_prepend_read
 * and on_append_read are provided to call arbitrary slots before and
 * after reading.
 */
class truncate
{
public:
    typedef signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;
    /**
     * For the truncate reader/writer, a subgroup is defined as the dataset
     * which contains the data to be read or written.
     * For the append reader/writer, a subgroup is defined as the group
     * comprising the samples, time and step datasets, where samples
     * contains the data to be read or written. Additional attributes
     * should always be attached to the subgroup, never the samples
     * dataset.
     * To give both readers the same API for convenient use in template
     * functions in C++ unit tests, we define a subgroup type.
     */
    typedef H5::DataSet subgroup_type;

    /** open reader group */
    truncate(
        H5::Group const& root
      , std::vector<std::string> const& location
    );
    /** connect data slot for reading */
    template <typename T>
    connection on_read(
        subgroup_type& dataset
      , std::function<T ()> const& slot
      , std::vector<std::string> const& location
    );
    /** connect slot called before reading */
    connection on_prepend_read(slot_function_type const& slot);
    /** connect slot called after reading */
    connection on_append_read(slot_function_type const& slot);
    /** read datasets */
    void read();
    /** Lua bindings */
    static void luaopen(lua_State* L);

    /**
     * returns reader group
     */
    H5::Group const& group() const
    {
        return group_;
    }

private:
    template <typename T>
    static void read_dataset(
        H5::DataSet dataset
      , std::function<T ()> const& slot
    );

    /** reader group */
    H5::Group group_;
    /** signal emitted for reading datasets */
    signal_type on_read_;
    /** signal emitted before reading datasets */
    signal_type on_prepend_read_;
    /** signal emitted before after datasets */
    signal_type on_append_read_;
};

} // namespace h5md
} // namespace readers
} // namespace io
} // namespace halmd

#endif /* ! HALMD_IO_READERS_H5MD_TRUNCATE_HPP */
