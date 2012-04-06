/*
 * Copyright © 2011-2012  Peter Colberg
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

#ifndef HALMD_IO_WRITERS_H5MD_TRUNCATE_HPP
#define HALMD_IO_WRITERS_H5MD_TRUNCATE_HPP

#include <boost/function.hpp>
#include <boost/multi_array.hpp>
#include <lua.hpp>

#include <h5xx/h5xx.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace io {
namespace writers {
namespace h5md {

/**
 * H5MD dataset writer (truncate mode)
 *
 * This module implements collective writing to one or multiple H5MD
 * datasets. Upon initialisation, the writer is assigned a collective
 * H5MD group. A dataset within this group is created by connecting a
 * data slot to the on_write signal.
 *
 * The writer provides a common write slot, which may be connected to
 * the sampler to write to the datasets at a fixed interval. Further
 * signals on_prepend_write and on_append_write are provided to call
 * arbitrary slots before and after writing.
 */
class truncate
{
private:
    typedef signal<void ()> signal_type;

public:
    typedef signal_type::slot_function_type slot_function_type;
    /**
     * For the truncate reader/writer, a subgroup is defined as the dataset
     * which contains the data to be read or written.
     * For the append reader/writer, a subgroup is defined as the group
     * comprising the samples, time and step datasets, where samples
     * contains the data to be read or written. Additional attributes
     * should always be attached to the subgroup, never the samples
     * dataset.
     * To give both writers the same API for convenient use in template
     * functions in C++ unit tests, we define a subgroup type.
     */
    typedef H5::DataSet subgroup_type;

    /** open writer group */
    truncate(
        H5::Group const& root
      , std::vector<std::string> const& location
    );
    /** connect data slot for writing */
    template <typename T>
    connection on_write(
        subgroup_type& dataset
      , boost::function<T ()> const& slot
      , std::vector<std::string> const& location
    );
    /** connect slot called before writing */
    connection on_prepend_write(slot_function_type const& slot);
    /** connect slot called after writing */
    connection on_append_write(slot_function_type const& slot);
    /** write datasets */
    void write();
    /** Lua bindings */
    static void luaopen(lua_State* L);

    /**
     * returns writer group
     */
    H5::Group const& group() const
    {
        return group_;
    }

private:
    /** writer group */
    H5::Group group_;
    /** signal emitted for writing datasets */
    signal_type on_write_;
    /** signal emitted before writing datasets */
    signal_type on_prepend_write_;
    /** signal emitted before after datasets */
    signal_type on_append_write_;
};

} // namespace h5md
} // namespace writers
} // namespace io
} // namespace halmd

#endif /* ! HALMD_IO_WRITERS_H5MD_TRUNCATE_HPP */
