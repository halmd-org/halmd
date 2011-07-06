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

#include <ctime>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/writers/h5md.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

namespace halmd {
namespace io {
namespace writers {

/**
 * This overwrites any existing file at the given path.
 */
h5md::h5md(string const& path)
{
    file_ = H5::H5File(path, H5F_ACC_TRUNC);

    H5::Group attr = file_.createGroup("h5md");
    h5xx::write_attribute(attr, "creation_time", time(NULL));
    h5xx::write_attribute(attr, "creator", PROGRAM_NAME);
    h5xx::write_attribute(attr, "creator_version", PROGRAM_VERSION);
    h5xx::write_attribute(attr, "version", h5md::version());
}

void h5md::flush()
{
    file_.flush(H5F_SCOPE_GLOBAL);
}

void h5md::close()
{
    file_.close();
}

h5md::version_type h5md::version()
{
    version_type version = {{ 0, 0 }};
    return version;
}

static signal<void ()>::slot_function_type
wrap_flush(shared_ptr<h5md> instance)
{
    return bind(&h5md::flush, instance);
}

void h5md::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("writers")
            [
                class_<h5md>("h5md")
                    .def(constructor<string const&>())
                    // wrap as slot to support periodic flushing of file
                    .property("flush", &wrap_flush)
                    .def("close", &h5md::close)
                    .property("file", &h5md::file)
                    .property("path", &h5md::path)
                    .scope
                    [
                        def("version", &h5md::version)
                    ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_writers_h5md(lua_State* L)
{
    h5md::luaopen(L);
    return 0;
}

} // namespace writers
} // namespace io
} // namespace halmd
