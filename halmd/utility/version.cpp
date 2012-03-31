/*
 * Copyright © 2012  Peter Colberg
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

#include <boost/algorithm/string/join.hpp>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/utility/hostname.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

namespace halmd {

/**
 * Log HALMD version, build flags, command line and host name.
 */
static void prologue(vector<string> const& arg)
{
    LOG(PROJECT_NAME " (" PROGRAM_DESC ") " PROGRAM_VERSION);
    LOG("variant: " << PROGRAM_VARIANT);
#ifndef NDEBUG
    LOG_WARNING("built with enabled debugging");
#endif
    LOG("command line: " << join(arg, " "));
    LOG("host name: " << host_name());
}

HALMD_LUA_API int luaopen_libhalmd_utility_version(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("utility")
        [
            namespace_("version")
            [
                def("prologue", &prologue)
            ]
        ]
    ];
    return 0;
}

} // namespace halmd
