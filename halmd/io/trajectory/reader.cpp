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

#include <halmd/io/trajectory/reader.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace io { namespace trajectory
{

/**
 * Assemble module options
 */
template <int dimension>
void reader<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("trajectory-file,J", po::value<string>(),
         "trajectory input file")
        ("trajectory-sample,S", po::value<ssize_t>()->default_value(-1),
         "trajectory sample for initial state")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    register_any_converter<string>();
    register_any_converter<ssize_t>();
}

template <typename T>
static void register_lua(lua_State* L, char const* class_name)
{
    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("io")
            [
                namespace_("trajectory")
                [
                    class_<T, shared_ptr<T> >(class_name)
                        .scope
                        [
                            def("options", &T::options)
                        ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        bind(&register_lua<reader<3> >, _1, "reader_3_")
    ]
    [
        bind(&register_lua<reader<2> >, _1, "reader_2_")
    ];
}

}} // namespace io::trajectory

} // namespace halmd
