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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/position.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Assemble module options
 */
template <int dimension>
void position<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("position",
         po::value<string>()->default_value("lattice"),
         "initial particle positions module")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    using namespace lua_wrapper;
    register_any_converter<string>();
}


template <typename T>
static void register_lua(char const* class_name)
{
    using namespace luabind;
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                class_<T, shared_ptr<T> >(class_name)
                    .scope
                    [
                        def("options", &T::options)
                    ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    register_lua<position<3> >("position_3_");
    register_lua<position<2> >("position_2_");
}

template class position<3>;
template class position<2>;

} // namespace mdsim

} // namespace halmd
