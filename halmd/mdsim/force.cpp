/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/force.hpp>
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
void force<dimension>::options(options_definition& options)
{
    options.add_options()
        ("force", po::value<string>()->default_value("lennard_jones"),
         "specify force module")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    register_any_converter<string>();
}

template <int dimension>
void force<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    string class_name("force_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                class_<force, shared_ptr<force> >(class_name.c_str())
                    .def("compute", &force::compute)
                    .scope
                    [
                        def("options", &force::options)
                    ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &force<3>::luaopen
    ]
    [
        &force<2>::luaopen
    ];
}

// explicit instantiation
template class force<3>;
template class force<2>;

} // namespace mdsim

} // namespace halmd
