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
#include <halmd/mdsim/integrator.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

template <int dimension>
double const integrator<dimension>::default_timestep = 0.001;

/**
 * Assemble module options
 */
template <int dimension>
void integrator<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("integrator", po::value<string>()->default_value("verlet"),
         "specify integration module")
        ("timestep,h", po::value<double>()->default_value(default_timestep),
         "integration timestep")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    using namespace lua_wrapper;
    register_any_converter<string>();
    register_any_converter<double>();
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
                    .property("timestep", (double (T::*)() const) &T::timestep)
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
    register_lua<integrator<3> >("integrator_3_");
    register_lua<integrator<2> >("integrator_2_");
}

// explicit instantiation
template class integrator<3>;
template class integrator<2>;

} // namespace mdsim

} // namespace halmd
