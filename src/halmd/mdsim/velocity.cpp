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
#include <halmd/mdsim/velocity.hpp>
#include <halmd/utility/lua.hpp>

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
void velocity<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("velocity",
         po::value<string>()->default_value("boltzmann"),
         "initial particle velocities module")
        ;
}

template <typename T>
static void register_lua(char const* class_name)
{
    using namespace luabind;
    lua_registry::get()->push_back
    ((
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
    ));
}

static __attribute__((constructor)) void register_lua()
{
    register_lua<velocity<3> >("velocity_3_");
    register_lua<velocity<2> >("velocity_2_");
}

template class velocity<3>;
template class velocity<2>;

} // namespace mdsim

} // namespace halmd
