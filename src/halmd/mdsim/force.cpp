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
void force<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("force", po::value<string>()->default_value("lj"),
         "specify force module")
        ;
}

/**
 * Resolve module dependencies
 */
template <int dimension>
void force<dimension>::depends()
{
    modules::depends<_Self, particle_type>::required();
}

template <int dimension>
force<dimension>::force(modules::factory& factory, po::options const& vm)
  // dependency injection
  : particle(modules::fetch<particle_type>(factory, vm))
{
}

template <typename T>
static luabind::scope register_lua(char const* class_name)
{
    using namespace luabind;
    return
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
        ];
}

static lua_registry::iterator dummy = (
    lua_registry::get()->push_back( register_lua<force<3> >("force_3_") )
  , lua_registry::get()->push_back( register_lua<force<2> >("force_2_") )
  , lua_registry::get()->end()
);

// explicit instantiation
template class force<3>;
template class force<2>;

} // namespace mdsim

} // namespace halmd
