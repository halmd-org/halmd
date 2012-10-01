/*
 * Copyright © 2011  Felix Höfling
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

#include <boost/make_shared.hpp>
#include <limits>

#include <halmd/observables/host/samples/binned_phase_space.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {
namespace samples {

template <int dimension, typename float_type>
binned_phase_space<dimension, float_type>::binned_phase_space(
    shared_ptr<data_sample_type const> data_sample
  , unsigned int species
  , fixed_vector<unsigned int, dimension> const& nbin
)
  // initialise public attributes
  : cell_array(nbin)
  , nbin(nbin)
  , step(numeric_limits<step_type>::max())
  // private module dependencies
  , data_sample_(data_sample)
  , species_(species)
{}

template <int dimension, typename float_type>
static int wrap_dimension(binned_phase_space<dimension, float_type> const&)
{
    return dimension;
}

template <int dimension, typename float_type>
void binned_phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("binned_phase_space_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>() + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                namespace_("samples")
                [
                    class_<binned_phase_space>(class_name.c_str())
                        .property("dimension", &wrap_dimension<dimension, float_type>)
                ]
            ]

          , namespace_("samples")
            [
                def("binned_phase_space", &make_shared<binned_phase_space
                  , shared_ptr<data_sample_type const>
                  , unsigned int
                  , fixed_vector<unsigned int, dimension>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_samples_binned_phase_space(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    binned_phase_space<3, double>::luaopen(L);
    binned_phase_space<2, double>::luaopen(L);
#endif
    binned_phase_space<3, float>::luaopen(L);
    binned_phase_space<2, float>::luaopen(L);
    return 0;
}

#ifndef USE_HOST_SINGLE_PRECISION
template class binned_phase_space<3, double>;
template class binned_phase_space<2, double>;
#endif
template class binned_phase_space<3, float>;
template class binned_phase_space<2, float>;

} // namespace observables
} // namespace host
} // namespace samples
} // namespace halmd
