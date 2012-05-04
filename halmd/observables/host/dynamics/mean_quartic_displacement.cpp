/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/host/dynamics/mean_quartic_displacement.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {
namespace dynamics {

template <int dimension, typename float_type>
mean_quartic_displacement<dimension, float_type>::mean_quartic_displacement(
    size_t type
)
  // member initialisation
  : type_(type)
{
    LOG("initialise mean-quartic displacement of " << string(1, 'A' + type) << " particles");
}

/**
 * Compute mean-quartic displacement of two position sample vectors.
 *
 * @param first particles positions of one species at time t1
 * @param second particles positions of one species at time t2
 * @returns accumulated mean-quartic displacement
 */
template <int dimension, typename float_type>
typename mean_quartic_displacement<dimension, float_type>::accumulator_type
mean_quartic_displacement<dimension, float_type>::compute(
    sample_type const& first
  , sample_type const& second
)
{
    accumulator_type acc;
    typename sample_type::sample_vector::const_iterator r1, r2, end = first.r[type_]->end();
    for (r1 = first.r[type_]->begin(), r2 = second.r[type_]->begin(); r1 != end; ++r1, ++r2) {
        // accumulate quartic displacement
        acc(correlate_function_type()(*r1, *r2));
    }
    return acc;
}

template <typename tcf_type>
static boost::shared_ptr<tcf_type>
wrap_tcf(size_t type, typename tcf_type::sample_type const&)
{
    return boost::make_shared<tcf_type>(type);
}

template <int dimension, typename float_type>
void mean_quartic_displacement<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("mean_quartic_displacement_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<mean_quartic_displacement>(class_name.c_str())

              , def("mean_quartic_displacement", &wrap_tcf<mean_quartic_displacement>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_dynamics_mean_quartic_displacement(lua_State* L)
{
    mean_quartic_displacement<3, double>::luaopen(L);
    mean_quartic_displacement<2, double>::luaopen(L);
    mean_quartic_displacement<3, float>::luaopen(L);
    mean_quartic_displacement<2, float>::luaopen(L);
    observables::dynamics::correlation<mean_quartic_displacement<3, double> >::luaopen(L);
    observables::dynamics::correlation<mean_quartic_displacement<2, double> >::luaopen(L);
    observables::dynamics::correlation<mean_quartic_displacement<3, float> >::luaopen(L);
    observables::dynamics::correlation<mean_quartic_displacement<2, float> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class mean_quartic_displacement<3, double>;
template class mean_quartic_displacement<2, double>;
template class mean_quartic_displacement<3, float>;
template class mean_quartic_displacement<2, float>;

} // namespace dynamics
} // namespace host

namespace dynamics
{

// explicit instantiation
template class correlation<host::dynamics::mean_quartic_displacement<3, double> >;
template class correlation<host::dynamics::mean_quartic_displacement<2, double> >;
template class correlation<host::dynamics::mean_quartic_displacement<3, float> >;
template class correlation<host::dynamics::mean_quartic_displacement<2, float> >;

} // namespace dynamics
} // namespace observables
} // namespace halmd
