/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/io/logger.hpp>
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/host/dynamics/mean_quartic_displacement.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace dynamics {

/**
 * Compute mean-quartic displacement of two position sample vectors.
 *
 * @param first particles positions of one species at time t1
 * @param second particles positions of one species at time t2
 * @returns accumulated mean-quartic displacement
 */
template <int dimension, typename float_type>
void mean_quartic_displacement<dimension, float_type>::operator() (
    sample_type const& first
  , sample_type const& second
  , accumulator<result_type>& result
)
{
    accumulator<result_type> acc;
    typename sample_type::array_type::const_iterator r1, r2, end = first.data().end();
    for (r1 = first.data().begin(), r2 = second.data().begin(); r1 != end; ++r1, ++r2) {
        // accumulate quartic displacement
        acc(correlate_function_type()(*r1, *r2));
    }
    result(acc);
}

template <typename tcf_type>
static std::shared_ptr<tcf_type>
select_tcf_by_acquire(std::function<std::shared_ptr<typename tcf_type::sample_type const> ()> const&)
{
    return std::make_shared<tcf_type>();
}

template <int dimension, typename float_type>
void mean_quartic_displacement<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<mean_quartic_displacement>()

              , def("mean_quartic_displacement", &select_tcf_by_acquire<mean_quartic_displacement>)
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
