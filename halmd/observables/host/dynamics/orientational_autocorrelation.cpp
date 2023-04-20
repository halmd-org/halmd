/*
 * Copyright © 2016       Manuel Dibak
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
#include <halmd/observables/host/dynamics/orientational_autocorrelation.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace dynamics {

/**
 * Compute orientation autocorrelation of two orientation sample vectors.
 *
 * @param first particles velocities of one species at time t1
 * @param second particles velocities of one species at time t2
 * @returns accumulated velocity autocorrelation
 */
template <int dimension, typename float_type>
void orientational_autocorrelation<dimension, float_type>::operator() (
    sample_type const& first
  , sample_type const& second
  , accumulator<result_type>& result
)
{
    accumulator<result_type> acc;
    typename sample_type::array_type::const_iterator u1, u2, end = first.data().end();

    for (u1 = first.data().begin(), u2 = second.data().begin(); u1 != end; ++u1, ++u2) {
        // accumulate orientational correlation
        acc(correlate_function_type()(*u1, *u2));
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
void orientational_autocorrelation<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<orientational_autocorrelation>()

              , def("orientational_autocorrelation", &select_tcf_by_acquire<orientational_autocorrelation>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_dynamics_orientational_autocorrelation(lua_State* L)
{
    orientational_autocorrelation<3, double>::luaopen(L);
    orientational_autocorrelation<2, double>::luaopen(L);
    orientational_autocorrelation<3, float>::luaopen(L);
    orientational_autocorrelation<2, float>::luaopen(L);
    observables::dynamics::correlation<orientational_autocorrelation<3, double>>::luaopen(L);
    observables::dynamics::correlation<orientational_autocorrelation<2, double>>::luaopen(L);
    observables::dynamics::correlation<orientational_autocorrelation<3, float>>::luaopen(L);
    observables::dynamics::correlation<orientational_autocorrelation<2, float>>::luaopen(L);
    return 0;
}

// explicit instantiation
template class orientational_autocorrelation<3, double>;
template class orientational_autocorrelation<2, double>;
template class orientational_autocorrelation<3, float>;
template class orientational_autocorrelation<2, float>;

} // namespace dynamics
} // namespace host

namespace dynamics {

// explicit instantiation
template class correlation<host::dynamics::orientational_autocorrelation<3, double>>;
template class correlation<host::dynamics::orientational_autocorrelation<2, double>>;
template class correlation<host::dynamics::orientational_autocorrelation<3, float>>;
template class correlation<host::dynamics::orientational_autocorrelation<2, float>>;

} // namespace dynamics
} // namespace observables
} // namespace halmd
