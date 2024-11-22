/*
 * Copyright © 2024 Felix Höfling
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

#include <memory>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/host/dynamics/self_intermediate_scattering_function.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace dynamics {

template <int dimension, typename float_type>
template <typename MultiArray>
void self_intermediate_scattering_function<dimension, float_type>::operator() (
    sample_type const& first
  , sample_type const& second
  , MultiArray&& result
) const
{
    // compute sum of exponentials over all particles,
    // F(q, t2 - t1) = sum_j exp(-i q · [R(t2) - R(t1)])
    // FIXME for simplicity, we return only the real part
    //
    // 1st loop: iterate over particles in the sample
    auto r2_it = second.data().begin();
    for (auto const& r1 : first.data()) {
        auto const& r2 = *r2_it++;
        // displacement of this particle: R(t2) - R(t1)
        vector_type dr = r2 - r1;

        // 2nd loop: iterate over wavevector shells,
        // encoded as index ranges to the wavevector array
        auto output = result.begin();
        for (auto idx_range : wavevector_->shell()) {
            // accumulate results for equal wavenumber
            accumulator<result_type> acc;

            // iterate over wavevectors in index range
            auto const& q = wavevector_->value();
            for (size_t i = idx_range.first; i != idx_range.second; ++i) {
                float_type q_dr = inner_prod(static_cast<vector_type>(q[i]), dr);
                acc(cos(q_dr)); // FIXME complex_type({{ cos(q_r), -sin(q_r) }});
            }

            // write to output iterator, which accumulates the result
            (*output++)(acc);
        }
    }
}

template <typename tcf_type>
static std::shared_ptr<tcf_type>
select_tcf_by_acquire(
    std::function<std::shared_ptr<typename tcf_type::sample_type const> ()> const&
  , std::shared_ptr<typename tcf_type::wavevector_type const> wavevector
)
{
    return std::make_shared<tcf_type>(wavevector);
}

template <int dimension, typename float_type>
void self_intermediate_scattering_function<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<self_intermediate_scattering_function>()

              , def("self_intermediate_scattering_function", &select_tcf_by_acquire<self_intermediate_scattering_function>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_dynamics_self_intermediate_scattering_function(lua_State* L)
{
    // provide module for both floating-point types since
    // it may also process GPU samples that are in single precision.
    self_intermediate_scattering_function<3, float>::luaopen(L);
    self_intermediate_scattering_function<2, float>::luaopen(L);
    self_intermediate_scattering_function<3, double>::luaopen(L);
    self_intermediate_scattering_function<2, double>::luaopen(L);
    observables::dynamics::correlation<self_intermediate_scattering_function<3, float>>::luaopen(L);
    observables::dynamics::correlation<self_intermediate_scattering_function<2, float>>::luaopen(L);
    observables::dynamics::correlation<self_intermediate_scattering_function<3, double>>::luaopen(L);
    observables::dynamics::correlation<self_intermediate_scattering_function<2, double>>::luaopen(L);

    return 0;
}

// explicit instantiation, for both floating-point types (see above)
template class self_intermediate_scattering_function<3, float>;
template class self_intermediate_scattering_function<2, float>;
template class self_intermediate_scattering_function<3, double>;
template class self_intermediate_scattering_function<2, double>;

} // namespace dynamics
} // namespace host

namespace dynamics {

// explicit instantiation, for both floating-point types (see above)
template class correlation<host::dynamics::self_intermediate_scattering_function<3, float>>;
template class correlation<host::dynamics::self_intermediate_scattering_function<2, float>>;
template class correlation<host::dynamics::self_intermediate_scattering_function<3, double>>;
template class correlation<host::dynamics::self_intermediate_scattering_function<2, double>>;

} // namespace dynamics
} // namespace observables
} // namespace halmd
