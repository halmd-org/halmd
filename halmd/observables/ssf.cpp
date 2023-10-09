/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#include <functional>
#include <stdexcept>

#include <halmd/observables/ssf.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/owner_equal.hpp>
#include <halmd/utility/signal.hpp>

using namespace std;

namespace halmd {
namespace observables {

template <int dimension>
ssf<dimension>::ssf(
    mode_acquisitor_type mode1
  , mode_acquisitor_type mode2
  , shared_ptr<wavevector_type const> wavevector
  , double norm
  , shared_ptr<logger> logger
)
  // dependency injection
  : mode1_(mode1)
  , mode2_(mode2)
  , wavevector_(wavevector)
  , logger_(logger)
  // initialise members
  , norm_(norm)
  , result_(wavevector_->wavenumber().size())
{
    LOG("normalisation factor: " << norm)
}

template <int dimension>
typename ssf<dimension>::result_type const&
ssf<dimension>::sample()
{
    // retrieve cached density modes stored within std::shared_ptr
    mode_type mode1 = mode1_();
    mode_type mode2 = mode2_();

    // track the validity of the caches via std::weak_ptr
    if (!owner_equal(mode1_observer_, mode1) || !owner_equal(mode2_observer_, mode2))
    {
        LOG_DEBUG("sampling");

        scoped_timer_type timer(runtime_.sample);

        assert(mode1->size() == mode2->size());
        assert(mode1->size() == wavevector_->value().size());

        // accumulate products of density modes with equal wavenumber,
        // iterate over wavevector shells encoded as index ranges to the wavevector array
        auto result = begin(result_);
        auto rho1_begin = begin(*mode1);
        auto rho2_begin = begin(*mode2);
        for (auto idx_range : wavevector_->shell()) {
            accumulator<double> acc;
            auto rho_q1 = rho1_begin + idx_range.first;
            auto rho_q2 = rho2_begin + idx_range.first;
            // iterate over wavevectors and density modes simultaneously
            for (size_t i = idx_range.first; i != idx_range.second; ++i, ++rho_q1, ++rho_q2) {
                // compute Re[rho_q1 (rho_q2)*]
                double re = ((*rho_q1)[0] * (*rho_q2)[0] + (*rho_q1)[1] * (*rho_q2)[1]);
                // accumulate results for this wavenumber
                acc(re / norm_);
            }
            // transform accumulator to array (mean, error_of_mean, count)
            // and write to output iterator
            *result++ = {{
                mean(acc)
              , count(acc) > 1 ? error_of_mean(acc) : 0
              , static_cast<double>(count(acc))
            }};
        }

        // update observers
        mode1_observer_ = mode1;
        mode2_observer_ = mode2;
    }

    return result_;
}

template <int dimension>
void ssf<dimension>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<ssf>()
                .property("sampler", &ssf::sampler)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("sample", &runtime::sample)
                ]
                .def_readonly("runtime", &ssf::runtime_)

          , def("ssf", &std::make_shared<ssf
              , mode_acquisitor_type
              , mode_acquisitor_type
              , std::shared_ptr<wavevector_type const>
              , double
              , std::shared_ptr<logger>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_ssf(lua_State* L)
{
    ssf<3>::luaopen(L);
    ssf<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class ssf<3>;
template class ssf<2>;

} // namespace observables
} // namespace halmd
