/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#include <functional>
#include <stdexcept>

#include <halmd/observables/ssf.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/signal.hpp>

using namespace std;

namespace halmd {
namespace observables {

template <int dimension>
ssf<dimension>::ssf(
    std::shared_ptr<wavevector_type const> wavevector
  , double norm
  , std::shared_ptr<clock_type const> clock
  , std::shared_ptr<logger_type> logger
)
  // dependency injection
  : wavevector_(wavevector)
  , norm_(norm)
  , clock_(clock)
  , logger_(logger)
  // initialise members
  , result_(wavevector_->wavenumber().size())
  , step_(numeric_limits<step_type>::max())
{
}

template <int dimension>
typename ssf<dimension>::result_type const&
ssf<dimension>::sample(density_mode_type const& mode1, density_mode_type const& mode2)
{
    scoped_timer_type timer(runtime_.sample);

    if (step_ == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return result_;
    }

    LOG_TRACE("sampling");

    typename rho_vector_type::const_iterator rho1_begin = mode1.rho().begin();
    typename rho_vector_type::const_iterator rho2_begin = mode2.rho().begin();
    typename result_type::iterator result = result_.begin();

    // accumulate products of density modes with equal wavenumber,
    // iterate over wavevector shells encoded as index ranges to the wavevector array
    for (auto idx_range : wavevector_->shell()) {
        accumulator<double> acc;
        auto rho_q1 = rho1_begin + idx_range.first;
        auto rho_q2 = rho2_begin + idx_range.first;
        // iterate over wavevectors and density modes simultaneously
        for (size_t i = idx_range.first; i != idx_range.second; ++i, ++rho_q1, ++rho_q2) {
            // compute Re[rho_q rho_q^*]
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

    // store simulation step as time stamp
    step_ = clock_->step();

    return result_;
}

template <typename result_type, typename ssf_type, typename slot_type>
static std::function<result_type const& ()>
wrap_sample(std::shared_ptr<ssf_type> ssf, slot_type const& mode1, slot_type const& mode2)
{
    return [=]() -> result_type const& {
        return ssf->sample(*mode1(), *mode2());
    };
}

template <int dimension>
void ssf<dimension>::luaopen(lua_State* L)
{
    typedef std::function<std::shared_ptr<density_mode_type const> ()> slot_type;

    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<ssf>()
                .def("sample", &wrap_sample<result_type, ssf, slot_type>)
                .property("wavevector", &ssf::wavevector)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("sample", &runtime::sample)
                ]
                .def_readonly("runtime", &ssf::runtime_)

          , def("ssf", &std::make_shared<ssf
              , std::shared_ptr<wavevector_type const>
              , double
              , std::shared_ptr<clock_type const>
              , std::shared_ptr<logger_type>
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
