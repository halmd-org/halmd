/*
 * Copyright © 2011-2012  Felix Höfling and Peter Colberg
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

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <stdexcept>

#include <halmd/observables/ssf.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {

template <int dimension>
ssf<dimension>::ssf(
    shared_ptr<wavevector_type const> wavevector
  , double norm
  , shared_ptr<clock_type const> clock
  , shared_ptr<logger_type> logger
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

    if (mode1.step != clock_->step()) {
        throw logic_error("first density modes sample was not updated");
    }
    if (mode2.step != clock_->step()) {
        throw logic_error("second density modes sample was not updated");
    }

    typename rho_vector_type::const_iterator rho_q1 = mode1.rho->begin();
    typename rho_vector_type::const_iterator rho_q2 = mode2.rho->begin();
    typename result_type::iterator result = result_.begin();

    typename wavevector_type::map_type const& wavevector = wavevector_->value();
    typename wavevector_type::map_type::const_iterator q = wavevector.begin();
    typename wavevector_type::map_type::const_iterator q_next = q; ++q_next;

    // accumulate products of density modes with equal wavenumber,
    // iterate over sorted list of wavevectors
    accumulator<double> acc;
    while (q != wavevector.end()) {
        // compute Re[rho_q rho_q^*], add result to output accumulator
        double re = (real(*rho_q1) * real(*rho_q2) + imag(*rho_q1) * imag(*rho_q2));
        acc(re / norm_);
        // find end of range with equal wavenumber
        if (q_next == wavevector.end() || q->first != q_next->first) {
            // transform accumulator to array (mean, error_of_mean, count)
            (*result)[0] = mean(acc);
            (*result)[1] = count(acc) > 1 ? error_of_mean(acc) : 0;
            (*result)[2] = count(acc);
            acc.reset();
            // next wavenumber: increment output iterator
            ++result;
        }
        ++q; ++q_next; ++rho_q1; ++rho_q2;
    }

    // store simulation step as time stamp
    step_ = clock_->step();

    return result_;
}

template <typename result_type, typename ssf_type, typename slot_type>
static result_type const&
sample(shared_ptr<ssf_type> ssf, slot_type const& mode1, slot_type const& mode2)
{
    return ssf->sample(*mode1(), *mode2());
}

template <typename result_type, typename ssf_type, typename slot_type>
static function<result_type const& ()>
wrap_sample(shared_ptr<ssf_type> ssf, slot_type const& mode1, slot_type const& mode2)
{
    return bind(&sample<result_type, ssf_type, slot_type>, ssf, mode1, mode2);
}

template <int dimension>
void ssf<dimension>::luaopen(lua_State* L)
{
    typedef function<shared_ptr<density_mode_type const> ()> slot_type;

    using namespace luabind;
    static string class_name("ssf_" + lexical_cast<string>(dimension));
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<ssf>(class_name.c_str())
                .def("sample", &wrap_sample<result_type, ssf, slot_type>)
                .property("wavevector", &ssf::wavevector)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("sample", &runtime::sample)
                ]
                .def_readonly("runtime", &ssf::runtime_)

          , def("ssf", &make_shared<ssf
              , shared_ptr<wavevector_type const>
              , double
              , shared_ptr<clock_type const>
              , shared_ptr<logger_type>
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
