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
#include <boost/foreach.hpp>

#include <halmd/observables/host/density_mode.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
density_mode<dimension, float_type>::density_mode(
    boost::shared_ptr<wavevector_type const> wavevector
  , boost::shared_ptr<clock_type const> clock
  , boost::shared_ptr<logger_type> logger
)
    // dependency injection
  : wavevector_(wavevector)
  , clock_(clock)
  , logger_(logger)
{
}

/**
 * Acquire sample of all density modes from phase space sample
 */
template <int dimension, typename float_type>
boost::shared_ptr<typename density_mode<dimension, float_type>::sample_type const>
density_mode<dimension, float_type>::acquire(phase_space_type const& phase_space)
{
    scoped_timer_type timer(runtime_.acquire);

    if (rho_sample_ && rho_sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return rho_sample_;
    }

    LOG_TRACE("acquire sample");

    if (phase_space.step() != clock_->step()) {
        throw logic_error("host phase space sample was not updated");
    }

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    rho_sample_ = boost::make_shared<sample_type>(wavevector_->value().size(), clock_->step());

    // compute density modes
    mode_array_type& rho_vector = rho_sample_->rho();
    // initialise result array
    fill(rho_vector.begin(), rho_vector.end(), 0);
    // compute sum of exponentials: rho_q = sum_r exp(-i q·r)
    // 1st loop: iterate over particles
    BOOST_FOREACH (vector_type const& r, phase_space.position()) {
        typename mode_array_type::iterator rho_q = rho_vector.begin();
        typedef pair<double, vector_type> map_value_type; // pair: (wavenumber, wavevector)
        // 2nd loop: iterate over wavevectors
        BOOST_FOREACH (map_value_type const& q_pair, wavevector_->value()) {
            float_type q_r = inner_prod(static_cast<vector_type>(q_pair.second), r);
            *rho_q++ += mode_type(cos(q_r), -sin(q_r));
        }
    }

    return rho_sample_;
}

template <typename sample_type, typename density_mode_type, typename slot_type>
static boost::shared_ptr<sample_type const>
acquire(boost::shared_ptr<density_mode_type> density_mode, slot_type const& phase_space)
{
    return density_mode->acquire(*phase_space());
}

template <typename sample_type, typename density_mode_type, typename slot_type>
static boost::function<boost::shared_ptr<sample_type const> ()>
wrap_acquire(boost::shared_ptr<density_mode_type> density_mode, slot_type const& phase_space)
{
    return bind(&acquire<sample_type, density_mode_type, slot_type>, density_mode, phase_space);
}

template <int dimension, typename float_type>
void density_mode<dimension, float_type>::luaopen(lua_State* L)
{
    typedef boost::function<boost::shared_ptr<phase_space_type const> ()> slot_type;

    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<density_mode>()
                .def("acquire", &wrap_acquire<sample_type, density_mode, slot_type>)
                .property("wavevector", &density_mode::wavevector)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("acquire", &runtime::acquire)
                ]
                .def_readonly("runtime", &density_mode::runtime_)

          , def("density_mode", &boost::make_shared<density_mode
              , boost::shared_ptr<wavevector_type const>
              , boost::shared_ptr<clock_type const>
              , boost::shared_ptr<logger_type>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_density_mode(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    density_mode<3, double>::luaopen(L);
    density_mode<2, double>::luaopen(L);
#else
    density_mode<3, float>::luaopen(L);
    density_mode<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class density_mode<3, double>;
template class density_mode<2, double>;
#else
template class density_mode<3, float>;
template class density_mode<2, float>;
#endif

}}  // namespace observables::host

}  // namespace halmd
