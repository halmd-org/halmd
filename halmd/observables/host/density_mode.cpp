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
    shared_ptr<phase_space_type const> phase_space
  , shared_ptr<wavevector_type const> wavevector
  , shared_ptr<clock_type const> clock
  , shared_ptr<logger_type> logger
)
    // dependency injection
  : phase_space_(phase_space)
  , wavevector_(wavevector)
  , clock_(clock)
  , logger_(logger)
    // memory allocation
  , rho_sample_(phase_space_->r.size(), wavevector_->value().size())
{
}

/**
 * Acquire sample of all density modes from phase space sample
 */
template <int dimension, typename float_type>
void density_mode<dimension, float_type>::acquire()
{
    scoped_timer_type timer(runtime_.acquire);

    if (rho_sample_.step == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return;
    }

    typedef typename phase_space_type::sample_vector_ptr positions_vector_ptr_type;
    typedef typename density_mode_sample_type::mode_vector_type mode_vector_type;

    // trigger update of phase space sample
    on_acquire_();

    LOG_TRACE("acquire sample");

    if (phase_space_->step != clock_->step()) {
        throw logic_error("host phase space sample was not updated");
    }

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    rho_sample_.reset();

    // compute density modes separately for each particle type
    // 1st loop: iterate over particle types
    unsigned int type = 0;
    BOOST_FOREACH (positions_vector_ptr_type const r_sample, phase_space_->r) {
        mode_vector_type& rho_vector = *rho_sample_.rho[type]; //< dereference shared_ptr
        // initialise result array
        fill(rho_vector.begin(), rho_vector.end(), 0);
        // compute sum of exponentials: rho_q = sum_r exp(-i q·r)
        // 2nd loop: iterate over particles of the same type
        BOOST_FOREACH (vector_type const& r, *r_sample) {
            typename mode_vector_type::iterator rho_q = rho_vector.begin();
            typedef pair<double, vector_type> map_value_type; // pair: (wavenumber, wavevector)
            // 3rd loop: iterate over wavevectors
            BOOST_FOREACH (map_value_type const& q_pair, wavevector_->value()) {
                float_type q_r = inner_prod(static_cast<vector_type>(q_pair.second), r);
                *rho_q++ += mode_type(cos(q_r), -sin(q_r));
            }
        }
        ++type;
    }
    rho_sample_.step = clock_->step();
}

template <int dimension, typename float_type>
void density_mode<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("density_mode_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<density_mode, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .def(constructor<
                        shared_ptr<phase_space_type const>
                      , shared_ptr<wavevector_type const>
                      , shared_ptr<clock_type const>
                      , shared_ptr<logger_type>
                    >())
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                    ]
                    .def_readonly("runtime", &density_mode::runtime_)
            ]
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
