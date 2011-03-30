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

#include <halmd/io/logger.hpp>
#include <halmd/observables/host/density_mode.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace host
{

template <int dimension, typename float_type>
density_mode<dimension, float_type>::density_mode(
    shared_ptr<phase_space_type> phase_space
  , shared_ptr<wavevector_type> wavevector
)
    // dependency injection
  : phase_space_(phase_space)
  , wavevector_(wavevector)
    // memory allocation
  , rho_sample_(phase_space_->sample->r.size(), wavevector_->value().size())
{
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void density_mode<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.sample, "sample", "computation of density modes");
}

/**
 * Acquire sample of all density modes from phase space sample
 */
template <int dimension, typename float_type>
void density_mode<dimension, float_type>::acquire(double time)
{
    scoped_timer<timer> timer_(runtime_.sample);

    if (rho_sample_.time == time) {
        LOG_TRACE("[density_mode] sample is up to date");
        return;
    }

    typedef typename phase_space_type::sample_type::sample_vector_ptr positions_vector_ptr_type;
    typedef typename density_mode_sample_type::mode_vector_type mode_vector_type;

    // trigger update of phase space sample
    on_acquire_(time);

    LOG_TRACE("[density_mode] acquire sample");

    if (phase_space_->sample->time != time) {
        throw logic_error("host phase space sample was not updated");
    }

    // compute density modes separately for each particle type
    // 1st loop: iterate over particle types
    unsigned int type = 0;
    BOOST_FOREACH (positions_vector_ptr_type const r_sample, phase_space_->sample->r) {
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
    rho_sample_.time = time;
}

template <int dimension, typename float_type>
void density_mode<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("density_mode_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("observables")
            [
                namespace_("host")
                [
                    class_<density_mode, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<phase_space_type>
                          , shared_ptr<wavevector_type>
                        >())
                        .def("register_runtimes", &density_mode::register_runtimes)
                ]
            ]
        ]
    ];
}

namespace  // limit symbols to translation unit
{

__attribute__ ((constructor)) void register_lua()
{
    lua_wrapper::register_(1)	//< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &density_mode<3, double>::luaopen
    ]
    [
        &density_mode<2, double>::luaopen
    ];
#else
    [
        &density_mode<3, float>::luaopen
    ]
    [
        &density_mode<2, float>::luaopen
    ];
#endif
}

}  // namespace

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
