/*
 * Copyright © 2008-2011  Peter Colberg
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
#include <boost/tuple/tuple.hpp>

#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace velocities {

template <int dimension, typename float_type>
boltzmann<dimension, float_type>::boltzmann(
    boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<random_type> random
  , double temperature
  , boost::shared_ptr<logger_type> logger
)
  : _Base(particle, logger)
  // dependency injection
  , particle_(particle)
  , random_(random)
  , logger_(logger)
  // set parameters
  , temp_(temperature)
{
    LOG("Boltzmann velocity distribution temperature: T = " << temp_);
}

/**
 * Initialise velocities from Maxwell-Boltzmann distribution
 */
template <int dimension, typename float_type>
void boltzmann<dimension, float_type>::set()
{
    scoped_timer_type timer(runtime_.set);

    // assuming equal (unit) mass for all particle types
    vector_type v_cm;
    float_type vv;
    tie(v_cm, vv) = gaussian(sqrt(temp_));

    // center velocities around origin, then rescale to exactly
    // match the desired temperature;
    // temp = vv / dimension
    // vv changes to vv - v_cm^2 after shifting
    float_type scale = sqrt(temp_ * dimension / (vv - inner_prod(v_cm, v_cm)));
    boltzmann::shift_rescale(-v_cm, scale);

    LOG_DEBUG("velocities rescaled by factor " << scale);
    LOG_DEBUG("assigned Boltzmann-distributed velocities");
}

/**
 * Assign new velocities from Gaussian distribution
 */
template <int dimension, typename float_type>
pair<typename boltzmann<dimension, float_type>::vector_type, float_type>
inline boltzmann<dimension, float_type>::gaussian(float_type sigma)
{
    vector_type v_cm = 0;
    float_type vv = 0;
    float_type r = 0;
    bool r_valid = false;

    BOOST_FOREACH (vector_type& v, particle_->velocity()) {
        // assign two components at a time
        for (unsigned i=0; i < dimension-1; i+=2) {
            tie(v[i], v[i+1]) = random_->normal(sigma);
        }
        // handle last component separately for odd dimensions
        if (dimension % 2 == 1) {
            if (r_valid) {
                v[dimension-1] = r;
            }
            else {
                tie(v[dimension-1], r) = random_->normal(sigma);
            }
            r_valid = !r_valid;
        }
        v_cm += v;
        vv += inner_prod(v, v);
    }

    v_cm /= particle_->velocity().size();
    vv /= particle_->velocity().size();
    return make_pair(v_cm, vv);
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(boltzmann<dimension, float_type> const&)
{
    return boltzmann<dimension, float_type>::module_name();
}

template <int dimension, typename float_type>
void boltzmann<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("velocities")
                [
                    class_<boltzmann, boost::shared_ptr<boltzmann> >(class_name.c_str())
                        .def(constructor<
                            boost::shared_ptr<particle_type>
                           , boost::shared_ptr<random_type>
                           , double
                           , boost::shared_ptr<logger_type>
                        >())
                        .property("temperature", &boltzmann::temperature)
                        .property("module_name", &module_name_wrapper<dimension, float_type>)
                        .def("set", &boltzmann::set)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("set", &runtime::set)
                        ]
                        .def_readonly("runtime", &boltzmann::runtime_)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_velocities_boltzmann(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    boltzmann<3, double>::luaopen(L);
    boltzmann<2, double>::luaopen(L);
#else
    boltzmann<3, float>::luaopen(L);
    boltzmann<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class boltzmann<3, double>;
template class boltzmann<2, double>;
#else
template class boltzmann<3, float>;
template class boltzmann<2, float>;
#endif

} // namespace mdsim
} // namespace host
} // namespace velocities
} // namespace halmd
