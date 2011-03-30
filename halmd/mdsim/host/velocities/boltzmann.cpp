/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace velocities
{

template <int dimension, typename float_type>
boltzmann<dimension, float_type>::boltzmann(
    shared_ptr<particle_type> particle
  , shared_ptr<random_type> random
  , double temperature
)
  : _Base(particle)
  // dependency injection
  , particle(particle)
  , random(random)
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

    BOOST_FOREACH (vector_type& v, particle->v) {
        // assign two components at a time
        for (unsigned i=0; i < dimension-1; i+=2) {
            tie(v[i], v[i+1]) = random->normal(sigma);
        }
        // handle last component separately for odd dimensions
        if (dimension % 2 == 1) {
            if (r_valid) {
                v[dimension-1] = r;
            }
            else {
                tie(v[dimension-1], r) = random->normal(sigma);
            }
            r_valid = !r_valid;
        }
        v_cm += v;
        vv += inner_prod(v, v);
    }

    v_cm /= particle->v.size();
    vv /= particle->v.size();
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
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    namespace_("velocities")
                    [
                        class_<boltzmann, shared_ptr<_Base_Base>, bases<_Base_Base, _Base> >(class_name.c_str())
                            .def(constructor<shared_ptr<particle_type>, shared_ptr<random_type>, double>())
                            .property("temperature", &boltzmann::temperature)
                            .property("module_name", &module_name_wrapper<dimension, float_type>)
                    ]
                ]
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &boltzmann<3, double>::luaopen
    ]
    [
        &boltzmann<2, double>::luaopen
    ];
#else
    [
        &boltzmann<3, float>::luaopen
    ]
    [
        &boltzmann<2, float>::luaopen
    ];
#endif
}

} // namespace

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class boltzmann<3, double>;
template class boltzmann<2, double>;
#else
template class boltzmann<3, float>;
template class boltzmann<2, float>;
#endif

}}} // namespace mdsim::host::velocities

} // namespace halmd
