/*
 * Copyright © 2010 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/mdsim/host/velocity.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace velocities {

template <int dimension, typename float_type>
boltzmann<dimension, float_type>::boltzmann(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<random_type> random
  , double temperature
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , random_(random)
  , logger_(logger)
  // set parameters
  , temp_(temperature)
{
    LOG("Boltzmann velocity distribution temperature: T = " << temp_);
}

template <int dimension, typename float_type>
void boltzmann<dimension, float_type>::set()
{
    scoped_timer_type timer(runtime_.set);

    LOG_TRACE("assign Boltzmann-distributed velocities");

    auto velocity = make_cache_mutable(particle_->velocity());
    mass_array_type const& mass = read_cache(particle_->mass());
    size_type nparticle = particle_->nparticle();

    float_type const sigma = std::sqrt(temp_);
    fixed_vector<double, dimension> mv = 0;
    double mv2 = 0;
    double m = 0;
    float_type r;
    bool r_valid = false;

    for (size_type i = 0; i < nparticle; ++i) {
        vector_type& v = (*velocity)[i];
        // assign two components at a time
        for (unsigned int j = 0; j < dimension - 1; j += 2) {
            boost::tie(v[j], v[j + 1]) = random_->normal(sigma);
        }
        // handle last component separately for odd dimensions
        if (dimension % 2 == 1) {
            if (r_valid) {
                v[dimension - 1] = r;
            }
            else {
                boost::tie(v[dimension - 1], r) = random_->normal(sigma);
            }
            r_valid = !r_valid;
        }
        double m_i = mass[i];
        v /= std::sqrt(m_i);
        mv += m_i * v;
        mv2 += m_i * inner_prod(v, v);
        m += m_i;
    }

    fixed_vector<double, dimension> v_cm = mv / m;
    double scale = std::sqrt(nparticle * temp_ * dimension / (mv2 - m * inner_prod(v_cm, v_cm)));
    LOG_TRACE("shift velocities by " << -v_cm << " and resacle by factor " << scale);
    shift_rescale_velocity(*particle_, -v_cm, scale);
}

template <int dimension, typename float_type>
void boltzmann<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("velocities")
            [
                class_<boltzmann>()
                    .property("temperature", &boltzmann::temperature)
                    .def("set", &boltzmann::set)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("set", &runtime::set)
                    ]
                    .def_readonly("runtime", &boltzmann::runtime_)

              , def("boltzmann", &std::make_shared<boltzmann
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<random_type>
                  , double
                  , std::shared_ptr<logger>
                >)
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

} // namespace velocities
} // namespace host
} // namespace mdsim
} // namespace halmd
