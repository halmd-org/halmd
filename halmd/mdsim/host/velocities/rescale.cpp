/*
 * Copyright © 2025 Giorgia Marcelli
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

#include <halmd/mdsim/host/velocities/rescale.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <cmath>

namespace halmd {
namespace mdsim {
namespace host {
namespace velocities {

template <int dimension, typename float_type>
rescale<dimension, float_type>::rescale(
    std::shared_ptr<particle_type> particle,
    double target_energy,
    std::shared_ptr<halmd::logger> logger
)
  : particle_(particle),
    logger_(logger)
{
    set_target_energy(target_energy);
}

template <int dimension, typename float_type>
void rescale<dimension, float_type>::set_target_energy(double energy)
{
    target_energy_ = energy;
    LOG("target energy: " << target_energy_);
}

template <int dimension, typename float_type>
void rescale<dimension, float_type>::set()
{
    auto const& en_pot_array = read_cache(particle_->potential_energy());
    auto const& mass = read_cache(particle_->mass());

    scoped_timer_type timer(runtime_.set);

    LOG_DEBUG("rescale particle velocities to match target energy, for each particle");

    auto& velocity = *make_cache_mutable(particle_->velocity());

    for (size_type i = 0; i < particle_->nparticle(); ++i)
    {
        auto& v = velocity[i];
        float_type mass_ = mass[i];
        float_type en_kin = mass_ * inner_prod(v, v) / 2;
        float_type en_pot = en_pot_array[i];

        float_type target_total_energy = static_cast<float_type>(target_energy_);
        float_type energy_diff = target_total_energy - en_pot;

        if (energy_diff < float_type(0)) {
            LOG_ERROR("target energy (" << target_total_energy
                      << ") is less than potential energy (" << en_pot
                      << ") of particle #" << i
                     );
            throw std::runtime_error("target energy can not be matched by velocity rescaling");
        }

        // safety guard: handle zero kinetic energy (i.e., v = 0)
        if (en_kin == float_type(0)) {
            LOG_WARNING_ONCE("kinetic energy is zero for some particle(s). Initialising first velocity component.");
            LOG_DEBUG("zero velocity of particle #" << i);

            // let v point along the first axis, set magnitude to match the desired kinetic energy
            v[0] = std::sqrt(2 * energy_diff / mass_);
        }
        else {
            // compute rescaling factor and apply to all velocities
            float_type scaling = std::sqrt(energy_diff / en_kin);
            v *= scaling;
        }
    }
}


template <int dimension, typename float_type>
void rescale<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("velocities")
            [
                class_<rescale>()
                    .def("set", &rescale::set)
                    .property("target_energy", &rescale::target_energy, &rescale::set_target_energy)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("set", &runtime::set)
                    ]
                    .def_readonly("runtime", &rescale::runtime_)

              , def("rescale", &std::make_shared<rescale
                  , std::shared_ptr<particle_type>
                  , double
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_velocities_rescale(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    rescale<3, double>::luaopen(L);
    rescale<2, double>::luaopen(L);
#else
    rescale<3, float>::luaopen(L);
    rescale<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class rescale<3, double>;
template class rescale<2, double>;
#else
template class rescale<3, float>;
template class rescale<2, float>;
#endif

} // namespace velocities
} // namespace host
} // namespace mdsim
} // namespace halmd
