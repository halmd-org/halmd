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
    std::shared_ptr<thermo_type> thermo,
    double target_energy,
    std::shared_ptr<halmd::logger> logger
)
  : particle_(particle),
    thermo_(thermo),
    target_energy_(target_energy),
    logger_(logger)
{}

template <int dimension, typename float_type>
void rescale<dimension, float_type>::set()
{
    auto const& en_pot = read_cache(particle_->potential_energy());

    scoped_timer_type timer(runtime_.set);

    LOG_DEBUG("Rescale particle velocities to match target energy, for each particle");

    auto velocity = make_cache_mutable(particle_->velocity());
    std::size_t n = particle_->nparticle();
    float_type target_total_energy = static_cast<float_type>(target_energy_);
    
    // Loop over all particles to compute current total energy and rescale velocities
    for (std::size_t i = 0; i < n; ++i) {
        float_type ke = 0.0;
        for (int d = 0; d < dimension; ++d)
            ke += 0.5 * velocity(i, d) * velocity(i, d);  // assuming mass = 1

        float_type total_energy = ke + en_pot(i);
        float_type scaling = std::sqrt(target_total_energy / total_energy);

        for (int d = 0; d < dimension; ++d)
            velocity(i, d) *= scaling;
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
                    .def("set_target_energy", &rescale::set_target_energy)
                    .def("target_energy", &rescale::target_energy)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("set", &runtime::set)
                    ]
                    .def_readonly("runtime", &rescale::runtime_)
              , def("rescale", &std::make_shared<rescale
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<thermo_type>
                  , double
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_velocities_rescale(lua_State* L)
{
#ifdef USE_HOST_SINGLE_PRECISION
    rescale<3, float>::luaopen(L);
    rescale<2, float>::luaopen(L);
#endif
#ifdef USE_HOST_DOUBLE_SINGLE_PRECISION
    rescale<3, dsfloat>::luaopen(L);
    rescale<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

#ifdef USE_HOST_SINGLE_PRECISION
template class rescale<3, float>;
template class rescale<2, float>;
#endif
#ifdef USE_HOST_DOUBLE_SINGLE_PRECISION
template class rescale<3, dsfloat>;
template class rescale<2, dsfloat>;
#endif

} // namespace velocities
} // namespace host
} // namespace mdsim
} // namespace halmd
