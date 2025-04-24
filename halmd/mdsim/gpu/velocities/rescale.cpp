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

#include <halmd/mdsim/gpu/velocities/rescale.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <cmath>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {

template <int dimension, typename float_type>
rescale<dimension, float_type>::rescale(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<thermo_type> thermo
  , double target_energy
  , std::shared_ptr<halmd::logger> logger
)
  : particle_(particle)
  , thermo_(thermo)
  , target_energy_(target_energy)
  , logger_(logger)
{}

template <int dimension, typename float_type>
void rescale<dimension, float_type>::set()
{
    // read access to potential energy of each particle
    auto const& en_pot = read_cache(particle_->potential_energy());

    scoped_timer_type timer(runtime_.set);

    LOG_DEBUG("Rescale particle velocities to match target energy, for each particle");

    // Access particle velocity buffer
    auto velocity = make_cache_mutable(particle_->velocity());

    // Configure and launch the rescale kernel
    // This performs for each particle: (1) calculate total energy and scaling factor, (2) apply velocity scaling
    configure_kernel(wrapper_type::kernel.rescale, particle_->dim(), true);

    wrapper_type::kernel.rescale(
        velocity->data(),             // velocities (device pointer)
        en_pot.data(),                // potential energy (device pointer)
        particle_->nparticle(),       // number of particles
        target_energy_                // target energy per particle
    );
    cuda::thread::synchronize();
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

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_velocities_rescale(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    rescale<3, float>::luaopen(L);
    rescale<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    rescale<3, dsfloat>::luaopen(L);
    rescale<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

#ifdef USE_GPU_SINGLE_PRECISION
template class rescale<3, float>;
template class rescale<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class rescale<3, dsfloat>;
template class rescale<2, dsfloat>;
#endif

} // namespace velocities
} // namespace gpu
} // namespace mdsim
} // namespace halmd
