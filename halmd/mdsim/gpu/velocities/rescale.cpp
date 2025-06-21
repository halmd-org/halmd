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
#include <halmd/utility/gpu/configure_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {

template <int dimension, typename float_type>
rescale<dimension, float_type>::rescale(
    std::shared_ptr<particle_type> particle
  , double target_energy
  , mode_selection mode
  , std::shared_ptr<halmd::logger> logger
)
  : particle_(particle)
  , mode_(mode)
  , logger_(logger)
{
    LOG("mode of operation: " << ((mode_ == nve) ? "nve" : "none"));
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
    // read access to potential energy of each particle
    auto const& en_pot = read_cache(particle_->potential_energy());

    scoped_timer_type timer(runtime_.set);

    LOG_DEBUG("rescale particle velocities to match target energy, for each particle");

    // access particle velocity array
    auto& velocity = *make_cache_mutable(particle_->velocity());
    cuda::memory::device::vector<int> retcode(1);

    if (mode_ == nve) {
        // configure and launch the rescale kernel
        // This performs for each particle:
        // (1) calculate total energy and scaling factor, (2) apply velocity scaling
        configure_kernel(wrapper_type::kernel.rescale_nve, particle_->dim(), true);

        wrapper_type::kernel.rescale_nve(
            velocity.data()               // velocities (device pointer)
          , en_pot.data()                 // potential energy (device pointer)
          , particle_->nparticle()        // number of particles
          , target_energy_                // target energy per particle
          , retcode.data()                // return code
        );
        // cuda::thread::synchronize();   // cuda::copy below is blocking
    }
    else {
        LOG_ERROR("mode of operation is not yet supported");
    }

    // obtain return code and test bits
    int r;
    cuda::copy(retcode.begin(), retcode.begin() + 1, &r);

    if (r & warning) {
        LOG_WARNING_ONCE("kinetic energy is zero for some particle(s). Initialising first velocity component");
    }
    else if (r & failure) {
        LOG_ERROR("target energy is less than potential energy");
        throw std::runtime_error("target energy can not be matched by velocity rescaling");
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
                  , mode_selection
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
