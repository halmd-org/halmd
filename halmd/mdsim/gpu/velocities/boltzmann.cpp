/*
 * Copyright © 2010-2014 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {

template <int dimension, typename float_type, typename RandomNumberGenerator>
boltzmann<dimension, float_type, RandomNumberGenerator>::boltzmann(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<random_type> random
  , double temperature
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , random_(random)
  , logger_(logger)
  // allocate GPU memory
  , g_mv_(random_->rng().dim.blocks_per_grid())
  , g_mv2_(random_->rng().dim.blocks_per_grid())
  , g_m_(random_->rng().dim.blocks_per_grid())
{
    set_temperature(temperature);
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void boltzmann<dimension, float_type, RandomNumberGenerator>::set_temperature(double temperature)
{
    temp_ =  temperature;
    LOG("temperature of Boltzmann distribution: " << float(temp_));
}

/**
 * Initialise velocities from Maxwell-Boltzmann distribution
 *
 * The particle velocities need to fullfill two constraints:
 *
 *  1. centre of mass velocity shall be zero
 *  2. temperature of the distribution shall equal exactly the given value
 *
 * The above order is chosen as shifting the centre of mass velocity
 * means altering the first moment of the velocity distribution, which
 * in consequence affects the second moment, i.e. the temperature.
 *
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void boltzmann<dimension, float_type, RandomNumberGenerator>::set()
{
    auto velocity = make_cache_mutable(particle_->velocity());

    scoped_timer_type timer(runtime_.set);

    // generate Maxwell-Boltzmann distributed velocities,
    // assuming equal (unit) mass for all particle types
    wrapper_type::kernel.gaussian.configure(random_->rng().dim.grid, random_->rng().dim.block);
    wrapper_type::kernel.gaussian(
        velocity->data()
      , particle_->nparticle()
      , particle_->dim().threads()
      , temp_
      , g_mv_
      , g_mv2_
      , g_m_
      , random_->rng().rng()
    );
    cuda::thread::synchronize();

    // set centre of mass velocity to zero and
    // rescale velocities to accurate temperature
    wrapper_type::kernel.shift.configure(particle_->dim().grid,
      particle_->dim().block, g_mv2_.size() * (2 + dimension) * sizeof(dsfloat));
    wrapper_type::kernel.shift(
        velocity->data()
      , particle_->nparticle()
      , particle_->dim().threads()
      , temp_
      , g_mv_
      , g_mv2_
      , g_m_
      , g_mv2_.size()
    );
    cuda::thread::synchronize();

    LOG_TRACE("assigned Boltzmann-distributed velocities");
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void boltzmann<dimension, float_type, RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("velocities")
            [
                class_<boltzmann>()
                    .property("temperature", &boltzmann::temperature, &boltzmann::set_temperature)
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

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_velocities_boltzmann(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    boltzmann<3, float, random::gpu::mrg32k3a>::luaopen(L);
    boltzmann<2, float, random::gpu::mrg32k3a>::luaopen(L);
    boltzmann<3, float, random::gpu::rand48>::luaopen(L);
    boltzmann<2, float, random::gpu::rand48>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    boltzmann<3, dsfloat, random::gpu::mrg32k3a>::luaopen(L);
    boltzmann<2, dsfloat, random::gpu::mrg32k3a>::luaopen(L);
    boltzmann<3, dsfloat, random::gpu::rand48>::luaopen(L);
    boltzmann<2, dsfloat, random::gpu::rand48>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class boltzmann<3, float, random::gpu::mrg32k3a>;
template class boltzmann<2, float, random::gpu::mrg32k3a>;
template class boltzmann<3, float, random::gpu::rand48>;
template class boltzmann<2, float, random::gpu::rand48>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class boltzmann<3, dsfloat, random::gpu::mrg32k3a>;
template class boltzmann<2, dsfloat, random::gpu::mrg32k3a>;
template class boltzmann<3, dsfloat, random::gpu::rand48>;
template class boltzmann<2, dsfloat, random::gpu::rand48>;
#endif

} // namespace velocities
} // namespace gpu
} // namespace mdsim
} // namespace halmd
