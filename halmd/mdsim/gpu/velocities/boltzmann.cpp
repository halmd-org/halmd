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
  , std::shared_ptr<logger_type> logger
)
  : _Base(particle, logger)
  // dependency injection
  , particle_(particle)
  , random_(random)
  , logger_(logger)
  // select thread-dependent implementation
  , gaussian_impl_(get_gaussian_impl(random_->rng().dim.threads_per_block()))
  // set parameters
  , temp_(temperature)
  // allocate GPU memory
  , g_vcm_(2 * random_->rng().dim.blocks_per_grid())
  , g_vv_(random_->rng().dim.blocks_per_grid())
{
    LOG("Boltzmann velocity distribution temperature: T = " << temp_);
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
typename boltzmann<dimension, float_type, RandomNumberGenerator>::gaussian_impl_type
boltzmann<dimension, float_type, RandomNumberGenerator>::get_gaussian_impl(int threads)
{
    switch (threads) {
      case 512:
        return wrapper_type::kernel.gaussian_impl_512;
      case 256:
        return wrapper_type::kernel.gaussian_impl_256;
      case 128:
        return wrapper_type::kernel.gaussian_impl_128;
      case 64:
        return wrapper_type::kernel.gaussian_impl_64;
      case 32:
        return wrapper_type::kernel.gaussian_impl_32;
      default:
        throw std::logic_error("invalid gaussian thread count");
    }
}

/**
 * Initialise velocities from Maxwell-Boltzmann distribution
 *
 * The particle velocities need to fullfill two constraints:
 *
 *  1. center of mass velocity shall be zero
 *  2. temperature of the distribution shall equal exactly the given value
 *
 * The above order is chosen as shifting the center of mass velocity
 * means altering the first moment of the velocity distribution, which
 * in consequence affects the second moment, i.e. the temperature.
 *
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void boltzmann<dimension, float_type, RandomNumberGenerator>::set()
{
    cache_proxy<velocity_array_type> velocity = particle_->velocity();

    scoped_timer_type timer(runtime_.set);

    // generate Maxwell-Boltzmann distributed velocities,
    // assuming equal (unit) mass for all particle types
    cuda::configure(
        random_->rng().dim.grid
      , random_->rng().dim.block
      , random_->rng().dim.threads_per_block() * (1 + dimension) * sizeof(dsfloat)
    );
    gaussian_impl_(
        &*velocity->begin()
      , particle_->nparticle()
      , particle_->dim.threads()
      , temp_
      , g_vcm_
      , g_vv_
      , random_->rng().rng()
    );
    cuda::thread::synchronize();

    // set center of mass velocity to zero and
    // rescale velocities to accurate temperature
    cuda::configure(
        particle_->dim.grid
      , particle_->dim.block
      , g_vv_.size() * (1 + dimension) * sizeof(dsfloat)
    );
    wrapper_type::kernel.shift_rescale(
        &*velocity->begin()
      , particle_->nparticle()
      , particle_->dim.threads()
      , temp_
      , g_vcm_
      , g_vv_
      , g_vv_.size()
    );
    cuda::thread::synchronize();

    LOG_DEBUG("assigned Boltzmann-distributed velocities");
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
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_velocities_boltzmann(lua_State* L)
{
    boltzmann<3, float, random::gpu::rand48>::luaopen(L);
    boltzmann<2, float, random::gpu::rand48>::luaopen(L);
    return 0;
}

// explicit instantiation
template class boltzmann<3, float, random::gpu::rand48>;
template class boltzmann<2, float, random::gpu::rand48>;

} // namespace velocities
} // namespace gpu
} // namespace mdsim
} // namespace halmd
