/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {

template <int dimension, typename float_type, typename RandomNumberGenerator>
boltzmann<dimension, float_type, RandomNumberGenerator>::boltzmann(
    shared_ptr<particle_type> particle
  , shared_ptr<random_type> random
  , double temperature
  , shared_ptr<logger_type> logger
)
  : _Base(particle, logger)
  // dependency injection
  , particle_(particle)
  , random_(random)
  , logger_(logger)
  // select thread-dependent implementation
  , gaussian_impl_(get_gaussian_impl(random_->rng.dim.threads_per_block()))
  // set parameters
  , temp_(temperature)
  // allocate GPU memory
  , g_vcm_(2 * random_->rng.dim.blocks_per_grid())
  , g_vv_(random_->rng.dim.blocks_per_grid())
{
    // copy random number generator parameters to GPU
    cuda::copy(random_->rng.rng(), wrapper_type::kernel.rng);

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
    // generate Maxwell-Boltzmann distributed velocities,
    // assuming equal (unit) mass for all particle types
    cuda::configure(
        random_->rng.dim.grid
      , random_->rng.dim.block
      , random_->rng.dim.threads_per_block() * (1 + dimension) * sizeof(dsfloat)
    );
    gaussian_impl_(
        particle_->g_v
      , particle_->nbox
      , particle_->dim.threads()
      , temp_
      , g_vcm_
      , g_vv_
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
        particle_->g_v
      , particle_->nbox
      , particle_->dim.threads()
      , temp_
      , g_vcm_
      , g_vv_
      , g_vv_.size()
    );
    cuda::thread::synchronize();

#ifdef USE_HILBERT_ORDER
    // make thermostat independent of neighbour list update frequency or skin
//     order_velocities(); boltzmann is not a thermostat!
#endif

    LOG_DEBUG("assigned Boltzmann-distributed velocities");
//    LOG_DEBUG("velocities rescaled by factor " << scale);
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
static char const* module_name_wrapper(boltzmann<dimension, float_type, RandomNumberGenerator> const&)
{
    return boltzmann<dimension, float_type, RandomNumberGenerator>::module_name();
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void boltzmann<dimension, float_type, RandomNumberGenerator>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("velocities")
                [
                    class_<boltzmann, shared_ptr<_Base_Base>, bases<_Base_Base, _Base> >(class_name.c_str())
                        .def(constructor<
                             shared_ptr<particle_type>
                           , shared_ptr<random_type>
                           , double
                           , shared_ptr<logger_type>
                         >())
                        .property("temperature", &boltzmann::temperature)
                        .property("module_name", &module_name_wrapper<dimension, float_type, RandomNumberGenerator>)
                ]
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

} // namespace mdsim
} // namespace gpu
} // namespace velocities
} // namespace halmd
