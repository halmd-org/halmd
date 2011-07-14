/*
 * Copyright © 2010  Peter Colberg
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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/velocity.hpp>
#include <halmd/mdsim/gpu/velocity_kernel.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension, typename float_type>
velocity<dimension, float_type>::velocity(
    shared_ptr<particle_type> particle
)
  // dependency injection
  : particle(particle)
  // set parameters
  , dim_(particle->dim) // FIXME not used?
{
    cuda::copy(particle->nbox, get_velocity_kernel<dimension>().nbox);
}

/**
 * Rescale magnitude of all velocities by 'factor'
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::rescale(double factor)
{
    LOG_TRACE("rescale particle velocities by a factor of " << factor);
    cuda::configure(dim_.grid, dim_.block);
    get_velocity_kernel<dimension>().rescale(
        particle->g_v
      , particle->dim.threads()
      , factor
    );
}

/**
 * Shift all velocities by 'delta'
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::shift(vector_type const& delta)
{
    LOG_TRACE("shift particle velocities by " << delta);
    cuda::configure(dim_.grid, dim_.block);
    get_velocity_kernel<dimension>().shift(
        particle->g_v
      , particle->dim.threads()
      , delta
    );
}

/**
 * First shift, then rescale all velocities
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::shift_rescale(vector_type const& delta, double factor)
{
    LOG_TRACE("shift particle velocities by " << delta << " and rescale by a factor of " << factor);
    cuda::configure(dim_.grid, dim_.block);
    get_velocity_kernel<dimension>().shift_rescale(
        particle->g_v
      , particle->dim.threads()
      , delta
      , factor
    );
}

template <int dimension, typename float_type>
void velocity<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("velocity_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                class_<velocity, shared_ptr<_Base>, _Base>(class_name.c_str())
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_velocity(lua_State* L)
{
    velocity<3, float>::luaopen(L);
    velocity<2, float>::luaopen(L);
    return 0;
}

template class velocity<3, float>;
template class velocity<2, float>;

} // namespace mdsim
} // namespace gpu
} // namespace halmd
