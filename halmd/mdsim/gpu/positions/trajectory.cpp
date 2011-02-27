/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <algorithm>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/positions/trajectory.hpp>
#include <halmd/mdsim/gpu/positions/trajectory_kernel.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace positions
{

using namespace boost;
using namespace std;

template <int dimension, typename float_type>
trajectory<dimension, float_type>::trajectory(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<sample_type> sample
)
  // dependency injection
  : particle(particle)
  , box(box)
  , sample(sample)
{
    try {
        cuda::copy(static_cast<vector_type>(box->length()), trajectory_wrapper<dimension>::kernel.box_length);
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("[trajectory] failed to copy box length to GPU");
        throw;
    }
}

/**
 * set particle positions
 */
template <int dimension, typename float_type>
void trajectory<dimension, float_type>::set()
{
    // assign particle coordinates
    for (size_t j = 0, i = 0; j < particle->ntype; i += particle->ntypes[j], ++j) {
        copy(sample->r[j]->begin(), sample->r[j]->end(), &particle->h_r[i]);
    }

#ifdef USE_VERLET_DSFUN
    // erase particle position vectors (double-single precision)
    cuda::memset(particle->g_r, 0, particle->g_r.capacity());
#endif

    try {
        cuda::copy(particle->h_r, particle->g_r);
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("[trajectory] failed to copy particle positions to GPU");
        throw;
    }

    // shift particle positions to range (-L/2, L/2)
    try {
        cuda::configure(particle->dim.grid, particle->dim.block);
        trajectory_wrapper<dimension>::kernel.reduce_periodic(particle->g_r);
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("[trajectory] failed to reduce particle positions on GPU");
        throw;
    }

    // assign particle image vectors
    cuda::memset(particle->g_image, 0, particle->g_image.capacity());

    LOG("set particle positions from trajectory sample");
}

template <int dimension, typename float_type>
void trajectory<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("trajectory_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("gpu")
                [
                    namespace_("positions")
                    [
                        class_<trajectory, shared_ptr<_Base>, _Base>(class_name.c_str())
                            .def(constructor<
                                 shared_ptr<particle_type>
                               , shared_ptr<box_type>
                               , shared_ptr<sample_type>
                            >())
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
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &trajectory<3, float>::luaopen
    ]
    [
        &trajectory<2, float>::luaopen
    ];
}

} // namespace

// explicit instantiation
template class trajectory<3, float>;
template class trajectory<2, float>;

}}} // namespace mdsim::gpu::positions

} // namespace halmd
