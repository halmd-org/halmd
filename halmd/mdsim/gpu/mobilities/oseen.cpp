/*
 * Copyright © 2011-2012  Michael Kopp
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
#include <halmd/mdsim/gpu/mobilities/oseen.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace mobilities {

template <int dimension, typename float_type>
oseen<dimension, float_type>::oseen(
    boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<box_type> box
  , float radius
  , float viscosity
  , int order
  , boost::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle(particle)
  , box(box)
  , logger_(logger)
  // set parameters
  , radius_(radius)
  , viscosity_(viscosity)
  , self_mobility_(1/(6*3.14159265358979312*viscosity*radius))
  , order_(order)
{
    // log
    LOG("Particle radii: a = " << radius_);
    LOG("Dynamic viscosity of fluid: eta = " << viscosity_);
    LOG("Order of accurancy of hydrodynamic interaction in (a/r): " << order_);
    if( order_ <= 2 ) LOG( "Using Oseen Tensor for hydrodynamic interaction");
    if( order_ >= 3 ) LOG( "Using Rotne-Prager Tensor for hydrodynamic interaction");
#ifdef USE_OSEEN_DSFUN
    LOG("Accumulating velocities in double-single precision");
#else
    LOG_WARNING("Accumulating velocities in single precision");
#endif
}

/**
 * compute velocities from forces
 */
template <int dimension, typename float_type>
void oseen<dimension, float_type>::compute_velocities()
{
    scoped_timer<timer> timer_(runtime_.compute_velocities); // measure time 'till destruction

    // call kernel
    try {
        // configure cuda parameter (place this immediately before wrapper)
        cuda::configure(
            particle->dim.grid
          , particle->dim.block
          , WARP_SIZE * (sizeof(float4) + sizeof(typename wrapper_type::gpu_vector_type))
          // Allocate shared memory for position (float4) and forces (float4 in 3D, float2 in 2D).
        );
        if (order_ <= 2) // oseen
            wrapper_type::wrapper.compute_velocities_oseen(
                    particle->g_r, particle->g_f, particle->g_v, particle->nbox, static_cast<vector_type>(box->length()), radius_, self_mobility_
                    );
        else // rotne
            wrapper_type::wrapper.compute_velocities_rotne(
                    particle->g_r, particle->g_f, particle->g_v, particle->nbox, static_cast<vector_type>(box->length()), radius_, self_mobility_
                    );
        cuda::thread::synchronize();

        // Since parameters such as radius, self mobility etc. are constants, one
        // could also store them on the GPU (cuda::symbol<>) instead of passing them
        // as arguments.

    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream computation of velocities via oseen on GPU");
        throw;
    }
}

//! Compute mobility matrix -- not implemented yet
template <int dimension, typename float_type>
void oseen<dimension, float_type>::compute()
{
    scoped_timer<timer> timer_(runtime_.compute); // measure time 'till destruction
}

// Wrapper to connect set with slots.
template <typename mobility_type>
typename signal<void ()>::slot_function_type
wrap_compute(shared_ptr<mobility_type> self)
{
    return bind(&mobility_type::compute, self);
}

template <typename mobility_type>
typename signal<void ()>::slot_function_type
wrap_compute_velocities(shared_ptr<mobility_type> self)
{
    return bind(&mobility_type::compute_velocities, self);
}

template <int dimension, typename float_type>
void oseen<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("mobilities")
                [
                    class_<oseen, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                            boost::shared_ptr<particle_type>
                          , boost::shared_ptr<box_type>
                          , float
                          , float
                          , int
                          , shared_ptr<logger_type>
                         >())
                        .property("compute", &wrap_compute<oseen>)
                        .property("compute_velocities", &wrap_compute_velocities<oseen>)
                        .property("radius", &oseen::radius)
                        .property("viscosity", &oseen::viscosity)
                        .property("order", &oseen::order)
                        .property("self_mobility", &oseen::self_mobility)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("compute", &runtime::compute)
                                .def_readonly("compute_velocities", &runtime::compute_velocities)
                        ]
                        .def_readonly("runtime", &oseen::runtime_)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_mobilities_oseen(lua_State* L)
{
    oseen<3, float>::luaopen(L);
    oseen<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class oseen<3, float>;
template class oseen<2, float>;

} // namespace mobilities
} // namespace gpu
} // namespace mdsim
} // namespace halmd
