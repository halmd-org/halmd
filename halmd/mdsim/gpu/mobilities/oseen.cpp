/*
 * Copyright © 2011-2012 Michael Kopp
 * Copyright © 2011-2014 Felix Höfling
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
  , boost::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
  // set parameters
  , radius_(radius)
  , viscosity_(viscosity)
  , order_(order)
{
    LOG("Particle radii: a = " << radius_);
    LOG("Dynamic viscosity of fluid: η = " << viscosity_);
    LOG("Order of accurancy of hydrodynamic interaction in (a/r): " << order_);
    if( order_ <= 2 ) LOG( "using Oseen tensor for hydrodynamic interaction");
    if( order_ >= 3 ) LOG( "using Rotne-Prager tensor for hydrodynamic interaction");
#ifdef USE_OSEEN_DSFUN
    LOG("accumulating velocities in double-single precision");
#else
    LOG_WARNING("accumulating velocities in single precision");
#endif
}

/**
 * compute velocities from forces
 */
template <int dimension, typename float_type>
void oseen<dimension, float_type>::compute_velocity()
{
    scoped_timer<timer> timer_(runtime_.compute_velocity);

    // Stokes mobility: μ = 1 / 6 π η a, a = radius, η = shear viscosity
    float_type self_mobility = 1 / (6 * 3.14159265358979312 * viscosity_ * radius_);

    try {
        cuda::configure(
            particle_->dim.grid
          , particle_->dim.block
          , WARP_SIZE * (sizeof(float4) + sizeof(typename wrapper_type::gpu_vector_type))
          // Allocate shared memory for position (float4) and forces (float4 in 3D, float2 in 2D).
        );
        if (order_ <= 2) // Oseen tensor
        {
            wrapper_type::wrapper.compute_velocity_oseen(
                particle_->g_r, particle_->g_f, particle_->g_v
              , particle_->nbox, static_cast<vector_type>(box_->length()), radius_, self_mobility
            );
        }
        else {
            // Rotne-Prager tensor
            wrapper_type::wrapper.compute_velocity_rotne(
                particle_->g_r, particle_->g_f, particle_->g_v
              , particle_->nbox, static_cast<vector_type>(box_->length()), radius_, self_mobility
            );
        }
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream computation of velocities via mobility tensor on GPU");
        throw;
    }
}

//! Compute mobility matrix -- not yet implemented
template <int dimension, typename float_type>
void oseen<dimension, float_type>::compute()
{
    scoped_timer<timer> timer_(runtime_.compute);
    LOG_ERROR("computation and storage of mobility matrix not yet implemented");
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
wrap_compute_velocity(shared_ptr<mobility_type> self)
{
    return bind(&mobility_type::compute_velocity, self);
}

template <int dimension, typename float_type>
void oseen<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("oseen_" + lexical_cast<string>(dimension) + "_");

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
                          , shared_ptr<logger>
                         >())
                        .property("compute", &wrap_compute<oseen>)
                        .property("compute_velocity", &wrap_compute_velocity<oseen>)
                        .property("radius", &oseen::radius)
                        .property("viscosity", &oseen::viscosity)
                        .property("order", &oseen::order)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("compute", &runtime::compute)
                                .def_readonly("compute_velocity", &runtime::compute_velocity)
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
