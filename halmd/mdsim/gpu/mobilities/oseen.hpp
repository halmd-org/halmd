/*
 * Copyright © 2011-2012 Michael Kopp
 * Copyright © 2011-2012 Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_MOBILITIES_OSEEN_HPP
#define HALMD_MDSIM_GPU_MOBILITIES_OSEEN_HPP

#include <cuda_wrapper/cuda_wrapper.hpp> // cuda::vector
#include <lua.hpp>
#include <utility>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/mobilities/oseen_kernel.hpp>
#include <halmd/mdsim/mobility.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace mobilities {

/**
 * \brief Class to compute hydrodynamic interactions.
 *
 * In lowest approximation the Oseen tensor is used, in higher approximation
 * Rotne-Prager. Currently only one particle type is supported (only one
 * radius).
 */
template <int dimension, typename float_type>
class oseen
  : public mdsim::mobility<dimension>
{
    typedef mdsim::mobility<dimension> _Base;

public:
    typedef mdsim::box<dimension> box_type;
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename particle_type::gpu_vector_type gpu_vector_type;
    typedef halmd::mdsim::gpu::mobilities::oseen_wrapper<dimension> wrapper_type;

    static void luaopen(lua_State* L);

    oseen(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , float radius
      , float viscosity
      , int order
      , boost::shared_ptr<halmd::logger> logger = boost::make_shared<halmd::logger>()
    );

    // inherited functions
    virtual void compute();
    virtual void compute_velocity();

    /** returns radius */
    float radius() const
    {
        return radius_;
    }

    /** returns dynamic viscosity of fluid */
    float viscosity() const
    {
        return viscosity_;
    }

    /** returns order of integration */
    int order() const
    {
        return order_;
    }

protected:
    /** particle instance */
    boost::shared_ptr<particle_type> particle_;
    /** box instance */
    boost::shared_ptr<box_type> box_;
    /** module logger */
    boost::shared_ptr<logger> logger_;

    /** hydrodynamic radius */
    float radius_;
    /** dynamic viscosity of fluid */
    float viscosity_;
    /** order of accuracy of hydrodynamic interaction in (a/r) */
    int order_;

    typedef utility::profiler::accumulator_type accumulator_type;

    struct runtime
    {
        accumulator_type compute_velocity;
        accumulator_type compute;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mobilities
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_MOBILITIES_OSEEN_HPP */
