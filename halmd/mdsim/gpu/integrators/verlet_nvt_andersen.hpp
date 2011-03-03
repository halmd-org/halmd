/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_INTEGRATORS_VERLET_NVT_ANDERSEN_HPP
#define HALMD_MDSIM_GPU_INTEGRATORS_VERLET_NVT_ANDERSEN_HPP

#include <lua.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/integrators/nvt.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{

template <int dimension, typename float_type, typename RandomNumberGenerator>
class verlet_nvt_andersen
  : public mdsim::integrators::nvt<dimension>
{
public:
    typedef mdsim::integrators::nvt<dimension> _Base;
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef random::gpu::random<RandomNumberGenerator> random_type;
    typedef utility::profiler profiler_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename random_type::rng_type rng_type;
    typedef verlet_nvt_andersen_wrapper<dimension, rng_type> wrapper_type;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type integrate;
        accumulator_type finalize;
    };

    static char const* module_name() { return "verlet_nvt_andersen"; }

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;
    boost::shared_ptr<random_type> random;

    static void luaopen(lua_State* L);

    verlet_nvt_andersen(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , boost::shared_ptr<random_type> random
      , float_type timestep
      , float_type temperature
      , float_type coll_rate
    );
    void register_runtimes(profiler_type& profiler);
    virtual void integrate();
    virtual void finalize();
    virtual void timestep(double timestep);
    virtual void temperature(double temperature);

    //! returns integration time-step
    virtual double timestep() const
    {
        return timestep_;
    }

    //! returns temperature of heat bath
    virtual double temperature() const
    {
        return temperature_;
    }

    //! returns collision rate with the heat bath
    float_type collision_rate() const
    {
        return coll_rate_;
    }

private:
    /** integration time-step */
    float_type timestep_;
    /** half time-step */
    float_type timestep_half_;
    /** temperature of the heat bath */
    float_type temperature_;
    /** square root of temperature */
    float_type sqrt_temperature_;
    /** collision rate with the heat bath */
    float_type coll_rate_;
    /** probability of a collision with the heat bath during a timestep */
    float_type coll_prob_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

}}} // namespace mdsim::gpu::integrators

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATORS_VERLET_NVT_ANDERSEN_HPP */
