/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP

#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_group.hpp>
#include <halmd/observables/host/samples/sample.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace gpu {

/**
 * phase space sampler abstraction
 *
 * abstract base class defining the interface for the actual sampling implementation
 * this interface is implemented for each data type and for GPU and CPU samples
 */
class phase_space_sampler {
public:
    virtual std::shared_ptr<sample_base> acquire(void) = 0;
    virtual void set(std::shared_ptr<sample_base const> sample) = 0;
    virtual luaponte::object acquire_lua(lua_State* L, std::shared_ptr<phase_space_sampler> self) = 0;
    virtual luaponte::object data_lua(lua_State* L, std::shared_ptr<phase_space_sampler> self) = 0;
    virtual void set_lua(luaponte::object sample) = 0;
};

class phase_space_host_cache;

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
class phase_space
{
public:
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef mdsim::box<dimension> box_type;
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef mdsim::clock clock_type;

    /**
     * Construct phase_space sampler from particle group.
     */
    phase_space(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<particle_group_type> particle_group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<clock_type const> clock
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    void set(std::string const& name, std::shared_ptr<sample_base const> sample) {
        get_sampler(name)->set(sample);
    }

    template<typename sample_type>
    std::shared_ptr<sample_type const> acquire(std::string const& name)
    {
        auto sample = get_sampler(name)->acquire();
        if (sample->type() != (sample_type::gpu ? typeid(gpu_sample<typename sample_type::data_type>) : typeid(typename sample_type::data_type))) {
            throw std::runtime_error("invalid sample data type");
        }
        return std::static_pointer_cast<sample_type const>(sample);
    }

    std::shared_ptr<phase_space_sampler> get_sampler(std::string const& name);

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** particle instance to particle group */
    std::shared_ptr<particle_type> particle_;
    /** particle group */
    std::shared_ptr<particle_group_type> particle_group_;
    /** simulation box */
    std::shared_ptr<box_type const> box_;
    /** simulation clock */
    std::shared_ptr<clock_type const> clock_;
    /** logger instance */
    std::shared_ptr<logger> logger_;

    std::unordered_map<std::string, std::shared_ptr<phase_space_sampler>> samplers_;

    std::map<mdsim::gpu::particle_array*, std::shared_ptr<phase_space_host_cache>> host_cache_;

    typedef halmd::utility::profiler::accumulator_type accumulator_type;
    typedef halmd::utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
        accumulator_type reset;
        accumulator_type set;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP */
