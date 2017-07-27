/*
 * Copyright © 2016       Daniel Kirchner
 * Copyright © 2008-2012  Felix Höfling
 * Copyright © 2008-2012  Peter Colberg
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
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_group.hpp>
#include <halmd/observables/host/samples/sample.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace gpu {

/**
 * phase space sampler abstraction (gpu)
 *
 * abstract base class defining the interface for the actual sampling implementation
 * this interface is implemented for each data type
 */
class phase_space_sampler_gpu
{
public:
    virtual std::shared_ptr<sample_base> acquire(void) = 0;
    virtual void set(std::shared_ptr<sample_base const> sample) = 0;
    virtual luaponte::object acquire_lua(lua_State* L, std::shared_ptr<phase_space_sampler_gpu> self) = 0;
    virtual luaponte::object data_lua(lua_State* L, std::shared_ptr<phase_space_sampler_gpu> self) = 0;
    virtual void set_lua(luaponte::object sample) = 0;
};

/**
 * phase space sampler abstraction (host)
 *
 * abstract base class defining the interface for the actual sampling implementation
 * this interface is implemented for each data type
 */
class phase_space_sampler_host
{
public:
    virtual std::shared_ptr<sample_base> acquire(void) = 0;
    virtual void set(std::shared_ptr<sample_base const> sample) = 0;
    virtual luaponte::object acquire_lua(lua_State* L, std::shared_ptr<phase_space_sampler_host> self) = 0;
    virtual luaponte::object data_lua(lua_State* L, std::shared_ptr<phase_space_sampler_host> self) = 0;
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
    typedef fixed_vector<float, dimension> vector_type;

    /**
     * Construct phase_space sampler from particle group.
     */
    phase_space(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<particle_group_type> particle_group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    void set(std::string const& name, std::shared_ptr<sample_base const> sample) {
        if (sample->gpu()) {
            get_sampler_gpu(name)->set(sample);
        } else {
            get_sampler_host(name)->set(sample);
        }
    }

    /**
     * Acquires a sample of the given type (template parameter) of
     * a particle array.
     *
     * @param name      identifier of the particle array to be sampled
     */
    template<typename sample_type>
    std::shared_ptr<sample_type const> acquire(std::string const& name)
    {
        if (sample_type::gpu_sample) {
            auto sample = get_sampler_gpu(name)->acquire();
            if (sample->type() != typeid(typename sample_type::data_type)) {
                throw std::runtime_error("invalid sample data type");
            }
            return std::static_pointer_cast<sample_type const>(sample);
        } else {
            auto sample = get_sampler_host(name)->acquire();
            if (sample->type() != typeid(typename sample_type::data_type)) {
                throw std::runtime_error("invalid sample data type");
            }
            return std::static_pointer_cast<sample_type const>(sample);
        }
    }

    std::shared_ptr<phase_space_sampler_gpu> get_sampler_gpu(std::string const& name);
    std::shared_ptr<phase_space_sampler_host> get_sampler_host(std::string const& name);

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
    /** logger instance */
    std::shared_ptr<logger> logger_;

    /** Associative container mapping identifiers of particle arrays to
        already created gpu sampler implementations */
    std::unordered_map<std::string, std::shared_ptr<phase_space_sampler_gpu>> gpu_samplers_;

    /** Associative container mapping identifiers of particle arrays to
    already created gpu sampler implementations */
    std::unordered_map<std::string, std::shared_ptr<phase_space_sampler_host>> host_samplers_;

    /** Associative container mapping GPU particle arrays to their host cache. */
    std::map<mdsim::gpu::particle_array_gpu_base*, std::shared_ptr<phase_space_host_cache>> host_cache_;

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
