/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#ifndef HALMD_OBSERVABLES_HOST_DENSITY_MODE_HPP
#define HALMD_OBSERVABLES_HOST_DENSITY_MODE_HPP

#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_group.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/owner_equal.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/raw_array.hpp>

namespace halmd {
namespace observables {
namespace host {

/**
 * Compute Fourier modes of the particle density.
 *
 * @f$ \rho_{\vec q} = \sum_{i=1}^N \exp(\textrm{i}\vec q \cdot \vec r_i) @f$
 *
 * The result is stored and returned within a std::shared_ptr, allowing
 * efficient copying, e.g., in dynamics::blocking_scheme.  Further, the result
 * may be tracked by std::weak_ptr providing a similar functionality as
 * halmd::cache
 */
template <int dimension, typename float_type>
class density_mode
{
public:
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::particle_group particle_group_type;
    typedef observables::utility::wavevector<dimension> wavevector_type;
    typedef raw_array<fixed_vector<double, 2>> result_type;

    density_mode(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<particle_group_type> particle_group
      , std::shared_ptr<wavevector_type const> wavevector
      , std::shared_ptr<logger> logger = std::make_shared<logger>("density_mode")
    );

    /** Compute density modes from particle group.
     *
     * The result is re-computed only if the particle positions have been
     * modified. In this case, the data array managed by std::shared_ptr is
     * re-allocated before.
     */
    std::shared_ptr<result_type const> acquire();

    /**
     * Functor wrapping acquire() for class instance stored within std::shared_ptr
     */
    static std::function<std::shared_ptr<result_type const> ()>
    acquisitor(std::shared_ptr<density_mode> self)
    {
        return [=]() {
           return self->acquire();
        };
    }

    /**
     * Return wavevector instance passed to constructor
     */
    std::shared_ptr<wavevector_type const> wavevector()
    {
        return wavevector_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef fixed_vector<float_type, dimension> vector_type;

    /** system state */
    std::shared_ptr<particle_type const> particle_;
    /** particle group */
    std::shared_ptr<particle_group_type> particle_group_;
    /** wavevector list */
    std::shared_ptr<wavevector_type const> wavevector_;
    /** logger instance */
    std::shared_ptr<logger> logger_;

    /** result for the density modes */
    std::shared_ptr<result_type> result_;
    /** cache observer for particle positions */
    cache<> position_cache_;
    /** cache observer for particle group */
    cache<> group_cache_;

    typedef halmd::utility::profiler::accumulator_type accumulator_type;
    typedef halmd::utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace observables
} // namespace host
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_DENSITY_MODE_HPP */
