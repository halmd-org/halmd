/*
 * Copyright © 2013      Nicolas Höft
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

#ifndef HALMD_MDSIM_GPU_FORCES_TABLULATED_EXTERNAL_HPP
#define HALMD_MDSIM_GPU_FORCES_TABLULATED_EXTERNAL_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/forces/tabulated_external_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/mdsim/forces/interpolation/linear.hpp>

#include <lua.hpp>

#include <memory>
#include <tuple>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * class template for modules implementing force based on pretabulated values
 */
template <int dimension, typename float_type, typename force_interpolation_type>
class tabulated_external
{
public:
    typedef float_type coefficient_value_type;
    typedef cuda::vector<coefficient_value_type> coefficient_array_type;

    typedef particle<dimension, float_type> particle_type;
    typedef box<dimension> box_type;
    typedef typename particle_type::position_type position_type;
    typedef typename mdsim::forces::interpolation::linear<dimension, float_type> virial_interpolation_type;
    typedef logger logger_type;

    tabulated_external(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<force_interpolation_type> force_interpolation
      , std::shared_ptr<virial_interpolation_type> virial_interpolation
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    /**
     * Test if the cache is up-to-date and if not, inform the particle
     * module about it (usually done at on_prepend_force())
     */
    void check_cache();

    /**
     * Apply the force to the particles
     */
    void apply();

    /**
     * Return the const reference of the interpolation coefficients
     */
    coefficient_array_type const& coefficients() const
    {
        return force_coefficients_;
    }

    /**
     * Return the reference of the interpolation coefficients
     */
    coefficient_array_type& coefficients()
    {
        return force_coefficients_;
    }

    /**
     * Return the const reference of the interpolation coefficients
     */
    coefficient_array_type const& virial_coefficients() const
    {
        return virial_coefficients_;
    }

    /**
     * Return the reference of the interpolation coefficients
     */
    coefficient_array_type& virial_coefficients()
    {
        return virial_coefficients_;
    }

    /**
     * Return total number of needed coeffcients for interpolation
     */
    size_t ncoefficients() const
    {
        return force_coefficients_.size();
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::force_type force_type;
    typedef tabulated_external_wrapper<dimension, force_interpolation_type, virial_interpolation_type> gpu_wrapper;
    typedef fixed_vector<unsigned int, dimension> index_type;

    /** compute forces */
    void compute_();
    /** compute forces with auxiliary variables */
    void compute_aux_();

    /** state of system */
    std::shared_ptr<particle_type> particle_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** edges of the precalculated grid */
    std::shared_ptr<position_type const> grid_edges_;
    /** force interpolation functor */
    std::shared_ptr<force_interpolation_type const> force_interpolation_;
    /** virial interpolation functor */
    std::shared_ptr<virial_interpolation_type const> virial_interpolation_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;

    /** coeffcients for the force interpolation scheme */
    coefficient_array_type force_coefficients_;
    /** coeffcients for the virial interpolation scheme */
    coefficient_array_type virial_coefficients_;

    /** cache observer of net force per particle */
    cache<> force_cache_;
    /** cache observer of auxiliary variables */
    cache<> aux_cache_;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type compute;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};


} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_TABLULATED_EXTERNAL_HPP */
