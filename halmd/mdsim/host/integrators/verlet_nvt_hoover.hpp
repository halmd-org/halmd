/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_INTEGRATORS_VERLET_NVT_HOOVER_HPP
#define HALMD_MDSIM_HOST_INTEGRATORS_VERLET_NVT_HOOVER_HPP

#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

/**
 * Nosé-Hoover-Chain thermostat with a chain of length 2
 *
 * The implementation follows
 * Frenkel & Smit: Understanding Molecular Simulation (Appendix E.2)
 */

template <int dimension, typename float_type>
class verlet_nvt_hoover
{
public:
    typedef particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef box<dimension> box_type;
    typedef logger logger_type;
    typedef fixed_vector<float_type, 2> chain_type;

    static void luaopen(lua_State* L);

    verlet_nvt_hoover(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
      , float_type timestep
      , float_type temperature
      , float_type resonance_frequency
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    void integrate();
    void finalize();
    void set_timestep(double timestep);
    void set_temperature(double temperature);
    void set_mass(chain_type const& mass);

    //! returns integration time-step
    double timestep() const
    {
        return timestep_;
    }

    //! returns temperature of heat bath
    double temperature() const
    {
        return temperature_;
    }

    //! returns resonance frequency of heat bath
    double resonance_frequency() const
    {
        return resonance_frequency_;
    }

    //! returns coupling parameters: `mass' of the heat bath variables
    chain_type const& mass() const
    {
        return mass_xi_;
    }

    //! returns energy per particle of the Nosé-Hoover chain
    double en_nhc() const
    {
        return en_nhc_;
    }

    /**
     * chain of heat bath variables
     *
     * In analogy with the particle positions and velocities, these variables are accessible to the public.
     */
    chain_type xi;
    chain_type v_xi;

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::image_array_type image_array_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_type::force_array_type force_array_type;
    typedef typename particle_type::mass_array_type mass_array_type;
    typedef typename particle_type::size_type size_type;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type integrate;
        accumulator_type finalize;
        accumulator_type propagate;
        accumulator_type rescale; //< for compatibility with GPU backend
    };

    // propagate chain of Nosé-Hoover variables
    void propagate_chain();

    /** system state */
    std::shared_ptr<particle_type> particle_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;

    /** integration time-step */
    float_type timestep_;
    /** fractions of the time-step */
    float_type timestep_half_;
    float_type timestep_4_;
    float_type timestep_8_;
    /** temperature of the heat bath */
    float_type temperature_;
    /** target value for twice the total kinetic energy */
    float_type en_kin_target_2_;
    /** energy of chain variables per particle */
    float_type en_nhc_;

    /** resonance frequency of heat bath, determines coupling parameters below */
    float_type resonance_frequency_;
    /** coupling parameters: `mass' of the heat bath variables */
    chain_type mass_xi_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace integrators
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_INTEGRATORS_VERLET_NVT_HOOVER_HPP */
