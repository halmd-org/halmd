/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/io/statevars/writer.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/integrators/nvt.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace integrators
{

/**
 * Nosé-Hoover-Chain thermostat with a chain of length 2
 *
 * The implementation follows
 * Frenkel & Smit: Understanding Molecular Simulation (Appendix E.2)
 */

template <int dimension, typename float_type>
class verlet_nvt_hoover
  : public mdsim::integrators::nvt<dimension>
{
public:
    typedef mdsim::integrators::nvt<dimension> _Base;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef io::statevars::writer<dimension> writer_type;
    typedef utility::profiler profiler_type;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type integrate;
        accumulator_type finalize;
        accumulator_type propagate;
    };

    static char const* module_name() { return "verlet_nvt_hoover"; }

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    static void luaopen(lua_State* L);

    verlet_nvt_hoover(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , float_type timestep
      , float_type temperature
      , float_type resonance_frequency
    );
    void register_observables(writer_type& writer);
    void register_runtimes(profiler_type& profiler);

    virtual void integrate();
    virtual void finalize();
    virtual void timestep(double timestep);
    virtual void temperature(double temperature);
    virtual void mass(fixed_vector<double, 2> const& mass);

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

    //! returns resonance frequency of heat bath
    virtual double resonance_frequency() const
    {
        return resonance_frequency_;
    }

    //! returns coupling parameters: `mass' of the heat bath variables
    fixed_vector<double, 2> const& mass() const
    {
        return mass_xi_;
    }

    //! returns energy per particle of the Nosé-Hoover chain
    double en_nhc() const
    {
        return en_nhc_;
    }

    /** chain of heat bath variables */
    fixed_vector<float_type, 2> xi;
    fixed_vector<float_type, 2> v_xi;

protected:
    // propagate chain of Nosé-Hoover variables
    void propagate_chain();

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
    fixed_vector<float_type, 2> mass_xi_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

}}} // namespace mdsim::host::integrators

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_INTEGRATORS_VERLET_NVT_HOOVER_HPP */
