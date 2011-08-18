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

#ifndef HALMD_MDSIM_GPU_INTEGRATORS_VERLET_NVT_HOOVER_HPP
#define HALMD_MDSIM_GPU_INTEGRATORS_VERLET_NVT_HOOVER_HPP

#include <boost/make_shared.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <lua.hpp>

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/integrators/verlet_nvt_hoover_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/integrators/nvt.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type>
class verlet_nvt_hoover
  : public mdsim::integrators::nvt<dimension>
{
public:
    typedef mdsim::integrators::nvt<dimension> _Base;
    typedef gpu::particle<dimension, float> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef logger logger_type;
    typedef fixed_vector<float_type, 2> chain_type;

    typedef typename particle_type::vector_type vector_type;
    typedef typename boost::mpl::if_<
        boost::is_same<float_type, double>, dsfloat, float_type
    >::type gpu_float_type;

    static char const* module_name() { return "verlet_nvt_hoover"; }

    static void luaopen(lua_State* L);

    verlet_nvt_hoover(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type const> box
      , float_type timestep
      , float_type temperature
      , float_type resonance_frequency
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    virtual void integrate();
    virtual void finalize();
    virtual void timestep(double timestep);
    virtual void temperature(double temperature);
    virtual void set_mass(chain_type const& mass);

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
    virtual chain_type const& mass() const
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
    typedef verlet_nvt_hoover_wrapper<dimension, gpu_float_type> wrapper_type;
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type integrate;
        accumulator_type finalize;
        accumulator_type propagate;
        accumulator_type rescale;
    };

    /** propagate chain of Nosé-Hoover variables */
    float_type propagate_chain();

    /** module dependencies */
    boost::shared_ptr<particle_type> particle_;
    boost::shared_ptr<box_type const> box_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;

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

    /** functor to compute actual value of total kinetic energy (multiplied by 2) */
    algorithm::gpu::reduce<
        algorithm::gpu::sum_                    // reduce_transform
      , fixed_vector<float, dimension>          // input_type
      , float4                                  // coalesced_input_type
      , dsfloat                                 // output_type
      , dsfloat                                 // coalesced_output_type
      , float_type                              // host_output_type
      , algorithm::gpu::square_                 // input_transform
    > compute_en_kin_2_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATORS_VERLET_NVT_HOOVER_HPP */
