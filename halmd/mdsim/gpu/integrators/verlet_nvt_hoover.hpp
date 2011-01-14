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

#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/integrators/verlet_nvt_hoover_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/integrators/nvt.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{

template <int dimension, typename float_type>
class verlet_nvt_hoover
  : public mdsim::integrators::nvt<dimension>
{
public:
    typedef mdsim::integrators::nvt<dimension> _Base;
    typedef gpu::particle<dimension, float> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef utility::profiler profiler_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename boost::mpl::if_<
        boost::is_same<float_type, double>, dsfloat, float_type
    >::type gpu_float_type;
    typedef verlet_nvt_hoover_wrapper<dimension, gpu_float_type> wrapper_type;

    static char const* module_name() { return "verlet_nvt_hoover"; }

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    static void luaopen(lua_State* L);

    verlet_nvt_hoover(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , float_type timestep
      , float_type temperature
      , fixed_vector<float_type, 2>  const& mass
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

    //! returns coupling parameters: `mass' of the heat bath variables
    fixed_vector<float_type, 2> const& mass() const
    {
        return mass_xi_;
    }

    /** chain of heat bath variables */
    fixed_vector<float_type, 2> xi;
    fixed_vector<float_type, 2> v_xi;

    // module runtime accumulator descriptions
    HALMD_PROFILING_TAG( integrate_, "first half-step of velocity-Verlet (+ Nosé-Hoover chain)" );
    HALMD_PROFILING_TAG( finalize_, "second half-step of velocity-Verlet (+ Nosé-Hoover chain)" );
    HALMD_PROFILING_TAG( propagate_, "propagate Nosé-Hoover chain" );
    HALMD_PROFILING_TAG( rescale_, "rescale velocities in Nosé-Hoover thermostat" );

private:
    // propagate chain of Nosé-Hoover variables
    float_type propagate_chain();
    // compute actual value of total kinetic energy (multiplied by 2)
    float_type compute_en_kin_2() const;


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

    /** coupling parameters: `mass' of the heat bath variables */
    fixed_vector<float_type, 2> mass_xi_;

    boost::fusion::map<
        boost::fusion::pair<integrate_, accumulator<double> >
      , boost::fusion::pair<finalize_, accumulator<double> >
      , boost::fusion::pair<propagate_, accumulator<double> >
      , boost::fusion::pair<rescale_, accumulator<double> >
    > runtime_;
};

}}} // namespace mdsim::gpu::integrators

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATORS_VERLET_NVT_HOOVER_HPP */
