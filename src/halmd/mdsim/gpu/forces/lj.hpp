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

#ifndef HALMD_MDSIM_GPU_FORCES_LJ_HPP
#define HALMD_MDSIM_GPU_FORCES_LJ_HPP

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/options.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{

template <int dimension, typename float_type>
class lj
  : public mdsim::gpu::force<dimension, float_type>
{
public:
    typedef mdsim::gpu::force<dimension, float_type> _Base;
    typedef typename _Base::matrix_type matrix_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_stress_tensor_type gpu_stress_tensor_type;
    typedef gpu::particle<dimension, float> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef utility::profiler profiler_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    lj(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , boost::array<float, 3> const& cutoff
      , boost::array<float, 3> const& epsilon
      , boost::array<float, 3> const& sigma
    );
    void register_runtimes(profiler_type& profiler);
    virtual void compute();

    //! returns potential cutoff distance
    virtual matrix_type const& cutoff()
    {
        return r_cut_;
    }

    //! returns potential energies of particles
    virtual cuda::vector<float> const& potential_energy()
    {
        return g_en_pot_;
    }

    /** potential part of stress tensors of particles */
    virtual cuda::vector<gpu_stress_tensor_type> const& potential_stress()
    {
        return g_stress_pot_;
    }

    // module runtime accumulator descriptions
    HALMD_PROFILE_TAG( compute_, "Lennard-Jones forces" );

private:
    /** potential well depths in MD units */
    matrix_type epsilon_;
    /** pair separation in MD units */
    matrix_type sigma_;
    /** cutoff length in units of sigma */
    matrix_type r_cut_sigma_;
    /** cutoff length in MD units */
    matrix_type r_cut_;
    /** square of cutoff length */
    matrix_type rr_cut_;
    /** square of pair separation */
    matrix_type sigma2_;
    /** potential energy at cutoff length in MD units */
    matrix_type en_cut_;
    /** Lennard-Jones potential parameters */
    cuda::vector<float4> g_ljparam_;
    /** potential energy for each particle */
    cuda::vector<float> g_en_pot_;
    /** potential part of stress tensor for each particle */
    cuda::vector<gpu_stress_tensor_type> g_stress_pot_;

    boost::fusion::map<
        boost::fusion::pair<compute_, accumulator<double> >
    > runtime_;
};

}}} // namespace mdsim::gpu::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_LJ_HPP */
