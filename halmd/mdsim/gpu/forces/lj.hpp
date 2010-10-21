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

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

// #include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/forces/pair_short_ranged.hpp>
#include <halmd/mdsim/gpu/forces/lj_kernel.hpp>
// #include <halmd/mdsim/gpu/particle.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{

/**
 * define Lennard-Jones potential and parameters
 */
template <int dimension, typename float_type>
class lj_potential
{
public:
    typedef typename mdsim::gpu::force<dimension, float_type>::matrix_type matrix_type;

    lj_potential(
        unsigned ntype
      , boost::array<float, 3> const& cutoff
      , boost::array<float, 3> const& epsilon
      , boost::array<float, 3> const& sigma
    );

    lj_wrapper<dimension> const& get_kernel() const
    {
        return lj_wrapper<dimension>::kernel;
    }

    cuda::vector<float4> const& g_param() const { return g_param_; }

    matrix_type const& r_cut() const { return r_cut_; }

    float_type r_cut(unsigned a, unsigned b) const
    {
        return r_cut_(a, b);
    }

    float_type rr_cut(unsigned a, unsigned b) const
    {
        return rr_cut_(a, b);
    }

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
    /** potential parameters at CUDA device */
    cuda::vector<float4> g_param_;
};

template <int dimension, typename float_type>
class lj
  : public pair_short_ranged<dimension, float_type, lj_potential<dimension, float_type> >
{
public:
    typedef lj_potential<dimension, float_type> potential_type;
    typedef mdsim::gpu::forces::pair_short_ranged<dimension, float_type, potential_type> _Base;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::particle_type particle_type;
    typedef typename _Base::box_type box_type;

    static void luaopen(lua_State* L);

    lj(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , boost::array<float, 3> const& cutoff
      , boost::array<float, 3> const& epsilon
      , boost::array<float, 3> const& sigma
    );
};

}}} // namespace mdsim::gpu::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_LJ_HPP */
