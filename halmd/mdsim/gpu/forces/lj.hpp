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

#include <boost/numeric/ublas/symmetric.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/forces/lj_kernel.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{

/**
 * define Lennard-Jones potential and parameters
 */
template <typename float_type>
class lj_potential
{
public:
    typedef lj_kernel::lj_potential gpu_potential_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;

    static char const* name() { return "Lennard-Jones"; }
    static char const* module_name() { return "lennard_jones"; }

    static void luaopen(lua_State* L);

    lj_potential(
        unsigned ntype
      , boost::array<float, 3> const& cutoff
      , boost::array<float, 3> const& epsilon
      , boost::array<float, 3> const& sigma
    );

    cuda::texture<float4> const& get_kernel_param() const
    {
        return lj_wrapper::param;
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

}}} // namespace mdsim::gpu::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_LJ_HPP */
