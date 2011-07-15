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

#ifndef HALMD_MDSIM_GPU_FORCES_LENNARD_JONES_SIMPLE_HPP
#define HALMD_MDSIM_GPU_FORCES_LENNARD_JONES_SIMPLE_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/forces/lennard_jones_simple_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * define Lennard-Jones potential and parameters
 * for a single species
 *
 * The usual LJ units are employed, only parameter is
 * the potential cutoff.
 */
template <typename float_type>
class lennard_jones_simple
{
public:
    typedef lennard_jones_simple_kernel::lennard_jones_simple gpu_potential_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;

    static char const* module_name() { return "lennard_jones_simple"; }

    static void luaopen(lua_State* L);

    lennard_jones_simple(float_type cutoff);

    void bind_textures() const {}

    matrix_type const& r_cut() const
    {
        return r_cut_;
    }

    matrix_type const& epsilon() const
    {
        return epsilon_;
    }

    matrix_type const& sigma() const
    {
        return sigma_;
    }

private:
    /** cutoff length in MD units, r_cut() must return a matrix */
    matrix_type r_cut_;
    /** square of cutoff length */
    float_type rr_cut_;
    /** potential energy at cutoff length in MD units */
    float_type en_cut_;
    /** potential well depths in MD units, for coherence with lennard_jones only */
    const matrix_type epsilon_;
    /** pair separation in MD units, for coherence with lennard_jones only */
    const matrix_type sigma_;
};

} // namespace mdsim
} // namespace gpu
} // namespace forces
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_LENNARD_JONES_SIMPLE_HPP */
