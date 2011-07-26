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

#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/forces/lennard_jones_simple_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * define Lennard-Jones potential and parameters
 * for a single species (consituting a "simple liquid")
 *
 * The usual LJ units are employed, the only parameter is
 * the potential cutoff.
 */
template <typename float_type>
class lennard_jones_simple
{
public:
    typedef lennard_jones_simple_kernel::lennard_jones_simple gpu_potential_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef logger logger_type;

    static char const* module_name() { return "lennard_jones_simple"; }

    static void luaopen(lua_State* L);

    lennard_jones_simple(
        float_type cutoff
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    void bind_textures() const {}

    matrix_type r_cut() const
    {
        // construct 1×1 matrix with cutoff
        matrix_type A(1, 1);
        A(0, 0) = r_cut_;
        return A;
    }

    matrix_type epsilon() const
    {
        // construct 1×1 matrix with ε=1
        matrix_type A(1, 1);
        A(0, 0) = 1;
        return A;
    }

    matrix_type sigma() const
    {
        // construct 1×1 matrix with σ=1
        matrix_type A(1, 1);
        A(0, 0) = 1;
        return A;
    }

private:
    /** cutoff length in MD units, r_cut() must return a matrix */
    float_type r_cut_;
    /** square of cutoff length */
    float_type rr_cut_;
    /** potential energy at cutoff length in MD units */
    float_type en_cut_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
};

} // namespace mdsim
} // namespace gpu
} // namespace forces
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_LENNARD_JONES_SIMPLE_HPP */
