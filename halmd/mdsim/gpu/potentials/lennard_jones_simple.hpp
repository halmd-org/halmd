/*
 * Copyright © 2010-2011 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_LENNARD_JONES_SIMPLE_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_LENNARD_JONES_SIMPLE_HPP

#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/lennard_jones_simple_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {

/**
 * Lennard-Jones potential for a single species (constituting a "simple liquid").
 *
 * The usual LJ units are employed, the only parameter is the potential cutoff.
 */
template <typename float_type>
class lennard_jones_simple
{
private:
    typedef boost::numeric::ublas::scalar_matrix<float_type> scalar_matrix_type;

public:
    typedef lennard_jones_simple_kernel::lennard_jones_simple gpu_potential_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef logger logger_type;

    lennard_jones_simple(
        unsigned int ntype1
      , unsigned int ntype2
      , float_type cutoff
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    void bind_textures() const {}

    matrix_type r_cut() const
    {
        return scalar_matrix_type(ntype1_, ntype2_, r_cut_);
    }

    matrix_type epsilon() const
    {
        return scalar_matrix_type(ntype1_, ntype2_, 1);
    }

    matrix_type sigma() const
    {
        return scalar_matrix_type(ntype1_, ntype2_, 1);
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** cutoff length in MD units, r_cut() must return a matrix */
    float_type r_cut_;
    /** square of cutoff length */
    float_type rr_cut_;
    /** potential energy at cutoff length in MD units */
    float_type en_cut_;
    /** number of first particle types */
    unsigned int ntype1_;
    /** number of second particle types */
    unsigned int ntype2_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
};

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_LENNARD_JONES_SIMPLE_HPP */
