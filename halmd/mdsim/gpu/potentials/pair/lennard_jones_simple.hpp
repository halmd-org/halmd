/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_SIMPLE_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_SIMPLE_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/pair/lennard_jones_simple_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

/**
 * Lennard-Jones potential for a single species (constituting a "simple liquid").
 *
 * The usual LJ units are employed
 */
template <typename float_type_>
class lennard_jones_simple
{
private:
    typedef boost::numeric::ublas::scalar_matrix<float_type_> scalar_matrix_type;

public:
    typedef float_type_ float_type;
    typedef lennard_jones_simple_kernel::lennard_jones_simple gpu_potential_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;

    lennard_jones_simple(
        std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    void bind_textures() const {}

    // FIXME are the following functions actually needed?
    matrix_type epsilon() const
    {
        return scalar_matrix_type(1, 1, 1);
    }

    matrix_type sigma() const
    {
        return scalar_matrix_type(1, 1, 1);
    }

    unsigned int size1() const
    {
        return 1U;
    }

    unsigned int size2() const
    {
        return 1U;
    }

    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        return lennard_jones_simple_kernel::compute(rr);
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_SIMPLE_HPP */
