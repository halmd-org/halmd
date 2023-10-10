/*
 * Copyright © 2023 Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_YUKAWA_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_YUKAWA_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/pair/yukawa_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

/**
 * define Yukawa potential and parameters
 */
template <typename float_type_>
class yukawa
{
public:
    typedef float_type_ float_type;
    typedef yukawa_kernel::yukawa gpu_potential_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;

    yukawa(
        matrix_type const& amplitude
      , matrix_type const& sigma
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /** return gpu potential with texture */
    gpu_potential_type get_gpu_potential()
    {
        t_param_ = cuda::texture<float2>(g_param_);
        return gpu_potential_type(t_param_);
    }

    matrix_type const& amplitude() const
    {
        return amplitude_;
    }

    matrix_type const& sigma() const
    {
        return sigma_;
    }

    unsigned int size1() const
    {
        return sigma_.size1();
    }

    unsigned int size2() const
    {
        return sigma_.size2();
    }

    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        return yukawa_kernel::compute(rr, amplitude_(a,b), sigma_(a,b));
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** interaction strength, in MD units of energy x length*/
    matrix_type amplitude_;
    /** screening length (interaction range parameter), in MD length units */
    matrix_type sigma_;
    /** potential parameters at CUDA device */
    cuda::memory::device::vector<float2> g_param_;
    /** array of potential parameters for all combinations of particle types */
    cuda::texture<float2> t_param_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_YUKAWA_HPP */
