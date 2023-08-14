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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_CUSTOM_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_CUSTOM_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/pair/custom_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

/**
 * define custom potential and parameters
 */
template <typename float_type_>
class custom
{
public:
    typedef float_type_ float_type;
    typedef custom_kernel::custom gpu_potential_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;

    custom(
        matrix_type const& sigma    // one parameter must be named sigma and is used as unit of length in the truncations
      , matrix_type const& param2   // FIXME rename param[2-3]
      , matrix_type const& param3
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /** return gpu potential with texture */
    gpu_potential_type get_gpu_potential()
    {
        t_param_ = cuda::texture<float4>(g_param_);
        return gpu_potential_type(t_param_);
    }

    matrix_type const& sigma() const
    {
        return sigma_;
    }

    // FIXME rename param2
    matrix_type const& param2() const
    {
        return param2_;
    }

    // FIXME rename param3
    matrix_type const& param3() const
    {
        return param3_;
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
        return custom_kernel::compute(rr, sigma_(a,b), param2_(a,b), param3_(a,b));
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** interaction range parameter, in MD units */
    matrix_type sigma_;
    /** FIXME second potential parameter, in MD units */
    matrix_type param2_;
    /** FIXME third potential parameter, in MD units */
    matrix_type param3_;
    /** potential parameters at CUDA device */
    cuda::memory::device::vector<float4> g_param_;
    /** array of potential parameters for all combinations of particle types */
    cuda::texture<float4> t_param_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_CUSTOM_HPP */
