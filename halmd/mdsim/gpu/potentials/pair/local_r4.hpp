/*
 * Copyright © 2016 Daniel Kirchner
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_LOCAL_R4_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_LOCAL_R4_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/pair/local_r4_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

/**
 * define Lennard-Jones potential and parameters
 */
template <typename potential_type>
class local_r4 : public potential_type
{
public:
    typedef typename potential_type::float_type float_type;
    typedef typename potential_type::gpu_potential_type parent_potential;
    typedef local_r4_kernel::local_r4<parent_potential> gpu_potential_type;
    typedef typename potential_type::matrix_type matrix_type;

    template<typename... Args>
    local_r4(float_type h, Args&&... args)
            : potential_type (std::forward<Args>(args)...), rri_smooth_(std::pow(h, -2)) {
    }

    /** bind textures before kernel invocation */
    void bind_textures() const
    {
        cuda::copy(rri_smooth_, local_r4_wrapper<parent_potential>::rri_smooth);
        potential_type::bind_textures();
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L) {
        using namespace luaponte;
        module(L, "libhalmd")
        [
                namespace_("mdsim")
                [
                        namespace_("gpu")
                        [
                                namespace_("potentials")
                                [
                                        namespace_("pair")
                                        [

                                                class_<local_r4, potential_type, std::shared_ptr<local_r4> >()
                                              , def("local_r4", &std::make_shared<local_r4, float_type, potential_type const&>)
                                        ]
                                ]
                        ]
                ]
        ];
    }
private:
    float_type rri_smooth_;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_LOCAL_R4_HPP */
