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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_HARD_CORE_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_HARD_CORE_HPP

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/pair/adapters/hard_core_kernel.hpp>
#include <halmd/utility/matrix_shape.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace adapters {

/**
 * define hard_core potential adapter
 */
template <typename potential_type>
class hard_core : public potential_type
{
public:
    typedef typename potential_type::float_type float_type;
    typedef typename potential_type::gpu_potential_type parent_potential;
    typedef hard_core_kernel::hard_core<parent_potential> gpu_potential_type;
    typedef typename potential_type::matrix_type matrix_type;

    template<typename... Args>
    hard_core(matrix_type const& core, Args&&... args)
            : potential_type (std::forward<Args>(args)...)
            , r_core_sigma_(check_shape(core, this->sigma()))
            , r_core_(element_prod(core, this->sigma()))
            , g_param_(this->size1() * this->size2())
    {
        LOG("core radius r_core/σ = " << r_core_sigma_);

        cuda::host::vector<float> param(g_param_.size());
        for (size_t i = 0; i < param.size(); ++i) {
            param[i] = r_core_.data()[i];
        }

        cuda::copy(param, g_param_);
    }

    /** bind textures before kernel invocation */
    void bind_textures() const
    {
        hard_core_wrapper<parent_potential>::param.bind(g_param_);
        potential_type::bind_textures();
    }

    matrix_type const& r_core_sigma() const
    {
        return r_core_sigma_;
    }

    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        float_type r = sqrtf(rr);
        float_type r_s = r - r_core_(a,b);
        float_type f_abs, en_pot;
        tie(f_abs, en_pot) = potential_type::operator()(r_s*r_s, a, b);
        f_abs *= r_s / r;
        return make_tuple(f_abs, en_pot);
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

                                                class_<hard_core, potential_type, std::shared_ptr<hard_core> >()
                                                    .property("r_core_sigma", &hard_core::r_core_sigma)
                                              , def("hard_core", &std::make_shared<hard_core
                                                                                 , matrix_type const&
                                                                                 , potential_type const&>)
                                        ]
                                ]
                        ]
                ]
        ];
    }
private:
    /** core radius in units of sigma */
    matrix_type r_core_sigma_;
    /** core radius in MD units */
    matrix_type r_core_;
    /** adapter parameters at CUDA device */
    cuda::vector<float> g_param_;
};

} // namespace adapters
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_HARD_CORE_HPP */
