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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SHIFTED_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SHIFTED_HPP

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations/shifted_kernel.hpp>
#include <halmd/utility/matrix_shape.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace truncations {

/**
 * define (energy) shifted potential truncation
 */
template <typename potential_type>
class shifted
  : public potential_type
{
public:
    typedef typename potential_type::float_type float_type;
    typedef typename potential_type::gpu_potential_type parent_potential;
    typedef shifted_kernel::shifted<parent_potential> gpu_potential_type;
    typedef typename potential_type::matrix_type matrix_type;

    template<typename... Args>
    shifted(matrix_type const& cutoff, Args&&... args)
            : potential_type (std::forward<Args>(args)...)
            , r_cut_sigma_(check_shape(cutoff, this->sigma()))
            , r_cut_(element_prod(this->sigma(), r_cut_sigma_))
            , rr_cut_(element_prod(r_cut_, r_cut_))
            , en_cut_(this->size1(), this->size2())
            , g_param_(this->size1() * this->size2())
    {

        for (size_t i = 0; i < this->size1(); ++i) {
            for (size_t j = 0; j < this->size2(); ++j) {
                std::tie(std::ignore, en_cut_(i,j)) = potential_type::operator()(rr_cut_(i, j), i, j);
            }
        }

        LOG("potential cutoff length: r_c = " << r_cut_sigma_);
        LOG("potential cutoff energy: U = " << en_cut_);

        cuda::host::vector<float2> param(g_param_.size());
        for (size_t i = 0; i < param.size(); ++i) {
            fixed_vector<float, 2> p;
            p[shifted_kernel::RR_CUT] = rr_cut_.data()[i];
            p[shifted_kernel::EN_CUT] = en_cut_.data()[i];
            param[i] = p;
        }

        cuda::copy(param, g_param_);
    }

    /** bind textures before kernel invocation */
    void bind_textures() const
    {
        shifted_wrapper<parent_potential>::param.bind(g_param_);
        potential_type::bind_textures();
    }

    matrix_type const& r_cut() const
    {
        return r_cut_;
    }

    float_type r_cut(unsigned a, unsigned b) const
    {
        return r_cut_(a, b);
    }

    float_type rr_cut(unsigned a, unsigned b) const
    {
        return rr_cut_(a, b);
    }

    matrix_type const& r_cut_sigma() const
    {
        return r_cut_sigma_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L)
    {
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
                            class_<shifted, potential_type, std::shared_ptr<shifted> >()
                                .property("r_cut", (matrix_type const& (shifted::*)() const) &shifted::r_cut)
                                .property("r_cut_sigma", &shifted::r_cut_sigma)

                          , def("shifted", &std::make_shared<shifted
                              , matrix_type const&
                              , potential_type const&>
                            )
                        ]
                    ]
                ]
            ]
        ];
    }

private:
    /** cutoff length in units of sigma */
    matrix_type r_cut_sigma_;
    /** cutoff length in MD units */
    matrix_type r_cut_;
    /** square of cutoff length */
    matrix_type rr_cut_;
    /** potential energy at cutoff length in MD units */
    matrix_type en_cut_;
    /** adapter parameters at CUDA device */
    cuda::vector<float2> g_param_;
};

} // namespace truncations
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SHIFTED_HPP */
