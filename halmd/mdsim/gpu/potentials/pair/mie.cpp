/*
 * Copyright © 2008-2013 Felix Höfling
 * Copyright © 2008-2010 Peter Colberg
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

#include <boost/numeric/ublas/io.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <cmath>
#include <stdexcept>
#include <string>

#include <halmd/mdsim/gpu/forces/pair_full.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/potentials/pair/mie.hpp>
#include <halmd/mdsim/gpu/potentials/pair/mie_kernel.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations/truncations.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost::numeric::ublas;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

/**
 * Initialise Lennard-Jones potential parameters
 */
template <typename float_type>
mie<float_type>::mie(
    matrix_type const& epsilon
  , matrix_type const& sigma
  , uint_matrix_type const& index_m
  , uint_matrix_type const& index_n
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , index_m_(check_shape(index_m, epsilon))
  , index_n_(check_shape(index_n, epsilon))
  , sigma2_(element_prod(sigma_, sigma_))
  , g_param_(size1() * size2())
  , t_param_(g_param_)
  , logger_(logger)
{
    LOG("potential well depths: ε = " << epsilon_);
    LOG("interaction range: σ = " << sigma_);
    LOG("index of repulsion: m = " << index_m_);
    LOG("index of attraction: n = " << index_n_);

    // check conditions on power low indices (after logging output)
    for (unsigned i = 0; i < index_m_.size1(); ++i) {
        for (unsigned j = 0; j < index_m_.size2(); ++j) {
            // indices must be even
            if (index_m_(i, j) & 1 || index_n_(i, j) & 1) {
                throw std::logic_error("power law indices of potential must be even");
            }
            if (index_m_(i, j) <= index_n_(i, j)) {
                throw std::logic_error("repulsive part of potential must be stronger than attraction");
            }
        }
    }

    cuda::memory::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p;
        p[mie_kernel::EPSILON] = epsilon_.data()[i];
        p[mie_kernel::SIGMA2] = sigma2_.data()[i];
        p[mie_kernel::INDEX_M_2] = index_m_.data()[i] / 2;
        p[mie_kernel::INDEX_N_2] = index_n_.data()[i] / 2;
        param[i] = p;
    }
    cuda::copy(param.begin(), param.end(), g_param_.begin());
}

template <typename float_type>
void mie<float_type>::luaopen(lua_State* L)
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
                        class_<mie, std::shared_ptr<mie> >("mie")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , uint_matrix_type const&
                              , uint_matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("epsilon", &mie::epsilon)
                            .property("sigma", &mie::sigma)
                            .property("index_m", &mie::index_m)
                            .property("index_n", &mie::index_n)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_pair_mie(lua_State* L)
{
    mie<float>::luaopen(L);
#ifdef USE_GPU_SINGLE_PRECISION
    forces::pair_full<3, float, mie<float> >::luaopen(L);
    forces::pair_full<2, float, mie<float> >::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    forces::pair_full<3, dsfloat, mie<float> >::luaopen(L);
    forces::pair_full<2, dsfloat, mie<float> >::luaopen(L);
#endif
    truncations::truncations_luaopen<mie<float> >(L);
    return 0;
}

// explicit instantiation
template class mie<float>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(mie<float>)

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifdef USE_GPU_SINGLE_PRECISION
template class pair_full<3, float, potentials::pair::mie<float> >;
template class pair_full<2, float, potentials::pair::mie<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::mie<float>)
#endif

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class pair_full<3, dsfloat, potentials::pair::mie<float> >;
template class pair_full<2, dsfloat, potentials::pair::mie<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(dsfloat, potentials::pair::mie<float>)
#endif

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
