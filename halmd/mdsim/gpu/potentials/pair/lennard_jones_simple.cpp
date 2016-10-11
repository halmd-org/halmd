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

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/forces/pair_full.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/potentials/pair/lennard_jones_simple.hpp>
#include <halmd/mdsim/gpu/potentials/pair/lennard_jones_simple_kernel.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

/**
 * Initialise Lennard-Jones potential parameters
 */
template <typename float_type>
lennard_jones_simple<float_type>::lennard_jones_simple(
    std::shared_ptr<logger> logger
)
  // initialise members
  : logger_(logger)
{
    LOG("using optimised version for a single species with ε = 1, σ = 1");
}

template <typename float_type>
void lennard_jones_simple<float_type>::luaopen(lua_State* L)
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
                        class_<lennard_jones_simple, std::shared_ptr<lennard_jones_simple> >("lennard_jones_simple")
                            .def(constructor<std::shared_ptr<logger> >())
                            .property("epsilon", &lennard_jones_simple::epsilon)
                            .property("sigma", &lennard_jones_simple::sigma)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_pair_lennard_jones_simple(lua_State* L)
{
    lennard_jones_simple<float>::luaopen(L);
    forces::pair_full<3, float, lennard_jones_simple<float> >::luaopen(L);
    forces::pair_full<2, float, lennard_jones_simple<float> >::luaopen(L);
    truncations_luaopen<float, lennard_jones_simple<float> >(L);
    return 0;
}

// explicit instantiation
template class lennard_jones_simple<float>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(lennard_jones_simple<float>);

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
template class pair_full<3, float, potentials::pair::lennard_jones_simple<float> >;
template class pair_full<2, float, potentials::pair::lennard_jones_simple<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::lennard_jones_simple<float>);

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
