/*
 * Copyright © 2010-2013 Felix Höfling
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

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/forces/trunc/local_r4.hpp>
#include <halmd/mdsim/gpu/forces/pair_full.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/potentials/pair/lennard_jones_simple.hpp>
#include <halmd/mdsim/gpu/potentials/pair/lennard_jones_simple_kernel.hpp>
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
    float_type cutoff
  , std::shared_ptr<logger> logger
)
  // initialise members
  : r_cut_(cutoff)
  , rr_cut_(cutoff * cutoff)
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    float_type rri_cut = 1 / rr_cut_;
    float_type r6i_cut = rri_cut * rri_cut * rri_cut;
    en_cut_ = 4 * r6i_cut * (r6i_cut - 1);

    LOG("using optimised version for a single species with ε = 1, σ = 1");
    LOG("potential cutoff length: r_c = " << r_cut_);
    LOG("potential cutoff energy: U = " << en_cut_);

    cuda::copy(rr_cut_, lennard_jones_simple_wrapper::rr_cut);
    cuda::copy(en_cut_, lennard_jones_simple_wrapper::en_cut);
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
                            .def(constructor<
                                float_type
                              , std::shared_ptr<logger>
                            >())
                            // provide Lua interface coherent with lennard_jones
                            .property("r_cut", &lennard_jones_simple::r_cut)
                            .property("r_cut_sigma", &lennard_jones_simple::r_cut) // note σ=1
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
    forces::pair_trunc<3, float, lennard_jones_simple<float> >::luaopen(L);
    forces::pair_trunc<2, float, lennard_jones_simple<float> >::luaopen(L);
    forces::pair_trunc<3, float, lennard_jones_simple<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
    forces::pair_trunc<2, float, lennard_jones_simple<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class lennard_jones_simple<float>;

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
template class pair_full<3, float, potentials::pair::lennard_jones_simple<float> >;
template class pair_full<2, float, potentials::pair::lennard_jones_simple<float> >;
template class pair_trunc<3, float, potentials::pair::lennard_jones_simple<float> >;
template class pair_trunc<2, float, potentials::pair::lennard_jones_simple<float> >;
template class pair_trunc<3, float, potentials::pair::lennard_jones_simple<float>, mdsim::forces::trunc::local_r4<float> >;
template class pair_trunc<2, float, potentials::pair::lennard_jones_simple<float>, mdsim::forces::trunc::local_r4<float> >;

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
