/*
 * Copyright © 2014 Sutapa Roy
 * Copyright © 2014 Felix Höfling
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

#include <boost/numeric/ublas/io.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <stdexcept>
#include <string>

#include <halmd/mdsim/gpu/forces/external.hpp>
#include <halmd/mdsim/gpu/potentials/external/slit.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {

/**
 * Initialise slit potential
 */
template <int dimension, typename float_type>
slit<dimension, float_type>::slit(
    scalar_container_type const& offset
  , vector_container_type const& surface_normal
  , matrix_container_type const& epsilon
  , matrix_container_type const& sigma
  , matrix_container_type const& wetting
  , shared_ptr<logger> logger
)
  // initialise attributes
  : offset_(offset)
  , surface_normal_(surface_normal)
  , epsilon_(epsilon)
  , sigma_(sigma)
  , wetting_(wetting)
  , g_param_geometry_(surface_normal_.size())
  , g_param_potential_(epsilon_.size1() * epsilon_.size2())
  , logger_(logger)
{
    unsigned int nwall = surface_normal_.size();
    unsigned int nspecies = epsilon_.size2();

    // check parameter size
    if (offset_.size() != nwall) {
        throw invalid_argument("geometry parameters have mismatching shapes");
    }

    if (epsilon_.size1() != nwall || sigma_.size1() != nwall || wetting_.size1() != nwall
        || sigma_.size2() != nspecies || wetting_.size2() != nspecies
       ) {
        throw invalid_argument("potential parameters have mismatching shapes");
    }
    
    LOG("number of walls: " << nwall);
    LOG("wall positions: d₀ = " << offset_);
    LOG("surface normals: n = " << surface_normal_);
    LOG("interaction strengths: epsilon = " << epsilon_);
    LOG("interaction ranges: sigma = " << sigma_);
    LOG("wetting paramters: c = " << wetting_);

    // merge geometry parameters in a single array and copy to device
    cuda::host::vector<float4> param_geometry(g_param_geometry_.size());
    for (size_t i = 0; i < nwall; ++i) {
            param_geometry[i] <<= tie(surface_normal_(i), offset_(i));
    }
    cuda::copy(param_geometry, g_param_geometry_);

    // merge potential parameters in a single array and copy to device
    cuda::host::vector<float4> param_potential(g_param_potential_.size());
    for (size_t i = 0; i < nwall; ++i) {
        for (size_t j = 0; j < nspecies; ++j) {
              fixed_vector<float, 4> p;
              p[0] = epsilon_(i, j);
              p[1] = sigma_(i, j);
              p[2] = wetting_(i, j);
              param_potential[j * nwall + i] = p;
        }
    }
    cuda::copy(param_potential, g_param_potential_);
}

template <int dimension, typename float_type>
void slit<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static string class_name("slit_" + to_string(dimension));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("potentials")
                [
                    namespace_("external")
                    [
                        class_<slit, shared_ptr<slit>>(class_name.c_str())
                            .def(constructor<
                                 scalar_container_type const& 
                               , vector_container_type const& 
                               , matrix_container_type const& 
                               , matrix_container_type const& 
                               , matrix_container_type const& 
                               , shared_ptr<logger>
                             >())
                            .property("offset", &slit::offset)
                    ]
                ]
            ]
        ]
    ];
}


HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_external_slit(lua_State* L)
{
    slit<3, float>::luaopen(L);
    slit<2, float>::luaopen(L);
    forces::external<3, float, slit<3, float>>::luaopen(L);
    forces::external<2, float, slit<2, float>>::luaopen(L);
}

// explicit instantiation
template class slit<3, float>;
template class slit<2, float>;

} // namespace external
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
using namespace potentials::external;

template class external<3, float, slit<3, float>>;
template class external<2, float, slit<2, float>>;

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
