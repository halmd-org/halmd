/*
 * Copyright © 2014-2015 Sutapa Roy
 * Copyright © 2014-2015 Felix Höfling
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
#include <halmd/mdsim/gpu/potentials/external/planar_wall.hpp>
#include <halmd/mdsim/gpu/potentials/external/planar_wall_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {

/**
 * Initialise planar_wall potential
 */
template <int dimension, typename float_type>
planar_wall<dimension, float_type>::planar_wall(
    scalar_container_type const& offset
  , vector_container_type const& surface_normal
  , matrix_container_type const& epsilon
  , matrix_container_type const& sigma
  , matrix_container_type const& wetting
  , matrix_container_type const& cutoff
  , float_type smoothing
  , shared_ptr<logger> logger
)
  // initialise attributes
  : offset_(offset)
  , surface_normal_(surface_normal)
  , epsilon_(epsilon)
  , sigma_(sigma)
  , wetting_(wetting)
  , cutoff_(cutoff)
  , smoothing_(smoothing)
  , g_param_geometry_(surface_normal_.size())
  , g_param_potential_(epsilon_.size1() * epsilon_.size2())
  , logger_(logger)
{
    unsigned int nwall = epsilon_.size1();
    unsigned int nspecies = epsilon_.size2();

    // check parameter size
    if (offset_.size() != nwall) {
        throw invalid_argument("geometry parameters have mismatching shapes");
    }

    if (epsilon_.size1() != nwall || sigma_.size1() != nwall || wetting_.size1() != nwall || cutoff_.size1() != nwall
        || sigma_.size2() != nspecies || wetting_.size2() != nspecies || cutoff_.size2() != nspecies
       ) {
        throw invalid_argument("potential parameters have mismatching shapes");
    }

    LOG("number of walls: " << nwall);
    LOG("wall position: d₀ = " << offset_);
    LOG("surface normal: n = " << surface_normal_);
    LOG("interaction strength: epsilon = " << epsilon_);
    LOG("interaction range: sigma = " << sigma_);
    LOG("wetting paramter: w = " << wetting_);
    LOG("cutoff length for wall potential: rc = " << cutoff_);
    LOG("smoothing parameter for wall potential: h = " << smoothing_);

    // impose normalisation of surface normals
    for (auto& n : surface_normal_) {
        n /= norm_2(n);
    }

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
              using namespace planar_wall_kernel;

              fixed_vector<float, 4> p;
              p[EPSILON] = epsilon_(i, j);
              p[SIGMA]   = sigma_(i, j);
              p[WETTING] = wetting_(i, j);
              p[CUTOFF]  = cutoff_(i, j);
              param_potential[j * nwall + i] = p;
        }
    }
    cuda::copy(param_potential, g_param_potential_);

    // copy CUDA symbols
    cuda::copy(nwall, planar_wall_wrapper::nwall);
    cuda::copy(smoothing_, planar_wall_wrapper::smoothing);
}

template <int dimension, typename float_type>
void planar_wall<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static string class_name("planar_wall_" + to_string(dimension));
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
                        class_<planar_wall, shared_ptr<planar_wall>>(class_name.c_str())
                            .def(constructor<
                                 scalar_container_type const&
                               , vector_container_type const&
                               , matrix_container_type const&
                               , matrix_container_type const&
                               , matrix_container_type const&
                               , matrix_container_type const&
                               , float_type
                               , shared_ptr<logger>
                             >())
                            .property("offset", &planar_wall::offset)
                            .property("surface_normal", &planar_wall::surface_normal)
                            .property("epsilon", &planar_wall::epsilon)
                            .property("sigma", &planar_wall::sigma)
                            .property("wetting", &planar_wall::wetting)
                            .property("cutoff", &planar_wall::cutoff)
                            .property("smoothing", &planar_wall::smoothing)
                    ]
                ]
            ]
        ]
    ];
}


HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_external_planar_wall(lua_State* L)
{
    planar_wall<3, float>::luaopen(L);
    planar_wall<2, float>::luaopen(L);
    forces::external<3, float, planar_wall<3, float>>::luaopen(L);
    forces::external<2, float, planar_wall<2, float>>::luaopen(L);
    return 0;
}

// explicit instantiation
template class planar_wall<3, float>;
template class planar_wall<2, float>;

} // namespace external
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
using namespace potentials::external;

template class external<3, float, planar_wall<3, float>>;
template class external<2, float, planar_wall<2, float>>;

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
