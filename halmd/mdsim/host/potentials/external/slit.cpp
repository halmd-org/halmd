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
#include <stdexcept>
#include <string>

#include <halmd/mdsim/host/forces/external.hpp>
#include <halmd/mdsim/host/potentials/external/slit.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
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
  , matrix_container_type const& cutoff
  , float_type smoothing
  , shared_ptr<logger> logger
)
  // allocate potential parameters
  : offset_(offset)
  , surface_normal_(surface_normal)
  , epsilon_(epsilon)
  , sigma_(sigma)
  , wetting_(wetting)
  , cutoff_(cutoff)
  , smoothing_(smoothing)
  , logger_(logger)
{

    // check parameter size
    if (offset_.size() != sigma_.size1()) {
        throw invalid_argument("geometry parameters have mismatching shapes");
    }

    if (epsilon_.size1() != sigma_.size1() || epsilon_.size1() != wetting_.size1() 
     || epsilon_.size2() != sigma_.size2() || epsilon_.size2() != wetting_.size2()
       ) {
        throw invalid_argument("potential parameters have mismatching shapes");
    }

    LOG("wall position: d₀ = " << offset_);
    LOG("surface normal: n = " << surface_normal_);
    LOG("interaction strength: epsilon = " << epsilon_);
    LOG("interaction range: sigma = " << sigma_);
    LOG("wetting paramter: w = " << wetting_);
    LOG("cutoff length: rc = " << cutoff_);
    LOG("smoothing parameter: h = " << smoothing_);

    // impose normalisation of surface normals
    for (auto& n : surface_normal_) {
        n /= norm_2(n);
    }
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
            namespace_("host")
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
                               , matrix_container_type const& 
                               , float_type
                               , shared_ptr<logger>
                             >())
                            .property("offset", &slit::offset)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_external_slit(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    slit<3, double>::luaopen(L);
    slit<2, double>::luaopen(L);
    forces::external<3, double, slit<3, double>>::luaopen(L);
    forces::external<2, double, slit<2, double>>::luaopen(L);
#else
    slit<3, float>::luaopen(L);
    slit<2, float>::luaopen(L);
    forces::external<3, float, slit<3, float>>::luaopen(L);
    forces::external<2, float, slit<2, float>>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class slit<3, double>;
template class slit<2, double>;
#else
template class slit<3, float>;
template class slit<2, float>;
#endif

} // namespace external
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
using namespace potentials::external;

#ifndef USE_HOST_SINGLE_PRECISION
template class external<3, double, slit<3, double>>;
template class external<2, double, slit<2, double>>;
#else
template class external<3, float, slit<3, float>>;
template class external<2, float, slit<2, float>>;
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
