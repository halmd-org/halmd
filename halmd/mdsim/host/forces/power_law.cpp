/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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
#include <cmath>
#include <string>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/mdsim/host/forces/power_law.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace forces {

/**
 * Initialise potential parameters
 */
template <typename float_type>
power_law<float_type>::power_law(
    unsigned int ntype
  , int index
  , array<float, 3> const& cutoff
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
  , shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : index_(index)
  , epsilon_(scalar_matrix<float_type>(ntype, ntype, 1))
  , sigma_(scalar_matrix<float_type>(ntype, ntype, 1))
  , r_cut_(ntype, ntype)
  , r_cut_sigma_(ntype, ntype)
  , rr_cut_(ntype, ntype)
  , en_cut_(scalar_matrix<float_type>(ntype, ntype, 0))
  , logger_(logger)
{
    // FIXME support any number of types
    for (unsigned i = 0; i < std::min(ntype, 2U); ++i) {
        for (unsigned j = i; j < std::min(ntype, 2U); ++j) {
            epsilon_(i, j) = epsilon[i + j];
            sigma_(i, j) = sigma[i + j];
            r_cut_sigma_(i, j) = cutoff[i + j];
        }
    }

    // precalculate derived parameters
    for (unsigned i = 0; i < ntype; ++i) {
        for (unsigned j = i; j < ntype; ++j) {
            r_cut_(i, j) = r_cut_sigma_(i, j) * sigma_(i, j);
            rr_cut_(i, j) = std::pow(r_cut_(i, j), 2);
            // energy shift due to truncation at cutoff length
            en_cut_(i, j) = (*this)(rr_cut_(i, j), i, j).get<1>();
        }
    }

    LOG("potential: power law index: n = " << index_);
    LOG("potential: interaction strength ε = " << epsilon_);
    LOG("potential: interaction range σ = " << sigma_);
    LOG("potential: cutoff length r_c = " << r_cut_sigma_);
    LOG("potential: cutoff energy U = " << en_cut_);
}

template <typename float_type>
void power_law<float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("forces")
                [
                    class_<power_law, shared_ptr<power_law> >(module_name())
                        .def(constructor<
                            unsigned int
                          , int
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , shared_ptr<logger_type>
                        >())
                        .property("index", &power_law::index)
                        .property("r_cut", (matrix_type const& (power_law::*)() const) &power_law::r_cut)
                        .property("r_cut_sigma", &power_law::r_cut_sigma)
                        .property("epsilon", &power_law::epsilon)
                        .property("sigma", &power_law::sigma)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_forces_power_law(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    typedef double float_type;
#else
    typedef float float_type;
#endif
    power_law<float_type>::luaopen(L);
    pair_trunc<3, float_type, power_law<float_type> >::luaopen(L);
    pair_trunc<2, float_type, power_law<float_type> >::luaopen(L);
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class power_law<double>;
template class pair_trunc<3, double, power_law<double> >;
template class pair_trunc<2, double, power_law<double> >;
#else
template class power_law<float>;
template class pair_trunc<3, float, power_law<float> >;
template class pair_trunc<2, float, power_law<float> >;
#endif

} // namespace mdsim
} // namespace host
} // namespace forces
} // namespace halmd
