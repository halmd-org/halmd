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

#include <halmd/mdsim/forces/trunc/local_r4.hpp>
#include <halmd/mdsim/host/forces/pair_full.hpp>
#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/potentials/morse.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::assign;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {

/**
 * Initialise Morse potential parameters
 */
template <typename float_type>
morse<float_type>::morse(
    unsigned ntype1
  , unsigned ntype2
  , matrix_type const& cutoff
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , matrix_type const& r_min
  , std::shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(sigma)
  , r_min_sigma_(r_min)
  , r_cut_sigma_(cutoff)
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , en_cut_(ntype1, ntype2)
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned i = 0; i < en_cut_.size1(); ++i) {
        for (unsigned j = 0; j < en_cut_.size2(); ++j) {
            en_cut_(i, j) = get<1>((*this)(rr_cut_(i, j), i, j));
        }
    }

    LOG("depth of potential well: ε = " << epsilon_);
    LOG("width of potential well: σ = " << sigma_);
    LOG("position of potential well: r_min / σ = " << r_min_sigma_);
    LOG("cutoff radius of potential: r_c / σ = " << r_cut_sigma_);
    LOG("potential energy at cutoff: U = " << en_cut_);
}

template <typename float_type>
void morse<float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("potentials")
                [
                    class_<morse, std::shared_ptr<morse> >(module_name())
                        .def(constructor<
                            unsigned int
                          , unsigned int
                          , matrix_type const&
                          , matrix_type const&
                          , matrix_type const&
                          , matrix_type const&
                          , std::shared_ptr<logger_type>
                        >())
                        .property("r_cut", (matrix_type const& (morse::*)() const) &morse::r_cut)
                        .property("r_cut_sigma", &morse::r_cut_sigma)
                        .property("epsilon", &morse::epsilon)
                        .property("sigma", &morse::sigma)
                        .property("r_min_sigma", &morse::r_min_sigma)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_morse(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    morse<double>::luaopen(L);
    forces::pair_full<3, double, morse<double> >::luaopen(L);
    forces::pair_full<2, double, morse<double> >::luaopen(L);
    forces::pair_trunc<3, double, morse<double> >::luaopen(L);
    forces::pair_trunc<2, double, morse<double> >::luaopen(L);
    forces::pair_trunc<3, double, morse<double>, mdsim::forces::trunc::local_r4<double> >::luaopen(L);
    forces::pair_trunc<2, double, morse<double>, mdsim::forces::trunc::local_r4<double> >::luaopen(L);
#else
    morse<float>::luaopen(L);
    forces::pair_full<3, float, morse<float> >::luaopen(L);
    forces::pair_full<2, float, morse<float> >::luaopen(L);
    forces::pair_trunc<3, float, morse<float> >::luaopen(L);
    forces::pair_trunc<2, float, morse<float> >::luaopen(L);
    forces::pair_trunc<3, float, morse<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
    forces::pair_trunc<2, float, morse<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class morse<double>;
#else
template class morse<float>;
#endif

} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::morse<double> >;
template class pair_full<2, double, potentials::morse<double> >;
template class pair_trunc<3, double, potentials::morse<double> >;
template class pair_trunc<2, double, potentials::morse<double> >;
template class pair_trunc<3, double, potentials::morse<double>, mdsim::forces::trunc::local_r4<double> >;
template class pair_trunc<2, double, potentials::morse<double>, mdsim::forces::trunc::local_r4<double> >;
#else
template class pair_full<3, float, potentials::morse<float> >;
template class pair_full<2, float, potentials::morse<float> >;
template class pair_trunc<3, float, potentials::morse<float> >;
template class pair_trunc<2, float, potentials::morse<float> >;
template class pair_trunc<3, float, potentials::morse<float>, mdsim::forces::trunc::local_r4<float> >;
template class pair_trunc<2, float, potentials::morse<float>, mdsim::forces::trunc::local_r4<float> >;
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
