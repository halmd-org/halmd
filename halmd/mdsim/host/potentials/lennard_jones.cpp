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

#include <halmd/mdsim/host/potentials/lennard_jones.hpp>
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
 * Initialise Lennard-Jones potential parameters
 */
template <typename float_type>
lennard_jones<float_type>::lennard_jones(
    unsigned int ntype1
  , unsigned int ntype2
  , matrix_type const& cutoff
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(sigma)
  , r_cut_sigma_(cutoff)
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , sigma2_(element_prod(sigma_, sigma_))
  , en_cut_(ntype1, ntype2)
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned i = 0; i < en_cut_.size1(); ++i) {
        for (unsigned j = 0; j < en_cut_.size2(); ++j) {
            en_cut_(i, j) = (*this)(rr_cut_(i, j), i, j).get<1>();
        }
    }

    LOG("potential well depths: ε = " << epsilon_);
    LOG("potential core width: σ = " << sigma_);
    LOG("potential cutoff length: r_c = " << r_cut_sigma_);
    LOG("potential cutoff energy: U = " << en_cut_);
}

template <typename float_type>
void lennard_jones<float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("potentials")
                [
                    class_<lennard_jones, shared_ptr<lennard_jones> >(module_name())
                        .def(constructor<
                            unsigned int
                          , unsigned int
                          , matrix_type const&
                          , matrix_type const&
                          , matrix_type const&
                          , shared_ptr<logger_type>
                        >())
                        .property("r_cut", (matrix_type const& (lennard_jones::*)() const) &lennard_jones::r_cut)
                        .property("r_cut_sigma", &lennard_jones::r_cut_sigma)
                        .property("epsilon", &lennard_jones::epsilon)
                        .property("sigma", &lennard_jones::sigma)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_lennard_jones(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    lennard_jones<double>::luaopen(L);
#else
    lennard_jones<float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class lennard_jones<double>;
#else
template class lennard_jones<float>;
#endif

} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd