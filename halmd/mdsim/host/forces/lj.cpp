/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/forces/lj.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace boost::assign;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**
 * Assemble module options
 */
template <typename float_type>
void lj_potential<float_type>::options(options_definition& options)
{
    options.add_options()
        ("cutoff", po::value<boost::array<float, 3> >()->default_value(default_cutoff()),
         "truncate potential at cutoff radius")
        ("epsilon", po::value<boost::array<float, 3> >()->default_value(default_epsilon()),
         "potential well depths AA,AB,BB")
        ("sigma", po::value<boost::array<float, 3> >()->default_value(default_sigma()),
         "collision diameters AA,AB,BB")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    register_any_converter<boost::array<float, 3> >();
}

/**
 * Initialise Lennard-Jones potential parameters
 */
template <typename float_type>
lj_potential<float_type>::lj_potential(
    unsigned ntype
  , array<float, 3> const& cutoff
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
)
  // allocate potential parameters
  : epsilon_(scalar_matrix<float_type>(ntype, ntype, 1))
  , sigma_(scalar_matrix<float_type>(ntype, ntype, 1))
  , r_cut_sigma_(ntype, ntype)
  , r_cut_(ntype, ntype)
  , rr_cut_(ntype, ntype)
  , sigma2_(ntype, ntype)
  , en_cut_(scalar_matrix<float_type>(ntype, ntype, 0))
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
            sigma2_(i, j) = std::pow(sigma_(i, j), 2);
            // energy shift due to truncation at cutoff length
            en_cut_(i, j) = (*this)(rr_cut_(i, j), i, j).second;
        }
    }

    LOG("potential well depths: ε = " << epsilon_);
    LOG("potential core width: σ = " << sigma_);
    LOG("potential cutoff length: r_c = " << r_cut_sigma_);
    LOG("potential cutoff energy: U = " << en_cut_);
}

template <typename float_type>
void lj_potential<float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    string class_name("lj_potential");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    namespace_("forces")
                    [
                        class_<lj_potential, shared_ptr<lj_potential> >(class_name.c_str())
                            .def(constructor<
                                unsigned
                              , array<float, 3> const&
                              , array<float, 3> const&
                              , array<float, 3> const&
                            >())
                            .scope
                            [
                                def("options", &lj_potential::options)
                            ]
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
#ifndef USE_HOST_SINGLE_PRECISION
    typedef double float_type;
#else
    typedef float float_type;
#endif

    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &lj_potential<float_type>::luaopen
    ];

    lua_wrapper::register_(2) //< distance of derived to base class
    [
        &pair_trunc<3, float_type, lj_potential<float_type> >::luaopen
    ]
    [
        &pair_trunc<2, float_type, lj_potential<float_type> >::luaopen
    ];
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class lj_potential<double>;
template class pair_trunc<3, double, lj_potential<double> >;
template class pair_trunc<2, double, lj_potential<double> >;
#else
template class lj_potential<float>;
template class pair_trunc<3, float, lj_potential<float> >;
template class pair_trunc<2, float, lj_potential<float> >;
#endif

}}} // namespace mdsim::host::forces

} // namespace halmd
