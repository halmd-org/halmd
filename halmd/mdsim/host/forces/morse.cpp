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
#include <halmd/mdsim/host/forces/morse.hpp>
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
void morse<float_type>::options(po::options_description& desc)
{
    desc.add_options()
//         ("cutoff", po::value<boost::array<float, 3> >()->default_value(default_cutoff()),
//          "truncate potential at cutoff radius")
//         ("epsilon", po::value<boost::array<float, 3> >()->default_value(default_epsilon()),
//          "depths of potential well for AA,AB,BB interactions")
//         ("sigma", po::value<boost::array<float, 3> >()->default_value(default_sigma()),
//          "widths of potential well for AA,AB,BB interactions")
        ("morse-minimum", po::value<boost::array<float, 3> >()->default_value(default_r_min()),
         "positions of potential minimum for AA,AB,BB interactions")
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
 * Initialise Morse potential parameters
 */
template <typename float_type>
morse<float_type>::morse(
    unsigned ntype
  , array<float, 3> const& cutoff
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
  , array<float, 3> const& r_min
)
  // allocate potential parameters
  : epsilon_(scalar_matrix<float_type>(ntype, ntype, 1))
  , sigma_(scalar_matrix<float_type>(ntype, ntype, 1))
  , r_min_(ntype, ntype)
  , r_min_sigma_(ntype, ntype)
  , en_cut_(scalar_matrix<float_type>(ntype, ntype, 0))
  , r_cut_(ntype, ntype)
  , r_cut_sigma_(ntype, ntype)
  , rr_cut_(ntype, ntype)
{
   // FIXME support any number of types
    for (unsigned i = 0; i < std::min(ntype, 2U); ++i) {
        for (unsigned j = i; j < std::min(ntype, 2U); ++j) {
            epsilon_(i, j) = epsilon[i + j];
            sigma_(i, j) = sigma[i + j];
            r_cut_sigma_(i, j) = cutoff[i + j];
            r_min_sigma_(i, j) = r_min[i + j];
        }
    }

    // precalculate derived parameters
    for (unsigned i = 0; i < ntype; ++i) {
        for (unsigned j = i; j < ntype; ++j) {
            r_min_(i, j) = r_min_sigma_(i, j) * sigma_(i, j);
            r_cut_(i, j) = r_cut_sigma_(i, j) * sigma_(i, j);
            rr_cut_(i, j) = std::pow(r_cut_(i, j), 2);
            // energy shift due to truncation at cutoff length
            en_cut_(i, j) = (*this)(rr_cut_(i, j), i, j).second;
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
    using namespace luabind;
    module(L, "halmd_wrapper")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("forces")
                [
                    class_<morse, shared_ptr<morse> >(module_name())
                        .def(constructor<
                            unsigned
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , array<float, 3> const&
                        >())
                        .scope
                        [
                            def("options", &morse::options)
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
        &morse<float_type>::luaopen
    ];

    lua_wrapper::register_(2) //< distance of derived to base class
    [
        &pair_trunc<3, float_type, morse<float_type> >::luaopen
    ]
    [
        &pair_trunc<2, float_type, morse<float_type> >::luaopen
    ];
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class morse<double>;
template class pair_trunc<3, double, morse<double> >;
template class pair_trunc<2, double, morse<double> >;
#else
template class morse<float>;
template class pair_trunc<3, float, morse<float> >;
template class pair_trunc<2, float, morse<float> >;
#endif

}}} // namespace mdsim::host::forces

} // namespace halmd
