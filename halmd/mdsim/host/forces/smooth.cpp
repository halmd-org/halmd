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

#include <boost/algorithm/string/predicate.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**
 * Initialise parameters
 */
template <int dimension, typename float_type>
smooth<dimension, float_type>::smooth(double r_smooth)
  // initialise parameters
  : r_smooth_(r_smooth)
  , rri_smooth_(std::pow(r_smooth_, -2))
{
    LOG("scale parameter for potential smoothing function: " << r_smooth_);
}

template <int dimension, typename float_type>
void smooth<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("smooth_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    namespace_("forces")
                    [
                        class_<smooth, shared_ptr<smooth> >(class_name.c_str())
                            .def(constructor<double>())
                    ]
                ]
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &smooth<3, double>::luaopen
    ]
    [
        &smooth<2, double>::luaopen
    ];
#else
    [
        &smooth<3, float>::luaopen
    ]
    [
        &smooth<2, float>::luaopen
    ];
#endif
}

} // namespace

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class smooth<3, double>;
template class smooth<2, double>;
#else
template class smooth<3, float>;
template class smooth<2, float>;
#endif

}}} // namespace mdsim::host::forces

} // namespace halmd
