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

#include <halmd/observables/host/mean_quartic_displacement.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace host
{

/**
 * Compute mean-quartic displacement of two position sample vectors.
 *
 * @param first particles positions of one species at time t1
 * @param second particles positions of one species at time t2
 * @returns accumulated mean-quartic displacement
 */
template <int dimension, typename float_type>
typename mean_quartic_displacement<dimension, float_type>::result_type
mean_quartic_displacement<dimension, float_type>::compute(
    sample_vector const& first
  , sample_vector const& second
)
{
    result_type acc;
    typename sample_vector::const_iterator r2, r1, end = first.end();
    for (r2 = first.begin(), r1 = second.begin(); r2 != end; ++r2, ++r1) {
        // displacement of particle
        vector_type dr = *r1 - *r2;
        // square displacement
        float_type rr = inner_prod(dr, dr);
        // accumulate quartic displacement
        acc(rr * rr);
    }
    return acc;
}

template <int dimension, typename float_type>
void mean_quartic_displacement<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("mean_quartic_displacement_" + lexical_cast<string>(dimension) + "_");
    module(L, "halmd_wrapper")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<mean_quartic_displacement>(class_name.c_str())
                    .def(constructor<>())
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &mean_quartic_displacement<3, double>::luaopen
    ]
    [
        &mean_quartic_displacement<2, double>::luaopen
    ];
#else
    [
        &mean_quartic_displacement<3, float>::luaopen
    ]
    [
        &mean_quartic_displacement<2, float>::luaopen
    ];
#endif
}

} // namespace

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class mean_quartic_displacement<3, double>;
template class mean_quartic_displacement<2, double>;
#else
template class mean_quartic_displacement<3, float>;
template class mean_quartic_displacement<2, float>;
#endif

}} // namespace observables::host

} // namespace halmd
