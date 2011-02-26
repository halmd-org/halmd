/*
 * Copyright © 2011  Felix Höfling
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

#include <boost/bind.hpp>
#include <cmath>
#include <iterator>
#include <sstream>

#include <halmd/algorithm/host/pick_lattice_points.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/observables/utility/wavevectors.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace utility
{

template <int dimension>
wavevectors<dimension>::wavevectors(
    vector<double> const& wavenumbers
  , vector_type const& box_length
  , double tolerance
  , unsigned int max_count
)
  // initialise members
  : wavenumbers_(wavenumbers)
  , tolerance_(tolerance)
  , max_count_(max_count)
{
    ostringstream s;
    copy(wavenumbers_.begin(), wavenumbers_.end(), ostream_iterator<double>(s, " "));
    LOG("wavenumber grid: " << s.str());
    LOG("tolerance on wavevector magnitude: " << tolerance_);
    LOG("maximum number of wavevectors per wavenumber: " << max_count_);

    // construct wavevectors and store as key/value pairs (wavenumber, wavevector)
    algorithm::host::pick_lattice_points_from_shell(
        wavenumbers.begin(), wavenumbers.end()
      , back_inserter(wavevectors_)
      , element_div(vector_type(2 * M_PI), box_length)
      , tolerance
      , max_count
    );

    // sort wavevector map according to keys (wavenumber)
    stable_sort(
        wavevectors_.begin(), wavevectors_.end()
      , bind(less<double>(), bind(&map_type::value_type::first, _1), bind(&map_type::value_type::first, _2))
    );

    // remove wavenumbers with no compatible wavevectors
    for (vector<double>::iterator q_it = wavenumbers_.begin(); q_it != wavenumbers_.end(); ++q_it) {
        // find wavevector q with |q| = *q_it
        typename map_type::const_iterator found = find_if(
            wavevectors_.begin(), wavevectors_.end()
          , bind(equal_to<double>(), bind(&map_type::value_type::first, _1), *q_it)
        );
        if (found == wavevectors_.end()) {
            LOG_WARNING("No wavevector compatible with |q| ≈ " << *q_it << ". Value discarded");
            wavenumbers_.erase(q_it--);   // post-decrement iterator, increment at end of loop
        }
    }

    LOG_DEBUG("total number of wavevectors: " << wavevectors_.size());
}

template <int dimension>
void wavevectors<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("wavevectors_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("observables")
            [
                namespace_("utility")
                [
                    class_<wavevectors, shared_ptr<wavevectors> >(class_name.c_str())
                        .def(constructor<
                             vector<double> const&
                           , vector_type const&
                           , double, unsigned int
                        >())
                        .property("wavenumbers", &wavevectors::wavenumbers)
                        .property("values", &wavevectors::values)
                        .property("tolerance", &wavevectors::tolerance)
                        .property("maximum_count", &wavevectors::maximum_count)
                ]
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &wavevectors<3>::luaopen
    ]
    [
        &wavevectors<2>::luaopen
    ];
}

} // namespace

// explicit instantiation
template class wavevectors<3>;
template class wavevectors<2>;

}} // namespace observables::utility

} // namespace halmd
