/*
 * Copyright © 2011-2013  Felix Höfling
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
#include <exception>
#include <iterator>
#include <sstream>

#include <halmd/algorithm/host/pick_lattice_points.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace observables {
namespace utility {

template <int dimension>
wavevector<dimension>::wavevector(
    vector<double> const& wavenumber
  , vector_type const& box_length
  , double tolerance
  , unsigned int max_count
)
  // initialise members
  : wavenumber_(wavenumber)
  , box_length_(box_length)
  , tolerance_(tolerance)
  , max_count_(max_count)
  , logger_(make_shared<logger>("wavevector"))
{
    auto first = begin(wavenumber_);
    auto last = end(wavenumber_);
    LOG("use " << last - first << " wavenumbers"
        << " from " << *min_element(first, last) << " to " << *max_element(first, last)
    );
#ifndef NDEBUG
    ostringstream s;
    copy(first, last, ostream_iterator<double>(s, " "));
    LOG_DEBUG("wavenumber grid: " << s.str());
#endif

    LOG("tolerance on magnitude: " << tolerance_);
    LOG("maximum shell size: " << max_count_);

    // construct wavevectors and store as key/value pairs (wavenumber, wavevector)
    algorithm::host::pick_lattice_points_from_shell(
        first, last
      , back_inserter(wavevector_)
      , element_div(vector_type(2 * M_PI), box_length_)
      , tolerance_
      , max_count_
    );

    // sort wavevector map according to keys (wavenumber)
    stable_sort(
        begin(wavevector_), end(wavevector_)
      , bind(less<double>(), bind(&map_type::value_type::first, _1), bind(&map_type::value_type::first, _2))
    );

    // remove wavenumbers with no compatible wavevectors
    for (vector<double>::iterator q_it = begin(wavenumber_); q_it != end(wavenumber_); ++q_it) {
        // find wavevector q with |q| = *q_it
        auto found = find_if(
            begin(wavevector_), end(wavevector_)
          , bind(equal_to<double>(), bind(&map_type::value_type::first, _1), *q_it)
        );
        if (found == end(wavevector_)) {
            LOG_WARNING("reciprocal lattice not compatible with |q| ≈ " << *q_it << ", value discarded");
            wavenumber_.erase(q_it--);   // post-decrement iterator, increment at end of loop
        }
    }

    if (wavenumber_.empty()) {
        LOG_WARNING("empty wavenumber grid");
        throw std::logic_error("Constraints on wavevectors are incompatible with geometry of simulation box.");
    }

    LOG_DEBUG("total number of wavevectors: " << wavevector_.size());
}

template <typename wavevector_type, typename wavenumber_array_type>
static std::function<wavenumber_array_type const& ()>
wrap_wavenumber(std::shared_ptr<wavevector_type const> self)
{
    return [=]() -> wavenumber_array_type const& {
        return self->wavenumber();
    };
}

template <int dimension>
static unsigned int wrap_dimension(wavevector<dimension> const&)
{
    return dimension;
}

template <int dimension>
void wavevector<dimension>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static string class_name("wavevector_" + boost::lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("utility")
            [
                class_<wavevector, std::shared_ptr<wavevector> >(class_name.c_str())
                    .def(constructor<
                         vector<double> const&
                       , vector_type const&
                       , double, unsigned int
                    >())
                    .property("wavenumber", &wrap_wavenumber<wavevector, wavenumber_array_type>)
                    .property("value", &wavevector::value)
                    .property("tolerance", &wavevector::tolerance)
                    .property("max_count", &wavevector::max_count)
                    .property("dimension", &wrap_dimension<dimension>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_utility_wavevector(lua_State* L)
{
    wavevector<3>::luaopen(L);
    wavevector<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class wavevector<3>;
template class wavevector<2>;

} // namespace observables
} // namespace utility
} // namespace halmd
