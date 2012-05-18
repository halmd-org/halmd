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
#include <exception>
#include <iterator>
#include <sstream>

#include <halmd/algorithm/host/pick_lattice_points.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
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
{
    ostringstream s;
    copy(wavenumber_.begin(), wavenumber_.end(), ostream_iterator<double>(s, " "));
    LOG("wavenumber grid: " << s.str());

    init_();
}

template <int dimension>
wavevector<dimension>::wavevector(
    double max_wavenumber
  , unsigned int decimation
  , vector_type const& box_length
  , double tolerance
  , unsigned int max_count
)
  // initialise members
  : box_length_(box_length)
  , tolerance_(tolerance)
  , max_count_(max_count)
{
    LOG("maximum wavenumber: " << max_wavenumber);

    // set up semi-linearly spaced wavenumber grid
    // determine q_min for the initial spacing
    double h = 2 * M_PI / norm_inf(box_length_); //< norm_inf returns the maximum coordinate
    unsigned int i = 0;
    for (double q = h; q < max_wavenumber; ) {
        wavenumber_.push_back(q);
        q += h;
        // double grid spacing after every 'decimation' number of points
        if (decimation > 0 && ++i % decimation == 0) {
            h *= 2;
        }
    }

    init_();
}

template <int dimension>
void wavevector<dimension>::init_()
{
    LOG("tolerance on wavevector magnitude: " << tolerance_);
    LOG("maximum number of wavevectors per wavenumber: " << max_count_);

    // construct wavevectors and store as key/value pairs (wavenumber, wavevector)
    algorithm::host::pick_lattice_points_from_shell(
        wavenumber_.begin(), wavenumber_.end()
      , back_inserter(wavevector_)
      , element_div(vector_type(2 * M_PI), box_length_)
      , tolerance_
      , max_count_
    );

    // sort wavevector map according to keys (wavenumber)
    stable_sort(
        wavevector_.begin(), wavevector_.end()
      , bind(less<double>(), bind(&map_type::value_type::first, _1), bind(&map_type::value_type::first, _2))
    );

    // remove wavenumbers with no compatible wavevectors
    for (vector<double>::iterator q_it = wavenumber_.begin(); q_it != wavenumber_.end(); ++q_it) {
        // find wavevector q with |q| = *q_it
        typename map_type::const_iterator found = find_if(
            wavevector_.begin(), wavevector_.end()
          , bind(equal_to<double>(), bind(&map_type::value_type::first, _1), *q_it)
        );
        if (found == wavevector_.end()) {
            LOG_WARNING("No wavevector compatible with |q| ≈ " << *q_it << ". Value discarded");
            wavenumber_.erase(q_it--);   // post-decrement iterator, increment at end of loop
        }
    }

    if (wavenumber_.empty()) {
        LOG_WARNING("Wavenumber grid is empty.");
        throw std::logic_error("Constraints on wavevectors are incompatible with geometry of simulation box.");
    }

    LOG_DEBUG("total number of wavevectors: " << wavevector_.size());
}

template <typename wavevector_type, typename wavenumber_array_type>
static boost::function<wavenumber_array_type const& ()>
wrap_wavenumber(boost::shared_ptr<wavevector_type const> wavevector)
{
    return bind(&wavevector_type::wavenumber, wavevector);
}

template <int dimension>
static unsigned int wrap_dimension(wavevector<dimension> const&)
{
    return dimension;
}

template <int dimension>
void wavevector<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("wavevector_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("utility")
            [
                class_<wavevector, boost::shared_ptr<wavevector> >(class_name.c_str())
                    .def(constructor<
                         vector<double> const&
                       , vector_type const&
                       , double, unsigned int
                    >())
                    .def(constructor<
                         double, unsigned int
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
