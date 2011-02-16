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

#include <cmath>
#include <iterator>
#include <sstream>

#include <halmd/algorithm/host/pick_lattice_points.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/observables/utility/wavevectors.hpp>

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

    algorithm::host::pick_lattice_points_from_shell(
        wavenumbers.begin(), wavenumbers.end()
      , inserter(wavevectors_, wavevectors_.begin())
      , element_div(vector_type(2 * M_PI), box_length)
      , tolerance
      , max_count
    );

    // remove wavenumbers with no compatible wavevectors
    for (vector<double>::iterator q_it = wavenumbers_.begin(); q_it != wavenumbers_.end(); ++q_it) {
        if (!wavevectors_.count(*q_it)) {
            LOG_WARNING("No wavevector compatible with |q| ≈ " << *q_it << ". Value discarded");
            wavenumbers_.erase(q_it--);   // post-decrement iterator, increment at end of loop
        }
    }

    LOG_DEBUG("total number of wavevectors: " << wavevectors_.size());
}

// explicit instantiation
template class wavevectors<3>;
template class wavevectors<2>;

}} // namespace observables::utility

} // namespace halmd
