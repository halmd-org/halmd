/*
 * Copyright © 2011-2020 Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <map>
#include <sstream>

#include <halmd/algorithm/host/pick_lattice_points.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/multi_index.hpp>

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
  , filter_type const& filter
)
  // initialise members
  : wavenumber_(wavenumber)
  , box_length_(box_length)
  , tolerance_(tolerance)
  , max_count_(max_count)
  , filter_(filter)
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
    if (filter_ != filter_type(1)) {
        LOG("apply filter on wavevectors: " << filter_);
        if (norm_inf(filter_) > 1 ) {
            throw std::invalid_argument("filter values must be 0 or 1");
        }
    }

    // construct wavevectors and store as key/value pairs (wavenumber iterator, wavevector).
    multimap<decltype(first), vector_type> wavevector_map;
    algorithm::host::pick_lattice_points_from_shell(
        first, last
      , inserter(wavevector_map, end(wavevector_map))
      , element_div(vector_type(2 * M_PI), box_length_)
      , tolerance_
      , max_count_
      , filter_
    );

    if (wavevector_map.empty()) {
        LOG_WARNING("no matching wavevectors");
        throw logic_error("Constraints on wavevectors are incompatible with geometry of simulation box.");
    }

    vector<decltype(first)> discarded;
    for (auto q_it = begin(wavenumber_); q_it != end(wavenumber_); ++q_it) {
        // find wavevector range with |q| = *q_it
        auto q_range = wavevector_map.equal_range(q_it);

        if (q_range.first != q_range.second) {
            // append wavevectors to wavevector list and store shell
            typename shell_array_type::value_type idx;
            idx.first = wavevector_.size();
            for (auto it = q_range.first; it != q_range.second; ++it) {
                assert(abs(*it->first - *q_it) < *q_it * tolerance_);
                assert(abs(norm_2(it->second) - *q_it) < 2 * *q_it * tolerance_);
                wavevector_.push_back(it->second);
            }
            idx.second = wavevector_.size();
            shell_.push_back(idx);
        }
        else {
            // remove wavenumbers with empty wavevector shells,
            // postpone deletion since we must not invalidate the iterators stored in wavevector_map
            LOG_WARNING("reciprocal lattice not compatible with |q| ≈ " << *q_it << ", value discarded");
            discarded.push_back(q_it);
        }
    }
    // remove discarded wavenumbers in reverse order
    reverse(begin(discarded), end(discarded));
    for (auto q_it : discarded) {
        wavenumber_.erase(q_it);
    }

    LOG_INFO("total number of wavevectors found: " << wavevector_.size());
}

template <int dimension>
wavevector<dimension>::wavevector(
    vector<double> const& wavenumber
  , vector_type const& box_length
  , filter_type const& filter
)
  // initialise members
  : wavenumber_(wavenumber)
  , box_length_(box_length)
  , tolerance_(0)
  , max_count_(0)
  , filter_(filter)
  , logger_(make_shared<logger>("wavevector"))
{
    typedef fixed_vector<int, dimension> index_type;

    // sort wavenumber array
    std::sort(begin(wavenumber_), end(wavenumber_));
    double q_max = wavenumber_.back();

    LOG("construct dense grid of wavevectors");
    LOG("maximum wavenumber: " << q_max);
    if (filter_ != filter_type(1)) {
        LOG("apply filter on wavevectors: " << filter_);
        if (norm_inf(filter_) > 1 ) {
            throw std::invalid_argument("filter values must be 0 or 1");
        }
    }

    // determine unit cell in reciprocal lattice, 2π / L[i]
    // and number of grid points per dimension up to q_max
    vector_type unit_cell = element_div(vector_type(2 * M_PI), box_length_);
    auto max_n = static_cast<index_type>(floor(element_div(vector_type(q_max), unit_cell)));

    // apply wavevector filter for each Cartesian component,
    // max_n[j] = 0 implies q[j] = 0 below.
    max_n = element_prod(max_n, static_cast<index_type>(filter_));

    LOG_INFO("grid points per dimension: " << (2 * max_n + index_type(1)));
#ifndef NDEBUG
    ostringstream s;
    copy(begin(wavenumber_), end(wavenumber_), ostream_iterator<double>(s, " "));
    LOG_DEBUG("wavenumber shells (upper bounds): " << s.str());
#endif

    // construct dense grid of wavevectors,
    // loop over multi-index: -max_n[i] ≤ n[i] ≤ max_n[i]
    for (index_type idx = -max_n; idx[dimension-1] <= max_n[dimension-1]; ) {
        // wavevector: q[i] = n[i] * 2 \pi / L[i]
        vector_type q = element_prod(unit_cell, static_cast<vector_type>(idx));
        // apply |q| < q_max
        if (inner_prod(q, q) <= q_max * q_max) {
            wavevector_.push_back(q);
        }

        // increment index tuple at end of loop,
        // obey -max_n[j] ≤ idx[j] ≤ max_n[j]
        ++idx[0];                            // always increment first 'digit' (element)
        for (unsigned int j = 0; j < dimension - 1; ++j) {
            if (idx[j] <= max_n[j]) {           // test upper bound
                break;                          // still within range, exit loop
            }
            idx[j] = -max_n[j];                 // reset this 'digit'
            ++idx[j+1];                         // increment next 'digit'
        }
    }
    LOG_DEBUG("total number of wavevectors: " << wavevector_.size());

    // partition into wavenumber shells
    double q_lower = 0;
    auto first = begin(wavevector_);
    auto last = end(wavevector_);
    auto start = first;
    for (auto q_it = begin(wavenumber_); q_it != end(wavenumber_); ) {
        // partition into wavevectors with norm smaller (<) and greater (≥) than q_bound,
        // returned iterator points at second partition
        double bound = (*q_it) * (*q_it);
        auto next = std::partition(start, last, [=](vector_type const& q) { return inner_prod(q, q) < bound; } );

        // if some wavevectors met the criterion
        if (next != start) {
            LOG_TRACE(next - start << " wavevectors in shell with " << q_lower << " ≤ |q| < " << *q_it);
            shell_.push_back(std::make_pair(start - first, next - first));
            start = next;
            q_lower = *q_it++;
        }
        else {
            // remove wavenumbers with empty wavevector shells,
            // postpone deletion since we must not invalidate our loop iterator
            LOG_WARNING("reciprocal lattice not compatible with " << q_lower << " ≤ |q| < " << *q_it << ", value discarded");
            q_it = wavenumber_.erase(q_it);
        }
    }
    assert(start == last);                        // no wavevectors left with |q| ≥ q_max
    assert(wavenumber_.size() == shell_.size());  // correspondence between wavenumbers and shells
}

template <typename wavevector_type>
static function<typename wavevector_type::wavevector_array_type const& ()>
wrap_value(shared_ptr<wavevector_type const> self)
{
    return [=]() -> typename wavevector_type::wavevector_array_type const& {
        return self->value();
    };
}

template <typename wavevector_type>
static function<typename wavevector_type::wavenumber_array_type const& ()>
wrap_wavenumber(shared_ptr<wavevector_type const> self)
{
    return [=]() -> typename wavevector_type::wavenumber_array_type const& {
        return self->wavenumber();
    };
}

template <typename T>
static bool equal(shared_ptr<T const> self, shared_ptr<T const> other)
{
    // compare pointers of managed objects. I could not get working
    // owner_equal() or owner_before() with shared_ptr's passed from Lua
    return self == other;
}


template <int dimension>
void wavevector<dimension>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("utility")
            [
                class_<wavevector>()
                    .property("wavenumber", &wrap_wavenumber<wavevector>)
                    .property("value", &wrap_value<wavevector>)
                    .def("shell", &wavevector::shell)
                    .def("__eq", &equal<wavevector>) // operator= in Lua
              , def("wavevector", &make_shared<wavevector
                  , vector<double> const&
                  , vector_type const&
                  , double
                  , unsigned int
                  , filter_type const&
                 >)
              , def("wavevector", &make_shared<wavevector
                  , vector<double> const&
                  , vector_type const&
                  , filter_type const&
                 >)

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
