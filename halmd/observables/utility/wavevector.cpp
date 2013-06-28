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

#include <cmath>
#include <exception>
#include <iterator>
#include <map>
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

    // construct wavevectors and store as key/value pairs (wavenumber, wavevector).
    //
    // The use of floating-point numbers as keys could result in rounding
    // issues (which is unlikely since they are not subject to arithmetic
    // operations), hence we consider keys equivalent if they differ by less
    // than the floating point precision.
    // x ~ y <=> !less_tol(x, y) && !less_tol(y, x) <=> |x – y| ≤ ε max(|x|, |y|)
    // assuming x, y ≥ 0:
    // !(x ~ y) <=> y – x > ε y or x – y > ε x <=> x < y(1 – ε) or y < x(1 – ε)
    auto less_tol = [=](double x, double y) { // comparison functor
        return x < y * (1 - numeric_limits<double>::epsilon());
    };
    // Note that the comparison predicate passed to std::multimap must define a
    // strict weak ordering, i.e. the equivalence relation must be transitive:
    // a ~ b, b ~ c ⇒ a ~ c. This is only the case for less_tol() if the
    // supplied keys are separated by more than ε.
    // (Counter example for ε = 0.1: 1 ~ 1.1, 1.1 ~ 1.2, but not 1 ~ 1.2).
    multimap<double, vector_type, decltype(less_tol)> wavevector_map(less_tol);
    algorithm::host::pick_lattice_points_from_shell(
        first, last
      , inserter(wavevector_map, end(wavevector_map))
      , element_div(vector_type(2 * M_PI), box_length_)
      , tolerance_
      , max_count_
    );

    if (wavevector_map.empty()) {
        LOG_WARNING("no matching wavevectors");
        throw logic_error("Constraints on wavevectors are incompatible with geometry of simulation box.");
    }

    for (auto q_it = begin(wavenumber_); q_it != end(wavenumber_); ++q_it) {
        // find wavevector range with |q| = *q_it
        auto q_range = wavevector_map.equal_range(*q_it);

        if (q_range.first != q_range.second) {
            // append wavevectors to wavevector list and store shell
            typename shell_array_type::value_type idx;
            idx.first = wavevector_.size();
            for (auto it = q_range.first; it != q_range.second; ++it) {
                assert(abs(it->first - *q_it) < *q_it * tolerance_);
                assert(abs(norm_2(it->second) - *q_it) < 2 * *q_it * tolerance_);
                wavevector_.push_back(it->second);
            }
            idx.second = wavevector_.size();
            shell_.push_back(idx);
        }
        else {
            // remove wavenumbers with empty wavevector shells
            LOG_WARNING("reciprocal lattice not compatible with |q| ≈ " << *q_it << ", value discarded");
            wavenumber_.erase(q_it--);   // post-decrement iterator, increment at end of loop
        }
    }

    LOG_DEBUG("total number of wavevectors found: " << wavevector_.size());
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
                    .def(constructor<
                         vector<double> const&
                       , vector_type const&
                       , double, unsigned int
                    >())
                    .property("wavenumber", &wrap_wavenumber<wavevector>)
                    .property("value", &wrap_value<wavevector>)
                    .def("__eq", &equal<wavevector>) // operator= in Lua

              , def("wavevector", &make_shared<wavevector
                  , vector<double> const&
                  , vector_type const&
                  , double
                  , unsigned int
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
