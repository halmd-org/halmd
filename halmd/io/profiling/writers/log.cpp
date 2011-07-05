/*
 * Copyright © 2010  Peter Colberg and Felix Höfling
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

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <iomanip>

#include <halmd/io/logger.hpp>
#include <halmd/io/profiling/writers/log.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace std { // needed for Boost.Log << 1.0

using namespace halmd;

/**
 * output accumulator results to stream,
 * use a suitable unit of time
 */
template <typename T>
static ostream& operator<<(ostream& os, accumulator<T> const& acc)
{
    T value = mean(acc);
    T error = count(acc) > 1 ? error_of_mean(acc) : 0;

    T const conversion[] = { 24 * 3600, 3600, 60, 1, 1e-3, 1e-6, 1e-9 };
    char const* const unit[] = { "d", "h", "min", "s", "ms", "µs", "ns" };
    unsigned const N = 7;

    unsigned i;
    for (i = 0; i < N - 1; ++i) {
        if (value > conversion[i]) {
            break;
        }
    }
    value /= conversion[i];
    error /= conversion[i];

    // let number of digits depend on the error value
    if (count(acc) > 1) {
        unsigned prec = static_cast<unsigned>(max(0., ceil(-log10(error)) + 1));
        os << fixed << setprecision(prec) << value << " " << unit[i];
        // at least 2 digits for the error
        prec = static_cast<unsigned>(max(2., ceil(log10(error))));
        os << resetiosflags(ios_base::floatfield) << setprecision(prec)
           << " (" << error << " " << unit[i] << ", " << count(acc) << " calls)";
    }
    else {
        os << setprecision(3) << value << " " << unit[i];
    }

    return os;
}

} // namespace std

namespace halmd {
namespace io {
namespace profiling {
namespace writers {

/**
 * compare total accumulated runtimes of acc_desc_pairs
 */
template <typename T>
bool less_total_runtime(T x, T y)
{
    return mean(*x.first) * count(*x.first) < mean(*y.first) * count(*y.first);
}

/**
 * write log entries for all runtime accumulators,
 * sorted by their total accumulated runtime
 */
void log::write()
{
    stable_sort(
        accumulators_.begin()
      , accumulators_.end()
      , bind(&less_total_runtime<acc_desc_pair_type>, _1, _2)
    );

    BOOST_FOREACH(acc_desc_pair_type const& x, accumulators_) {
        LOG(x.second << ": " << *x.first);
    }
}

void log::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("profiling")
            [
                namespace_("writers")
                [
                    class_<log, shared_ptr<_Base>, _Base>("log")
                        .def(constructor<>())
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_profiling_writers_log(lua_State* L)
{
    log::luaopen(L);
    return 0;
}

} // namespace io
} // namespace profiling
} // namespace writers
} // namespace halmd
