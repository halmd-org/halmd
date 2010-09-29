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
#include <halmd/io/profile/writers/log.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace std // needed for Boost.Log << 1.0
{

using namespace halmd;

/**
 * output accumulator results to stream
 */
template <typename T>
static ostream& operator<<(ostream& os, accumulator<T> const& acc)
{
    os << fixed << setprecision(4) << mean(acc) * 1E3 << " ms";
    if (count(acc) > 1) {
        os << " (" << error_of_mean(acc) * 1E3 << " ms, " << count(acc) << " calls)";
    }
    return os;
}

} // namespace std

namespace halmd
{
namespace io { namespace profile { namespace writers
{

/**
 * compare total accumulated runtimes of acc_desc_pairs
 */
template <typename T>
bool greater_total_runtime(T x, T y)
{
    return mean(*x.first) * count(*x.first) > mean(*y.first) * count(*y.first);
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
      , bind(&greater_total_runtime<acc_desc_pair_type>, _1, _2)
    );

    BOOST_FOREACH(acc_desc_pair_type const& x, accumulators_) {
        LOG(x.second << ": " << *x.first);
    }
}

static __attribute__((constructor)) void register_lua()
{
    typedef log::_Base _Base;

    using namespace luabind;
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        namespace_("halmd_wrapper")
        [
            namespace_("io")
            [
                namespace_("profile")
                [
                    namespace_("writers")
                    [
                        class_<log, shared_ptr<log>, _Base>("log")
                            .def(constructor<>())
                    ]
                ]
            ]
        ]
    ];
}

}}} // namespace io::profile::writers

} // namespace halmd
