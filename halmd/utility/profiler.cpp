/*
 * Copyright © 2010-2024 Felix Höfling
 * Copyright © 2010-2011 Peter Colberg
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
#include <boost/foreach.hpp>
#include <iomanip>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {

// put this operator<< into the namespace of accumulator, in order to make ADL
// find the implementation
namespace numeric {
namespace detail {
/**
 * output accumulator results to stream,
 * use a suitable unit of time
 */
template <typename value_type>
static std::ostream& operator<<(std::ostream& os, accumulator<value_type> const& acc)
{
    value_type value = mean(acc);
    value_type error = count(acc) > 1 ? error_of_mean(acc) : 0;

    value_type const conversion[] = { 24 * 3600, 3600, 60, 1, 1e-3, 1e-6, 1e-9 };
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
        unsigned prec = static_cast<unsigned>(std::max(0., ceil(-log10(error)) + 1));
        os << std::fixed << std::setprecision(prec) << value << " " << unit[i];
        // at least 2 digits for the error
        prec = static_cast<unsigned>(std::max(2., ceil(log10(error))));
        os << resetiosflags(std::ios_base::floatfield) << std::setprecision(prec)
           << " (" << error << " " << unit[i] << ", " << count(acc) << " calls)";
    }
    else {
        os << std::setprecision(3) << value << " " << unit[i];
    }

    return os;
}
} // namespace detail
} // namespace numeric

namespace utility {

profiler::profiler()
{
    LOG("profiler timer resolution: " << 1.E9 * timer_type::elapsed_min() << " ns");
}

connection profiler::on_profile(std::shared_ptr<accumulator_type> acc, std::string const& desc)
{
    return accumulators_.connect(make_pair(acc, desc));
}

connection profiler::on_prepend_profile(slot_function_type const& slot)
{
    return on_prepend_profile_.connect(slot);
}

connection profiler::on_append_profile(slot_function_type const& slot)
{
    return on_append_profile_.connect(slot);
}

void profiler::profile()
{
    on_prepend_profile_();
    log();
    for (slots_const_iterator acc = accumulators_.begin(); acc != accumulators_.end(); ++acc) {
        acc->first->reset();
    }
    on_append_profile_();
}


/**
 * return total runtime from runtime accumulator:
 * average time per call × number of calls
 */
template <typename accumulator_type>
inline typename accumulator_type::value_type total_runtime(accumulator_type const& acc)
{
    return mean(acc) * count(acc);
}

/**
 * compare total accumulated runtimes of accumulators
 */
struct less_total_runtime
{
    template <typename accumulator_pair_type>
    bool operator()(accumulator_pair_type const& acc1, accumulator_pair_type const& acc2) const
    {
        return total_runtime(*acc1.first) < total_runtime(*acc2.first);
    }
};

/**
 * write log entries for all runtime accumulators,
 * sorted by their total accumulated runtime
 */
void profiler::log() const
{
    if (accumulators_.empty()) return;

    std::vector<accumulator_pair_type> accumulators(accumulators_.begin(), accumulators_.end());

    stable_sort(accumulators.begin(), accumulators.end(), less_total_runtime());

    double maximum_runtime = total_runtime(*accumulators.back().first);
    BOOST_FOREACH(accumulator_pair_type const& acc, accumulators) {
        double fraction = total_runtime(*acc.first) / maximum_runtime;
        logging::severity_level log_level = logging::info;  // default logging level: INFO
        if (fraction >= 0.20) {                             // increase level of profilers that cover more than 20% of the total runtime
            log_level = logging::message;
        }
        if (count(*acc.first) == 0) {                       // show empty profilers only at DEBUG level
            log_level = logging::debug;
        }

        HALMD_LOG(
            log_level
          , "[" << std::setw(5) << std::fixed << std::setprecision(1) << fraction * 100 << "%] "
                << acc.second << ": " << *acc.first
        );
    }
}

void profiler::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("utility")
        [
            class_<profiler, std::shared_ptr<profiler> >("profiler")
                .def(constructor<>())
                .def("on_profile", &profiler::on_profile)
                .def("on_prepend_profile", &profiler::on_prepend_profile)
                .def("on_append_profile", &profiler::on_append_profile)
                .def("profile", &profiler::profile)
                .scope
                [
                    class_<accumulator_type>("accumulator")
                ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_utility_profiler(lua_State* L)
{
    profiler::luaopen(L);
    return 0;
}

} // namespace utility
} // namespace halmd
