/*
 * Copyright © 2013  Nicolas Höft
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

#ifndef HALMD_OBSERVABLES_UTILITY_ACCUMULATOR_HPP
#define HALMD_OBSERVABLES_UTILITY_ACCUMULATOR_HPP

#include <halmd/io/logger.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <functional>

namespace halmd {
namespace observables {
namespace utility {

/**
 * Calculate a sum based on an acquire function. This class basically wraps
 * numeric/accumulator and binds its methods to Lua.
 */
template <typename sample_type>
class accumulator
{
public:
    typedef sample_type value_type;
    typedef numeric::detail::accumulator<value_type> accumulator_type;
    typedef std::function<value_type ()> sample_function_type;

    /**
     * Initialise accumulator
     */
    accumulator(
        sample_function_type const& sample
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    )
      : sample_(sample)
      , logger_(logger)
    { }

    /** Aquire a new sample and add this to the total sum */
    void sample()
    {
        LOG_TRACE("acquire sample");
        acc_(sample_());
    }

    /** return the total sum of samples aquired */
    value_type sum() const
    {
        return numeric::detail::sum(acc_);
    }

    /** mean value */
    value_type mean() const
    {
        return numeric::detail::mean(acc_);
    }

    /** variance */
    value_type variance() const
    {
        return numeric::detail::variance(acc_);
    }


    /** standard error of mean value */
    value_type error_of_mean() const
    {
        return numeric::detail::error_of_mean(acc_);
    }

    /** return the number of samples aquired since last reset() */
    typename accumulator_type::size_type count() const
    {
        return numeric::detail::count(acc_);
    }

    /** reset the sum and number of samples to zero. */
    void reset()
    {
        LOG_TRACE("reset accumulator");
        acc_.reset();
    }

    /** Lua bindings */
    static void luaopen(lua_State* L);

protected:
    /** acquire function yielding the data to be accumulated */
    sample_function_type sample_;
    /** accumulator instance */
    accumulator_type acc_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};


} // namespace observables
} // namespace utility
} // namespace halmd


#endif /* ! HALMD_OBSERVABLES_UTILITY_ACCUMULATOR_HPP */
