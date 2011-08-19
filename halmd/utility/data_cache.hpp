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

#ifndef HALMD_UTILITY_DATA_CACHE_HPP
#define HALMD_UTILITY_DATA_CACHE_HPP

#include <boost/shared_ptr.hpp>
#include <exception>
#include <limits>

#include <halmd/mdsim/clock.hpp>

namespace halmd {

/**
 * Simple clock-based data cache for a single value.
 *
 * Assignment updates the data cache and the internal time stamp, which is
 * derived from the simulation clock. Cache coherency is tested with valid()
 * and the conversion operator yields the cached data, if valid.
 */
template <typename T>
class data_cache
{
public:
    typedef mdsim::clock clock_type;

    /** Construct cache and store simulation clock */
    data_cache(boost::shared_ptr<clock_type const> clock)
      : clock_(clock), step_(std::numeric_limits<step_type>::max())
    {}

    /** Returns true if internal time stamp and simulation clock agree */
    bool valid() const
    {
        return step_ == clock_->step();
    }

    /** Cache passed data by copy and update time stamp from simulation clock. */
    template <typename S>
    T const& operator=(S const& data)
    {
        data_ = data;
        step_ = clock_->step();
        return data_;
    }

    /** Returns cached data. Throws on invalid time stamp. */
    operator T const&() const
    {
        if (!valid()) {
            throw std::logic_error("Reading data from outdated cache.");
        }
        return data_;
    }

private:
    typedef clock_type::step_type step_type;

    boost::shared_ptr<clock_type const> clock_;

    T data_;
    step_type step_;
};

} // namespace halmd

#endif /* ! HALMD_UTILITY_DATA_CACHE_HPP */
