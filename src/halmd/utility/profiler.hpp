/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_UTILITY_PROFILER_HPP
#define HALMD_UTILITY_PROFILER_HPP

#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/map.hpp>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{

namespace utility
{

/**
 * This module delegates registrations of runtime accumulator
 * maps of other modules to available profile writer modules.
 */
class profiler
{
public:
    // module definitions
    typedef profiler _Self;
    static void options(po::options_description& desc) {}
    static void depends() {}
    static void select(po::options const& vm) {}

    profiler(modules::factory& factory, po::options const& vm) {}

    template <typename AccumulatorMap>
    void register_map(AccumulatorMap const& map)
    {
        // FIXME
    }
};

} // namespace utility

} // namespace halmd

#endif /* ! HALMD_UTILITY_PROFILER_HPP */
