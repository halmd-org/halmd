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
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/map.hpp>
#include <lua.hpp>

#include <halmd/numeric/accumulator.hpp>

namespace halmd
{
namespace io { namespace profiling
{

// forward declaration
class writer;

}} // namespace io::profiling

namespace utility
{

namespace detail { namespace profiler
{

// forward declaration
struct visitor;

}} // namespace detail::profiler


/**
 * This module delegates registrations of runtime accumulator
 * maps of other modules to available profiling writer modules.
 */
class profiler
{
public:
    typedef io::profiling::writer profiling_writer_type;
    typedef accumulator<double> accumulator_type;

    std::vector<boost::shared_ptr<profiling_writer_type> > profiling_writers;

    static void luaopen(lua_State* L);

    profiler(std::vector<boost::shared_ptr<profiling_writer_type> > profiling_writers);

    template <typename AccumulatorMap>
    void register_map(AccumulatorMap const& map)
    {
        boost::fusion::for_each(map, detail::profiler::visitor(*this));
    }

private:
    friend class detail::profiler::visitor;
    void register_accumulator(
        std::type_info const& tag
      , accumulator_type const& acc
      , std::string const& desc
    ) const;
};

/**
 * Define tag for runtime accumulator within module.
 */
#define HALMD_PROFILING_TAG(__tag__, __desc__)    \
    struct __tag__ {                            \
        static std::string desc() {             \
            return __desc__;                    \
        }                                       \
    }

namespace detail { namespace profiler
{

/**
 * Extract tag type and accumulator reference from boost::fusion::pair.
 */
struct visitor
{
    template <typename TaggedAccumulator>
    void operator()(TaggedAccumulator const& acc) const
    {
        p.register_accumulator(
            typeid(typename TaggedAccumulator::first_type)
          , acc.second
          , TaggedAccumulator::first_type::desc()
        );
    }
    visitor(utility::profiler const& p) : p(p) {}
    utility::profiler const& p;
};

}} // namespace detail::profiler

} // namespace utility

} // namespace halmd

#endif /* ! HALMD_UTILITY_PROFILER_HPP */
