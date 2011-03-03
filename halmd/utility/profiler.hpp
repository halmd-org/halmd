/*
 * Copyright © 2010-2011  Peter Colberg
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

#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/io/profiling/writer.hpp>
#include <halmd/numeric/accumulator.hpp>

namespace halmd
{
namespace utility
{

class profiler
{
public:
    typedef io::profiling::writer writer_type;
    typedef writer_type::tag_type tag_type;
    typedef writer_type::accumulator_type accumulator_type;
    typedef std::vector<boost::shared_ptr<writer_type> > writers_type;

    static void luaopen(lua_State* L);
    profiler(writers_type writer, tag_type const& tag);
    void register_runtime(accumulator_type const& runtime, std::string const& desc);

private:
    writers_type writer_;
    tag_type tag_;
};

} // namespace utility

} // namespace halmd

#endif /* ! HALMD_UTILITY_PROFILER_HPP */
