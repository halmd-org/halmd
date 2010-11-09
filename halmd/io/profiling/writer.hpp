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

#ifndef HALMD_IO_PROFILING_WRITER_HPP
#define HALMD_IO_PROFILING_WRITER_HPP

#include <lua.hpp>
#include <string>
#include <vector>

#include <halmd/numeric/accumulator.hpp>

namespace halmd
{
namespace utility
{

// forward declaration
class profiler;

} // namespace utility

namespace io { namespace profiling
{

/**
 * Abstract base class of a profiler writer.
 */
class writer
{
public:
    typedef accumulator<double> accumulator_type;

    static void luaopen(lua_State* L);

    writer() {}
    virtual ~writer() {}
    virtual void write() = 0;

protected:
    friend class utility::profiler;
    virtual void register_accumulator(
        std::vector<std::string> const& tag
      , accumulator<double> const& acc
      , std::string const& desc
    ) = 0;
};

}} // namespace io::profiling

} // namespace halmd

#endif /* ! HALMD_IO_PROFILING_WRITER_HPP */
