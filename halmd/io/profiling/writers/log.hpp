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

#ifndef HALMD_IO_PROFILING_LOG_HPP
#define HALMD_IO_PROFILING_LOG_HPP

#include <lua.hpp>
#include <utility>
#include <vector>

#include <halmd/io/profiling/writer.hpp>

namespace halmd
{
namespace io { namespace profiling { namespace writers
{

/**
 * This module writes runtime accumulator results to the log.
 */
class log
  : public profiling::writer
{
public:
    typedef profiling::writer _Base;
    typedef _Base::tag_type tag_type;
    typedef _Base::accumulator_type accumulator_type;
    typedef std::pair<accumulator_type const*, std::string> acc_desc_pair_type;

    static void luaopen(lua_State* L);

    log() {}
    virtual void write();

    /**
    * register runtime accumulator
    */
    virtual void register_accumulator(
        tag_type const& tag
      , accumulator_type const& acc
      , std::string const& desc
    )
    {
        accumulators_.push_back(make_pair(&acc, desc));
    }

private:
    /** list of registered accumulators and their descriptions */
    std::vector<acc_desc_pair_type> accumulators_;
};

}}} // namespace io::profiling::writers

} // namespace halmd

#endif /* ! HALMD_IO_PROFILING_LOG_HPP */
