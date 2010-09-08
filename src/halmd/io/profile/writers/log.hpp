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

#ifndef HALMD_IO_PROFILE_LOG_HPP
#define HALMD_IO_PROFILE_LOG_HPP

#include <string>
#include <utility>
#include <vector>

#include <halmd/io/profile/writer.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace io { namespace profile { namespace writers
{

/**
 * This module writes runtime accumulator results to the log.
 */
class log
  : public profile::writer
{
public:
    // module definitions
    typedef log _Self;
    typedef profile::writer _Base;
    static void options(po::options_description& desc) {}
    static void depends() {}
    static void select(po::options const& vm) {}

    typedef _Base::accumulator_type accumulator_type;
    typedef std::pair<accumulator_type const*, std::string> acc_desc_pair_type;

    log(modules::factory& factory, po::options const& vm);
    void write();

private:
    /**
    * register runtime accumulator
    */
    void register_accumulator(
        std::vector<std::string> const& tag
      , accumulator_type const& acc
      , std::string const& desc
    )
    {
        accumulators_.push_back(make_pair(&acc, desc));
    }

    /** list of registered accumulators and their descriptions */
    std::vector<acc_desc_pair_type> accumulators_;
};

}}} // namespace io::profile::writers

} // namespace halmd

#endif /* ! HALMD_IO_PROFILE_LOG_HPP */
