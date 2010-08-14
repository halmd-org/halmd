/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_SCRIPT_HPP
#define HALMD_SCRIPT_HPP

#include <halmd/io/profile/writer.hpp>
#include <halmd/main.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/utility/module.hpp>

namespace halmd
{

/**
 * NVE ensemble run
 */
template <int dimension>
class script
  : public halmd::main
{
public:
    // module definitions
    typedef script _Self;
    typedef halmd::main _Base;
    static void depends();
    static void options(po::options_description& desc);
    static void select(po::options const& vm) {}

    typedef mdsim::core<dimension> core_type;
    typedef io::profile::writer profile_writer_type;

    script(modules::factory& factory, po::options const& vm);
    virtual ~script() {}
    virtual void run();
    uint64_t steps() { return steps_; }
    double time() { return time_; }

    shared_ptr<core_type> core;
    std::vector<shared_ptr<profile_writer_type> > profile_writers;

protected:
    /** number of integration steps */
    uint64_t steps_;
    /** integration time in MD units */
    double time_;
};

} // namespace halmd

#endif /* ! HALMD_SCRIPT_HPP */
