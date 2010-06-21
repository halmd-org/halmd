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

#ifndef HALMD_IO_TRAJECTORY_WRITER_HPP
#define HALMD_IO_TRAJECTORY_WRITER_HPP

#include <halmd/mdsim/samples/host/trajectory.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace io { namespace trajectory
{

template <int dimension>
class writer
{
public:
    // module definitions
    typedef writer _Self;
    static void options(po::options_description& desc) {}
    static void depends() {}
    static void select(po::options const& vm) {}

    writer(modules::factory& factory, po::options const& vm) {}
    virtual ~writer() {}
    virtual std::string extension() const = 0;
    virtual void append() = 0;
    virtual void flush() = 0;
};

}} // namespace io::trajectory

} // namespace halmd

#endif /* ! HALMD_IO_TRAJECTORY_WRITER_HPP */
