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

#ifndef HALMD_IO_TRAJECTORY_READER_HPP
#define HALMD_IO_TRAJECTORY_READER_HPP

#include <halmd/mdsim/samples/host/trajectory.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/options.hpp>

namespace halmd
{
namespace io { namespace trajectory
{

template <int dimension>
class reader
{
public:
    // module definitions
    typedef reader _Self;
    static void options(po::options_description& desc);
    static void depends() {}
    static void select(po::variables_map const& vm);

    reader(modules::factory& factory, po::variables_map const& vm);
    virtual ~reader() {}

protected:
    /** absolute path to HDF5 trajectory file */
    std::string const path_;
    /** offset relative to start (non-negative) or end (negative) of dataset */
    ssize_t const offset_;
};

}} // namespace io::trajectory

} // namespace halmd

#endif /* ! HALMD_IO_TRAJECTORY_READER_HPP */
