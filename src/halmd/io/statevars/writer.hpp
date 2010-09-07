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

#ifndef HALMD_IO_STATEVARS_WRITER_HPP
#define HALMD_IO_STATEVARS_WRITER_HPP

#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim
{

// forward declaration
template <int dimension>
class thermodynamics;

} // namespace mdsim

namespace io { namespace statevars
{

/**
 * Abstract base class of a writer of macroscopic state variables.
 */
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
    virtual void write() = 0;

protected:
    template <int dimension>
    friend class mdsim::thermodynamics;

    virtual void register_observable(
        std::string const& tag
      , double const* value_ptr
      , std::string const& desc
    ) = 0;
};

}} // namespace io::statevars

} // namespace halmd

#endif /* ! HALMD_IO_STATEVARS_WRITER_HPP */
