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
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace observables
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
template <int dimension>
class writer
{
public:
    // module definitions
    typedef writer _Self;
    typedef typename observables::thermodynamics<dimension>::vector_type vector_type;

    static void options(po::options_description& desc) {}
    static void depends() {}
    static void select(po::options const& vm) {}

    writer(modules::factory& factory, po::options const& vm) {}
    virtual ~writer() {}
    virtual void write() = 0;

protected:
    friend class observables::thermodynamics<dimension>;

    /**
     * register an observable for output in HDF5 file
     * as dataset named by 'tag' with description 'desc',
     * data are read from *value_ptr at each invocation of write()
     */
    virtual void register_observable(
        std::string const& tag
      , double const* value_ptr
      , std::string const& desc
    ) = 0;

    virtual void register_observable(
        std::string const& tag
      , vector_type const* value_ptr
      , std::string const& desc
    ) = 0;

    virtual void register_observable(
        std::string const& tag
      , std::vector<double> const* value_ptr
      , std::string const& desc
    ) = 0;

    virtual void register_observable(
        std::string const& tag
      , std::vector<vector_type> const* value_ptr
      , std::string const& desc
    ) = 0;
};

}} // namespace io::statevars

} // namespace halmd

#endif /* ! HALMD_IO_STATEVARS_WRITER_HPP */
