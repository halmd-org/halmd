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

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/options.hpp>

namespace halmd
{
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
    typedef typename mdsim::type_traits<dimension, double>::vector_type vector_type;

    static void options(po::options_description& desc) {}
    static void depends() {}
    static void select(po::variables_map const& vm) {}

    writer(modules::factory& factory, po::variables_map const& vm) {}
    virtual ~writer() {}
    virtual void write() = 0;

    /**
     * Register an observable for output as dataset named by 'tag'
     * with description 'desc'. Data are read from *value_ptr at
     * each invocation of write().
     *
     * convenience template for the protected virtual function register_observable() */
    template <typename T>
    void register_observable(std::string const& tag, T const* value_ptr, std::string const& desc)
    {
        register_observable(tag, value_ptr, typeid(T), desc);
    }

    /**
     * Write single dataset 'value' named by 'tag' with description 'desc'
     *
     * convenience template for the protected virtual function write_dataset() */
    template <typename T>
    void write_dataset(std::string const& tag, T const& value, std::string const& desc)
    {
        write_dataset(tag, &value, typeid(T), desc);
    }

protected:
    /**
     * Register an observable for output as dataset named by 'tag'
     * with description 'desc'. Data are read from *value_ptr at
     * each invocation of write(). The value type is passed
     * separately to this pure virtual function
     */
    virtual void register_observable(
        std::string const& tag
      , void const* value_ptr
      , std::type_info const& value_type
      , std::string const& desc
    ) = 0;

    /**
     * Write single dataset named by 'tag' with description 'desc'.
     * Data are read from *value_ptr, the value type is passed
     * separately to this pure virtual function
     */
    virtual void write_dataset(
        std::string const& tag
      , void const* value_ptr
      , std::type_info const& value_type
      , std::string const& desc
    ) = 0;
};

}} // namespace io::statevars

} // namespace halmd

#endif /* ! HALMD_IO_STATEVARS_WRITER_HPP */
