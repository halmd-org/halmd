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

#ifndef HALMD_IO_STATEVARS_HDF5_HPP
#define HALMD_IO_STATEVARS_HDF5_HPP

#include <H5xx.hpp>

#include <halmd/io/statevars/writer.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace io { namespace statevars { namespace writers
{

/**
 * Write results for macroscopic state variables to an HDF5 file.
 */
template <int dimension>
class hdf5
  : public statevars::writer<dimension>
{
public:
    // module definitions
    typedef hdf5 _Self;
    typedef statevars::writer<dimension> _Base;
    typedef typename _Base::vector_type vector_type;

    static void options(po::options_description& desc) {}
    static void depends() {}
    static void select(po::options const& vm) {}

    typedef boost::function<void ()> writer_functor;

    hdf5(modules::factory& factory, po::options const& vm);
    void write();

private:
    // virtual register function, only called by base class
    virtual void register_observable(
        std::string const& tag
      , void const* value_ptr
      , std::type_info const& value_type
      , std::string const& desc
    );

    // templates for register functions
    template <typename T>
    void register_observable(
        std::string const& tag
      , T const* value_ptr
      , std::string const& desc
    );

    template <typename T>
    void register_observable(
        std::string const& tag
      , std::vector<T> const* value_ptr
      , std::string const& desc
    );

    // virtual write_dataset function, only called by base class
    virtual void write_dataset(
        std::string const& tag
      , void const* value
      , std::type_info const& value_type
      , std::string const& desc
    );

    // templates for write_dataset function
    template <typename T>
    void write_dataset(
        std::string const& tag
      , T const& value
      , std::string const& desc
    );

    template <typename T>
    void write_dataset(
        std::string const& tag
      , std::vector<T> const& value
      , std::string const& desc
    );

    /** HDF5 file object */
    H5::H5File file_;
    /** list of registered writer functors */
    std::vector<writer_functor> writer_;
};

}}} // namespace io::statevars::writers

} // namespace halmd

#endif /* ! HALMD_IO_STATEVARS_HDF5_HPP */
