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

#include <halmd/deprecated/util/H5xx.hpp>
#include <halmd/io/statevars/writer.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace io { namespace statevars { namespace writers
{

/**
 * Write results for macroscopic state variables to a HDF5 file.
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
    void register_scalar_observable(
        std::string const& tag
      , double const* value_ptr
      , std::string const& desc
    );

    static void write_scalar_observable(
        H5::DataSet const& dset
      , double const* value_ptr
    );

    void register_vector_observable(
        std::string const& tag
      , vector_type const* value_ptr
      , std::string const& desc
    );

    static void write_vector_observable(
        H5::DataSet const& dset
      , vector_type const* value_ptr
    );

    H5::H5File file_;
    std::vector<writer_functor> writer_;
};

}}} // namespace io::statevars::writers

} // namespace halmd

#endif /* ! HALMD_IO_STATEVARS_HDF5_HPP */
