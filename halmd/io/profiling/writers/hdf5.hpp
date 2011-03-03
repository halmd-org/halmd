/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_IO_PROFILING_HDF5_HPP
#define HALMD_IO_PROFILING_HDF5_HPP

#include <lua.hpp>

#include <h5xx/h5xx.hpp>
#include <halmd/io/profiling/writer.hpp>

namespace halmd
{
namespace io { namespace profiling { namespace writers
{

/**
 * This module writes runtime accumulator results to an HDF5 file.
 */
class hdf5
  : public profiling::writer
{
public:
    typedef profiling::writer _Base;
    typedef _Base::tag_type tag_type;
    typedef _Base::accumulator_type accumulator_type;
    typedef boost::function<void ()> writer_functor;

    static void luaopen(lua_State* L);

    hdf5(std::string const& file_name);
    virtual void write();

    H5::H5File const& file() const
    {
        return file_;
    }

    virtual void register_accumulator(
        tag_type const& tag
      , accumulator_type const& acc
      , std::string const& desc
    );

private:
    static void write_accumulator(
        H5::DataSet const& dset
      , accumulator_type const& acc
    );

    H5::H5File file_;
    std::vector<writer_functor> writer_;
};

}}} // namespace io::profiling::writers

} // namespace halmd

#endif /* ! HALMD_IO_PROFILING_HDF5_HPP */
