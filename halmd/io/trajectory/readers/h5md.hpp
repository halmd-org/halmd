/*
 * Copyright © 2008-2011  Peter Colberg
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

#ifndef HALMD_IO_TRAJECTORY_READERS_H5MD_HPP
#define HALMD_IO_TRAJECTORY_READERS_H5MD_HPP

#include <boost/optional.hpp>
#include <lua.hpp>

#include <halmd/io/trajectory/reader.hpp>
#include <halmd/observables/host/samples/trajectory.hpp>

namespace halmd
{
namespace io { namespace trajectory { namespace readers
{

/**
 * H5MD trajectory file reader
 */
template <int dimension, typename float_type>
class h5md
  : public trajectory::reader<dimension>
{
public:
    typedef trajectory::reader<dimension> _Base;
    typedef observables::host::samples::trajectory<dimension, float_type> sample_type;
    typedef typename sample_type::sample_vector sample_vector_type;
    typedef typename sample_type::sample_vector_ptr sample_vector_ptr;

    static void luaopen(lua_State* L);
    static boost::optional<H5::H5File> format(std::string const& file_name);

    h5md(
        boost::shared_ptr<sample_type> sample
      , std::string const& file_name
      , ssize_t offset
    );

    boost::shared_ptr<sample_type> sample;

private:
    /** absolute path to HDF5 trajectory file */
    std::string const path_;
    /** offset relative to start (non-negative) or end (negative) of dataset */
    ssize_t const offset_;
};

}}} // namespace io::trajectory::readers

} // namespace halmd

#endif /* ! HALMD_IO_TRAJECTORY_READERS_H5MD_HPP */
