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

#ifndef HALMD_IO_TRAJECTORY_HDF5_READER_HPP
#define HALMD_IO_TRAJECTORY_HDF5_READER_HPP

#include <halmd/io/trajectory/reader.hpp>
#include <halmd/mdsim/samples/host/trajectory.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/util/H5xx.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace io { namespace trajectory { namespace readers
{

template <int dimension, typename float_type>
class hdf5
  : public trajectory::reader<dimension>
{
public:
    typedef trajectory::reader<dimension> _Base;
    typedef mdsim::samples::host::trajectory<dimension, float_type> sample_type;
    typedef typename sample_type::sample_vector sample_vector;
    typedef typename sample_type::sample_vector_ptr sample_vector_ptr;

    static void resolve(po::options const& vm);
    hdf5(po::options const& vm);

    shared_ptr<sample_type> sample;

private:
    size_t read(H5::DataSet dset, sample_vector_ptr sample);
    size_t read(H5::DataSet dset, float_type& sample);

    using _Base::path_;
    using _Base::offset_;
};

}}} // namespace io::trajectory::readers

} // namespace halmd

#endif /* ! HALMD_IO_TRAJECTORY_HDF5_READER_HPP */
