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

#include <boost/foreach.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/io/trajectory/writers/hdf5.hpp>

using namespace boost;
using namespace boost::filesystem;
using namespace std;

namespace halmd
{
namespace io { namespace trajectory { namespace writers
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::resolve(po::options const& vm)
{
    module<sample_type>::required(vm);
}

/**
 * read sample from HDF5 trajectory file
 */
template <int dimension, typename float_type>
hdf5<dimension, float_type>::hdf5(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , sample(module<sample_type>::fetch(vm))
  // initialize parameters
  , path_(initial_path() / (vm["output"].as<string>() + extension()))
  , file_(path_.file_string(), H5F_ACC_TRUNC)
{
    LOG("write trajectories to file: " << path_.file_string());

    H5::Group root = file_.createGroup("trajectory");

    for (size_t i = 0; i < sample->r.size(); ++i) {
        H5::Group type;
        if (sample->r.size() > 1) {
            type = root.createGroup(string(1, 'A' + i));
        }
        else {
            type = root;
        }
        H5::DataSet r = create_vector_dataset(type, "position", sample->r[i]);
        H5::DataSet v = create_vector_dataset(type, "velocity", sample->v[i]);

        // We bind the functions to write the datasets, using a
        // *reference* to the sample vector pointer so it remains
        // valid after reallocation, and a *copy* of the HDF5
        // dataset instance which goes out of scope at the end of
        // this loop.

        // particle positions
        writer_.insert(make_pair(
            H5xx::path(r)
          , bind(&hdf5<dimension, float_type>::write_vector_dataset, this, r, ref(sample->r[i]))
        ));
        // particle velocities
        writer_.insert(make_pair(
            H5xx::path(v)
          , bind(&hdf5<dimension, float_type>::write_vector_dataset, this, v, ref(sample->v[i]))
        ));
    }

    // simulation time in reduced units
    double* _FIXME_time = new double(42);
    H5::DataSet t = create_scalar_dataset(root, "time", *_FIXME_time);
    writer_.insert(make_pair(
        H5xx::path(t)
      , bind(&hdf5<dimension, float_type>::write_scalar_dataset, this, t, ref(*_FIXME_time))
    ));
}

/**
 * append samples to datasets
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::append()
{
    BOOST_FOREACH( typename writer_map::value_type const& writer, writer_ ) {
        LOG_DEBUG("writing dataset " << writer.first);
        writer.second();
    }
}

/**
 * flush HDF5 file
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::flush()
{
    file_.flush(H5F_SCOPE_GLOBAL);
}

/**
 * create vector sample dataset
 */
template <int dimension, typename float_type>
H5::DataSet hdf5<dimension, float_type>::create_vector_dataset(
    H5::Group where
  , std::string const& name
  , sample_vector_ptr sample
  )
{
    H5::DataType const tid = H5xx::ctype<float_type>();
    size_t const size = sample->size();

    // vector sample file dataspace
    hsize_t dim[3] = { 0, size, dimension };
    hsize_t max_dim[3] = { H5S_UNLIMITED, size, dimension };
    H5::DataSpace ds(3, dim, max_dim);

    H5::DSetCreatPropList cparms;
    hsize_t chunk_dim[3] = { 1, size, dimension };
    cparms.setChunk(3, chunk_dim);
    // enable GZIP compression
    cparms.setDeflate(6);

    return where.createDataSet(name, tid, ds, cparms);
}

/**
 * create scalar sample dataset
 */
template <int dimension, typename float_type>
H5::DataSet hdf5<dimension, float_type>::create_scalar_dataset(
    H5::Group where
  , std::string const& name
  , float_type sample
  )
{
    H5::DataType const tid = H5xx::ctype<double>();

    // scalar sample file dataspace
    hsize_t dim[1] = { 0 };
    hsize_t max_dim[1] = { H5S_UNLIMITED };
    H5::DataSpace ds(1, dim, max_dim);

    H5::DSetCreatPropList cparms;
    hsize_t chunk_dim[1] = { 1 };
    cparms.setChunk(1, chunk_dim);

    return where.createDataSet(name, tid, ds, cparms);
}

/**
 * write vector sample dataset
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::write_vector_dataset(
    H5::DataSet dset
  , sample_vector_ptr sample
  )
{
    H5::DataType const tid = H5xx::ctype<float_type>();
    size_t const size = sample->size();

    // extend vector sample file dataspace
    H5::DataSpace ds(dset.getSpace());
    hsize_t dim[3];
    ds.getSimpleExtentDims(dim);
    hsize_t count[3]  = { 1, 1, 1 };
    hsize_t start[3]  = { dim[0], 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, size, dimension };
    dim[0]++;
    ds.setExtentSimple(3, dim);
    ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // vector sample memory dataspace
    hsize_t dim_sample[2] = { size, dimension };
    H5::DataSpace ds_sample(2, dim_sample);

    dset.extend(dim);
    dset.write(sample->data(), tid, ds_sample, ds);
}

/**
 * write scalar sample dataset
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::write_scalar_dataset(
    H5::DataSet dset
  , double sample
  )
{
    H5::DataType const tid = H5xx::ctype<double>();

    // extend scalar sample file dataspace
    H5::DataSpace ds(dset.getSpace());
    hsize_t dim[1];
    ds.getSimpleExtentDims(dim);
    hsize_t count[1]  = { 1 };
    hsize_t start[1]  = { dim[0] };
    hsize_t stride[1] = { 1 };
    hsize_t block[1]  = { 1 };
    dim[0]++;
    ds.setExtentSimple(1, dim);
    ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    dset.extend(dim);
    dset.write(&sample, tid, H5S_SCALAR, ds);
}

// explicit instantiation
template class hdf5<3, double>;
template class hdf5<3, float>;
template class hdf5<2, double>;
template class hdf5<2, float>;

}}} // namespace io::trajectory::writers

template class module<io::trajectory::writers::hdf5<3, double> >;
template class module<io::trajectory::writers::hdf5<3, float> >;
template class module<io::trajectory::writers::hdf5<2, double> >;
template class module<io::trajectory::writers::hdf5<2, float> >;

} // namespace halmd
