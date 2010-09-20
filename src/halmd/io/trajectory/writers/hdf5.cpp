/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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
#include <H5xx.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/io/trajectory/writers/hdf5.hpp>

using namespace boost;
using namespace boost::filesystem;
using namespace std;
using namespace H5;

namespace halmd
{
namespace io { namespace trajectory { namespace writers
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::depends()
{
    modules::depends<_Self, sample_type>::required();
}

/**
 * read sample from HDF5 trajectory file
 */
template <int dimension, typename float_type>
hdf5<dimension, float_type>::hdf5(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  // dependency injection
  , sample(modules::fetch<sample_type>(factory, vm))
  // initialize parameters
  , path_(initial_path() / (vm["output"].as<string>() + file_extension()))
  , file_(path_.file_string(), H5F_ACC_TRUNC)
{
    LOG("write trajectory to file: " << path_.file_string());

    // create parameter group
    H5::Group param = open_group(file_, "/param");

    // store file version
    array<unsigned char, 2> version = {{ 1, 0 }};
    attribute(param, "file_version") = version;

    // create trajectory group
    H5::Group root = file_.createGroup("trajectory");

    for (size_t i = 0; i < sample->r.size(); ++i) {
        H5::Group type;
        if (sample->r.size() > 1) {
            type = root.createGroup(string(1, 'A' + i));
        }
        else {
            type = root;
        }
        // FIXME fixed_vector<> is not converted to boost::array
        typedef boost::array<
            typename sample_vector_type::value_type::value_type
          , sample_vector_type::value_type::static_size
        > array_type;
        typedef std::vector<array_type> output_type;
        size_t size = sample->r[i]->size();
        H5::DataSet r = create_dataset<output_type>(type, "position", size);
        H5::DataSet v = create_dataset<output_type>(type, "velocity", size);

        // particle positions
        writers_.push_back(
            make_dataset_writer(r, reinterpret_cast<output_type*>(&*sample->r[i]))
        );
        // particle velocities
        writers_.push_back(
            make_dataset_writer(v, reinterpret_cast<output_type*>(&*sample->v[i]))
        );
    }

    // simulation time in reduced units
    H5::DataSet t = create_dataset<double>(root, "time");
    writers_.push_back(make_dataset_writer(t, &sample->time));
}

/**
 * append samples to datasets
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::append()
{
    sample->acquire();
    BOOST_FOREACH (boost::function<void ()> const& writer, writers_) {
        writer();
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
