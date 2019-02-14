/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE trajectory
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <memory>
#include <vector>

#include <halmd/io/readers/h5md/append.hpp>
#include <halmd/io/readers/h5md/file.hpp>
#include <halmd/io/writers/h5md/append.hpp>
#include <halmd/io/writers/h5md/file.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/observables/host/samples/sample.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>

boost::array<std::string, 3> const types = {{ "A", "B", "C" }};

template <typename position_sample_type, typename velocity_sample_type, typename writer_type>
void on_write_sample(std::vector<std::shared_ptr<position_sample_type> > const& position_samples
                   , std::vector<std::shared_ptr<velocity_sample_type> > const& velocity_samples
                   , std::shared_ptr<writer_type> writer)
{
    typedef typename position_sample_type::array_type position_array_type;
    typedef typename velocity_sample_type::array_type velocity_array_type;
    typedef typename writer_type::subgroup_type subgroup_type;

    for (unsigned int type = 0; type < position_samples.size(); ++type) {
        std::shared_ptr<position_sample_type> sample = position_samples[type];
        {
            subgroup_type group;
            writer->template on_write<position_array_type const&>(
                group
              , [=]() -> position_array_type const& {
                    return sample->data();
                }
              , {types[type], "position"}
            );
            BOOST_CHECK_EQUAL(h5xx::path(group), "/trajectory/" + types[type] + "/position");
        }
    }
    for (unsigned int type = 0; type < velocity_samples.size(); ++type) {
        std::shared_ptr<velocity_sample_type> sample = velocity_samples[type];
        {
            subgroup_type group;
            writer->template on_write<velocity_array_type const&>(
                group
              , [=]() -> velocity_array_type const& {
                    return sample->data();
                }
              , {types[type], "velocity"}
            );
            BOOST_CHECK_EQUAL(h5xx::path(group), "/trajectory/" + types[type] + "/velocity");
        }
    }
}

template <typename position_sample_type, typename velocity_sample_type, typename reader_type>
void on_read_sample(std::vector<std::shared_ptr<position_sample_type> > const& position_samples
                  , std::vector<std::shared_ptr<velocity_sample_type> > const& velocity_samples
                  , std::shared_ptr<reader_type> reader)
{
    typedef typename position_sample_type::data_type position_type;
    typedef typename velocity_sample_type::data_type velocity_type;
    typedef typename reader_type::subgroup_type subgroup_type;

    for (unsigned int type = 0; type < position_samples.size(); ++type) {
        std::shared_ptr<position_sample_type> position_sample = position_samples[type];
        {
            std::shared_ptr<std::vector<position_type>> array = std::make_shared<std::vector<position_type>>();
            subgroup_type group;
            reader->template on_read<std::vector<position_type>&>(
                group
              , [=]() -> std::vector<position_type>& {
                    return *array;
                }
              , {types[type], "position"}
            );
            reader->on_append_read([=]() {
                std::copy(
                    array->begin()
                  , array->end()
                  , position_sample->data().begin()
                );
            });
            BOOST_CHECK_EQUAL(h5xx::path(group), "/trajectory/" + types[type] + "/position");
        }
    }
    for (unsigned int type = 0; type < velocity_samples.size(); ++type) {
        std::shared_ptr<velocity_sample_type> velocity_sample = velocity_samples[type];

        {
            std::shared_ptr<std::vector<velocity_type>> array = std::make_shared<std::vector<velocity_type>>();
            subgroup_type group;
            reader->template on_read<std::vector<velocity_type>&>(
                group
              , [=]() -> std::vector<velocity_type>& {
                    return *array;
                }
              , {types[type], "velocity"}
            );
            reader->on_append_read([=]() {
                std::copy(
                    array->begin()
                  , array->end()
                  , velocity_sample->data().begin()
                );
            });
            BOOST_CHECK_EQUAL(h5xx::path(group), "/trajectory/" + types[type] + "/velocity");
        }
    }
}

template <int dimension>
void h5md(std::vector<unsigned int> const& ntypes)
{
    typedef halmd::observables::host::samples::sample<dimension, float> float_position_sample_type;
    typedef halmd::observables::host::samples::sample<dimension, float> float_velocity_sample_type;
    typedef halmd::observables::host::samples::sample<dimension, double> double_position_sample_type;
    typedef halmd::observables::host::samples::sample<dimension, double> double_velocity_sample_type;
    typedef typename double_position_sample_type::array_type double_position_array_type;
    typedef typename double_velocity_sample_type::array_type double_velocity_array_type;

    typedef halmd::fixed_vector<float, dimension> float_vector_type;
    typedef halmd::fixed_vector<double, dimension> double_vector_type;

    std::string filename("test_io_h5md_trajectory_" + std::to_string(dimension) + "d_single" + std::to_string (ntypes.size()) + ".trj");

    BOOST_TEST_MESSAGE("Testing " << ntypes.size() << " particle types");

    // construct phase space sample and fill with positions and velocities
    std::vector<std::shared_ptr<double_position_sample_type> > double_position_sample;
    std::vector<std::shared_ptr<double_velocity_sample_type> > double_velocity_sample;
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        double_position_sample.push_back(std::make_shared<double_position_sample_type>(ntypes[type]));
        double_velocity_sample.push_back(std::make_shared<double_velocity_sample_type>(ntypes[type]));
        double_position_array_type& r_sample = double_position_sample[type]->data();
        double_velocity_array_type& v_sample = double_velocity_sample[type]->data();
        for (unsigned int i = 0; i < ntypes[type]; ++i) {
            double_vector_type& r = r_sample[i];
            r[0] = type;
            r[1] = 1. / (i + 1);
            if (dimension > 2) {
                r[2] = -i;
            }

            double_vector_type& v = v_sample[i];
            v[0] = i + 1;
            v[1] = sqrt(i + 1);
            if (dimension > 2) {
                v[2] = 1L << (i % 64);
            }
        }
    }

    // copy sample to single precision
    std::vector<std::shared_ptr<float_position_sample_type> > float_position_sample;
    std::vector<std::shared_ptr<float_velocity_sample_type> > float_velocity_sample;
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        float_position_sample.push_back(std::make_shared<float_position_sample_type>(ntypes[type]));
        float_velocity_sample.push_back(std::make_shared<float_velocity_sample_type>(ntypes[type]));
        std::transform(
            double_position_sample[type]->data().begin()
          , double_position_sample[type]->data().end()
          , float_position_sample[type]->data().begin()
          , [](double_vector_type const& r) {
                return float_vector_type(r);
            }
        );
        std::transform(
            double_velocity_sample[type]->data().begin()
          , double_velocity_sample[type]->data().end()
          , float_velocity_sample[type]->data().begin()
          , [](double_vector_type const& v) {
                return float_vector_type(v);
            }
        );
    }

    // write single-precision sample to file
    // use time-step not exactly representable as float-point value
    std::shared_ptr<halmd::mdsim::clock> clock = std::make_shared<halmd::mdsim::clock>();
    std::shared_ptr<halmd::io::writers::h5md::file> writer_file =
        std::make_shared<halmd::io::writers::h5md::file>(filename, "", "", true);
    std::shared_ptr<halmd::io::writers::h5md::append> writer =
        std::make_shared<halmd::io::writers::h5md::append>(writer_file->root(), std::vector<std::string>{"trajectory"}, clock);

    on_write_sample(float_position_sample, float_velocity_sample, writer);

    writer->write();
    writer_file->flush();

    // overwrite with double-precision data,
    // resetting the shared_ptr first closes the HDF5 file
    writer.reset();
    writer_file.reset();

    filename = std::string ("test_io_h5md_trajectory_" + std::to_string(dimension) + "d_double" + std::to_string (ntypes.size()) + ".trj");
    writer_file = std::make_shared<halmd::io::writers::h5md::file>(filename, "", "", true);
    writer = std::make_shared<halmd::io::writers::h5md::append>(writer_file->root(), std::vector<std::string>{"trajectory"}, clock);

    on_write_sample(double_position_sample, double_velocity_sample, writer);

    writer->write();
    writer_file->flush();

    // simulate an integration step for the very first particle
    clock->set_timestep(1 / 6.);
    clock->advance();
    double_position_sample[0]->data()[0] += double_velocity_sample[0]->data()[0];
    double_velocity_sample[0]->data()[0] = double_vector_type(sqrt(2));

    // deconstruct the file module before writing
    // the HDF5 library will keep the file open as long as groups
    // or datasets are open, which is the case for the writer
    writer_file.reset();

    writer->write();

    // deconstruct the writer to flush and close the HDF5 file
    writer.reset();

    // test integrity of H5MD file
    BOOST_CHECK(halmd::io::readers::h5md::file::check(filename));

    // read phase space sample #1 from file in double precision
    // reading is done upon construction, so we use an unnamed, temporary reader object
    // allocate memory for reading back the phase space sample
    std::vector<std::shared_ptr<double_position_sample_type> > double_position_sample_;
    std::vector<std::shared_ptr<double_velocity_sample_type> > double_velocity_sample_;
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        double_position_sample_.push_back(std::make_shared<double_position_sample_type>(ntypes[type]));
        double_velocity_sample_.push_back(std::make_shared<double_velocity_sample_type>(ntypes[type]));
    }

    std::shared_ptr<halmd::io::readers::h5md::file> reader_file =
        std::make_shared<halmd::io::readers::h5md::file>(filename);
    std::shared_ptr<halmd::io::readers::h5md::append> reader =
        std::make_shared<halmd::io::readers::h5md::append>(reader_file->root(), std::vector<std::string>{"trajectory"});

    on_read_sample(double_position_sample_, double_velocity_sample_, reader);

    // read at time 1/6. with maximum tolerated rounding error of 100 × 1/6. × epsilon
    reader->read_at_time(0.16666666666667);

    // check binary equality of written and read data
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        BOOST_CHECK_EQUAL_COLLECTIONS(
            double_position_sample_[type]->data().begin()
          , double_position_sample_[type]->data().end()
          , double_position_sample[type]->data().begin()
          , double_position_sample[type]->data().end()
        );
        BOOST_CHECK_EQUAL_COLLECTIONS(
            double_velocity_sample_[type]->data().begin()
          , double_velocity_sample_[type]->data().end()
          , double_velocity_sample[type]->data().begin()
          , double_velocity_sample[type]->data().end()
        );
    }

    // read phase space sample #0 from file in single precision
    std::vector<std::shared_ptr<float_position_sample_type> > float_position_sample_;
    std::vector<std::shared_ptr<float_velocity_sample_type> > float_velocity_sample_;
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        float_position_sample_.push_back(std::make_shared<float_position_sample_type>(ntypes[type]));
        float_velocity_sample_.push_back(std::make_shared<float_velocity_sample_type>(ntypes[type]));
    }

    // reconstruct the reader to replace slots to double with float sample
    reader.reset();
    reader = std::make_shared<halmd::io::readers::h5md::append>(reader_file->root(), std::vector<std::string>{"trajectory"});

    // deconstruct file module to check that the HDF5 library
    // keeps the file open if reader module still exists
    reader_file.reset();

    on_read_sample(float_position_sample_, float_velocity_sample_, reader);

    reader->read_at_time(0);

    // check binary equality of written and read data,
    // note that float_sample was not modified and thus corresponds to #0
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        BOOST_CHECK_EQUAL_COLLECTIONS(
            float_position_sample_[type]->data().begin()
          , float_position_sample_[type]->data().end()
          , float_position_sample[type]->data().begin()
          , float_position_sample[type]->data().end()
        );
        BOOST_CHECK_EQUAL_COLLECTIONS(
            float_velocity_sample_[type]->data().begin()
          , float_velocity_sample_[type]->data().end()
          , float_velocity_sample[type]->data().begin()
          , float_velocity_sample[type]->data().end()
        );
    }

    // close and remove file
    reader.reset();
#ifdef NDEBUG
    remove(filename.c_str());
#endif
}

HALMD_TEST_INIT( trajectory )
{
    using namespace boost::unit_test;

    std::vector<std::vector<unsigned int> > ntypes = {{1}, {1, 10}, {1, 10, 100}};

    test_suite* ts1 = BOOST_TEST_SUITE( "2d" );
    ts1->add( BOOST_PARAM_TEST_CASE( &h5md<2>, ntypes.begin(), ntypes.end() ) );

    test_suite* ts2 = BOOST_TEST_SUITE( "3d" );
    ts2->add( BOOST_PARAM_TEST_CASE( &h5md<3>, ntypes.begin(), ntypes.end() ) );

    framework::master_test_suite().add( ts1 );
    framework::master_test_suite().add( ts2 );
}
