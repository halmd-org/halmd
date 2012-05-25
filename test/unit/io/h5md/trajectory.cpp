/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE trajectory
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include <halmd/io/readers/h5md/append.hpp>
#include <halmd/io/readers/h5md/file.hpp>
#include <halmd/io/writers/h5md/append.hpp>
#include <halmd/io/writers/h5md/file.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>

using namespace boost;
using namespace boost::assign; // list_of
using namespace halmd;
using namespace halmd::io; // avoid ambiguity of io:: between halmd::io and boost::io
using namespace std;

boost::array<string, 3> const types = {{ "A", "B", "C" }};

template <typename sample_type, typename writer_type>
void on_write_sample(vector<boost::shared_ptr<sample_type> > const& samples, boost::shared_ptr<writer_type> writer)
{
    typedef typename sample_type::position_array_type position_array_type;
    typedef typename sample_type::velocity_array_type velocity_array_type;
    typedef typename writer_type::subgroup_type subgroup_type;

    for (unsigned int type = 0; type < samples.size(); ++type) {
        boost::shared_ptr<sample_type> sample = samples[type];
        {
            subgroup_type group;
            writer->template on_write<position_array_type const&>(
                group
              , [=]() -> position_array_type const& {
                    return sample->position();
                }
              , {types[type], "position"}
            );
            BOOST_CHECK_EQUAL(h5xx::path(group), "/trajectory/" + types[type] + "/position");
        }
        {
            subgroup_type group;
            writer->template on_write<velocity_array_type const&>(
                group
              , [=]() -> velocity_array_type const& {
                    return sample->velocity();
                }
              , {types[type], "velocity"}
            );
            BOOST_CHECK_EQUAL(h5xx::path(group), "/trajectory/" + types[type] + "/velocity");
        }
    }
}

template <typename sample_type, typename reader_type>
void on_read_sample(vector<boost::shared_ptr<sample_type> > const& samples, boost::shared_ptr<reader_type> reader)
{
    typedef typename sample_type::position_array_type position_array_type;
    typedef typename sample_type::velocity_array_type velocity_array_type;
    typedef typename sample_type::position_array_type::value_type position_type;
    typedef typename sample_type::velocity_array_type::value_type velocity_type;
    typedef typename reader_type::subgroup_type subgroup_type;

    for (unsigned int type = 0; type < samples.size(); ++type) {
        boost::shared_ptr<sample_type> sample = samples[type];
        {
            boost::shared_ptr<std::vector<position_type>> array = boost::make_shared<std::vector<position_type>>();
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
                  , sample->position().begin()
                );
            });
            BOOST_CHECK_EQUAL(h5xx::path(group), "/trajectory/" + types[type] + "/position");
        }
        {
            boost::shared_ptr<std::vector<velocity_type>> array = boost::make_shared<std::vector<velocity_type>>();
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
                  , sample->velocity().begin()
                );
            });
            BOOST_CHECK_EQUAL(h5xx::path(group), "/trajectory/" + types[type] + "/velocity");
        }
    }
}

template <int dimension>
void h5md(vector<unsigned int> const& ntypes)
{
    typedef observables::host::samples::phase_space<dimension, float> float_sample_type;
    typedef observables::host::samples::phase_space<dimension, double> double_sample_type;
    typedef typename float_sample_type::position_array_type float_position_array_type;
    typedef typename float_sample_type::velocity_array_type float_velocity_array_type;
    typedef typename double_sample_type::position_array_type double_position_array_type;
    typedef typename double_sample_type::velocity_array_type double_velocity_array_type;

    typedef fixed_vector<float, dimension> float_vector_type;
    typedef fixed_vector<double, dimension> double_vector_type;

    string const filename("test_io_h5md_trajectory_" + lexical_cast<string>(dimension) + "d.trj");

    BOOST_MESSAGE("Testing " << ntypes.size() << " particle types");

    // construct phase space sample and fill with positions and velocities
    vector<boost::shared_ptr<double_sample_type> > double_sample;
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        double_sample.push_back(boost::make_shared<double_sample_type>(ntypes[type]));
        double_position_array_type& r_sample = double_sample[type]->position();
        double_velocity_array_type& v_sample = double_sample[type]->velocity();
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
                v[2] = static_cast<double>(1L << (i % 64));
            }
        }
    }

    // copy sample to single precision
    vector<boost::shared_ptr<float_sample_type> > float_sample;
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        float_sample.push_back(boost::make_shared<float_sample_type>(ntypes[type]));
        transform(
            double_sample[type]->position().begin()
          , double_sample[type]->position().end()
          , float_sample[type]->position().begin()
          , [](double_vector_type const& r) {
                return float_vector_type(r);
            }
        );
        transform(
            double_sample[type]->velocity().begin()
          , double_sample[type]->velocity().end()
          , float_sample[type]->velocity().begin()
          , [](double_vector_type const& v) {
                return float_vector_type(v);
            }
        );
    }

    // write single-precision sample to file
    // use time-step not exactly representable as float-point value
    boost::shared_ptr<mdsim::clock> clock = make_shared<mdsim::clock>();
    boost::shared_ptr<writers::h5md::file> writer_file =
        boost::make_shared<writers::h5md::file>(filename);
    boost::shared_ptr<writers::h5md::append> writer =
        boost::make_shared<writers::h5md::append>(writer_file->root(), list_of("trajectory"), clock);

    on_write_sample(float_sample, writer);

    writer->write();
    writer_file->flush();

    // overwrite with double-precision data,
    // resetting the shared_ptr first closes the HDF5 file
    writer.reset();
    writer_file.reset();
    writer_file = boost::make_shared<writers::h5md::file>(filename);
    writer = boost::make_shared<writers::h5md::append>(writer_file->root(), list_of("trajectory"), clock);

    on_write_sample(double_sample, writer);

    writer->write();
    writer_file->flush();

    // simulate an integration step for the very first particle
    clock->set_timestep(1 / 6.);
    clock->advance();
    double_sample[0]->position()[0] += double_sample[0]->velocity()[0];
    double_sample[0]->velocity()[0] = double_vector_type(sqrt(2));

    // deconstruct the file module before writing
    // the HDF5 library will keep the file open as long as groups
    // or datasets are open, which is the case for the writer
    writer_file.reset();

    writer->write();

    // deconstruct the writer to flush and close the HDF5 file
    writer.reset();

    // test integrity of H5MD file
    BOOST_CHECK(readers::h5md::file::check(filename));

    // read phase space sample #1 from file in double precision
    // reading is done upon construction, so we use an unnamed, temporary reader object
    // allocate memory for reading back the phase space sample
    vector<boost::shared_ptr<double_sample_type> > double_sample_;
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        double_sample_.push_back(boost::make_shared<double_sample_type>(ntypes[type]));
    }

    boost::shared_ptr<readers::h5md::file> reader_file =
        boost::make_shared<readers::h5md::file>(filename);
    boost::shared_ptr<readers::h5md::append> reader =
        boost::make_shared<readers::h5md::append>(reader_file->root(), list_of("trajectory"));

    on_read_sample(double_sample_, reader);

    // read at time 1/6. with maximum tolerated rounding error of 100 × 1/6. × epsilon
    reader->read_at_time(0.16666666666667);

    // check binary equality of written and read data
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        BOOST_CHECK_EQUAL_COLLECTIONS(
            double_sample_[type]->position().begin()
          , double_sample_[type]->position().end()
          , double_sample[type]->position().begin()
          , double_sample[type]->position().end()
        );
        BOOST_CHECK_EQUAL_COLLECTIONS(
            double_sample_[type]->velocity().begin()
          , double_sample_[type]->velocity().end()
          , double_sample[type]->velocity().begin()
          , double_sample[type]->velocity().end()
        );
    }

    // read phase space sample #0 from file in single precision
    vector<boost::shared_ptr<float_sample_type> > float_sample_;
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        float_sample_.push_back(boost::make_shared<float_sample_type>(ntypes[type]));
    }

    // reconstruct the reader to replace slots to double with float sample
    reader.reset();
    reader = boost::make_shared<readers::h5md::append>(reader_file->root(), list_of("trajectory"));

    // deconstruct file module to check that the HDF5 library
    // keeps the file open if reader module still exists
    reader_file.reset();

    on_read_sample(float_sample_, reader);

    reader->read_at_time(0);

    // check binary equality of written and read data,
    // note that float_sample was not modified and thus corresponds to #0
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        BOOST_CHECK_EQUAL_COLLECTIONS(
            float_sample_[type]->position().begin()
          , float_sample_[type]->position().end()
          , float_sample[type]->position().begin()
          , float_sample[type]->position().end()
        );
        BOOST_CHECK_EQUAL_COLLECTIONS(
            float_sample_[type]->velocity().begin()
          , float_sample_[type]->velocity().end()
          , float_sample[type]->velocity().begin()
          , float_sample[type]->velocity().end()
        );
    }

    // close and remove file
    reader.reset();
#ifdef NDEBUG
    remove(filename.c_str());
#endif
}

HALMD_TEST_INIT( init_unit_test_suite )
{
    using namespace boost::unit_test;
    using namespace boost::unit_test::framework;

    vector<vector<unsigned int> > ntypes;
    ntypes += list_of(1), list_of(1)(10), list_of(1)(10)(100);

    test_suite* ts1 = BOOST_TEST_SUITE( "2d" );
    ts1->add( BOOST_PARAM_TEST_CASE( &h5md<2>, ntypes.begin(), ntypes.end() ) );

    test_suite* ts2 = BOOST_TEST_SUITE( "3d" );
    ts2->add( BOOST_PARAM_TEST_CASE( &h5md<3>, ntypes.begin(), ntypes.end() ) );

    master_test_suite().add( ts1 );
    master_test_suite().add( ts2 );
}
