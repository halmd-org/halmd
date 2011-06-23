/*
 * Copyright © 2011  Felix Höfling
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

#define BOOST_TEST_MODULE trajectory_h5md
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/casts.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <vector>
#include <stdint.h>

#include <halmd/io/trajectory/readers/h5md.hpp>
#include <halmd/io/trajectory/writers/h5md.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>
#include <test/tools/init.hpp>

using namespace boost;
using namespace halmd;
using namespace halmd::io; // avoid ambiguity of io:: between halmd::io and boost::io
using namespace std;

template <int dimension>
void h5md(std::vector<unsigned int> const& ntypes)
{
    typedef observables::host::samples::phase_space<dimension, float> float_sample_type;
    typedef observables::host::samples::phase_space<dimension, double> double_sample_type;

    typedef fixed_vector<float, dimension> float_vector_type;
    typedef fixed_vector<double, dimension> double_vector_type;

    string const filename("test_io_trajectory_h5md.trj");

    BOOST_MESSAGE("Testing " << ntypes.size() << " particle types");

    // construct phase space sample and fill with positions and velocities
    shared_ptr<double_sample_type> double_sample = make_shared<double_sample_type>(ntypes);
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        typename double_sample_type::sample_vector& r_sample = *double_sample->r[type];
        typename double_sample_type::sample_vector& v_sample = *double_sample->v[type];
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
    double_sample->step = 0;

    // copy sample to single precision
    shared_ptr<float_sample_type> float_sample = make_shared<float_sample_type>(ntypes);
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        transform(
            double_sample->r[type]->begin()
          , double_sample->r[type]->end()
          , float_sample->r[type]->begin()
          , lambda::ll_static_cast<float_vector_type>(lambda::_1)
        );
        transform(
            double_sample->v[type]->begin()
          , double_sample->v[type]->end()
          , float_sample->v[type]->begin()
          , lambda::ll_static_cast<float_vector_type>(lambda::_1)
        );
    }
    float_sample->step = double_sample->step;

    // write single-precision sample to file
    shared_ptr<mdsim::clock> clock = make_shared<mdsim::clock>(1); // bogus time-step
    shared_ptr<trajectory::writer<dimension> > writer =
        make_shared<trajectory::writers::h5md<dimension, float> >(float_sample, clock, filename);

    writer->append(0);
    writer->flush();

    // overwrite with double-precision data,
    // resetting the shared_ptr first closes the HDF5 file
    writer.reset();
    writer =
        make_shared<trajectory::writers::h5md<dimension, double> >(double_sample, clock, filename);

    writer->append(clock->step());
    writer->flush();

    // simulate an integration step for the very first particle
    clock->advance();
    (*double_sample->r[0])[0] += (*double_sample->v[0])[0];
    (*double_sample->v[0])[0] = double_vector_type(sqrt(2));
    double_sample->step = clock->step();

    writer->append(clock->step());
    writer->flush();

    // test integrity of H5MD file
    bool is_h5md = trajectory::readers::h5md<dimension, double>::format(filename);
    BOOST_CHECK(is_h5md);

    // read phase space sample #1 from file in double precision
    // reading is done upon construction, so we use an unnamed, temporary reader object
    // allocate memory for reading back the phase space sample
    shared_ptr<double_sample_type> double_sample_ = make_shared<double_sample_type>(ntypes);
    trajectory::readers::h5md<dimension, double>(double_sample_, filename, 1);

    // check binary equality of written and read data
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        BOOST_CHECK(equal(
            double_sample_->r[type]->begin()
          , double_sample_->r[type]->end()
          , double_sample->r[type]->begin()
        ));
        BOOST_CHECK(equal(
            double_sample_->v[type]->begin()
          , double_sample_->v[type]->end()
          , double_sample->v[type]->begin()
        ));
    }


    // read phase space sample #0 from file in single precision
    shared_ptr<float_sample_type> float_sample_ = make_shared<float_sample_type>(ntypes);
    trajectory::readers::h5md<dimension, float>(float_sample_, filename, 0);

    // check binary equality of written and read data,
    // note that float_sample was not modified and thus corresponds to #0
    for (unsigned int type = 0; type < ntypes.size(); ++type) {
        BOOST_CHECK(equal(
            float_sample_->r[type]->begin()
          , float_sample_->r[type]->end()
          , float_sample->r[type]->begin()
        ));
        BOOST_CHECK(equal(
            float_sample_->v[type]->begin()
          , float_sample_->v[type]->end()
          , float_sample->v[type]->begin()
        ));
    }

    // close and remove file
    writer.reset();
    remove(filename.c_str());
}

HALMD_TEST_INIT( init_unit_test_suite )
{
    using namespace boost::assign;
    using namespace boost::unit_test;
    using namespace boost::unit_test::framework;

    std::vector<std::vector<unsigned int> > ntypes;
    ntypes += list_of(1), list_of(1)(10), list_of(1)(10)(100);

    test_suite* ts = BOOST_TEST_SUITE( "trajectory_h5md" );

    test_suite* ts1 = BOOST_TEST_SUITE( "2d" );
    ts1->add( BOOST_PARAM_TEST_CASE( &h5md<2>, ntypes.begin(), ntypes.end() ) );

    test_suite* ts2 = BOOST_TEST_SUITE( "3d" );
    ts2->add( BOOST_PARAM_TEST_CASE( &h5md<3>, ntypes.begin(), ntypes.end() ) );

    ts->add( ts1 );
    ts->add( ts2 );

    master_test_suite().add( ts );
}
