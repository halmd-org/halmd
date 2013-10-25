/*
 * Copyright © 2013 Felix Höfling
 * Copyright © 2012 Peter Colberg
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

#define BOOST_TEST_MODULE continuation
#include <boost/test/unit_test.hpp>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

#include <boost/lexical_cast.hpp>
#include <h5xx/h5xx.hpp>

#include <vector>

void compare_datasets(H5::DataSet const& ds1, H5::DataSet const& ds2, hssize_t offset1, hssize_t offset2)
{
    std::vector<hsize_t> dims(ds1.getSpace().getSimpleExtentNdims());
    ds1.getSpace().getSimpleExtentDims(&*dims.begin());

    BOOST_REQUIRE( dims.size() == 2 || dims.size() == 3 );
    if (dims.size() == 2) {
        // compare scalar mass, species, …
        // assumes lossless conversion from integer to floating-point
        std::vector<double> array1, array2;
        h5xx::read_chunked_dataset(ds1, array1, offset1);
        h5xx::read_chunked_dataset(ds2, array2, offset2);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            array1.begin()
          , array1.end()
          , array2.begin()
          , array2.end()
        );
    }
    else if (dims.size() == 3) {
        // compare 2 or 3-dimensional positions, velocities, …
        BOOST_REQUIRE( dims[2] == 2 || dims[2] == 3 );
        if (dims[2] == 3) {
            std::vector<halmd::fixed_vector<double, 3> > array1, array2;
            h5xx::read_chunked_dataset(ds1, array1, offset1);
            h5xx::read_chunked_dataset(ds2, array2, offset2);
            BOOST_CHECK_EQUAL_COLLECTIONS(
                array1.begin()
              , array1.end()
              , array2.begin()
              , array2.end()
            );
        }
        else if (dims[2] == 2) {
            std::vector<halmd::fixed_vector<double, 2> > array1, array2;
            h5xx::read_chunked_dataset(ds1, array1, offset1);
            h5xx::read_chunked_dataset(ds2, array2, offset2);
            BOOST_CHECK_EQUAL_COLLECTIONS(
                array1.begin()
              , array1.end()
              , array2.begin()
              , array2.end()
            );
        }
    }
}

/**
 * This test case compares two H5MD particles groups given a pair of
 * H5MD files and dataset offsets. The comparison is binary exact, which
 * assumes reproducible floating-point operations used in the periodic
 * extension of particle positions for writing.
 */
BOOST_AUTO_TEST_CASE( compare_particles )
{
    using namespace boost::unit_test::framework;
    std::vector<std::string> arg(
        master_test_suite().argv
      , master_test_suite().argv + master_test_suite().argc
    );

    BOOST_REQUIRE_EQUAL( arg.size(), 5u );
    H5::H5File file1(arg[1], H5F_ACC_RDONLY);
    H5::H5File file2(arg[2], H5F_ACC_RDONLY);
    hssize_t offset1 = boost::lexical_cast<hssize_t>(arg[3]);
    hssize_t offset2 = boost::lexical_cast<hssize_t>(arg[4]);
    H5::Group trajectory1 = file1.openGroup("particles");
    H5::Group trajectory2 = file2.openGroup("particles");

    unsigned int ngroup = trajectory1.getNumObjs();
    for (unsigned int i = 0; i < ngroup; ++i) {
        H5::Group group1 = trajectory1.openGroup(trajectory1.getObjnameByIdx(i));
        BOOST_TEST_MESSAGE( "group " << h5xx::path(group1) );
        H5::Group group2 = trajectory2.openGroup(trajectory1.getObjnameByIdx(i));

        unsigned int narray = group1.getNumObjs();
        for (unsigned int j = 0; j < narray; ++j) {
            H5::Group subgroup1 = group1.openGroup(group1.getObjnameByIdx(j));
            H5::Group subgroup2 = group2.openGroup(group2.getObjnameByIdx(j));
            if (group1.getObjnameByIdx(j) == "box") {
                H5::DataSet value1 = subgroup1.openGroup("edges").openDataSet("value");
                BOOST_TEST_MESSAGE( "dataset " << h5xx::path(value1) );
                H5::DataSet value2 = subgroup2.openGroup("edges").openDataSet("value");
                compare_datasets(value1, value2, offset1, offset2);

                value1 = subgroup1.openGroup("offset").openDataSet("value");
                BOOST_TEST_MESSAGE( "dataset " << h5xx::path(value1) );
                value2 = subgroup2.openGroup("offset").openDataSet("value");
                compare_datasets(value1, value2, offset1, offset2);
            }
            else {
                H5::DataSet value1 = subgroup1.openDataSet("value");
                BOOST_TEST_MESSAGE( "dataset " << h5xx::path(value1) );
                H5::DataSet value2 = subgroup2.openDataSet("value");
                compare_datasets(value1, value2, offset1, offset2);
            }
        }
    }
}
