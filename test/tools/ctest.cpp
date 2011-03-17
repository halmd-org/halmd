/*
 * Copyright © 2011  Peter Colberg
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

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>

/**
 * Print CTEST_FULL_OUTPUT to avoid ctest truncation of output.
 *
 * Requires --log_level=message or --log_level=test_suite.
 */
struct ctest_full_output
{
    ctest_full_output()
    {
        BOOST_TEST_MESSAGE( "Avoid ctest truncation of output: CTEST_FULL_OUTPUT" );
    }
};

/**
 * To avoid static initialization order fiasco, add the CTEST_FULL_OUTPUT
 * printer as a global fixture instead of an __attribute__((constructor))
 * function. The global fixture is instantiated before the test run.
 */
BOOST_GLOBAL_FIXTURE( ctest_full_output );
