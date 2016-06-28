/*
 * Copyright © 2011-2012  Peter Colberg
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

#ifndef TEST_TOOLS_CTEST_HPP
#define TEST_TOOLS_CTEST_HPP

#include <boost/test/unit_test_suite.hpp> // BOOST_GLOBAL_FIXTURE
#include <boost/version.hpp>
// avoid warnings about superfluous semicola
// after BOOST_GLOBAL_FIXTURE for boost < 1.59
#if BOOST_VERSION < 105900 && defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wpedantic"
#endif

/**
 * Print CTEST_FULL_OUTPUT to avoid ctest truncation of output.
 *
 * Requires --log_level=message or --log_level=test_suite.
 */
struct ctest_full_output
{
    ctest_full_output();
};

/**
 * To avoid static initialization order fiasco, add the CTEST_FULL_OUTPUT
 * printer as a global fixture instead of an __attribute__((constructor))
 * function. The global fixture is instantiated before the test run.
 */
BOOST_GLOBAL_FIXTURE( ctest_full_output );

#ifndef HALMD_TEST_NO_LOGGING
/**
 * Enable logging to console.
 *
 * If built without debugging (NDEBUG), log with severity info.
 *
 * Otherwise, log with severity debug.
 */
struct ctest_logging
{
    ctest_logging();
};

BOOST_GLOBAL_FIXTURE( ctest_logging );
#endif

#endif /* ! TEST_TOOLS_CTEST_HPP */
