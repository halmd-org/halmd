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

#include <boost/test/test_tools.hpp> // BOOST_TEST_MESSAGE

#include <halmd/io/logger.hpp>
#include <halmd/version.h>

using namespace halmd;

// We may not include <test/tools/ctest.hpp> in this file, since this would
// add the global fixtures ctest_full_output, and in turn cause duplicated
// output in a test linking against this library.
// Therefore repeat the class definition of ctest_full_output here.

struct ctest_full_output
{
    ctest_full_output();
};

ctest_full_output::ctest_full_output()
{
    BOOST_TEST_MESSAGE( "Avoid ctest truncation of output: CTEST_FULL_OUTPUT" );
    BOOST_TEST_MESSAGE( PROJECT_NAME " " PROGRAM_VERSION );
}

struct ctest_logging
{
    ctest_logging();
};

ctest_logging::ctest_logging()
{
#ifdef NDEBUG
    logging::get().open_console(logging::info);
#else
    logging::get().open_console(logging::info);     // logging::debug creates huge logging files (~1 GB)
#endif
}
