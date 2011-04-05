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

#include <halmd/io/logger.hpp>

/**
 * enable logging to console
 */
struct log_to_console
  : halmd::logger
{
    log_to_console()
    {
#ifdef NDEBUG
        halmd::logger::log_to_console(logger::warning);
#else
        halmd::logger::log_to_console(logger::debug);
#endif
    }
};

BOOST_GLOBAL_FIXTURE( log_to_console );
