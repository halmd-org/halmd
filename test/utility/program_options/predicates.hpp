/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_TEST_UTILITY_PROGRAM_OPTIONS_PREDICATES_HPP
#define HALMD_TEST_UTILITY_PROGRAM_OPTIONS_PREDICATES_HPP

#include <boost/test/test_tools.hpp>

#include <halmd/utility/program_options/program_options.hpp>

/**
 * Returns true test predicate if po::store throws an exception of given
 * type and has given option name, otherwise returns false test predicate.
 */
template <typename error_type>
boost::test_tools::predicate_result
store_throws_option(
    halmd::po::parsed_options const& parsed, halmd::po::variables_map& vm
  , std::string const& option_name)
{
    try {
        halmd::po::store(parsed, vm);
    }
    catch (error_type const& e) {
        if (e.get_option_name() != option_name) {
            boost::test_tools::predicate_result res(false);
            res.message()
                << "Exception has wrong option [" << e.get_option_name()
                << " != " << option_name << "]"
                ;
            return res;
        }
        return true;
    }
    boost::test_tools::predicate_result res(false);
    res.message() << "No exception thrown";
    return res;
}

#endif /* ! HALMD_TEST_UTILITY_PROGRAM_OPTIONS_PREDICATES_HPP */
