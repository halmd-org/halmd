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

#ifndef HALMD_UTILITY_PROGRAM_OPTIONS_PROGRAM_OPTIONS_HPP
#define HALMD_UTILITY_PROGRAM_OPTIONS_PROGRAM_OPTIONS_HPP

#include <boost/program_options.hpp>

#include <halmd/utility/program_options/accumulating_value.hpp>
#include <halmd/utility/program_options/array.hpp>
#include <halmd/utility/program_options/errors.hpp>
#include <halmd/utility/program_options/typed_value.hpp>
#include <halmd/utility/program_options/variables_map.hpp>

namespace halmd {
namespace po {

// import Boost Program Options into this namespace for convenience
using namespace boost::program_options;

} // namespace po
} // namespace halmd

#endif /* ! HALMD_UTILITY_PROGRAM_OPTIONS_PROGRAM_OPTIONS_HPP */
