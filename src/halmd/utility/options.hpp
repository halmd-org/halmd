/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_UTILITY_OPTIONS_HPP
#define HALMD_UTILITY_OPTIONS_HPP

#include <boost/program_options.hpp>

#include <halmd/options.hpp> // FIXME transition
#include <halmd/utility/detail/options.hpp>

namespace halmd
{
namespace po
{

using boost::program_options::options_description;
using boost::program_options::value;

using halmd::options; // FIXME transition

} // namespace po

} // namespace halmd

#endif /* ! HALMD_UTILITY_OPTIONS_HPP */
