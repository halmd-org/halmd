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

#ifndef HALMD_UTILITY_FILESYSTEM_HPP
#define HALMD_UTILITY_FILESYSTEM_HPP

#include <boost/filesystem.hpp>

namespace halmd
{

/**
 * returns absolute filesystem path
 *
 * @param path relative or absolute path
 */
std::string absolute_path(std::string const& path)
{
#if defined(BOOST_FILESYSTEM_VERSION) && BOOST_FILESYSTEM_VERSION >= 3
    return boost::filesystem::absolute(path).string();
#else
    return boost::filesystem::complete(path).string();
#endif
}

} // namespace halmd

#endif /* HALMD_UTILITY_FILESYSTEM_HPP */
