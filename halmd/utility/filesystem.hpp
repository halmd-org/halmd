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
#include <cassert>

namespace halmd {

/**
 * returns absolute filesystem path
 *
 * @param path relative or absolute path
 */
inline std::string absolute_path(std::string const& path)
{
    return boost::filesystem::absolute(path).string();
}

/**
 * returns true if first path contains second path
 */
inline bool contains_path(boost::filesystem::path const& p1, boost::filesystem::path const& p2)
{
    assert( p1.is_absolute() );
    assert( p2.is_absolute() );
    boost::filesystem::path p(p2);
    while (!p.empty()) {
        if (boost::filesystem::equivalent(p1, p)) {
            return true;
        }
        p = p.parent_path();
    }
    return false;
}

} // namespace halmd

#endif /* HALMD_UTILITY_FILESYSTEM_HPP */
