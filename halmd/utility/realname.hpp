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

#ifndef HALMD_UTILITY_REALNAME_HPP
#define HALMD_UTILITY_REALNAME_HPP

#include <boost/system/error_code.hpp>
#include <boost/system/system_error.hpp>
#include <pwd.h>
#include <string>
#include <stdexcept> // std::runtime_error
#include <sys/types.h>
#include <unistd.h>
#include <vector>

namespace halmd {

/**
 * returns real name of the user of the calling process
 *
 * This function retrieves the password file entry for the real user id of
 * the calling process (or optional argument uid). If the GECOS field is
 * empty, the login name of the user is returned. If the GECOS field is
 * non-empty, either the part up to the first comma, or otherwise the
 * entire field, is returned.
 *
 * For details, see passwd(5) and getpwnam(3).
 *
 * See this bug report for rationale on splitting upon first comma:
 * https://bugzilla.samba.org/show_bug.cgi?id=5198
 *
 * Notice that in contrast to Samba, we always split upon first comma, i.e.
 * not only if at least three commas are detected, as the department of the
 * author of this function uses the GECOS format "Full Name,group".
 */
inline std::string realname(uid_t uid = getuid())
{
    passwd pwd;
    std::size_t buflen = sysconf(_SC_GETPW_R_SIZE_MAX); // e.g. 1024
    std::vector<char> buf(buflen);
    passwd* result;
    boost::system::error_code ec(
        getpwuid_r(uid, &pwd, &*buf.begin(), buf.size(), &result)
      , boost::system::get_posix_category()
    );
    if (ec != boost::system::errc::success) {
        throw boost::system::system_error(ec);
    }
    if (!result) {
        throw std::runtime_error("no password record found");
    }
    std::string gecos = result->pw_gecos;
    if (gecos.empty()) {
        return result->pw_name;
    }
    std::size_t pos = gecos.find_first_of(',');
    if (pos == 0) {
        return result->pw_name;
    }
    if (pos != std::string::npos) {
        return gecos.substr(0, pos);
    }
    return gecos;
}

} // namespace halmd

#endif /* HALMD_UTILITY_REALNAME_HPP */
