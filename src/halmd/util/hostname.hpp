/* hostname.hpp
 *
 * Copyright © 2009  Peter Colberg
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

#ifndef HALMD_UTIL_HOSTNAME_HPP
#define HALMD_UTIL_HOSTNAME_HPP

#include <cstring>
#include <exception>
#include <netdb.h>
#include <string>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/utsname.h>

namespace halmd
{

/**
 * get hostname of this machine
 */
std::string get_hostname()
{
    utsname name;

    if (uname(&name) != 0) {
        throw std::runtime_error(strerror(errno));
    }
    return name.nodename;
}

/**
 * get fully qualified hostname of given node
 */
std::string get_fqdn_hostname(std::string const& node)
{
    struct addrinfo hints, *result;
    int gai_result;

    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_UNSPEC; // either IPV4 or IPV6
    hints.ai_flags = AI_CANONNAME;
    hints.ai_socktype = SOCK_STREAM;

    if ((gai_result = getaddrinfo(node.c_str(), "domain", &hints, &result)) != 0) {
        throw std::runtime_error(gai_strerror(gai_result));
    }

    // return first resolved hostname
    std::string hostname(result->ai_canonname);
    freeaddrinfo(result);
    return hostname;
}

} // namespace halmd

#endif /* ! HALMD_UTIL_HOSTNAME_HPP */
