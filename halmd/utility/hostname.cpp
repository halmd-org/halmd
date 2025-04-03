/*
 * Copyright © 2010  Peter Colberg
 * Copyright © 2014  Nicolas Höft
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

#include <boost/asio.hpp>

#include <halmd/utility/hostname.hpp>

namespace halmd {

/**
 * resolve canonical host name
 */
std::string host_name()
{
    boost::asio::io_context ios;
    boost::system::error_code ec;
    namespace bai = boost::asio::ip;
    bai::tcp::resolver resolver(ios);
    bai::tcp::resolver::query query(bai::host_name(), "", bai::tcp::resolver::query::canonical_name);
    // use non-throwing variant of resolve()
    bai::tcp::resolver::iterator addr_iter = resolver.resolve(query, ec);
    bai::tcp::resolver::iterator end;
    if (addr_iter != end && !ec) {
        return addr_iter->host_name();
    }
    return query.host_name(); // failed to resolve canonical host name
}

} // namespace halmd
