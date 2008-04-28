/* xdr - library routines for external data representation
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef XDR_ERROR_HPP
#define XDR_ERROR_HPP

#include <exception>


namespace xdr
{

#define XDR_SAFE_CALL(func, ...) if (!func(__VA_ARGS__)) throw xdr::exception(#func " filter primitive failed")


/**
 * XDR procedure exception
 */
class exception: public std::exception
{
public:
    exception(char const* str): str_(str) {}

    char const* what() const throw()
    {
	return str_;
    }

private:
    char const* str_;
};

} // namespace xdr

#endif /* ! XDR_ERROR_HPP */
