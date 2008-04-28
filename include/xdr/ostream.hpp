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

#ifndef XDR_OSTREAM_HPP
#define XDR_OSTREAM_HPP

#include <rpc/xdr.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <string>
#include "error.hpp"


namespace xdr {

class ostream
{
public:
    /**
     * initialize the XDR stream object
     */
    ostream(FILE* file)
    {
	xdrstdio_create(&xdrs_, file, XDR_ENCODE);
    }

    /**
     * destroy the XDR stream object
     */
    ~ostream()
    {
	xdr_destroy(&xdrs_);
    }

    /**
     * translate from C character to XDR
     */
    ostream& operator<<(char const& data)
    {
	XDR_SAFE_CALL(xdr_char, &xdrs_, const_cast<char*>(&data));
	return *this;
    }

    /**
     * translate from C unsigned character to XDR
     */
    ostream& operator<<(unsigned char const& data)
    {
	XDR_SAFE_CALL(xdr_u_char, &xdrs_, const_cast<unsigned char*>(&data));
	return *this;
    }

    /**
     * translate from C short integer to XDR
     */
    ostream& operator<<(short const& data)
    {
	XDR_SAFE_CALL(xdr_short, &xdrs_, const_cast<short*>(&data));
	return *this;
    }

    /**
     * translate from C unsigned short integer to XDR
     */
    ostream& operator<<(unsigned short const& data)
    {
	XDR_SAFE_CALL(xdr_u_short, &xdrs_, const_cast<unsigned short*>(&data));
	return *this;
    }

    /**
     * translate from C integer to XDR
     */
    ostream& operator<<(int const& data)
    {
	XDR_SAFE_CALL(xdr_int, &xdrs_, const_cast<int*>(&data));
	return *this;
    }

    /**
     * translate from C unsigned integer to XDR
     */
    ostream& operator<<(unsigned int const& data)
    {
	XDR_SAFE_CALL(xdr_u_int, &xdrs_, const_cast<unsigned int*>(&data));
	return *this;
    }

    /**
     * translate from C long integer to XDR
     */
    ostream& operator<<(long const& data)
    {
	XDR_SAFE_CALL(xdr_long, &xdrs_, const_cast<long*>(&data));
	return *this;
    }

    /**
     * translate from C unsigned long integer to XDR
     */
    ostream& operator<<(unsigned long const& data)
    {
	XDR_SAFE_CALL(xdr_u_long, &xdrs_, const_cast<unsigned long*>(&data));
	return *this;
    }

    /**
     * translate from C float to XDR
     */
    ostream& operator<<(float const& data)
    {
	XDR_SAFE_CALL(xdr_float, &xdrs_, const_cast<float*>(&data));
	return *this;
    }

    /**
     * translate from C double to XDR
     */
    ostream& operator<<(double const& data)
    {
	XDR_SAFE_CALL(xdr_double, &xdrs_, const_cast<double*>(&data));
	return *this;
    }

    /**
     * translate from C string to XDR
     */
    ostream& operator<<(char const* data)
    {
	char *string = strdup(data);
	XDR_SAFE_CALL(xdr_wrapstring, &xdrs_, &string);
	free(string);
	return *this;
    }

    /**
     * translate from C++ string to XDR
     */
    ostream& operator<<(std::string const& data)
    {
	char *string = strdup(data.c_str());
	XDR_SAFE_CALL(xdr_wrapstring, &xdrs_, &string);
	free(string);
	return *this;
    }

private:
    XDR xdrs_;
};


/**
 * translate from STL sequence container to XDR array
 */
template <typename T>
ostream& operator<<(ostream& xdrs, T const& data)
{
    xdrs << u_int(data.size());
    for (typename T::const_iterator it = data.begin(); it != data.end(); ++it) {
	xdrs << *it;
    }
    return xdrs;
}

} // namespace xdr

#endif /* ! XDR_OSTREAM_HPP */
