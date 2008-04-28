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

#ifndef XDR_ISTREAM_HPP
#define XDR_ISTREAM_HPP

#include <rpc/xdr.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <string>
#include "error.hpp"


namespace xdr {

class istream
{
public:
    /**
     * initialize the XDR stream object
     */
    istream(FILE* file)
    {
	xdrstdio_create(&xdrs_, file, XDR_DECODE);
    }

    /**
     * destroy the XDR stream object
     */
    ~istream()
    {
	xdr_destroy(&xdrs_);
    }

    /**
     * translate to C character from XDR
     */
    istream& operator>>(char& data)
    {
	XDR_SAFE_CALL(xdr_char, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C unsigned character from XDR
     */
    istream& operator>>(unsigned char& data)
    {
	XDR_SAFE_CALL(xdr_u_char, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C short integer from XDR
     */
    istream& operator>>(short& data)
    {
	XDR_SAFE_CALL(xdr_short, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C unsigned short integer from XDR
     */
    istream& operator>>(unsigned short& data)
    {
	XDR_SAFE_CALL(xdr_u_short, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C integer from XDR
     */
    istream& operator>>(int& data)
    {
	XDR_SAFE_CALL(xdr_int, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C unsigned integer from XDR
     */
    istream& operator>>(unsigned int& data)
    {
	XDR_SAFE_CALL(xdr_u_int, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C long integer from XDR
     */
    istream& operator>>(long& data)
    {
	XDR_SAFE_CALL(xdr_long, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C unsigned long integer from XDR
     */
    istream& operator>>(unsigned long& data)
    {
	XDR_SAFE_CALL(xdr_u_long, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C float from XDR
     */
    istream& operator>>(float& data)
    {
	XDR_SAFE_CALL(xdr_float, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C double from XDR
     */
    istream& operator>>(double& data)
    {
	XDR_SAFE_CALL(xdr_double, &xdrs_, &data);
	return *this;
    }

    /**
     * translate to C++ string from XDR
     */
    istream& operator>>(std::string& data)
    {
	char *string = NULL;
	XDR_SAFE_CALL(xdr_wrapstring, &xdrs_, &string);
	data = string;
	free(string);
	return *this;
    }

private:
    XDR xdrs_;
};


/**
 * translate to STL sequence container from XDR array
 */
template <typename T>
istream& operator>>(istream& xdrs, T& data)
{
    u_int size;
    xdrs >> size;
    data.resize(size);
    for (typename T::iterator it = data.begin(); it != data.end(); ++it) {
	xdrs >> *it;
    }
    return xdrs;
}

} // namespace xdr

#endif /* ! XDR_ISTREAM_HPP */
