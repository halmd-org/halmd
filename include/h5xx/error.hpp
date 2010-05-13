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

#ifndef H5XX_ERROR_HPP
#define H5XX_ERROR_HPP

#include <exception>
#include <string>

#include <hdf5.h>
#if H5_VERS_MAJOR <= 1 && H5_VERS_MINOR < 8
# error "h5xx wrapper requires HDF5 >= 1.8"
#endif

namespace h5xx
{

// The C++ wrapper shipped with HDF5 implements exceptions in an
// unusual way. H5::Exception does not derive from std::exception.
// Instead of holding an error description in the exception
// object, the contents of an error stack are dumped to stderr,
// regardless of whether the exception is caught or not.
//
// Disabling the error stack printing and catching H5::Exception
// to print the error message manually is not a solution either,
// as this hides the throw location of the exception, which impedes
// debugging with gdb.
//
// The solution is to implement a custom error handler that
// throws a standards-compliant exception, h5xx:error, which
// returns an error message via the what() function.
//

// For a detailed analysis of the shortcomings of the HDF5 C++ library,
// see "HDF5 C++ API and programming model" by Manoj Rajagopalan.
//
// http://mail.hdfgroup.org/pipermail/hdf-forum_hdfgroup.org/2010-April/003048.html
//

/**
 * HDF5 exception
 */
class error
  : virtual public std::exception
{
public:
    /**
     * set HDF5 library error description
     */
    error(H5E_error_t const* err)
      : desc_(err->func_name + std::string(": ") + err->desc) {}

    /**
     * set custom error description
     */
    error(std::string const& desc)
      : desc_(desc) {}

    virtual ~error() throw() {}

    /**
     * returns error description
     */
    char const* what() const throw()
    {
        return desc_.c_str();
    }

private:
    std::string desc_;
};

/**
 * wrap HDF5 C API calls with this macro for error interception
 */
#define H5XX_CALL(expr) \
    do { \
        h5xx::_error_handler _no_print; \
        if ((expr) < 0) { \
            H5E_error_t const* _err; \
            H5Ewalk(H5E_DEFAULT, H5E_WALK_DOWNWARD, h5xx::_walk_stack, &_err); \
            throw h5xx::error(_err); \
        } \
    } while(0)

/**
 * retrieve error description on top of error stack
 */
inline herr_t _walk_stack(unsigned int n, H5E_error_t const* err, void* err_ptr)
{
    if (n == 0) {
        *reinterpret_cast<void const**>(err_ptr) = err;
    }
    // value of HDF5-internal macro SUCCESS
    return 0;
}

/**
 * silence HDF5 error stack printing
 */
struct _error_handler
{
#ifndef HDF5_DEBUG
    _error_handler()
    {
        // We do not use H5Eget_auto to save and later restore the
        // current error handler, as for HDF5 1.8.x the type of
        // the returned function pointer may vary depending on the
        // compile time option --with-default-api-version.

        H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    }

    ~_error_handler()
    {
        H5Eset_auto(H5E_DEFAULT, reinterpret_cast<H5E_auto_t>(H5Eprint), NULL);
    }
#endif /* ! HDF5_DEBUG */
};

} // namespace h5xx

#endif /* ! H5XX_ERROR_HPP */
