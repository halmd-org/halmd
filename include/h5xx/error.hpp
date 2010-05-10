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

// ensure compatibility with HDF5 1.6.x and 1.8.x
#define H5_USE_16_API_DEFAULT
#include <hdf5.h>

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
 * wrap HDF5 C API calls with this macro for error handling
 */
#define H5XX_CALL(expr) \
    do { \
        h5xx::error::_register register_; \
        if ((expr) < 0) { \
            throw error(); \
        } \
    } while(0)

/**
 * HDF5 exception
 */
class error
  : virtual public std::exception
{
public:
    explicit error(std::string err) : err_(err) {}
    virtual ~error() throw() {}

    /**
     * retrieve error description from HDF5 error stack
     */
    error()
    {
        H5Ewalk(H5E_WALK_DOWNWARD, reinterpret_cast<H5E_walk_t>(walk), &err_);
    }

    /**
     * returns error description
     */
    char const* what() const throw()
    {
        return err_.c_str();
    }

    /**
     * scoped error handler
     */
    class _register
    {
    public:
        /**
         * set error handler
         */
        _register()
        {
            // We do not use H5Eget_auto to save and later restore the
            // current error handler, as for HDF5 1.8.x the type of
            // the returned function pointer may vary depending on the
            // compile time option --with-default-api-version, even if
            // we set a compatibility macro to explicitly demand the
            // HDF 1.6 API.

            H5Eset_auto(NULL, NULL);
        }

        /**
         * unset error handler
         */
        ~_register()
        {
            H5Eset_auto(reinterpret_cast<H5E_auto_t>(H5Eprint), NULL);
        }
    };

private:
    /**
     * retrieve function name and error description of top stack entry
     */
    static herr_t walk(int n, H5E_error_t* err_desc, std::string* err)
    {
        if (n == 0) {
            *err = err_desc->func_name;
            if (err_desc->desc) {
                *err += std::string(": ") + err_desc->desc;
            }
        }
        // value of HDF5-internal macro SUCCESS
        return 0;
    }

    /** error description */
    std::string err_;
};

} // namespace h5xx

#endif /* ! H5XX_ERROR_HPP */
