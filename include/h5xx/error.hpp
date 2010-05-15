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

#ifndef H5XX_ERROR_HPP
#define H5XX_ERROR_HPP

#include <exception>
#include <string>

#include <h5xx/compat.hpp>

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
            h5xx::throw_exception(); \
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
 * throw HDF5 library exception
 */
inline void throw_exception()
{
    H5E_error_t const* _err;
#ifndef H5XX_USE_16_API
    H5Ewalk(H5E_DEFAULT, H5E_WALK_DOWNWARD, reinterpret_cast<H5E_walk_t>(_walk_stack), &_err);
#else
    H5Ewalk(H5E_WALK_DOWNWARD, reinterpret_cast<H5E_walk_t>(_walk_stack), &_err);
#endif
    throw h5xx::error(_err);
}

/**
 * silence HDF5 error stack printing
 *
 * This is an exception-safe adaptation of the HDF5 library macros
 * H5E_BEGIN_TRY and H5E_END_TRY for all supported HDF5 versions.
 */
struct _error_handler
{
#ifndef HDF5_DEBUG
# ifndef H5XX_USE_16_API
#  ifndef H5_NO_DEPRECATED_SYMBOLS
    unsigned saved_is_v2;
    union {
        H5E_auto1_t efunc1;
        H5E_auto2_t efunc2;
    } saved;
    void *saved_edata;

    _error_handler()
    {
        H5Eauto_is_v2(H5E_DEFAULT, &saved_is_v2);
        if (saved_is_v2) {
            H5Eget_auto2(H5E_DEFAULT, &saved.efunc2, &saved_edata);
            H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
        }
        else {
            H5Eget_auto1(&saved.efunc1, &saved_edata);
            H5Eset_auto1(NULL, NULL);
        }
    }

    ~_error_handler()
    {
        if(saved_is_v2) {
            H5Eset_auto2(H5E_DEFAULT, saved.efunc2, saved_edata);
        }
        else {
            H5Eset_auto1(saved.efunc1, saved_edata);
        }
    }
#  else /* H5_NO_DEPRECATED_SYMBOLS */
    H5E_auto_t saved_efunc;
    void* saved_edata;

    _error_handler()
    {
        H5Eget_auto(H5E_DEFAULT, &saved_efunc, &saved_edata);
        H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    }

    ~_error_handler()
    {
        H5Eset_auto(H5E_DEFAULT, saved_efunc, saved_edata);
    }
#  endif /* H5_NO_DEPRECATED_SYMBOLS */
# else /* H5XX_USE_16_API */
    H5E_auto_t saved_efunc;
    void* saved_edata;

    _error_handler()
    {
        H5Eget_auto(&saved_efunc, &saved_edata);
        H5Eset_auto(NULL, NULL);
    }

    ~_error_handler()
    {
        H5Eset_auto(saved_efunc, saved_edata);
    }
# endif /* H5XX_USE_16_API */
#endif /* ! HDF5_DEBUG */
};

} // namespace h5xx

#endif /* ! H5XX_ERROR_HPP */
