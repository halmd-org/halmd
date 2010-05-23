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

#include <boost/operators.hpp>
#include <deque>
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

#define H5XX_ERROR_SCOPE_BEGIN {
#define H5XX_ERROR_SCOPE_END   } // satisfy Vim syntax highlighting

/**
 * wrap HDF5 C API calls with this macro for error interception
 */
#define H5XX_CHECK(expr) \
    H5E_BEGIN_TRY { \
        try { \
            if ((expr) < 0) { \
                throw h5xx::error(); \
            } \
        } \
        catch (...) { \
            H5XX_ERROR_SCOPE_BEGIN H5E_END_TRY; \
            throw; \
        } \
    } H5E_END_TRY

/**
 * HDF5 major and minor error number
 */
#ifdef H5XX_USE_16_API
typedef std::pair<H5E_major_t, H5E_minor_t> error_type;
#else
typedef std::pair<hid_t, hid_t> error_type;
#endif

/**
 * HDF5 error description
 */
struct error_description
  : boost::equality_comparable<
        error_description
      , error_type
    >
{
    /** major and minor error number */
    error_type const type;
    /** error summary */
    std::string const desc;

    /**
     * set HDF5 library error description
     */
    explicit error_description(H5E_error_t const& e)
      : type(e.maj_num, e.min_num)
      , desc(e.func_name + (e.desc ? std::string("(): ") + e.desc : std::string("()")))
    {}

    /**
     * check error type for equality
     */
    bool operator==(error_type const& type) const
    {
        return this->type == type;
    }
};

/**
 * HDF5 exception
 */
class error
  : virtual public std::exception
{
public:
    typedef std::deque<error_description> stack_type;

    /** error stack */
    stack_type const stack;

    /**
     * set HDF5 library error description
     */
    error()
    {
#ifdef H5XX_USE_16_API
        H5Ewalk(
            H5E_WALK_DOWNWARD
          , reinterpret_cast<H5E_walk_t>(walk_cb)
          , const_cast<stack_type*>(&stack)
        );
#else
        H5Ewalk(
            H5E_DEFAULT
          , H5E_WALK_DOWNWARD
          , reinterpret_cast<H5E_walk_t>(walk_cb)
          , const_cast<stack_type*>(&stack)
        );
#endif
    }

    virtual ~error() throw() {}

    /**
     * returns error description
     */
    char const* what() const throw()
    {
        if (stack.empty()) {
            return "HDF5 error stack is empty";
        }
        return stack.front().desc.c_str();
    }

    /**
     * returns number of equal error types in stack
     */
    unsigned int count(error_type const& type) const
    {
        return std::count(stack.begin(), stack.end(), type);
    }

private:
    /**
     * retrieve error description in error stack
     */
    static herr_t walk_cb(unsigned int n, H5E_error_t const* err, stack_type* stack)
    {
        stack->push_back(error_description(*err));
        return 0; // indicate SUCCESS
    }
};

} // namespace h5xx

#endif /* ! H5XX_ERROR_HPP */
