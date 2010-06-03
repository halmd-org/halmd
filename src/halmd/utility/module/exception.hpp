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

#ifndef HALMD_UTILITY_MODULE_EXCEPTION_HPP
#define HALMD_UTILITY_MODULE_EXCEPTION_HPP

#include <exception>
#include <string>

#include <halmd/utility/module/demangle.hpp>

namespace halmd
{
namespace utility { namespace module
{

// This is the base class of all module exceptions.
class module_error
  : public virtual std::exception // virtual inheritance avoids ambiguity
{
public:
    module_error(std::string const& msg)
      : msg_(msg)
    {}

    virtual ~module_error() throw () {}

    virtual const char* what() const throw()
    {
        return msg_.c_str();
    }

private:
    std::string msg_;
};

// This exception is thrown in the resolve function of a module
// if it is unsuitable due to missing program option(s).
class unsuitable_module
  : public module_error
{
public:
    unsuitable_module(std::string const& msg)
      : module_error(msg)
    {}
};

// This exception is thrown in the resolve function of the
// factory class, if no module is available for a given type.
class unresolvable_dependency
  : public module_error
{
public:
    unresolvable_dependency(std::string const& msg)
      : module_error(msg)
    {}
};

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_EXCEPTION_HPP */
