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

namespace halmd
{
namespace utility { namespace module
{

/**
 * Module exception
 */
class module_exception
  : public std::exception
{
public:
    module_exception(std::string const& what_) : what_(what_) {}
    virtual ~module_exception() throw () {}
    virtual const char* what() const throw()
    {
        return what_.c_str();
    }
private:
    std::string what_;
};

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_EXCEPTION_HPP */
