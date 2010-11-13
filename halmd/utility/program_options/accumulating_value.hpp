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

#ifndef HALMD_UTILITY_PROGRAM_OPTIONS_ACCUMULATING_VALUE_HPP
#define HALMD_UTILITY_PROGRAM_OPTIONS_ACCUMULATING_VALUE_HPP

#include <boost/program_options.hpp>

#include <halmd/utility/program_options/typed_value.hpp>

namespace halmd
{
namespace po
{

/**
 * Accumulating program option value
 */
template <typename T, typename charT = char>
class accumulating_value
  : public extended_typed_value_base<accumulating_value<T, charT>, T, charT>
{
    //
    // Originally written by Bryan Green <bgreen@nas.nasa.gov>
    // http://article.gmane.org/gmane.comp.lib.boost.user/29084
    //

public:
    accumulating_value(T* store_to)
      : extended_typed_value_base<accumulating_value<T, charT>, T, charT>(store_to)
      , origin(0)
    {
        extended_typed_value_base<accumulating_value<T, charT>, T, charT>::zero_tokens();
    }

    accumulating_value* default_value(T const& v)
    {
        // setting a default value sets the origin to that value
        origin = v;
        extended_typed_value_base<accumulating_value<T, charT>, T, charT>::default_value(v);
        return this;
    }

    accumulating_value* default_value(T const& v, std::string const& textual)
    {
        // setting a default value sets the origin to that value
        origin = v;
        extended_typed_value_base<accumulating_value<T, charT>, T, charT>::default_value(v, textual);
        return this;
    }

public:
    void xparse(boost::any& value_store, std::vector<std::basic_string<charT> > const& new_tokens) const
    {
        // if this is the first occurrence of the option, initialize it to the origin
        if (value_store.empty()) {
            value_store = boost::any(origin);
        }
        ++boost::any_cast<T&>(value_store);
    }

private:
    T origin; //< the numeric origin from which to increment upward
};

template <typename T>
accumulating_value<T>* accum_value(T* v)
{
    accumulating_value<T>* r = new accumulating_value<T>(v);
    return r;
}

template <typename T>
accumulating_value<T>* accum_value()
{
    return accum_value<T>(0);
}

} // namespace po

} // namespace halmd

#endif /* ! HALMD_UTILITY_PROGRAM_OPTIONS_ACCUMULATING_VALUE_HPP */
