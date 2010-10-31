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

#ifndef HALMD_UTILITY_PROGRAM_OPTIONS_TYPED_VALUE_HPP
#define HALMD_UTILITY_PROGRAM_OPTIONS_TYPED_VALUE_HPP

#include <boost/program_options.hpp>

namespace halmd
{
namespace po
{

class value_semantic
{
public:
    virtual ~value_semantic() {}

    virtual std::vector<std::string> const& conflicting() const = 0;

    virtual std::vector<std::string> const& dependent() const = 0;
};

template <typename typed_value, typename T, typename charT = char>
class typed_value_base
  : public boost::program_options::typed_value<T, charT>
  , public value_semantic
{
public:
    typed_value_base(T* store_to)
      : boost::program_options::typed_value<T, charT>(store_to) {}

    typed_value* default_value(const T& v)
    {
        boost::program_options::typed_value<T, charT>::default_value(v);
        return static_cast<typed_value*>(this);
    }

    typed_value* default_value(const T& v, const std::string& textual)
    {
        boost::program_options::typed_value<T, charT>::default_value(v, textual);
        return static_cast<typed_value*>(this);
    }

    typed_value* implicit_value(const T &v)
    {
        boost::program_options::typed_value<T, charT>::implicit_value(v);
        return static_cast<typed_value*>(this);
    }

    typed_value* implicit_value(const T &v, const std::string& textual)
    {
        boost::program_options::typed_value<T, charT>::implicit_value(v, textual);
        return static_cast<typed_value*>(this);
    }

    typed_value* notifier(boost::function1<void, const T&> f)
    {
        boost::program_options::typed_value<T, charT>::notifier(f);
        return static_cast<typed_value*>(this);
    }

    typed_value* composing()
    {
        boost::program_options::typed_value<T, charT>::composing();
        return static_cast<typed_value*>(this);
    }

    typed_value* multitoken()
    {
        boost::program_options::typed_value<T, charT>::multitoken();
        return static_cast<typed_value*>(this);
    }

    typed_value* zero_tokens()
    {
        boost::program_options::typed_value<T, charT>::zero_tokens();
        return static_cast<typed_value*>(this);
    }

    typed_value* required()
    {
        boost::program_options::typed_value<T, charT>::required();
        return static_cast<typed_value*>(this);
    }

    typed_value* conflicts(std::string const& name)
    {
        conflicting_.push_back(name);
        return static_cast<typed_value*>(this);
    }

    typed_value* depends(std::string const& name)
    {
        dependent_.push_back(name);
        return static_cast<typed_value*>(this);
    }

public:
    virtual std::vector<std::string> const& conflicting() const
    {
        return conflicting_;
    }

    virtual std::vector<std::string> const& dependent() const
    {
        return dependent_;
    }

private:
    std::vector<std::string> conflicting_;
    std::vector<std::string> dependent_;
};

template <typename T, typename charT = char>
class typed_value
  : public typed_value_base<typed_value<T, charT>, T, charT>
{
public:
    typed_value(T* store_to)
      : typed_value_base<typed_value<T, charT>, T, charT>(store_to) {}
};

template <typename T>
typed_value<T>* value(T* v)
{
    typed_value<T>* r = new typed_value<T>(v);
    return r;
}

template <typename T>
typed_value<T>* value()
{
    return value<T>(0);
}

inline typed_value<bool>* bool_switch(bool* v)
{
    typed_value<bool>* r = new typed_value<bool>(v);
    r->default_value(0);
    r->zero_tokens();
    return r;
}

inline typed_value<bool>* bool_switch()
{
    return bool_switch(0);
}

} // namespace po

} // namespace halmd

#endif /* ! HALMD_UTILITY_PROGRAM_OPTIONS_TYPED_VALUE_HPP */
