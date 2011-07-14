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
#include <boost/version.hpp>

namespace halmd {
namespace po {

class extended_value_semantic
{
public:
    virtual ~extended_value_semantic() {}

    virtual std::vector<std::string> const& conflicting() const = 0;

    virtual std::vector<std::string> const& dependent() const = 0;
};

template <typename extended_typed_value, typename T, typename charT = char>
class extended_typed_value_base
  : public boost::program_options::typed_value<T, charT>
  , public extended_value_semantic
{
private:
    typedef boost::program_options::typed_value<T, charT> _Base;

public:
    extended_typed_value_base(T* store_to)
      : _Base(store_to) {}

    extended_typed_value* default_value(const T& v)
    {
        _Base::default_value(v);
        return static_cast<extended_typed_value*>(this);
    }

    extended_typed_value* default_value(const T& v, const std::string& textual)
    {
        _Base::default_value(v, textual);
        return static_cast<extended_typed_value*>(this);
    }

    extended_typed_value* implicit_value(const T &v)
    {
        _Base::implicit_value(v);
        return static_cast<extended_typed_value*>(this);
    }

    extended_typed_value* implicit_value(const T &v, const std::string& textual)
    {
        _Base::implicit_value(v, textual);
        return static_cast<extended_typed_value*>(this);
    }

    extended_typed_value* notifier(boost::function1<void, const T&> f)
    {
        _Base::notifier(f);
        return static_cast<extended_typed_value*>(this);
    }

    extended_typed_value* composing()
    {
        _Base::composing();
        return static_cast<extended_typed_value*>(this);
    }

    extended_typed_value* multitoken()
    {
        _Base::multitoken();
        return static_cast<extended_typed_value*>(this);
    }

    extended_typed_value* zero_tokens()
    {
        _Base::zero_tokens();
        return static_cast<extended_typed_value*>(this);
    }

#if BOOST_VERSION >= 104200
    extended_typed_value* required()
    {
        _Base::required();
        return static_cast<extended_typed_value*>(this);
    }
#endif /* BOOST_VERSION >= 104200 */

    extended_typed_value* conflicts(std::string const& name)
    {
        conflicting_.push_back(name);
        return static_cast<extended_typed_value*>(this);
    }

    extended_typed_value* depends(std::string const& name)
    {
        dependent_.push_back(name);
        return static_cast<extended_typed_value*>(this);
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
class extended_typed_value
  : public extended_typed_value_base<extended_typed_value<T, charT>, T, charT>
{
public:
    extended_typed_value(T* store_to)
      : extended_typed_value_base<extended_typed_value<T, charT>, T, charT>(store_to) {}
};

template <typename T>
extended_typed_value<T>* value(T* v)
{
    extended_typed_value<T>* r = new extended_typed_value<T>(v);
    return r;
}

template <typename T>
extended_typed_value<T>* value()
{
    return value<T>(0);
}

inline extended_typed_value<bool>* bool_switch(bool* v)
{
    extended_typed_value<bool>* r = new extended_typed_value<bool>(v);
    r->default_value(0);
    r->zero_tokens();
    return r;
}

inline extended_typed_value<bool>* bool_switch()
{
    return bool_switch(0);
}

} // namespace po
} // namespace halmd

#endif /* ! HALMD_UTILITY_PROGRAM_OPTIONS_TYPED_VALUE_HPP */
