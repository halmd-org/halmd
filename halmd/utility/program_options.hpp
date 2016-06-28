/*
 * Copyright © 2010-2012  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_UTILITY_PROGRAM_OPTIONS_HPP
#define HALMD_UTILITY_PROGRAM_OPTIONS_HPP

#include <boost/program_options/option.hpp>
#include <boost/program_options/value_semantic.hpp>

namespace halmd {

/**
 * This functor is an extra_style_parser for boost::program_options.
 *
 * If a positional argument is encountered, i.e. a string not beginning
 * with '-', a '--' terminator is prepended to the list of arguments.
 *
 * If the argument equals '-', the argument is substituted with '',
 * and a '--' terminator is prepended to the list of arguments.
 */
class inject_option_terminator
{
private:
    typedef std::vector<boost::program_options::option> result_type;

public:
    result_type operator()(std::vector<std::string>& args) const
    {
        std::vector<std::string>::iterator i(args.begin());
        if ((*i).substr(0, 1) != "-" || *i == "-") {
            if (*i == "-") {
                (*i).clear();
            }
            args.insert(i, "--");
        }
        return result_type();
    }
};

/**
 * Accumulating program option value
 */
template <typename T, typename charT = char>
class accumulating_value
  : public boost::program_options::typed_value<T, charT>
{
private:
    typedef boost::program_options::typed_value<T, charT> _Base;

    // Originally written by Bryan Green <bgreen@nas.nasa.gov>
    // http://article.gmane.org/gmane.comp.lib.boost.user/29084

public:
    accumulating_value(T* store_to)
      : _Base(store_to)
      , origin_(0)
    {
        _Base::zero_tokens();
    }

    accumulating_value* default_value(T const& v)
    {
        // setting a default value sets the origin to that value
        origin_ = v;
        _Base::default_value(v);
        return this;
    }

    accumulating_value* default_value(T const& v, std::string const& textual)
    {
        // setting a default value sets the origin to that value
        origin_ = v;
        _Base::default_value(v, textual);
        return this;
    }

public:
    void xparse(boost::any& v, std::vector<std::basic_string<charT> > const&) const
    {
        // if this is the first occurrence of the option, initialize it to the origin
        if (v.empty()) {
            v = boost::any(origin_);
        }
        ++boost::any_cast<T&>(v);
    }

private:
    /** numeric origin from which to increment upward */
    T origin_;
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

} // namespace halmd

#endif /* ! HALMD_UTILITY_PROGRAM_OPTIONS_HPP */
