/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   functional.hpp
 * \author Andrey Semashev
 * \date   30.03.2008
 *
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_DETAIL_FUNCTIONAL_HPP_INCLUDED_
#define BOOST_LOG_DETAIL_FUNCTIONAL_HPP_INCLUDED_

#include <algorithm>
#include <functional>
#include <boost/regex.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_unsigned.hpp>
#include <boost/log/detail/prologue.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace aux {

//! Comparison implementation for integral types
template< template< typename > class PredT, typename T, typename U >
inline bool compare_integral(T const& left, U const& right)
{
    // We do these tricks mostly to silence compiler warnings like 'signed/unsigned mismatch'
    typedef typename mpl::eval_if_c<
        (sizeof(T) > sizeof(U)),
        mpl::identity< T >,
        mpl::eval_if_c<
            (sizeof(U) > sizeof(T)),
            mpl::identity< U >,
            mpl::if_<
                is_unsigned< T >,
                T,
                U
            >
        >
    >::type actual_type;
    return PredT< actual_type >()(left, right);
}

//! Equality predicate
struct equal_to
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& left, U const& right) const
    {
        return op(left, right, typename mpl::and_< is_integral< T >, is_integral< U > >::type());
    }

private:
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::false_ const&)
    {
        return (left == right);
    }
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::true_ const&)
    {
        return compare_integral< std::equal_to >(left, right);
    }
};

//! Inequality predicate
struct not_equal_to
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& left, U const& right) const
    {
        return op(left, right, typename mpl::and_< is_integral< T >, is_integral< U > >::type());
    }

private:
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::false_ const&)
    {
        return (left != right);
    }
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::true_ const&)
    {
        return compare_integral< std::not_equal_to >(left, right);
    }
};

//! Less predicate
struct less
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& left, U const& right) const
    {
        return op(left, right, typename mpl::and_< is_integral< T >, is_integral< U > >::type());
    }

private:
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::false_ const&)
    {
        return (left < right);
    }
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::true_ const&)
    {
        return compare_integral< std::less >(left, right);
    }
};

//! Greater predicate
struct greater
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& left, U const& right) const
    {
        return op(left, right, typename mpl::and_< is_integral< T >, is_integral< U > >::type());
    }

private:
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::false_ const&)
    {
        return (left > right);
    }
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::true_ const&)
    {
        return compare_integral< std::greater >(left, right);
    }
};

//! Less or equal predicate
struct less_equal
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& left, U const& right) const
    {
        return op(left, right, typename mpl::and_< is_integral< T >, is_integral< U > >::type());
    }

private:
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::false_ const&)
    {
        return (left <= right);
    }
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::true_ const&)
    {
        return compare_integral< std::less_equal >(left, right);
    }
};

//! Greater or equal predicate
struct greater_equal
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& left, U const& right) const
    {
        return op(left, right, typename mpl::and_< is_integral< T >, is_integral< U > >::type());
    }

private:
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::false_ const&)
    {
        return (left >= right);
    }
    template< typename T, typename U >
    static bool op(T const& left, U const& right, mpl::true_ const&)
    {
        return compare_integral< std::greater_equal >(left, right);
    }
};

//! The in_range functor
struct in_range_fun
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& value, std::pair< U, U > const& rng) const
    {
        return op(value, rng, typename mpl::and_< is_integral< T >, is_integral< U > >::type());
    }

private:
    template< typename T, typename U >
    static bool op(T const& value, std::pair< U, U > const& rng, mpl::false_ const&)
    {
        return (value >= rng.first && value < rng.second);
    }
    template< typename T, typename U >
    static bool op(T const& value, std::pair< U, U > const& rng, mpl::true_ const&)
    {
        return (compare_integral< std::greater_equal >(value, rng.first)
            && compare_integral< std::less >(value, rng.second));
    }
};

//! The begins_with functor
struct begins_with_fun
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& left, U const& right) const
    {
        if (left.size() >= right.size())
            return std::equal(right.begin(), right.end(), left.begin());
        else
            return false;
    }
};

//! The ends_with functor
struct ends_with_fun
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& left, U const& right) const
    {
        if (left.size() >= right.size())
            return std::equal(right.begin(), right.end(), left.end() - right.size());
        else
            return false;
    }
};

//! The contains functor
struct contains_fun
{
    typedef bool result_type;

    template< typename T, typename U >
    bool operator() (T const& left, U const& right) const
    {
        if (left.size() >= right.size())
        {
            bool result = false;
            typename T::const_iterator search_end = left.end() - right.size() + 1;
            for (typename T::const_iterator it = left.begin(); it != search_end && !result; ++it)
                result = std::equal(right.begin(), right.end(), it);

            return result;
        }
        else
            return false;
    }
};

//! The regex matching functor
struct matches_fun
{
    typedef bool result_type;

private:
    match_flag_type m_Flags;

public:
    explicit matches_fun(match_flag_type flags = match_default) : m_Flags(flags)
    {
    }

    template< typename T, typename CharT, typename RegexTraitsT >
    bool operator() (T const& str, basic_regex< CharT, RegexTraitsT > const& expr) const
    {
        return regex_match(str.begin(), str.end(), expr, m_Flags);
    }
};

//! Second argument binder
template< typename FunT, typename SecondArgT >
struct binder2nd :
    private FunT
{
    typedef typename FunT::result_type result_type;

    binder2nd(FunT const& fun, SecondArgT const& arg) : FunT(fun), m_Arg(arg) {}

    template< typename T >
    result_type operator() (T const& arg) const
    {
        return FunT::operator() (arg, m_Arg);
    }

private:
    SecondArgT m_Arg;
};

template< typename FunT, typename SecondArgT >
inline binder2nd< FunT, SecondArgT > bind2nd(FunT const& fun, SecondArgT const& arg)
{
    return binder2nd< FunT, SecondArgT >(fun, arg);
}

} // namespace aux

} // namespace log

} // namespace boost

#endif // BOOST_LOG_DETAIL_FUNCTIONAL_HPP_INCLUDED_
