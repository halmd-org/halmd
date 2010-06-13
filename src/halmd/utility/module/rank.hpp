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

#ifndef HALMD_UTILITY_MODULE_RESOLVER_HPP
#define HALMD_UTILITY_MODULE_RESOLVER_HPP

#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_object.hpp>
#include <boost/utility/enable_if.hpp>
#include <string>

#include <halmd/utility/module/demangle.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace utility { namespace module
{

// import into namespace
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::enable_if;
using boost::is_object;

// forward declaration
class untyped_rank_base;

/**
 * The rank template models a forest of disjunct class hierarchy
 * trees with a strict weak ordered STL container. This allows
 * lookups with logarithmic complexity in a linear sequence.
 *
 * Two rank classes should be ordered thus that a derived
 * class comes before its base class, and any subtree of
 * a class hierarchy tree is a consecutive subsequence of the
 * sequence.
 *
 * The corresponding rules for ordering are, with descending
 * priority, (1) a derived class comes before its base class,
 * (2) classes with a common base class are ordered according
 * to the collation order of the direct descendants of this
 * base class, and (3) classes without a common base class are
 * ordered according to the collation order of their respective
 * base classes.
 */

typedef shared_ptr<untyped_rank_base> rank;

/**
 * A type-independent base class, which is used in the factory to
 * store pointers of arbitrary rank instances in a container.
 */
class untyped_rank_base
{
public:
    // C++0x will support scoped and strongly typed enums.
    // http://www2.research.att.com/~bs/C++0xFAQ.html#enum
    struct result
    {
        enum type
        {
            left,      // left is before right rank in collation order
            right,
            left_base, // left is base of right rank
            right_base,
        };
    };

    virtual ~untyped_rank_base() {}
    virtual result::type operator<(rank const& other) const = 0;
    virtual result::type operator>(rank const& other) const = 0;
    virtual std::string name() const = 0;
};

/**
 * A base template, which corresponds to the base class of a
 * module (e.g. force, particle).
 */
template <typename T, typename Enable = void>
class typed_rank_base
  : public untyped_rank_base
{
public:
    typedef untyped_rank_base _Base;
    typedef typename _Base::result result;
    typedef T _Base_type;

    virtual typename result::type operator<(rank const& other) const
    {
        if (dynamic_pointer_cast<typed_rank_base<T> >(other)) {
            return result::left_base;
        }
        return other->operator>(rank(new typed_rank_base<T>));
    }

    virtual typename result::type operator>(rank const& other) const
    {
        if (dynamic_pointer_cast<typed_rank_base<T> >(other)) {
            return result::right_base;
        }
        if (typeid(*other).before(typeid(typed_rank_base<T>))) {
            return result::left;
        }
        return result::right;
    }

    virtual std::string name() const
    {
        return demangled_name<T>();
    }
};

/**
 * A deriving template to the base template.
 */
template <typename T>
class typed_rank_base<T, typename enable_if<is_object<typename T::_Base> >::type>
  : public typed_rank_base<typename T::_Base>
{
public:
    typedef typed_rank_base<typename T::_Base> _Base;
    typedef typename _Base::result result;

    /**
     * The less-than operator is called by a rank ordering class
     * and marks the first phase of strict weak ordering.
     */
    virtual typename result::type operator<(rank const& other) const
    {
        // We cast the other rank to every base rank that the
        // rank on the *left* side of the inequality (directly or
        // indirectly) derives from, starting with the top-most
        // base rank and going downward to derived base ranks.
        //
        // Any common base class will trigger a rank base
        // exception, which is caught by the directly derived
        // rank. When the most derived base rank is reached, a
        // rank of the next derived rank is passed to the
        // greater-than operator of the other rank.

        typename result::type const result = _Base::operator<(other);
        if (result == result::left_base) {
            if (dynamic_pointer_cast<typed_rank_base<T> >(other)) {
                return result::left_base;
            }
            return other->operator>(rank(new typed_rank_base<T>));
        }
        return result;
    }

    /**
     * The greater-than operator is called by the less-than
     * operator of the other rank and marks the second phase.
     */
    virtual typename result::type operator>(rank const& other) const
    {
        // This essentially does the same as the less-than
        // operator, but for the rank on the *right* side of the
        // inequality. If the most derived base rank is reached,
        // we order the two ranks by their collation order, which
        // is compiler-dependent but deterministic.

        typename result::type const result = _Base::operator>(other);
        if (result == result::right_base) {
            if (dynamic_pointer_cast<typed_rank_base<T> >(other)) {
                return result::right_base;
            }
            if (typeid(*other).before(typeid(typed_rank_base<T>))) {
                return result::left;
            }
            return result::right;
        }
        return result;
    }

    virtual std::string name() const
    {
        return demangled_name<T>();
    }
};

template <typename T>
class typed_rank
  : public typed_rank_base<T>
{
public:
    typedef typed_rank_base<T> _Base;
};

/**
 * A class to implement strict weak ordering of rank instances
 * in STL sequences, following the rules described above.
 */
struct compare_rank
{
    typedef untyped_rank_base::result result;

    bool operator()(rank left, rank right) const
    {
        if (typeid(*left) != typeid(*right)) {
            switch (left->operator<(right)) {
              case result::left:
                return false;
              case result::right:
                return true;
              case result::left_base:
                return false;
              case result::right_base:
                return true;
              // compiler emits warning if enumeration value is unhandled
            }
        }
        return false;
    }
};

/**
 * A class to implement strict weak ordering of rank instances
 * in STL sequences, following the rules described above, with
 * the exception of defining a base and derived ranks as equal.
 */
struct compare_rank_base
{
    typedef untyped_rank_base::result result;

    bool operator()(rank left, rank right) const
    {
        if (typeid(*left) != typeid(*right)) {
            switch (left->operator<(right)) {
              case result::left:
                return false;
              case result::right:
                return true;
              case result::left_base:
                return false;
              case result::right_base:
                return false;
              // compiler emits warning if enumeration value is unhandled
            }
        }
        return false;
    }
};

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_RESOLVER_HPP */
