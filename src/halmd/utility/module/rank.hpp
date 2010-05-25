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

#ifndef HALMD_UTILITY_MODULE_RANK_HPP
#define HALMD_UTILITY_MODULE_RANK_HPP

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

/**
 * Exceptions for base-to-derived rank signaling.
 */
struct _rank_left_is_base {};
struct _rank_right_is_base {};

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
template <typename T = void, typename Enable = void>
class rank;

/**
 * A type-independent base class, which is used in the factory to
 * store pointers of arbitrary rank instances in a container.
 */
template <>
class rank<>
{
public:
    typedef shared_ptr<rank<> > _Rank_ptr;

    virtual ~rank() {}
    virtual bool operator<(_Rank_ptr const& other) const = 0;
    virtual std::string name() const = 0;

private:
    template <typename T, typename Enable> friend class rank;
    virtual bool operator>(_Rank_ptr const& other) const = 0;
};

/**
 * A base template, which corresponds to the base class of a
 * module (e.g. force, particle).
 */
template <typename T, typename Enable>
class rank
  : public rank<>
{
public:
    typedef T _Module_base;
    typedef shared_ptr<rank<> > _Rank_ptr;

    virtual bool operator<(_Rank_ptr const& other) const
    {
        if (dynamic_pointer_cast<rank<T> >(other)) {
            throw _rank_left_is_base();
        }
        return !(other->operator>(_Rank_ptr(new rank<T>)));
    }

    virtual bool operator>(_Rank_ptr const& other) const
    {
        if (dynamic_pointer_cast<rank<T> >(other)) {
            throw _rank_right_is_base();
        }
        return typeid(*other).before(typeid(rank<T>));
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
class rank<T, typename enable_if<is_object<typename T::_Base> >::type>
  : public rank<typename T::_Base>
{
public:
    typedef typename T::_Base _Base;
    typedef shared_ptr<rank<> > _Rank_ptr;

    /**
     * The less-than operator is called by a rank ordering class
     * and marks the first phase of strict weak ordering.
     */
    virtual bool operator<(_Rank_ptr const& other) const
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

        try {
            return rank<_Base>::operator<(other);
        }
        catch (_rank_left_is_base const&)
        {
            if (dynamic_pointer_cast<rank<T> >(other)) {
                throw;
            }
            return !(other->operator>(_Rank_ptr(new rank<T>)));
        }
    }

    /**
     * The greater-than operator is called by the less-than
     * operator of the other rank and marks the second phase.
     */
    virtual bool operator>(_Rank_ptr const& other) const
    {
        // This essentially does the same as the less-than
        // operator, but for the rank on the *right* side of the
        // inequality. If the most derived base rank is reached,
        // we order the two ranks by their collation order, which
        // is compiler-dependent but deterministic.

        try {
            return rank<_Base>::operator>(other);
        }
        catch (_rank_right_is_base const&)
        {
            if (dynamic_pointer_cast<rank<T> >(other)) {
                throw;
            }
            return typeid(*other).before(typeid(rank<T>));
        }
    }

    virtual std::string name() const
    {
        return demangled_name<T>();
    }
};

/**
 * A class to implement strict weak ordering of rank instances
 * in STL sequences, following the rules described above.
 */
struct rank_order
{
    typedef shared_ptr<rank<> > _Rank_ptr;

    bool operator()(_Rank_ptr left, _Rank_ptr right) const
    {
        if (typeid(*left) != typeid(*right)) {
            try {
                return left->operator<(right);
            }
            catch (_rank_left_is_base const&) {
                return false;
            }
            catch (_rank_right_is_base const&) {
                return true;
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
struct rank_order_equal_base
{
    typedef shared_ptr<rank<> > _Rank_ptr;

    bool operator()(_Rank_ptr left, _Rank_ptr right) const
    {
        if (typeid(*left) != typeid(*right)) {
            try {
                return left->operator<(right);
            }
            catch (_rank_left_is_base const&) {
                return false;
            }
            catch (_rank_right_is_base const&) {
                return false;
            }
        }
        return false;
    }
};

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_RANK_HPP */
