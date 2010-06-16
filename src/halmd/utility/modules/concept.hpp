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

#ifndef HALMD_UTILITY_MODULES_CONCEPT_HPP
#define HALMD_UTILITY_MODULES_CONCEPT_HPP

#include <boost/concept_check.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/utility/options.hpp>

namespace halmd
{
namespace modules
{

/**
 * The ModuleConcept is a template to check at compile time
 * whether a module class fulfills the requirements of a module.
 *
 * It is intended to aid the programmer in detecting undefined or
 * misdefined types, or missing or unintentionally inherited
 * functions.
 */
template <typename T, typename Enable = void>
struct ModuleConcept
{
    typedef typename T::_Self _Self;

    // check function existence and signature
    template <typename Signature, Signature* Base>
    struct Function {};

    BOOST_CONCEPT_USAGE(ModuleConcept)
    {
        BOOST_MPL_ASSERT(( boost::is_same<T, _Self> )); // _Self points to the class itself
        Function<void (), &T::depends>();
        Function<void (po::options_description&), &T::options>();
        Function<void (po::options const&), &T::select>();
    }
};

template <typename T>
struct ModuleConcept<T, typename boost::enable_if<boost::is_object<typename T::_Base> >::type>
  : ModuleConcept<typename T::_Base>
{
    typedef typename T::_Self _Self;
    typedef typename T::_Base _Base;

    // check for different functions in derived and base class
    template <typename Signature, Signature* Derived, Signature* Base>
    struct FunctionNotFromBase {};

    // declared but not defined, to cause error in case of identical functions
    template <typename Signature, Signature* Base>
    struct FunctionNotFromBase<Signature, Base, Base>;

    BOOST_CONCEPT_USAGE(ModuleConcept)
    {
        BOOST_MPL_ASSERT(( boost::is_same<_Self, T> ));
        BOOST_MPL_ASSERT(( boost::is_base_of<_Base, T> )); // _Base is (any) base of class
        FunctionNotFromBase<void (), &T::depends, &_Base::depends>();
        FunctionNotFromBase<void (po::options_description&), &T::options, &_Base::options>();
        FunctionNotFromBase<void (po::options const&), &T::select, &_Base::select>();
    }
};

} // namespace modules

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_CONCEPT_HPP */
