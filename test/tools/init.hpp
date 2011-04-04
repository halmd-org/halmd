/*
 * Copyright © 2011  Peter Colberg
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

#ifndef HALMD_TEST_INIT_HPP
#define HALMD_TEST_INIT_HPP

/**
   Portable alternative to __attribute__((constructor)).

   This macro may be used to define a function without arguments to be called
   at program startup, *before* main(). This may for example be useful to
   register parametrized unit tests with the master test suite.

   Example usage:

   HALMD_TEST_INIT( init_unit_test_suite )
   {
        // ...
   }

   Note that the order of static initialization is undefined. For consequences
   refer to the C++ FAQ entry about “static initialization order fiasco”.

   http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.14
 */
#define HALMD_TEST_INIT( name ) namespace { struct name { name(); } _ ## name; } name::name()

#endif /* ! HALMD_TEST_INIT_HPP */
