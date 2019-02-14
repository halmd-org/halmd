/*
 * Copyright © 2011  Peter Colberg
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

#ifndef HALMD_TEST_PERFORMANCE_FUNCTION_CALLS_EXTERN_HPP
#define HALMD_TEST_PERFORMANCE_FUNCTION_CALLS_EXTERN_HPP

#include <functional>

#include <halmd/utility/signal.hpp>

std::function<void (double)> bind_noop();
void bind_noop(halmd::signal<void (double)>& sig);
void noop(double);

#endif /* ! HALMD_TEST_PERFORMANCE_FUNCTION_CALLS_EXTERN_HPP */
