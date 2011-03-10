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

#include <boost/bind.hpp>

#include <test/performance/function_calls_noop.hpp>

using namespace boost;
using namespace std;

/**
 * Make function to noop using boost::bind
 */
function<void (double)> bind_noop()
{
    return bind(&noop, _1);
}

/**
 * Add noop function to signal
 */
void bind_noop(halmd::signal<void (double)>& sig)
{
    sig.connect(bind_noop());
}

/**
 * An empty function
 */
void noop(double) {}
