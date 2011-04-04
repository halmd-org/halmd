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

#include <cmath>
#include <vector>

#include <test/tools/init.hpp>

// defined in init.cpp
std::vector<double>& vector();

// This function is defined in this separate translation unit, to
// ensure there are no symbol conflicts between translation units
// for initialization functions of the same name.
HALMD_TEST_INIT( append_vector )
{
    vector().push_back(M_E);
    vector().push_back(M_LN2);
}
