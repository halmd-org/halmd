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

#define BOOST_TEST_MODULE init
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include <halmd/utility/init.hpp>

int& scalar()
{
    static int i = 1;
    return i;
}

std::vector<double>& vector()
{
    static std::vector<double> v;
    return v;
}

BOOST_AUTO_TEST_CASE( init_scalar )
{
    BOOST_CHECK_EQUAL( scalar(), 42 );
}

template <typename T>
typename std::vector<T>::const_iterator find(std::vector<T> const& v, T const& t)
{
    return std::find(v.begin(), v.end(), t);
}

BOOST_AUTO_TEST_CASE( init_vector )
{
    BOOST_CHECK_EQUAL( vector().size(), 6LU );
    BOOST_CHECK( find(vector(), M_PI_2) == (find(vector(), M_PI) + 1) );
    BOOST_CHECK( find(vector(), M_PI_4) == (find(vector(), M_PI) + 2) );
    BOOST_CHECK( find(vector(), M_PI)   != vector().end() );
    BOOST_CHECK( find(vector(), M_E)    != vector().end() );
    BOOST_CHECK( find(vector(), M_LN2)  != vector().end() );
    BOOST_CHECK( find(vector(), M_1_PI) != vector().end() );
}

HALMD_INIT( set_scalar )
{
    scalar() = 42;
}

HALMD_INIT( append_vector1 )
{
    vector().push_back(M_PI);
    vector().push_back(M_PI_2);
    vector().push_back(M_PI_4);
}

HALMD_INIT( append_vector2 )
{
    vector().push_back(M_E);
    vector().push_back(M_LN2);
}

HALMD_INIT( append_vector3 )
{
    vector().push_back(M_1_PI);
}
