/*
 * Copyright © 2013 Nicolas Höft
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

#include <cstdlib>
#include <ctime>
#include <limits>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/banded.hpp>

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE linear
#include <boost/test/unit_test.hpp>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/multi_index.hpp>
#include <halmd/utility/tuple.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/forces/interpolation/linear.hpp>
#ifdef HALMD_WITH_GPU
# include <test/tools/cuda.hpp>
#endif

#include <test/tools/ctest.hpp>


template<int dimension, typename float_type>
struct cubic_function
{
    typedef halmd::fixed_vector<float_type, dimension> vector_type;

    std::vector<float_type> coefficients(vector_type r) {
        std::vector<float_type> c;
        float_type x(0), y(0), z(0);

        x = r[0];
        if(dimension > 1)
            y = r[1];
        if(dimension > 2)
            z = r[2];

        c.push_back(f(x, y, z));
        c.push_back(fdx(x, y, z));
        if (dimension > 1) {
            c.push_back(fdy(x, y, z));
            c.push_back(fdxdy(x, y, z));
        }
        if (dimension > 2) {
            c.push_back(fdz(x, y, z));
            c.push_back(fdxdz(x, y, z));
            c.push_back(fdydz(x, y, z));
            c.push_back(fdxdydz(x, y, z));
        }

        return c;
    }

    float_type f(float_type x, float_type y, float_type z)
    {
        return 5*x*x*x - 3*y*y*y + z*z*z + 2*x*x*y - y*y*z - 2*z*z + x*y*z + 35;
    }

    float_type fdx(float_type x, float_type y, float_type z)
    {
        return 15*x*x + 4*x*y + y*z;
    }

    float_type fdy(float_type x, float_type y, float_type z)
    {
        return -9*y*y + 2*x*x - 2*y*z + x*z;
    }

    float_type fdz(float_type x, float_type y, float_type z)
    {
        return 3*z*z - y*y - 4*z + x*y;
    }

    float_type fdxdy(float_type x, float_type y, float_type z)
    {
        return 4*x + z;
    }

    float_type fdxdz(float_type x, float_type y, float_type z)
    {
        return y;
    }

    float_type fdydz(float_type x, float_type y, float_type z)
    {
        return -2*y + x;
    }

    float_type fdxdydz(float_type x, float_type y, float_type z)
    {
        return 1;
    }
};

BOOST_AUTO_TEST_SUITE( host )

BOOST_AUTO_TEST_CASE( linear )
{
    enum { dimension = 2 };
    typedef double float_type;
    typedef halmd::mdsim::forces::interpolation::linear<dimension, float_type> interpolation_type;
    cubic_function<dimension, float_type> fn;
    typedef typename interpolation_type::vector_type vector_type;
    typedef typename interpolation_type::vector_type force_vector_type;
    typedef typename interpolation_type::index_type index_type;
    typedef typename interpolation_type::size_type size_type;
    typedef halmd::mdsim::box<dimension> box_type;
    std::vector<float_type> coefficients;

    // set up a box with a length of 3x3x3 and nodes at (0,0,0) to (3,3,3)
    std::shared_ptr<box_type> box;
    float box_length = 3;
    index_type ngrid_points(4);

    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = box_length;
    }
    box = std::make_shared<box_type>(edges);

    interpolation_type interpolation_(box->length(), box->origin(), ngrid_points);
    vector_type grid_basis = interpolation_.grid_basis();
    size_type total_grid_points = interpolation_.total_knots();

    coefficients.resize(total_grid_points * interpolation_type::coefficients_per_knot);

    for (size_type i = 0; i < total_grid_points; ++i) {

        index_type nb_idx = halmd::offset_to_multi_index(i, ngrid_points);
        vector_type nb_pos = element_prod(static_cast<vector_type>(nb_idx), grid_basis);

        int const cidx = interpolation_type::coefficients_per_knot * halmd::multi_index_to_offset(nb_idx, ngrid_points);
        // shift the neighbour position into the box
        nb_pos += box->origin();

        std::vector<float_type> c = fn.coefficients(nb_pos);
        for (unsigned int j = 0; j < interpolation_type::coefficients_per_knot; ++j) {
            coefficients[cidx + j] = c[j];
        }
    }

    float_type const eps = std::numeric_limits<float_type>::epsilon();

    //std::srand(std::time(0));
    float_type max_eps = 0;
    int max_points = 50;
    int total_points = pow(max_points, dimension);
    index_type max_points_vec(max_points);

    index_type old_ridx(0);
    for (int i = 0; i < total_points; ++i) {
        vector_type r;
        index_type ridx = halmd::offset_to_multi_index(i, max_points_vec);

        for (int d = 0; d < dimension; ++d) {
            r[d] = (ridx[d]) * box->length()[d] / max_points_vec[d] - box->length()[d]/2;
        }
        if(old_ridx[dimension-1] != ridx[dimension-1] && dimension > 1)
            std::cout << "\n";
        old_ridx = ridx;
        //std::cout << ridx << " " << r << std::endl;
        // put the particle back into the box
        box->reduce_periodic(r);

        float_type x(0), y(0), z(0);
        x = r[0];
        if (dimension > 1)
            y = r[1];
        if (dimension > 2)
            z = r[2];
        float_type epot_analytic = fn.f(x, y, z);
        /*vector_type force_analytic;
        force_analytic[0] = -fn.fdx(x, y, z);
        if (dimension > 1)
            force_analytic[1] = -fn.fdy(x, y, z);
        if (dimension > 2)
            force_analytic[2] = -fn.fdz(x, y, z);*/

        float_type epot_interpol;
        force_vector_type force_interpol;
        halmd::tie(epot_interpol, force_interpol) = interpolation_.operator()<force_vector_type>(r, &*coefficients.begin());

        //BOOST_CHECK_SMALL(epot_interpol - epot_analytic, eps);
        max_eps = std::max(epot_interpol - epot_analytic, max_eps);
        std::cout << r << " " << epot_analytic << " " << epot_interpol << std::endl;

        /*for (int d = 0; d < dimension; ++d) {
            max_eps = std::max(force_interpol[d] - force_analytic[d], max_eps);
            //BOOST_CHECK_SMALL(force_interpol[d] - force_analytic[d], eps);
//             if(force_interpol[d] - force_analytic[d] > eps)
//             {
//                 std::cout << r << std::endl;
//             }
        }*/
    }
    //std::cout << "max eps: " << max_eps << std::endl;
    //std::cout << "factor : " << max_eps / eps << std::endl;

}

BOOST_AUTO_TEST_SUITE_END() // host

