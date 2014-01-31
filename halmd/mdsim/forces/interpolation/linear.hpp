/*
 * Copyright © 2014  Nicolas Höft
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

#ifndef HALMD_MDSIM_FORCES_INTERPOLATION_LINEAR_HPP
#define HALMD_MDSIM_FORCES_INTERPOLATION_LINEAR_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/multi_index.hpp>
#include <halmd/utility/tuple.hpp>

#ifndef __CUDACC__
# include <lua.hpp>
#endif

namespace halmd {
namespace mdsim {
namespace forces {
namespace interpolation {

/**
 * This is the linear interpolation scheme.
 * It assumes and enforces a uniform grid.
 */
template <int dimension, typename float_type>
class linear
{
public:
    enum {
        coefficients_per_knot = 1
        // each dimension doubles the number of next neighbours
      , neighbours = (1 << dimension)
    };

    typedef unsigned int size_type;
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef fixed_vector<size_type, dimension> index_type;

    /**
     * Calculate the interpolated value at position r
     * using the provided grid knot coefficients.
     */
    template <typename force_vector_type, typename position_type, typename coefficient_type>
    HALMD_GPU_ENABLED
    tuple<float_type, force_vector_type>
    operator()(
        position_type const& r
      , coefficient_type* coefficients
    ) const;

#ifndef __CUDACC__
    linear(
        vector_type length
      , vector_type origin
      , index_type nknots
    );

    /**
     * Bind class to Lua
     */
    static void luaopen(lua_State* L);
#endif // !__CUDACC__

    /**
     * Total number of grid knots.
     */
    size_type total_knots() const
    {
        return total_knots_;
    }

    /**
     * Number of knots in each dimension.
     */
    index_type nknots() const
    {
        return nknots_;
    }

    /**
     * Return the basis vector of the grid cell.
     */
    vector_type grid_basis() const
    {
        return grid_basis_;
    }

private:
    template<int neighbour, typename coefficient_type>
    HALMD_GPU_ENABLED
    float_type linear_factor(vector_type const r, index_type const first_neighbour_index, coefficient_type* coefficients) const;

    /** number of grid points in each spatial direction */
    index_type nknots_;
    /** origin of the box */
    vector_type origin_;
    /** distances between nodes */
    vector_type grid_basis_;
    /** edge length of the domain */
    vector_type cell_length_;
    /** total number of nodes */
    size_type total_knots_;
};

template <int dimension, typename float_type>
template <int neighbour, typename coefficient_type>
HALMD_GPU_ENABLED
float_type linear<dimension, float_type>::linear_factor(vector_type const r, index_type const first_neighbour_index, coefficient_type* coefficients) const
{
    index_type nb_idx;
    float_type factor = 1;
    for (int d = 0; d < dimension; ++d) {
        nb_idx[d] = (first_neighbour_index[d] + ((neighbour >> d) & 1));
        // if d'th bit is set, the current neighbours position
        // is larger than the point to be interpolated
        if(neighbour & (1 << d)) {
            factor *= r[d];
        }
        else {
            factor *= (1 - r[d]);
        }
    }
    int coefficient_idx = coefficients_per_knot * multi_index_to_offset(nb_idx, nknots_);
    return factor*coefficients[coefficient_idx];
}

//
// The linear scheme is a first order interpolation scheme
// that works the following
//
//    E_pot(x) = (1 - x) * E_{pot}(x_{i}) + x * E_{pot}(x_{i+1})
//
// constraints are sorted in the following way for each neighbour point
//  - function value
//
template <int dimension, typename float_type>
template <typename force_vector_type, typename position_type, typename coefficient_type>
HALMD_GPU_ENABLED tuple<float_type, force_vector_type>
linear<dimension, float_type>::operator()(
    position_type const& r
  , coefficient_type* coefficients
) const
{
    float_type epot(0);
    force_vector_type f(0);

    // reduced distance
    vector_type rd;
    index_type first_nb_idx;

    // Calculate the 'first' neighbour, i.e. the neighbour in the lower left edge.
    // All other neighbours have at least in one spatial direction a higher index.
    for (int d = 0; d < dimension; ++d) {
        // Shift the particle position that all coordinates are in (0, L)
        // and not in (-L/2, L/2) as given by the force module.
        // The conversation from position_type to vector_type is necessary because it can be
        // of a different type than the float_type used in the interpolation functor
        // e.g. position_type is of type fixed_vector<float> while the interpolation uses
        // fixed_vector<double>
        float_type r_ = static_cast<float_type>(r[d]) - origin_[d];
        // If the particle is outside the unit cell of the grid, shift it back into it.
        // This allows having multiple, periodic grid boxes.
        r_ -= cell_length_[d] * static_cast<unsigned int>(r_ / cell_length_[d]);
        first_nb_idx[d] = static_cast<size_type>(r_ / grid_basis_[d]);

        // reduced distance, the parametrisation for the spline between the knots
        // is in [0,1]
        rd[d] = (r_ - first_nb_idx[d] * grid_basis_[d]) / grid_basis_[d];
    }

    epot += linear_factor<0>(rd, first_nb_idx, coefficients);
    epot += linear_factor<1>(rd, first_nb_idx, coefficients);
    if (dimension > 1) {
        epot += linear_factor<2>(rd, first_nb_idx, coefficients);
        epot += linear_factor<3>(rd, first_nb_idx, coefficients);
    }
    if (dimension > 2) {
        epot += linear_factor<4>(rd, first_nb_idx, coefficients);
        epot += linear_factor<5>(rd, first_nb_idx, coefficients);
        epot += linear_factor<6>(rd, first_nb_idx, coefficients);
        epot += linear_factor<7>(rd, first_nb_idx, coefficients);
    }

    return make_tuple(epot, f);
}

} // namespace interpolation
} // namespace forces
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_FORCES_INTERPOLATION_LINEAR_HPP */
