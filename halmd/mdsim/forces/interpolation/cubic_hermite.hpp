/*
 * Copyright © 2013  Nicolas Höft
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

#ifndef HALMD_MDSIM_FORCES_INTERPOLATION_CUBIC_HERMITE_HPP
#define HALMD_MDSIM_FORCES_INTERPOLATION_CUBIC_HERMITE_HPP

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
 * This is the cubic hermite interpolation scheme.
 * It assumes and enforces a uniform grid.
 */
template <int dimension, typename float_type>
class cubic_hermite
{
public:
    enum {
        coefficients_per_knot = (1 << dimension)
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
    cubic_hermite(
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

    // prefactor for f(x) where x_i > x
    HALMD_GPU_ENABLED
    float_type h10(float_type x) const
    {
        return x * x * (3 - 2*x);
    }

    // prefactor for f'(x_i), where x_i > x
    HALMD_GPU_ENABLED
    float_type h11(float_type x, float_type dx) const
    {
        return dx * x * x * (x - 1);
    }

    HALMD_GPU_ENABLED
    float_type hd10(float_type x) const
    {
        return 6*x*(1 - x);
    }

    HALMD_GPU_ENABLED
    float_type hd11(float_type x, float_type dx) const
    {
        return dx * x * (3*x - 2);
    }

    template <int bitmask, typename force_vector_type>
    HALMD_GPU_ENABLED
    void hermite_factors(vector_type const rd
                       , unsigned int const nb
                       , float_type const coefficient
                       , float_type& epot
                       , force_vector_type& f
    ) const;

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
template <int bitmask, typename force_vector_type>
HALMD_GPU_ENABLED
void cubic_hermite<dimension, float_type>::hermite_factors(vector_type const rd
                  , unsigned int const nb
                  , float_type const coefficient
                  , float_type& epot
                  , force_vector_type& f
) const
{
    vector_type h;
    vector_type hd;

    // The bitmask gives the information about whether the coefficient
    // is the derivative of the interpolated function in x, y or z
    // direction (or any combination of them), indicated by a set bit
    // (first bit for x derivate, second for y and so on).
    //
    // This dertermines the prefactors for the coefficient (either h_j^0 or
    // h_j^1, the latter in case it is a derivative).
    // To determine the second parameter j, the neighbour index (nb) is needed.
    // nb's function is identical to the bitmask and if a bit is set, it
    // tells this function whether the neighbour position is larger or smaller
    // than the coordinate to interpolate.

    for (int d = 0; d < dimension; ++d) {
        // test if the d-th bit is set
        if (bitmask & (1 << d)) {
            // coefficient is a derivative in dimension 'd'
            if((nb & (1 << d)) == 0) {
                // neighbour position is smaller than r
                h[d] = h11(1 - rd[d], -grid_basis_[d]);
                hd[d] = -hd11(1 - rd[d], -grid_basis_[d]);
            }
            else {
                // neighbour position is larger than r
                h[d] = h11(rd[d], grid_basis_[d]);
                hd[d] = hd11(rd[d], grid_basis_[d]);
            }
        }
        else {
            // coefficient is not a derivative in dimension 'd'
            if((nb & (1 << d)) == 0) {
                // neighbour position is smaller than r
                h[d]  = h10(1 - rd[d]);
                hd[d] = -hd10(1 - rd[d]);
            }
            else {
                // neighbour position is larger than r
                h[d]  = h10(rd[d]);
                hd[d] = hd10(rd[d]);
            }
        }
    }

    if(dimension == 1) {
        epot += h[0] * coefficient;
        f[0] += hd[0] * coefficient;
    }
    else if(dimension == 2) {
        epot += h[0]*h[1] * coefficient;
        f[0] += hd[0]*h[1] * coefficient;
        f[1] += h[0]*hd[1] * coefficient;
    }
    else if(dimension == 3) {
        epot += h[0]*h[1]*h[2] * coefficient;
        f[0] += hd[0]*h[1]*h[2] * coefficient;
        f[1] += h[0]*hd[1]*h[2] * coefficient;
        f[2] += h[0]*h[1]*hd[2] * coefficient;
    }
}

//
// The cubic hermite scheme is a 3rd order interpolation scheme
// that
//
//    E_pot(r) =  h_0^0(r) * E_{pot}(r_{i}) + h_1^0(r) * E_{pot}(r_{i+1})
//              + h_0^1(r) * F(r_{i})       + h_1^1(r) * F(r_{i+1})
//
// with
//
//   h01(r) = r² * (3 - 2*r)
//   h11(r) = Δx * r² (r - 1)
//
// constraints are sorted in the following way for each neighbour point
//  - function value
//  - d/dx
//  - d/dy
//  - d²/dxdy
//  - d/dz
//  - d²/dxdz
//  - d²/dydz
//  - d³/dxdydz
//
template <int dimension, typename float_type>
template <typename force_vector_type, typename position_type, typename coefficient_type>
HALMD_GPU_ENABLED tuple<float_type, force_vector_type>
cubic_hermite<dimension, float_type>::operator()(
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

    for (int i = 0; i < neighbours; ++i) {
        index_type nb_idx;
        for (int d = 0; d < dimension; ++d) {
            nb_idx[d] = (first_nb_idx[d] + ((i >> d) & 1));
        }
        int cidx = coefficients_per_knot * multi_index_to_offset(nb_idx, nknots_);
        float_type c0, c1, c2, c3;
#ifdef __CUDACC__
        float4 cc = *((float4*)&coefficients[cidx]);
        c0 = cc.x;
        c1 = cc.y;
        c2 = cc.z;
        c3 = cc.w;
#else
        c0 = coefficients[cidx + 0];
        c1 = coefficients[cidx + 1];
        c2 = coefficients[cidx + 2];
        c3 = coefficients[cidx + 3];
#endif

        hermite_factors<0>(rd, i, c0, epot, f);
        hermite_factors<1>(rd, i, c1, epot, f);
        if (dimension > 1) {
            hermite_factors<2>(rd, i, c2, epot, f);
            hermite_factors<3>(rd, i, c3, epot, f);
        }
        if (dimension > 2) {
#ifdef __CUDACC__
            cc = *((float4*)&coefficients[cidx + 4]);
            c0 = cc.x;
            c1 = cc.y;
            c2 = cc.z;
            c3 = cc.w;
#else
            c0 = coefficients[cidx + 4];
            c1 = coefficients[cidx + 5];
            c2 = coefficients[cidx + 6];
            c3 = coefficients[cidx + 7];
#endif
            hermite_factors<4>(rd, i, c0, epot, f);
            hermite_factors<5>(rd, i, c1, epot, f);
            hermite_factors<6>(rd, i, c2, epot, f);
            hermite_factors<7>(rd, i, c3, epot, f);
        }
    }

    // the force is  F(x) = -∇V(x)
    for (int d = 0; d < dimension; ++d) {
        f[d] = -f[d]/grid_basis_[d];
    }

    return make_tuple(epot, f);
}

} // namespace interpolation
} // namespace forces
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_FORCES_INTERPOLATION_CUBIC_HERMITE_HPP */
