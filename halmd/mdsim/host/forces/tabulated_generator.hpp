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

#ifndef HALMD_MDSIM_HOST_FORCES_TABLULATED_GENERATOR_HPP
#define HALMD_MDSIM_HOST_FORCES_TABLULATED_GENERATOR_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/algorithm/multi_range.hpp>
#include <halmd/utility/multi_index.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <memory>

namespace halmd {
namespace mdsim {
namespace host {
namespace forces {

/**
 * This class is actually a helper class that generates the interpolation
 * coefficients for a given potential at the grid edges to be used with
 * with the tabulated force module.
 * The generated coefficients will include all possible mixed derivatives
 * (i.e. in 2D up to d²/dxdy and in 3D up to d³/dxdydz)
 */
template <int dimension, typename float_type, typename potential_type>
class tabulated_generator
{
public:
    typedef particle<dimension, float_type> particle_type;
    typedef box<dimension> box_type;
    typedef logger logger_type;

    typedef typename box_type::vector_type vector_type;
    typedef float_type coefficient_value_type;
    typedef raw_array<coefficient_value_type> coefficient_array_type;
    typedef unsigned int size_type;
    typedef fixed_vector<size_type, dimension> grid_size_type;

    tabulated_generator(
        std::shared_ptr<potential_type const> potential
      , std::shared_ptr<particle_type const> particle
      , std::shared_ptr<box_type const> box
      , grid_size_type const& grid_size
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    /**
     * Return the const reference of the interpolation coefficients
     */
    coefficient_array_type const& coefficients() const
    {
        return coefficients_;
    }

    /**
     * Return the reference of the interpolation coefficients
     */
    coefficient_array_type& coefficients()
    {
        return coefficients_;
    }

    /**
     * Calculate the interpolation coefficients for the given potenial
     */
    void compute();

    /**
     * Bind class to Lua
     */
    static void luaopen(lua_State* L);

private:
    enum { coefficients_per_knot = 1 << dimension };
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::position_type position_type;
    typedef typename particle_type::species_array_type species_array_type;
    typedef typename particle_type::species_type species_type;

    /** pair potential */
    std::shared_ptr<potential_type const> potential_;
    /** system state */
    std::shared_ptr<particle_type const> particle_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** number of grid points in each dimension */
    grid_size_type const grid_size_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;

    /** cache observer of force per particle */
    std::tuple<cache<>, cache<>> coefficient_cache_;

    /** constraints for the interpolation scheme */
    coefficient_array_type coefficients_;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type compute;
    };
    /** profiling runtime accumulators */
    runtime runtime_;
};

template <int dimension, typename float_type, typename potential_type>
tabulated_generator<dimension, float_type, potential_type>::tabulated_generator(
    std::shared_ptr<potential_type const> potential
  , std::shared_ptr<particle_type const> particle
  , std::shared_ptr<box_type const> box
  , grid_size_type const& grid_size
  , std::shared_ptr<logger_type> logger
)
  : potential_(potential)
  , particle_(particle)
  , box_(box)
  , grid_size_(grid_size)
  , logger_(logger)
{
    size_type total_knots = std::accumulate(grid_size.begin(), grid_size.end(), 1,  std::multiplies<size_type>());
    coefficients_ = coefficient_array_type(coefficients_per_knot * total_knots);
}

template <int dimension, typename float_type, typename potential_type>
void
tabulated_generator<dimension, float_type, potential_type>::compute()
{
    cache<position_array_type> const& position_cache = particle_->position();
    cache<species_array_type> const& species_cache = particle_->species();

    auto current_state = std::tie(position_cache, species_cache);

    if (coefficient_cache_ == current_state) {
        return;
    }

    position_type grid_basis;
    position_type box_length = box_->length();
    for (int i = 0; i < dimension; ++i) {
        // substract one, as the grid knot lattice contains both
        // edge points at the beginning and at the end
        grid_basis[i] = box_length[i] / (grid_size_[i] - 1);
    }

    position_array_type const& position = read_cache(particle_->position());
    species_array_type const& species = read_cache(particle_->species());
    size_type nparticle = particle_->nparticle();

    LOG_TRACE("compute coefficients for tabulated force");

    scoped_timer_type timer(runtime_.compute);

    multi_range_for_each(
        grid_size_type(0)
      , grid_size_
      , [&](grid_size_type const& index) {
          size_type c_offset = coefficients_per_knot * multi_index_to_offset(index, grid_size_);
          position_type r1 = element_prod(static_cast<position_type>(index), grid_basis);

          // use a fixed species, always assume that the test particle
          // is the first species, ie. species 0
          species_type const a = 0;
          float_type en_pot(0);
          vector_type first_der(0);
          vector_type sec_der(0);
          float_type third_der(0);

          for (size_type i = 0; i < nparticle; ++i) {
              species_type const b = species[i];
              // particle distance vector
              // in contrast to the interpolation box, the simulation box does not
              // start at zero-origin, so shift the position of the particle accordingly
              position_type r = r1 - (position[i] - box_->origin());
              box_->reduce_periodic(r);

              float_type rr = inner_prod(r, r);
              // truncate potential at cutoff length
              if (rr >= potential_->rr_cut(a, b))
                  continue;

              float_type fval, pot;
              boost::tie(fval, pot) = (*potential_)(rr, a, b);
              en_pot += pot;
              // The potential functor returns the force, but we need the first
              // derivative of the potential, therefore apply minus sign
              first_der -= r * fval;

              if (dimension > 1) {
                  float_type second_derivative, third_derivative;

                  boost::tie(second_derivative, third_derivative) = potential_->derivatives(rr, a, b);
                  // d²/dxdy
                  sec_der[0] += second_derivative * r[0] * r[1];
                  if (dimension > 2) {
                    // d²/dxdz
                    sec_der[1] += second_derivative * r[0] * r[2];
                    // d²/dydz
                    sec_der[2] += second_derivative * r[1] * r[2];
                    third_der += third_derivative * r[0] * r[1] * r[2];
                  }
              }
          }
          // set the coefficients as needed
          coefficients_[c_offset + 0] = en_pot;
          coefficients_[c_offset + 1] = first_der[0]; // d/dx
          if(dimension > 1) {
              coefficients_[c_offset + 2] = first_der[1]; // d/dy
              coefficients_[c_offset + 3] = sec_der[0];   // d²/dxdy
          }
          if (dimension > 2) {
            coefficients_[c_offset + 4] = first_der[2]; // d/dz
            coefficients_[c_offset + 5] = sec_der[1];   // d²/dxdz
            coefficients_[c_offset + 6] = sec_der[2];   // d²/dydz
            coefficients_[c_offset + 7] = third_der;    // d³/dxdydz
          }
      }
    );

    coefficient_cache_ = current_state;
}

/**
 * Copy coefficients array to given array.
 */
template <typename force_type, typename iterator_type>
inline iterator_type
get_coefficients(force_type& force, iterator_type const& first)
{
    force.compute();
    return std::copy(force.coefficients().begin(), force.coefficients().end(), first);
}

template <typename force_type>
static std::function<std::vector<typename force_type::coefficient_value_type> ()>
wrap_get_coefficients(std::shared_ptr<force_type> self)
{
    return [=]() -> std::vector<typename force_type::coefficient_value_type> {
        std::vector<typename force_type::coefficient_value_type> output;
        {
            output.reserve(self->coefficients().size());
        }
        get_coefficients(*self, std::back_inserter(output));
        return std::move(output);
    };
}

template <int dimension, typename float_type, typename potential_type>
void tabulated_generator<dimension, float_type, potential_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("forces")
            [
                class_<tabulated_generator>()
                    .property("get_coefficients", &wrap_get_coefficients<tabulated_generator>)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                    ]
                    .def_readonly("runtime", &tabulated_generator::runtime_)

              , def("tabulated_generator", &std::make_shared<tabulated_generator,
                    std::shared_ptr<potential_type const>
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<box_type const>
                  , grid_size_type const&
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}


} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* HALMD_MDSIM_HOST_FORCES_TABLULATED_GENERATOR_HPP */
