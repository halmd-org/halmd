/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP

#include <lua.hpp>
#include <limits>
#include <stdexcept> // std::logic_error
#include <vector>

#include <halmd/mdsim/clock.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/raw_array.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace samples {

template <int dimension, typename float_type>
class phase_space
{
public:
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef raw_array<vector_type> position_array_type;
    typedef raw_array<vector_type> velocity_array_type;
    typedef raw_array<unsigned int> species_array_type;
    typedef typename mdsim::clock::step_type step_type;

    /**
     * Construct phase space sample.
     *
     * @param nparticle number of particles
     * @param step simulation step the sample is taken (optional)
     */
    phase_space(std::size_t nparticle, step_type step = std::numeric_limits<step_type>::max());

    /**
     * Returns const reference to particle positions.
     *
     * The positions are extended with their periodic image vectors.
     */
    position_array_type const& position() const
    {
        return position_;
    }

    /**
     * Returns non-const reference to particle positions.
     *
     * The positions are extended with their periodic image vectors.
     */
    position_array_type& position()
    {
        return position_;
    }

    /**
     * Returns const reference to particle velocities.
     */
    velocity_array_type const& velocity() const
    {
        return velocity_;
    }

    /**
     * Returns non-const reference to particle velocities.
     */
    velocity_array_type& velocity()
    {
        return velocity_;
    }

    /**
     * Returns const reference to particle species.
     */
    species_array_type const& species() const
    {
        return species_;
    }

    /**
     * Returns non-const reference to particle species.
     */
    species_array_type& species()
    {
        return species_;
    }

    /**
     * Returns simulation step when the sample was taken.
     */
    step_type step() const
    {
        if (step_ == std::numeric_limits<step_type>::max()) {
            throw std::logic_error("step not set in phase space sample");
        }
        return step_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** periodically extended particle positions */
    position_array_type position_;
    /** particle velocities */
    velocity_array_type velocity_;
    /** particle species */
    species_array_type species_;
    /** simulation step when sample was taken */
    step_type step_;
};

template <int dimension, typename float_type>
inline phase_space<dimension, float_type>::phase_space(std::size_t nparticle, step_type step)
  : position_(nparticle)
  , velocity_(nparticle)
  , species_(nparticle)
  , step_(step)
{
}

} // namespace samples
} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP */
