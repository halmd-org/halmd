/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP

#include <lua.hpp>
#include <limits>
#include <stdexcept> // std::logic_error
#include <vector>

#include <halmd/mdsim/clock.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/host/samples/sample.hpp>
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
    typedef raw_array<float_type> mass_array_type;
    typedef typename mdsim::clock::step_type step_type;

    typedef sample<dimension, float_type> position_sample_type;
    typedef sample<dimension, float_type> velocity_sample_type;
    typedef sample<1, unsigned int> species_sample_type;
    typedef sample<1, float_type> mass_sample_type;

    /**
     * Construct phase space sample.
     *
     * @param nparticle number of particles
     * @param step simulation step the sample is taken (optional)
     */
    phase_space(std::size_t nparticle, step_type step = std::numeric_limits<step_type>::max());

    phase_space(std::shared_ptr<sample<dimension, float_type> const> position
             ,  std::shared_ptr<sample<dimension, float_type> const> velocity
             ,  std::shared_ptr<sample<1, unsigned int> const> species
             ,  std::shared_ptr<sample<1, float_type> const> mass
             ,  step_type step = std::numeric_limits<step_type>::max())
            : position_(position), velocity_(velocity), species_(species), mass_(mass), step_ (step)
    {
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

    std::shared_ptr<position_sample_type const> position_sample() const {
        return position_;
    };

    void set_position_sample(std::shared_ptr<position_sample_type const> sample) {
        position_ = sample;
    };

    std::shared_ptr<velocity_sample_type const> velocity_sample() const {
        return velocity_;
    };

    void set_velocity_sample(std::shared_ptr<velocity_sample_type const> sample) {
        velocity_ = sample;
    };

    std::shared_ptr<species_sample_type const> species_sample() const {
        return species_;
    };

    void set_species_sample(std::shared_ptr<species_sample_type const> sample) {
        species_ = sample;
    };

    std::shared_ptr<mass_sample_type const> mass_sample() const {
        return mass_;
    };

    void set_mass_sample(std::shared_ptr<mass_sample_type const> sample) {
        mass_ = sample;
    };


    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** periodically extended particle positions */
    std::shared_ptr<position_sample_type const> position_;
    /** particle velocities */
    std::shared_ptr<velocity_sample_type const> velocity_;
    /** particle species */
    std::shared_ptr<species_sample_type const> species_;
    /** particle mass */
    std::shared_ptr<mass_sample_type const> mass_;
    /** simulation step when sample was taken */
    step_type step_;
};

template <int dimension, typename float_type>
inline phase_space<dimension, float_type>::phase_space(std::size_t nparticle, step_type step)
  : position_(new sample<dimension, float_type>(nparticle))
  , velocity_(new sample<dimension, float_type>(nparticle))
  , species_(new sample<1, unsigned int>(nparticle))
  , mass_(new sample<1, float_type>(nparticle))
  , step_(step)
{
}

} // namespace samples
} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP */
