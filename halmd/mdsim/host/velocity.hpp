/*
 * Copyright © 2010, 2012 Peter Colberg
 * Copyright © 2010 Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_VELOCITY_HPP
#define HALMD_MDSIM_HOST_VELOCITY_HPP

#include <halmd/mdsim/host/particle_group.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace mdsim {
namespace host {

/**
 * Rescale magnitude of all velocities by 'factor'
 */
template <typename particle_type>
inline void rescale_velocity(particle_type& particle, double factor)
{
    auto velocity = make_cache_mutable(particle.velocity());
    for (unsigned int i = 0; i < particle.nparticle(); i++) {
        (*velocity)[i] *= factor;
    }
}

/**
 * Rescale magnitude of velocities of group by 'factor'
 */
template <typename particle_type>
inline void rescale_velocity_group(particle_type& particle, particle_group& group, double factor)
{
    auto velocity = make_cache_mutable(particle.velocity());
    for (typename particle_group::size_type i : *group.unordered()) {
        (*velocity)[i] *= factor;
    }
}

/**
 * Shift all velocities by 'delta'
 */
template <typename particle_type>
inline void shift_velocity(particle_type& particle, fixed_vector<double, particle_type::velocity_type::static_size> const& delta)
{
    auto velocity = make_cache_mutable(particle.velocity());
    for (unsigned int i = 0; i < particle.nparticle(); i++) {
        (*velocity)[i] += delta;
    }
}

/**
 * Shift velocities of group by 'delta'
 */
template <typename particle_type>
inline void shift_velocity_group(particle_type& particle, particle_group& group, fixed_vector<double, particle_type::velocity_type::static_size> const& delta)
{
    auto velocity = make_cache_mutable(particle.velocity());
    for (typename particle_group::size_type i : *group.unordered()) {
        (*velocity)[i] += delta;
    }
}

/**
 * First shift, then rescale all velocities
 */
template <typename particle_type>
inline void shift_rescale_velocity(particle_type& particle, fixed_vector<double, particle_type::velocity_type::static_size> const& delta, double factor)
{
    auto velocity = make_cache_mutable(particle.velocity());
    for (unsigned int i = 0; i < particle.nparticle(); i++) {
        (*velocity)[i] += delta;
        (*velocity)[i] *= factor;
    }
}

/**
 * First shift, then rescale velocities of group
 */
template <typename particle_type>
inline void shift_rescale_velocity_group(particle_type& particle, particle_group& group, fixed_vector<double, particle_type::velocity_type::static_size> const& delta, double factor)
{
    auto velocity = make_cache_mutable(particle.velocity());
    for (typename particle_group::size_type i : *group.unordered()) {
        auto& v = (*velocity)[i];
        v += delta;
        v *= factor;
    }
}

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_VELOCITY_HPP */
