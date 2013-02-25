/*
 * Copyright © 2012-2013 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2012      Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_PARTICLE_GROUP_HPP
#define HALMD_MDSIM_HOST_PARTICLE_GROUP_HPP

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/raw_array.hpp>

#include <lua.hpp>

#include <algorithm>
#include <tuple>

namespace halmd {
namespace mdsim {
namespace host {

/**
 * A particle group represents a subset of particles, which is defined
 * by an instance of particle together with a sequence of indices.
 */
class particle_group
{
public:
    typedef raw_array<unsigned int> array_type;
    typedef typename array_type::value_type size_type;

    /**
     * Returns ordered sequence of particle indices.
     */
    virtual cache<array_type> const& ordered() = 0;

    /**
     * Returns unordered sequence of particle indices.
     */
    virtual cache<array_type> const& unordered() = 0;

    /**
     * Returns number of particles.
     */
    virtual cache<size_type> const& size() = 0;

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);
};

/**
 * Copy ordered sequence of particle indices to array in host memory.
 */
template <typename iterator_type>
inline iterator_type
get_ordered(particle_group& group, iterator_type const& first)
{
    auto const& ordered = read_cache(group.ordered());
    return std::copy(ordered.begin(), ordered.end(), first);
}

/**
 * Copy unordered sequence of particle indices to array in host memory.
 */
template <typename iterator_type>
inline iterator_type
get_unordered(particle_group& group, iterator_type const& first)
{
    auto const& unordered = read_cache(group.unordered());
    return std::copy(unordered.begin(), unordered.end(), first);
}

/**
 * Compute mean kinetic energy per particle.
 */
template <typename particle_type>
double get_mean_en_kin(particle_type const& particle, particle_group& group)
{
    typename particle_group::array_type const& unordered        = *group.unordered();
    auto const& velocity = read_cache(particle.velocity());
    typename particle_type::mass_array_type const& mass         = *particle.mass();

    double mv2 = 0;
    for (typename particle_group::size_type i : unordered) {
        mv2 += mass[i] * inner_prod(velocity[i], velocity[i]);
    }
    return  0.5 * mv2 / unordered.size();
}

/**
 * Compute centre of mass.
 */
template <typename particle_type, typename box_type>
fixed_vector<double, particle_type::velocity_type::static_size>
get_r_cm(particle_type const& particle, particle_group& group, box_type const& box)
{
    auto const& unordered = read_cache(group.unordered());
    auto const& position = read_cache(particle.position());
    auto const& image = read_cache(particle.image());
    auto const& mass = read_cache(particle.mass());

    fixed_vector<double, particle_type::velocity_type::static_size> mr = 0;
    double m = 0;
    for (typename particle_group::size_type i : unordered) {
        auto r = position[i];
        box.extend_periodic(r, image[i]);
        mr += mass[i] * r;
        m += mass[i];
    }
    return mr / m;
}

/**
 * Compute velocity of centre of mass, and mean mass
 */
template <typename particle_type>
std::tuple<fixed_vector<double, particle_type::velocity_type::static_size>, double>
get_v_cm_and_mean_mass(particle_type const& particle, particle_group& group)
{
    auto const& unordered = read_cache(group.unordered());
    auto const& velocity = read_cache(particle.velocity());
    auto const& mass = read_cache(particle.mass());

    fixed_vector<double, particle_type::velocity_type::static_size> mv = 0;
    double m = 0;
    for (typename particle_group::size_type i : unordered) {
        mv += mass[i] * velocity[i];
        m += mass[i];
    }
    return std::make_tuple(mv / m, m / unordered.size());
}

/**
 * Compute velocity of centre of mass.
 */
template <typename particle_type>
fixed_vector<double, particle_type::velocity_type::static_size>
get_v_cm(particle_type const& particle, particle_group& group)
{
    return std::get<0>(get_v_cm_and_mean_mass(particle, group));
}

/**
 * Compute mean potential energy per particle.
 */
template <typename force_type>
double get_mean_en_pot(force_type& force, particle_group& group)
{
    auto const& unordered = read_cache(group.unordered());
    auto const& en_pot = read_cache(force.en_pot());

    double sum = 0;
    for (typename particle_group::size_type i : unordered) {
        sum += en_pot[i];
    }
    return sum / unordered.size();
}

/**
 * Compute mean virial per particle.
 */
template <typename force_type>
double get_mean_virial(force_type& force, particle_group& group)
{
    enum { dimension = force_type::net_force_type::static_size };

    auto const& unordered = read_cache(group.unordered());
    auto const& stress_pot = read_cache(force.stress_pot());

    double sum = 0;
    for (typename particle_group::size_type i : unordered) {
        // compute trace of the stress tensor
        for (int j = 0; j < dimension; ++j) {
            sum += stress_pot[i][j];
        }
    }
    return sum / unordered.size();
}

/**
 * Compute mean hypervirial per particle.
 */
template <typename force_type>
double get_mean_hypervirial(force_type& force, particle_group& group)
{
    auto const& unordered = read_cache(group.unordered());
    auto const& hypervirial = read_cache(force.hypervirial());

    double sum = 0;
    for (typename particle_group::size_type i : unordered) {
        sum += hypervirial[i];
    }
    return sum / unordered.size();
}

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_GROUP_HPP */
