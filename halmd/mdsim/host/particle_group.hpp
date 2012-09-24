/*
 * Copyright © 2012 Peter Colberg
 * Copyright © 2012 Felix Höfling
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
    typedef typename particle_group::array_type array_type;
    cache_proxy<array_type const> ordered = group.ordered();
    return std::copy(ordered->begin(), ordered->end(), first);
}

/**
 * Copy unordered sequence of particle indices to array in host memory.
 */
template <typename iterator_type>
inline iterator_type
get_unordered(particle_group& group, iterator_type const& first)
{
    typedef typename particle_group::array_type array_type;
    cache_proxy<array_type const> unordered = group.unordered();
    return std::copy(unordered->begin(), unordered->end(), first);
}

/**
 * Compute mean kinetic energy per particle.
 */
template <typename particle_type>
double get_mean_en_kin(particle_type const& particle, particle_group& group)
{
    typedef typename particle_group::array_type group_array_type;
    typedef typename particle_group::size_type size_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_type::mass_array_type mass_array_type;

    cache_proxy<group_array_type const> unordered = group.unordered();
    cache_proxy<velocity_array_type const> velocity = particle.velocity();
    cache_proxy<mass_array_type const> mass = particle.mass();

    double mv2 = 0;
    for (size_type i : *unordered) {
        mv2 += (*mass)[i] * inner_prod((*velocity)[i], (*velocity)[i]);
    }
    return  0.5 * mv2 / unordered->size();
}

/**
 * Compute velocity of centre of mass, and mean mass.
 */
template <typename particle_type>
std::tuple<fixed_vector<double, particle_type::velocity_type::static_size>, double>
get_v_cm_and_mean_mass(particle_type const& particle, particle_group& group)
{
    typedef typename particle_group::array_type group_array_type;
    typedef typename particle_group::size_type size_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_type::mass_array_type mass_array_type;

    cache_proxy<group_array_type const> unordered = group.unordered();
    cache_proxy<velocity_array_type const> velocity = particle.velocity();
    cache_proxy<mass_array_type const> mass = particle.mass();

    fixed_vector<double, particle_type::velocity_type::static_size> mv = 0;
    double m = 0;
    for (size_type i : *unordered) {
        mv += (*mass)[i] * (*velocity)[i];
        m += (*mass)[i];
    }
    return std::make_tuple(mv / m, m / unordered->size());
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
    typedef typename particle_group::size_type size_type;
    typedef typename particle_group::array_type group_array_type;
    typedef typename force_type::en_pot_array_type en_pot_array_type;

    cache_proxy<group_array_type const> unordered = group.unordered();
    cache_proxy<en_pot_array_type const> en_pot = force.en_pot();

    double sum = 0;
    for (size_type i : *unordered) {
        sum += (*en_pot)[i];
    }
    return sum / unordered->size();
}

/**
 * Compute mean virial per particle.
 */
template <typename force_type>
double get_mean_virial(force_type& force, particle_group& group)
{
    typedef typename particle_group::size_type size_type;
    typedef typename particle_group::array_type group_array_type;
    typedef typename force_type::stress_pot_array_type stress_pot_array_type;

    cache_proxy<group_array_type const> unordered = group.unordered();
    cache_proxy<stress_pot_array_type const> stress_pot = force.stress_pot();

    double sum = 0;
    for (size_type i : *unordered) {
        sum += (*stress_pot)[i][0];
    }
    return sum / unordered->size();
}

/**
 * Compute mean hypervirial per particle.
 */
template <typename force_type>
double get_mean_hypervirial(force_type& force, particle_group& group)
{
    typedef typename particle_group::size_type size_type;
    typedef typename particle_group::array_type group_array_type;
    typedef typename force_type::hypervirial_array_type hypervirial_array_type;

    cache_proxy<group_array_type const> unordered = group.unordered();
    cache_proxy<hypervirial_array_type const> hypervirial = force.hypervirial();

    double sum = 0;
    for (size_type i : *unordered) {
        sum += (*hypervirial)[i];
    }
    return sum / unordered->size();
}

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_GROUP_HPP */
