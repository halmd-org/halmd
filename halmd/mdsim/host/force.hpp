/*
 * Copyright © 2012 Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_FORCE_HPP
#define HALMD_MDSIM_HOST_FORCE_HPP

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/raw_array.hpp>

#include <lua.hpp>

#include <algorithm>
#include <memory>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class force
{
public:
    typedef fixed_vector<float_type, dimension> net_force_type;
    typedef double en_pot_type;
    typedef typename type_traits<dimension, float_type>::stress_tensor_type stress_pot_type;
    typedef double hypervirial_type;

    typedef raw_array<net_force_type> net_force_array_type;
    typedef raw_array<en_pot_type> en_pot_array_type;
    typedef raw_array<stress_pot_type> stress_pot_array_type;
    typedef raw_array<hypervirial_type> hypervirial_array_type;

    virtual ~force() {}

    /**
     * Returns const reference to net force per particle.
     */
    virtual cache<net_force_array_type> const& net_force() = 0;

    /**
     * Returns const reference to potential energy per particle.
     */
    virtual cache<en_pot_array_type> const& en_pot() = 0;

    /**
     * Returns const reference to potential part of stress tensor of each particle.
     */
    virtual cache<stress_pot_array_type> const& stress_pot() = 0;

    /**
     * Returns const reference to hypervirial per particle.
     */
    virtual cache<hypervirial_array_type> const& hypervirial() = 0;

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);
};

/**
 * Copy net force per particle to given array.
 */
template <typename force_type, typename iterator_type>
inline iterator_type
get_net_force(force_type& force, iterator_type const& first)
{
    typedef typename force_type::net_force_array_type net_force_array_type;
    cache_proxy<net_force_array_type const> net_force = force.net_force();
    return std::copy(net_force->begin(), net_force->end(), first);
}

/**
 * Copy potential energy per particle to given array.
 */
template <typename force_type, typename iterator_type>
inline iterator_type
get_en_pot(force_type& force, iterator_type const& first)
{
    typedef typename force_type::en_pot_array_type en_pot_array_type;
    cache_proxy<en_pot_array_type const> en_pot = force.en_pot();
    return std::copy(en_pot->begin(), en_pot->end(), first);
}

/**
 * Copy potential part of stress tensor per particle to given array.
 */
template <typename force_type, typename iterator_type>
inline iterator_type
get_stress_pot(force_type& force, iterator_type const& first)
{
    typedef typename force_type::stress_pot_array_type stress_pot_array_type;
    cache_proxy<stress_pot_array_type const> stress_pot = force.stress_pot();
    return std::copy(stress_pot->begin(), stress_pot->end(), first);
}

/**
 * Copy hypervirial per particle to given array.
 */
template <typename force_type, typename iterator_type>
inline iterator_type
get_hypervirial(force_type& force, iterator_type const& first)
{
    typedef typename force_type::hypervirial_array_type hypervirial_array_type;
    cache_proxy<hypervirial_array_type const> hypervirial = force.hypervirial();
    return std::copy(hypervirial->begin(), hypervirial->end(), first);
}

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCE_HPP */
