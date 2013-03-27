/*
 * Copyright © 2013 Felix Höfling
 * Copyright © 2013 Nicolas Höft
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

#ifndef HALMD_MDSIM_GPU_FORCE_HPP
#define HALMD_MDSIM_GPU_FORCE_HPP

#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/cache.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <algorithm>
#include <memory>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension, typename float_type>
class force
{
public:
    typedef fixed_vector<float_type, dimension> net_force_type;
    typedef float_type en_pot_type;
    typedef typename type_traits<dimension, float_type>::stress_tensor_type stress_pot_type;
    typedef float_type hypervirial_type;

    typedef cuda::vector<typename type_traits<dimension, float_type>::gpu::coalesced_vector_type> net_force_array_type;
    typedef cuda::vector<en_pot_type> en_pot_array_type;
    typedef cuda::vector<typename stress_pot_type::value_type> stress_pot_array_type;
    typedef cuda::vector<hypervirial_type> hypervirial_array_type;

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
     * Returns const reference to potential part of stress tensor per particle.
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
    cache_proxy<net_force_array_type const> g_net_force = force.net_force();
    cuda::host::vector<typename net_force_array_type::value_type> h_net_force(g_net_force->size());
    cuda::copy(g_net_force->begin(), g_net_force->end(), h_net_force.begin());
    return std::copy(h_net_force.begin(), h_net_force.end(), first);
}

/**
 * Copy potential energy per particle to given array.
 */
template <typename force_type, typename iterator_type>
inline iterator_type
get_en_pot(force_type& force, iterator_type const& first)
{
    typedef typename force_type::en_pot_array_type en_pot_array_type;
    cache_proxy<en_pot_array_type const> g_en_pot = force.en_pot();
    cuda::host::vector<typename en_pot_array_type::value_type> h_en_pot(g_en_pot->size());
    cuda::copy(g_en_pot->begin(), g_en_pot->end(), h_en_pot.begin());
    return std::copy(h_en_pot.begin(), h_en_pot.end(), first);
}

/**
 * Copy potential part of stress tensor per particle to given array.
 */
template <typename force_type, typename iterator_type>
inline iterator_type
get_stress_pot(force_type& force, iterator_type const& first)
{
    // copy data from GPU to host
    typedef typename force_type::stress_pot_array_type stress_pot_array_type;
    cache_proxy<stress_pot_array_type const> g_stress_pot = force.stress_pot();
    cuda::host::vector<typename stress_pot_array_type::value_type> h_stress_pot(g_stress_pot->size());
    h_stress_pot.reserve(g_stress_pot->capacity());
    cuda::copy(g_stress_pot->begin(), g_stress_pot->begin() + g_stress_pot->capacity(), h_stress_pot.begin());

    // convert from column-major to row-major layout
    typedef typename force_type::stress_pot_type stress_pot_type;
    unsigned int stride = h_stress_pot.capacity() / stress_pot_type::static_size;
    iterator_type output = first;
    for (auto const& stress : h_stress_pot) {
        *output++ = read_stress_tensor<stress_pot_type>(&stress, stride);
    }
    return output;
}

/**
 * Copy hypervirial per particle to given array.
 */
template <typename force_type, typename iterator_type>
inline iterator_type
get_hypervirial(force_type& force, iterator_type const& first)
{
    typedef typename force_type::hypervirial_array_type hypervirial_array_type;
    cache_proxy<hypervirial_array_type const> g_hypervirial = force.hypervirial();
    cuda::host::vector<typename hypervirial_array_type::value_type> h_hypervirial(g_hypervirial->size());
    cuda::copy(g_hypervirial->begin(), g_hypervirial->end(), h_hypervirial.begin());
    return std::copy(h_hypervirial.begin(), h_hypervirial.end(), first);
}

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCE_HPP */
