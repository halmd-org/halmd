/*
 * Copyright © 2012-2016 Felix Höfling
 * Copyright © 2016      Daniel Kirchner
 * Copyright © 2012      Peter Colberg
 * Copyright © 2013-2015 Nicolas Höft
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

#ifndef HALMD_MDSIM_GPU_PARTICLE_GROUP_HPP
#define HALMD_MDSIM_GPU_PARTICLE_GROUP_HPP

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/gpu/thermodynamics_kernel.hpp>
#include <halmd/mdsim/gpu/particle_group_kernel.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <algorithm>
#include <tuple>
#include <stdexcept>

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * A particle group represents a subset of particles, which is defined
 * by an instance of particle together with a sequence of indices.
 */
class particle_group
{
public:
    typedef cuda::vector<unsigned int> array_type;
    typedef array_type::value_type size_type;
    typedef cuda::host::vector<size_type> host_array_type;

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
     * Returns ordered sequence of particle indices in host memory.
     * If the stored cached copy of the indices is no longer valid,
     * it is automatically updated from the GPU storage.
     */
    host_array_type const& ordered_host_cached()
    {
        if (!(ordered_observer_ == ordered())) {
            auto const& indices = read_cache(ordered());
            cached_ordered_.resize(indices.size());
            cuda::copy(indices.begin(), indices.end(), cached_ordered_.begin());
            ordered_observer_ = ordered();
        }
        return cached_ordered_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    host_array_type cached_ordered_;
    cache<> ordered_observer_;
};

/**
 * Copy ordered sequence of particle indices to array in host memory.
 */
template <typename iterator_type>
inline iterator_type
get_ordered(particle_group& group, iterator_type const& first)
{
    typedef particle_group::array_type array_type;
    array_type const& g_ordered = read_cache(group.ordered());
    cuda::host::vector<array_type::value_type> h_ordered(g_ordered.size());
    cuda::copy(g_ordered.begin(), g_ordered.end(), h_ordered.begin());
    return std::copy(h_ordered.begin(), h_ordered.end(), first);
}

/**
 * Copy unordered sequence of particle indices to array in host memory.
 */
template <typename iterator_type>
inline iterator_type
get_unordered(particle_group& group, iterator_type const& first)
{
    typedef particle_group::array_type array_type;
    array_type const& g_unordered = read_cache(group.unordered());
    cuda::host::vector<array_type::value_type> h_unordered(g_unordered.size());
    cuda::copy(g_unordered.begin(), g_unordered.end(), h_unordered.begin());
    return std::copy(h_unordered.begin(), h_unordered.end(), first);
}

/**
 * Compute mean kinetic energy per particle.
 */
template <typename particle_type>
double get_mean_en_kin(particle_type const& particle, particle_group& group)
{
    typedef particle_group::array_type group_array_type;
    unsigned int constexpr dimension = particle_type::velocity_type::static_size;
    typedef observables::gpu::kinetic_energy<dimension, dsfloat> accumulator_type;

    group_array_type const& unordered = read_cache(group.unordered());
    accumulator_type::get().bind(*particle.velocity());
    return reduce(&*unordered.begin(), &*unordered.end(), accumulator_type())() / unordered.size();
}

/**
 * Compute centre of mass.
 */
template <typename particle_type, typename box_type>
fixed_vector<double, particle_type::velocity_type::static_size>
get_r_cm(particle_type const& particle, particle_group& group, box_type const& box)
{
    typedef particle_group::array_type group_array_type;
    typedef typename particle_type::position_type position_type;
    unsigned int constexpr dimension = particle_type::position_type::static_size;
    typedef observables::gpu::centre_of_mass<dimension, dsfloat> accumulator_type;

    group_array_type const& unordered = read_cache(group.unordered());

    accumulator_type::get_position().bind(*particle.position());
    accumulator_type::get_image().bind(*particle.image());
    accumulator_type::get_velocity().bind(*particle.velocity());
    return reduce(
        std::make_tuple(&*unordered.begin(), static_cast<position_type>(box.length()))
      , std::make_tuple(&*unordered.end())
      , accumulator_type()
   )();
}

/**
 * Compute velocity of centre of mass.
 */
template <typename particle_type>
std::tuple<fixed_vector<double, particle_type::velocity_type::static_size>, double>
get_v_cm_and_mean_mass(particle_type const& particle, particle_group& group)
{
    typedef particle_group::array_type group_array_type;
    unsigned int constexpr dimension = particle_type::velocity_type::static_size;
    typedef observables::gpu::velocity_of_centre_of_mass<dimension, dsfloat> accumulator_type;

    group_array_type const& unordered = read_cache(group.unordered());

    accumulator_type::get().bind(*particle.velocity());
    fixed_vector<double, dimension> mv;
    double m;
    std::tie(mv, m) = reduce(&*unordered.begin(), &*unordered.end(), accumulator_type())();
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
 * Compute total force on all particles.
 */
template <typename particle_type>
fixed_vector<double, particle_type::force_type::static_size>
get_total_force(particle_type& particle, particle_group& group)
{
    typedef particle_group::array_type group_array_type;
    unsigned int constexpr dimension = particle_type::force_type::static_size;
    typedef observables::gpu::total_force<dimension, dsfloat> accumulator_type;

    group_array_type const& unordered = read_cache(group.unordered());

    accumulator_type::get().bind(*particle.force());
    return reduce(&*unordered.begin(), &*unordered.end(), accumulator_type())();
}

/**
 * Compute mean potential energy per particle.
 */
template <typename particle_type>
double get_mean_en_pot(particle_type& particle, particle_group& group)
{
    typedef particle_group::array_type group_array_type;
    typedef observables::gpu::potential_energy<dsfloat> accumulator_type;

    group_array_type const& unordered = read_cache(group.unordered());

    accumulator_type::get().bind(*particle.potential_energy());
    return reduce(&*unordered.begin(), &*unordered.end(), accumulator_type())() / unordered.size();
}

/**
 * Compute mean virial per particle.
 */
template <typename particle_type>
double get_mean_virial(particle_type& particle, particle_group& group)
{
    typedef particle_group::array_type group_array_type;
    typedef typename particle_type::stress_pot_type stress_pot_type;
    unsigned int constexpr dimension = particle_type::force_type::static_size;
    typedef observables::gpu::virial<dimension, dsfloat> accumulator_type;

    group_array_type const& unordered = read_cache(group.unordered());

    unsigned int stride = particle.stress_pot()->capacity() / stress_pot_type::static_size;
    accumulator_type::get().bind(*particle.stress_pot());
    return reduce(&*unordered.begin(), &*unordered.end(), accumulator_type(stride))() / unordered.size();
}

/**
 * Compute stress tensor.
 */
template <typename particle_type>
typename type_traits<particle_type::force_type::static_size, double>::stress_tensor_type
get_stress_tensor(particle_type& particle, particle_group& group)
{
    typedef particle_group::array_type group_array_type;
    typedef typename particle_type::stress_pot_array_type stress_pot_array_type;
    typedef typename particle_type::stress_pot_type stress_pot_type;

    enum { dimension = particle_type::force_type::static_size };

    typedef typename type_traits<dimension, double>::stress_tensor_type stress_tensor_type;
    typedef observables::gpu::stress_tensor<dimension, dsfloat> accumulator_type;

    group_array_type const& unordered = read_cache(group.unordered());
    stress_pot_array_type const& stress_pot = read_cache(particle.stress_pot());

    unsigned int stride = stress_pot.capacity() / stress_pot_type::static_size;
    accumulator_type::get_stress_pot().bind(stress_pot);
    accumulator_type::get_velocity().bind(*particle.velocity());

    return stress_tensor_type(reduce(&*unordered.begin(), &*unordered.end(), accumulator_type(stride))());
}

/**
 * Copy all particles from a group into a given particle instance of the
 * same size as the group
 */
template <typename particle_type>
void particle_group_to_particle(particle_type const& particle_src, particle_group& group, particle_type& particle_dst)
{
    typedef typename particle_type::float_type float_type;
    enum { dimension = particle_type::force_type::static_size };

    if(*group.size() != particle_dst.nparticle()) {
        LOG_TRACE("group size: " << *group.size() << ", destination particle size: " << particle_dst.nparticle());
        throw std::logic_error("source group size and destination particle size must match!");
    }

    auto const& ordered = read_cache(group.ordered());
    auto position = make_cache_mutable(particle_dst.position());
    auto image    = make_cache_mutable(particle_dst.image());
    auto velocity = make_cache_mutable(particle_dst.velocity());

    particle_group_wrapper<dimension, float_type>::kernel.r.bind(read_cache(particle_src.position()));
    particle_group_wrapper<dimension, float_type>::kernel.image.bind(read_cache(particle_src.image()));
    particle_group_wrapper<dimension, float_type>::kernel.v.bind(read_cache(particle_src.velocity()));

    configure_kernel(particle_group_wrapper<dimension, float_type>::kernel.particle_group_to_particle, ordered.size());
    particle_group_wrapper<dimension, float_type>::kernel.particle_group_to_particle(
        &*ordered.begin()
      , position->data()
      , &*image->begin()
      , velocity->data()
      , ordered.size()
    );
}

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_GROUP_HPP */
