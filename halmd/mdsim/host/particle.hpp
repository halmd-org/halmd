/*
 * Copyright © 2008-2012 Peter Colberg
 * Copyright © 2010-2012 Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_PARTICLE_HPP
#define HALMD_MDSIM_HOST_PARTICLE_HPP

#include <halmd/utility/profiler.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/raw_array.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

#include <lua.hpp>

#include <algorithm>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class particle
{
public:
    typedef fixed_vector<float_type, dimension> vector_type;

    typedef unsigned int size_type;
    typedef vector_type position_type;
    typedef vector_type image_type;
    typedef vector_type velocity_type;
    typedef unsigned int tag_type;
    typedef unsigned int reverse_tag_type;
    typedef unsigned int species_type;
    typedef double mass_type;

    typedef raw_array<position_type> position_array_type;
    typedef raw_array<image_type> image_array_type;
    typedef raw_array<velocity_type> velocity_array_type;
    typedef raw_array<tag_type> tag_array_type;
    typedef raw_array<reverse_tag_type> reverse_tag_array_type;
    typedef raw_array<species_type> species_array_type;
    typedef raw_array<mass_type> mass_array_type;

    void rearrange(std::vector<unsigned int> const& index);

    /**
     * Allocate particle arrays in host memory.
     *
     * @param nparticle number of particles
     * @param nspecies number of particle species
     *
     * All particle arrays, except the masses, are initialised to zero.
     * The particle masses are initialised to unit mass.
     */
    particle(size_type nparticle, unsigned int nspecies);

    /**
     * Returns number of particles.
     *
     * Currently the number of particles is fixed at construction of
     * particle. This may change in the future, to allow for chemical
     * reactions that do not conserve the number of particles, or to
     * transfer particles between domains of different processors.
     */
    size_type nparticle() const
    {
        return nparticle_;
    }

    /**
     * Returns number of species.
     */
    unsigned int nspecies() const
    {
        return nspecies_;
    }

    /**
     * Returns const reference to particle positions.
     */
    cache<position_array_type> const& position() const
    {
        return position_;
    }

    /**
     * Returns non-const reference to particle positions.
     */
    cache<position_array_type>& position()
    {
        return position_;
    }

    /**
     * Returns non-const reference to particle images.
     */
    cache<image_array_type> const& image() const
    {
        return image_;
    }

    /**
     * Returns const reference to particle images.
     */
    cache<image_array_type>& image()
    {
        return image_;
    }

    /**
     * Returns const reference to particle velocities.
     */
    cache<velocity_array_type> const& velocity() const
    {
        return velocity_;
    }

    /**
     * Returns non-const reference to particle velocities.
     */
    cache<velocity_array_type>& velocity()
    {
        return velocity_;
    }

    /**
     * Returns const reference to particle tags.
     */
    cache<tag_array_type> const& tag() const
    {
        return tag_;
    }

    /**
     * Returns non-const reference to particle tags.
     */
    cache<tag_array_type>& tag()
    {
        return tag_;
    }

    /**
     * Returns const reference to particle reverse tags.
     */
    cache<reverse_tag_array_type> const& reverse_tag() const
    {
        return reverse_tag_;
    }

    /**
     * Returns non-const reference to particle reverse tags.
     */
    cache<reverse_tag_array_type>& reverse_tag()
    {
        return reverse_tag_;
    }

    /**
     * Returns const reference to particle species.
     */
    cache<species_array_type> const& species() const
    {
        return species_;
    }

    /**
     * Returns non-const reference to particle species.
     */
    cache<species_array_type>& species()
    {
        return species_;
    }

    /**
     * Returns const reference to particle masses.
     */
    cache<mass_array_type> const& mass() const
    {
        return mass_;
    }

    /**
     * Returns non-const reference to particle masses.
     */
    cache<mass_array_type>& mass()
    {
        return mass_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** number of particles */
    unsigned int nparticle_;
    /** number of particle species */
    unsigned int nspecies_;
    /** positions, reduced to extended domain box */
    cache<position_array_type> position_;
    /** minimum image vectors */
    cache<image_array_type> image_;
    /** velocities */
    cache<velocity_array_type> velocity_;
    /** particle tags */
    cache<tag_array_type> tag_;
    /** reverse particle tags */
    cache<reverse_tag_array_type> reverse_tag_;
    /** particle species */
    cache<species_array_type> species_;
    /** particle masses */
    cache<mass_array_type> mass_;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type rearrange;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

/**
 * Copy particle positions to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_position(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::position_array_type position_array_type;
    cache_proxy<position_array_type const> position = particle.position();
    return std::copy(position->begin(), position->end(), first);
}

/**
 * Copy particle positions from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_position(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::position_array_type position_array_type;
    cache_proxy<position_array_type> position = particle.position();
    iterator_type input = first;
    for (typename particle_type::position_type& value : *position) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle images to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_image(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::image_array_type image_array_type;
    cache_proxy<image_array_type const> image = particle.image();
    return std::copy(image->begin(), image->end(), first);
}

/**
 * Copy particle images from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_image(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::image_array_type image_array_type;
    cache_proxy<image_array_type> image = particle.image();
    iterator_type input = first;
    for (typename particle_type::image_type& value : *image) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle velocities to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_velocity(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::velocity_array_type velocity_array_type;
    cache_proxy<velocity_array_type const> velocity = particle.velocity();
    return std::copy(velocity->begin(), velocity->end(), first);
}

/**
 * Copy particle velocities from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_velocity(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::velocity_array_type velocity_array_type;
    cache_proxy<velocity_array_type> velocity = particle.velocity();
    iterator_type input = first;
    for (typename particle_type::velocity_type& value : *velocity) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle tags to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_tag(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::tag_array_type tag_array_type;
    cache_proxy<tag_array_type const> tag = particle.tag();
    return std::copy(tag->begin(), tag->end(), first);
}

/**
 * Copy particle tags from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_tag(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::tag_array_type tag_array_type;
    cache_proxy<tag_array_type> tag = particle.tag();
    iterator_type input = first;
    for (typename particle_type::tag_type& value : *tag) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle reverse tags to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_reverse_tag(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::reverse_tag_array_type reverse_tag_array_type;
    cache_proxy<reverse_tag_array_type const> reverse_tag = particle.reverse_tag();
    return std::copy(reverse_tag->begin(), reverse_tag->end(), first);
}

/**
 * Copy particle reverse tags from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_reverse_tag(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::reverse_tag_array_type reverse_tag_array_type;
    cache_proxy<reverse_tag_array_type> reverse_tag = particle.reverse_tag();
    iterator_type input = first;
    for (typename particle_type::reverse_tag_type& value : *reverse_tag) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle species to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_species(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::species_array_type species_array_type;
    cache_proxy<species_array_type const> species = particle.species();
    return std::copy(species->begin(), species->end(), first);
}

/**
 * Copy particle species from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_species(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::species_array_type species_array_type;
    cache_proxy<species_array_type> species = particle.species();
    iterator_type input = first;
    for (typename particle_type::species_type& value : *species) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle masses to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_mass(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::mass_array_type mass_array_type;
    cache_proxy<mass_array_type const> mass = particle.mass();
    return std::copy(mass->begin(), mass->end(), first);
}

/**
 * Copy particle masses from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_mass(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::mass_array_type mass_array_type;
    cache_proxy<mass_array_type> mass = particle.mass();
    iterator_type input = first;
    for (typename particle_type::mass_type& value : *mass) {
        value = *input++;
    }
    return input;
}


} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_HPP */
