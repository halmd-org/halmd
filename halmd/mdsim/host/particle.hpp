/*
 * Copyright © 2016      Daniel Kirchner
 * Copyright © 2010-2012 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_PARTICLE_HPP
#define HALMD_MDSIM_HOST_PARTICLE_HPP

#include <halmd/mdsim/host/particle_array.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/raw_array.hpp>

#include <lua.hpp>

#include <unordered_map>
#include <algorithm>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type_>
class particle
{
public:
    typedef float_type_ float_type;
    typedef halmd::signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;

    typedef fixed_vector<float_type, dimension> vector_type;

    typedef unsigned int size_type;
    typedef vector_type position_type;
    typedef vector_type image_type;
    typedef vector_type orientation_type;
    typedef vector_type velocity_type;
    typedef unsigned int id_type;
    typedef unsigned int reverse_id_type;
    typedef unsigned int species_type;
    typedef float_type mass_type;
    typedef fixed_vector<float_type, dimension> force_type;
    typedef typename type_traits<dimension, float_type>::pseudo_vector_type torque_type;
    typedef float_type en_pot_type;
    typedef typename type_traits<dimension, float_type>::stress_tensor_type stress_pot_type;

    typedef raw_array<position_type> position_array_type;
    typedef raw_array<orientation_type> orientation_array_type;
    typedef raw_array<image_type> image_array_type;
    typedef raw_array<velocity_type> velocity_array_type;
    typedef raw_array<id_type> id_array_type;
    typedef raw_array<reverse_id_type> reverse_id_array_type;
    typedef raw_array<species_type> species_array_type;
    typedef raw_array<mass_type> mass_array_type;
    typedef raw_array<force_type> force_array_type;
    typedef raw_array<torque_type> torque_array_type;
    typedef raw_array<en_pot_type> en_pot_array_type;
    typedef raw_array<stress_pot_type> stress_pot_array_type;

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
     * register typed particle data
     *
     * @param name identifier for the particle data
     * @param update_function optional update function
     * @return the newly created particle array
     *
     * throws an exception if a particle array with the same name already exists
     */
    template<typename T>
    std::shared_ptr<particle_array_typed<T>>
    register_data(std::string const& name, std::function<void()> update_function = std::function<void()>())
    {
        auto ptr = particle_array::create<T>(nparticle_, capacity_, update_function);
        if (!data_.insert(std::make_pair(name, ptr)).second) {
            throw std::runtime_error("a particle array named \"" + name + "\" already exists");
        }
        return ptr;
    }

    /**
     * get data from named particle array with iterator
     *
     * @param name identifier of the particle array
     * @param first output iterator
     * @return output iterator
     *
     * throws an exception if the array does not exist or has an invalid type
     */
    template<typename T, typename iterator_type>
    iterator_type get_data(std::string const& name, iterator_type const& first) const
    {
        return particle_array::cast<T>(get_array(name))->get_data(first);
    }

    /**
     * set data in named particle array with iterator
     *
     * @param name identifier of the particle array
     * @param first input iterator
     * @return input iterator
     *
     * throws an exception if the array does not exist or has an invalid type
     */
    template <typename T, typename iterator_type>
    iterator_type set_data(std::string const& name, iterator_type const& first)
    {
        return particle_array::cast<T>(get_array(name))->set_data(first);
    }

    /**
     * Returns const reference to the data of a named particle array.
     *
     * @param name identifier of the particle array
     * @return const reference to the data
     *
     * throws an exception if the array does not exist or has an invalid type
     */
    template<typename T>
    cache<raw_array<T>> const& data(std::string const& name) const
    {
        return particle_array::cast<T>(get_array(name))->data();
    }

    /**
     * Returns non-const reference to the data of a named particle array.
     *
     * @param name identifier of the particle array
     * @return non-const reference to the data
     *
     * throws an exception if the array does not exist or has an invalid type
     */
    template<typename T>
    cache<raw_array<T>>& mutable_data(std::string const& name)
    {
        return particle_array::cast<T>(get_array(name))->mutable_data();
    }

    /**
     * Returns const reference to particle positions.
     */
    cache<position_array_type> const& position() const
    {
        return data<position_type>("position");
    }

    /**
     * Returns non-const reference to particle positions.
     */
    cache<position_array_type>& position()
    {
        return mutable_data<position_type>("position");
    }

    /**
     * Returns non-const reference to particle images.
     */
    cache<image_array_type> const& image() const
    {
        return data<image_type>("image");
    }

    /**
     * Returns const reference to particle images.
     */
    cache<image_array_type>& image()
    {
        return mutable_data<image_type>("image");
    }

    /**
     * Returns const reference to particle orientations.
     */
    cache<orientation_array_type> const& orientation() const
    {
        return data<orientation_type>("orientation");
    }

    /**
     * Returns non-const reference to particle orientations.
     */
    cache<orientation_array_type>& orientation()
    {
        return mutable_data<orientation_type>("orientation");
    }

    /**
     * Returns const reference to particle velocities.
     */
    cache<velocity_array_type> const& velocity() const
    {
        return data<velocity_type>("velocity");
    }

    /**
     * Returns non-const reference to particle velocities.
     */
    cache<velocity_array_type>& velocity()
    {
        return mutable_data<velocity_type>("velocity");
    }

    /**
     * Returns const reference to particle IDs.
     */
    cache<id_array_type> const& id() const
    {
        return data<id_type>("id");
    }

    /**
     * Returns non-const reference to particle IDs.
     */
    cache<id_array_type>& id()
    {
        return mutable_data<id_type>("id");
    }

    /**
     * Returns const reference to particle reverse IDs.
     */
    cache<reverse_id_array_type> const& reverse_id() const
    {
        return data<reverse_id_type>("reverse_id");
    }

    /**
     * Returns non-const reference to particle reverse IDs.
     */
    cache<reverse_id_array_type>& reverse_id()
    {
        return mutable_data<reverse_id_type>("reverse_id");
    }

    /**
     * Returns const reference to particle species.
     */
    cache<species_array_type> const& species() const
    {
        return data<species_type>("species");
    }

    /**
     * Returns non-const reference to particle species.
     */
    cache<species_array_type>& species()
    {
        return mutable_data<species_type>("species");
    }

    /**
     * Returns const reference to particle masses.
     */
    cache<mass_array_type> const& mass() const
    {
        return data<mass_type>("mass");
    }

    /**
     * Returns non-const reference to particle masses.
     */
    cache<mass_array_type>& mass()
    {
        return mutable_data<mass_type>("mass");
    }
    /**
     * Returns const reference to particle force.
     */
    cache<force_array_type> const& force()
    {
        return data<force_type>("force");
    }

    /**
     * Returns non-const reference to particle force.
     */
    cache<force_array_type>& mutable_force()
    {
        return mutable_data<force_type>("force");
    }

    /**
     * Returns const reference to particle torque.
     */
    cache<torque_array_type> const& torque()
    {
        return data<torque_type>("torque");
    }

    /**
     * Returns non-const reference to particle torque.
     */
    cache<torque_array_type>& mutable_torque()
    {
        return mutable_data<torque_type>("torque");
    }

    /**
     * Returns const reference to potential energy of particles.
     */
    cache<en_pot_array_type> const& potential_energy()
    {
        return data<en_pot_type>("potential_energy");
    }

    /**
     * Returns non-const reference to potential energy of particles.
     */
    cache<en_pot_array_type>& mutable_potential_energy()
    {
        return mutable_data<en_pot_type>("potential_energy");
    }

    /**
     * Returns const reference to potential part of stress tensor.
     */
    cache<stress_pot_array_type> const& stress_pot()
    {
        return data<stress_pot_type>("potential_stress_tensor");
    }

    /**
     * Returns non-const reference to potential part of stress tensor.
     */
    cache<stress_pot_array_type>& mutable_stress_pot()
    {
        return mutable_data<stress_pot_type>("potential_stress_tensor");
    }

    /**
     * Enable computation of auxiliary variables.
     *
     * The flag is reset after the next trigger of on_force_().
     */
    void aux_enable();

    /**
     * Returns true if computation of auxiliary variables is enabled.
     */
    bool aux_enabled() const
    {
        return aux_enabled_;
    }

    /**
     * Returns true if the force has to be reset to zero prior to reading.
     */
    bool force_zero()
    {
        return force_zero_;
    }

    /**
     * Disable a reset of the force to zero upon reading.
     *
     * Must be called after computation of the first force contribution.
     */
    void force_zero_disable()
    {
        force_zero_ = false;
    }

    /**
     * Indicate that a force update (ie. triggering on_force_()) is required
     */
    void mark_force_dirty()
    {
        force_dirty_ = true;
    }

    /**
     * Indicate that an auxiliary update (ie. triggering on_force_() with
     * aux_enabled) is required
     */
    void mark_aux_dirty()
    {
        aux_dirty_ = true;
    }

    connection on_prepend_force(slot_function_type const& slot)
    {
        return on_prepend_force_.connect(slot);
    }

    connection on_force(slot_function_type const& slot)
    {
        return on_force_.connect(slot);
    }

    connection on_append_force(slot_function_type const& slot)
    {
        return on_append_force_.connect(slot);
    }

    std::shared_ptr<particle_array> const& get_array(std::string const& name) const
    {
        auto it = data_.find(name);
        if(it == data_.end()) {
            throw std::invalid_argument("particle array \"" + name + "\" not registered");
        }
        return it->second;
    }

    bool has_array(std::string const& name) const
    {
        return data_.find(name) != data_.end();
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** number of particles */
    unsigned int nparticle_;
    /** array size */
    unsigned int capacity_;
    /** number of particle species */
    unsigned int nspecies_;

    /** map of the stored particle arrays */
    std::unordered_map<std::string, std::shared_ptr<particle_array>> data_;

    /** flag that the force has to be reset to zero prior to reading */
    bool force_zero_;
    /** flag that the force cache is dirty (not up to date) */
    bool force_dirty_;
    /** flag that the caches of the auxiliary variables are dirty (not up to date) */
    bool aux_dirty_;
    /** flag that the computation of auxiliary variables is requested */
    bool aux_enabled_;

    /**
     * Update all forces and auxiliary variables if needed. The auxiliary
     * variables are guaranteed to be up-to-date upon return if with_aux was
     * set to true.
     *
     * Auxiliary variables are computed only if they are out of date
     * (aux_dirty_ == true) and if either with_aux or aux_enabled_ is true.
     *
     * Emit a warning if the force update would be necessary solely to compute
     * the auxiliary variables, which indicates a performance problem.
     */
    void update_force_(bool with_aux=false);

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type rearrange;
    };

    /** profiling runtime accumulators */
    runtime runtime_;

    signal_type on_prepend_force_;
    signal_type on_force_;
    signal_type on_append_force_;
};

/**
 * Copy particle positions to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_position(particle_type const& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::position_type>("position", first);
}

/**
 * Copy particle positions from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_position(particle_type& particle, iterator_type const& first)
{
    return particle.template set_data<typename particle_type::position_type>("position", first);
}

/**
 * Copy particle images to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_image(particle_type const& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::image_type>("image", first);
}

/**
 * Copy particle images from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_image(particle_type& particle, iterator_type const& first)
{
    return particle.template set_data<typename particle_type::image_type>("image", first);
}

/**
 * Copy particle orientations to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_orientation(particle_type const& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::orientation_type>("orientation", first);
}

/**
 * Copy particle orientations from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_orientation(particle_type& particle, iterator_type const& first)
{
    return particle.template set_data<typename particle_type::orientation_type>("orientation", first);
}

/**
 * Copy particle velocities to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_velocity(particle_type const& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::velocity_type>("velocity", first);
}

/**
 * Copy particle velocities from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_velocity(particle_type& particle, iterator_type const& first)
{
    return particle.template set_data<typename particle_type::velocity_type>("velocity", first);
}

/**
 * Copy particle IDs to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_id(particle_type const& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::id_type>("id", first);
}

/**
 * Copy particle IDs from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_id(particle_type& particle, iterator_type const& first)
{
    return particle.template set_data<typename particle_type::id_type>("id", first);
}

/**
 * Copy particle reverse IDs to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_reverse_id(particle_type const& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::reverse_id_type>("reverse_id", first);
}

/**
 * Copy particle reverse IDs from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_reverse_id(particle_type& particle, iterator_type const& first)
{
    return particle.template set_data<typename particle_type::reverse_id_type>("reverse_id", first);
}

/**
 * Copy particle species to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_species(particle_type const& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::species_type>("species", first);
}

/**
 * Copy particle species from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_species(particle_type& particle, iterator_type const& first)
{
    return particle.template set_data<typename particle_type::species_type>("species", first);
}

/**
 * Copy particle masses to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_mass(particle_type const& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::mass_type>("mass", first);
}

/**
 * Copy particle masses from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_mass(particle_type& particle, iterator_type const& first)
{
    return particle.template set_data<typename particle_type::mass_type>("mass", first);
}

/**
 * Copy force per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_force(particle_type& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::force_type>("force", first);
}

/**
 * Copy torque per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_torque(particle_type& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::torque_type>("torque", first);
}

/**
 * Copy potential energy per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_potential_energy(particle_type& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::en_pot_type>("potential_energy", first);
}

/**
 * Copy potential part of stress tensor per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_stress_pot(particle_type& particle, iterator_type const& first)
{
    return particle.template get_data<typename particle_type::stress_pot_type>("potential_stress_tensor", first);
}

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_HPP */
