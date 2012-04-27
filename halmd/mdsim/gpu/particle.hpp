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

#ifndef HALMD_MDSIM_GPU_PARTICLE_HPP
#define HALMD_MDSIM_GPU_PARTICLE_HPP

#include <lua.hpp>
#include <vector>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension, typename float_type>
class particle
{
public:
    typedef typename type_traits<dimension, float_type>::vector_type vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type gpu_vector_type;
    struct defaults;

    typedef vector_type position_type;
    typedef vector_type image_type;
    typedef vector_type velocity_type;
    typedef unsigned int tag_type;
    typedef unsigned int reverse_tag_type;
    typedef unsigned int species_type;
    typedef float mass_type;
    typedef vector_type force_type;
    typedef float_type en_pot_type;
    typedef typename type_traits<dimension, float_type>::stress_tensor_type stress_pot_type;
    typedef float_type hypervirial_type;

    typedef cuda::vector<float4> position_array_type;
    typedef cuda::vector<gpu_vector_type> image_array_type;
    typedef cuda::vector<float4> velocity_array_type;
    typedef cuda::vector<tag_type> tag_array_type;
    typedef cuda::vector<reverse_tag_type> reverse_tag_array_type;
    typedef cuda::vector<gpu_vector_type> force_array_type;
    typedef cuda::vector<en_pot_type> en_pot_array_type;
    typedef cuda::vector<typename type_traits<dimension, float_type>::gpu::stress_tensor_type> stress_pot_array_type;
    typedef cuda::vector<hypervirial_type> hypervirial_array_type;

    void set();
    void rearrange(cuda::vector<unsigned int> const& g_index);

    /** grid and block dimensions for CUDA calls */
    cuda::config const dim;

    /**
     * Allocate particle arrays in GPU memory.
     */
    particle(std::size_t nparticle, unsigned int threads = defaults::threads());

    /**
     * Returns number of particles.
     *
     * Currently the number of particles is fixed at construction of
     * particle. This may change in the future, to allow for chemical
     * reactions that do not conserve the number of particles, or to
     * transfer particles between domains of different processors.
     */
    std::size_t nparticle() const
    {
        return g_tag_.size();
    }

    /**
     * Returns number of particle placeholders.
     *
     * Currently the number of placeholders, i.e. the element count of the
     * particle arrays in memory, is equal to the total number of kernel
     * threads, which is a multiple of the number of threads per block,
     * and greater or equal than the number of particles.
     */
    std::size_t nplaceholder() const
    {
        return g_tag_.capacity();
    }

    /**
     * Returns number of species.
     */
    unsigned int nspecies() const
    {
        return nspecies_;
    }

    /**
     * Returns non-const reference to particle positions and species.
     */
    position_array_type const& position() const
    {
        return g_position_;
    }

    /**
     * Returns const reference to particle positions and species.
     */
    position_array_type& position()
    {
        return g_position_;
    }

    /**
     * Copy particle positions to given array.
     */
    void get_position(std::vector<position_type>& position);

    /**
     * Copy particle positions from given array.
     */
    void set_position(std::vector<position_type> const& position);

    /**
     * Copy particle species to given array.
     */
    void get_species(std::vector<species_type>& species);

    /**
     * Copy particle species from given array.
     */
    void set_species(std::vector<species_type> const& species);

    /**
     * Returns non-const reference to particle images.
     */
    image_array_type const& image() const
    {
        return g_image_;
    }

    /**
     * Returns const reference to particle images.
     */
    image_array_type& image()
    {
        return g_image_;
    }

    /**
     * Copy particle images to given array.
     */
    void get_image(std::vector<image_type>& image);

    /**
     * Copy particle images from given array.
     */
    void set_image(std::vector<image_type> const& image);

    /**
     * Returns non-const reference to particle velocities and masses.
     */
    velocity_array_type const& velocity() const
    {
        return g_velocity_;
    }

    /**
     * Returns const reference to particle velocities and masses.
     */
    velocity_array_type& velocity()
    {
        return g_velocity_;
    }

    /**
     * Copy particle velocities to given array.
     */
    void get_velocity(std::vector<velocity_type>& velocity);

    /**
     * Copy particle velocities from given array.
     */
    void set_velocity(std::vector<velocity_type> const& velocity);

    /**
     * Copy particle masses to given array.
     */
    void get_mass(std::vector<mass_type>& mass);

    /**
     * Copy particle masses from given array.
     */
    void set_mass(std::vector<mass_type> const& mass);

    /**
     * Set particle masses to scalar.
     *
     * This includes the masses of particle placeholders.
     */
    void set_mass(float_type mass);

    /**
     * Returns non-const reference to particle tags.
     */
    tag_array_type const& tag() const
    {
        return g_tag_;
    }

    /**
     * Returns const reference to particle tags.
     */
    tag_array_type& tag()
    {
        return g_tag_;
    }

    /**
     * Copy particle tags to given array.
     */
    void get_tag(std::vector<tag_type>& tag);

    /**
     * Copy particle tags from given array.
     */
    void set_tag(std::vector<tag_type> const& tag);

    /**
     * Returns non-const reference to particle reverse tags.
     */
    reverse_tag_array_type const& reverse_tag() const
    {
        return g_reverse_tag_;
    }

    /**
     * Returns const reference to particle reverse tags.
     */
    reverse_tag_array_type& reverse_tag()
    {
        return g_reverse_tag_;
    }

    /**
     * Copy particle reverse tags to given array.
     */
    void get_reverse_tag(std::vector<reverse_tag_type>& reverse_tag);

    /**
     * Copy particle reverse tags from given array.
     */
    void set_reverse_tag(std::vector<reverse_tag_type> const& reverse_tag);

    /**
     * Returns non-const reference to force per particle.
     */
    force_array_type const& force() const
    {
        return g_force_;
    }

    /**
     * Returns const reference to force per particle.
     */
    force_array_type& force()
    {
        return g_force_;
    }

    /**
     * Copy force per particle to given array.
     */
    void get_force(std::vector<force_type>& force);

    /**
     * Copy force per particle from given array.
     */
    void set_force(std::vector<force_type> const& force);

    /**
     * Returns const reference to potential energy per particle.
     *
     * This method checks that the computation of auxiliary variables was enabled.
     */
    en_pot_array_type const& en_pot() const
    {
        assert_aux_valid();
        return g_en_pot_;
    }

    /**
     * Returns non-const reference to potential energy per particle.
     */
    en_pot_array_type& en_pot()
    {
        return g_en_pot_;
    }

    /**
     * Copy potential energy per particle to given array.
     */
    void get_en_pot(std::vector<en_pot_type>& en_pot);

    /**
     * Copy potential energy per particle from given array.
     */
    void set_en_pot(std::vector<en_pot_type> const& en_pot);

    /**
     * Returns const reference to potential part of stress tensor per particle.
     *
     * This method checks that the computation of auxiliary variables was enabled.
     */
    stress_pot_array_type const& stress_pot() const
    {
        assert_aux_valid();
        return g_stress_pot_;
    }

    /**
     * Returns non-const reference to potential part of stress tensor per particle.
     */
    stress_pot_array_type& stress_pot()
    {
        return g_stress_pot_;
    }

    /**
     * Copy potential part of stress tensor per particle to given array.
     */
    void get_stress_pot(std::vector<stress_pot_type>& stress_pot);

    /**
     * Copy potential part of stress tensor per particle from given array.
     */
    void set_stress_pot(std::vector<stress_pot_type> const& stress_pot);

    /**
     * Returns const reference to hypervirial per particle.
     *
     * This method checks that the computation of auxiliary variables was enabled.
     */
    hypervirial_array_type const& hypervirial() const
    {
        assert_aux_valid();
        return g_hypervirial_;
    }

    /**
     * Returns non-const reference to hypervirial per particle.
     */
    hypervirial_array_type& hypervirial()
    {
        return g_hypervirial_;
    }

    /**
     * Copy hypervirial per particle to given array.
     */
    void get_hypervirial(std::vector<hypervirial_type>& hypervirial);

    /**
     * Copy hypervirial per particle from given array.
     */
    void set_hypervirial(std::vector<hypervirial_type> const& hypervirial);

    /**
     * Enable computation of auxiliary variables.
     *
     * The flag is reset by the next call to prepare().
     */
    void aux_enable();

    /**
     * Returns true if computation of auxiliary variables is enabled.
     */
    bool aux_valid() const
    {
        return aux_valid_;
    }

    /**
     * Reset forces, and optionally auxiliary variables, to zero.
     */
    void prepare();

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** number of particle species */
    unsigned int nspecies_;
    /** positions, species */
    position_array_type g_position_;
    /** minimum image vectors */
    image_array_type g_image_;
    /** velocities, masses */
    velocity_array_type g_velocity_;
    /** particle tags */
    tag_array_type g_tag_;
    /** reverse particle tags */
    reverse_tag_array_type g_reverse_tag_;
    /** force per particle */
    force_array_type g_force_;
    /** potential energy per particle */
    en_pot_array_type g_en_pot_;
    /** potential part of stress tensor per particle */
    stress_pot_array_type g_stress_pot_;
    /** hypervirial per particle */
    hypervirial_array_type g_hypervirial_;

    /** flag for enabling the computation of auxiliary variables this step */
    bool aux_flag_;
    /** flag that indicates the auxiliary variables are computed this step */
    bool aux_valid_;

    void assert_aux_valid() const
    {
        if (!aux_valid_) {
            throw std::logic_error("auxiliary variables were not enabled in particle");
        }
    }

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

template <int dimension, typename float_type>
struct particle<dimension, float_type>::defaults
{
    static unsigned int threads();
};

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_HPP */
