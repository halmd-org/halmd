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

#ifndef HALMD_MDSIM_HOST_PARTICLE_HPP
#define HALMD_MDSIM_HOST_PARTICLE_HPP

#include <lua.hpp>
#include <vector>

#include <halmd/mdsim/particle.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/raw_allocator.hpp>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class particle
  : public mdsim::particle<dimension>
{
    typedef mdsim::particle<dimension> _Base;

public:
    typedef typename type_traits<dimension, float_type>::vector_type vector_type;

    typedef vector_type position_type;
    typedef vector_type image_type;
    typedef vector_type velocity_type;
    typedef unsigned int tag_type;
    typedef unsigned int reverse_tag_type;
    typedef unsigned int species_type;
    typedef double mass_type;
    typedef vector_type force_type;
    typedef double en_pot_type;
    typedef typename type_traits<dimension, float_type>::stress_tensor_type stress_pot_type;
    typedef double hypervirial_type;

    typedef std::vector<position_type> position_array_type;
    typedef std::vector<image_type> image_array_type;
    typedef std::vector<velocity_type> velocity_array_type;
    typedef std::vector<tag_type> tag_array_type;
    typedef std::vector<reverse_tag_type> reverse_tag_array_type;
    typedef std::vector<species_type> species_array_type;
    typedef std::vector<mass_type> mass_array_type;
    typedef std::vector<force_type, raw_allocator<force_type> > force_array_type;
    typedef std::vector<en_pot_type, raw_allocator<en_pot_type> > en_pot_array_type;
    typedef std::vector<stress_pot_type, raw_allocator<stress_pot_type> > stress_pot_array_type;
    typedef std::vector<hypervirial_type, raw_allocator<hypervirial_type> > hypervirial_array_type;

    static void luaopen(lua_State* L);

    particle(
        std::vector<unsigned int> const& particles
      , std::vector<double> const& mass
    );
    virtual void set();
    void rearrange(std::vector<unsigned int> const& index);

    /** positions, reduced to extended domain box */
    std::vector<vector_type> r;
    /** minimum image vectors */
    std::vector<vector_type> image;
    /** velocities */
    std::vector<vector_type> v;
    /** globally unique particle numbers */
    std::vector<unsigned int> tag;
    /** reverse particle tags */
    std::vector<unsigned int> reverse_tag;
    /** types */
    std::vector<unsigned int> type;

private:
    /** number of particles in simulation box */
    using _Base::nbox;
    /** number of particle types */
    using _Base::ntype;

public:
    /** number of particles per type */
    using _Base::ntypes;
    /** mass per type */
    using _Base::mass;

    /**
     * Enable computation of auxiliary variables.
     *
     * The flag is reset by the next call to prepare().
     */
    void aux_enable();

    /**
     * Reset forces, and optionally auxiliary variables, to zero.
     */
    void prepare();

    /**
     * Returns number of particles.
     */
    std::size_t nparticle() const
    {
        return nbox;
    }

    /**
     * Returns number of species.
     */
    unsigned int nspecies() const
    {
        return ntype;
    }

    /**
     * Returns non-const reference to force per particle.
     */
    force_array_type const& force() const
    {
       return force_;
    }

    /**
     * Returns const reference to force per particle.
     */
    force_array_type& force()
    {
       return force_;
    }

    /**
     * Returns const reference to potential energy per particle.
     *
     * This method checks that the computation of auxiliary variables was enabled.
     */
    en_pot_array_type const& en_pot() const
    {
        assert_aux_valid();
        return en_pot_;
    }

    /**
     * Returns non-const reference to potential energy per particle.
     */
    en_pot_array_type& en_pot()
    {
        return en_pot_;
    }

    /**
     * Returns const reference to potential part of stress tensor per particle.
     *
     * This method checks that the computation of auxiliary variables was enabled.
     */
    stress_pot_array_type const& stress_pot() const
    {
        assert_aux_valid();
        return stress_pot_;
    }

    /**
     * Returns non-const reference to potential part of stress tensor per particle.
     */
    stress_pot_array_type& stress_pot()
    {
        return stress_pot_;
    }

    /**
     * Returns const reference to hypervirial per particle.
     *
     * This method checks that the computation of auxiliary variables was enabled.
     */
    hypervirial_array_type const& hypervirial() const
    {
        assert_aux_valid();
        return hypervirial_;
    }

    /**
     * Returns non-const reference to hypervirial per particle.
     */
    hypervirial_array_type& hypervirial()
    {
        return hypervirial_;
    }

    /**
     * Returns true if computation of auxiliary variables is enabled.
     */
    bool aux_valid() const
    {
        return aux_valid_;
    }

private:
    /** force per particle */
    force_array_type force_;
    /** potential energy per particle */
    en_pot_array_type en_pot_;
    /** potential part of stress tensor per particle */
    stress_pot_array_type stress_pot_;
    /** hypervirial per particle */
    hypervirial_array_type hypervirial_;

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

} // namespace mdsim
} // namespace host
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_HPP */
