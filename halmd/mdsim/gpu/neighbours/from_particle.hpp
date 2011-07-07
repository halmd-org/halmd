/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_NEIGHBOURS_FROM_PARTICLE_HPP
#define HALMD_MDSIM_GPU_NEIGHBOURS_FROM_PARTICLE_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {

/**
 * Neighbour lists without binning
 *
 * This module builds neighbour lists by direct iteration over all
 * particles, which scales as O(N^2) with the number of particles N.
 * The GPU kernel uses shared memory to share particle coordinates between
 * threads of a block, as described in the section “N-squared MD” in
 *
 * J. A. van Meel, A. Arnold, D. Frenkel, S.F. Portegies Zwart and R. G.
 * Belleman, Harvesting graphics power for MD simulations, Molecular
 * Simulation, 34 (3) 259-266 (2008)
 * http://dx.doi.org/10.1080/08927020701744295
 * http://arxiv.org/abs/0709.3225
 */
template <int dimension, typename float_type>
class from_particle
  : public gpu::neighbour
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef mdsim::box<dimension> box_type;
    typedef utility::profiler profiler_type;
    struct defaults;
    typedef typename from_particle::signal_type signal_type; // import type from base class
    typedef typename from_particle::slot_function_type slot_function_type; // import type from base class

    static void luaopen(lua_State* L);

    from_particle(
        boost::shared_ptr<particle_type const> particle1
      , boost::shared_ptr<particle_type const> particle2
      , boost::shared_ptr<box_type const> box
      , matrix_type const& r_cut
      , double skin
      , double cell_occupancy = defaults::occupancy()
    );
    void register_runtimes(profiler_type& profiler);
    virtual void update();

    virtual void on_update(slot_function_type const& slot)
    {
        on_update_.connect(slot);
    }

    //! returns neighbour list skin in MD units
    float_type r_skin() const
    {
        return r_skin_;
    }

    //! returns neighbour list skin in MD units
    float_type cell_occupancy() const
    {
        return nu_cell_;
    }

    /**
     * neighbour lists
     */
    virtual cuda::vector<unsigned int> const& g_neighbour() const
    {
       return g_neighbour_;
    }

    /**
     * number of placeholders per neighbour list
     */
    virtual unsigned int size() const
    {
        return size_;
    }

    /**
     * neighbour list stride
     */
    virtual unsigned int stride() const
    {
        return stride_;
    }

private:
    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type update;
    };

    boost::shared_ptr<particle_type const> particle1_;
    boost::shared_ptr<particle_type const> particle2_;
    boost::shared_ptr<box_type const> box_;

    /** neighbour list skin in MD units */
    float_type r_skin_;
    /** (cutoff lengths + neighbour list skin)² */
    matrix_type rr_cut_skin_;
    /** (cutoff lengths + neighbour list skin)² */
    cuda::vector<float_type> g_rr_cut_skin_;
    /** FIXME average desired cell occupancy */
    float_type nu_cell_;
    /** neighbour lists */
    cuda::vector<unsigned int> g_neighbour_;
    /** number of placeholders per neighbour list */
    unsigned int size_;
    /** neighbour list stride */
    unsigned int stride_;

    /** profiling runtime accumulators */
    runtime runtime_;
    /** signal emitted before neighbour list update */
    signal_type on_update_;
};

template <int dimension, typename float_type>
struct from_particle<dimension, float_type>::defaults
{
    static float_type occupancy();
};

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOURS_FROM_PARTICLE_HPP */
