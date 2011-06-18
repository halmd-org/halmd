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

#ifndef HALMD_MDSIM_GPU_NEIGHBOUR_WITH_BINNING_HPP
#define HALMD_MDSIM_GPU_NEIGHBOUR_WITH_BINNING_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/binning.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/neighbour_with_binning_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension, typename float_type>
class neighbour_with_binning
  : public gpu::neighbour
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef mdsim::box<dimension> box_type;
    typedef gpu::binning<dimension, float_type> binning_type;
    typedef utility::profiler profiler_type;

    static void luaopen(lua_State* L);

    neighbour_with_binning(
        boost::shared_ptr<particle_type const> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<binning_type const> binning
      , matrix_type const& r_cut
      , double skin
      , double cell_occupancy
    );
    void register_runtimes(profiler_type& profiler);
    void update();

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
    typedef gpu::neighbour _Base;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type update;
    };

    boost::shared_ptr<particle_type const> particle_;
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<binning_type const> binning_;

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
};

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOUR_WITH_BINNING_HPP */
