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

#ifndef HALMD_MDSIM_GPU_NEIGHBOUR_HPP
#define HALMD_MDSIM_GPU_NEIGHBOUR_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/binning.hpp>
#include <halmd/mdsim/gpu/neighbour_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/neighbour.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

template <int dimension, typename float_type>
class neighbour
  : public mdsim::neighbour<dimension>
{
public:
    typedef mdsim::neighbour<dimension> _Base;
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef mdsim::box<dimension> box_type;
    typedef gpu::binning<dimension, float_type> binning_type;
    typedef utility::profiler profiler_type;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type update;
    };

    static void luaopen(lua_State* L);

    neighbour(
        boost::shared_ptr<particle_type const> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<binning_type const> binning
      , matrix_type const& r_cut
      , double skin
      , double cell_occupancy
    );
    void register_runtimes(profiler_type& profiler);
    virtual void update();

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
    cuda::vector<unsigned int> const& g_neighbour() const
    {
       return g_neighbour_;
    }

    /**
     * number of placeholders per neighbour list
     */
    unsigned int size() const
    {
        return size_;
    }

    /**
     * neighbour list stride
     */
    unsigned int stride() const
    {
        return stride_;
    }

private:
    boost::shared_ptr<particle_type const> particle_;
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<binning_type const> binning_;

    /** neighbour list skin in MD units */
    float_type r_skin_;
    /** half neighbour list skin */
    float_type rr_skin_half_;
    /** (cutoff lengths + neighbour list skin)² */
    matrix_type rr_cut_skin_;
    /** (cutoff lengths + neighbour list skin)² */
    cuda::vector<float_type> g_rr_cut_skin_;
    /** FIXME average desired cell occupancy */
    float_type nu_cell_;
    /** particle positions at last neighbour list update */
    cuda::vector<float4> g_r0_;
    /** block-reduced squared particle distances */
    cuda::vector<float> g_rr_;
    /** block-reduced squared particle distances */
    cuda::host::vector<float> h_rr_;

    /** neighbour lists */
    cuda::vector<unsigned int> g_neighbour_;
    /** number of placeholders per neighbour list */
    unsigned int size_;
    /** neighbour list stride */
    unsigned int stride_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOUR_HPP */
