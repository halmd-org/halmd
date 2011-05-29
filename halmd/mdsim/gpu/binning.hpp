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

#ifndef HALMD_MDSIM_GPU_BINNING_HPP
#define HALMD_MDSIM_GPU_BINNING_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/binning_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

template <int dimension, typename float_type>
class binning
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef mdsim::box<dimension> box_type;
    typedef utility::profiler profiler_type;

    typedef fixed_vector<unsigned int, dimension> cell_size_type;
    typedef fixed_vector<int, dimension> cell_diff_type;

    static void luaopen(lua_State* L);

    binning(
        boost::shared_ptr<particle_type const> particle
      , boost::shared_ptr<box_type const> box
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

    //! returns average desired cell occupancy
    float_type cell_occupancy() const
    {
        return nu_cell_;
    }

    //! returns average effective cell occupancy
    float_type effective_cell_occupancy() const
    {
        return nu_cell_eff_;
    }

    //! returns cell size
    unsigned int cell_size() const
    {
        return cell_size_;
    }

    /**
     * number of cells per dimension
     */
    cell_size_type ncell() const
    {
        return ncell_;
    }

    /**
     * cell kernel dimensions
     */
    cuda::config const& dim_cell() const
    {
        return dim_cell_;
    }

    /**
     * cell lists
     */
    cuda::vector<unsigned int> const& g_cell() const
    {
        return g_cell_;
    }

private:
    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type update;
    };

    boost::shared_ptr<particle_type const> particle_;
    boost::shared_ptr<box_type const> box_;

    /** neighbour list skin in MD units */
    float_type r_skin_;
    /** half neighbour list skin */
    float_type rr_skin_half_;
    /** (cutoff lengths + neighbour list skin)² */
    matrix_type rr_cut_skin_;
    /** average desired cell occupancy */
    float_type nu_cell_;
    /** average effective cell occupancy */
    float_type nu_cell_eff_;
    /** number of cells per dimension */
    cell_size_type ncell_;
    /** number of placeholders per cell */
    unsigned int cell_size_;
    /** cell edge lengths */
    vector_type cell_length_;
    /** CUDA cell kernel execution configuration */
    cuda::config dim_cell_;
    /** cell lists in global device memory */
    cuda::vector<unsigned int> g_cell_;

    /** GPU radix sort */
    algorithm::gpu::radix_sort<unsigned int> sort_;
    /** cell indices for particles */
    cuda::vector<unsigned int> g_cell_index_;
    /** particle permutation */
    cuda::vector<unsigned int> g_cell_permutation_;
    /** cell offsets in sorted particle list */
    cuda::vector<unsigned int> g_cell_offset_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_BINNING_HPP */
