/*
 * Copyright © 2008-2014  Felix Höfling
 * Copyright © 2014       Nicolas Höft
 * Copyright © 2008-2012  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_BINNING_HPP
#define HALMD_MDSIM_GPU_BINNING_HPP

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/algorithm/multi_range.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/binning_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/multi_index.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension, typename float_type>
class binning
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::matrix<float> matrix_type;
    typedef mdsim::box<dimension> box_type;

    typedef cuda::vector<unsigned int> array_type;
    typedef fixed_vector<unsigned int, dimension> cell_size_type;
    typedef fixed_vector<int, dimension> cell_diff_type;

    static void luaopen(lua_State* L);

    binning(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<box_type const> box
      , matrix_type const& r_cut
      , double skin
      , double cell_occupancy = 0.5
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    //! returns neighbour list skin in MD units
    float r_skin() const
    {
        return r_skin_;
    }

    /**
     * cell size, i.e. number of particle placeholders per cell
     */
    unsigned int cell_size() const
    {
        return cell_size_;
    }

    /**
     * cells edge lengths
     */
    vector_type cell_length() const
    {
        return cell_length_;
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
    cache<array_type> const& g_cell();

private:
    /** update cell lists */
    void update();
    /** set number of placeholders per cell and reallocate memory */
    void set_cell_size(size_t cell_size);

    std::shared_ptr<particle_type const> particle_;
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger> logger_;

    /** neighbour list skin in MD units */
    float r_skin_;
    /** maximum cutoff */
    float r_cut_max_;
    /** CUDA device properties */
    cuda::device::properties device_properties_;

    /** number of cells per dimension */
    cell_size_type ncell_;
    /** number of placeholders per cell */
    unsigned int cell_size_;
    /** cell edge lengths */
    vector_type cell_length_;
    /** CUDA cell kernel execution configuration */
    cuda::config dim_cell_;

    /** cell lists in global device memory */
    cache<array_type> g_cell_;
    /** cache observer for cell list update */
    cache<> cell_cache_;

    /** cell indices for particles */
    array_type g_cell_index_;
    /** particle permutation */
    array_type g_cell_permutation_;
    /** cell offsets in sorted particle list */
    array_type g_cell_offset_;

    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        utility::profiler::accumulator_type update;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

/**
 * Copy cells to multi-range output iterator.
 *
 * @param output multi-range output iterator
 *
 * A multi-range iterator is a functor that accepts a multi-dimensional
 * index of array type, and returns an output iterator for the given
 * index. The particle indices in the cell of the given index are
 * then copied to the returned output iterator.
 */
template <typename binning_type, typename output_iterator>
inline void
get_cell(binning_type& binning, output_iterator output)
{
    typedef typename binning_type::cell_size_type cell_size_type;
    typedef typename binning_type::array_type array_type;
    array_type const& g_cell = read_cache(binning.g_cell());
    cell_size_type ncell = binning.ncell();
    unsigned int cell_size = binning.cell_size();
    cuda::host::vector<unsigned int> h_cell(g_cell.size());
    cuda::copy(g_cell.begin(), g_cell.end(), h_cell.begin());
    multi_range_for_each(
        cell_size_type(0)
      , ncell
      , [&](cell_size_type const& index) {
            unsigned int offset = multi_index_to_offset(index, ncell);
            auto begin = h_cell.begin() + cell_size * offset;
            std::remove_copy(begin, begin + cell_size, output(index), -1u);
        }
    );
}

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_BINNING_HPP */
