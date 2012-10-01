/*
 * Copyright © 2011  Felix Höfling
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

#ifndef HALMD_OBSERVABLES_GPU_SAMPLES_BINNED_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_GPU_SAMPLES_BINNED_PHASE_SPACE_HPP

#include <boost/make_shared.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/gpu/samples/phase_space.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace samples {

/**
 * store binning data for a phase space sample
 */
template <int dimension, typename float_type>
class binned_phase_space
{
public:
    typedef gpu::samples::phase_space<dimension, float_type> data_sample_type;
    typedef typename data_sample_type::vector_type vector_type;
    typedef mdsim::clock::step_type step_type;
    typedef logger logger_type;

    typedef cuda::vector<unsigned int> cell_array_type;
    typedef fixed_vector<unsigned int, dimension> cell_size_type;
    typedef fixed_vector<int, dimension> cell_diff_type;


    /** flattened array of index cells, each cell holds a list of particle indices */
    cell_array_type g_cell;
    /** number of placeholders per cell */
    unsigned int cell_size;
    /** number of bins per dimension */
    cell_size_type nbin;
    /** configuration for kernels operating on g_cell: one CUDA block per cell */
    cuda::config dim;

    /** edge lengths of a binning cell */
    vector_type cell_length;
    /** origin of lower, left cell (0, … 0) */
    vector_type cell_origin;
    /** simulation step when binning was done */
    step_type step;

    static void luaopen(lua_State* L);

    /**
     * construct binned phase space sample of given size
     *
     * @param data_sample phase space sample underlying the binning
     * @param species particle species selected from sample
     * @param nbin number of bins per space dimension
     * @param occupancy average occupancy of binning cells
     *
     * The cell size is computed as [[#particles / #bins] × occupancy]
     * and increased to the next multiple of warp size; [.] denotes rounding up.
     */
    binned_phase_space(
        boost::shared_ptr<data_sample_type const> data_sample
      , unsigned int species
      , fixed_vector<unsigned int, dimension> const& nbin
      , double occupancy
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    typename data_sample_type::sample_vector const& position_data() const
    {
        return *data_sample_->r[species_];
    }

    typename data_sample_type::sample_vector const& velocity_data() const
    {
        return *data_sample_->v[species_];
    }

private:
    /** corresponding phase space sample */
    boost::shared_ptr<data_sample_type const> data_sample_;
    /** particle species selected from phase space sample */
    unsigned int species_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
};

} // namespace samples
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_SAMPLES_BINNED_PHASE_SPACE_HPP */
