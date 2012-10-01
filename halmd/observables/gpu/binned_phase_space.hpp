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

#ifndef HALMD_OBSERVABLES_GPU_BINNED_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_GPU_BINNED_PHASE_SPACE_HPP

#include <boost/make_shared.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/observables/binned_phase_space.hpp>
#include <halmd/observables/gpu/samples/binned_phase_space.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace observables {
namespace gpu {

/**
 * bin phase space sample into spatial cells
 *
 * compute index cell lists from a given phase space sample
 */
template <int dimension, typename float_type>
class binned_phase_space
  : public observables::binned_phase_space<dimension>
{
    typedef observables::binned_phase_space<dimension> _Base;
    typedef typename _Base::signal_type signal_type;

public:
    typedef gpu::samples::binned_phase_space<dimension, float_type> sample_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::clock clock_type;
    typedef logger logger_type;
    typedef typename _Base::slot_function_type slot_function_type;

    static void luaopen(lua_State* L);

    binned_phase_space(
        boost::shared_ptr<sample_type> sample
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<clock_type const> clock
      , unsigned int threads
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    virtual void acquire();

    virtual connection on_acquire(slot_function_type const& slot)
    {
        return on_acquire_.connect(slot);
    }

    /** returns midpoint positions of cell grid along specified axis */
    virtual std::vector<double> const& position(unsigned int i) const
    {
        return position_[i];
    }

private:
    typedef typename sample_type::data_sample_type data_sample_type;
    typedef typename sample_type::cell_array_type cell_array_type;
    typedef typename sample_type::cell_size_type cell_size_type;
    typedef typename sample_type::vector_type vector_type;
    typedef boost::array<std::vector<double>, dimension> position_type;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
    };

    /** module dependencies */
    boost::shared_ptr<sample_type> sample_;
    boost::shared_ptr<clock_type const> clock_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;

    /**
     * midpoint positions of cell grid
     *
     * contains one list of grid points per axis
     */
    position_type position_;

    /** number of particles */
    unsigned int npart_;
    /** CUDA execution configuration for kernels */
    cuda::config dim_;
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

    signal_type on_acquire_;
};

} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_BINNED_PHASE_SPACE_HPP */
