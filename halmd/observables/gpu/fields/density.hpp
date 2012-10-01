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

#ifndef HALMD_OBSERVABLES_GPU_FIELDS_DENSITY_HPP
#define HALMD_OBSERVABLES_GPU_FIELDS_DENSITY_HPP

#include <boost/make_shared.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/observables/fields/density.hpp>
#include <halmd/observables/gpu/samples/binned_phase_space.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace fields {

/**
 * compute number density field from binned phase space sample
 */
template <int dimension, typename float_type>
class density
  : public observables::fields::density<dimension>
{
    typedef observables::fields::density<dimension> _Base;
    typedef typename _Base::signal_type signal_type;

public:
    typedef typename _Base::result_type result_type;
    typedef typename _Base::slot_function_type slot_function_type;
    typedef gpu::samples::binned_phase_space<dimension, float_type> sample_type;
    typedef mdsim::clock clock_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    density(
        boost::shared_ptr<sample_type const> sample
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    virtual void sample();

    virtual connection on_sample(slot_function_type const& slot)
    {
        return on_sample_.connect(slot);
    }

    virtual result_type const& value() const
    {
        return density_;
    }

private:
    typedef clock_type::step_type step_type;
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type sample;
    };

    /** module dependencies */
    boost::shared_ptr<sample_type const> sample_;
    boost::shared_ptr<clock_type const> clock_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;

    /** number density field */
    result_type density_;
    /** number density field on device
     *
     *  Use a memory layout that corresponds to density_.
     */
    cuda::vector<typename result_type::element> g_density_;
    /** timestamp of result */
    step_type step_;

    /** profiling runtime accumulators */
    runtime runtime_;

    signal_type on_sample_;
};

} // namespace fields
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_FIELDS_DENSITY_HPP */
