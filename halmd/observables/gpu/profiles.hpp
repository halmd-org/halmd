/*
 * Copyright © 2010-2011  Felix Höfling
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

#ifndef HALMD_OBSERVABLES_GPU_PROFILES_HPP
#define HALMD_OBSERVABLES_GPU_PROFILES_HPP

#include <boost/make_shared.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/observables/profiles.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace gpu {

/**
 * compute profiles for density, stress tensor, ...
 *
 * the potential part of the stress tensor is
 * computed and stored by the force modules
 */

template <int dimension, typename float_type>
class profiles
    : public observables::profiles<dimension>
{
public:
    // module definitions
    typedef observables::profiles<dimension> _Base;

    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef typename _Base::box_type box_type;
    typedef typename _Base::clock_type clock_type;
    typedef mdsim::gpu::force<dimension, float_type> force_type;
    typedef logger logger_type;

    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::stress_tensor_type stress_tensor_type;
    typedef mdsim::type_traits<dimension, float_type> _type_traits;
    typedef typename _type_traits::gpu::coalesced_vector_type coalesced_vector_type;
    typedef typename _type_traits::vector_type gpu_vector_type;
    typedef fixed_vector<unsigned, dimension> index_type;

    static void luaopen(lua_State* L);

    profiles(
        boost::shared_ptr<particle_type const> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<force_type const> force
      , fixed_vector<unsigned, dimension> const& ngrid
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type sample;
        accumulator_type bins;
        accumulator_type sort;
        accumulator_type boundaries;
        accumulator_type reduce;
        accumulator_type stress_tensor;
        accumulator_type copy;
    };

    /** module dependencies */
    boost::shared_ptr<particle_type const> particle_;
    boost::shared_ptr<force_type const> force_;
    using _Base::box_;
    using _Base::clock_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;

    /** compute all profiles */
    void compute_profiles();

    /** GPU radix sort */
    algorithm::gpu::radix_sort<unsigned> sort_;
    /** CUDA execution configuration of bin reduction kernel */
    cuda::config dim_bins_;

    /** bin ID of each particle according to its position
      *
      * The simulation box is divided in prod(ngrid_) bins of extent
      * given by spacing_. The bin ID is computed from
      * ID = x + ngrid_[0] * (y + ngrid_[1] * z)
      * where (x, y, z) is the bin index of the three-dimensional grid
      */
    cuda::vector<uint> g_bins_;
    /** permutation of particles indices after sort */
    cuda::vector<uint> g_permutation_;
    /** bin boundaries as a pair of begin/end pointers in g_permutation_
      *
      * negative entries indicate empty bins */
    cuda::vector<int2> g_boundaries_;
    /** bin boundaries in page-locked host memory
      *
      * This variable yields the density profiles finally.
      */
    cuda::host::vector<int2> h_boundaries_;

    /** diagonal of total stress tensor per bin in device memory */
    cuda::vector<coalesced_vector_type> g_stress_diag_;
    /** diagonal of total stress tensor per bin in page-locked host memory */
    cuda::host::vector<coalesced_vector_type> h_stress_diag_;

    /** density profiles along each axis */
    using _Base::density_profile_;
    /** profiles of the stress tensor diagonal elements along each axis */
    using _Base::stress_tensor_profile_;
    /** number of bins for each axis */
    using _Base::ngrid_;
    /** grid spacing for each axis */
    using _Base::spacing_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_PROFILES_HPP */
