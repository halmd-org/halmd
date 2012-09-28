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

#ifndef HALMD_OBSERVABLES_HOST_PROFILES_HPP
#define HALMD_OBSERVABLES_HOST_PROFILES_HPP

#include <boost/make_shared.hpp>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/observables/profiles.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace host {

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
    typedef observables::profiles<dimension> _Base;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef typename _Base::box_type box_type;
    typedef typename _Base::clock_type clock_type;
    typedef mdsim::host::force<dimension, float_type> force_type;
    typedef logger logger_type;

    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::stress_tensor_type stress_tensor_type;

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
    };

    /** compute all profiles */
    void compute_profiles();

    /** module dependencies */
    boost::shared_ptr<particle_type const> particle_;
    boost::shared_ptr<force_type const> force_;
    using _Base::box_;
    using _Base::clock_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;

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

} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_PROFILES_HPP */
