/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_HOST_DYNAMICS_VELOCITY_AUTOCORRELATION_HPP
#define HALMD_OBSERVABLES_HOST_DYNAMICS_VELOCITY_AUTOCORRELATION_HPP

#include <lua.hpp>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/dynamics/velocity_autocorrelation.hpp>
#include <halmd/observables/host/samples/sample.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace dynamics {

/**
 * Velocity autocorrelation.
 */
template <int dimension, typename float_type>
class velocity_autocorrelation
{
public:
    typedef host::samples::sample<dimension, float_type> sample_type;
    typedef typename sample_type::data_type vector_type;
    typedef double result_type;

    static void luaopen(lua_State* L);

    void operator() (sample_type const& first, sample_type const& second, accumulator<result_type>& result);

private:
    typedef observables::dynamics::velocity_autocorrelation<dimension, float_type> correlate_function_type;
};

} // namespace dynamics
} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_DYNAMICS_VELOCITY_AUTOCORRELATION_HPP */
