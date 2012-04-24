/*
 * Copyright © 2012  Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_DYNAMICS_TAGGED_PARTICLE_HPP
#define HALMD_OBSERVABLES_GPU_DYNAMICS_TAGGED_PARTICLE_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/accumulator.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

template <typename correlation_function, typename output_type>
class tagged_particle
{
private:
    typedef typename correlation_function::vector_type vector_type;

public:
    typedef float4 first_argument_type;
    typedef float4 second_argument_type;

    HALMD_GPU_ENABLED void operator()(float4 const& first, float4 const& second)
    {
        acc_(correlation_function()(vector_type(first), vector_type(second)));
    }

    HALMD_GPU_ENABLED void operator()(tagged_particle const& acc)
    {
        acc_(acc.acc_);
    }

    accumulator<output_type> const& operator()() const
    {
        return acc_;
    }

private:
    accumulator<output_type> acc_;
};

} // namespace dynamics
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_DYNAMICS_TAGGED_PARTICLE_HPP */
