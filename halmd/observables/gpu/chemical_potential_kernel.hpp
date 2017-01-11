/*
 * Copyright © 2017 Felix Höfling
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

#ifndef HALMD_OBSERVABLES_GPU_CHEMICAL_POTENTIAL_KERNEL_HPP
#define HALMD_OBSERVABLES_GPU_CHEMICAL_POTENTIAL_KERNEL_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/accumulator.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

#ifndef __CUDACC__
# include <cmath>
#endif

namespace halmd {
namespace observables {
namespace gpu {

#ifndef __CUDACC__
  using std::isfinite;
#endif

/**
 * Compute canonical partition sum, @f$ \sum_i exp(-U_i / kT) @f$.
 */
template <typename float_type>
class partition_sum
{
    typedef unsigned int size_type;

public:
    /** element pointer type of input array */
    typedef float const* iterator;

    /**
     * Initialise canonical partition sum, store temperature
     */
    partition_sum(float temperature) : Z_(), temperature_(temperature) {}

    /**
     * Accumulate canonical partition sum of a particle.
     */
    HALMD_GPU_ENABLED void operator()(float en_pot)
    {
        if (isfinite(en_pot)) {                           // catch NaN/Inf etc.
            Z_(exp(-en_pot / temperature_));              // prefer exp() over __expf() since
                                                          // the argument can vary strongly TODO check
        }
        else {
            Z_(0);
        }
    }

    /**
     * Accumulate canonical partition sum of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(partition_sum const& acc)
    {
        Z_(acc.Z_);
    }

    /**
     * Returns canonical partition sum.
     */
    HALMD_GPU_ENABLED accumulator<float_type> operator()() const
    {
        return Z_;
    }

private:
    /** canonical partition sum */
    accumulator<float_type> Z_;
    /** assumed system temperature */
    float temperature_;
};

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_CHEMICAL_POTENTIAL_KERNEL_HPP */
