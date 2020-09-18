/*
 * Copyright 2016 Daniel Kirchner
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

#ifndef HALMD_OBSERVABLES_GPU_SAMPLES_SAMPLE_HPP
#define HALMD_OBSERVABLES_GPU_SAMPLES_SAMPLE_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/sample.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace samples {

/**
 * GPU phase_space sample.
 *
 * In order to be able to distinguish between two component data packed in a float4
 * and two component data packed in a float2 (resp. between four component data and
 * three component data packed into a float4), the dimension as well as the _GPU_
 * data type are needed as template parameters.
 */
template<int dimension_, typename data_type_>
class sample
  : public sample_base
{
public:
    static constexpr int dimension = dimension_;
    static constexpr bool gpu_sample = true;

    virtual bool gpu() const {
        return true;
    }

    typedef data_type_ data_type;
    typedef cuda::memory::device::vector<data_type> array_type;

    sample(std::size_t nparticles) : data_(nparticles) {}

    virtual std::type_info const& type() const
    {
        return typeid(data_type);
    }

    array_type const& data() const
    {
        return data_;
    }

    array_type& data()
    {
        return data_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** actual GPU data */
    array_type data_;
};

} // namespace samples
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* !defined HALMD_OBSERVABLES_GPU_SAMPLES_SAMPLE_HPP */
