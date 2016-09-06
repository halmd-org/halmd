/*
 * Copyright 2016 Daniel Kirchner
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HALMD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HALMD.  If not, see <http://www.gnu.org/licenses/>.
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

template<int dimension_, typename data_type_>
class sample : public sample_base {
public:
    static constexpr int dimension = dimension_;
    static constexpr bool gpu = true;

    typedef data_type_ data_type;
    typedef cuda::vector<data_type> array_type;

    sample(std::size_t nparticles) : data_(nparticles)
    {
    }

    virtual std::type_info const& type() const {
        return typeid(gpu_sample<data_type>);
    }

    array_type const& data() const {
        return data_;
    }

    array_type& data() {
        return data_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    array_type data_;
};

} // namespace samples
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* !defined HALMD_OBSERVABLES_GPU_SAMPLES_SAMPLE_HPP */
