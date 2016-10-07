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

#ifndef HALMD_OBSERVABLES_HOST_SAMPLES_SAMPLE_HPP
#define HALMD_OBSERVABLES_HOST_SAMPLES_SAMPLE_HPP

#include <lua.hpp>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/sample.hpp>
#include <halmd/utility/raw_array.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace samples {

template<int dimension_, typename scalar_type>
class sample
  : public sample_base
{
public:
    static constexpr int dimension = dimension_;
    static constexpr bool gpu = false;

    typedef typename std::conditional<dimension == 1, scalar_type, fixed_vector<scalar_type, dimension>>::type data_type;
    typedef raw_array<data_type> array_type;

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

    data_type const& maximum() const
    {
        return *std::max_element(data_.begin(), data_.end());
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    array_type data_;
};

} // namespace samples
} // namespace host
} // namespace observables
} // namespace halmd

#endif /* !defined HALMD_OBSERVABLES_HOST_SAMPLES_SAMPLE_HPP */
