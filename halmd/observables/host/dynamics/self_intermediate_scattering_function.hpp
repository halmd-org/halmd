/*
 * Copyright © 2024 Felix Höfling
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
 * <http://www.gnu.org/licenses/.
 */

#ifndef HALMD_OBSERVABLES_HOST_DYNAMICS_SELF_INTERMEDIATE_SCATTERING_FUNCTION_HPP
#define HALMD_OBSERVABLES_HOST_DYNAMICS_SELF_INTERMEDIATE_SCATTERING_FUNCTION_HPP

#include <lua.hpp>
#include <memory>

#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/observables/host/samples/sample.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace dynamics {

/**
 * Intermediate scattering function
 */
template <int dimension, typename float_type>
class self_intermediate_scattering_function
{
public:
    typedef observables::host::samples::sample<dimension, float_type> sample_type;
    typedef typename sample_type::data_type vector_type;
    typedef double result_type;
    enum { result_rank = 1 };

    typedef observables::utility::wavevector<dimension> wavevector_type;

    self_intermediate_scattering_function(
        std::shared_ptr<wavevector_type const> wavevector
    )
      : wavevector_(wavevector)
      , result_shape_(wavevector->shell().size())
    {}

    static void luaopen(lua_State* L);

    /**
     * Compute time correlation from two phase space samples
     *
     * @param first phase space sample at initial time t1
     * @param second phase space sample at later time t2
     * @param result returns F(k,t) for lag time t = t2 - t1
     *
     * If instantiated with dynamics::correlation<>, the template parameter is
     */
    template <typename MultiArray>
    void operator() (sample_type const& first, sample_type const& second, MultiArray&& result) const;

    /**
     * Return shape of result array.
     */
    unsigned int const* result_shape() const
    {
        return &result_shape_;
    }

private:
    /** wavevector module */
    std::shared_ptr<wavevector_type const> wavevector_;
    /** shape of result array */
    unsigned int result_shape_;
};

} // namespace dynamics
} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_DYNAMICS_SELF_INTERMEDIATE_SCATTERING_FUNCTION_HPP */
