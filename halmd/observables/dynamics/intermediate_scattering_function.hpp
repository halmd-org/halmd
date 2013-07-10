/*
 * Copyright © 2013 Felix Höfling
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

#ifndef HALMD_OBSERVABLES_DYNAMICS_INTERMEDIATE_SCATTERING_FUNCTION_HPP
#define HALMD_OBSERVABLES_DYNAMICS_INTERMEDIATE_SCATTERING_FUNCTION_HPP

#include <lua.hpp>
#include <memory>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/raw_array.hpp>

namespace halmd {
namespace observables {
namespace dynamics {

/**
 * Intermediate scattering function
 */
template <int dimension>
class intermediate_scattering_function
{
public:
    typedef raw_array<fixed_vector<double, 2>> sample_type;
    typedef double result_type;
    enum { result_rank = 1 };

    typedef observables::utility::wavevector<dimension> wavevector_type;

    intermediate_scattering_function(
        std::shared_ptr<wavevector_type const> wavevector
      , double norm
    )
      : wavevector_(wavevector)
      , norm_(norm)
      , result_shape_(wavevector->shell().size())
    {}

    static void luaopen(lua_State* L);

    /**
     * Compute time correlation from two density mode samples
     *
     * @param rho1 density mode sample at initial time t1
     * @param rho2 density mode sample at later time t2
     * @param result returns S(k,t) for lag time t = t2 - t1
     *
     * If instantiated with dynamics::correlation<>, the template argument is
     * boost::multi_array<accumulator<double>, 3>::subarray<1>::value
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
    /** normalisation factor */
    double norm_;
    /** shape of result array */
    unsigned int result_shape_;
};

template <int dimension> template <typename MultiArray>
void intermediate_scattering_function<dimension>::operator() (
    sample_type const& first
  , sample_type const& second
  , MultiArray&& result
) const
{
    // accumulate products of density modes with equal wavenumber,
    // iterate over wavevector shells encoded as index ranges to the wavevector array
    auto output = result.begin();
    for (auto idx_range : wavevector_->shell()) {
        accumulator<result_type> acc;
        auto rho1 = begin(first) + idx_range.first;
        auto rho2 = begin(second) + idx_range.first;
        // iterate over wavevectors and density modes simultaneously
        for (size_t i = idx_range.first; i != idx_range.second; ++i, ++rho1, ++rho2) {
            // compute Re[rho1 (rho2)*]
            double re = ((*rho1)[0] * (*rho2)[0] + (*rho1)[1] * (*rho2)[1]);
            // accumulate results for this wavenumber
            acc(re / norm_);
        }
        // write to output iterator, which accumulates the result
        (*output++)(acc);
    }
}

} // namespace dynamics
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_INTERMEDIATE_SCATTERING_FUNCTION_HPP */
