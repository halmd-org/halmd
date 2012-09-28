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

#ifndef HALMD_OBSERVABLES_PROFILES_HPP
#define HALMD_OBSERVABLES_PROFILES_HPP

#include <boost/make_shared.hpp>
#include <lua.hpp>
#include <stdint.h>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace observables {

/**
 * compute profiles for density, stress tensor, ...
 *
 * the potential part of the stress tensor is
 * computed and stored by the force modules
 */

template <int dimension>
class profiles
{
public:
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::clock clock_type;
    typedef logger logger_type;

    typedef typename mdsim::type_traits<dimension, double>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, double>::stress_tensor_type stress_tensor_type;
    typedef typename signal<void ()>::slot_function_type slot_function_type;

    static void luaopen(lua_State* L);

    profiles(
        boost::shared_ptr<box_type const> box
      , boost::shared_ptr<clock_type const> clock
      , fixed_vector<unsigned, dimension> const& ngrid
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    // sample profiles and store with given time stamp
    virtual void sample();

    /** return density profile along given axis */
    std::vector<double> const& density_profile(unsigned axis) const
    {
        return density_profile_[axis];
    }

    /** return profile of the stress tensor diagonal elements along given axis */
    std::vector<vector_type> const& stress_tensor_profile(unsigned axis) const
    {
        return stress_tensor_profile_[axis];
    }

    /** return position grid along given axis */
    std::vector<double> const& position(unsigned axis) const
    {
        return position_[axis];
    }

protected:
    /** compute all profiles */
    virtual void compute_profiles() = 0;

    /** module dependencies */
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<clock_type const> clock_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;

    /** density profiles along each axis */
    boost::array<std::vector<double>, dimension> density_profile_;
    /** profiles of the stress tensor diagonal elements along each axis */
    boost::array<std::vector<vector_type>, dimension> stress_tensor_profile_;
    /** corresponding positions, grid points for each axis
      *
      * grid points specify the centres of the histogram bins
      */
    boost::array<std::vector<double>, dimension> position_;

    /** number of bins for each axis */
    fixed_vector<unsigned, dimension> ngrid_;
    /** grid spacing for each axis */
    vector_type spacing_;
    /** simulation step when data where obtained */
    uint64_t step_;
};

} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_PROFILES_HPP */
