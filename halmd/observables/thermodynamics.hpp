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

#ifndef HALMD_OBSERVABLES_THERMODYNAMICS_HPP
#define HALMD_OBSERVABLES_THERMODYNAMICS_HPP

#include <boost/foreach.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace observables {

/**
 * compute thermodynamic state variables such as pressure,
 * temperature, potential energy, total energy
 *
 * potential energy and the potential part of the stress tensor
 * are computed and stored by the force modules
 */

template <int dimension>
class thermodynamics
{
public:
    typedef mdsim::box<dimension> box_type;
    typedef typename mdsim::type_traits<dimension, double>::vector_type vector_type;
    typedef typename signal<void ()>::slot_function_type slot_function_type;

    static void luaopen(lua_State* L);

    thermodynamics(
        boost::shared_ptr<box_type const> box
    );

    // basic quantities with backend-specific evaluation

    /** potential energy per particle */
    virtual double en_pot() = 0;
    /** kinetic energy per particle */
    virtual double en_kin() = 0;
    /** mean velocity per particle */
    virtual vector_type const& v_cm() = 0;
    /** virial sum */
    virtual double virial() = 0;
    /** hypervirial sum */
    virtual double hypervirial() = 0;

    /** clear all data caches */
    virtual void clear_cache() = 0;

    // compute derived quantities on the fly

    /** total pressure */
    double pressure()
    {
        return density() * (temp() + virial() / dimension);
    }

    /** system temperature */
    double temp() { return 2 * en_kin() / dimension; }
    /** particle density */
    double density() { return box_->density(); }
    /** total energy per particle */
    double en_tot() { return en_pot() + en_kin(); }

private:
    /** module dependencies */
    boost::shared_ptr<box_type const> box_;
};

} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_HPP */
