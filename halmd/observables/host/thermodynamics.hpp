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

#ifndef HALMD_OBSERVABLES_HOST_THERMODYNAMICS_HPP
#define HALMD_OBSERVABLES_HOST_THERMODYNAMICS_HPP

#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/observables/thermodynamics.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/particle.hpp>

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
class thermodynamics
    : public observables::thermodynamics<dimension>
{
private:
    typedef observables::thermodynamics<dimension> _Base;

public:
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef typename _Base::box_type box_type;
    typedef typename _Base::clock_type clock_type;
    typedef mdsim::host::force<dimension, float_type> force_type;
    typedef typename _Base::logger_type logger_type;
    typedef typename particle_type::vector_type vector_type;

    static void luaopen(lua_State* L);

    thermodynamics(
        boost::shared_ptr<particle_type const> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<force_type> force
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    virtual void prepare();
    virtual void sample();

    virtual double en_kin();
    virtual vector_type v_cm();

    virtual double en_pot()
    {
        if (!force_->aux_flag()) {
            throw std::logic_error("Potential energy not enabled in force module");
        }
        return force_->potential_energy();
    }

    virtual double virial()
    {
        if (!force_->aux_flag()) {
            throw std::logic_error("Stress tensor not enabled in force module");
        }
        return force_->stress_tensor_pot()[0];
    }

    virtual double hypervirial()
    {
        if (!force_->aux_flag()) {
            throw std::logic_error("Hypervirial not enabled in force module");
        }
        return force_->hypervirial();
    }

private:
    boost::shared_ptr<particle_type const> particle_;
    boost::shared_ptr<force_type> force_;
};

} // namespace observables
} // namespace host
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_HPP */
