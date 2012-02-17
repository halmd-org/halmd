/*
 * Copyright © 2011 Michael Kopp
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

#ifndef HALMD_MDSIM_HOST_INTEGRATORS_EULER_HPP
#define HALMD_MDSIM_HOST_INTEGRATORS_EULER_HPP

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

template <int dimension, typename float_type>
class euler
  : public mdsim::integrator<dimension>
{
public:
    typedef mdsim::integrator<dimension> _Base;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef logger logger_type;

    static char const* module_name() { return "euler"; }

    static void luaopen(lua_State* L);

    euler(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type const> box
      , double timestep
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    virtual void integrate();
    virtual void finalize() {} // the Euler integrator consists of a single step only

    //! set integration timestep
    virtual void timestep(double timestep);

    //! returns integration timestep
    virtual double timestep() const
    {
        return timestep_;
    }

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type integrate;
    };

    boost::shared_ptr<particle_type> particle_;
    boost::shared_ptr<box_type const> box_;

    /** integration time-step */
    float_type timestep_;
    /** profiling runtime accumulators */
    runtime runtime_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
};

} // namespace integrators
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_INTEGRATORS_EULER_HPP */
