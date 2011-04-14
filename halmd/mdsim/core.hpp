/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_CORE_HPP
#define HALMD_MDSIM_CORE_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/mdsim/neighbour.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/mdsim/position.hpp>
#include <halmd/mdsim/sort.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/utility/profiler.hpp>

/** HAL’s MD package */
namespace halmd
{
/** Molecular Dynamics simulation modules */
namespace mdsim
{

template <int dimension>
class core
{
public:
    typedef mdsim::particle<dimension> particle_type;
    typedef mdsim::force<dimension> force_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::neighbour<dimension> neighbour_type;
    typedef mdsim::sort<dimension> sort_type;
    typedef mdsim::integrator<dimension> integrator_type;
    typedef mdsim::position<dimension> position_type;
    typedef mdsim::velocity<dimension> velocity_type;
    typedef mdsim::clock clock_type;
    typedef utility::profiler profiler_type;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type prepare;
        accumulator_type mdstep;
    };

    static void luaopen(lua_State* L);

    core();
    void register_runtimes(profiler_type& profiler);
    void prepare();
    void mdstep();

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;
    boost::shared_ptr<force_type> force;
    boost::shared_ptr<neighbour_type> neighbour;
    boost::shared_ptr<sort_type> sort;
    boost::shared_ptr<integrator_type> integrator;
    boost::shared_ptr<position_type> position;
    boost::shared_ptr<velocity_type> velocity;
    const boost::shared_ptr<clock_type> clock;

    //! return current simulation step
    uint64_t step() const {
        return clock->step();
    }

private:
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_CORE_HPP */
