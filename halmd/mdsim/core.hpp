/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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
#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/mdsim/neighbour.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/mdsim/position.hpp>
#include <halmd/mdsim/sort.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/options.hpp>
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
    typedef utility::profiler profiler_type;

    static void options(po::options_description& desc);
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

    // module runtime accumulator descriptions
    HALMD_PROFILE_TAG( prepare_, "microscopic state preparation" );
    HALMD_PROFILE_TAG( mdstep_, "MD integration step" );

    uint64_t step_counter() const { return step_counter_; }

private:
    // list of profiling timers
    boost::fusion::map<
        boost::fusion::pair<prepare_, accumulator<double> >
      , boost::fusion::pair<mdstep_, accumulator<double> >
    > runtime_;

    // MD step counter
    uint64_t step_counter_;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_CORE_HPP */
