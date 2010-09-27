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

#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/mdsim/neighbour.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/mdsim/position.hpp>
#include <halmd/mdsim/sort.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/options.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace mdsim
{

template <int dimension>
class core
{
public:
    // module definitions
    typedef core _Self;
    static void options(po::options_description& desc);
    static void depends();
    static void select(po::variables_map const& vm);

    typedef mdsim::force<dimension> force_type;
    typedef mdsim::neighbour<dimension> neighbour_type;
    typedef mdsim::sort<dimension> sort_type;
    typedef mdsim::integrator<dimension> integrator_type;
    typedef mdsim::particle<dimension> particle_type;
    typedef mdsim::position<dimension> position_type;
    typedef mdsim::velocity<dimension> velocity_type;
    typedef utility::profiler profiler_type;

    core(modules::factory& factory, po::variables_map const& vm);
    void register_runtimes(profiler_type& profiler);
    void prepare();
    void mdstep();

    shared_ptr<profiler_type> profiler;
    shared_ptr<force_type> force;
    shared_ptr<neighbour_type> neighbour;
    shared_ptr<sort_type> sort;
    shared_ptr<integrator_type> integrator;
    shared_ptr<particle_type> particle;
    shared_ptr<position_type> position;
    shared_ptr<velocity_type> velocity;

    // module runtime accumulator descriptions
    HALMD_PROFILE_TAG( prepare_, "microscopic state preparation" );
    HALMD_PROFILE_TAG( mdstep_, "MD integration step" );

    uint64_t step_counter() const { return step_counter_; }
    double time() const { return step_counter_ * integrator->timestep(); }

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
