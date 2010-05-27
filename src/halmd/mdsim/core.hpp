/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <halmd/core.hpp>
#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/mdsim/neighbor.hpp>
#include <halmd/mdsim/position.hpp>
#include <halmd/mdsim/sort.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/utility/module.hpp>

namespace halmd
{
namespace mdsim
{

template <int dimension>
class core
  : public halmd::core
{
public:
    typedef halmd::core _Base;
    typedef mdsim::force<dimension> force_type;
    typedef mdsim::neighbor<dimension> neighbor_type;
    typedef mdsim::sort<dimension> sort_type;
    typedef mdsim::integrator<dimension> integrator_type;
    typedef mdsim::position<dimension> position_type;
    typedef mdsim::velocity<dimension> velocity_type;

    static void options(po::options_description& desc);
    static void resolve(po::options const& vm);
    core(po::options const& vm);
    /** Run complete simulation */
    void run();
    /** Initialise simulation */
    void init();
    /** Perform a single MD integration step */
    void mdstep();
    uint64_t steps() { return steps_; }
    double time() { return time_; }

    shared_ptr<force_type> force;
    shared_ptr<neighbor_type> neighbor;
    shared_ptr<sort_type> sort;
    shared_ptr<integrator_type> integrator;
    shared_ptr<position_type> position;
    shared_ptr<velocity_type> velocity;

private:
    /** number of integration steps */
    uint64_t steps_;
    /** integration time in MD units */
    double time_;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_CORE_HPP */
