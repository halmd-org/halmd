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

#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/mdsim/neighbor.hpp>
#include <halmd/mdsim/position.hpp>
#include <halmd/mdsim/sort.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim
{

template <int dimension>
class core
{
public:
    core(options const& vm);
    void run();
    uint64_t steps() { return steps_; }
    double time() { return time_; }

    boost::shared_ptr<mdsim::force<dimension> > force;
    boost::shared_ptr<mdsim::neighbor<dimension> > neighbor;
    boost::shared_ptr<mdsim::sort<dimension> > sort;
    boost::shared_ptr<mdsim::integrator<dimension> > integrator;
    boost::shared_ptr<mdsim::position<dimension> > position;
    boost::shared_ptr<mdsim::velocity<dimension> > velocity;

private:
    /** number of integration steps */
    uint64_t steps_;
    /** integration time in MD units */
    double time_;
};

}} // namespace halmd::mdsim

#endif /* ! HALMD_MDSIM_CORE_HPP */
