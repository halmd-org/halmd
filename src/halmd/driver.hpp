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

#ifndef HALMD_DRIVER_HPP
#define HALMD_DRIVER_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/core.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/utility/module.hpp>

namespace halmd
{

/**
 * This class drives the MD integration and the evaluation.
 */
template <int dimension>
class driver
  : public halmd::core
{
public:
    typedef halmd::core _Base;
    typedef mdsim::core<dimension> core_type;

    static void options(po::options_description& desc);
    static void resolve(po::options const& vm);
    driver(po::options const& vm);
    void run();
    uint64_t steps() { return steps_; }
    double time() { return time_; }

    shared_ptr<core_type> core;

protected:
    /** number of integration steps */
    uint64_t steps_;
    /** integration time in MD units */
    double time_;
};

} // namespace halmd

#endif /* ! HALMD_DRIVER_HPP */
