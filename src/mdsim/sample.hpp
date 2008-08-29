/* Phase space sample
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef MDSIM_SAMPLE_HPP
#define MDSIM_SAMPLE_HPP

#include <boost/function.hpp>
#include <vector>
#include "config.hpp"

namespace mdsim {

/**
 * MD simulation sample
 */
struct trajectory_sample
{
    /** trajectory sample visitor type */
    typedef boost::function<void (std::vector<hvector>&, std::vector<hvector>&)> visitor;

    /** periodically reduced particle positions */
    std::vector<hvector> r;
    /** periodically extended particle positions */
    std::vector<hvector> R;
    /** particle velocities */
    std::vector<hvector> v;
    /** potential energy per particle */
    double en_pot;
    /** virial equation sum per particle */
    double virial;
};

} // namespace mdsim

#endif /* ! MDSIM_SAMPLE_HPP */
