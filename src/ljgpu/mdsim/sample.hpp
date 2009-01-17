/* Phase space sample
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_MDSIM_SAMPLE_HPP
#define LJGPU_MDSIM_SAMPLE_HPP

#include <boost/array.hpp>
#include <boost/function.hpp>
#include <ljgpu/mdsim/impl.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <vector>

namespace ljgpu {

/**
 * phase space sample for A or B particles
 */
template <typename float_type, int dimension>
struct phase_space_sample
{
    typedef vector<double, dimension> position_vector;
    typedef vector<float_type, dimension> velocity_vector;
    typedef std::vector<position_vector> position_sample_vector;
    typedef std::vector<velocity_vector> velocity_sample_vector;

    /** periodically extended particle positions */
    position_sample_vector r;
    /** particle velocities */
    velocity_sample_vector v;
};

/**
 * MD simulation sample for A and B particles
 */
template <typename float_type, int dimension>
struct mdsim_sample
: public boost::array<phase_space_sample<float_type, dimension>, 2>
{
    typedef boost::function<void (mdsim_sample<float_type, dimension>&)> sample_visitor;
    typedef phase_space_sample<float_type, dimension> uniform_sample;
    typedef typename uniform_sample::position_vector position_vector;
    typedef typename uniform_sample::velocity_vector velocity_vector;
    typedef typename uniform_sample::position_sample_vector position_sample_vector;
    typedef typename uniform_sample::velocity_sample_vector velocity_sample_vector;

    /** potential energy per particle */
    double en_pot;
    /** virial equation sum per particle */
    double virial;
    /** simulation box length */
    float box;
};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_SAMPLE_HPP */
