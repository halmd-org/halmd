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

#ifndef HALMD_MDSIM_BOX_HPP
#define HALMD_MDSIM_BOX_HPP

#include <vector>

#include <halmd/utility/module.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/numeric/host/blas/vector.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim
{

template <int dimension>
class box
{
public:
    typedef box module_type;
    typedef numeric::host::blas::vector<double, dimension> vector_type;
    typedef mdsim::particle<dimension> particle_type;

    boost::shared_ptr<particle_type> particle;

    static void options(po::options_description& desc);
    static void resolve(po::options const& vm);
    box(po::options const& vm);
    virtual ~box() {}
    void length(vector_type const& value_type);
    vector_type const& length() { return length_; }
    void density(double value_type);
    double density() { return density_; }

protected:
    /** edge lengths of cuboid */
    vector_type length_;
    /** edge lengths of cuboid relative to maximum edge length */
    vector_type scale_;
    /** number density */
    double density_;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_BOX_HPP */
