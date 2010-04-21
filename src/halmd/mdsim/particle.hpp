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

#ifndef HALMD_MDSIM_PARTICLE_HPP
#define HALMD_MDSIM_PARTICLE_HPP

#include <vector>

#include <halmd/mdsim/module.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim
{

template <int dimension>
class particle
{
public:
    particle(options const& vm);
    virtual ~particle() {}

    /** number of particles in simulation box */
    unsigned int nbox;
    /** number of particle types */
    unsigned int ntype;
    /** number of particles per type */
    std::vector<unsigned int> ntypes;
};

template <int dimension>
class module<particle<dimension> >
{
public:
    typedef boost::shared_ptr<particle<dimension> > pointer;
    static pointer fetch(options const& vm);

private:
    static pointer singleton_;
};

}} // namespace halmd::mdsim

#endif /* ! HALMD_MDSIM_PARTICLE_HPP */
