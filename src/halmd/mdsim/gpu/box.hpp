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

#ifndef HALMD_MDSIM_GPU_BOX_HPP
#define HALMD_MDSIM_GPU_BOX_HPP

#include <vector>

#include <halmd/mdsim/box.hpp>
// #include <halmd/utility/module.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

template <int dimension>
class box
  : public mdsim::box<dimension>
{
public:
    // module definitions
    typedef box _Self;
    typedef mdsim::box<dimension> _Base;
    static void depends();
    static void select(po::options const& vm) {}
    static void options(po::options_description& desc) {}

    typedef typename _Base::vector_type vector_type;
    typedef mdsim::particle<dimension> particle_type;

    shared_ptr<particle_type> particle;

    box(modules::factory& factory, po::options const& vm);
    virtual ~box() {}

protected:
    /** edge lengths of cuboid */
    using _Base::length_;
};

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_BOX_HPP */
