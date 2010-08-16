/*
 * Copyright © 2010  Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_THERMODYNAMICS_HPP
#define HALMD_MDSIM_HOST_THERMODYNAMICS_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <vector>

#include <halmd/mdsim/thermodynamics.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace host
{

template <int dimension, typename float_type>
class thermodynamics
    : public mdsim::thermodynamics<dimension>
{
public:
    // module definitions
    typedef thermodynamics _Self;
    typedef mdsim::thermodynamics<dimension> _Base;
    static void depends();
    static void options(po::options_description& desc) {}
    static void select(po::options const& vm) {}

    typedef host::particle<dimension, float_type> particle_type;
    typedef host::force<dimension, float_type> force_type;

    typedef typename particle_type::vector_type vector_type;
    typedef typename _Base::virial_type virial_type;

    shared_ptr<particle_type> particle;
    shared_ptr<force_type> force;

    thermodynamics(modules::factory& factory, po::options const& vm);
    virtual ~thermodynamics() {}

    double en_kin() const;
    vector_type v_cm() const;
    double en_pot() const { return force->en_pot_; }
    std::vector<virial_type> const& virial() { return force->virial_; }
};

}} // namespace mdsim::host

} // namespace halmd

#endif /* ! HALMD_MDSIM_THERMODYNAMICS_HPP */
