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
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace host
{

template <int dimension>
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

    typedef mdsim::host::particle<dimension, double> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename _Base::virial_type virial_type;

    shared_ptr<particle_type> particle;

    thermodynamics(modules::factory& factory, po::options const& vm);
    virtual ~thermodynamics() {}

    double en_kin() const;
    vector_type v_cm() const;
    double en_pot() const { return en_pot_; }
    std::vector<virial_type> const& virial() const { return virial_; }

    template <typename T>
    static virial_type virial_tensor(T rr, vector_type const& r);

    /** average potential energy per particle */
    double en_pot_;
    /** virial tensor for each particle type */
    std::vector<virial_type> virial_;
};

/**
 * Trace and off-diagonal elements of distance tensor
 */

template <> template <typename T>
inline typename thermodynamics<3>::virial_type
thermodynamics<3>::virial_tensor(T rr, typename thermodynamics<3>::vector_type const& r)
{
    thermodynamics<3>::virial_type v;
    v[0] = rr;
    v[1] = r[1] * r[2];
    v[2] = r[2] * r[0];
    v[3] = r[0] * r[1];
    return v;
}

template <> template <typename T>
inline typename thermodynamics<2>::virial_type
thermodynamics<2>::virial_tensor(T rr, typename thermodynamics<2>::vector_type const& r)
{
    thermodynamics<2>::virial_type v;
    v[0] = rr;
    v[1] = r[0] * r[1];
    return v;
}

}} // namespace mdsim::host

} // namespace halmd

#endif /* ! HALMD_MDSIM_THERMODYNAMICS_HPP */
