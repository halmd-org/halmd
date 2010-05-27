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

#ifndef HALMD_MDSIM_THERMODYNAMICS_HPP
#define HALMD_MDSIM_THERMODYNAMICS_HPP

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <vector>

#include <halmd/mdsim/box.hpp>
#include <halmd/numeric/host/blas/vector.hpp>
// #include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim
{

/**
 * compute and store thermodynamic state variables such as pressure,
 * temperature, potential energy, total energy
 */

template <int dimension>
class thermodynamics
{
public:
    typedef mdsim::box<dimension> box_type;
    typedef numeric::host::blas::vector<double, dimension> vector_type;
    typedef numeric::host::blas::vector<double, 1 + (dimension - 1) * dimension / 2> virial_type;

    shared_ptr<box_type> box;

    static void options(po::options_description& desc) {}
    static void resolve(po::options const& vm);
    thermodynamics(po::options const& vm);
    virtual ~thermodynamics() {}

    /** potential energy per particle */
    virtual double en_pot() const = 0;
    /** kinetic energy per particle */
    virtual double en_kin() const = 0;
    /** mean velocity per particle */
    virtual vector_type v_cm() const = 0;
    /** virial tensor per particle for each particle type */
    virtual std::vector<virial_type> const& virial() const = 0;

    /** total pressure */
    double pressure() const;
    /** system temperature */
    double temp() const { return 2 * en_kin() / dimension; }
    /** particle density */
    double density() const { return box->density(); }
    /** total energy per particle */
    double en_tot() const { return en_pot() + en_kin(); }
};

template <int dimension>
inline double thermodynamics<dimension>::pressure() const
{
    double virial_sum = 0;
    // store reference to avoid redundant calls to virtual function
    std::vector<virial_type> const& vir = virial();
    BOOST_FOREACH(virial_type const& v, vir) {
        virial_sum += v[0];
    }
    virial_sum /= dimension;
    // TODO: check '-' sign
    return box->density() * (temp() - virial_sum);
}

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_THERMODYNAMICS_HPP */
