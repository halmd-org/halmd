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

#include <halmd/io/statevars/writer.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace mdsim
{

/**
 * compute thermodynamic state variables such as pressure,
 * temperature, potential energy, total energy
 *
 * potential energy and the potential part of the stress tensor
 * are computed and stored by the force modules
 */

template <int dimension>
class thermodynamics
{
public:
    // module definitions
    typedef thermodynamics _Self;
    static void options(po::options_description& desc);
    static void depends();
    static void select(po::options const& vm);

    typedef mdsim::box<dimension> box_type;
    typedef io::statevars::writer<dimension> writer_type;
    typedef utility::profiler profiler_type;
    typedef fixed_vector<double, dimension> vector_type;

    shared_ptr<box_type> box;
    shared_ptr<writer_type> writer;
    shared_ptr<profiler_type> profiler;

    thermodynamics(modules::factory& factory, po::options const& vm);
    virtual ~thermodynamics() {}

    // sample macroscopic state variables and store with given time
    void sample(double time);

    /** potential energy per particle */
    virtual double en_pot() const = 0;
    /** kinetic energy per particle */
    virtual double en_kin() const = 0;
    /** mean velocity per particle */
    virtual vector_type v_cm() const = 0;
    /** virial sum */
    virtual double virial() const = 0;

    /** total pressure */
    double pressure() const
    {
        return box->density() * (temp() + virial() / dimension);
    }
    /** system temperature */
    double temp() const { return 2 * en_kin() / dimension; }
    /** particle density */
    double density() const { return box->density(); }
    /** total energy per particle */
    double en_tot() const { return en_pot() + en_kin(); }

private:
    // sample() passes values to HDF5 writer via a fixed location in memory
    double en_pot_;
    double en_kin_;
    double en_tot_;
    vector_type v_cm_;
    double pressure_;
    double temp_;
    double density_;
    double time_;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_THERMODYNAMICS_HPP */
