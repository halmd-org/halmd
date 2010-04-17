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

#ifndef HALMD_MDSIM_FORCE_HPP
#define HALMD_MDSIM_FORCE_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <vector>

#include <halmd/mdsim/module.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/numeric/host/blas/vector.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim
{

template <int dimension>
class force
{
public:
    typedef numeric::host::blas::vector<double, dimension> vector_type;
    typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> matrix_type;
    typedef numeric::host::blas::vector<double, 1 + (dimension - 1) * dimension / 2> virial_type;
    typedef mdsim::particle<dimension> particle_type;

    force(options const& vm);
    virtual ~force() {}
    virtual void compute() = 0;
    virtual matrix_type const& cutoff() = 0;

    double en_pot() { return en_pot_; }
    std::vector<virial_type> const& virial() { return virial_; }

    boost::shared_ptr<particle_type> particle;

protected:
    /** average potential energy per particle */
    double en_pot_;
    /** average virial per particle for each particle type */
    std::vector<virial_type> virial_;
};

template <int dimension>
class module<force<dimension> >
{
public:
    typedef boost::shared_ptr<force<dimension> > pointer;
    static pointer fetch(options const& vm);

private:
    static pointer singleton_;
};

}} // namespace halmd::mdsim

#endif /* ! HALMD_MDSIM_FORCE_HPP */
