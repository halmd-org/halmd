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

#ifndef HALMD_MDSIM_HOST_FORCE_HPP
#define HALMD_MDSIM_HOST_FORCE_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>

#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace observables { namespace host
{

// forward declaration
template <int dimension, typename float_type>
class thermodynamics;

}} // namespace observables::host

namespace mdsim { namespace host
{

template <int dimension, typename float_type>
class force
  : public mdsim::force<dimension>
{
public:
    // module definitions
    typedef force _Self;
    typedef mdsim::force<dimension> _Base;
    static void depends();
    static void options(po::options_description& desc) {}
    static void select(po::options const& vm) {}

    typedef type_traits<dimension, float_type> _type_traits;
    typedef typename _type_traits::vector_type vector_type;
    typedef typename _type_traits::stress_tensor_type stress_tensor_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;

    typedef host::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef host::forces::smooth<dimension, float_type> smooth_type;

    shared_ptr<particle_type> particle;
    shared_ptr<box_type> box;
    shared_ptr<smooth_type> smooth;

    force(modules::factory& factory, po::options const& vm);
    virtual ~force() {};
    virtual void compute() = 0;
    virtual matrix_type const& cutoff() = 0;

protected:
    /** average potential energy per particle */
    double en_pot_;
    /** potential part of stress tensor */
    stress_tensor_type stress_pot_;

    friend class observables::host::thermodynamics<dimension, float_type>;
};

/**
 * Trace and off-diagonal elements of distance tensor
 */
template <typename float_type>
typename force<3, float_type>::stress_tensor_type
make_stress_tensor(float_type rr, fixed_vector<float_type, 3> const& r)
{
    typename force<3, float_type>::stress_tensor_type v;
    v[0] = rr;
    v[1] = r[1] * r[2];
    v[2] = r[2] * r[0];
    v[3] = r[0] * r[1];
    return v;
}

template <typename float_type>
typename force<2, float_type>::stress_tensor_type
make_stress_tensor(float_type rr, fixed_vector<float_type, 2> const& r)
{
    typename force<2, float_type>::stress_tensor_type v;
    v[0] = rr;
    v[1] = r[0] * r[1];
    return v;
}

}} // namespace mdsim::host

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCE_HPP */
