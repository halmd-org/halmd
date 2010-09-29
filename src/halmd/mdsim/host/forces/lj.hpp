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

#ifndef HALMD_MDSIM_HOST_FORCES_LJ_HPP
#define HALMD_MDSIM_HOST_FORCES_LJ_HPP

#include <boost/assign.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/options.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

template <int dimension, typename float_type>
class lj
  : public mdsim::host::force<dimension, float_type>
{
public:
    static void options(po::options_description& desc);

    typedef mdsim::host::force<dimension, float_type> _Base;
    typedef typename _Base::matrix_type matrix_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::stress_tensor_type stress_tensor_type;

    typedef host::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef host::forces::smooth<dimension, float_type> smooth_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;
    boost::shared_ptr<smooth_type> smooth;

    static boost::array<float, 3> default_cutoff()
    {
        return boost::assign::list_of(2.5f)(2.5f)(2.5f);
    }

    static boost::array<float, 3> default_epsilon()
    {
        return boost::assign::list_of(1.0f)(1.5f)(0.5f);
    }

    static boost::array<float, 3> default_sigma()
    {
        return boost::assign::list_of(1.0f)(0.8f)(0.88f);
    }

    lj(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      // FIXME , boost::shared_ptr<smooth_type> smooth
      , boost::array<float, 3> const& cutoff
      , boost::array<float, 3> const& epsilon
      , boost::array<float, 3> const& sigma
    );
    virtual void compute();

    //! returns potential cutoff distance
    virtual matrix_type const& cutoff()
    {
        return r_cut_;
    }

    //! returns average potential energy per particle
    virtual double potential_energy()
    {
        return en_pot_;
    }

    //! potential part of stress tensor
    virtual stress_tensor_type potential_stress()
    {
        return stress_pot_;
    }

private:
    /** potential well depths in MD units */
    matrix_type epsilon_;
    /** pair separation in MD units */
    matrix_type sigma_;
    /** cutoff length in units of sigma */
    matrix_type r_cut_sigma_;
    /** cutoff length in MD units */
    matrix_type r_cut_;
    /** square of cutoff length */
    matrix_type rr_cut_;
    /** square of pair separation */
    matrix_type sigma2_;
    /** potential energy at cutoff length in MD units */
    matrix_type en_cut_;
    /** average potential energy per particle */
    double en_pot_;
    /** potential part of stress tensor */
    stress_tensor_type stress_pot_;
};

}}} // namespace mdsim::host::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_LJ_HPP */
