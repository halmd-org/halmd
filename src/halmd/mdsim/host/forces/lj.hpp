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

#ifndef HALMD_MDSIM_HOST_FORCES_LJ_HPP
#define HALMD_MDSIM_HOST_FORCES_LJ_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/host/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

template <int dimension, typename float_type>
class lj
  : public mdsim::force<dimension>
{
public:
    typedef mdsim::force<dimension> _Base;
    typedef typename _Base::matrix_type matrix_type;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef host::box<dimension> box_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    static void resolve(po::options const& vm);
    static po::options_description options();
    lj(po::options const& vm);
    virtual ~lj() {}
    virtual void compute();
    matrix_type const& cutoff() { return r_cut_; }

protected:
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
    using _Base::en_pot_;
    /** average virial per particle for each particle types */
    using _Base::virial_;
};

}}} // namespace mdsim::host::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_LJ_HPP */
