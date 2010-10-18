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

#ifndef HALMD_MDSIM_HOST_FORCES_POWER_LAW_HPP
#define HALMD_MDSIM_HOST_FORCES_POWER_LAW_HPP

#include <lua.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/force.hpp>
// FIXME #include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/options.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**
 * A power-law potential @f$r^{-n}@f$ is often used for
 * repulsive smooth spheres. A big advantage is
 * its scale invariance (in the absence of a cutoff).
 */

template <int dimension, typename float_type>
class power_law
  : public mdsim::host::force<dimension, float_type>
{
public:
    typedef mdsim::host::force<dimension, float_type> _Base;
    typedef typename _Base::matrix_type matrix_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::stress_tensor_type stress_tensor_type;

    typedef host::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    // FIXME typedef host::forces::smooth<dimension, float_type> smooth_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;
    // FIXME boost::shared_ptr<smooth_type> smooth;

    //! default power-law index
    static int const default_index;

    static void options(po::options_description& desc);
    static void luaopen(lua_State* L);

    power_law(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      // FIXME , boost::shared_ptr<smooth_type> smooth
      , int index
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
    /** power law index */
    int index_;
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

    /** optimise pow() function by providing the index at compile time */
    template <int index>
    void compute_impl();
};

}}} // namespace mdsim::host::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_POWER_LAW_HPP */
