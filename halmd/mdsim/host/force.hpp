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
#include <lua.hpp>

#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/program_options/program_options.hpp>

namespace halmd
{
namespace mdsim { namespace host
{

template <int dimension, typename float_type>
class force
  : public mdsim::force<dimension>
{
public:
    typedef mdsim::force<dimension> _Base;
    typedef type_traits<dimension, float_type> _type_traits;
    typedef typename _type_traits::vector_type vector_type;
    typedef typename _type_traits::stress_tensor_type stress_tensor_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;

    static void luaopen(lua_State* L);

    force() {}
    virtual matrix_type const& cutoff() = 0;
    virtual double potential_energy() = 0;
    virtual stress_tensor_type stress_tensor_pot() = 0;
};

}} // namespace mdsim::host

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCE_HPP */
