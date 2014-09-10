/*
 * Copyright © 2014 Sutapa Roy
 * Copyright © 2014 Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_EXTERNAL_SLIT_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_EXTERNAL_SLIT_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace external {

/**
 * define slit potential and parameters
 */
template <int dimension, typename float_type>
class slit
{
public:
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_container_type;

    slit(
        float_type const& width
      , vector_type const& offset
      , vector_type const& surface_normal
      , matrix_container_type const& epsilon
      , matrix_container_type const& sigma
      , matrix_container_type const& wetting
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Compute force and potential energy due to slit walls.
     * Form of the potential is given here:
     * u_i(d)=epsilon_i*[(2/15)*(sigma_i/d)**9-wetting_i*(sigma_i/d)**3].
     */
    boost::tuple<vector_type, float_type> operator()(vector_type const& r, unsigned int species) const
    {
        // compute distance to bottom wall
        float_type d = width_2_ + inner_prod(r, surface_normal_) - offset_dot_normal_;
        // energy and force due to bottom wall
        float_type fval1;
        float_type en_pot1;
        tie(fval1, en_pot1) = lj_wall_(d, epsilon_(species, 0), sigma_(species, 0), wetting_(species, 0));

        // compute distance to top wall
        d = width_ - d;
        // energy and force due to top wall
        float_type fval2;
        float_type en_pot2;
        tie(fval2, en_pot2) = lj_wall_(d, epsilon_(species, 1), sigma_(species, 1), wetting_(species, 1));

        vector_type force = (fval1 - fval2 ) * surface_normal_;

        return make_tuple(force, en_pot1 + en_pot2);
    }

    float_type const& width() const
    {
        return width_;
    }

    vector_type const& offset() const
    {
        return offset_;
    }
    
    vector_type const& surface_normal() const
    {
        return surface_normal_;
    }

    matrix_container_type const& epsilon() const
    {
        return epsilon_;
    }
    
    matrix_container_type const& sigma() const
    {
        return sigma_;
    }
    
    matrix_container_type const& wetting() const
    {
        return wetting_;
    }
    
    unsigned int size() const
    {
        return surface_normal_.size();
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** ... */
    boost::tuple<float_type, float_type> lj_wall_(float_type d, float_type epsilon, float_type sigma, float_type wetting) const
    {
      /*
       * energy and force due to top wall
       */
        float_type sigma_d_3 = std::pow(sigma / d, 3);
        float_type en_pot = epsilon * sigma_d_3 * (float_type(2) / 15 * sigma_d_3 * sigma_d_3 - wetting);
        float_type fval = float_type(3) * epsilon * sigma_d_3 * (float_type(2) / 5 * sigma_d_3 * sigma_d_3 - wetting)/d;
        return make_tuple(fval, en_pot);
    }
    
    /** slitwidth in MD units */
    float_type width_;
    /** position of slit centre in MD units */
    vector_type offset_;
    /** outward normal vector for the top wall in MD units */
    vector_type surface_normal_;
    /** interaction strengths for wall potential in MD units */
    matrix_container_type epsilon_;
    /** interaction ranges for wall potential in MD units */
    matrix_container_type sigma_;
    /** wetting parameters for wall potential in MD units */
    matrix_container_type wetting_;
    /** dot product of slit centre and surface normal vector */
    float_type offset_dot_normal_;
    /** half of slit width*/
    float_type width_2_;
    
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace external
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_EXTERNAL_SLIT_HPP */
