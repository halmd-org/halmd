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
#include <cmath>
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
    typedef boost::numeric::ublas::vector<float_type> scalar_container_type;
    typedef boost::numeric::ublas::vector<vector_type> vector_container_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_container_type;

    slit(
        scalar_container_type const& offset
      , vector_container_type const& surface_normal
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
        float_type en_pot = 0;
        vector_type force = 0;

        // loop over walls
        for (unsigned int i = 0; i < surface_normal_.size(); ++i) {

          // compute absolute distance to wall i
          float_type d = inner_prod(r, surface_normal_(i)) - offset_(i);
          
          float_type epsilon = epsilon_(i, species);
          float_type sigma = sigma_(i, species);
          float_type c = wetting_(i, species);

          // energy and force due to wall i
          float_type d3i = std::pow(sigma / d, 3);
          float_type d6i = float_type(2) / 15 * d3i * d3i;
          float_type eps_d3i = std::copysign(1, d) * epsilon * d3i;
          float_type fval = 3 * eps_d3i * (3 * d6i - c) / d;
          force += fval * surface_normal_(i);
          en_pot += eps_d3i * (d6i - c);
        }

        return make_tuple(force, en_pot);
    }

    scalar_container_type const& offset() const
    {
        return offset_;
    }
    
    vector_container_type const& surface_normal() const
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
        return offset_.size();
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** wall positions in MD units */
    scalar_container_type offset_;
    /** wall normal vectors in MD units */
    vector_container_type surface_normal_;
    /** interaction strengths for wall potential in MD units */
    matrix_container_type epsilon_;
    /** interaction ranges for wall potential in MD units */
    matrix_container_type sigma_;
    /** wetting parameters for wall potential in MD units */
    matrix_container_type wetting_;
    
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace external
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_EXTERNAL_SLIT_HPP */
