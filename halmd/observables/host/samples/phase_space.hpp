/*
 * Copyright © 2008-2011  Peter Colberg
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

#ifndef HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <stdint.h>
#include <vector>

#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace samples {

template <int dimension, typename float_type>
class phase_space
{
public:
    typedef fixed_vector<float_type, dimension> vector_type;

    /** sample vector type for all particles of a species */
    typedef std::vector<vector_type> sample_vector;
    /** sample pointer type for all particle of a species */
    typedef boost::shared_ptr<sample_vector> sample_vector_ptr;
    /** sample pointer type for all species */
    typedef std::vector<sample_vector_ptr> sample_vector_ptr_vector;

    /** periodically extended particle positions */
    sample_vector_ptr_vector r;
    /** particle velocities */
    sample_vector_ptr_vector v;
    /** simulation step when sample was taken */
    uint64_t step;

    static void luaopen(lua_State* L);

    phase_space(std::vector<unsigned int> ntypes);
    /** get read-only particle positions of given type */
    sample_vector const& position(unsigned int type) const;
    /** get read-only particle velocities of given type */
    sample_vector const& velocity(unsigned int type) const;
    /** get read-write particle positions of given type */
    sample_vector& position(unsigned int type);
    /** get read-write particle velocities of given type */
    sample_vector& velocity(unsigned int type);
};

} // namespace observables
} // namespace host
} // namespace samples
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP */
