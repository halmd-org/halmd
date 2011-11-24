/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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
#include <limits>
#include <lua.hpp>
#include <stdint.h>
#include <vector>

#include <halmd/mdsim/clock.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/raw_allocator.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace samples {

template <int dimension, typename float_type>
class phase_space
{
private:
    typedef mdsim::clock clock_type;

public:
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef typename clock_type::step_type step_type;

    /** sample vector type for all particles of a species */
    typedef std::vector<vector_type, raw_allocator<vector_type> > sample_vector;
    /** sample vector type for all particles of a species */
    typedef std::vector<unsigned int, raw_allocator<unsigned int> > type_vector;

    /** periodically extended particle positions */
    boost::shared_ptr<sample_vector> r;
    /** particle velocities */
    boost::shared_ptr<sample_vector> v;
    /** particle types */
    boost::shared_ptr<type_vector> type;
    /** simulation step when sample was taken */
    step_type step;

    static void luaopen(lua_State* L);

    phase_space(unsigned int npart);

    /**
     * Free shared pointers and re-allocate memory
     * if containers are shared with some other object.
     *
     * Values are not initialised.
     *
     * @param force if true then enforce reallocation
     */
    void reset(bool force=false);

    /** get read-only particle positions of given type */
    sample_vector const& position() const;
    /** get read-only particle velocities of given type */
    sample_vector const& velocity() const;
    /** get read-only particle types of given type */
    type_vector const& types() const;
    /** get read-write particle positions of given type */
    sample_vector& position();
    /** get read-write particle velocities of given type */
    sample_vector& velocity();
    /** get read-write particle types of given type */
    type_vector& types();
};

template <int dimension, typename float_type>
inline phase_space<dimension, float_type>::phase_space(unsigned int npart)
  // allocate sample pointers
  : r(new sample_vector(npart))
  , v(new sample_vector(npart))
  , type(new type_vector(npart))
  // initialise attributes
  , step(std::numeric_limits<step_type>::max())
{
}

template <int dimension, typename float_type>
inline void phase_space<dimension, float_type>::reset(bool force)
{
    // free shared pointers and re-allocate memory
    if (force || !r.unique()) {
        r.reset(new sample_vector(r->size()));
    }
    if (force || !v.unique()) {
        v.reset(new sample_vector(v->size()));
    }
    if (force || !type.unique()) {
        type.reset(new type_vector(type->size()));
    }
    // make time stamp invalid
    step = std::numeric_limits<step_type>::max();
}

} // namespace samples
} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_SAMPLES_PHASE_SPACE_HPP */
