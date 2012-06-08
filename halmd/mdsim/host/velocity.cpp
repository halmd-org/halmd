/*
 * Copyright © 2010  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/host/velocity.hpp>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
velocity<dimension, float_type>::velocity(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , logger_(logger)
{
}

/**
 * Rescale magnitude of all velocities by 'factor'
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::rescale(double factor)
{
    cache_proxy<velocity_array_type> velocity = particle_->velocity();
    for (vector_type& v : *velocity) {
        v *= factor;
    }
    LOG("velocities rescaled by factor " << factor);
}

/**
 * Shift all velocities by 'delta'
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::shift(vector_type const& delta)
{
    cache_proxy<velocity_array_type> velocity = particle_->velocity();
    for (vector_type& v : *velocity) {
        v += delta;
    }
}

/**
 * First shift, then rescale all velocities
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::shift_rescale(vector_type const& delta, double factor)
{
    cache_proxy<velocity_array_type> velocity = particle_->velocity();
    for (vector_type& v : *velocity) {
        v += delta;
        v *= factor;
    }
}

#ifndef USE_HOST_SINGLE_PRECISION
template class velocity<3, double>;
template class velocity<2, double>;
#else
template class velocity<3, float>;
template class velocity<2, float>;
#endif

} // namespace mdsim
} // namespace host
} // namespace halmd
