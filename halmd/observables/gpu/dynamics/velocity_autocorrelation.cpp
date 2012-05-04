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

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <cassert>
#include <stdexcept>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/gpu/dynamics/velocity_autocorrelation.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

template <int dimension, typename float_type>
velocity_autocorrelation<dimension, float_type>::velocity_autocorrelation(
    size_t type
  , unsigned int blocks
  , unsigned int threads
)
  // member initialisation
  : type_(type)
  , blocks_(blocks)
  , threads_(threads)
  , compute_(select_compute(threads_))
  , g_acc_(blocks_)
  , h_acc_(blocks_)
{
    LOG("initialise velocity autocorrelation of " << string(1, 'A' + type) << " particles");
}

template <int dimension, typename float_type>
unsigned int velocity_autocorrelation<dimension, float_type>::defaults::blocks() {
    return 32;
}

template <int dimension, typename float_type>
unsigned int velocity_autocorrelation<dimension, float_type>::defaults::threads() {
    return 32 << DEVICE_SCALE;
}

template <int dimension, typename float_type>
typename velocity_autocorrelation<dimension, float_type>::accumulator_type
velocity_autocorrelation<dimension, float_type>::compute(
    sample_type const& first
  , sample_type const& second
)
{
    cuda::configure(blocks_, threads_);
    sample_vector_type const& v1 = *first.v[type_];
    sample_vector_type const& v2 = *second.v[type_];
    assert(v1.size() == v2.size());
    compute_(v1, v2, v1.size(), g_acc_);
    cuda::copy(g_acc_, h_acc_); // implicit synchronize
    return for_each(h_acc_.begin(), h_acc_.end(), accumulator_type());
}

template <int dimension, typename float_type>
typename velocity_autocorrelation<dimension, float_type>::compute_function_type
velocity_autocorrelation<dimension, float_type>::select_compute(unsigned int threads)
{
    switch (threads) {
        case 512:
            return velocity_autocorrelation_wrapper<dimension, 512>::wrapper.compute;
        case 256:
            return velocity_autocorrelation_wrapper<dimension, 256>::wrapper.compute;
        case 128:
            return velocity_autocorrelation_wrapper<dimension, 128>::wrapper.compute;
        case 64:
            return velocity_autocorrelation_wrapper<dimension, 64>::wrapper.compute;
        case 32:
            return velocity_autocorrelation_wrapper<dimension, 32>::wrapper.compute;
        default:
            throw invalid_argument("number of threads");
    }
}


template <typename tcf_type>
static boost::shared_ptr<tcf_type>
wrap_tcf(size_t type, typename tcf_type::sample_type const&)
{
    return make_shared<tcf_type>(type);
}

template <int dimension, typename float_type>
void velocity_autocorrelation<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name(demangled_name<velocity_autocorrelation>());
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<velocity_autocorrelation>(class_name.c_str())

              , def("velocity_autocorrelation", &wrap_tcf<velocity_autocorrelation>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_dynamics_velocity_autocorrelation(lua_State* L)
{
    velocity_autocorrelation<3, float>::luaopen(L);
    velocity_autocorrelation<2, float>::luaopen(L);
    observables::dynamics::correlation<velocity_autocorrelation<3, float> >::luaopen(L);
    observables::dynamics::correlation<velocity_autocorrelation<2, float> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class velocity_autocorrelation<3, float>;
template class velocity_autocorrelation<2, float>;

} // namespace dynamics
} // namespace gpu

namespace dynamics
{

// explicit instantiation
template class correlation<gpu::dynamics::velocity_autocorrelation<3, float> >;
template class correlation<gpu::dynamics::velocity_autocorrelation<2, float> >;

} // namespace dynamics
} // namespace observables
} // namespace halmd
