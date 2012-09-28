/*
 * Copyright © 2010-2011  Felix Höfling
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
#include <boost/foreach.hpp>
#include <functional>

#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/observables/host/profiles.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
profiles<dimension, float_type>::profiles(
    shared_ptr<particle_type const> particle
  , shared_ptr<box_type const> box
  , shared_ptr<clock_type const> clock
  , shared_ptr<force_type const> force
  , fixed_vector<unsigned, dimension> const& ngrid
  , shared_ptr<logger_type> logger
)
  : _Base(box, clock, ngrid, logger)
  , particle_(particle)
  , force_(force)
  , logger_(logger)
{
}

/**
 * compute all profiles
 */
template <int dimension, typename float_type>
void profiles<dimension, float_type>::compute_profiles()
{
    scoped_timer_type timer_(runtime_.sample);

    // initialise accumulators
    BOOST_FOREACH (vector<double>& profile, density_profile_) {
        fill(profile.begin(), profile.end(), 0);
    }
    BOOST_FOREACH (vector<vector_type>& profile, stress_tensor_profile_) {
        fill(profile.begin(), profile.end(), 0);
    }

    // bin particle positions appropriately
    for (size_t i = 0; i < particle_->r.size(); ++i) {
        vector_type const& r = particle_->r[i];
        vector_type origin = box_->origin();
        // stress tensor due to particle 'i'
        stress_tensor_type stress = force_->stress_tensor_pot()[i];
        stress += mdsim::make_stress_tensor(particle_->v[i]);
        for (unsigned axis = 0; axis < dimension; ++axis) {
            size_t index = (r[axis] - origin[axis]) / spacing_[axis];
            index %= density_profile_[axis].size();
            // density profile is simply a histogram
            density_profile_[axis][index] += 1;
            // accumulate stress tensors, but keep diagonal only
            for (unsigned k=0; k < dimension; ++k) {
                stress_tensor_profile_[axis][index][k] += stress[k];
            }
        }
    }

    // normalisation of density and stress profiles by slab volume
    double box_volume = box_->volume();
    for (unsigned axis = 0; axis < dimension; ++axis) {
        double slab_vol = box_volume / box_->length()[axis] * spacing_[axis];
        for (size_t index = 0 ; index < ngrid_[axis]; ++index) {
            stress_tensor_profile_[axis][index] /= slab_vol;
            density_profile_[axis][index] /= slab_vol;
        }
    }
}

template <int dimension, typename float_type>
void profiles<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("profiles_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<profiles, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .def(constructor<
                        shared_ptr<particle_type const>
                      , shared_ptr<box_type const>
                      , shared_ptr<clock_type const>
                      , shared_ptr<force_type const>
                      , fixed_vector<unsigned, dimension> const&
                      , shared_ptr<logger_type>
                    >())
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("sample", &runtime::sample)
                    ]
                    .def_readonly("runtime", &profiles::runtime_)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_profiles(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    profiles<3, double>::luaopen(L);
    profiles<2, double>::luaopen(L);
#else
    profiles<3, float>::luaopen(L);
    profiles<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class profiles<3, double>;
template class profiles<2, double>;
#else
template class profiles<3, float>;
template class profiles<2, float>;
#endif

} // namespace host
} // namespace observables
} // namespace halmd
