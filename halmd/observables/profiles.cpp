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

#include <limits>

#include <halmd/observables/profiles.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {

template <int dimension>
profiles<dimension>::profiles(
    shared_ptr<box_type const> box
  , shared_ptr<clock_type const> clock
  , fixed_vector<unsigned, dimension> const& ngrid
  , shared_ptr<logger_type> logger
)
  // dependency injection
  : box_(box)
  , clock_(clock)
  , logger_(logger)
  // initialise members
  , ngrid_(ngrid)
  , spacing_(element_div(box_->length(), static_cast<vector_type>(ngrid_)))
  , step_(numeric_limits<uint64_t>::max())
{
    // setup result vectors, no need to initialise values
    for (unsigned axis=0; axis < dimension; ++axis) {
        unsigned n = ngrid_[axis];
        position_[axis].resize(n);
        density_profile_[axis].resize(n);
        stress_tensor_profile_[axis].resize(n);
    }

    // setup position grid for each axis
    vector_type origin = box_->origin();
    for (unsigned axis=0; axis < dimension; ++axis) {
        for (unsigned i=0; i < ngrid_[axis]; ++i) {
            position_[axis][i] = (i + 0.5) * spacing_[axis] + origin[axis];
        }
    }
}

/**
 * Sample all profiles
 */
template <int dimension>
void profiles<dimension>::sample()
{
    if (step_ == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return;
    }

    LOG_TRACE("acquire sample");

    compute_profiles();
    step_ = clock_->step();
}

template <typename profiles_type>
static typename profiles_type::slot_function_type
sample_wrapper(shared_ptr<profiles_type> profiles)
{
    return bind(&profiles_type::sample, profiles);
}

template <typename profiles_type>
static function<vector<double> const& ()>
wrap_density_profile(shared_ptr<profiles_type> profiles, unsigned axis)
{
    return bind(&profiles_type::density_profile, profiles, axis);
}

template <typename profiles_type>
static function<vector<typename profiles_type::vector_type> const& ()>
wrap_stress_tensor_profile(shared_ptr<profiles_type> profiles, unsigned axis)
{
    return bind(&profiles_type::stress_tensor_profile, profiles, axis);
}

template <typename profiles_type>
static function<vector<double> const& ()>
wrap_position(shared_ptr<profiles_type> profiles, unsigned axis)
{
    return bind(&profiles_type::position, profiles, axis);
}

template <int dimension>
void profiles<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("profiles_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<profiles, shared_ptr<profiles> >(class_name.c_str())
                .property("sample", &sample_wrapper<profiles>)
                .def("density_profile", &wrap_density_profile<profiles>)
                .def("stress_tensor_profile", &wrap_stress_tensor_profile<profiles>)
                .def("position", &wrap_position<profiles>)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_profiles(lua_State* L)
{
    profiles<3>::luaopen(L);
    profiles<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class profiles<3>;
template class profiles<2>;

} // namespace observables
} // namespace halmd
