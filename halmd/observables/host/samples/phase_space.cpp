/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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
#include <functional>
#include <luaponte/luaponte.hpp>
#include <luaponte/out_value_policy.hpp>
#include <memory>

#include <halmd/observables/host/samples/phase_space.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace samples {

template <typename phase_space_type>
static std::size_t wrap_nparticle(phase_space_type const& self)
{
    return self.position().size();
}

template <typename phase_space_type>
static std::size_t wrap_nspecies(phase_space_type const& self)
{
    return 1 + *std::max_element(self.species().begin(), self.species().end());
}

template <typename phase_space_type>
static std::size_t wrap_dimension(phase_space_type const&)
{
    return phase_space_type::position_array_type::value_type::static_size;
}

template <typename phase_space_type>
static std::function<std::vector<typename phase_space_type::position_array_type::value_type>& ()>
wrap_position(std::shared_ptr<phase_space_type> self, std::function<void ()>& array_to_sample)
{
    typedef std::vector<typename phase_space_type::position_array_type::value_type> array_type;
    std::shared_ptr<array_type> array = std::make_shared<array_type>();
    array_to_sample = [=]() {
        if (self->position().size() != array->size()) {
            throw std::runtime_error("phase space sample has mismatching size");
        }
        std::copy(
            array->begin()
          , array->end()
          , self->position().begin()
        );
        array->clear();
    };
    return [=]() -> array_type& {
        return *array;
    };
}

template <typename phase_space_type>
static std::function<std::vector<typename phase_space_type::velocity_array_type::value_type>& ()>
wrap_velocity(std::shared_ptr<phase_space_type> self, std::function<void ()>& array_to_sample)
{
    typedef std::vector<typename phase_space_type::velocity_array_type::value_type> array_type;
    std::shared_ptr<array_type> array = std::make_shared<array_type>();
    array_to_sample = [=]() {
        if (self->velocity().size() != array->size()) {
            throw std::runtime_error("phase space sample has mismatching size");
        }
        std::copy(
            array->begin()
          , array->end()
          , self->velocity().begin()
        );
        array->clear();
    };
    return [=]() -> array_type& {
        return *array;
    };
}

template <typename phase_space_type>
static std::function<std::vector<typename phase_space_type::species_array_type::value_type>& ()>
wrap_species(std::shared_ptr<phase_space_type> self, std::function<void ()>& array_to_sample)
{
    typedef std::vector<typename phase_space_type::species_array_type::value_type> array_type;
    std::shared_ptr<array_type> array = std::make_shared<array_type>();
    array_to_sample = [=]() {
        if (self->species().size() != array->size()) {
            throw std::runtime_error("phase space sample has mismatching size");
        }
        std::transform(
            array->begin()
          , array->end()
          , self->species().begin()
          , [](typename phase_space_type::species_array_type::value_type s) {
                return s - 1;
            }
        );
        array->clear();
    };
    return [=]() -> array_type& {
        return *array;
    };
}

template <typename phase_space_type>
static std::function<std::vector<typename phase_space_type::mass_array_type::value_type>& ()>
wrap_mass(std::shared_ptr<phase_space_type> self, std::function<void ()>& array_to_sample)
{
    typedef std::vector<typename phase_space_type::mass_array_type::value_type> array_type;
    std::shared_ptr<array_type> array = std::make_shared<array_type>();
    array_to_sample = [=]() {
        if (self->mass().size() != array->size()) {
            throw std::runtime_error("phase space sample has mismatching size");
        }
        std::copy(
            array->begin()
          , array->end()
          , self->mass().begin()
        );
        array->clear();
    };
    return [=]() -> array_type& {
        return *array;
    };
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string const class_name = "phase_space_" + std::to_string(dimension) + "_" + demangled_name<float_type>();
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                namespace_("samples")
                [
                    class_<phase_space, std::shared_ptr<phase_space> >(class_name.c_str())
                        .def(constructor<std::size_t>())
                        .property("nparticle", &wrap_nparticle<phase_space>)
                        .property("nspecies", &wrap_nspecies<phase_space>)
                        .property("dimension", &wrap_dimension<phase_space>)
                        .def("position", &wrap_position<phase_space>, pure_out_value(_2))
                        .def("velocity", &wrap_velocity<phase_space>, pure_out_value(_2))
                        .def("species", &wrap_species<phase_space>, pure_out_value(_2))
                        .def("mass", &wrap_mass<phase_space>, pure_out_value(_2))
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_samples_phase_space(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    phase_space<3, double>::luaopen(L);
    phase_space<2, double>::luaopen(L);
    observables::samples::blocking_scheme<phase_space<3, double> >::luaopen(L);
    observables::samples::blocking_scheme<phase_space<2, double> >::luaopen(L);
#endif
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
    observables::samples::blocking_scheme<phase_space<3, float> >::luaopen(L);
    observables::samples::blocking_scheme<phase_space<2, float> >::luaopen(L);
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class phase_space<3, double>;
template class phase_space<2, double>;
#endif
template class phase_space<3, float>;
template class phase_space<2, float>;

} // namespace samples
} // namespace host

namespace samples
{

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class blocking_scheme<host::samples::phase_space<3, double> >;
template class blocking_scheme<host::samples::phase_space<2, double> >;
#endif
template class blocking_scheme<host::samples::phase_space<3, float> >;
template class blocking_scheme<host::samples::phase_space<2, float> >;

} // namespace samples
} // namespace observables
} // namespace halmd
