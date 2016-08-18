/*
 * Copyright 2016 Daniel Kirchner
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HALMD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HALMD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <luaponte/luaponte.hpp>
#include <memory>

#include <halmd/observables/host/samples/sample.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace samples {

template <typename sample_type>
static std::size_t wrap_nparticle(sample_type const& self)
{
    return self.data().size();
}

template <typename sample_type>
static std::size_t wrap_dimension(sample_type const&)
{
    return sample_type::dimension;
}

template <typename sample_type>
static std::function<void(std::vector<typename sample_type::data_type> const&)>
wrap_data_setter(std::shared_ptr<sample_type> self)
{
    return [self](std::vector<typename sample_type::data_type> const& data) {
        if (self->data().size() != data.size()) {
            throw std::runtime_error("phase space sample has mismatching size");
        }
        std::copy(data.begin(), data.end(), self->data().begin());
    };
}

template <typename sample_type>
static std::function<typename sample_type::array_type const&()>
wrap_data_getter(std::shared_ptr<sample_type> self)
{
    return [self]() -> typename sample_type::array_type const& {
        return self->data();
    };
}

template <typename sample_type>
static typename sample_type::data_type wrap_maximum(sample_type const& self) {
    return *std::max_element(self.data().begin(), self.data().end());
}

static std::string space_to_underscore(std::string const& input)
{
    std::string result(input);
    std::replace(result.begin(), result.end(), ' ', '_');
    return result;
}

template <int dimension, typename scalar_type>
void sample<dimension, scalar_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string const class_name = "sample_" + std::to_string(dimension) + "_" + space_to_underscore(demangled_name<scalar_type>());

    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                namespace_("samples")
                [
                    class_<sample, std::shared_ptr<sample> >(class_name.c_str())
                        .def(constructor<std::size_t>())
                        .property("nparticle", &wrap_nparticle<sample>)
                        .property("dimension", &wrap_dimension<sample>)
                        .def("data_setter", &wrap_data_setter<sample>)
                        .def("maximum", &wrap_maximum<sample>)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_samples_sample(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    sample<4, double>::luaopen(L);
    sample<3, double>::luaopen(L);
    sample<2, double>::luaopen(L);
    sample<1, double>::luaopen(L);
#endif
    sample<4, float>::luaopen(L);
    sample<3, float>::luaopen(L);
    sample<2, float>::luaopen(L);
    sample<1, float>::luaopen(L);

    sample<4, int>::luaopen(L);
    sample<3, int>::luaopen(L);
    sample<2, int>::luaopen(L);
    sample<1, int>::luaopen(L);

    sample<4, unsigned int>::luaopen(L);
    sample<3, unsigned int>::luaopen(L);
    sample<2, unsigned int>::luaopen(L);
    sample<1, unsigned int>::luaopen(L);
    return 0;
}


} // namespace samples
} // namespace host
} // namespace observables
} // namespace halmd
