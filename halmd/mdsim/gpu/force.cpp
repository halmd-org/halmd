/*
 * Copyright © 2013 Felix Höfling
 * Copyright © 2013 Nicolas Höft
 * Copyright © 2012 Peter Colberg
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

#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <functional>
#include <iterator>
#include <utility>
#include <vector>

namespace halmd {
namespace mdsim {
namespace gpu {

template <typename force_type>
static std::function<std::vector<typename force_type::net_force_type> ()>
wrap_get_net_force(std::shared_ptr<force_type> self)
{
    return [=]() -> std::vector<typename force_type::net_force_type> {
        std::vector<typename force_type::net_force_type> output;
        {
            output.reserve(self->net_force()->size());
        }
        get_net_force(*self, std::back_inserter(output));
        return std::move(output);
    };
}

template <typename force_type>
static std::function<std::vector<typename force_type::en_pot_type> ()>
wrap_get_en_pot(std::shared_ptr<force_type> self)
{
    return [=]() -> std::vector<typename force_type::en_pot_type> {
        std::vector<typename force_type::en_pot_type> output;
        {
            output.reserve(self->en_pot()->size());
        }
        get_en_pot(*self, std::back_inserter(output));
        return std::move(output);
    };
}

template <typename force_type>
static std::function<std::vector<typename force_type::stress_pot_type> ()>
wrap_get_stress_pot(std::shared_ptr<force_type> self)
{
    return [=]() -> std::vector<typename force_type::stress_pot_type> {
        std::vector<typename force_type::stress_pot_type> output;
        {
            output.reserve(self->stress_pot()->size());
        }
        get_stress_pot(*self, std::back_inserter(output));
        return std::move(output);
    };
}

template <typename force_type>
static std::function<std::vector<typename force_type::hypervirial_type> ()>
wrap_get_hypervirial(std::shared_ptr<force_type> self)
{
    return [=]() -> std::vector<typename force_type::hypervirial_type> {
        std::vector<typename force_type::hypervirial_type> output;
        {
            output.reserve(self->hypervirial()->size());
        }
        get_hypervirial(*self, std::back_inserter(output));
        return std::move(output);
    };
}

template <int dimension, typename float_type>
void force<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L)
    [
        class_<force>()
            .property("get_net_force", &wrap_get_net_force<force>)
            .property("get_en_pot", &wrap_get_en_pot<force>)
            .property("get_stress_pot", &wrap_get_stress_pot<force>)
            .property("get_hypervirial", &wrap_get_hypervirial<force>)
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_force(lua_State* L)
{
    force<3, float>::luaopen(L);
    force<2, float>::luaopen(L);
    return 0;
}

} // namespace gpu
} // namespace mdsim
} // namespace halmd
