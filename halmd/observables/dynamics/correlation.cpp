/*
 * Copyright © 2011 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#include <halmd/observables/dynamics/correlation.hpp>

namespace halmd {
namespace observables {
namespace dynamics {

/**
 * Adapt Lua function as correlation function.
 */
class correlation_adaptor
{
public:
    typedef luaponte::object sample_type;
    typedef double result_type;
    typedef accumulator<double> accumulator_type;
    typedef std::function<result_type (sample_type, sample_type)> function_type;

    correlation_adaptor(function_type const& function) : function_(function) {}

    void operator() (sample_type const& first, sample_type const& second, accumulator_type& result) const
    {
        result(function_(first, second));
    }

private:
    function_type function_;
};

void correlation_base::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<correlation_base>()

              , class_<correlation_adaptor, std::shared_ptr<correlation_adaptor>>("correlation_adaptor")
                    .def(constructor<typename correlation_adaptor::function_type>())
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_dynamics_correlation(lua_State* L)
{
    correlation_base::luaopen(L);
    correlation<correlation_adaptor>::luaopen(L);
    return 0;
}

// explicit instantiation
template class correlation<correlation_adaptor>;

} // namespace dynamics
} // namespace observables
} // namespace halmd
