/*
 * Copyright © 2013  Nicolas Höft
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

#include <halmd/io/logger.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/observables/utility/accumulator.hpp>

namespace halmd {
namespace observables {
namespace utility {

template <typename accumulator_type>
static std::function<void ()>
wrap_sample(std::shared_ptr<accumulator_type> self)
{
    return [=]() {
        self->sample();
    };
}

template <typename accumulator_type>
static std::function<void ()>
wrap_reset(std::shared_ptr<accumulator_type> self)
{
    return [=]() {
        self->reset();
    };
}

template <typename accumulator_type>
static std::function<typename accumulator_type::value_type ()>
wrap_sum(std::shared_ptr<accumulator_type> self)
{
    return [=]() {
        return self->sum();
    };
}

template <typename accumulator_type>
static std::function<typename accumulator_type::value_type ()>
wrap_mean(std::shared_ptr<accumulator_type> self)
{
    return [=]() {
        return self->mean();
    };
}

template <typename accumulator_type>
static std::function<typename accumulator_type::value_type ()>
wrap_error_of_mean(std::shared_ptr<accumulator_type> self)
{
    return [=]() {
        return self->error_of_mean();
    };
}

template <typename accumulator_type>
static std::function<typename accumulator_type::value_type ()>
wrap_variance(std::shared_ptr<accumulator_type> self)
{
    return [=]() {
        return self->variance();
    };
}

template <typename accumulator_type>
static std::function<typename accumulator_type::accumulator_type::size_type ()>
wrap_count(std::shared_ptr<accumulator_type> self)
{
    return [=]() {
        return self->count();
    };
}


template <typename sample_type>
void accumulator<sample_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("utility")
            [
                class_<accumulator>()
                    .property("sample", &wrap_sample<accumulator<sample_type> >)
                    .property("sum", &wrap_sum<accumulator<sample_type> >)
                    .property("mean", &wrap_mean<accumulator<sample_type> >)
                    .property("error_of_mean", &wrap_error_of_mean<accumulator<sample_type> >)
                    .property("variance", &wrap_variance<accumulator<sample_type> >)
                    .property("count", &wrap_count<accumulator<sample_type> >)
                    .property("reset", &wrap_reset<accumulator<sample_type> >)
              , def("accumulator", &std::make_shared<accumulator,
                    sample_function_type
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_utility_accumulator(lua_State* L)
{
    accumulator<double>::luaopen(L);
    accumulator<fixed_vector<double, 2>>::luaopen(L);
    accumulator<fixed_vector<double, 3>>::luaopen(L);
    accumulator<fixed_vector<double, 4>>::luaopen(L);
    accumulator<fixed_vector<double, 6>>::luaopen(L);
    return 0;
}

// explicit instantiation
template class accumulator<double>;
template class accumulator<fixed_vector<double, 2>>;
template class accumulator<fixed_vector<double, 3>>;
template class accumulator<fixed_vector<double, 4>>;
template class accumulator<fixed_vector<double, 6>>;

} // namespace observables
} // namespace utility
} // namespace halmd
