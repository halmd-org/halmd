/*
 * Copyright © 2011-2014 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <functional>

#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/utility/raw_array.hpp>

namespace halmd {
namespace observables {
namespace dynamics {

/**
 * Adapt Lua function as correlation function.
 *
 * General case for tensor-valued correlation functions of rank N.
 */
template <unsigned int N>
class correlation_adaptor
{
public:
    typedef luaponte::object sample_type;
    typedef double result_type;
    enum { result_rank = N };

    typedef std::function<boost::multi_array<result_type, N> (sample_type, sample_type)> function_type;
    typedef boost::array<unsigned int, N> shape_type;

    correlation_adaptor(function_type const& function, shape_type const& result_shape)
        : function_(function), result_shape_(result_shape) {}

    template <typename MultiArray>
    void operator() (sample_type const& first, sample_type const& second, MultiArray&& result) const
    {
        // evaluate correlation function
        auto const& value = function_(first, second);
        // check for size to prevent memory access violoations
        if (value.num_elements() != result.num_elements()) {
            throw std::logic_error("result of correlation function has mismatching shape");
        }

        // accumulate result element-wise
        auto output = result.origin();
        for (unsigned int i = 0; i < value.num_elements(); ++i) {
            (*output++)(value.origin()[i]);
        }
    }

    /**
     * Return shape of result array.
     */
    unsigned int const* result_shape() const
    {
        return result_shape_.data();
    }

private:
    function_type function_;
    shape_type result_shape_;
};

/**
 * Specialisation for correlation functions that yield a scalar
 */
template <>
class correlation_adaptor<0>
{
public:
    typedef luaponte::object sample_type;
    typedef double result_type;

    typedef std::function<result_type (sample_type, sample_type)> function_type;
    typedef boost::array<unsigned int, 0> shape_type;

    correlation_adaptor(function_type const& function, shape_type const& shape = shape_type())
      : function_(function) {}

    void operator() (sample_type const& first, sample_type const& second, accumulator<result_type>& result) const
    {
        result(function_(first, second));
    }

private:
    function_type function_;
};

HALMD_LUA_API int luaopen_libhalmd_observables_dynamics_correlation_adaptor(lua_State* L)
{
    using namespace luaponte;

    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<correlation_adaptor<0>>()
              , def("correlation_adaptor", &std::make_shared<correlation_adaptor<0>,
                    typename correlation_adaptor<0>::function_type const&
                  , typename correlation_adaptor<0>::shape_type const&
                >)

              , class_<correlation_adaptor<1>>()
              , def("correlation_adaptor", &std::make_shared<correlation_adaptor<1>,
                    typename correlation_adaptor<1>::function_type const&
                  , typename correlation_adaptor<1>::shape_type const&
                >)
            ]
        ]
    ];

    correlation<correlation_adaptor<0>>::luaopen(L);
    correlation<correlation_adaptor<1>>::luaopen(L);

    return 0;
}

// explicit instantiation
template class correlation<correlation_adaptor<0>>;
template class correlation<correlation_adaptor<1>>;

} // namespace dynamics
} // namespace observables
} // namespace halmd
