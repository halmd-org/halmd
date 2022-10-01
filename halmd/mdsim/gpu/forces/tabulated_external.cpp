/*
 * Copyright © 2013 Nicolas Höft
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

#include <halmd/mdsim/gpu/forces/tabulated_external.hpp>
#include <halmd/mdsim/forces/interpolation/cubic_hermite.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <luaponte/out_value_policy.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

using namespace halmd::mdsim::forces::interpolation;

template <int dimension, typename float_type, typename force_interpolation_type>
tabulated_external<dimension, float_type, force_interpolation_type>::tabulated_external(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<force_interpolation_type> force_interpolation
  , std::shared_ptr<virial_interpolation_type> virial_interpolation
  , std::shared_ptr<logger_type> logger
)
  : particle_(particle)
  , box_(box)
  , force_interpolation_(force_interpolation)
  , virial_interpolation_(virial_interpolation)
  , logger_(logger)
  , force_coefficients_(force_interpolation_->total_knots() * force_interpolation_type::coefficients_per_knot)
  , virial_coefficients_(virial_interpolation_->total_knots() * virial_interpolation_type::coefficients_per_knot)
{
    cuda::memset(force_coefficients_.begin(), force_coefficients_.end(), 0);
    cuda::memset(virial_coefficients_.begin(), virial_coefficients_.end(), 0);
}

template <int dimension, typename float_type, typename force_interpolation_type>
void tabulated_external<dimension, float_type, force_interpolation_type>::check_cache()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (force_cache_ != position_cache) {
        particle_->mark_force_dirty();
    }

    if (aux_cache_ != position_cache) {
        particle_->mark_aux_dirty();
    }
}

template <int dimension, typename float_type, typename force_interpolation_type>
void tabulated_external<dimension, float_type, force_interpolation_type>::apply()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (particle_->aux_enabled()) {
        compute_aux_();
        force_cache_ = position_cache;
        aux_cache_ = force_cache_;
    }
    else {
        compute_();
        force_cache_ = position_cache;
    }
    particle_->force_zero_disable();
}


template <int dimension, typename float_type, typename force_interpolation_type>
inline void tabulated_external<dimension, float_type, force_interpolation_type>::compute_()
{
    position_array_type const& position = read_cache(particle_->position());
    auto force = make_cache_mutable(particle_->mutable_force());

    LOG_TRACE("compute forces");

    scoped_timer_type timer(runtime_.compute);

    configure_kernel(gpu_wrapper::kernel.compute, particle_->dim(), true);
    gpu_wrapper::kernel.compute(
        position.data()
      , &*force->begin()
      , nullptr
      , nullptr
      , particle_->nspecies()
      , static_cast<position_type>(box_->length())
      , &*force_coefficients_.begin()
      , nullptr
      , *force_interpolation_
      , *virial_interpolation_
      , particle_->force_zero()
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename force_interpolation_type>
inline void tabulated_external<dimension, float_type, force_interpolation_type>::compute_aux_()
{
    position_array_type const& position = read_cache(particle_->position());
    auto force = make_cache_mutable(particle_->mutable_force());
    auto en_pot = make_cache_mutable(particle_->mutable_potential_energy());
    auto stress_pot = make_cache_mutable(particle_->mutable_stress_pot());

    LOG_TRACE("compute forces with auxiliary variables");

    scoped_timer_type timer(runtime_.compute);

    configure_kernel(gpu_wrapper::kernel.compute_aux, particle_->dim(), true);
    gpu_wrapper::kernel.compute_aux(
        position.data()
      , &*force->begin()
      , &*en_pot->begin()
      , &*stress_pot->begin()
      , particle_->nspecies()
      , static_cast<position_type>(box_->length())
      , &*force_coefficients_.begin()
      , &*virial_coefficients_.begin()
      , *force_interpolation_
      , *virial_interpolation_
      , particle_->force_zero()
    );
    cuda::thread::synchronize();
}

template<typename force_type, typename iterator_type>
inline iterator_type
set_coefficients(force_type& tabulated, iterator_type const& first)
{
    typedef typename force_type::coefficient_array_type coefficient_array_type;
    typedef typename force_type::coefficient_value_type value_type;

    coefficient_array_type& g_coefficient = tabulated.coefficients();
    cuda::memory::host::vector<value_type> h_coefficient(g_coefficient.size());
    iterator_type input = first;
    for (value_type& constraint : h_coefficient) {
        constraint = *input++;
    }
    cuda::copy(h_coefficient.begin(), h_coefficient.end(), g_coefficient.begin());
    return input;
}

/**
 * Copy interpolation coefficients per particle to given array.
 */
template <typename force_type, typename iterator_type>
inline iterator_type
get_coefficients(force_type& tabulated, iterator_type const& first)
{
    typedef typename force_type::coefficient_array_type coefficient_array_type;
    coefficient_array_type const& g_coefficient = tabulated.coefficients();
    cuda::memory::host::vector<typename coefficient_array_type::value_type> h_coefficient(g_coefficient.size());
    cuda::copy(g_coefficient.begin(), g_coefficient.end(), h_coefficient.begin());
    return std::copy(h_coefficient.begin(), h_coefficient.end(), first);
}

template<typename force_type>
static void
wrap_set_coefficients(std::shared_ptr<force_type> self, std::vector<typename force_type::coefficient_value_type> const& input)
{
    if(input.size() != self->ncoefficients()) {
        throw std::invalid_argument("input array size not equal to number needed constraints");
    }
    set_coefficients(*self, input.begin());
}

template<typename force_type>
static std::function<std::vector<typename force_type::coefficient_value_type> ()>
wrap_get_coefficients(std::shared_ptr<force_type> self)
{
    return [=]() -> std::vector<typename force_type::coefficient_value_type> {
        std::vector<typename force_type::coefficient_value_type> output;
        {
            output.reserve(self->ncoefficients());
        }
        get_coefficients(*self, std::back_inserter(output));
        return output;
    };
}

template <typename force_type>
static std::function<std::vector<typename force_type::coefficient_value_type>& ()>
wrap_coefficients(std::shared_ptr<force_type> self, std::function<void ()>& array_to_sample)
{
    typedef std::vector<typename force_type::coefficient_value_type> array_type;
    std::shared_ptr<array_type> array = std::make_shared<array_type>();

    array_to_sample = [=]() {
        if (self->coefficients().size() != array->size()) {
            throw std::runtime_error("input array size not equal to number needed constraints");
        }
        set_coefficients(*self, array->begin());
        array->clear();
    };
    return [=]() -> array_type& {
        return *array;
    };
}

template<typename tabulated_type, typename iterator_type>
inline iterator_type
set_virial_coefficients(tabulated_type& tabulated, iterator_type const& first)
{
    typedef typename tabulated_type::coefficient_array_type coefficient_array_type;
    typedef typename tabulated_type::coefficient_value_type value_type;

    coefficient_array_type& g_coefficient = tabulated.virial_coefficients();
    cuda::memory::host::vector<value_type> h_coefficient(g_coefficient.size());
    iterator_type input = first;
    for (value_type& constraint : h_coefficient) {
        constraint = *input++;
    }
    cuda::copy(h_coefficient.begin(), h_coefficient.end(), g_coefficient.begin());
    return input;
}

/**
 * Copy interpolation coefficients per particle to given array.
 */
template <typename tabulated_type, typename iterator_type>
inline iterator_type
get_virial_coefficients(tabulated_type& tabulated, iterator_type const& first)
{
    typedef typename tabulated_type::coefficient_array_type coefficient_array_type;
    coefficient_array_type const& g_coefficient = tabulated.virial_coefficients();
    cuda::memory::host::vector<typename coefficient_array_type::value_type> h_coefficient(g_coefficient.size());
    cuda::copy(g_coefficient.begin(), g_coefficient.end(), h_coefficient.begin());
    return std::copy(h_coefficient.begin(), h_coefficient.end(), first);
}

template<typename tabulated_type>
static void
wrap_set_virial_coefficients(std::shared_ptr<tabulated_type> self, std::vector<typename tabulated_type::coefficient_value_type> const& input)
{
    if(input.size() != self->ncoefficients()) {
        throw std::invalid_argument("input array size not equal to number needed constraints");
    }
    set_virial_coefficients(*self, input.begin());
}


template<typename tabulated_type>
static std::function<std::vector<typename tabulated_type::coefficient_value_type> ()>
wrap_get_virial_coefficients(std::shared_ptr<tabulated_type> self)
{
    return [=]() -> std::vector<typename tabulated_type::coefficient_value_type> {
        std::vector<typename tabulated_type::coefficient_value_type> output;
        {
            output.reserve(self->ncoefficients());
        }
        get_virial_coefficients(*self, std::back_inserter(output));
        return output;
    };
}

template <typename tabulated_type>
static std::function<std::vector<typename tabulated_type::coefficient_value_type>& ()>
wrap_virial_coefficients(std::shared_ptr<tabulated_type> self, std::function<void ()>& array_to_sample)
{
    typedef std::vector<typename tabulated_type::coefficient_value_type> array_type;
    std::shared_ptr<array_type> array = std::make_shared<array_type>();

    array_to_sample = [=]() {
        if (self->virial_coefficients().size() != array->size()) {
            throw std::runtime_error("input array size not equal to number needed constraints");
        }
        set_virial_coefficients(*self, array->begin());
        array->clear();
    };
    return [=]() -> array_type& {
        return *array;
    };
}

template <int dimension, typename float_type, typename force_interpolation_type>
void tabulated_external<dimension, float_type, force_interpolation_type>::luaopen(lua_State* L)
{
    using namespace boost::placeholders;
    using namespace luaponte;

    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("forces")
            [
                class_<tabulated_external>()
                    .def("set_coefficients", &wrap_set_coefficients<tabulated_external>)
                    .property("get_coefficients", &wrap_get_coefficients<tabulated_external>)
                    .def("set_virial_coefficients", &wrap_set_virial_coefficients<tabulated_external>)
                    .property("get_virial_coefficients", &wrap_get_virial_coefficients<tabulated_external>)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                    ]
                    .def_readonly("runtime", &tabulated_external::runtime_)
                    .def("coefficients", &wrap_coefficients<tabulated_external>, pure_out_value(_2))
                    .def("virial_coefficients", &wrap_virial_coefficients<tabulated_external>, pure_out_value(_2))
                    .def("check_cache", &tabulated_external::check_cache)
                    .def("apply", &tabulated_external::apply)
              , def("tabulated_external", &std::make_shared<tabulated_external,
                    std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<force_interpolation_type>
                  , std::shared_ptr<virial_interpolation_type>
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_forces_tabulated_external(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    tabulated_external<2, float, cubic_hermite<2, float> >::luaopen(L);
    tabulated_external<3, float, cubic_hermite<3, float> >::luaopen(L);
#endif

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    tabulated_external<2, dsfloat, cubic_hermite<2, float> >::luaopen(L);
    tabulated_external<3, dsfloat, cubic_hermite<3, float> >::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class tabulated_external<2, float, cubic_hermite<2, float> >;
template class tabulated_external<3, float, cubic_hermite<3, float> >;
#endif

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class tabulated_external<2, dsfloat, cubic_hermite<2, float> >;
template class tabulated_external<3, dsfloat, cubic_hermite<3, float> >;
#endif
} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd