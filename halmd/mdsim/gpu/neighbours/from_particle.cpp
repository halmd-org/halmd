/*
 * Copyright © 2008-2011  Felix Höfling
 * Copyright © 2014       Nicolas Höft
 * Copyright © 2008-2011  Peter Colberg
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

#include <halmd/mdsim/gpu/neighbours/from_particle.hpp>
#include <halmd/mdsim/gpu/neighbours/from_particle_kernel.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {

/**
 * construct neighbour list module
 *
 * @param particle mdsim::gpu::particle instance
 * @param box mdsim::box instance
 * @param cutoff force cutoff radius
 * @param skin neighbour list skin
 * @param cell_occupancy desired average cell occupancy
 */
template <int dimension, typename float_type>
from_particle<dimension, float_type>::from_particle(
    std::shared_ptr<particle_type const> particle1
  , std::shared_ptr<particle_type const> particle2
  , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>> displacement
  , std::shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , double skin
  , double cell_occupancy
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle1_(particle1)
  , particle2_(particle2)
  , displacement1_(displacement.first)
  , displacement2_(displacement.second)
  , box_(box)
  , logger_(logger)
  // allocate parameters
  , r_skin_(skin)
  , r_cut_max_(*std::max_element(r_cut.data().begin(), r_cut.data().end()))
  , rr_cut_skin_(particle1_->nspecies(), particle2_->nspecies())
  , g_rr_cut_skin_(rr_cut_skin_.data().size())
  , nu_cell_(cell_occupancy) // FIXME neighbour list occupancy
{
    for (size_t i = 0; i < particle1_->nspecies(); ++i) {
        for (size_t j = 0; j < particle2_->nspecies(); ++j) {
            rr_cut_skin_(i, j) = std::pow(r_cut(i, j) + r_skin_, 2);
        }
    }
    cuda::copy(rr_cut_skin_.data().begin(), rr_cut_skin_.data().end(),
        g_rr_cut_skin_.begin());
    LOG("neighbour list skin: " << r_skin_);

    set_occupancy(cell_occupancy);
}

template <int dimension, typename float_type>
void from_particle<dimension, float_type>::set_occupancy(double cell_occupancy)
{
    nu_cell_ = cell_occupancy;
    // volume of n-dimensional sphere with neighbour list radius
    // volume of unit sphere: V_d = π^(d/2) / Γ(1+d/2), Γ(1) = 1, Γ(1/2) = √π
    float unit_sphere[5] = {0, 2, M_PI, 4 * M_PI / 3, M_PI * M_PI / 2 };
    assert(dimension <= 4);
    float neighbour_sphere = unit_sphere[dimension] * std::pow(r_cut_max_ + r_skin_, dimension);
    // partial number density
    float density = particle1_->nparticle() / box_->volume();
    // number of placeholders per neighbour list
    size_ = static_cast<size_t>(ceil(neighbour_sphere * (density / nu_cell_)));
    // at least cell_size (or warp_size?) placeholders
    // FIXME what is a sensible lower bound?
    // size_ = max(size_, binning_->cell_size());
    // number of neighbour lists
    stride_ = particle1_->dim().threads();
    // allocate neighbour lists
    auto g_neighbour = make_cache_mutable(g_neighbour_);
    g_neighbour->resize(stride_ * size_);
    LOG("number of placeholders per neighbour list: " << size_);
}

template <int dimension, typename float_type>
cache<typename from_particle<dimension, float_type>::array_type> const&
from_particle<dimension, float_type>::g_neighbour()
{
    cache<reverse_id_array_type> const& reverse_id_cache1 = particle1_->reverse_id();
    cache<reverse_id_array_type> const& reverse_id_cache2 = particle2_->reverse_id();

    auto current_cache = std::tie(reverse_id_cache1, reverse_id_cache2);

    if (neighbour_cache_ != current_cache|| float(displacement1_->compute()) > float(r_skin_ / 2)
        || float(displacement2_->compute()) > float(r_skin_ / 2)) {
        on_prepend_update_();
        update();
        displacement1_->zero();
        displacement2_->zero();
        neighbour_cache_ = current_cache;
        on_append_update_();
    }
    return g_neighbour_;
}

/**
 * Update neighbour lists
 */
template <int dimension, typename float_type>
void from_particle<dimension, float_type>::update()
{
    position_array_type const& position1 = read_cache(particle1_->position());
    position_array_type const& position2 = read_cache(particle2_->position());
    auto g_neighbour = make_cache_mutable(g_neighbour_);

    LOG_TRACE("update neighbour lists");


    bool overcrowded = false;
    do {
        scoped_timer_type timer(runtime_.update);

        // mark neighbour list placeholders as virtual particles
        cuda::memset(g_neighbour->begin(), g_neighbour->end(), 0xFF);
        // build neighbour lists
        cuda::vector<int> g_overflow(1);
        cuda::host::vector<int> h_overflow(1);
        cuda::memset(g_overflow.begin(), g_overflow.end(), 0);

        cuda::texture<float> rr_cut_skin(g_rr_cut_skin_);

        configure_kernel(
            get_from_particle_kernel<dimension>().update
          , particle1_->dim()
          , true
          , sizeof(unsigned int) + sizeof(vector_type)
        );
        get_from_particle_kernel<dimension>().update(
            rr_cut_skin
          , position1.data()
          , particle1_->nparticle()
          , position2.data()
          , particle2_->nparticle()
          , particle1_->nspecies()
          , particle2_->nspecies()
          , static_cast<vector_type>(box_->length())
          , &*g_neighbour->begin()
          , size_
          , stride_
          , g_overflow
        );

        cuda::copy(g_overflow.begin(), g_overflow.end(), h_overflow.begin());
        cuda::thread::synchronize();

        overcrowded = h_overflow.front() > 0;
        if (overcrowded) {
            LOG("failed to bin " << h_overflow.front() << " particles, reducing occupancy");
            set_occupancy(nu_cell_ / 2);
        }
    } while (overcrowded);
}

template <int dimension, typename float_type>
float from_particle<dimension, float_type>::defaults::occupancy() {
    return 0.4;
}

template<typename float_type>
struct variant_name;

template<>
struct variant_name<float>
{
    static constexpr const char *name = "float";
};

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template<>
struct variant_name<dsfloat>
{
    static constexpr const char *name = "dsfloat";
};
#endif

template <int dimension, typename float_type>
void from_particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    std::string const defaults_name("defaults_" + std::to_string(dimension) + "_" + std::string(variant_name<float_type>::name));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("neighbours")
            [
                class_<from_particle, _Base>()
                    .property("r_skin", &from_particle::r_skin)
                    .property("cell_occupancy", &from_particle::cell_occupancy)
                    .def("on_prepend_update", &from_particle::on_prepend_update)
                    .def("on_append_update", &from_particle::on_append_update)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("update", &runtime::update)
                    ]
                    .def_readonly("runtime", &from_particle::runtime_)
              , def("from_particle", &std::make_shared<from_particle
                    , std::shared_ptr<particle_type const>
                    , std::shared_ptr<particle_type const>
                    , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>>
                    , std::shared_ptr<box_type const>
                    , matrix_type const&
                    , double
                    , double
                    , std::shared_ptr<logger>
                  >)
            ]
          , namespace_(defaults_name.c_str())
            [
                namespace_("from_particle")
                [
                    def("occupancy", &defaults::occupancy)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_neighbours_from_particle(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    from_particle<3, float>::luaopen(L);
    from_particle<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    from_particle<3, dsfloat>::luaopen(L);
    from_particle<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class from_particle<3, float>;
template class from_particle<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class from_particle<3, dsfloat>;
template class from_particle<2, dsfloat>;
#endif

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd
