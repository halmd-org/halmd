/*
 * Copyright © 2019 Roya Ebrahimi Viand
 * Copyright © 2014-2015 Nicolas Höft
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

#include <halmd/algorithm/gpu/copy_if.hpp>
#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/mdsim/geometries/cuboid.hpp>
#include <halmd/mdsim/geometries/cylinder.hpp>
#include <halmd/mdsim/geometries/sphere.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_groups/region_species.hpp>
#include <halmd/mdsim/gpu/particle_groups/region_species_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <exception>
#include <stdexcept>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_groups {

template <int dimension, typename float_type, typename geometry_type>
region_species<dimension, float_type, geometry_type>::region_species(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<geometry_type const> geometry
  , geometry_selection geometry_sel
  , unsigned int species
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , geometry_(geometry)
  , geometry_selection_(geometry_sel)
  , species_(species)
  , logger_(logger)
{
    geometry_->log(logger_);
    try {
        auto mask = make_cache_mutable(mask_);
        auto selection = make_cache_mutable(selection_);
        mask->reserve(particle_->dim().threads());
        mask->resize(particle_->nparticle());
        selection->reserve(particle_->dim().threads());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate global device memory");
        throw;
    }
}
template <int dimension, typename float_type, typename geometry_type>
cache<typename region_species<dimension, float_type, geometry_type>::array_type> const&
region_species<dimension, float_type, geometry_type>::selection()
{
    update_selection_();
    return selection_;
}

template <int dimension, typename float_type, typename geometry_type>
cache<typename region_species<dimension, float_type, geometry_type>::array_type> const&
region_species<dimension, float_type, geometry_type>::mask()
{
    update_mask_();
    return mask_;
}

/**
 * update the particle lists for the region
 */
template <int dimension, typename float_type, typename geometry_type>
void region_species<dimension, float_type, geometry_type>::update_mask_()
{
    cache<position_array_type> const& position_cache = particle_->position();
    if (position_cache != mask_cache_) {
        auto const& position = read_cache(position_cache);

        LOG_DEBUG("update selection mask for region and species");
        scoped_timer_type timer(runtime_.update_mask);

        auto mask = make_cache_mutable(mask_);

        auto& kernel = region_species_wrapper<dimension, geometry_type>::kernel;
        // calculate "bin", ie. inside/outside the region
        cuda::memset((*mask).begin(), (*mask).end(), 0xFF);
        kernel.compute_mask.configure(particle_->dim().grid, particle_->dim().block);
        kernel.compute_mask(
            position.data()
          , particle_->nparticle()
          , mask->data()
          , *geometry_
          , geometry_selection_ == excluded ? particle_groups::excluded : particle_groups::included
          , species_
        );
        mask_cache_ = position_cache;
    }
}

/**
 * update the particle lists for the region
 */
template <int dimension, typename float_type, typename geometry_type>
void region_species<dimension, float_type, geometry_type>::update_selection_()
{
    auto const& position_cache = particle_->position();

    if(position_cache != selection_cache_) {
        auto const& position = read_cache(position_cache);
        unsigned int nparticle = particle_->nparticle();

        LOG_DEBUG("update particle selection for region and species");
        scoped_timer_type timer(runtime_.update_selection);

        auto selection = make_cache_mutable(selection_);
        auto size = make_cache_mutable(size_);

        // allocate memory for maximum selection (i.e., all particles)
        selection->resize(nparticle);

        auto const& kernel = region_species_wrapper<dimension, geometry_type>::kernel;
        unsigned int s = kernel.copy_selection(
            position.data()
          , nparticle
          , selection->data()
          , *geometry_
          , geometry_selection_ == excluded ? particle_groups::excluded : particle_groups::included
          , species_
        );
        // shrink array to actual size of selection
        selection->resize(s);
        *size = s;

        selection_cache_ = position_cache;
    }
}

template <int dimension, typename float_type, typename geometry_type>
cache<typename region_species<dimension, float_type, geometry_type>::array_type> const&
region_species<dimension, float_type, geometry_type>::ordered()  // ID order
{
    auto const& s = selection();
    if (s != ordered_cache_) {
//        auto const& id = read_cache(particle_->id());

        LOG_WARNING("sorting selection of particle indices not yet implemented");
//        LOG_DEBUG("sorting selection of particle indices");
        scoped_timer_type timer(runtime_.sort_selection);

        auto ordered = make_cache_mutable(ordered_);
        ordered->resize(s->size());
        cuda::copy(s->begin(), s->end(), ordered->begin());

        // TODO: bring selection in ID order, sort 'ordered' by key 'id'
        // 1) copy IDs of the selection, 2) perform in-place sort
        // radix_sort(key.begin(), key.end(), ordered->begin());

        ordered_cache_ = s;
    }
    return ordered_;
}

template <int dimension, typename float_type, typename geometry_type>
cache<typename region_species<dimension, float_type, geometry_type>::array_type> const&
region_species<dimension, float_type, geometry_type>::unordered()  // memory order
{
    return selection();
}

template <int dimension, typename float_type, typename geometry_type>
cache<typename region_species<dimension, float_type, geometry_type>::size_type> const&
region_species<dimension, float_type, geometry_type>::size()
{
    update_selection_();
    return size_;
}

template <typename particle_group_type, typename particle_type>
static void wrap_to_particle(
    std::shared_ptr<particle_group_type> self
  , std::shared_ptr<particle_type> particle_src
  , std::shared_ptr<particle_type> particle_dst
)
{
    particle_group_to_particle(*particle_src, *self, *particle_dst);
}

template <int dimension, typename float_type, typename geometry_type>
void region_species<dimension, float_type, geometry_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("particle_groups")
            [
                class_<region_species, particle_group>()
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("update_mask", &runtime::update_mask)
                            .def_readonly("update_selection", &runtime::update_selection)
                            .def_readonly("sort_selection", &runtime::sort_selection)
                    ]
                    .def_readonly("runtime", &region_species::runtime_)
                    .def("to_particle", &wrap_to_particle<region_species<dimension, float_type, geometry_type>, particle_type>)

              , def("region_species", &std::make_shared<region_species<dimension, float_type, geometry_type>
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<geometry_type const>
                  , geometry_selection
                  , unsigned int
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_particle_groups_region_species(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    region_species<3, float, mdsim::geometries::cuboid<3, float>>::luaopen(L);
    region_species<2, float, mdsim::geometries::cuboid<2, float>>::luaopen(L);
    region_species<3, float, mdsim::geometries::cylinder<3, float>>::luaopen(L);
    region_species<2, float, mdsim::geometries::cylinder<2, float>>::luaopen(L);
    region_species<3, float, mdsim::geometries::sphere<3, float>>::luaopen(L);
    region_species<2, float, mdsim::geometries::sphere<2, float>>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    region_species<3, dsfloat, mdsim::geometries::cuboid<3, float>>::luaopen(L);
    region_species<2, dsfloat, mdsim::geometries::cuboid<2, float>>::luaopen(L);
    region_species<3, dsfloat, mdsim::geometries::cylinder<3, float>>::luaopen(L);
    region_species<2, dsfloat, mdsim::geometries::cylinder<2, float>>::luaopen(L);
    region_species<3, dsfloat, mdsim::geometries::sphere<3, float>>::luaopen(L);
    region_species<2, dsfloat, mdsim::geometries::sphere<2, float>>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class region_species<3, float, mdsim::geometries::cuboid<3, float>>;
template class region_species<2, float, mdsim::geometries::cuboid<2, float>>;
template class region_species<3, float, mdsim::geometries::cylinder<3, float>>;
template class region_species<2, float, mdsim::geometries::cylinder<2, float>>;
template class region_species<3, float, mdsim::geometries::sphere<3, float>>;
template class region_species<2, float, mdsim::geometries::sphere<2, float>>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class region_species<3, dsfloat, mdsim::geometries::cuboid<3, float>>;
template class region_species<2, dsfloat, mdsim::geometries::cuboid<2, float>>;
template class region_species<3, dsfloat, mdsim::geometries::cylinder<3, float>>;
template class region_species<2, dsfloat, mdsim::geometries::cylinder<2, float>>;
template class region_species<3, dsfloat, mdsim::geometries::sphere<3, float>>;
template class region_species<2, dsfloat, mdsim::geometries::sphere<2, float>>;
#endif

} // namespace particle_groups
} // namespace gpu
} // namespace mdsim
} // namespace halmd
