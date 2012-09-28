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

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <functional>

#include <halmd/observables/gpu/profiles.hpp>
#include <halmd/observables/gpu/profiles_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
profiles<dimension, float_type>::profiles(
    shared_ptr<particle_type const> particle
  , shared_ptr<box_type const> box
  , shared_ptr<clock_type const> clock
  , shared_ptr<force_type const> force
  , fixed_vector<unsigned, dimension> const& ngrid
  , shared_ptr<logger_type> logger
)
  // dependency injection
  : _Base(box, clock, ngrid, logger)
  , particle_(particle)
  , force_(force)
  , logger_(logger)
  // member initialisation
  , sort_(particle_->nbox, particle_->dim.threads_per_block())
  , dim_bins_(
        accumulate(ngrid_.begin(), ngrid_.end()-1, 1, multiplies<unsigned>())
      , ngrid_.back()
    )   //< (blocks, threads per block)
  // memory allocation
  , g_bins_(particle_->nbox)
  , g_permutation_(particle_->nbox)
  , g_boundaries_(dim_bins_.threads())
  , h_boundaries_(dim_bins_.threads())
  , g_stress_diag_(dim_bins_.threads())
  , h_stress_diag_(dim_bins_.threads())
{
    // initialise device constants
    try {
        cuda::copy(particle_->nbox, get_profiles_kernel<dimension>().nbox);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("failed to copy cell parameters to device symbols");
        throw;
    }
}

/**
 * compute all profiles
 */
template <int dimension, typename float_type>
void profiles<dimension, float_type>::compute_profiles()
{
    scoped_timer_type timer_(runtime_.sample);

    unsigned nbins = accumulate(ngrid_.begin(), ngrid_.end(), 1, multiplies<unsigned>());

    // compute bin ID for each particle
    // and store serial particle index in g_permutation_
    {
        typedef typename profiles_wrapper<dimension>::vector_type gpu_vector_type;

        LOG_TRACE("compute profile bins");
        scoped_timer<timer> timer_(runtime_.bins);
        cuda::configure(particle_->dim.grid, particle_->dim.block);
        get_profiles_kernel<dimension>().compute_bins(
            particle_->g_r, g_bins_, g_permutation_
          , ngrid_
          , static_cast<gpu_vector_type>(spacing_)
          , static_cast<gpu_vector_type>(box_->origin())
        );
    }

    // sort bin IDs and generate permutation of particle indices
    {
        LOG_TRACE("sort bin IDs");
        scoped_timer<timer> timer_(runtime_.sort);
        cuda::configure(particle_->dim.grid, particle_->dim.block);
        sort_(g_bins_, g_permutation_);
    }

    // find bin boundaries in sorted list of bin IDs
    {
        LOG_TRACE("find boundaries");
        scoped_timer<timer> timer_(runtime_.boundaries);
        cuda::memset(g_boundaries_, 0xFF);                          //< initialise with -1
        cuda::configure(particle_->dim.grid, particle_->dim.block);
        get_profiles_kernel<dimension>().find_boundaries(g_bins_, g_boundaries_);
    }

    // collect stress tensor for each bin
    {
        LOG_TRACE("collect stress tensors");
        scoped_timer<timer> timer_(runtime_.stress_tensor);
        cuda::configure(dim_bins_.grid, dim_bins_.block);
        get_profiles_kernel<dimension>().collect_stress_tensor(
            g_boundaries_, g_permutation_
          , get<0>(force_->stress_tensor_pot())
          , particle_->g_v
          , g_stress_diag_
          , nbins
        );
    }

    // copy profiles
    try {
        LOG_TRACE("copy profiles from device to host");
        scoped_timer<timer> timer_(runtime_.copy);
        cuda::copy(g_boundaries_, h_boundaries_);
        cuda::copy(g_stress_diag_, h_stress_diag_);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("failed to copy profiles from GPU to host");
        throw;
    }

    // compute number of particles per slab
    // TODO move to CUDA kernel
    // access permuations & stress tensor via texture fetches?
    {
        LOG_TRACE("reduce multi-dimensional histograms");
        scoped_timer<timer> timer_(runtime_.reduce);
        // initialise accumulators
        BOOST_FOREACH (std::vector<double>& profile, density_profile_) {
            fill(profile.begin(), profile.end(), 0);
        }
        BOOST_FOREACH (std::vector<vector_type>& profile, stress_tensor_profile_) {
            fill(profile.begin(), profile.end(), 0);
        }
        for (unsigned id = 0; id < nbins; ++id) {
            // number of particles in bin #id
            unsigned N = h_boundaries_[id].y - h_boundaries_[id].x;
            // total stress tensor in bin #id
            vector_type stress_diag = static_cast<vector_type>(
                static_cast<gpu_vector_type>(h_stress_diag_[id])
            );
            // compute index from bin ID by inverting ID = x + ngrid_[0] * (y + ngrid_[1] * z)
            // and assign to profile bins
            unsigned id_ = id;
            for (unsigned axis = 0; axis < dimension; ++axis) {
                unsigned n = ngrid_[axis];
                density_profile_[axis][id_ % n] += N;
                stress_tensor_profile_[axis][id_ % n] += stress_diag;
                id_ /= n;
            }
        }
    }
#ifndef NDEBUG
    // check that nobody has got lost
    for (unsigned axis = 0; axis < dimension; ++axis) {
        assert(std::accumulate(
                density_profile_[axis].begin()
              , density_profile_[axis].end()
              , 0., plus<double>()
            ) == particle_->nbox);
    }
#endif

    // normalisation of density and stress profiles by slab volume
    LOG_TRACE("normalisation");
    double box_volume = box_->volume();
    for (unsigned axis = 0; axis < dimension; ++axis) {
        double slab_vol = box_volume / box_->length()[axis] * spacing_[axis];
        for (unsigned index = 0; index < ngrid_[axis]; ++index) {
            stress_tensor_profile_[axis][index] /= slab_vol;
            density_profile_[axis][index] /= slab_vol;
        }
    }
}

template <int dimension, typename float_type>
void profiles<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("profiles_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                class_<profiles, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .def(constructor<
                        shared_ptr<particle_type const>
                      , shared_ptr<box_type const>
                      , shared_ptr<clock_type const>
                      , shared_ptr<force_type const>
                      , fixed_vector<unsigned, dimension> const&
                      , shared_ptr<logger_type>
                    >())
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("sample", &runtime::sample)
                            .def_readonly("bins", &runtime::bins)
                            .def_readonly("sort", &runtime::sort)
                            .def_readonly("boundaries", &runtime::boundaries)
                            .def_readonly("stress_tensor", &runtime::stress_tensor)
                            .def_readonly("reduce", &runtime::reduce)
                            .def_readonly("copy", &runtime::copy)
                    ]
                    .def_readonly("runtime", &profiles::runtime_)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_profiles(lua_State* L)
{
    profiles<3, float>::luaopen(L);
    profiles<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class profiles<3, float>;
template class profiles<2, float>;

} // namespace gpu
} // namespace observables
} // namespace halmd
