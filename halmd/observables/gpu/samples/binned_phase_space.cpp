/*
 * Copyright © 2011  Felix Höfling
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

#include <boost/make_shared.hpp>
#include <limits>

#include <halmd/io/logger.hpp>
#include <halmd/observables/gpu/samples/binned_phase_space.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace gpu {
namespace samples {

template <int dimension, typename float_type>
binned_phase_space<dimension, float_type>::binned_phase_space(
    shared_ptr<data_sample_type const> data_sample
  , unsigned int species
  , fixed_vector<unsigned int, dimension> const& nbin
  , double occupancy
  , shared_ptr<logger_type> logger
)
  // initialise public attributes
  : nbin(nbin)
  , step(numeric_limits<step_type>::max())
  // private module dependencies
  , data_sample_(data_sample)
  , species_(species)
  , logger_(logger)
{
    // compute cell size as [[#particles / #bins] / occupancy]
    // and make it a multiple of warp size
    size_t total_bins = accumulate(nbin.begin(), nbin.end(), 1u, multiplies<unsigned int>());
    size_t npart = position_data().size(); // number of particles
    cuda::device::properties device_prop(cuda::device::get());
    size_t warp_size = device_prop.warp_size();

    cell_size = ceil(((npart + total_bins - 1) / total_bins) / occupancy);
    cell_size = ((cell_size + warp_size - 1) / warp_size) * warp_size;

    // allocate memory
    try {
        g_cell.resize(total_bins * cell_size);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate binning cells in global device memory");
        throw;
    }

    LOG("number of placeholders per cell: " << cell_size);
    LOG("desired average cell occupancy: " << occupancy);
    LOG("effective average cell occupancy: " << static_cast<double>(npart) / g_cell.size());

    // as CUDA compute capability ≤ 1.3 supports only 2-dimensional block grids
    // we redistribute the number of bins per dimension
    unsigned int maximum = *max_element(nbin.begin(), nbin.end());
    dim = device::validate(cuda::config(
        dim3(total_bins / maximum, maximum)
      , min(cell_size, device_prop.max_threads_per_block())
    ));
    LOG("using " << dim.threads_per_block() << " CUDA threads in "
        << dim.grid.x << "×" << dim.grid.y << "×" << dim.grid.z << " blocks"
    );
}

template <int dimension, typename float_type>
static int wrap_dimension(binned_phase_space<dimension, float_type> const&)
{
    return dimension;
}

template <int dimension, typename float_type>
void binned_phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("binned_phase_space_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>() + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                namespace_("samples")
                [
                    class_<binned_phase_space>(class_name.c_str())
                        .property("dimension", &wrap_dimension<dimension, float_type>)
                ]
            ]

          , namespace_("samples")
            [
                def("binned_phase_space", &make_shared<binned_phase_space
                  , shared_ptr<data_sample_type const>
                  , unsigned int
                  , fixed_vector<unsigned int, dimension>
                  , double
                  , shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_samples_binned_phase_space(lua_State* L)
{
    binned_phase_space<3, float>::luaopen(L);
    binned_phase_space<2, float>::luaopen(L);
    return 0;
}

template class binned_phase_space<3, float>;
template class binned_phase_space<2, float>;

} // namespace samples
} // namespace gpu
} // namespace observables
} // namespace halmd
