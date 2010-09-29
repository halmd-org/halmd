/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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
#include <boost/array.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <cmath>
#include <limits>
#include <numeric>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/position/lattice_kernel.hpp>
#include <halmd/mdsim/gpu/position/lattice.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu { namespace position
{

template <int dimension, typename float_type, typename RandomNumberGenerator>
lattice<dimension, float_type, RandomNumberGenerator>::lattice(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<random_type> random
)
  // dependency injection
  : particle(particle)
  , box(box)
  , random(random)
{
}

/**
 * Place particles on a face-centered cubic (fcc) lattice
 *
 * The task is to determine the minimum lattice distance of an fcc lattice
 * that fits into a rectangular parallelepiped (the simulation box) and has
 * equally many or more lattice sites than the number of particles.
 *
 * The number of lattice unit cells is defined as
 *
 *    s = floor(L_x / a) floor(L_y / a) floor(L_z / a) ≥ ceil(N / u)
 *
 *    N       number of particles
 *    L_i     box edge lengths for i ∈ {x, y, z}
 *    a       lattice distance
 *    u       number of particle per unit cell (4 in 3D, 2 in 2D)
 *
 * The upper bound for the lattice distance is given by
 *
 *    a ≤ (L_x L_y L_z / ceil(N / u))^(1/3)
 *
 * This yields lower bounds for the number of unit cells per dimension
 *
 *    n_i ≥ floor(L_i / a)
 *
 * The minimum a is then determined by iteratively increasing any of the
 * n_i to yield the nearest smaller value for the lattice distance until
 *
 *    s * u ≥ N
 *
 * is satisfied.
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void lattice<dimension, float_type, RandomNumberGenerator>::set()
{
    LOG("randomly permuting particle types");
    random->shuffle(particle->g_r);

    // TODO: handle non-cubic boxes
    for (unsigned i=1; i < dimension; i++) {
        assert(box->length()[0] == box->length()[i]);
    }

    LOG("placing particles on face-centered cubic (fcc) lattice");

    // particles per 2- or 3-dimensional unit cell
    unsigned int const m = 2 * (dimension - 1);
    // lower boundary for number of particles per lattice dimension
    unsigned int n = static_cast<unsigned int>(pow(particle->nbox / m, 1. / dimension));
    // lower boundary for total number of lattice sites
    unsigned int N = m * static_cast<unsigned int>(pow(static_cast<double>(n), dimension));

    if (N < particle->nbox) {
        n += 1;
        N = m * static_cast<unsigned int>(pow(static_cast<double>(n), dimension));
    }
    if (N > particle->nbox) {
        LOG_WARNING("lattice not fully occupied (" << N << " sites)");
    }

    // minimum distance in 2- or 3-dimensional fcc lattice
    LOG("minimum lattice distance: " << (box->length()[0] / n) / sqrt(2));

#ifdef USE_VERLET_DSFUN
    cuda::memset(particle->g_r, 0, particle->g_r.capacity());
#endif

//     boost::array<high_resolution_timer, 2> timer;
    cuda::thread::synchronize();
    try {
//         timer[0].record();
        cuda::configure(particle->dim.grid, particle->dim.block);
        get_lattice_kernel<dimension>().fcc(particle->g_r, n, box->length()[0]);
        cuda::thread::synchronize();
//         timer[1].record();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw runtime_error("failed to compute particle lattice positions on GPU");
    }
//     m_times["lattice"] += timer[1] - timer[0];

    // ??? shift particle positions to range (-L/2, L/2)
//    box->reduce_periodic(r);

//    ??? assign_positions();

    // reset particle image vectors
    cuda::memset(particle->g_image, 0, particle->g_image.capacity());
}

template <typename T>
static void register_lua(char const* class_name)
{
    typedef typename T::_Base _Base;
    typedef typename T::particle_type particle_type;
    typedef typename T::box_type box_type;
    typedef typename T::random_type random_type;

    using namespace luabind;
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("gpu")
                [
                    namespace_("position")
                    [
                        class_<T, shared_ptr<_Base>, bases<_Base> >(class_name)
                            .def(constructor<shared_ptr<particle_type>, shared_ptr<box_type>, shared_ptr<random_type> >())
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    register_lua<lattice<3, float, random::gpu::rand48> >("lattice_3_");
    register_lua<lattice<2, float, random::gpu::rand48> >("lattice_2_");
}

// explicit instantiation
template class lattice<3, float, random::gpu::rand48>;
template class lattice<2, float, random::gpu::rand48>;

}}} // namespace mdsim::gpu::position

} // namespace halmd
