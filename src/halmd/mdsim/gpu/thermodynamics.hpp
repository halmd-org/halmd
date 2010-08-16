/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_THERMODYNAMICS_HPP
#define HALMD_MDSIM_GPU_THERMODYNAMICS_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <vector>

#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/thermodynamics.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

template <int dimension, typename float_type>
class force;

namespace forces {
template <int dimension, typename float_type>
class lj;
}

template <int dimension, typename float_type>
class thermodynamics
    : public mdsim::thermodynamics<dimension>
{
public:
    // module definitions
    typedef thermodynamics _Self;
    typedef mdsim::thermodynamics<dimension> _Base;
    static void depends();
    static void options(po::options_description& desc) {}
    static void select(po::options const& vm) {}

    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename particle_type::gpu_vector_type gpu_vector_type;
    typedef typename _Base::virial_type virial_type;
    typedef typename boost::mpl::if_c<virial_type::static_size >= 3,
                         float4, float2>::type gpu_virial_type;

    shared_ptr<particle_type> particle;

    thermodynamics(modules::factory& factory, po::options const& vm);
    virtual ~thermodynamics() {}

    double en_kin() const;
    vector_type v_cm() const;
    double en_pot() const;
    std::vector<virial_type> const& virial();

private:
    /** virial for each particle type */
    /** the GPU implementation returns the total virial sum only at index 0 */
    std::vector<virial_type> virial_;
    /** potential energy for each particle */
    cuda::vector<float> g_en_pot_;
    /** virial for each particle */
    cuda::vector<gpu_virial_type> g_virial_;

    // FIXME can we inherit the friendship from force to its derived classes?
    friend class gpu::force<dimension, float_type>;
    friend class gpu::forces::lj<dimension, float_type>;
};


}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_THERMODYNAMICS_HPP */
