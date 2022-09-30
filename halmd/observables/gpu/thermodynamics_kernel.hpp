/*
 * Copyright © 2020  Roya Ebrahimi Viand
 * Copyright © 2013  Nicolas Höft
 * Copyright © 2012  Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_THERMODYNAMICS_KERNEL_HPP
#define HALMD_OBSERVABLES_GPU_THERMODYNAMICS_KERNEL_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/iterator.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>


#ifndef __CUDACC__
# include <tuple>
#endif

namespace halmd {
namespace observables {
namespace gpu {

/**
 * Compute total kinetic energy.
 */
template <int dimension, typename float_type>
class kinetic_energy
{
private:
    typedef unsigned int size_type;

public:
    /** element pointer type of input array */
    typedef size_type const* iterator;

    /**
     * Initialise kinetic energy to zero.
     */
    kinetic_energy(cudaTextureObject_t texture) : mv2_(0), texture_(texture) {}

    /**
     * Accumulate kinetic energy of a particle.
     */
    inline HALMD_GPU_ENABLED void operator()(size_type i);

    /**
     * Accumulate kinetic energy of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(kinetic_energy const& acc)
    {
        mv2_ += acc.mv2_;
    }

    /**
     * Returns total kinetic energy.
     */
    HALMD_GPU_ENABLED double operator()() const
    {
        return 0.5 * mv2_;
    }

private:
    /** sum over mass × square of velocity vector */
    float_type mv2_;
    /** texture with velocities and masses */
    cudaTextureObject_t texture_;
};

/**
 * Compute centre of mass.
 */
template <int dimension, typename float_type>
class centre_of_mass
{
private:
    typedef unsigned int size_type;
    typedef fixed_vector<float_type, dimension> vector_type;

public:
    /** element pointer type of input array */
    typedef zip_iterator<size_type const*, constant_iterator<fixed_vector<float, dimension> > > iterator;
    typedef typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;

    /**
     * Initialise momentum and mass to zero.
     */
    centre_of_mass(cudaTextureObject_t t_position, cudaTextureObject_t t_image
      , cudaTextureObject_t t_velocity) :
        mr_(0), m_(0), t_position_(t_position), t_image_(t_image)
      , t_velocity_(t_velocity) {}

    /**
     * Accumulate momentum and mass of a particle.
     */
    inline HALMD_GPU_ENABLED void operator()(typename iterator::value_type const& value);

    /**
     * Accumulate centre of mass of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(centre_of_mass const& acc)
    {
        mr_ += acc.mr_;
        m_ += acc.m_;
    }

    /**
     * Returns centre of mass.
     */
    HALMD_GPU_ENABLED fixed_vector<double, dimension> operator()() const
    {
        return fixed_vector<double, dimension>(mr_ / m_);
    }

private:
    /** sum over momentum vector */
    vector_type mr_;
    /** sum over mass */
    float_type m_;
    /** texture with positions and species */
    cudaTextureObject_t t_position_;
    /** texture with images */
    cudaTextureObject_t t_image_;
    /** texture with velocities and masses */
    cudaTextureObject_t t_velocity_;
};

/**
 * Compute velocity of centre of mass.
 */
template <int dimension, typename float_type>
class velocity_of_centre_of_mass
{
private:
    typedef unsigned int size_type;
    typedef fixed_vector<float_type, dimension> vector_type;

public:
    /** element pointer type of input array */
    typedef size_type const* iterator;

    /**
     * Initialise momentum and mass to zero.
     */
    velocity_of_centre_of_mass(cudaTextureObject_t texture) : mv_(0), m_(0), texture_(texture) {}

    /**
     * Accumulate momentum and mass of a particle.
     */
    inline HALMD_GPU_ENABLED void operator()(size_type i);

    /**
     * Accumulate velocity centre of mass of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(velocity_of_centre_of_mass const& acc)
    {
        mv_ += acc.mv_;
        m_ += acc.m_;
    }

#ifndef __CUDACC__
    /**
     * Returns centre of mass momentum and total mass.
     */
    HALMD_GPU_ENABLED std::tuple<fixed_vector<double, dimension>, double> operator()() const
    {
        return std::make_tuple(fixed_vector<double, dimension>(mv_), double(m_));
    }
#endif

private:
    /** sum over momentum vector */
    vector_type mv_;
    /** sum over mass */
    float_type m_;
    /** texture with velocities and masses */
    cudaTextureObject_t texture_;
};

/**
 * Compute total force.
 */
template <int dimension, typename float_type>
class total_force
{
private:
    typedef unsigned int size_type;
    typedef fixed_vector<float_type, dimension> vector_type;

public:
    /** element pointer type of input array */
    typedef size_type const* iterator;
    typedef typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type gpu_force_type;

    /**
     * Initialise force to zero.
     */
    total_force(cudaTextureObject_t texture) : force_(0), texture_(texture) {}

    /**
     * Add force of particle i.
     */
    inline HALMD_GPU_ENABLED void operator()(size_type i);

    /**
     * Add forces of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(total_force const& acc)
    {
        force_ += acc.force_;
    }

#ifndef __CUDACC__
    /**
     * Returns total force.
     */
    HALMD_GPU_ENABLED fixed_vector<double, dimension> operator()() const
    {
        return fixed_vector<double, dimension>(force_);
    }
#endif

private:
    /** total force */
    vector_type force_;
    /** texture with forces */
    cudaTextureObject_t texture_;
};

/**
 * Compute total potential energy.
 */
template <typename float_type>
class potential_energy
{
private:
    typedef unsigned int size_type;

public:
    typedef size_type const* iterator;

    /**
     * Initialise potential energy to zero
     */
    potential_energy(cudaTextureObject_t texture) : en_pot_(0), texture_(texture) {}

    /**
     * Accumulate potential energy of a particle.
     */
    inline HALMD_GPU_ENABLED void operator()(size_type i);

    /**
     * Accumulate potential energy of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(potential_energy const& acc)
    {
        en_pot_ += acc.en_pot_;
    }

    /**
     * Returns total potential energy.
     */
    HALMD_GPU_ENABLED double operator()() const
    {
        return en_pot_;
    }

private:
    /** total potential energy */
    float_type en_pot_;
    /** texture with potential energies */
    cudaTextureObject_t texture_;
};

/**
 * Compute total virial sum.
 */
template <int dimension, typename float_type>
class virial
{
private:
    typedef unsigned int size_type;

public:
    /** element pointer type of input array */
    typedef size_type const* iterator;

    /**
     * Initialise virial sum to zero and store number of strides
     */
    virial(cudaTextureObject_t texture, unsigned int stride) : virial_(0), stride_(stride), texture_(texture) {}

    /**
     * Accumulate stress tensor diagonal of a particle.
     */
    inline HALMD_GPU_ENABLED void operator()(size_type i);

    /**
     * Accumulate virial sum of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(virial const& acc)
    {
        virial_ += acc.virial_;
    }

    /**
     * Returns total virial sum.
     */
    HALMD_GPU_ENABLED double operator()() const
    {
        return virial_;
    }

private:
    /** total virial sum */
    float_type virial_;
    /**
     * stride of the stress tensor array in device memory
     *
     * Note that the stride is defined by GTDIM in the force kernel, which may
     * be different in the reduce kernel. Thus we pass the value inside the
     * accumulation functor.
     **/
    unsigned int stride_;
    /** texture with stress tensors */
    cudaTextureObject_t texture_;
};

/**
 * Compute (full) stress tensor sum from potential part of the stress tensor.
 */
template <int dimension, typename float_type>
class stress_tensor
{
private:
    typedef unsigned int size_type;
    typedef typename mdsim::type_traits<dimension, float_type>::stress_tensor_type stress_tensor_type;

public:
    /** element pointer type of input array */
    typedef size_type const* iterator;

    /**
     * Initialise stress tensor sum to zero and store number of strides
     */
    stress_tensor(cudaTextureObject_t t_velocity, cudaTextureObject_t t_stress_pot, unsigned int stride) :
        stride_(stride), stress_tensor_(0), t_velocity_(t_velocity), t_stress_pot_(t_stress_pot) {}

    /**
     * Accumulate stress tensor diagonal of a particle.
     */
    inline HALMD_GPU_ENABLED void operator()(size_type i);

    /**
     * Accumulate stress tensor sum of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(stress_tensor const& acc)
    {
        stress_tensor_ += acc.stress_tensor_;
    }

    /**
     * Returns total stress tensor sum.
     */
    HALMD_GPU_ENABLED stress_tensor_type operator()() const
    {
        return stress_tensor_;
    }

private:
    /**
     * stride of the stress tensor array in device memory
     *
     * Note that the stride is defined by GTDIM in the force kernel, which may
     * be different in the reduce kernel. Thus we pass the value inside the
     * accumulation functor.
     **/
    unsigned int stride_;
    /** sum of stress tensors */
    stress_tensor_type stress_tensor_;
    /** texture with velocities */
    cudaTextureObject_t t_velocity_;
    /** texture with stress tensors */
    cudaTextureObject_t t_stress_pot_;
};

/**
 * Compute mean heat flux.
 */
template <int dimension, typename float_type>
class heat_flux
{
private:
    typedef unsigned int size_type;
    typedef fixed_vector<float_type, dimension> vector_type;

public:
    /** element pointer type of input array */
    typedef size_type const* iterator;
    /**
     * Initialise heat flux to zero.
     */
    heat_flux(unsigned int stride) : hf_(0), stride_(stride) {}
    /**
     * Accumulate heat flux of a particle.
     */
    inline HALMD_GPU_ENABLED void operator()(size_type i);

    /**
     * Accumulate heat flux of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(heat_flux const& acc)
    {
        hf_ += acc.hf_;
    }

    /**
     * Returns total heat flux.
     */
    HALMD_GPU_ENABLED fixed_vector<double, dimension> operator()() const
    {
        return fixed_vector<double, dimension>(hf_);
    }

    /**
     * Returns reference to texture with velocities and masses.
     */
    static cuda::texture<float4> const& get()
    {
        return texture_;
    }
    /**
     * Returns reference to texture with potential part of stress tensors.
     */
    static cuda::texture<float> const& get_stress_pot()
    {
        return stress_pot_texture_;
    }

    /**
     * Returns reference to texture with potential energy.
     */
    static cuda::texture<float> const& get_potential_energy()
    {
        return potential_texture_;
    }

private:
    /** sum over heat flux vector */
    vector_type hf_;
    /** texture with velocities and masses*/
    static cuda::texture<float4> const texture_;
    /** texture with potential energies*/
    static cuda::texture<float> const potential_texture_;
    /** texture with stress tensors */
    static cuda::texture<float> const stress_pot_texture_;
    /**
     * stride of the stress tensor array in device memory
     *
     * Note that the stride is defined by GTDIM in the force kernel, which may
     * be different in the reduce kernel. Thus we pass the value inside the
     * accumulation functor.
     **/
    unsigned int stride_;
};

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_KERNEL_HPP */
