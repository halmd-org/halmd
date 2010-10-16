/* Time correlation functions
 *
 * Copyright © 2008-2010  Peter Colberg
 *                        Felix Höfling
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

#ifndef HALMD_SAMPLE_TCF_BASE_HPP
#define HALMD_SAMPLE_TCF_BASE_HPP

#include <algorithm>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <vector>

#include <halmd/math/accum.hpp>
#include <halmd/util/H5xx.hpp>

namespace halmd {

template <int dimension>
struct tcf_sample
{
    /** q vector in momentum space */
    typedef vector<double, dimension> vector_type;
    /** |q| value vector */
    typedef std::vector<double> q_value_vector;
    /** vector of q vectors for different |q| values */
    typedef std::vector<std::vector<vector_type> > q_vector_vector;
    /** real and imaginary components of Fourier transformed density rho(q) */
    typedef std::pair<double, double> density_pair;
    /** vector of Fourier transformed densities for |q| value */
    typedef std::vector<density_pair> density_vector;
    /** vector of density vectors */
    typedef std::vector<density_vector> density_vector_vector;
    /** self-intermediate scattering function */
    typedef std::vector<std::vector<double> > isf_vector_vector;
    /** off-diagonal elements of virial stress tensor */
    typedef boost::array<double, (dimension - 1) * dimension / 2> virial_tensor;

    /** Fourier transformed density for different |q| values and vectors */
    boost::shared_ptr<density_vector_vector> rho;
    /** self-intermediate scattering function for different |q| values and vectors */
    boost::shared_ptr<isf_vector_vector> isf;
    /** off-diagonal elements of virial stress tensor */
    boost::shared_ptr<virial_tensor> virial;
    /** time-integral of virial stress tensor */
    boost::shared_ptr<virial_tensor> helfand;
};

/** correlation function result types */
typedef boost::multi_array<accumulator<double>, 2> tcf_unary_result_type;
typedef boost::multi_array<accumulator<double>, 3> tcf_binary_result_type;

template <template <int> class sample_type>
struct correlation_function
{
    /** HDF5 dataset */
    H5::DataSet dataset;
    /** particle type */
    size_t type;
};

/**
 * mean-square displacement
 */
template <template <int> class sample_type>
struct mean_square_displacement;

/**
 * mean-quartic displacement
 */
template <template <int> class sample_type>
struct mean_quartic_displacement;

/**
 * velocity autocorrelation
 */
template <template <int> class sample_type>
struct velocity_autocorrelation;

/**
 * correlation functions sorted by squared displacements
 */
template <template <int> class sample_type>
struct sorted_by_msd;

/**
 * mean squared displacement for most mobile particles given lower boundary
 */
template <template <int> class sample_type>
struct mean_square_displacement_mobile;

/**
 * mean squared displacement for most immobile particles given upper boundary
 */
template <template <int> class sample_type>
struct mean_square_displacement_immobile;

/**
 * mean quartic displacement for most mobile particles given lower boundary
 */
template <template <int> class sample_type>
struct mean_quartic_displacement_mobile;

/**
 * mean quartic displacement for most immobile particles given upper boundary
 */
template <template <int> class sample_type>
struct mean_quartic_displacement_immobile;

/**
 * velocity autocorrelation for most mobile particles given lower boundary
 */
template <template <int> class sample_type>
struct velocity_autocorrelation_mobile;

/**
 * velocity autocorrelation for most immobile particles given upper boundary
 */
template <template <int> class sample_type>
struct velocity_autocorrelation_immobile;

/**
 * intermediate scattering function
 */
template <template <int> class sample_type>
struct intermediate_scattering_function : correlation_function<sample_type>
{
    /** block sample results */
    tcf_binary_result_type result;

    using correlation_function<sample_type>::type;

    char const* name() const { return "ISF"; }

    /**
     * autocorrelate samples in block
     */
    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
        typedef typename input_iterator::first_type sample_iterator;
        typedef typename sample_iterator::value_type::value_type sample_type;
        typedef typename sample_type::density_vector_vector density_vector_vector;
        typedef typename output_iterator::value_type result_vector;

        sample_iterator sample;
        typename density_vector_vector::const_iterator rho0, rho0_0;
        typename density_vector_vector::value_type::const_iterator rho1, rho1_0;
        typename result_vector::iterator result0;

        // accumulate intermediate scattering functions on host
        for (sample = first.first; sample != last.first; ++sample, ++result) {
            for (rho0 = (*sample)[type].rho->begin(), rho0_0 = (*first.first)[type].rho->begin(), result0 = result->begin(); rho0 != (*sample)[type].rho->end(); ++rho0, ++rho0_0, ++result0) {
                for (rho1 = rho0->begin(), rho1_0 = rho0_0->begin(); rho1 != rho0->end(); ++rho1, ++rho1_0) {
                    *result0 += (rho1->first * rho1_0->first + rho1->second * rho1_0->second) / (*sample)[type].r->size();
                }
            }
        }
    }
};

/**
 * self-intermediate scattering function
 */
template <template <int> class sample_type>
struct self_intermediate_scattering_function;

/**
 * squared self-intermediate scattering function
 */
template <template <int> class sample_type>
struct squared_self_intermediate_scattering_function : correlation_function<sample_type>
{
    /** block sample results */
    tcf_binary_result_type result;

    using correlation_function<sample_type>::type;

    char const* name() const { return "SISF2"; }

    /**
     * autocorrelate samples in block
     */
    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
        typedef typename input_iterator::first_type sample_iterator;
        typedef typename sample_iterator::value_type::value_type sample_type;
        typedef typename sample_type::q_vector_vector q_vector_vector;
        typedef typename sample_type::isf_vector_vector isf_vector_vector;
        typedef typename output_iterator::value_type result_vector;

        sample_iterator sample;
        typename isf_vector_vector::iterator isf0;
        typename isf_vector_vector::value_type::iterator isf;
        typename result_vector::iterator result0;

        for (sample = first.first; sample != last.first; ++sample, ++result) {
            for (isf0 = (*sample)[type].isf->begin(), result0 = result->begin(); isf0 != (*sample)[type].isf->end(); ++isf0, ++result0) {
                for (isf = isf0->begin(); isf != isf0->end(); ++isf) {
                    *result0 += (*isf) * (*isf);
                }
            }
        }
    }
};

/**
 * virial stress
 */
template <template <int> class sample_type>
struct virial_stress : correlation_function<sample_type>
{
    /** block sample results */
    tcf_unary_result_type result;

    using correlation_function<sample_type>::type;

    char const* name() const { return "STRESS"; }

    /**
     * autocorrelate samples in block
     */
    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
        typedef typename input_iterator::first_type sample_iterator;
        typedef typename sample_iterator::value_type::value_type sample_type;
        typedef typename sample_type::vector_type vector_type;
        typedef typename sample_type::virial_tensor::const_iterator virial_iterator;
        enum { dimension = vector_type::static_size };

        for (sample_iterator sample = first.first; sample != last.first; ++sample, ++result) {
            for (virial_iterator vir = (*sample)[type].virial->begin(), vir0 = (*first.first)[type].virial->begin(); vir != (*sample)[type].virial->end(); ++vir, ++vir0) {
                *result += ((*vir) * (*vir0)) * (*sample)[type].r->size();
            }
        }
    }
};

/**
 * Helfand moment for virial stress
 *
 * for background see doi:10.1103/PhysRevB.60.3169 and doi:10.1063/1.2724820
 */
template <template <int> class sample_type>
struct helfand_moment : correlation_function<sample_type>
{
    /** block sample results */
    tcf_unary_result_type result;

    using correlation_function<sample_type>::type;

    char const* name() const { return "HELFAND"; }

    /**
     * autocorrelate samples in block
     */
    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
        typedef typename input_iterator::first_type sample_iterator;
        typedef typename sample_iterator::value_type::value_type sample_type;
        typedef typename sample_type::vector_type vector_type;
        typedef typename sample_type::virial_tensor::const_iterator virial_iterator;
        typedef typename sample_type::virial_tensor::value_type virial_value;
        enum { dimension = vector_type::static_size };

        for (sample_iterator sample = first.first; sample != last.first; ++sample, ++result) {
            for (virial_iterator h = (*sample)[type].helfand->begin(),
                 h0 = (*first.first)[type].helfand->begin();
                 h != (*sample)[type].helfand->end(); ++h, ++h0) {
                virial_value dh = *h - *h0;
                *result += (dh * dh) * (*sample)[type].r->size();
            }
        }
    }
};

} // namespace halmd

#endif /* ! HALMD_SAMPLE_TCF_BASE_HPP */
