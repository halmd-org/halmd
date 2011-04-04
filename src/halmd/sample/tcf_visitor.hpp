/* Time correlation function visitors
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

#ifndef HALMD_SAMPLE_TCF_VISITOR_HPP
#define HALMD_SAMPLE_TCF_VISITOR_HPP

#include <algorithm>
#include <boost/array.hpp>
// requires boost 1.37.0 or patch from http://svn.boost.org/trac/boost/ticket/1852
#include <boost/circular_buffer.hpp>
#include <boost/multi_array.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/variant.hpp>
#include <vector>

#ifdef WITH_CUDA
# include <halmd/sample/tcf_gpu.hpp>
#endif
#include <halmd/sample/tcf_host.hpp>
#include <halmd/util/H5xx.hpp>

namespace halmd {

class tcf_sample_phase_space : public boost::static_visitor<>
{
public:
    /**
     * sample from host memory to host memory
     */
    template <typename T, int dimension>
    void operator()(std::vector<tcf_host_sample<dimension> >& sample, std::vector<trajectory_host_sample<T, dimension> > const& sample_) const
    {
        typedef std::vector<trajectory_host_sample<T, dimension> > trajectory_sample_type;
        typedef typename trajectory_sample_type::const_iterator sample_iterator;
        typedef tcf_host_sample<dimension> sample_type;
        typedef typename sample_type::sample_vector sample_vector;
        typedef boost::shared_ptr<sample_vector> sample_ptr;

        for (sample_iterator s = sample_.begin(); s != sample_.end(); ++s) {
            sample_ptr r(new sample_vector(s->r->begin(), s->r->end()));
            sample_ptr v(new sample_vector(s->v->begin(), s->v->end()));
            sample.push_back(sample_type(r, v));
        }
    }

#ifdef WITH_CUDA
    template <int dimension>
    void operator()(std::vector<tcf_host_sample<dimension> >& sample, std::vector<trajectory_gpu_sample<dimension> > const& sample_) const
    {
        throw std::runtime_error("correlation host backend does not support GPU samples");
    }

    /**
     * sample from global device memory to global device memory
     */
    template <int dimension>
    void operator()(std::vector<tcf_gpu_sample<dimension> >& sample, std::vector<trajectory_gpu_sample<dimension> > const& sample_) const
    {
        typedef std::vector<trajectory_gpu_sample<dimension> > trajectory_sample_type;
        typedef typename trajectory_sample_type::const_iterator sample_iterator;

        for (sample_iterator s = sample_.begin(); s != sample_.end(); ++s) {
            // copy shared global device memory pointers
            sample.push_back(tcf_gpu_sample<dimension>(s->r, s->v));
        }
    }

    /**
     * sample from host memory to global device memory
     */
    template <typename T, int dimension>
    void operator()(std::vector<tcf_gpu_sample<dimension> >& sample, std::vector<trajectory_host_sample<T, dimension> > const& sample_) const
    {
        typedef std::vector<trajectory_host_sample<T, dimension> > trajectory_sample_type;
        typedef typename trajectory_sample_type::const_iterator sample_iterator;
        typedef tcf_gpu_sample<dimension> sample_type;
        typedef typename sample_type::gpu_sample_vector gpu_sample_vector;
        typedef boost::shared_ptr<gpu_sample_vector> gpu_sample_ptr;
        typedef typename sample_type::gpu_vector_type gpu_vector_type;

        for (sample_iterator s = sample_.begin(); s != sample_.end(); ++s) {
            // copy sample to page-locked host memory
            cuda::host::vector<gpu_vector_type> h_r(s->r->begin(), s->r->end());
            cuda::host::vector<gpu_vector_type> h_v(s->v->begin(), s->v->end());
            // copy from host to GPU
            gpu_sample_ptr g_r(new gpu_sample_vector(h_r.size()));
            gpu_sample_ptr g_v(new gpu_sample_vector(h_v.size()));
            sample.push_back(sample_type(g_r, g_v));
            cuda::copy(h_r, *g_r);
            cuda::copy(h_v, *g_v);
        }
    }
#endif /* WITH_CUDA */
};

/**
 * compute Fourier-transformed densities
 */
template <typename U>
class _tcf_fourier_transform_sample : public boost::static_visitor<>
{
public:
    _tcf_fourier_transform_sample(U const& q_vector) : q_vector(q_vector) {}

    template <typename T>
    void operator()(T& sample) const
    {
        for (typename T::iterator s = sample.begin(); s != sample.end(); ++s) {
            (*s)(q_vector);
        }
    }

private:
    U const& q_vector;
};

template <typename U>
_tcf_fourier_transform_sample<U> tcf_fourier_transform_sample(U const& q_vector)
{
   return _tcf_fourier_transform_sample<U>(q_vector);
}

/**
 * copy off-diagonal elements of virial tensor
 */
template <typename U>
class _tcf_sample_virial : public boost::static_visitor<>
{
public:
    _tcf_sample_virial(U const& virial) : virial(virial) {}

    template <typename T>
    void operator()(T& sample) const
    {
        typename U::const_iterator vir = virial.begin();
        for (typename T::iterator s = sample.begin(); s != sample.end(); ++s, ++vir) {
            s->virial.reset(new typename T::value_type::virial_tensor);
            std::copy(vir->begin() + 1, vir->end(), s->virial->begin());
        }
    }

private:
    U const& virial;
};

template <typename U>
_tcf_sample_virial<U> tcf_sample_virial(U const& virial)
{
   return _tcf_sample_virial<U>(virial);
}

/**
 * copy off-diagonal elements of Helfand moment
 */
template <typename U>
class _tcf_sample_helfand : public boost::static_visitor<>
{
public:
    _tcf_sample_helfand(U const& helfand) : helfand(helfand) {}

    template <typename T>
    void operator()(T& sample) const
    {
        assert(helfand.size() == sample.size());
        typename U::const_iterator h = helfand.begin();
        for (typename T::iterator s = sample.begin(); s != sample.end(); ++s, ++h) {
            s->helfand.reset(new typename T::value_type::virial_tensor);
            std::copy(h->begin() + 1, h->end(), s->helfand->begin());
        }
    }

private:
    U const& helfand;
};

template <typename U>
_tcf_sample_helfand<U> tcf_sample_helfand(U const& helfand)
{
    return _tcf_sample_helfand<U>(helfand);
}

class tcf_block_add_sample : public boost::static_visitor<>
{
public:
    template <typename T, typename U>
    void operator()(T const&, U const&) const
    {
        throw std::runtime_error("block sample mismatch");
    }

    template <typename T>
    void operator()(boost::circular_buffer<T>& block, T const& sample) const
    {
        block.push_back(sample);
    }
};

class tcf_block_is_full : public boost::static_visitor<bool>
{
public:
    template <typename T>
    bool operator()(boost::circular_buffer<T> const& block) const
    {
        return block.full();
    }
};

class tcf_block_clear : public boost::static_visitor<>
{
public:
    template <typename T>
    void operator()(boost::circular_buffer<T>& block) const
    {
        block.clear();
    }
};

class tcf_block_pop_front : public boost::static_visitor<>
{
public:
    template <typename T>
    void operator()(boost::circular_buffer<T>& block) const
    {
        block.pop_front();
    }
};

class tcf_block_size : public boost::static_visitor<size_t>
{
public:
    template <typename T>
    size_t operator()(boost::circular_buffer<T> const& block) const
    {
        return block.size();
    }
};

/**
 * set particle type for a correlation function
 */
class tcf_set_type : public boost::static_visitor<>
{
public:
    tcf_set_type(size_t type) : type(type) {}

    template <typename T>
    void operator()(T& tcf) const
    {
        tcf.type = type;
    }

private:
    size_t type;
};

class tcf_add_mobile_particle_filter : public boost::static_visitor<>
{
public:
    tcf_add_mobile_particle_filter(float fraction) : fraction(fraction) {}

    template <template <int> class sample_type>
    void operator()(mean_square_displacement_mobile<sample_type>& tcf) const
    {
        tcf.mobile_fraction.push_back(fraction);
    }

    template <template <int> class sample_type>
    void operator()(mean_quartic_displacement_mobile<sample_type>& tcf) const
    {
        tcf.mobile_fraction.push_back(fraction);
    }

    template <template <int> class sample_type>
    void operator()(velocity_autocorrelation_mobile<sample_type>& tcf) const
    {
        tcf.mobile_fraction.push_back(fraction);
    }

    template <typename T>
    void operator()(T& tcf) const
    {
        // noop
    }

private:
    float fraction;
};

class tcf_add_immobile_particle_filter : public boost::static_visitor<>
{
public:
    tcf_add_immobile_particle_filter(float fraction) : fraction(fraction) {}

    template <template <int> class sample_type>
    void operator()(mean_square_displacement_immobile<sample_type>& tcf) const
    {
        tcf.immobile_fraction.push_back(fraction);
    }

    template <template <int> class sample_type>
    void operator()(mean_quartic_displacement_immobile<sample_type>& tcf) const
    {
        tcf.immobile_fraction.push_back(fraction);
    }

    template <template <int> class sample_type>
    void operator()(velocity_autocorrelation_immobile<sample_type>& tcf) const
    {
        tcf.immobile_fraction.push_back(fraction);
    }

    template <typename T>
    void operator()(T& tcf) const
    {
        // noop
    }

private:
    float fraction;
};

#ifdef WITH_CUDA
class tcf_get_sorted_by_msd
  : public boost::static_visitor<
        boost::multi_array<cuda::vector<float4>, 2> const*
    >
{
public:
    tcf_get_sorted_by_msd(size_t type) : type(type) {}

    template <template <int> class sample_type>
    boost::multi_array<cuda::vector<float4>, 2> const*
        operator()(sorted_by_msd<sample_type> const& tcf) const
    {
        if (tcf.type == type) {
            return &tcf.result;
        }
        return (boost::multi_array<cuda::vector<float4>, 2> const*) 0; // noop
    }

    template <typename T>
    boost::multi_array<cuda::vector<float4>, 2> const*
        operator()(T const& tcf) const
    {
        return (boost::multi_array<cuda::vector<float4>, 2> const*) 0; // noop
    }

private:
    size_t const type;
};
#endif /* WITH_CUDA */

/**
 * apply correlation function to block of phase space samples
 */
template <typename V, typename tcf_vector_type>
class _tcf_correlate_block : public boost::static_visitor<>
{
public:
    _tcf_correlate_block(unsigned int block, V const& q_vector, tcf_vector_type const& tcf_vector) : block(block), q_vector(q_vector), tcf_vector(tcf_vector) {}

    template <typename T, typename U>
    void operator()(T&, U&) const
    {
        throw std::runtime_error("correlation function mismatch");
    }

    template <typename T, template <int> class sample_type, int dimension>
    typename boost::enable_if<boost::is_base_of<correlation_function<sample_type>, T>, void>::type
    operator()(T& tcf, boost::circular_buffer<std::vector<sample_type<dimension> > >& sample) const
    {
        if (tcf.result.num_elements()) {
            tcf(
                std::make_pair(sample.begin(), q_vector.begin())
              , std::make_pair(sample.end(), q_vector.end())
              , tcf.result[block].begin()
            );
        }
    }

#ifdef WITH_CUDA
    template <template <int> class sample_type, int dimension>
    typename boost::enable_if<
        boost::is_base_of<
            correlation_function<sample_type>
          , mean_square_displacement_mobile<sample_type>
        >, void>::type
    operator()(
        mean_square_displacement_mobile<sample_type>& tcf
      , boost::circular_buffer<std::vector<sample_type<dimension> > >& sample
    ) const
    {
        if (tcf.result.num_elements()) {
            typename tcf_vector_type::const_iterator it, end = tcf_vector.end();
            for (it = tcf_vector.begin(); it != end; ++it) {
                if (boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it)) {
                    tcf(
                        /* ignore zero MSD */ ++(*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].begin()
                      , (*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].end()
                      , /* ignore zero MSD */ ++(tcf.result[block].begin())
                    );
                    return;
                }
            }
            throw std::logic_error("no pointer to sorted mean square displacements");
        }
    }

    template <template <int> class sample_type, int dimension>
    typename boost::enable_if<
        boost::is_base_of<
            correlation_function<sample_type>
          , mean_square_displacement_immobile<sample_type>
        >, void>::type
    operator()(
        mean_square_displacement_immobile<sample_type>& tcf
      , boost::circular_buffer<std::vector<sample_type<dimension> > >& sample
    ) const
    {
        if (tcf.result.num_elements()) {
            typename tcf_vector_type::const_iterator it, end = tcf_vector.end();
            for (it = tcf_vector.begin(); it != end; ++it) {
                if (boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it)) {
                    tcf(
                        /* ignore zero MSD */ ++(*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].begin()
                      , (*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].end()
                      , /* ignore zero MSD */ ++(tcf.result[block].begin())
                    );
                    return;
                }
            }
            throw std::logic_error("no pointer to sorted mean square displacements");
        }
    }

    template <template <int> class sample_type, int dimension>
    typename boost::enable_if<
        boost::is_base_of<
            correlation_function<sample_type>
          , mean_quartic_displacement_mobile<sample_type>
        >, void>::type
    operator()(
        mean_quartic_displacement_mobile<sample_type>& tcf
      , boost::circular_buffer<std::vector<sample_type<dimension> > >& sample
    ) const
    {
        if (tcf.result.num_elements()) {
            typename tcf_vector_type::const_iterator it, end = tcf_vector.end();
            for (it = tcf_vector.begin(); it != end; ++it) {
                if (boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it)) {
                    tcf(
                        /* ignore zero MSD */ ++(*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].begin()
                      , (*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].end()
                      , /* ignore zero MSD */ ++(tcf.result[block].begin())
                    );
                    return;
                }
            }
            throw std::logic_error("no pointer to sorted mean quartic displacements");
        }
    }

    template <template <int> class sample_type, int dimension>
    typename boost::enable_if<
        boost::is_base_of<
            correlation_function<sample_type>
          , mean_quartic_displacement_immobile<sample_type>
        >, void>::type
    operator()(
        mean_quartic_displacement_immobile<sample_type>& tcf
      , boost::circular_buffer<std::vector<sample_type<dimension> > >& sample
    ) const
    {
        if (tcf.result.num_elements()) {
            typename tcf_vector_type::const_iterator it, end = tcf_vector.end();
            for (it = tcf_vector.begin(); it != end; ++it) {
                if (boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it)) {
                    tcf(
                        /* ignore zero MSD */ ++(*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].begin()
                      , (*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].end()
                      , /* ignore zero MSD */ ++(tcf.result[block].begin())
                    );
                    return;
                }
            }
            throw std::logic_error("no pointer to sorted mean quartic displacements");
        }
    }

    template <template <int> class sample_type, int dimension>
    typename boost::enable_if<
        boost::is_base_of<
            correlation_function<sample_type>
          , velocity_autocorrelation_mobile<sample_type>
        >, void>::type
    operator()(
        velocity_autocorrelation_mobile<sample_type>& tcf
      , boost::circular_buffer<std::vector<sample_type<dimension> > >& sample
    ) const
    {
        if (tcf.result.num_elements()) {
            typename tcf_vector_type::const_iterator it, end = tcf_vector.end();
            for (it = tcf_vector.begin(); it != end; ++it) {
                if (boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it)) {
                    tcf(
                        /* ignore zero MSD */ ++(*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].begin()
                      , (*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].end()
                      , /* ignore zero MSD */ ++(tcf.result[block].begin())
                    );
                    return;
                }
            }
            throw std::logic_error("no pointer to sorted velocity autocorrelations");
        }
    }

    template <template <int> class sample_type, int dimension>
    typename boost::enable_if<
        boost::is_base_of<
            correlation_function<sample_type>
          , velocity_autocorrelation_immobile<sample_type>
        >, void>::type
    operator()(
        velocity_autocorrelation_immobile<sample_type>& tcf
      , boost::circular_buffer<std::vector<sample_type<dimension> > >& sample
    ) const
    {
        if (tcf.result.num_elements()) {
            typename tcf_vector_type::const_iterator it, end = tcf_vector.end();
            for (it = tcf_vector.begin(); it != end; ++it) {
                if (boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it)) {
                    tcf(
                        /* ignore zero MSD */ ++(*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].begin()
                      , (*boost::apply_visitor(tcf_get_sorted_by_msd(tcf.type), *it))[block].end()
                      , /* ignore zero MSD */ ++(tcf.result[block].begin())
                    );
                    return;
                }
            }
            throw std::logic_error("no pointer to sorted velocity autocorrelations");
        }
    }
#endif /* WITH_CUDA */

private:
    unsigned int block;
    V const& q_vector;
    tcf_vector_type const& tcf_vector;
};

template <typename T, typename U>
_tcf_correlate_block<T, U> tcf_correlate_block(unsigned int block, T const& q_vector, U const& tcf_vector)
{
    return _tcf_correlate_block<T, U>(block, q_vector, tcf_vector);
}

/**
 * retrieve name of a correlation function
 */
class tcf_name : public boost::static_visitor<char const*>
{
public:
    template <typename T>
    char const* operator()(T const& tcf) const
    {
        return tcf.name();
    }
};

/**
 * create correlation function HDF5 dataset
 */
class tcf_create_dataset : public boost::static_visitor<>
{
public:
    tcf_create_dataset(H5::H5File& file, size_t types) : file(file), types(types) {}

    template <typename T>
    void operator()(T& tcf) const
    {
        if (tcf.result.num_elements()) {
            H5::Group root(file.openGroup("/"));
            if (types > 1) {
                std::string name;
                // AA, BB, ... correlation function
                name.push_back('A' + tcf.type);
                name.push_back('A' + tcf.type);
                try {
                    H5XX_NO_AUTO_PRINT(H5::GroupIException);
                    root = root.createGroup(name);
                }
                catch (H5::GroupIException const&) {
                    root = root.openGroup(name);
                }
            }
            tcf.dataset = create_dataset(root, tcf.name(), tcf.result);
        }
    }

    template <template <int> class sample_type>
    void operator()(sorted_by_msd<sample_type>& tcf) const
    {
        // noop
    }

    static H5::DataSet create_dataset(H5::Group const& node, char const* name, tcf_unary_result_type const& result)
    {
        // extensible dataspace for unary correlation function results
        hsize_t dim[3] = { 0, result.shape()[1], 5 };
        hsize_t max_dim[3] = { H5S_UNLIMITED, result.shape()[1], 5 };
        hsize_t chunk_dim[3] = { 1, result.shape()[1], 3 };
        H5::DataSpace ds(3, dim, max_dim);
        H5::DSetCreatPropList cparms;
        cparms.setChunk(3, chunk_dim);

        return node.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds, cparms);
    }

    static H5::DataSet create_dataset(H5::Group const& node, char const* name, tcf_binary_result_type const& result)
    {
        // extensible dataspace for binary correlation function results
        hsize_t dim[4] = { result.shape()[2], 0, result.shape()[1], 6 };
        hsize_t max_dim[4] = { result.shape()[2], H5S_UNLIMITED, result.shape()[1], 6 };
        hsize_t chunk_dim[4] = { result.shape()[2], 1, result.shape()[1], 4 };
        H5::DataSpace ds(4, dim, max_dim);
        H5::DSetCreatPropList cparms;
        cparms.setChunk(4, chunk_dim);

        return node.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds, cparms);
    }

private:
    H5::H5File& file;
    size_t types;
};

/**
 * allocate correlation functions results
 */
class tcf_allocate_results : public boost::static_visitor<>
{
public:
    tcf_allocate_results(unsigned int block_count, unsigned int block_size, unsigned int q_values)
        : block_count(block_count), block_size(block_size), q_values(q_values) {}

    template <typename T>
    void operator()(T& tcf) const
    {
        resize(tcf.result, q_values);
    }

    template <template <int> class sample_type>
    void operator()(sorted_by_msd<sample_type>& tcf) const
    {
        tcf.result.resize(boost::extents[block_count][block_size]);
    }

    template <template <int> class sample_type>
    void operator()(mean_square_displacement_mobile<sample_type>& tcf) const
    {
        resize(tcf.result, tcf.mobile_fraction.size());
    }

    template <template <int> class sample_type>
    void operator()(mean_square_displacement_immobile<sample_type>& tcf) const
    {
        resize(tcf.result, tcf.immobile_fraction.size());
    }

    template <template <int> class sample_type>
    void operator()(mean_quartic_displacement_mobile<sample_type>& tcf) const
    {
        resize(tcf.result, tcf.mobile_fraction.size());
    }

    template <template <int> class sample_type>
    void operator()(mean_quartic_displacement_immobile<sample_type>& tcf) const
    {
        resize(tcf.result, tcf.immobile_fraction.size());
    }

    template <template <int> class sample_type>
    void operator()(velocity_autocorrelation_mobile<sample_type>& tcf) const
    {
        resize(tcf.result, tcf.mobile_fraction.size());
    }

    template <template <int> class sample_type>
    void operator()(velocity_autocorrelation_immobile<sample_type>& tcf) const
    {
        resize(tcf.result, tcf.immobile_fraction.size());
    }

    void resize(tcf_unary_result_type& result, unsigned int) const
    {
        result.resize(boost::extents[block_count][block_size]);
    }

    void resize(tcf_binary_result_type& result, unsigned int count) const
    {
        result.resize(boost::extents[block_count][block_size][count]);
    }

private:
    unsigned int block_count;
    unsigned int block_size;
    unsigned int q_values;
};

/**
 * write correlation function results to HDF5 file
 */
class tcf_write_results : public boost::static_visitor<>
{
public:
    typedef boost::multi_array<double, 2> block_time_type;
    typedef std::vector<double> q_value_vector;

public:
    tcf_write_results(block_time_type const& block_time, q_value_vector const& q_value, unsigned int max_blocks)
        : block_time(block_time), q_value(q_value), max_blocks(max_blocks) {}

    template <typename T>
    void operator()(T& tcf) const
    {
        if (tcf.result.num_elements()) {
            write(tcf.dataset, tcf.result, q_value);
        }
    }

    template <template <int> class sample_type>
    void operator()(sorted_by_msd<sample_type>& tcf) const
    {
        // noop
    }

    template <template <int> class sample_type>
    void operator()(mean_square_displacement_mobile<sample_type>& tcf) const
    {
        if (tcf.result.num_elements()) {
            write(tcf.dataset, tcf.result, tcf.mobile_fraction);
        }
    }

    template <template <int> class sample_type>
    void operator()(mean_square_displacement_immobile<sample_type>& tcf) const
    {
        if (tcf.result.num_elements()) {
            write(tcf.dataset, tcf.result, tcf.immobile_fraction);
        }
    }

    template <template <int> class sample_type>
    void operator()(mean_quartic_displacement_mobile<sample_type>& tcf) const
    {
        if (tcf.result.num_elements()) {
            write(tcf.dataset, tcf.result, tcf.mobile_fraction);
        }
    }

    template <template <int> class sample_type>
    void operator()(mean_quartic_displacement_immobile<sample_type>& tcf) const
    {
        if (tcf.result.num_elements()) {
            write(tcf.dataset, tcf.result, tcf.immobile_fraction);
        }
    }

    template <template <int> class sample_type>
    void operator()(velocity_autocorrelation_mobile<sample_type>& tcf) const
    {
        if (tcf.result.num_elements()) {
            write(tcf.dataset, tcf.result, tcf.mobile_fraction);
        }
    }

    template <template <int> class sample_type>
    void operator()(velocity_autocorrelation_immobile<sample_type>& tcf) const
    {
        if (tcf.result.num_elements()) {
            write(tcf.dataset, tcf.result, tcf.immobile_fraction);
        }
    }

    template <typename vector_type>
    void write(H5::DataSet& dataset, tcf_unary_result_type const& result, vector_type const&) const
    {
        // dataset dimensions
        boost::array<hsize_t, 3> dim = {{ max_blocks, result.shape()[1], 5 }};
        // memory buffer for results
        boost::multi_array<double, 3> data(dim);

        for (unsigned int j = 0; j < dim[0]; ++j) {
            for (unsigned int k = 0; k < dim[1]; ++k) {
                // time interval
                data[j][k][0] = block_time[j][k];
                // mean average
                data[j][k][1] = result[j][k].mean();
                // standard error of mean
                data[j][k][2] = result[j][k].err();
                // variance
                data[j][k][3] = result[j][k].var();
                // count
                data[j][k][4] = result[j][k].count();
            }
        }
        dataset.extend(dim.c_array());
        // write results to HDF5 file
        dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    }

    template <typename vector_type>
    void write(H5::DataSet& dataset, tcf_binary_result_type const& result, vector_type const& vector) const
    {
        // dataset dimensions
        boost::array<hsize_t, 4> dim = {{ result.shape()[2], max_blocks, result.shape()[1], 6 }};
        // memory buffer for results
        boost::multi_array<double, 4> data(dim);

        for (unsigned int j = 0; j < dim[0]; ++j) {
            for (unsigned int k = 0; k < dim[1]; ++k) {
                for (unsigned int l = 0; l < dim[2]; ++l) {
                    // q-value
                    data[j][k][l][0] = vector[j];
                    // time interval
                    data[j][k][l][1] = block_time[k][l];
                    // mean average
                    data[j][k][l][2] = result[k][l][j].mean();
                    // standard error of mean
                    data[j][k][l][3] = result[k][l][j].err();
                    // variance
                    data[j][k][l][4] = result[k][l][j].var();
                    // count
                    data[j][k][l][5] = result[k][l][j].count();
                }
            }
        }
        dataset.extend(dim.c_array());
        // write results to HDF5 file
        dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    }

private:
    block_time_type const& block_time;
    q_value_vector const& q_value;
    unsigned int max_blocks;
};

} // namespace halmd

#endif /* ! HALMD_SAMPLE_TCF_VISITOR_HPP */
