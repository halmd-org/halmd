/* Molecular Dynamics simulation
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_MDSIM_BASE_HPP
#define HALMD_MDSIM_BASE_HPP

#include <boost/assign.hpp>
#include <boost/noncopyable.hpp>
#include <boost/type_traits.hpp>
#include <boost/variant.hpp>
#include <cmath>
#include <limits>

#include <halmd/mdsim/exception.hpp>
#include <halmd/mdsim/impl.hpp>
#include <halmd/mdsim/traits.hpp>
#include <halmd/mdsim/variant.hpp>
#include <halmd/sample/H5param.hpp>
#include <halmd/sample/perf.hpp>
#include <halmd/util/exception.hpp>
#include <halmd/util/log.hpp>

#define foreach BOOST_FOREACH

namespace halmd
{

template <typename mdsim_impl, int dimension_>
class mdsim_base : boost::noncopyable
{
public:
    typedef mdsim_impl impl_type;
    typedef mdsim_traits<impl_type, dimension_> traits_type;
    typedef typename traits_type::float_type float_type;
    typedef typename traits_type::vector_type vector_type;
    typedef typename traits_type::host_sample_type host_sample_type;
    typedef host_sample_type trajectory_sample_type;
    typedef boost::variant<host_sample_type> trajectory_sample_variant;
    typedef typename traits_type::energy_sample_type energy_sample_type;
    typedef typename energy_sample_type::virial_tensor virial_tensor;
    enum { dimension = dimension_ };

public:
    mdsim_base() :
        mixture_(UNARY) {}

    /** set number of particles */
    void particles(unsigned int value);
    /** set number of A and B particles in binary mixture */
    void particles(boost::array<unsigned int, 2> const& value);
    /** set particle density */
    void density(float_type value);
    /** set periodic box length */
    void box(float_type value);
    /** set system state from phase space sample */
    void state(host_sample_type& sample, float_type box);

    /** returns number of particles */
    unsigned int particles() const { return npart; }
    /** returns particle density */
    float_type density() const { return density_; }
    /** returns periodic box length */
    float_type box() const { return box_; }
    /** returns and resets CPU or GPU time accumulators */
    perf::counters times();

    mixture_type mixture() const { return mixture_; }
    bool is_binary () const { return mixture_ == BINARY; }

protected:
    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

protected:
    /** number of particles */
    unsigned int npart;
    /** number of A and B particles in binary mixture */
    boost::array<unsigned int, 2> mpart;
    /** particle density */
    float_type density_;
    /** periodic box length */
    float_type box_;
    /** GPU time accumulators */
    perf::counters mutable m_times;

    /** uniform fluid or binary mixture */
    mixture_type mixture_;
};

template <typename mdsim_impl, int dimension_>
void mdsim_base<mdsim_impl, dimension_>::particles(unsigned int value)
{
    if (value < 1) {
        throw exception("invalid number of particles");
    }
    npart = value;
    mpart = boost::assign::list_of(npart)(0);
    mixture_ = UNARY;
    LOG("number of particles: " << npart);
}

template <typename mdsim_impl, int dimension_>
void mdsim_base<mdsim_impl, dimension_>::particles(boost::array<unsigned int, 2> const& value)
{
    if (*std::min_element(value.begin(), value.end()) < 1) {
        throw exception("invalid number of A or B particles");
    }
    mpart = value;
    npart = std::accumulate(mpart.begin(), mpart.end(), 0);
    mixture_ = BINARY;
    LOG("binary mixture with " << mpart[0] << " A particles and " << mpart[1] << " B particles");
}

template <typename mdsim_impl, int dimension_>
void mdsim_base<mdsim_impl, dimension_>::density(float_type value)
{
    // set particle density
    density_ = value;
    LOG("particle density: " << density_);

    // compute periodic box length
    box_ = std::pow(npart / density_, (float_type) 1 / dimension);
    LOG("periodic simulation box length: " << box_);
}

template <typename mdsim_impl, int dimension_>
void mdsim_base<mdsim_impl, dimension_>::box(float_type value)
{
    // set periodic box length
    box_ = value;
    LOG("periodic simulation box length: " << box_);

    // compute particle density
    density_ = npart / std::pow(box_, (float_type) dimension);
    LOG("particle density: " << density_);
}

template <typename mdsim_impl, int dimension_>
void mdsim_base<mdsim_impl, dimension_>::state(host_sample_type& sample, float_type box)
{
    typedef typename host_sample_type::value_type sample_type;
    typedef typename sample_type::position_vector position_vector;
    typedef typename position_vector::value_type position_value;

    for (size_t i = 0; i < sample.size(); ++i) {
        if (sample[i].r->size() != mpart[i]) {
            throw exception("mismatching number of particles in phase space sample");
        }
    }

    foreach (sample_type& s, sample) {
        foreach (position_vector &r, *s.r) {
            // apply periodic boundary conditions to positions
            r -= floor(r / static_cast<position_value>(box)) * static_cast<position_value>(box);
        }
    }

    if (std::fabs(box - box_) > (box_ * std::numeric_limits<float>::epsilon())) {
        position_value const scale = box_ / box;
        foreach (sample_type& s, sample) {
            foreach (position_vector &r, *s.r) {
                r *= scale;
            }
        }
        LOG("rescaled periodic simulation box length from " << box);
    }
}

template <typename mdsim_impl, int dimension_>
perf::counters mdsim_base<mdsim_impl, dimension_>::times()
{
    perf::counters times(m_times);
    foreach (perf::counter& i, m_times) {
        // reset performance counter
        i.second.clear();
    }
    return times;
}

template <typename mdsim_impl, int dimension_>
void mdsim_base<mdsim_impl, dimension_>::param(H5param& param) const
{
    H5xx::group node(param["mdsim"]);
    node["box_length"] = box_;
    node["density"] = density_;
    if (mixture_ == BINARY) {
        node["particles"] = mpart;
    }
    else {
        node["particles"] = npart;
    }
}

} // namespace halmd

#undef foreach

#endif /* ! HALMD_MDSIM_BASE_HPP */
