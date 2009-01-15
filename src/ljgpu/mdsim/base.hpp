/* Molecular Dynamics simulation
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef LJGPU_MDSIM_BASE_HPP
#define LJGPU_MDSIM_BASE_HPP

#include <cmath>
#include <limits>
#include <ljgpu/mdsim/impl.hpp>
#include <ljgpu/mdsim/traits.hpp>
#include <ljgpu/sample/H5param.hpp>
#include <ljgpu/sample/perf.hpp>
#include <ljgpu/util/exception.hpp>
#include <ljgpu/util/log.hpp>

namespace ljgpu
{

template <typename mdsim_impl>
class mdsim_base
{
public:
    typedef mdsim_impl impl_type;
    typedef mdsim_traits<impl_type> traits_type;
    typedef typename traits_type::float_type float_type;
    typedef typename traits_type::vector_type vector_type;
    typedef typename traits_type::sample_type sample_type;
    typedef typename sample_type::sample_visitor sample_visitor;
    enum { dimension = traits_type::dimension };

public:
    /** set number of particles */
    void particles(unsigned int value);
    /** set particle density */
    void density(float_type value);
    /** set periodic box length */
    void box(float_type value);
    /** set system state from phase space sample */
    void sample(sample_visitor read);

    /** returns number of particles */
    unsigned int particles() const { return npart; }
    /** returns particle density */
    float_type density() const { return density_; }
    /** returns periodic box length */
    float_type box() const { return box_; }
    /** returns trajectory sample */
    sample_type const& sample() const { return m_sample; }
    /** returns and resets CPU or GPU time accumulators */
    perf::counters times();

protected:
    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

protected:
    /** number of particles */
    unsigned int npart;
    /** particle density */
    float_type density_;
    /** periodic box length */
    float_type box_;
    /** trajectory sample in swappable host memory */
    sample_type m_sample;
    /** GPU time accumulators */
    perf::counters m_times;
};

template <typename mdsim_impl>
void mdsim_base<mdsim_impl>::particles(unsigned int value)
{
    // validate particle number
    if (value < 1) {
	throw exception("invalid number of particles");
    }
    // set particle number
    npart = value;
    LOG("number of particles: " << npart);
}

template <typename mdsim_impl>
void mdsim_base<mdsim_impl>::density(float_type value)
{
    // set particle density
    density_ = value;
    LOG("particle density: " << density_);

    // compute periodic box length
    box_ = std::pow(npart / density_, (float_type) 1 / dimension);
    LOG("periodic simulation box length: " << box_);
}

template <typename mdsim_impl>
void mdsim_base<mdsim_impl>::box(float_type value)
{
    // set periodic box length
    box_ = value;
    LOG("periodic simulation box length: " << box_);

    // compute particle density
    density_ = npart / std::pow(box_, (float_type) dimension);
    LOG("particle density: " << density_);
}

template <typename mdsim_impl>
void mdsim_base<mdsim_impl>::sample(mdsim_base<mdsim_impl>::sample_visitor read)
{
    typedef typename sample_type::position_vector position_vector;
    typedef typename position_vector::value_type position_value;

    read(m_sample);

    if (m_sample.r.size() != npart) {
	throw exception("mismatching number of particles in phase space sample");
    }

    position_value const box = m_sample.box;
    BOOST_FOREACH(position_vector &r, m_sample.r) {
	// apply periodic boundary conditions to positions
	r -= floor(r / box) * box;
    }

    if (std::fabs(box - box_) > (box_ * std::numeric_limits<float>::epsilon())) {
	LOG("rescaling periodic simulation box length from " << box);
	position_value const scale = box_ / box;
	BOOST_FOREACH(position_vector &r, m_sample.r) {
	    r *= scale;
	}
	m_sample.box = box_;
    }
}

template <typename mdsim_impl>
perf::counters mdsim_base<mdsim_impl>::times()
{
    perf::counters times(m_times);
    BOOST_FOREACH(perf::counter& i, m_times) {
	// reset performance counter
	i.second.clear();
    }
    return times;
}

template <typename mdsim_impl>
void mdsim_base<mdsim_impl>::param(H5param& param) const
{
    H5xx::group node(param["mdsim"]);
    node["box_length"] = box_;
    node["density"] = density_;
    node["particles"] = npart;
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_BASE_HPP */
