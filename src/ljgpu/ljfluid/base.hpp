/* Lennard-Jones fluid simulation
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

#ifndef LJGPU_LJFLUID_BASE_HPP
#define LJGPU_LJFLUID_BASE_HPP

#include <boost/foreach.hpp>
#include <ljgpu/ljfluid/traits.hpp>
#include <ljgpu/sample/perf.hpp>
#include <ljgpu/sample/sample.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <ljgpu/util/exception.hpp>
#include <ljgpu/util/log.hpp>

namespace ljgpu
{

/**
 * Lennard-Jones fluid interface
 */
template <typename ljfluid_impl>
class ljfluid_base
{
public:
    typedef ljfluid_traits<ljfluid_impl> traits_type;
    typedef typename traits_type::float_type float_type;
    typedef typename traits_type::vector_type vector_type;
    typedef typename traits_type::sample_type sample_type;
    enum { dimension = traits_type::dimension };

public:
    /** set number of particles */
    void particles(unsigned int value);
    /** set particle density */
    void density(float_type value);
    /** set periodic box length */
    void box(float_type value);
    /** set simulation timestep */
    void timestep(float_type value);
    /** set potential cutoff radius */
    void cutoff_radius(float_type value);
#ifdef USE_POTENTIAL_SMOOTHING
    /** set potential smoothing function scale parameter */
    void potential_smoothing(float_type value);
#endif

    /** returns number of particles */
    unsigned int particles() const { return npart; }
    /** returns particle density */
    float_type density() const { return density_; }
    /** returns periodic box length */
    float_type box() const { return box_; }
    /** returns simulation timestep */
    float_type timestep() const { return timestep_; }
    /** returns potential cutoff radius */
    float_type cutoff_radius() const { return r_cut; }
#ifdef USE_POTENTIAL_SMOOTHING
    /** returns potential smoothing function scale parameter */
    float_type potential_smoothing() const { return r_smooth; }
#endif

    /** returns trajectory sample */
    sample_type const& sample() const { return m_sample; }
    /** returns and resets CPU or GPU time accumulators */
    perf_counters times();

    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

protected:
    /** number of particles in system */
    unsigned int npart;
    /** particle density */
    float_type density_;
    /** periodic box length */
    float_type box_;
    /** simulation timestep */
    float_type timestep_;
    /** cutoff radius for shifted Lennard-Jones potential */
    float_type r_cut;
    /** squared cutoff radius */
    float_type rr_cut;
    /** potential energy at cutoff radius */
    float_type en_cut;
#ifdef USE_POTENTIAL_SMOOTHING
    /** potential smoothing function scale parameter */
    float_type r_smooth;
    /** squared inverse potential smoothing function scale parameter */
    float_type rri_smooth;
#endif

    /** trajectory sample in swappable host memory */
    sample_type m_sample;
    /** GPU time accumulators */
    perf_counters m_times;
};

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::particles(unsigned int value)
{
    // validate particle number
    if (value < 1) {
	throw exception("invalid number of particles");
    }
    // set particle number
    npart = value;
    LOG("number of particles: " << npart);
}

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::cutoff_radius(float_type value)
{
    r_cut = value;
    LOG("potential cutoff radius: " << r_cut);

    // squared cutoff radius
    rr_cut = std::pow(r_cut, 2);
    // potential energy at cutoff radius
    float_type rri_cut = 1 / rr_cut;
    float_type r6i_cut = rri_cut * rri_cut * rri_cut;
    en_cut = 4 * r6i_cut * (r6i_cut - 1);

    LOG("potential cutoff energy: " << en_cut);
}

#ifdef USE_POTENTIAL_SMOOTHING
template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::potential_smoothing(float_type value)
{
    r_smooth = value;
    LOG("potential smoothing function scale parameter: " << r_smooth);

    // squared inverse potential smoothing function scale parameter
    rri_smooth = std::pow(r_smooth, -2);
}
#endif /* USE_POTENTIAL_SMOOTHING */

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::density(float_type value)
{
    // set particle density
    density_ = value;
    LOG("particle density: " << density_);

    // compute periodic box length
    box_ = std::pow(npart / density_, (float_type) 1 / dimension);
    LOG("periodic simulation box length: " << box_);
}

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::box(float_type value)
{
    // set periodic box length
    box_ = value;
    LOG("periodic simulation box length: " << box_);

    // compute particle density
    density_ = npart / std::pow(box_, (float_type) dimension);
    LOG("particle density: " << density_);
}

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::timestep(float_type value)
{
    // set simulation timestep
    timestep_ = value;
    LOG("simulation timestep: " << timestep_);
}

template <typename ljfluid_impl>
perf_counters ljfluid_base<ljfluid_impl>::times()
{
    perf_counters times(m_times);
    BOOST_FOREACH(perf_counters::value_type& i, m_times) {
	// reset performance counter
	i.second.clear();
    }
    return times;
}

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::attrs(H5::Group const& param) const
{
    H5xx::group node(param.createGroup("mdsim"));
    node["box_length"] = box();
    node["cutoff_radius"] = cutoff_radius();
    node["density"] = density();
    node["dimension"] = (unsigned int) dimension;
    node["particles"] = particles();
#ifdef USE_POTENTIAL_SMOOTHING
    node["potential_smoothing"] = potential_smoothing();
#endif
    node["timestep"] = timestep();
}

} // namespace ljgpu

#endif /* ! LJGPU_LJFLUID_BASE_HPP */
