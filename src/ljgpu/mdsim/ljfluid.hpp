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

#ifndef LJGPU_MDSIM_LJFLUID_HPP
#define LJGPU_MDSIM_LJFLUID_HPP

#include <ljgpu/mdsim/ljfluid_gpu.hpp>
#include <ljgpu/mdsim/ljfluid_host.hpp>

namespace ljgpu
{

template <template<int> class ljfluid_impl, int dimension>
class ljfluid : public ljfluid_impl<dimension>
{
public:
    typedef typename ljfluid_impl<dimension>::float_type float_type;
    typedef typename ljfluid_impl<dimension>::vector_type vector_type;
    typedef typename ljfluid_impl<dimension>::trajectory_sample trajectory_sample;

public:
    /** returns and resets CPU or GPU time accumulators */
    perf_counters times();
    /** get trajectory sample */
    trajectory_sample const& trajectory() const { return m_sample; }

    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

protected:
    using ljfluid_impl<dimension>::m_sample;
    using ljfluid_impl<dimension>::m_times;

    using ljfluid_impl<dimension>::npart;
    using ljfluid_impl<dimension>::density_;
    using ljfluid_impl<dimension>::box_;
    using ljfluid_impl<dimension>::timestep_;
    using ljfluid_impl<dimension>::r_cut;
#ifdef USE_POTENTIAL_SMOOTHING
    using ljfluid_impl<dimension>::r_smooth;
#endif
};

template <template<int> class ljfluid_impl, int dimension>
void ljfluid<ljfluid_impl, dimension>::attrs(H5::Group const& param) const
{
    H5xx::group node(param.createGroup("mdsim"));
    node["dimension"] = dimension;
    node["particles"] = npart;
    node["density"] = density_;
    node["box_length"] = box_;
    node["timestep"] = timestep_;
    node["cutoff_radius"] = r_cut;
#ifdef USE_POTENTIAL_SMOOTHING
    node["potential_smoothing"] = r_smooth;
#endif

    // implementation-dependent attributes
    ljfluid_impl<dimension>::attrs(param);
}

template <template<int> class ljfluid_impl, int dimension>
perf_counters ljfluid<ljfluid_impl, dimension>::times()
{
    perf_counters times(m_times);
    BOOST_FOREACH(perf_counters::value_type& i, m_times) {
	// reset performance counter
	i.second.clear();
    }
    return times;
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_HPP */
