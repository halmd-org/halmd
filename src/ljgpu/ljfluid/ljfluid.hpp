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

#ifndef LJGPU_LJFLUID_LJFLUID_HPP
#define LJGPU_LJFLUID_LJFLUID_HPP

#include <ljgpu/ljfluid/ljfluid_gpu.hpp>
#include <ljgpu/ljfluid/ljfluid_host.hpp>

namespace ljgpu
{

/**
 * Lennard-Jones fluid interface
 */
template <template<int> class ljfluid_impl, int dimension>
class ljfluid_base : public ljfluid_impl<dimension>
{
public:
    typedef ljfluid_impl<dimension> _Base;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::trajectory_sample trajectory_sample;

public:
    ljfluid_base(options const& opt);

    using _Base::box;
    using _Base::cutoff_radius;
    using _Base::density;
    using _Base::particles;
#ifdef USE_POTENTIAL_SMOOTHING
    using _Base::potential_smoothing;
#endif
    using _Base::timestep;

    /** returns and resets CPU or GPU time accumulators */
    perf_counters times();
    /** get trajectory sample */
    trajectory_sample const& trajectory() const { return m_sample; }
    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

protected:
    using _Base::m_sample;
    using _Base::m_times;
};

template <template<int> class ljfluid_impl, int dimension>
ljfluid_base<ljfluid_impl, dimension>::ljfluid_base(options const& opt)
{
    LOG("positional coordinates dimension: " << dimension);

    particles(opt["particles"].as<unsigned int>());
    if (opt["density"].defaulted() && !opt["box-length"].empty()) {
	box(opt["box-length"].as<float>());
    }
    else {
	density(opt["density"].as<float>());
    }
    cutoff_radius(opt["cutoff"].as<float>());
#ifdef USE_POTENTIAL_SMOOTHING
    potential_smoothing(opt["smoothing"].as<float>());
#endif
    timestep(opt["timestep"].as<float>());
}

template <template<int> class ljfluid_impl, int dimension>
void ljfluid_base<ljfluid_impl, dimension>::attrs(H5::Group const& param) const
{
    H5xx::group node(param.createGroup("mdsim"));
    node["box_length"] = box();
    node["cutoff_radius"] = cutoff_radius();
    node["density"] = density();
    node["dimension"] = dimension;
    node["particles"] = particles();
#ifdef USE_POTENTIAL_SMOOTHING
    node["potential_smoothing"] = potential_smoothing();
#endif
    node["timestep"] = timestep();

    // implementation-dependent attributes
    ljfluid_impl<dimension>::attrs(param);
}

template <template<int> class ljfluid_impl, int dimension>
perf_counters ljfluid_base<ljfluid_impl, dimension>::times()
{
    perf_counters times(m_times);
    BOOST_FOREACH(perf_counters::value_type& i, m_times) {
	// reset performance counter
	i.second.clear();
    }
    return times;
}

/**
 * specialised Lennard-Jones fluid interface
 */
template<template<int> class ljfluid_impl, int dimension>
class ljfluid;

template<int dimension>
class ljfluid<ljfluid_host, dimension> : public ljfluid_base<ljfluid_host, dimension>
{
public:
    typedef ljfluid_base<ljfluid_host, dimension> _Base;

public:
    ljfluid(options const& opt) : _Base(opt)
    {
	init_cell();
    }

    void mdstep() {}

    void sample() {}

    void synchronize()
    {
	_Base::mdstep();
    }

    using _Base::init_cell;
};

template<int dimension>
class ljfluid<ljfluid_gpu_square, dimension> : public ljfluid_base<ljfluid_gpu_square, dimension>
{
public:
    typedef ljfluid_base<ljfluid_gpu_square, dimension> _Base;

public:
    ljfluid(options const& opt) : _Base(opt)
    {
	threads(opt["threads"].as<unsigned int>());
    }

    using _Base::threads;
};

template<int dimension>
class ljfluid<ljfluid_gpu_cell, dimension> : public ljfluid_base<ljfluid_gpu_cell, dimension>
{
public:
    typedef ljfluid_base<ljfluid_gpu_cell, dimension> _Base;

public:
    ljfluid(options const& opt) : _Base(opt)
    {
	cell_occupancy(opt["cell-occupancy"].as<float>());
	threads(opt["threads"].as<unsigned int>());
    }

    using _Base::cell_occupancy;
    using _Base::threads;
};

template<int dimension>
class ljfluid<ljfluid_gpu_neighbour, dimension> : public ljfluid_base<ljfluid_gpu_neighbour, dimension>
{
public:
    typedef ljfluid_base<ljfluid_gpu_neighbour, dimension> _Base;

public:
    ljfluid(options const& opt) : _Base(opt)
    {
	cell_occupancy(opt["cell-occupancy"].as<float>());
	threads(opt["threads"].as<unsigned int>());
    }

    using _Base::cell_occupancy;
    using _Base::threads;
};

} // namespace ljgpu

#endif /* ! LJGPU_LJFLUID_LJFLUID_HPP */
