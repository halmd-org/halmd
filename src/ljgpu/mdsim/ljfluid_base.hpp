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

#ifndef LJGPU_MDSIM_LJFLUID_BASE_HPP
#define LJGPU_MDSIM_LJFLUID_BASE_HPP

#include <boost/foreach.hpp>
#include <ljgpu/mdsim/base.hpp>
#include <ljgpu/mdsim/traits.hpp>

namespace ljgpu
{

/**
 * Lennard-Jones fluid interface
 */
template <typename ljfluid_impl>
class ljfluid_base : public mdsim_base<ljfluid_impl>
{
public:
    typedef mdsim_base<ljfluid_impl> _Base;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::sample_type sample_type;
    typedef typename sample_type::sample_visitor sample_visitor;
    enum { dimension = _Base::dimension };

public:
    ljfluid_base();
    /** set simulation timestep */
    void timestep(float_type value);
    /** set potential cutoff radius */
    void cutoff_radius(float_type value);
    /** set potential smoothing function scale parameter */
    void potential_smoothing(float_type value);
    /** set heat bath collision probability and temperature */
    void thermostat(float_type nu, float_type temp);

    /** returns simulation timestep */
    float_type timestep() const { return timestep_; }
    /** returns potential cutoff radius */
    float_type cutoff_radius() const { return r_cut; }
    /** returns potential smoothing function scale parameter */
    float_type potential_smoothing() const { return r_smooth; }

    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

protected:
    using _Base::npart;
    using _Base::box_;
    using _Base::density_;

    /** simulation timestep */
    float_type timestep_;
    /** cutoff radius for shifted Lennard-Jones potential */
    float_type r_cut;
    /** squared cutoff radius */
    float_type rr_cut;
    /** potential energy at cutoff radius */
    float_type en_cut;
    /** potential smoothing function scale parameter */
    float_type r_smooth;
    /** squared inverse potential smoothing function scale parameter */
    float_type rri_smooth;
    /** heat bath collision probability */
    float_type thermostat_nu;
    /** heat bath temperature */
    float_type thermostat_temp;
};

template <typename ljfluid_impl>
ljfluid_base<ljfluid_impl>::ljfluid_base() : r_smooth(0), thermostat_nu(0)
{
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

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::potential_smoothing(float_type value)
{
    r_smooth = value;
    LOG("potential smoothing function scale parameter: " << r_smooth);

    // squared inverse potential smoothing function scale parameter
    rri_smooth = std::pow(r_smooth, -2);
}

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::timestep(float_type value)
{
    // set simulation timestep
    timestep_ = value;
    LOG("simulation timestep: " << timestep_);
}

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::thermostat(float_type nu, float_type temp)
{
    thermostat_nu = nu;
    LOG("heat bath collision probability: " << thermostat_nu);

    thermostat_temp = temp;
    LOG("heat bath temperature: " << thermostat_temp);
}

template <typename ljfluid_impl>
void ljfluid_base<ljfluid_impl>::param(H5param& param) const
{
    _Base::param(param);

    H5xx::group node(param["mdsim"]);
    node["cutoff_radius"] = r_cut;
    node["timestep"] = timestep_;
    if (r_smooth > 0) {
	node["potential_smoothing"] = r_smooth;
    }
    if (thermostat_nu > 0) {
	node["thermostat_nu"] = thermostat_nu;
	node["thermostat_temp"] = thermostat_temp;
    }
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_BASE_HPP */
