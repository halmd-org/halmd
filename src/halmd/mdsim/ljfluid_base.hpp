/* Lennard-Jones fluid simulation
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

#ifndef HALMD_MDSIM_LJFLUID_BASE_HPP
#define HALMD_MDSIM_LJFLUID_BASE_HPP

#include <boost/foreach.hpp>
#include <halmd/mdsim/base.hpp>
#include <halmd/mdsim/traits.hpp>

namespace halmd
{

/**
 * Lennard-Jones fluid interface
 */
template <typename ljfluid_impl, int dimension>
class ljfluid_base : public mdsim_base<ljfluid_impl, dimension>
{
public:
    typedef mdsim_base<ljfluid_impl, dimension> _Base;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::host_sample_type host_sample_type;
    typedef typename _Base::energy_sample_type energy_sample_type;

public:
    ljfluid_base() :
        epsilon_(boost::assign::list_of(1)(0)(0)),
        sigma_(boost::assign::list_of(1)(0)(0)),
        r_smooth(0),
        thermostat_nu(0),
        thermostat_steps(0),
        thermostat_count(0),
        potential_(C0POT) {}

    /** set simulation timestep */
    void timestep(double value);
    /** set potential cutoff radius */
    void cutoff_radius(boost::array<float, 3> const& value);
    /** set potential smoothing function scale parameter */
    void potential_smoothing(float_type value);
    /** set heat bath collision probability and temperature */
    void thermostat(float_type nu, float_type temp);
    /** set potential well depths */
    void epsilon(boost::array<float, 3> const& value);
    /** set collision diameters */
    void sigma(boost::array<float, 3> const& value);

    /** returns simulation timestep */
    double timestep() const { return timestep_; }
    /** returns potential cutoff radius */
    boost::array<float_type, 3> const& cutoff_radius() const { return r_cut_sigma; }
    /** returns potential smoothing function scale parameter */
    float_type potential_smoothing() const { return r_smooth; }

    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

    potential_type potential() const { return potential_ ; }

protected:
    using _Base::npart;
    using _Base::mpart;
    using _Base::box_;
    using _Base::density_;
    using _Base::mixture_;

    /** cutoff radius in units of sigma */
    boost::array<float_type, 3> r_cut_sigma;
    /** cutoff radii in binary mixture */
    boost::array<float_type, 3> r_cut;
    /** squared cutoff radii */
    boost::array<float_type, 3> rr_cut;
    /** potential well depths */
    boost::array<float_type, 3> epsilon_;
    /** collision diameters */
    boost::array<float_type, 3> sigma_;
    /** squared collision diameters */
    boost::array<float_type, 3> sigma2_;
    /** Lennard-Jones potential at cutoff radius in units of epsilon */
    boost::array<float_type, 3> en_cut;
    /** potential smoothing function scale parameter */
    float_type r_smooth;
    /** squared inverse potential smoothing function scale parameter */
    float_type rri_smooth;
    /** simulation timestep */
    double timestep_;
    /** heat bath collision probability */
    float_type thermostat_nu;
    /** heat bath coupling frequency */
    unsigned int thermostat_steps;
    /** MD steps since last heat bath coupling */
    unsigned int thermostat_count;
    /** heat bath temperature */
    float_type thermostat_temp;

    /** C⁰ or C²-smooth potential */
    potential_type potential_;
};

template <typename ljfluid_impl, int dimension>
void ljfluid_base<ljfluid_impl, dimension>::epsilon(boost::array<float, 3> const& value)
{
    epsilon_ = value;
    LOG("potential well depths: ε(AA) = " << epsilon_[0] << ", ε(AB) = " << epsilon_[1] << ", ε(BB) = " << epsilon_[2]);
}

template <typename ljfluid_impl, int dimension>
void ljfluid_base<ljfluid_impl, dimension>::sigma(boost::array<float, 3> const& value)
{
    sigma_ = value;
    LOG("collision diameters: σ(AA) = " << sigma_[0] << ", σ(AB) = " << sigma_[1] << ", σ(BB) = " << sigma_[2]);

    for (size_t i = 0; i < sigma_.size(); ++i) {
        sigma2_[i] = std::pow(sigma_[i], 2);
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_base<ljfluid_impl, dimension>::cutoff_radius(boost::array<float, 3> const& value)
{
    r_cut_sigma = value;

    for (size_t i = 0; i < r_cut_sigma.size(); ++i) {
        float_type rri_cut = 1 / std::pow(r_cut_sigma[i], 2);
        float_type r6i_cut = rri_cut * rri_cut * rri_cut;
        en_cut[i] = 4 * r6i_cut * (r6i_cut - 1);
    }

    if (mixture_ == BINARY) {
        LOG("potential cutoff radii: r(AA) = " << r_cut_sigma[0] << ", r(AB) = " << r_cut_sigma[1] << ", r(BB) = " << r_cut_sigma[2]);
        LOG("potential cutoff energies: E(AA) = " << en_cut[0] << ", E(AB) = " << en_cut[1] << ", E(BB) = " << en_cut[2]);
    }
    else {
        LOG("potential cutoff radius: " << r_cut_sigma[0]);
        LOG("potential cutoff energy: " << en_cut[0]);
    }

    for (size_t i = 0; i < sigma_.size(); ++i) {
        r_cut[i] = r_cut_sigma[i] * sigma_[i];
        rr_cut[i] = std::pow(r_cut[i], 2);
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_base<ljfluid_impl, dimension>::potential_smoothing(float_type value)
{
    r_smooth = value;
    potential_ = C2POT;
    LOG("potential smoothing function scale parameter: " << r_smooth);

    // squared inverse potential smoothing function scale parameter
    rri_smooth = std::pow(r_smooth, -2);
}

template <typename ljfluid_impl, int dimension>
void ljfluid_base<ljfluid_impl, dimension>::timestep(double value)
{
    // set simulation timestep
    timestep_ = value;
    LOG("simulation timestep: " << timestep_);
}

template <typename ljfluid_impl, int dimension>
void ljfluid_base<ljfluid_impl, dimension>::thermostat(float_type nu, float_type temp)
{
    thermostat_nu = nu;
    LOG("heat bath collision probability: " << thermostat_nu);
    thermostat_steps = static_cast<unsigned int>(std::max(round(1 / (nu * timestep_)), 1.));
    LOG("heat bath coupling frequency: " << thermostat_steps);
    thermostat_temp = temp;
    LOG("heat bath temperature: " << thermostat_temp);
}

template <typename ljfluid_impl, int dimension>
void ljfluid_base<ljfluid_impl, dimension>::param(H5param& param) const
{
    _Base::param(param);

    H5xx::group node(param["mdsim"]);
    if (mixture_ == BINARY) {
        node["potential_epsilon"] = epsilon_;
        node["potential_sigma"] = sigma_;
    }
    node["cutoff_radius"] = r_cut_sigma;
    node["timestep"] = timestep_;
    if (potential_ == C2POT) {
        node["potential_smoothing"] = r_smooth;
    }
    if (thermostat_steps) {
        node["thermostat_nu"] = thermostat_nu;
        node["thermostat_steps"] = thermostat_steps;
        node["thermostat_temp"] = thermostat_temp;
    }
}

} // namespace halmd

#endif /* ! HALMD_MDSIM_LJFLUID_BASE_HPP */
