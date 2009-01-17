/* Thermodynamic equilibrium properties
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

#ifndef LJGPU_SAMPLE_ENERGY_HPP
#define LJGPU_SAMPLE_ENERGY_HPP

#include <H5Cpp.h>
#include <algorithm>
#include <boost/foreach.hpp>
#include <ljgpu/math/stat.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/sample/H5param.hpp>
#include <string>
#include <vector>

#define foreach BOOST_FOREACH

namespace ljgpu
{

/**
 * Thermodynamic equilibrium properties
 */
template <int dimension>
class energy
{
public:
    /** time and scalar property */
    typedef std::pair<double, double> scalar_pair;
    /** time and vector property */
    typedef std::pair<double, vector<double, dimension> > vector_pair;

    enum {
	/** HDF5 dataset chunk size */
	CHUNK_SIZE = 5000,
    };

public:
    energy() : m_samples(0), m_samples_buffer(0), m_samples_file(0) {}
    /** create HDF5 thermodynamic equilibrium properties output file */
    void open(std::string const& filename);
    /** write thermodynamic equilibrium properties to HDF5 file */
    void flush(bool force=true);
    /** close HDF5 file */
    void close();

    /** sample thermodynamic equilibrium properties */
    template <typename sample_type>
    void sample(sample_type const& sample, float density, double time);

    /** returns HDF5 parameter group */
    operator H5param() { return m_file; }
private:
    /** number of aquired samples */
    unsigned int m_samples;
    /** number of samples in memory */
    uint64_t m_samples_buffer;
    /** number of samples committed to file */
    uint64_t m_samples_file;

    /** thermodynamic equilibrium properties */
    std::vector<scalar_pair> m_en_pot;
    std::vector<scalar_pair> m_en_kin;
    std::vector<scalar_pair> m_en_tot;
    std::vector<scalar_pair> m_temp;
    std::vector<scalar_pair> m_press;
    std::vector<vector_pair> m_v_cm;
    /** HDF5 thermodynamic equilibrium properties output file */
    H5::H5File m_file;
    /** HDF5 datasets */
    boost::array<H5::DataSet, 6> m_dataset;
    /** HDF5 floating-point data type */
    H5::DataType m_tid;
};

/**
 * sample thermodynamic equilibrium properties
 */
template <int dimension>
template <typename sample_type>
void energy<dimension>::sample(sample_type const& sample, float density, double time)
{
    typedef typename sample_type::uniform_sample uniform_sample;
    // ensure double-precision floating point arithmetic
    typedef vector<double, dimension> velocity_vector;

    vector<double, dimension> v_cm = 0;
    double vv = 0;
    size_t i = 0;
    foreach (uniform_sample const& sample_, sample) {
	foreach (velocity_vector v, sample_.v) {
	    // center of mass velocity
	    v_cm += (v - v_cm) / ++i;
	    // mean squared velocity
	    vv += (v * v - vv) / i;
	}
    }

    // mean potential energy per particle
    m_en_pot.push_back(scalar_pair(time, sample.en_pot));
    // mean kinetic energy per particle
    m_en_kin.push_back(scalar_pair(time, vv / 2));
    // mean total energy per particle
    m_en_tot.push_back(scalar_pair(time, m_en_pot.back().second + m_en_kin.back().second));
    // temperature
    m_temp.push_back(scalar_pair(time, vv / dimension));
    // pressure
    m_press.push_back(scalar_pair(time, density / dimension * (vv + sample.virial)));
    // velocity center of mass
    m_v_cm.push_back(vector_pair(time, v_cm));

    m_samples++;
    m_samples_buffer++;

    if (m_samples_buffer >= CHUNK_SIZE) {
	// commit full buffers to file
	flush(false);
    }
}

} // namespace ljgpu

#undef foreach

#endif /* ! LJGPU_SAMPLE_ENERGY_HPP */
