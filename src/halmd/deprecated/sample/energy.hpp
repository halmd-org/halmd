/* Thermodynamic equilibrium properties
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

#ifndef HALMD_SAMPLE_ENERGY_HPP
#define HALMD_SAMPLE_ENERGY_HPP

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/noncopyable.hpp>
#include <string>
#include <vector>

#include <halmd/math/accum.hpp>
#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/mdsim/backend/sample.hpp>
#include <halmd/sample/H5param.hpp>
#include <halmd/util/H5xx.hpp>

#define foreach BOOST_FOREACH

namespace halmd
{

/**
 * Thermodynamic equilibrium properties
 */
template <int dimension>
class energy : boost::noncopyable
{
public:
    /** io flags */
    enum openmode {
        in = 0x1,
        out = 0x2,
    };

    /** time and scalar property */
    typedef std::pair<double, double> scalar_pair;
    /** time and vector property */
    typedef std::pair<double, vector<double, dimension> > vector_pair;

    enum {
        /** HDF5 dataset chunk size */
        CHUNK_SIZE = 5000,
    };

public:
    energy() : m_samples(0), m_samples_buffer(0), m_samples_file(0), m_is_open(false) {}
    /** create HDF5 thermodynamic equilibrium properties output file */
    void open(std::string const& filename, openmode mode = in);
    /** write thermodynamic equilibrium properties to HDF5 file */
    void flush(bool force=true);
    /** close HDF5 file */
    void close();
    /** returns true iff associated with HDF5 file */
    bool is_open() const { return m_is_open; }

    /** sample thermodynamic equilibrium properties */
    void sample(energy_sample<dimension> const& sample, float density, double time);
    /** calculate temperature from given start time */
    void temperature(accumulator<double>& en, double start_time = 0);

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
    /** true iff associated with HDF5 file */
    bool m_is_open;
};

template <int dimension>
void energy<dimension>::sample(energy_sample<dimension> const& sample, float density, double time)
{
    typedef typename energy_sample<dimension>::virial_tensor virial_tensor;

    // mean potential energy per particle
    m_en_pot.push_back(scalar_pair(time, sample.en_pot));
    // mean kinetic energy per particle
    m_en_kin.push_back(scalar_pair(time, sample.vv / 2));
    // mean total energy per particle
    m_en_tot.push_back(scalar_pair(time, m_en_pot.back().second + m_en_kin.back().second));
    // temperature
    m_temp.push_back(scalar_pair(time, sample.vv / dimension));
    // pressure
    double virial = 0;
    foreach(virial_tensor const& vir, sample.virial) {
        virial += vir.front();
    }
    m_press.push_back(scalar_pair(time, density / dimension * virial));
    // velocity center of mass
    m_v_cm.push_back(vector_pair(time, sample.v_cm));

    m_samples++;
    m_samples_buffer++;

    if (m_samples_buffer >= CHUNK_SIZE) {
        // commit full buffers to file
        flush(false);
    }
}

} // namespace halmd

#undef foreach

#endif /* ! HALMD_SAMPLE_ENERGY_HPP */
