/* Thermodynamic equilibrium properties
 *
 * Copyright (C) 2008  Peter Colberg
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

#ifndef MDSIM_ENERGY_HPP
#define MDSIM_ENERGY_HPP

#include <H5Cpp.h>
#include <algorithm>
#include <string>
#include <vector>
#include "H5param.hpp"
#include "config.hpp"

namespace mdsim
{

/**
 * Thermodynamic equilibrium properties
 */
class energy
{
public:
    /** time and scalar property */
    typedef std::pair<double, double> scalar_pair;
    /** time and vector property */
    typedef std::pair<double, hvector> vector_pair;

    enum {
	/** HDF5 dataset chunk size */
	CHUNK_SIZE = 2500,
    };

public:
    energy() : m_samples(0), m_samples_buffer(0), m_samples_file(0) {}
    /** create HDF5 thermodynamic equilibrium properties output file */
    void open(std::string const& filename);
    /** returns HDF5 parameter group */
    H5param attrs();
    /** sample thermodynamic equilibrium properties */
    void sample(std::vector<hvector> const& v, double const& en_pot, double const& virial, double const& density, double const& time);
    /** write thermodynamic equilibrium properties to HDF5 file */
    void flush(bool force=true);
    /** close HDF5 thermodynamic equilibrium properties output file */
    void close();

private:
    /** number of aquired samples */
    uint64_t m_samples;
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

} // namespace mdsim

#endif /* ! MDSIM_ENERGY_HPP */
