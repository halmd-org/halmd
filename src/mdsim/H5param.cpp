/* Molecular Dynamics simulation parameters
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

#include "H5param.hpp"
#include "exception.hpp"
#include "version.h"


namespace mdsim
{

/**
 * read parameters from HDF5 file
 */
void H5param::read(H5::Group const& root)
{
    try {
	H5::Attribute attr;
	H5::Group node;

	// hard spheres simulation parameters
	node = root.openGroup("mdsim");

	// positional coordinates dimension
	attr = node.openAttribute("dimension");
	attr.read(H5::PredType::NATIVE_UINT, &dimension_);
	// number of particles
	attr = node.openAttribute("particles");
	attr.read(H5::PredType::NATIVE_UINT, &particles_);
	// pair separation at which particle collision occurs
	attr = node.openAttribute("pair_separation");
	attr.read(H5::PredType::NATIVE_DOUBLE, &pair_separation_);
	// particle density
	attr = node.openAttribute("density");
	attr.read(H5::PredType::NATIVE_DOUBLE, &density_);
	// periodic box length
	attr = node.openAttribute("box_length");
	attr.read(H5::PredType::NATIVE_DOUBLE, &box_length_);
	// sample timestep
	attr = node.openAttribute("timestep");
	attr.read(H5::PredType::NATIVE_DOUBLE, &timestep_);

	// correlation function parameters
	node = root.openGroup("autocorrelation");

	// number of sample steps
	attr = node.openAttribute("steps");
	attr.read(H5::PredType::NATIVE_UINT64, &steps_);
	// total sample time
	attr = node.openAttribute("time");
	attr.read(H5::PredType::NATIVE_DOUBLE, &time_);
	// block size
	attr = node.openAttribute("block_size");
	attr.read(H5::PredType::NATIVE_UINT, &block_size_);
	// block shift
	attr = node.openAttribute("block_shift");
	attr.read(H5::PredType::NATIVE_UINT, &block_shift_);
	// block count
	attr = node.openAttribute("block_count");
	attr.read(H5::PredType::NATIVE_UINT, &block_count_);
	// maximum number of samples per block
	attr = node.openAttribute("max_samples");
	attr.read(H5::PredType::NATIVE_UINT64, &max_samples_);
	// number of k-valuess
	attr = node.openAttribute("k_values");
	attr.read(H5::PredType::NATIVE_UINT, &k_values_);
    }
    catch (H5::Exception const& e) {
	throw exception("failed to read parameters from HDF5 input file");
    }
}

/**
 * write parameters to HDF5 file
 */
void H5param::write(H5::Group root) const
{
    try {
	H5::Attribute attr;
	H5::Group node;

	// hard spheres simulation parameters
	node = root.createGroup("mdsim");

	// positional coordinates dimension
	attr = node.createAttribute("dimension", H5::PredType::NATIVE_UINT, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_UINT, &dimension_);
	// number of particles
	attr = node.createAttribute("particles", H5::PredType::NATIVE_UINT, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_UINT, &particles_);
	// pair separation at which particle collision occurs
	attr = node.createAttribute("pair_separation", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_DOUBLE, &pair_separation_);
	// particle density
	attr = node.createAttribute("density", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_DOUBLE, &density_);
	// periodic box length
	attr = node.createAttribute("box_length", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_DOUBLE, &box_length_);
	// number of cells per dimension
	attr = node.createAttribute("cells", H5::PredType::NATIVE_UINT, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_UINT, &cells_);
	// cell length
	attr = node.createAttribute("cell_length", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_DOUBLE, &cell_length_);
	// sample timestep
	attr = node.createAttribute("timestep", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_DOUBLE, &timestep_);

	// correlation function parameters
	node = root.createGroup("autocorrelation");

	// number of sample steps
	attr = node.createAttribute("steps", H5::PredType::NATIVE_UINT64, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_UINT64, &steps_);
	// total sample time
	attr = node.createAttribute("time", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_DOUBLE, &time_);
	// block size
	attr = node.createAttribute("block_size", H5::PredType::NATIVE_UINT, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_UINT, &block_size_);
	// block shift
	attr = node.createAttribute("block_shift", H5::PredType::NATIVE_UINT, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_UINT, &block_shift_);
	// block count
	attr = node.createAttribute("block_count", H5::PredType::NATIVE_UINT, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_UINT, &block_count_);
	// maximum number of samples per block
	attr = node.createAttribute("max_samples", H5::PredType::NATIVE_UINT64, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_UINT64, &max_samples_);
	// number of k-values
	attr = node.createAttribute("k_values", H5::PredType::NATIVE_UINT, H5S_SCALAR);
	attr.write(H5::PredType::NATIVE_UINT, &k_values_);

	// program info
	node = root.createGroup("program");

	// program name
	attr = node.createAttribute("name", H5::StrType(H5::PredType::C_S1, 256), H5S_SCALAR);
	attr.write(H5::StrType(H5::PredType::C_S1, 256), PROGRAM_NAME);
	// program version
	attr = node.createAttribute("version", H5::StrType(H5::PredType::C_S1, 256), H5S_SCALAR);
	attr.write(H5::StrType(H5::PredType::C_S1, 256), PROGRAM_VERSION);
	// program variant
	attr = node.createAttribute("variant", H5::StrType(H5::PredType::C_S1, 256), H5S_SCALAR);
	attr.write(H5::StrType(H5::PredType::C_S1, 256), PROGRAM_VARIANT);
    }
    catch (H5::Exception const& e) {
	throw exception("failed to write parameters to HDF5 input file");
    }
}

} // namespace mdsim
