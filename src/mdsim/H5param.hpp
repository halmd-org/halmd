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

#ifndef MDSIM_H5PARAM_HPP
#define MDSIM_H5PARAM_HPP

#include <H5Cpp.h>
#include <stdint.h>
#include <string>


namespace mdsim
{

/**
 * Molecular Dynamics simulation parameters
 */
class H5param
{
public:
    /** read parameters from HDF5 file */
    void read(H5::Group const& root);
    /** write parameters to HDF5 file */
    void write(H5::Group root) const;

    /**
     * returns positional coordinates dimension
     */
    unsigned int const& dimension() const
    {
	return dimension_;
    }

    /**
     * set positional coordinates dimension
     */
    void dimension(unsigned int const& value)
    {
	dimension_ = value;
    }

    /**
     * returns number of particles
     */
    unsigned int const& particles() const
    {
	return particles_;
    }

    /**
     * set number of particles
     */
    void particles(unsigned int const& value)
    {
	particles_ = value;
    }

    /**
     * set number of cells per dimension
     */
    void cells(unsigned int const& value)
    {
	cells_ = value;
    }

    /**
     * returns number of CUDA execution threads
     */
    unsigned int const& threads() const
    {
	return threads_;
    }

    /**
     * set number of CUDA execution threads
     */
    void threads(unsigned int const& value)
    {
	threads_ = value;
    }

    /**
     * returns number of CUDA execution blocks
     */
    unsigned int const& blocks() const
    {
	return blocks_;
    }

    /**
     * set number of CUDA execution blocks
     */
    void blocks(unsigned int const& value)
    {
	blocks_ = value;
    }

    /**
     * returns particle density
     */
    double const& density() const
    {
	return density_;
    }

    /**
     * set particle density
     */
    void density(double const& value)
    {
	density_ = value;
    }

    /**
     * returns periodic box length
     */
    double const& box_length() const
    {
	return box_length_;
    }

    /**
     * set periodic box length
     */
    void box_length(double const& value)
    {
	box_length_ = value;
    }

    /**
     * set cell length
     */
    void cell_length(double const& value)
    {
	cell_length_ = value;
    }

    /**
     * returns simulation timestep
     */
    double const& timestep() const
    {
	return timestep_;
    }

    /**
     * returns simulation timestep
     */
    void timestep(double const& value)
    {
	timestep_ = value;
    }

    /**
     * returns number of simulation steps
     */
    uint64_t const& steps() const
    {
	return steps_;
    }

    /**
     * set number of simulation steps
     */
    void steps(uint64_t const& value)
    {
	steps_ = value;
    }

    /**
     * returns total simulation time
     */
    double const& time() const
    {
	return time_;
    }

    /**
     * returns total simulation time
     */
    void time(double const& value)
    {
	time_ = value;
    }

    /**
     * returns block size
     */
    unsigned int const& block_size() const
    {
	return block_size_;
    }

    /**
     * set block size
     */
    void block_size(unsigned int const& value)
    {
	block_size_ = value;
    }

    /**
     * returns block shift
     */
    unsigned int const& block_shift() const
    {
	return block_shift_;
    }

    /**
     * set block shift
     */
    void block_shift(unsigned int const& value)
    {
	block_shift_ = value;
    }

    /**
     * returns block count
     */
    unsigned int const& block_count() const
    {
	return block_count_;
    }

    /**
     * set block count
     */
    void block_count(unsigned int const& value)
    {
	block_count_ = value;
    }

    /**
     * returns maximum number of samples per block
     */
    uint64_t const& max_samples() const
    {
	return max_samples_;
    }

    /**
     * set maximum number of samples per block
     */
    void max_samples(uint64_t const& value)
    {
	max_samples_ = value;
    }

    /**
     * returns number of q-values
     */
    unsigned int const& q_values() const
    {
	return q_values_;
    }

    /**
     * set number of q-values
     */
    void q_values(unsigned int const& value)
    {
	q_values_ = value;
    }

    /**
     * returns Lennard-Jones potential cutoff distance
     */
    double const& cutoff_distance() const
    {
	return cutoff_distance_;
    }

    /**
     * set Lennard-Jones potential cutoff distance
     */
    void cutoff_distance(double const& value)
    {
	cutoff_distance_ = value;
    }

private:
    /** positional coordinates dimension */
    unsigned int dimension_;
    /** number of particles */
    unsigned int particles_;
    /** number of cells per dimension */
    unsigned int cells_;
    /** number of CUDA execution threads */
    unsigned int threads_;
    /** number of CUDA execution blocks */
    unsigned int blocks_;
    /** particle density */
    double density_;
    /** periodic box length */
    double box_length_;
    /** cell length */
    double cell_length_;
    /** simulation timestep */
    double timestep_;
    /** number of simulation steps */
    uint64_t steps_;
    /** total simulation time */
    double time_;
    /** block size */
    unsigned int block_size_;
    /** block shift */
    unsigned int block_shift_;
    /** block count */
    unsigned int block_count_;
    /** maximum number of samples per block */
    uint64_t max_samples_;
    /** number of q-values */
    unsigned int q_values_;
    /** Lennard-Jones potential cutoff distance */
    double cutoff_distance_;
};

} // namespace mdsim

#endif /* ! MDSIM_H5PARAM_HPP */
