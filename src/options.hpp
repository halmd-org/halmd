/* Molecular Dynamics simulation program options
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

#ifndef MDSIM_OPTIONS_HPP
#define MDSIM_OPTIONS_HPP

#include <stdint.h>
#include <string>


namespace mdsim
{

/**
 * Molecular Dynamics simulation program options
 */
class options
{
public:
    class exit_exception
    {
    public:
	exit_exception(int status) : status_(status) {}

	int status() const
	{
	    return status_;
	}

    private:
	int status_;
    };

public:
    options();

    void parse(int argc, char** argv);

    uint64_t npart() const
    {
	return npart_;
    }

    float density() const
    {
	return density_;
    }

    float timestep() const
    {
	return timestep_;
    }

    float temp() const
    {
	return temp_;
    }

    uint64_t steps() const
    {
	return steps_;
    }

    uint64_t avgsteps() const
    {
	return avgsteps_;
    }

    unsigned int block_count() const
    {
	return block_count_;
    }

    unsigned int block_size() const
    {
	return block_size_;
    }

    unsigned int block_shift() const
    {
	return block_shift_;
    }

    uint64_t max_samples() const
    {
	return max_samples_;
    }

    std::string correlations_output_file() const
    {
	return output_file_prefix_ + ".tcf";
    }

    unsigned int device() const
    {
	return device_;
    }

    unsigned int blocks() const
    {
	return blocks_;
    }

    unsigned int threads() const
    {
	return threads_;
    }

    unsigned int rngseed() const
    {
	return rngseed_;
    }

    std::string trajectory_output_file() const
    {
	return output_file_prefix_ + ".trj";
    }

    bool quiet() const
    {
	return quiet_;
    }

private:
    /** number of particles */
    uint64_t npart_;
    /** density */
    float density_;
    /** simulation timestep */
    float timestep_;
    /** initial temperature */
    float temp_;
    /** number of simulation steps */
    uint64_t steps_;
    /** number of steps for average accumulation */
    uint64_t avgsteps_;

    /** block count */
    unsigned int block_count_;
    /** block size */
    unsigned int block_size_;
    /** block shift */
    unsigned int block_shift_;
    /** maximum number of samples per block */
    uint64_t max_samples_;

    /** CUDA device */
    unsigned int device_;
    /** number of blocks per grid */
    unsigned int blocks_;
    /** number of threads per block */
    unsigned int threads_;

    /** random number generator seed */
    unsigned int rngseed_;
    /** output file prefix */
    std::string output_file_prefix_;
    /** suppress status output */
    bool quiet_;
};

} // namespace mdsim

#endif /* ! MDSIM_OPTIONS_HPP */
