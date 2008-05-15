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

    double density() const
    {
	return density_;
    }

    double timestep() const
    {
	return timestep_;
    }

    double temp() const
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

    char const* tcf_output() const
    {
	return tcf_output_.c_str();
    }

    unsigned int rngseed() const
    {
	return rngseed_;
    }

    char const* output() const
    {
	return output_.c_str();
    }

private:
    /** number of particles */
    uint64_t npart_;
    /** density */
    double density_;
    /** simulation timestep */
    double timestep_;
    /** initial temperature */
    double temp_;
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
    /** trajectory correlation functions output file */
    std::string tcf_output_;

    /** random number generator seed */
    unsigned int rngseed_;
    /** trajectory output file */
    std::string output_;
};

} // namespace mdsim

#endif /* ! MDSIM_OPTIONS_HPP */
