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

#include <boost/program_options.hpp>
#include <cuda_wrapper.hpp>
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

    uint64_t const& npart() const
    {
	return npart_;
    }

    float const& density() const
    {
	return density_;
    }

    float const& timestep() const
    {
	return timestep_;
    }

    float const& temp() const
    {
	return temp_;
    }

    uint64_t const& steps() const
    {
	return steps_;
    }

    unsigned int const& block_size() const
    {
	return block_size_;
    }

    uint64_t const& max_samples() const
    {
	return max_samples_;
    }

    std::string correlations_output_file() const
    {
	return output_file_prefix_ + ".tcf";
    }

    unsigned short const& device() const
    {
	return device_;
    }

    cuda::config dim() const
    {
	return cuda::config(dim3((npart_ + threads_ - 1) / threads_), dim3(threads_));
    }

    unsigned int const& rngseed() const
    {
	return rngseed_;
    }

    std::string trajectory_output_file() const
    {
	return output_file_prefix_ + ".trj";
    }

    std::string energy_output_file() const
    {
	return output_file_prefix_ + ".tep";
    }

    std::string logfile() const
    {
	return output_file_prefix_ + ".log";
    }

    /**
     * returns verbosity
     */
    int verbosity() const
    {
	if (vm.count("verbose")) {
	    return vm["verbose"].as<int>();
	}
	return 0;
    }

    std::string const& trajectory_input_file() const
    {
	return trajectory_input_file_;
    }

    int const& sample() const
    {
	return sample_;
    }

private:
    /** parsed program options */
    boost::program_options::variables_map vm;

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

    /** block size */
    unsigned int block_size_;
    /** maximum number of samples per block */
    uint64_t max_samples_;

    /** CUDA device */
    unsigned short device_;
    /** number of threads per block */
    unsigned int threads_;

    /** random number generator seed */
    unsigned int rngseed_;
    /** output file prefix */
    std::string output_file_prefix_;
    /** trajectory input file */
    std::string trajectory_input_file_;
    /** sample in trajectory input file */
    int sample_;
};

} // namespace mdsim

#endif /* ! MDSIM_OPTIONS_HPP */
