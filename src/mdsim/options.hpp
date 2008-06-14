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
#include <stdint.h>
#include <string>


namespace mdsim
{

/**
 * option value
 */
template <typename T>
class option_value : public boost::program_options::variable_value
{
public:
    option_value(boost::program_options::variable_value const& vv) : boost::program_options::variable_value(vv) {}

    /**
     * returns the contained value
     */
    T const& value() const
    {
	return boost::program_options::variable_value::as<T>();
    }

    T& value()
    {
	return boost::program_options::variable_value::as<T>();
    }
};

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
    options() {}
    /** parse program option values */
    void parse(int argc, char** argv);

    /**
     * returns number of particles
     */
    option_value<unsigned int> particles() const
    {
	return vm["particles"];
    }

    /**
     * returns particle density
     */
    option_value<double> density() const
    {
	return vm["density"];
    }

    /**
     * returns simulation box length
     */
    option_value<double> box_length() const
    {
	return vm["box-length"];
    }

    /**
     * returns simulation timestep
     */
    option_value<double> timestep() const
    {
	return vm["timestep"];
    }

    /**
     * returns initial system temperature
     */
    option_value<double> temperature() const
    {
	return vm["temperature"];
    }

    /**
     * returns number of simulation steps
     */
    option_value<uint64_t> steps() const
    {
	return vm["steps"];
    }

#ifndef BENCHMARK
    /**
     * returns total simulation time
     */
    option_value<double> time() const
    {
	return vm["time"];
    }
#endif

#ifndef BENCHMARK
    /**
     * returns block size
     */
    option_value<unsigned int> block_size() const
    {
	return vm["block-size"];
    }

    /**
     * returns maximum number of samples per block
     */
    option_value<uint64_t> max_samples() const
    {
	return vm["max-samples"];
    }
#endif

    /**
     * returns output file prefix
     */
    option_value<std::string> output_file_prefix() const
    {
	return vm["output"];
    }

    /**
     * returns random number generator seed
     */
    option_value<unsigned int> rng_seed() const
    {
	return vm["rng-seed"];
    }

    /**
     * returns verbosity
     */
    option_value<int> verbosity() const
    {
	return vm["verbose"];
    }

    /**
     * returns trajectory input file
     */
    option_value<std::string> trajectory_input_file() const
    {
	return vm["trajectory"];
    }

    /**
     * returns sample in trajectory input file
     */
    option_value<int64_t> trajectory_sample() const
    {
	return vm["sample"];
    }

    /**
     * returns whether to perform a trial run without simulation
     */
    option_value<bool> dry_run() const
    {
	return vm["dry-run"];
    }

private:
    /** parsed program options */
    boost::program_options::variables_map vm;
};

} // namespace mdsim

#endif /* ! MDSIM_OPTIONS_HPP */
