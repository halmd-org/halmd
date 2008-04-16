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

#include <sys/types.h>


namespace mdsim
{

/**
 * Molecular Dynamics simulation program options
 */
class options
{
public:
    class exception
    {
    public:
	exception(int status) : status_(status) {}

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

    size_t npart() const
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

    size_t steps() const
    {
	return steps_;
    }

    size_t avgsteps() const
    {
	return avgsteps_;
    }

    unsigned int rngseed() const
    {
	return rngseed_;
    }

private:
    /** number of particles */
    size_t npart_;
    /** density */
    double density_;
    /** simulation timestep */
    double timestep_;
    /** initial temperature */
    double temp_;
    /** number of simulation steps */
    size_t steps_;
    /** number of steps for average accumulation */
    size_t avgsteps_;
    /** random number generator seed */
    unsigned int rngseed_;
};

} // namespace mdsim

#endif /* ! MDSIM_OPTIONS_HPP */
