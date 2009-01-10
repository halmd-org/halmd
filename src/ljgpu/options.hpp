/* Molecular Dynamics simulation program options
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

#ifndef LJGPU_OPTIONS_HPP
#define LJGPU_OPTIONS_HPP

#include <boost/program_options.hpp>
#include <ljgpu/ljfluid/impl.hpp>
#include <stdint.h>
#include <string>

namespace ljgpu
{

/**
 * Molecular Dynamics simulation program options
 */
class options
{
public:
    typedef boost::program_options::options_description description;

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
    /** parse program option values */
    void parse(int argc, char** argv);
    /** parse backend option values */
    void parse(description const& opt);

    /**
     * return option value
     */
    boost::program_options::variable_value const& operator[](std::string const& vv) const
    {
	return vm[vv];
    }

private:
    /** parsed program options */
    boost::program_options::variables_map vm;
    /** unrecognised program options */
    std::vector<std::string> unparsed;
};

struct mdsim_options_description : public boost::program_options::options_description
{
    typedef boost::program_options::options_description _Base;
    mdsim_options_description();
};

template <template <int> class backend>
struct options_description : public mdsim_options_description
{
    options_description();
};

} // namespace ljgpu

#endif /* ! LJGPU_OPTIONS_HPP */
