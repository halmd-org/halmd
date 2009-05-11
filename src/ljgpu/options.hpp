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

#include <boost/mpl/filter_view.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/noncopyable.hpp>
#include <boost/program_options.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <ljgpu/mdsim/impl.hpp>
#include <stdint.h>
#include <string>

namespace ljgpu
{

/**
 * Molecular Dynamics simulation program options
 */
class options : boost::noncopyable
{
public:
    /**
     * build program options for an implementation using its base classes
     */
    template <typename mdsim_impl>
    class description : public boost::program_options::options_description
    {
    public:
	typedef boost::program_options::options_description _Base;

	description() : _Base("MD simulation options")
	{
	    boost::mpl::for_each<
		boost::mpl::transform_view<
		    boost::mpl::filter_view<
			boost::mpl::vector<
			    mdsim_impl_base,
			    ljfluid_impl_base,
			    ljfluid_impl_gpu_base,
			    ljfluid_impl_gpu_square,
			    ljfluid_impl_gpu_cell,
			    ljfluid_impl_gpu_neighbour,
			    ljfluid_impl_host,
			    hardsphere_impl>,
			boost::is_base_of<boost::mpl::_, mdsim_impl> >,
		    options::add<boost::mpl::_> > >(boost::ref(*this));
	}

	template <typename functor>
	void operator()(functor& f)
	{
	    f(*this);
	}
    };

private:
    /**
     * functor to add implementation-specific options
     */
    template <typename mdsim_impl>
    struct add
    {
	void operator()(boost::program_options::options_description& desc);
    };

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
    void parse(boost::program_options::options_description const& opt);

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

} // namespace ljgpu

#endif /* ! LJGPU_OPTIONS_HPP */
