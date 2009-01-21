/* Molecular Dynamics simulation of a Lennard-Jones fluid
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

#include <ljgpu/options.hpp>
#include <ljgpu/util/dl.hpp>

namespace ljgpu
{

struct mdlib : public dl::library
{
    void open(std::string const& name)
    {
	dl::library::open(name);
	options.set(*this, "mdlib_options");
	mdsim.set(*this, "mdlib_mdsim");
	version.set(*this, "mdlib_version");
    }

    dl::symbol<options::description ()> options;
    dl::symbol<void (ljgpu::options const&)> mdsim;
    dl::symbol<std::string ()> version;
};

} // namespace ljgpu
