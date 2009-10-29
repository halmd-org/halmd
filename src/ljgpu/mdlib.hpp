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

#ifndef LJGPU_MDLIB_HPP
#define LJGPU_MDLIB_HPP

#include <ljgpu/options.hpp>
#include <ljgpu/util/dl.hpp>

extern "C" int mdlib_mdsim(ljgpu::options const& opt);
extern "C" boost::program_options::options_description mdlib_options();
extern "C" std::string mdlib_backend();
extern "C" std::string mdlib_variant();
extern "C" std::string mdlib_version();

namespace ljgpu
{

#ifndef BACKEND_EXECUTABLES

struct mdlib : public dl::library
{
    void open(boost::filesystem::path const& path)
    {
        dl::library::open(path);
        options.set(*this, "mdlib_options");
        mdsim.set(*this, "mdlib_mdsim");
        backend.set(*this, "mdlib_backend");
        variant.set(*this, "mdlib_variant");
        version.set(*this, "mdlib_version");
    }

    dl::symbol<boost::program_options::options_description ()> options;
    dl::symbol<int (ljgpu::options const&)> mdsim;
    dl::symbol<std::string ()> backend;
    dl::symbol<std::string ()> variant;
    dl::symbol<std::string ()> version;
};

#else /* ! BACKEND_EXECUTABLES */

struct mdlib
{
    mdlib() : options(mdlib_options), mdsim(mdlib_mdsim), backend(mdlib_backend), variant(mdlib_variant) {}

    boost::function<boost::program_options::options_description ()> options;
    boost::function<int (ljgpu::options const&)> mdsim;
    boost::function<std::string ()> backend;
    boost::function<std::string ()> variant;
};

#endif /* ! BACKEND_EXECUTABLES */

} // namespace ljgpu

#endif /* ! LJGPU_MDLIB_HPP */
