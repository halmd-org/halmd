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

#include <ljgpu/mdsim/hardsphere.hpp>
#ifdef WITH_CUDA
# include <ljgpu/mdsim/ljfluid_gpu_cell.hpp>
# include <ljgpu/mdsim/ljfluid_gpu_nbr.hpp>
# include <ljgpu/mdsim/ljfluid_gpu_square.hpp>
#endif
#include <ljgpu/mdsim/ljfluid_host.hpp>
#include <ljgpu/mdsim/mdsim.hpp>
#include <ljgpu/options.hpp>
#include <ljgpu/version.h>
#include <string>
using namespace ljgpu;

extern "C" void mdlib_mdsim(options const& opt)
{
    using namespace std;
    using namespace boost;

    int const dimension = opt["dimension"].as<int>();
    if (dimension == 3) {
	mdsim<MDSIM_BACKEND<MDSIM_IMPL<3> > > md(opt);
	md();
    }
    else if (dimension == 2) {
	mdsim<MDSIM_BACKEND<MDSIM_IMPL<2> > > md(opt);
	md();
    }
    else {
	throw logic_error("invalid dimension: " + lexical_cast<string>(dimension));
    }
}

extern "C" options::description mdlib_options()
{
    return options_description<MDSIM_IMPL>();
}

extern "C" std::string mdlib_version()
{
    return PROGRAM_VERSION;
}
