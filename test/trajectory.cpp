/* MD simulation trajectory writer
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

#include <stdlib.h>
#include <cuda_wrapper.hpp>
#include "../src/trajectory.hpp"

#define ATOMS 100
#define STEPS 1000

int main(int argc, char** argv)
{
#ifdef DIM_3D
    mdsim::trajectory<3, cuda::host::vector<float3> > traj("ljfluid3d.trj", ATOMS, STEPS);
    cuda::host::vector<float3> r(0), v(0);
#else
    mdsim::trajectory<2, cuda::host::vector<float2> > traj("ljfluid2d.trj", ATOMS, STEPS);
    cuda::host::vector<float2> r(0), v(0);
#endif

    srand48(42);
    for (size_t i = 0; i < STEPS; ++i) {
	for (size_t j = 0; j < ATOMS; ++j) {
#ifdef DIM_3D
	    r.push_back(make_float3(drand48() * 10., drand48() * 0.01 + 3.14, drand48() * 0.01 + 42.));
	    v.push_back(make_float3(drand48() - 0.5, (drand48() - 0.5) * 0.01 - 3.14, (drand48() - 0.5) - 42.));
#else
	    r.push_back(make_float2(drand48() * 10., drand48() * 0.01 + 42.));
	    v.push_back(make_float2(drand48() - 0.5, (drand48() - 0.5) - 42.));
#endif
	}
	traj.write(r, v);
	r.clear();
	v.clear();
    }

    return EXIT_SUCCESS;
}
