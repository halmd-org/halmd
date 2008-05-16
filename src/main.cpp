/* Molecular Dynamics Simulation of a Lennard-Jones fluid
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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <cuda_wrapper.hpp>
#include <iostream>
#include "autocorrelation.hpp"
#include "exception.hpp"
#include "ljfluid.hpp"
#include "mdsim.hpp"
#include "options.hpp"
#include "rand48.hpp"
#include "trajectory.hpp"
#include "vector2d.hpp"
#include "vector3d.hpp"
#include "version.h"
#include "time.hpp"
using namespace std;


int main(int argc, char **argv)
{
    mdsim::options opts;

    try {
	opts.parse(argc, argv);
    }
    catch (mdsim::options::exit_exception const& e) {
	return e.status();
    }

    cuda::device::set(opts.device());

    const cuda::config dim(dim3(opts.blocks()), dim3(opts.threads()));
    mdsim::rand48 rng(dim);
#ifdef DIM_3D
    mdsim::ljfluid<3, vector3d<float> > fluid(opts.npart(), dim);
    mdsim::mdsim<3, vector3d<float> > sim;
    mdsim::trajectory<3, cuda::host::vector<vector3d<float> > > traj(opts);
    mdsim::autocorrelation<3, vector3d<float> > tcf(opts);
#else
    mdsim::ljfluid<2, vector2d<float> > fluid(opts.npart(), dim);
    mdsim::mdsim<2, vector2d<float> > sim;
    mdsim::trajectory<2, cuda::host::vector<vector2d<float> > > traj(opts);
    mdsim::autocorrelation<2, vector2d<float> > tcf(opts);
#endif

    if (opts.steps() < tcf.min_samples()) {
	throw mdsim::exception("less simulation steps than minimum required number of samples");
    }

    rng.set(opts.rngseed());

    try {
	fluid.density(opts.density());
	fluid.timestep(opts.timestep());
	fluid.temperature(opts.temp(), rng);
    }
    catch (string const& e) {
	cerr << PROGRAM_NAME ": " << e << endl;
	return EXIT_FAILURE;
    }

    cout << "### particles(" << fluid.particles() << ")" << endl;
    cout << "### density(" << fluid.density() << ")" << endl;
    cout << "### box(" << fluid.box() << ")" << endl;
    cout << "### timestep(" << fluid.timestep() << ")" << endl;
    cout << endl;

    mdsim::timer timer;
    boost::posix_time::ptime time_start = boost::posix_time::second_clock::local_time();
    mdsim::accumulator<float> time_estimated;

    timer.start();

    for (uint64_t i = 1; i <= opts.steps(); i++) {
	try {
	    sim.step(fluid);
	}
	catch (cuda::error const& e) {
	    fprintf(stderr, PROGRAM_NAME ": CUDA ERROR: %s\n", e.what());
	    return EXIT_FAILURE;
	}

	fluid.sample(tcf);
	fluid.sample(traj);

	if (i % opts.avgsteps())
	    continue;

	// compute elapsed time since start of simulation
	float time_elapsed = (boost::posix_time::second_clock::local_time() - time_start).total_seconds();
	// accumulate estimated finish time
	time_estimated += time_elapsed * opts.steps() / i;
	// print mean average of estimated finish time
	cerr << "\r" << "Estimated finish time:\t" << time_start + boost::posix_time::seconds(time_estimated.mean());

	cout << "## steps(" << i << ")" << endl;

	cout << "# en_pot(" << sim.en_pot().mean() << ")" << endl;
	cout << "# sigma_en_pot(" << sim.en_pot().std() << ")" << endl;
	cout << "# en_kin(" << sim.en_kin().mean() << ")" << endl;
	cout << "# sigma_en_kin(" << sim.en_kin().std() << ")" << endl;
	cout << "# en_tot(" << sim.en_tot().mean() << ")" << endl;
	cout << "# sigma_en_tot(" << sim.en_tot().std() << ")" << endl;
	cout << "# temp(" << sim.temp().mean() << ")" << endl;
	cout << "# sigma_temp(" << sim.temp().std() << ")" << endl;
	cout << "# pressure(" << sim.pressure().mean() << ")" << endl;
	cout << "# sigma_pressure(" << sim.pressure().std() << ")" << endl;
	cout << "# vel_cm(" << sim.vel_cm().mean() << ")" << endl;
	cout << "# sigma_vel_cm(" << sim.vel_cm().std() << ")" << endl;

	sim.clear();
    }

    tcf.flush();
    tcf.write(opts.correlations_output_file(), opts.timestep());

    timer.stop();

    cerr << "\n\n" << "GPU time: " << fluid.gputime() << "s" << endl;
    cerr << "Device memory transfer time: " << fluid.memtime() << "s" << endl;
    cerr << "Elapsed time: " << timer.elapsed() << "s" << endl;

    return EXIT_SUCCESS;
}
