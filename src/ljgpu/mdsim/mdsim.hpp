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

#ifndef LJGPU_MDSIM_HPP
#define LJGPU_MDSIM_HPP

#include <boost/bind.hpp>
#include <boost/multi_array.hpp>
#ifdef WITH_CUDA
# include <cuda_wrapper.hpp>
#endif
#include <fstream>
#include <iostream>
#include <ljgpu/mdsim/impl.hpp>
#include <ljgpu/mdsim/sample.hpp>
#include <ljgpu/mdsim/traits.hpp>
#include <ljgpu/options.hpp>
#include <ljgpu/sample/H5param.hpp>
#include <ljgpu/sample/correlation.hpp>
#include <ljgpu/sample/energy.hpp>
#include <ljgpu/sample/perf.hpp>
#include <ljgpu/sample/trajectory.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <ljgpu/util/log.hpp>
#include <ljgpu/util/signal.hpp>
#include <ljgpu/util/timer.hpp>
#include <ljgpu/version.h>
#include <sched.h>
#include <stdint.h>
#include <vector>
#include <unistd.h>

namespace ljgpu
{

/**
 * Molecular Dynamics simulation of a Lennard-Jones fluid
 */
template <typename mdsim_backend>
class mdsim
{
public:
    typedef typename mdsim_backend::impl_type impl_type;
    typedef typename mdsim_backend::traits_type traits_type;
    typedef typename traits_type::float_type float_type;
    typedef typename traits_type::vector_type vector_type;
    typedef typename traits_type::sample_type sample_type;
    enum { dimension = traits_type::dimension };

public:
    enum {
	/** HDF5 buffers flush to disk interval in seconds */
	FLUSH_TO_DISK_INTERVAL = 900,
	/** waiting time in seconds before runtime estimate after block completion */
	TIME_ESTIMATE_WAIT_AFTER_BLOCK = 300,
	/** runtime estimate interval in seconds */
	TIME_ESTIMATE_INTERVAL = 1800,
    };

public:
    /** initialize MD simulation program */
    mdsim(options const& opt);
    /** run MD simulation program */
    void operator()();
    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

private:
    /** backend-specific API wrappers */
    void epsilon(boost::true_type const&);
    void sigma(boost::true_type const&);
    void cutoff_radius(boost::true_type const&);
    void potential_smoothing(boost::true_type const&);
    void pair_separation(boost::true_type const&);
    void cell_occupancy(boost::true_type const&);
    void nbl_skin(boost::true_type const&);
    void init_cells(boost::true_type const&);
    void init_event_list(boost::true_type const&);
    void threads(boost::true_type const&);
    void thermostat(boost::true_type const&);
    void cuda_device(boost::true_type const&);
    void cuda_allocated_memory(boost::true_type const&);
    void stream(boost::true_type const&);
    void copy(boost::true_type const&);

    void epsilon(boost::false_type const&) {}
    void sigma(boost::false_type const&) {}
    void cutoff_radius(boost::false_type const&) {}
    void potential_smoothing(boost::false_type const&) {}
    void pair_separation(boost::false_type const&) {}
    void cell_occupancy(boost::false_type const&) {}
    void nbl_skin(boost::false_type const&) {}
    void init_cells(boost::false_type const&) {}
    void init_event_list(boost::false_type const&) {}
    void threads(boost::false_type const&) {}
    void thermostat(boost::false_type const&) {}
    void cuda_device(boost::false_type const&) {}
    void cuda_allocated_memory(boost::false_type const&) {}
    void stream(boost::false_type const&) {}
    void copy(boost::false_type const&) {}

    /** bind process to CPU core(s) */
    void cpu_set(std::vector<int> const& cpu_set);

private:
    /** program options */
    options const& opt;
    /** Lennard-Jones fluid simulation */
    mdsim_backend fluid;
    /** block correlations */
    correlation<dimension> tcf;
    /**  trajectory file writer */
    trajectory traj;
    /** thermodynamic equilibrium properties */
    energy<dimension> tep;
    /** performance data */
    perf prf;
};

/**
 * initialize MD simulation program
 */
template <typename mdsim_backend>
mdsim<mdsim_backend>::mdsim(options const& opt) : opt(opt)
{
    // set CPU core(s)
    if (!opt["processor"].empty()) {
	cpu_set(opt["processor"].as<std::vector<int> >());
    }
    // set CUDA device
    cuda_device(boost::is_base_of<ljfluid_impl_gpu_base<dimension>, impl_type>());

    LOG("positional coordinates dimension: " << dimension);

    // set Lennard-Jones or hard-sphere system parameters
    if (!opt["binary"].empty() && opt["particles"].defaulted()) {
	fluid.particles(opt["binary"].as<boost::array<unsigned int, 2> >());
    }
    else {
	fluid.particles(opt["particles"].as<unsigned int>());
    }
    if (opt["density"].defaulted() && !opt["box-length"].empty()) {
	fluid.box(opt["box-length"].as<float>());
    }
    else {
	fluid.density(opt["density"].as<float>());
    }
    fluid.timestep(opt["timestep"].as<float>());

    epsilon(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());
    sigma(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());
    cutoff_radius(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());
    potential_smoothing(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());
    pair_separation(boost::is_base_of<hardsphere_impl<dimension>, impl_type>());

    cell_occupancy(boost::is_base_of<ljfluid_impl_gpu_cell<dimension>, impl_type>());
    cell_occupancy(boost::is_base_of<ljfluid_impl_gpu_neighbour<dimension>, impl_type>());
    nbl_skin(boost::is_base_of<ljfluid_impl_gpu_neighbour<dimension>, impl_type>());
    threads(boost::is_base_of<ljfluid_impl_gpu_base<dimension>, impl_type>());

    nbl_skin(boost::is_base_of<ljfluid_impl_host<dimension>, impl_type>());
    init_cells(boost::is_base_of<hardsphere_impl<dimension>, impl_type>());

    thermostat(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());

    // initialize random number generator with seed
    if (opt["random-seed"].empty()) {
	LOG("obtaining 32-bit integer seed from /dev/random");
	unsigned int seed;
	try {
	    std::ifstream rand;
	    rand.exceptions(std::ifstream::eofbit|std::ifstream::failbit|std::ifstream::badbit);
	    rand.open("/dev/random");
	    rand.read(reinterpret_cast<char*>(&seed), sizeof(seed));
	    rand.close();
	}
	catch (std::ifstream::failure const& e) {
	    throw std::logic_error(std::string("failed to read from /dev/random: ") + e.what());
	}
	fluid.rng(seed);
    }
    else {
	fluid.rng(opt["random-seed"].as<unsigned int>());
    }

    if (!opt["trajectory-sample"].empty()) {
	int64_t const index = opt["trajectory-sample"].as<int64_t>();
	// open trajectory input file
	traj.open(opt["trajectory"].as<std::string>(), trajectory::in);
	// read trajectory sample and restore system state
	fluid.sample(boost::bind(&trajectory::read<sample_type>, boost::ref(traj), _1, index));
	// close trajectory input file
	traj.close();
    }
    else {
	// arrange particles on a face-centered cubic (fcc) lattice
	fluid.lattice();
    }

    if (opt["trajectory-sample"].empty() || opt["discard-velocities"].as<bool>()) {
	// set system temperature according to Maxwell-Boltzmann distribution
	fluid.temperature(opt["temperature"].as<float>());
    }
    // initialize event list
    init_event_list(boost::is_base_of<hardsphere_impl<dimension>, impl_type>());

    // print GPU memory usage
    cuda_allocated_memory(boost::is_base_of<ljfluid_impl_gpu_base<dimension>, impl_type>());

    if (!opt["time"].empty()) {
	// set total simulation time
	tcf.time(opt["time"].as<float>(), fluid.timestep());
    }
    else {
	// set total number of simulation steps
	tcf.steps(opt["steps"].as<uint64_t>(), fluid.timestep());
    }
    // set sample rate for lowest block level
    tcf.sample_rate(opt["sample-rate"].as<unsigned int>());
    // set block size
    tcf.block_size(opt["block-size"].as<unsigned int>());
    // set maximum number of samples per block
    tcf.max_samples(opt["max-samples"].as<uint64_t>());

    std::vector<float> q;
    if (!opt["q-values"].empty()) {
	boost::multi_array<float, 1> v = opt["q-values"].as<boost::multi_array<float, 1> >();
	q.assign(v.begin(), v.end());
    }
    else {
	// static structure factor peak at q ~ 2pi/sigma
	q.push_back(2 * M_PI);
    }
    tcf.q_values(q, opt["q-error"].as<float>(), fluid.box());
}

/**
 * run MD simulation program
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::operator()()
{
    if (opt["dry-run"].as<bool>()) {
	// test parameters only
	return;
    }
    if (opt["daemon"].as<bool>()) {
	// run program in background
	daemon(0, 0);
    }

    // measure elapsed realtime
    real_timer timer;

    // performance data
    prf.open(opt["output"].as<std::string>() + ".prf");
    H5param(prf) << *this << fluid << tcf;

    // time correlation functions
    if (!opt["disable-correlation"].as<bool>()) {
	tcf.open(opt["output"].as<std::string>() + ".tcf", fluid.mixture() == BINARY);
	H5param(tcf) << *this << fluid << tcf;
    }
    // trajectory file writer
    traj.open(opt["output"].as<std::string>() + ".trj", trajectory::out);
    H5param(traj) << *this << fluid << tcf;
    // thermodynamic equilibrium properties
    if (!opt["disable-energy"].as<bool>()) {
	tep.open(opt["output"].as<std::string>() + ".tep");
	H5param(tep) << *this << fluid << tcf;
    }

    // schedule first disk flush
    alarm(FLUSH_TO_DISK_INTERVAL);

    LOG("starting MD simulation");
    timer.start();

    for (count_timer<uint64_t> step(0); step < tcf.steps(); ++step) {
	// check if sample is acquired for given simulation step
	if (tcf.sample(step)) {
	    // copy previous MD simulation state from GPU to host
	    copy(boost::is_base_of<ljfluid_impl_gpu_base<dimension>, impl_type>());
	    // sample phase space
	    copy(boost::is_base_of<hardsphere_impl<dimension>, impl_type>());
	}
	else {
	    // cell lists implementation requires GPU to host copy every step
	    copy(boost::is_base_of<ljfluid_impl_gpu_cell<dimension>, impl_type>());
	}

	// stream next MD simulation program step on GPU
	stream(boost::is_base_of<ljfluid_impl_gpu_base<dimension>, impl_type>());

	// check if sample is acquired for given simulation step
	if (tcf.sample(step)) {
	    bool flush = false;
	    // simulation time
	    double time = step * (double)fluid.timestep();
	    // sample time correlation functions
	    if (!opt["disable-correlation"].as<bool>()) {
		tcf.sample(fluid.sample(), step, flush);
	    }
	    // sample thermodynamic equilibrium properties
	    if (!opt["disable-energy"].as<bool>()) {
		tep.sample(fluid.sample(), fluid.density(), time);
	    }
	    // sample trajectory
	    if (opt["enable-trajectory"].as<bool>() || step == 0) {
		traj.write(fluid.sample(), time);
		if (step == 0) {
		    traj.flush();
		}
	    }
	    // acquired maximum number of samples for a block level
	    if (flush) {
		// sample performance counters
		prf.sample(fluid.times());
		// write partial results to HDF5 files and flush to disk
		if (!opt["disable-correlation"].as<bool>()) {
		    tcf.flush();
		}
		if (opt["enable-trajectory"].as<bool>()) {
		    traj.flush();
		}
		if (!opt["disable-energy"].as<bool>()) {
		    tep.flush();
		}
		prf.flush();
		LOG("flushed HDF5 buffers to disk");
		// schedule remaining runtime estimate
		step.clear();
		step.set(TIME_ESTIMATE_WAIT_AFTER_BLOCK);
		// schedule next disk flush
		alarm(FLUSH_TO_DISK_INTERVAL);
	    }
	}
	// synchronize MD simulation program step on GPU
	fluid.mdstep();

	// check whether a runtime estimate has finished
	if (step.estimate() > 0) {
	    LOG("estimated remaining runtime: " << real_timer::format(step.estimate()));
	    step.clear();
	    // schedule next remaining runtime estimate
	    step.set(TIME_ESTIMATE_INTERVAL);
	}

	// process signal event
	if (signal::poll()) {
	    if (signal::signal != SIGALRM) {
		LOG_WARNING("trapped signal " << signal::signal << " at simulation step " << step);
	    }
	    if (signal::signal == SIGUSR1) {
		// schedule runtime estimate now
		step.set(0);
	    }
	    else if (signal::signal == SIGHUP || signal::signal == SIGALRM) {
		// sample performance counters
		prf.sample(fluid.times());
		// write partial results to HDF5 files and flush to disk
		if (!opt["disable-correlation"].as<bool>()) {
		    tcf.flush();
		}
		if (opt["enable-trajectory"].as<bool>()) {
		    traj.flush();
		}
		if (!opt["disable-energy"].as<bool>()) {
		    tep.flush();
		}
		prf.flush();
		LOG("flushed HDF5 buffers to disk");
		// schedule next disk flush
		alarm(FLUSH_TO_DISK_INTERVAL);
	    }
	    else if (signal::signal == SIGINT || signal::signal == SIGTERM) {
		LOG_WARNING("aborting simulation");
		break;
	    }
	}
    }

    // copy previous MD simulation state from GPU to host
    copy(boost::is_base_of<ljfluid_impl_gpu_base<dimension>, impl_type>());
    // sample phase space
    copy(boost::is_base_of<hardsphere_impl<dimension>, impl_type>());
    // save last phase space sample
    traj.write(fluid.sample(), tcf.steps() * fluid.timestep());
    // sample performance counters
    prf.sample(fluid.times());

    // print performance statistics
    std::stringstream is;
    std::string str;
    is << prf.times();
    while (std::getline(is, str)) {
	LOG(str);
    }

    timer.stop();
    LOG("finished MD simulation");
    LOG("total MD simulation runtime: " << timer);

    // cancel previously scheduled disk flush
    alarm(0);
    // close HDF5 output files
    if (!opt["disable-correlation"].as<bool>()) {
	tcf.close();
    }
    traj.close();
    if (!opt["disable-energy"].as<bool>()) {
	tep.close();
    }
    prf.close();
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::epsilon(boost::true_type const&)
{
    if (!opt["binary"].empty() && opt["particles"].defaulted()) {
	fluid.epsilon(opt["epsilon"].as<boost::array<float, 3> >());
    }
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::sigma(boost::true_type const&)
{
    if (!opt["binary"].empty() && opt["particles"].defaulted()) {
	fluid.sigma(opt["sigma"].as<boost::array<float, 3> >());
    }
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::cutoff_radius(boost::true_type const&)
{
    fluid.cutoff_radius(opt["cutoff"].as<float>());
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::potential_smoothing(boost::true_type const&)
{
    if (!opt["smooth"].empty()) {
	fluid.potential_smoothing(opt["smooth"].as<float>());
    }
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::pair_separation(boost::true_type const&)
{
    fluid.pair_separation(opt["pair-separation"].as<float>());
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::cell_occupancy(boost::true_type const&)
{
    fluid.cell_occupancy(opt["cell-occupancy"].as<float>());
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::nbl_skin(boost::true_type const&)
{
    fluid.nbl_skin(opt["skin"].as<float>());
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::threads(boost::true_type const&)
{
    fluid.threads(opt["threads"].as<unsigned int>());
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::init_cells(boost::true_type const&)
{
    fluid.init_cells();
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::init_event_list(boost::true_type const&)
{
    fluid.init_event_list();
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::thermostat(boost::true_type const&)
{
    if (!opt["thermostat"].empty()) {
	float const nu = opt["thermostat"].as<float>();
	float const temp = opt["temperature"].as<float>();
	fluid.thermostat(nu, temp);
    }
}

/**
 * set CUDA device for host context
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::cuda_device(boost::true_type const&)
{
#ifdef WITH_CUDA
    int dev = opt["device"].as<int>();
    if (opt["device"].defaulted()) {
	char const* env = getenv("CUDA_DEVICE");
	if (env != NULL && *env != '\0') {
	    char* endptr;
	    int i = strtol(env, &endptr, 10);
	    if (*endptr != '\0' || i < 0) {
		throw std::logic_error(std::string("CUDA_DEVICE environment variable invalid: ") + env);
	    }
	    dev = i;
	}
    }
    cuda::device::set(dev);
    LOG("CUDA device: " << cuda::device::get());

    // query CUDA device properties
    cuda::device::properties prop(cuda::device::get());
    LOG("CUDA device name: " << prop.name());
    LOG("CUDA device total global memory: " << prop.total_global_mem() << " bytes");
    LOG("CUDA device shared memory per block: " << prop.shared_mem_per_block() << " bytes");
    LOG("CUDA device registers per block: " << prop.regs_per_block());
    LOG("CUDA device warp size: " << prop.warp_size());
    LOG("CUDA device maximum number of threads per block: " << prop.max_threads_per_block());
    LOG("CUDA device total constant memory: " << prop.total_const_mem());
    LOG("CUDA device major revision: " << prop.major());
    LOG("CUDA device minor revision: " << prop.minor());
    LOG("CUDA device clock frequency: " << prop.clock_rate() << " kHz");
#endif /* WITH_CUDA */
}

/**
 * print GPU memory usage
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::cuda_allocated_memory(boost::true_type const&)
{
#ifdef WITH_CUDA
    int const dev = cuda::device::get();
    LOG("GPU allocated global device memory: " << cuda::device::mem_get_used(dev) << " bytes");
    LOG("GPU available global device memory: " << cuda::device::mem_get_free(dev) << " bytes");
    LOG("GPU total global device memory: " << cuda::device::mem_get_total(dev) << " bytes");
#endif /* WITH_CUDA */
}

/**
 * bind process to CPU core(s)
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::cpu_set(std::vector<int> const& cpu_set)
{
    cpu_set_t mask;
    CPU_ZERO(&mask);
    BOOST_FOREACH(int cpu, cpu_set) {
	LOG("adding CPU core " << cpu << " to process CPU affinity mask");
	CPU_SET(cpu, &mask);
    }
    if (0 != sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask)) {
	throw std::logic_error("failed to set process CPU affinity mask");
    }
}

/**
 * write parameters to HDF5 parameter group
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::param(H5param& param) const
{
    H5xx::group node(param["mdsim"]);
    node["backend"] = opt["backend"].as<std::string>();
    node["dimension"] = (unsigned int) dimension;

    node = param["program"];
    node["name"] = PROGRAM_NAME;
    node["version"] = PROGRAM_VERSION;
    node["variant"] = PROGRAM_VARIANT;
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::stream(boost::true_type const&)
{
    fluid.stream();
}

template <typename mdsim_backend>
void mdsim<mdsim_backend>::copy(boost::true_type const&)
{
    fluid.copy();
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_HPP */
