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

#include <boost/assign.hpp>
#include <boost/multi_array.hpp>
#include <boost/noncopyable.hpp>
#include <boost/unordered_map.hpp>
#include <fstream>
#include <iostream>
#include <ljgpu/mdsim/impl.hpp>
#include <ljgpu/mdsim/sample.hpp>
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
#include <stdint.h>
#include <vector>
#include <unistd.h>

namespace ljgpu
{

/**
 * Molecular Dynamics simulation of a Lennard-Jones fluid
 */
template <typename mdsim_backend>
class mdsim : boost::noncopyable
{
public:
    typedef typename mdsim_backend::impl_type impl_type;
    typedef typename mdsim_backend::float_type float_type;
    typedef typename mdsim_backend::vector_type vector_type;
    typedef typename mdsim_backend::host_sample_type host_sample_type;
    typedef typename mdsim_backend::trajectory_sample_type trajectory_sample_type;
    typedef typename mdsim_backend::trajectory_sample_variant trajectory_sample_variant;
    typedef typename mdsim_backend::energy_sample_type energy_sample_type;
    enum { dimension = mdsim_backend::dimension };

public:
    /** initialise MD simulation */
    mdsim(options const& opt);
    /** run MD simulation */
    void operator()();
    /** write parameters to HDF5 parameter group */
    void param(H5param param) const;

private:
    /** aquire fluid sample */
    void sample_fluid(uint64_t step, bool dump);
    /** process fluid sample */
    bool sample_properties(uint64_t step, bool perf);
    /** open HDF5 output files */
    void open();
    /** write partial results to HDF5 files and flush to disk */
    void flush();
    /** close HDF5 output files */
    void close();
    /** obtain integer seed from file */
    static unsigned int read_random_seed(std::string const& fn);

    /** backend-specific API wrappers */
    void epsilon(boost::false_type const&) {}
    void epsilon(boost::true_type const&)
    {
	if (!m_opt["binary"].empty() && m_opt["particles"].defaulted()) {
	    m_fluid.epsilon(m_opt["epsilon"].as<boost::array<float, 3> >());
	}
    }

    void sigma(boost::false_type const&) {}
    void sigma(boost::true_type const&)
    {
	if (!m_opt["binary"].empty() && m_opt["particles"].defaulted()) {
	    m_fluid.sigma(m_opt["sigma"].as<boost::array<float, 3> >());
	}
    }

    void cutoff_radius(boost::false_type const&) {}
    void cutoff_radius(boost::true_type const&)
    {
	m_fluid.cutoff_radius(m_opt["cutoff"].as<float>());
    }

    void potential_smoothing(boost::false_type const&) {}
    void potential_smoothing(boost::true_type const&)
    {
	if (!m_opt["smooth"].empty()) {
	    m_fluid.potential_smoothing(m_opt["smooth"].as<float>());
	}
    }

    void pair_separation(boost::false_type const&) {}
    void pair_separation(boost::true_type const&)
    {
	m_fluid.pair_separation(m_opt["pair-separation"].as<float>());
    }

    void cell_occupancy(boost::false_type const&) {}
    void cell_occupancy(boost::true_type const&)
    {
	m_fluid.cell_occupancy(m_opt["cell-occupancy"].as<float>());
    }

    void nbl_skin(boost::false_type const&) {}
    void nbl_skin(boost::true_type const&)
    {
	m_fluid.nbl_skin(m_opt["skin"].as<float>());
    }

    void threads(boost::false_type const&) {}
    void threads(boost::true_type const&)
    {
	m_fluid.threads(m_opt["threads"].as<unsigned int>());
    }

    void init_cells(boost::false_type const&) {}
    void init_cells(boost::true_type const&)
    {
	m_fluid.init_cells();
    }

    void init_event_list(boost::false_type const&) {}
    void init_event_list(boost::true_type const&)
    {
	m_fluid.init_event_list();
    }

    void rescale_energy(boost::false_type const&) {}
    void rescale_energy(boost::true_type const&)
    {
	if (!m_opt["energy"].empty()) {
	    m_fluid.rescale_energy(m_opt["energy"].as<float>());
	}
    }

    void thermostat(boost::false_type const&) {}
    void thermostat(boost::true_type const&)
    {
	if (!m_opt["thermostat"].empty()) {
	    float const nu = m_opt["thermostat"].as<float>();
	    float const temp = m_opt["temperature"].as<float>();
	    m_fluid.thermostat(nu, temp);
	}
    }

    void stream(boost::false_type const&) {}
    void stream(boost::true_type const&)
    {
	m_fluid.stream();
    }

private:
    /** program options */
    options const& m_opt;
    /** Lennard-Jones fluid simulation */
    mdsim_backend m_fluid;

    /** block correlations */
    correlation<dimension> m_corr;
    /**  trajectory file writer */
    trajectory m_traj;
    /** thermodynamic equilibrium properties */
    energy<dimension> m_en;
    /** performance data */
    perf m_perf;

    /** trajectory sample */
    host_sample_type m_traj_sample;
    /** trajectory sample for correlation functions */
    trajectory_sample_variant m_corr_sample;
    /** thermodynamic equilibrium properties sample */
    energy_sample_type m_en_sample;

    bool m_is_traj_step;
    bool m_is_corr_step;
    bool m_is_en_step;
    bool m_is_corr_sample_gpu;
};

/**
 * initialise MD simulation
 */
template <typename mdsim_backend>
mdsim<mdsim_backend>::mdsim(options const& opt) : m_opt(opt)
{
    LOG("positional coordinates dimension: " << dimension);

    // number of particles in periodic simulation box
    if (!m_opt["binary"].empty() && m_opt["particles"].defaulted()) {
	m_fluid.particles(m_opt["binary"].as<boost::array<unsigned int, 2> >());
    }
    else {
	m_fluid.particles(m_opt["particles"].as<unsigned int>());
    }
    if (m_opt["density"].defaulted() && !m_opt["box-length"].empty()) {
	// periodic simulation box length
	m_fluid.box(m_opt["box-length"].as<float>());
    }
    else {
	// number density
	m_fluid.density(m_opt["density"].as<float>());
    }
    // simulation timestep
    m_fluid.timestep(m_opt["timestep"].as<double>());

    // potential well depths
    epsilon(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());
    // collision diameters
    sigma(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());
    // potential cutoff radius
    cutoff_radius(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());
    // potential smoothing function scale parameter
    potential_smoothing(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());
    // pair separation at which particle collision occurs
    pair_separation(boost::is_base_of<hardsphere_impl<dimension>, impl_type>());

    // desired average cell occupancy
    cell_occupancy(boost::is_base_of<ljfluid_impl_gpu_cell<dimension>, impl_type>());
    cell_occupancy(boost::is_base_of<ljfluid_impl_gpu_neighbour<dimension>, impl_type>());
    // neighbour list skin
    nbl_skin(boost::is_base_of<ljfluid_impl_gpu_neighbour<dimension>, impl_type>());
    nbl_skin(boost::is_base_of<ljfluid_impl_host<dimension>, impl_type>());
    // number of CUDA execution threads
    threads(boost::is_base_of<ljfluid_impl_gpu_base<dimension>, impl_type>());
    // initialise hard-sphere cell lists
    init_cells(boost::is_base_of<hardsphere_impl<dimension>, impl_type>());

    // heat bath collision probability and temperature
    thermostat(typename mdsim_backend::has_thermostat());

    if (m_opt["random-seed"].empty()) {
	m_fluid.rng(read_random_seed("/dev/random"));
    }
    else {
	m_fluid.rng(m_opt["random-seed"].as<unsigned int>());
    }

    if (!m_opt["trajectory-sample"].empty()) {
	// read trajectory sample and periodic simulation box length
	host_sample_type sample;
	m_traj.open(m_opt["trajectory"].as<std::string>(), trajectory::in);
	m_traj.read(sample, m_opt["trajectory-sample"].as<int64_t>());
	float box = H5param(m_traj)["mdsim"]["box_length"].as<float>();
	m_traj.close();
	// restore system state
	m_fluid.state(sample, box);
    }
    else {
	// arrange particles on a face-centered cubic (fcc) lattice
	m_fluid.lattice();
    }

    // rescale mean particle energy
    rescale_energy(boost::is_base_of<ljfluid_impl_base<dimension>, impl_type>());

    if (m_opt["trajectory-sample"].empty() || !m_opt["temperature"].defaulted()) {
	// initialise velocities from Maxwell-Boltzmann distribution
	m_fluid.temperature(m_opt["temperature"].as<float>());
    }
    // initialise hard-sphere event list
    init_event_list(boost::is_base_of<hardsphere_impl<dimension>, impl_type>());

    if (m_opt["steps"].defaulted() && !m_opt["time"].empty()) {
	// total simulation time
	m_corr.time(m_opt["time"].as<double>(), m_fluid.timestep());
    }
    else {
	// total number of simulation steps
	m_corr.steps(m_opt["steps"].as<uint64_t>(), m_fluid.timestep());
    }
    // sample rate for lowest block level
    m_corr.sample_rate(m_opt["sample-rate"].as<unsigned int>());
    // minimum number of trajectory samples
    m_corr.min_samples(m_opt["min-samples"].as<uint64_t>());
    // maximum number of samples per block
    m_corr.max_samples(m_opt["max-samples"].as<uint64_t>());
    // block size
    m_corr.block_size(m_opt["block-size"].as<unsigned int>());

    if (!m_opt["disable-correlation"].as<bool>()) {
	std::vector<float> q;
	if (m_opt["q-values"].empty()) {
	    // static structure factor peak at q ~ 2pi/sigma
	    q.push_back(2 * M_PI);
	}
	else {
	    typedef boost::multi_array<float, 1> q_value_vector;
	    q_value_vector v = m_opt["q-values"].as<q_value_vector>();
	    q.assign(v.begin(), v.end());
	}
	m_corr.q_values(q, m_opt["q-error"].as<float>(), m_fluid.box());

	std::string const backend(m_opt["tcf-backend"].as<std::string>());
	if (backend == "host") {
	    m_corr.add_host_correlation_functions(m_fluid.is_binary() ? 2 : 1);
	    m_is_corr_sample_gpu = false;
	}
#if WITH_CUDA
	else if (backend == "gpu") {
	    m_corr.add_gpu_correlation_functions(m_fluid.is_binary() ? 2 : 1);
	    m_is_corr_sample_gpu = mdsim_backend::has_trajectory_gpu_sample::value;
	}
#endif
	else {
	    throw std::logic_error("unknown correlation function backend: " + backend);
	}
    }
}

/**
 * run MD simulation
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::operator()()
{
    /** HDF5 buffers flush to disk interval in seconds */
    enum { FLUSH_TO_DISK_INTERVAL = 900 };
    /** runtime estimate interval in seconds */
    enum { TIME_ESTIMATE_INTERVAL = 1800 };
    /** waiting time in seconds before runtime estimate after block completion */
    enum { TIME_ESTIMATE_WAIT_AFTER_BLOCK = 300 };

    boost::unordered_map<int, std::string> const signals = boost::assign::map_list_of
	(SIGINT, "INT")
	(SIGTERM, "TERM")
	(SIGHUP, "HUP")
	(SIGUSR1, "USR1")
	(SIGUSR2, "USR2");

    // open HDF5 output files
    open();
    // schedule first disk flush
    alarm(FLUSH_TO_DISK_INTERVAL);

    LOG("starting MD simulation");
    real_timer timer;
    timer.start();

    count_timer<uint64_t> step;

    for (step = 0; step < m_corr.steps(); ++step) {
	sample_fluid(step, !step);
	// stream next MD simulation program step on GPU
	stream(boost::is_base_of<ljfluid_impl_gpu_base<dimension>, impl_type>());

	if (sample_properties(step, false)) {
	    // acquired maximum number of samples for a block level
	    flush();
	    step.set(TIME_ESTIMATE_WAIT_AFTER_BLOCK);
	    alarm(FLUSH_TO_DISK_INTERVAL);
	}

	try {
	    // synchronize MD simulation program step on GPU
	    m_fluid.mdstep();
	}
	catch (potential_energy_divergence const& e) {
	    step++;
	    LOG_ERROR(e.what() << " at step " << step);
	    break;
	}

	if (step.estimate() > 0) {
	    LOG("estimated remaining runtime: " << real_timer::format(step.estimate()));
	    step.set(TIME_ESTIMATE_INTERVAL);
	}

	if (signal::poll()) {
	    if (signals.find(signal::signal) != signals.end()) {
		std::string const& name = signals.at(signal::signal);
		LOG_WARNING("trapped signal " + name + " at simulation step " << step);
	    }
	    if (signal::signal == SIGUSR2) {
		// block process until further signal is received
		LOG("pausing simulation");
		unsigned int seconds = alarm(0);
		signal::wait();
		alarm(seconds);
		LOG("resuming simulation");
	    }
	    if (signal::signal == SIGINT || signal::signal == SIGTERM) {
		step++;
		LOG_WARNING("aborting simulation at step " << step);
		break;
	    }
	    else if (signal::signal == SIGHUP || signal::signal == SIGALRM) {
		flush();
		alarm(FLUSH_TO_DISK_INTERVAL);
	    }
	    else if (signal::signal == SIGUSR1) {
		step.set(0);
	    }
	}
    }
    sample_fluid(step, true);
    sample_properties(step, true);

    timer.stop();
    LOG("finished MD simulation");

    // print performance statistics
    std::stringstream is;
    std::string str;
    is << m_perf.times();
    while (std::getline(is, str)) {
	LOG(str);
    }
    LOG("total MD simulation runtime: " << timer);

    // close HDF5 output files
    close();
}

/**
 * aquire fluid sample
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::sample_fluid(uint64_t step, bool dump)
{
    m_is_corr_step = (m_corr.is_sample_step(step) && m_corr.is_open());

    // sample trajectory on GPU
    if (m_is_corr_sample_gpu && m_is_corr_step) {
	trajectory_sample_type sample;
	m_fluid.sample(sample);
	// copy pointers to GPU phase space samples in binary mixture
	m_corr_sample = sample;
    }

    m_is_traj_step = dump || (m_corr.is_trajectory_step(step) && !m_opt["disable-trajectory"].as<bool>());

    // sample trajectory on host
    if ((!m_is_corr_sample_gpu && m_is_corr_step) || m_is_traj_step) {
	host_sample_type sample;
	m_fluid.sample(sample);
	// copy pointers to host phase space samples in binary mixture
	if (!m_is_corr_sample_gpu) {
	    m_corr_sample = sample;
	}
	m_traj_sample = sample;
    }

    m_is_en_step = (m_en.is_open() && m_corr.is_sample_step(step));

    // sample thermodynamic equilibrium properties
    if (m_is_corr_step || m_is_en_step) {
	m_fluid.sample(m_en_sample);
    }
}

/**
 * process fluid sample
 */
template <typename mdsim_backend>
bool mdsim<mdsim_backend>::sample_properties(uint64_t step, bool perf)
{
    double time = step * static_cast<double>(m_fluid.timestep());
    bool flush = false;

    if (m_is_corr_step) {
	m_corr.sample(m_corr_sample, m_en_sample, step, flush);
    }
    if (m_is_traj_step) {
	m_traj.write(m_traj_sample, time);
    }
    if (m_is_en_step) {
	m_en.sample(m_en_sample, m_fluid.density(), time);
    }
    if (perf) {
	m_perf.sample(m_fluid.times());
    }
    return flush;
}

/**
 * open HDF5 output files
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::open()
{
    std::string fn = m_opt["output"].as<std::string>();

    if (!m_opt["disable-correlation"].as<bool>()) {
	m_corr.open(fn + ".tcf", m_fluid.is_binary() ? 2 : 1);
	param(m_corr);
    }

    m_traj.open(fn + ".trj", trajectory::out);
    param(m_traj);

    if (!m_opt["disable-energy"].as<bool>()) {
	m_en.open(fn + ".tep");
	param(m_en);
    }

    m_perf.open(fn + ".prf");
    param(m_perf);
}

/**
 * write partial results to HDF5 files and flush to disk
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::flush()
{
    if (m_corr.is_open()) {
	m_corr.flush();
    }
    if (m_traj.is_open()) {
	m_traj.flush();
    }
    if (m_en.is_open()) {
	m_en.flush();
    }
    if (m_perf.is_open()) {
	m_perf.sample(m_fluid.times());
	m_perf.flush();
    }
    LOG("flushed HDF5 buffers to disk");
}

/**
 * close HDF5 output files
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::close()
{
    if (m_corr.is_open()) {
	m_corr.close();
    }
    if (m_traj.is_open()) {
	m_traj.close();
    }
    if (m_en.is_open()) {
	m_en.close();
    }
    if (m_perf.is_open()) {
	m_perf.close();
    }
}

/**
 * obtain integer seed from file
 */
template <typename mdsim_backend>
unsigned int mdsim<mdsim_backend>::read_random_seed(std::string const& fn)
{
    using namespace std;
    unsigned int seed;

    LOG("obtaining 32-bit integer seed from " + fn);
    try {
	ifstream rand;
	rand.exceptions(ifstream::eofbit|ifstream::failbit|ifstream::badbit);
	rand.open(fn.c_str());
	rand.read(reinterpret_cast<char*>(&seed), sizeof(seed));
	rand.close();
    }
    catch (ifstream::failure const& e) {
	throw logic_error("failed to read from " + fn + e.what());
    }
    return seed;
}

/**
 * write parameters to HDF5 parameter group
 */
template <typename mdsim_backend>
void mdsim<mdsim_backend>::param(H5param param) const
{
    H5xx::group node(param["mdsim"]);
    node["backend"] = m_opt["backend"].as<std::string>();
    if (!m_opt["disable-correlation"].as<bool>()) {
	node["tcf_backend"] = m_opt["tcf-backend"].as<std::string>();
    }
    node["dimension"] = (unsigned int) dimension;

    node = param["program"];
    node["name"] = PROGRAM_NAME;
    node["version"] = PROGRAM_VERSION;
    node["variant"] = PROGRAM_VARIANT;

    param << m_fluid << m_corr;
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_HPP */
