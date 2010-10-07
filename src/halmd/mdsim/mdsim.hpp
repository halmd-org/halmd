/* Molecular Dynamics simulation of a Lennard-Jones fluid
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#ifndef HALMD_MDSIM_HPP
#define HALMD_MDSIM_HPP

#include <boost/assign.hpp>
#include <boost/multi_array.hpp>
#include <boost/noncopyable.hpp>
#include <boost/unordered_map.hpp>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <vector>
#include <unistd.h>

#include <halmd/mdsim/impl.hpp>
#include <halmd/mdsim/sample.hpp>
#include <halmd/options.hpp>
#include <halmd/sample/H5param.hpp>
#include <halmd/sample/correlation.hpp>
#include <halmd/sample/energy.hpp>
#include <halmd/sample/perf.hpp>
#include <halmd/sample/trajectory.hpp>
#include <halmd/util/H5xx.hpp>
#include <halmd/util/log.hpp>
#include <halmd/util/signal.hpp>
#include <halmd/util/timer.hpp>
#include <halmd/version.h>

namespace halmd
{

#define IMPL(x) typename mdsim_backend::impl_type::impl_##x()

/**
 * Molecular Dynamics simulation of a Lennard-Jones fluid
 */
template <typename mdsim_backend>
class mdsim : boost::noncopyable
{
public:
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
    int operator()();
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
        boost::array<float, 3> r_cut;
        try {
            r_cut = m_opt["cutoff"].as<boost::array<float, 3> >();
        }
        catch (boost::bad_any_cast const&) {
            // backwards compatibility
            r_cut[0] = r_cut[1] = r_cut[2] = m_opt["cutoff"].as<float>();
        }
        m_fluid.cutoff_radius(r_cut);
    }

    void potential_smoothing(boost::false_type const&) {}
    void potential_smoothing(boost::true_type const&)
    {
        if (!m_opt["smooth"].empty() && m_opt["smooth"].as<float>() > 0) {
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
        if (!m_opt["measure-temperature-after-time"].empty()) {
            m_en.open(m_opt["energy"].as<std::string>(), energy<dimension>::in);
            accumulator<double> en;
            m_en.temperature(en, m_opt["measure-temperature-after-time"].as<double>());
            LOG("temperature: " << en.mean() << " (" << en.std() << ", " << en.count() << " samples)");
            m_en.close();
            // rescale velocities to yield desired temperature
            m_fluid.rescale_velocities(std::sqrt(m_opt["temperature"].as<float>() / en.mean()));
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

    void pause(boost::false_type const&)
    {
        signal::wait();
    }

#ifdef WITH_CUDA
    void pause(boost::true_type const&)
    {
#ifndef __DEVICE_EMULATION__
        // pop current context from context stack
        cuda::driver::context::floating ctx;
#endif
        // pause simulation
        signal::wait();
        // push floating context onto context stack
    }
#endif /* WITH_CUDA */

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

    /** current MD step */
    count_timer<uint64_t> step;
};

/**
 * initialise MD simulation
 */
template <typename mdsim_backend>
mdsim<mdsim_backend>::mdsim(options const& opt) : m_opt(opt)
{
    LOG("positional coordinates dimension: " << static_cast<unsigned>(dimension));

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
    epsilon(IMPL(lennard_jones_potential));
    // collision diameters
    sigma(IMPL(lennard_jones_potential));
    // potential cutoff radius
    cutoff_radius(IMPL(lennard_jones_potential));
    // potential smoothing function scale parameter
    potential_smoothing(IMPL(lennard_jones_potential));
    // pair separation at which particle collision occurs
    pair_separation(IMPL(hardsphere_potential));

    // desired average cell occupancy
    cell_occupancy(IMPL(fixed_size_cell_lists));
    // neighbour list skin
    nbl_skin(IMPL(neighbour_lists));
    // number of CUDA execution threads
    threads(IMPL(gpu));
    // initialise hard-sphere cell lists
    init_cells(IMPL(hardsphere_cell_lists));

    // heat bath collision probability and temperature
    thermostat(IMPL(thermostat));

    if (m_opt["random-seed"].empty()) {
        m_fluid.rng(read_random_seed(m_opt["random-device"].as<std::string>()));
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
    rescale_energy(IMPL(lennard_jones_potential));

    if (m_opt["trajectory-sample"].empty() || (m_opt["measure-temperature-after-time"].empty() && !m_opt["temperature"].defaulted())) {
        // initialise velocities from Maxwell-Boltzmann distribution
        m_fluid.temperature(m_opt["temperature"].as<float>());
    }
    // initialise hard-sphere event list
    init_event_list(IMPL(hardsphere_event_lists));

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
    if (!m_opt["max-samples"].empty()) {
        boost::multi_array<uint64_t, 1> max_samples(m_opt["max-samples"].as<boost::multi_array<uint64_t, 1> >());
        m_corr.max_samples(std::vector<uint64_t>(max_samples.begin(), max_samples.end()));
    }
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
            m_is_corr_sample_gpu = IMPL(trajectory_gpu_sample);
        }
#endif
        else {
            throw std::logic_error("unknown correlation function backend: " + backend);
        }

        if (!m_opt["fastest-particle-fraction"].empty()) {
            boost::multi_array<float, 1> fastest_particle_fraction(
                m_opt["fastest-particle-fraction"].as<boost::multi_array<float, 1> >()
            );
            std::for_each(
                fastest_particle_fraction.begin()
              , fastest_particle_fraction.end()
              , boost::bind(&correlation<dimension>::add_fastest_particle_vacf_filter, &m_corr, _1)
            );
        }
        if (!m_opt["slowest-particle-fraction"].empty())
        {
            boost::multi_array<float, 1> slowest_particle_fraction(
                m_opt["slowest-particle-fraction"].as<boost::multi_array<float, 1> >()
            );
            std::for_each(
                slowest_particle_fraction.begin()
              , slowest_particle_fraction.end()
              , boost::bind(&correlation<dimension>::add_slowest_particle_vacf_filter, &m_corr, _1)
            );
        }
    }
}

/**
 * run MD simulation
 */
template <typename mdsim_backend>
int mdsim<mdsim_backend>::operator()()
{
    /** HDF5 buffers flush to disk interval in seconds */
    enum { FLUSH_TO_DISK_INTERVAL = 900 };
    /** runtime estimate interval in seconds */
    enum { TIME_ESTIMATE_INTERVAL = 1800 };
    /** waiting time in seconds before runtime estimate after block completion */
    enum { TIME_ESTIMATE_WAIT_AFTER_BLOCK = 300 };

    /** exit code */
    int status_ = HALMD_EXIT_SUCCESS;

    boost::unordered_map<int, std::string> const signals = boost::assign::map_list_of
        (SIGINT, "INT")
        (SIGTERM, "TERM")
        (SIGHUP, "HUP")
        (SIGUSR1, "USR1")
        (SIGTSTP, "TSTP");

    // open HDF5 output files
    open();
    // schedule first disk flush
    alarm(FLUSH_TO_DISK_INTERVAL);

    LOG("starting MD simulation");
    real_timer timer;
    timer.start();

    for (this->step = 0; this->step < m_corr.steps(); ++this->step) {
        sample_fluid(this->step, !this->step);

        if (sample_properties(this->step, false)) {
            // acquired maximum number of samples for a block level
            flush();
            this->step.set(TIME_ESTIMATE_WAIT_AFTER_BLOCK);
            alarm(FLUSH_TO_DISK_INTERVAL);
        }

        try {
            // MD integration step
            m_fluid.mdstep();
        }
        catch (potential_energy_divergence const& e) {
            this->step++;
            LOG_ERROR(e.what() << " at step " << this->step);
            status_ = HALMD_EXIT_POTENTIAL_ENERGY_DIVERGENCE;
            break;
        }

        if (this->step.estimate() > 0) {
            LOG(real_timer::format(this->step.estimate()) << " estimated remaining runtime at step " << this->step);
            this->step.set(TIME_ESTIMATE_INTERVAL);
        }

        if (signal::poll()) {
            if (signals.find(signal::signal) != signals.end()) {
                std::string const& name = signals.at(signal::signal);
                LOG_WARNING("trapped signal " + name + " at simulation step " << this->step);
            }
            if (signal::signal == SIGTSTP) {
                // block process until further signal is received
                LOG("pausing simulation");
                unsigned int seconds = alarm(0);
                pause(IMPL(gpu));
                alarm(seconds);
                LOG("resuming simulation");
            }
            if (signal::signal == SIGINT || signal::signal == SIGTERM) {
                this->step++;
                LOG_WARNING("aborting simulation at step " << this->step);
                status_ = HALMD_EXIT_TERM;
                break;
            }
            else if (signal::signal == SIGHUP || signal::signal == SIGALRM) {
                flush();
                alarm(FLUSH_TO_DISK_INTERVAL);
            }
            else if (signal::signal == SIGUSR1) {
                this->step.set(0);
            }
        }
    }
    sample_fluid(this->step, true);
    sample_properties(this->step, true);

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

    return status_;
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

    m_is_en_step = (m_corr.is_sample_step(step) || m_corr.is_trajectory_step(step)) && m_en.is_open();

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
        m_en.open(fn + ".tep", energy<dimension>::out);
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
        param(m_corr);
        m_corr.close();
    }
    if (m_traj.is_open()) {
        param(m_traj);
        m_traj.close();
    }
    if (m_en.is_open()) {
        param(m_en);
        m_en.close();
    }
    if (m_perf.is_open()) {
        param(m_perf);
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
        throw logic_error("failed to read from " + fn);
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
    node["backend"] = MDSIM_BACKEND;
    if (!m_opt["disable-correlation"].as<bool>()) {
        node["tcf_backend"] = m_opt["tcf-backend"].as<std::string>();
    }
    node["dimension"] = static_cast<unsigned>(dimension);
    if (this->step > 0) {
        node["effective_steps"] = static_cast<uint64_t>(this->step);
    }

    node = param["program"];
    node["name"] = PROGRAM_NAME;
    node["version"] = PROGRAM_VERSION;
    node["variant"] = PROGRAM_VARIANT;

    param << m_fluid << m_corr;
}

} // namespace halmd

#undef IMPL

#endif /* ! HALMD_MDSIM_HPP */
