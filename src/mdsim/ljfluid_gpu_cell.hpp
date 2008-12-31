/* Lennard-Jones fluid simulation using CUDA
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

#ifndef MDSIM_LJFLUID_CELL_GPU_HPP
#define MDSIM_LJFLUID_CELL_GPU_HPP

#include <boost/array.hpp>
#include <cuda_wrapper.hpp>
#include "config.hpp"
#include "perf.hpp"
#include "radix.hpp"
#include "rand48.hpp"
#include "sample.hpp"

namespace mdsim
{

/**
 * Lennard-Jones fluid simulation using CUDA
 */
class ljfluid
{
public:
    /** set number of particles in system */
    void particles(unsigned int value);
    /** set particle density */
    void density(float value);
    /** set periodic box length */
    void box(float value);
    /** set potential cutoff radius */
    void cutoff_radius(float value);
#ifdef USE_POTENTIAL_SMOOTHING
    /** set potential smoothing function scale parameter */
    void potential_smoothing(float value);
#endif
#ifdef USE_CELL
    /** set desired average cell occupancy */
    void cell_occupancy(float value);
#endif
    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /** restore system state from phase space sample */
    void restore(trajectory_sample::visitor visitor);

    /** seed random number generator */
    void rng(unsigned int seed);
    /** restore random number generator from state */
    void rng(mdsim::rand48::state_type const& state);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(float temp);
    /** set simulation timestep */
    void timestep(float value);

    /** get number of particles */
    unsigned int const& particles() const { return npart; }
#ifdef USE_CELL
    /** get number of cells per dimension */
    unsigned int const& cells() const { return ncell; }
    /** get total number of cell placeholders */
    unsigned int const& placeholders() const { return nplace; }
    /** get cell length */
    float const& cell_length() const { return cell_length_; }
    /** get effective average cell occupancy */
    float const& cell_occupancy() const { return cell_occupancy_; }
    /** get number of placeholders per cell */
    unsigned int const& cell_size() const { return cell_size_; }
#endif
    /** get number of CUDA execution blocks */
    unsigned int blocks() const { return dim_.blocks_per_grid(); }
    /** get number of CUDA execution threads */
    unsigned int threads() const { return dim_.threads_per_block(); }
    /** get particle density */
    float const& density() const { return density_; }
    /** get periodic box length */
    float const& box() const { return box_; }
    /** get potential cutoff radius */
    float const& cutoff_radius() const { return r_cut; }
#ifdef USE_POTENTIAL_SMOOTHING
    /** get potential smoothing function scale parameter */
    float const& potential_smoothing() const { return r_smooth; }
#endif
    /** get simulation timestep */
    float const& timestep() const { return timestep_; }
    /** returns and resets CUDA time statistics */
    perf_counters times();

    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

    /** stream MD simulation step on GPU */
    void mdstep();
    /** synchronize MD simulation step on GPU */
    void synchronize();
    /** copy MD simulation step results from GPU to host */
    void sample();
    /** get trajectory sample */
    trajectory_sample const& trajectory() const { return h_sample; }

private:
    /** number of particles in system */
    unsigned int npart;
#ifdef USE_CELL
    /** number of cells per dimension */
    unsigned int ncell;
    /** total number of cell placeholders */
    unsigned int nplace;
    /** cell length */
    float cell_length_;
    /** effective average cell occupancy */
    float cell_occupancy_;
    /** number of placeholders per cell */
    unsigned int cell_size_;
#endif
    /** particle density */
    float density_;
    /** periodic box length */
    float box_;
    /** simulation timestep */
    float timestep_;
    /** cutoff radius for shifted Lennard-Jones potential */
    float r_cut;
    /** maximum velocity magnitude after last MD step */
    float v_max;
#ifdef USE_CELL
    /** cell skin */
    float r_skin;
    /** sum over maximum velocity magnitudes since last cell lists update */
    float v_max_sum;
#endif
#ifdef USE_POTENTIAL_SMOOTHING
    /** potential smoothing function scale parameter */
    float r_smooth;
#endif

    /** trajectory sample in swappable host memory */
    trajectory_sample h_sample;

#ifdef USE_CELL
    /** cell placeholders in page-locked host memory */
    struct {
	/** periodically reduced particle positions */
	cuda::host::vector<gvector> r;
	/** periodically extended particle positions */
	cuda::host::vector<gvector> R;
	/** particle velocities */
	cuda::host::vector<gvector> v;
	/** particle number tags */
	cuda::host::vector<int> n;
	/** potential energies per particle */
	cuda::host::vector<float> en;
	/** virial equation sums per particle */
	cuda::host::vector<float> virial;
    } h_cell;
#else
    /** system state in page-locked host memory */
    struct {
	/** periodically reduced particle positions */
	cuda::host::vector<gvector> r;
	/** periodically extended particle positions */
	cuda::host::vector<gvector> R;
	/** particle velocities */
	cuda::host::vector<gvector> v;
	/** potential energies per particle */
	cuda::host::vector<float> en;
	/** virial equation sums per particle */
	cuda::host::vector<float> virial;
    } h_part;
#endif

#ifdef USE_CELL
    /** system state in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<gvector> r;
	/** periodically extended particle positions */
	cuda::vector<gvector> R;
	/** particle velocities */
	cuda::vector<gvector> v;
	/** particle number tags */
	cuda::vector<int> n;
	/** particle forces */
	cuda::vector<gvector> f;
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_cell;

    /** system state double buffer in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<gvector> r;
	/** periodically extended particle positions */
	cuda::vector<gvector> R;
	/** particle velocities */
	cuda::vector<gvector> v;
	/** particle number tags */
	cuda::vector<int> n;
    } g_cell2;
#else
    /** system state in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<gvector> r;
	/** periodically extended particle positions */
	cuda::vector<gvector> R;
	/** particle velocities */
	cuda::vector<gvector> v;
	/** particle forces */
	cuda::vector<gvector> f;
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_part;
#endif

    /** GPU random number generator */
    mdsim::rand48 rng_;
    /** GPU mdsim::radix sort for particle positions */
    mdsim::radix_sort<gvector> radix_sort;
    /** CUDA execution dimensions */
    cuda::config dim_;
#ifdef USE_CELL
    /** CUDA execution dimensions for cell-specific kernels */
    cuda::config dim_cell_;
#endif
    /** CUDA asynchronous execution */
    cuda::stream stream_;
    /** CUDA events for kernel timing */
#ifdef USE_CELL
    boost::array<cuda::event, 5> event_;
#else
    boost::array<cuda::event, 3> event_;
#endif
    /** CUDA time statistics */
    perf_counters m_times;
};

} // namespace mdsim

#endif /* ! MDSIM_LJFLUID_CELL_GPU_HPP */
