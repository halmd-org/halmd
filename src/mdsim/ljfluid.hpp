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

#ifndef MDSIM_LJFLUID_HPP
#define MDSIM_LJFLUID_HPP

#include <boost/array.hpp>
#include <cuda_wrapper.hpp>
#include <stdint.h>
#include "H5xx.hpp"
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
    enum {
	REDUCE_BLOCKS = 16,
	REDUCE_THREADS = 512,
    };

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
    /** get total number of placeholders per neighbour list */
    unsigned int const& neighbours() const { return nbl_size; }
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
    /** returns and resets GPU time statistics */
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
    /** first leapfrog step of integration of differential equations of motion */
    void velocity_verlet(cuda::stream& stream);
    /** Lennard-Jones force calculation */
    void update_forces(cuda::stream& stream);
    /** potential energy sum calculation */
    void potential_energy(cuda::stream& stream);
    /** virial equation sum calculation */
    void virial_sum(cuda::stream& stream);
#ifdef USE_CELL
    /** maximum velocity calculation */
    void maximum_velocity(cuda::stream& stream);
    /** assign particles to cells */
    void assign_cells(cuda::stream& stream);
    /** update neighbour lists */
    void update_neighbours(cuda::stream& stream);
#ifdef USE_HILBERT_ORDER
    /** order particles after Hilbert space-filling curve */
    void hilbert_order(cuda::stream& stream);
#endif
#endif

private:
    /** number of particles in system */
    unsigned int npart;
#ifdef USE_CELL
    /** number of cells per dimension */
    unsigned int ncell;
    /** total number of cell placeholders */
    unsigned int nplace;
    /** number of placeholders per neighbour list */
    unsigned int nbl_size;
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

    /** system state in page-locked host memory */
    struct {
	/** periodically reduced particle positions */
	cuda::host::vector<gvector> r;
	/** periodically extended particle positions */
	cuda::host::vector<gvector> R;
	/** particle velocities */
	cuda::host::vector<gvector> v;
	/** particle tags */
	cuda::host::vector<int> tag;
	/** blockwise maximum particles velocity magnitudes */
	cuda::host::vector<float> v_max;
	/** blockwise potential energies sum per particle */
	cuda::host::vector<float2> en_sum;
	/** blockwise virial equation sum per particle */
	cuda::host::vector<float2> virial_sum;
    } h_part;

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
	/** particle tags */
	cuda::vector<int> tag;
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
	/** blockwise maximum particles velocity magnitudes */
	cuda::vector<float> v_max;
	/** blockwise potential energies sum per particle */
	cuda::vector<float2> en_sum;
	/** blockwise virial equation sum per particle */
	cuda::vector<float2> virial_sum;
    } g_part;

#ifdef USE_CELL
    /** double buffers for particle sorting */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<gvector> r;
	/** periodically extended particle positions */
	cuda::vector<gvector> R;
	/** particle velocities */
	cuda::vector<gvector> v;
	/** particle tags */
	cuda::vector<int> tag;
    } g_sort;

    /** auxiliary device memory arrays for particle sorting */
    struct {
	/** particle cells */
	cuda::vector<uint> cell;
	/** cell offsets in sorted particle list */
	cuda::vector<int> offset;
	/** permutation indices */
	cuda::vector<int> idx;
    } g_aux;

    /** cell lists in global device memory */
    cuda::vector<int> g_cell;
    /** neighbour lists in global device memory */
    cuda::vector<int> g_nbl;
#endif

    /** GPU random number generator */
    mdsim::rand48 rng_;
#ifdef USE_CELL
    /** GPU mdsim::radix sort */
    mdsim::radix_sort<int> radix_sort;
#endif

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
    boost::array<cuda::event, 9> event_;
#else
    boost::array<cuda::event, 5> event_;
#endif
    /** GPU time statistics */
    perf_counters m_times;
};

} // namespace mdsim

#endif /* ! MDSIM_LJFLUID_HPP */
