#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2011  Felix Höfling
#
# This file is part of HALMD.
#
# HALMD is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#
# Script for computing 1D density profile from simulation data of wetting process
#

import numpy as np
from tqdm import tqdm


def compute_density_modes_1D(wavevector: np.ndarray, positions: np.ndarray, coord_index: int) -> np.ndarray[complex]:
    """
    Compute time-averaged density modes in 1D.

    Args:
    - coord_index: index of coordinate of interest (coord_index ∈ {0, 1, 2})
    """

    density_modes =  np.zeros(wavevector.shape[0], dtype = complex)
    for i, k_i in enumerate(tqdm(wavevector)):
        density_modes[i] = np.sum(np.exp(1j * k_i[coord_index] * positions[:,:, coord_index]))

    return density_modes / positions.shape[0]


def wavevector_to_grid_indices(wavevector: np.ndarray[float], box_edges: np.ndarray[float]) -> np.ndarray[float]:
    """
    Map wave indices to the appropriate indices in a discrete Fourier grid.

    Args:
    - wavevector: the wave indices to map, given as an array of shape (n_vectors, 3)
    - box_edges: the physical dimensions of the real-space grid, given as (Lx, Ly, Lz)

    Returns:
    - grid_indices: an array containing the grid indices for all wavevectors
    """

    # Compute the grid spacing in reciprocal space
    dkx, dky, dkz =2 * np.pi / box_edges

    # Map the wavevector components to grid indices
    wx = np.array(np.round(wavevector[:,0] / dkx), dtype =int)
    wy = np.array(np.round(wavevector[:,1] / dky), dtype =int)
    wz = np.array(np.round(wavevector[:,2] / dkz), dtype =int)
    assert np.min(wx)==-np.max(wx) , "Density-modes need to be on a symmetric grid around 0, i.e. k_max = -k_min in x- direction"
    assert np.min(wy)==-np.max(wy) , "Density-modes need to be on a symmetric grid around 0, i.e. k_max = -k_min in y- direction"
    assert np.min(wz)==-np.max(wz) , "Density-modes need to be on a symmetric grid around 0, i.e. k_max = -k_min in z- direction"

    # Combine wx, wy, wz into a single array of shape (199, 3)
    grid_indices = np.column_stack((wx, wy, wz))

    return grid_indices


def compute_1D_density(density_modes: np.ndarray[complex], box_edges: np.ndarray[float], wavevector: np.ndarray[float], smoothing: bool = False, sigma: float = None) -> tuple:
    """
    Compute 1D density given the wave vectors and corresponding density modes.
    It's applicable when filter = {1,0,0} is set for wavevector in HALMD simulation.

    Args:
    - smoothing: if True, a Gaussian filter with standard deviation of sigma will be applied to the obtained density

    Returns:
    - real_grid: grid in real space that hold values of density
    - density: density normalized by volume for each point in real_grid
    """

    # Grid indices of a 1D grid
    grid_indices = wavevector_to_grid_indices(wavevector, box_edges)[:,0]

    if smoothing:
        # smooth out rapid oscillations with Gaussians
        k_max = np.max(np.linalg.norm(wavevector, axis=1))
        if not sigma:
            sigma = 4*np.pi/k_max
        gaussian = np.exp(-0.5 * sigma**2 * (np.linalg.norm(wavevector, axis=1))**2)
        density_modes = gaussian * (np.real(density_modes) + 1j*np.imag(density_modes))

    # shift density modes such that first entry contains the zero frequency term
    density_modes = np.concatenate((density_modes[-len(density_modes)//2:], density_modes[:-len(density_modes)//2]))

    # Apply inverse Fourier Transform to obtain density
    density_1D = np.fft.fftshift(np.fft.ifft(density_modes)).real

    # normalise with bulk density
    indices_left, indices_right = int(len(density_1D)//2 - (1/10)*len(density_1D)), int(len(density_1D)//2 + (1/10)*len(density_1D))
    bulk_density = np.average(density_1D[indices_left:indices_right])
    density_1D /= bulk_density

    # enforce boundary conditions
    boundary_average = 0.5 * (density_1D[0] + density_1D[-1])
    density_1D[-1] = boundary_average
    density_1D[0] = boundary_average

    # generate the corresponding positions in real space
    real_grid = box_edges[0] / len(grid_indices)* np.linspace(np.min(grid_indices), np.max(grid_indices), len(grid_indices))

    return real_grid, density_1D