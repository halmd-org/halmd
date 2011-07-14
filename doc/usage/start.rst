Getting started
===============

Program parameters
------------------

HALMD has three ways of accepting program parameters:

* pass directly as command line options
* read from parameter config file [RC format]
* read from HDF5 data file, optionally resuming from a prior trajectory

Options are described in the command line help::

  halmd --help

The group of global options must precede module-specific options, which are
introduced by their respective module name.


A first example: Lennard-Jones fluids
-------------------------------------


Many simulation scenarios can be controlled by command line options, there is
no need to setup a separate simulation script. Let us consider a simple fluid
with 20,000 Lennard-Jones particles at density :math:`\rho^*=0.8`.
Equilibration is done with the Andersen thermostat at a temperature
:math:`T^*=2` over 10,000 steps ::

  halmd -v \
    --force=lennard_jones --integrator=verlet_nvt_andersen \
    --timestep 0.005 --time 50 \
    box --density 0.8 --particles 20000 \
    verlet_nvt_andersen --temperature 2 \
    sampler --trajectory 10000

Many parameters have sensible default values and may be omitted, e.g, the
collision rate of the Andersen thermostat (10), or the cutoff radius of the
potential (:math:`r_c=2.5\sigma`). The option ``-v`` makes the output more
verbose, check that your CUDA device has beed detected properly.

The equivalent parameter configuration file would look like this ::

  # file example.rc
  force=lennard_jones
  integrator=verlet_nvt_andersen
  timestep=0.005
  time 50

  [box]
  density=0.8
  particles=20000

  [verlet_nvt_andersen]
  temperature=2

  [sampler]
  trajectory=10000

It is passed directly to HALMD, while parameters may be overwritten: ::

  halmd -v -c example.rc --time 10

The trajectory is written every 10,000 steps, i.e., at the beginning and the
end of the simulation.  The initial configuration of the particles is a fcc
lattice. The default output settings yield an HDF5 file with a time stamp in
its name, ``halmd_%Y%m%d_%H%M%S.trj`` and a corresponding log file.

We may now continue the simulation at constant energy by resuming from the
trajectory file and selecting a different integrator ::

  halmd -v \
    --integrator=verlet --timestep 0.001 --time 100 \
    trajectory --file halmd_20110715_160545.trj
    sampler --state-vars 100

This will continue the simulation over 10⁵ steps and output thermodynamic
(macroscopic) state variables every 100 steps (potential energy, instantaneous
"temperature", pressure, ...) to a file with the extension ``msv``.

If the HDF5 tools are properly installed, we may have a quick overview of the output file ::

  h5ls halmd_20110715_160920.trj

or look at a specific data set ::

  h5dump -d EPOT halmd_20110713_161511.msv | less

Very convenient and powerful access to the HDF5 output files is provided by the
Python packages `h5py <http://alfven.org/wp/hdf5-for-python>`_ together with
`PyLab <http://www.scipy.org/PyLab>`_. An examplary Python script for
accessing, post-processing and plotting the output data is provided in
the sources at ``examples/plotting/plot_h5md.py``.
