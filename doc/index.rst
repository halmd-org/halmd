Introduction
============

HALMD is a high-precision molecular dynamics package for the large-scale
simulation of simple and complex liquids. HALMD supports acceleration through
CUDA-enabled graphics processors. HALMD lays emphasis on long-time numerical
stability by using double-single precision arithmetic for numerically crucial
parts of the MD integration step. HALMD is freely available and licensed under
the GPL.

HALMD is written and developed by Peter Colberg and Felix Höfling.

.. note::

   A description of the implementation, performance tests, numerical stability
   tests and an application to the glass transition of the Kob-Andersen mixture
   is found in the article by Peter H. Colberg and Felix Höfling, `Highly
   accelerated simulations of glassy dynamics using GPUs: Caveats on limited
   floating-point precision <http://dx.doi.org/10.1016/j.cpc.2011.01.009>`_,
   Computer Physics Communications (2011), doi:10.1016/j.cpc.2011.01.009
   [`arXiv:0912.3824 <http://arxiv.org/abs/0912.3824>`_].

The name *HALMD* was chosen in honour of the machine `HAL
<http://en.wikipedia.org/wiki/HAL_9000>`_ at the Arnold Sommerfeld Center
for Theoretical Physics of the Ludwigs-Maximilians-Universität München.
HAL has been the project's first GPGPU machine, equipped initially with two
NVIDIA GeForce 8800 Ultra. HAL survived a critical air condition failure and
to this day houses two NVIDIA GeForce GTX 280.

Features
========

HAL's MD package features

* GPU acceleration by up to a factor of 80 (with GeForce GTX 280),

* double-single arithmetic for excellent numerical long-time stability,

* potential smoothing for optimal energy conservation,

* online evaluation of various dynamic correlation functions,

* simulations at constant energy in the NVE ensemble,

* Andersen thermostat for NVT simulations.


Physics applications
---------------------

HAL's MD package is designed to study

* the dynamics of simple liquids and binary mixtures,

* particles interacting via attractive or purely repulsive potentials
  (truncated and shifted Lennard-Jones potentials),

* both two- and three-dimensional systems.


Technical features
------------------

HAL's MD package employs

* simulation parameters passed via command line options, INI-style
  configuration files and HDF5 output files,

* compressed and structured output of trajectories, thermodynamic variables and
  correlation functions in `HDF5 <http://hdfgroup.org/HDF5/>`_ format,

* a modular design with various backends loaded as dynamic libraries or
  separate executables,

* easily extensible template-based C++ code.


Getting HAL's MD package
========================

Please refer to the :ref:`Download section <download>`.

