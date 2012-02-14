FAQ
***

**With HALMD, I could do simulations of breath-taking quality and have obtained new scientific insight. How may I thank you?**
  Please acknowledge the use of HALMD in your publications by citing our article:

   P. H. Colberg and F. Höfling,
   *Highly accelerated simulations of glassy dynamics using GPUs: Caveats on limited floating-point precision,*
   Comp. Phys. Comm. *182*, 1120 (2011) `[Link] <http://dx.doi.org/10.1016/j.cpc.2011.01.009>`_.

**Why does HALMD abort with “[ERROR] overcrowded placeholders in cell/neighbour lists update”?**
  The GPU backends of HALMD employ fixed-size cell and neighbour lists to
  accommodate for memory access pattern constraints.
  By default, a cell or neighbour list is allocated twice as many placeholders
  as the average expected number of particles. This may not suffice in case
  of large local density fluctuations, e.g. for unequilibrated systems.

  Try lowering the ``--cell-occupancy`` value.

**Why does HALMD abort with “[ERROR] potential energy diverged”?**
  An infinite potential energy sum of one or more particles indicates that the
  integration time-step is too large.

  Try lowering the ``--timestep`` value.

**nvcc fails with error: inline function ‘__signbit’ cannot be declared weak**
  CUDA 2.3 (or less) is not compatible with GCC 4.4.
  As a work around install GCC 4.3 and place symlinks in the default CUDA
  compiler directory, e.g. if the CUDA toolkit is located in ``/opt/cuda``,
  symlink ``/opt/cuda/bin/{gcc,g++}`` to ``/usr/bin/{gcc,g++}-4.3``, respectively.
  The compiler directory may be overriden with the ``--compiler-bindir`` option.

