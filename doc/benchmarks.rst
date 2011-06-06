Benchmarks
**********

Simple Lennard-Jones fluid in 3 dimensions
==========================================

Parameters:

    * 64,000 particles, number density :math:`\rho = 0.2\sigma^3`
    * force: lennard_jones_simple (:math:`r_c = 2.5\sigma, r_\text{skin} = 0.5\sigma`)
    * integrator: verlet (NVE, double-single precision, :math:`\delta t^* = 0.002`)

+--------------------+-------------------------------+------------------+
| Hardware           | time per MD step and particle | steps per second |
+====================+===============================+==================+
| Intel Xeon E5520   |                               |                  |
+--------------------+-------------------------------+------------------+
| NVIDIA Tesla C1060 | 41.9 ns                       | 373              |
+--------------------+-------------------------------+------------------+
| NVIDIA Tesla C2050 | 31.5 ns                       | 491              |
+--------------------+-------------------------------+------------------+

Results were obtained from 1 independent measurement and are based on commit
31fc650, using CUDA 3.2. Each run consisted of NVT equilibration at
:math:`T^*=1.5` over :math:`\Delta t^*=100` (10⁴ steps), followed by
benchmarking 5 times 10⁴ NVE steps in a row.
