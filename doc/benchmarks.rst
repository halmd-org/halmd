Benchmarks
**********

Simple Lennard-Jones fluid in 3 dimensions
==========================================

Parameters:

    * 64,000 particles, number density :math:`\rho = 0.2\sigma^3`
    * force: lennard_jones (:math:`r_c = 2.5\sigma, r_skin = 0.5\sigma`)
    * integrator: verlet (NVE, double-single precision, :math:`\delta t^* = 0.002`)

+--------------------+-------------------------------+------------------+
| Hardware           | time per MD step and particle | steps per second |
+====================+===============================+==================+
| Intel Xeon E5520   |                               |                  |
+--------------------+-------------------------------+------------------+
| NVIDIA Tesla C1060 | 46.4 ± 1 ns                   | 337 ± 8          |
+--------------------+-------------------------------+------------------+
| NVIDIA Tesla C2050 | 33.0 ± 0.7 ns                 | 474 ± 10         |
+--------------------+-------------------------------+------------------+

Results were obtained from 3 independent measurements and are based on commit
0410a09, using CUDA 3.2. Each run consisted of NVT equilibration at
:math:`T^*=1.5` over :math:`\Delta t^*=100` (10⁴ steps), followed by
benchmarking 5 times 10⁴ NVE steps in a row.
