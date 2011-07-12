Benchmarks
**********

Simple Lennard-Jones fluid in 3 dimensions
==========================================

Parameters:

    * 64,000 particles, number density :math:`\rho = 0.2\sigma^3`
    * force: lennard_jones_simple (:math:`r_c = 2.5\sigma, r_\text{skin} = 0.5\sigma`)
    * integrator: verlet (NVE, :math:`\delta t^* = 0.002`)

+--------------------+-------------------------------+------------------+---------------+-----------------------+
| Hardware           | time per MD step and particle | steps per second | FP precision  | compilation details   |
+====================+===============================+==================+===============+=======================+
| Intel Xeon E5540   | 626 ns                        | 24.9             | double        | GCC 4.3.4, -O3        |
+--------------------+-------------------------------+------------------+---------------+-----------------------+
| NVIDIA Tesla S1070 | 38.9 ns                       | 402              | double-single | CUDA 3.2, -arch sm_12 |
|                    +-------------------------------+------------------+---------------+-----------------------+
|                    | 33.6 ns                       | 465              | single        | CUDA 3.2, -arch sm_12 |
+--------------------+-------------------------------+------------------+---------------+-----------------------+
| NVIDIA Tesla S2050 | 30.4 ns                       | 514              | double-single | CUDA 3.2, -arch sm_12 |
|                    +-------------------------------+------------------+---------------+-----------------------+
|                    | 27.1 ns                       | 577              | single        | CUDA 3.2, -arch sm_12 |
+--------------------+-------------------------------+------------------+---------------+-----------------------+

Results were obtained from 1 independent measurement and are based on commit
8eb50a7. Each run consisted of NVT equilibration at :math:`T^*=1.5` over
:math:`\Delta t^*=100` (10⁴ steps), followed by benchmarking 5 times 10⁴ NVE
steps in a row.

Supercooled binary mixture (Kob-Andersen)
=========================================

Parameters:

    * 64,000 particles, number density :math:`\rho = 1.2\sigma^3`
    * force: lennard_jones with 2 particle spezies :math:`A, B`

      (:math:`\epsilon_{AA}=1`, :math:`\epsilon_{AB}=.5`, :math:`\epsilon_{BB}=1.5`,
      :math:`\sigma_{AA}=1`, :math:`\sigma_{AB}=.88`, :math:`\sigma_{BB}=.8`,
      :math:`r_c = 2.5\sigma`, :math:`r_\text{skin} = 0.5\sigma`)

    * integrator: verlet (NVE, :math:`\delta t^* = 0.001`)

+--------------------+-------------------------------+------------------+---------------+-----------------------+
| Hardware           | time per MD step and particle | steps per second | FP precision  | compilation details   |
+====================+===============================+==================+===============+=======================+
| Intel Xeon E5520   | 3035 ns                       | 5.15             | double        | GCC 4.4.1, -O3        |
+--------------------+-------------------------------+------------------+---------------+-----------------------+
| NVIDIA Tesla S1070 | 75.3 ns                       | 208              | double-single | CUDA 3.2, -arch sm_12 |
+--------------------+-------------------------------+------------------+---------------+-----------------------+

Results were obtained from 1 independent measurement and are based on commit
94928a1. Each run consisted of NVT equilibration at :math:`T^*=1.0` over
:math:`\Delta t^*=50` (10⁴ steps), followed by benchmarking 5 times 10⁴ NVE
steps in a row.
