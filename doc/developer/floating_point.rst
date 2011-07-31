Floating-point precision
========================

Most of the classes in halmd/mdsim carry a template parameter ``float_type``,
which allows to choose between single and double floating-point precision.
(This is not fully implemented yet.) In most cases, single precision will be
sufficient, notable exceptions are situations where many values are
accumulated. In particular, this is the case for the position and velocity
variables in the velocity-Verlet integrator, which a higher precision (double
or double-single) should always be used for.

The distinction of different floating-point precision ends at the level of the
phase space data, i.e., after the particle positions etc. have been copied to
``observables::samples::phase_space``. All derived observables are usually
accumulated quantities (e.g. mean kinetic energy or density modes) and shall be
computed with double precision. Although the relative statistical fluctuations
of these quantities will be much larger than 10⁻⁷, we shall stick to standard
usage and employ double precision for the subsequent analysis. Note that still
many classes in halmd/observables will carry the template parameter
``float_type``, which is, however, only used to specify the precision of the
input data, the results shall be in double precision in all cases.

At the level of file output (halmd/io), one may consider to write the final
results in single precision as an optimisation measure to save disk space. (The
high precision bits are mostly unneeded for questions in science, and such a
random noise is not compressed very well.) Again, all intermediate computations
from the particle coordinates down to the final results shall be done in double
precision to avoid potential artifacts.
