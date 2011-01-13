/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#ifndef HALMD_TEST_MODULES_HPP
#define HALMD_TEST_MODULES_HPP

#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <exception>
#include <string>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/forces/lennard_jones.hpp>
#include <halmd/mdsim/host/forces/zero.hpp>
#include <halmd/mdsim/host/integrators/verlet.hpp>
#include <halmd/mdsim/host/integrators/verlet_nvt_andersen.hpp>
#include <halmd/mdsim/host/integrators/verlet_nvt_hoover.hpp>
#include <halmd/mdsim/host/neighbour.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/observables/host/trajectory.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/read_integer.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/forces/lennard_jones.hpp>
# include <halmd/mdsim/gpu/forces/zero.hpp>
# include <halmd/mdsim/gpu/integrators/verlet.hpp>
# include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen.hpp>
# include <halmd/mdsim/gpu/integrators/verlet_nvt_hoover.hpp>
# include <halmd/mdsim/gpu/neighbour.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/observables/gpu/thermodynamics.hpp>
# include <halmd/observables/gpu/trajectory.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif

namespace halmd
{
namespace test
{

#ifdef WITH_CUDA
boost::shared_ptr<utility::gpu::device> make_device(std::string const& backend)
{
    if (backend == "gpu") {
        static boost::weak_ptr<utility::gpu::device> device;
        boost::shared_ptr<utility::gpu::device> device_(device.lock());
        if (!device_) {
            device_ = boost::make_shared<utility::gpu::device>(
                std::vector<int>()   // devices
              , 128                  // threads
            );
            device = device_;
        }
        return device_;
    }
    if (backend == "host") {
        return boost::shared_ptr<utility::gpu::device>(); // null pointer
    }
    throw std::runtime_error("unknown backend: " + backend);
}
#endif /* WITH_CUDA */

template <int dimension>
boost::shared_ptr<mdsim::particle<dimension> > make_particle(
    std::string const& backend
  , unsigned int npart
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return boost::make_shared<mdsim::gpu::particle<dimension, float> >(
            make_device(backend)
          , std::vector<unsigned int>(1, npart)
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return boost::make_shared<mdsim::host::particle<dimension, double> >(
            std::vector<unsigned int>(1, npart)
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

boost::shared_ptr<halmd::random::random> make_random(
    std::string const& backend
  , unsigned int seed
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
        return boost::make_shared<random_type>(
            make_device(backend)
          , seed
          , 32                  // blocks
          , 32 << DEVICE_SCALE  // threads
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return boost::make_shared<halmd::random::host::random>(
            seed
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

boost::shared_ptr<halmd::random::random> make_random(
    std::string const& backend
  , std::string const& filename
)
{
    return make_random(backend, read_integer<unsigned int>(filename));
}

template <int dimension>
boost::shared_ptr<mdsim::box<dimension> > make_box(
    boost::shared_ptr<mdsim::particle<dimension> > particle
  , double density
  , fixed_vector<double, dimension> ratios = fixed_vector<double, dimension>(1)
)
{
    return boost::make_shared<mdsim::box<dimension> >(particle, density, ratios);
}


template <int dimension, typename vector_type>
boost::shared_ptr<mdsim::position<dimension> > make_lattice(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
  , boost::shared_ptr<halmd::random::random> random
  , vector_type const& slab
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return boost::make_shared<mdsim::gpu::positions::lattice<dimension, float, halmd::random::gpu::rand48> >(
            boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , boost::dynamic_pointer_cast<halmd::random::gpu::random<halmd::random::gpu::rand48> >(random)
          , slab
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return boost::make_shared<mdsim::host::positions::lattice<dimension, double> >(
            boost::dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , boost::dynamic_pointer_cast<halmd::random::host::random>(random)
          , slab
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

template <int dimension>
boost::shared_ptr<mdsim::position<dimension> > make_lattice(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
  , boost::shared_ptr<halmd::random::random> random
)
{
    typedef typename mdsim::box<dimension>::vector_type vector_type;
    return make_lattice(backend, particle, box, random, vector_type(1));
}

template <int dimension, typename float_type>
boost::shared_ptr<observables::trajectory<dimension> > make_trajectory_host(
    boost::shared_ptr<observables::host::samples::trajectory<dimension, float_type> > sample
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
)
{
    return boost::make_shared<observables::host::trajectory<dimension, float_type> >(
        sample
      , boost::dynamic_pointer_cast<mdsim::host::particle<dimension, float_type> >(particle)
      , box
    );
}

#ifdef WITH_CUDA
template <int dimension, typename float_type>
boost::shared_ptr<observables::trajectory<dimension> > make_trajectory_gpu(
    boost::shared_ptr<observables::host::samples::trajectory<dimension, float_type> > sample
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
)
{
    return boost::make_shared<observables::gpu::trajectory<observables::host::samples::trajectory<dimension, float_type> > >(
        sample
      , boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float_type> >(particle)
      , box
    );
}
#endif

template <int dimension>
boost::shared_ptr<mdsim::integrator<dimension> > make_verlet_integrator(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
  , double timestep
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return boost::make_shared<mdsim::gpu::integrators::verlet<dimension, float> >(
            boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , timestep
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return boost::make_shared<mdsim::host::integrators::verlet<dimension, double> >(
            boost::dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , timestep
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

template <int dimension>
boost::shared_ptr<mdsim::integrator<dimension> > make_verlet_nvt_andersen_integrator(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
  , boost::shared_ptr<halmd::random::random> random
  , double timestep
  , double temperature
  , double collision_rate
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return boost::make_shared<mdsim::gpu::integrators::verlet_nvt_andersen<dimension, float, halmd::random::gpu::rand48> >(
            boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , boost::dynamic_pointer_cast<halmd::random::gpu::random<halmd::random::gpu::rand48> >(random)
          , timestep, temperature, collision_rate
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return boost::make_shared<mdsim::host::integrators::verlet_nvt_andersen<dimension, double> >(
            boost::dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , boost::dynamic_pointer_cast<halmd::random::host::random>(random)
          , timestep, temperature, collision_rate
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

template <int dimension, typename float_type>
boost::shared_ptr<mdsim::integrator<dimension> > make_verlet_nvt_hoover_integrator(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
  , double timestep
  , double temperature
  , fixed_vector<double, 2> mass
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return boost::make_shared<mdsim::gpu::integrators::verlet_nvt_hoover<dimension, float_type> >(
            boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , timestep, temperature, mass
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return boost::make_shared<mdsim::host::integrators::verlet_nvt_hoover<dimension, float_type> >(
            boost::dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , timestep, temperature, mass
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

template <int dimension>
boost::shared_ptr<mdsim::force<dimension> > make_lennard_jones_force(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
  , boost::array<float, 3> cutoff
  , boost::array<float, 3> epsilon
  , boost::array<float, 3> sigma
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        typedef mdsim::gpu::forces::lennard_jones<float> potential_type;
        typedef mdsim::gpu::forces::pair_trunc<dimension, float, potential_type> force_type;
        return boost::make_shared<force_type>(
            boost::make_shared<potential_type>(particle->ntype, cutoff, epsilon, sigma)
          , boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        typedef mdsim::host::forces::lennard_jones<double> potential_type;
        typedef mdsim::host::forces::pair_trunc<dimension, double, potential_type> force_type;
        return boost::make_shared<force_type>(
            boost::make_shared<potential_type>(particle->ntype, cutoff, epsilon, sigma)
          , boost::dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

template <int dimension>
boost::shared_ptr<mdsim::force<dimension> > make_zero_force(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return boost::make_shared<mdsim::gpu::forces::zero<dimension, float> >(
            boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return boost::make_shared<mdsim::host::forces::zero<dimension, double> >(
            boost::dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

template <int dimension>
boost::shared_ptr<mdsim::neighbour<dimension> > make_neighbour(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
  , boost::shared_ptr<mdsim::force<dimension> > force
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        typedef mdsim::gpu::forces::lennard_jones<float> potential_type;
        typedef mdsim::gpu::forces::pair_trunc<dimension, float, potential_type> force_type;
        return boost::make_shared<mdsim::gpu::neighbour<dimension, float> >(
            boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , boost::dynamic_pointer_cast<force_type>(force)->r_cut()
          , 0.5         // skin
          , 0.4         // cell occupancy
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        typedef mdsim::host::forces::lennard_jones<double> potential_type;
        typedef mdsim::host::forces::pair_trunc<dimension, double, potential_type> force_type;
        return boost::make_shared<mdsim::host::neighbour<dimension, double> >(
            boost::dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , boost::dynamic_pointer_cast<force_type>(force)->r_cut()
          , 0.5         // skin
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

template <int dimension>
boost::shared_ptr<mdsim::velocity<dimension> > make_boltzmann(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<halmd::random::random> random
  , double temperature
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return boost::make_shared<mdsim::gpu::velocities::boltzmann<dimension, float, halmd::random::gpu::rand48> >(
            boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , boost::dynamic_pointer_cast<halmd::random::gpu::random<halmd::random::gpu::rand48> >(random)
          , temperature
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return boost::make_shared<mdsim::host::velocities::boltzmann<dimension, double> >(
            boost::dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , boost::dynamic_pointer_cast<halmd::random::host::random>(random)
          , temperature
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

template <int dimension>
boost::shared_ptr<observables::thermodynamics<dimension> > make_thermodynamics(
    std::string const& backend
  , boost::shared_ptr<mdsim::particle<dimension> > particle
  , boost::shared_ptr<mdsim::box<dimension> > box
  , boost::shared_ptr<mdsim::force<dimension> > force
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return boost::make_shared<observables::gpu::thermodynamics<dimension, float> >(
            boost::dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , boost::dynamic_pointer_cast<mdsim::gpu::force<dimension, float> >(force)
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return boost::make_shared<observables::host::thermodynamics<dimension, double> >(
            boost::dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , boost::dynamic_pointer_cast<mdsim::host::force<dimension, double> >(force)
        );
    }
    throw std::runtime_error("unknown backend: " + backend);
}

}  // namespace test

}  // namespace halmd

#endif  /* ! HALMD_TEST_MODULES_HPP */

