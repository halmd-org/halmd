#!/usr/bin/env halmd
--
-- Copyright © 2022 Felix Höfling
--
-- This file is part of HALMD.
--
-- HALMD is free software: you can redistribute it and/or modify
-- it under the terms of the GNU Lesser General Public License as
-- published by the Free Software Foundation, either version 3 of
-- the License, or (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU Lesser General Public License for more details.
--
-- You should have received a copy of the GNU Lesser General
-- Public License along with this program.  If not, see
-- <http://www.gnu.org/licenses/>.
--

-- grab modules
local mdsim = halmd.mdsim
local observables = halmd.observables
local writers = halmd.io.writers
local utility = halmd.utility

--
-- Setup and run simulation
--
function main(args)
    local nparticle = 32000   -- total number of particles
    local density = 0.9       -- number density
    local temperature = 1.2   -- temperature of hot reservoir
    local cutoff = 2.5        -- potential cutoff

    local timestep = 0.002    -- integration timestep
    local steps = 100000      -- number of integration steps

    local thermo_width = 10   -- width of thermostatted region
    local delta_temp = 0.4    -- temperature difference

    -- derive edge length of three-dimensional box from particle number and density
    local length = math.pow(nparticle / density, 1 / 3)

    -- create simulation box
    local box = mdsim.box({length = {length, length, length}})

    -- create system state with 1 particle species
    local particle = mdsim.particle({dimension = 3, particles = nparticle})

    -- define group of all particles
    local all_group = mdsim.particle_groups.all({particle = particle})

    -- set initial particle positions
    mdsim.positions.lattice({box = box, particle = particle}):set()

    -- set initial particle velocities
    mdsim.velocities.boltzmann({particle = particle, group = all_group, temperature = temperature}):set()

    -- define interaction using truncated Lennard-Jones potential
    local potential = mdsim.potentials.pair.lennard_jones():truncate({cutoff = cutoff})
    -- compute forces
    local force = mdsim.forces.pair({
        box = box, particle = particle, potential = potential,
    })

    -- melt initial fcc crystal
    local integrator = mdsim.integrators.verlet_nvt_boltzmann_value({
        box = box, particle = particle, group = all_group, mean_velocity = {0, 0, 0}
      , timestep = timestep, temperature = temperature, rate = 5
    })
    observables.sampler:run(steps / 10)
    integrator:disconnect()

    -- define regions and particle groups for integrators
    local lowest_corner = box:lowest_corner()
    local cold_group = mdsim.particle_groups.region({
        particle = particle, box = box, selection = "included", label = "cold"
      , geometry = mdsim.geometries.cuboid({
            lowest_corner = lowest_corner
          , length = {thermo_width, length, length}
        })
    })

    lowest_corner[1] = lowest_corner[1] + thermo_width
    local nve_group = mdsim.particle_groups.region({
        particle = particle, box = box, selection = "included", label = "NVE"
      , geometry = mdsim.geometries.cuboid({
            lowest_corner = lowest_corner
          , length = {length - 2 * thermo_width, length, length}
        })
    })

    lowest_corner[1] = length / 2 - thermo_width
    local hot_group = mdsim.particle_groups.region({
        particle = particle, box = box, selection = "included", label = "hot"
      , geometry = mdsim.geometries.cuboid({
            lowest_corner = lowest_corner
          , length = {thermo_width, length, length}
        })
    })

    -- define velocity-Verlet integrators, with Boltzmann thermostat in hot and cold regions
    mdsim.integrators.verlet({
        box = box, particle = particle, group = nve_group, timestep = timestep
    })

    mdsim.integrators.verlet_nvt_boltzmann_value({
        box = box, particle = particle, group = cold_group, mean_velocity = {0, 0, 0}
      , timestep = timestep, temperature = temperature - delta_temp, rate = 5
    })

    mdsim.integrators.verlet_nvt_boltzmann_value({
        box = box, particle = particle, group = hot_group, mean_velocity = {0, 0, 0}
      , timestep = timestep, temperature = temperature, rate = 5
    })

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output)})

    -- connect H5MD writer, output only first and last step
    observables.phase_space({box = box, group = all_group}):writer({
        file = file
      , fields = {"position", "velocity", "species", "mass"}
    })

    -- sample basic thermodynamical properties
    observables.thermodynamics({box = box, group = all_group})
       :writer({file = file, every = 100, fields = {"potential_energy", "pressure", "temperature", "kinetic_energy"}})
    observables.thermodynamics({box = box, group = nve_group})
       :writer({file = file, every = 100})
    observables.thermodynamics({box = box, group = cold_group})
       :writer({file = file, every = 100})
    observables.thermodynamics({box = box, group = hot_group})
       :writer({file = file, every = 100})

    -- sample temperature profile using kinetic energy modes

    -- define wave number grid and write density modes to H5MD file
    local kmax = 2 * math.pi / (.25)    -- spatial resolution: σ / 4
    local density_wavevector = observables.utility.wavevector({
        box = box, wavenumber = {kmax}, filter = {1, 0, 0}, dense = true
    })

    observables.density_mode({
        group = all_group, wavevector = density_wavevector
    }):writer({file = file, every = 250})

    observables.kinetic_energy_density_mode({
        group = all_group, wavevector = density_wavevector
    }):writer({file = file, every = 250})

    -- estimate remaining runtime
    observables.runtime_estimate({steps = steps, first = 10, interval = 900, sample = 60})

    -- sample initial state
    observables.sampler:sample()

    -- run simulation
    observables.sampler:run(steps)

    -- log profiler results
    utility.profiler:profile()
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "temperature_gradient_%Y%m%d_%H%M%S", help = "prefix of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})
end
