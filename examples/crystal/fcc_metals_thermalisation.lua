#!/usr/bin/env halmd
--
-- Copyright © 2023 Felix Höfling
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
local constants = halmd.utility.constants
local log = halmd.io.log
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local writers = halmd.io.writers
local utility = halmd.utility

-- The next line is not needed if the definition files are located
-- in the same folder as the simulation script. The regular expression
-- is used to construct a path relative to the current script.
package.path = arg[0]:match("@?(.*/)") .. "../?.lua;" .. package.path
local definitions = {
    fcc_metals = require("definitions/fcc_metals_morse_MacDonaldMacDonald_1981")
}

--
-- Setup and run simulation
--
function main(args)
    -- get potential parameters
    local substance = args.substance
    local parameters = definitions.fcc_metals.parameters[substance]
    local units = definitions.fcc_metals.units
    log.message(("using potential parameters for %s (%s)"):format(substance, parameters[1]))

    -- query (nominal) lattice constant from potential parameters.
    -- r_min (or r₀) is the distance between nearest lattice points,
    -- so for the fcc lattice, we have a_lat = √2 r₀.
    local lattice_constant = math.sqrt(2) * parameters.r_min
    log.info(("lattice constant: %g Å"):format(lattice_constant))

    -- linear extend of cubic simulation box in multiples of the lattice constant
    local ncells = args.ncells
    local box_length = lattice_constant * ncells
    local nparticle = math.pow(ncells, 3) * 4       -- 4 atoms per fcc unit cell

    -- create cubic simulation domain with periodic boundary conditions
    local box = mdsim.box({length = {box_length, box_length, box_length}})

    -- create system state
    local particle = mdsim.particle({dimension = box.dimension, particles = nparticle})

    -- set particle masses
    local mass = parameters.mass
    log.info(("atomic mass: %g u"):format(mass * units.mass / constants.Da))

    -- set initial particle positions
    mdsim.positions.lattice({box = box, particle = particle})
       :set()

    -- set initial particle velocities
    --
    -- use twice the value of the target temperature since we start with a
    -- perfect lattice configuration which is force free (the potential energy
    -- fluctuation is zero)
    local temperature = constants.kB * args.temperature / units.energy    -- convert from Kelvin to simulation units (eV)
    log.message(("heat bath temperature: %g K = %g eV / kB"):format(args.temperature, temperature))
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle, temperature = 2 * temperature
    }):set()

    -- define Lennard-Jones pair potential (with parameters ε=1 and σ=1 for a single species)
    -- and register computation of pair forces
    definitions.fcc_metals.create_pair_force({
        substance = substance, box = box, particle = particle
      , cutoff = args.cutoff, smoothing = args.smoothing
      , neighbour = { disable_binning = true, disable_sorting = true, unroll_force_loop = true }
    })

    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- select all particles
    local all_group = mdsim.particle_groups.all({particle = particle})

    -- sample phase space
    local phase_space = observables.phase_space({box = box, group = all_group})
    -- write trajectory of particle groups to H5MD file
    local interval = args.sampling.trajectory or steps
    if interval > 0 then
        phase_space:writer({
            file = file, fields = {"position", "velocity"}, every = interval
        })
    end

    -- Sample macroscopic state variables.
    local interval = args.sampling.state_vars
    if interval > 0 then
        local msv = observables.thermodynamics({box = box, group = all_group})
        msv:writer({file = file, every = interval})
    end

    -- sample initial state
    observables.sampler:sample()

    -- add velocity-Verlet integrator with Boltzmann distribution
    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box, particle = particle, timestep = args.timestep
      , temperature = temperature, rate = args.rate
    })

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({steps = steps})

    -- run simulation
    observables.sampler:run(steps)
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time
      , default = "copper_thermalisation_rho{density:g}_T{temperature:.2f}_%Y%m%d_%H%M%S", help = "basename of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})

    parser:add_argument("random-seed", {type = "integer", action = parser.action.random_seed
      , help = "seed for random number generator"})

    -- build choices from definition file
    local choices = {}
    for k,v in pairs(definitions.fcc_metals.parameters) do
        choices[k] = v[1]
    end
    parser:add_argument("substance", {type = "string", required = true, choices = choices, help = "chemical symbol of the metallic substance"})
    parser:add_argument("ncells", {type = "integer", default = 10, help = "number of elementary fcc cells along each axis of the simulation box"})

    parser:add_argument("cutoff", {type = "float32", help = "potential cutoff radius"})
    parser:add_argument("smoothing", {type = "number", default = 0.01, help = "cutoff smoothing parameter"})
    parser:add_argument("temperature", {type = "number", default = 700, help = "thermostat temperature (Kelvin)"})
    parser:add_argument("rate", {type = "number", default = 2, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default = 100, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
end
