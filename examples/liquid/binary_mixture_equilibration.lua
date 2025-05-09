#!/usr/bin/env halmd
--
-- Copyright © 2011-2023 Felix Höfling
-- Copyright © 2010-2012 Peter Colberg
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
local log = halmd.io.log
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local random = halmd.random
local writers = halmd.io.writers
local utility = halmd.utility

-- search definition files in the top-level path relative to the simulation script
package.path = utility.abspath("../?.lua;") .. package.path
local definitions = { kob_andersen = require("definitions/kob_andersen") }

--
-- Setup and run simulation
--
function main(args)
    -- total number of particles from sum of particles per species
    local nspecies = #args.particles
    local nparticle = numeric.sum(args.particles)

    -- derive edge lengths from number of particles, density and edge ratios
    local volume = nparticle / args.density
    local ratios = args.ratios
    local dimension = #ratios
    local unit_edge = math.pow(volume / numeric.prod(ratios), 1 / dimension)
    local length = {}
    for i = 1, #ratios do
        length[i] = unit_edge * ratios[i]
    end
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle, species = nspecies})
    -- set particle species, with continuous range of IDs per species
    local species = {}
    local groups = {}
    local offset = 0
    for s = 0, nspecies - 1 do
        local nparticle = assert(args.particles[s + 1])
        for i = 1, nparticle do
            table.insert(species, s)
        end
        local label = string.char(string.byte("A") + s)
        groups[label] = mdsim.particle_groups.id_range({
            particle = particle
          , range = {offset + 1, offset + nparticle}
          , label = label
        })
        offset = offset + nparticle
    end
    particle.data["species"] = species
    -- set initial particle positions sequentially on an fcc lattice
    local lattice = mdsim.positions.lattice({box = box, particle = particle})
    lattice:set()
    -- randomly shuffle the positions
    particle.data["position"] = random.generator({memory = "host"}):shuffle(particle.data["position"])
    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle
      , temperature = args.initial_temperature
    })
    boltzmann:set()

    -- define pair potential of Kob-Andersen binary mixture
    -- and register computation of pair forces
    definitions.kob_andersen.create_pair_force({
        box = box, particle = particle
    })

    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- sample macroscopic state variables for all particles.
    if args.sampling.state_vars > 0 then
        local msv = observables.thermodynamics({
            box = box
          , group = mdsim.particle_groups.all({particle = particle})
        })
        msv:writer({file = file, every = args.sampling.state_vars})
    end

    -- sample each particle group separately
    for label, group in utility.sorted(groups) do
        -- sample phase space
        local phase_space = observables.phase_space({box = box, group = group})
        -- write trajectory of particle groups to H5MD file
        phase_space:writer({
            file = file
          , fields = {"position", "velocity", "species", "mass"}
          , every = args.sampling.trajectory or steps
        })

        -- Sample macroscopic state variables per particle group.
        local msv = observables.thermodynamics({box = box, group = group})
        msv:writer({file = file, every = args.sampling.state_vars})
    end

    -- sample initial state
    observables.sampler:sample()

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({steps = steps})

    -- add velocity-Verlet integrator with Boltzmann thermostat (NVT)
    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box, particle = particle
      , timestep = args.timestep, temperature = args.initial_temperature, rate = args.rate
    })

    -- run first 10% of the simulation in NVT ensemble at elevated temperature
    -- in order to melt the initial fcc crystal
    observables.sampler:run(steps / 10)

    -- run remaining first half of the simulation in NVT ensemble at the target temperature
    integrator:set_temperature(args.temperature)
    observables.sampler:run(steps / 2 - steps / 10)

    -- log intermediate profiler results and reset accumulators
    utility.profiler:profile()

    -- disconnect NVT integrator from sampler and profiler
    integrator.disconnect()

    -- add velocity-Verlet integrator (NVE)
    integrator = mdsim.integrators.verlet({
        box = box, particle = particle, timestep = args.timestep
    })

    -- run remaining part of the simulation in NVE ensemble
    -- to prepare for the NVE production run
    observables.sampler:run(steps - steps / 2)
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "binary_mixture_equilibration_%Y%m%d_%H%M%S", help = "basename of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})

    parser:add_argument("particles", {type = "vector", dtype = "integer", default = {4000, 1000}, help = "number of particles"})
    parser:add_argument("density", {type = "number", default = 1.2, help = "particle number density"})
    parser:add_argument("ratios", {type = "vector", dtype = "number", action = function(args, key, value)
        if #value ~= 2 and #value ~= 3 then
            error(("box ratios has invalid dimension '%d'"):format(#value), 0)
        end
        args[key] = value
    end, default = {1, 1, 1}, help = "relative aspect ratios of simulation box"})
    parser:add_argument("masses", {type = "vector", dtype = "number", default = {1}, help = "particle masses"})
    parser:add_argument("initial-temperature", {type = "number", default = 1.5, help = "initial temperature"})
    parser:add_argument("temperature", {type = "number", default = 0.7, help = "target temperature"})
    parser:add_argument("rate", {type = "number", default = 2, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default = 1000, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.002, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
end
