#!/usr/bin/env halmd
--
-- Copyright © 2013      Nicolas Höft
-- Copyright © 2010-2015 Felix Höfling
-- Copyright © 2010-2012 Peter Colberg
--
-- This file is part of HALMD.
--
-- HALMD is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with this program.  If not, see <http://www.gnu.org/licenses/>.
--

local halmd = require("halmd")
local rescale_velocity = require("rescale_velocity")

-- grab modules
local log = halmd.io.log
local mdsim = halmd.mdsim
local observables = halmd.observables
local dynamics = halmd.observables.dynamics
local readers = halmd.io.readers
local writers = halmd.io.writers

--
-- Setup and run simulation
--
local function shear_viscosity(args)
    local nparticle = 10000   -- total number of particles
    local equi_steps = 5e5    -- steps used to equilibrate the system

    -- create simulation domain with periodic boundary conditions
    -- derive edge lengths from number of particles and density
    local length = math.pow(nparticle / args.density, 1 / 3)
    local box = mdsim.box({length = {length, length, length}})

    -- create system state
    local particle = mdsim.particle({dimension = 3, particles = nparticle})

    -- set initial particle positions
    local lattice = mdsim.positions.lattice({box = box, particle = particle})
    lattice:set()
    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle
      , temperature = args.temperature
    })
    boltzmann:set()

    -- Lennard-Jones potential
    local potential = mdsim.potentials.pair.lennard_jones({cutoff = args.cutoff})

    -- compute forces
    local force = mdsim.forces.pair_trunc({
        box = box, particle = particle, potential = potential
    })

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output)})

    -- select all particles
    local particle_group = mdsim.particle_groups.all({particle = particle})

    -- sample phase space
    local phase_space = observables.phase_space({box = box, group = particle_group})

    -- Sample macroscopic state variables.
    local msv = observables.thermodynamics({box = box, group = particle_group})
    local interval = args.sampling.state_vars
    if interval > 0 then
        msv:writer({file = file, every = interval})
    end

    local internal_energy = observables.utility.accumulator({
        aquire = msv.internal_energy
      , every = 500
      , start = math.floor(equi_steps / 2)
      , desc = "averaged internal energy"
      , aux_enable = {particle}
    })

    -- sample initial state
    observables.sampler:sample()

    -- add velocity-Verlet integrator with Boltzmann thermostat
    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box
      , particle = particle
      , timestep = args.timestep
      , temperature = args.temperature
      , rate = args.rate
    })

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({
        steps = equi_steps, first = 10, interval = 900, sample = 60
    })

    -- run equilibration
    observables.sampler:run(equi_steps)
    -- log profiler results
    halmd.utility.profiler:profile()

    integrator:disconnect()
    runtime:disconnect()

    --
    -- After equilibration, we can measure the correlation functions now
    --

    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- determine the mean energy from equilibration run
    -- and rescale velocities to match the mean energy of the equilibration
    rescale_velocity({msv = msv, internal_energy = internal_energy:mean()})
    internal_energy:disconnect()

    local temperature = observables.utility.accumulator({
        aquire = msv.temperature
      , every = 200
      , desc = "averaged temperature"
      , aux_enable = {particle}
    })

    temperature:writer({
        file = file
      , location = {"observables", msv.group.label, "averaged_temperature"}
      , every = math.floor(steps / 1000)
      , reset = true
    })

    -- monitor internal energy in NVE run in addition to the fields registered above
    interval = args.sampling.state_vars
    if interval > 0 then
        msv:writer({file = file, fields = {"internal_energy"}, every = interval})
    end

    -- replace thermostat integrator by NVE velocity-Verlet
    integrator = mdsim.integrators.verlet({
        box = box
      , particle = particle
      , timestep = args.timestep
    })

    -- setup blocking scheme
    local max_lag = steps * integrator.timestep / 10
    local interval = 100
    local blocking_scheme = dynamics.blocking_scheme({
        max_lag = max_lag
      , every = interval
      , size = 10
    })

    local helfand_moment = dynamics.helfand_moment({thermodynamics = msv, interval = 5})
    blocking_scheme:correlation({tcf = helfand_moment, file = file})

    local stress_tensor_autocorrelation = dynamics.stress_tensor_autocorrelation({thermodynamics = msv})
    blocking_scheme:correlation({tcf = stress_tensor_autocorrelation, file = file})

    runtime = observables.runtime_estimate({
        steps = steps, first = 10, interval = 900, sample = 60
    })

    -- run NVE simulation
    observables.sampler:run(steps)

    -- log profiler results
    halmd.utility.profiler:profile()
end

--
-- Parse command-line arguments.
--
local function parse_args()
    local parser = halmd.utility.program_options.argument_parser()

    parser:add_argument("output,o", {type = "string", action = function(args, key, value)
        -- substitute current time
        args[key] = os.date(value)
    end, default = "shear_viscosity_%Y%m%d", help = "prefix of output files"})
    -- _%Y%m%d_%H%M%S

    parser:add_argument("verbose,v", {type = "accumulate", action = function(args, key, value)
        local level = {
            -- console, file
            {"warning", "info" },
            {"info"   , "info" },
            {"debug"  , "debug"},
            {"trace"  , "trace"},
        }
        args[key] = level[value] or level[#level]
    end, default = 1, help = "increase logging verbosity"})

    parser:add_argument("density", {type = "number", default = 0.8442, help = "particle number density"})
    parser:add_argument("cutoff", {type = "float32", default = 2.5, help = "potential cutoff radius"})
    parser:add_argument("masses", {type = "vector", dtype = "number", default = {1}, help = "particle masses"})
    parser:add_argument("temperature", {type = "number", default = 0.722, help = "initial system temperature"})
    parser:add_argument("rate", {type = "number", default = 5, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default = 1e4, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})

    return parser:parse_args()
end

local args = parse_args()

-- log to console
halmd.io.log.open_console({severity = args.verbose[1]})
-- log to file
halmd.io.log.open_file(("%s.log"):format(args.output), {severity = args.verbose[2]})
-- log version
halmd.utility.version.prologue()

-- run simulation
shear_viscosity(args)
