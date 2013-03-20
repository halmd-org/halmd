#!/usr/bin/env halmd
--
-- Copyright © 2010-2013 Felix Höfling
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
local function liquid(args)
    -- open H5MD trajectory file for reading
    local file = readers.h5md({path = args.trajectory})

    -- construct a phase space reader and sample
    local reader, sample = observables.phase_space.reader({
        file = file, location = {"trajectory", "all"}, fields = {"position", "velocity", "species", "mass"}
    })
    -- read phase space sample at last step in file
    reader:read_at_step(-1)
    -- determine system parameters from phase space sample
    local nparticle = assert(sample.nparticle)
    local nspecies = assert(sample.nspecies)
    local dimension = assert(sample.dimension)

    -- read edge vectors of simulation domain from file
    local edges = mdsim.box.reader(file)
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({edges = edges})

    -- close H5MD trajectory file
    file:close()

    -- create system state
    local particle = mdsim.particle({box = box, particles = nparticle, species = nspecies})

    -- smoothly truncated Lennard-Jones potential
    local potential = mdsim.potentials.lennard_jones({particle = particle, cutoff = args.cutoff})
    -- smooth truncation
    local trunc
    if args.smoothing > 0 then
        trunc = mdsim.forces.trunc.local_r4({h = args.smoothing})
    end
    -- compute forces
    local force = mdsim.forces.pair_trunc({
        box = box, particle = particle, potential = potential, trunc = trunc
    })

    -- add velocity-Verlet integrator
    local integrator = mdsim.integrators.verlet({
        box = box
      , particle = particle
      , force = force
      , timestep = args.timestep
    })

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output)})
    -- write box specification to H5MD file
    box:writer(file)

    -- select all particles
    local particle_group = mdsim.particle_groups.all({particle = particle})

    -- sample phase space
    local phase_space = observables.phase_space({box = box, group = particle_group})
    -- set particle positions, velocities, species
    phase_space:set(sample)
    -- write trajectory of particle groups to H5MD file
    local interval = args.sampling.trajectory or args.steps
    if interval > 0 then
        phase_space:writer({file = file, fields = {"position", "velocity", "species", "mass"}, every = interval})
    end

    -- sample macroscopic state variables
    local msv = observables.thermodynamics({box = box, group = particle_group, force = force})
    local interval = args.sampling.state_vars
    if interval > 0 then
        msv:writer(file, {every = interval})
    end

    -- time correlation functions
    local interval = args.sampling.correlation
    if interval > 0 then
        -- setup blocking scheme
        local max_lag = args.steps * integrator.timestep / 10
        local blocking_scheme = dynamics.blocking_scheme({
            max_lag = max_lag
          , every = interval
          , size = 10
        })

        -- compute mean-square displacement
        local msd = dynamics.mean_square_displacement({phase_space = phase_space})
        blocking_scheme:correlation(msd, file)
        -- compute mean-quartic displacement
        local mqd = dynamics.mean_quartic_displacement({phase_space = phase_space})
        blocking_scheme:correlation(mqd, file)
        -- compute velocity autocorrelation function
        local vacf = dynamics.velocity_autocorrelation({phase_space = phase_space})
        blocking_scheme:correlation(vacf, file)

        -- compute interdiffusion coefficient
        local selfdiffusion = dynamics.correlation({
            -- acquire centre of mass
            acquire = function()
                return msv:r_cm()
            end
            -- correlate centre of mass at first and second point in time
          , correlate = function(first, second)
                local result = 0
                for i = 1, #first do
                    result = result + math.pow(second[i] - first[i], 2)
                end
                return result
            end
            -- H5MD group names
          , group = {particle_group.label, "selfdiffusion"}
            -- profiling description
          , desc = "selfdiffusion coefficient of A particles"
        })
        blocking_scheme:correlation(selfdiffusion, file)
    end

    -- sample initial state
    observables.sampler:sample()

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({
        steps = args.steps
      , first = 10
      , interval = 900
      , sample = 60
    })

    -- run simulation
    observables.sampler:run(args.steps)

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
    end, default = "lennard_jones_%Y%m%d_%H%M%S", help = "prefix of output files"})

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

    parser:add_argument("trajectory", {type = "string", required = true, action = function(args, key, value)
        if not readers.h5md.check(value) then
            error(("not an H5MD file: %s"):format(value), 0)
        end
        args[key] = value
    end, help = "trajectory file name"})

    parser:add_argument("cutoff", {type = "number", default = math.pow(2, 1 / 6), help = "potential cutoff radius"})
    parser:add_argument("smoothing", {type = "number", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("steps", {type = "integer", default = 10000, help = "number of simulation steps"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
    sampling:add_argument("correlation", {type = "integer", default = 100, help = "for correlation functions"})

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
liquid(args)
