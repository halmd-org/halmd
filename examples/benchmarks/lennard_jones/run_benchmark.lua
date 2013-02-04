#!/usr/bin/env halmd
--
-- Copyright © 2012 Felix Höfling
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
local readers = halmd.io.readers
local writers = halmd.io.writers
local utility = halmd.utility

--
-- Setup and run simulation
--
local function lennard_jones(args)
    local timestep = 0.002   -- integration timestep
    local steps = 10000      -- number of integration steps
    local cutoff = 3.0       -- potential cutoff
    local count = args.count -- number of repetitions

    -- open H5MD trajectory file and read phase space data for A and B particles
    local file = readers.h5md({path = args.trajectory})

    local trajectory = file.root:open_group("trajectory")
    -- construct a phase space reader and sample
    local reader, sample = observables.phase_space.reader(file, {group = "all"})
    -- read phase space sample at last step in file
    log.info("number of particles: %d", sample.nparticle)
    reader:read_at_step(-1)

    -- read edge vectors of simulation domain from file
    local edges = mdsim.box.reader(file)
    -- create simulation box
    local box = mdsim.box({edges = edges})

    -- close H5MD trajectory file
    file:close()

    -- open H5MD file writer
    local writer = writers.h5md({path = ("%s.h5"):format(args.output)})
    -- write box specification to H5MD file
    box:writer(writer)

    -- create system state
    local particle = mdsim.particle({box = box, particles = sample.nparticle, species = 1})

    -- setup particles from trajectory sample
    local all_group = mdsim.particle_groups.all({particle = particle, label = "all"})
    -- construct phase space instance
    local phase_space = observables.phase_space({box = box, group = all_group})
    -- set particle positions, velocities, species
    phase_space:set(sample)
    -- write phase space data of group to H5MD file, but only first and last step
    phase_space:writer(writer)

    -- define interaction of Kob-Andersen mixture using truncated Lennard-Jones potential
    local potential = mdsim.potentials.lennard_jones({particle = particle, cutoff = cutoff})
    -- compute forces
    local force = mdsim.forces.pair_trunc({
        box = box, particle = particle, potential = potential, -- neighbour_skin = 0.7
    })

    -- define velocity-Verlet integrator
    local integrator = mdsim.integrators.verlet({
        box = box, particle = particle, force = force, timestep = timestep
    })

    -- estimate remaining runtime
    observables.runtime_estimate({steps = count * steps, first = 10, interval = 900, sample = 60})

    -- sample initial state
    observables.sampler:sample()

    -- run simulation several times and output profiling data
    for i = 1, count do
        observables.sampler:run(steps)
        utility.profiler:profile()
    end
end

--
-- Parse command-line arguments.
--
local function parse_args()
    local parser = utility.program_options.argument_parser()

    parser:add_argument("output,o", {type = "string", action = function(args, key, value)
        -- substitute current time
        args[key] = os.date(value)
    end, default = "lennard_jones_equilibration_%Y%m%d_%H%M%S", help = "prefix of output files"})

    parser:add_argument("trajectory", {type = "string", required = true, action = function(args, key, value)
        if not readers.h5md.check(value) then
            error(("not an H5MD file: %s"):format(value), 0)
        end
        args[key] = value
    end, help = "trajectory file name"})

    parser:add_argument("count", {type = "number", default = 5, help = "number of repetitions"})

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

    return parser:parse_args()
end

local args = parse_args()

-- log to console
log.open_console({severity = args.verbose[1]})
-- log to file
log.open_file(("%s.log"):format(args.output), {severity = args.verbose[2]})
-- log version
utility.version.prologue()

-- run simulation
lennard_jones(args)
