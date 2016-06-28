#!/usr/bin/env halmd
--
-- Copyright © 2012 Felix Höfling
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

local halmd = require("halmd")

-- grab modules
local mdsim = halmd.mdsim
local observables = halmd.observables
local writers = halmd.io.writers
local utility = halmd.utility

--
-- Setup and run simulation
--
local function lennard_jones(args)
    local nparticle = 64000   -- total number of particles
    local density = 0.4       -- number density
    local temperature = 1.2   -- heat bath temperature
    local cutoff = 3.0        -- potential cutoff

    local timestep = 0.01     -- integration timestep
    local steps = 10000       -- number of integration steps

    -- derive edge length of three-dimensional box from particle number and density
    local length = math.pow(nparticle / density, 1 / 3)

    -- create simulation box
    local box = mdsim.box({length = {length, length, length}})

    -- create system state with 1 particle species
    local particle = mdsim.particle({dimension = 3, particles = nparticle})

    -- set initial particle positions
    mdsim.positions.lattice({box = box, particle = particle}):set()

    -- set initial particle velocities
    mdsim.velocities.boltzmann({particle = particle, temperature = temperature}):set()

    -- define interaction using truncated Lennard-Jones potential
    local potential = mdsim.potentials.pair.lennard_jones({cutoff = cutoff})
    -- compute forces
    local force = mdsim.forces.pair_trunc({
        box = box, particle = particle, potential = potential, -- neighbour_skin = 0.7
    })

    -- define velocity-Verlet integrator with Andersen thermostat
    local integrator = mdsim.integrators.verlet_nvt_andersen({
        box = box, particle = particle
      , timestep = timestep, temperature = temperature, rate = 0.1
    })

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output)})

    -- sample group of all particles
    local all_group = mdsim.particle_groups.all({particle = particle})

    -- connect H5MD writer, output only first and last step
    observables.phase_space({box = box, group = all_group}):writer({
        file = file
      , fields = {"position", "velocity", "species", "mass"}
    })

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
local function parse_args()
    local parser = utility.program_options.argument_parser()

    parser:add_argument("output,o", {type = "string", action = function(args, key, value)
        -- substitute current time
        args[key] = os.date(value)
    end, default = "lennard_jones_benchmark_configuration_%Y%m%d_%H%M%S", help = "prefix of output files"})

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
halmd.io.log.open_console({severity = args.verbose[1]})
-- log to file
halmd.io.log.open_file(("%s.log"):format(args.output), {severity = args.verbose[2]})
-- log version
utility.version.prologue()

-- run simulation
lennard_jones(args)
