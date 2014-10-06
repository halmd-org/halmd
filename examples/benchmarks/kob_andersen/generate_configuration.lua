#!/usr/bin/env halmd
--
-- Copyright © 2011-2014 Felix Höfling
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
local mdsim = halmd.mdsim
local observables = halmd.observables
local random = halmd.random
local writers = halmd.io.writers
local utility = halmd.utility

--
-- Setup and run simulation
--
local function kob_andersen(args)
    local nparticle = 256000  -- total number of particles
    local concentration = 0.8 -- concentration of A particles
    local density = 1.2       -- number density
    local temperature = 0.7   -- heat bath temperature

    local timestep = 0.005    -- integration timestep
    local steps = 20000       -- number of integration steps

    -- number of particles in groups A and B
    local ngroup = { math.floor(nparticle * concentration) }
    ngroup[2] = nparticle - ngroup[1]

    -- derive edge length of three-dimensional box from particle number and density
    local length = math.pow(nparticle / density, 1 / 3)

    -- create simulation box
    local box = mdsim.box({length = {length, length, length}})

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle, species = 2})

    -- set particle species, with continuous range of tags per species:
    -- construct array with particle species: (0, 0, … 0, 1, 1, … 1)
    local species = {}
    for i = 1, ngroup[1] do table.insert(species, 0) end
    for i = 1, ngroup[2] do table.insert(species, 1) end
    particle:set_species(species)

    -- set initial particle positions, randomise the particle species
    mdsim.positions.lattice({box = box, particle = particle}):set()
    -- randomly shuffle the positions
    particle:set_position(random.generator({memory = "host"}):shuffle(particle:get_position()))

    -- set initial particle velocities
    mdsim.velocities.boltzmann({particle = particle, temperature = temperature}):set()

    -- define interaction of Kob-Andersen mixture using truncated Lennard-Jones potential
    local potential = mdsim.potentials.pair.lennard_jones({
        epsilon = {{1, 1.5}, {1.5, 0.5}} -- ((AA, AB), (BA, BB))
      , sigma = {{1, 0.8}, {0.8, 0.88}} -- ((AA, AB), (BA, BB))
      , cutoff = 2.5
    })
    -- compute forces
    local force = mdsim.forces.pair_trunc({box = box, particle = particle, potential = potential})

    -- define velocity-Verlet integrator with Andersen thermostat
    local integrator = mdsim.integrators.verlet_nvt_andersen({
        box = box, particle = particle,
        timestep = timestep, temperature = temperature, rate = 0.1
    })

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output)})

    -- sample each particle group separately and write to H5MD file
    local offset = 0
    for s = 0, 1 do
        local group = mdsim.particle_groups.from_range({
            particle = particle
          , range = {offset + 1, offset + ngroup[s + 1]}
          , label = string.char(string.byte('A') + s)
        })
        offset = offset + ngroup[s + 1]

        -- connect H5MD writer, output only first and last step
        observables.phase_space({box = box, group = group}):writer({
            file = file
          , fields = {"position", "velocity", "species", "mass"}
        })
    end

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
    end, default = "kob_andersen_benchmark_configuration_%Y%m%d_%H%M%S", help = "prefix of output files"})

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
kob_andersen(args)
