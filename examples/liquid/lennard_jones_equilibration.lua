#!/usr/bin/env halmd
--
-- Copyright © 2010-2014 Felix Höfling
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
local numeric = halmd.numeric
local observables = halmd.observables
local writers = halmd.io.writers

--
-- Setup and run simulation
--
local function liquid(args)
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
    -- set particle species, with continuous range of tags per species
    local species = {}
    for s = 0, nspecies - 1 do
        local nparticle = assert(args.particles[s + 1])
        for i = 1, nparticle do table.insert(species, s) end
    end
    particle:set_species(species)

    local origin = {}
    for i = 1, dimension do
        origin[i] = 0
        length[i] = length[i]/2
    end
    local geometry = mdsim.geometries.cuboid({origin = origin, length = length})
    local region = mdsim.region({particle = particle, label = "upper quadrant", geometry = geometry, box = box})

    -- set initial particle positions
    local lattice = mdsim.positions.lattice({box = box, particle = particle})
    lattice:set()
    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle
      , temperature = args.temperature
    })
    boltzmann:set()

    -- smoothly truncated Lennard-Jones potential
    local potential = mdsim.potentials.pair.lennard_jones({cutoff = args.cutoff, species = particle.nspecies})
    -- smooth truncation
    local trunc = nil
    if args.smoothing > 0 then
        trunc = mdsim.forces.trunc.local_r4({h = args.smoothing})
    end
    -- compute forces
    local force = mdsim.forces.pair_trunc({
        box = box, particle = particle, potential = potential, trunc = trunc
    })

    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output)})

    -- select all particles
    local particle_group = mdsim.particle_groups.all({particle = particle})
    local group_included = mdsim.particle_groups.from_region({particle = particle, region = region, selection = "included", label = "included"})

    local msv_local = observables.thermodynamics({box = box, group = group_included})
    observables.sampler:on_sample(
        function()
            halmd.io.log.info("particles in region: %d, temp: %f", group_included.size, msv_local:temperature())
        end
      , 200
      , 0
    )

    -- sample phase space
    local phase_space = observables.phase_space({box = box, group = particle_group})
    -- write trajectory of particle groups to H5MD file
    local interval = args.sampling.trajectory or steps
    if interval > 0 then
        phase_space:writer({
            file = file, fields = {"position", "velocity", "species", "mass"}, every = interval
        })
    end

    -- Sample macroscopic state variables.
    local msv
    local interval = args.sampling.state_vars
    if interval > 0 then
        msv = observables.thermodynamics({box = box, group = particle_group})
        msv:writer({file = file, every = interval})
    end

    local accumulator = observables.utility.accumulator({
         aquire = msv.internal_energy
       , every = 10
       , desc = "Averaged internal energy"
       , aux_enable = {particle}
     })
     accumulator:writer({
         file = file
       , location = {"observables", "averaged_internal_energy"}
       , every = 200
       , reset = true
     })

    -- sample initial state
    observables.sampler:sample()

    -- add velocity-Verlet integrator with Boltzmann distribution
    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box
      , particle = particle
      , timestep = args.timestep
      , temperature = args.temperature
      , rate = args.rate
    })

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({
        steps = steps, first = 10, interval = 900, sample = 60
    })

    -- run simulation
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
    end, default = "lennard_jones_equilibration_%Y%m%d_%H%M%S", help = "prefix of output files"})

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
    parser:add_argument("random-seed", {type = "integer", help = "seed for random number generator"})

    parser:add_argument("particles", {type = "vector", dtype = "integer", default = {10000}, help = "number of particles"})
    parser:add_argument("density", {type = "number", default = 0.75, help = "particle number density"})
    parser:add_argument("ratios", {type = "vector", dtype = "number", action = function(args, key, value)
        if #value ~= 2 and #value ~= 3 then
            error(("box ratios has invalid dimension '%d'"):format(#value), 0)
        end
        args[key] = value
    end, default = {1, 1, 1}, help = "relative aspect ratios of simulation box"})
    parser:add_argument("cutoff", {type = "float32", default = math.pow(2, 1 / 6), help = "potential cutoff radius"})
    parser:add_argument("smoothing", {type = "number", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("masses", {type = "vector", dtype = "number", default = {1}, help = "particle masses"})
    parser:add_argument("temperature", {type = "number", default = 1.5, help = "initial system temperature"})
    parser:add_argument("rate", {type = "number", default = 0.1, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default = 100, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.005, help = "integration time step"})

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

-- seed the random number generator
if args.random_seed then
    halmd.random.generator({seed = args.random_seed})
end

-- run simulation
liquid(args)
