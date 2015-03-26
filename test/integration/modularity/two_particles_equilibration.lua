#!/usr/bin/env halmd
--
-- Copyright © 2013-2014 Felix Höfling
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
local random = halmd.random
local writers = halmd.io.writers
local utility = halmd.utility

--
-- Setup and equilibrate mixture using two distinct instances of the particle module
--
local function setup(args)

    local dimension = args.dimension      -- dimension of space
    local density = args.density          -- number density
    local np = args.particles             -- number of particles per species
    assert(#np == 2)
    local temp0 = 2 * args.temperature    -- initial temperature

    -- derive edge lengths from number of particles and density
    local length = { math.pow((np[1] + np[2]) / density, 1 / dimension) }
    for i = 2, dimension do
        length[i] = length[1]
    end
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})

    -- create system state for all particles first
    local particle = mdsim.particle({dimension = dimension, particles = np[1] + np[2], species = 2})
    -- set particle species, with continuous range of tags per species
    local species = {}
    local groups = {}
    local offset = 0
    for s = 0, 1 do
        local n = np[s + 1]
        for i = 1, n do
            table.insert(species, s)
        end
        local label = string.char(string.byte("A") + s)
        groups[label] = mdsim.particle_groups.from_range({
            particle = particle
          , range = {offset + 1, offset + n}
          , label = label
        })
        offset = offset + n
    end
    particle:set_species(species)
    -- set initial particle positions
    mdsim.positions.lattice({box = box, particle = particle}):set()
    -- randomly shuffle the positions
    particle:set_position(random.generator({memory = "host"}):shuffle(particle:get_position()))
    -- set initial particle velocities
    mdsim.velocities.boltzmann({particle = particle, temperature = temp0}):set()

    -- split global particle instance in two particle instances, one per species
    local particle_ = {
        A = groups["A"]:to_particle({label = "A"})
      , B = groups["B"]:to_particle({label = "B"})
    }
    particle.disconnect() -- needed?
    particle = particle_; particle_ = nil

    -- truncated Lennard-Jones potential
    -- FIXME move cutoff to pair_trunc
    local potential = mdsim.potentials.pair.lennard_jones({
        epsilon = {
            {1  , 1.5} -- AA, AB
          , {1.5, 0.5} -- BA, BB
        }
      , sigma = {
            {1  , 0.8 } -- AA, AB
          , {0.8, 0.88} -- BA, BB
        }
      , cutoff = 2.5
    })

    -- create binning modules explicitly and therefore only once for each particle instance
    local binning = {
        A = mdsim.binning({
            box = box
          , particle = particle["A"]
          , r_cut = potential.r_cut
        })
      , B = mdsim.binning({
            box = box
          , particle = particle["B"]
          , r_cut = potential.r_cut
        })
    }
    -- define interaction forces with smoothly truncated potential
    local force = {}
    local trunc = mdsim.forces.trunc.local_r4({h = 0.005})
    for label1, p1 in pairs(particle) do
        for label2, p2 in pairs(particle) do
            local neighbour = mdsim.neighbour({
                box = box
              , particle = { p1, p2 }
              , r_cut = potential.r_cut
              , binning = { binning[label1], binning[label2] }
            })
            force[label1 .. label2] = mdsim.forces.pair_trunc({
                box = box
              , particle = { p1, p2 }
              , potential = potential, trunc = trunc
              , neighbour = neighbour
            })
        end
    end

    return box, particle, args
end

--
-- equilibrate system of several instances of the particle module
--
local function equilibrate(box, particle, args)

    local temp = args.temperature         -- target temperature
    local temp0 = 2 * temp                -- initial temperature
    local collision_rate = 0.1            -- "collision" rate with heat bath
    local timestep = args.timestep        -- integration timestep
    local steps = math.ceil(args.time / timestep) -- number of integration steps

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output)})

    -- sample macroscopic state variables for all particles.
    -- FIXME provide user-defined thermodynamic variables and
    -- compute the sum from all particle instances by hand
--     if args.sampling.state_vars > 0 then
--         local msv = observables.thermodynamics({
--             box = box
--           , group = mdsim.particle_groups.all({particle = particle, global = false})
--         })
--         msv:writer({file = file, every = args.sampling.state_vars})
--     end

    -- sample each particle instance separately
    for label, p in utility.sorted(particle) do
        local group = mdsim.particle_groups.all({particle = p, global = false})
        -- write phase space trajectory to H5MD file
        observables.phase_space({box = box, group = group})
          : writer({
                file = file
              , fields = {"position", "velocity", "species", "mass"}
              , every = args.sampling.trajectory or steps
            })

        -- thermodynamic variables
        observables.thermodynamics({box = box, group = group})
          : writer({file = file, every = args.sampling.state_vars})
    end

    -- sample initial state
    observables.sampler:sample()

    -- add velocity-Verlet integrators with Boltzmann thermostat (NVT)
    local integrator = {}
    for k,v in pairs(particle) do
        integrator[k] = mdsim.integrators.verlet_nvt_boltzmann({
            box = box, particle = v
          , timestep = timestep, temperature = temp0, rate = collision_rate
        })
    end

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({
        steps = steps, first = 10, interval = 900, sample = 60
    })

    -- run first 10% of the simulation in NVT ensemble at elevated temperature
    -- in order to melt the initial fcc crystal
    observables.sampler:run(math.floor(steps / 10))

    -- switch to actual target temperature
    for k,v in pairs(integrator) do v:set_temperature(temp) end

    -- run remaining first half of the simulation in NVT ensemble at the target temperature
    observables.sampler:run(math.floor(steps / 2) - math.floor(steps / 10))

    -- log intermediate profiler results and reset accumulators
    halmd.utility.profiler:profile()

    -- disconnect NVT integrator from sampler and profiler
    -- and replace by velocity-Verlet integrator (NVE)
    for k,v in pairs(integrator) do
        v:disconnect()
        integrator[k] = mdsim.integrators.verlet({
            box = box, particle = particle[k], timestep = timestep
        })
    end

    -- run remaining part of the simulation in NVE ensemble
    -- to prepare for the NVE production run
    observables.sampler:run(steps - math.floor(steps / 2))

    -- log profiler results
    halmd.utility.profiler:profile()
end

--
-- Parse command-line arguments.
--
local function parse_args()
    local parser = halmd.utility.program_options.argument_parser()

    parser:add_argument("output,o",
        {type = "string", default = "two_particles_equilibration", help = "prefix of output files"})

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

    parser:add_argument("particles", {type = "vector", dtype = "integer", default = {4000, 1000}, help = "number of particles"})
    parser:add_argument("density", {type = "number", default = 1.2, help = "particle number density"})
    parser:add_argument("dimension", {type = "number", default = 3, help = "dimension of space"})
    parser:add_argument("temperature", {type = "number", default = 0.7, help = "target temperature"})
    parser:add_argument("time", {type = "number", default = 1000, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.002, help = "integration time step"})

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

-- set up system and perform equilibration run
equilibrate(setup(args))
