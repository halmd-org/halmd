#!/usr/bin/env halmd
--
-- Copyright © 2010-2012  Peter Colberg and Felix Höfling
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

--
-- Setup and run simulation
--
local function liquid(args)
    -- total number of particles from sum of particles per species
    local nspecies = #args.particles
    local nparticle = 0
    for i = 1, nspecies do
        nparticle = nparticle + args.particles[i]
    end
    -- derive edge lengths from number of particles, density and edge ratios
    local volume = nparticle / args.density
    local dimension = #args.ratios
    local det = 1
    for i = 1, #args.ratios do
        det = det * args.ratios[i]
    end
    local length = {}
    for i = 1, #args.ratios do
        length[i] = args.ratios[i] * math.pow(volume / det, 1 / dimension)
    end
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})

    -- label particles A, B, …

    -- create system state
    local particle = mdsim.particle({box = box, particles = nparticle, species = nspecies})
    -- add velocity-Verlet integrator
    local integrator = mdsim.integrators.verlet({box = box, particle = particle, timestep = args.timestep})
    -- pair potential
    local potential = mdsim.potentials.lennard_jones({particle = particle})
    -- add force
    local force = mdsim.forces.pair_trunc{box = box, particle = particle, potential = potential}
    -- set initial particle positions
    local lattice = mdsim.positions.lattice{box = box, particle = particle}
    lattice:set()
    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann{particle = particle}
    boltzmann:set()

    -- H5MD file writer
    local writer = writers.h5md({path = ("%s.trj"):format(args.output)})
    -- write box specification to H5MD file
    box:writer(writer)

    -- Set particle species, with continuous range of tags per species
    local species = {}
    for s = 1, nspecies do
        local nparticle = assert(args.particles[s])
        for i = 1, nparticle do table.insert(species, s) end
    end
    particle:set_species(species)

    local particle_group = mdsim.particle_groups.all{
        particle = particle -- FIXME , species = species
    }
    local phase_space = observables.phase_space({box = box, group = particle_group})
    -- write trajectory of particle groups to H5MD file
    phase_space:writer(writer, {every = args.sampling.trajectory})

    -- H5MD file writer
    local writer = writers.h5md({path = ("%s.obs"):format(args.output)})
    -- Sample macroscopic state variables.
    local msv = observables.thermodynamics({box = box, particle = particle})
    msv:writer(writer, {every = args.sampling.state_vars})

    -- Sample static structure factors, construct density modes before.
    -- FIXME local density_mode = observables.density_mode{
    -- FIXME     phase_space = phase_space, max_wavevector = 15
    -- FIXME }
    -- FIXME observables.ssf{density_mode = density_mode, every = args.sampling.structure}

    -- setup blocking scheme for correlation functions
    local max_lag = args.steps * integrator.timestep / 10
    local blocking_scheme = observables.dynamics.blocking_scheme({max_lag = max_lag, every = 100, size = 10})

    -- compute mean-square displacement
    local msd = observables.dynamics.mean_square_displacement({phase_space = phase_space})
    blocking_scheme:correlation(msd, writer)
    -- compute mean-quartic displacement
    local mqd = observables.dynamics.mean_quartic_displacement({phase_space = phase_space})
    blocking_scheme:correlation(mqd, writer)
    -- compute velocity autocorrelation function
    local vacf = observables.dynamics.velocity_autocorrelation({phase_space = phase_space})
    blocking_scheme:correlation(vacf, writer)
    -- compute intermediate scattering function from density modes different than those used for ssf computation
    -- FIXME density_mode = observables.density_mode{
    -- FIXME     phase_space = phase_space, max_wavevector = 12, decimation = 2
    -- FIXME }
    -- FIXME observables.dynamics.correlation{sampler = density_mode, correlation = "intermediate_scattering_function"}

    -- setup simulation box
    observables.sampler:setup()

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({steps = args.steps, first = 10, interval = 900, sample = 60})

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
    end, default = "liquid_%Y%m%d_%H%M%S", help = "prefix of output files"})

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

    parser:add_argument("particles", {type = "vector", dtype = "integer", default = {1000}, help = "number of particles"})
    parser:add_argument("density", {type = "number", default = 0.75, help = "particle number density"})
    parser:add_argument("ratios", {type = "vector", dtype = "number", action = function(args, key, value)
        if #value ~= 2 and #value ~= 3 then
            error(("box ratios has invalid dimension '%d'"):format(#value), 0)
        end
        args[key] = value
    end, default = {1, 1, 1}, help = "relative aspect ratios of simulation box"})
    parser:add_argument("masses", {type = "vector", dtype = "number", default = {1}, help = "particle masses"})

    parser:add_argument("ensemble", {type = "string", choices = {
        nve = "Constant NVE",
        nvt = "Constant NVT",
    }, default = "nve", help = "statistical ensemble"})

    parser:add_argument("steps", {type = "integer", default = 10000, help = "number of simulation steps"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals"})
    sampling:add_argument("trajectory", {type = "integer", default = 1000, help = "sampling interval for trajectory"})
    sampling:add_argument("structure", {type = "integer", default = 1000, help = "sampling interval for structural properties"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "sampling interval for state variables"})

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
