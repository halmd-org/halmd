#!/usr/bin/env halmd
--
-- Copyright © 2024 Felix Höfling
-- Copyright © 2024 Jung Nguyen
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
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local dynamics = halmd.observables.dynamics
local readers = halmd.io.readers
local writers = halmd.io.writers
local utility = halmd.utility

-- search definition files in the top-level path relative to the simulation script
package.path = utility.abspath("../?.lua;") .. package.path
local definitions = {
    lennard_jones = require("definitions/lennard_jones")
  , slit_pore = require("definitions/slit_pore")
}

--
-- read particle data and box dimensions from input file
--
function read_sample(args)
    local input = utility.assert_kwarg(args, "input")
    local location = utility.assert_kwarg(args, "location")
    local fields = args.fields      -- optional, see phase_space.reader()
    local step = args.step or -1    -- by default, use the last step stored in the file

    -- open H5MD file for reading
    local file = readers.h5md({path = input})

    -- read edge vectors of simulation domain from file,
    -- only cuboid boxes are supported, so 'edges' must be a diagonal matrix
    local edges = mdsim.box.reader({file = file, location = location})
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({edges = edges})

    -- construct a phase space reader and sample
    local reader, sample = observables.phase_space.reader({
        file = file, location = location, fields = fields
    })
    -- read phase space sample at given step in file
    reader:read_at_step(step)

    -- determine system parameters from phase space sample
    local nparticle = assert(sample.nparticle)
    local dimension = assert(sample.dimension)
    if (assert(sample.nspecies) ~= 1) then
        log.error("Simulation script expects input file with particle data for a single species.")
    end

    -- close H5MD file
    file:close()

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle})

    -- use phase space "sampler" to set 'all' particle data (e.g., positions
    -- and velocities) by converting the sample obtained from the file
    local all_group = mdsim.particle_groups.all({particle = particle})
    observables.phase_space({box = box, group = all_group}):set(sample)

    return box, particle
end

--
-- setup and run simulation
--
function main(args)
    -- read particle data and box dimensions from input file
    local box, particle = read_sample({
        input = args.input, location = {"particles", "all"}, fields = {"position", "velocity"}
    })

    -- define Lennard-Jones pair potential (with parameters ε=1 and σ=1 for a single species)
    -- and register computation of pair forces
    definitions.lennard_jones.create_pair_force({
        box = box, particle = particle, cutoff = args.pair_cutoff, smoothing = args.smoothing
    })

    -- define Lennard-Jones wall potential to form a slit pore
    -- and register computation of this external force
    definitions.slit_pore.create_wall_force({
        box = box, particle = particle
      , pore_width = args.pore_width
      , epsilon = args.wall_epsilon, wetting = args.wetting
      , cutoff = args.wall_cutoff, smoothing = args.smoothing
    })

    -- add velocity-Verlet integrator
    local integrator = mdsim.integrators.verlet({
        box = box, particle = particle, timestep = args.timestep
    })

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- select all particles to be used by the observables
    local all_group = mdsim.particle_groups.all({particle = particle})

    -- write trajectory of particle groups to H5MD file
    local steps = math.ceil(args.time / args.timestep)
    local traj_interval = args.sampling.trajectory or steps
    if traj_interval > 0 then
        -- reuse instance of phase space sampler from above
        observables.phase_space({box = box, group = all_group})
           :writer({file = file, fields = {"position", "velocity"}, every = traj_interval})
    end

    -- sample macroscopic state variables.
    local msv = observables.thermodynamics({box = box, group = all_group})
    if args.sampling.state_vars > 0 then
        msv:writer({file = file,
                    fields = {
                        "potential_energy", "pressure", "temperature"  -- fluctuating quantities
                      , "internal_energy", "center_of_mass_velocity"   -- conserved quantities
                    },
                    every = args.sampling.state_vars})
    end

    -- set up wavevectors, compute density modes and static structure factor
    if args.sampling.structure > 0 or args.sampling.correlation > 0 then
        -- set up dense wavevector grid compatible with the periodic simulation box,
        -- filter wave vectors to point perpendicular to the wall surfaces
        wavevector = observables.utility.wavevector({box = box, wavenumber = {50}, dense = true, filter = {1,0,0}})

        -- compute density modes and output their time series
        local density_mode = observables.density_mode({group = all_group, wavevector = wavevector})
        if args.sampling.structure > 0 then
            density_mode:writer({file = file, every = args.sampling.structure})
        end
    end

    -- sample initial state
    observables.sampler:sample()

    -- estimate remaining runtime
    observables.runtime_estimate({steps = steps})

    -- run simulation
    observables.sampler:run(steps)
end


--
-- parse command-line arguments
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time
        , default ="slit_pore_planar_walls_production_D{pore_width:g}_rho{density:g}_T{temperature:.2f}_%Y%m%d_%H%M%S"
        , help = "basename of output files"
    })
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})

    parser:add_argument("input", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD input file"})

    parser:add_argument("pore-width", {type = "float32", default = 20, help = "width of the slit pore"})
    parser:add_argument("wetting", {type = "float32", default = 1, help = "wetting parameter of the wall potential"})
    parser:add_argument("wall-epsilon", {type = "float32", default = 1, help = "interaction strengths of the wall potential"})
    parser:add_argument("wall-cutoff", {type = "float32", default = 2.5, help = "cutoff distance of the wall potential"})
    parser:add_argument("pair-cutoff", {type = "float32", default = 2.5, help = "cutoff radius of pair potential"})
    parser:add_argument("smoothing", {type = "float32", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("time", {type = "number", default = 100, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
    sampling:add_argument("structure", {type = "integer", default = 1000, help = "for density modes, static structure factor"})
    sampling:add_argument("correlation", {type = "integer", default = 100, help = "for correlation functions"})
end
