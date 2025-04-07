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
local writers = halmd.io.writers
local utility = halmd.utility

-- search definition files in the top-level path relative to the simulation script
package.path = utility.abspath("../?.lua;") .. package.path
local definitions = {
    lennard_jones = require("definitions/lennard_jones")
  , slit_pore = require("definitions/slit_pore")
}

--
-- set up particle system and define interactions
--
function setup(args)
-- the following code is part of the documentation in doc/recipes/create_site_pore.rst.in,
-- line numbers must not be changed
-- define geometry
local pore_width = args.pore_width
local box_length = args.box_length
local length = {pore_width + 10, box_length, box_length}   -- add extra space "behind" the walls
local slab = {pore_width / length[1], 1, 1}
-- compute particle number from pore volume and mean density
local nparticle = math.floor(args.density * numeric.prod(slab) * numeric.prod(length))

-- create system state
local box = mdsim.box({length = length})
local particle = mdsim.particle({dimension = #length, particles = nparticle})

-- set initial particle positions
mdsim.positions.lattice({box = box, particle = particle, slab = slab}):set()

-- set initial particle velocities
mdsim.velocities.boltzmann({particle = particle, temperature = args.temperature}):set()

-- define Lennard-Jones pair potential (with parameters ε=1 and σ=1 for a single species)
-- and register computation of pair forces
definitions.lennard_jones.create_pair_force({
    box = box, particle = particle, cutoff = args.pair_cutoff, smoothing = args.smoothing
})

-- define Lennard-Jones wall potential to form a slit pore
-- and register computation of this external force
definitions.slit_pore.create_wall_force({
    box = box, particle = particle
  , pore_width = pore_width
  , epsilon = args.wall_epsilon, wetting = args.wetting
  , cutoff = args.wall_cutoff, smoothing = args.smoothing
})
-- end of usage in doc/recipes/create_slit_pore.rst.in

    return box, particle
end

--
-- setup and run simulation
--
function main(args)
    -- set up system and interactions
    local box, particle = setup(args)

    -- H5MD file writer for output
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- select all particles
    local all_group = mdsim.particle_groups.all({particle = particle})

    -- sample phase space
    local phase_space = observables.phase_space({box = box, group = all_group})

    -- write trajectory of particle groups to H5MD file
    local steps = math.ceil(args.time / args.timestep)
    local traj_interval = args.sampling.trajectory or steps
    if traj_interval > 0 then
        phase_space:writer({file = file, fields = {"position", "image", "velocity"}, every = traj_interval})
    end

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
    local runtime = observables.runtime_estimate({steps = steps})

    -- run simulation
    observables.sampler:run(steps)
end

--
-- parse command-line arguments
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time
        , default ="slit_pore_planar_walls_equilibration_D{pore_width:g}_rho{density:g}_T{temperature:.2f}_%Y%m%d_%H%M%S"
        , help = "basename of output files"
    })
    parser:add_argument("overwrite", {type = "boolean", default = true, help = "overwrite output file"})
    parser:add_argument("random-seed", {type = "integer", action = parser.action.random_seed
        , help = "seed for random number generator"
    })

    parser:add_argument("box-length", {type = "float32", default = 50, help = "box size parallel to the pore walls"})
    parser:add_argument("pore-width", {type = "float32", default = 20, help = "width of the slit pore"})
    parser:add_argument("wetting", {type = "float32", default = 1, help = "wetting parameter of the wall potential"})
    parser:add_argument("wall-epsilon", {type = "float32", default = 1, help = "interaction strengths of the wall potential"})
    parser:add_argument("wall-cutoff", {type = "float32", default = 2.5, help = "cutoff distance of the wall potential"})
    parser:add_argument("pair-cutoff", {type = "float32", default = 2.5, help = "cutoff radius of pair potential"})
    parser:add_argument("smoothing", {type = "float32", default = 0.005, help = "cutoff smoothing parameter"})

    parser:add_argument("density", {type = "number", default = 0.8, help = "particle number density"})
    parser:add_argument("temperature", {type = "float32", default = 1.5, help = "initial system temperature"})
    parser:add_argument("rate", {type = "number", default = 4, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default = 100, help = "integration time"})
    parser:add_argument("timestep", {type = "float32", default = 0.001, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", default = 1000, help = "save snapshot for every $trajectory time step"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "calculate state variables for every $state-vars time step"})
end
