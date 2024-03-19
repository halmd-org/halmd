#!/usr/bin/env halmd
--
-- Copyright © 2010-2023 Felix Höfling
-- Copyright © 2010-2012 Peter Colberg
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
local SetForces = require("potential") -- dynamic loading


--
-- setup and run simulation
--
function main(args)

    -- H5MD file writer for output
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite =
    args.overwrite})

	-- create simulation box with PBCs centered at (0,0,0)
    local length = {50, 30, 30}
    local box = mdsim.box({length = length})
    local slab = {args.pore_width / length[1],1,1}

    -- create system state
    local particle = mdsim.particle({dimension = #length, particles = math.floor(args.density * numeric.prod(slab) * numeric.prod(length))})
    -- set initial particle positions
    mdsim.positions.lattice({box = box, particle = particle, slab = slab}):set()

    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({particle = particle, temperature = args.temperature}):set()

	-- select particles
    local all_group = mdsim.particle_groups.all({particle = particle})

	-- set fluid-fluid forces
    local ff_forces = SetForces.create_pair_force(box, particle, args)

	-- set fluid-wall forces
    local wall_forces = SetForces.create_wall_force(box, particle, args)

    -- sample phase space
	local phase_space = observables.phase_space({box = box, group = all_group})

    -- write trajectory of particle groups to H5MD file
    local steps = math.ceil(args.time / args.timestep)
    local traj_interval = args.sampling.trajectory or steps
    if traj_interval > 0 then
        phase_space:writer({file = file, fields = {"position", "velocity"}, every = traj_interval})
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
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time, default ="slit_1d_2walls_equilibration", help = "basename of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = true, help = "overwrite output file"})
    parser:add_argument("random-seed", {type = "integer", action = parser.action.random_seed, help = "seed for random number generator"})
    parser:add_argument("pore-width", {type = "float32", default = 20, help = "pore width"})
    parser:add_argument("density", {type = "number", default = 0.8, help = "particle number density"})
    parser:add_argument("cutoff", {type = "float32", default = 2.5, help = "potential cutoff radius"})
    parser:add_argument("epsilon", {type = "float32", default = 1, help = "interaction strengths in wall potential"})
    parser:add_argument("wetting", {type = "float32", default = 1, help = "wetting parameters in wall potential"})
    parser:add_argument("smoothing", {type = "float32", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("temperature", {type = "float32", default = 1.5, help = "initial system temperature"})
    parser:add_argument("rate", {type = "number", default = 4, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default = 100, help = "integration time"})
    parser:add_argument("timestep", {type = "float32", default = 0.001, help = "integration time step"})
    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", default = 1000, help = "save snapshot for every $trajectory time step"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "calculate state variables for every $state-vars time step"})
end