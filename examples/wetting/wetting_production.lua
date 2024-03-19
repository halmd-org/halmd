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
local dynamics = halmd.observables.dynamics
local readers = halmd.io.readers
local writers = halmd.io.writers
local utility = halmd.utility
local SetForces = require("potential") -- dynamic loading


--
-- setup and run simulation
--
function main(args)

    -- open H5MD file for reading
    local file = readers.h5md({path = args.input})

    -- construct a phase space reader and sample
    local reader, sample = observables.phase_space.reader({
        file = file, location = {"particles", "all"}, fields = {"position", "velocity"}
    })

    -- read phase space sample at last step in file
    reader:read_at_step(-1)
    -- determine system parameters from phase space sample
    local nparticle = assert(sample.nparticle)
    local dimension = assert(sample.dimension)
    if (assert(sample.nspecies) ~= 1) then
        log.error("Simulation script expects input file with particle data for a single species.")
    end

    -- read edge vectors of simulation domain from file
    local edges = mdsim.box.reader({file = file, location = {"particles", "all"}})
    local length = {edges[1][1], edges[2][2], edges[3][3]}
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({edges = edges})

    -- close H5MD file
    file:close()

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle})

    -- select particles
    local all_group = mdsim.particle_groups.all({particle = particle})

	-- set fluid-fluid forces
    local ff_forces = SetForces.create_pair_force(box, particle, args)

    -- set fluid-wall forces
    local wall_forces = SetForces.create_wall_force(box, particle, args)

    -- add velocity-Verlet integrator
    local integrator = mdsim.integrators.verlet({
        box = box, particle = particle, timestep = args.timestep
    })

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- use phase space "sampler" to set particle positions and velocities
    -- by converting the sample obtained from the file,
    -- the sampler will be reused below for observables and file output
    local phase_space = observables.phase_space({box = box, group = all_group})
    phase_space:set(sample)

    -- write trajectory of particle groups to H5MD file
    local steps = math.ceil(args.time / args.timestep)
    local traj_interval = args.sampling.trajectory or steps
    if traj_interval > 0 then
        -- reuse instance of phase space sampler from above
        phase_space:writer({file = file, fields = {"position", "velocity"}, every = traj_interval})
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
        -- set up wavevector grid compatible with the periodic simulation box
        wavevector = observables.utility.wavevector({box = box, wavenumber = {50}, dense =true, filter = {1,0,0}})

        -- compute density modes and output their time series
        local density_mode = observables.density_mode({group = all_group, wavevector = wavevector})
        if args.sampling.structure > 0 then
            density_mode:writer({file = file, every = args.sampling.structure})
        end

        -- compute static structure factor from density modes
        local ssf = observables.ssf({density_mode = density_mode, norm = nparticle})
        -- output averages over a certain number of configurations each
        ssf:writer({file = file, every = args.sampling.structure})
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
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "slit_1d_2walls_production", help = "basename of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = true, help = "overwrite output file"})
    parser:add_argument("pore-width", {type = "number", default = 20, help = "pore_width"})
    parser:add_argument("input", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD input file"})
    parser:add_argument("rescale-to-energy", {type = "number", help = "rescale velocities to match given internal energy"})
    parser:add_argument("cutoff", {type = "float32", default = math.pow(2, 1 / 6), help = "potential cutoff radius"})
    parser:add_argument("epsilon", {type = "float32", default = 1, help = "interaction strengths in wall potential"})
    parser:add_argument("wetting", {type = "float32", default = 1, help = "wetting parameters in wall potential"})
    parser:add_argument("smoothing", {type = "number", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("time", {type = "number", default = 100, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})
    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
    sampling:add_argument("structure", {type = "integer", default = 1000, help = "for density modes, static structure factor"})
    sampling:add_argument("correlation", {type = "integer", default = 100, help = "for correlation functions"})
end