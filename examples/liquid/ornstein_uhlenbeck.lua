#!/usr/bin/env halmd
--
-- Copyright © 2010-2014 Felix Höfling
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
local writers = halmd.io.writers
local utility = halmd.utility

--
-- Setup and run simulation
--
function main(args)
    -- total number of particles from sum of particles per species
    local nparticle = args.particles

    -- derive edge lengths from number of particles, density and edge ratios
    local volume = nparticle / args.density
    local ratios = args.ratios
    local dimension = #ratios
    local unit_edge = math.pow(volume / numeric.prod(ratios), 1 / dimension)
    local length = {}

    if dimension == 2 then
        length = {100, 100} 
    else
        length = {100, 100, 100} 
    end
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle})

    local position_table = particle.data["position"]
    for i = 1, #position_table do
        if dimension == 2 then
            position_table[i] = {1,1}
        else 
            position_table[i] = {1,1,1}
        end
        particle.data["position"] = position_table
    end

    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle
      , temperature = args.temperature
    }):set()

    -- define harmonic potential
    local potential_args = {
        stiffness = args.stiffness,
        offset = args.offset
    }
    local potential = mdsim.potentials.external.harmonic(potential_args)

    -- register computation of external forces
    mdsim.forces.external({
        box = box, particle = particle, potential = potential
    })

    -- convert integration time to number of steps
    local init_steps = math.ceil(args.init_time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- select all particles
    local all_group = mdsim.particle_groups.all({particle = particle})

    -- sample phase space
    local phase_space = observables.phase_space({box = box, group = all_group})
    -- write trajectory of particle groups to H5MD file
    local interval = args.sampling.trajectory or init_steps
    if interval > 0 then
        phase_space:writer({file = file, fields = {"position"}, every = interval})
    end

    -- Sample macroscopic state variables.
    local interval = args.sampling.state_vars
    if interval > 0 then
        observables.thermodynamics({box = box, group = all_group})
          :writer({file = file, fields = {"potential_energy", "center_of_mass"}, every = interval})
    end

    -- sample initial state
    observables.sampler:sample()

    -- add Brownian integrator
    local integrator = mdsim.integrators.brownian({
        box = box, particle = particle, timestep = args.timestep, temperature = args.temperature, diff_const = args.diff_const
    })

    -- run simulation
    observables.sampler:run(init_steps)


    -- -- simulation after initialization

    -- -- convert integration time to number of steps
    local steps = math.ceil((args.time - args.init_time) / args.timestep)

    -- time correlation functions
    local interval = args.sampling.correlation
    if interval > 0 then
        -- setup blocking scheme
        local max_lag = steps * integrator.timestep / 10
        local blocking_scheme = dynamics.blocking_scheme({max_lag = max_lag, every = interval, size = 10})

        -- compute mean-square displacement
        local msd = dynamics.mean_square_displacement({phase_space = phase_space})
        blocking_scheme:correlation({tcf = msd, file = file})
        -- compute mean-quartic displacement
        local mqd = dynamics.mean_quartic_displacement({phase_space = phase_space})
        blocking_scheme:correlation({tcf = mqd, file = file})
    end

    -- run simulation
    observables.sampler:run(steps)

end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "ornstein_uhlenbeck_D{diff_const:.2f}_T{temperature:.2f}_%Y%m%d_%H%M%S", help = "basename of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})

    parser:add_argument("random-seed", {type = "integer", action = parser.action.random_seed,
        help = "seed for random number generator"})

    parser:add_argument("particles", {type = "integer", default = 5000, help = "number of particles"})
    parser:add_argument("density", {type = "number", default = 0.75, help = "particle number density"})
    parser:add_argument("ratios", {type = "vector", dtype = "number", action = function(args, key, value)
        if #value ~= 2 and #value ~= 3 then
            error(("box ratios has invalid dimension '%d'"):format(#value), 0)
        end
        args[key] = value
    end, default = {1, 1, 1}, help = "relative aspect ratios of simulation box"})

    parser:add_argument("masses", {type = "vector", dtype = "number", default = {1}, help = "particle masses"})
    parser:add_argument("stiffness", {type = "vector", dtype = "number", default = {2}, help = "stiffness matrix/constant for harmonic potential"})
    parser:add_argument("offset", {type = "vector", dtype = "number", default = {0,0,0}, help = "offset for harmonic potential"})
    parser:add_argument("diff_const", {type = "number", default = 0.3, help = "diffusion constant"})
    parser:add_argument("temperature", {type = "number", default = 3, help = "initial system temperature"})
    parser:add_argument("rate", {type = "number", default = 2, help = "heat bath collision rate"})
    parser:add_argument("init_time", {type = "number", default = 100, help = "initialization time"})
    parser:add_argument("time", {type = "number", default = 200, help = "total integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.01, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
    sampling:add_argument("correlation", {type = "integer", default = 10, help = "for correlation functions"})
end
