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
local device = halmd.device
local mdsim = halmd.mdsim
local observables = halmd.observables
local readers = halmd.io.readers
local writers = halmd.io.writers
-- grab C++ library
local po = libhalmd.po

--
-- Simple liquid simulation
--
local liquid = {_NAME = "liquid"} -- FIXME implement halmd.module("name")

halmd.modules.register(liquid)

function liquid.new(args)
    -- load the device module to log (optional) GPU properties
    device{} -- singleton
    -- FIXME support reading multiple species groups into single particle
    local reader = readers.trajectory{group = "A"}

    -- label particles A, B, …

    -- create system state
    local particle = mdsim.particle{
        particles = assert(args.particles)
      , masses = assert(args.masses)
      , dimension = assert(args.dimension)
      , label = "all" -- FIXME make optional
    }
    -- create simulation box
    mdsim.box{particles = {particle}}
    -- add integrator
    mdsim.integrator{particle = particle}
    -- add force
    local force = mdsim.force{particle = particle}
    -- set initial particle positions (optionally from reader)
    mdsim.position{reader = reader, particle = particle}
    -- set initial particle velocities (optionally from reader)
    mdsim.velocity{reader = reader, particle = particle}

    -- Construct sampler.
    local sampler = observables.sampler{}

    -- Construct particle groups and phase space samplers by species (species are numbered 0, 1, 2, ...)
    local species = {} for i = 1, #args.particles do species[i] = i - 1 end -- FIXME avoid explicit for-loop!?
    local particle_group = observables.samples.particle_group{
        particle = particle, species = species
    }
    local phase_space = observables.phase_space{particle_group = particle_group}

    -- Sample macroscopic state variables.
    observables.thermodynamics{particle_group = { particle }, force = { force }, every = args.state_vars or 10 }
    observables.thermodynamics{particle_group = { particle_group[1] }, force = { force }, every = args.state_vars or 10 }

    -- Write trajectory to H5MD file.
    writers.trajectory{particle_group = particle_group, every = args.trajectory or 10}

    -- Sample static structure factors, construct density modes before.
    local density_mode = observables.density_mode{
        phase_space = phase_space, max_wavevector = 15
    }
    observables.ssf{density_mode = density_mode, every = args.structure or 10}

    -- compute mean-square displacement
    observables.dynamics.correlation{sampler = phase_space, correlation = "mean_square_displacement"}
    -- compute mean-quartic displacement
    observables.dynamics.correlation{sampler = phase_space, correlation = "mean_quartic_displacement"}
    -- compute velocity autocorrelation function
    observables.dynamics.correlation{sampler = phase_space, correlation = "velocity_autocorrelation"}
    -- compute intermediate scattering function from density modes different than those used for ssf computation
    density_mode = observables.density_mode{
        phase_space = phase_space, max_wavevector = 12, decimation = 2
    }
--     observables.dynamics.correlation{sampler = density_mode, correlation = "intermediate_scattering_function"}

    -- yield sampler.setup slot from Lua to C++ to setup simulation box
    coroutine.yield(sampler:setup())

    -- yield sampler.run slot from Lua to C++ to run simulation
    coroutine.yield(sampler:run())
end

function liquid.options(desc, globals)
    globals:add("particles", po.uint_array():default({1000}), "number of particles")
    globals:add("masses", po.uint_array():default({1}), "particle masses")
    globals:add("dimension", po.uint():default(3):notifier(function(value)
        if value ~= 2 and value ~= 3 then
            error(("invalid dimension '%d'"):format(value), 0)
        end
    end), "dimension of positional coordinates")
--    globals:add("trajectory", po.uint64():default(0), "sampling interval for trajectory") -- FIXME boost::any_cast
--    globals:add("structure", po.uint64():default(0), "sampling interval for structural properties") -- FIXME boost::any_cast
--    globals:add("state-vars", po.uint64():default(0), "sampling interval for state variables") -- FIXME boost::any_cast
end

-- FIXME explicitly load modules to collect deprecated module options
require("halmd.option")
require("halmd.log")
require("halmd.device")
require("halmd.io.readers.trajectory")
require("halmd.io.writers.trajectory")
require("halmd.mdsim.box")
require("halmd.mdsim.force")
require("halmd.mdsim.integrator")
require("halmd.mdsim.particle")
require("halmd.mdsim.position")
require("halmd.mdsim.velocity")
require("halmd.observables.dynamics.correlation")
require("halmd.observables.sampler")
require("halmd.observables.ssf")
require("halmd.observables.thermodynamics")

return liquid
