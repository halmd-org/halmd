--
-- Copyright © 2010-2011  Peter Colberg and Felix Höfling
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
      , label = (function()
          -- generate labels "A", "B", "C", … according to number of species
          local label = {}
          for i = 1, #args.particles do
              label[i] = string.char(string.byte("A") + i - 1)
          end
          return label
      end)()
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

    -- Write trajectory to H5MD file.
    -- FIXME support filtering of particles by species
    writers.trajectory{particle = particle, group = "A"}
    -- Sample macroscopic state variables.
    observables.thermodynamics{particle = particle, force = force}
    -- Sample static structure factor.
    observables.ssf{particle = particle}
    -- compute mean-square displacement
    observables.dynamics.correlation{particle = particle, correlation = "mean_square_displacement"}
    -- compute mean-quartic displacement
    observables.dynamics.correlation{particle = particle, correlation = "mean_quartic_displacement"}
    -- compute velocity autocorrelation
    observables.dynamics.correlation{particle = particle, correlation = "velocity_autocorrelation"}

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
end

return liquid
