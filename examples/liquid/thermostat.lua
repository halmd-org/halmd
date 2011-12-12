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
-- Simple liquid simulation with thermostat, equilibration and production phase
--
local thermostat = {_NAME = "thermostat"} -- FIXME implement halmd.module("name")

halmd.modules.register(thermostat)

function thermostat.new(args)
    -- load the device module to log (optional) GPU properties
    device{} -- singleton
    -- open (optional) H5MD file and read simulation parameters
    local reader = readers.trajectory{group = "liquid"}

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
    mdsim.box{
        particles = {particle}
    }
    -- add integrator
    local thermostat = mdsim.integrators.verlet_nvt_andersen{
        particle = particle
    }
    -- add force
    local force = mdsim.force{
        particle = particle
    }
    -- set initial particle positions (optionally from reader)
    mdsim.position{
        reader = reader
      , particle = particle
    }
    -- set initial particle velocities (optionally from reader)
    mdsim.velocity{
        reader = reader
      , particle = particle
    }

    -- Construct sampler.
    local sampler = observables.sampler{}

    -- Write trajectory to H5MD file.
    writers.trajectory{particle = particle, group = "liquid"}
    -- Sample macroscopic state variables.
    observables.thermodynamics{particle = particle, force = force}

    -- yield sampler.setup slot from Lua to C++ to setup simulation box
    coroutine.yield(sampler:setup())

    -- thermostat system in NVT ensemble
    coroutine.yield(sampler:run(assert(args.thermostat)))

    -- equilibrate system in NVE ensemble
    thermostat:disconnect()
    mdsim.integrators.verlet{particle = particle}
    coroutine.yield(sampler:run(assert(args.equilibrate)))

    -- Sample static structure factor.
    observables.ssf{particle = particle}

    -- yield sampler.run slot from Lua to C++ to run simulation
    coroutine.yield(sampler:run())
end

function thermostat.options(desc, globals)
    globals:add("particles", po.uint_array():default({1000}), "number of particles")
    globals:add("masses", po.uint_array():default({1}), "particle masses")
    globals:add("dimension", po.uint():default(3):notifier(function(value)
        if value ~= 2 and value ~= 3 then
            error(("invalid dimension '%d'"):format(value), 0)
        end
    end), "dimension of positional coordinates")
    globals:add("thermostat", po.uint():default(200), "thermostat steps")
    globals:add("equilibrate", po.uint():default(200), "equilibation steps")
end

return thermostat
