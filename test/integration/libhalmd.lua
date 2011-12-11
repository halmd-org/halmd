--
-- Copyright © 2011  Peter Colberg
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

--
-- This test ensures that libhalmd.so fulfills the requirement of a binary
-- Lua extension module, i.e. is loadable from the upstream Lua interpreter.
--

local halmd = require("halmd")

-- load the device module to log (optional) GPU properties
halmd.device{} -- singleton
-- create system state
local particle = halmd.mdsim.particle{particles = {1000}, masses = {1}, label = {"A"}, dimension = 3}
-- create simulation box with particles
halmd.mdsim.box{particles = {particle}, dimension = 3, length = {10, 10, 10}}
-- add integrator
halmd.mdsim.integrator{particle = particle}
-- add force
halmd.mdsim.force{particle = particle}
-- set initial particle positions
halmd.mdsim.position{particle = particle}
-- set initial particle velocities
halmd.mdsim.velocity{particle = particle}

-- Construct sampler.
local sampler = halmd.observables.sampler{steps = 100}

-- setup simulation box
sampler:setup()()
-- run simulation
sampler:run()()
