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

--
-- Setup and run simulation
--
local function script()
    -- load the device module to log (optional) GPU properties
    device() -- singleton
    -- open (optional) H5MD file and read simulation parameters
    local reader = readers.trajectory()

    -- create simulation box with particles
    mdsim.box()
    -- add integrator
    mdsim.integrator()
    -- add force
    local force = mdsim.force()
    -- set initial particle positions (optionally from reader)
    mdsim.position{reader = reader}
    -- set initial particle velocities (optionally from reader)
    mdsim.velocity{reader = reader}

    -- Construct sampler.
    local sampler = observables.sampler()

    -- Write trajectory to H5MD file.
    writers.trajectory{}
    -- Sample macroscopic state variables.
    observables.thermodynamics{force = force}
    -- Sample static structure factor.
    observables.ssf{}
    -- compute mean-square displacement
    observables.dynamics.correlation{correlation = "mean_square_displacement"}
    -- compute mean-quartic displacement
    observables.dynamics.correlation{correlation = "mean_quartic_displacement"}
    -- compute velocity autocorrelation
    observables.dynamics.correlation{correlation = "velocity_autocorrelation"}

    -- yield sampler.setup slot from Lua to C++ to setup simulation box
    coroutine.yield(sampler:setup())

    -- yield sampler.run slot from Lua to C++ to run simulation
    coroutine.yield(sampler:run())
end

return script
