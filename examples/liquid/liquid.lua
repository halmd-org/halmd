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
local log = halmd.io.log
local mdsim = halmd.mdsim
local observables = halmd.observables
local readers = halmd.io.readers
local writers = halmd.io.writers

--
-- Setup and run simulation
--
local function liquid(args)
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

--
-- Parse command-line arguments.
--
local function parse_args()
    local parser = halmd.utility.program_options.argument_parser()

    parser:add_argument("output,o", {type = "string", action = function(args, key, value)
        -- substitute current time
        args[key] = os.date(value)
    end, default = "liquid_%Y%m%d_%H%M%S", help = "prefix of output files"})

    parser:add_argument("verbose,v", {type = "accumulate", action = function(args, key, value)
        local level = {
            -- console, file
            {"warning", "info" },
            {"info"   , "info" },
            {"debug"  , "debug"},
            {"trace"  , "trace"},
        }
        args[key] = level[value] or level[#level]
    end, default = 1, help = "increase logging verbosity"})

    parser:add_argument("disable-gpu", {type = "boolean", help = "disable GPU acceleration"})
    parser:add_argument("devices", {type = "vector", dtype = "integer", help = "CUDA device(s)"})

    return parser:parse_args()
end

local args = parse_args()

-- log to console
halmd.io.log.open_console({severity = args.verbose[1]})
-- log to file
halmd.io.log.open_file(("%s.log"):format(args.output), {severity = args.verbose[2]})
-- initialise or disable GPU
halmd.utility.device({disable_gpu = args.disable_gpu, devices = args.devices})

-- run simulation
liquid(args)
