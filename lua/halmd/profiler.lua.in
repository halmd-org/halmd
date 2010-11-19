--
-- Copyright © 2010  Peter Colberg
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

require("halmd.modules")

-- grab environment
local profiler_wrapper = halmd_wrapper.utility.profiler
local profiling_writers = require("halmd.io.profiling.writers")
local hooks = require("halmd.hooks")
local pairs = pairs
local table = table

module("halmd.profiler", halmd.modules.register)

--
-- construct profiler module
--
function new()
    local writers = profiling_writers()
    local profiler = profiler_wrapper(writers)

    hooks.register_object_hook(function(profilable)
        if profilable.register_runtimes then
            profilable:register_runtimes(profiler)
        end
    end)

    return profiler
end
