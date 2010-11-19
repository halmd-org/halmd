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
local statevars_writers = {
    hdf5 = require("halmd.io.statevars.writers.hdf5")
}
local hooks = require("halmd.hooks")
local pairs = pairs
local table = table

module("halmd.io.statevars.writers", halmd.modules.register)

--
-- construct statevars module
--
function new()
    local writer = statevars_writers.hdf5()

    hooks.register_object_hook(function(observable)
        if observable.register_observables then
            observable:register_observables(writer)
        end
    end)

    return writer
end
