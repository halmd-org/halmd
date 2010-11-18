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
local hooks = require("halmd.hooks")
local assert = assert

module("halmd.parameter", halmd.modules.register)

--
-- register HDF5 writer as parameter writer
--
-- @param writer HDF5 writer object
--
function register_writer(writer)
    hooks.register_module_hook(function(module, object)
        if module.write_parameters then
            -- HDF5 file with write access
            local file = writer:file()
            -- open parameter/namespace group
            local namespace = assert(module.namespace)
            local param = file:open_group("param")
            local group = file:open_group(namespace)

            module.write_parameters(object, group, param)
        end
    end)
end
