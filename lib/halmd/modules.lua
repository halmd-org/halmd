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

-- grab environment
local table = table
local ipairs = ipairs

--
-- This module keeps a registry of all HALMD Lua modules.
--
-- Modules register themselves by supplying the modules.register
-- function as an additional parameter to the module() command.
--
module("halmd.modules")

local registry = {} -- ordered list of registered HALMD modules

--
-- register a module
--
function register(module)
    table.insert(registry, module)
end

--
-- query options of registered modules
--
-- @param desc Boost.Program_options options description
--
function options(desc)
    for i, module in ipairs(registry) do
        if module.options then
            module.options(desc)
        end
    end
end
