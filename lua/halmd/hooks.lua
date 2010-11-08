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
local pairs = pairs
local ipairs = ipairs
local setmetatable = setmetatable

module("halmd.hooks")

local hooks = {}

local objects = setmetatable({}, { __mode = "k" }) -- table with weak keys

function register_hook(func)
    for object, module in pairs(objects) do
        func(object, module)
    end
    table.insert(hooks, func)
end

function register_object(object, module)
    if not objects[object] then
        for i, func in ipairs(hooks) do
            func(object, module)
        end
        objects[object] = module
    end
end
