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

local objects = setmetatable({}, { __mode = "k" }) -- table with weak keys

local object_hooks = {}

--
-- Register object hook, and apply to registered objects.
--
-- An object hook is applied once to each C++ object.
--
-- @param hook function, receives C++ object as argument
--
function register_object_hook(hook)
    for object, _ in pairs(objects) do
        hook(object)
    end
    table.insert(object_hooks, hook)
end

local module_hooks = {}

--
-- Register module hook, and apply to registered objects.
--
-- A module hook is applied once to each pair of Lua module and C++ object.
--
-- @param hook function, receives Lua module and C++ object as arguments
--
function register_module_hook(hook)
    for _, modules in pairs(objects) do
        for module, object in pairs(modules) do
            hook(module, object)
        end
    end
    table.insert(module_hooks, hook)
end

--
-- Register C++ object, and apply registered object and module hooks.
--
-- @param object C++ object
-- @param module Lua module
--
function register_object(object, module)
    if not objects[object] then
        for _, hook in ipairs(object_hooks) do
            hook(object)
        end
        objects[object] = setmetatable({}, { __mode = "v" }) -- table with weak values
    end
    local modules = objects[object]
    if not modules[module] then
        for _, hook in ipairs(module_hooks) do
            hook(module, object)
        end
        modules[module] = object
    end
end
