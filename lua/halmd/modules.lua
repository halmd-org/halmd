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
local hooks = require("halmd.hooks")
local po = halmd_wrapper.po
local ipairs = ipairs
local pairs = pairs
local rawget = rawget
local rawset = rawset
local setmetatable = setmetatable
local string = string
local table = table

--
-- This module keeps a registry of all HALMD Lua modules.
--
-- Modules register themselves by supplying the function register
-- as an additional parameter to the module() command.
--
module("halmd.modules")

-- ordered list of registered HALMD modules
local modules = {}

-- parsed module options
local vm = {}

--
-- register a module
--
-- @param module Lua module
--
function register(module)
    table.insert(modules, module)

    --
    -- construct C++ module
    --
    -- @param module Lua module
    -- @param args script parameter table
    -- @returns C++ module object
    --
    local new = function(module, args)
        --
        -- get module parameter
        --
        -- @param self parameter table
        -- @param key parameter name
        -- @returns parameter value, or nil if not found
        --
        local param = function(self, key)
            -- command-line option value
            local value = rawget(vm, key)

            -- script parameter
            if not value and args then
                value = rawget(args, key)
            end

            -- cache parameter
            rawset(self, key, value)

            return value
        end

        local args = setmetatable({}, { __index = param })
        local object = module.new(args)
        if object then
            hooks.register_object(object, module)
        end
        return object
    end

    local defaults = {
        -- option/parameter namespace
        namespace = module._NAME:match("[^.]+$")
    }

    setmetatable(module, { __call = new, __index = defaults })
end

--
-- query options of registered modules
--
-- @param desc Boost.Program_options options description
--
function options(desc)
    local list = {}
    local map = setmetatable({}, {
        __index = function(self, key)
            local desc = po.options_description(key)
            table.insert(list, desc)
            rawset(self, key, desc)
            return desc
        end
    })
    for i, module in ipairs(modules) do
        if module.options then
            module.options(map[module.namespace])
        end
    end
    for i, v in ipairs(list) do
        desc:add(v)
    end
end

--
-- Set parsed command line options
--
-- @param args Boost.Program_options variables_map
--
-- This function is called by halmd::script.
--
function parsed(args)
    for arg, value in pairs(args) do
        local key = string.gsub(arg, "-", "_") -- e.g. options.power_law_index
        vm[key] = value
    end
end
