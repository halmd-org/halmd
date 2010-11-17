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
local assert = assert
local ipairs = ipairs
local pairs = pairs
local rawget = rawget
local rawset = rawset
local setmetatable = setmetatable
local string = string
local table = table
local type = type

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
            local value
            if module.options then
                local vm = vm[module.namespace]
                if type(vm) == "table" then
                    value = rawget(vm, key)
                end
            end
            if not value then
                value = rawget(vm, key)
            end

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
-- @param parser options_parser instance
--
function options(parser)
    for _, module in ipairs(modules) do
        if module.options then
            local desc = po.options_description()
            local globals = po.options_description()
            module.options(desc, globals)
            parser:add(desc, module.namespace)
            parser:add(globals)
        end
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
    for k, v in pairs(args) do
        vm[k] = v
    end
end
