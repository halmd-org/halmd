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
local modules = require("halmd.modules")
local setmetatable = setmetatable
local pairs = pairs
local string = string

module("halmd.options", modules.register)

setmetatable(_M, {
    --
    -- Set parsed command line options
    --
    -- @param self Module halmd.options
    -- @param parsed Boost.Program_options variables_map
    --
    -- This function is called by halmd::script.
    --
    __call = function(self, parsed)
        for k, v in pairs(parsed) do
            option = string.gsub(k, "-", "_") -- e.g. options.power_law_index
            self[option] = v:value()
        end
    end
})
