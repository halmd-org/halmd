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

require("halmd.io.trajectory.readers")

-- grab environment
local readers = halmd.io.trajectory.readers
local po = halmd_wrapper.po
local error = error
local io = io
local pairs = pairs

module("halmd.io.trajectory.reader", halmd.modules.register)

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    desc:add("trajectory", po.string():notifier(function(value)

        -- check whether file exists and is readable
        local file, message = io.open(value)
        if not file then
            error(message, 0)
        end
        file:close()

        -- check for reader that can handle file format
        local reader
        for _, module in pairs(readers) do
            local format = module.format
            if format and format(value) then
                reader = module
                break
            end
        end
        if not reader then
            error(value .. ": unknown trajectory file format", 0)
        end

    end), "trajectory input file")

    desc:add("sample", po.int64():depends("trajectory"), "trajectory sample for initial state")
end
