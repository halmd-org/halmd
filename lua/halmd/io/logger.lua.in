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
local logger_wrapper = halmd_wrapper.io.logger
local setmetatable = setmetatable
local rawset = rawset
local rawget = rawget

module("halmd.io.logger")

--
--
-- @param level Logging severity level
--
local function log(self, level)
    local log = logger_wrapper[level]
    if log then
        rawset(self, level, log)
    else
        rawset(self, level, function() end)
    end
    return rawget(self, level)
end

setmetatable(_M, { __index = log })
