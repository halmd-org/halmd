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
local core_wrapper = {
    [2] = halmd_wrapper.mdsim.core_2_
  , [3] = halmd_wrapper.mdsim.core_3_
}
local h5 = halmd_wrapper.h5
local po = halmd_wrapper.po
local assert = assert

module("halmd.mdsim.core", halmd.modules.register)

local core -- singleton instance

--
-- construct core module
--
function new(args)
    local dimension = args.dimension or 3 -- default value
    if not core then
        core = core_wrapper[dimension]()
        assert(core.dimension == dimension)
    end
    return core
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    desc:add("dimension", po.uint(), "dimension of positional coordinates")
end

--
-- write module parameters to HDF5 group
--
-- @param core module instance
-- @param group HDF5 group
--
function write_parameters(core, group)
    group:write_attribute("dimension", h5.uint(), core.dimension)
end
