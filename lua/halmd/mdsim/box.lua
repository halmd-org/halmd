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
local box_wrapper = {
    [2] = halmd_wrapper.mdsim.box_2_
  , [3] = halmd_wrapper.mdsim.box_3_
}
local mdsim = {
  core = require("halmd.mdsim.core")
}
local h5 = halmd_wrapper.h5
local po = halmd_wrapper.po
local assert = assert

module("halmd.mdsim.box", halmd.modules.register)

--
-- construct box module instance
--
-- @param args parameter table
-- @returns box module instance
--
function new(args)
    local density = args.density or 0.75 -- default value
    local ratios = args.ratios or {}
    local length = args.length -- optional

    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local particle = assert(core.particle)

    local box = box_wrapper[dimension]
    if length then
        -- complete missing values by repeating the last entry
        for i = #length + 1, dimension do
            length[i] = length[#length]
        end
        return box(particle, length)
    else
        -- fill up missing values with 1
        for i = #ratios + 1, dimension do
            ratios[i] = 1 -- `neutral' aspect ratio
        end
        return box(particle, density, ratios)
    end
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    desc:add("density,d", po.float(), "particle density")
    desc:add("length,L", po.float_array():conflicts("density"), "edge lengths of simulation box")
    desc:add("ratios", po.float_array():conflicts("length"), "aspect ratios of simulation box (specify relative edge lengths)")
end

--
-- write module parameters to HDF5 group
--
-- @param box module instance
-- @param group HDF5 group
--
function write_parameters(box, group)
    group:write_attribute("length", h5.float_array(), box.length)
    group:write_attribute("density", h5.float(), box.density)
end
