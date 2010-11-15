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
    local box_ratios = args.box_ratios or {}
    local box_length = args.box_length -- optional

    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local particle = assert(core.particle)

    local box = box_wrapper[dimension]
    if box_length then
        -- complete missing values by repeating the last entry
        for i = #box_length + 1, dimension do
            box_length[i] = box_length[#box_length]
        end
        return box(particle, box_length)
    else
        -- fill up missing values with 1
        for i = #box_ratios + 1, dimension do
            box_ratios[i] = 1 -- `neutral' aspect ratio
        end
        return box(particle, density, box_ratios)
    end
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    desc:add("density,d", po.float(), "particle density")
    desc:add("box-length,L", po.array_float():conflicts("density"), "edge lengths of simulation box")
    desc:add("box-ratios", po.array_float():conflicts("box-length"), "aspect ratios of simulation box (specify relative edge lengths)")
end
