--
-- Copyright © 2010  Peter Colberg and Felix Höfling
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
local lattice_wrapper = {
    host = {
        [2] = halmd_wrapper.mdsim.host.positions.lattice_2_
      , [3] = halmd_wrapper.mdsim.host.positions.lattice_3_
    }
}
if halmd_wrapper.mdsim.gpu then
    lattice_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.positions.lattice_2_
      , [3] = halmd_wrapper.mdsim.gpu.positions.lattice_3_
    }
end
local mdsim = {
    core = require("halmd.mdsim.core")
}
local random = require("halmd.random")
local device = require("halmd.device")
local h5 = halmd_wrapper.h5
local po = halmd_wrapper.po
local assert = assert

module("halmd.mdsim.positions.lattice", halmd.modules.register)

--
-- construct lattice module
--
function new(args)
    local slab = args.slab or {} -- optional

    -- dependency injection
    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local particle = assert(core.particle)
    local box = assert(core.box)
    local random = assert(random())

    -- fill up missing values with 1
    for i = #slab + 1, dimension do
        slab[i] = 1 -- full box size
    end

    local lattice
    if device() then
        lattice = lattice_wrapper.gpu[dimension]
    else
        lattice = lattice_wrapper.host[dimension]
    end
    return lattice(particle, box, random, slab)
end

--
-- returns module description
--
function name()
    return "Face-centered cubic lattice"
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    desc:add("slab", po.float_array(), "slab extents as a fraction of simulation box")
end

--
-- write module parameters to HDF5 group
--
-- @param instance of lattice module
-- @param group HDF5 group
--
function write_parameters(lattice, group)
    group:write_attribute("slab", h5.float_array(), lattice.slab)
end
