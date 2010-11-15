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
local random_wrapper = {
    read_integer = halmd_wrapper.random.random.read_integer
}
if halmd_wrapper.gpu then
    random_wrapper.gpu = halmd_wrapper.gpu.random
end
local device = require("halmd.device")
local po = halmd_wrapper.po
local assert = assert

module("halmd.gpu.random", halmd.modules.register)

local random -- singleton instance

--
-- construct random module
--
function new(args)
    local file = args.random_file or "/dev/random" -- default value
    local seed = args.random_seed -- optional
    local blocks = args.random_blocks or 32 -- default value
    local threads = args.random_threads or 256 -- default value FIXME DEVICE_SCALE

    if not random then
        if not seed then
            seed = random_wrapper.read_integer(file)
        end
        random = random_wrapper.gpu.rand48(device(), seed, blocks, threads)
    end
    return random
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    if random_wrapper.gpu then
        desc:add("random-seed", po.uint(), "random number generator integer seed")
        desc:add("random-file", po.string(), "read random seed from file")
        desc:add("random-blocks", po.uint(), "number of CUDA blocks")
        desc:add("random-threads", po.uint(), "number of CUDA threads per block")
    end
end
