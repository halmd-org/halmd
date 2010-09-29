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
local args = require("halmd.options")
local assert = assert

module("halmd.gpu.random", halmd.modules.register)

local random -- singleton instance

--
-- construct random module
--
function new()
    local file = assert(args.random_file)
    local seed = args.random_seed -- optional
    local blocks = assert(args.random_blocks)
    local threads = assert(args.random_threads)

    if not random then
        if not seed then
            seed = random_wrapper.read_integer(file)
        end
        random = random_wrapper.gpu.rand48(device(), seed, blocks, threads)
    end
    return random
end

if random_wrapper.gpu then
    options = random_wrapper.gpu.options
end
