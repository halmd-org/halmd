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
    host = halmd_wrapper.host.random
  , read_integer = halmd_wrapper.random.random.read_integer
}
local po = halmd_wrapper.po
local assert = assert

module("halmd.host.random", halmd.modules.register)

function options(desc)
    desc:add("random-seed", po.uint(), "random number generator integer seed")
    desc:add("random-file", po.string(), "read random seed from file")
end

local random -- singleton instance

--
-- construct random module
--
function new(args)
    local file = args.random_file or "/dev/random" -- default value
    local seed = args.random_seed -- optional

    if not random then
        if not seed then
            seed = random_wrapper.read_integer(file)
        end
        random = random_wrapper.host.gfsr4(seed)
    end
    return random
end
