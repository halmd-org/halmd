--
-- Copyright © 2014 Felix Höfling
-- Copyright © 2011 Peter Colberg
--
-- This file is part of HALMD.
--
-- HALMD is free software: you can redistribute it and/or modify
-- it under the terms of the GNU Lesser General Public License as
-- published by the Free Software Foundation, either version 3 of
-- the License, or (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU Lesser General Public License for more details.
--
-- You should have received a copy of the GNU Lesser General
-- Public License along with this program.  If not, see
-- <http://www.gnu.org/licenses/>.
--

local halmd = halmd

function test()
    local box = halmd.mdsim.box{length = {20, 20, 20}}

    local excluded = halmd.mdsim.positions.excluded_volume({box = box, cell_length = 10})
    excluded:exclude_sphere({3, 3, 3}, 1)
    assert(not excluded:place_sphere({3, 3, 3}, 1))
    assert(excluded:place_sphere({4, 4, 4}, 1))
    assert(excluded:place_sphere({4, 4, 4}, 2 * math.sqrt(3) - 1))
    assert(not excluded:place_sphere({4, 4, 4}, 2.00001 * math.sqrt(3) - 1))
    excluded:exclude_sphere({4, 4, 4}, 1)
    assert(not excluded:place_sphere({4, 4, 4}, 1))
    excluded:exclude_sphere({3, 3, 3}, 10)
end

-- the following function body is part of the documentation in
-- lua/mdsim/positions/excluded_volume.lua.in and shall not be indented
function example()
local random = math.random
math.randomseed(os.time())

local edge = 50
local box = halmd.mdsim.box{length = {edge, edge, edge}}
local excluded = halmd.mdsim.positions.excluded_volume{box = box, cell_length = 1}

local obstacles = {}
local diameter = 1
local repeats = 50
for i = 1, 1000 do
    for j = 1, repeats do
        local r = {edge * random(), edge * random(), edge * random()}
        if excluded:place_sphere(r, diameter) then
            obstacles[i] = r
            excluded:exclude_sphere(r, diameter)
            break
        end
    end
    if not obstacles[i] then
        error(("cannot place obstacle %d after %d repeats"):format(i, repeats))
    end
end

local particle = halmd.mdsim.particle{dimension = box.dimension, particles = #obstacles}
particle:set_position(obstacles)
end

function run()
    test()
    example()
end
