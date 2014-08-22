--
-- Copyright © 2014-2015 Nicolas Höft
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

local halmd = require("halmd")
local mdsim = halmd.mdsim
local observables = halmd.observables
local writers = halmd.io.writers

halmd.io.log.open_console()

local function setup(args)
    local dimension = args.dimension      -- dimension of space
    local np = args.particles             -- number of particles

    local length = {}
    for d = 1, dimension do
        length[d] = 10
    end
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})

    -- create system state for all particles first
    local particle = mdsim.particle({dimension = dimension, particles = np, species = 1})

    -- set initial particle positions
    mdsim.positions.lattice({box = box, particle = particle}):set()

    if particle.memory == "gpu" then
      -- sort particles *once*, as they do not move
      local sort = mdsim.sorts.hilbert({box = box, particle = particle})
      sort:order()
    end

    local origin = {}
    for i = 1, dimension do
        origin[i] = 0
        length[i] = length[i]/2
    end
    local geometry = mdsim.geometries.cuboid({origin = origin, length = length})
    local region = {}
    region["included"] = mdsim.region({particle = particle, label = "upper quadrant (included)", geometry = geometry, selection = "included", box = box})
    region["excluded"] = mdsim.region({particle = particle, label = "upper quadrant (excluded)", geometry = geometry, selection = "excluded", box = box})

    return box, particle, region, {origin = origin, length = length},  args
end

local function test(box, particle, region, cuboid, args)
    -- construct included/excluded particle groups
    local group_included = mdsim.particle_groups.from_region({particle = particle, region = region["included"], label = "included"})
    local group_excluded = mdsim.particle_groups.from_region({particle = particle, region = region["excluded"], label = "excluded"})
    -- check if the total number of particles is correct
    assert(group_excluded.size + group_included.size == args.particles)

    local dimension = particle.dimension

    -- convert to a new particle instance in order to access the positions
    local particle_exc = group_excluded:to_particle()
    local particle_inc = group_included:to_particle()

    -- for included/excluded make sure that the particles have been sorted
    -- into the respecting group correctly
    local positions_inc = particle_inc:get_position()
    for i = 1, group_included.size do
        local p = positions_inc[i]
        for d = 1, dimension do
            local l = p[d] - cuboid.origin[d]
            assert(l < cuboid.length[d] and l > 0)
        end
    end
    local positions_exc = particle_exc:get_position()
    for i = 1, group_excluded.size do
        local p = positions_exc[i]
        local outside = false
        for d = 1, dimension do
            local l = p[d] - cuboid.origin[d]
            if l > cuboid.length[d] or l < 0 then
                outside = true
            end
        end
        assert(outside)
    end
end

--
-- Parse command-line arguments.
--
local function parse_args()
    local parser = halmd.utility.program_options.argument_parser()

    parser:add_argument("output,o",
        {type = "string", default = "from_region_test", help = "prefix of output files"})
    parser:add_argument("verbose,v", {type = "accumulate", action = function(args, key, value)
        local level = {
            -- console, file
            {"warning", "info" },
            {"info"   , "info" },
            {"debug"  , "debug"},
            {"trace"  , "trace"},
        }
        args[key] = level[value] or level[#level]
    end, default = 1, help = "increase logging verbosity"})

    parser:add_argument("particles", {type = "number", default = 10000, help = "number of particles"})
    parser:add_argument("dimension", {type = "number", default = 3, help = "dimension of space"})

    return parser:parse_args()
end

local args = parse_args()

-- log to console
halmd.io.log.open_console({severity = args.verbose[1]})
-- log to file
halmd.io.log.open_file(("%s.log"):format(args.output), {severity = args.verbose[2]})
-- log version
halmd.utility.version.prologue()

-- set up system and perform test
test(setup(args))

-- log profiler results
halmd.utility.profiler:profile()
