--
-- Copyright © 2014-2015 Nicolas Höft
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

local mdsim = halmd.mdsim
local observables = halmd.observables
local writers = halmd.io.writers

local function setup(args)
    local dimension = args.dimension      -- dimension of space
    local np = args.particles             -- number of particles

    local length = {}
    for i = 1, dimension do
        length[i] = 10
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

    -- select particles within/not within upper quadrant:
    -- define geometry first
    local lowest_corner = {}
    for i = 1, dimension do
        lowest_corner[i] = 0
        length[i] = length[i] / 2
    end
    local cuboid = mdsim.geometries.cuboid({lowest_corner = lowest_corner, length = length})

    -- construct included/excluded particle groups, label is inherited from 'region'
    local group = {}
    group["included"] = mdsim.particle_groups.region({
        particle = particle, geometry = cuboid, selection = "included"
      , label = "upper quadrant (included)"
    })
    group["excluded"] = mdsim.particle_groups.region({
        particle = particle, geometry = cuboid, selection = "excluded"
      , label = "upper quadrant (excluded)"
    })

    return group, cuboid, args
end

local function test(group, cuboid, args)
    -- check if the total number of particles is correct
    assert(group["excluded"].size + group["included"].size == args.particles)
    local lowest_corner = cuboid.lowest_corner
    local length = cuboid.length

    -- convert to a new particle instance in order to access the positions
    local particle_exc = group["excluded"]:to_particle()
    local particle_inc = group["included"]:to_particle()

    -- for included/excluded check that the particles have been sorted
    -- into the respective group correctly
    local positions_inc = particle_inc.data.position
    for i = 1, group["included"].size do
        local r = positions_inc[i]
        for j = 1, #r do
            local dr = r[j] - lowest_corner[j]
            assert(dr < length[j] and dr > 0, ("particle #%d not included in selection"):format(i))
        end
    end

    local positions_exc = particle_exc.data.position
    for i = 1, group["excluded"].size do
        local r = positions_exc[i]
        local outside = false
        for j = 1, #r do
            local dr = r[j] - lowest_corner[j]
            if dr > length[j] or dr < 0 then
                outside = true
            end
        end
        assert(outside, ("particle #%d not excluded from selection"):format(i))
    end
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", default = "region_test", help = "prefix of output files"})

    parser:add_argument("particles", {type = "number", default = 10000, help = "number of particles"})
    parser:add_argument("dimension", {type = "number", default = 3, help = "dimension of space"})
end

--
-- set up system and perform test
--
function main(args)
    test(setup(args))
end
