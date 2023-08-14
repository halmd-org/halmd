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

local log = halmd.io.log
local mdsim = halmd.mdsim
local observables = halmd.observables
local writers = halmd.io.writers

local function setup(geometry, args)
    local dimension = args.dimension      -- dimension of space
    local np = args.particles             -- number of particles
    local L = args.box_length             -- edge length of cubic box

    local length = {}
    for i = 1, dimension do
        length[i] = L
    end
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})

    -- create system state for all particles first
    local particle = mdsim.particle({dimension = dimension, particles = np, species = 1})

    -- set initial particle positions
    mdsim.positions.lattice({box = box, particle = particle}):set()

--    if particle.memory == "gpu" then  -- FIXME why not for host, missing 'binning' - skip sorting at all?
--      -- sort particles *once*, as they do not move
--      mdsim.sorts.hilbert({box = box, particle = particle})
--        : order()
--    end

    -- construct included/excluded particle groups, label is inherited from 'region'
    local group = {}
    group["included"] = mdsim.particle_groups.region({
        particle = particle, geometry = geometry, selection = "included"
      , label = "included"
    })
    group["excluded"] = mdsim.particle_groups.region({
        particle = particle, geometry = geometry, selection = "excluded"
      , label = "excluded"
    })

    return group
end

-- named tabled of test functions
test = {}

--
-- test cuboid geometry: select upper quadrant of the box
--
test["cuboid"] = function(args)
    local dimension = args.dimension      -- dimension of space
    local L = args.box_length             -- edge length of cubic box

    -- select particles within/not within upper quadrant:
    -- define geometry first
    local length = {}
    local lowest_corner = {}
    for i = 1, dimension do
        lowest_corner[i] = 0
        length[i] = L / 2
    end
    local cuboid = mdsim.geometries.cuboid({lowest_corner = lowest_corner, length = length})

    local group = setup(cuboid, args)

    -- check if the total number of particles is correct
    log.info(("%d particles included in selection"):format(group["included"].size))
    log.info(("%d particles excluded from selection"):format(group["excluded"].size))
    assert(group["excluded"].size + group["included"].size == args.particles)

    -- convert groups to new particle instance in order to access the positions
    local positions = {}
    for label, g in pairs(group) do
        positions[label] = g:to_particle().data.position
    end

    -- for included/excluded check that the particles have been sorted
    -- into the respective group correctly
    for i, r in ipairs(positions["included"]) do
        for j = 1, #r do
            local dr = r[j] - lowest_corner[j]
            assert(dr <= length[j] and dr >= 0, ("particle #%d included in selection, but it should not"):format(i))
        end
    end

    for i, r in ipairs(positions["excluded"]) do
        local outside = false
        for j = 1, #r do
            local dr = r[j] - lowest_corner[j]
            if dr > length[j] or dr < 0 then
                outside = true
            end
        end
        assert(outside, ("particle #%d excluded from selection, but it should not"):format(i))
    end
end

--
-- test spherical geometry: half-sphere placed at the face with x=L/2 of the box
--
test["sphere"] = function(args)
    local dimension = args.dimension      -- dimension of space
    local L = args.box_length             -- edge length of cubic box

    -- select particles within/not within a sphere that sticks out of the simulation box
    -- define geometry first
    local radius = L / 2
    local centre = { L / 2 }
    for i = #centre + 1, dimension do
        centre[i] = 0
    end
    local sphere = mdsim.geometries.sphere({centre = centre, radius = radius})

    local group = setup(sphere, args)

    -- check if the total number of particles is correct
    log.info(("%d particles included in selection"):format(group["included"].size))
    log.info(("%d particles excluded from selection"):format(group["excluded"].size))
    assert(group["excluded"].size + group["included"].size == args.particles)

    -- convert groups to new particle instance in order to access the positions
    local positions = {}
    for label, g in pairs(group) do
        positions[label] = g:to_particle().data.position
    end

    -- for included/excluded check that the particles have been sorted
    -- into the respective group correctly
    for i, r in ipairs(positions["included"]) do
        local distance = 0
        for j = 1, #r do
            local dr = r[j] - centre[j]
            distance = distance + dr * dr
        end
        distance = math.sqrt(distance)
        assert(distance <= radius, ("particle #%d included in selection, but it should not"):format(i))
    end

    for i, r in ipairs(positions["excluded"]) do
        local distance = 0
        for j = 1, #r do
            local dr = r[j] - centre[j]
            distance = distance + dr * dr
        end
        distance = math.sqrt(distance)
        assert(distance > radius, ("particle #%d excluded from selection, but it should not"):format(i))
    end
end

--
-- test spherical geometry: half-sphere placed at the face with x=L/2 of the box
--
test["cylinder"] = function(args)
    local dimension = args.dimension      -- dimension of space
    local L = args.box_length             -- edge length of cubic box
    local tolerance = 1e-12               -- tolerance when verifying selection criteria

    -- select particles within/not within a cylinder that runs parallel to the
    -- z=0, x=y diagonal of the simulation box. The cylinder has a fixed width
    -- of 1.5σ.
    -- define geometry first
    local radius = 0.75
    local length = L
    local axis = { 1, 1 }
    for i = #axis + 1, dimension do
        axis[i] = 0
    end
    local centre = { math.sqrt(2) * radius }
    for i = #centre + 1, dimension do
        centre[i] = 0
    end

    local cylinder = mdsim.geometries.cylinder({axis = axis, centre = centre, radius = radius, length = length})

    local group = setup(cylinder, args)

    -- check if the total number of particles is correct
    log.info(("%d particles included in selection"):format(group["included"].size))
    log.info(("%d particles excluded from selection"):format(group["excluded"].size))
    assert(group["excluded"].size + group["included"].size == args.particles)

    -- convert groups to new particle instance in order to access the positions
    local positions = {}
    for label, g in pairs(group) do
        positions[label] = g:to_particle().data.position
    end

    -- normalise axis
    local norm = 0
    for j = 1, #axis do
        norm = norm + axis[j] * axis[j]
    end
    norm = math.sqrt(norm)
    for j = 1, #axis do
        axis[j] = axis[j] / norm
    end

    -- for included/excluded check that the particles have been sorted
    -- into the respective group correctly
    for i, r in ipairs(positions["included"]) do
        local dr2 = 0
        local dr_n = 0
        for j = 1, #r do
            local dr = r[j] - centre[j]
            dr2 = dr2 + dr * dr
            dr_n = dr_n + dr * axis[j]
        end
        local distance = math.sqrt(dr2 - dr_n * dr_n)
        assert(distance <= radius * (1 + tolerance) and math.abs(dr_n) <= length / 2 * (1 + tolerance),
               ("particle #%d included in selection, but it should not"):format(i))
    end

    for i, r in ipairs(positions["excluded"]) do
        local dr2 = 0
        local dr_n = 0
        for j = 1, #r do
            local dr = r[j] - centre[j]
            dr2 = dr2 + dr * dr
            dr_n = dr_n + dr * axis[j]
        end
        local distance = math.sqrt(dr2 - dr_n * dr_n)
        assert(distance > radius * (1 - tolerance) or math.abs(dr_n) > length / 2 * (1 - tolerance),
               ("particle #%d excluded from selection, but it should not"):format(i))
    end
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", default = "region_test", help = "prefix of output files"})

    parser:add_argument("run_test", {type = "string", help = "run only selected test case"})
    parser:add_argument("particles", {type = "number", default = 10000, help = "number of particles"})
    parser:add_argument("dimension", {type = "number", default = 3, help = "dimension of space"})
    parser:add_argument("box-length", {type = "number", default = 10, help = "edge length of cubic box"})
end

--
-- set up system and perform test
--
function main(args)
    -- run selected test case or, by default, all tests
    local test_case = args.run_test
    local cases = test_case and { test_case } or {"cuboid", "cylinder", "sphere"}

    for i,case in ipairs(cases) do
        log.message(("Running test case '%s' ..."):format(case))
        assert(test[case], ("test case '%s' is not registered"):format(case))

        -- call test function
        test[case](args)
        log.message(("Test case '%s' finished."):format(case))
        log.message("")
    end
end
