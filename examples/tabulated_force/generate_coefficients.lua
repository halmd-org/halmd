#!/usr/bin/env halmd
--
-- Copyright © 2013  Nicolas Höft
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

-- grab modules
local log = halmd.io.log
local mdsim = halmd.mdsim
local numeric = halmd.utility.numeric
local observables = halmd.observables
local writers = halmd.io.writers

local function flatten(list)
  if type(list) ~= "table" then return {list} end
  local flat_list = {}
  for _, elem in ipairs(list) do
    for _, val in ipairs(flatten(elem)) do
      flat_list[#flat_list + 1] = val
    end
  end
  return flat_list
end

local function generate_configuration(args)
    local particles = args.particles
    local box       = args.box
    local sigma     = args.sigma
    local excluded  = args.excluded_volume
    --local obstacles = args.obstacles

    local edges = box.length
    local nspecies = #particles
    local dimension = #edges

    if type(sigma) == "number" then
        sigma = numeric.scalar_matrix(nspecies, nspecies, sigma)
    end

    local cell_length = numeric.max(flatten(sigma))
    print(cell_length)

    log.info("Placing particles at random positions..")
    local random = math.random
    math.randomseed(os.time())

    if not excluded then
        excluded = halmd.mdsim.positions.excluded_volume{box = box, cell_length = cell_length}
    end

    local obstacles = {}

    local diameter
    local repeats = 1000
    cnt = 1
    for s = 1, nspecies do
        diameter = 1.2*sigma[s][s]
        for i = 1, particles[s] do
            local placed = false
            for j = 1, repeats do
                -- generate position
                local r = {}
                for d = 1, dimension do
                    r[d] = edges[d] * random()
                end
                -- test if position is valid
                if excluded:place_sphere(r, diameter) then
                    -- place particle and mark it as excluded volume
                    obstacles[cnt] = r
                    cnt = cnt + 1
                    placed = true
                    excluded:exclude_sphere(r, diameter)
                    break
                end
            end
            if not placed then
                error(("cannot place obstacle %d of species %d after %d repeats"):format(i, s, repeats))
            end
        end
    end
    log.info("Random position placement done. Set %d particles.", #obstacles)
    return {positions = obstacles, excluded_volume = excluded}
end

local function main(args)
    local box = mdsim.box({length = {10,5,1}})
    local nparticle_ext = 10
    local particle_external = mdsim.particle({box = box, particles = nparticle_ext, species = 2, label="E", memory = "host"})
    local species_ext = {}
    for i = 1, nparticle_ext do
        -- assign species 1 to all particles forming the external potential
        -- species 0 would be our test particle
        table.insert(species_ext, 1)
    end
    particle_external:set_species(species_ext)

    -- set postions of particle matrix
    local obstacle_config_ext = generate_configuration({
        particles = {nparticle_ext}
      , sigma = 1
      , box = box
    })
    particle_external:set_position(obstacle_config_ext.positions)


    local nknots = {51, 26, 6} -- separation of 0.2

    -- truncated Lennard-Jones potential
    local potential = mdsim.potentials.lennard_jones({
      , epsilon = {
            {1  , 1.5} -- AA, AB
          , {1.5, 0.5} -- BA, BB
        }
      , sigma = {
            {1  , 0.8 } -- AA, AB
          , {0.8, 0.88} -- BA, BB
        }
      , cutoff = 2.5
      , memory = "host"
    })

    local tabulated_generator = mdsim.forces.tabulated_generator({potential = potential, particle = particle_external, box = box, grid_size = nknots})

    local file = writers.h5md({path = "coefficients.h5"})
    --local coefficients = tabulated_generator:get_coefficients()
    tabulated_generator:write_coefficients({file = file})

end

--
-- Parse command-line arguments.
--
local function parse_args()
    local parser = halmd.utility.program_options.argument_parser()

    parser:add_argument("output,o", {type = "string", action = function(args, key, value)
        -- substitute current time
        args[key] = os.date(value)
    end, default = "test_%Y%m%d", help = "prefix of output files"})

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

    return parser:parse_args()
end

local args = parse_args()

-- log to console
halmd.io.log.open_console({severity = args.verbose[1]})
-- log to file
-- halmd.io.log.open_file(("%s.log"):format(args.output), {severity = args.verbose[2]})
-- log version
halmd.utility.version.prologue()

-- run simulation
main(args)
