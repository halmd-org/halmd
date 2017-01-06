--
-- Copyright © 2016 Felix Höfling
-- Copyright © 2016 Arthur Straube
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

local halmd = require("halmd")
local mdsim = require("halmd.mdsim")
local observables = require("halmd.observables")
local sampler = require("halmd.observables.sampler")

function rel_error(a, b)
    return math.abs(a - b) / math.abs(b)
end

function test_construction(args)
    -- construct prerequisites
    local box = mdsim.box({length = args.box_length})
    local particle = mdsim.particle({particles = args.particles, species = 2, dimension = box.dimension})
    mdsim.positions.lattice({particle = particle, box = box}):set()

    -- construct chemical potential module
    local chemical_potential = observables.chemical_potential({
        box = box
      , particle = particle
      , temperature = args.temperature
      , test_particles = args.test_particles
      , position = "random"
    })

    -- define interaction: Kob-Andersen mixture
    local potential = mdsim.potentials.pair.lennard_jones({
        epsilon = {
            {1  , 1.5} -- AA, AB
          , {1.5, 0.5} -- BA, BB
        }
      , sigma = {
            {1  , 0.8 } -- AA, AB
          , {0.8, 0.88} -- BA, BB
        }
    })
    potential = potential:truncate({cutoff = 2.5})
    chemical_potential:add_force({"pair"
      , particle = particle
      , potential = potential
    })

    -- check parameter passing
    assert(rel_error(chemical_potential.temperature, args.temperature) < 1e-15)

    return chemical_potential
end

function test_methods(chemical_potential)
    chemical_potential.temperature = 1.5    -- set temperature
    assert(rel_error(chemical_potential.temperature, 1.5) < 1e-6)

    chemical_potential:set_position()
    chemical_potential:sample(0)(nil) -- FIXME nil shouldn't be needed
end

function test_writer(chemical_potential, args)
    if args.output then
        local file = halmd.io.writers.h5md({path = ("%s.h5"):format(args.output)})
        local writer = chemical_potential:writer({file = file, every = 1, species = {"A", "B"}})

        -- write positions of test particles
        local box = mdsim.box({length = args.box_length})
        local test_group = mdsim.particle_groups.all({particle = chemical_potential.test_particle})
        observables.phase_space({box = box, group = test_group})
            :writer({file = file, fields = {"position"}, every = 1})

        sampler:sample() -- sample current state
        writer:disconnect()
    end
end

function test_single_particle(args)
    -- construct prerequisites
    local box = mdsim.box({length = {5, 5, 5}})
    local particle = mdsim.particle({particles = 1, species = 1, dimension = box.dimension})
    particle:set_position({{0,0,0}})

    -- construct chemical potential module
    local chemical_potential = observables.chemical_potential({
        box = box
      , particle = particle
      , temperature = 2
      , test_particles = {8e6} -- maximum particle number supported for CUDA compute capability 2.0
      , position = "lattice"
    })

    -- define interaction
    chemical_potential:add_force({"pair"
      , particle = particle
      , potential = mdsim.potentials.pair.lennard_jones({})
    })

    -- compute chemical potential and check with result
    -- obtained from quadrature of the Boltzmann factor over the box
    local result = chemical_potential:sample(0)(nil)
    print(("chemical potential of a single Lennard-Jones particle in the box: %.4g ± %.1g")
            :format(result:mean(), result:error_of_mean()))
    local reference = -0.0365878169
    if not (math.abs(result:mean() - reference) < 5 * result:error_of_mean()) then
        error(("mismatch in chemical potential [%f != %f]"):format(result:mean(), reference))
    end
end

-- define command line arguments
function define_args(parser)
    parser:add_argument("output,o", {type = "string", help = "prefix of output files"})
    parser:add_argument("particles", {type = "integer", default = 1600, help = "number of particles"})

    parser:add_argument("box-length", {type = "vector", dtype = "number", default = {9, 10, 11}
      , help = "edge lengths of simulation box"
    })

    parser:add_argument("temperature", {type = "number", default = 3.0, help = "temperature"})

    parser:add_argument("test-particles", {type = "vector", dtype = "number", default = {100, 200}
      , help = "number of test particles"
    })
end

-- start tests
function main(args)
    local chemical_potential = test_construction(args)
    test_methods(chemical_potential, args)
    test_writer(chemical_potential, args)

    chemical_potential:disconnect()

    if #args.box_length == 3 then
        test_single_particle(args)
    end
end

