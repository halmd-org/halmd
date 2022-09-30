#!/usr/bin/env halmd
--
-- Copyright © 2011-2021 Felix Höfling
-- Copyright © 2010-2012 Peter Colberg
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

-- grab modules
local mdsim = halmd.mdsim
local observables = halmd.observables
local random = halmd.random
local writers = halmd.io.writers
local utility = halmd.utility

--
-- Setup and run simulation
--
function main(args)
    local nparticle = args.tiny and 4096 or 256000  -- total number of particles
    local concentration = 0.8 -- concentration of A particles
    local density = 1.2       -- number density
    local temperature = 0.7   -- heat bath temperature

    local timestep = 0.005    -- integration timestep
    local steps = 20000       -- number of integration steps

    -- number of particles in groups A and B
    local ngroup = { math.floor(nparticle * concentration) }
    ngroup[2] = nparticle - ngroup[1]

    -- derive edge length of three-dimensional box from particle number and density
    local length = math.pow(nparticle / density, 1 / 3)

    -- create simulation box
    local box = mdsim.box({length = {length, length, length}})

    -- create system state
    local particle = mdsim.particle({dimension = 3, particles = nparticle, species = 2})

    -- set particle species, with continuous range of IDs per species:
    -- construct array with particle species: (0, 0, … 0, 1, 1, … 1)
    local species = {}
    for i = 1, ngroup[1] do table.insert(species, 0) end
    for i = 1, ngroup[2] do table.insert(species, 1) end
    particle.data["species"] = species

    -- set initial particle positions, randomise the particle species
    mdsim.positions.lattice({box = box, particle = particle}):set()
    -- randomly shuffle the positions
    particle.data["position"] = random.generator({memory = "host"}):shuffle(particle.data["position"])

    -- set initial particle velocities
    mdsim.velocities.boltzmann({particle = particle, temperature = temperature}):set()

    -- define interaction of Kob-Andersen mixture using truncated Lennard-Jones potential
    local potential = mdsim.potentials.pair.lennard_jones({
        epsilon = {{1, 1.5}, {1.5, 0.5}} -- ((AA, AB), (BA, BB))
      , sigma = {{1, 0.8}, {0.8, 0.88}} -- ((AA, AB), (BA, BB))
    }):truncate({cutoff = 2.5})
    -- compute forces
    local force = mdsim.forces.pair({
        box = box, particle = particle, potential = potential, neighbour = { unroll_force_loop = args.tiny }
    })

    -- define velocity-Verlet integrator with Andersen thermostat
    local integrator = mdsim.integrators.verlet_nvt_andersen({
        box = box, particle = particle,
        timestep = timestep, temperature = temperature, rate = 2
    })

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output)})

    -- sample each particle group separately and write to H5MD file
    local offset = 0
    for s = 0, 1 do
        local group = mdsim.particle_groups.id_range({
            particle = particle
          , range = {offset + 1, offset + ngroup[s + 1]}
          , label = string.char(string.byte('A') + s)
        })
        offset = offset + ngroup[s + 1]

        -- connect H5MD writer, output only first and last step
        observables.phase_space({box = box, group = group}):writer({
            file = file
          , fields = {"position", "velocity", "species", "mass"}
        })
    end

    -- estimate remaining runtime
    observables.runtime_estimate({steps = steps, first = 10, interval = 900, sample = 60})

    -- sample initial state
    observables.sampler:sample()

    -- run simulation
    observables.sampler:run(steps)

    -- log profiler results
    utility.profiler:profile()
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "kob_andersen_benchmark_configuration_%Y%m%d_%H%M%S", help = "prefix of output files"})
    parser:add_argument("tiny", {type = "boolean", help = "generate a tiny system (4096 particles)"})
end
