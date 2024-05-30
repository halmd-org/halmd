#!/usr/bin/env halmd
--
-- Copyright © 2011-2012 Felix Höfling
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
local log = halmd.io.log
local mdsim = halmd.mdsim
local observables = halmd.observables
local readers = halmd.io.readers
local writers = halmd.io.writers
local utility = halmd.utility

--
-- Setup and run simulation
--
function main(args)
    local timestep = 0.001   -- integration timestep
    local steps = 10000      -- number of integration steps
    local count = args.count -- number of repetitions

    -- open H5MD trajectory file and read phase space data for A and B particles
    local file = readers.h5md({path = args.trajectory})

    local samples = {}
    local nparticle = 0
    local edges
    for i, label in ipairs({'A', 'B'}) do
        -- construct a phase space reader and sample
        local reader, sample = observables.phase_space.reader({
            file = file
          , location = {"particles", label}
          , fields = {"position", "velocity", "species", "mass"}
        })
        samples[label] = sample
        -- read phase space sample at last step in file
        log.info("number of %s particles: %d", label, sample.nparticle)
        reader:read_at_step(-1)
        -- read edge vectors of simulation domain from particle group
        edges = mdsim.box.reader({file = file, location = {"particles", label}})
        -- determine system parameters from phase space sample
        nparticle = nparticle + sample.nparticle
    end

    -- close H5MD trajectory file
    file:close()

    -- create simulation box
    local box = mdsim.box({edges = edges})

    -- open H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output)})

    -- create system state
    local particle = mdsim.particle({dimension = 3, particles = nparticle, species = 2, precision = args.precision})

    -- setup and sample each particle group separately
    local offset = 0
    for label, sample in utility.sorted(samples) do
        local group = mdsim.particle_groups.id_range({
            particle = particle
          , range = {offset + 1, offset + sample.nparticle}
          , label = label
        })
        offset = offset + sample.nparticle

        -- construct phase space instance
        local phase_space = observables.phase_space({box = box, group = group})
        -- set particle positions, velocities, species
        phase_space:set(sample)
        -- write phase space data of group to H5MD file, but only first and last step
        phase_space:writer({file = file, fields = {"position", "velocity", "species", "mass"}, every = steps})
    end

    -- define interaction of Kob-Andersen mixture using truncated Lennard-Jones potential
    local potential = mdsim.potentials.pair.lennard_jones({
        epsilon = {{1, 1.5}, {1.5, 0.5}} -- ((AA, AB), (BA, BB))
      , sigma = {{1, 0.8}, {0.8, 0.88}} -- ((AA, AB), (BA, BB))
    }):truncate({cutoff = 2.5})
    -- compute pair forces
    --
    -- use smaller skin width and increased neighbour list occupancy as the
    -- high density of the system makes diffusion slow and suppresses density
    -- fluctuations.
    local force = mdsim.forces.pair({
        box = box, particle = particle, potential = potential
      , neighbour = { skin = 0.3, occupancy = 0.7 }
    })

    -- define velocity-Verlet integrator
    local integrator = mdsim.integrators.verlet({
        box = box, particle = particle, timestep = timestep
    })

    -- estimate remaining runtime
    observables.runtime_estimate({steps = count * steps, first = 10, interval = 900, sample = 60})

    -- sample initial state
    observables.sampler:sample()

    -- run simulation several times and output profiling data
    for i = 1, count do
        observables.sampler:run(steps)
        utility.profiler:profile()
    end
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.substitute_date_time_action,
        default = "kob_andersen_benchmark_%Y%m%d_%H%M%S", help = "prefix of output files"})

    parser:add_argument("trajectory", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD trajectory file"})

    parser:add_argument("count", {type = "number", default = 5, help = "number of repetitions"})
    parser:add_argument("precision", {type = "string", help = "floating-point precision"})
end
