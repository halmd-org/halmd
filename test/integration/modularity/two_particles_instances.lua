#!/usr/bin/env halmd
--
-- Copyright © 2013 Felix Höfling
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
local dynamics = halmd.observables.dynamics
local readers = halmd.io.readers
local writers = halmd.io.writers

--
-- restore phase space point from file, using distinct instances of the
-- particle module
--
local function restore(args)
    -- open H5MD file for reading
    local file = readers.h5md({path = args.input})

    local samples = {}
    local nspecies = 0
    local label = "A"
    local group = file.root:open_group("particles")
    while group:exists_group(label) do
        -- construct a phase space reader and sample
        local reader, sample = observables.phase_space.reader({
            file = file
          , location = {"particles", label}
          , fields = {"position", "velocity", "species", "mass"}
        })
        samples[label] = sample
        -- read phase space sample at last step in file
        log.message("number of %s particles: %d", label, sample.nparticle)
        reader:read_at_step(-1)
        -- increment label and species count
        label = string.char(string.byte(label) + 1)
        nspecies = nspecies + 1
    end
    local group = nil -- let garbage collector close the HDF5 group (hopefully)

    -- read edge vectors of simulation domain from file and recreate box with
    -- periodic boundary conditions
    local box = mdsim.box({edges = mdsim.box.reader({file = file, location = {"particles", "A"}})})

    -- close H5MD file
    file:close()

    -- create system state, one particle instance per species
    local particle = {}
    for label, sample in pairs(samples) do
        local p = mdsim.particle({
            dimension = box.dimension, particles = sample.nparticle, species = nspecies, label = label
        })
        observables.phase_space({
            box = box
          , group = mdsim.particle_groups.all({particle = p, global = false})
        }):set(sample)
        particle[label] = p
    end

    -- truncated Lennard-Jones potential
    local potential = mdsim.potentials.pair.lennard_jones({
        epsilon = {
            {1  , 1.5} -- AA, AB
          , {1.5, 0.5} -- BA, BB
        }
      , sigma = {
            {1  , 0.8 } -- AA, AB
          , {0.8, 0.88} -- BA, BB
        }
    }):truncate({"smooth_r4", cutoff = 2.5, h = 0.005})
    --[[
    -- create binning modules explicitly and therefore only once for each particle instance
    local binning = {
        A = mdsim.binning({
            box = box
          , particle = particle["A"]
          , r_cut = potential.r_cut
        })
      , B = mdsim.binning({
            box = box
          , particle = particle["B"]
          , r_cut = potential.r_cut
        })
    }
    --]]
    -- define interaction forces with smoothly truncated potential
    local force = {}
    for label1, p1 in pairs(particle) do
        for label2, p2 in pairs(particle) do
            local neighbour = mdsim.neighbour({
                box = box
              , particle = { p1, p2 }
              , r_cut = potential.r_cut
              --, binning = { binning[label1], binning[label2] }
            })
            force[label1 .. label2] = mdsim.forces.pair({
                box = box
              , particle = { p1, p2 }
              , potential = potential
              , label = label1 .. label2  -- FIXME do not infer logger from potential
              , neighbour = { disable_binning = true }
              --, neighbour = neighbour
            })
        end
    end
    
    return box, particle, args
end

local function production(box, particle, args)
    local timestep = args.timestep                -- integration timestep
    local steps = math.ceil(args.time / timestep) -- number of integration steps

    -- add velocity-Verlet integrators
    local integrator = {}
    for k,v in pairs(particle) do
        integrator[k] = mdsim.integrators.verlet({
            box = box, particle = v, timestep = timestep
        })
    end

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = true})

    -- total number of particles
    local nparticle = 0
    for k,v in pairs(particle) do
        nparticle = nparticle + v.nparticle
    end
  
    -- sample each particle instance (i.e., species) separately
    for label, p in pairs(particle) do
        -- select particles of species
        local group = mdsim.particle_groups.all({particle = p, global = false})

        -- write phase space trajectory to H5MD file
        local phase_space = observables.phase_space({box = box, group = group})
        phase_space:writer({
                file = file
              , fields = {"position", "velocity", "species", "mass", "potential_energy"}
              , every = args.sampling.trajectory
            })

        -- sample macroscopic state variables
        local thermo = observables.thermodynamics({box = box, group = group})
          : writer({
                file = file
              , fields = {
                    "potential_energy", "pressure", "temperature"  -- fluctuating quantities
                  , "internal_energy", "center_of_mass_velocity"   -- conserved quantities
                }
              , every = args.sampling.state_vars
            })
    
    end
       
    -- sample initial state
    observables.sampler:sample()

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({ steps = steps })

    -- run simulation
    observables.sampler:run(steps)
end

function main(args)
    -- restore simulation and run production
    production(restore(args))
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", default = "two_particles", help = "prefix of output files"})
    parser:add_argument("input", {type = "string", required = true, action = function(args, key, value)
        if not readers.h5md.check(value) then
            error(("not an H5MD file: %s"):format(value), 0)
        end
        args[key] = value
    end, help = "H5MD input file"})
    parser:add_argument("time", {type = "number", default = 1000, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})
    parser:add_argument("random-seed", {type = "integer", action = parser.action.random_seed, help = "seed for random number generator"})
    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
  
end
