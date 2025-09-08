#!/usr/bin/env halmd
--
-- Copyright © 2011-2023 Felix Höfling
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
local dynamics = halmd.observables.dynamics
local readers = halmd.io.readers
local writers = halmd.io.writers
local utility = halmd.utility

-- search definition files in the top-level path relative to the simulation script
package.path = utility.abspath("../?.lua;") .. package.path

--
-- restore phase space point from file, using a unique instance of the
-- particle module
--
local function restore(args)
    -- open H5MD file for reading
    local file = readers.h5md({path = args.input})

    local samples = {}
    local nparticle = 0
    local nspecies = 0
    local edges
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
        log.info("number of %s particles: %d", label, sample.nparticle)
        reader:read_at_step(-1)
        -- read edge vectors of simulation domain from particle group
        edges = mdsim.box.reader({file = file, location = {"particles", label}})
        -- determine system parameters from phase space sample
        nparticle = nparticle + sample.nparticle
        nspecies = nspecies + 1
        label = string.char(string.byte(label) + 1)
    end
    local group = nil -- let garbage collector close the HDF5 group (hopefully)
    local dimension = assert(samples.A.dimension)

    -- close H5MD file
    file:close()

    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({edges = edges})

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle, species = nspecies})

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
    
    -- register computation of pair forces
    local force = mdsim.forces.pair({ box = box, particle = particle, potential = potential})

    return box, particle, samples, args
  
end

local function production(box, particle, samples, args)
  
  -- convert integration time to number of steps
  local steps = math.ceil(args.time / args.timestep)

  -- add velocity-Verlet integrator
  local integrator = mdsim.integrators.verlet({
    box = box, particle = particle, timestep = args.timestep
  })

  -- H5MD file writer
  local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

  -- sample each particle group separately
  local offset = 0
  for label, sample in utility.sorted(samples) do
    -- select particles of species
      local group = mdsim.particle_groups.id_range({
          particle = particle
        , range = {offset + 1, offset + sample.nparticle}
        , label = label
      })
      offset = offset + sample.nparticle

      -- sample phase space
      local phase_space = observables.phase_space({box = box, group = group})

      -- set particle positions, velocities, species
      phase_space:set(sample)

      -- write phase space trajectory to H5MD file
      phase_space:writer({
          file = file
        , fields = {"position", "velocity", "species", "mass"}
        , every = args.sampling.trajectory
      })

      -- sample macroscopic state variables
      local thermo = observables.thermodynamics({box = box, group = group}):writer({
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
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "one_particle_instance_%Y%m%d_%H%M%S", help = "basename of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = true, help = "overwrite output file"})

    parser:add_argument("input", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD input file"})

    parser:add_argument("time", {type = "number", default = 1000, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})
    parser:add_argument("random-seed", {type = "integer", action = parser.action.random_seed, help = "seed for random number generator"})
    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
      
end
