--
-- Copyright © 2019 Abbas Gholami
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

-- grab modules
local device = halmd.utility.device
local log = halmd.io.log
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local writers = halmd.io.writers
local utility = halmd.observables.utility

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

    log.info("Placing particles at random positions..")

    if not excluded then
        excluded = halmd.mdsim.positions.excluded_volume{box = box, cell_length = cell_length}
    end
    local A = 0
    local B = 0
    local obstacles = {}

    local diameter
    local repeats = 1000
    local cnt = 1
    for s = 1, nspecies do
        diameter = sigma[s][s]
        for i = 1, particles[s] do
            local placed = false
            for j = 1, repeats do
                -- generate position
                local r = {}

                r[1] = edges[1] / (math.floor (edges[1] / diameter)) * (i-1)

                r[2] = edges[2] / (math.floor (edges[2] / diameter)) * A
                if i%edges[1] == 0 then
                     A = A + 1
                end
                r[3] = edges[3] / (math.floor (edges[3] / diameter)) * B
                if i%(edges[1]*edges[2]) == 0 then
                     B = B + 1
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


-------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------

function main()

--PARAMETERS
    local box = mdsim.box({length = {10, 10, 10}})
    local nparticles = device.gpu and 10000 or 100
    local nknots = {101, 101, 101} --with 0.1 spacing
    local sigma = 1
    local temperature = 1
    local steps = device.gpu and 10000 or 100
    local timestep = 1e-1
    local amplitude = 1.5
    local period = 2.5

--INTERNAL PARTICLES
    local particle = mdsim.particle({dimension = box.dimension, particles = nparticles, species = 1})

--Initial Configuration
    local configuration = generate_configuration({particles = {nparticles}, sigma = sigma, box = box})
    particle.data["position"] = configuration.positions

--Initial Species
    local species = {}
    for i = 1, nparticles do table.insert(species, 0) end
    particle.data["species"] = species

--Initial Velocities
    local boltzmann = mdsim.velocities.boltzmann({particle = particle, temperature = temperature})
    boltzmann:set()

--POTENTIAL
     local potential = math.cos
     local potential_derivatives = math.sin
     local coefficients = {}
     local dx = box.length[1] / (nknots[1]-1)
     local T_1 = 2 * 3.14 / period

     for i = 1, 8*nknots[1]*nknots[2]*nknots[3], 8 do

         local j = (i-1)/8 + 1
         local residue = j-math.floor(j/nknots[1])*nknots[1]
         if residue == 0 then
            residue = 1
         end

         coefficients[i]   = amplitude*potential( T_1 * dx * ( residue - 1 ) )
         coefficients[i+1] = -amplitude * T_1 * potential_derivatives( T_1 * dx * ( residue - 1) )
         coefficients[i+2] = 0
         coefficients[i+3] = 0
         coefficients[i+4] = 0
         coefficients[i+5] = 0
         coefficients[i+6] = 0
         coefficients[i+7] = 0
     end

--INTERPOLATION
    local interpolation = mdsim.forces.interpolation.cubic_hermite({box = box, nknots = nknots, precision = "double"})
    local virial_interpolation = mdsim.forces.interpolation.linear({box = box, nknots = nknots, precision = particle.precision})

--FORCE IMPLEMENTATION
    local tabulated_forces = mdsim.forces.tabulated_external({particle = particle, box = box, interpolation = interpolation, virial_interpolation = virial_interpolation})

    --Reading Coefficients
    --local file = halmd.io.readers.h5md({path = "coefficients.h5"})
    tabulated_forces:set_coefficients(coefficients)

--WRITER
    -- H5MD file writer
    local file = halmd.io.writers.h5md({path = "tabulated_external.h5", overwrite = true})

    -- select all particles
    local particle_group = halmd.mdsim.particle_groups.all({particle = particle})

    -- sample phase space
    local phase_space = halmd.observables.phase_space({box = box, group = particle_group})

    -- write trajectory of particle groups to H5MD file
    phase_space:writer({file = file, fields = {"position", "velocity", "species", "mass"}, every = 1})

    -- Sample macroscopic state variables
    local msv = halmd.observables.thermodynamics({box = box, group = particle_group})
    msv:writer({file = file, every = 100})

    local accumulator = halmd.observables.utility.accumulator({acquire = msv.internal_energy, every = 1, desc = "Averaged internal energy", aux_enable = {particle}})
    accumulator:writer({file = file, location = {"observables", "averaged_internal_energy"}, every = 1, reset = true})

--RUN
    -- sample initial state
    halmd.observables.sampler:sample()

    local integrator = halmd.mdsim.integrators.verlet_nvt_boltzmann({box = box, particle = particle, timestep = timestep, temperature = temperature, rate = 2})

    -- estimate remaining runtime
    local runtime = halmd.observables.runtime_estimate({steps = steps, first = 10, interval = 900, sample = 60})

    -- run simulation
    halmd.observables.sampler:run(steps)

    halmd.utility.profiler:profile()
end

main()
