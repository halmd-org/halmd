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
local definitions = { kob_andersen = require("definitions/kob_andersen") }

--
-- Setup and run simulation
--
function main(args)
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

    -- define pair potential of Kob-Andersen binary mixture
    -- and register computation of pair forces
    definitions.kob_andersen.create_pair_force({
        box = box, particle = particle
    })

    -- add velocity-Verlet integrator
    local integrator = mdsim.integrators.verlet({
        box = box, particle = particle, timestep = args.timestep
    })

    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- set up wavevector grid compatible with the periodic simulation box
    -- if the computation of structural information is requested
    local wavevector
    if args.sampling.structure > 0 or args.sampling.correlation > 0 then
        local grid = args.wavevector.wavenumbers
        if not grid then
            grid = observables.utility.semilog_grid({
                start = 2 * math.pi / halmd.numeric.max(box.length)
              , stop = args.wavevector.maximum
              , decimation = args.wavevector.decimation
            }).value
        end
        wavevector = observables.utility.wavevector({
            box = box, wavenumber = grid
          , tolerance = args.wavevector.tolerance, max_count = args.wavevector.max_count
        })
    end

    -- sample each particle group separately
    local density_mode = {}
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
        -- write trajectory of particle groups to H5MD file
        phase_space:writer({
            file = file
          , fields = {"position", "velocity", "species", "mass"}
          , every = args.sampling.trajectory
        })

        -- sample macroscopic state variables
        local msv = observables.thermodynamics({box = box, group = group})
        msv:writer({
            file = file
          , fields = {
                "potential_energy", "pressure", "temperature"  -- fluctuating quantities
              , "internal_energy", "center_of_mass_velocity"   -- conserved quantities
            }
          , every = args.sampling.state_vars
        })

        -- sample structural data
        if wavevector then
            local interval = args.sampling.structure
            -- compute density modes and output their time series
            local density_mode_ = observables.density_mode({group = group, wavevector = wavevector})
            density_mode_:writer({file = file, every = interval})
            density_mode[label] = density_mode_

            -- compute static structure factor from correlation of density
            -- modes of previous particle groups and of this one
            for label, rho in pairs(density_mode) do
                local ssf = observables.ssf({density_mode = {rho, density_mode_}, norm = nparticle})
                -- output averages over a certain number of configurations each
                local average = args.sampling.average
                if average and average > 0 then
                    halmd.io.log.warning("Averaging of static structure factors not yet supported")
--                    local total_ssf = observables.utility.accumulator({
--                        acquire = ssf.acquire, every = interval, desc = "ssf " .. ssf.label
--                    })
--                    total_ssf:writer({
--                        file = file
--                      , location = {"structure", ssf.label, "static_structure_factor"}
--                      , every = average * interval
--                      , reset = true})
                else
                    ssf:writer({file = file, every = interval})
                end
            end
        end

        -- setup blocking scheme for correlation functions
        local max_lag = steps * integrator.timestep / 10
        local blocking_scheme = dynamics.blocking_scheme({
            max_lag = max_lag
          , every = args.sampling.correlation
          , size = 10
        })

        -- compute mean-square displacement
        local msd = dynamics.mean_square_displacement({phase_space = phase_space})
        blocking_scheme:correlation({tcf = msd, file = file})
        -- compute mean-quartic displacement
        local mqd = dynamics.mean_quartic_displacement({phase_space = phase_space})
        blocking_scheme:correlation({tcf = mqd, file = file})
        -- compute velocity autocorrelation function
        local vacf = dynamics.velocity_autocorrelation({phase_space = phase_space})
        blocking_scheme:correlation({tcf = vacf, file = file})

        -- compute intermediate scattering function from correlation of density
        -- modes of previous particle groups and of this one
        for l, rho in pairs(density_mode) do
            local isf = dynamics.intermediate_scattering_function({
                density_mode = (l == label) and rho or {rho, density_mode[label]}
              , norm = nparticle
            })
            blocking_scheme:correlation({tcf = isf, file = file})
        end
   end

    -- sample initial state
    observables.sampler:sample()

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({steps = steps})

    -- run simulation
    observables.sampler:run(steps)
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "binary_mixture_%Y%m%d_%H%M%S", help = "basename of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})

    parser:add_argument("input", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD input file"})

    parser:add_argument("time", {type = "number", default = 1000, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
    sampling:add_argument("structure", {type = "integer", default = 1000, help = "for density modes, static structure factor"})
    sampling:add_argument("correlation", {type = "integer", default = 100, help = "for correlation functions"})
    sampling:add_argument("average", {type = "integer", help = "output averages of given number of samples"})

    local wavevector = parser:add_argument_group("wavevector", {help = "wavevector shells in reciprocal space"})
    observables.utility.wavevector.add_options(wavevector, {tolerance = 0.01, max_count = 7})
    observables.utility.semilog_grid.add_options(wavevector, {maximum = 25, decimation = 0})
end
