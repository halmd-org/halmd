#!/usr/bin/env halmd
--
-- Copyright © 2010-2023 Felix Höfling
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
local numeric = halmd.numeric
local observables = halmd.observables
local dynamics = halmd.observables.dynamics
local log = halmd.io.log
local readers = halmd.io.readers
local writers = halmd.io.writers
local utility = halmd.utility

-- search definition files in the top-level path relative to the simulation script
package.path = utility.abspath("../?.lua;") .. package.path
local definitions = { lennard_jones = require("definitions/lennard_jones") }

--
-- Setup and run simulation
--
function main(args)
    -- open H5MD file for reading
    local file = readers.h5md({path = args.input})

    -- construct a phase space reader and sample
    local reader, sample = observables.phase_space.reader({
        file = file, location = {"particles", "all"}, fields = {"position", "velocity"}
    })
    -- read phase space sample at last step in file
    reader:read_at_step(-1)
    -- determine system parameters from phase space sample
    local nparticle = assert(sample.nparticle)
    local dimension = assert(sample.dimension)
    if (assert(sample.nspecies) ~= 1) then
        log.error("Simulation script expects input file with particle data for a single species.")
    end

    -- read edge vectors of simulation domain from file
    local edges = mdsim.box.reader({file = file, location = {"particles", "all"}})
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({edges = edges})

    -- close H5MD file
    file:close()

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle})

    -- select all particles
    local all_group = mdsim.particle_groups.all({particle = particle})

    -- use phase space "sampler" to set particle positions and velocities
    -- by converting the sample obtained from the file,
    -- the sampler will be reused below for observables and file output
    local phase_space = observables.phase_space({box = box, group = all_group})
    phase_space:set(sample)

    -- define Lennard-Jones pair potential (with parameters ε=1 and σ=1 for a single species)
    -- and register computation of pair forces
    definitions.lennard_jones.create_pair_force({
        box = box, particle = particle, cutoff = args.cutoff, smoothing = args.smoothing
    })

    -- add velocity-Verlet integrator
    local integrator = mdsim.integrators.verlet({
        box = box, particle = particle, timestep = args.timestep
    })

    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- write trajectory of particle groups to H5MD file
    local interval = args.sampling.trajectory or steps
    if interval > 0 then
        -- reuse instance of phase space sampler from above
        phase_space:writer({file = file, fields = {"position", "velocity"}, every = interval})
    end

    -- sample macroscopic state variables
    local msv = observables.thermodynamics({box = box, group = all_group})
    local interval = args.sampling.state_vars
    if interval > 0 then
        msv:writer({
            file = file
          , fields = {
                "potential_energy", "pressure", "temperature"  -- fluctuating quantities
              , "internal_energy", "center_of_mass_velocity"   -- conserved quantities
            }
          , every = args.sampling.state_vars
        })
    end

    -- set up wavevectors, compute density modes and static structure factor
    local density_mode
    local interval = args.sampling.structure
    if interval > 0 or args.sampling.correlation > 0 then
        -- set up wavevector grid compatible with the periodic simulation box
        local grid = args.wavevector.wavenumbers
        if not grid then
            grid = observables.utility.semilog_grid({
                start = 2 * math.pi / numeric.max(box.length)
              , stop = args.wavevector.maximum
              , decimation = args.wavevector.decimation
            }).value
        end
        wavevector = observables.utility.wavevector({
            box = box, wavenumber = grid
          , tolerance = args.wavevector.tolerance, max_count = args.wavevector.max_count
        })

        -- compute density modes and output their time series,
        density_mode = observables.density_mode({group = all_group, wavevector = wavevector})
        if interval > 0 then
            density_mode:writer({file = file, every = interval})
        end

        -- compute static structure factor from density modes
        local ssf = observables.ssf({density_mode = density_mode, norm = nparticle})

        -- output averages over a certain number of configurations each
        local average = args.sampling.average
        if average and average > 0 then
            halmd.io.log.warning("Averaging of static structure factors not yet supported")
--            local total_ssf = observables.utility.accumulator({
--                acquire = ssf.acquire, every = interval, desc = "ssf"
--            })
--            total_ssf:writer({
--                file = file
--              , location = {"structure", ssf.label, "static_structure_factor"}
--              , every = average * interval
--              , reset = true})
        else
            ssf:writer({file = file, every = interval})
        end
    end

    -- time correlation functions
    local interval = args.sampling.correlation
    if interval > 0 then
        -- setup blocking scheme
        local max_lag = steps * integrator.timestep / 10
        local blocking_scheme = dynamics.blocking_scheme({
            max_lag = max_lag
          , every = interval
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

        if density_mode then
            -- compute intermediate scattering function
            local isf = dynamics.intermediate_scattering_function({density_mode = density_mode, norm = nparticle})
            blocking_scheme:correlation({tcf = isf, file = file })
        end
    end

    -- rescale velocities of all particles
    if args.rescale_to_energy then
        local rescale_velocity = require("rescale_velocity")
        rescale_velocity({msv = msv, internal_energy = args.rescale_to_energy})
    end

    -- sample initial state
    observables.sampler:sample()

    -- estimate remaining runtime
    observables.runtime_estimate({steps = steps})

    -- run simulation
    observables.sampler:run(steps)
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "lennard_jones_rc{cutoff:g}_%Y%m%d_%H%M%S", help = "basename of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})

    parser:add_argument("input", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD input file"})
    parser:add_argument("rescale-to-energy", {type = "number", help = "rescale velocities to match given internal energy"})

    parser:add_argument("cutoff", {type = "float32", default = math.pow(2, 1 / 6), help = "potential cutoff radius"})
    parser:add_argument("smoothing", {type = "number", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("time", {type = "number", default = 100, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
    sampling:add_argument("structure", {type = "integer", default = 1000, help = "for density modes, static structure factor"})
    sampling:add_argument("correlation", {type = "integer", default = 100, help = "for correlation functions"})
    sampling:add_argument("average", {type = "integer", help = "output averages of given number of samples"})

    local wavevector = parser:add_argument_group("wavevector", {help = "wavevector shells in reciprocal space"})
    observables.utility.wavevector.add_options(wavevector, {tolerance = 0.01, max_count = 7})
    observables.utility.semilog_grid.add_options(wavevector, {maximum = 15, decimation = 0})
end
