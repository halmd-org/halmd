--
-- Copyright © 2014 Felix Höfling
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

local mdsim = halmd.mdsim
local sampler = halmd.observables.sampler

function rel_error(a, b)
    return math.abs(a - b) / math.abs(b)
end

function test_construction(args)
    -- construct prerequisites
    local box = mdsim.box({length = args.box_length})
    local particle = mdsim.particle({particles = args.particles, dimension = box.dimension})

    -- construct integrator
    local integrator = mdsim.integrators.verlet_nvt_hoover({
        particle = particle, box = box
      , timestep = args.timestep
      , temperature = args.temperature
      , resonance_frequency = args.frequency
    })

    -- check parameter passing
    assert(rel_error(integrator.timestep, args.timestep) < args.parameter_tolerance)
    assert(rel_error(integrator.temperature, args.temperature) < args.parameter_tolerance)

    -- check integrator masses
    local dimension = box.dimension
    local omega = 2 * math.pi * args.frequency
    assert(rel_error(integrator.mass[1], dimension * particle.nparticle * args.temperature / (omega * omega)) < 1e-6)
    assert(rel_error(integrator.mass[2], args.temperature / (omega * omega)) < 1e-6)

    return integrator
end

function test_methods(integrator)
    integrator:set_timestep(0.002)
    assert(rel_error(integrator.timestep, 0.002) < 1e-6)

    integrator:set_temperature(1.0)
    assert(rel_error(integrator.temperature, 1.0) < 1e-6)

    local mass = {10, 1}
    integrator:set_mass(mass)
    for i,m in ipairs(integrator.mass) do
        assert(rel_error(m, mass[i]) < 1e-6)
    end

    integrator:integrate()
    integrator:finalize()
end

function test_writer(integrator, args)
    if args.output then
        local file = halmd.io.writers.h5md({path = ("%s.h5"):format(args.output)})
        local writer = integrator:writer({
              file = file
            , fields = {"position", "velocity", "internal_energy"}
            , every = 1
        })
        sampler:sample() -- sample current state
        sampler:run(1)   -- step integrator and sample again
        writer:disconnect()
    end
end

-- define command line arguments
function define_args(parser)
    parser:add_argument("particles", {type = "integer", default=1600, help = "number of particles"})
    parser:add_argument("box-length", {type = "vector", dtype = "number", default = {20, 40, 20}
      , help = "edge lengths of simulation box"
    })
    parser:add_argument("timestep", {type = "number", default=0.002, help = "integration time step"})
    parser:add_argument("temperature", {type = "number", default=3.0, help = "temperature"})
    parser:add_argument("frequency", {type = "number", default=5.0, help = "resonance frequency"})
    parser:add_argument("output,o", {type = "string", help = "prefix of output files"})
    parser:add_argument("parameter-tolerance", {type = "number", default=1e-15
      , help = "parameter passing tolerance"
    })
end

-- start tests
function main(args)
    local integrator = test_construction(args)
    test_methods(integrator, args)
    test_writer(integrator, args)

    integrator:disconnect()
end
