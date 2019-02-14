#!/usr/bin/env halmd
--
-- Copyright © 2018 Felix Höfling
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

local log = require("halmd.io.log")
local mdsim = require("halmd.mdsim")
local utility = require("halmd.observables.utility")

function test_dense_grid_2d()
    -- parameters
    local pi = math.pi
    local q = { 4, 3, 1, 1.3 }  -- unordered
    local box = mdsim.box({length = {.5 * pi, 4 * pi}})

    -- construct wavevectors
    local wavevector = utility.wavevector({wavenumber=q, box=box, dense=true})

    log.debug("Wavevectors:")
    local q_vec = wavevector:value()
    for i,v in ipairs(q_vec) do
        log.debug("  %d:\tq = (%g, %g),\t|q| = %g", i-1, v[1], v[2], math.sqrt(v[1]*v[1] + v[2]*v[2]))
    end

    log.debug("Wavevector shells:")
    q = wavevector:wavenumber()     -- wavenumbers are sorted
    local shell = wavevector:shell()
    for i,v in ipairs(shell) do
        log.debug("  %g:\t[%d, %d)", q[i], v[1], v[2])
    end

    local shell_upper = {3, 5, 11, 15}
    assert(#q == #shell)
    assert(#q == #shell_upper)
    assert(#q_vec == shell_upper[#shell_upper])
    for i,v in ipairs(shell) do
        assert(shell_upper[i] == v[2], ("shell_upper[%d] != v[2] (%d != %d)"):format(i, shell_upper[i], v[2]))
    end
end

function main()
    test_dense_grid_2d()
end
