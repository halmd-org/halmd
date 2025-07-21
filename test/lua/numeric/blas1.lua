#!/usr/bin/env halmd
--
-- Copyright © 2023 Viktor Skoblin
-- Copyright © 2023 Felix Höfling
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
local numeric = halmd.numeric

function test_vector_ops(vector1, vector2, scalar)
    local assert_close = numeric.assert_close

    local sum = numeric.add(vector1, vector2)
    for i, x in ipairs(sum) do
        assert_close(x, vector1[i] + vector2[i])
    end

    local sub1 = numeric.subtract(sum, vector1)
    for i, x in ipairs(sub1) do
        assert_close(x, vector2[i])
    end

    local sum_plus_scalar = numeric.add_scalar(sum, scalar)
    for i, x in ipairs(sum_plus_scalar) do
        assert_close(x, vector1[i] + vector2[i] + scalar)
    end

    local subtract_scalar = numeric.subtract_scalar(sum_plus_scalar, scalar)
    for i, x in ipairs(subtract_scalar) do
        assert_close(x, sum[i])
    end

    local product = numeric.multiply(vector1, vector2)
    for i, x in ipairs(product) do
        assert_close(x, vector1[i] * vector2[i])
    end

    local quotient = numeric.divide(product, vector1)
    for i, x in ipairs(quotient) do
        assert_close(x, vector2[i])
    end

    local product_scalar = numeric.multiply_scalar(scalar, vector1)
    for i, x in ipairs(product_scalar) do
        assert_close(x, scalar * vector1[i])
    end

    local divide_by_scalar = numeric.divide_scalar(product_scalar, scalar)
    for i, x in ipairs(divide_by_scalar) do
        assert_close(x, vector1[i])
    end
end

function test_norm()
    local vector = {-2, 5, 0, 3}
    local norm1 = numeric.norm_1(vector)
    assert(norm1 == 10)

    local norm2 = numeric.norm_2(vector)
    numeric.assert_close(norm2, math.sqrt(38))
end

function main()
    test_vector_ops({1, 1. / 4, 6}, {1./ 2, 6, 2}, 4)
    test_vector_ops({4, 1. / 7, -1}, {9, 2, 1./ 3}, -1. / 7)
    test_norm()
end
