#!/usr/bin/env halmd
--
-- Copyright © 2023 Viktor Skoblin
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
local utility = halmd.utility

function test_repeat(element, n)
    log.info(("Repeating element '%s' %d times ."):format(element, n))
    t = utility.repeat_element(element, n)
    assert(#t == n, ("Length of result list does not match: %d != %d"):format(#t, n))
    -- if the element is a table, we want to have a copy of this object but not the reference.
    if type(element) == "table" then
        for i = 1, #t do
            assert(t[i] ~= element)
            for j = 1, #t[i] do
                assert(t[i][j] == element[j])
            end
        end
    else
        for i = 1, #t do
            assert(t[i] == element)
        end
    end
end

function test_concat()
    local check = { 1, 4, 2, "spoon", 3, 78, 45, 3, "table", 2, 8, 235, 342 }
    local t1 = { 1, 4, 2, "spoon", 3 }
    local t2 = { 78, 45, 3, "table", 2, 8, 235, 342 }

    local t12 = utility.concat(t1, t2)
    assert(#t1 == 5 and #t2 == 8)    -- check that t1 and t2 were not modified
    for k, v in pairs(t12) do
        assert(check[k] == v, ("test failed for element %s at key %s, got %s"):format(check[k], k, v))
    end

    local t10 = utility.concat(t1, {})
    for k, v in pairs(t10) do
        assert(t1[k] == v, ("test failed for element %s at key %s, got %s"):format(t1[k], k, v))
    end
end

function main()
    test_repeat({ 1, 4, 2 }, 0)
    test_repeat({ 1, 4, 2 }, 1)
    test_repeat({ 1, 4, 2 }, 4)
    test_repeat(3, 0)
    test_repeat(3, 1)
    test_repeat("potential", 2)
    test_repeat(nil, 0) -- nil can't be repeated!

    test_concat()
end
