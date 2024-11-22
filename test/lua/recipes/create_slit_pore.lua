-- The next line adds the example tree of the source directory to the Lua
-- search path. The regular expression is used to construct a path relative to
-- the current script.
package.path = arg[0]:match("@?(.*/)") .. "../../../examples/?.lua;" .. package.path
require("wetting/wetting_equilibration")   -- defines setup()

function main()
    local args = {
        box_length = 50, pore_width = 20
      , wall_epsilon = 1, wetting = 1.5, smoothing = 0.005
      , density = 0.75, temperature  = 0.9
    }

    -- call function from example script
    setup(args)
end
