# This is meant to be run as a script.
# It loads the required operators from file and saves the output again.

using MultivariateCreativeTelescoping
using JSON


infile = ARGS[1]
outfile = ARGS[2]

# Load input
data = JSON.parsefile(infile)


W = OreAlg(
    order = data["order"],
    ratdiffvars = (String.(data["ratdiffvars"][1]), String.(data["ratdiffvars"][2])),
    poldiffvars = (String.(data["poldiffvars"][1]), String.(data["poldiffvars"][2]))
)

# Set up annihilator of 1/Ft
ann = [parse_OrePoly(String(s), W) for s in data["ann"]] 
numerator = parse_OrePoly(String(data["numerator"]), W)

# Compute the PF operator
init = weyl_closure_init(W)
gb = weyl_closure(ann, W, init)
LDE = MCT(numerator, gb, W)

# Write output
open(outfile, "w") do io
    redirect_stdout(io) do 
        prettyprint(LDE, W)
    end
end
