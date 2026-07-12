# Interface to HypersurfaceRegions.jl
# This is meant to be run as a script.
# It loads the required operators from file and saves the output again.

using HypersurfaceRegions
using JSON


infile = ARGS[1];
outfile = ARGS[2];

# Load input
data = JSON.parsefile(infile);


# Create each variable individually
for var_name in data["variables"]
    eval(Meta.parse("@var " * var_name))
end


# Load system from list of strings
system = [eval(Meta.parse(poly_str)) for poly_str in data["system"]];

# Compute regions
regs = regions(system);

# Write output
open(outfile, "w") do io
    redirect_stdout(io) do 
        println([region.critical_points[1] for region in regs.region_list])
    end
end
