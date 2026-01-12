# This is based on the example code provided in HypersurfaceRegions.jl

using HypersurfaceRegions
@var x0 x1 x2 x3;
fdef = x0^4 + 2*x0^2*x1^2 + x1^4 + 2*x0^2*x2^2 + 2*x1^2*x2^2 + x2^4 + 2*x0^2*x3^2 + 2*x1^2*x3^2 + 2*x2^2*x3^2 + x3^4 - 2*x0^3 - 2*x0*x1^2 - 2*x0*x2^2 - 2*x0*x3^2 - x0^2 - x1^2 - x2^2 - x3^2 + 2*x0 - 1/1000;
critValPoly = x0^6 - 3*x0^5 + 1251/1000*x0^4 + 1249/500*x0^3 - 563/250*x0^2 + 503/1000*x0 - 251/1000000;
system = [fdef;critValPoly];
regs = regions(system);
println("[")
for region in regs.region_list
  println(string(region.critical_points[1]) * ",")
end
println("]")

