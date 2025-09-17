# Created: 2025-09-17
# Author: Nicolas Weiss

# Instructions to run this: 
# --> I recommend to install the julia add-on for VSCode.
# --> Then the below script can be executed line by line in a julia REPL 
#       using the command "Julia: Send current line or selection to REPL".
# --> You can set a shortcut for this, e.g. ALT + ENTER.

using HomotopyContinuation
const HC = HomotopyContinuation 

# Writing an example of a deformed system, 
# intersection of 3 lp balls in R_2

const n = 2
const p = 4

mus = [[0,0], [1,1], [1,0]]

@var x[1:n] t[1:1]

# The intersection then defined by fs[i] \geq 0. for all i:
fs = [1 - ((x[1]-mu[1])^p + (x[2]-mu[2])^p) for mu in mus]

F = prod(fs)
Ft = F - t[1]

# We now consider the "deformed" setting of t[1] = epsilon

epsilon = 0.001
Feps = Ft(t[1] => epsilon)

# The projection onto the $x[2]$ axis has branchpoints defined by the vanishing of
criticalPointsSystem = [Feps, differentiate(Feps, x[1])]

# The real valued critical points are:
res = HC.solve(criticalPointsSystem) # This might take a bit.
criticalPoints = HC.real_solutions(res)

# TODO: Add a plot of the deformed setup with marked critical points.

# Let us now continue the critical points back to eps = 0.

criticalPointsSystemParametric = [F, differentiate(F, x[1])]

HCSystem = System(criticalPointsSystemParametric, variables = x[1:2], parameters = t[1:1]) 
res = HC.solve(HCSystem, criticalPoints; start_parameters = [epsilon], target_parameters = [0]) 

# Note that the continuation of these solutions might well move to "singular" points of the system for t=0.
# We are only interested in the endpoints:

# Really only want to get a real number for the further computations. 
# TODO: See if we can also proceed with the "complex values" in a reasonable way.
HCcriticalPoints = [[result.solution[1].re, result.solution[2].re] for result in res]

# Let us now identify those, which upto accuracy acc are lying on the boundary 

# of our set of interest \bicap \{f_i \leq 0 \}:

acc0 = 0.0001 # might have to modify adaptively to identify correctly the points in the branch ideal.

# Need to check that point is close to satisfying all fs[i](pt) is greater equal 0 (upto acc):

on_boundary(point, acc) = all([f(x[1] => point[1], x[2] => point[2]) > -acc for f in fs])

crit_points_on_boundary(points, acc) = [pt for pt in points if on_boundary(pt,acc)==true]

HConBoundary = crit_points_on_boundary(HCcriticalPoints, acc0)

# TODO: plot them and check whether it makes sense.

xs = [point[1] for point in HConBoundary]
ys = [point[2] for point in HConBoundary]
fig = plot();
scatter!(fig, xs, ys, markersize=10, markercolor=:blue);

display(fig)

r = range(-1,1, length=100)
Xs = [X for X in r for Y in r]
Ys = [Y for X in r for Y in r]
Zs = [[F(x[1]=> X, x[2]=>Y) for Y in r] for X in r]
contour(r, r, Zs, linecolor=:red, linewidth=2)

# TODO: Plot doesn't work :()

# TODO: Now identify from this which are the boundary points of the original setting. 
#           Better: Check is real and on boundary in one go, to keep a list of the same length before and after HomotopyContinuation.