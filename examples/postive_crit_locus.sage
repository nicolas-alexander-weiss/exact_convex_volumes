"""
@author: Nicolas Weiss
@date: 2026-07-12

The following is an example case, where the critical locus of the deformed intersection is positive dimensional.
"""

from exact_convex_volumes.volumes import *
from exact_convex_volumes.tools import *

R =  PolynomialRing(QQ, "x", 4)
x = R.gens()

fs = [tools.shifted_lp_poly(R, 2, [0,0,0,0]), tools.shifted_lp_poly(R, 2, [1,0,0,0])]

tdef = 1/100
x1val = 1/100
proj_var = x[0]

# This will use the interface to HypersurfaceRegions.jl to compute the critical values corresponding to the deformed intersection.
project_deformed_intersection(fs, def_value=tdef, var_value_pairs={x[Integer(1)]:x1val}, proj_var=proj_var, prec=Integer(200))