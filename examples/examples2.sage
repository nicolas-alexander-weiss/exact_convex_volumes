#
# Add path to the package.
import sys
sys.path.append("../src/volumes")

from volumes import *
from tools import *



def last_variable_proj_var_strategy(fs, deform_value, var_value_pairs):
    evaluated_poly = partial_eval_poly(eval_poly(deformed_product(fs), [deform_value]), var_value_pairs)

    return evaluated_poly.parent().gens()[-1]

def example2():

    n = 3
    p = 2

    prec = 200

    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    mus = [[0,0,0], [1,0,0]]

    fs = [shifted_lp_poly(S, p, mu) for mu in mus]

    strategy = {"proj_var":last_variable_proj_var_strategy}

    #slice_vols = dict((val, volume1(fs, deform_value, var_value_pairs={x[1]:val}, prec=prec, debug_level=3, strategy=strategy)) for val in initial_points)
    vol = volume2(fs, prec, strategy=None)
    print("The volumes of the slices: {}".format(vol))
