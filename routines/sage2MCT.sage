# Author: Nicolas Weiss
# Goal: Provide an "interface" to the multivariate creative telescoping implementation by Hadrien Brochet
#       -> [Done] Generate Julia code 
#       -> [TODO] Run and read the result
# 
#   Link to Hadrien Brochet's Julia package:


#
# Remark: Currently this needs to be run from within the src folder. Or pasted into a sage idle opened there.
#



######################

# With installation:
# import exact_convex_volumes

# Without installation:
from context import exact_convex_volumes

# In any case:
from exact_convex_volumes.volumes import *
from exact_convex_volumes.tools import *
from exact_convex_volumes.msolve_interface import *

######################

def get_x_from_Dx(D):
    """ D is a generator in the rational Weyl algebra C[x0...xn]
    returns the xi, such that D*xi = xi*D + 1.
    """
    candidates = [x for x in D.parent().base_ring().gens() if D * x == x * D + 1]
    assert(len(candidates) == 1)

    return candidates[0]
    
def get_Dx_from_x(x, W):
    """ For x and a corresponding rational Weyl algebra over C[x0,...,xn]
    returns the generator D of W satisfying D*x = x*D + 1.
    """ 
    candidates = [D for D in W.gens() if D * x == x * D + 1]
    assert(len(candidates) == 1)

    return candidates[0]

def define_WA_code(A, proj_var):
    """ Given the rational function A, it writes the code to set up
    the Weyl algebra 

        QQ(t)[x]<dt, dx>
    
    where proj_var = t, and x stands for the remaining variables in A.
    """

    W = rational_weyl_algebra(A.parent().ring())

    # Identify Dt:
    
    Dt = get_Dx_from_x(proj_var, W)

    # Set up the order
    order = "lex " + str(Dt) + " > grevlex " + " ".join([str(get_x_from_Dx(D)) for D in W.gens() if not D == Dt]) + " " + " ".join([str(D) for D in W.gens() if not D == Dt])

    print('W = OreAlg(order = "' + order + '",ratdiffvars=(["' + str(proj_var) + '"],["' + str(Dt) + '"]) , poldiffvars=( [' + ", ".join([ '"' + str(get_x_from_Dx(D)) + '"' for D in W.gens() if not D == Dt]) + '],[' + ", ".join([ '"' + str(D) + '"' for D in W.gens() if not D == Dt]) + '])' + ')')

def set_up_annihilator_code(A):
    """ Given the rational function A(t,x), it writes the code to set up
    a definite ideal that annihilates A.
    """

    # TODO: Consider factoring out the numerator!

    W = rational_weyl_algebra(A.parent().ring())

    annA = W.ideal([A*D-D(A) for D in W.gens()])

    print('ann = [' + ", ".join(['parse_OrePoly("' + str(P) + '",W)' for P in annA.gens()]) + ']')


def set_up_MCT_julia_code(A, proj_var):
    """ Given the rational function A and a chosen projection variable, it
    sets up the julia code for the creative telescoping.
    """

    num = A.numerator() # Factor out the numerator.

    define_WA_code(A, proj_var)
    set_up_annihilator_code(A/num)

    print('init = weyl_closure_init(W)')
    print('gb = weyl_closure(ann,W,init)')

    print('LDE = MCT(parse_OrePoly("' + str(num) + '",W), gb, W)')
    

def example1():
    """ Code for 1 L4 ball centered at the origin.
    """
    n = 2
    p = 4
    R = PolynomialRing(QQ, "x", n)
    x = R.gens()

    fs = [shifted_lp_poly(R, p, [0,0])]
    
    def_value = 0
    var_value_pairs = {}
    proj_var = x[0]
    A = construct_integrand(fs, def_value, var_value_pairs,proj_var)

    set_up_MCT_julia_code(A,proj_var)

def example2():
    """ Code for the intersection of 2 L4 balls centered at [0,0] and [1,0]
    """
    n = 2
    p = 4
    R = PolynomialRing(QQ, "x", n)
    x = R.gens()

    fs = [shifted_lp_poly(R, p, [0,0]), shifted_lp_poly(R, p, [1,0])]
    
    def_var_name = "t"
    At = construct_integrand_t(fs, def_var_name)

    proj_var = At.parent().ring()(def_var_name)

    set_up_MCT_julia_code(At,proj_var)


def example3():
    """ Code for the intersection of 3 L4 balls centered at [0,0], [1,0], and [0,1].
    """
    n = 2
    p = 4
    R = PolynomialRing(QQ, "x", n)
    x = R.gens()

    fs = [shifted_lp_poly(R, p, [0,0]), shifted_lp_poly(R, p, [1,0]), shifted_lp_poly(R, p, [0,1])]
    
    def_var_name = "t"
    At = construct_integrand_t(fs, def_var_name)

    proj_var = At.parent().ring()(def_var_name)

    set_up_MCT_julia_code(At, proj_var)


def twoL2inR4():
    """ Intersection of two L2 balls in R4. """
    n = 4
    p = 2
    R = PolynomialRing(QQ, "x", n)
    x = R.gens()

    fs = [shifted_lp_poly(R, p, [0,0,0,0]), shifted_lp_poly(R, p, [1,0,0,0])]
    
    def_var_name = "t"
    At = construct_integrand_t(fs, def_var_name)

    proj_var = At.parent().ring()(def_var_name)

    set_up_MCT_julia_code(At, proj_var)

    # Similarly for the deformations, set up the Julia code:
    initial_points_t = [1/50, 1/40, 3/100]

    proj_var = x[3]
    for def_value in initial_points_t:
        A = construct_integrand(fs, def_value, {}, proj_var)
        print()
        set_up_MCT_julia_code(A, proj_var)
        print()

def twoL4inR4():
    n = 4
    p = 4
    R = PolynomialRing(QQ, "x", n)
    x = R.gens()

    fs = [shifted_lp_poly(R, p, [0,0,0,0]), shifted_lp_poly(R, p, [1,0,0,0])]
    
    def_var_name = "t"
    At = construct_integrand_t(fs, def_var_name)

    proj_var = At.parent().ring()(def_var_name)

    set_up_MCT_julia_code(At, proj_var)





if __name__ == "__main__":
    twoL4inR4()