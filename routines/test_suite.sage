# To run this test suite, make sure to have sage-preparsed routine.sage first,
# so that we can import the defined methods.

load("routine.sage")
load("sageM2.sage")

from scipy.optimize import *
import numpy as np

def test0():
    """ Build the lp poly centered at zero and test that the unit vectors are contained
        in the lp ball.   
    """
    R = ZZ
    S = PolynomialRing(R, "x", 3)
    f = S.gens()[0]^2 + S.gens()[1]^2 + S.gens()[2]^2 + 10
    #f = shifted_lp_poly(R, 4, [0,0,0])

    assert(eval_poly(f, [0,0,0])==10)
    assert(eval_poly(f, [1,1,2])==16)

def test1():
    """ Test that the Lp poly is the one we expect it to be in specific cases:
    """
    R = ZZ
    S = PolynomialRing(R, "x", 3)
    f = 1 - (S.gens()[0] - 5)^2 - S.gens()[1]^2 - S.gens()[2]^2
    assert(shifted_lp_poly(S, 2, [5,0,0]) == f)

def test2():
    """ Test that if we evaluate the deformed poly at 0 we get the original poly.
    """
    R = ZZ
    S = PolynomialRing(R, "x", 3)
    lp = shifted_lp_poly(S, 4, [1,1,1])
    def_lp = deformed_poly(lp, "t")
    assert(lp == eval_poly(def_lp, [0]))

def test3():
    """ Test elimination function. Check that the projection of x^2 = y to the x axis gives everything.
    """
    S = PolynomialRing(QQ, "x", 3)
    I = S.ideal([S.gens()[0]^2-S.gens()[1], S.gens()[2]-S.gens()[0]])
    eliminatedIdeal = eliminateM2(I, "x0,x1")
    assert(len(eliminatedIdeal.gens()) == 1 and eliminatedIdeal.gens()[0] == 0)

def test4():
    """ Test the construction of the deformed lp_ball intersection polynomial.
    """
    n = 2
    p = 4
    mus = [[0,0], [1,0], [0,1]]
    R = PolynomialRing(QQ, "x", 2)
    x = R.gens()

    fs = [shifted_lp_poly(R, p, mu) for mu in mus]
    
    deformed_prod = deformed_product(fs, "t")

    Rt = deformed_prod.parent()
    t = Rt.gens()[0]

    assert(prod([1- sum([(x[i] - mu[i])^p for i in range(0,n)]) for mu in mus]) - t == deformed_prod)
    
def test5():
    """ Test the partial evaluation map.
    """
    R = PolynomialRing(QQ, "x", 3)
    x = R.gens()

    f = (x[0]-2)^2*x[1] + x[2]
    var_value_pairs = {x[0]:1, x[2]:5}
    
    evaluation = partial_eval_poly(f, var_value_pairs)

    assert(evaluation == x[1] + 5)
    assert(len(evaluation.parent().gens()) == 1)

def test6():
    """ Test identifying the real roots.
    """
    R = QQ["x"]
    x = R.gens()[0]

    f = x^3 + QQ(2^(-15))

    prec_bits = 1
    real_roots = identify_real_roots(f,prec_bits)
    assert(len(real_roots) == 3) # Since precision is too low.

    prec_bits = 10 # Might be interesting to try this out with more fine grained precision.
    real_roots = identify_real_roots(f,prec_bits)
    assert(len(real_roots) == 1)

def test7():
    """ Test the 1 dimensional volume.

    TODO: Compare with the exact value in the 2 dimensional case.
    """
    n = 2
    p = 4
    mus = [[0,0], [1,0], [0,1]]
    R = PolynomialRing(QQ, "x", 2)
    x = R.gens()

    fs = [shifted_lp_poly(R, p, mu) for mu in mus]

    vol = get_1_dim_volume(fs, var_value_pairs={x[0]:QQ(0.5)}, def_value=QQ(0.01), prec=20)

    assert(vol > 0 and vol < sqrt(2)) # Check the obvious bounds for the line width could be.

def test8():
    """ Test branch point computation.
    """
    R = PolynomialRing(QQ, "x", 2)
    x = R.gens()
    f = x[0]^2 + x[1]^2 - 1
    proj_var = x[1]
    # Projection onto x[1] axis, with branchpoints (0, +- 1)
    bpoints = branch_points(f, proj_var) 

    assert(len(bpoints) == 2)
    assert(all([point[x[0]]==0 and ( point[x[1]] == 1 or point[x[1]] == -1) for point in bpoints]))

def test9():
    """ Test the computation of real branch points of the deformed intersection. 
    Should compute the real branch points, or rather their value, when the branchpoint lies
    within the intersection of all lp balls.

    # TODO: No assertion so far, but to check that it runs without errors.
    """
    n = 2
    p = 4

    R = PolynomialRing(QQ, "x", n)
    mus = [[QQ(0), QQ(0)], [QQ(1),QQ(0)], [QQ(0), QQ(1)]]
    fs = [shifted_lp_poly(R, p, mu) for mu in mus]

    def_value = QQ(0.1)
    proj_var = R.gens()[1]

    f = eval_poly(deformed_product(fs), [def_value])

    # Below doesn't really terminate.

    # crit_val = project_deformed_intersection(fs, def_value, proj_var, var_value_pairs={})

    # print(crit_val)

    # RMK: 2025-10-18: Doesn't terminate in reasonable time. Should solve just numerically, with fixed precision. (the precision we need.)
    # Would be much better to do this via something similar to 

    # Numerical instead:
    # --> For example with Lasser's "Hierarchy of relaxation"

    # Without thinking, too much, start with the scipy optimizer first.
    # So we want to optimize proj_var given that f >= 0 and f_1, ..., f_m >=0. I.e. this way we stay in the intersection.

    # Can provide all the non-linear constraints:
    symbolic_functions = [poly.base_extend(SR) for poly in [f]+fs]
    constraints = [NonlinearConstraint(lambda x, y=fnc:y(x0=x[0], x1=x[1]), lb=0, ub=np.inf, jac=lambda x,y=fnc:np.array([y.gradient()[0](x0=x[0],x1=x[1]),y.gradient()[1](x0=x[0],x1=x[1])])) for fnc in symbolic_functions]

    start_point = (0.5,0.5) # Should already be in the semi-algebraic set!
    #obj_fun = lambda x : x[0]

    # And use BFGS, though need to add the correct gradient as information.

    res_min = minimize(lambda x : x[0], start_point, constraints=constraints)
    # And if we want to maximize, we simply take "-minimize(-obj_fun,...)"
    # print(res_min.success, res_min.fun, res_min.message)

    res_max = minimize(lambda x : -x[0], start_point, constraints=constraints)
    # print(res_max.success, res_max.fun, res_max.message)

    # SOMETHING BREAKS ABOVE IF I LET BOTH MINIMIZERS RUN...

    # obj_fun = lambda x : -x[0]
    # res_max = minimize(obj_fun, start_point, constraints=[NonlinearConstraint(lambda x, y=fnc:y(x0=x[0], x1=x[1]), lb=0, ub=np.inf) for fnc in symbolic_functions])
    # print(res_max)
    # Issue: Will have to sample a point in the region of interest.


    # Should talk with someone who has more experience: Lorenzo Baldi !!
    # (This would be much more effective than trying to set something up myself.)

    # TODO: Next: Reach out to Lorenzo. And besides that, implement the other steps involving the ore_algebra and solving.
    # Then run it step by step through our basic example that had worked out correctly.


def test10():
    """ Test the construction of the integrand in a basic example.
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
    assert(A == (-4*x[1]^4)/(-x[0]^4 - x[1]^4 + 1))

def test11():
    """ Test the construction of the rational Weyl algebra.
    """
    S = PolynomialRing(QQ, "x", 2)
    x = S.gens()
    W = rational_weyl_algebra(S)
    d = W.gens()

    assert(d[0] * x[0] == x[0]*d[0] + 1)
    # print(W)
    # print(d[0] * x[0])

def test12():
    """ Test the get Picard Fuchs:
    """
    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [shifted_lp_poly(S, p, [0,0])]
    deform_value = 0
    var_value_pairs = {}
    proj_var = S.gens()[0]

    
    W = rational_weyl_algebra(S)
    d = W.gens()

    P = get_picard_fuchs(fs, deform_value, var_value_pairs, proj_var, strategy=None)
    
    # Computed this before:
    assert(P == (-x[0]^4 + 1)*d[0] + (x[0]^3))

if __name__ == "__main__":
    test0()
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
    test9()
    test10()
    test11()
    test12()