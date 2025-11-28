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
    print("[TEST0]")
    R = ZZ
    S = PolynomialRing(R, "x", 3)
    f = S.gens()[0]^2 + S.gens()[1]^2 + S.gens()[2]^2 + 10
    #f = shifted_lp_poly(R, 4, [0,0,0])

    assert(eval_poly(f, [0,0,0])==10)
    assert(eval_poly(f, [1,1,2])==16)

def test1():
    """ Test that the Lp poly is the one we expect it to be in specific cases:
    """
    print("[TEST1]")
    R = ZZ
    S = PolynomialRing(R, "x", 3)
    f = 1 - (S.gens()[0] - 5)^2 - S.gens()[1]^2 - S.gens()[2]^2
    assert(shifted_lp_poly(S, 2, [5,0,0]) == f)

def test2():
    """ Test that if we evaluate the deformed poly at 0 we get the original poly.
    """
    print("[TEST2]")
    R = ZZ
    S = PolynomialRing(R, "x", 3)
    lp = shifted_lp_poly(S, 4, [1,1,1])
    def_lp = deformed_poly(lp, "t")
    assert(lp == eval_poly(def_lp, [0]))

def test3():
    """ Test elimination function. Check that the projection of x^2 = y to the x axis gives everything.
    """
    print("[TEST3]")
    S = PolynomialRing(QQ, "x", 3)
    I = S.ideal([S.gens()[0]^2-S.gens()[1], S.gens()[2]-S.gens()[0]])
    eliminatedIdeal = eliminateM2(I, "x0,x1")
    assert(len(eliminatedIdeal.gens()) == 1 and eliminatedIdeal.gens()[0] == 0)

def test4():
    """ Test the construction of the deformed lp_ball intersection polynomial.
    """
    print("[TEST4]")
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
    print("[TEST5]")
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
    print("[TEST6]")
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
    print("[TEST7]")
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
    print("[TEST8]")
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

    REMARK
    ------
        [TODO] This test is no longer needed. Recycle slot for another test.
    """
    print("[TEST9]")
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
    print("[TEST10]")
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
    print("[TEST11]")
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
    print("[TEST12]")
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

def test13():
    """ Test construct_integrand_t()
    """
    print("[TEST13]")
    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [shifted_lp_poly(S, p, [0,0])]

    def_var_name = "t"

    At = construct_integrand_t(fs, def_var_name)

    assert(At == At.parent()("(-4*x0^4)/(-x0^4 - x1^4 - t + 1)"))

def test14():
    """ Test get_picard_fuchs_t()
    # TODO: Add a condition here. So far only a termination test of the code.
    """
    print("[TEST14]")

    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [shifted_lp_poly(S, p, [0,0])]

    def_var_name = "t"
    ft = deformed_product(fs, def_var_name)

    ft_flattened = ft.parent().flattening_morphism()(ft) 

    At = construct_integrand_t(fs, def_var_name, strategy=None)

    # print("Construct the integrand:")
    # print(At)
    # print(At.parent())

    # print("\nConstruct the rational Weyl algebra:")
    Wt = rational_weyl_algebra(At.parent().ring())

    # print(Wt)

    annAt = Wt.ideal([At*D-D(At) for D in Wt.gens()]) # construct the annihilating ideal
    # print(annAt)
    intIdeal_t = creative_telescoping(annAt, At.parent().ring()(def_var_name))

    # print(intIdeal_t)
    
def test15():
    """ Set up an example where the CT might not terminate (i.e. for comparison with e.g. Magma or other CT systems.)

    """
    if ONLY_FAST_TESTS and not 15 in MUST_DO:
        print("[TEST15] -- Skipped slow test")
        return
    print("[TEST15] Long test!")

    n = 2
    p = 4

    R = PolynomialRing(QQ, "x", n)
    mus = [[QQ(0), QQ(0)], [QQ(1),QQ(0)], [QQ(0), QQ(1)]]
    fs = [shifted_lp_poly(R, p, mu) for mu in mus]

    deform_value = QQ(0.1)
    proj_var = R.gens()[1]

    var_value_pairs = {}

    print("Computing Picard Fuchs operator for the deformed slice. Might take a moment.")
    P = get_picard_fuchs(fs, deform_value, var_value_pairs, proj_var, strategy=None)
    # Finishes, but takes a while!

def test16():
    """ A test of the computation of the "inside branch_points"
    """
    print("[TEST16]")

    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    mus = [(1,0), (1,1), (0,1)]
    prec_bits = 400

    fs = [shifted_lp_poly(S, p, mu) for mu in mus]

    def_value = 0.1

    fdeft = deformed_product(fs, "t")
    fdef = eval_poly(fdeft, [def_value])

    proj_var  = x[1]

    # Classical way of solving for branch points using msolve interface:

    # branch_system = branch_points_system(fdef, proj_var)
    # bpts = variety_msolve(branch_system, prec_bits)
    # print(bpts)
    # plt = contour_plot(fdef, (0,2), (0,2), contours=(0,), fill=False)
    # plt.save_image('Contour.png')

    target_accuracy = 10^(-100)
    target_branch_points = [
        {   x[0]: 0.520089189436583572099386548315109453709892711432221927121966608506140183711041309893220576274838347848254425230357792400,
            x[1]: 0.115612795644764287999445860503433810430423048493184113046313751725931712460186211614019685391545667827399700045578468316},
        
        {   x[0]: 0.65986137169719393343842666348786166256121414985538656002676237617740667458440343820447246037894511649916856006920436507,
            x[1]: 0.96342885339173973353718256114649404783334788940698883387078736518040784929246713957594138522332784033658300019879997926},
    ]

    inside_bpts = get_inside_branch_points(fs, def_value, proj_var, {}, prec_bits)

    RBF = RealBallField(prec_bits - 2)

    for pt in inside_bpts:
        # Check that each of the computed inside branch_points is found
        assert(len([bp for bp in target_branch_points if  abs(RBF(bp[x[0]])-RBF(pt[x[0]])) < target_accuracy and  abs(RBF(bp[x[1]])-RBF(pt[x[1]])) < target_accuracy]) == 1)
    
    # Also check that its all, i.e. two.
    assert(len(inside_bpts) == 2)


def test17():
    """ Basic Volume computation: Volume of deformed slice for mu = [0,0] and [0,1]

    # TODO: Numerically compute the volume by sampling to verify.
    --> So far mainly a test that it runs through.

    """ 

    if ONLY_FAST_TESTS and not 17 in MUST_DO:
        print("[TEST17] -- Skipped slow test")
        return
    print("[TEST17]")

    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [shifted_lp_poly(S, p, [0,0]), shifted_lp_poly(S, p, [1,0])]

    t0 = QQ(0.3)

    prec = 200
    volume = volume1(fs, t0, {}, prec)
    target_vol = 1.16484606433490063654808600776036941031049501636400071712
    CBF = ComplexBallField(prec - 2)

    assert(abs(CBF(volume) - CBF(target_vol)) < 10^(-10)) # Just a crude comparison
    

def test18():
    """ Bernd's hour volume example: n=2, p=4, mus = [(0,0), (1,0)]
    """
    
    if ONLY_FAST_TESTS and not 18 in MUST_DO:
        print("[TEST18] -- Skipped slow test")
        return
    
    print("[TEST18]")

    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    mus = [(0,0), (1,0)]
    prec_bits = 400

    fs = [shifted_lp_poly(S, p, mu) for mu in mus]

    vol = volume2(fs, prec=prec_bits)

    CBF = ComplexBallField(prec_bits - 2)
    target_vol = 1.714482859044855236241617456389863266850907891292196002387016141008363973273758352945472686246402767128126792
    target_accuracy = 10^(-100)

    assert(abs(CBF(vol) - CBF(target_vol)) < target_accuracy)

def test19():
    """ The volume is translation invariant
    """
    print("[TEST19] Not yet implemented!")
    pass

def test20():
    """ Computation of volume by rescaling? --> Adapt to non-unit balls.
    """
    print("[TEST20] Not yet implemented!")
    pass

def test21():
    """ Compute the PicardFuchsT for deformation with 2 pts. (For time taking.)
    """

    if ONLY_FAST_TESTS and not 21 in MUST_DO:
        print("[TEST21] -- Skipped slow test")
        return
    
    print("[TEST21] -- Slow test!")

    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    mus = [(0,0), (1,0)]
    prec_bits = 400

    fs = [shifted_lp_poly(S, p, mu) for mu in mus]

    Pt = get_picard_fuchs_t(fs)

def test22():
    """ Test that volume1 also runs correctly, when giving it only 1 polynomial, without deforming. 
    """
    print("[TEST22]")

    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [shifted_lp_poly(S, p, [0,0])]

    prec_bits = 400
    volume = volume1(fs, 0, {}, prec=400)

    # TODO Add assertion on the volume.

def test23():
    """ Test that volume2 also runs correctly when only
    providing 1 polynomial is input.
    """
    print("[TEST23]")

    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [shifted_lp_poly(S, p, [0,0])]

    prec = 400

    vol1 = volume1(fs, 0, {}, prec)
    vol2 = volume2(fs, prec)

    target_accuracy = 10^(-100)

    CBF = ComplexBallField(prec - 2)
    assert(abs(CBF(vol1) - CBF(vol2)) < target_accuracy)

def test24():
    """ Volume of L2 ball in RR^4.

    Remark: In the implementation of Mark Mezzaroba, this is mentioned to fail due to the CT algorithm.
    https://src.koda.cnrs.fr/marc.mezzarobba.3/volumes/-/blob/main/volume.py?ref_type=heads
    """
    print("[TEST24]")

    n = 4
    p = 2
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [shifted_lp_poly(S, p, [0,0,0,0])]

    prec = 400
    ##
    P = get_picard_fuchs(fs, 0, {}, x[0])
    # print("The Picard Fuchs operator: {}".format(P))
    assert(P == 1)



if __name__ == "__main__":
    global ONLY_FAST_TESTS, MUST_DO
    ONLY_FAST_TESTS = True
    MUST_DO = []


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
    test13()
    test14()
    test15()
    test16()
    test17()
    test18()
    test19()
    test20()
    test21()
    test22()
    test23()
    test24()