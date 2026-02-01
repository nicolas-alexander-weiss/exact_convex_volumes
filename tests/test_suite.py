# Run with
#           sage --python test_suite.py
#

# Load the package from context if not globally installed
from context import exact_convex_volumes

# Load the sub-modules of the volumes package
from exact_convex_volumes import volumes, tools
from exact_convex_volumes import m2_interface as m2
from exact_convex_volumes import msolve_interface as msolve

# Sage imports
from sage.all import (

    ZZ, QQ, PolynomialRing,
    RealBallField, ComplexBallField,
    prod, sqrt,
    diff,
)

# Python imports
import numpy as np


#
# msolve interface tests:
#

def msolve_test1():
    """ Compute the intersection of the unit ball and the x-axis.
    """
    pass

msolve_interface_tests = [msolve_test1]


#
# m2_interface tests:
#

def m2_test_1():
    """ Test elimination function. Check that the projection of x^2 = y to the x axis gives everything.
    """
    S = PolynomialRing(QQ, "x", 3)
    I = S.ideal([S.gens()[0]**2-S.gens()[1], S.gens()[2]-S.gens()[0]])
    eliminatedIdeal = m2.eliminate(I, "x0,x1")
    assert(len(eliminatedIdeal.gens()) == 1 and eliminatedIdeal.gens()[0] == 0)

m2_interface_tests = [m2_test_1]

#
# Tools tests:
#

def tools_test_1():
    """ Build the lp poly centered at zero and test that the unit vectors are contained
        in the lp ball.   
    """
    R = ZZ
    S = PolynomialRing(R, "x", 3)
    f = S.gens()[0]**2 + S.gens()[1]**2 + S.gens()[2]**2 + 10

    pt = [0,0,0]
    val = tools.eval_poly(f, pt)
    expected_val = 10
    if not val == expected_val:
        raise ValueError(f"Expected that {f} evaluates at {pt} to {expected_val}, but got {val}.")
    
    pt = [1,1,2]
    val = tools.eval_poly(f, pt)
    expected_val = 16
    if not val == expected_val:
        raise ValueError(f"Expected that {f} evaluates at {pt} to {expected_val}, but got {val}.")

def tools_test_2():
    """ Test that the Lp poly is the one we expect it to be in specific cases:
    """
    R = ZZ
    S = PolynomialRing(R, "x", 3)
    f = 1 - (S.gens()[0] - 5)**2 - S.gens()[1]**2 - S.gens()[2]**2
    assert(tools.shifted_lp_poly(S, 2, [5,0,0]) == f)

def tools_test_3():
    """ Test that if we evaluate the deformed poly at 0 we get the original poly.
    """
    R = ZZ
    S = PolynomialRing(R, "x", 3)
    lp = tools.shifted_lp_poly(S, 4, [1,1,1])
    def_lp = tools.deformed_poly(lp, "t")
    assert(lp == tools.eval_poly(def_lp, [0]))

def tools_test_4():
    """ Test the construction of the deformed lp_ball intersection polynomial.
    """
    n = 2
    p = 4
    mus = [[0,0], [1,0], [0,1]]
    R = PolynomialRing(QQ, "x", 2)
    x = R.gens()

    fs = [tools.shifted_lp_poly(R, p, mu) for mu in mus]
    
    deformed_prod = tools.deformed_product(fs, "t")

    Rt = deformed_prod.parent()
    t = Rt.gens()[0]

    assert(prod([1- sum([(x[i] - mu[i])**p for i in range(0,n)]) for mu in mus]) - t == deformed_prod)

def tools_test_5():
    """ Test the partial evaluation map.
    """
    R = PolynomialRing(QQ, "x", 3)
    x = R.gens()

    f = (x[0]-2)**2 * x[1] + x[2]
    var_value_pairs = {x[0]:1, x[2]:5}
    
    evaluation = tools.partial_eval_poly(f, var_value_pairs)

    assert(evaluation == x[1] + 5)
    assert(len(evaluation.parent().gens()) == 1)

tools_tests = [tools_test_1, tools_test_2, tools_test_3, tools_test_4, tools_test_5]



#
# Volume tests:
#

def volume_test_1():
    """ Test the 1 dimensional volume. 
    [TODO: Add more precise value comparison.]
    """
    n = 2
    p = 2
    mus = [[0,0],[1,0], [0,1]]
    R = PolynomialRing(QQ, "x", 2)
    x = R.gens()

    fs = [tools.shifted_lp_poly(R, p, mu) for mu in mus]

    vol = volumes.get_1_dim_volume(fs, var_value_pairs={x[0]:QQ("1/2")}, def_value=QQ("1/100"), prec=20)
    assert(vol > 0 and vol < sqrt(2)) # Check what would be the obvious bounds for the length of the line.


def volume_test_2():
    """ Test critical locus computation.
    """
    R = PolynomialRing(QQ, "x", 2)
    x = R.gens()
    f = x[0]**2 + x[1]**2 - 1
    proj_var = x[1]
    # Projection onto x[1] axis, with critical points (0, +- 1)
    crit_points = volumes.get_critical_points(f, proj_var) 

    assert(len(crit_points) == 2)
    assert(all([point[x[0]]==0 and ( point[x[1]] == 1 or point[x[1]] == -1) for point in crit_points]))

def volume_test_3():
    """
    Check that the critical points computation raises an error if the ideal is positive dimensional.
    """
    n = 4
    p = 2

    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    mus = [[0,0,0,0], [1,0,0,0]]

    fs = [tools.shifted_lp_poly(S, p, mu) for mu in mus]
    def_value = QQ("1/1000")
    F = tools.eval_poly(tools.deformed_product(fs, "t"), def_value)
    
    proj_var = x[0]

    try:
        crit_points = volumes.get_critical_points(F, proj_var) 
    except volumes.PosDimCritLocusError:
        crit_locus_pos_dim = True

    assert(crit_locus_pos_dim)


def volume_test_4():
    """ test_ the construction of the rational integrand in a basic example.
    """
    n = 2
    p = 4
    R = PolynomialRing(QQ, "x", n)
    x = R.gens()

    fs = [tools.shifted_lp_poly(R, p, [0,0])]
    
    def_value = 0
    var_value_pairs = {}
    proj_var = x[0]
    A = volumes.construct_integrand(fs, def_value, var_value_pairs, proj_var)

    assert(A == diff(fs[0], x[1]) * x[1] / fs[0])

def volume_test_5():
    """ Test the construction of the rational Weyl algebra.
    """
    S = PolynomialRing(QQ, "x", 2)
    x = S.gens()
    W = volumes.rational_weyl_algebra(S)
    d = W.gens()
    assert(d[0] * x[0] == x[0]*d[0] + 1)


def volume_test_6():
    """ Test the get Picard Fuchs:
    """
    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [tools.shifted_lp_poly(S, p, [0,0])]
    deform_value = 0
    var_value_pairs = {}
    proj_var = S.gens()[0]
    
    W = volumes.rational_weyl_algebra(S)
    d = W.gens()

    P = volumes.get_picard_fuchs(fs, deform_value, var_value_pairs, proj_var, strategy=None)
    
    # Computed this before:
    assert(P == (-x[0]**4 + 1) * d[0] + (x[0]**3))


def volume_test_7():
    """ Test construct_integrand_t()
    """
    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [tools.shifted_lp_poly(S, p, [0,0])]

    def_var_name = "t"

    At = volumes.construct_integrand_t(fs, def_var_name)

    assert(At == At.parent()("(-4*x0^4)/(-x0^4 - x1^4 - t + 1)"))


def volume_test_8():
    """ Test get_picard_fuchs_t()
    # TODO: Add a condition here. So far only a termination test of the code.
    """
    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)

    fs = [tools.shifted_lp_poly(S, p, [0,0])]

    def_var_name = "t"
   
    At = volumes.construct_integrand_t(fs, def_var_name, strategy=None)
    Wt = volumes.rational_weyl_algebra(At.parent().ring())

    annAt = Wt.ideal([At*D-D(At) for D in Wt.gens()]) # construct the annihilating ideal

    volumes.creative_telescoping(annAt, At.parent().ring()(def_var_name))


def volume_test_9():
    """ A test of the computation of the "inside branch_points"
    """
    n = 2
    p = 4
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    mus = [(1,0), (1,1), (0,1)]
    fs = [tools.shifted_lp_poly(S, p, mu) for mu in mus]

    def_value = QQ("1/10")
    proj_var  = x[1]

    prec_bits = 400
    RBF = RealBallField(prec_bits - 2)

    inside_crit_pts = volumes.get_inside_critical_points(fs, def_value, proj_var, {}, prec_bits)

    target_accuracy = RBF(10**(-100))

    target_crit_points = [
        {   x[0]: RBF("0.520089189436583572099386548315109453709892711432221927121966608506140183711041309893220576274838347848254425230357792400"),
            x[1]: RBF("0.115612795644764287999445860503433810430423048493184113046313751725931712460186211614019685391545667827399700045578468316")},
        
        {   x[0]: RBF("0.65986137169719393343842666348786166256121414985538656002676237617740667458440343820447246037894511649916856006920436507"),
            x[1]: RBF("0.96342885339173973353718256114649404783334788940698883387078736518040784929246713957594138522332784033658300019879997926")},
    ]

    for pt in inside_crit_pts:
        # Check that each of the computed inside_crit_points is found
        assert(len([cp for cp in target_crit_points if  abs(cp[x[0]]-pt[x[0]]) < target_accuracy and  abs(cp[x[1]]-pt[x[1]]) < target_accuracy]) == 1)
    
    # Also check that its all, i.e. two.
    assert(len(inside_crit_pts) == 2)

def volume_test_10():
    """ Chyzak's algorithmmight not result in the correct operators. For example if the rational function 
    decomposes as a sum of derivatives.

    Mentioned in:
    https://src.koda.cnrs.fr/marc.mezzarobba.3/volumes/-/blob/main/volume.py?ref_type=heads
    """
    n = 4
    p = 2
    S = PolynomialRing(QQ, "x", n)
    x = S.gens()

    fs = [tools.shifted_lp_poly(S, p, [0,0,0,0])]

    prec = 400
    ##
    P = volumes.get_picard_fuchs(fs, 0, {}, x[0])
    # print("The Picard Fuchs operator: {}".format(P))
    assert(P == 1)


volume_tests = [volume_test_1, volume_test_2, volume_test_3, volume_test_4, volume_test_5, volume_test_6, volume_test_7,
                volume_test_8, volume_test_9, volume_test_10]




#
# Long tests, i.e. integration tests of the various components.
# In particular, tests on volume1 and volume2.
#


def long_test_1():
    """ Compare the volume of the L2, L4 and L6 ball in RR^3 with the output of the computation.
    
    Note that the volume of an Lp ball of radius r in RR^n is described by the formula (see Wikipedia, https://en.wikipedia.org/wiki/Volume_of_an_n-ball)

    r^n * (2 Gamma(1/p + 1))^n / Gamma(n/p + 1) 

    """
    print("[Long test 1] Compare closed form for vol of L2, L4 and L6 ball in RR^3 with output of volume1.")

    prec_bits = 200

    n = 3
    p = 2
    r = 1

    S = PolynomialRing(QQ, "x", n)

    use_complex = True

    for p in [2,4,6]:
        vol_p_n_r_closed = volumes.vol_lp_ball_closed_formula(n,p,r, prec_bits, use_complex)

        fs = [tools.shifted_lp_poly(S, p, np.zeros((n,), int))]
        vol_p_n = volumes.volume1(fs, 0, {}, prec_bits)

        print(f"Computed volume: {vol_p_n}")
        print(f"Closed form volume: {vol_p_n_r_closed}")
        print(f"Difference: {vol_p_n - vol_p_n_r_closed}")

        # Both outputs will be in the complex ball field of the same precision. 
        # So we can compare if the mid point of the closed formula output lives in the output
        # of our computation:
        assert(vol_p_n.parent() == vol_p_n_r_closed.parent())
        assert(vol_p_n_r_closed.mid() in vol_p_n)


long_tests = [long_test_1]


    
# def test15():
#     """ Set up an example where the CT might not terminate (i.e. for comparison with e.g. Magma or other CT systems.)

#     """
#     if ONLY_FAST_TESTS and not 15 in MUST_DO:
#         print("[TEST15] -- Skipped slow test")
#         return
#     print("[TEST15] Long test!")

#     n = 2
#     p = 4

#     R = PolynomialRing(QQ, "x", n)
#     mus = [[QQ(0), QQ(0)], [QQ(1),QQ(0)], [QQ(0), QQ(1)]]
#     fs = [shifted_lp_poly(R, p, mu) for mu in mus]

#     deform_value = QQ(0.1)
#     proj_var = R.gens()[1]

#     var_value_pairs = {}

#     print("Computing Picard Fuchs operator for the deformed slice. Might take a moment.")
#     P = get_picard_fuchs(fs, deform_value, var_value_pairs, proj_var, strategy=None)
#     # Finishes, but takes a while!


# def test17():
#     """ Basic Volume computation: Volume of deformed slice for mu = [0,0] and [0,1]

#     # TODO: Numerically compute the volume by sampling to verify.
#     --> So far mainly a test that it runs through.

#     """ 

#     if ONLY_FAST_TESTS and not 17 in MUST_DO:
#         print("[TEST17] -- Skipped slow test")
#         return
#     print("[TEST17]")

#     n = 2
#     p = 4
#     S = PolynomialRing(QQ, "x", n)
#     x = S.gens()

#     fs = [shifted_lp_poly(S, p, [0,0]), shifted_lp_poly(S, p, [1,0])]

#     t0 = QQ(0.3)

#     prec = 200
#     volume = volume1(fs, t0, {}, prec)
#     target_vol = 1.16484606433490063654808600776036941031049501636400071712
#     CBF = ComplexBallField(prec - 2)

#     assert(abs(CBF(volume) - CBF(target_vol)) < 10^(-10)) # Just a crude comparison
    

# def test18():
#     """ Bernd's hour volume example: n=2, p=4, mus = [(0,0), (1,0)]
#     """
    
#     if ONLY_FAST_TESTS and not 18 in MUST_DO:
#         print("[TEST18] -- Skipped slow test")
#         return
    
#     print("[TEST18]")

#     n = 2
#     p = 4
#     S = PolynomialRing(QQ, "x", n)
#     x = S.gens()

#     mus = [(0,0), (1,0)]
#     prec_bits = 400

#     fs = [shifted_lp_poly(S, p, mu) for mu in mus]

#     vol = volume2(fs, prec=prec_bits)

#     CBF = ComplexBallField(prec_bits - 2)
#     target_vol = 1.714482859044855236241617456389863266850907891292196002387016141008363973273758352945472686246402767128126792
#     target_accuracy = 10^(-100)

#     assert(abs(CBF(vol) - CBF(target_vol)) < target_accuracy)

# def test19():
#     """ The volume is translation invariant
#     """
#     print("[TEST19] Not yet implemented!")
#     pass

# def test20():
#     """ Computation of volume by rescaling? --> Adapt to non-unit balls.
#     """
#     print("[TEST20] Not yet implemented!")
#     pass


# def test21():
#     """ Compute the PicardFuchsT for deformation with 2 pts. (For time taking.)
#     """

#     if ONLY_FAST_TESTS and not 21 in MUST_DO:
#         print("[TEST21] -- Skipped slow test")
#         return
    
#     print("[TEST21] -- Slow test!")

#     n = 2
#     p = 4
#     S = PolynomialRing(QQ, "x", n)
#     x = S.gens()

#     mus = [(0,0), (1,0)]
#     prec_bits = 400

#     fs = [shifted_lp_poly(S, p, mu) for mu in mus]

#     Pt = get_picard_fuchs_t(fs)


# def test22():
#     """ Test that volume1 also runs correctly, when giving it only 1 polynomial, without deforming. 
#     """
#     print("[TEST22]")

#     n = 2
#     p = 4
#     S = PolynomialRing(QQ, "x", n)
#     x = S.gens()

#     fs = [shifted_lp_poly(S, p, [0,0])]

#     prec_bits = 400
#     volume = volume1(fs, 0, {}, prec=400)
#     # TODO Add assertion on the volume.

# def test23():
#     """ Test that volume2 also runs correctly when only
#     providing 1 polynomial is input.
#     """
#     print("[TEST23]")

#     n = 2
#     p = 4
#     S = PolynomialRing(QQ, "x", n)
#     x = S.gens()

#     fs = [shifted_lp_poly(S, p, [0,0])]

#     prec = 400

#     vol1 = volume1(fs, 0, {}, prec)
#     vol2 = volume2(fs, prec)

#     target_accuracy = 10^(-100)

#     CBF = ComplexBallField(prec - 2)
#     assert(abs(CBF(vol1) - CBF(vol2)) < target_accuracy)



# def test25():
#     """ Check that error is raised when the initial points are not good, i.e. the system is not good.
#     """
#     # TODO
#     print("[TEST25] Not yet implemented!")
#     pass

# def test26():
#     """ Make sure that the output of volume 1 is actually an element of a complex ball field of expected precision.
#     """


#     prec_bits = 200

#     n = 3
#     p = 4
#     r = 1

#     S = PolynomialRing(QQ, "x", n)
#     x = S.gens()


#     fs = [shifted_lp_poly(S, p, np.zeros((n,), int))]
#     vol_p_n = volume1(fs, 0, {}, prec_bits)

#     # print("Parent of the vol1 output expression: {}".format(vol_p_n.parent()))
#     assert((vol_p_n).parent() == ComplexBallField(prec_bits))
    




if __name__ == "__main__":
    global ONLY_FAST_TESTS, MUST_DO
    ONLY_FAST_TESTS = True
    MUST_DO = [] # Enter the number of the test here.


    # test0()
    # test1()
    # test2()
    # test3()
    # test4()
    # test5()
    # test6()
    # test7()
    # test8()
    # test9()
    # test10()
    # test11()
    # test12()
    # test13()
    # test14()
    # test15()
    # test16()
    # test17()
    # test18()
    # test19()
    # test20()
    # test21()
    # test22()
    # test23()
    # test24()
    # test25()
    # test26()
    # test27()

    print("Testing the msolve interface.")

    failures = []
    for test in msolve_interface_tests:
        try:
            test()
        except AssertionError:
            failures.append((str(test), "Failed with assertion error."))

    print("Testing the m2 interface.")
    for test in m2_interface_tests:
        try:
            test()
        except AssertionError:
            failures.append((str(test), "Failed with assertion error."))

    print("Testing the tools.")
    for test in tools_tests:
        try:
            test()
        except AssertionError:
            failures.append((str(test), "Failed with assertion error."))

    
    print("Testing the volume methods.")
    for test in volume_tests:
        try:
            test()
        except AssertionError:
            failures.append((str(test), "Failed with assertion error."))


    do_long_tests = input("Do long tests as well? y / [n]")
    if do_long_tests == "y":
        print("Running computational tests.")
        for test in long_tests:
            try:
                test()
            except AssertionError:
                failures.append((str(test), "Failed with assertion error."))

    if len(failures) == 0:
        print("All good.")
    else:
        print("Encountered failures: ")
        for failure in failures:
            print(f"{failure[0]}: {failure[1]}")