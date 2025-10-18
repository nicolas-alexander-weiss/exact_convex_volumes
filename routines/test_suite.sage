# To run this test suite, make sure to have sage-preparsed routine.sage first,
# so that we can import the defined methods.

load("routine.sage")
load("sageM2.sage")

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

if __name__ == "__main__":
    test0()
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()