load("routine.sage")

def last_variable_proj_var_strategy(fs, deform_value, var_value_pairs):
    evaluated_poly = partial_eval_poly(eval_poly(deformed_product(fs), [deform_value]), var_value_pairs)

    return evaluated_poly.parent().gens()[-1]

def example_2d_triangle():
    prec = 200
    S = PolynomialRing(QQ, "x", 2)
    x = S.gens()

    f1 = -x[0] - x[1] +1/2
    f2 = x[0] - x[1] + 1/2
    f3 = x[1] 

    fs = [f1,f2,f3]

    strategy = {"proj_var":last_variable_proj_var_strategy}

    #slice_vols = dict((val, volume1(fs, deform_value, var_value_pairs={x[1]:val}, prec=prec, debug_level=3, strategy=strategy)) for val in initial_points)
    vol = volume2(fs, prec, strategy=None)
    print("The volumes of the slices: {}".format(vol))

def example_4d_square():
    prec = 200
    S = PolynomialRing(QQ, "x", 4)
    x = S.gens()

    f1 = x[0]
    f2 = x[1]
    f3 = x[2]
    f4 = x[3]
    f5 = 1 - x[0]
    f6 = 1 - x[1]
    f7 = 1 - x[2]
    f8 = 1 - x[3]

    fs = [f1,f2,f3,f4,f5,f6,f7,f8]

    strategy = {"proj_var":last_variable_proj_var_strategy}

    #slice_vols = dict((val, volume1(fs, deform_value, var_value_pairs={x[1]:val}, prec=prec, debug_level=3, strategy=strategy)) for val in initial_points)
    vol = volume2(fs, prec, strategy=None)
    print("The volumes of the slices: {}".format(vol))
