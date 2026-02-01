# sys.path.insert(0, '/Users/lakshmiramesh/Desktop/ore_algebra/src')

# Helper Modules

from . import tools
from . import msolve_interface as msolve
from . import m2_interface as m2

# Ore Algebra Package

from ore_algebra import OreAlgebra

# Python Imports

import numpy as np

# Sage Imports

from sage.all import (
    ZZ, QQ, AA,
    PolynomialRing,
    RealBallField, ComplexBallField,
    diff, log, power,
    vector, matrix
)

# Custom exceptions:

class BadPointsError(Exception):
    pass

class PosDimCritLocusError(Exception):
    pass

class CertificateError(Exception):
    pass

# Global parameters

NUM_BITS_PRECISION = 200 # i.e. Precision of 1 / 2^NUM_BITS_PRECISION.

# Picard-Fuchs-Operator-Caching
class PicardFuchsCache:
    """Cache object to store computed Picard-Fuchs-Operators.

    Note that this only builds a string representation, so this is sensitive to any kind of reordering or renaming of inputs.
    
    This way, computations can be re-run to increase precision.
    """

    def __init__(self):
        self.cache = {}     # for the picard-fuchs operators of the deformed slices
        self.cache_t = {}   # for the picard-fuchs operator of the deformation family
    
    ### PF versions (after deformation)

    def PF_input_repr(fs, def_value, var_value_pairs, proj_var):
        """ Flattens the inputs together as a tuple pairs denoting the input. 
        Dicts will be flattened to a tuple of key value pairs.      
        """
        input_params = (("fs",tuple(fs)), ("def_value", def_value), ("var_value_pairs", tuple(var_value_pairs.items())), ("proj_var",proj_var))

        return input_params

    def contains_PF(self, fs, def_value, var_value_pairs, proj_var):
        """ A computed Picard-Fuchs operator should be characterized by the cohomology class that it integrates.
        In our case, since the Picard-Fuchs operator will be returned from the from the "get_picard_fuchs()" and
        "get_picard_fuchs_t()" function, the input will be either the fs, or additionally the deformation and slice
        slice parameters (var_value_pairs).

        [TODO] Extend this doc string.

        Input
        ------
            fs : Polynomials over QQ of the same polynomial ring, common positivity locus defining the region of interest.
            proj_var : The axis onto which the region will be projected. The variable in which the PF operator is defined.
            def_value : prod(fs) - def_value will define the deformed slices.
            var_value_pairs : At which the above expression will be evaluated.
        """
        return PicardFuchsCache.PF_input_repr(fs, def_value, var_value_pairs, proj_var) in self.cache.keys()

    def add_PF(self, P, fs, def_value, var_value_pairs, proj_var):
        """ Builds a representation and then stores the operator P in the dictionary.
        """
        representation = PicardFuchsCache.PF_input_repr(fs, def_value, var_value_pairs, proj_var)
        self.cache[representation] = P


    def retrieve_PF(self, fs, def_value, var_value_pairs, proj_var):
        """ Builds a representation and retrieves the PF operator from the cache (dictionary).
        """
        representation = PicardFuchsCache.PF_input_repr(fs, def_value, var_value_pairs, proj_var)
        return self.cache[representation]

    ### PF_t versions:

    def PF_t_input_repr(fs, def_var_name):
        """
        [TODO] In principal should also have a way to go backwards.
        """
        input_params = (("fs", tuple(fs)), ("def_var_name",def_var_name))
        return input_params

    def contains_PF_t(self, fs, def_var_name):
        """
        See doc string of contains_PF.
        """
        return PicardFuchsCache.PF_t_input_repr(fs, def_var_name) in self.cache_t.keys()

    def add_PF_t(self, P, fs, def_var_name):
        """ Builds a representation and then stores the operator P in the dictionary.
        """
        representation = PicardFuchsCache.PF_t_input_repr(fs, def_var_name)
        self.cache_t[representation] = P

    def retrieve_PF_t(self, fs, def_var_name):
        """ Builds a representation and retrieves the PF operator from the cache (dictionary).
        """
        representation = PicardFuchsCache.PF_t_input_repr(fs, def_var_name)
        return self.cache_t[representation]

def get_inside_points(fs, var_value_pairs, points):
    """
    Returns those pts in points that lies in the restricted intersection of the fs,
    i.e. they satisfy 
        f(var_value_pairs | point) > 0 for all f in fs.
    
    Input
    --------
    fs : Polynomials over QQ of the same polynomial ring, common positivity locus defining the region of interest.
    var_value_pairs : At which the above expression will be evaluated.
    points : Points in the remaining variables that shall be checked for containment.
    """
    fs_restricted = [tools.partial_eval_poly(f, var_value_pairs) for f in fs]
    # for pt in points:
    #     print([tools.partial_eval_poly(f, pt, infer_target_base_ring = True)  for f in fs_restricted ])
    return [pt for pt in points if all([tools.partial_eval_poly(f, pt, infer_target_base_ring = True) > 0 for f in fs_restricted ])]
    

def get_1_dim_volume(fs, var_value_pairs, def_value, prec=NUM_BITS_PRECISION):
    """ In the setting of shifted lp-balls for p even, this return volume of the 
    deformed intersection of lp-balls ((Prod fs) - t) after intersection with the line defined by the 
    coordinate values for all but one coordinate.

    Input
    ------
        fs : List of polynomials in the same polynomial ring over QQ. 
        var_value_pairs : List of pairs (variable, value)
        def_value : Value, by which the product of the fs will be deformed.
        prec : Natural number, defining precision up to 2^(-prec).

    Output
    ------
        Length of the line defined by var_value_pairs intersecting with the deformed intersection of the {f >= 0}.

    Caveat:
    ------ 
        Assumes that the line actually intersects the deformed object and 
        that $t$ lies between 0 and the first critival value of the projection 
        from {(prod(fs) - t) = 0} to the t axis.

    TODO
    ------
    # TODO: Raise error when variable appears twice. 
    # TODO: Also need to make the caveat more precise.

    [TODO] Rephrase both in terms of convexity, vs. just the two points within all the other {f > 0} for f in fs.
            and then check that there are only two in the intersection.
    """
    # Deform the product to prod(fs) - t by the specified value for t.
    def_poly = tools.eval_poly(tools.deformed_product(fs), def_value)

    # Partially evaluate def_poly to land in a univariate polynomial ring.
    univariate_poly = tools.partial_eval_poly(def_poly, var_value_pairs) 
    var_name = univariate_poly.parent().gens()[0]

    # Solutions
    real_variety = msolve.variety_msolve([univariate_poly], prec) 
    num_roots = len(real_variety)
    
    if not (num_roots >= 2):
        raise ValueError(f"[1DimVol] Expected slice at to intersect in at least two points, but got {num_roots} for def_value ={def_value} and var_value_pairs = {var_value_pairs}.")

    real_values = []
    if len(fs) == 1 and def_value == 0:
        # There is a single undeformed concave polynomial. In that case, there exist
        # just two intersection points:
        if not (num_roots == 2):
            raise ValueError(f"[1DimVol] For a non-deformed single polynomial, expected exactly two points in the intersection, but got {num_roots} for def_value ={def_value} and var_value_pairs = {var_value_pairs}.")
        real_values = [pt[var_name] for pt in real_variety]
    else:
        # Extract the pts that lie inside f > 0 for all f in fs:
        pts_inside_intersection = get_inside_points(fs, var_value_pairs, real_variety)
        if not (len(pts_inside_intersection) == 2):
            raise ValueError(f"[1DimVol] In the concave case, lines can intersect the deformed intersection only in 2 pts, but got {pts_inside_intersection} for def_value ={def_value} and var_value_pairs = {var_value_pairs}.")
        real_values = [pt[var_name] for pt in pts_inside_intersection]

    real_values.sort()

    # It is ensured by the above that there are exactly 2 points:
    return real_values[1] - real_values[0]


    # OLD HEURISTIC: The relevant line segment is bounded by the middle two real roots of univariate_poly,
    # see also the respective proposition in our paper. Will have evenly many roots in this setting.

    # Sort the roots in increasing order:
    real_values = [pt[var_name] for pt in real_variety]
    real_values.sort()
    # print("Real values: {}".format(real_values))

    # The following difference will necesarily be positive.
    # [TODO] This uses the fact that in the lp ball setting,
    #           the relevant intersection points will be the middle two.
    #           Should change this to the msolve base approach. 
    #           (Since we just need to check what leads to positive values.)
    return real_values[num_roots // 2] - real_values[num_roots // 2 - 1]


def get_inside_critical_points(fs, def_value, proj_var, var_value_pairs, prec=NUM_BITS_PRECISION, debug_level=0):
    """ Return the crit points of prod(fs)-t, evaluated at def_value and var_value_pairs, that lie
    inside {f > 0} for all f in fs.

    Input
    ------
        fs              : A list of multivariate polynomials
        def_value       : In QQ, value by which the product is deformed.
        proj_var        : The variable, such that we project onto the proj_var-axis.
        var_value_pairs : The values already restricted to.
    
    Output
    ------
        inside_crit_points : List of points of the form {xi:valxi ... for all i}, where i 
                                goes over all variables that remain after the value substitutions.

    Remark
    ------
        The resulting computed critical points are those that lie on boundary of the smooth deformation of the intersection.
    """
    assert(def_value > 0) # The branch points would satisfy f==0 for some of the f in fs, and hence this case is not allowed.

    # Deform and compute branch points of the projection.
    fdef = tools.partial_eval_poly(tools.eval_poly(tools.deformed_product(fs), [def_value]), var_value_pairs) 

    #fs_restricted = [tools.partial_eval_poly(f, var_value_pairs) for f in fs]

    try:
        crit_points = get_critical_points(fdef, proj_var, prec) 

        if debug_level > 0:
            print("Computed critical points for the proj_var = {} are: {}".format(proj_var, crit_points))

        # Now identify the critical points satisfying f > 0 for all f in fs:
        inside_crit_points = get_inside_points(fs, var_value_pairs, crit_points)
        
        #[point for point in crit_points if all([tools.partial_eval_poly(f, point, infer_target_base_ring = True) > 0 for f in fs_restricted ])]
    
    except PosDimCritLocusError as error:
        print("Critical locus of projection is positive dimensional: {}".format(error))
        raise error
    
    return inside_crit_points

def get_critical_points(f, proj_var, prec=NUM_BITS_PRECISION):
    """ Computes all the real critical points, i.e. points on V(f) relative
        to the projection onto the proj_var-axis.

        The critical points of the projection are the solutions to the ideal 
            I = (f) + ( diff(f, x_i) | x_i != proj_var).
        
        The solutions are computed using groebner bases and real root isolation in msolve.
    
    Input
    ------
        f : Multivariate polynomial over QQ
        proj_var : A variable out of f.parent().gens()

    Output
    ------
        A representation of the critical points.

    Caveat
    ------
        Raises a PosDimCritLocusError if the critical locus is not 0 dimensional.
    """
    R = f.parent()

    # critical locus system.
    system = [f] + [diff(f, x) for x in R.gens() if not (x == proj_var)]

    # Now solve the ideal, assuming it is just a collection of points.
    if R.ideal(system).dimension() != 0:
        raise PosDimCritLocusError("The system {} is not 0 dimensional! It has dim = {}".format(system, R.ideal(system).dimension()))

    variety = msolve.variety_msolve(system, prec)

    return variety

def project_deformed_intersection(fs, def_value, proj_var, var_value_pairs, prec=NUM_BITS_PRECISION, debug_level=0):
    """ Computes the minimum and maximum value that proj_var takes on the deformed volume_intersection,
    after restricting to the chosen var_value_pairs.

    Input
    -------
        fs              : A list of multivariate polynomials
        def_value       : In QQ, value by which the product is deformed.
        proj_var        : The variable, such that we project onto the proj_var-axis.
        var_value_pairs : The values already restricted to.

    Output
    -------
        An interval representing the [min, max] value of proj_var on the deformed intersection.
    
        This can be done in terms of computing real branch_points lying on the boundary 
        of the deformed intersection and then projecting onto proj_var.

    Caveat
    -------
        The semi-algebraic set "intersection" refers to the intersection of {f_i >= 0}.
        The "deformed semi-algebraic set" refers to the connected component of {prod(f) - t >= 0}
        that lies within the "intersection".

    """
    # If only 1 lp ball, then no need to deform!
    if def_value == 0:
        # If not deformed, only 1 poly (smooth bdry!) supported
        assert(len(fs) == 1)
        return [{proj_var:point[proj_var]} for point in get_critical_points(fs[0], proj_var, prec)]    

    # Now identify the branchpoints satisfying f >= 0 for all f in fs:
    try:

        # TODO: Rename this to inside critical points!
        inside_branch_points = get_inside_critical_points(fs, def_value, proj_var, var_value_pairs, prec, debug_level=debug_level)
        # Now return only the proj_var values 
        return [{proj_var:point[proj_var]} for point in inside_branch_points]

    except PosDimCritLocusError:
        # In this case we are computing points in the complement instead:
        print("Reverting to computing the relevant critical values using HypersurfaceRegions.jl")

        fdef = tools.partial_eval_poly(tools.eval_poly(tools.deformed_product(fs), [def_value]), var_value_pairs)
        
        # Compute the polynomial defining the critical locus:
        crit_val_ideal = m2.branchIdeal(fdef, proj_var)
        if len(crit_val_ideal.gens()) != 1:
            raise ValueError("Expected ideal of critical values to have only one generator. Instead got: {}".format(crit_val_ideal))
        
        # Extract the polynomial. Ideally, we should map it to the same ring as fdef, but this is not necessary here.
        crit_val_poly = crit_val_ideal.gens()[0]
        pot_crit_vals = crit_val_poly.roots(AA, multiplicities=False) 
        
        if debug_level > 0:
            print("pot_crit_vals: ", pot_crit_vals)

        # Setting up the input to julia:
        print("Please execute the below in Julia: ")
        print("using HypersurfaceRegions")
        print("@var " + " ".join([str(xi) for xi in fdef.parent().gens()]) + ";")
        print("fdef = {};".format(str(fdef)))
        print("critValPoly = {};".format(str(crit_val_poly)))
        print("system = [fdef;critValPoly];")
        print("regs = regions(system);")
        #print('println("[")')
        print('println("["); for region in regs.region_list\n  println(string(region.critical_points[1]) * ",")\nend; println("]")')
        #print('println("]")')

        # TODO: Should actually make sure the precision is high enough!
        sampled_points = sage_eval(input("Please input the sampled points from Julia as one line:"))
        
        # Todo do with a join
        sampled_points = [dict( zip(fdef.parent().gens(), pnt) )  for pnt in sampled_points]
        
        print("\nReceived the sampled points: {}".format(sampled_points))

        # Identify the relevant critical values (in the roots of crit_val_poly)

        # Check which of the sample points lie in our region: fdef > 0 and fi > 0 for all fi in fs:
        inside_points = []
        for point in sampled_points:
            # TODO: Can also merge the pnt with the var_value_pairs instead of restricting the polynomials.
            # TODO: var_value pairs have QQ values and the point has real values.
            if all([tools.partial_eval_poly(fun, point | var_value_pairs ) > 0 for fun in [fdef] + fs]):
                inside_points.append(point) 

        if debug_level > 0:
            print("The polys: ",[fdef] + fs)
            print("Print inside_points: ", inside_points)
        
        minval_inside_points = min([pt[proj_var] for pt in inside_points])
        maxval_inside_points = max([pt[proj_var] for pt in inside_points])

        if debug_level > 0:
            print("Min and max inside value:", minval_inside_points, maxval_inside_points)

        relevant_crit_values_numerical = [ max([pt for pt in pot_crit_vals if pt < minval_inside_points]),
                                 min([pt for pt in pot_crit_vals if pt > maxval_inside_points])
        ]
        
        return [{proj_var:val} for val in relevant_crit_values_numerical]

    


def get_picard_fuchs(fs, def_value, var_value_pairs, proj_var, strategy=None, debug_level=0):
    """ Computes the Picard Fuchs operator for 
    Vol(proj_var) = Vol( p^{-1}(proj_var) \cap {fs \geq 0 forall s} \cap slice(var_value_pairs)).

    Input
    ------
    fs              : Polynomials defining semi-alg set (by fs >= 0)
    def_value    : Value in QQ, by which to smooth prod(fs).
    var_value_pairs : defining slice to restrict to.

    Output
    ------
    P : Picard Fuchs operator, in WeylAlgebra D_{proj_var}

    TODO
    -----
    Make the procedure by which to compute the intersection ideal more informed.
    """

    # Construct the integrand
    A = construct_integrand(fs, def_value, var_value_pairs, proj_var, debug_level=debug_level)
    
    W = rational_weyl_algebra(A.parent().ring())

    annA = W.ideal([A*D-D(A) for D in W.gens()]) # construct the annihilating ideal

    # To be precise, below we simply construct some subset of the integration ideal, 
    # but it suffices to be non-empty.
    allowed_pole = A.denominator().change_ring(QQ)
    intIdeal = creative_telescoping(annA, proj_var, allowed_pole=allowed_pole, debug_level=debug_level)

    return intIdeal.gens()[0]

def annihilator_deformed_intersection(fs, def_value, var_value_pairs, proj_var, debug_level=0):
    """ Returns an annihilating d-finite ideal in the rational Weyl algebra 
        for the rational integrand 
        
        d_i(F) * x_i / F 
        where x_i not = proj_var
        and 
        F(x,t) = prod(fs) - t

        after evaluating at t=def_value and var_value_pairs.
    """

    A = construct_integrand(fs, def_value, var_value_pairs, proj_var, debug_level=debug_level)
    
    W = rational_weyl_algebra(A.parent().ring())

    annA = W.ideal([A*D-D(A) for D in W.gens()]) # construct the annihilating ideal
    
    return annA

def construct_integrand(fs, def_value, var_value_pairs, proj_var, prim_var_name=None, debug_level=0):
    """ By standard considerations, see for example our paper, the function 
    vol(proj_var) can be expressed as period of a rational function A. 
    We construct this function rational function here, considering  closely the order of indeterminate variables
    and the proj var.

    Input
    ------
    See get_picard_fuchs()

    Output
    -----
    A : Element of a FractionField. (TODO: Specify more clearly the variables.)
    """

    fdef = tools.partial_eval_poly(tools.eval_poly(tools.deformed_product(fs), [def_value]), var_value_pairs)
    
    if debug_level > 0:
        print("[PF] fdef = {}".format(fdef))

    # Multiply out the denominators and change to integer ring:
    fdef_ZZ = fdef.numerator().change_ring(ZZ)
    proj_var_ZZ = fdef_ZZ.parent()(proj_var)

    prim_var = [var for var in fdef_ZZ.parent().gens() if not var == proj_var_ZZ][0] # TODO: Here taking the first, should make this more strategic here.
    
    # TODO: Test this before committing.
    if prim_var_name != None:
        prim_var = fdef_ZZ.parent()(prim_var_name)

    sgn = (-1)**0 # TODO: Should depend on the choice of prim_var, to account for ordering.

    A = sgn * diff(fdef_ZZ, prim_var) * prim_var / fdef_ZZ # Automatically constructs the fraction field.

    return A


def get_picard_fuchs_t(fs, strategy=None, debug_level=0):
    """ Computes the Picard Fuchs operator for 
    Vol(t) = Vol( prod(fs) - t >= 0) \cap {fs \geq 0 forall s} 

    Input
    ------
    fs              : Polynomials defining semi-alg set (by fs >= 0)

    Output
    ------
    P : Picard Fuchs operator, in WeylAlgebra D_{t}

    TODO
    -----
    Make the procedure by which to compute the intersection ideal more informed.
    """

    def_var_name = "t"

    # Construct the integrand, note that here this is over the ring QQ[x..., t], rather than QQ[x..][t]
    At = construct_integrand_t(fs, def_var_name, strategy)

    Wt = rational_weyl_algebra(At.parent().ring())

    annAt = Wt.ideal([At*D-D(At) for D in Wt.gens()]) # construct the annihilating ideal

    # To be precise, below we simply construct some subset of the integration ideal, 
    # but it suffices to be non-empty.
    allowed_pole = At.denominator().change_ring(QQ)
    intIdeal_t = creative_telescoping(annAt, At.parent().ring()(def_var_name), allowed_pole=allowed_pole, debug_level=debug_level)

    return intIdeal_t.gens()[0]


def construct_integrand_t(fs, def_var_name, strategy=None, prim_var_name=None):
    """ See construct_integrand(). Constructs the picard fuchs operator for the deformed slice
    vol(t) = vol({prod(fs)-t >= 0} \cap {fs >= 0})

    Input
    -----
        See get_picard_fuchs_t()
    """

    ft = tools.deformed_product(fs, def_var_name)

    ft_flattened = ft.parent().flattening_morphism()(ft) # flatten the ring.

    # Multiply out the denominators and change to integer ring:
    ft_flattened_ZZ = ft_flattened.numerator().change_ring(ZZ)


    proj_var_ZZ = ft_flattened_ZZ.parent()(def_var_name)

    prim_var = [var for var in ft_flattened_ZZ.parent().gens() if not var == proj_var_ZZ][0] # TODO: Here taking the first, should make this more strategic here.
    
    # TODO: Test this before committing.
    if prim_var_name != None:
        prim_var = ft_flattened_ZZ.parent()(prim_var_name)



    sgn = (-1)**0 # TODO: Should depend on the choice of prim_var, to account for ordering. (Or maybe it doesn't matter? Either does annihilate it.)

    At = sgn * diff(ft_flattened_ZZ, prim_var) * prim_var / ft_flattened_ZZ # Automatically constructs the fraction field.

    return At

def rational_weyl_algebra(polyRing):
    """ Constructs the rational Weyl algebra for a specified polynomial ring.
    """
    fracField = polyRing.fraction_field()

    return OreAlgebra(fracField, *[("D" + str(var), {}, {var : polyRing(1)}) for var in polyRing.gens()])

def creative_telescoping(I, proj_var, allowed_pole=None, strategy=None, debug_level=0):
    """ For an ideal in the rational Weyl algebra W, it carries out creative telescoping
    sequentially to eliminate all but proj_var.

    The result is an ideal contained in:     (I + dx_1 * D + ...+ dx_n *D) \cap D_x0
    (In this example we assume x0 to be the proj_var)

    Output
    ------
    ct_ : Returns a list of of univariate differential operators of length 1 (!).

    TODO
    ------
    - Should extend this, so that at least temporarily we store / output for debuginfo the certificates, too.
    - Should include some basic strategy.
    """

    if debug_level > 0:
        print("[CT] Entering creative telescoping:")

    W = I.ring()
    ct_system = list(I.gens())

    if debug_level > 2:
        print("\nThe ideal is generated by: {}".format(list(I.gens())))

    for Dvar in [Dvar for Dvar in W.gens() if not Dvar == W("D"+str(proj_var))]:
        ct_ideal = ct_system[0].parent().ideal(ct_system)
        
        if debug_level > 0:
            print("\n[CT] Integrating out the variable {}".format(str(Dvar)))

        # ct_system, certificates = ct_ideal.ct(Dvar, certificates=True)
        ct_system = ct_ideal.ct(Dvar, certificates=False)

        # Verify that the certificates only contain allowed poles. 
        # Factor first and then check for the factors.
        # TODO: This can likely be done better. (i.e. check)

        # if allowed_pole!= None:
        #     for cert in certificates:
        #         cert_denom = (W.base_ring()(cert)).denominator()
        #         cert_denom = cert_denom.change_ring(QQ)

        #         if cert_denom in QQ: # Skip if its a constant
        #             continue

        #         # This has to be looked into again.
        #         if not (cert_denom.radical() % allowed_pole.radical() == 0 and allowed_pole.radical() % cert_denom.radical() == 0):
        #             raise CertificateError("[CT] The certificate contains a pole that is not appearing in fdef:\ncert_denom.radical().factor() = {}\nallowed_poly.radical().factor() = {}".format(list(cert_denom.radical().factor()), list(allowed_pole.radical().factor())))
                        
        # if debug_level > 2:
        #     print("\n[CT] The resulting system: {}".format(ct_system))
        #     print("\n[CT] The corresponding certificates: {}".format(certificates))


    return ct_system[0].parent().ideal(ct_system)

def solve_diff_op(P, initial_conditions, evaluation_condition, prec, apparent_singular_points=[], debug_level=0):
    """ Given a linear differential operator in 1 variable, i.e. an ODE, solve it given the provided initial conditions
    and output the value of the solution at the requested point.

    The initial conditions specify the coefficients to the local series solutions, given the exponents of the starting monomials.
    In our case, since our solutions are bounded in the specified range, there will be no log part in the starting monomials.

    The initial condition solver is used in two scenarios:
    - Vol1:
        - Here the operator to be considered is P*dx, i.e. we solve for the value of the integral \int_a^x vol(x0)dx0.
        - At a, the corresponding starting monomial is 1 with value 0. (RMK. "a" will be a branchvalue, i.e. a singular point of the operator!)
        - At the other points the initial data is the volume of the slice at a given value of x. Thus, the corresponding starting monomial is (x-x_0)^1. 
        - The evaluation point will be the next branch_value, i.e. again a singular point of the picard fuchs operator.
    - Vol2: 
        - Directly solve the picard fuchs operator in t, i.e. for the slices of the deformed intersection and then analytically continue to t=0.
        - The initial conditions all are taken at smooth points, smaller than the first branch-value of the projection onto t. 
        - The corresponding starting monomials will always be (t-t_0)^0 = 1.


    Input
    ------
    P : Element of Weyl algebra over polynomial or rational function ring in 1 variable.
    initial_conditions : { x_val_1:{"exponent":expon, "coef":coef }, x_val_2:...}
    evaluation_condition : Dictionary containing evaluation point and the exponent 
                            of the starting monomial, whose coef shall be read of, e.g. {"pt":0, "exponent":0}

    Example
    -------
    initial_condition = {1/10:{"exponent":1, "coef":3.12223532...}, }

    Caveat
    ------
    - [TODO] Since we really depend on the order of the initial conditions to be always the same, we should store them as a list.
    - [TODO] Consider the singular points of the differential operator.
    - This is still buggy in general.
    - [TODO] Need proof that this will work in our considered scenarios.
    - [TODO] Need to check that the initial conditions are good (and can do this before evaluating I think! So we don't waste the computations of the slice volumes.)


    - [TODO] Assure that the accuracy is in the correct ball precision. (Rather than providing a float as accuracy).
            --> Or all together change the accuracy! (e.g. to be input rather as decimal digits. e.g. as 1e-100)

    """

    if debug_level > 0:
        print("\n[ICS] Entering the initial condition solver. ")

    # Assert that univariate differential operator
    assert(len(P.parent().gens()) == 1 and len(P.parent().base_ring().gens()) == 1)

    # Extract the variable:
    t = P.parent().base_ring().gens()[0]

    # Assert that the starting monomials corresponding to the exponents do exist in the standard monomials
    assert(all([(t-pt)**condition["exponent"] in P.local_basis_monomials(pt) for pt,condition in initial_conditions.items()]))

    # Assert that we provide as many initial conditions as the order of the operator
    assert(P.order() == len(initial_conditions.keys()))

    # Set up the linear system of rk = P.order() 

    # Choose a base point for the linear system, e.g. the evaluation point
    eval_point = evaluation_condition["pt"]

    # Set up our variables as elements in polynomial ring over the complex ball field:
    CBF = ComplexBallField(prec)
    R = PolynomialRing(CBF, "a", P.order())
    a = R.gens()

    # coef_vector:
    ini_eval_point = list(a) # [a0, a1,..., a{rk-1}]

    # Transition to the other points and yields linear equations in the a's.
    # More precisely: Given the unkown coefs "a" at the base point, yields an affine linear expression in "a"
    #   for the coef of the std_monomial provided in the initial data.
    lin_eqs = [(P.numerical_transition_matrix([eval_point] + [xi for xi in reversed(sorted(apparent_singular_points))] + [pt], eps=2**(-prec)) * vector(ini_eval_point))[P.local_basis_monomials(pt).index((t-pt)**condition["exponent"])] for pt,condition in initial_conditions.items()]

    # Extract the constant terms in the affine linear equations.
    constant_terms = vector([tools.eval_poly(lin_eq, [CBF.zero() for ai in a]) for lin_eq in lin_eqs])

    # Express in terms of a matrix (over the complex ball field defined above):
    M = matrix([[CBF(lin_eq.coefficient(ai)) for ai in a ] for lin_eq in lin_eqs])

    # Check that determinant is non-zero (i.e. the complex ball of the determinant doesn't contain 0.)
    if not (M.determinant() != 0):
        raise BadPointsError("Determinant of the system cannot be distinguished from zero.")

    # Set up the initial data as a vector
    initial_data_vector = vector([condition["coef"] for pt,condition in initial_conditions.items()]).change_ring(CBF)

    # Solve
    eval_point_coefs = M.change_ring(CBF).inverse() * (initial_data_vector - constant_terms)

    # Read of the coefficient that we wanted to evaluate at:
    return eval_point_coefs[P.local_basis_monomials(eval_point).index((t-eval_point)**evaluation_condition["exponent"])]

def volume1(fs, def_value, var_value_pairs, prec=NUM_BITS_PRECISION, strategy=None, debug_level=0):
    """ Computes the smooth volume of the deformed intersection of lp balls. 

    Input
    -----
    fs : list of multivariate polynomials over QQ
    def_value : non-negative element in QQ
    var_value_pairs : Dictionary with  variable:value   key-value-pairs, where value is in QQ.

    Output
    ------
    The volume of the deformation of the intersection of all {f>=0 | f in fs}
    to prod(fs)-def_value. In both case only the points in the slice specified by var_value_pairs are considered.

    The output is an element of a complex ball field!

    Caveat
    ------
    The output will be a number in the complex ball field with prec many bits of precision.



    [TODO] Add tests to test-suite.
    [TODO] Add assertions.
    [TODO] Implement good initial conditions.
    [TODO] Check if branch_points computation also works to t=0 case (so that don't need separate case.)

    [TODO] Make the choice of a projection variable strategic.

    [TODO] Assert that there is actually only a unique root in the singular locus corresponding to either numerical branch_point.

    """
    # The evaluated polynomial (for reference)
    evaluated_poly = tools.partial_eval_poly(tools.eval_poly(tools.deformed_product(fs), [def_value]), var_value_pairs)
    if debug_level > 0:
        print("Deformation value t: {}".format(def_value))
        print("Slice taken at: {}".format(var_value_pairs))
        print("Restricted deformed product: {}".format(evaluated_poly))

    # Early exit if already univariate, then return 1-dim volume.
    if len(evaluated_poly.parent().gens()) == 1:
        if debug_level > 0:
            print("The restricted polynomial is univariate, hence going into the 1-dim-volume routine.")
        return get_1_dim_volume(fs, var_value_pairs, def_value, prec)

    # Fix a projection variable (here just the first undetermined variable): 

    if strategy == None:
        proj_var = evaluated_poly.parent().gens()[0]
    else:
        proj_var = strategy["proj_var"](fs, def_value, var_value_pairs)

    if debug_level > 0:
        print("proj_var: {}".format(proj_var))

    # Get the Picard-Fuchs operator and define the operator to be solved.
    P = get_picard_fuchs(fs, def_value, var_value_pairs, proj_var, strategy, debug_level=debug_level)
    Pdx = P * P.parent().gens()[0]

    lead_coef = P.leading_coefficient().numerator()
    singular_locus = lead_coef.roots(AA, multiplicities=False) 

    if debug_level > 0:
        print("\nOrder Pdx:", Pdx.order())
        print("Leading coef:", Pdx.leading_coefficient())
        print("Singular locus:", singular_locus)
        if debug_level > 2:
            print("Pdx = ", Pdx, "\n")

    # Determine the branch points and corresponding critical values bounding the deformed set.
    
    # Numerical critical values:
    relevant_crit_val = [pt[proj_var] for pt in project_deformed_intersection(fs, def_value, proj_var, var_value_pairs, prec)]
    if debug_level > 0:
        print("Relevant crit vals:", relevant_crit_val)

    if len(relevant_crit_val) != 2:
        raise ValueError("Computation returned too many (#={}) critical values: {}".format(len(relevant_crit_val), relevant_crit_val))
    
    min_crit_val = min(relevant_crit_val)
    max_crit_val = max(relevant_crit_val)

    if debug_level > 0:
        print("Relevant critical val: ".format(relevant_crit_val))

    # Set the interval in which we want to sample points, identify our branch points among
    # TODO: Make sure that this uniquely identifies points of the singular locus.
    xmin = singular_locus[np.argmin([np.abs(root - min_crit_val) for root in singular_locus])]
    xmax = singular_locus[np.argmin([np.abs(root - max_crit_val) for root in singular_locus])]
    
    if debug_level > 0:
        print("Identified the relevant critical values in the singular locus as: \nxmin = {}\nxmax = {}".format(xmin, xmax))

    # apparent singularities (that are not singularities of our solution:)
    xapparent = sorted([xi for xi in singular_locus if not xi in [xmin, xmax] and (xmin < xi) and (xi < xmax)])

    if debug_level > 0:
        print("Apparent singular points in the interval: ", xapparent)

    # Determine points for initial data, such that transition matrix becomes invertible.
    CBF = ComplexBallField(prec)
    lower_bound = CBF(xmin).real()
    upper_bound = CBF(min(xapparent)).real() if len(xapparent) > 0 else CBF(xmax).real()
    initial_points = get_good_initial_points(P, lower_bound, upper_bound, debug_level=debug_level)
    
    if debug_level > 0:
        print("sampled initial points:", initial_points)

    # Determine initial conditions (TODO: in parallel)
    #{"pt":0, "exponent":0} { x_val_1:{"exponent":mon, "coef":coef }, x_val_2:...}
    initial_conditions = {xmin:{"exponent":0, "coef":0}} | {proj_var_val:{"exponent":1, "coef":volume1(fs, def_value, var_value_pairs=var_value_pairs | {proj_var:proj_var_val}, prec=prec, strategy=strategy, debug_level=debug_level)} for proj_var_val in initial_points}
    
    if debug_level > 0:
        print("Initial conditions", initial_conditions)

    evaluation_condition = {"pt":xmax, "exponent":0}

    if debug_level > 0:
        print("Eval condition: ", evaluation_condition)

    # Solve the initial value problem
    volume = solve_diff_op(Pdx, initial_conditions, evaluation_condition, prec, xapparent, debug_level=debug_level)

    # Return the volume
    return volume

def get_good_initial_points(P, x0, x1, debug_level=0):
    """ Provides n rational points strictly between x0 and x1, 
    where n is the order of the ordinary differential operator to be solved.

    Input
    -----
    [TODO]

    Output
    ------
    [TODO]

    Remark
    ------
    [TODO] Actually provide the invertibility check of the initial points.
    """
    n = P.order()
    q = sample_n_rational_points(x0,x1,n, debug_level=debug_level)

    # TODO: Actually check invertibility of the linear system.
    # Rmk: Note that with probability 1, the system is invertible. Ideally this check while solving the differential operator
    #   and then raising an error.

    return q


def sample_n_rational_points(x0,x1,n, base=10, debug_level=0):
    """ Linearly samples n rational points strictly between x0 and x1.
    Determines q_start and delta, accordingly, so that

    qk = qstart + k*delta (for k = 0,..., n-1)

    yields n rational points between x0 and x1.

    In case x0 and x1 are given as algebraic numbers, it first truncates the imaginary pary.

    Determines a suitable interval in base 10.

    Remark
    ------
    [TODO] - Print delta on demand. 
    [TODO] - Which log function is this? Is it accurate?
    """

    if debug_level > 0:
        print("Sampling n = {} many points between x0 = {} and x1 = {}.".format(n, x0, x1))
        print("Chosen base = {}.".format(base))

    xmax = max(x0,x1)
    xmin = min(x0,x1)

    # Determine integer N, such that 10^N < xmax-xmin 
    N = np.floor(log(xmax-xmin)/log(base))

    # Choose qstart in 2*10^(N-1) neigborhood of xmin.
    # Use QQ, to avoid floats!
    qstart = QQ((np.floor(xmin / power(QQ(10), N-1)) + 2) * power(QQ(10), N-1))

    # Generate the samples linearly
    if n == 1:
        return [qstart]

    # n != 1:
    delta = QQ(10**(N-1) / (n-1))
    q = [qstart + k*delta for k in range(0,n)]

    return q


def volume2(fs, prec, strategy=None, debug_level=0):
    """ Computes the volume of semi-algebraic convex bodies.

    Assumes that convex body is then given as
        C = {x in RR^n | f(x) > 0 for all f in fs}
    where 
        f is concave, for all f in fs.
    
    Parameters
    ----------
        fs : list of polys in QQ[x_1,..., x_n]
            Concave polynomials defining the convex body.
        prec : int
            Number of binary digits of precision.
        debug_level : int, optional
            Amount of intermediate results to be printed.
            
    Returns
    -------
    sage.rings.complex_arb.ComplexBall
        The volume of the convex body C with precision 2^{-prec}.

    Caveat
    ------
    It is assumed that the convex body C is full-dimensional.
    Concavity is NOT explicitly checked for.
    """

    if debug_level > 0:
        print("[Vol2] Start of volume2:")
        print("[Vol2] fs = {}".format(fs))

    # Early exit, if only 1 function:
    if len(fs) == 1:
        return volume1(fs,0, {}, prec, strategy)

    # TODO: Assertions
    assert(tools.all_from_same_parent_ring(fs))

    # Get Picard Fuchs
    Pt = get_picard_fuchs_t(fs, strategy, debug_level=debug_level)
    t = Pt.parent().base_ring().gens()[0]

    if debug_level > 0:
        print("[Vol2] Pt = ".format(Pt))
        print("[Vol2] order(Pt) = {}".format(Pt.order()))


    # Choose initial points between 0 and the smallest singular point larger than 0.
    lead_coef = Pt.leading_coefficient().numerator()
    sols = msolve.variety_msolve([lead_coef], prec)    # Compute real solutions with msolve
    
    # TODO: Define t here. Why does this work without?
    singular_locus = [sol[t] for sol in sols]

    if debug_level > 0:
        print("[Vol2] Lead coef = {}".format(lead_coef))
        print("[Vol2] Singular Locus = {}".format(singular_locus))

    assert(0 in singular_locus)
    index_zero = sorted(singular_locus).index(0)
    
    smallest_positive_singular_point = sorted(singular_locus)[index_zero + 1] if index_zero + 1 < len(singular_locus) else QQ(1/10)

    CBF = ComplexBallField(prec)
    initial_points = get_good_initial_points(Pt, 0, CBF(smallest_positive_singular_point).real())

    if debug_level > 0:
        print("[Vol2] Initial points for Pt: {}".format(initial_points))

    # Compute initial conditions (volumes of the t slices)
    initial_conditions = {ti:{"exponent":0, "coef":volume1(fs, ti, var_value_pairs={}, prec=prec, strategy=strategy, debug_level=debug_level)} for ti in initial_points}

    if debug_level > 0:
        print("[Vol2] Initial conditions for Pt: {}".format(initial_conditions))

    # Set evaluation condition at 0
    evaluation_condition = {"pt":0, "exponent":0}

    # Solve the initial value problem
    volume = solve_diff_op(Pt, initial_conditions, evaluation_condition, prec, [])
    
    return volume


def vol_lp_ball_closed_formula(n,p,r, prec, use_complex=True):
    """ The closed formula for the volume of an Lp ball in R^n of radius r.

    Uses the implementation of the gamma function in the complex resp real ball 
    field. The closed formula, see (en.wikipedia.org/wiki/Volume_of_an_n-ball),
    is given by:
                Vol_{n,p,r} = r^n * (2 Gamma(1/p + 1))^n / Gamma(n/p + 1) 

    Parameters
    ----------
    n : int
        Dimension of RR^n.
    p : int
        Defining lp norm.
    prec : int
        Number of binary digits of precision.
    use_complex : bool
        Whether or not to produce output as complex ball (or real).
    
    Returns
    -------
    sage.rings.real_arb.RealBall
        Volume of the Lp ball in R^n of radius r with precision 2^{-prec}.
    """
    if use_complex:
        BF = ComplexBallField(prec) 
    else:
        BF = RealBallField(prec)

    return  r**n * (2 * BF(QQ(1)/QQ(p) + 1).gamma())**n / BF(QQ(n)/QQ(p) + 1).gamma()