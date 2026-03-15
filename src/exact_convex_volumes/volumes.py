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
    vector, matrix,
    sage_eval
)

#
# Volumes as objects:
#

class SmoothVolume:
    """ SmoothVolume is the object which holds all intermediate results of the 
    smooth volume computation. In particular, they usually hold all the parameters
    with which the smooth volume function were computed.

    Note that the previous functions are still available, but they are wrappers to the 
    methods of the SmoothVolume class. Also, all helper functions such as for the computation
    for the PF operators remain separate.

    Features:
    - [TODO] Resume computations when they were aborted.
    - [TODO] Obtain insights into the computation afterwards
    - [TODO] Allows for recomputation with increased precision.

    Other:
    - [TODO] Have suitable getter and setter functions?

    """

    debug_level = 0 # Indicating amount of extra info printed during computation.

    def __init__(self, fs, def_value, var_value_pairs, prec, strategy=None, debug_level=0):
        """
        Docstring for __init__
        
        fs : list of multivariate polynomials over QQ
        def_value : non-negative element in QQ
        var_value_pairs : Dictionary with  variable:value   key-value-pairs, where value is in QQ.
        prec : Number of binary digits of precision, i.e. an accuracy of 2^{-prec}.

        strategy : dict (see volume computation for details)
        debug_level : int (to regulate the amount of printed details.)
        """
        self.fs = fs
        self.def_value = def_value
        self.var_value_pairs = var_value_pairs

        self.prec = prec
        self.strategy = strategy
        self.debug_level = debug_level

        # Initialize data for the computation:
        self.vol = None # in CBF(prec)
        self.PF_slice_vol = None # To hold the computed PF operator which annihilates the slice volume.
        self.proj_var = None # To hold the variable onto which we project in this step.
        self.slice_volumes = None # dict indexed by elements of QQ with values SmoothVolume
        
    def get_fdef(self):
        return tools.partial_eval_poly(tools.eval_poly(tools.deformed_product(self.fs), [self.def_value]), self.var_value_pairs) 

    def is_one_dim(self):
        """
        Returns true if only one free variable remains after 
        fixing all values in self.var_value_pairs.
        """
        return self.get_fdef().parent().ngens() == 1


    def start_computation(self):
        """ Runs the computation based on the provided data.

        [TODO] Allow for resumption of computation at later point.

        -- Refactoring into Objects:
        [Done] Go through the below and make sure that it runs based on the input to the object.
        [TODO] Adapt the debug messages to reflect the structure of the paper in fact.
        [TODO] Make 1dim volume a method.
        [TODO] Consider making the other helper functions a method, too.

        [TODO] Make initial conditions SmoothVolume objects too. 
            - First create dict of slice objects
            - Then turn into dict of initial conditions

        [TODO] Allow for parallel computation of the slice volumes.

        """
        # The evaluated polynomial
        evaluated_poly = self.get_fdef()

        if self.debug_level > 0:
            print("[Vol1] Deformation value t: {}".format(self.def_value))
            print("[Vol1] Slice taken at: {}".format(self.var_value_pairs))
            print("[Vol1] Restricted deformed product: {}".format(evaluated_poly))

        # Early exit if already univariate, then return 1-dim volume.
        if len(evaluated_poly.parent().gens()) == 1:
            if self.debug_level > 0:
                print("[Vol1] The restricted polynomial is univariate, hence going into the 1-dim-volume routine.")
            self.vol = get_1_dim_volume(self.fs, self.var_value_pairs, self.def_value, self.prec)
            return

        # Fix a projection variable (here just the first undetermined variable): 

        if self.strategy == None:
            self.proj_var = evaluated_poly.parent().gens()[0]
        else:
            self.proj_var = self.strategy["proj_var"](self.fs, self.def_value, self.var_value_pairs)

        if self.debug_level > 0:
            print("[Vol1] proj_var: {}".format(self.proj_var))

        # Get the Picard-Fuchs operator and define the operator to be solved.
        self.PF_slice_vol = get_picard_fuchs(self.fs, self.def_value, self.var_value_pairs, self.proj_var, self.strategy, debug_level=self.debug_level)
        
        Pdx = self.PF_slice_vol * self.PF_slice_vol.parent().gens()[0]

        lead_coef = self.PF_slice_vol.leading_coefficient().numerator()
        singular_locus = lead_coef.roots(AA, multiplicities=False) 

        if self.debug_level > 0:
            print("\n[Vol1] Order Pdx:", Pdx.order())
            print("[Vol1] Leading coef:", Pdx.leading_coefficient())
            print("[Vol1] Singular locus:", singular_locus)
            if self.debug_level > 2:
                print("[Vol1] Pdx = ", Pdx, "\n")

        # Determine the branch points and corresponding critical values bounding the deformed set.
        
        # Numerical critical values:
        relevant_crit_val = [pt[self.proj_var] for pt in project_deformed_intersection(self.fs, self.def_value, self.proj_var, self.var_value_pairs, self.prec)]
        if self.debug_level > 0:
            print("[Vol1] Relevant CritVals:", relevant_crit_val)

        if len(relevant_crit_val) != 2:
            raise ValueError("Computation returned too many (#={}) critical values: {}".format(len(relevant_crit_val), relevant_crit_val))
        
        min_crit_val = min(relevant_crit_val)
        max_crit_val = max(relevant_crit_val)

        # Set the interval in which we want to sample points, identify our branch points among
        # TODO: Make sure that this uniquely identifies points of the singular locus.
        xmin = singular_locus[np.argmin([np.abs(root - min_crit_val) for root in singular_locus])]
        xmax = singular_locus[np.argmin([np.abs(root - max_crit_val) for root in singular_locus])]
        
        if self.debug_level > 0:
            print("[Vol1] Identified the relevant critical values in the singular locus as: \nxmin = {}\nxmax = {}".format(xmin, xmax))

        # apparent singularities (that are not singularities of our solution:)
        xapparent = sorted([xi for xi in singular_locus if not xi in [xmin, xmax] and (xmin < xi) and (xi < xmax)])

        if self.debug_level > 0:
            print("[Vol1] Apparent singular points in the interval: ", xapparent)

        # Determine points for initial data, such that transition matrix becomes invertible.
        CBF = ComplexBallField(self.prec)
        lower_bound = CBF(xmin).real()
        upper_bound = CBF(min(xapparent)).real() if len(xapparent) > 0 else CBF(xmax).real()
        initial_points = get_good_initial_points(self.PF_slice_vol, lower_bound, upper_bound, self.prec, debug_level=self.debug_level)
        
        if self.debug_level > 0:
            print("[Vol1] Sampled initial points:", initial_points)


        # Compute slice volumes first:
        self.slice_volumes = {}
        for pt_val in initial_points:
            self.slice_volumes[pt_val] = volume1(self.fs, self.def_value, var_value_pairs=self.var_value_pairs | {self.proj_var:pt_val}, prec=self.prec, strategy=self.strategy, debug_level=self.debug_level)

        # Construct initial conditions dict:
        # {"pt":0, "exponent":0} { x_val_1:{"exponent":mon, "coef":coef }, x_val_2:...}
        self.initial_conditions = {xmin:{"exponent":0, "coef":0}} | {proj_var_val:{"exponent":1, "coef":self.slice_volumes[proj_var_val]} for proj_var_val in initial_points}
        
        if self.debug_level > 0:
            print("[Vol1] Initial conditions", self.initial_conditions)

        evaluation_condition = {"pt":xmax, "exponent":0}

        if self.debug_level > 0:
            print("[Vol1] Eval condition: ", evaluation_condition)

        # Solve the initial value problem
        self.vol = solve_diff_op(Pdx, self.initial_conditions, evaluation_condition, self.prec, xapparent, debug_level=self.debug_level)


    def get_volume(self):
        """
        Returns the volume of the smooth semi-algebraic set defined by
        {prod(fs)-def_val>0} \cap C \cap {var_value_pairs}.
        """
        if self.vol == None:
            self.start_computation()
        
        return self.vol
    
    def __repr__(self):
        """Returns a representation of the str."""
        if self.vol == None:
            return "SmoothVolume: None (use .get_volume() or .start_computation() to initiate computation)."
            
        return "SmoothVolume: " + str(self.vol)


class Volume:
    """ This class encapsulates the computation of volume2.
    In particular, it only considers a list of fs and then 
    carries out the necessary computations.

    Input
    -----
    fs : List of concave polynomials in QQ[x_1,..,x_n].
    prec : Target number of precision bits.

    Features:
    ------
    - [TODO]: Ensure that precision is really the output precision.
    
    """
    pass

# 
#
#
        
# Custom exceptions:

class BadPointsError(Exception):
    pass

class PosDimCritLocusError(Exception):
    pass

class CertificateError(Exception):
    pass

# Global parameters

NUM_BITS_PRECISION = 200 # i.e. Precision of 1 / 2^NUM_BITS_PRECISION.

# Methods

def get_inside_points(fs, var_value_pairs, points):
    """ Returns all points that satisfy f(var_value_pairs | point) > 0 for all f in fs.

    Parameters
    ----------
    fs : list of polys in QQ[x_1,..., x_n]
        Concave polynomials defining the convex body.
    var_value_pairs : dict of (x_i : val_x_i) pairs.
        The QQ values of the x_i to which the fs should be restricted.
    points : list of dict
        Of pts, given each by {var_1 : val_1,... }.

    Returns
    -------
    list of dict
        A points that satisfy f(var_value_pairs | point) > 0 for all f in fs.

    """
    fs_restricted = [tools.partial_eval_poly(f, var_value_pairs) for f in fs]
    # for pt in points:
    #     print([tools.partial_eval_poly(f, pt, infer_target_base_ring = True)  for f in fs_restricted ])
    return [pt for pt in points if all([tools.partial_eval_poly(f, pt, infer_target_base_ring = True) > 0 for f in fs_restricted ])]
    

def get_1_dim_volume(fs, var_value_pairs, def_value, prec=NUM_BITS_PRECISION):
    """ In the setting of shifted lp-balls for p even, this return volume of the 
    deformed intersection of lp-balls ((Prod fs) - t) after intersection with the line defined by the 
    coordinate values for all but one coordinate.

    Parameters
    ----------
    fs : list
        List of polynomials in the same polynomial ring over QQ. 
    var_value_pairs : dict of (x_i : val_x_i) pairs.
        The QQ values of the x_i to which the fs should be restricted.
    def_value : QQ
        Value, by which the product of the fs will be deformed.
    prec : int
        Number of binary digits of precision.

    Returns
    ------
    RealBall
        Length of the line segment defined by var_value_pairs intersecting with the deformed intersection of the {f >= 0}.

    Caveat
    -------
        Assumes that the line actually intersects the deformed object and 
        that $t$ lies between 0 and the first critival value of the projection 
        from {(prod(fs) - t) = 0} to the t axis.

    TODOs
    -----
    - TODO: Raise error when variable appears twice. 
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


def get_inside_critical_points(fs, def_value, proj_var, var_value_pairs, prec=NUM_BITS_PRECISION, debug_level=0):
    """ Return the crit points of prod(fs)-t, evaluated at def_value and var_value_pairs, that lie inside {f > 0} for all f in fs.

    Parameters
    ---------
    fs              : list
        A list of multivariate polynomials.
    def_value       : QQ
        Value by which the product is deformed.
    proj_var        : variable
        The variable, such that we project onto the proj_var-axis.
    var_value_pairs : dict
        The values already restricted to.

    Returns
    -------
    inside_crit_points : list
        List of points of the form {xi:valxi ... for all i}, where i 
        goes over all variables that remain after the value substitutions.

    Remark
    ------
        The resulting computed critical points are those that lie on boundary of the smooth deformation of the intersection.
    """
    assert(def_value > 0) # The branch points would satisfy f==0 for some of the f in fs, and hence this case is not allowed.

    # Deform and compute branch points of the projection.
    fdef = tools.partial_eval_poly(tools.eval_poly(tools.deformed_product(fs), [def_value]), var_value_pairs) 

    try:
        crit_points = get_critical_points(fdef, proj_var, prec) 

        if debug_level > 0:
            print("Computed critical points for the proj_var = {} are: {}".format(proj_var, crit_points))

        # Now identify the critical points satisfying f > 0 for all f in fs:
        inside_crit_points = get_inside_points(fs, var_value_pairs, crit_points)
        
    except PosDimCritLocusError as error:
        print("Critical locus of projection is positive dimensional: {}".format(error))
        raise error
    
    return inside_crit_points

def get_critical_points(f, proj_var, prec=NUM_BITS_PRECISION):
    """ Computes all the real critical points.

    The critical points of the projection are the solutions to the ideal 
        I = (f) + ( diff(f, x_i) | x_i != proj_var).
    
    The solutions are computed using groebner bases and real root isolation in msolve.
    
    Parameters
    ----------
    f : Polynomial in QQ[x_1,...,x_n]
    proj_var : A variable out of f.parent().gens()

    Returns
    -------
    critical_points : list of dict
        The list of critical points represented by [{proj_var:value1}, {proj_var:value2},...]

    Raises
    ------
        PosDimCritLocusError if the critical locus is not 0 dimensional.
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

    Parameters
    ----------
    fs : list of polys in QQ[x_1,..., x_n]
        Concave polynomials defining the convex body.
    def_value : QQ > 0
        The deformation value.
    proj_var : in fs[0].parent().gens()
        Variable onto which Cdef shall be projected next.
    var_value_pairs : dict of (x_i : val_x_i) pairs.
        The QQ values of the x_i to which the fs should be restricted.
    prec : int
        Number of binary digits of precision.
    debug_level : int, optional
        Amount of intermediate results to be printed.

    Returns
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

    Parameters
    ----------
    fs : list of polys in QQ[x_1,..., x_n]
        Concave polynomials defining the convex body.
    def_value : QQ > 0
        The deformation value.
    var_value_pairs : dict of (x_i : val_x_i) pairs.
        The QQ values of the x_i to which the fs should be restricted.
    proj_var : in QQ[x_1,..., x_n]
        Variable that shall be projected onto next, i.e. which will not be used
        to construct the integrand.
    strategy : dict, optional
        A dictionary providing strategies, such as order of projections.
    debug_level : int, optional
        Amount of intermediate results to be printed.

    Returns
    -------
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
    intIdeal = creative_telescoping(annA, proj_var, debug_level=debug_level)

    return intIdeal.gens()[0]

def annihilator_deformed_intersection(fs, def_value, var_value_pairs, proj_var, debug_level=0):
    """ Returns an annihilating d-finite ideal for the rational integrand.

    The rational integrand is under standard considerations given by
    A = diff(F, x_i) * x_i / F      where x_i not = proj_var
    and  
    F = prod(fs) - def_value        after restricting to var_value_pairs.
    
    Normally, the integrand is constructed with construct_integrand().

    Parameters
    ----------
    fs : list of polys in QQ[x_1,..., x_n]
        Concave polynomials defining the convex body.
    def_value : QQ > 0
        The deformation value.
    var_value_pairs : dict of (x_i : val_x_i) pairs.
        The QQ values of the x_i to which the fs should be restricted.
    proj_var : in QQ[x_1,..., x_n]
        Variable that shall be projected onto next, i.e. which will not be used
        to construct the integrand.
    debug_level : int, optional
        Amount of intermediate results to be printed.
            
    Returns
    -------
    ideal
        A d-finite ideal in the rational Weyl algebra annihilating A.
    """

    A = construct_integrand(fs, def_value, var_value_pairs, proj_var, debug_level=debug_level)
    
    W = rational_weyl_algebra(A.parent().ring())

    annA = W.ideal([A*D-D(A) for D in W.gens()]) # construct the annihilating ideal
    
    return annA

def construct_integrand(fs, def_value, var_value_pairs, proj_var, prim_var_name=None, debug_level=0):
    """ Returns the rational function for which Vol(Cdef cap proj_var) is a period.
    
    By standard considerations, see for example our paper, the function 
    vol(proj_var) can be expressed as period of a rational function A. 
    We construct this function rational function here, considering 
    closely the order of indeterminate variables and the proj var.

    Parameters
    ----------
    fs : list of polys in QQ[x_1,..., x_n]
            Concave polynomials defining the convex body.
    def_value : QQ > 0
        The deformation value.
    var_value_pairs : dict of (x_i : val_x_i) pairs.
        The QQ values of the x_i to which the fs should be restricted.
    proj_var : in fs[0].parent().gens()
        Variable onto which Cdef shall be projected next.
    prim_var_name : str, optional
        Can specify a variable w.r.t. which the rational function is constructed.
    debug_level : int, optional
        Amount of intermediate results to be printed.

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

    # TODO: Here taking the first, should make this more strategic here.
    prim_var = [var for var in fdef_ZZ.parent().gens() if not var == proj_var_ZZ][0] 

    if prim_var_name != None:
        prim_var = fdef_ZZ.parent()(prim_var_name)

    sgn = (-1)**0 # TODO: Should depend on the choice of prim_var, to account for ordering.

    A = sgn * diff(fdef_ZZ, prim_var) * prim_var / fdef_ZZ # Automatically constructs the fraction field.

    return A


def get_picard_fuchs_t(fs, strategy=None, debug_level=0):
    """ Computes the Picard Fuchs operator for Vol(Ct).

    The functions in t>0 
    Vol(t) = Vol( prod(fs) - t >= 0) \cap {fs \geq 0 forall s} 
    is a period of a rational function and so there exists 
    a PF operator that annihilates. 

    Parameters
    ----------
    fs : list in QQ[x_1,..., x_n]
        List of concave polynomials defining semi algebraic set C.
    strategy : dict, optional
        [TODO] Not yet implemented
    debug_level : int, optional
            Amount of intermediate results to be printed.
    
    Returns
    -------
    UnivariateDifferentialOperatorOverUnivariateRing
        The PF operator annihilating Vol(t).

    TODOs
    -----
    - [TODO] Make the procedure by which to compute the intersection ideal more informed.
    """

    def_var_name = "t"

    # Construct the integrand, note that here this is over the ring QQ[x..., t], rather than QQ[x..][t]
    At = construct_integrand_t(fs, def_var_name, strategy)

    Wt = rational_weyl_algebra(At.parent().ring())

    annAt = Wt.ideal([At*D-D(At) for D in Wt.gens()]) # construct the annihilating ideal

    # To be precise, below we simply construct some subset of the integration ideal, 
    # but it suffices to be non-empty.
    intIdeal_t = creative_telescoping(annAt, At.parent().ring()(def_var_name), debug_level=debug_level)

    return intIdeal_t.gens()[0]


def construct_integrand_t(fs, def_var_name, strategy=None, prim_var_name=None):
    """Returns the rational function At for which Vol(Ct) is a period.
    
    By standard considerations, see for example our paper, the function
        vol(t) = vol({prod(fs)-t >= 0} \cap {fs >= 0}) 
    can be expressed as period of a rational function At.


    Parameters
    ----------
    fs : list of polys in QQ[x_1,..., x_n]
        Concave polynomials defining the convex body.
    def_var_name : str
        Variable name of the deformation parameters, e.g. "t".
    strategy : dict, optional
        [TODO] Not implemented here.
    prim_var_name : str, optional
        Can specify a variable w.r.t. which the rational function is constructed.

    Output
    -----
    At : Element of a FractionField.
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

def creative_telescoping(I, proj_var, strategy=None, debug_level=0):
    """ Integrates out all but the variable proj_var from the R or D ideal I.

    The result is an ideal contained in:     (I + dx_1 * D + ...+ dx_n *D) \cap D_x0
    (In this example we assume x0 to be the proj_var)

    Parameters
    ----------
    I : ideal of RationalWeylAlgebra
        Ideal in Weyl algebra over field of rational functions.
    proj_var : in I.parent().base_ring().gens()
        The variable that should remains after creative telescoping.
    strategy : dict, optional
        [TODO] Not implemented yet.
    debug_level : int, optional
            Amount of intermediate results to be printed.
    
    Output
    ------
    ideal
        Returns a list of of univariate differential operators of length 1 (!).

    TODOs
    -----
    - [TODO] Should extend this, so that at least temporarily we store / output for debuginfo the certificates, too.
    - [TODO] Should include some basic strategy.
    """

    if debug_level > 0:
        print("[CT] Entering creative telescoping:")

    W = I.ring()
    ct_system = list(I.gens())

    if debug_level > 2:
        print("\n[CT] The ideal is generated by: {}".format(list(I.gens())))

    for Dvar in [Dvar for Dvar in W.gens() if not Dvar == W("D"+str(proj_var))]:
        ct_ideal = ct_system[0].parent().ideal(ct_system)
        
        if debug_level > 0:
            print("\n[CT] Integrating out the variable {}".format(str(Dvar)))

        ct_system = ct_ideal.ct(Dvar, certificates=False)

    return ct_system[0].parent().ideal(ct_system)


def solve_diff_op(P, initial_conditions, evaluation_condition, prec, apparent_singular_points=[], debug_level=0):
    """ Solves a linear univariate differential operator with specified initial and evaluation conditions.

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

    Parameters
    ----------
    P : UnivariateDifferentialOperatorOverUnivariateRing
        PF operator to be solved.
    initial_conditions : dict
        E.g. { x_val_1:{"exponent":expon, "coef":coef }, x_val_2:...}
    evaluation_condition : dict
        Dictionary containing evaluation point and the exponent of the 
        starting monomial, whose coef shall be read of, e.g. {"pt":0, "exponent":0}
    prec : int
        Number of binary digits of precision.
    apparent_singular_points : list, optional
        List of points through which analytic continuation from the initial points 
        to the evaluation pt should always go through (i.e. they are only apparent
        singular points of P)
    debug_level : int, optional
            Amount of intermediate results to be printed.
           
    Returns
    -------
    sage.rings.complex_arb.ComplexBall
        The solution of P according to initial_conditions and eval_conditions up to precision 2^{-prec}.

    Example
    -------
    initial_condition = {1/10:{"exponent":1, "coef":3.12223532...}, }
    evaluation_condition = {"pt":0, "exponent":0}
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
    """ Computes the volume of the smooth deformations of semi-algebraic convex bodies.

    Assumes that convex body is given as
        C = {x in RR^n | f(x) > 0 for all f in fs}
    where 
        f is concave, for all f in fs.

    The deformation then is
        Cdef =  C \cap {x \in RR^n | prod(fs) - def_value > 0}

    The var_value_pairs, additionally restrict this to
        { x_i = val | (x_i, val) in var_value_pairs}
    
    Parameters
    ----------
    fs : list of polys in QQ[x_1,..., x_n]
        Concave polynomials defining the convex body.
    def_value : QQ > 0
        The deformation value.
    var_value_pairs : dict of (x_i : val_x_i) pairs.
        The QQ values of the x_i to which the fs should be restricted.
    prec : int
        Number of binary digits of precision.
    strategy : dict, optional
        A dictionary providing strategies, such as order of projections.
    debug_level : int, optional
        Amount of intermediate results to be printed.
            
    Returns
    -------
    sage.rings.complex_arb.ComplexBall
        The volume of the deformed convex body Cdef with precision 2^{-prec}.

    Caveat
    ------
    It is assumed that the convex body C, respectively the slice defined by 
    var_value_paris is full-dimensional.

    Concavity is NOT explicitly checked for.
    
    Remark
    ------
    This is a wrapper for the SmoothVolume class.
    It constructs a SmoothVolume object here and return output of get_volume().


    TODOs
    ------
    [TODO] Check that the initial points are actually good before evaluating (or add option to do so.)
    [TODO] Assert that there is actually only a unique root in the singular locus corresponding to either numerical branch_point.
    [TODO] Parallelize computation of slices.
    """

    volObject = SmoothVolume(fs=fs,def_value=def_value,var_value_pairs=var_value_pairs,prec=prec, strategy=strategy,debug_level=debug_level)

    # Return the volume
    return volObject.get_volume()


def get_good_initial_points(P, x0, x1, prec, debug_level=0):
    """ Provides order(P) many rational points strictly between x0 and x1.

    Whether the conditions are good is checked up to prec. (i.e. bad if 0 \in CBF(determinant))

    Parameters
    ----------
    P : UnivariateDifferentialOperatorOverUnivariateRing
        PF for which good initial points shall be found.
    x0 : RealBall or int or QQ
    x1 : RealBall or int or QQ
    prec : int
            Number of binary digits of precision.
    debug_level : int, optional
        Amount of intermediate results to be printed.
    
    Returns
    ----------
    list of QQ
        List of good initial points for P.

    TODOs
    ------
    [TODO] Actually provide the invertibility check of the initial points.  Perturb if needed.
    """
    
    n = P.order()
    q = sample_n_rational_points(x0,x1,n, debug_level=debug_level)

    # TODO: Actually check invertibility of the linear system.
    # and then slightly perturb the points if not invertible.

    return q

def sample_n_rational_points(x0,x1,n, base=10, debug_level=0):
    """ Linearly samples n rational points strictly between x0 and x1.

    Determines q_start and delta, accordingly, so that 
        qk = qstart + k*delta (for k = 0,..., n-1)
    yields n rational points between x0 and x1.

    Parameters
    ----------
    x0 : RealBall or int or QQ
    x1 : RealBall or int or QQ
    n : int
        Number of points to be sampled between x0 and x1.
    base : int, optional
        Base in which the logarithm is computed.
    debug_level : int, optional
        Amount of intermediate results to be printed.
    
    Returns
    -------
    list of QQ
        List of n rational points strictly between x0 and x1.

    TODOs
    -----
    [TODO] Make sure log and floor always work as expected.
    """

    if debug_level > 0:
        print("[SamplePoints] Sampling n = {} many points between x0 = {} and x1 = {}.".format(n, x0, x1))
        print("[SamplePoints] Chosen base = {}.".format(base))

    xmax = max(x0,x1)
    xmin = min(x0,x1)

    # Determine integer N, such that 10^N < xmax-xmin 
    N = np.floor(log(xmax-xmin)/log(base))

    # Choose qstart in 2*10^(N-1) neigborhood of xmin.
    qstart = QQ((np.floor(xmin / power(QQ(10), N-1)) + 2) * power(QQ(10), N-1))

    # Generate the samples linearly
    if n == 1:
        return [qstart]

    # n != 1:
    delta = QQ(10**(N-1) / (n-1))
    q = [qstart + k*delta for k in range(0,n)]

    if debug_level > 1:
        print("[SamplePoints] qstart = {}".format(qstart))
        print("[SamplePoints] delta = {}".format(delta))

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
    strategy : dict, optional
        A dictionary providing strategies, such as order of projections.
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
    initial_points = get_good_initial_points(Pt, 0, CBF(smallest_positive_singular_point).real(), prec)

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