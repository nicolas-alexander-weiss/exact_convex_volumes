"""
@author: Nicolas Weiss
@goal: Providing a little wrapper to access elimination in M2 from sage (based on the sage M2 interface).

"""

# Python Imports

# Sage Imports

from sage.all import (
    macaulay2,
    QQ
)



def eliminate(I, varNames, useM2=False):
    """
    Input
    -----
    I : ideal in QQ[x0, ..., xn]
    varNames : string of variables to be eliminated

    Output
    ------
    Eliminate varNames from I using M2 and returns the result within the same ring as I.

    Example
    -----

    var('x,y,z')
    S.<x,y,z> = PolynomialRing(QQ)

    I = S.ideal([x^2 + y^2 + z^2, 2*y, 2*z])

    eliminate(I, 'y,z')

    """

    if not type(varNames) == str:
        raise TypeError("varNames excepted to be of type str")

    if useM2:
        m2ideal = macaulay2(I)
        eliminatedIdeal = m2ideal.eliminate("{" + varNames + "}") 
        return eliminatedIdeal.sage()
    else:
        # Write code to do elimination in Sage.
        return I.elimination_ideal([I.ring()(str(var_name.strip())) for var_name in varNames.split(",")])

# TODO: If f lives in a ring of the form $QQ(t)[x0,...,x_n]$ we have to manually construct the M2 object.
# TODO: It would be best to specifiy if we want to keep certain variables as parametric
def branchIdeal(f, projectedVar):
    """
    Input
    -----
    f : Polynomial in QQ[x0, ..., xn]
    projectedVar : variable onto which to project

    Output
    ------
    Returns the ideal of branch points (critical values of the projection) 
    in QQ[projectedVar].

    Example
    -------
    var('x,y,z')
    S.<x,y,z> = PolynomialRing(QQ)

    f = x^2 + y^2 + z^2

    I = branchIdeal(f,x); I 
    """
    R = f.parent()

    # Currently, this handles R = QQ[x0,..., xn], but not base_ring being a fraction field.
    if not R.base_ring() == QQ:
        raise TypeError("Currently can handle only polynomial f over QQ")

    fiberVars = [g for g in R.gens() if g != projectedVar]; 

    # Construct ideal of branch points in V(f)
    I = R.ideal([f] + [f.derivative(v) for v in fiberVars])

    eliminatedIdeal = eliminate(I, ",".join(map(str, fiberVars)))

    return (R.base_ring()[projectedVar]).ideal(eliminatedIdeal)