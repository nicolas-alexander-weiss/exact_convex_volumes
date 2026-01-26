#
# @author: Nicolas Weiss
# @created: 2025-10-31
#
# @goal: Interact with msolve from sage, retrieving interval arithmetic results with guaranteed precision.
#
# @python-version: > 3.4
#
# TODO: The assertions should really be turned into errors that are more informative!
#

# Python Imports

import subprocess
import os
from datetime import datetime

# Sage Imports

from sage.all import (
                        
    QQ,
    diff,
    sage_eval,
    RealBallField, 
    RealIntervalField
)



def branch_points_system(f, proj_var):
    """ The branchpoints of the projection onto the proj_var-axis
    are the solutions to the system:

    f = 0,
    diff(f, xi) for xi not proj_var

    Input
    ------
    f : Element of polynomial ring over QQ.
    proj_var : Element of f.parent().gens()

    Output
    ------
    The list of generators to the ideal.

    Caveat
    ------
    - Requires that f lives in a polynomial ring over QQ. For example (QQ[x])[t] is not allowed, 
    but has to be flattened first.
    """

    assert(f.parent().base_ring() == QQ)

    return [f] + [diff(f, var) for var in f.parent().gens() if not (var == proj_var)]


def write_msolve_file(system, file_name, characteristic=0):
    """ Writes an msolve input file for the system of polynomials in specified in system.

    Input
    ------
    system : List of polynomials over the same polynomial ring over QQ
    file_name : String
    characteristic : Natural Number >= 0, this is used as input to msolve, too.

    """

    # Check that base_ring is QQ and that all polynomials lie in the same ring.
    assert(system[0].base_ring() == QQ)
    assert(all([poly.base_ring() == system[0].base_ring() for poly in system]))

    poly_ring = system[0].parent()

    # Write to the output file
    with open(file_name, "w") as out_file:
        # Write the variable names as the first row
        out_file.write(",".join([str(var) for var in poly_ring.gens()]))
        out_file.write("\n")

        # Write the characteristic of the field as the second row
        out_file.write(str(characteristic))
        out_file.write("\n")

        # Next write all the polynomials, 1 per row, ended by a "," (except last)
        out_file.write(",\n".join([str(poly) for poly in system]))


def read_msolve_output(file_name):
    """ Reading in the output file of msolve.

    Note: The file ends with a ":" (then we can stop)

    Successful Example File structure: (annotated here by the line number)
    line 1:         [0, [1,
    line 2:         [[[-273881754523108621558138698965051022467 / 2^128, -136940877261554310779069349482525511233 / 2^127], [127701611743269739908586015828574437907 / 2^126, 510806446973078959634344063314297751629 / 2^128], [12507467247990595072232567913918321782687 / 2^127, 25014934495981190144465135827836643565375 / 2^128]], [[127701611743269739908586015828574437907 / 2^126, 510806446973078959634344063314297751629 / 2^128], [-273881754523108621558138698965051022467 / 2^128, -136940877261554310779069349482525511233 / 2^127], [12507467247990595072232567913918321782687 / 2^127, 25014934495981190144465135827836643565375 / 2^128]], [[290016790287279909900914490826241735255 / 2^127, 580033580574559819801828981652483470511 / 2^128], [290016790287279909900914490826241735255 / 2^127, 580033580574559819801828981652483470511 / 2^128], [10201656202392956683560662369959977360999 / 2^128, 1275207025299119585445082796244997170125 / 2^125]], [[170141183460469231731687303715884105727 / 2^128, 170141183460469231731687303715884105729 / 2^128], [21097451832125581058166014184115227427 / 2^123, 675118458628018593861312453891687277665 / 2^128], [-460820905525367823493297946583722440723 / 2^382, 460820905525367823493297946441493312999 / 2^382]], [[188851325292295987614616248873271680739 / 2^128, 47212831323073996903654062218317920185 / 2^126], [188851325292295987614616248873271680739 / 2^128, 47212831323073996903654062218317920185 / 2^126], [117565802873507594341607564620084590397 / 2^127, 235131605747015188683215129240169180795 / 2^128]], [[1, 1], [1, 1], [0, 0]], [[21097451832125581058166014184115227427 / 2^123, 675118458628018593861312453891687277665 / 2^128], [170141183460469231731687303715884105727 / 2^128, 170141183460469231731687303715884105729 / 2^128], [-690527393910386501608815396611790428336503 / 2^381, 5524219151283092012870523172851497613948983 / 2^384]], [[5446275213858333065436088620978774959 / 2^128, 340392200866145816589839582669969721 / 2^124], [10633823966279326983230393290174834703 / 2^124, 10633823966279326983230519674310678513 / 2^124], [-662737094745247 / 2^128, 662737094745247 / 2^128]], [[85070591730234615865843651856684484641 / 2^127, 85070591730234615865843651859199621087 / 2^127], [85098050216536454147449390164437713 / 2^122, 5446275213858333065436760973174276863 / 2^128], [-1313364001 / 2^128, 1313364001 / 2^128]], [[16395153586771595646955130527489217495 / 2^127, 16395153586771595646955130527489250265 / 2^127], [32790307173543191293910261054978439451 / 2^128, 8197576793385797823477565263744624017 / 2^126], [-12596525515723031973012167686913679381 / 2^128, -12596525515723031973012167686913632969 / 2^128]], [[0, 0], [0, 0], [0, 0]], [[83444594588840401951197726375664784343 / 2^127, 166889189177680803902395452751329568687 / 2^128], [559376840048630554274471107371728468471 / 2^128, 69922105006078819284308888421466058559 / 2^125], [-1270106343953908409925400429977976753977 / 2^128, -158763292994238551240675053747247094247 / 2^125]], [[559376840048630554274471107371728468471 / 2^128, 69922105006078819284308888421466058559 / 2^125], [83444594588840401951197726375664784343 / 2^127, 166889189177680803902395452751329568687 / 2^128], [-1270106343953908409925400429977976753977 / 2^128, -158763292994238551240675053747247094247 / 2^125]]]
    line 3:         ]]:

    joined (colon removed):     [0, [1, [[[-273881754523108621558138698965051022467 / 2^128, -136940877261554310779069349482525511233 / 2^127], [127701611743269739908586015828574437907 / 2^126, 510806446973078959634344063314297751629 / 2^128], [12507467247990595072232567913918321782687 / 2^127, 25014934495981190144465135827836643565375 / 2^128]], [[127701611743269739908586015828574437907 / 2^126, 510806446973078959634344063314297751629 / 2^128], [-273881754523108621558138698965051022467 / 2^128, -136940877261554310779069349482525511233 / 2^127], [12507467247990595072232567913918321782687 / 2^127, 25014934495981190144465135827836643565375 / 2^128]], [[290016790287279909900914490826241735255 / 2^127, 580033580574559819801828981652483470511 / 2^128], [290016790287279909900914490826241735255 / 2^127, 580033580574559819801828981652483470511 / 2^128], [10201656202392956683560662369959977360999 / 2^128, 1275207025299119585445082796244997170125 / 2^125]], [[170141183460469231731687303715884105727 / 2^128, 170141183460469231731687303715884105729 / 2^128], [21097451832125581058166014184115227427 / 2^123, 675118458628018593861312453891687277665 / 2^128], [-460820905525367823493297946583722440723 / 2^382, 460820905525367823493297946441493312999 / 2^382]], [[188851325292295987614616248873271680739 / 2^128, 47212831323073996903654062218317920185 / 2^126], [188851325292295987614616248873271680739 / 2^128, 47212831323073996903654062218317920185 / 2^126], [117565802873507594341607564620084590397 / 2^127, 235131605747015188683215129240169180795 / 2^128]], [[1, 1], [1, 1], [0, 0]], [[21097451832125581058166014184115227427 / 2^123, 675118458628018593861312453891687277665 / 2^128], [170141183460469231731687303715884105727 / 2^128, 170141183460469231731687303715884105729 / 2^128], [-690527393910386501608815396611790428336503 / 2^381, 5524219151283092012870523172851497613948983 / 2^384]], [[5446275213858333065436088620978774959 / 2^128, 340392200866145816589839582669969721 / 2^124], [10633823966279326983230393290174834703 / 2^124, 10633823966279326983230519674310678513 / 2^124], [-662737094745247 / 2^128, 662737094745247 / 2^128]], [[85070591730234615865843651856684484641 / 2^127, 85070591730234615865843651859199621087 / 2^127], [85098050216536454147449390164437713 / 2^122, 5446275213858333065436760973174276863 / 2^128], [-1313364001 / 2^128, 1313364001 / 2^128]], [[16395153586771595646955130527489217495 / 2^127, 16395153586771595646955130527489250265 / 2^127], [32790307173543191293910261054978439451 / 2^128, 8197576793385797823477565263744624017 / 2^126], [-12596525515723031973012167686913679381 / 2^128, -12596525515723031973012167686913632969 / 2^128]], [[0, 0], [0, 0], [0, 0]], [[83444594588840401951197726375664784343 / 2^127, 166889189177680803902395452751329568687 / 2^128], [559376840048630554274471107371728468471 / 2^128, 69922105006078819284308888421466058559 / 2^125], [-1270106343953908409925400429977976753977 / 2^128, -158763292994238551240675053747247094247 / 2^125]], [[559376840048630554274471107371728468471 / 2^128, 69922105006078819284308888421466058559 / 2^125], [83444594588840401951197726375664784343 / 2^127, 166889189177680803902395452751329568687 / 2^128], [-1270106343953908409925400429977976753977 / 2^128, -158763292994238551240675053747247094247 / 2^125]]]]]

    Failed Example (i.e. positive dimensional solution)
    [TODO]
    """

    with open(file_name, "r") as in_file:
        # Read line by line and concatenate.  EllipticCurve
        # Note: If at end of file, readline will return an empty string.

        line1 = in_file.readline().strip()
        line2 = in_file.readline().strip()
        line3 = in_file.readline().strip()

        # Assert the common structure
        assert(line1[0] == "[")
        assert(line1[-1] == ",")
        assert(line2[:2] == "[[")
        assert(line3[:3] == "]]:")

        # Concatenate and evaluate as sage literal. The coordinates, should be correctly parsed in QQ.

        return sage_eval(line1 + line2 + line3[:-1])
        
def extract_msolve_solutions(msolve_out_object, vars, real_ball_field):
    """ Read the solution out of the list of lists.

    Input
    -----
    msolve_out_object : The return of read_msolve_output
    vars : The variables to which the computed coordinates will be matched.
    real_ball_field : Interval arithmetic field in which output should be provided.

    Output
    ------
    points : list of dicts {x0:124.2523..., x1:-2.023..., x2:0.22211123.... }, indexed by the variables.
    """
    # Create a shorter handle
    ms_out = msolve_out_object

    # Check that the output is a 0-dimensional variety.
    assert(ms_out[0] == 0)
    
    # Check that the next list has two elements
    assert(len(ms_out[1]) == 2)
    num_solution_lists = ms_out[1][0]
    solution_lists = ms_out[1][1:]

    # Assert that as many coordinates as variables are provided
    assert(len(vars) == len(solution_lists[0][0]))

  
    points = []
    RI = RealIntervalField(real_ball_field.prec())
    for solution_list in solution_lists:
        for solution in solution_list:
            # Turn the intervals into elements of the real ball field.
            sol = {
                    var:real_ball_field(RI(*interval))  for var,interval in zip(vars,solution)
            }
            points.append(sol)

    return points

def variety_msolve(system, prec):
    """ Compute the real solutions of a 0-dimensional system via msolve with a required precision.

    Input
    -----
    system : a list of polynomials in a polynomial ring over QQ. 
    prec : number of binary digits of precision.

    Output
    -----
    A variety encode as a list of points [pt, pt, pt], where each point is a dictionary {x:.., y:.., z:..} with 
    keys being the variable names and values in the real ball field with 200 bits precision.
    """
    # Assert that all polynomials are from the same ring and over QQ.
    R = system[0].parent()

    assert(all([R == poly.parent() for poly in system]))
    assert(R.base_ring() == QQ)

    # Create a msolve_tmp folder if doesn't exist.
    msolve_tmp_folder_name = "msolve_tmp"
    if not os.path.isdir(msolve_tmp_folder_name):
        os.mkdir(msolve_tmp_folder_name)

    # Determine a filename (like timestamp)
    current_time = datetime.now() # TODO: Also should check if there already is a file with that time stamp...
    ms_prefix = current_time.strftime("%Y%m%d%H%M%S")

    # Write the msolve input file 
    in_file_path = os.path.join(msolve_tmp_folder_name, ms_prefix + "-in.ms")
    write_msolve_file(system, in_file_path)
    # Execute msolve
    out_file_path = os.path.join(msolve_tmp_folder_name, ms_prefix + "-out.ms")
    subprocess.run(["msolve", "-f", in_file_path, "-o", out_file_path, "-p", str(prec)])

    # Read the result
    msolve_obj = read_msolve_output(out_file_path)
    
    # Extract solutions
    RB = RealBallField(prec)
    sols = extract_msolve_solutions(msolve_obj, R.gens(), RB)
    
    return sols


# if __name__ == "__main__":
#     S = PolynomialRing(QQ, "x", 2)
#     p = 4
#     mus = [(1,0), (0,1), (1,1)]

#     fs = [shifted_lp_poly(S, p, mu) for mu in mus]

#     ft = deformed_product(fs, "t")

#     ft_flattened = ft.parent().flattening_morphism()(ft)

#     S.inject_variables()
#     ft.parent().inject_variables()

#     proj_var = t

#     branch_points_system = branch_point_system(ft_flattened, t)

#     in_file_name = "ms-files/branch_system_t.msolve"

#     write_msolve_file(branch_points_system, in_file_name)

#     out_file_name = "ms-files/out.msolve"
#     # Execute the system

#     prec = 200
#     subprocess.run(["msolve", "-f", in_file_name, "-o", out_file_name, "-p", str(prec)])

#     msolve_obj = read_msolve_output(out_file_name)

#     RB = RealBallField(prec)

#     sols = extract_msolve_solutions(msolve_obj, ft_flattened.parent().gens(), RB)

