# sys.path.insert(0, '/Users/lakshmiramesh/Desktop/ore_algebra/src')
import ore_algebra
from ore_algebra import *


def shifted_lp_poly(Rt, p, mu):

    lp_poly = Rt(1)
    for i in range(0,len(mu)):
        lp_poly = lp_poly - (Rt.gens()[i] - mu[i])^(p)
    
    return lp_poly
        


def volume_intersection(dim, p, translation_vectors, precis): #assume that dim > 1, p>1 even, translation vectors give non-empty intersection!
    """ Computes the volume of the intersection of L_p balls shifted by translation_vectors in RR^dim. 

    Input
    -----
    dim : integer dimension of ambient space
    p : integer, assumed even.
    translation_vectors : list of vectors in QQ^dim
    precis : integer precision digits

    Output
    ------
    The volume of the intersection of the translated L_p balls as a fixed precision real number.

    Example
    -----
    
    dim = 2
    p = 4
    translation_vectors = [(0,0), (1,0)]
    
    vol = volume_intersection(dim,p,translation_vectors)  # 1.71448285904485523624162985240570841916093...

    """
    
    k = len(translation_vectors)
    Rt = PolynomialRing(ZZ, 'x', dim+1) # note that xdim+1 is t, i.e. Rt.gens()[-1] or Rt.gens()[xdim]
    Ft = Rt.fraction_field()    
    Wt = OreAlgebra(Ft, *[("Dx" + str(i), {}, {"x" + str(i) : 1}) for i in range(0, dim+1, 1)])
    
    # create complex ring for creative telescoping and real ring for precision 
    # 

    # make the polynomial ft = f1*...*fk - t
    ft = Rt(1)
    for mu in translation_vectors: 
        ft = ft * shifted_lp_poly(Rt, p, mu)

    ft = ft - Rt.gens()[-1]

    At = Ft(y*ft.derivative(Rt.gens()[0])/ft) # differentiat by x0. this means we always telescope wrt x0 first.  
    annt = Wt.ideal([At*D-D(At) for D in Wt.gens()])

    # telescope to get Picard-Fuchs in xdim+1  (i.e. eliminate all variables but the last one)
    ct0 = annt.ct(Wt.gens()[0], certificates=False)
    for i in range(1, dim): 
        ct0 = ct0[0].parent().ideal(ct0).ct(Wt.gens()[i], certificates=False)
    order0 = ct0[0].order()

    # CHOOSE SLICES 
    max0 = QQ[Rt.gens()[dim+1]](ct0[0].leading_coefficient()).roots(RB, multiplicities=False)
    # TODO : Find the t0 smallest positive, real critical value and choose order0-1 many points between 0 and  t0
    # L0 = [...] values of x0 for slices!

    # now we compute order0 - 1 volumes of x0-slices. 
    x0_volumes = []
    for i in range(0,7):
        # smooth algorithm begins here
        j=0
        something = smooth_volume(j,L[i],f,dim)
        x0_volumes.append(something)
    # analytically continue to get volume at x0 = 0
    return final_volume


def smooth_volume(j, point, f1, dim, RB):
    R = PolynomialRing(ZZ, dim-j, 'x') #Note that now, x0 = x, x1 = y, ...
    F = R.fraction_field()    
    W = OreAlgebra(F, *[("Dx" + str(i), {}, {"x" + str(i) : 1}) for i in range(0, dim-j, 1)])
    # TODO
    # f = None
    f = F(f1(R.gens()[dim - j] = QQ(point))) 
    A = F(y*f.derivative(R.gens()[0])/f)
    ann = W.ideal([A*D-D(A) for D in W.gens()])
    ct = ann.ct(W.gens()[0], certificates=False)
    if j == dim-2:
        # do nothing
        pass
    else:
        for k in range(1,dim-j-1):
            ct = ct[0].parent().ideal(ct).ct(W.gens()[k], certificates=False)
    L = QQ[R.gens()[dim-j]](ct[0].leading_coefficient()).roots(RB, multiplicities=False)
    # TODO kick out unwanted critical points
    initial_conditions = []
    for i in len(L):
        if j == dim-2: 
            # TODO actually compute a volume 
            # here we must use our conjecture.
            f = f(R.gens()[dim - j] = L[i])
            # we select the middle two points
            roots = (QQ[R.gens()[0]](f)).roots(RB, multiplicities=False)
            len_roots = len(roots) + 1 # this should be even
            slice_vol = roots[len_roots/2] - roots[len_roots/2 - 1]
        else:
            slice_vol = smooth_volume(j+1, L[i], f, dim)
            initial_conditions.append(slice_vol)

    # TODO use the initial conditions to compute the volume of this slice. 
    return volume        
    