sys.path.insert(0, '/Users/lakshmiramesh/Desktop/ore_algebra/src')
import ore_algebra
from ore_algebra import *
def volume_intersection(dim, p, translation_vectors): #assume that dim > 1, p>1 even, translation vectors give non-empty intersection!
    k = len(translation_vectors)
    Rt = PolynomialRing(ZZ, dim+1, 'x') #note that xdim+1 is t
    f = Rt(1)
    Ft = Rt.fraction_field()    
    Wt = OreAlgebra(Ft, *[("Dx" + str(i), {}, {"x" + str(i) : 1}) for i in range(0, dim+1, 1)])
    
    for i in range(0,k+1): #make the polynomials
        this_poly = 1
        for j in range(0,dim):
            this_poly = this_poly - Rt.gens()[i]^(p)
        f = ft*this_poly
    f = f - Rt.gens()[dim+1]
    At = Ft(y*ft.derivative(Rt.gens()[0])/ft) #differentiat by x0. this means we always telescope wrt x1 first.  
    annt = Wt.ideal([At*D-D(At) for D in Wt.gens()])
    ct0 = annt.ct(Wt.gens()[0], certificates=False)
    for i in range(1:dim): #telescope to get Picard-Fuchs in xdim+1 
        ct0 = ct0[0].parent().ideal(ct0).ct(Wt.gens()[i], certificates=False)
    order0 = ct0[0].order()
    #CHOOSE SLICES 
    max0 = QQ[Rt.gens()[dim+1]](ct0[0].leading_coefficient()).roots(QQbar, multiplicities=False)
    #TODO : Find the t0 smallest positive, real critical value and choose order0-1 many points between 0 and  t0
    #L0 = [...] values of x0 for slices!

    #now we compute order0 - 1 volumes of x0-slices. 
    x0_volumes = []
    for i in range(0,7):
        #smooth algorithm begins here
        j=0
        something = smooth_volume(j,L[i],f,dim)
        x0_volumes.append(something)
    #analytically continue to get volume at x0 = 0
    return final_volume


def smooth_volume(j, point, f1, dim):
    R = PolynomialRing(ZZ, dim-j, 'x') #Note that now, x0 = x, x1 = y, ...
    F = R.fraction_field()    
    W = OreAlgebra(F, *[("Dx" + str(i), {}, {"x" + str(i) : 1}) for i in range(0, dim-j, 1)])
    f = F(f1(R.gens()[j] = QQ(point))) #i hope this is correct
    A = F(y*f.derivative(R.gens()[0])/f)
    ann = W.ideal([A*D-D(A) for D in W.gens()])
    ct = ann.ct(W.gens()[0], certificates=False)
    if j == dim-1:
        #do nothing
    else:
        for k in range(1,dim-j):
            ct = ct[0].parent().ideal(ct).ct(W.gens()[k], certificates=False)
    L = QQ[R.gens()[dim-j]](ct[0].leading_coefficient()).roots(QQbar, multiplicities=False)
    #kick out unwanted critical points
    initial_conditions = []
    for i in len(L):
        if j == dim-1: 
            #actually compute a volume 
        else:
            somthing = smooth_volume(j+1, L[i], f, dim)
            initial_conditions.append(somthing)
    #use the initial conditions to compute the volume of this slice. 
    return volume        
    