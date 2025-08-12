-- Computing branchpoints 

--- Our deformed polynomial

R = QQ[t,x,y,z,w]

f1 = 1 - (x^4 + y^4 + z^4 + w^4)
f2 = 1 - ((x-1)^4 + y^4 + z^4 + w^4)

ft = f1*f2 -t


-- Compute the branchpoints for the projection on z, for fixed t

S = QQ[x,y,w]

t0 = 0 -- 1/1000
z0 = 1/10
-- w0 = 1/10
phi = map(S,R,{t0, x, y, z0, w})

-- Compute the deformed intersection
ft0 = phi(ft) 



-- Compute the branch points of the projection onto the y-axis (i.e. in the V(ft0))

-- Cond 1: ft0 vanishes
-- Cond 2: dx, dw  derivative of ft0 vanishing.

I = ideal (ft0, diff(x, ft0), diff(w, ft0))

-- Now project onto the y-axis, i.e. eliminate

branchideal = eliminate({x,w}, I)

branchpoly = sub((gens branchideal)_0_0, QQ[x])

branchpoints = roots(branchpoly, Precision=>50, Unique=>true)