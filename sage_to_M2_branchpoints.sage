#
# Given the input of generators of an ideal in sage, 
# compute the elimination with regards to specific variables in M2
# 

var('x,y')

f = x^2 + y^2

dyf = derivative(f, y)

R = macaulay2('QQ[x,y]')

I = macaulay2.ideal( ('y^2 - x^3', 'x - y') ); I 


