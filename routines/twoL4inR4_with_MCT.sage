# TODO: Compute the volume of the intersection of two L4 balls in R4 with centers 0 and (1,0,0,0) 
#   But compute the telescopers usign Julia.

# Output: [4.37856654871924288558479945106166817027285180219207050298477722007062424113388554779398481745775 +/- 9.11e-96] + [+/- 2.72e-96]*I
# Time: 15h
# Remark: Only worked with 500 bits target/working precision.

#
# The numerical experiment with is here:


from fractions import Fraction
import numpy as np

# Setup
n = 4
p = 4

# Number of samples and batch_size
num_samples = 200*1000*1000
batch_size =  1000000

# Region in which to sample:
radius_box = 1
center_box = 0

# Indicator function of the lp unit ball:
ind_funs = [lambda x_vec, mu=mu : np.sum(np.power(x_vec - np.array(mu), p), axis=1) <= 1 for mu in [[0,0,0,0],[1,0,0,0]]]

# Setup a random number generator:
seed = int(12345)
rng = np.random.default_rng(seed)

num_within_support = 0
for i in range(0, num_samples, batch_size):
    num_samples_batch = batch_size if (i + batch_size < num_samples) else num_samples - i 
    
    # Sample points
    sampled_points_batch = rng.uniform(-radius_box, radius_box, (num_samples_batch, n)) + center_box

    # Cancel out function by function the samples that are not in their support
    in_the_support_batch = np.full((num_samples_batch,), True)
    for ind_fun in ind_funs:
        # Check entries if previously or now not in the support of ind_fun
        in_the_support_batch &= ind_fun(sampled_points_batch)

    # apply indicator function to each entry and count 
    num_within_support_batch = np.count_nonzero(in_the_support_batch) 
    num_within_support += num_within_support_batch

# Return ratio of volume of box that corresponds to the support.
volume_mc =  num_within_support / num_samples * np.power(2 * radius_box, n)

print("Monte-Carlo approximated volume of Lp cap Lp (p={}) unit ball in R^n (n={}): {:.100}".format(p, n, float(volume_mc)))






######################

# With installation:
# import exact_convex_volumes

# Without installation:
from context import exact_convex_volumes

# In any case:
from exact_convex_volumes.volumes import *
from exact_convex_volumes.tools import *
from exact_convex_volumes.msolve_interface import *

######################

debug_level = 3

n = 4
p = 4
R = PolynomialRing(QQ, "x", n)
x = R.gens()

fs = [shifted_lp_poly(R, p, [0,0,0,0]), shifted_lp_poly(R, p, [1,0,0,0])]

def last_variable_proj_var_strategy(fs, deform_value, var_value_pairs):
    evaluated_poly = partial_eval_poly(eval_poly(deformed_product(fs), [deform_value]), var_value_pairs)

    return evaluated_poly.parent().gens()[-1]

strategy={"proj_var":last_variable_proj_var_strategy}


prec = 500  # Important. Accuracy should be high enough! (last time it faded!)

# Now essentially manually go through the volume 2 function.
vol = Volume(fs, prec, strategy=strategy, debug_level=3, use_julia_for_CT=True)

%time vol.start_computation()

print(vol.vol)