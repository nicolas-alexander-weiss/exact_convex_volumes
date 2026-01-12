""" 
@created: 2025-04-07
@author: Nicolas Weiss

@goal: Implement various Monte Carlo based methods to compute volumes of semi-algebraic convex bodies.

@usage: Activate the environment (if setup) (in Unix / Apple) with "source .venv/bin/activate". 
        The files can then be run with "python3 montecarlo_approach.py"
 
"""

from fractions import Fraction
import numpy as np




MAX_BATCH_SIZE = 100000

STD_NUM_SAMPLES = 1000*1000
STD_BATCH_SIZE = 1000*100
STD_SEED = 12345


#
# (vectorized) helper functions:
#
def p_norm_squared(x_vec, p):
    """
    Computes the square of the p-norm for each entry in parallel.

    Assumes that p is even.
    Assumes that x_vec a (num_points, n) shaped numpy array.
    """
    assert(p % 2 == 0)
    return np.sum(np.power(x_vec, p), axis=1)

#
# Key function:
#

# For now just writing the batched function, which can be used to define the other ones.
def montecarlo_support(n, ind_funs, center_box, radius_box, num_samples, batch_size, seed):
    """ Computes the volume of the intersection of the box and the support
    of the provided indicator function.

    Input:
    ------------------------------------------------------------
    n           : dimension of ambient space (natural number)

    ind_funs    : list of parallelized functions, returnining True on points within the set.
    
    center_box  : center point of box
    radius_box  : half the side length of the box

    num_samples : the number of points to sample within [-1, 1]^n
    batch_size  : number of samples per batch
    seed        : seed for randomness

    Output:
    ------------------------------------------------------------
    v       : approximately    Vol({x in R^n|  (ind_fun(x) = 1) and (x in box) })
    """

    if num_samples <= 0: raise TypeError("num_samples needs to be positive (provided: {})".format(num_samples))

    # Initialize random number generator with seed
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
    return num_within_support / num_samples * np.power(2 * radius_box, n)

#
# Application to unit ball:
#

def compute_volume_unit_ball_batched(n, p, num_samples=STD_NUM_SAMPLES, batch_size=STD_BATCH_SIZE, seed=STD_SEED):
    """ Computes the volume of the L_p unit ball in n-dim space.
    Uniformly samples points in the [-1, +1]^n hyper cube and returns the 
    ratio of points inside the ball multiplied by the volume of the hypercube (2^n)

    Input:
    ------------------------------------------------------------
    n       : dimension of ambient space (natural number)
    p       : integer defining the norm (even natural number)
    num_samples
            : the number of points to sample within [-1, 1]^n
    seed    : seed for randomness
    
    Output:
    ------------------------------------------------------------
    Vol({x  |  |x|_p <= 1})
    """

    if num_samples <= 0: raise TypeError("num_samples needs to be positive (provided: {})".format(num_samples))

    if not p % 2 == 0: raise TypeError("Expected an even integer p (provided: {})".format(p))

    ind_fun = lambda x_vec : np.sum(np.power(x_vec, p), axis=1) <= 1

    return montecarlo_support(n, ind_funs=[ind_fun], center_box=0, radius_box=1, num_samples=num_samples, batch_size=batch_size, seed=seed)

# Handles to specific cases of batching.
def compute_volume_unit_ball_parallel(n, p, num_samples, seed):
    return compute_volume_unit_ball_batched(n,p, num_samples, batch_size=num_samples, seed=seed)
    
def compute_volume_unit_ball_sequential(n, p, num_samples, seed):
    return compute_volume_unit_ball_batched(n,p, num_samples, batch_size=1, seed=seed)


####### Particular Setups

def intersection_of_k_Lp_n_balls(n,p,centers, center_box, radius_box, num_samples, batch_size, seed):
    # centers should be a list of $n$-vectors
    assert(centers.shape[1] == n)

    # Set up the indicator functions
    ind_funcs = [
        lambda xvec, shift=mu : p_norm_squared(xvec - shift, p) <= 1 for mu in centers
    ]

    # Compute volume
    return montecarlo_support(n, ind_funcs, center_box, radius_box, num_samples, batch_size, seed)




###### Example computations ############

def intersection_L_p_n_balls_shifted_by_basis_vectors(n,p,k, num_samples=STD_NUM_SAMPLES, batch_size=STD_BATCH_SIZE, seed = STD_SEED):
    # Takes L_p balls centered at 0 and the first k basis vectors.
    
    assert(n >= k)

    centers = np.vstack((np.zeros(n,),np.identity(n)[:k])) # The origin and the first k basis vectors in n dims.

    # Set up sampling region: Since one of the balls is centered at 0, that box suffices.    --- (This is just [0,1]^n with the center point at (..,1/2,..))
    center_box = np.zeros((n,))
    radius_box = 1

    vol = intersection_of_k_Lp_n_balls(n, p, centers, center_box, radius_box, num_samples, batch_size, seed)
            

    print("The volume of {} L(p={},n={}) shifted by 0 and the first {} basis vectors is:\n {}".format(
        k, p, n, k, vol
    ))


def previous_examples(num_samples, batch_size, seed):
    n = 2 
    p = 4

    #
    # Computation of the volume of the unit p ball in R^n
    #

    # exact_volume = 4/3 * np.pi
    
    # approximated_volume = compute_volume_unit_ball_batched(n, p, num_samples, batch_size, seed) 
    # #    compute_volume_unit_ball_parallel(n,p,num_samples,seed)
        
    # print("(n,p) = ({},{})".format(n,p))
    # # print("{} <- exact volume".format(exact_volume))
    # print("{} <- approximated volume of unit L_p ball in R^n".format(approximated_volume))
    
    # After shifting the box (should get e.g. half the volume)

    center_box = np.array([0,0])
    ind_fun_unit_ball = lambda x_vec : np.sum(np.power(x_vec, p), axis=1) <= 1
    radius_box = 1
    #volume_shifted_box = montecarlo_support(n, [ind_fun_unit_ball], center_box, radius_box, num_samples, batch_size, seed)

    #print("")
    #print("Vol after shifting box: {}".format(volume_shifted_box))

    # When increasing the box radius, the accuracy should go down.

    radius_box = 2
    #volume_dilated_box = montecarlo_support(n, [ind_fun_unit_ball], 0, radius_box, num_samples, batch_size, seed)
    
    #print("")
    #print("Volume after dilating the box: {}".format(volume_dilated_box))

    #
    # Computation of intersection of balls:  B \cap B + tx
    #
    
    
    t = 1
    x = np.array([1,0])

    ind_fun_shifted_ball = lambda x_vec : p_norm_squared(x_vec - t*x, p) <= 1

    radius_box = 1 # Note the intersection of unit ball and another ball will still be contained in the centered cube of radius 1.
    volume_intersection_two_balls = montecarlo_support(n, [ind_fun_unit_ball, ind_fun_shifted_ball], 0, radius_box, num_samples, batch_size, seed)

    print("")
    print("Volume of B intersected with B + t*x: {}".format(volume_intersection_two_balls))
    print("--> t = {} and x = {}".format(t,x))

    print("Similarly: {}".format(montecarlo_support(n, [lambda xvec : p_norm_squared(xvec, p) <= 1,lambda xvec : p_norm_squared(xvec - x, p) <= 1], 0, radius_box, num_samples, batch_size, seed)))

    print("")
    # print("When taking the box around the center of the shifted ball, it shouldnt change: {}".format(
    #     montecarlo_support(n, 
    #                        [ind_fun_unit_ball, ind_fun_shifted_ball], 
    #                        center_box=t*x,                
    #                        radius_box=radius_box, num_samples=num_samples, batch_size=batch_size, seed=seed)
    # ))


    #
    # Next: Could visualize & plot how the volume of the intersection changes for changing t.
    #


if __name__ == "__main__":
    # intersection_L_p_n_balls_shifted_by_basis_vectors(n=2,p=4,k=1)
    # previous_examples(1000000, 10000, 12345)

    n=4
    p=2


    mus = [np.array([0,0,0,0]), np.array([1,0,0,0])]

    ind_fun_unit_ball = lambda x_vec : p_norm_squared(x_vec - mus[0], p) <= 1

    ind_fun_shifted_ball = lambda x_vec : p_norm_squared(x_vec - mus[1], p) <= 1

    center_box = np.array([0,0,0,0])

    radius_box = 1 # Note the intersection of unit ball and another ball will still be contained in the centered cube of radius 1.
    num_samples = 100*1000*1000
    batch_size = 1000
    seed=1235612

    volume_intersection_two_balls = montecarlo_support(n, [ind_fun_unit_ball, ind_fun_shifted_ball], 0, radius_box, num_samples, batch_size, seed)

    print("The volume approximation: {}".format(volume_intersection_two_balls))

