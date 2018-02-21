## utility.py ##
'''
Code containing auxilliary functions

'''
import numpy as np

### KL Divergence ###
def KL_div(N_0,N_1):
    """
        Calculates the KL divergence D(N_0||N_1) between two normal distributions
    ---
    Inputs:
        N_0: list of [mean,variance] (Integration wrt)
        N_1: list of [mean,variance]

    ---
    Output: 
        D = distance measure

    """
    m0 = N_0[0]; v0 = N_0[1]
    m1 = N_1[0]; v1 = N_1[1]
    
    D = 0.5*(v0/v1 + (m1-m0)*(1/v1)*(m1-m0)-1+np.log(v1/v0))
    return D