###################################################################################
# R_STATS
###################################################################################
"""





"""


import numpy as np

#######################################################
# 2D correlation
#######################################################



def correlation(var1,var2):
    '''
    Normalized cross correlation between 2 same size arrays
    
    Uses only values where both arrays are finite
    
    '''
    var1_t = var1[np.logical_and(np.isfinite(var1),np.isfinite(var2))]
    var2_t = var2[np.logical_and(np.isfinite(var1),np.isfinite(var2))]

    var1_n = (var1_t - np.mean(var1_t))/np.std(var1_t)
    var2_n = (var2_t - np.mean(var2_t))/np.std(var2_t)
    
    return np.correlate(var1_n,var2_n)/len(var1_n)


    



