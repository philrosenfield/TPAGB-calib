import numpy as np

def count_uncert_ratio(numerator, denominator):
    '''
    combining poisson error for taking a ratio
    '''
    n = float(numerator)
    d = float(denominator)
    return (n / d) * (1./np.sqrt(n) + 1./np.sqrt(d))

