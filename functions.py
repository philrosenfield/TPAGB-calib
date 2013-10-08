import numpy as np

def gaussian(x, p):
    '''
    gaussian(arr,p): p[0] = norm, p[1] = mean, p[2]=sigma
    '''
    return p[0] * np.exp( -1 * (x - p[1])**2 / (2 * p[2]**2))

def double_gaussian(x, p):
    '''
    gaussian(arr,p): p[0] = norm1, p[1] = mean1, p[2]=sigma1
                     p[3] = norm2, p[4] = mean2, p[5]=sigma2
    '''
    return gaussian(x, p[:3]) + gaussian(x, p[3:])
   
def mp_double_gauss(p, fjac=None, x=None, y=None, err=None):
    '''
    double gaussian for mpfit
    '''
    model = double_gaussian(x, p)
    status = 0
    return [status, (y - model) / err]

def mp_gaussian(p, fjac=None, x=None, y=None, err=None):
    '''
    single gaussian for mpfit
    '''
    model = gaussian(x, p)
    status = 0
    return [status, (y - model) / err]
