import numpy as np
import scipy.stats as sts

def quantile_normalization(multidim_array):
    '''perform quantile normalization on columns of 
    a multidimensional array.
    to use:
    > import quantile_normalization as qn
    > normalized_array = qn.quantile_normalization(unnormalized_array)'''
    m = []
    for i in multidim_array.T:
        m.append(sts.rankdata(i)-1)
    
    b = []
    for i in multidim_array.T:
        f=  i.copy()
        f.sort()
        b.append(f)
    b = np.array(b).T
    avs = np.average(b, axis=1)
    #print avs
    c = []
    for i in xrange(len(m)):
        q = m[i]
        q = q.astype(int)
        c.append(avs[q])
    return np.array(c).T


