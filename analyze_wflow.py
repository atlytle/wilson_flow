"""Analysis of Wilson flow data generated on the cray at TIFR.

"""

import sys
sys.path.append('../tifr')
import numpy as np
import pylab as p

from parse_wflow import file2numpy
from resample import JK_block, JKsigma

inroot = "/Users/atlytle/Documents/wflow_out/"

def inname(L, B, m):
    return inroot+"output_{0}.{0}.{0}.{0}_{1}_{2}.txt".format(L, B, m)

def dfb(ns):
    """Symmetric difference of a series.
    Forward and backward difference used for endpoints.
    """
    result = np.zeros((len(ns)))
    # Endpoints.
    result[0] = ns[1]-ns[0]
    result[-1] = ns[-1] - ns[-2]
    # Interior.
    result[1:-1] = (ns[2:] - ns[:-2])/2.
    
    return result

def plot_items(items):
    """Generic plot routine.
    """
    p.xlabel('step')
    x = range(len(items[0]))
    for y in items:
        p.errorbar(x,y)
    p.show()

def main(argv):
    dat1 = file2numpy(inname(16,5.2875,0.025))
    dat2 = file2numpy(inname(16, 5.4, 0.025))
    nconf, nstep, nmeas = dat1.shape
    #print dat2.shape
    #print dat1[0,:,1]
    #print dat1[0,:,2]
    plot_items([dat1[x,:,2] for x in range(1,nconf)])#, dat1[0,:,2], dat2[0,:,1], dat2[0,:,2]])

    #JKblock
    t2E = JK_block(dat1[1:,:,2])
    print t2E[0]
    print JKsigma(t2E)
    t = np.linspace(0,1.0,nstep)
    plot_items([t2E[0], 100*t*dfb(t2E[0])])
    
    nconf = dat2.shape[0]
    plot_items([dat2[x,:,2] for x in range(1,nconf)])
    t2E_2 = JK_block(dat2[:,:,2])
    plot_items([t2E_2[0], 100*t*dfb(t2E_2[0])])

    #plot
    #extract t0
    #extract w0
    #dfb(np.array([1,3,3,4,7,6]))
    return 0
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
