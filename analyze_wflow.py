"""Analysis of Wilson flow data generated on the cray at TIFR.

"""

import sys
sys.path.append('../tifr')
import numpy as np
import pylab as p
from scipy.interpolate import interp1d

from parse_wflow import file2numpy
from resample import JK_block, JKsigma

np.set_printoptions(precision=4)

def inname(L, B, m, id=''):
    root = "/Users/atlytle/Documents/wflow_out/"
    file = "output_{0}.{0}.{0}.{0}_{1}_{2}{3}.txt".format(L, B, m, id)
    return root + file

def dfb(ns):
    """Symmetric difference of a series.
    Forward and backward difference used for endpoints.
    ... notation makes it work with 1d time series, or a stack of time series.
    """
    result = np.zeros(ns.shape)
    # Endpoints.
    result[...,0] = ns[...,1]-ns[...,0]
    result[...,-1] = ns[...,-1] - ns[...,-2]
    # Interior.
    result[...,1:-1] = (ns[...,2:] - ns[...,:-2])/2.
    
    return result

def dfb2(f):
    """More accurate estimate of the derivative -- O(dt^4) in the interior.
    """
    df = np.zeros(f.shape)
    # Endpoints.
    df[...,0] = f[...,1]-f[...,0]  
    df[...,-1] = f[...,-1] - f[...,-2]  # O(dt).  
    df[...,1] = (f[...,2]-f[...,0])/2.
    df[...,-2] = (f[...,-1] - f[...,-3])/2.  # O(dt^2).
    # Interior.
    df[...,2:-2] = -f[...,4:]/12. +\
                  (2/3.)*(f[...,3:-1] - f[...,1:-3]) + f[...,:-4]/12.

    return df


def broot(f, xi, xf, tol=0.0001):
    '''Find a root of f(x) in [xi,xf] via bisection.
    '''
    s = np.sign(f(xi))
    # Make sure f changes sign on the interval.
    assert s != np.sign(f(xf))
    # Check if the endpoints are roots.
    if abs(f(xi)) < tol:
        return xi
    if abs(f(xf)) < tol:
        return xf
    while True:
        mp = (xi+xf)/2
        if abs(f(mp)) < tol:
            return mp
        if np.sign(f(mp)) == s:
            xi = mp
        else:
            xf = mp

def plot_items(items):
    """Generic plot routine.
    """
    p.xlabel('step')
    x = range(len(items[0]))
    for y in items:
        p.errorbar(x,y)
    p.show()

def main(argv):
    dt = 0.01
    dat1 = file2numpy(inname(16,5.2875,0.025), offset=8) # erratic 1.
    dat2 = file2numpy(inname(16, 5.4, 0.025), offset=8)
    dat3 = file2numpy(inname(24, 5.6, 0.025), offset=8)
    dat4 = file2numpy(inname(16, 5.4, 0.025, '_s'), offset=2)
    dat5 = file2numpy(inname(16, 5.4, 0.025, 'b'), offset=10)
    dat6 = file2numpy(inname(16, 5.4, 0.015, 'b'), offset=9)
    dat7 = file2numpy(inname(16, 5.4, 0.05, 'b'), offset=9)  # erratic 1-4.
    dat8 = file2numpy(inname(16, 5.5, 0.025, 'b'), offset=9)
    dat9 = file2numpy(inname(16, 5.5, 0.05, 'b'), offset=9)
    dat10 = file2numpy(inname(16, 5.6, 0.025, 'b'), offset=9)
    dat11 = file2numpy(inname(16,5.2875,0.025, 'b'), offset=9) # erratic 1?


    for d in dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10, dat11:
        print d.shape

    for d in dat10, dat11:
        nconf, nstep, nmeas = d.shape
        plot_items([d[x,:,3] for x in range(nconf)])

    print 1.5*dat1[1,:,2]
    print dat11[1,:,1]

    #plot_items([dat7[x,:,1] for x in range(4)])
    #plot_items([dat7[x,:,1] for x in range(4,16)])


    # JKblock.
    # nconf, nstep, nmeas = dat2.shape
    # t2E = 1.5*JK_block(dat2[:16,:,2])
    # t = np.linspace(dt, dt*nstep, nstep)
    # p.xlim([0,2])
    # p.errorbar(t, t2E[0])
    # p.errorbar(t, t*dfb2(t2E[0])/dt, fmt='k')
    # p.errorbar(t, 1.5*dat2[0,:,2], fmt='--')

    # nconf, nstep, nmeas = dat4.shape
    # t2E = JK_block(dat4[:,:,2])
    # t = np.linspace(0,dt*nstep,nstep)
    # p.errorbar(t, t2E[0])
    # p.errorbar(t, t*dfb2(t2E[0])/dt)
    # p.errorbar(t, dat4[0,:,2], fmt='--')


    
    # p.show()
    
    #t2E = JK_block(dat3[1:,:,2])
    #t = np.linspace(0,dt*nstep,nstep)
    #plot_items([t2E[0], 100*t*dfb(t2E[0])])
    #plot_items([x for x in t2E])

    #error = JKsigma(t2E)
    #error_dt = t*JKsigma(dfb2(t2E))/dt

    # Plot.

    # p.errorbar(range(nstep),t2E[0]+error)
    # p.errorbar(range(nstep),t2E[0]-error)
    # p.errorbar(range(nstep),t*dfb2(t2E[0])/dt+error_dt)
    # p.errorbar(range(nstep),t*dfb2(t2E[0])/dt-error_dt)
    # p.show()

    # Extract t0, w0.

    #f=interp1d(t,t2E[0]-0.3, kind='cubic')
    #g=interp1d(t, 100*t*dfb(t2E[0])-0.3, kind='cubic')
    #print broot(f,0.01,0.99)
    #print broot(g,0.01,0.99)
    return 0
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
