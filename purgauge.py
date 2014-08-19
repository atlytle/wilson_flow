"""Analysis of Saumen's pure gauge Wilson flow observables.
"""

import sys
import numpy as np
from scipy.interpolate import interp1d

from parse_wflow import convert_block as cb
from tifr.resample import JK_block, JKsigma
from analyze_wflow import plot_items, dfb2, broot

def inname(B, conf):
    root = "/Users/atlytle/Documents/wflow_out"
    return root + '/bt{0}/flow.{1}'.format(B, conf)

def main(argv):

    # Read in data.
    blocks = []
    for conf in 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500:
        with open(inname(6.17, conf), 'r') as f:
            blocks.append(f.readlines())
    d617 = np.array(map(cb, blocks))
    
    blocks = []
    for conf in 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000:
        with open(inname(5.96, conf), 'r') as f:
            blocks.append(f.readlines())
    d596 = np.array(map(cb, blocks))

    dt = 0.01

    # Analyze Beta = 5.96
    nconf, nstep, nmeas = d596.shape
    t = np.linspace(0,dt*nstep,nstep)
    #plot_items([d596[x,:,4] for x in range(nconf)])
    t2Ep = JK_block(d596[:,:,2])  # Ep = E from plaquette.
    t2Ec = JK_block(d596[:,:,4])  # Ec = E from clover term.
    dt2Ep = t*dfb2(t2Ep[0])/dt  # Logarithmic derivative.
    dt2Ec = t*dfb2(t2Ec[0])/dt
    plot_items([t2Ep[0], t2Ec[0], dt2Ep, dt2Ec])
    for s in t2Ep[0], t2Ec[0], dt2Ep, dt2Ec:
        f = interp1d(t, s-0.3, kind='cubic')
        print np.sqrt(broot(f, 0.0, nstep*dt))

    # Analyze Beta = 6.17
    dt = 0.01
    nconf, nstep, nmeas = d617.shape
    t = np.linspace(0,dt*nstep,nstep)
    #plot_items([d617[x,:,4] for x in range(nconf)])
    t2Ep = JK_block(d617[:,:,2])  
    t2Ec = JK_block(d617[:,:,4])  
    dt2Ep = t*dfb2(t2Ep[0])/dt  
    dt2Ec = t*dfb2(t2Ec[0])/dt
    plot_items([t2Ep[0], t2Ec[0], dt2Ep, dt2Ec])
    for s in t2Ep[0], t2Ec[0], dt2Ep, dt2Ec:
        f = interp1d(t, s-0.3, kind='cubic')
        print np.sqrt(broot(f, 0.0, nstep*dt))

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))