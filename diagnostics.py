"""Diagnostics of wflow integration.

Dependence of observables on step size, integration scheme, etc.
"""

import sys
import numpy as np
import pylab as p

from parse_wflow import strip_blocks

inroot = "/Users/atlytle/Documents/wflow_out/"
inname = inroot+"output_dt.txt"

def convert(ts,ys):
    """Index a time-series by the time.
    Used to manipulate time-series of different granularities."""
    d = {}
    for t,y in zip(ts,ys):
        d[str(t)] = y
    return d

def difference(d, dref):
    """Take difference of two time series created using 'convert()'.
    Reference dref should have finer granularity to 
    ensure the difference is possible.
    """
    times = d.keys()
    times.sort()
    delta = [d[t]-dref[t] for t in times]
    return delta

def main(argv):
    with open(inname, 'r') as f:
        dat = f.readlines()

    blocks = strip_blocks(dat, 8)
    blocks = map(np.array, blocks)

    # Index the series by time.
    ds = []
    for b in blocks:
        t = b[:,0]
        y = b[:,2]
        #p.plot(t,y)
        ds.append(convert(t,y))
    #p.show()
    
    # Take differences of time series.
    ref = ds[-1]  # Most accurate.
    deltas = []
    for d in ds:
        tmp = d.keys()
        tmp.sort()  # Ascending time order.
        delta = []
        for t in tmp:
            delta.append(d[t]-ref[t])
        deltas.append(delta)

    # Plot differences. 
    p.xlabel('t')
    p.ylabel('$\Delta$')
    color = ['b', 'r', 'm', 'k', 'g']
    dt = [.1,.05,.025,.01,.005]
    for i in range(len(deltas)-1):
        t = blocks[i][:,0]
        if i % 2 == 0:
            f = '--'+color[i/2]
            l = 'dt={0}'.format(dt[i/2])
        else:
            f = '-'+color[i/2]
            l = None
        p.errorbar(t, deltas[i], fmt=f, label=l)
    p.legend(loc=0)
    p.ylim([-.001,.001])
    p.show()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))