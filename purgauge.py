"""Analysis of Saumen's pure gauge Wilson flow observables.
"""

import sys
import numpy as np

from parse_wflow import convert_block as cb
from tifr.resample import JKsigma

def inname(B, conf):
    root = "/Users/atlytle/Documents/wflow_out"
    return root + '/bt{0}/flow.{1}'.format(B, conf)

def main(argv):

    # Read in data.
    blocks = []
    for conf in 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500:
        with open(inname(6.17, 5500), 'r') as f:
            blocks.append(f.readlines())

    blocks = np.array(map(cb, blocks))
    
    print blocks.shape

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))