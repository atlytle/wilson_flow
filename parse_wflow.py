"""Read raw data from Sourendu's wflow (fortran) output.

Convert to numpy arrays.
"""

import sys
import numpy as np

offset = 8
nstep = 100
def convert_block(block):
    "Convert a block of wflow data (one configuration) to floats."
    block = [x.strip().split('  ') for x in block]
    block = [map(float, x) for x in block]
    return block

def file2numpy(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Extract blocks of data from the output file.
    i = 0
    blocks = []
    for x in lines:
        if x.strip() == 'Wilson flow: program version    1.000':
            blocks.append(lines[i+offset:i+offset+nstep])
        i+=1
    del i
    del lines

    # Convert to numerical data.
    blocks = map(convert_block, blocks)
    blocks = np.array(blocks)
    return blocks

def main(argv):
    blocks = file2numpy(argv[0])
    
    print blocks[:,:,1].shape
    print blocks[0,:,1]

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
