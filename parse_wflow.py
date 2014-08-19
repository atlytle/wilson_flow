"""Read raw data from Sourendu's wflow (fortran) output.

Convert to numpy arrays.
"""

import sys
import numpy as np

def convert_block(block):
    "Convert a block of wflow data (one configuration) to floats."
    block = [x.strip().split('  ') for x in block]
    block = [map(float, x) for x in block]
    return block

def strip_blocks(lines, offset):
    """Pull out numerical blocks of data from wflow output.
    
    offset gives the linelength of the header parts.
    """
    i = 0
    blocks = []
    for x in lines:
        if x.strip() == 'Wilson flow: program version    1.000':  # Begin block.
            block = []
            for y in lines[i+offset:]:
                if y.strip() == '':  # End block.
                    break
                block.append(y.strip())
            blocks.append(block)
        i+=1
    # Convert to numerical data.
    blocks = map(convert_block, blocks)
    
    return blocks

def file2numpy(filename):
    """Convert wflow output files to numpy data.
    """
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    offset = 8
    blocks = strip_blocks(lines, offset)
    
    return np.array(blocks)

def file2gnuplot(filename):
    "Convert wflow output files to gnuplot files."
    pass

def main(argv):
    blocks = file2numpy(argv[0])
    
    print blocks[:,:,1].shape  # (nconfig, nstep, nmeas)
    print blocks[0,:,2]

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
