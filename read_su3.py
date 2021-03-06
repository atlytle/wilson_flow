"""Read binary configuration files and reconstruct SU(3) matrices.
"""

import sys
import itertools
import numpy as np
from os.path import getsize

nx, ny, nz, nt = 16, 16, 16, 16
Nsite = nx*ny*nz*nt
nc = 3

def readhead(filename):
    "Read header of fortran binaries."
    f = open(filename, "rb")
    bsize = np.fromfile(f, dtype='>i', count=1)  # Bytes of actual data.
    beta = np.fromfile(f, dtype='>d', count=1)  # Beta.
    nf = np.fromfile(f, dtype='>i', count=1)  # Nflav.
    # Haven't parsed other stuff.
    return beta, nf

def readfort(filename):
    '''Read gauge configuration from Fortran90 binary.
    
    There are six complex entries per su3 matrix as only two rows are stored.
    The first and last entries are 4byte integers specifying the size in bytes
    of the rest of the file.  This is a Fortran convention.
    '''
    csize = 16  # Complex number size in bytes.
    hsize = 96  # Header size in bytes.
    assert getsize(filename) == 4 + hsize + 4*Nsite*6*csize + 4
    f = open(filename, "rb")
    bsize = np.fromfile(f, dtype='>i', count=1)  # Bytes of actual data.
    assert bsize == hsize + 4*Nsite*6*csize
    f.seek(hsize, 1)  # Ignore header.
    tmp = np.fromfile(f, dtype='>c16', count=4*Nsite*6)  # Two rows stored.
    # Last entry gives bytesize again.
    assert bsize == np.fromfile(f, dtype='>i', count=1)

    f.close()
    return tmp

def fort_id(j,i,x):
    '''Index structure of binary configurations. See lat.F.'''
    return j*2*(Nsite*4) + i*(Nsite*4) + x

# Put into a more familiar matrix form.
new_ids = [fort_id(j,i,x) for x,i,j in itertools.product(
                                       range(4*Nsite), range(2), range(3))]

def matrec(r1, r2):
    "Construct SU(3) matrix from first two rows."
    r3 = np.conj(np.cross(r1, r2))
    return np.array([r1,r2,r3])
    
def reconstruct(config):
    "Reorder data and reconstruct configurations."
    dat = config[new_ids]
    dat = dat.reshape(4*Nsite,2,3)
    dat = np.array(map(lambda x:matrec(x[0],x[1]), dat))
    assert dat.shape == (4*Nsite,3,3)
    
    return dat
    
def main(files):

    tmp = readfort(files[0])
    v = tmp[0], tmp[2*4*Nsite], tmp[4*4*Nsite]
    print v
    print np.inner(v, np.conj(v))
    dat = tmp[new_ids]
    print dat[0], dat[1], dat[2]
    v2 = dat[3], dat[4], dat[5]  # Second row of first matrx.
    print np.inner(v2, np.conj(v2))
    print np.inner(v, np.conj(v2))

    x = matrec(v,v2)
    assert x.shape == (3,3)
    print np.dot(x, np.conj(np.transpose(x)))

    # Reconstruct matrices.
    dat = dat.reshape(4*Nsite,2,3)
    dat = np.array(map(lambda x:matrec(x[0],x[1]), dat))
    assert dat.shape == (4*Nsite,3,3)
    print np.dot(dat[42], np.conj(np.transpose(dat[42])))
    print ''
    print dat[0]
    print dat[1]

    return 0
    

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

