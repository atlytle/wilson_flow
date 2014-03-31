import sys
import numpy as np
from os.path import getsize

nx, ny, nz, nt = 16, 16, 16, 16
Nsite = nx*ny*nz*nt
nc = 3
ns = 4

def readfort(file):
    '''Read gauge configuration from Fortran90 binary.
    
    There are six complex entries per su3 matrix as only two rows are stored.
    The first and last entries are 4byte integers specifying the size in bytes
    of the rest of the file.  This is a Fortran convention.
    '''
    csize = 16  # Complex number size in bytes.
    hsize = 96  # Header size in bytes.
    assert getsize(file) == 4 + hsize + 4*Nsite*6*csize + 4
    f = open(file, "rb")
    bsize = np.fromfile(f, dtype='>i', count=1)  # Bytes of actual data.
    assert bsize == hsize + 4*Nsite*6*csize
    f.seek(hsize, 1)  # Ignore header.
    tmp = np.fromfile(f, dtype='>c16', count=4*Nsite*6)  # Two rows stored.
    # Last entry gives bytesize again.
    assert bsize == np.fromfile(f, dtype='>i', count=1)

    f.close()
    return tmp


def main(files):

    tmp = readfort(files[0])
    v = tmp[0], tmp[2*4*Nsite], tmp[4*4*Nsite]
    print np.inner(v, np.conj(v))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
