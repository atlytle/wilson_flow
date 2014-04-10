import sys
import numpy as np

from read_su3 import readfort, reconstruct
from su3 import plaq, plaqt, plaqs, test_su3, F, Fsq

nx,ny,nz,nt = 16, 16, 16, 16
Nsite = nx*ny*nz*nt
nc = 3

unit = np.array([np.identity(3, dtype=complex) for i in range(4*Nsite)])

def main(files):
    dat = readfort(files[0])
    dat = reconstruct(dat)
    print dat.shape
    #print test_su3(dat)
    print 'plaq', plaq(dat)
    print 'plaqt', plaqt(dat)
    print 'plaqs', plaqs(dat)
    #tmp = F(dat, (0,0,0,0), 0, 1)
    #print tmp
    #print np.dot(tmp, tmp)
    #print Fsq(dat)
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
