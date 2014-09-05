import sys
import numpy as np

from read_su3 import readfort, reconstruct
from su3 import plaq, plaqt, plaqs, test_su3, F, Fsq, Z, cmult
from expsu3 import expsu3

nx,ny,nz,nt = 16, 16, 16, 16
Nsite = nx*ny*nz*nt
nc = 3

unit = np.array([np.identity(3, dtype=complex) for i in range(4*Nsite)])

def main(files):
    print 'Setting up configuration..'
    dat = readfort(files[0])
    dat = reconstruct(dat)
    #print dat.shape
    #print test_su3(dat)
    
    print 'Calculating plaq..'
    print 'plaqt:', plaqt(dat)
    
    #print 'plaqt:', plaqt(dat)
    #print 'plaqs:', plaqs(dat)
    
    #tmp1 = F(dat, (0,0,0,0), 0, 1)
    #print tmp1
    #print ''
    
    #print Fsq(dat)

    dt = .01
    # print 'Calculating force term..'
    # z = Z(dat)
    # print 'Exponentiating force term..'
    # expZ = map(expsu3, -1j*dt*z)
    # #print test_su3(expZ)
    # print 'Updating gauge field..'
    # dat = cmult(expZ,dat)

    # print 'Calculating plaq..'
    # print 'plaq:', plaq(dat)

    print 'Calculating force term..'
    z = Z(dat)
    #print "z[0]"
    #print z[0]
    #print ""
    print 'Exponentiating force term..'
    expZ = map(expsu3, -1j*dt*z)
    #print expZ[0]
    # print test_su3(expZ)
    print 'Updating gauge field..'
    dat = cmult(expZ,dat)

    print 'Calculating plaq..'
    print 'plaqt:', plaqt(dat)
    #print expsu3(-1j*z[0])


    print 'done.'
    return 0
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
