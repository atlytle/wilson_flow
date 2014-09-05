from random import random
import numpy as np
from numpy import sin, cos, arccos, sqrt, exp

adj = lambda x: np.transpose(np.conjugate(x))
norm = lambda x: np.trace(np.dot(x,adj(x)))

def expsu3(M):
    "Exponentiate anti-Hermitian traceless matrix ala hep-lat/0311018."
    assert M.shape == (3,3)
    #print np.trace(M)
    #print adj(M)
    
    c0 = np.trace(reduce(np.dot,[M,M,M]))/3.  # Tr(M^3)/3.
    cflag = False
    if c0 < 0 :
        c0 = -c0
        cflag = True  # Use symmetry properties of f_i.
    c1 = np.trace(np.dot(M,M))/2.  # Tr(M^2)/2.
    
    # Derived quantities.
    c0max = 2*(c1/3)**(1.5)
    th = arccos(c0/c0max)
    u = sqrt(c1/3)*cos(th/3)
    uu = u*u
    w = sqrt(c1)*sin(th/3)
    ww = w*w
    d = 9*uu - ww
    e2iu = exp(2j*u)
    emiu = exp(-1j*u)
    cw = cos(w)
    
    # Numerical evaluation of sin(w)/w.
    if abs(w) < 0.05:
        xi0 = 1-(ww/6)*(1-(ww/20)*(1-ww/42))
    else:
        xi0 = sin(w)/w
        
    h0 = (uu-ww)*e2iu + emiu*(8*uu*cw + 2j*u*(3*uu+ww)*xi0)
    h1 = 2*u*e2iu - emiu*(2*u*cw-1j*(3*uu-ww)*xi0)
    h2 = e2iu - emiu*(cw+3j*u*xi0)
    
    f0, f1, f2 = h0/d, h1/d, h2/d
    if cflag:
        f0 = np.conjugate(f0)
        f1 = -np.conjugate(f1)
        f2 = np.conjugate(f2)
    
    return f0*np.identity(3) + f1*M + f2*np.dot(M,M)
    
def test_xi0(w):
    ww = w*w
    return 1-(ww/6)*(1-(ww/20)*(1-ww/42)), sin(w)/w
    
def randherm():
    '''A "random" traceless hermitian matrix (3x3).'''
    c = map(complex, 
            [random() for i in range(9)], [random() for i in range(9)])
    r = np.array(c).reshape((3,3))
    tmp = r + adj(r)
    return tmp - np.trace(tmp)/3
    
def test():
    tol = np.power(10.,-13)  # 10^-15 typical, 10^-14 rarely exceeded.
    n = 1000000
    print "Running", n, "iterations with tol=", tol
    for i in range(n):
        h = randherm()
        u = expsu3(randherm())
        tmp = norm(u)-3.
        if abs(tmp) > tol:
            print "Iteration", i, "failed:"
            print "tol exceeded:", abs(tmp)
            print "Hermitian:", h
            print "Unitary?:", u
            return 1
    print "tol never exceeded."
    return 0
    
if __name__ == "__main__":
    test()
    
    
