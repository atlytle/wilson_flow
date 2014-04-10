import itertools
import numpy as np

nx,ny,nz,nt = 16, 16, 16, 16
Nsite = nx*ny*nz*nt
nc = 3

adj = lambda x: np.transpose(np.conj(x))  # Adjoint.


# Could just store all indexing info and avoid function calls?
def index(x,y,z,t,mu):
    "Linear index corresponding to link (x,y,z,t,mu)."
    x,y,z,t = x%nx, y%ny, z%nz, t%nt  # Periodic BCs.
    
    # Data stored according to "checkerboard".
    # You might reorder this data earlier in read_su3.py.
    if (x+y+z+t) % 2 == 0:
        return mu*nt*nz*ny*nx + (t*nz*ny*nx + z*ny*nx + y*nx + x)/2
    else:
        return mu*nt*nz*ny*nx + Nsite/2 + (t*nz*ny*nx + z*ny*nx + y*nx + x)/2
    
def ix(xvec, mu):
    "Linear index corresponding to link (x,y,z,t,mu)."
    x,y,z,t = xvec[0]%nx, xvec[1]%ny, xvec[2]%nz, xvec[3]%nt  # Periodic BCs.
    return index(x,y,z,t,mu)
    
def fx(i):
    "Link coordinates (x,y,z,t), mu corresponding to linear index i."
    mu = i/Nsite
    i = i-mu*Nsite
    eo, i = i/(Nsite/2), i%(Nsite/2)  # Which checkerboard.
    t,i = i/(nz*ny*nx/2), i%(nz*ny*nx/2)
    z, i = i/(ny*nx/2), i%(ny*nx/2)
    y, i = i/(nx/2), i%(nx/2)  # try divmod.
    x = 2*i 
    if not eo and (t+z+y) % 2 == 1: x +=1
    if eo and (t+z+y) % 2 == 0: x +=1

    return (x,y,z,t), mu
    
def test_index():
    "Consistency of ix() and fx()."
    for i in range(4*Nsite):
        assert i == ix(*fx(i))

def plaq(config):
    "Calculate plaquette values on a configuration."
    # Slow way of doing things, doesn't use numpy arrays well.
    # Instead you should allocate indice arrays, reorder and multiply.
    plaq = np.zeros((3,3), dtype=np.complex128)
    for mu in range(4):
        for nu in range(mu):
            print mu, nu
            for x, y, z, t in itertools.product(
                       range(nx), range(ny), range(nz), range(nt)):
                xvec = np.array([x,y,z,t])

                mh = muhat(mu)
                nh = muhat(nu)
                
                m1 = config[ix(xvec, mu)]
                m2 = config[ix(xvec+mh, nu)]
                
                m3 = config[ix(xvec, nu)]
                m4 = config[ix(xvec + nh, mu)]
                
                tmp = np.dot(np.dot(m1,m2), adj(np.dot(m3,m4)))
                plaq += tmp
            #print (plaq/(6*Nsite)).real
    return np.trace(plaq/(6*Nsite)).real
    
def plaqt(config):
    "Calculate temporal plaquette values on a configuration."
    # Slow way of doing things, doesn't use numpy arrays well.
    # Instead you should allocate indice arrays, reorder and multiply.
    plaq = np.zeros((3,3), dtype=np.complex128)
    for mu in 3,:
        for nu in 0,1,2:
            print mu, nu
            for x, y, z, t in itertools.product(
                       range(nx), range(ny), range(nz), range(nt)):
                xvec = np.array([x,y,z,t])

                mh = muhat(mu)
                nh = muhat(nu)
                
                m1 = config[ix(xvec, mu)]
                m2 = config[ix(xvec+mh, nu)]
                
                m3 = config[ix(xvec, nu)]
                m4 = config[ix(xvec + nh, mu)]
                
                tmp = np.dot(np.dot(m1,m2), adj(np.dot(m3,m4)))
                plaq += tmp
            #print (plaq/(3*Nsite)).real
    return np.trace(plaq/(3*Nsite)).real
    
def plaqs(config):
    "Calculate spatial plaquette values on a configuration."
    # Slow way of doing things, doesn't use numpy arrays well.
    # Instead you should allocate indice arrays, reorder and multiply.
    plaq = np.zeros((3,3), dtype=np.complex128)
    for mu, nu in (1,2), (1,0), (2,0):
        print mu, nu
        for x, y, z, t in itertools.product(
                   range(nx), range(ny), range(nz), range(nt)):
            xvec = np.array([x,y,z,t])

            mh = muhat(mu)
            nh = muhat(nu)
            
            m1 = config[ix(xvec, mu)]
            m2 = config[ix(xvec+mh, nu)]
            
            m3 = config[ix(xvec, nu)]
            m4 = config[ix(xvec + nh, mu)]
            
            tmp = np.dot(np.dot(m1,m2), adj(np.dot(m3,m4)))
            plaq += tmp
        #print (plaq/(3*Nsite)).real
    return np.trace(plaq/(3*Nsite)).real
    
def F(config, xvec, mu, nu):
    "Lattice field strength tensor F_{mu nu}."
    F = np.zeros((3,3), dtype=complex)
    mh = muhat(mu)
    nh = muhat(nu)
    
    # Clover leaves.
    m1 = config[ix(xvec, mu)]
    m2 = config[ix(xvec+mh, nu)]
    m3 = adj(config[ix(xvec+nh, mu)])
    m4 = adj(config[ix(xvec, nu)])
    F += reduce(np.dot, [m1, m2, m3, m4])
    
    m1 = config[ix(xvec, nu)]
    m2 = adj(config[ix(xvec-mh+nh, mu)])
    m3 = adj(config[ix(xvec-mh, nu)])
    m4 = config[ix(xvec-mh, mu)]
    F += reduce(np.dot, [m1, m2, m3, m4])
    
    m1 = adj(config[ix(xvec-mh, mu)])
    m2 = adj(config[ix(xvec-mh-nh, nu)])
    m3 = config[ix(xvec-mh-nh, mu)]
    m4 = config[ix(xvec-nh, nu)]
    F += reduce(np.dot, [m1, m2, m3, m4])
    
    m1 = adj(config[ix(xvec-nh, nu)])
    m2 = config[ix(xvec-nh, mu)]
    m3 = config[ix(xvec+mh-nh, nu)]
    m4 = adj(config[ix(xvec, mu)])
    F += reduce(np.dot, [m1, m2, m3, m4])
    
    F =(F - adj(F))/8
    return F
    
def Fsq(config):
    "Trace of F^2."
    Fsq = np.zeros((3,3), dtype=complex)
    for x, y, z, t in itertools.product(
                       range(nx), range(ny), range(nz), range(nt)):
        xvec = np.array([x,y,z,t])
        for mu in range(4):
            for nu in range(mu+1,4):
                tmp = F(config, xvec, mu, nu)
                Fsq += np.dot(tmp, tmp)
                
    return np.trace(Fsq)/Nsite/6
    
def Z(config):
    '''Force matrix used to update configurations in Wilson Flow.
    
    Z(x, mu) = \partial_{x,mu} \sum_p plaq.'''
    
    # Initialize result.
    result = np.zeros(config.shape, dtype=complex)
    
    for i in range(config.shape[0]):
        m1 = config[i]  # U_mu(x).
        xvec, mu = fx(i)
        mh = muhat(mu)
        d = [0,1,2,3]
        d.remove(mu)  # nu != mu
        Zxmu = np.zeros((3,3), dtype=complex)
        for nu in d:
            nh = muhat(nu)
            # Staple up.
            m2 = config[ix(xvec+mh, nu)]
            m3 = adj(config[ix(xvec+nh, mu)])
            m4 = adj(config[ix(xvec, nu)])
            su = reduce(np.dot, [m1,m2,m3,m4])
            # Staple down.
            m5 = adj(config[ix(xvec-nh+mh, nu)])
            m6 = adj(config[ix(xvec-nh, mu)])
            m7 = config[ix(xvec-nh, nu)]
            sd = reduce(np.dot, [m1,m5,m6,m7])
            
            Zxmu -= Hproject(su + sd)
            
        result[i] = Zxmu
    return result
    
def Hproject(M):
    "Convert M to a traceless anti-Hermitian matrix."
    tmp = M-adj(M)
    return tmp/2. - np.trace(tmp)/6.
            
    
    
def test_su3(config):
    "Check all matrices are SU(3)."
    res = np.zeros((3,3), dtype=complex)
    for x, y, z, t in itertools.product(
                       range(nx), range(ny), range(nz), range(nt)):
        xvec = np.array([x,y,z,t])
        for mu in range(4):
            tmp = config[ix(xvec, mu)]
            res += np.dot(tmp, adj(tmp))
    return res/(4*Nsite)
    
def muhat(mu):
    if mu == 0:
        return np.array([1,0,0,0])
    if mu == 1:
        return np.array([0,1,0,0])
    if mu == 2:
        return np.array([0,0,1,0])
    if mu == 3:
        return np.array([0,0,0,1])
    
#def plaq2(config, xvec, mu, nu):
#    #ix = index
#    mh = muhat(mu)
#    nh = muhat(nu)
#    
#    # mult these
#    ix(xvec, mu)
#    ix(xvec+mh, nu)
#    
#    # conjugate mult these
#    ix(xvec, nu)
#    ix(xvec+nh, mu)
    
def cmult(c1, c2):
    '''Multiply color matrices of configurations.'''
    pass
    
