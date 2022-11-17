import numpy as np
from smallexpimv import smallexpimv

class expimv_tridiag:
  def __init__(self,m):
    self.ldz = m
    self.z = np.zeros([m,m], dtype=np.double, order='F')
    self.w = np.zeros([m], dtype=np.double, order='F')
    self.isuppz = np.zeros([2*m], dtype=np.intc, order='F')
    self.work = np.zeros([20*m], dtype=np.double, order='F')
    self.iwork = np.zeros([10*m], dtype=np.intc, order='F')
    self.cwork = np.zeros([m], dtype=np.cdouble, order='F')
    self.pmone = -1

  def precompute(self,td,tsd,n):
    self.n = n
    info = smallexpimv.expimv_tridiag_precompute(td, tsd, self.w, self.z, n, self.isuppz, self.work, self.iwork)
    return info

  def eval(self, dt, beta):
    x = np.zeros([self.n], dtype=np.cdouble)
    smallexpimv.expimv_tridiag_apply(dt, self.pmone, beta, self.w, self.z, x, self.cwork)
    return x

class expimv:
  def __init__(self, kmax, ldh):
    self.Y = np.zeros([ldh*(kmax+1)], dtype=np.cdouble, order='F')
    self.pmone = -1
    self.kmax = kmax
    self.tol = 10**-15
    self.maxrestart = 200

  def apply(self,dt,H,x,m,nest=1.0):
    self.Y[:m] = x
    (dtdone, info ) = smallexpimv.expimv( dt, H, m, self.Y, self.pmone, self.kmax, self.tol, nest, self.maxrestart)
    y=np.zeros([m], dtype=np.cdouble)
    y[:] = self.Y[:m]
    return y






