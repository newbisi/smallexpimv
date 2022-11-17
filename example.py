# python3 code
import smallexpimv_pyclasses as expm
import numpy as np
import numpy.linalg

mmax=10

n=4
td = -2*np.ones([mmax], dtype=np.double)
tsd = np.ones([mmax-1],dtype=np.double)

beta = 1.0
dt = 0.01

expT = expm.expimv_tridiag(mmax)
expT.precompute(td,tsd,n)
y=expT.eval(dt,beta)

print(len(y))

print(np.abs(y))
expH = expm.expimv(40,mmax)

H = np.zeros([mmax,mmax])
for j in range(n): 
  H[j,j]=-2.0
for j in range(n-1): 
  H[j+1,j]=1.0
  H[j,j+1]=1.0

x = np.zeros([n])
x[0]=beta

y2 = expH.apply(dt,H,x,n)
print(y-y2)

print(numpy.linalg.norm(y)) 
print(numpy.linalg.norm(y2)) 

# python3 -m numpy.f2py -c F_smallexpimv.F90 -m smallexpimv -lblas
# python3 -m numpy.f2py -h smexp.pyf -m smallexpimv F_smallexpimv.F90
