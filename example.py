# python3 code
import smallexpimv_pyclasses as expm
import numpy as np
import numpy.linalg

# define a maximal dimension for matrix H,
# required to setup working memory in case of multiple H with different dimensions
nmax=40

# setup example matrix
n=30

H = np.zeros([nmax,nmax])
for j in range(n): 
  H[j,j]=-2.0
for j in range(n-1): 
  H[j+1,j]=1.0
  H[j,j+1]=1.0
# for tridiagonal setup
td = -2*np.ones([nmax], dtype=np.double)
tsd = np.ones([nmax-1],dtype=np.double)

# scaling factor and time-step
beta = 1.0
dt = 1.0

# setup solver for imaginary exponential for tridiagonal matirx
expT = expm.expimv_tridiag(nmax)
# solve
expT.precompute(td,tsd,n)
y=expT.eval(dt,beta)

# setup solver for imaginary exponential for general matrix
ntaylor = 40
expH = expm.expimv(ntaylor,nmax)
#solve
x = np.zeros([n])
x[0] = beta
y2 = expH.apply(dt,H,x,n)

# compare results
erry=numpy.linalg.norm(y-y2)
print(f"comparing both methods, |y-y2| = {erry}")

# imaginary exponential conserves mass, check
ern1=abs(numpy.linalg.norm(y) - beta)
ern2=abs(numpy.linalg.norm(y2) - beta)
print(f"error in norm of solution, ||y|-b| = {ern1} and ||y2|-b|2 = {ern2}")
