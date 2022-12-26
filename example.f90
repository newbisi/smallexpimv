program example
  use smallexpimv

  implicit none
  integer :: n, nmax, ntaylor, k

  double precision :: beta, dt

  double precision, pointer, dimension(:) :: td, tsd
  complex*16, pointer, dimension(:,:) :: Hmat
  complex*16, pointer, dimension(:) :: cwork, Y, v, u, uer

  integer, pointer, dimension(:) :: iwork, isuppz

  double precision, pointer, dimension(:) :: work, w, z

  double precision :: tol, nrmest, dtdone, unrm, x1,x2
  integer :: pmone, maxrestart , info

  double precision, external :: dznrm2

  pmone = -1

  nmax = 40
  n = 30
  ntaylor = 40
  tol = 10E-15
  maxrestart = 200
  nrmest = 1.0

  beta = 1.0
  dt = 1.0

  allocate(Hmat(nmax,nmax), td(nmax), tsd(nmax-1), u(nmax), v(nmax), uer(nmax))
  allocate(work(20*nmax), w(nmax), z(nmax*nmax))
  allocate(cwork(nmax), iwork(10*nmax), isuppz(2*nmax))
  allocate(Y(nmax*(ntaylor+1)))

  ! define tridiagonal test matrix
  do k=1, nmax
    td(k) = -2.0
  end do
  do k=1, nmax-1
    tsd(k) = 1.0
  end do

  call expimv_tridiag_precompute(td, tsd, w, z, n, nmax, isuppz, work, iwork, info)
  call expimv_tridiag_apply(dt, pmone, beta, w, z, v, n, nmax, u)

  ! define test matrix
  do k=1, nmax
    Hmat(k,k) = dcmplx(-2.0,0.0)
  end do
  do k=1, nmax-1
    Hmat(k+1,k) = dcmplx(1.0,0.0)
    Hmat(k,k+1) = dcmplx(1.0,0.0)
  end do

  ! define initial vector
  do k=1, nmax
    Y(k) = 0.0
  end do
  Y(1) = beta

  call expimv( dt, Hmat, nmax, n, Y, pmone, ntaylor, tol, nrmest, maxrestart, dtdone, info )
  ! solution is in Y(1),..,Y(n)

  ! comparing both solutions and norm of solutions
  do k=1, n
    uer(k) = Y(k)-v(k)
  end do

  unrm =dznrm2(n, uer(1),1)
  write(*,*) "comparing different methods, error in solutions =", unrm

  x1 = abs(dznrm2(n, v(1),1)-beta)
  x2 = abs(dznrm2(n, Y(1),1)-beta)
  write(*,*) "error in norm x1 =", x1, "x2=", x2

  deallocate(Hmat, u, v, uer)
  deallocate(work, w, z, cwork, iwork,isuppz)
  deallocate(Y)

end program example

