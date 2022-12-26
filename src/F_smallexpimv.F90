module smallexpimv

implicit none

private
public expimv, expimv_tridiag_precompute, expimv_tridiag_apply
contains

subroutine expimv( dt, H, ldh, m, Y, pmone, kmax, tol, nrmest, maxrestart, dtdone, info )
    ! compute exp(sig * dt * h) * Y(1:muse) -> Y(1:muse)
    ! Y is of size 'muse x (kmax + 1)'
    ! for a small, dense matrix h, dim muse (submatrix of ks%mmax)
    ! norm(y) = 1
    ! norm(H) approx nrmest
    double precision, intent(in) :: dt, tol, nrmest
    integer, intent(in) :: pmone, ldh, m, kmax, maxrestart
    complex*16, intent(in), dimension(ldh, ldh) :: H
    complex*16, intent(inout), dimension(ldh * (kmax+1)) :: Y

    double precision, intent(out) :: dtdone
    integer, intent(out) :: info

    integer :: inorm
    integer :: k, j, taylorrun
    double precision :: tnow, tend, dtguess, dt1, dt2, dtx, nrmy
    double precision :: qnew, qold, normtest, normtest1, rtol
    integer :: kuse
    logical :: stoptaylor
    complex*16 :: ytemp

    external :: zgemv, zaxpy, zscal
    double precision, external :: dznrm2
    integer, external :: izamax

    complex*16 :: sig

    sig = dcmplx(0.0d0, pmone*1.0d0)
    info = 1
    dtdone = dt
    
    nrmy = dznrm2(m, Y(1),1)

    rtol = tol*nrmy*nrmest    ! scale tolerance by size of entries and beta
    
    stoptaylor = .false.

    tnow = 0.0
    tend = tnow + dt

    taylorrun = 0

    ! scalar case, exponential of matrix with dimension 1
    if (m .eq. 1) then
        ytemp = Y(1)
        Y(1) = exp(sig * dt * h(1,1)) * ytemp
    else

outloop: do while (tnow < tend)
    taylorrun = taylorrun + 1

    dtguess = tend-tnow

    kuse = kmax
    qnew = 1
precomputemv: do k = 1, kmax
        call zgemv('N', m, m, sig * dcmplx(dtguess/k, 0.0), H, &
        ldh, Y(1 + (k-1)*m), 1, dcmplx(0.0, 0.0), Y(1 + k*m), 1)

        qold = qnew
        inorm = izamax(m, Y(1 + k*m), 1)
        qnew = abs(Y(k*m + inorm))

        if ( qold + qnew < rtol ) then
            kuse = k
            stoptaylor = .true.
            exit precomputemv
        endif
    end do precomputemv

    if (stoptaylor) then
        dtx=1.0
    else
        dt1=(rtol/(2*qnew))**(1./(kuse))
        dt2=(rtol/(2*qold))**(1./(kuse-1))
        dtx = min(min(dt1,dt2),1.0)
        call zscal(m,dcmplx(dtx**kuse, 0.0), Y(1+kuse*m),1)! first scaling step for summationtaylor loop
    endif

summationtaylor: do k = 1, kuse-1
        call zaxpy(m, dcmplx(dtx**(kuse-k), 0.0), Y(1+(kuse-k)*m), 1, Y(1+kuse*m), 1)
    end do summationtaylor
    call zaxpy(m, dcmplx(1.0, 0.0), Y(1+kuse*m), 1, Y(1), 1)

    tnow = tnow + dtguess*dtx

    if ((taylorrun > maxrestart)) then
        write(*, *) "warning, extensive number of iterations for exponential of small matrix, iteration",taylorrun
    endif

end do outloop

endif


end subroutine

subroutine expimv_tridiag_precompute(td, tsd, w, z, n, ldz, isuppz, work, iwork, info)
    ! compute eigenvalues and eigenvectors of tridiagonal jacobi matrix with entries saved in workd
    ! as a result, eigenbasis z and eigenvalues w are also saved in workd
    implicit none

    integer :: infoeigs
    integer :: m, lwork, liwork

    integer, intent(in) :: ldz, n
    double precision, intent(inout), dimension(ldz, ldz) :: Z
    double precision, intent(inout), dimension(ldz) :: w
    double precision, intent(in), dimension(ldz) :: td
    double precision, intent(in), dimension(ldz-1) :: tsd

    double precision, intent(in), dimension(20 * ldz) :: work
    integer, intent(in), dimension(10 * ldz) :: iwork
    integer, intent(in), dimension(2*ldz) :: isuppz
    integer, intent(out) :: info

    ! not sure if tsd(n) is overwritten/used as workspace ???

! for krylov, ldz = mmax, n = m
! Z is matrix of dimension (ldz, n) <-> (m, mmax) for krylov

    external :: dstevr
    !double precision, external :: dznrm2

    lwork = 20 * ldz
    liwork = 10 * ldz

    if (n > 1) then
        call dstevr('V', 'A', n, td, tsd, 0, 0, 0, 0, 0, &
        m, w, z, ldz, isuppz, work, lwork, iwork, liwork, infoeigs)
        if ( infoeigs .ne. 0 ) then
            write(*,*) " dstevr returns error !!!! (error on eigendecomposition for small exponential) "
        endif
        info = infoeigs
    else
        w(1) = td(1)
        z(1,1) = 1.0
        info = 1
    endif

end subroutine

subroutine expimv_tridiag_apply(dt, pmone, beta, w, z, x, n, ldz, cwork)
    ! use eigenbasis of a tridiagonal real symmetric matrix J (previously computed by expimv_tridiag_precompute)
    ! to compute beta*exp(-i*dt*J)*e1
    ! J = Z * diag(w) * Z.t
    ! w is vector of eigenvalues
    ! x .. output x = beta*exp(pm*i*dt*J)*e1
    ! cwork .. work vector
    ! pm = 1 or -1
    implicit none

    double precision, dimension(ldz, ldz), intent(in) :: Z

    double precision, dimension(ldz), intent(in) :: w 
    complex*16, dimension(n), intent(inout) :: x
    complex*16, dimension(ldz), intent(inout) :: cwork

    double precision, intent(in) :: dt, beta
    integer, intent(in) :: pmone, n, ldz
    integer :: k

    complex*16 sig

    external :: zgemv

    sig = dcmplx(0.0d0, pmone*1.0d0)
    if (n.eq.1) then
        ! z(1,1) = 1
        x(1) = dcmplx(beta, 0.0) * exp( sig * dcmplx(dt * w(1), 0.0))
    else
        do k = 1, n
            cwork(k) = exp( sig * dcmplx(dt * w(k), 0.0)) * dcmplx(Z(1,k), 0.0d0) ! Z(1 + (k - 1) * ldz)
        end do
        call zgemv('n', n, n, dcmplx(beta, 0.0), dcmplx(Z, 0.0), ldz, cwork, 1, dcmplx(0.0, 0.0), x, 1)
    end if

end subroutine

end module smallexpimv

