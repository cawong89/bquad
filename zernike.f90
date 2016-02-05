! =================================================================================================
! =================================================================================================
! ==========================                                           ============================
! ==========================      Zernike Polynomial Evaluation        ============================
! ==========================                                           ============================
! =================================================================================================
! =================================================================================================
!
! Author: Christopher A. Wong
! Date: 26 January 2016
! Source format: Fortran 2008
!
! This module contains procedures to evaluate the Zernike polynomials on the unit disk.
!
! Dependencies: LAPACK
!
! Externally called procedures:
!   zern_eval       --- Main Zernike polynomial evaluation routine. Call this externally.
!   zern_buildpars  --- Builds all parameters needed to call the main routine. Call first.
!   zern_bquad      --- Computes Zernike matrix for bilinear quadrature computations.
!
! Internal procedures:
!   getnumrad       --- Computes the number of radial polynomials.
!   Kconst          --- Computes recursive constants for q-recursive algorithm.
!   buildK          --- Builds all recursive constants.
!   idxwrap         --- Computes linear index of R(m,n)
!   idxwrap2        --- Computes linear index of Z(m,n)
!   shiftdown       --- Computes linear index of R(n,n+2)
!   cart_to_polar   --- Conversion from cartesian to polar.
!   Radial_eval     --- Evaluates radial polynomials.
!
! =================================================================================================
!
! Usage:
!
! To evaluate Zernike polynomials up to a maximum degree d, first use
!
!       call zern_buildpars(d)
!
! to precompute the recursion coefficients needed to compute the polynomials. Then, given any sets
! of points x, you can use
!
!       call zern_eval(x,numpt,d,output)
!
! where 'numpt' is the number of cartesian points in x, to calculate the output array of polynomial
! evaluations. If you wish to calculate the Zernike polynomials up to a new maximum degree, you must
! run zern_buildpars again with the new value of d before you start evaluating.
! =================================================================================================

module zernike

implicit none
private

integer,    parameter                       :: dp = kind(1.0d0)
integer                                     :: numrad, numfunc, maxdeg
real(dp),   dimension(:), allocatable       :: K1, K2, K3
save

public  :: zern_eval, zern_buildpars, zern_bquad

contains

! =================================================================================================

! This outer subroutine evaluates all Zernike polynomials up to a degree 'maxdeg' at a set of
! points x. This routine automatically evaluates the polynomials if the radial component is zero,
! and passes to the main subroutine Zern_eval for nonzero radial components.

subroutine zern_eval(x,numpt,func_out)

implicit none

! I/O Variables
real(dp),   dimension(:,:),             intent(in)      :: x
integer,                                intent(in)      :: numpt
real(dp),   dimension(:,:),             intent(out)     :: func_out

! Local variables
integer                                 :: m,n,p
real(dp),   dimension(1:numpt)          :: r, theta
real(dp),   dimension(1:numpt,1:numrad) :: Rad  ! Radial polynomial array
real(dp),   dimension(1:numpt,1:maxdeg) :: cosarray, sinarray ! Stores values of cos(m*theta), etc

! Convert to polar coordinates
call cart_to_polar(x,r,theta)

! Compute radial polynomials
Rad = 0.0d0
call Radial_eval(r,numpt,Rad)

! Precompute sin(m*theta), cos(m*theta) using elemental functions

do m = 1,maxdeg
    cosarray(:,m) = cos(m*theta)
    sinarray(:,m) = sin(m*theta)
end do

! Compute Zernike polynomials

do n = 0,maxdeg
    do p = 0,2*n,2
        m = p - n
        if (m < 0) then
            func_out(:,idxwrap2(m,n)) = Rad(:,idxwrap(-m,n)) * sinarray(:,-m)
        else if (m > 0) then
            func_out(:,idxwrap2(m,n)) = Rad(:,idxwrap(m,n)) * cosarray(:,m)
        else if (m == 0) then
            func_out(:,idxwrap2(m,n)) = Rad(:,idxwrap(m,n))
        end if
    end do
end do



end subroutine zern_eval



! =================================================================================================

subroutine zern_buildpars(maxdeg_in)

! Computes and stores the following parameters within the modules scope:
! numrad, numfunc, maxdeg, K1, K2, K3
! This is used to rapidly evaluate at many different input points while under the same parameters

    implicit none

    integer,    intent(in)  :: maxdeg_in
    integer                 :: ksize

 ! Precompute parameters

    maxdeg = maxdeg_in
    numrad = getnumrad(maxdeg)
    numfunc = (maxdeg+1) * (maxdeg+2)/2

    if (mod(maxdeg,2) == 0) then
        ksize = (maxdeg/2 + 1)**2
    else
        ksize = (maxdeg+3) * (maxdeg+1) / 4
    end if

! Allocate q-recursive array

    if (allocated(K1) .eqv. .true.) then
        deallocate(K1, K2, K3)
    end if

    allocate(K1(1:ksize),K2(1:ksize),K3(1:ksize))

    call buildK(maxdeg)


end subroutine zern_buildpars

! =================================================================================================

! This special routine is used only for computing bilinear quadratures. Given degree "maxdeg", this
! calculates a matrix A(x) = inv(F(x)) * G(x), where F(x) is the Zernike matrix up to degree
! 'maxdeg' - 1 and [F(x) G(x)] is the Zernike matrix up to degree 'maxdeg'.

subroutine zern_bquad(x,numpt,A_out)

    implicit none

    real(dp),   dimension(:),                   intent(in)  :: x        ! input array of points
    integer,                                    intent(in)  :: numpt    ! number of points
    real(dp),   dimension(:,:),                 intent(out) :: A_out    ! output matrix

    real(dp),   dimension(1:numpt,1:numfunc)                :: output
    real(dp),   dimension(1:numpt,1:(numfunc - maxdeg - 1)) :: F
    real(dp),   dimension(1:numpt,1:(maxdeg+1))             :: G
    integer,    dimension(1:numpt)                          :: IPIV     ! LAPACK pivot info
    integer                                                 :: INFO     ! LAPACK error value

! ---------

    ! Check that a symmetric bilinear quadrature has been requested
    if (numpt /= (numfunc - maxdeg - 1)) then
        print *, 'Error: Non-symmetric bilinear quadrature requested.'
        stop
    end if


    ! Evaluate matrix functions
    call zern_eval(reshape(x, (/ 2, numpt /) ),numpt,output)
    F = output(:,1:(numfunc - maxdeg -1))
    G = output(:,(numfunc - maxdeg):numfunc)

    ! Compute objective function A(x) = inv(F(x)) * G(x) using LAPACK routines
    call dgetrf(numpt,numpt,F,numpt,IPIV,INFO)
    call dgetrs('N',numpt,maxdeg+1,F,numpt,IPIV,G,numpt,INFO)

    A_out = G

end subroutine zern_bquad


! =================================================================================================

pure integer function getnumrad(n)

! Computes the number of radial polynomials up to degree n

    implicit none

    integer, intent(in) :: n

    if (mod(n,2) == 1) then
        getnumrad = (n+3)*(n+1)/4
    else
        getnumrad = (n/2 + 1)**2
    end if

end function getnumrad

! =================================================================================================

subroutine Kconst(K1i,K2i,K3i,m,n)

! Computes special constant parameters in recursive Zernike evaluation
! Inputs:
!   m,n         --- Integer, m < n, n - m even
! Outputs:
!   K1i,K2i,K3i    --- Double, special recursive constants in Q-Recursive Algorithm

implicit none

integer,    intent(in)  :: m,n              ! Input values
real(dp),   intent(out) :: K1i, K2i, K3i    ! Output coefficients


K3i = - 4 * real((m+2)*(m+1),dp) / real((n+m+2)*(n-m),dp)

K2i = K3i * real((n+m+4)*(n-m-2),dp) / real(4*(m+3),dp) + m + 2

K1i = real((m+4)*(m+3)/2,dp) - (m+4) * K2i + K3i/8 * real((n+m+6)*(n -m - 4),dp)

end subroutine Kconst
! =================================================================================================

subroutine buildK(maxdeg_in)

! Precomputes the matrix of special coefficients K(m,n) in the recursive Zernike evaluation
! If n is odd then each K has (n+3)(n+1)/4 elements
! If n is even then each K has (n/2 + 1)^2 elements

implicit none

integer,                        intent(in)  :: maxdeg_in       ! Input, highest order

! Local variables
integer                                     :: m,n,p,j

do n = 0,maxdeg_in
    do p = 2,n,2
        m = n - p
        j = idxwrap(m,n)
        call Kconst(K1(j),K2(j),K3(j),m,n)
    end do
end do

end subroutine

! =================================================================================================

pure integer function idxwrap(m,n)

! Computes the linear index of R(m,n) and K(m,n)

implicit none

integer, intent(in) :: m,n

if (mod(n,2) == 1) then
    idxwrap = (n+3) * (n+1)/4 - (n-m)/2
else
    idxwrap = (n/2 + 1)**2 - (n-m)/2
end if

end function idxwrap

! =================================================================================================

pure integer function idxwrap2(m,n)

! Computes the linear index of Z(m,n); n - m even

implicit none

integer, intent(in) :: m,n

idxwrap2 =  (n*(n+1) + (n + m + 2))/2

end function idxwrap2

! =================================================================================================

pure integer function shiftdown(n)

! Computes the linear index of R(n,n+2)

implicit none

integer, intent(in) :: n

if (mod(n,2) == 1) then
    shiftdown = n + 2 + (n + 3)/2 * (n + 1)/2
else
    shiftdown = n + 2 + (n/2 + 1)**2
end if

end function shiftdown

! =================================================================================================

subroutine cart_to_polar(x,r,theta)

implicit none

real(dp),   dimension(:,:), intent(in)  :: x        ! Input value in R^2 cartesian coordinates
real(dp),   dimension(:),   intent(out) :: r,theta  ! Output in polar coordinates

r = sqrt(x(1,:)**2 + x(2,:)**2)

theta = atan2(x(2,:),x(1,:))

end subroutine cart_to_polar


! =================================================================================================
! ================================== Radial Polynomial Subroutine =================================
! =================================================================================================

! This subroutine evaluates all radial polynomials up to degree 'maxdeg' at a specified point x.
! For n - m even, the Zernike polynomials are given by
!
!       Z(m,n,r,theta) = Rad(m,n,r) * cos(m * theta)    m >= 0
!       Z(m,n,r,theta) = Rad(m,n,r) * sin(m * theta)    m < 0
!
! The so-called "q-recursion" on Rad(m,n,r) is given by
!
!       Rad(n,n,r) = r^n
!       Rad(n-2,n,r) = n * Rad(n,n,r) + (n-1)* Rad(n-2,n-2,r)             n >= 2
!       Rad(m,n,r) = K_1 * Rad(m+4,n,r) + (K_2 + K_3/r^2) * Rad(m+2,n,r)  n - m > 2
!
! where the coefficients K_1, K_2, K_3 are given by
!
!       K_1 = (m+4)(m+3)/2 - (m+4) K_2 + K_3 (n + m + 6)(n - m - 4)/8
!       K_2 = K_3 (n + m + 4) (n - m - 2)/(4(m+3)) + m + 2
!       K_3 = - 4(m + 2)(m+1)/ ((n + m + 2)(n - m))

! This method *might* fail for extremely small r. Therefore for sufficiently small r we use the
! Prata recursive method to estimate the radial polynomial.

subroutine Radial_eval(r, numpt, Rad)

implicit none

! I/O Variables
integer,                                intent(in)      :: numpt        ! Number of input points
real(dp),   dimension(1:numpt),         intent(inout)   :: r            ! input radial values
real(dp),   dimension(:,:),             intent(out)     :: Rad          ! k x numrad output matrix

! Local Variables
real(dp),   dimension(1:numpt,1:numrad)             :: Radtemp      ! Sorted evaluation matrix
logical,    dimension(1:numpt)                      :: rnz          ! Is TRUE for a nonzero radius
real(dp)                                            :: tolr2        ! Tolerance for r^2
integer                                             :: j,n,m,p      ! loop variables
real(dp)                                            :: rtemp
integer                                             :: nnz          ! # of nonzero radial values

! Determine nonzero radial values and permute near-zero ones to the back
tolr2 = 10.0d-30
rnz = .true.
nnz = numpt

do j = numpt,1,-1 ! Permutes from back to front; this is optimal
    if (r(j)**2 <= tolr2) then
        rnz(j) = .false.
        nnz = nnz - 1
        rtemp = r(j)
        r(j:(numpt-1)) = r((j+1):numpt)
        r(numpt) = rtemp
    end if
end do

! Compute Rad(n,n,r)
do n = 0,maxdeg
    Radtemp(:,idxwrap(n,n)) = r**n
end do

! ------ Prata method for near-zero radial values ------

! Compute Rad(0,n,r) for n even, using a two-term expansion (since r << 1)
do n = 0,maxdeg,2
    j = n/2
    Radtemp(nnz+1:numpt,idxwrap(0,n)) = (-1)**(mod(j,2)) * &
                                    (1 - ((j+1)* r(nnz+1:numpt))*(j * r(nnz+1:numpt)))
end do

! Prata recursion for Rad(m,n,r), 0 < m < n
do n = 1,maxdeg
    do p = 2,n,2
        m = n - p
        j = idxwrap(m,n)
        Radtemp(nnz+1:numpt,j) = real(2*n,dp)/(n + m) * r(nnz+1:numpt) * &
                                Radtemp(nnz+1:numpt,idxwrap(m-1,n-1)) + &
                                real(m - n,dp)/(n + m) * Radtemp(nnz+1:numpt,j - n)
    end do
end do



! ------ Q-Recursive algorithm for nonzero radial values ------



! Compute Rad(n-2,n,r)
do n = 2,maxdeg
    j = shiftdown(n-2)
    Radtemp(1:nnz,j) = (1 - n) * Radtemp(1:nnz,j - n) + n * Radtemp(1:nnz,j+1)
    !print *, Radtemp(1:nnz,j - n)
    !print *, Radtemp(1:nnz,j+1)
end do

! Compute Rad(m,n,r) for n - m > 2
do n = 4, maxdeg
    do p = 4,n,2
        m = n - p
        j = idxwrap(m,n)
        Radtemp(1:nnz,j) = K1(idxwrap(m,n)) * Radtemp(1:nnz,j+2) + &
                        K2(idxwrap(m,n)) * Radtemp(1:nnz,j+1) + &
                        K3(idxwrap(m,n)) * Radtemp(1:nnz,j+1)/(r(1:nnz)**2)
    end do
end do

! ------ Re-sort output matrix ------
do j = 1,numpt
    if (rnz(j) .eqv. .false.) then
        Rad(j,:) = Radtemp(numpt,:)
        Radtemp((j+1):numpt,:) = Radtemp(j:(numpt-1),:)
    else
        Rad(j,:) = Radtemp(j,:)
    end if
end do

end subroutine Radial_eval

end module zernike
