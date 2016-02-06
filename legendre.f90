! =================================================================================================
! =================================================================================================
! ==========================                                           ============================
! ==========================     Legendre Polynomial Evaluation        ============================
! ==========================                                           ============================
! =================================================================================================
! =================================================================================================
!
! Author: Christopher A. Wong
! Date: 6 February 2016
! Format: Fortran 2008
!
! This module contains the procedures evaluate the Legendre polynomials and their derivatives up to
! a given degree on the d-dimensional unit cube.
!
! Dependencies: LAPACK
!
! Externally called procedures:
!   leg_buildpars   --- Initial subroutine called to form Legendre polynomial parameters
!   leg_eval        --- Main subroutine for evaluating Legendre polynomials.
!   leg_quad2       --- Computes Legendre matrix for 2D bilinear quadrature computations.
!
! Internal procedures:
!   choose          --- Computes binomial coefficients.


module legendre

    implicit none
    private
    integer,    parameter                       :: dp = kind(1.0d0)
    integer                                     :: maxdeg, d, numfuncs

    public  :: leg_eval, leg_quad2

contains

! =================================================================================================

! Binomial coefficient function

pure integer function choose (n, k)

    implicit none
    integer,    intent(in)  :: n
    integer,    intent(in)  :: k
    integer                 :: i,kk

    if (n-k < k) then
        kk = n-k
    else
        kk = k
    end if

    choose = 1

    if (kk >=1) then
        do i = 1, kk
            choose = (choose * (n-kk+i))/i
        end do
    end if

end function choose

! =================================================================================================

subroutine leg_buildpars(deg_in,d_in)

    implicit none

    integer,    intent(in)  :: deg_in, d_in

    maxdeg = deg_in
    d = d_in
    numfuncs = choose(maxdeg+d,d)

end subroutine leg_buildpars

! =================================================================================================
!
! Legendre polynomial evaluation
!
! Inputs:
!   x           --- input point as a dx1 array
!   need_deriv  --- optional, set to .TRUE. to compute derivatives, default is .FALSE.
!
! Outputs:
!   func_out    --- evaluations of multivariate orthogonal polynomials at x
!   deriv_out   --- values of the 1st order partial derivatives of the multivariate orthogonal
!       polynomials at x, given as a matrix with rows as functions
!   flag        --- set to .true. if an error occurs

subroutine leg_eval(x, need_deriv, func_out,deriv_out,flag)

implicit none

real(dp),   dimension(:),   intent(in)  :: x
logical,    optional,       intent(in)  :: need_deriv
real(dp),   dimension(:),   intent(out) :: func_out
real(dp),   dimension(:,:), intent(out) :: deriv_out
logical,                    intent(out) :: flag

logical                                 :: deriv_flag   ! Determines if derivatives are computed
integer                                 :: i, ii, j, jj, k, kk ! loop variables

integer,    dimension(1:choose(maxdeg+d,d),1:d) :: degidx ! stores degrees of multivar. polys
real(dp),   dimension(0:maxdeg,1:d)             :: leg1D, legder1D ! stores 1d Leg poly evaluations

! Initialize
func_out    = 1.0d0
deriv_out   = 1.0d0
flag        = .false.

if (present(need_deriv)) then
    deriv_flag = need_deriv
else
    deriv_flag = .false.
end if

! Check whether inputs have the correct sizes
if (size(func_out) /= numfuncs) then
    flag = .true.
else if ( (deriv_flag) .and. (size(deriv_out,1) /= numfuncs)) then
    flag = .true.
else if ( (deriv_flag) .and. (size(deriv_out,2) /=d) ) then
    flag = .true.
end if

! form the list of polynomial degrees in the tensor product of legendre polys
! this is equivalent to forming the set of all ordered partitions of positive
! integers up to maxdeg

if (d == 1) then
    degidx(1:numfuncs,1) = [ (j, j = 0, numfuncs-1) ]
else if (d == 2) then
    k = 1                               ! index counter
    do i = 0, maxdeg                    ! total degree
        do j = 0,i                      ! degree in first variable
            degidx(k,1:2) = [ j,i - j ]
            k = k + 1
        end do
    end do
else if (d == 3) then
    kk = 1                              ! index counter
    do i = 0, maxdeg                    ! total degree
        do j = 0, i                     ! total degree of the first two terms
            do k = 0, j                 ! degree in the first variable
                degidx(kk,1:3) = [ k, j-k, i-j ]
                kk = kk + 1
            end do
        end do
    end do

else
    flag = .true.
end if

! Evaluate 1D Legendre polynomials and derivatives at the d coordinates of x
! using the three-term recurrence formula for orthogonal polynomials

! set values for degrees 0 and 1
leg1D(0,1:d)    = 1.0d0
leg1D(1,1:d)    = x


if (deriv_flag) then
    legder1D(0,1:d) = 0.0d0
    legder1D(1,1:d) = 1.0d0
end if

do concurrent (i = 1:d) ! looping through d dimensions
    do j = 1, maxdeg -1 ! looping through degrees

        leg1D(j+1,i) =( (2*j+1) * x(i) * leg1D(j,i) - j*leg1D(j-1,i) ) / (j + 1.0d0)

        if (deriv_flag) then
            legder1D(j+1,i) = ( (2*j+1) * leg1D(j,i) + (2*j+1) * x(i) * legder1D(j,i) - &
                j* legder1D(j-1,i) ) / (j + 1.0d0)
        end if

    end do
end do

! renormalization
if (deriv_flag) then
    do concurrent (i = 0:maxdeg)
        leg1D(i,1:d)    = leg1D(i,1:d) * sqrt((2.0d0*i +1)/2)
        legder1D(i,1:d) = legder1D(i,1:d) * sqrt((2.0d0*i +1)/2)
    end do
else
    do concurrent (i = 0:maxdeg)
        leg1D(i,1:d)    = leg1D(i,1:d) * sqrt((2.0d0*i +1)/2)
    end do
end if


! Evaluate d-dimensional Legendre polynomials

if (d == 1) then
    func_out(1:numfuncs) = leg1D(0:numfuncs-1,1)

    if (deriv_flag) then
        deriv_out = legder1D
    end if
else
    do i = 1, numfuncs ! loop over all functions
        do j = 1,d ! loop over all dimensions
            func_out(i) = func_out(i) * leg1D(degidx(i,j),j)

            if (deriv_flag) then
                do k = 1,d ! loop over all partial derivatives

                    if (j /= k) then
                        deriv_out(i,k) = deriv_out(i,k) * leg1D(degidx(i,j),j)
                    else if (j == k) then
                        deriv_out(i,k) = deriv_out(i,k) * legder1D(degidx(i,j),j)
                    end if

                end do
            end if
        end do
    end do
end if

end subroutine leg_eval

! =================================================================================================
!
! This special subroutine is only used for computing bilinear quadrature that is exact on a space
! of polynomials and minimal on polynomials of the next degree. It is designed to be called with
! the NLopt optimization library.
!
subroutine leg_quad2(x,numpt,A_out)

    implicit none
    real(dp),   dimension(:),                   intent(in)  :: x        ! input array of points
    integer,                                    intent(in)  :: numpt    ! number of points
    real(dp),   dimension(:,:),                 intent(out) :: A_out    ! output matrix

    real(dp),   dimension(1:numpt,1:numfuncs)               :: output
    real(dp),   dimension(1:numpt,1:choose(maxdeg+1,2))     :: F
    real(dp),   dimension(1:numpt,1:(maxdeg+1))             :: G
    integer,    dimension(1:numpt)                          :: IPIV     ! LAPACK pivot info
    integer                                                 :: INFO     ! LAPACK error value
    logical                                                 :: flag

    ! ---------

    ! Check that a symmetric bilinear quadrature has been requested
    if (numpt /= (numfuncs - maxdeg - 1)) then
        print *, 'Error: Non-symmetric bilinear quadrature requested.'
        stop
    end if


    ! Evaluate matrix functions
    call leg_eval(reshape(x, (/ 2, numpt /) ), need_deriv = .false., output, 0.0d0, flag)
    F = output(:,1:(numfuncs - maxdeg -1))
    G = output(:,(numfuncs - maxdeg):numfuncs)

    ! Compute objective function A(x) = inv(F(x)) * G(x) using LAPACK routines
    call dgetrf(numpt,numpt,F,numpt,IPIV,INFO)
    call dgetrs('N',numpt,maxdeg+1,F,numpt,IPIV,G,numpt,INFO)

    A_out = G

end subroutine leg_quad2

! =================================================================================================



end module legendre
