! =================================================================================================
! =================================================================================================
! ========================                                               ==========================
! ========================   SYMMETRIC INVERTIBLE BILINEAR QUADRATURE    ==========================
! ========================                                               ==========================
! =================================================================================================
! =================================================================================================
!
!
! Author: Christopher A. Wong
! Date: 7 March 2016
! Source: Fortran 2008
!
! Dependencies: LAPACK, NLopt
!
! This module contains the necessary routines to compute a stable and accurate symmetric bilinear
! quadrature rule that is exact for a continuous inner product on a prescribed finite-dimensional
! subspace of continuous functions.
!
! Externally called procedures:
!   bquad_buildpars
!   bquad_construct
!
! Internal procedures:
!   sigmamax
!	ptperturb
!

! =================================================================================================
!                                       INTERFACE MODULE
! =================================================================================================
!
! This module specifies the interface for the subroutine that computes the matrix function
!                               A(x) = inv(F(x)) * G(x)
! where F(x), G(x) are the collocation matrices for the function space to be integrated exactly and
! the orthogonal function space to be minimized, respectively.

module bquad_interface

    implicit none
    private
    integer,    parameter                       :: dp = kind(1.0d0)


    abstract interface
        subroutine eval(x,numpt,output)
            import dp
            real(dp),   dimension(:),   intent(in)  :: x
            integer,                    intent(in)  :: numpt
            real(dp),   dimension(:,:), intent(out) :: output
        end subroutine eval

    end interface

    public  :: eval

end module bquad_interface



module bquad

    use bquad_interface
    implicit none
    private
    save
    integer,    parameter                       :: dp = kind(1.0d0)
    integer                                     :: numpt, dm, numg, mink, LWORK

    public  :: bquad_buildpars, bquad_construct

    contains

! =================================================================================================

! Sets the following module-global parameters
!   numpt   --- number of points = dimension of function space on which quadrature is exact
!   dm      --- dimension of ambient space
!   numg    --- Dimension of orthogonal function space on which bilinear quadrature is minimal
!	mink	--- Minimum between numpt and numg
!	LWORK	--- Optimal amount of memory required for SVD calculation.

subroutine bquad_buildpars(numpt_in, dm_in, numg_in)

    implicit none
    integer,    intent(in)  :: numpt_in, dm_in, numg_in

    numpt = numpt_in
    dm = dm_in
    numg = numg_in
    mink = min(numpt,numg)
    LWORK = max(1, 3*mink + max(numpt,numg), 5*mink)

end subroutine bquad_buildpars

! =================================================================================================

subroutine bquad_construct(n, A_in, x, W, sigma_out, time)

	implicit none
	integer,					intent(in)		:: n
	procedure(eval)								:: A_in
	real(dp),	dimension(:),	intent(inout)	:: x
	real(dp),	dimension(:,:),	intent(out)		:: W
	real(dp),					intent(out)		:: sigma_out
	real(dp),	optional,		intent(out)		:: time
	
	real(dp)									:: eps
	integer										:: AlgIter, OptIter, MaxOptIter
	integer,	dimension(1:4)					:: alg_list
	logical										:: time_flag
	
	! 12 - Praxis, 25 - BOBYQA, 28 - Nelder Mead, 29 - SBPLX
	alg_list = (/ 28, 29, 25, 12 /)
	
	! Set parameters
	MaxOptIter = 10
	eps = 5.0d-3
	if (present(time)) then
		time_flag = .TRUE.
	else
		time_flag = .FALSE.
	end if
	
	! Iterate over different algorithms
	do OptIter = 1,MaxOptIter	
	
		call ptperturb(x,eps)
		
		do AlgIter = 1,4
			call optimize(n, x, A_in, time_flag, sigma_out, alg_list(AlgIter))
		end do
		
	end do
	
	! Construct weight matrix
	W = 1	! Temporary
	
	print *, 'Final min. sigma: ', sigma_out
	

end subroutine bquad_construct

! =================================================================================================

! Calls an external optimization library to minimize the bilinear quadrature objective function
!
!   Inputs:
!       n           --- number of variables; equals numpt * dm
!       x           --- initial guess value, DIM(1:n)
!       A_in        --- subroutine that computes A(x) = inv(F(x) * G(x))
!		time 		--- logical, if TRUE indicates that runtime should be output
!       alg         --- optional, integer identifier for NLopt optimization algorithm
!   Outputs:
!       x           --- bilinear quadrature evaluation points, DIM(1:n)
!		sigma		--- Minimized singular value


subroutine optimize(n, x, A_in, time_flag, sigma, alg)

    implicit none

    integer,                    intent(in)      :: n
    procedure(eval)                             :: A_in
    real(dp),   dimension(:),   intent(inout)   :: x
    integer,    optional,       intent(in)      :: alg
    logical,					intent(in)		:: time_flag
	real(dp),					intent(out)		:: sigma
    

    real(dp),   parameter                       :: tol_rel = 1.0d-6
    real(dp),   parameter                       :: tol_abs = 1.0d-6
    integer                                     :: ires, NL_alg
    integer*8                                   :: opt
    
    real(dp)									:: time
    integer 									:: clock0, clock1, clockrate, ticks

	! Initialize clock
	if (time_flag) then
    	call system_clock(clock0, count_rate = clockrate)
	end if

    ! Set algorithm
    if ( time_flag ) then
        NL_alg = alg
    else
        NL_alg = 12
    end if

    ! Build NLopt library optimization options
    ! 12 - Praxis, 25 - BOBYQA, 28 - Nelder Mead, 29 - SBPLX
    opt = 0
    call nlo_create(opt, NL_alg, n)  ! Creates option object with desired algorithm
    call nlo_set_min_objective(ires,opt,sigmamax,A_in) ! Sets objective function to sigmamax
    call nlo_set_ftol_rel(ires, opt, tol_rel)   ! Specifies relative tolerance
    call nlo_set_ftol_abs(ires, opt, tol_abs)   ! Specifies absolute tolerance
    call nlo_set_maxeval(ires,opt,0)
!     call nlo_set_maxtime(ires,opt,0.0d0)

    ! Invoke NLopt
    call nlo_optimize(ires, opt, x, sigma)
    call nlo_destroy(opt)

    print *, 'Min. sigma: ', sigma
    
    ! Compute run time
	if (time_flag) then
    	call system_clock(clock1)
    	time = real(clock1-clock0,dp) / real(clockrate,dp)
    	print *, "Algorithm ", NL_alg, " runtime: ", time
	end if


end subroutine optimize

! =================================================================================================

! Computes || inv(F(x)) * G(x) ||_2, the objective function in the bilinear quadrature formulation
! Input variables "grad" and "need_grad" are only present to conform to NLopt subroutine format

subroutine sigmamax(output, n, x, grad, need_grad, A_in)

    implicit none

    real(dp),                   intent(out) :: output
    integer,                    intent(in)  :: n
    real(dp),   dimension(1:n), intent(in)  :: x, grad
    integer,                    intent(in)  :: need_grad
    procedure(eval)                         :: A_in

    real(dp),   dimension(1:numpt,1:numg)   :: eval_out
    real(dp),   dimension(1:mink)           :: sigma
    real(dp),   dimension(1:LWORK)          :: WORK
    integer                                 :: INFO


    call A_in(x,numpt,eval_out)

    ! Compute singular values with LAPACK
    call dgesvd('N','N',numpt,numg,eval_out,numpt,sigma,0.0d0,1,0.0d0,1,WORK,LWORK,INFO)

    output = sigma(1)

end subroutine sigmamax

! =================================================================================================

! Slightly perturbs an input point by a small random value. Used to force optimization procedure to
! escape from local but nonglobal minima.

subroutine ptperturb(x,eps)

	implicit none
	
	real(dp),	dimension(:),	intent(inout)	:: x
	real(dp),					intent(in)		:: eps
	
	integer										:: loop, seed
	
	call system_clock(count = seed)
    call srand(seed)
    
    do loop = 1, size(x,1)
    	x(loop) = x(loop) + eps * (2*rand() - 1)
    end do
	

end subroutine ptperturb


! =================================================================================================


end module bquad
