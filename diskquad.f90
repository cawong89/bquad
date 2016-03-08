program diskquad

    use zernike
    use bquad

    implicit none

    integer,    parameter                       :: dp = kind(1.0d0)
    integer,    parameter                       :: pi = 4*atan(1.0d0)
    integer										:: seed, numrand
    integer                                     :: maxdeg, numpt, numg, dm
    real(dp),   dimension(:),   allocatable     :: x0, y0
    real(dp),   dimension(:,:), allocatable     :: W
    real(dp)                                    :: sigma_out

    integer                                     :: k
    real(dp)                                    :: r, theta
    logical                                     :: input_flag

    real(dp)									:: time

!   include 'nlopt.f'

! ================================================================================================

    input_flag = .true.
    dm = 2										! Dimension of space
    maxdeg = 8									! Bquad evaluates IP for degrees up to maxdeg - 1
    numrand = 1000								! Number of random initial pts to generate
    numpt = maxdeg * (maxdeg + 1)/2				! Number of eval pts in quadrature
    numg = maxdeg + 1							! Dimension of minimizing space

    allocate(x0(1:numpt *dm), W(1:numpt,1:numpt), y0(1:numpt *dm))
   
	! Build parameters

    call zern_buildpars(maxdeg)
    call bquad_buildpars(numpt,dm,numg)

    ! Read input points
    if (input_flag) then
        open(unit = 1, file = "input.txt", action = "read", status = "old")

        do k = 1,numpt
            read(1,*) x0(2*k-1), x0(2*k)
        end do

    else
    	call system_clock(count = seed)
        call srand(seed)
        do k = 1,numpt
            r = sqrt(rand())
            theta = 2*pi*rand()
            x0(2*k - 1) = r*cos(theta)
            x0(2*k) = r*sin(theta)
            !print *, 'Error: Input flag set to false.'
            !stop
        end do

    end if


    ! Construct quadrature

	call bquad_construct(dm*numpt,zern_bquad,x0,W,sigma_out,time)


    ! Write output points

    open(unit = 2, file = "output.txt", action = "write", status = "replace")

    do k = 1, numpt
        write(2,*) x0(2*k -1), x0(2*k)
    end do


end program diskquad
