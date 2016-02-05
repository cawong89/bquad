program diskquad

    use zernike
    use bquad

    implicit none

    integer,    parameter                       :: dp = kind(1.0d0)
    integer,    parameter                       :: pi = 4*atan(1.0d0)
    integer										:: seed
    integer                                     :: maxdeg, numpt, numg, dm
    real(dp),   dimension(:),   allocatable     :: x0, y0
    real(dp),   dimension(:,:), allocatable     :: W
    real(dp)                                    :: sigma_out

    real(dp),   dimension(:,:), allocatable     :: A_debug
    integer                                     :: k
    real(dp)                                    :: r, theta
    logical                                     :: input_flag

    integer										:: j,jj
    integer,	dimension(1:4)					:: alg_list

!   include 'nlopt.f'

! ================================================================================================

    input_flag = .false.
    dm = 2
    maxdeg = 3
    numpt = maxdeg * (maxdeg + 1)/2
    numg = maxdeg + 1

    allocate(x0(1:numpt *dm), W(1:numpt,1:numpt), y0(1:numpt *dm))
    allocate(A_debug(1:numpt, 1:numg))

    alg_list = (/ 28, 29, 25, 12 /)

    ! Read input points
    if (input_flag .eqv. .true.) then
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

    ! Build parameters

    call zern_buildpars(maxdeg)
    call bquad_buildpars(numpt,dm,numg)

    ! Construct quadrature

!     call zern_bquad(x0,numpt,A_debug)
!     print *, A_debug

	do j = 1,3
		do jj = 1,4
			call bquad_construct(dm*numpt,zern_bquad,x0,W,sigma_out, alg = alg_list(jj))
		end do
	end do

    ! Write output points

    open(unit = 2, file = "output.txt", action = "write", status = "replace")

    do k = 1, numpt
        write(2,*) x0(2*k -1), x0(2*k)
    end do


end program diskquad
