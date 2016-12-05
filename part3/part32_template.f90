!-------------------------------
module syncmodule
	implicit none
	complex(kind=16), parameter :: ii=cmplx(0.0,1.0) !ii = sqrt(-1)
    integer :: ntotal, a !total number of oscillators, 
	real(kind=8) :: c,mu,sigma !coupling coefficient, mean, std
	save
end module syncmodule
!-------------------------------

program sync_mpi
    use mpi
    use syncmodule
    implicit none
    integer :: i1,j1
    integer :: nt !number of time steps, number of oscillators
    real(kind=8) :: dt,pi !time step
    integer :: myid, numprocs, ierr
    real(kind=8), allocatable, dimension(:) :: f0,w,f ! initial condition, frequencies, solution
    real(kind=8), allocatable, dimension(:) :: order !order parameter

 ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

!gather input
    open(unit=10,file='data.in')
        read(10,*) ntotal
        read(10,*) nt
        read(10,*) dt
        read(10,*) c
        read(10,*) sigma
        read(10,*) a
    close(10)

    allocate(f0(ntotal),f(ntotal),w(ntotal),order(nt))
	

!generate initial condition
    pi = acos(-1.d0)
    call random_number(f0)
    f0 = f0*2.d0*pi


!generate frequencies
    mu = 1.d0       
    call random_normal(ntotal,w)
    w = sigma*w+mu    
    
!compute solution
    call euler_mpi(MPI_COMM_WORLD,numprocs,ntotal,0.d0,f0,w,dt,nt,f,order)


!output solution (after completion of gather in euler_mpi)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       if (myid==0) then
        open(unit=11,file='theta.dat')
        do i1=1,ntotal
            write(11,*) f(i1)
        end do
        close(11)
        
        open(unit=12,file='order.dat')
        do i1=1,nt
	    write(12,*) order(i1)
	end do
	close(12)
    end if
    !can be loaded in python with: f=np.loadtxt('theta.dat')
   
    call MPI_FINALIZE(ierr)
end program sync_mpi


subroutine output



subroutine euler_mpi(comm,numprocs,n,t0,y0,w,dt,nt,y,order)
    !explicit Euler method, parallelized with mpi
    !input: 
    !comm: MPI communicator
    !numprocs: total number of processes
    !n: number of oscillators
    !t0: initial time
    !y0: initial phases of oscillators
    !w: array of frequencies, omega_i
    !dt: time step
    !nt: number of time steps
    !output: y, final solution
    !order: order at each time step
    use mpi
    use syncmodule
    implicit none
    integer, intent (in) :: n,nt
    real(kind=8), dimension(n), intent(in) :: y0,w
    real(kind=8), intent(in) :: t0,dt,order(nt)
    real(kind=8), dimension(n), intent(out) :: y
    real(kind=8) :: t
    integer :: i1,k,istart,iend
    integer :: comm,myid,ierr,numprocs
  

    call MPI_COMM_RANK(comm, myid, ierr)
    print *, 'start euler_mpi, myid=',myid

    !set initial conditions
    y = y0
    t = t0

    !generate decomposition and allocate sub-domain variables
    call mpe_decomp1d(size(y),numprocs,myid,istart,iend)
    print *, 'istart,iend,threadID=',istart,iend,myid

    
 
 
    !time marching
    do k = 1,nt
       

        call RHS_mpi(!add code here)

        ylocal= ylocal + dt*Rpart !ylocal must be declared and defined, Rpart must be declared, and 
                                  !should be returned by RHS_mpi


	!compute order, and store on myid==0
    
    end do
 
 
    print *, 'before collection',myid, maxval(abs(ylocal))
  

    call MPI_GATHER(!add code here &
								,MPI_INT,0,comm,ierr)
      !collect ylocal from each processor onto myid=0

    if (myid==0) then
        !compute disps
    end if

    !collect ylocal from each processor onto myid=0
	call MPI_GATHERV(ylocal,!add code here &
                disps,MPI_DOUBLE_PRECISION,0,comm,ierr)

    if (myid==0) print *, 'finished',maxval(abs(y))


end subroutine euler_mpi
!-------------------------
subroutine RHS_mpi(nn,t,w,f,rhs)
    !called by euler_mpi
    !rhs = df/dt
    use syncmodule
    implicit none
    integer, intent(in) :: nn
    real(kind=8), intent(in) :: t
!dimensions of variables below must be added    
    real(kind=8), dimension( ), intent(in) :: w 
    real(kind=8), dimension( ), intent(in) :: f
    real(kind=8), dimension( ), intent(out) :: rhs


!Add code to compute rhs


end subroutine RHS_mpi


!--------------------------------------------------------------------
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in online MPE documentation.
!  This file contains a routine for producing a decomposition of a 1-d array
!  when given a number of processors.  It may be used in "direct" product
!  decomposition.  The values returned assume a "global" domain in [1:n]
!
subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
    implicit none
    integer :: n, numprocs, myid, s, e
    integer :: nlocal
    integer :: deficit

    nlocal  = n / numprocs
    s       = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s       = s + min(myid,deficit)
    if (myid .lt. deficit) then
        nlocal = nlocal + 1
    endif
    e = s + nlocal - 1
    if (e .gt. n .or. myid .eq. numprocs-1) e = n

end subroutine MPE_DECOMP1D

!--------------------------------------------------------------------

subroutine random_normal(n,rn)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

IMPLICIT NONE
integer, intent(in) :: n
real(kind=8), intent(out) :: rn(n)
!     Local variables
integer :: i1
REAL(kind=8)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region
do i1=1,n

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156d0 * (v - 0.5d0)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.d0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
rn(i1) = v/u
end do
RETURN


END subroutine random_normal

