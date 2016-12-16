!Project part 3 by Mathilde Duverger CID: 00978498

!-------------------------------
module syncmodule
    implicit none
    complex(kind=8), parameter :: ii=cmplx(0.0,1.0) !ii = sqrt(-1)
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
    integer, intent (in) :: comm,numprocs,n,nt
    real(kind=8), dimension(n), intent(in) :: y0,w
    real(kind=8), intent(in) :: t0,dt
    real(kind=8), dimension(n), intent(out) :: y
    real(kind=8), dimension(nt), intent(out) :: order
    real(kind=8), allocatable, dimension(:) :: ylocal
    real(kind=8), allocatable, dimension(:) :: f
    real(kind=8), allocatable, dimension(:) :: Rpart
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, allocatable, dimension(:) :: Nper_proc, disps
    real(kind=8) :: t
    complex(kind=8), dimension(nt) :: R
    complex(kind=8), dimension(nt) :: Rtot
    integer :: i1,k,istart,iend,nn,receiver,sender
    integer :: myid,ierr
  
    call MPI_COMM_RANK(comm, myid, ierr)
    print *, 'start euler_mpi, myid=',myid

    !set initial conditions
    y = y0
    t = t0

    !generate decomposition and allocate sub-domain variables
    call mpe_decomp1d(size(y),numprocs,myid,istart,iend)
    print *, 'istart,iend,threadID=',istart,iend,myid
    
    !allocate and define variables
    allocate(ylocal(iend-istart+1))
    ylocal = y0(istart:iend)
    nn = iend-istart+2*a+1
    allocate(f(nn))
    allocate(RPart(nn-2*a))
    
    !time marching
    do k = 1,nt
        if (myid<numprocs-1) then
            receiver = myid+1
        else
            receiver = 0
        end if
    
        if (myid>0) then
            sender = myid-1
        else 
            sender = numprocs-1
        end if
        
        !send and receive the previous a oscillator's thetas
        call MPI_SEND(ylocal(size(ylocal)-a+1:size(ylocal)),a,MPI_DOUBLE_PRECISION,receiver,0,comm,ierr)
        call MPI_RECV(f(1:a),a,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,comm,status,ierr)
        
        if (myid>0) then
            receiver = myid-1
        else
            receiver = numprocs-1
        end if
        
        if (myid<numprocs-1) then
            sender = myid+1
        else
            sender=0
        end if
        
        !send and receive the next a oscillator's thetas
        call MPI_SEND(ylocal(1:a),a,MPI_DOUBLE_PRECISION,receiver,0,comm,ierr)
        call MPI_RECV(f(nn-a+1:nn),a,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,comm,status,ierr)
      
        f(a+1:nn-a)= ylocal  
            
        call RHS_mpi(nn,t,w(istart:iend),f,Rpart)

        ylocal= ylocal + dt*Rpart !ylocal must be declared and defined, Rpart must be declared, and 
                                  !should be returned by RHS_mpi

        t = t+dt
        
	!compute partial order R
        R(k) = sum(exp(ii*ylocal))
        
        !collect R from each processor onto myid=0 and compute order
        call MPI_REDUCE(R(k),Rtot(k),1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
        if (myid==0) then 
            order = abs(Rtot)/dble(n)
        end if
        
    end do
 
    print *, 'before collection',myid, maxval(abs(ylocal))
    
    allocate(Nper_proc(numprocs),disps(numprocs))
    
    call MPI_GATHER(iend-istart+1,1,MPI_INT,Nper_proc,1,MPI_INT,0,comm,ierr)
      !collect ylocal from each processor onto myid=0

    if (myid==0) then
        !compute disps
        disps(1)=0
        do i1=2,numprocs
            disps(i1)=disps(i1-1) + Nper_proc(i1-1)
        end do
    end if
    
    !collect ylocal from each processor onto myid=0
    call MPI_GATHERV(ylocal,iend-istart+1,MPI_DOUBLE_PRECISION,y,Nper_proc, &
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
    real(kind=8), dimension(nn-2*a), intent(in) :: w 
    real(kind=8), dimension(nn), intent(in) :: f
    real(kind=8), dimension(nn-2*a), intent(out) :: rhs
    real(kind=8) :: sumsin
    integer :: i,j
    
    !Add code to compute rhs
    do i=1,nn-2*a
        rhs(i) = w(i) - c*sum(sin(f(i+a)-f(i:i+2*a)))/dble(ntotal)
    end do

end subroutine RHS_mpi

!--------------------------------------------------------------------

!Part 3 question 6
!Consider the weakly-coupled model on a recursive network

 
!The OMP approach might be prefered for ease of setting up since 
!it uses shared memory, so we wouldn't have to think about 
!the communication between processors. This would obviously create an 
!issue for large problem sizes where all the information could
!not be stored on a unique shared memory.
!Suppose we have a small problem, we could use OMP do around 
!the loop inside RHS quite easily but it would not provide that 
!much of a speedup. If we want to use OMP do around the main time 
!marching loop  it would be much more complicated because information  
!about the previous loop is needed for each marching time. 
!In this case, we might need to use OMP barriers, which
!would waste time since a processor would need to wait 
!for the other one to have finished.
!Using OMP, we wouldn't need to choose a partition since
!OMP runs regardless of partitioning. 

!Given this problem, we would probably prefer using MPI
!since it will allows us to subdivise the problem into
!subparts that we can compute without having to wait for 
!the other processors to finish. However, using MPI
!implies needing to exchange information between processors.
!This involves partitioning the network in such a way that
!the minimum communication is needed. Moreover, the set up
!would be quite complicated as for efficiency we would want
!to call information only about the thetas for which Aij is one.
!In terms of partition, we would need to choose a partition 
!such that there is minimal communication between processors.
!Moreover, we might want to take a look at the degree of each node
!in the processors. Indeed, the highest the degree of a node, the 
!more information we need since we collect the thetas for all
!nodes connected to our node. We might want to strategically
!partition the nodes in such a way that the sum of degrees is
!approximately the same in each processor so there is not processor
!that is done before the others. MPI also allows working
!with bigger problems since it uses distributed memory.



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

