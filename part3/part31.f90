module sync
	implicit none
	complex(kind=8), parameter :: ii=cmplx(0.0,1.0) !ii = sqrt(-1)
	integer :: ntotal !total number of nodes
	real(kind=8) :: c !coupling coefficient
        real(kind=8), allocatable, dimension(:) :: w !array of frequencies
	save
    contains

!Project part 3 by Mathilde Duverger CID: 00978498

!---------------------------
subroutine rk4(t0,y0,dt,nt,y,order)
    !4th order RK method
    !input:
    !t0: initial time
    !y0: initial condition (array)
    !dt: time step
    !nt: number of time steps
    !output:
    !y: solution at each time step (array)
    !order: order parameter
    implicit none
    real(kind=8), dimension(:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
    integer, intent (in) :: nt
    real(kind=8), dimension(size(y0),nt+1), intent(out) :: y
    real(kind=8), dimension(nt), intent(out) :: order
    real(kind=8), dimension(size(y0)) :: f1, f2, f3, f4
    real(kind=8) :: t,halfdt,fac
    integer:: k
    
    halfdt = 0.5d0*dt
    fac = 1.d0/6.d0

    y(:,1) = y0 !initial condition
    t = t0 !initial time
    
    do k = 1, nt !advance nt time steps
    
        f1 = dt*RHS(t, y(:,k))

        f2 = dt*RHS(t + halfdt, y(:,k) + 0.5d0*f1)

        f3 = dt*RHS(t + halfdt, y(:,k) + 0.5d0*f2)

        f4 = dt*RHS(t + dt, y(:,k) + f3)
        
        y(:,k+1) = y(:,k) + (f1 + 2*f2  + 2*f3 + f4)*fac

        t = t + dt
        
        !add code to compute order
        order(k) = abs(sum(exp(ii*y(:,k+1))))/dble(ntotal)
        
    end do
        
        
end subroutine rk4

!---------------------------
function RHS(t,f)
    !RHS sync
    !f is the array of phases at time, t
    implicit none
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(:), intent(in) :: f
    real(kind=8), dimension(size(f)) :: RHS
    integer :: i

    !Add code to compute RHS
    do i=1,ntotal
        RHS(i) = w(i) - c*sum(sin(f(i)-f))/dble(ntotal)
    end do
 
end function RHS
!---------------------------
end module sync


