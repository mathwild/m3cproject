!Part 2 by Mathilde Duverger CID: 00978498

module flunet
        use omp_lib
	implicit none
        !add variables as needed
	save
	contains

subroutine rhs(P,n,y,t,a,b0,b1,g,k,w,dy)
    implicit none
    !Return RHS of network flu model
    !input: 
    !n: total number of nodes
    !y: S,E,C
    !t: time
    !a,b0,b1,g,k,w: model parameters
    !output: dy
    !I have added P the transport matrix as an input because
    !it is used when computing dy/dt
    integer, intent(in) :: n
    real(kind=8), dimension(n*3),intent(in) :: y
    real(kind=8), intent(in) :: t,a,b0,b1,g,k,w
    real(kind=8), dimension(n,n), intent(in) :: P
    real(kind=8), dimension(n*3), intent(out) :: dy
    real(kind=8), dimension(n) :: S,E,C,Prow
    real(kind=8) :: b,Pi
    integer :: i
    Pi = acos(-1.d0)
    S = y(1:n)
    E = y(n+1:2*n)
    C = y(2*n+1:3*n)
    b = b0 + b1*(1+cos(2*Pi*t))
    do i=1,n
        Prow = P(i,:)
        dy(i)= k*(1-S(i))-b*C(i)*S(i)+w*dot_product(Prow,S)-w*S(i)
        dy(n+i)= b*C(i)*S(i)-(k+a)*E(i)+w*dot_product(Prow,E)-w*E(i)
        dy(2*n+i)= a*E(i)-(g+k)*C(i)+w*dot_product(Prow,C)-w*C(i)
    end do

end subroutine rhs



subroutine rhs_omp(P,n,y,t,a,b0,b1,g,k,w,numthreads,dy)
    implicit none
    !Return RHS of network flu model, parallelized with OpenMP
    !input: 
    !n: total number of nodes
    !y: S,E,C
    !t: time
    !a,b0,b1,g,k,w: model parameters
    !numthreads: the parallel regions should use numthreads threads
    !output: dy, RHS
    !I have added P the transport matrix as an input because
    !it is used when computing dy/dt
    integer, intent(in) :: n, numthreads
    real(kind=8), dimension(n*3),intent(in) :: y
    real(kind=8), intent(in) :: t,a,b0,b1,g,k,w
    real(kind=8), dimension(n,n), intent(in) :: P
    real(kind=8), dimension(n*3), intent(out) :: dy
    real(kind=8), dimension(n) :: S,E,C
    real(kind=8) :: b,Pi 
    integer :: i
    Pi = acos(-1.d0)
    S = y(1:n)
    E = y(n+1:2*n)
    C = y(2*n+1:3*n)
    b = b0 + b1*(1+cos(2*Pi*t))
    !$call omp_set_num_threads(numthreads)
    !$OMP parallel do
    do i=1,n
        dy(i)= k*(1-S(i))-b*C(i)*S(i)+w*dot_product(P(i,:),S)-w*S(i)
        dy(n+i)= b*C(i)*S(i)-(k+a)*E(i)+w*dot_product(P(i,:),E)-w*E(i)
        dy(2*n+i)= a*E(i)-(g+k)*C(i)+w*dot_product(P(i,:),C)-w*C(i)
    end do
    !$OMP end parallel do
end subroutine rhs_omp



end module flunet
