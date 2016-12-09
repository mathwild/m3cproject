module flunet
        use omp_lib
	implicit none
        !add variables as needed
	save
	contains


subroutine rhs(y,t,n,a,b0,b1,g,k,w,dy)
    implicit none
    !Return RHS of network flu model
    !input: 
    !n: total number of nodes
    !y: S,E,C
    !t: time
    !a,b0,b1,g,k,w: model parameters
    !output: dy
    integer, intent(in) :: n
    real(kind=8), dimension(n*3),intent(in) :: y
    real(kind=8), intent(in) :: t,a,b0,b1,g,k,w
    real(kind=8), dimension(n*3), intent(out) :: dy
    real(kind=8), dimension(n) :: S,E,C
    S = y(1:n)
    E = y(n:2*n)
    C = y(2*n:3*n)
    b = b0 + b1*(1+cos(2*Pi*t))
    do i=1,n
        dy(i)= k*(1-S(i))-b*C(i)*S(i)+w*dot_product(P(i,:),S)-w*S(i)
        dy(n+i)= b*C(i)*S(i)-(k+a)*E(i)+w*dot_product(P(i,:),E)-w*E(i)
        dy(2*n+i)= a*E(i)-(g+k)*C(i)+w*dot_product(P(i,:),C)-w*C(i)
    end do

end subroutine rhs



subroutine rhs_omp(n,y,t,a,b0,b1,g,k,w,numthreads,dy)
    implicit none
    !Return RHS of network flu model, parallelized with OpenMP
    !input: 
    !n: total number of nodes
    !y: S,E,C
    !t: time
    !a,b0,b1,g,k,w: model parameters
    !numthreads: the parallel regions should use numthreads threads
    !output: dy, RHS
    integer, intent(in) :: n
    real(kind=8), dimension(n*3),intent(in) :: y
    real(kind=8), intent(in) :: t,a,b0,b1,g,k,w
    real(kind=8), dimension(n*3), intent(out) :: dy

end subroutine rhs_omp



end subroutine rhs_omp
