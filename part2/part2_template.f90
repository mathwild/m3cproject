module flunet
        use omp_lib
	implicit none
        !add variables as needed
	save
	contains


subroutine rhs(n,y,t,a,b0,b1,g,k,w,dy)
    implicit none
    !Return RHS of network flu model
    !input: 
    !n: total number of nodes
    !y: S,E,C
    !t: time
    !a,b0,b1,g,k,w: model parameters
    !output: dy, RHS
    integer, intent(in) :: n
    real(kind=8), dimension(n*3),intent(in) :: y
    real(kind=8), intent(in) :: t,a,b0,b1,g,k,w
    real(kind=8), dimension(n*3), intent(out) :: dy

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
