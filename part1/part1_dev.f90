!Project part 1
module rwmodule
    use network
!f2py -llapack -c part1_dev.f90 network.f90 --f90flags='-fopenmp' -lgomp -m rw
contains



subroutine rwnet(Ntime,m,X0,N0,L,Nt,isample,X,XM)
    !random walks on a recursive network
    !Input variables:
    !Ntime: number of time steps
    !m: number of walks
    !X0: initial node, (node with maximum degree if X0=0)
    !N0,L,Nt: recursive network parameters
    !Output: X, m random walks on network each with Ntime steps from initial node, X0
    !XM: fraction of walks at initial node, computed every isample timesteps
    implicit none
    integer, intent(in) :: Ntime,m,N0,L,Nt,X0,isample
    integer, dimension(Ntime+1,m), intent(out) :: X
    real(kind=8), dimension(Ntime/(isample)), intent(out) :: XM
    integer :: qmax
    integer, dimension(N0+Nt) :: qnet
    integer, dimension(N0+L*Nt,2) :: enet
    integer, dimension(2*(N0+L*Nt)) :: alist1
    integer, dimension(N0+Nt+1) :: alist2
    real(kind=8), dimension(N0+Nt,N0+Nt) :: S
    integer :: i,j,k,mcount,ncount,icount,wcount
    real(kind=8) :: y
    integer, dimension(1) :: u
    call generate(N0,L,Nt,qmax,qnet,enet)
    call adjacency_list(qnet,enet,alist1,alist2)
    !compute the stochastic matrix S
    do i=1,N0+Nt
        do j =1,N0+Nt
            ncount = 0
            do k =alist2(i),alist2(i+1)-1
                if (alist1(k)==j) then
                    ncount = ncount + 1
                end if
            end do
            S(i,j) = ncount/dble(qnet(i))
        end do 
    end do
    !transform each row of the stochastic matrix into
    !cumulative probabilities
    do i=1,N0+Nt
        do j=2,N0+Nt
            S(i,j) = S(i,j-1) + S(i,j)
        end do
    end do
    do mcount=1,m
        !set the initial node for each walker
        if (X0==0) then
            u = maxloc(qnet)
            X(1,mcount) = u(1)
        else
            X(1,mcount) = X0
        end if 
        !compute the node we move to at each timestep
        !to do so, we call a random number from a uniform 
        !distribution U(0,1). we then look at the row corresponding
        !to our present node and look at the cdf to assign the random 
        !number to one of the nodes.
        do i=2,Ntime+1
            call random_number(y)
            do j = 1,N0+Nt
                if (y<=S(X(i-1,mcount),j)) then
                    X(i,mcount) = j
                    exit  
                end if
            end do
        end do
    end do
    !we compute Xm, the fraction of walkers that are
    !at the initial node at each timestep
    icount = 1+isample
    i = 1
    do while (icount<Ntime+2) 
        wcount = 0
        do j=1,m
            if (X(icount,j)==X(1,j)) then
                wcount = wcount+1
            end if
        end do
        XM(i) = wcount/dble(m)
        icount = icount + isample
        i = i +1
    end do
    
end subroutine rwnet


subroutine rwnet_omp(Ntime,m,X0,N0,L,Nt,isample,numthreads,X,XM)
    !parallelized version of rwnet, parallel regions should
    !use numthreads threads
    implicit none
    integer, intent(in) :: Ntime,m,N0,L,Nt,X0,isample,numthreads
    integer, dimension(Ntime+1,m), intent(out) :: X
    real(kind=8), dimension(Ntime/(isample)), intent(out) :: XM
    integer :: qmax
    integer, dimension(N0+Nt) :: qnet
    integer, dimension(N0+L*Nt,2) :: enet
    integer, dimension(2*(N0+L*Nt)) :: alist1
    integer, dimension(N0+Nt+1) :: alist2
    real(kind=8), dimension(N0+Nt,N0+Nt) :: S
    integer :: i,j,k,mcount,ncount,icount,wcount
    real(kind=8) :: y
    integer, dimension(1) :: u
    call generate(N0,L,Nt,qmax,qnet,enet)
    call adjacency_list(qnet,enet,alist1,alist2)
    !$call omp_set_num_threads(numthreads)
    !$OMP parallel do
    !compute the stochastic matrix S
    do i=1,N0+Nt
        do j =1,N0+Nt
            ncount = 0
            do k =alist2(i),alist2(i+1)-1
                if (alist1(k)==j) then
                    ncount = ncount + 1
                end if
            end do
            S(i,j) = ncount/dble(qnet(i))
        end do 
    end do
    !$OMP end parallel do
    !$OMP parallel do
    !transform each row of the stochastic matrix into
    !cumulative probabilities
    do j=2,N0+Nt
        do i=1,N0+Nt
            S(i,j) = S(i,j-1) + S(i,j)
        end do
    end do
    !$OMP end parallel do
    !$OMP parallel do
    do mcount=1,m
        !set the initial node for each walker
        if (X0==0) then
            u = maxloc(qnet)
            X(1,mcount) = u(1)
        else
            X(1,mcount) = X0
        end if 
        !compute the node we move to at each timestep
        !to do so, we call a random number from a uniform 
        !distribution U(0,1). we then look at the row corresponding
        !to our present node and look at the cdf to assign the random 
        !number to one of the nodes.
        do i=2,Ntime+1
            call random_number(y)
            do j = 1,N0+Nt
                if (y<=S(X(i-1,mcount),j)) then
                    X(i,mcount) = j
                    exit  
                end if
            end do
        end do
    end do
    !$OMP end parallel do
    !we compute Xm, the fraction of walkers that are
    !at the initial node at each timestep
    icount = 1+isample
    i = 1
    do while (icount<Ntime+2) 
        wcount = 0
        do j=1,m
            if (X(icount,j)==X(1,j)) then
                wcount = wcount+1
            end if
        end do
        XM(i) = wcount/dble(m)
        icount = icount + isample
        i = i +1
    end do

end subroutine rwnet_omp


end module rwmodule
