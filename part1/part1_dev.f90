!Project part 1
module rwmodule
    use network

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
    integer, dimension(N0), intent(out) :: XM
    integer :: qmax
    integer, dimension(N0+Nt) :: qnet
    integer, dimension(N0+L*Nt,2) :: enet
    integer, dimension(2*(N0+L*Nt)) :: alist1
    integer, dimension(N0+Nt+1) :: alist2
    real(kind=8), dimension(N0+Nt,N0+Nt) :: S
    integer :: i,j,k,mcount,ncount
    real(kind=8) :: y
    integer, dimension(1) :: u
    call generate(N0,L,Nt,qmax,qnet,enet)
    call adjacency_list(qnet,enet,alist1,alist2)
    do mcount=1,m
        if (X0==0) then
            u = maxloc(qnet)
            X(1,mcount) = u(1)
        end if 
        do i=1,size(qnet)
            do j =1,size(qnet)
                ncount = 0
                do k =alist2(i),alist2(i+1)-1
                    if (alist1(k)==j) then
                        ncount = ncount + 1
                    end if
                end do
                S(i,j) = ncount/dble(qnet(i))
            end do 
        end do
        do i=1,size(qnet)
            do j=2,size(qnet)
                S(i,j) = S(i,j-1) + S(i,j)
            end do
        end do
        do i=2,Ntime+1
            do j = 1,size(qnet)
                call random_number(y)
                if (y<=S(X(i-1,mcount),j)) then
                    X(i,mcount) = j
                    exit  
                end if
            end do
        end do
    end do
end subroutine rwnet


subroutine rwnet_omp(Ntime,m,X0,N0,L,Nt,isample,numthreads,X,XM)
    !parallelized version of rwnet, parallel regions should
    !use numthreads threads
	implicit none
    integer, intent(in) :: Ntime,m,N0,L,Nt,X0,isample,numthreads
    integer, dimension(Ntime+1,m), intent(out) :: X
    integer, dimension(N0), intent(out) :: XM



end subroutine rwnet_omp


end module rwmodule
