!Recursive network model
module network

	contains


subroutine generate(N0,L,Nt,qmax,qnet,enet)
	!Generate recursive network corresponding
	!to model parameters provided as input
	!Output: max degree, degree list, edge list
	implicit none
	integer, intent(in) :: N0,L,Nt
	integer, intent(out) :: qmax
	integer, dimension(N0+Nt), intent(out) :: qnet
	integer, dimension(N0+L*Nt,2), intent(out) :: enet 
	integer :: i1,j1,k1,N,M,i_add,ntotal
	real(kind=8), allocatable, dimension(:) :: CP,P
	real(kind=8), allocatable, dimension(:,:) :: R
		
	
	!Initialize network, allocate variables
	ntotal = N0 + Nt
	
	qnet(1:N0) = 2
	do i1 = 1,N0-1
		enet(i1,:) = (/i1,i1+1/)
	end do
	enet(N0,:) = (/N0,1/)
	
	
	!generate random numbers
	allocate(R(Nt,L))

	call random_number(R)
	allocate(CP(ntotal),P(ntotal)) 
	
	N = N0 !total (initial) nodes
	M = N0 !total (initial) links
	!Iterate over Nt steps
	do i1=1,Nt

		!convert degrees into probabilities
		P(1:N) = dble(qnet(1:N))/sum(qnet(1:N))
		!Convert P into cumulative probability
		CP(1) = P(1)
		CP(N) = 1.d0
		do j1=2,N-1
			CP(j1) = CP(j1-1) + P(j1)
		end do		
	
		!assign links
		do j1=1,L
			i_add = minloc(abs(CP(1:N)-R(i1,j1)),1) !locate R in CP
			if (R(i1,j1)>CP(i_add)) i_add = i_add + 1 

			qnet(i_add) = qnet(i_add)+1
			enet(M+j1,:) = (/N+1,i_add/)
		end do
		qnet(N+1) = L	
		N = N+1 !update number of nodes
		M = M + L !update number of links
		
	end do
	!set qmax
	qmax = maxval(qnet)	

	deallocate(R,CP,P)
end subroutine generate

subroutine adjacency_list(qnet,enet,alist1,alist2)
	!Create adjacency lists corresponding
	!to degree and edge lists provided as input
	implicit none
	integer, dimension(:), intent(in) :: qnet
	integer, dimension(:,:), intent(in) :: enet
	integer, dimension(size(enet,1)*2), intent(out) :: alist1
	integer, dimension(size(qnet)), intent(out) :: alist2
	integer :: i,j,n
        alist2(1) = 1
        do i=2,size(alist2)
            alist2(i) = alist2(i-1) + qnet(i-1)
        end do
         n = 1
        do j=1,size(qnet)
            do i=1,size(enet,1)
                if (enet(i,1)==j) then
                    alist1(n) = enet(i,2)
                    n = n+1
                else if (enet(i,2)==j) then
                    alist1(n) = enet(i,1)
                    n = n+1
                end if
            end do
        end do
        
	
end subroutine adjacency_list



subroutine adjacency_matrix(ntotal,enet,anet)
	!Create adjacency matrix corresponding 
	!to edge list provided as input
	implicit none
	integer, intent(in) :: ntotal
	integer, dimension(:,:), intent(in) :: enet
	integer, dimension(ntotal,ntotal), intent(out) :: anet
	integer :: i1

	anet = 0
	do i1=1,size(enet,1)
		anet(enet(i1,1),enet(i1,2)) = anet(enet(i1,1),enet(i1,2)) + 1
	end do
	
	anet = anet + transpose(anet)
	
end subroutine adjacency_matrix

subroutine connectivity(qnet,anet,c)
	!Compute connectivity of pre-computed network
	!corresponding to degree list qnet and adjacency matrix
	!anet, both provided as input
	implicit none
	integer, dimension(:), intent(in) :: qnet
	integer, dimension(:,:), intent(in) ::  anet
	real(kind=8), intent(out) :: c
	integer :: i1,j1
	real(kind=8), dimension(size(anet,1),size(anet,2)) :: M !Laplacian matrix
	integer :: LWORK,INFO !variables for DSYEV
	real(kind=8) :: lmin
	real(kind=8), allocatable, dimension(:) :: W,WORK
		
	!generate Laplacian matrix, M
	M = 0.d0
	do i1=1,size(anet,1)
		M(i1,i1) = qnet(i1)
	end do
	M = M-anet

	!compute eigenvalues of M, select connectivity
	LWORK = 6*size(anet,1)
	allocate(W(size(anet,1)),WORK(LWORK))
	call DSYEV('N','U',size(M,1),M,size(M,1),W,WORK,LWORK,INFO)
	c = W(2)
	if (INFO .ne. 0) print *,'INFO=',INFO
	deallocate(W,WORK)
end subroutine connectivity


subroutine vary_connectivity(N0,Lmax,Nt,carray)
	!Compute connectivity for networks with L varying
	!between 1 and Lmax
	implicit none
	integer, intent(in) :: N0,Lmax,Nt
	real(kind=8), dimension(Lmax),intent(out) :: carray
	integer :: ntotal,ltotal,L,qmax
	integer, allocatable, dimension(:) :: qnet
	integer, allocatable, dimension(:,:) :: enet,anet

	ntotal = N0+Nt
	
	allocate(qnet(ntotal),anet(ntotal,ntotal))
	
	do L = 1,Lmax
		ltotal = N0+L*Nt
		allocate(enet(ltotal,2))
		call generate(N0,L,Nt,qmax,qnet,enet)
		call adjacency_matrix(ntotal,enet,anet)
		call connectivity(qnet,anet,carray(L))
		print *, L,carray(L)
		deallocate(enet)
	end do
	deallocate(qnet,anet)
	
	
end subroutine vary_connectivity


end module network
