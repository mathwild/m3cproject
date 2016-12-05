program test
    use network
    implicit none
    integer, parameter :: N0 = 3
    integer, parameter :: L = 1
    integer, parameter :: Nt = 4
    integer :: qmax
    integer, dimension(N0+Nt) :: qnet
    integer, dimension(N0+L*Nt,2) :: enet 
    integer, dimension(2*(N0+L*Nt)) :: alist1
    integer, dimension(N0+Nt) :: alist2
    integer, dimension(1) :: start
    call generate(N0,L,Nt,qmax,qnet,enet)
    print *, qnet
    call adjacency_list(qnet,enet,alist1,alist2)
    start = maxloc(qnet)
    print *, start
end program