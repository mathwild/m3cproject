program test
    use rwmodule
    implicit none
    integer, parameter :: N0 = 3
    integer, parameter :: L = 1
    integer, parameter :: Nt = 3
    integer :: Ntime,m,X0,isample
    integer, dimension(4,2) :: X
    integer, dimension(3) :: XM
    Ntime = 3
    m = 2
    X0 = 0
    isample = 1
    call rwnet(Ntime,m,X0,N0,L,Nt,isample,X,XM)
    print *, X
end program