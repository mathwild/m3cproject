!gfortran network.f90 part1_dev.f90 -o test.exe test.f90 -llapack

program test
    use rwmodule
    implicit none
    integer, parameter :: N0 = 3
    integer, parameter :: L = 1
    integer, parameter :: Nt = 2
    integer :: Ntime,m,X0,isample
    integer, dimension(11,4) :: X
    real(kind=8), dimension(10) :: XM
    real(kind=8) :: walltime
    integer(kind=8) :: start,finish,clockrate
    walltime = 0.d0
    start = 0.d0
    finish = 0.d0
    Ntime = 10
    m = 4
    X0 = 0
    isample = 1
    call rwnet(Ntime,m,X0,N0,L,Nt,isample,X,XM)
    print *, X
    !call system_clock(start)
    !call rwnet(Ntime,m,X0,N0,L,Nt,isample,X,XM)
    !call system_clock(finish,clockrate)
    !walltime = dble(finish-start)/dble(clockrate)
    !print *, 'w1', walltime
    !call system_clock(start)
    !call rwnet_omp(Ntime,m,X0,N0,L,Nt,isample,2,X,XM)
    !call system_clock(finish,clockrate)
    !walltime = dble(finish-start)/dble(clockrate)
    !print *, 'w2', walltime
end program