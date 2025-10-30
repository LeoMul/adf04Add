program main 
    implicit none 
    character*300 :: line 
    character*29  :: line2 
    real*8 :: energy
    integer :: lev ,ii , jj,numread 
    integer :: upp_old,low_old 
    integer :: upp_new,low_new
    integer,allocatable :: map(:)
    real*8,allocatable :: shifted(:)
    character*300,allocatable :: leveldata(:)
    integer :: numlev 
    integer :: io 
!
    open(1,file='adf04_ups',action='read')
!
    open(2,file='adf04_ups_shuf',action='write')
!
    open(3,file='indexShuffle',action='read')!
    read(3,*) numread 
    allocate(map(numread))
    do lev = 1,numread 
        read(3,*) ii ,jj 
        map(ii) = jj 
    end do 
    close(3)
!
    open(4,file='shiftedEnergies',action='read')
    allocate(shifted(numread))
    do lev = 1,numread
        read(4,*) shifted(lev)
    end do 
    close(4)
!
    read (1,'(A300)') line 
    write(2,'(A300)') line
    !read in levels 
!
    allocate(leveldata(numread))
    lev = 0 
    numlev = 1
    read (1,'(I5,A29,F20.4)') lev,line2,energy
    do while (lev .ne. -1 )
        leveldata(lev) = line2
        numlev = lev 
        read (1,'(I5,A29,F20.4)',iostat=io) lev,line2,energy
    end do 
!
    do lev = 1, numlev
        write(2,'(I5,A29,F20.4)') lev,leveldata(map(lev)),shifted(lev)
    end do 
    write(2,'(I5,A300)') -1,line 
!
    read(1,'(A300)') line 
    write(2,'(A300)') line
!
    read (1,'(2I5,A300)') upp_old,low_old, line
    do while (upp_old .ne. -1 )
        upp_new = map(upp_old)
        low_new = map(low_old) 
        write(2,'(2I5,A300)') max(upp_new,low_new), min(upp_new,low_new)&
        & ,line 
        read (1,'(2I5,A300)') upp_old,low_old, line
    end do
    write(2,'(I5)')  -1
!
!
    read(1,'(a300)',iostat=io) line
    do while (.not. IS_IOSTAT_END(io))
        write(2,'(a300)') line 
        read(1,'(a300)',iostat=io) line
    end do 
    write(2,*) 'C Shuffled by the adf04Shuffle.f90 of https://github.com/LeoMul'

    close(1)
    close(2)
end program 
