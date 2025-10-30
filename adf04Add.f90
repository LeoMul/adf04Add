program adf04_add 
    !Leo Patrick Mulholland
    !Queen's University Belfast (2025)
    
    !This program adds two adf04 files. 
    !Reads in all the transitions from file1 (call these t1)
    !Reads in each transition in file2 (call these t2)
    !If t2 corresponds to a transitiion t1 - add them and print out.
    !If not - ignore and move on.

    !This (slightly) accounts for the potentially not exactly matching of the two files.
    !This said, the program assumes there is the same number of levels in each file
    !And at this moment the same temperatures.
    !In principle, there could be a fit or interpolation done to get around this.
    !This would require a large 

    !Warning, this code is full of hacks.
    !While the adf04 format is standard, it does not make use of standard fortran formats.
    !This leads to the (probably incredibly) inefficient string manipulations required here.
    !This said, typically an adf04 is probably not more than 1 GB - so perhaps we can live
    !with a speed decrease here.

    !Additionally - this is a dumb code. 
    !It knows precisely nothing about your atomic stucture, 
    !A-values, shifting, angular identifications or political affiliations.
    !It is at the discretion of the user to correctly add files that correspond to the 
    !same level indices, and importantly the same ion.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Parameters - probably will never need to be changed
    integer,parameter   :: maxIter = 10000
    integer,parameter   :: maxNumTemps = 320
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Local Variables
    integer             :: yy,zz,ii,jj,kk,offset
    integer             :: lower,upper
    integer             :: numLevels,numLevelsSecondFile
    integer             :: maxNumTransitions
    integer             :: currentpointer
    integer             :: numTemps,numTemps2
    real*8, allocatable :: avalues(:)
    real*8, allocatable :: ups1(:,:) , ups2(:),upsinf(:)
    real*8              :: temperatures (maxNumTemps)
    real*8              :: temperatures2(maxNumTemps)
    integer,allocatable :: pointerarray(:)
    integer,allocatable :: pointer1(:)
    integer,allocatable :: pointer2(:)
    integer             :: pp ,iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Strings for manipulation
    character*8         :: dummyUpsilon
    character*9         :: dummyUpsilon2
    character*320       :: templine
    character*16        :: dummy
    character*1         :: dummy2
    !File names.      
    character*8         :: file1 = 'adf04_1'
    character*8         :: file2 = 'adf04_2'
    character*8         :: file3 = 'adf04_3'
    !File existence
    logical :: ex1,ex2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Program start
    !Check for relevant files.
    inquire(file=file1,exist=ex1)
    if (.not.ex1) stop 'file 1 not found.'
    inquire(file=file2,exist=ex2)
    if (.not.ex2) stop 'file 2 not found.'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Write to stream 2.
    open(2,file=file3,form='formatted')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Open file 1.
    open (1,file=file1,form='formatted')
    read (1,'(A320)') templine
    write(2,*) templine
    !first find number of levels.
    do ii=1,maxIter
        read (1,'(I5,A320)') yy,templine
        write(2,'(I5,A320)') yy,templine
        if (yy.eq.-1) exit
    end do 
    numLevels = ii-1
    maxNumTransitions = (numLevels * (numLevels + 1)) / 2
    print*,'Found ',ii-1,' atomic levels in file 1.'    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Allocate
    allocate(pointerarray(maxNumTransitions))
    allocate(avalues(maxNumTransitions))
    allocate(pointer1(maxNumTransitions))
    allocate(pointer2(maxNumTransitions))
    pointer1  = 0 
    pointer2  = 0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Read temps
    read(1,'(A16,A320)') dummy,templine
    call StripSpaces(templine)
    offset = 0
    do ii = 1,maxNumTemps
        dummy2 = templine(ii+offset:ii+offset)
        temperatures(ii) = a7toFloat(templine(ii+offset:ii+offset+6))
        offset = offset +6
        if (dummy2.eq.' ') then 
            exit 
        end if 
    end do
    print*,'I have found ',ii-1, 'temperatures.'
    numTemps = ii-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Allocate
    allocate(ups1(numTemps,maxNumTransitions))
    allocate(ups2(numTemps))
    allocate(upsinf(maxNumTransitions))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Read transitions.
    do ii = 1,maxNumTransitions
        !not sure if this will work with all compilre==
        read(1,'(2I5,A300)') yy,zz,templine 
        if (yy.eq.-1) exit
        lower = min(zz,yy)
        upper = max(zz,yy)
        !pp = upper-1 + numLevels * (lower-1)
        pp = upperTriangleIndexing(lower,upper,numLevels)
        pointer1(pp) = 1
        offset = 0
        avalues(pp) = a8toFloat(templine(1+offset:1+offset+7))
        offset = offset + 7
        do jj = 2,numTemps+1
            ups1(jj-1,pp) = a8toFloat(templine(jj+offset:jj+offset+7))
            offset = offset +7
        end do
        upsinf(pp) = a8toFloat(templine(jj+offset:jj+offset+7))
    end do 
    print*,'Found ',ii-1,' atomic transitions in file 1.'    
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Read second file
    !some repeated code, could be cleaned up later.
    open(1,file=file2)
    read(1,*)
    !first find number of levels.
    do ii=1,maxIter
        read(1,*) yy
        if (yy.eq.-1) exit
    end do 

    numLevelsSecondFile = ii-1
    print*,'Found ',ii-1,' atomic levels in file 2.'    
    if (numLevels.ne.numLevelsSecondFile) then 
        write(0,*) 'number of levels not same'
        stop
    end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Read temps
    read(1,'(A16,A320)') dummy,templine
    write(2,'(A16,A320)')dummy,templine
    call StripSpaces(templine)
    offset = 0
    do ii = 1,maxNumTemps
        dummy2 = templine(ii+offset:ii+offset)
        temperatures2(ii) = a7toFloat(templine(ii+offset:ii+offset+6))
        offset = offset +6
        if (dummy2.eq.' ') then 
            exit 
        end if 
    end do
    print*,'I have found ',ii-1, 'temperatures in file 2.'
    numTemps2 = ii-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Consistency check
    if (numTemps.ne.numTemps2) stop 'ntemps are not same.'
    do ii = 1,numTemps
        if(temperatures(ii).ne.temperatures2(ii)) stop 'temps not same'
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Read transitions, cross check with file 1.
    !If found in file1, add and output
    do ii = 1,maxNumTransitions
        !not sure if this will work with all compilers.
        read(1,'(2I4,A300)') yy,zz,templine 
        !read(1,'(2I4)',iostat=iostat) yy,zz
        if (yy.eq.-1) exit
        lower = min(zz,yy)
        upper = max(zz,yy)
        !pp = upper-1 + numLevels * (lower-1)
        !pp = upper + (lower * (lower + 1))/2 
        pp = upperTriangleIndexing(lower,upper,numLevels)
!
        offset = 0
        if (pointer1(pp) > 0 ) then 
            !if this transition is in the first file, read it and add it.
            offset = offset + 7
!
            do jj = 2,numTemps+1
                ups2(jj-1) = a8toFloat(templine(jj+offset:jj+offset+7))
                offset = offset +7
            end do
! 
            ups1(:,pp)  = ups1(:,pp) + ups2(:)  
        end if 

    end do 
    print*,'uhhh'
    do upper = 2,numLevels 

    do lower = 1,upper-1

!        ii = upper-1 + numLevels * (lower-1)
!       lower tria
!        ii = upper + (lower * (lower + 1))/2 
!        
        ii = upperTriangleIndexing(lower,upper,numLevels)
        print*,upper,lower,'--->',ii
    !do ii = 1 ,maxNumTransitions 
        pp = pointer1(ii)
        if (pp > 0) then 
            !pp = upper-1 + numLevels * (lower-1) 
            !upper = mod(ii,numLevels) + 1 
            !lower= (ii+1 - upper) / numLevels + 1 
            !print*,upper,lower,ii
            write(2,'(2I4)',advance='no') upper,lower
            write(dummyUpsilon,31) avalues(ii)
            write(2,33,advance='no') dummyUpsilon(1:4)
            write(2,34,advance='no') dummyUpsilon(6:8)
            do kk = 1,numtemps
               write(dummyUpsilon,31) ups1(kk,ii)
               write(2,33,advance='no') dummyUpsilon(1:4)
               write(2,34,advance='no') dummyUpsilon(6:8)
            end do 
            write(dummyUpsilon2,32) upsinf(ii)
            write(2,'(A5,A3)') dummyUpsilon2(1:5),dummyUpsilon2(7:9)

        end if 
    !end do 
    end do 
end do 

!        !ignore avalue in second file.
        !avalues(ii) = a8toFloat(templine(1+offset:1+offset+7))
        

        !ignore infinite energy point in this file.
        !upsinf(ii) = a8toFloat(templine(jj+offset:jj+offset+7))

        !do jj = 1,maxNumTransitions 
        !    if (pointerarray(jj) .eq. currentpointer) then 
        !        !We have found the transition, add and output.
        !        ups2(:) = ups2(:) + ups1(:,jj)
        !        !write(2,10) upper,lower, avalues(jj), ups2(:) , upsinf(jj)
        !        write(2,'(2I4)',advance='no') upper,lower
        !        !write(78,32,advance='no') ii,jj, '1.00-30'
        !        write(dummyUpsilon,31) avalues(jj)
        !        write(2,33,advance='no') dummyUpsilon(1:4)
        !        write(2,34,advance='no') dummyUpsilon(6:8)              
        !        !Dumb string manipulation. There must be a better way
        !        !However - it works?
        !        do kk = 1,numtemps
        !           write(dummyUpsilon,31) ups2(kk)
        !           write(2,33,advance='no') dummyUpsilon(1:4)
        !           write(2,34,advance='no') dummyUpsilon(6:8)
        !        end do 
        !        write(dummyUpsilon2,32) upsinf(jj)
        !        write(2,'(A5,A3)') dummyUpsilon2(1:5),dummyUpsilon2(7:9)
        !        !write(2,33,advance='no') dummyUpsilon2(1:6)
        !        !write(2,34)              dummyUpsilon2(7:9)
        !    end if 
        !end do 

    31  FORMAT(ES8.2)
    32  FORMAT(ES9.2)
    33 FORMAT(1X,A4)
    34 FORMAT(A3)

    write(2,'( I4)') -1
    write(2,'(2I4)') -1,-1

    print*,'Found ',ii-1,' atomic transitions in file 2.'    
    10 FORMAT(2I4,ES10.2,14ES10.2,ES10.2)
    close(1)

    11 FORMAT('C Resultant adf04 file, produced by https://github.com/LeoMul/adf04Add.')
    write(2,11)
    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    stop
    contains 


    function a7toFloat(a7string) result(float)
        !converts A.BC+DE 
        ! or      A.BC-DE
        !to A.BC \times 10 ^{DE}
        real*8 :: mantessa 
        integer :: exponent 
        real*8 :: float
        character*7 :: a7string
        character*4 :: a4string
        character*3 :: a3string 
        integer :: ierr

        a4string = a7string(1:4)
        a3string = a7string(5:7)

        read(a4string,1) mantessa 
        read(a3string,2,iostat=ierr) exponent
        if (ierr.ne.0) then 
            write(0,*) a3string
            stop 'error in a7tofloat'
        end if 

        float = mantessa * (10.0d0**exponent)

        1 format(F4.1)
        2 format(I3)
    end function

    function a8toFloat(a8string) result(float)
        !converts A.BC+DE 
        ! or      A.BC-DE
        !to A.BC \times 10 ^{DE}
        !or:
        !converts
        !converts -A.BC+DE 
        ! or      -A.BC-DE
        !to -A.BC \times 10 ^{DE}
        real*8 :: mantessa 
        integer :: exponent 
        real*8 :: float
        character*8 :: a8string
        character*4 :: a4string
        character*3 :: a3string 
        character*1 :: a1string
        integer :: ierr
        
        a1string = a8string(1:1)
        a4string = a8string(2:5)
        a3string = a8string(6:8)
        


        read(a4string,1) mantessa 
        read(a3string,2,iostat=ierr) exponent
        if (ierr.ne.0) then 
            write(0,*) a3string 
            write(0,*) a8string

            stop 'error in a8tofloat'
        end if 

        float = mantessa * (10.0d0**exponent)

        if (a1string .eq. '-') then 
            float = -1.0d0 * float 
        end if 


        1 format(F4.1)
        2 format(I3)
    end function

    subroutine StripSpaces(string)
    !    https://stackoverflow.com/questions/27179549/removing-whitespace-in-string
        character(len=*) :: string
        integer :: stringLen 
        integer :: last, actual

        stringLen = len (string)
        last = 1
        actual = 1

        do while (actual < stringLen)
            if (string(last:last) == ' ') then
                actual = actual + 1
                string(last:last) = string(actual:actual)
                string(actual:actual) = ' '
            else
                last = last + 1
                if (actual < last) &
                    actual = last
            endif
        end do

    end subroutine

    function upperTriangleIndexing(lower,upper,rowsize)
        integer :: lower, upper,rowsize 
        integer :: upperTriangleIndexing 
        if(lower.ge.upper) stop 'error'
        upperTriangleIndexing = (lower - 1) * rowsize  + &
            upper - lower - (lower*(lower-1))/2
    end function 

end program adf04_add