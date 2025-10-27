program adf04_red 
    !adf04_red - requires the .f90 file name extension.
    !reduces an adf04 file to a lower number of levels
    !it is the user's responsibility to know what they are 
    !reducing. 
    !Maintained at: github link.
    implicit none 
    character*320 :: templine
    logical       :: iamMinusOne = .false.
    integer       :: ii,jj ,upp,low
    integer       :: rstat,wstat
    integer       :: max_i       = 101 
    integer       :: counter
    integer       :: originalNumberLevels

    open(1,file='adf04') 
    open(2,file='adf04_reduced')
    read (1,'(A320)') templine
    write(2,'(A320)') templine
    originalNumberLevels = 0
    do while (.not. iamMinusOne) 
        read (1,'(I5,A320)',iostat = rstat)      ii,templine
        if (ii.le.max_i) then
            write(2,'(I5,A320)',iostat = wstat)  ii,templine
        end if 
        originalNumberLevels = originalNumberLevels + 1 
        if (ii.eq.-1) iamMinusOne = .true.
    end do

    write(0,*) 'levels finished'
    read (1,'(A320)',iostat = rstat) templine
    write(2,'(A320)',iostat = wstat) templine

    iamMinusOne = .false.
    do while (.not. iamMinusOne) 
        read (1,'(2I5,A320)',iostat = rstat)     jj,ii,templine
        if (rstat .lt.0) stop 'read error - emergency stop.'

        upp = max(jj,ii)
        low = min(jj,ii)
        if (ii.eq.-1 .or. jj.eq.-1) then 
            iamMinusOne = .true.
            write(2,'(I4)') -1 
        else if (upp.le.max_i) then
            write(2,'(2I5,A320)',iostat = wstat) jj,ii,templine
        end if 

        if (wstat .lt.0) stop 'write error - emergency stop.'


    end do
    read  (1,'(A320)',iostat = rstat) templine
    write (2,'(A320)',iostat = wstat) templine
    counter = 0 
    write (2, '(A)',advance = 'no') 'C'
    write (2,'(A,I5,A,I5,A)') ' Reduced from ',originalNumberLevels,   & 
    ' levels to ',max_i,&
    ' levels by the code of Leo Mulholland (github link)' 

    do while (rstat .ge. 0) 
        read  (1,'(A320)',iostat = rstat) templine
        write (2,'(A320)',iostat = wstat) templine
        counter = counter + 1
    end do 

    close(1)
    close(2)

end program adf04_red