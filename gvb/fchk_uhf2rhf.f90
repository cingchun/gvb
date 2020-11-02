! written by jxzou at 20180215
! This is a program/subroutine used to transform a Gaussian UHF type .fchk
! into a RHF or ROHF type .fchk file.

! Note: the Spin SCF Density section will be deleted.

! modified by jxzou at 20180328: input file can be a RHF type .fchk file
! modified by jxzou at 20180424: output file can be a ROHF type .fchk file

program main
 implicit none
 integer i
 integer, parameter :: iout = 6
 character(len=240) fchkname

 i = iargc()
 if(i == 0) then
  write(iout,'(A)') 'ERROR in subroutine fchk_uhf2rhf: you must provide a .fch(k) file.'
  write(iout,'(A)') 'Example: fchk_uhf2rhf a.fchk'
  stop
 end if

 fchkname = ' '
 call getarg(1,fchkname)

 call fchk_uhf2rhf(fchkname)
 stop
end program main

subroutine fchk_uhf2rhf(fchkname)
 implicit none
 integer k, nalpha, nbeta
 integer RENAME
 integer, parameter :: fchkid = 11, fchkid1 = 12
 character(len=*),intent(in) :: fchkname
 character(len=240) fchkname1, buffer
 logical rhf

 k = 0
 buffer = ' '
 fchkname1 = ' '
 rhf = .false.

 ! step 1: get nalpha and nbeta
 ! step 2: modify the second line: UHF->RHF or ROHF
 ! step 3: modify the Route section
 ! step 4: modify the 1st value of ILSW into 0
 ! step 5: change IOpCl value to 0; and change IROHF value to 1 (if ROHF)
 ! step 6: skip the Beta Orbital Energies
 ! step 7: copy the Alpha MO coefficients and skip the Beta MO coefficients
 ! step 8: skip the Spin SCF Density

 k = index(fchkname,'.fch')
 fchkname1 = fchkname(1:k-1)//'_r.fchk'
 open(unit=fchkid,file=TRIM(fchkname),status='old',position='rewind')
 open(unit=fchkid1,file=TRIM(fchkname1),status='replace')

 ! step 1: get nalpha and nbeta
 do while(.true.)
  read(fchkid,'(A)') buffer
  if(buffer(1:25) == 'Number of alpha electrons') exit
 end do
 BACKSPACE(fchkid)
 read(fchkid,'(A51,I10)') buffer, nalpha
 read(fchkid,'(A51,I10)') buffer, nbeta
 rewind(fchkid)

 ! step 2: modify the second line: UHF->RHF or ROHF
 read(fchkid,'(A)') buffer
 write(fchkid1,'(A)') TRIM(buffer)
 read(fchkid,'(A)') buffer
 if(buffer(11:13)=='UHF' .and. nalpha==nbeta) then
  buffer(11:13) = 'RHF'
 else if(buffer(11:13)=='UHF' .and. nalpha/=nbeta) then
  buffer(11:14) = 'ROHF'
 end if
 write(fchkid1,'(A)') TRIM(buffer)

 ! step 3: modify the Route section
 if(nalpha /= nbeta) then
  do while(.true.)
   read(fchkid,'(A)') buffer
   write(fchkid1,'(A)') TRIM(buffer)
   if(buffer(1:5) == 'Route') exit
  end do
  read(fchkid,'(A)') buffer
  k = index(buffer,'UHF')
  if(k == 0) k = index(buffer,'uhf')
  if(k /= 0) then
   buffer(k+4:240) = buffer(k+3:239)
   buffer(k:k+3) = 'ROHF'
  end if
  write(fchkid1,'(A)') TRIM(buffer)
 end if

 ! step 4: modify the 1st value of ILSW into 0
 do while(.true.)
  read(fchkid,'(A)') buffer
  write(fchkid1,'(A)') TRIM(buffer)
  if(buffer(1:4) == 'ILSW') exit
 end do
 read(fchkid,'(A)') buffer
 buffer(12:12) = '0'
 write(fchkid1,'(A)') TRIM(buffer)

 ! step 5: change IOpCl value to 0
 do while(.true.)
  read(fchkid,'(A)') buffer
  if(buffer(1:5) == 'IOpCl') exit
  write(fchkid1,'(A)') TRIM(buffer)
 end do
 write(fchkid1,'(A5,38X,A1,16X,A1)') 'IOpCl','I','0'
 ! change IROHF value to 1 (if ROHF)
 if(nalpha /= nbeta) then
  write(fchkid1,'(A5,38X,A1,16X,A1)') 'IROHF','I','1'
  read(fchkid,'(A)') buffer
 end if

 ! step 6: skip the Beta Orbital Energies
 do while(.true.)
  read(fchkid,'(A)') buffer
  if(buffer(1:12) == 'Beta Orbital') exit
  if(buffer(1:8) == 'Alpha MO') exit
  write(fchkid1,'(A)') TRIM(buffer)
 end do
 if(buffer(1:8) == 'Alpha MO') rhf = .true.
 if(.not. rhf) then
  do while(.true.)
   read(fchkid,'(A)') buffer
   if(buffer(1:8) == 'Alpha MO') exit
  end do
 end if
 BACKSPACE(fchkid)

 ! step 7: copy the Alpha MO coefficients and skip the Beta MO coefficients
 do while(.true.)
  read(fchkid,'(A)') buffer
  if(buffer(1:7) == 'Beta MO') exit
  if(buffer(1:11) == 'Orthonormal') exit
  if(buffer(1:9) == 'Total SCF') exit
  write(fchkid1,'(A)') TRIM(buffer)
 end do
 if(.not. rhf) then
  do while(.true.)
   read(fchkid,'(A)') buffer
   if(buffer(1:11) == 'Orthonormal') exit
   if(buffer(1:9) == 'Total SCF') exit
  end do
 end if
 BACKSPACE(fchkid)

 ! step 8: skip the Spin SCF Density
 do while(.true.)
  read(fchkid,'(A)') buffer
  if(buffer(1:8) == 'Spin SCF') exit
  if(buffer(1:16) == 'Mulliken Charges') exit
  write(fchkid1,'(A)') TRIM(buffer)
 end do
 if(.not. rhf) then
  do while(.true.)
   read(fchkid,'(A)') buffer
   if(buffer(1:16) == 'Mulliken Charges') exit
  end do
 end if
 BACKSPACE(fchkid)

 ! copy the rest content
 do while(.true.)
  read(fchkid,'(A)',iostat=k) buffer
  if(k /= 0) exit
  write(fchkid1,'(A)') TRIM(buffer)
 end do
 ! copy done

 close(fchkid)
 close(fchkid1)
 return
end subroutine fchk_uhf2rhf

