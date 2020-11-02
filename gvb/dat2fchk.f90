! written by jxzou at 20171108
! modified by jxzou at 20180222: support UHF type MOs

! Copy the MOs in .dat or .inp file (Gamess) to .fch(k) file (Gaussian).
! The orders of Cartesian f,g and h functions in .dat file will be permuted.
! Note: an initial .fch(k) file must be provided, and MOs in it will be replaced.

! The order of Cartesian functions can be known by adding keyword 'pop=reg' in .gjf file.
! In Gamess output files(.out, .gms), the order of Cartesian functions is printed without
!  extra keyword needed.
! Limitation: only deal with Cartesian basis functions

program main
 implicit none
 integer i, npair
 integer, parameter :: iout = 6
 character(len=4) gvb_or_uhf,string
 character(len=240) datname,fchkname

 npair = 0
 gvb_or_uhf = ' '
 string = ' '
 datname = ' '
 fchkname = ' '

 i = iargc()
 if(.not. (i==2 .or. i==3 .or. i==4)) then
  write(iout,'(A)') 'ERROR in subroutine dat2fchk: the number of arguments in command line is wrong!'
  write(iout,'(A)') 'Example 1 (for RHF, GVB and CASSCF): datfchk a.dat a.fchk'
  write(iout,'(A)') 'Example 2 (for GVB): datfchk a.dat a.fchk -gvb 4, where 4 means 4 pairs'
  write(iout,'(A)') ' and this subroutine will transform the order of MOs into that of Gaussian.'
  write(iout,'(A)') 'Example 3 (for UHF): datfchk a.dat a.fchk -uhf'
  stop
 end if

 call getarg(1,datname)
 call getarg(2,fchkname)

 if(i==3 .or. i==4) then
  call getarg(3,gvb_or_uhf)
  if(.not. (gvb_or_uhf=='-gvb' .or. gvb_or_uhf=='-uhf')) then
   write(iout,'(A)') 'ERROR in subroutine dat2fchk: the 3rd argument in command line is wrong!'
   write(iout,'(A)') "It must be '-gvb' or '-uhf'."
   stop
  end if
  if(gvb_or_uhf=='-gvb' .and. i==4) then
   call getarg(4,string)
   read(string,*) npair
  end if
 end if

 call dat2fchk(datname, fchkname, gvb_or_uhf, npair)
 stop
end program main

subroutine dat2fchk(datname,fchkname,gvb_or_uhf,npair)
 implicit none
 integer i,j,k
 integer RENAME
 integer ncoeff, nbf, nif, nocc
 ! nbf: the number of basis functions
 ! nif: the number of independent functions
 ! nocc: the number of occupied orbitals
 integer nalpha, nbeta
 integer nline, nleft, nstart
 integer n10fmark, n15gmark, n21hmark
 integer,allocatable :: f_mark(:), g_mark(:), h_mark(:) ! mark the index where f,g,h functions begin
 integer,parameter :: datid = 11, fchkid = 12, fchkid1 = 13
 integer,intent(in) :: npair
 integer,allocatable :: order(:)
 real(kind=8),allocatable :: alpha_coeff(:), alpha_coeff2(:,:)
 real(kind=8),allocatable :: beta_coeff(:), beta_coeff2(:,:)
 real(kind=8),allocatable :: temp_coeff(:,:)
 character(len=4),intent(in) :: gvb_or_uhf
 character(len=5) str1
 character(len=30) str2
 character(len=240),intent(in) :: datname, fchkname
 character(len=240) fchkname1, buffer

 str1 = ' '   ! initialization
 str2 = ' '
 buffer = ' '
 ncoeff = 0
 nalpha = 0
 nbeta = 0
 nline = 0
 nleft = 0
 nstart = 0

 fchkname1 = TRIM(fchkname)//'.d2f.tmp'
 open(unit=fchkid,file=TRIM(fchkname),status='old',position='rewind')

 ! find nalpha and nbeta
 do while(.true.)
  read(fchkid,'(A)',iostat=i) buffer
  if(buffer(1:25) == 'Number of alpha electrons') exit
 end do
 if(i /= 0) then
  write(*,'(A)') "ERROR in subroutine dat2fchk. Number of alpha electrons not found in the input file:"
  close(fchkid)
  return
 end if
 BACKSPACE(fchkid)
 read(fchkid,'(A49,2X,I10)') buffer, nalpha
 read(fchkid,'(A49,2X,I10)') buffer, nbeta
 nocc = nalpha
 ! found nalpha and nbeta

 ! find nbf and nif
 do while(.true.)
  read(fchkid,'(A)',iostat=i) buffer
  if(buffer(1:25) == 'Number of basis functions') exit
 end do
 if(i /= 0) then
  write(*,'(A)') "ERROR in subroutine dat2fchk: Number of basis functions not found in the input file:"
  write(*,'(A)') TRIM(fchkname1)
  close(fchkid)
  return
 end if
 BACKSPACE(fchkid)
 read(fchkid,'(A49,2X,I10)') buffer, nbf
 read(fchkid,'(A49,2X,I10)') buffer, nif
 ! nbf and nif found

 ! find strings 'Alpha MO'
 do while(.true.)
  read(fchkid,'(A)') buffer
  if(buffer(1:8) == 'Alpha MO') exit
 end do

 ! find ncoeff
 BACKSPACE(fchkid)
 read(fchkid,'(A49,2X,I10)') buffer, ncoeff
 ! ncoeff found

 if(nbf*nif /= ncoeff) then
  write(*,'(A)') 'ERROR in subroutine dat2fchk: nbf*nif/=ncoeff in '//TRIM(fchkname)//'. Cannot deal with that.'
  close(fchkid)
  return
 end if

 open(unit=datid,file=TRIM(datname),status='old',position='rewind')
 do while(.true.)
  read(datid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  k = index(buffer,'$VEC')
  if(k == 0) k = index(buffer,'$vec')
  if(k == 0) k = index(buffer,'$Vec')
  if(k /= 0) exit
 end do
 if(i /= 0) then
  write(*,'(A)') "ERROR in subroutine dat2fchk. No '$VEC' section in the input file:"
  write(*,'(A)') TRIM(datname)
  close(datid)
  close(fchkid)
  close(fchkid1,status='delete')
  return
 end if

 ! read Alpha MO
 allocate(alpha_coeff2(nbf,nif))
 alpha_coeff2 = 0.0d0
 nline = nbf/5
 nleft = nbf - nline*5
 do i = 1,nif,1
  k = 1
  do j = 1,nline,1
   read(datid,'(A)') buffer
   buffer = buffer(6:)
   read(buffer,'(5ES15.8)') alpha_coeff2(k:k+4,i)
   k = k + 5
  end do
  if(nleft > 0) then
   read(datid,'(A)') buffer
   buffer = buffer(6:)
   str1 = ' '
   write(str1,'(I5)') nleft
   str1 = ADJUSTL(str1)
   str2 = '('//TRIM(str1)//'ES15.8)'
   read(buffer,TRIM(str2)) alpha_coeff2(k:nbf,i)
  end if
 end do

 ! if '-uhf' is specified, read Beta MO
 if(gvb_or_uhf == '-uhf') then
  allocate(beta_coeff2(nbf,nif))
  beta_coeff2 = 0.0d0
  do i = 1,nif,1
   k = 1
   do j = 1,nline,1
    read(datid,'(A)') buffer
    buffer = buffer(6:)
    read(buffer,'(5ES15.8)') beta_coeff2(k:k+4,i)
    k = k + 5
   end do
   if(nleft > 0) then
    read(datid,'(A)') buffer
    buffer = buffer(6:)
    str1 = ' '
    write(str1,'(I5)') nleft
    str1 = ADJUSTL(str1)
    str2 = '('//TRIM(str1)//'ES15.8)'
    read(buffer,TRIM(str2)) beta_coeff2(k:nbf,i)
   end if
  end do
 end if

 ! find the $DATA section and record the indices of Cartesian f, g and h functions
 rewind(datid)
 do while(.true.)
  read(datid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  k = index(buffer,'$DATA')
  if(k == 0) k = index(buffer,'$data')
  if(k == 0) k = index(buffer,'$Data')
  if(k /= 0) exit
 end do
 if(i /= 0) then
  write(*,'(A)') "ERROR in subroutine dat2fchk: no '$DATA' section in the input file:"
  write(*,'(A)') TRIM(datname)
  close(datid)
  close(fchkid)
  return
 end if

 ! skip the point group and Title Card lines
 do i = 1,2
  read(datid,'(A)') buffer
 end do

 n10fmark = 0
 n15gmark = 0
 n21hmark = 0
 allocate(f_mark(nbf), g_mark(nbf), h_mark(nbf))
 f_mark = 0
 g_mark = 0
 h_mark = 0
 nbf = 0
 do while(.true.)
  read(datid,'(A)') buffer
  buffer = ADJUSTL(buffer)
  if(buffer(1:1)=='$' .or. buffer(1:1)==' ') exit ! meet '$END' or blank line
  do while(.true.)
   read(datid,'(A)') buffer
   buffer = ADJUSTL(buffer)
   if(buffer(1:1)=='$' .or. buffer(1:1)==' ') exit ! meet '$END' or blank line
   if(buffer(1:1) == 'S') then
    nbf = nbf + 1
   else if(buffer(1:1) == 'P') then
    nbf = nbf + 3
   else if(buffer(1:1) == 'L') then ! 'L' is 'SP' in Gaussian format
    nbf = nbf + 4
   else if(buffer(1:1) == 'D') then
    nbf = nbf + 6
   else if(buffer(1:1) == 'F') then
    n10fmark = n10fmark + 1
    f_mark(n10fmark) = nbf + 1
    nbf = nbf + 10
   else if(buffer(1:1) == 'G') then
    n15gmark = n15gmark + 1
    g_mark(n15gmark) = nbf + 1
    nbf = nbf + 15
   else if(buffer(1:1) == 'H') then
    n21hmark = n21hmark + 1
    h_mark(n21hmark) = nbf + 1
    nbf = nbf + 21
   end if
   call get_int_after_flag(buffer,' ',nline)
   do i = 1,nline,1  ! skip n lines
    read(datid,'(A)') buffer
   end do
  end do
 end do
 close(datid)
 ! done recording

 ! adjust the order of Cartesian f, g, h functions
 do i = 1,n10fmark,1
  call dat2fchk_permute_10f(nif,alpha_coeff2(f_mark(i)+3:f_mark(i)+8,:))
 end do
 do i = 1,n15gmark,1
  call dat2fchk_permute_15g(nif,alpha_coeff2(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1,n21hmark,1
  call dat2fchk_permute_21h(nif,alpha_coeff2(h_mark(i):h_mark(i)+20,:))
 end do
 if(gvb_or_uhf == '-uhf') then
  do i = 1,n10fmark,1
   call dat2fchk_permute_10f(nif,beta_coeff2(f_mark(i)+3:f_mark(i)+8,:))
  end do
  do i = 1,n15gmark,1
   call dat2fchk_permute_15g(nif,beta_coeff2(g_mark(i):g_mark(i)+14,:))
  end do
  do i = 1,n21hmark,1
   call dat2fchk_permute_21h(nif,beta_coeff2(h_mark(i):h_mark(i)+20,:))
  end do
 end if
 deallocate(f_mark, g_mark ,h_mark)
 ! adjust Cartesian functions finished

 ! if orbitals in GVB order are required, permute the active orbitals
 if(gvb_or_uhf=='-gvb' .and. npair>1) then
  allocate(order(2*npair))
  allocate(temp_coeff(nbf,2*npair))
  order = 0
  temp_coeff = 0.0d0
  forall(i = 1:npair)
   order(i) = 2*i - 1
  end forall
  forall(i = npair+1:2*npair)
   order(i) = 2*(2*npair-i+1)
  end forall
  nstart = nocc - npair + 1
  forall(i = 1:2*npair)
   order(i) = order(i) + nstart - 1
  end forall
  forall(i = 1:2*npair)
   temp_coeff(1:nbf,i) = alpha_coeff2(1:nbf,order(i))
  end forall
  forall(i = 1:2*npair)
   alpha_coeff2(1:nbf,i+nstart-1) = temp_coeff(1:nbf,i)
  end forall
  deallocate(order,temp_coeff)
 end if
 ! permute done

 ! output the MOs to .fch(k) file
 open(unit=fchkid1,file=TRIM(fchkname1),status='replace')
 rewind(fchkid)
 do while(.true.)
  read(fchkid,'(A)') buffer
  write(fchkid1,'(A)') TRIM(buffer)
  if(buffer(1:8) == 'Alpha MO') exit
 end do

 allocate(alpha_coeff(ncoeff))
 alpha_coeff = 0.0d0
 alpha_coeff = RESHAPE(alpha_coeff2,(/ncoeff/))
 write(fchkid1,'(5(1X,ES15.8))') (alpha_coeff(i),i=1,ncoeff)
 deallocate(alpha_coeff, alpha_coeff2)
 ! write Beta MOs, if any
 if(gvb_or_uhf == '-uhf') then
  do while(.true.) ! skip the original Beta MO in .fchk file
   read(fchkid,'(A)') buffer
   if(buffer(1:7) == 'Beta MO') exit
  end do
  write(fchkid1,'(A)') TRIM(buffer)
  allocate(beta_coeff(ncoeff))
  beta_coeff = 0.0d0
  beta_coeff = RESHAPE(beta_coeff2,(/ncoeff/))
  write(fchkid1,'(5(1X,ES15.8))') (beta_coeff(i),i=1,ncoeff)
  deallocate(beta_coeff, beta_coeff2)
 end if
 ! output MOs finished

 ! copy the rest of the .fch(k) file
 do while(.true.)
  read(fchkid,'(A)') buffer
  if(index(buffer,'=') /= 0) exit
 end do
 BACKSPACE(fchkid)

 do while(.true.)
  read(fchkid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  write(fchkid1,'(A)') TRIM(buffer)
 end do
 ! copy finished

 close(fchkid,status='delete')
 close(fchkid1)
 i = RENAME(TRIM(fchkname1),TRIM(fchkname))
 return
end subroutine dat2fchk

subroutine dat2fchk_permute_10f(nif,coeff)
 implicit none
 integer i,nif
 integer order(6)
 real(kind=8),intent(inout) :: coeff(6,nif)
 real(kind=8) coeff2(6,nif)
 data order /3, 1, 2, 5, 6, 4/
! From: the order of Cartesian f functions in Gamess
! To: the order of Cartesian f functions in Gaussian
!                     1    2    3    4    5    6
! From: XXX,YYY,ZZZ, XXY, XXZ, XYY, YYZ, XZZ, YZZ, XYZ
! To:   XXX,YYY,ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ

 coeff2 = 0.0d0
 forall (i = 1:6)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine dat2fchk_permute_10f

subroutine dat2fchk_permute_15g(nif,coeff)
 implicit none
 integer i
 integer,intent(in) :: nif
 integer order(15)
 real(kind=8),intent(inout) :: coeff(15,nif)
 real(kind=8) coeff2(15,nif)
 data order /3, 9, 12, 7, 2, 8, 15, 14, 6, 11, 13, 10, 5, 4, 1/
! From: the order of Cartesian g functions in Gamess
! To: the order of Cartesian g functions in Gaussian
!       1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! From: XXXX,YYYY,ZZZZ,XXXY,XXXZ,XYYY,YYYZ,XZZZ,YZZZ,XXYY,XXZZ,YYZZ,XXYZ,XYYZ,XYZZ
! To:   ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX

 coeff2 = 0.0d0
 forall(i = 1:15)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine dat2fchk_permute_15g

subroutine dat2fchk_permute_21h(nif,coeff)
 implicit none
 integer i,nif
 integer order(21)
 real(kind=8),intent(inout) :: coeff(21,nif)
 real(kind=8) coeff2(21,nif)
 data order /3, 9, 15, 13, 7, 2, 8, 18, 21, 17, 6, 14, 20, 19, 12, 11, 16, 10, 5, 4, 1/
! From: the order of Cartesian h functions in Gamess
! To: the order of Cartesian h functions in Gaussian
!       1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! From: XXXXX,YYYYY,ZZZZZ,XXXXY,XXXXZ,XYYYY,YYYYZ,XZZZZ,YZZZZ,XXXYY,XXXZZ,XXYYY,YYYZZ,XXZZZ,YYZZZ,XXXYZ,XYYYZ,XYZZZ,XXYYZ,XXYZZ,XYYZZ
! To:   ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX

 coeff2 = 0.0d0
 forall(i = 1:21)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine dat2fchk_permute_21h

subroutine get_int_after_flag(buffer,flag,k)
 implicit none
 integer,intent(out) :: k
 character(len=1),intent(in) :: flag
 character(len=240),intent(inout) :: buffer

 k = index(buffer,flag)
 if(k == 0) then
  write(*,'(A)') 'ERROR in subroutine get_int_after_flag: sth wrong in this line:'
  write(*,'(A)') TRIM(buffer)
  stop
 end if
 buffer(1:k) = ' '
 buffer = ADJUSTL(buffer)
 read(buffer,*) k
 return
end subroutine get_int_after_flag

