! written by jxzou at 20171120
! modified by jxzou at 20171224: use '-gvb npair' in command line to permute the active orbtials
! modified by jxzou at 20180211: support UHF type MOs
! modified by jxzou at 20180620: support open shell orbitals (doublet, triplet, ...)

! copy the MOs in .fch(k) file (Gaussian) to .inp or .dat file (GAMESS)
! Note: an initial .inp or .dat file must be provided, and MOs in it will be replaced.

! If '-gvb [npair]' is specified, the orbitals in active space will be permuted to the
!  order of Gamess. In this case, you must specify the argument [npair] in command line.

! After '-gvb [npair]' is specified, you can also append '-open [nopen]' if your system
!  is truely open shell, i.e., doublet, triplet and so on.

! If '-uhf' is specified, the .fchk file must include the Beta MO part

! The order of Cartesian functions in Gaussian can be acquired by adding keyword 'pop=reg' in .gjf file.
! In GAMESS output files(.out, .gms), the order of Cartesian functions is printed without
!  extra keyword needed.

! Limitation: only deal with Cartesian basis functions.

program main
 implicit none
 integer i, npair, nopen
 integer, parameter :: iout = 6
 character(len=4) gvb_or_uhf
 character(len=4) string
 character(len=240) fchkname,inpname

 npair = 0
 nopen = 0
 gvb_or_uhf = ' '
 string = ' '
 i = iargc()
 if(.not. (i==2 .or. i==3 .or. i==4 .or. i==6)) then
  write(iout,'(A)') 'ERROR in subroutine fchk2vec: the number of arguments in command line is wrong!'
  write(iout,'(A)') 'Example 1 (for RHF, GVB and CASSCF): fchk2vec a.fchk a.inp'
  write(iout,'(A)') 'Example 2 (for GVB): fchk2vec a.fchk a.inp -gvb 4, where 4 means 4 pairs'
  write(iout,'(A)') ' and this subroutine will transform the order of MOs into that of GAMESS.'
  write(iout,'(A)') 'Example 3 (for UHF): fchk2vec a.fchk a.inp -uhf'
  stop
 end if

 call getarg(1,fchkname)
 call getarg(2,inpname)

 if(i==3 .or. i==4 .or. i==6) then
  call getarg(3,gvb_or_uhf)
  if(.not. (gvb_or_uhf=='-gvb' .or. gvb_or_uhf== '-uhf')) then
   write(iout,'(A)') 'ERROR in subroutine fchk2vec: the 3rd argument in command line is wrong!'
   write(iout,'(A)') "It must be '-gvb' or '-uhf'."
   stop
  end if
  if(gvb_or_uhf=='-gvb' .and. i==4) then
   call getarg(4,string)
   read(string,*) npair
  else if(gvb_or_uhf=='-gvb' .and. i==6) then
   call getarg(4,string)
   read(string,*) npair
   call getarg(6,string)
   read(string,*) nopen
  end if
 end if

 call fchk2vec(fchkname, inpname, gvb_or_uhf, npair, nopen)
 stop
end program main

subroutine fchk2vec(fchkname, inpname, gvb_or_uhf, npair, nopen)
 implicit none
 integer i, j, k
 integer nline, nleft
 integer RENAME
 integer ncoeff, nbf, nif
 ! nbf: number of basis functions
 ! nif: number of independent functions, i.e., the number of MOs
 integer nalpha, nbeta, ncore
 integer n10fmark, n15gmark, n21hmark
 integer, allocatable :: f_mark(:), g_mark(:), h_mark(:) ! mark the index where f,g,h functions begin
 integer, parameter :: fchkid = 11, inpid = 12, inpid1 = 13
 integer, parameter :: iout = 6
 integer, intent(in) :: npair, nopen
 integer,allocatable :: order(:), shell_type(:)
 real(kind=8), allocatable :: alpha_coeff(:), alpha_coeff2(:,:)
 real(kind=8), allocatable :: beta_coeff(:), beta_coeff2(:,:)
 real(kind=8), allocatable :: temp_coeff(:,:), open_coeff(:,:)
 character(len=4), intent(in) :: gvb_or_uhf
 character(len=240), intent(in) :: fchkname, inpname
 character(len=240) inpname1, buffer
 logical vec_exist

 buffer = ' ' ! initialization
 ncoeff = 0
 nbf = 0
 nif = 0
 nalpha = 0
 nbeta = 0
 ncore = 0
 vec_exist = .true.

 inpname1 = TRIM(inpname)//'.f2v.tmp'
 open(unit=fchkid,file=TRIM(fchkname),status='old',position='rewind')

 ! find nalpha and nbeta
 do while(.true.)
  read(fchkid,'(A)',iostat=i) buffer
  if(buffer(1:25) == 'Number of alpha electrons') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine fchk2vec: Number of alpha electrons not found in the input file:"
  write(iout,'(A)') TRIM(fchkname)
  close(fchkid)
  return
 end if
 BACKSPACE(fchkid)
 read(fchkid,'(A49,2X,I10)') buffer, nalpha
 read(fchkid,'(A49,2X,I10)') buffer, nbeta
 if(nopen /= nalpha - nbeta) then
  write(iout,'(A)') 'Warning! Input parameter [nopen] is not equal to nalpha-nbeta, still continue.'
  write(iout,'(A)') TRIM(fchkname)
 end if
 ! nalpha and nbeta found

 ! find nbf and nif
 do while(.true.)
  read(fchkid,'(A)',iostat=i) buffer
  if(buffer(1:25) == 'Number of basis functions') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine fchk2vec: number of basis functions not found in the input file:"
  write(iout,'(A)') TRIM(fchkname)
  close(fchkid)
  return
 end if
 BACKSPACE(fchkid)
 read(fchkid,'(A49,2X,I10)') buffer, nbf
 read(fchkid,'(A49,2X,I10)') buffer, nif
 ! nbf and nif found

 ! find the 'Shell types' section
 do while(.true.)
  read(fchkid,'(A)',iostat=i) buffer
  if(buffer(1:11) == 'Shell types') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine fchk2vec: missing the 'Shell types' section in "//TRIM(fchkname)
  close(fchkid)
  return
 end if

 BACKSPACE(fchkid)
 read(fchkid,'(A49,2X,I10)') buffer, k
 allocate(shell_type(k))
 shell_type = 0
 read(fchkid,'(6(6X,I6))') (shell_type(i),i=1,k)
 ! read 'Shell types' done

 ! check if any spherical functions
 if( ANY(shell_type<-1) ) then
  write(iout,'(A)') 'ERROR in subroutine fchk2vec: spherical functions detected in '//TRIM(fchkname)
  close(fchkid)
  return
 end if
 ! check done

 ! record the indices of Cartesian f, g and h functions
 n10fmark = 0
 n15gmark = 0
 n21hmark = 0
 allocate(f_mark(k), g_mark(k), h_mark(k))
 f_mark = 0
 g_mark = 0
 h_mark = 0
 nbf = 0
 do i = 1,k,1
  if(shell_type(i) == 0) then       !'S'
   nbf = nbf + 1
  else if(shell_type(i) == 1) then  !'P'
   nbf = nbf + 3
  else if(shell_type(i) == -1) then !'SP' or 'L'
   nbf = nbf + 4
  else if(shell_type(i) == 2) then  ! 6D
   nbf = nbf + 6
  else if(shell_type(i) == 3) then  ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  else if(shell_type(i) == 4) then  ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf + 1
   nbf = nbf + 15
  else if(shell_type(i) == 5) then  ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf + 1
   nbf = nbf + 21
  end if
 end do
 deallocate(shell_type)
 ! done recording

 ! find strings 'Alpha MO'
 do while(.true.)
  read(fchkid,'(A)') buffer
  if(buffer(1:8) == 'Alpha MO') exit
 end do

 ! read ncoeff
 BACKSPACE(fchkid)
 read(fchkid,'(A49,2X,I10)') buffer, ncoeff
 ! ncoeff get

 if(nbf*nif /= ncoeff) then
  write(iout,'(A)') 'ERROR in subroutine fchk2vec: nbf*nif/=ncoeff in '//TRIM(fchkname)//', cannot deal with that.'
  close(fchkid)
  return
 end if

 ! read Alpha MO
 allocate(alpha_coeff(ncoeff))
 alpha_coeff = 0.0d0
 read(fchkid,'(5(1X,ES15.8))') (alpha_coeff(i),i=1,ncoeff)
 allocate(alpha_coeff2(nbf,nif))
 alpha_coeff2 = RESHAPE(alpha_coeff,(/nbf,nif/))
 deallocate(alpha_coeff)

 ! if '-uhf' is specified, read Beta MO
 if(gvb_or_uhf == '-uhf') then
  allocate(beta_coeff(ncoeff))
  beta_coeff = 0.0d0
  read(fchkid,'(A)') buffer   ! skip the Beta MO line
  read(fchkid,'(5(1X,ES15.8))') (beta_coeff(i),i=1,ncoeff)
  allocate(beta_coeff2(nbf,nif))
  beta_coeff2 = RESHAPE(beta_coeff,(/nbf,nif/))
  deallocate(beta_coeff)
 end if

 ! all MO in .fchk file read done
 close(fchkid)

 ! adjust the order of Cartesian f, g, etc functions
 do i = 1,n10fmark,1
  call fchk2vec_permute_10f(nif,alpha_coeff2(f_mark(i)+3:f_mark(i)+8,:))
 end do
 do i = 1,n15gmark,1
  call fchk2vec_permute_15g(nif,alpha_coeff2(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1,n21hmark,1
  call fchk2vec_permute_21h(nif,alpha_coeff2(h_mark(i):h_mark(i)+20,:))
 end do
 if(gvb_or_uhf == '-uhf') then
  do i = 1,n10fmark,1
   call fchk2vec_permute_10f(nif,beta_coeff2(f_mark(i)+3:f_mark(i)+8,:))
  end do
  do i = 1,n15gmark,1
   call fchk2vec_permute_15g(nif,beta_coeff2(g_mark(i):g_mark(i)+14,:))
  end do
  do i = 1,n21hmark,1
   call fchk2vec_permute_21h(nif,beta_coeff2(h_mark(i):h_mark(i)+20,:))
  end do
 end if
 deallocate(f_mark, g_mark, h_mark)
 ! adjustment finished

 ! if active orbitals in GAMESS order are required, permute the active orbitals
 if(gvb_or_uhf=='-gvb' .and. npair>1) then
  ncore = nalpha - npair - nopen
  allocate(order(2*npair))
  allocate(temp_coeff(nbf,2*npair))
  order = 0
  temp_coeff = 0.0d0

  if(nopen > 0) then
   allocate(open_coeff(nbf,nopen))  
   open_coeff(1:nbf,1:nopen) = alpha_coeff2(1:nbf,nalpha-nopen+1:nalpha)
  end if
  forall(i = 1:npair)
   order(2*i-1) = i
  end forall
  forall(i = 1:npair)
   order(2*i) = 2*npair + 1 - i + nopen
  end forall
  forall(i = 1:2*npair)
   temp_coeff(1:nbf,i) = alpha_coeff2(1:nbf,order(i)+ncore)
  end forall
  forall(i = 1:2*npair)
   alpha_coeff2(1:nbf,i+ncore+nopen) = temp_coeff(1:nbf,i)
  end forall
  deallocate(order, temp_coeff)
  if(nopen > 0) then
   alpha_coeff2(1:nbf,ncore+1:ncore+nopen) = open_coeff(1:nbf,1:nopen)
   deallocate(open_coeff)
  end if
 end if
 ! done permute

 ! output the MOs to the .inp (or .dat) file
 open(unit=inpid,file=TRIM(inpname),status='old',position='rewind')
 open(unit=inpid1,file=TRIM(inpname1),status='replace')
 do while(.true.)
  read(inpid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  write(inpid1,'(A)') TRIM(buffer)
  k = index(buffer,'$VEC')
  if(k == 0) k = index(buffer,'$vec')
  if(k == 0) k = index(buffer,'$Vec')
  if(k /= 0) exit
 end do
 if(i /= 0) then
  vec_exist = .false.
  write(iout,'(A)') "Warning! No '$VEC' section in the input file "//TRIM(inpname)//', still continue.'
 end if
 if(.not. vec_exist) write(inpid1,'(1X,A4)') '$VEC'

 ! write Alpha MO
 nline = nbf/5
 nleft = nbf - 5*nline
 do i = 1,nif, 1
  k = MOD(i,100)
  do j = 1, nline, 1
   write(inpid1,'(I2,I3,5ES15.8)') k, MOD(j,1000), alpha_coeff2(5*j-4:5*j,i)
  end do
  if(nleft > 0) then
   write(inpid1,'(I2,I3,5ES15.8)') k, MOD(j,1000), alpha_coeff2(5*j-4:nbf,i)
  end if
 end do
 ! write Beta MO, if any
 if(gvb_or_uhf == '-uhf') then
  do i = 1,nif, 1
   k = MOD(i,100)
   do j = 1, nline, 1
    write(inpid1,'(I2,I3,5ES15.8)') k, MOD(j,1000), beta_coeff2(5*j-4:5*j,i)
   end do
   if(nleft > 0) then
    write(inpid1,'(I2,I3,5ES15.8)') k, MOD(j,1000), beta_coeff2(5*j-4:nbf,i)
   end if
  end do
 end if
 write(inpid1,'(1X,A4)') '$END'
 ! output MOs finished

 deallocate(alpha_coeff2)
 if(allocated(beta_coeff2)) deallocate(beta_coeff2)

 ! skip the MOs in the file inpname
 if(vec_exist) then
  do while(.true.)
   read(inpid,'(A)') buffer
   i = index(buffer,'$END')
   j = index(buffer,'$end')
   k = index(buffer,'$End')
   if(i/=0 .or. j/=0 .or. k/=0) exit
  end do
 end if
 ! skip done

 ! copy the rest of the file inpname
 do while(.true.)
  read(inpid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  write(inpid1,'(A)') TRIM(buffer)
 end do
 ! copy finished

 close(inpid,status='delete')
 close(inpid1)
 i = RENAME(TRIM(inpname1),TRIM(inpname))
 return
end subroutine fchk2vec

subroutine fchk2vec_permute_10f(nif,coeff)
 implicit none
 integer i
 integer,intent(in) :: nif
 integer order(6)
 real(kind=8),intent(inout) :: coeff(6,nif)
 real(kind=8) coeff2(6,nif)
 data order /2, 3, 1, 6, 4, 5/
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in Gamess
!                     1    2    3    4    5    6
! From: XXX,YYY,ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ,XYZ
! To:   XXX,YYY,ZZZ, XXY, XXZ, XYY, YYZ, XZZ, YZZ,XYZ

 coeff2 = 0.0d0
 forall(i = 1:6)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2vec_permute_10f

subroutine fchk2vec_permute_15g(nif,coeff)
 implicit none
 integer i
 integer,intent(in) :: nif
 integer order(15)
 real(kind=8),intent(inout) :: coeff(15,nif)
 real(kind=8) coeff2(15,nif)
 data order /15, 5, 1, 14, 13, 9, 4, 6, 2, 12, 10, 3, 11, 8, 7/
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in Gamess
!       1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! From: ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! To:   XXXX,YYYY,ZZZZ,XXXY,XXXZ,XYYY,YYYZ,XZZZ,YZZZ,XXYY,XXZZ,YYZZ,XXYZ,XYYZ,XYZZ

 coeff2 = 0.0d0
 forall(i = 1:15)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2vec_permute_15g

subroutine fchk2vec_permute_21h(nif,coeff)
 implicit none
 integer i
 integer,intent(in) :: nif
 integer order(21)
 real(kind=8),intent(inout) :: coeff(21,nif)
 real(kind=8) coeff2(21,nif)
 data order /21, 6, 1, 20, 19, 11, 5, 7, 2, 18, 16, 15, 4, 12, 3, 17, 10, 8, 14, 13, 9/
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in Gamess
!       1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! From: ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! To:   XXXXX,YYYYY,ZZZZZ,XXXXY,XXXXZ,XYYYY,YYYYZ,XZZZZ,YZZZZ,XXXYY,XXXZZ,XXYYY,YYYZZ,XXZZZ,YYZZZ,XXXYZ,XYYYZ,XYZZZ,XXYYZ,XXYZZ,XYYZZ

 coeff2 = 0.0d0
 forall(i = 1:21)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2vec_permute_21h

