! written by jxzou at 20171203: adjust the orders of d,f,g, etc functions in
!  .fch(k) file, to the orders in PySCF

! modified by jxzou at 20180314: add a input parameter 'a' or 'b' to read Alpha or Beta MO in .fchk
! modified by jxzou at 20180323: support the case that D functions preceding L functions
! modified by jxzou at 20180406: code optimization
! modified by jxzou at 20180520: support linear dependence

! This subroutine is designed to be called as a module by Python.
!  For INTEL compiler, use
! -----------------------------------------------------------------------
!  f2py -m fchk2py -c fchk2py.f90 --fcompiler=intelem --compiler=intelem
! -----------------------------------------------------------------------
!  For GNU compiler, use
! --------------------------------
!  f2py -m fchk2py -c fchk2py.f90
! --------------------------------

!  to compile this file (a fchk2py.so file will be generated). Then in
!  Python you can import the fchk2py module.

!program main
! implicit none
! integer nbf
! real(kind=8),allocatable :: coeff2(:,:)
! character(len=240) fname
!
! fname = ' '
! coeff2 = 0.0d0
! call getarg(1,fname)
! nbf = 38
! allocate(coeff2(nbf,nbf))
! call fchk2py(fname,nbf,coeff2,'a')
! stop
!end program

! read the MOs in .fch(k) file and adjust its d,f,g, etc. functions order
!  of Gaussian to that of PySCF
subroutine fchk2py(fchkname, nbf, nif, coeff2, ab)
 implicit none
 integer i, k
 integer length
 integer ncoeff,nbf,nif
!f2py intent(in) :: nbf, nif
 integer n6dmark,n10fmark,n15gmark,n21hmark
 integer n5dmark,n7fmark, n9gmark, n11hmark
 integer,allocatable :: shell_type(:), shell_to_atom_map(:)
 integer,parameter :: fchkid = 11
 ! mark the index where d, f, g, h functions begin
 integer,allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 character(len=1) ab
!f2py intent(in) :: ab
 character(len=8) key
 character(len=8), parameter :: key1 = 'Alpha MO'
 character(len=7), parameter :: key2 = 'Beta MO'
 character(len=240) fchkname
!f2py intent(in) :: fchkname
 character(len=240) buffer
 real(kind=8),allocatable :: coeff(:)
 real(kind=8) coeff2(nbf,nif)
!f2py intent(out) :: coeff2

 key = ' '
 buffer = ' '
 ncoeff = 0

 key = key1
 if(ab/='a' .and. ab/='A') then
  key = key2//' '
 end if

 open(unit=fchkid,file=TRIM(fchkname),status='old',position='rewind')
 do while(.true.)
  read(fchkid,'(A)') buffer
  if(buffer(1:8) == key) exit
 end do
 BACKSPACE(fchkid)
 read(fchkid,'(A49,2X,I10)') buffer, ncoeff

! nif = nbf
 if(ncoeff /= nbf*nif) then
  write(*,'(A)') 'ERROR in subroutine fchk2py: the MOs in fchk is not square!'
  write(*,'(A)') 'The number of basis functions is not equal to that of MOs.'
  write(*,'(A)') TRIM(fchkname)
  close(fchkid)
  return
 end if

 ! read Alpha MO or Beta MO
 allocate(coeff(ncoeff))
 coeff = 0.0d0
 read(fchkid,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
 coeff2 = RESHAPE(coeff,(/nbf,nif/))
 deallocate(coeff)

 rewind(fchkid)
 ! find and read Shell types
 do while(.true.)
  read(fchkid,'(A)',iostat=i) buffer
  if(buffer(1:11) == 'Shell types') exit
 end do
 if(i /= 0) then
  write(*,'(A)') "ERROR in subroutine fchk2py: missing the 'Shell types' section in .fchk file!"
  write(*,'(A)') TRIM(fchkname)
  close(fchkid)
  return
 end if

 BACKSPACE(fchkid)
 read(fchkid,'(A49,2X,I10)') buffer, k
 allocate(shell_type(2*k))
 shell_type = 0
 read(fchkid,'(6(6X,I6))') (shell_type(i),i=1,k)
 ! read Shell types done

 ! find and read Shell to atom map
 do while(.true.)
  read(fchkid,'(A)',iostat=i) buffer
  if(buffer(1:13) == 'Shell to atom') exit
 end do
 if(i /= 0) then
  write(*,'(A)') "ERROR in subroutine fchk2py: missing the 'Shell to atom map' section in .fchk file!"
  write(*,'(A)') TRIM(fchkname)
  close(fchkid)
  return
 end if
 allocate(shell_to_atom_map(2*k))
 shell_to_atom_map = 0
 read(fchkid,'(6(6X,I6))') (shell_to_atom_map(i),i=1,k)
 ! read Shell to atom map done

 ! all information in .fchk file read done
 close(fchkid)

! first we adjust the basis functions in each MO according to the Shell to atom map
! this is to ensure that D comes after L functions
 ! split the 'L' into 'S' and 'P'
 call split_L_func(k, shell_type, shell_to_atom_map, length)

 ! sort the shell_type, shell_to_atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo(length, shell_type, shell_to_atom_map, nbf, nif, coeff2)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 n6dmark = 0
 n10fmark = 0
 n15gmark = 0
 n21hmark = 0
 n5dmark = 0
 n7fmark = 0
 n9gmark = 0
 n11hmark = 0
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k))
 d_mark = 0
 f_mark = 0
 g_mark = 0
 h_mark = 0
 nbf = 0
 do i = 1,k,1
  select case(shell_type(i))
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-1)   ! SP or L
   nbf = nbf + 4
  case(-2)   ! 5D
   n5dmark = n5dmark + 1
   d_mark(n5dmark) = nbf + 1
   nbf = nbf + 5
  case( 2)   ! 6D
   n6dmark = n6dmark + 1
   d_mark(n6dmark) = nbf + 1
   nbf = nbf + 6
  case(-3)   ! 7F
   n7fmark = n7fmark + 1
   f_mark(n7fmark) = nbf + 1
   nbf = nbf + 7
  case( 3)   ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  case(-4)   ! 9G
   n9gmark = n9gmark + 1
   g_mark(n9gmark) = nbf + 1
   nbf = nbf + 9
  case( 4)   ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf + 1
   nbf = nbf + 15
  case(-5)   ! 11H
   n11hmark = n11hmark + 1
   h_mark(n11hmark) = nbf + 1
   nbf = nbf + 11
  case( 5)   ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf + 1
   nbf = nbf + 21
  end select
 end do
 deallocate(shell_type)

 ! adjust the order of d, f, etc. functions
 do i = 1,n5dmark,1
  call fchk2py_permute_5d(nif,coeff2(d_mark(i):d_mark(i)+4,:))
 end do
 do i = 1,n6dmark,1
  call fchk2py_permute_6d(nif,coeff2(d_mark(i):d_mark(i)+5,:))
 end do
 do i = 1,n7fmark,1
  call fchk2py_permute_7f(nif,coeff2(f_mark(i):f_mark(i)+6,:))
 end do
 do i = 1,n10fmark,1
  call fchk2py_permute_10f(nif,coeff2(f_mark(i):f_mark(i)+9,:))
 end do
 do i = 1,n9gmark,1
  call fchk2py_permute_9g(nif,coeff2(g_mark(i):g_mark(i)+8,:))
 end do
 do i = 1,n15gmark,1
  call fchk2py_permute_15g(nif,coeff2(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1,n11hmark,1
  call fchk2py_permute_11h(nif,coeff2(h_mark(i):h_mark(i)+10,:))
 end do
 do i = 1,n21hmark,1
  call fchk2py_permute_21h(nif,coeff2(h_mark(i):h_mark(i)+20,:))
 end do
! adjustment finished

 deallocate(d_mark, f_mark, g_mark, h_mark)
 return
end subroutine fchk2py

! split the 'L' into 'S' and 'P'
subroutine split_L_func(k, shell_type, shell_to_atom_map, length)
 implicit none
 integer i, k0
 integer,intent(in) :: k
 integer,intent(inout) :: shell_type(2*k), shell_to_atom_map(2*k)
 integer,intent(out) :: length
 integer,allocatable :: temp1(:), temp2(:)

 k0 = 2*k
 length = k
 ! set initial values for arrays shell_type, assume 15 will not be used
 shell_type(k+1:k0) = 15
 i = 1
 do while(shell_type(i) /= 15)
  if(shell_type(i) /= -1) then
   i = i + 1
   cycle
  end if
  shell_type(i) = 0
  allocate( temp1(i+1 : k0-1), temp2(i+1 : k0-1) )
  temp1(i+1 : k0-1) = shell_type(i+1 : k0-1)
  shell_type(i+2 : k0) = temp1(i+1 : k0-1)
  temp2(i+1 : k0-1) = shell_to_atom_map(i+1 : k0-1)
  shell_to_atom_map(i+2 : k0) = temp2(i+1 : k0-1)
  deallocate(temp1, temp2)
  shell_type(i+1) = 1
  shell_to_atom_map(i+1) = shell_to_atom_map(i)
  i = i + 2
 end do

 length = i - 1
 shell_type(i : k0) = 0
 return
end subroutine split_L_func

! sort the shell_type, shell_to_atom_map by ascending order
! MOs will be adjusted accordingly
subroutine sort_shell_and_mo(ilen, shell_type, shell_to_atom_map, nbf, nif, coeff2)
 implicit none
 integer i, j, k
 integer ibegin, iend, natom
 integer jbegin, jend
 integer,parameter :: ntype = 10
 integer num0(ntype), num1(ntype), num(ntype)
 data num0 /0, 1, -2, 2, -3, 3, -4, 4, -5, 5/
 data num1 /1, 3, 5, 6, 7, 10, 9, 15, 11, 21/
 !          S  P  5D 6D 7F 10F 9G 15G 11H 21H

 integer,intent(in) :: ilen, nbf, nif
 integer,intent(inout) :: shell_type(ilen), shell_to_atom_map(ilen)
 integer,allocatable :: ith(:), ith_bas(:), tmp_type(:)
 real(kind=8),intent(inout) :: coeff2(nbf,nif)

 ! find the number of atoms
 natom = shell_to_atom_map(ilen)

 allocate(ith(0:natom), ith_bas(0:natom))
 ith = 0
 ith_bas = 0

 ! find the end position of each atom in array shell_to_atom_map
 do i = 1, natom, 1
  ith(i) = count(shell_to_atom_map==i) + ith(i-1)
 end do

 ! find the end position of basis functions between two atoms
 do i = 1, natom, 1
  ibegin = ith(i-1) + 1
  iend = ith(i)
  allocate(tmp_type(ibegin:iend))
  tmp_type = 0
  tmp_type = shell_type(ibegin:iend)
  num = 0
  do j = 1, ntype, 1
   num(j) = count(tmp_type == num0(j))
  end do
  k = 0
  do j = 1, ntype, 1
   k = k + num(j)*num1(j)
  end do
  ith_bas(i) = ith_bas(i-1) + k
  deallocate(tmp_type)
 end do

 ! adjust the MOs in each atom
 do i = 1, natom, 1
  ibegin = ith(i-1) + 1
  iend = ith(i)
  jbegin = ith_bas(i-1) + 1
  jend = ith_bas(i)
  call sort_shell_and_mo_in_each_atom(iend-ibegin+1, shell_type(ibegin:iend), &
  & jend-jbegin+1, nif, coeff2(jbegin:jend,1:nif))
 end do
 ! adjust the MOs in each atom done
 deallocate(ith, ith_bas)
 return
end subroutine sort_shell_and_mo

subroutine sort_shell_and_mo_in_each_atom(ilen1, shell_type, ilen2, nif, coeff2)
 implicit none
 integer i, tmp_type
 integer ibegin, iend, jbegin, jend
 integer,parameter :: ntype = 10
 integer,parameter :: num0(ntype) = (/0, 1, -2, 2, -3, 3, -4, 4, -5, 5/)
 integer,parameter :: num1(ntype) = (/1, 3, 5, 6, 7, 10, 9, 15, 11, 21/)
 !                                    S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer,parameter :: rnum(-5:5) = (/9, 7, 5, 3, 0, 1, 2, 4, 6, 8, 10/)

 integer,intent(in) :: nif, ilen1, ilen2
 integer,intent(inout) :: shell_type(ilen1)
 integer,allocatable :: ith_bas(:)
 real(kind=8),intent(inout) :: coeff2(ilen2,nif)
 real(kind=8),allocatable :: tmp_coeff1(:,:), tmp_coeff2(:,:)
 logical sort_done

 tmp_type = 0

 ! find the end position of basis functions within an atom
 allocate(ith_bas(0:ilen1))
 ith_bas = 0
 do i = 1, ilen1, 1
   ith_bas(i) = ith_bas(i-1) + num1(rnum(shell_type(i)))
 end do

 sort_done = .false.
 do while(.not. sort_done)
  sort_done = .true.
  do i = 1, ilen1-1, 1
    if(shell_type(i) == 0) cycle
    if(ABS(shell_type(i+1)) >= ABS(shell_type(i))) cycle
    sort_done = .false.
    tmp_type = shell_type(i+1)
    shell_type(i+1) = shell_type(i)
    shell_type(i) = tmp_type
    ibegin = ith_bas(i-1) + 1
    iend = ith_bas(i)
    jbegin = ith_bas(i) + 1
    jend = ith_bas(i+1)
    allocate(tmp_coeff1(ibegin:iend,nif), tmp_coeff2(jbegin:jend,nif))
    tmp_coeff1 = 0.0d0
    tmp_coeff2 = 0.0d0
    tmp_coeff1 = coeff2(ibegin:iend,:)
    tmp_coeff2 = coeff2(jbegin:jend,:)
    ith_bas(i) = ibegin+jend-jbegin
    coeff2(ibegin: ith_bas(i),:) = tmp_coeff2
    ith_bas(i+1) = jend+iend-jbegin+1
    coeff2(ibegin+jend-jbegin+1: ith_bas(i+1),:) = tmp_coeff1
    deallocate(tmp_coeff1, tmp_coeff2)
  end do
 end do
 deallocate(ith_bas)
 return
end subroutine sort_shell_and_mo_in_each_atom

subroutine fchk2py_permute_5d(nif,coeff)
 implicit none
 integer i
 integer order(5)
 integer,intent(in) :: nif
 real(kind=8),intent(inout) :: coeff(5,nif)
 real(kind=8) coeff2(5,nif)
 data order /5, 3, 1, 2, 4/
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in PySCF
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d-2, d-1, d0 , d+1, d+2

 coeff2 = 0.0d0
 forall(i = 1:5)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2py_permute_5d

subroutine fchk2py_permute_6d(nif,coeff)
 implicit none
 integer i
 integer order(6)
 integer,intent(in) :: nif
 real(kind=8),intent(inout) :: coeff(6,nif)
 real(kind=8) coeff2(6,nif)
 data order /1, 4, 5, 2, 6, 3/
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in PySCF
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,XZ,YY,YZ,ZZ

 coeff2 = 0.0d0
 forall(i = 1:6)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2py_permute_6d

subroutine fchk2py_permute_7f(nif,coeff)
 implicit none
 integer i
 integer order(7)
 integer,intent(in) :: nif
 real(kind=8),intent(inout) :: coeff(7,nif)
 real(kind=8) coeff2(7,nif)
 data order /7, 5, 3, 1, 2, 4, 6/
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in PySCF
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f-3, f-2, f-1, f0 , f+1, f+2, f+3

 coeff2 = 0.0d0
 forall(i = 1:7)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2py_permute_7f

subroutine fchk2py_permute_10f(nif,coeff)
 implicit none
 integer i
 integer order(10)
 integer,intent(in) :: nif
 real(kind=8),intent(inout) :: coeff(10,nif)
 real(kind=8) coeff2(10,nif)
 data order /1, 5, 6, 4, 10, 7, 2, 9, 8, 3/
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in PySCF
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ

 coeff2 = 0.0d0
 forall(i = 1:10)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2py_permute_10f

subroutine fchk2py_permute_9g(nif,coeff)
 implicit none
 integer i
 integer order(9)
 integer,intent(in) :: nif
 real(kind=8),intent(inout) :: coeff(9,nif)
 real(kind=8) coeff2(9,nif)
 data order /9, 7, 5, 3, 1, 2, 4, 6, 8/
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in PySCF
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4

 coeff2 = 0.0d0
 forall(i = 1:9)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2py_permute_9g

subroutine fchk2py_permute_15g(nif,coeff)
 implicit none
 integer i
 integer,intent(in) :: nif
 real(kind=8),intent(inout) :: coeff(15,nif)
 real(kind=8) coeff2(15,nif)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in PySCF
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz

 coeff2 = 0.0d0
 forall(i = 1:15)
  coeff2(i,:) = coeff(16-i,:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2py_permute_15g

subroutine fchk2py_permute_11h(nif,coeff)
 implicit none
 integer i
 integer order(11)
 integer,intent(in) :: nif
 real(kind=8),intent(inout) :: coeff(11,nif)
 real(kind=8) coeff2(11,nif)
 data order /11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10/
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in PySCF
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5

 coeff2 = 0.0d0
 forall(i = 1:11)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2py_permute_11h

subroutine fchk2py_permute_21h(nif,coeff)
 implicit none
 integer i
 integer,intent(in) :: nif
 real(kind=8),intent(inout) :: coeff(21,nif)
 real(kind=8) coeff2(21,nif)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in PySCF
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz

 coeff2 = 0.0d0
 forall(i = 1:21)
  coeff2(i,:) = coeff(22-i,:)
 end forall
 coeff = coeff2
 return
end subroutine fchk2py_permute_21h

