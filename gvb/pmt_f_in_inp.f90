! written by jxzou at 20171109: adjust the orders of 10 Cartesian f functions in .inp file

! f functions of MOs in the .inp file will be replaced by their new order, which corresponds
! to the order in Gamess
! the order of Cartesian f functions in Gaussian: XXX,YYY,ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ
! the order of Cartesian f functions in Gamess:   XXX,YYY,ZZZ, XXY, XXZ, XYY, YYZ, XZZ, YZZ, XYZ
! The order of Cartesian f functions can be known by adding keyword 'pop=full' in .gjf file.
! In Gamess output files(.out, .gms), the order of Cartesian f functions is printed without
!  extra keyword needed.
! Limitation1: only deal with MOs with no linear dependence and Cartesian basis functions
! Limitation2: the number of atoms which have f functions must be less than 1000
! Limitation3: nbf must be less than 1000

program main
 implicit none
 character(len=240) old_inp

 call getarg(1,old_inp)
 call pmt_f_in_inp(old_inp)
 stop
end program main

subroutine pmt_f_in_inp(old_inp)
 implicit none
 integer i,j,k,m
 integer RENAME
 integer ncoeff,nbf,nif,nmark
 ! nbf: number of basis functions
 ! nif: number of independent functions
 ! ham: highest angular momentum
 ! natom: the total number of atoms
 integer nline,remain
 integer f_mark(999) ! mark the index where f functions begin
 integer,parameter :: fid1 = 11, fid2 = 12
 real(kind=8),allocatable :: coeff(:,:),coeff2(:)
 character(len=5) str1
 character(len=30) str2
 character(len=240),intent(in) :: old_inp
 character(len=240) new_inp
 character(len=240) buffer

 str1 = ' '   ! initialization
 str2 = ' '
 buffer = ' '
 ncoeff = 0
 nmark = 0
 f_mark = 0
 nline = 0
 remain = 0
 k = index(old_inp,'.inp')
 if(k == 0) then
  write(*,'(A)') "ERROR in subroutine pmt_f_in_inp: the input file does not end with '.inp'!"
  return
 else
  new_inp = old_inp(1:k-1)//'_new.inp'
 end if

 open(unit=fid1,file=TRIM(old_inp),status='old',position='rewind')
 open(unit=fid2,file=TRIM(new_inp),status='replace')
 ! find natom and how many basis functions for each atom
 do while(.true.)
  read(fid1,'(A)',iostat=i) buffer
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buffer)
  k = index(buffer,'$data')
  if(k == 0) k = index(buffer,'$DATA')
  if(k /= 0) exit
 end do
 if(i /= 0) then
  write(*,'(A)') "ERROR in subroutine pmt_f_in_inp: no '$DATA' section in the input file:"
  write(*,'(A)') TRIM(old_inp)
  close(fid1)
  close(fid2,status='delete')
  return
 end if

 do i = 1,2,1  ! skip 2 line
  read(fid1,'(A)') buffer
  write(fid2,'(A)') TRIM(buffer)
 end do

 nbf = 0
 do while(.true.)
  read(fid1,'(A)') buffer
  write(fid2,'(A)') TRIM(buffer)
  buffer = ADJUSTL(buffer)
  if(buffer(1:1)=='$' .or. buffer(1:1)==' ') exit ! meet '$END' or blank line
  do while(.true.)
   read(fid1,'(A)') buffer
   write(fid2,'(A)') TRIM(buffer)
   buffer = ADJUSTL(buffer)
   if(buffer(1:1)=='$' .or. buffer(1:1)==' ') exit ! meet '$END' or blank line
   if(buffer(1:1) == 'S') then
    nbf = nbf + 1
   else if(buffer(1:1) == 'P') then
    nbf = nbf + 3
   else if(buffer(1:1) == 'D') then
    nbf = nbf + 6
   else if(buffer(1:1) == 'F') then
    nmark = nmark + 1
    f_mark(nmark) = nbf + 1
    nbf = nbf + 10
   end if
   call get_int_after_flag(buffer,' ',nline)
   do i = 1,nline,1  ! skip n lines
    read(fid1,'(A)') buffer
    write(fid2,'(A)') TRIM(buffer)
   end do
  end do
 end do
 ! natom and atom_basis found

 do while(.true.)
  read(fid1,'(A)',iostat=i) buffer
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buffer)
  k = index(buffer,'$vec')
  if(k == 0) k = index(buffer,'$VEC')
  if(k /= 0) exit
 end do
 if(i /= 0) then
  write(*,'(A)') "No '$VEC' section in the input file, still continue."
  write(*,'(A)') TRIM(old_inp)
  close(fid1)
  close(fid2,status='delete')
  return
 end if

 nif = nbf
 allocate(coeff(nbf,nif))
 coeff = 0.0d0
 remain = MOD(nbf,5)
 if(remain > 0) then
  nline = (nbf - remain)/5
 else
  nline = nbf/5
 end if
 ! read MOs from .inp file
 do i = 1,nif,1
  k = 1
  do j = 1,nline,1
   read(fid1,'(A)') buffer
   buffer = buffer(6:)
   read(buffer,'(5ES15.8)') coeff(k:k+4,i)
   k = k + 5
  end do
  if(remain > 0) then
   read(fid1,'(A)') buffer
   buffer = buffer(6:)
   str1 = ' '
   write(str1,'(I5)') remain
   str1 = ADJUSTL(str1)
   str2 = '('//TRIM(str1)//'ES15.8)'
   read(buffer,TRIM(str2)) coeff(k:nbf,i)
  end if
 end do
 !write(*,*) 'nline=',nline
 !write(*,*) 'remain=',remain
 !write(*,*) 'ncoeff=',ncoeff
 !write(*,*) 'k=',k

 ! if nmark > 0, adjust the order of Cartesian functions
 if(nmark > 0) then
  do i = 1,nmark,1
   call permute_f2(nif,coeff(f_mark(i)+3:f_mark(i)+8,:))
  end do
 end if
 ! adjustment finished

 ncoeff = nbf*nif
 ! output the MOs to new_inp file
 do i = 1,nif,1
  m = i/100
  m = i - m*100
  k = 1
  do j = 1,nline,1
   write(fid2,'(I2,I3,5ES15.8)') m,j,coeff(k:k+4,i)
   k = k + 5
  end do
  if(remain > 0) then
   str1 = ' '
   write(str1,'(I5)') remain
   str1 = ADJUSTL(str1)
   str2 = '(I2,I3,'//TRIM(str1)//'ES15.8)'
   write(fid2,TRIM(str2)) m,j,coeff(k:nbf,i)
  end if
 end do
 deallocate(coeff)
 ! output MOs finished

 ! copy the rest of old_inp file
 do while(.true.)
  read(fid1,'(A)',iostat=i) buffer
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buffer)
 end do
 ! copy finished

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(new_inp),TRIM(old_inp))
 return
end subroutine pmt_f_in_inp

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

subroutine permute_f2(nif,coeff)
 implicit none
 integer i,nif
 real(kind=8),intent(inout) :: coeff(6,nif)
 real(kind=8) temp_coeff(nif)

 temp_coeff = 0.0d0
 ! From: the order of Cartesian f functions in Gaussian
 ! To: the order of Cartesian f functions in Gamess
 ! From: XXX,YYY,ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ,XYZ
 ! To:   XXX,YYY,ZZZ, XXY, XXZ, XYY, YYZ, XZZ, YZZ,XYZ
  temp_coeff = coeff(1,:) ! XYY,XXY,XXZ,XZZ,YZZ,YYZ
  coeff(1,:) = coeff(2,:)
  coeff(2,:) = temp_coeff ! XXY,XYY,XXZ,XZZ,YZZ,YYZ
  temp_coeff = coeff(2,:)
  coeff(2,:) = coeff(3,:)
  coeff(3,:) = temp_coeff ! XXY,XXZ,XYY,XZZ,YZZ,YYZ
  temp_coeff = coeff(4,:)
  coeff(4,:) = coeff(5,:)
  coeff(5,:) = temp_coeff ! XXY,XXZ,XYY,YZZ,XZZ,YYZ
  temp_coeff = coeff(4,:)
  coeff(4,:) = coeff(6,:)
  coeff(6,:) = temp_coeff ! XXY,XXZ,XYY,YYZ,XZZ,YZZ
 return
end subroutine permute_f2

