! This is a program used to solve a (up to) quartic equation (only one variable x)
!  or to find the extrema of a quartic trigonometric function (only one variable theta)
! To solve a polynominal euqation, call general_quartic_solver(a,b,c,d,e,nroot,root)
! written by ZJX 20170328, modified at 20170505

module polyn_info
 implicit none
 real(kind=8),parameter :: zero = 1.0d-9
end module polyn_info

!program main
! implicit none
! integer i,k,fid1,fid2,nroot
! integer,external :: freeunit
! real(kind=8) theta,min_val
! real(kind=8) a,b,c,d,e,root
! character(len=150) fname,outfname
! character(len=250) buffer
! logical fopened
!
! buffer = ' '
! fname = 'C2H4_sto-3g-old.out'
! k = index(fname,'.')
! outfname = fname(1:k-1)//'_root.txt'
! fid1 = freeunit()   ! find a usable file id
! open(unit=fid1,file=fname,position='rewind')
! fid2 = freeunit()   ! find a usable file id
! open(unit=fid2,file=TRIM(outfname),status='replace')
!
! do while(.true.)
!  read(fid1,'(A)') buffer
!  buffer = ADJUSTL(buffer)
!  if(buffer(1:2) == 'k=') exit
! end do
! do while(.true.)
!  read(unit=fid1,fmt=*,iostat=i) k,k,k,k,theta,a,b,c,d
!  if(i /= 0) exit
!  call quartic_trig_minima(a,b,c,d,nroot,root,min_val)
!  if(nroot == 0) then
!   write(fid2,'(A)') 'No root!!!'
!  end if
!  write(fid2,"(4I3,2X,6F12.8)") k,k,k,k,root,a,b,c,d,min_val
! end do
! close(fid1)
! close(fid2)
! stop
!end program main

! function freeunit: finds an unopened Fortran I/O
! unit and returns its numerical value from 1 to 9999; 
!function freeunit()
! implicit none
! integer ::  input=5,iout=6
! integer freeunit
! logical used
! freeunit = 9  ! try each logical unit until an unopened one is found
! used = .true.
! do while (used)
!  freeunit = freeunit + 1
!  if (freeunit .gt. 9999) then
!   write (iout,10)
!10 format (/,' FREEUNIT  --  No Available Fortran',' I/O Units')
!  end if
!  inquire (unit=freeunit,opened=used)
! end do
! return
!end function freeunit

subroutine quartic_trig_minima(d,c,b,a,nroot,root0,min_val)
 use polyn_info, only: zero
 implicit none
 integer i,k,loc_min_val
 integer,intent(out) :: nroot
 real(kind=8),intent(in) :: a,b,c,d
 real(kind=8) a_p,b_p,c_p,d_p,e_p
 real(kind=8) f(4),d2f_dtheta(4)
 real(kind=8),intent(out) :: min_val,root0
 real(kind=8) root(4)
 logical minima(4)
 
 root = 0.0d0   ! initialization
 root0 = 0.0d0
 min_val = 0.0d0
 f(4) = 0.0d0
 d2f_dtheta = 0.0d0
 loc_min_val = 0
 ! coefficient transformation
 a_p = b
 b_p = 2.0d0*(c - 2.0d0*a)
 c_p = 3.0d0*(d - b)
 d_p = -2.0d0*c
 e_p = -d
 ! transformation done
 call general_quartic_solver(a_p,b_p,c_p,d_p,e_p,nroot,root)
 if(nroot == 0) return
 minima = .false.
 k = 0
 do i = 1,nroot,1
  root(i) = DATAN(root(i))
  f(i) = a*(DSIN(root(i))**4) + b*(DSIN(root(i))**3)*DCOS(root(i)) + &
         c*(DSIN(root(i))**2)*(DCOS(root(i))**2) + d*DSIN(root(i))*(DCOS(root(i))**3)
  d2f_dtheta(i) = 2.0d0*(b-d)*DSIN(4.0d0*root(i)) + 2.0d0*(c-a)*DCOS(4.0d0*root(i))- &
                  (b+d)*DSIN(2.0d0*root(i)) + 2.0d0*a*DCOS(2.0d0*root(i))
  if(d2f_dtheta(i) > zero) then
   minima(i) = .true.
   if(k == 0) then
    min_val = f(i)
    loc_min_val = i
    k = k + 1
   else
    if(f(i) < min_val) then
     loc_min_val = i ; min_val = f(i)
    end if
   end if
  end if
 end do
 if(ALL(minima .eqv. .false.)) then
  nroot = 0
  root = 0.0d0
  return
 end if
 nroot = 1
 root(1) = root(loc_min_val)
 root(2:4) = 0.0d0
 root0 = root(1)
 return
end subroutine quartic_trig_minima

subroutine general_quartic_solver(a,b,c,d,e,nroot,root)
 use polyn_info, only: zero
 implicit none
 integer i,j
 integer,intent(out) :: nroot
 integer nroot_c,nroot_q1,nroot_q2
 real(kind=8),intent(in) :: a,b,c,d,e
 real(kind=8),intent(out) :: root(4)
 real(kind=8) b_a,c_a,d_a,e_a
 real(kind=8) a_c,b_c,c_c,d_c,root_c(3),y
 real(kind=8) a_p,b_p,c_p,a_q1,b_q1,b_q2,c_q1,c_q2,root_q1(2),root_q2(2)
 real(kind=8) temp_value(4),tmpv
 logical alive(4)

 root = 0.0d0   ! initialization
 alive = .false.
 if(a<zero .and. a>-zero) alive(1) = .true.
 if(alive(1)) then   ! a = 0, degrade to a cubic equation
  call general_cubic_solver(b,c,d,e,nroot,root(1:3))
  return
 end if
 if(b<zero .and. b>-zero) alive(2) = .true.
 if(d<zero .and. d>-zero) alive(3) = .true.
 if(alive(2) .and. alive(3)) then ! b = d = 0, degrade to a biquadratic equation
  call general_biquadratic_solver(a,c,e,nroot,root)
  return
 end if
 ! coefficient transformation and solve a cubic equation
 b_a = b/a ; c_a = c/a ; d_a = d/a ; e_a = e/a
 a_c = 1.0d0
 b_c = -c_a
 c_c = b_a*d_a - 4.0d0*e_a
 d_c = 4.0d0*c_a*e_a - d_a*d_a - b_a*b_a*e_a
 call general_cubic_solver(a_c,b_c,c_c,d_c,nroot_c,root_c)
 y = root_c(1)
 ! done transform and solve cubic
 ! coefficient transformation and solve two quadratic equations
 a_p = 0.25d0*b_a*b_a - c_a + y
 b_p = 0.5d0*b_a*y - d_a
 alive = .false.
 if(a_p<zero .and. a_p>-zero) alive(1) = .true.
 if(b_p<zero .and. b_p>-zero) alive(2) = .true.
 if(alive(1) .and. alive(2)) then ! a_p = b_p = 0
  c_p = 0.25d0*y*y-e
  if(c_p < -zero) then
   nroot = 0
   return
  end if
  nroot = 4
  tmpv = b_a*b_a - 8.0d0*(y + 2.0d0*c_p)
  if(tmpv < 0.0d0) then
   root(1) = 0.0d0
  else
   root(1) = DSQRT(tmpv)
  end if
  root(2) = -root(1)
  tmpv = b_a*b_a - 8.0d0*(y - 2.0d0*c_p)
  if(tmpv < 0.0d0) then
   root(3) = 0.0d0
  else
   root(3) = DSQRT(tmpv)
  end if
  root(4) = -root(3)
  root = 0.25d0*(root - b_a)
  goto 100
 end if
 ! a_p /= 0 and b_p /= 0
 a_q1 = 1.0d0
 temp_value(1) = DSQRT(DABS(a_p))
 temp_value(2) = 0.5d0*b_p*temp_value(1)/a_p
 temp_value(3) = 0.5d0*b_a
 temp_value(4) = 0.5d0*y
 b_q1 = temp_value(3) - temp_value(1)
 b_q2 = temp_value(3) + temp_value(1)
 c_q1 = temp_value(4) - temp_value(2)
 c_q2 = temp_value(4) + temp_value(2)
 call general_quadratic_solver(a_q1,b_q1,c_q1,nroot_q1,root_q1)
 call general_quadratic_solver(a_q1,b_q2,c_q2,nroot_q2,root_q2)
 ! done transform and solve quadratic
 ! combine all roots
 nroot = nroot_q1 + nroot_q2
 if(nroot_q1 == 2) then
  root(1:2) = root_q1
  root(3:4) = root_q2
 else if(nroot_q1 == 1) then
  root(1:1) = root_q1(1:1)
  root(2:3) = root_q2
 else
  root(1:2) = root_q2
 end if
 ! delete identical roots
100 if(nroot > 1) then
  alive = .true.
  do i = 1,nroot,1
   if(.not. alive(i)) cycle
   do j = i+1,nroot,1
    if(.not. alive(j)) cycle
    if(DABS(root(j) - root(i)) < zero) then
     root(j) = 0.0d0
     alive(j) = .false.
     nroot = nroot - 1
    end if
   end do
  end do
  temp_value = 0.0d0
  j = 1
  do i = 1,nroot
   if(.not. alive(i)) cycle
   temp_value(j) = root(i)
   j = j + 1
  end do
  root = temp_value
 end if   ! delete done
 return
end subroutine general_quartic_solver

subroutine general_cubic_solver(a,b,c,d,nroot,root)
 use polyn_info, only: zero
 implicit none
 integer,intent(out) :: nroot
 real(kind=8),intent(in) :: a,b,c,d
 real(kind=8),intent(out) :: root(3)
 real(kind=8) p,q,a2,b2,delta,sqrt_delta
 real(kind=8) theta,cos_theta,temp_value
 real(kind=8) tmpv1,tmpv2
 real(kind=8),parameter :: PI = 3.1415926535898d0

 root = 0.0d0   ! initialization
 if(a<zero .and. a>-zero) then   ! a = 0, degrade to a quadratic equation
  call general_quadratic_solver(b,c,d,nroot,root(1:2))
  return
 end if
 ! coefficient transformation
 a2 = a**2
 b2 = b**2
 temp_value = b/(3.0d0*a)
 p = (3.0d0*a*c - b2)/(3.0d0*a2)
 q = (2.0d0*b2*b - 9.0d0*a*b*c + 27.0d0*a2*d)/(27.0d0*a2*a)
 ! transformation done
 delta = 0.25d0*(q*q) + (p**3)/27.0d0
 if(delta > zero) then
  nroot = 1
  sqrt_delta = DSQRT(delta)
  tmpv1 = - 0.5d0*q + sqrt_delta
  tmpv2 = - 0.5d0*q - sqrt_delta
  call dcurt(tmpv1)
  call dcurt(tmpv2)
  root(1) = tmpv1 + tmpv2
  root(1) = root(1) - temp_value
 else if(delta < -zero) then
  if(p > 0.0d0) p = 0.0d0
  nroot = 3
  cos_theta = -0.5d0*q*DSQRT(-27.0d0*p)/(p*p)
  cos_theta = min(max(-1.0d0,cos_theta),1.0d0) ! in case of value >1 or <-1
  theta = DACOS(cos_theta)/3
  q = 2.0d0*DSQRT(-p/3.0d0)
  cos_theta = 2.0d0*PI/3.0d0
  root(1) = q*DCOS(theta) - temp_value
  root(2) = q*DCOS(theta + cos_theta) - temp_value
  root(3) = q*DCOS(theta - cos_theta) - temp_value
 else
  nroot = 2
  tmpv1 = 0.5d0*q
  call dcurt(tmpv1)
  root(2) = tmpv1
  !root(2) = DSQRT(-p/3.0d0) ! this sometimes cause errors
  root(1) = -2.0d0*root(2) - temp_value
  root(2) = root(2) - temp_value
  ! delete identical roots
  if(DABS(root(1) - root(2)) <zero) then
   nroot = 1
   root(2) = 0.0d0
  end if
 end if
 return
end subroutine general_cubic_solver

subroutine general_biquadratic_solver(a,b,c,nroot,root)
 use polyn_info, only: zero
 implicit none
 integer,intent(out) :: nroot
 real(kind=8),intent(in) :: a,b,c
 real(kind=8),intent(out) :: root(4)
 real(kind=8) d,delta

 root = 0.0d0   ! initialization
 if(a<zero .and. a>-zero) then   ! a = 0, degrade to a quadratic equation
  d = 0.0d0
  call general_quadratic_solver(b,d,c,nroot,root(1:2))
  return
 end if
 call general_quadratic_solver(a,b,c,nroot,root(1:2))
 if(nroot == 2) then
  if(root(1) > zero) then
   nroot = 3
   root(3) = root(2)
   root(1) = DSQRT(root(1))
   root(2) = -root(1)
   if(root(3) > zero) then
    nroot = 4
    root(3) = DSQRT(root(3))
    root(4) = -root(3)
   else if(root(3) < -zero) then
    nroot = 2
    root(3) = 0.0d0
   else
    nroot = 3
   end if
  else if(root(1) < -zero) then
   nroot = 1
   root(1) = root(2)
   root(2) = 0.0d0
   if(root(1) > zero) then
    nroot = 2
    root(1) = DSQRT(root(1))
    root(2) = -root(1)
   else if(root(1) < -zero) then
    nroot = 0
    root(1) = 0.0d0
   else
    nroot = 1
   end if
  else
   if(root(2) > zero) then
    nroot = 3
    root(2) = DSQRT(root(2))
    root(3) = -root(2)
   else if(root(1) < -zero) then
    nroot = 1
    root(2) = 0.0d0
   else
    nroot = 1
   end if
  end if
 else if(nroot == 1) then
  if(root(1) > zero) then
   nroot = 2
   root(1) = DSQRT(root(1))
   root(2) = -root(1)
  else if(root(1) < -zero) then
   nroot = 0
   root = 0.0d0
  end if
 end if
 return
end subroutine general_biquadratic_solver

subroutine general_quadratic_solver(a,b,c,nroot,root)
 use polyn_info, only: zero
 implicit none
 integer,intent(out) :: nroot
 real(kind=8),intent(in) :: a,b,c
 real(kind=8),intent(out) :: root(2)
 real(kind=8) delta
 
 root = 0.0d0   ! initialization
 if(a<zero .and. a>-zero) then   ! a = 0
  if(b<zero .and. b>-zero) then   ! b = 0
   if(c<zero .and. c>-zero) then   ! c = 0
    nroot = 1
   else   ! c /= 0
    nroot = 0
    return
   end if
  else   ! b /= 0
   nroot = 1
   root(1) = -c/b
  end if
 else   ! a /= 0
  delta = b*b - 4.0d0*a*c
  if(delta > zero) then
   nroot = 2
   delta = DSQRT(delta)
   root(1) = 0.5d0*(-b + delta)/a
   root(2) = 0.5d0*(-b - delta)/a
  else if(delta < -zero) then
   nroot = 0
   return
  else
   nroot = 1
   root(1) = -0.5d0*b/a
  end if
 end if
 return
end subroutine general_quadratic_solver

! double precision for cubic root
! In expression 'a**(1.0d0/3.0d0)', if a<0.0d0, it may cause errors
! in gnu compiler, so call subroutine dcurt is safer
subroutine dcurt(num)
 implicit none
 real(kind=8),intent(inout) :: num
 if(num > 0.0d0) then
  num = num**(1.0d0/3.0d0)
 else
  num = (-num)**(1.0d0/3.0d0)
  num = -num
 end if
 return
end subroutine dcurt

