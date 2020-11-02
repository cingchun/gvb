! convert .fchk to .vec, which can be used for GAMESS $VEC section
! download from http://blog.sina.com.cn/s/blog_6594d5d30101jd9a.html
! this version has been modified to be more user-friendly by zjx
! the most recent modification: 20170209
Program Read_CS
        implicit none

!parameter
        integer natom
!      integer :: file_id
        character*150 :: file_in=' '
        character*150 :: file_out=' '
        character*8 :: flag='Alpha MO'

        integer :: error
        character*8 :: buffer

        double precision,allocatable :: NBO(:)
        character*10 :: short_string = ' '
        character*100 :: string
!count
        integer :: i,j,k

        call getarg(1,file_in)
        call getarg(2,short_string)
        read(short_string,*) natom
        k = index(file_in,'.fch')
        file_out = file_in(1:k-1)//'.vec'
        allocate (NBO(natom*natom),stat=error)
        if(error /= 0) then
          stop
        end if

        open(unit=101,file=file_in,status='OLD',iostat=error)
        if(error/=0) then
          stop
        end if

        do while(.TRUE.)
          read(101,"(A8)") buffer
          if(buffer==flag) then
            do i=1,natom**2
              read(101,"(E16.8)",advance='NO') NBO(i)
              if(mod(i,5)==0) then
                read(101,*)
              end if
            end do
            exit
          end if
        end do

        close(101)

        open(unit=102,file=file_out,status="REPLACE",iostat=error)
        if(error/=0) then
          stop
        end if

        write(102,"(A5)",advance='NO') ' $VEC'
        do i=1,natom
          do j=1,natom
            if(mod(j,5)==1) then
              write(102,*)
              k=i
              do while (k > 99)
                k=k-100
              end do
              write(102,"(I2,I3)",advance='NO'),k,ceiling(j/.5d1)
            end if
            write(102,"(E15.8)",advance='NO') NBO((i-1)*natom+j)
          end do
        end do

        write(102,*)
        write(102,"(A5)") ' $END'
        close(102)

        !open(unit=105,file="VEC.txt")
        !  do i=1,natom
        !    do j=1,natom
        !      write(105,*)i,j,NBO((i-1)*natom+j)
        !    end do
        !  end do
        !close(105)

        deallocate (NBO)

      End
