module gen_dat
   save
   integer,parameter :: lchar=1024
   !integer, parameter :: iout=0
   integer, parameter :: iout=17
   !integer :: nmodes_norm
   integer :: icoord_inp
   integer, parameter :: ifit_auto=0      ! TODO: Make it as input variable
   !character(len=lchar) :: cpotfile
   character(len=:), allocatable :: cpotfile
   character(len=:), allocatable :: input
   character(len=:), allocatable :: outstem
   character(len=:), allocatable :: outname
   integer :: icomp

   ! General Informations about system
   contains
     ! subroutine init_general_stuff
     !    implicit none
     !    nmodes_norm = n3_s - 6  !So far only for non linear molecules
     !    allocate(refcoord(n3_s))
     !    
     !    write(ci_form,'(A)')'(A50,1x,A1,1x,i10)'
     ! end subroutine

! ----------------------------------------------------------------------

      !subroutine dealloc_general_stuff
      !   implicit none
      !   deallocate(refcoord)
      !end subroutine
      
! ----------------------------------------------------------------------

      subroutine open_outfile
         implicit none
         open(iout,file=outstem//'.out',action='write')
      end subroutine

! ----------------------------------------------------------------------

      subroutine close_outfile
         implicit none
         close(iout)
      end subroutine

! ----------------------------------------------------------------------

      subroutine flushout()
         implicit none
         call flush(iout)
      end subroutine
end module

module time_management
   implicit none
   contains
      subroutine get_date(year,month,day,hour,minute,second)
         implicit none
         integer, dimension(8) :: date_time
         character(len=4) :: year
         character(len=3) :: month
         character(len=2) :: day,hour,minute,second
         character(len=10),dimension(3) :: dummy
         call date_and_time(dummy(1),dummy(2),dummy(3),date_time)
         call convert_i4(date_time(1),year)
         call get_month(date_time(2),month)
         call convert_i2(date_time(3),day)
         call convert_i2(date_time(5),hour)
         call convert_i2(date_time(6),minute)
         call convert_i2(date_time(7),second)
      end subroutine

! ----------------------------------------------------------------------

      double precision function wallcl() result(cl)
         integer :: c, cr
         call system_clock(count=c, count_rate=cr)
         cl = dble(c)/dble(cr)
      end function wallcl

! ----------------------------------------------------------------------

      subroutine get_month(im,cm)
         integer :: im
         character(len=3) :: cm
         if(im.eq.1) then
            cm = 'Jan'
         elseif(im.eq.2) then
            cm = 'Feb'
         elseif(im.eq.3) then
            cm = 'Mar'
         elseif(im.eq.4) then
            cm = 'Apr'
         elseif(im.eq.5) then
            cm = 'May'
         elseif(im.eq.6) then
            cm = 'Jun'
         elseif(im.eq.7) then
            cm = 'Jul'
         elseif(im.eq.8) then
            cm = 'Aug'
         elseif(im.eq.9) then
            cm = 'Sep'
         elseif(im.eq.10) then
            cm = 'Oct'
         elseif(im.eq.11) then
            cm = 'Nov'
         elseif(im.eq.12) then
            cm = 'Dec'
         endif
      end subroutine

! ----------------------------------------------------------------------

      subroutine convert_i2(i2,c2)
         implicit none
         integer :: i2
         character(len=2) :: c2
         if(i2.ge.10) then
            write(c2,'(i2)')i2
         else
            write(c2,'(a1,i1)')'0',i2
         endif
      end subroutine

! ----------------------------------------------------------------------

      subroutine convert_i4(i4,c4)
         implicit none
         integer :: i4
         character(len=4) :: c4
         if(i4.ge.1000) then
            write(c4,'(i4)')i4
         elseif(i4.ge.100.AND.i4.lt.1000) then
            write(c4,'(a1,i3)')'0',i4
         elseif(i4.ge.10.AND.i4.lt.100) then
            write(c4,'(a2,i2)')'00',i4
         else
            write(c4,'(a4,i1)')'000',i4
         endif
      end subroutine

end module time_management
