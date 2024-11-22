
! ---- HEADER ----------------------------------------------------------

subroutine write_out_header
   use gen_dat
   use time_management
   implicit none
   integer, dimension(8) :: date_time
   character(len=3) :: month
   character(len=2) :: day,hour,minute,second
   character(len=4) :: year
   character(len=10),dimension(3) :: dummy
   call get_date(year,month,day,hour,minute,second)
   write(iout,'(/,1x,38x,a24)')'*** MOLPRO INTERFACE ***'
   write(iout,'(1x,31x,a38,/)')'Reading and Fitting of Molpro Surfaces'
   call print_line(iout)
   write(iout,'(1x,a5,90x,a5)')'Date:','Time:'
   write(iout,'(1x,a3,a1,a2,a1,a4,81x,a2,a1,a2,a1,a2)')month,'-',day,'-',year,&
                                                hour,':',minute,':',second
   call print_line(iout)
end subroutine

! ----------------------------------------------------------------------

subroutine print_end
   use gen_dat
   implicit none
   write(iout,'(1x,2x,a18)')'Interface finished'
end subroutine

! ----------------------------------------------------------------------

subroutine print_line(iout)
   implicit none
   integer :: iout
   write(iout,'(1x,a100)')'************************************************************&
                                       **************************************************'
end subroutine

! ---- SURFACE ---------------------------------------------------------

subroutine write_head_surf
   use gen_dat
   implicit none
   write(iout,'(/,1x,2x,a43,/)')'STEP 1 * GRID REPRESENTATION OF THE SURFACE'
end subroutine

! ----------------------------------------------------------------------

subroutine write_end_surf(t1,t2)
   use gen_dat
   implicit none
   double precision :: t1,t2
   write(iout,'(/,1x,2x,a29,f9.2)')'Time for Reading the Surface:',t2-t1
   write(iout,'(/,1x,2x,a44)')'Finished: Grid Representation of the Surface'
   call print_line(iout)
end subroutine

! ---- POLY ------------------------------------------------------------

subroutine print_poly_headr
   use gen_dat
   implicit none
   write(iout,'(/,1x,2x,a28,/)')'STEP 2 * FITTING THE SURFACE'
end subroutine

! ----------------------------------------------------------------------

subroutine print_poly_information
   use gen_dat, only: iout
   use poly
   implicit none
   integer :: idim,isurf,itl,icheck
   double precision :: chi_m,chi_a
   double precision, dimension(:), pointer :: chi_all
   double precision, dimension(:,:), pointer :: chi_type
   integer, dimension(:,:), pointer :: modes
   character(len=30) :: text30
   character(len=7) :: text7,text7_2
   character(len=5) :: ct
   write(iout,'(/,1x,5x,a29)')'Information about the fitting'
   icheck = 0
   modes => p_modes(1)%p
   do isurf=1,ivnsurf_p(1)
      if(ifit_fctp(modes(1,isurf),1).ne.ifit_fctp(modes(1,1),1)) then
         icheck = 1
         exit
      endif
   enddo
   if(icheck.eq.1) then
      write(text30,'(A)')'Type of basis function:'
      write(iout,'(1x,7x,A30)')text30
      write(iout,'(1x,9x,a6,4x,a6,5x,a7)')'Coord.','Energy','Dip/Pol'
      modes => p_modes(1)%p
      do isurf=1,ivnsurf_p(1)
         call sget_char(ifit_fctp(modes(1,isurf),1),text7)
         call sget_char(ifit_fctp(modes(1,isurf),2),text7_2)
         write(iout,'(1x,9x,i6,4x,a7,4x,a7)')modes(1,isurf),text7,text7_2
      enddo
   else
      write(text30,'(A)')'Type of basis function:'
      call sget_char(ifit_fctp(modes(1,1),1),text7)
      write(iout,'(1x,7x,A30,a1,a7)')text30,'=',text7
   endif
   write(text30,'(A)')'Number of basis functions'
   write(iout,'(1x,7x,A30,a1,i3)')text30,'=',numpol_p
   write(iout,'(/,1x,5x,8x,a4,4x,a7,4x,a8)')'Dim.','Chi Av.','Chi Max'
   do idim=1,ncoup_p
      chi_all => p_chi(idim)%p
      chi_a = 0.0d0
      chi_m = 0.0d0
      do isurf=1,ivnsurf_p(idim)
         chi_a = chi_a + chi_all(isurf)
         if(chi_all(isurf).gt.chi_m) chi_m = chi_all(isurf)
      enddo
      chi_a = chi_a / dble(ivnsurf_p(idim))
      if(idim.eq.1) then
         write(iout,'(/,1x,5x,a6,2x,i2,a1,2x,es10.2,2x,es10.2)')'Energy',idim,'D',chi_a,chi_m
      else
         write(iout,'(1x,5x,6x,2x,i2,a1,2x,es10.2,2x,es10.2)')idim,'D',chi_a,chi_m
      endif
   enddo
   do idim=1,ncoupdip_p
      chi_type => p_chi_dip(idim)%p
      chi_a = 0.0d0
      chi_m = 0.0d0
      do isurf=1,ivnsurf_p(idim)
         do itl=1,3
            chi_a = chi_a + chi_type(itl,isurf)
            if(chi_type(itl,isurf).gt.chi_m) chi_m = chi_type(itl,isurf)
         enddo
      enddo
      chi_a = chi_a / dble(3*ivnsurf_p(idim))
      if(idim.eq.1) then
         write(iout,'(/,1x,5x,a6,2x,i2,a1,2x,es10.2,2x,es10.2)')'Dipole',idim,'D',chi_a,chi_m
      else
         write(iout,'(1x,5x,6x,2x,i2,a1,2x,es10.2,2x,es10.2)')idim,'D',chi_a,chi_m
      endif
   enddo
   do idim=1,ncouppol_p
      chi_type => p_chi_pol(idim)%p
      chi_a = 0.0d0
      chi_m = 0.0d0
      do isurf=1,ivnsurf_p(idim)
         do itl=1,6
            chi_a = chi_a + chi_type(itl,isurf)
            if(chi_type(itl,isurf).gt.chi_m) chi_m = chi_type(itl,isurf)
         enddo
      enddo
      chi_a = chi_a / dble(6*ivnsurf_p(idim))
      if(idim.eq.1) then
         write(iout,'(/,1x,5x,a6,2x,i2,a1,2x,es10.2,2x,es10.2)')'Polar',idim,'D',chi_a,chi_m
      else
         write(iout,'(1x,5x,6x,2x,i2,a1,2x,es10.2,2x,es10.2)')idim,'D',chi_a,chi_m
      endif
   enddo

end subroutine

! ----------------------------------------------------------------------

subroutine sget_char(idc,text7)
   implicit none
   integer :: idc
   character(len=7) :: text7
   if(idc.eq.1) then
      text7='Gauss'
   elseif(idc.eq.2) then
      text7='Bspline'
   elseif(idc.eq.3) then
      text7='Poly'
   elseif(idc.eq.4) then
      text7='Morse'
   elseif(idc.eq.5) then
      text7='Trigo'
   endif
end subroutine

! ----------------------------------------------------------------------

subroutine print_poly_end(t1,t2)
   use gen_dat
   implicit none
   double precision :: t1,t2
   write(iout,'(/,1x,2x,a29,f9.2)')'Time for Fitting the Surface:',t2-t1
   write(iout,'(/,1x,2x,a29)')'Finished: Fitting the Surface'
   call print_line(iout)
end subroutine
