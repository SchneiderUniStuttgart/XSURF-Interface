program molpro_interface
   use gen_dat, only: icoord_inp
   use surf_dat
   use surf_read
   use intc_info
   use poly
   use time_management, only: wallcl
   implicit none
   integer :: num_args,ix,idim,ires_read,isurf,ifit_inp
   character(len=lchar), dimension(:), allocatable :: args
   double precision :: t1,t2
   double precision, dimension(:), allocatable :: x1,x2
   integer, dimension(:,:), pointer :: modes_1
   double precision, dimension(9) :: coord
   double precision, dimension(1) :: val
   !TODO: Write better input and better output file
   !Add Input options: Choosing fit functions, the way of using it

   num_args = command_argument_count()
   allocate(args(num_args))
   do ix=1,num_args
      call get_command_argument(ix,args(ix))
   enddo
   
   input=trim(adjustl(args(1)))
   call get_stems(args(1))
   deallocate(args)
   call read_input
   call open_outfile
   call write_out_header
   call write_head_surf
   t1 = wallcl()
   call open_potfile
 
   !! Initialise reading potfile
   call flushout()
   call read_from_potfile_HEAD
   if(icoord_inp.eq.1) call read_coord_info_from_input
   call print_basic_surface_information
   call init_surf_dat_basic
   call flushout()
   do idim=1,ncoup_s_pot
      call read_from_potfile_NECK(idim)
      call print_dim_specific_headr(idim)
      call init_dim_spec_stuff(idim,ires_all_s(idim))
      call read_from_potfile_BODY(idim)
      call fit_xdim_surface(idim,ires_all_s(idim))
      call dealloc_igrid_nabi(idim)
   enddo
   call close_potfile
   t2 = wallcl()
   call write_end_surf(t1,t2)
   call flushout()
 
   call init_poly_basic(ncoup_s_pot,ncoupdip_s_pot,ncouppol_s_pot)
   call surf_to_poly
   call init_poly_general
   call print_poly_headr
   call flushout()
   t1 = wallcl()
   call calc_coef_xpoly(ires_read)
   call print_poly_information
   t2 = wallcl()
   call release_surf_dat

   call init_poly
   call flushout()
   call get_outname_from_potfile
   call print_coeffs(outname)
   call print_poly_end(t1,t2)
   call release_poly_after_fit
   call release_itypeint_1_2
   call print_end
   call close_outfile
   deallocate(input)
   deallocate(outstem)
end program

! -----------------------------------------------------------------------------------------------------

subroutine read_input
   ! Reads the input file and evaluates the needed values
   use gen_dat
   use surf_dat
   !use print_fitvals, only: niter
   implicit none
   integer :: iunit_inp,ios,idummy
   character(len=lchar) :: line,cdum
   iunit_inp = 12
   open(iunit_inp,file=input)
   icoord_inp = 0
   do
      read(iunit_inp,'(a)',iostat=ios) line
      if(ios.ne.0) exit
      call to_upper_case(lchar,line)
      if(line(1:3).eq.'POT') then
         !cpotfile=''
         cdum = ''
         call get_string_name(lchar,line,cdum)
         cpotfile = cdum(1:len_trim(cdum))
      elseif(line(1:5).eq.'COORD') then
         icoord_inp = 1
      !elseif(line(1:4).eq.'NDIM') then
      !   call get_int_from_input(lchar,line,ncoup_s)
      !elseif(line(1:4).eq.'NDIP') then
      !   call get_int_from_input(lchar,line,ncoupdip_s)
      !elseif(line(1:4).eq.'NPOL') then
      !   call get_int_from_input(lchar,line,ncouppol_s)
      endif
   enddo
   close(iunit_inp)
end subroutine

! -----------------------------------------------------------------------------------------------------

subroutine read_coord_info_from_input
   ! Reads the input file and evaluates the needed values
   use gen_dat
   use surf_dat
   !use print_fitvals, only: niter
   implicit none
   integer :: iunit_inp,ios,idummy,i,icrd,intfit
   character(len=lchar) :: line,cdum,keyword,chr
   iunit_inp = 12
   open(iunit_inp,file=input)
   do
      read(iunit_inp,'(a)',iostat=ios) line
      if(ios.ne.0) exit
      call to_upper_case(lchar,line)
      if(line(1:5).eq.'COORD') then
         cdum = line
         keyword = 'COORD'
         call get_int_from_keyword(lchar,cdum,icrd,keyword)
         do i=1,lchar
            if(line(i:i).eq.',') then
               if(line(i+1:i+6).eq.'FITFCT') then
                  keyword = 'FITFCT'
                  cdum = line
                  call get_char_from_keyword(lchar,cdum,chr,keyword)
                  call convert_chr_2_fitfct(lchar,chr,intfit)
                  ifit_fcts(icrd,1) = intfit
               endif
            endif
         enddo
      endif
   enddo
   close(iunit_inp)
end subroutine

! -----------------------------------------------------------------------------------------------------

subroutine convert_chr_2_fitfct(lchar,chr,intfit)
   implicit none
   integer :: intfit,lchar
   character(len=lchar) :: chr
   if(trim(chr).eq.'GAUSS') then
      intfit = 1
   elseif(trim(chr).eq.'BSPLINE') then
      intfit = 2
   elseif(trim(chr).eq.'POLY') then
      intfit = 3
   elseif(trim(chr).eq.'MORSE') then
      intfit = 4
   elseif(trim(chr).eq.'TRIGO') then
      intfit = 5
   endif
end subroutine

! -----------------------------------------------------------------------------------------------------

subroutine get_char_from_keyword(lchar,line,chrout,key)
   use gen_dat, only:iout
   implicit none
   integer :: lchar,ilenk,i,j,jsave,k
   character(len=lchar) :: line,key,chrout
   ilenk = len_trim(key)
   do i=1,lchar
      if(line(i:i+ilenk-1).eq.trim(key)) then
         do j=i,lchar
            if(line(j:j).eq.'=') then
               line(j:j) = ' '
               jsave = j
               exit
            endif
            line(j:j) = ' '
         enddo
         j2loop: do j=jsave+1,lchar
            if(line(j:j).eq.',') then
               do k=j,lchar
                  line(k:k) = ' '
               enddo
               exit j2loop
            endif
         enddo j2loop
         read(line,*)chrout
         return
      endif
      line(i:i) = ' '
   enddo
   write(iout,'(a)')'Usage of the keyword is wrong. Keyword:'
   write(iout,'(a)')trim(key)
   write(iout,'(a)')'Make sure that multiple keywords are sepearted by a&
                     comma, and that every keywords includes an &
                     equation mark'
   stop
end subroutine

! -----------------------------------------------------------------------------------------------------

subroutine get_int_from_keyword(lchar,line,intout,key)
   use gen_dat,only:iout
   implicit none
   integer :: lchar,intout,ilenk,i,j,jsave,k
   character(len=lchar) :: line,key
   ilenk = len_trim(key)
   do i=1,lchar
      if(line(i:i+ilenk-1).eq.trim(key)) then
         do j=i,lchar
            if(line(j:j).eq.'=') then
               line(j:j) = ' '
               jsave = j
               exit
            endif
            line(j:j) = ' '
         enddo
         j2loop: do j=jsave+1,lchar
            if(line(j:j).eq.',') then
               do k=j,lchar
                  line(k:k) = ' '
               enddo
               exit j2loop
            endif
         enddo j2loop
         read(line,*)intout
         return
      endif
      line(i:i) = ' '
   enddo
   write(iout,'(a)')'Usage of the keyword is wrong. Keyword:'
   write(iout,'(a)')trim(key)
   write(iout,'(a)')'Make sure that multiple keywords are sepearted by a&
                     comma, and that every keywords includes an &
                     equation mark'
   stop
end subroutine

! -----------------------------------------------------------------------------------------------------

subroutine get_int_from_input(l,line,ival)
   ! gets the integer from an input-line
   implicit none
   integer :: l,ival,i
   character(len=l) :: line

   do i=1,l
      if(line(i:i).eq.'=') then
         line(i:i) = ' '
         exit
      endif
      line(i:i) = ' '
   enddo
   line = trim(adjustl(line))
   read(line,*)ival
end subroutine


! -----------------------------------------------------------------------------------------------------

subroutine get_string_name(l,line,cnam)
   ! Get the name placed in ''
   implicit none
   integer :: l,idx,idx2,icnt
   character(len=l) :: line,cnam

   icnt = 1
   loop1: do idx=1,l
      if(line(idx:idx).eq."'".or.line(idx:idx).eq.'"') then
         do idx2=idx+1,l
            if(line(idx2:idx2).eq."'".or.line(idx2:idx2).eq.'"') exit loop1
            cnam(icnt:icnt) = line(idx2:idx2)
            icnt = icnt + 1
         enddo
      endif
   enddo loop1

end subroutine

! -----------------------------------------------------------------------------------------------------

subroutine to_upper_case(l,str_in)
   ! Transform lower letters to upper letters 
   implicit none

   integer :: i,l,icnt
   character(len=l) :: str_in
   character(26) :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(26) :: low = 'abcdefghijklmnopqrstuvwxyz'
   integer, dimension(1) :: ic

   icnt = 0
   do i=1,l
      if(str_in(i:i).eq.'"'.or.str_in(i:i).eq."'") icnt = icnt + 1
      if(mod(icnt,2).eq.1)cycle
      ic = INDEX(low, str_in(i:i))
      if (ic(1) > 0) str_in(i:i) = cap(ic(1):ic(1))
   end do

end subroutine

! -----------------------------------------------------------------------------------------------------

subroutine get_stems(finp)
   use gen_dat
   implicit none
   integer :: i,icnt
   character(len=lchar) :: finp
   icnt = 0
   do i=1,lchar
      if(finp(i:i).eq.'.')exit
      icnt = icnt + 1
   enddo
   outstem = finp(1:icnt)
end subroutine

! -----------------------------------------------------------------------------------------------------

subroutine get_outname_from_potfile
   use gen_dat
   implicit none
   integer :: i,iuse
   do i=1,len_trim(cpotfile)
      if(cpotfile(i:i).eq.'.') then
         iuse = i-1
      endif
   enddo
   outname = cpotfile(1:iuse)//'.apot'
end subroutine
