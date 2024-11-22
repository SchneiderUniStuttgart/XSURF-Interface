! Short documentation:
! Routines in mod_surf_extern:
!
! Routines for reading and writing to potfile and restartfile
!
!  1) read_from_potfile_HEAD: Reading the information out of the head of the potfile
!     Remark: some of the information is only for the user and is not read again
!  2) read_from_potfile_HEAD_mppx: interface for 3) for mppx mode
!  3) read_from_potfile_BODY: Reading surface information
!  4) read_List_restart_HEAD: reading the head of the restart file
!  5) read_List_restart_HEAD_mppx: interface for mppx mode
!  6) read_List_restart_BODY: reading body of restart file
!  7) read_List_restart_BODY_mppx: interface for mppx mode
!  8) close_potfile: closing the potfile
module surf_read
   use gen_dat, only: iout
   use surf_dat

   implicit none
   save
      integer :: iunitpot                                               ! Unit for potfile
   contains
!========================================================================
!=== Handling Potfile ===================================================
!========================================================================
      subroutine open_potfile
         implicit none
         write(iout,'(1x,5x,a32,2x,a)')'The potential is read from file:',trim(cpotfile)
         iunitpot = 13
         open(iunitpot,file=trim(adjustl(cpotfile)),action='read')
      end subroutine

!---------------------------------------------------------------------

      subroutine close_potfile
         implicit none
         close(iunitpot)
      end subroutine
!========================================================================
!=== Reading Potential =============================================
!========================================================================
      subroutine read_from_potfile_HEAD
         use intc_info
         use poly, only: ncart_p,pi_p
         implicit none

         integer :: ios,ios2
         integer :: i,j,icnt

         double precision :: d_intermed
         double precision, dimension(:),allocatable :: qval_primitives

         character(len=16384) :: line!,form
         character(1) :: text1
         character(16) :: text16

         read1loop: do
            read(iunitpot,'(A)',iostat=ios) line
            if(ios/=0) return
            ! General parameters like dimensino and stuff
            if(trim(line)=='### PROGRAM PARAMETERS')then
               read2loop: do
                  read(iunitpot,'(A)',iostat=ios) line
                  if(ios/=0) return
                  if(trim(line)=='') exit read2loop
                  if(trim(line(3:64+2))=='Order of the potential energy expansion:')then
                     read(line(2+64+16+1+1:2+64+16+1+4),'(I4)',iostat=ios2) ncoup_s_pot
                  elseif(trim(line(3:64+2))=='Order of the dipole surface expansion:')then
                     read(line(2+64+16+1+1:2+64+16+1+4),'(I4)',iostat=ios2) ncoupdip_s_pot
                  elseif(trim(line(3:64+2))=='Order of the polarizability surface expansion:')then
                     read(line(2+64+16+1+1:2+64+16+1+4),'(I4)',iostat=ios2) ncouppol_s_pot
                  elseif(trim(line(3:64+2))=='Type of coordinate:')then
                     read(line(2+64+16+1+1:2+64+16+1+4),'(I4)',iostat=ios2) itypeint
                  endif
               enddo read2loop
            ! refcoord
            elseif(trim(line)=='### GEOMETRY')then
               !refcoord=0.0d0
               icnt = 0
               do 
                  read(iunitpot,'(A)',iostat=ios) line
                  if(trim(line).eq.'')exit
                  icnt = icnt + 1
               enddo
               ncart_p = icnt * 3
            elseif(trim(line)=='### FIT FUNCTIONS')then
               call init_fitting
               read(iunitpot,'(A)')line
               do i=1,ncrd_s
                  read(iunitpot,'(A)')line
                  read(line,'(7x,I7,I8)')ifit_fcts(i,1),ifit_fcts(i,2)
               enddo
            !Internal Coordinates - Parameters
            elseif(trim(line)=='### PARAMETER INTERNAL COORDINATES') then
               !TODO: Add internal coordinate definitions!
               read(iunitpot,'(A)',iostat=ios) line
               read(line,'(4x,A16,A1,I6)')text16,text1,nstretch
               read(iunitpot,'(A)',iostat=ios) line
               read(line,'(4x,A16,A1,I6)')text16,text1,nbending
               read(iunitpot,'(A)',iostat=ios) line
               read(line,'(4x,A16,A1,I6)')text16,text1,ntorsion
               read(iunitpot,'(A)',iostat=ios) line
               call init_itypeint_1_2(ncrd_s)
            !Internal Coordinates - Definitions
            elseif(trim(line)=='### DEFINITION INTERNAL COORDINATES') then
               allocate(qval_primitives(ncrd_s))
               qval_primitives = 0.0d0
 
               read3loop: do
                  read(iunitpot,'(A)',iostat=ios) line
                  if(trim(adjustl(line))=='Stretching') then
                     do i=1,nstretch
                        read(iunitpot,'(A)',iostat=ios) line
                        read(line,'(4x,I5,I8,I8,I8,F18.12)')atstretch(3,i),atstretch(1,i),atstretch(2,i),&
                                                                    atstretch(4,i),d_intermed
                        qval_primitives(atstretch(3,i)) = d_intermed
                     enddo
                     !call get_connectivity
                  elseif(trim(adjustl(line))=='Bending') then
                     do i=1,nbending
                        read(iunitpot,'(A)',iostat=ios) line
                        read(line,'(4x,I5,I8,I8,I8,F18.12)')atbending(4,i),atbending(1,i),atbending(2,i),&
                                                               atbending(3,i),d_intermed
                        qval_primitives(atbending(4,i)) = d_intermed
                    enddo
                  elseif(trim(adjustl(line))=='Torsion') then
                     do i=1,ntorsion
                        read(iunitpot,'(A)',iostat=ios) line
                        read(line,'(4x,I5,I8,I8,I8,I8,F18.12)')attorsion(5,i),attorsion(1,i),attorsion(2,i),&
                                                          attorsion(3,i),attorsion(4,i),d_intermed
                        qval_primitives(attorsion(5,i)) = d_intermed
                     enddo
                     exit read3loop
                  endif
               enddo read3loop
               call determine_pf_intcs(pi_p)
               call allocate_refqvals(qval_primitives,ncrd_s)
               !call generate_refgeom(qval_primitives)
               deallocate(qval_primitives)
            elseif(trim(line)=='### REFERENCE')then
               exit ! last entry of HEAD
            endif
         enddo read1loop

      end subroutine

!---------------------------------------------------------------------

      subroutine read_from_potfile_NECK(idim)
         implicit none
         integer :: idim
         integer :: ios

         character(16384) :: line
         character(32) :: text32

         ! reading surfaces
         write(text32,'(A,I1,A)') '### ',idim,'D SURFACES'
         do
            read(iunitpot,'(A)',iostat=ios) line
            if(ios/=0) return
            if(trim(line)==trim(text32)) exit
         enddo
         !isthere=1
         ! reading surfaces info
         read(iunitpot,'(A)',iostat=ios) line ! ires
         read(line(45:52),'(I8)') ires_all_s(idim)
         read(iunitpot,'(A)',iostat=ios) line ! igrid_all(idim)
         read(line(45:52),'(I8)') igrid_all(idim)
         read(iunitpot,'(A)',iostat=ios) line ! ivnsurf(idim)
         read(line(45:52),'(I8)') ivnsurf_s(idim)
         call alloc_i_labxd_s(idim)
      end subroutine

!   !---------------------------------------------------------------------

      subroutine read_from_potfile_BODY(idim)
         implicit none
         integer :: idim
         integer :: i,ii,j,ios,jj,idummy
         integer, dimension(idim) :: modes

         character(16384) :: line, form

         write(form,'(A,I1,A,I2,A)') '(2x,',idim,'F16.10,',ires_all_s(idim),'F18.12)'
         read(iunitpot,'(A)',iostat=ios) line ! empty line
         ! ab initio points
         surfloop : do i=1,ivnsurf_s(idim)
            read(iunitpot,'(A)',iostat=ios) line
            read(line,'(22x,8I4)') (modes(ii),ii=1,idim)
            call iset_lab_grid_s(idim,i,modes)
            do
               read(iunitpot,'(A)',iostat=ios) line
               if(line(15:18)=='Grid')then
                  exit
               elseif(line(5:32)=='Number of ab initio points: ')then
                  read(line,'(4x,28x,I5,2x,I5,1x)') idummy,igrid_total(i)
               endif
            enddo
            read(iunitpot,'(A)',iostat=ios) line ! headline
            ! setting isurface_adr
            do j=i+1,ivnsurf_s(idim)+1
               isurface_adr(j)=isurface_adr(j)+igrid_total(i)
            enddo
            jj=1
            do
               read(iunitpot,'(A)',iostat=ios) line
               if(ios/=0) exit
               if(trim(line)=='') exit
                 read(line,form) (grid_abi_ptr(ii,jj+isurface_adr(i)),ii=1,idim+ires_all_s(idim))
               jj=jj+1
            enddo
            ! after reading grid, setting some grid information
            ! like the number of grid points for each direction
            do j=1,idim
               igrid_nabi(j,i)=1
               do ii=1,igrid_total(i)-1
                  if(grid_s(j,ii+isurface_adr(i))<grid_s(j,ii+1+isurface_adr(i)))then
                     if(abs(grid_s(j,ii+isurface_adr(i))+1.0d-6-grid_s(j,ii+1+isurface_adr(i)))>1.0d-6)then
                        igrid_nabi(j,i)=igrid_nabi(j,i)+1
                     endif
                  else
                     if(abs(grid_s(j,ii+(isurface_adr(i)))-grid_s(j,ii+1+isurface_adr(i)))>1.0d-6)then
                        exit
                     endif
                  endif
               enddo
            enddo
         enddo surfloop
         if(idim==1)then
            do i=1,ncrd_s
               xmin_s(i)=grid_s(1,1+isurface_adr(i))
               xmax_s(i)=grid_s(1,igrid_total(i)+isurface_adr(i))
            enddo
         endif
      end subroutine
end module surf_read
