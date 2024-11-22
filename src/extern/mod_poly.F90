! This module is the heart of the evaluation of the potential. But additionally you need the 
! files fittools.F90 and mod_intc_info.F90 in your source code, to make everything running.
! Keep in mind, that you also need LAPACK and BLAS routines for running the code.
! There are three Routines to use this module:
!
! 1. set_iout_p:
!     This Routine sets the general unitnumber for the error outputs of the code
! 2. get_surface_information:
!     With this routine, the surface is evaluated at a certain point
! 3. read_fit_coefs:
!     To read the file with the coefficients, this routine has to be called.



module poly
   implicit none


   save

!----------------------------------------------------------------------
   ! objects
   type iptr_array_1
      integer, dimension(:), pointer, contiguous :: p
   end type iptr_array_1

   type iptr_array_2
      integer, dimension(:,:), pointer, contiguous :: p
   end type iptr_array_2
!
   type ptr_array_1
      double precision, dimension(:), pointer, contiguous :: p
   end type ptr_array_1
!
   type ptr_array_2
      double precision, dimension(:,:), pointer, contiguous :: p
   end type ptr_array_2

   type ptr_array_3
      double precision, dimension(:,:,:), pointer, contiguous :: p
   end type ptr_array_3

!----------------------------------------------------------------------
   ! public variables
   ! some dimensions and number of surfaces
   double precision, parameter :: pi_p= 3.141592653589793238462643383279502884197169399375105820974944d0
   integer, parameter :: lch = 1024
   integer :: iout_p !-> TODO: Set this to the right value
   character(len=lch) :: ci_form
   character(len=lch) :: form2
   character(len=lch) :: formc
   integer :: ncart_p !Number of cartesian coordinates
   integer :: iunit_apot
   integer :: ncrd_p
   integer :: ngrid_max_p
   integer :: ncoup_p
   integer :: ncoupdip_p
   integer :: ncouppol_p
   integer :: ibas_min ! --> Make it as Input
   integer :: ibas_max ! --> Make it as Input
   integer :: numpol_p
   ! fit functions
   integer, dimension(:,:), allocatable :: ifit_fctp

   ! fine grid points
   double precision, dimension(:,:), allocatable :: fine_grid

   ! mode combinations
   type(iptr_array_2), dimension(:), allocatable :: p_modes

   ! coefficients and additional information
   type(ptr_array_2), dimension(:), allocatable :: p_coef_ene
   type(ptr_array_3), dimension(:), allocatable :: p_coef_dip
   type(ptr_array_3), dimension(:), allocatable :: p_coef_pol
   type(ptr_array_1), dimension(:), allocatable :: p_chi
   type(ptr_array_2), dimension(:), allocatable :: p_chi_dip
   type(ptr_array_2), dimension(:), allocatable :: p_chi_pol
   double precision, dimension(:), allocatable :: chi_max,chi_av

   integer,dimension(:), allocatable :: ires_all_p                      ! number of properties computed
   integer, dimension(:), allocatable :: ivnsurf_p                      ! number of surfaces per dimension
   double precision, dimension(:), allocatable :: xmin_p                ! Min value of elongation
   double precision, dimension(:), allocatable :: xmax_p                ! Max value of elongation
   contains
!========================================================================
!=== External Routines ==================================================
!========================================================================
! The following 3 Routines are needed to evaluate the potential surface
      subroutine set_iout_p(iout)
         ! Set the desired file, where the output should be printed
         ! Input:
         !     iout: unitnumber for output
         implicit none
         integer :: iout
         iout_p = iout
      end subroutine
      
! --------------------------------------------------------------

      subroutine get_surface_information(coord,ires,values)
         ! Returns the desired values of the potential for an given set of internal coordinates
         ! Input:
         !     coord: cartesian coordinates of the desired point on PES
         !     ires: Information which properties are needed
         !           1: energies
         !           4: energies+dipoles
         !          10: energies+dipoles+polarisabilites
         ! Output:
         !     values: desired values
         !                values(1): energy
         !              values(2:4): dipoles (x,y,z)
         !             values(2:10): polarisabilites (xx,yy,zz,xy,xz,yz)
         implicit none
         integer :: ires
         double precision, dimension(ires) :: values
         double precision, dimension(ncart_p) :: coord
         double precision, dimension(ncrd_p) :: qvalp
         call get_current_elongation_intc(coord,qvalp)
         call get_values_from_qvalp(qvalp,values,ires)
      end subroutine

! ----------------------------------------------------------------------

      subroutine read_fit_coefs(fname)
         ! Reading of the fitting coefficients from the given potential file
         ! Input:
         !     fname: Name of the coefficient file
         use intc_info
         implicit none
         integer :: ios,idim,isurf,idum1,idum2,i,ipol
         integer, dimension(:,:), pointer :: coords
         character(len=30) :: text30
         character(len=16) :: text16
         character(len=1) :: text1
         double precision, dimension(:,:), allocatable :: coefs_all
         character(len=*) :: fname
         character(len=lch) :: line,cdum1,cdum2,ccheck,cdum3,cdum4,cdum5
         integer, dimension(:), pointer :: ifite,ifitd,ifitp
         double precision :: d_intermed
         double precision, dimension(:),allocatable :: qval_primitives
         call set_ci_from
         open(iunit_apot,file=trim(fname),action='read')
      
         floop: do
            read(iunit_apot,'(A)',iostat=ios) line
            if(ios.ne.0) exit floop

            ! Read the relevant parameters
            if(trim(line).eq.'### PARAMETERS') then
               do 
                  read(iunit_apot,'(A)',iostat=ios) line ! highest coupling of surface
                  if(trim(line).eq.'') exit
                  if(ios.ne.0) then
                     ! Todo: Write error warning in outfile
                     exit
                  endif
                  read(line,ci_form)text30,cdum2,idum1
                  if(trim(text30).eq.'Number of cartesian Coordinates') then
                     ncart_p = idum1
                  elseif(trim(text30).eq.'Order of energy coupling') then
                     ncoup_p = idum1
                  elseif(trim(text30).eq.'Order of dipole coupling') then
                     ncoupdip_p = idum1
                  elseif(trim(cdum1).eq.'Order of polar coupling') then
                     ncouppol_p = idum1
                  elseif(trim(cdum1).eq.'Number of coordinates') then
                     ncrd_p = idum1
                  elseif(trim(text30).eq.'Number of basis functions') then
                     numpol_p = idum1
                  elseif(trim(text30).eq.'Type of int. coord.') then
                     itypeint = idum1
                  endif
               enddo
               call init_poly_basic_for_reading
            !Internal Coordinates - Parameters
            elseif(trim(line)=='### PARAMETER INTERNAL COORDINATES') then
               !TODO: Add internal coordinate definitions!
               read(iunit_apot,'(A)',iostat=ios) line
               read(line,'(4x,A16,A1,I6)')text16,text1,nstretch
               read(iunit_apot,'(A)',iostat=ios) line
               read(line,'(4x,A16,A1,I6)')text16,text1,nbending
               read(iunit_apot,'(A)',iostat=ios) line
               read(line,'(4x,A16,A1,I6)')text16,text1,ntorsion
               read(iunit_apot,'(A)',iostat=ios) line
               call init_itypeint_1_2(ncrd_p)
               !call allocate_refcoord_ic(refcoord)
            !!Internal Coordinates - Definitions
            elseif(trim(line)=='### DEFINITION INTERNAL COORDINATES') then
               allocate(qval_primitives(ncrd_p))
               qval_primitives = 0.0d0
 
               read3loop: do
                  read(iunit_apot,'(A)',iostat=ios) line
                  if(trim(adjustl(line))=='Stretching') then
                     do i=1,nstretch
                        read(iunit_apot,'(A)',iostat=ios) line
                        read(line,'(4x,I5,I8,I8,I8,F18.12)')atstretch(3,i),atstretch(1,i),atstretch(2,i),&
                                                                    atstretch(4,i),d_intermed
                        qval_primitives(atstretch(3,i)) = d_intermed
                     enddo
                     !call get_connectivity
                  elseif(trim(adjustl(line))=='Bending') then
                     do i=1,nbending
                        read(iunit_apot,'(A)',iostat=ios) line
                        read(line,'(4x,I5,I8,I8,I8,F18.12)')atbending(4,i),atbending(1,i),atbending(2,i),&
                                                               atbending(3,i),d_intermed
                        qval_primitives(atbending(4,i)) = d_intermed
                    enddo
                  elseif(trim(adjustl(line))=='Torsion') then
                     do i=1,ntorsion
                        read(iunit_apot,'(A)',iostat=ios) line
                        read(line,'(4x,I5,I8,I8,I8,I8,F18.12)')attorsion(5,i),attorsion(1,i),attorsion(2,i),&
                                                          attorsion(3,i),attorsion(4,i),d_intermed
                        qval_primitives(attorsion(5,i)) = d_intermed
                     enddo
                     exit read3loop
                  endif
               enddo read3loop
               call determine_pf_intcs(pi_p)
               call allocate_refqvals(qval_primitives,ncrd_p)
               !call generate_refgeom(qval_primitives)
               deallocate(qval_primitives)
            ! Read in the fit functions
            elseif(trim(line)=='### FIT FUNCTIONS') then
               read(iunit_apot,'(A)',iostat=ios) line
               do i=1,ncrd_p
                  read(iunit_apot,'(A)',iostat=ios) line
                  read(line,'(4x,I5,I5,I5)')idum1,ifit_fctp(i,1),ifit_fctp(i,2)
               enddo
            ! Read the coefficients of the fitfunctions
            elseif(trim(line).eq.'### COEFFICIENTS') then
               do idim=1,ncoup_p
                  write(ccheck,'(A,i1,A)')'##  ',idim,'D SURFACES'
                  sloop: do
                     read(iunit_apot,'(A)',iostat=ios) line
                     if(ios.ne.0) exit
                     
                     ! Position, where the idim-Surfaces start
                     if(trim(line).eq.(trim(ccheck))) then
                        read(iunit_apot,'(A)',iostat=ios) line ! number of surfaces for coupling idim
                        read(line,ci_form)text30,text1,ivnsurf_p(idim)
                        read(iunit_apot,'(A)',iostat=ios) line ! number of properties
                        read(line,ci_form)text30,text1,ires_all_p(idim)
                        call set_dim_spec_forms(idim)
                        call init_dim_spec_stuff_poly(idim)
                        allocate(coefs_all(ires_all_p(idim),numpol_p**idim))
                        coords => p_modes(idim)%p
                        
 
                        ! Read every single surface of coupling idim
                        do isurf=1,ivnsurf_p(idim)
                           !read(iunit_apot,'(A)',iostat=ios) line ! empty line
                           read(iunit_apot,'(A)',iostat=ios) line ! Coords included in surface
                           read(line,form2)cdum1,idum1,cdum2,cdum3,(coords(i,isurf),i=1,idim),cdum4,idum2,cdum5
                           read(iunit_apot,'(A)',iostat=ios) line ! basis function energy
                           ! Coefficients
                           do ipol=1,numpol_p**idim
                              read(iunit_apot,formc)idum1,coefs_all(:,ipol)
                           enddo
                           call sort_and_store_coefs(idim,isurf,coefs_all)
                        enddo
                        exit sloop
                     endif
                  enddo sloop
                  deallocate(coefs_all)
               enddo
            endif
         enddo floop
         close(iunit_apot)
      end subroutine
!========================================================================
!=== Initialisation and Allocation ======================================
!========================================================================
      subroutine init_poly
         implicit none
         write(ci_form,'(A)')'(A50,A1,i5)'
      end subroutine

! ----------------------------------------------------------------------

      subroutine release_poly_after_fit
         implicit none
         call dealloc_poly_fitfct_and_modes
         call dealloc_poly_chi
         call dealloc_poly_general
         call dealloc_poly_basic
      end subroutine

! ----------------------------------------------------------------------

      subroutine dealloc_poly_for_reading
         implicit none
         call dealloc_dim_spec_stuff_poly
         deallocate(p_modes)
         deallocate(p_coef_ene)
         deallocate(p_coef_dip)
         deallocate(p_coef_pol)
         deallocate(ires_all_p)
         deallocate(ivnsurf_p)
      end subroutine

! ----------------------------------------------------------------------

      subroutine init_poly_basic(ncoup_s_pot,ncoupdip_s_pot,ncouppol_s_pot)
         implicit none
         integer :: ncoup_s_pot,ncoupdip_s_pot,ncouppol_s_pot
         ncoup_p = ncoup_s_pot
         ncoupdip_p = ncoupdip_s_pot
         ncouppol_p = ncouppol_s_pot
         allocate(p_modes(ncoup_p))
         allocate(p_coef_ene(ncoup_p))
         allocate(p_coef_dip(ncoupdip_p))
         allocate(p_coef_pol(ncouppol_p))
         allocate(ires_all_p(ncoup_p))
         allocate(ivnsurf_p(ncoup_p))
         allocate(p_chi(ncoup_p))
         allocate(p_chi_dip(ncoupdip_p))
         allocate(p_chi_pol(ncouppol_p))
      end subroutine

! ----------------------------------------------------------------------

      subroutine dealloc_poly_basic
         implicit none
         deallocate(p_modes)
         deallocate(p_coef_ene)
         deallocate(p_coef_dip)
         deallocate(p_coef_pol)
         deallocate(ires_all_p)
         deallocate(ivnsurf_p)
         deallocate(p_chi)
         deallocate(p_chi_dip)
         deallocate(p_chi_pol)
      end subroutine
! ---- Chi values ------------------------------------------------------

      subroutine alloc_poly_chi
         implicit none
         allocate(chi_max(ncoup_p))
         allocate(chi_av(ncoup_p))
      end subroutine

! ----------------------------------------------------------------------

      subroutine dealloc_poly_chi
         implicit none
         deallocate(chi_max)
         deallocate(chi_av)
      end subroutine

! ----------------------------------------------------------------------

      subroutine init_poly_general
         implicit none
         integer :: idim,isurf
         integer, dimension(:), pointer :: grid_total_p
         integer, dimension(:,:), pointer :: ifit,modes
         integer, dimension(:), pointer :: ifite,ifitd,ifitp
         double precision, dimension(:,:,:), pointer :: fine_ene
         ibas_min = 7
         ibas_max = 10
         allocate(fine_grid(ngrid_max_p,ncrd_p))
         call get_fine_grid
         call alloc_poly_chi
         chi_max(1) = 1.0d-10
         chi_av(1) = 1.0d-10
         if(ncoup_p.ge.2) then
            chi_av(2) = 1.0d-10
            chi_max(2) = 1.0d-8
         endif
         if(ncoup_p.ge.3) then
            chi_av(3) = 1.0d-10
            chi_max(3) = 1.0d-6
         endif
         do idim=4,ncoup_p
            chi_max(idim) = 1.0d-5
            chi_av(idim) = 1.0d-5
         enddo
      end subroutine

! ----------------------------------------------------------------------

      subroutine dealloc_poly_general
         implicit none
         deallocate(fine_grid)
      end subroutine

! ---- Fit function ----------------------------------------------------

      subroutine init_poly_fitfct_and_modes
         implicit none
         integer :: idim
         allocate(xmin_p(ncrd_p))
         allocate(xmax_p(ncrd_p))
         allocate(ifit_fctp(ncrd_p,2))
         do idim=1,ncoup_p
            allocate(p_modes(idim)%p(idim,ivnsurf_p(idim)))
            allocate(p_chi(idim)%p(ivnsurf_p(idim)))
         enddo
         do idim=1,ncoupdip_p
            allocate(p_chi_dip(idim)%p(3,ivnsurf_p(idim)))
         enddo
         do idim=1,ncouppol_p
            allocate(p_chi_pol(idim)%p(6,ivnsurf_p(idim)))
         enddo
      end subroutine
      
! ----------------------------------------------------------------------

      subroutine dealloc_poly_fitfct_and_modes
         implicit none
         integer :: idim
         deallocate(xmin_p)
         deallocate(xmax_p)
         deallocate(ifit_fctp)
         do idim=1,ncoup_p
            deallocate(p_coef_ene(idim)%p)
            deallocate(p_modes(idim)%p)
            deallocate(p_chi(idim)%p)
         enddo
         do idim=1,ncoupdip_p
            deallocate(p_coef_dip(idim)%p)
            deallocate(p_chi_dip(idim)%p)
         enddo
         do idim=1,ncouppol_p
            deallocate(p_coef_pol(idim)%p)
            deallocate(p_chi_pol(idim)%p)
         enddo
      end subroutine

! ----------------------------------------------------------------------

      subroutine get_fine_grid
         implicit none
         integer :: imode,igrid
         double precision :: hgrid
         do imode=1,ncrd_p
            hgrid=abs(xmax_p(imode)-xmin_p(imode))/dble(ngrid_max_p-1)
            do igrid=1,ngrid_max_p
               fine_grid(igrid,imode) = xmin_p(imode) + dble(igrid-1) * hgrid
            enddo
         enddo
      end subroutine

!---- basic for reading ------------------------------------------------

      subroutine init_poly_basic_for_reading
      ! Todo: Need to be deallocated as well
         implicit none
         allocate(p_modes(ncoup_p))
         allocate(p_coef_ene(ncoup_p))
         allocate(p_coef_dip(ncoupdip_p))
         allocate(p_coef_pol(ncouppol_p))
         allocate(ifit_fctp(ncrd_p,2))
         allocate(ires_all_p(ncoup_p))
         allocate(ivnsurf_p(ncoup_p))
      end subroutine

! ----------------------------------------------------------------------

      subroutine dealloc_poly_basic_for_reading
      ! Todo: Need to be deallocated as well
         implicit none
         deallocate(p_modes)
         deallocate(p_coef_ene)
         deallocate(p_coef_dip)
         deallocate(p_coef_pol)
         deallocate(ifit_fctp)
         deallocate(ires_all_p)
         deallocate(ivnsurf_p)
      end subroutine

!---------------------------------------------------------------------

      subroutine init_dim_spec_stuff_poly(idim)
         implicit none
         integer :: idim
         allocate(p_coef_ene(idim)%p(numpol_p**idim,ivnsurf_p(idim)))
         allocate(p_modes(idim)%p(idim,ivnsurf_p(idim)))
         if(idim.le.ncoupdip_p) then
            ! Basis functions
            allocate(p_coef_dip(idim)%p(numpol_p**idim,3,ivnsurf_p(idim)))
         endif
         if(idim.le.ncouppol_p) then
            ! Basis functions
            allocate(p_coef_pol(idim)%p(numpol_p**idim,6,ivnsurf_p(idim)))
         endif
      end subroutine

!---------------------------------------------------------------------

      subroutine dealloc_dim_spec_stuff_poly
         implicit none
         integer :: idim
         do idim=1,ncoup_p
            ! Basis functions
            deallocate(p_coef_ene(idim)%p)
            deallocate(p_modes(idim)%p)
         enddo
         do idim=1,ncoupdip_p
            ! Basis functions
            deallocate(p_coef_dip(idim)%p)
         enddo
         do idim=1,ncouppol_p
            ! Basis functions
            deallocate(p_coef_pol(idim)%p)
         enddo
      end subroutine

!---------------------------------------------------------------------

      subroutine set_ci_from
         implicit none
         write(ci_form,'(A)')'(4x,A30,1x,A1,1x,i10)'
      end subroutine

!---------------------------------------------------------------------

      subroutine set_dim_spec_forms(idim)
         implicit none
         integer :: idim
         write(form2,'(A,I1,A)') '(A4,I1,A10,A7,',idim,'I4,3x,A1,I5,A1)'
         write(formc,'(a,i2,a)') '(I5,',ires_all_p(idim),'es20.10)'
      end subroutine
!========================================================================
!=== Potfile for fit ====================================================
!========================================================================
      subroutine sort_and_store_coefs(idim,isurf,coefs_all)
         implicit none
         integer :: idim,ipol,itype,icnt,isurf
         double precision, dimension(ires_all_p(idim),numpol_p**idim) :: coefs_all
         double precision, dimension(:,:), pointer :: coef_ene
         double precision, dimension(:,:,:), pointer :: coef_dip
         double precision, dimension(:,:,:), pointer :: coef_pol

         ! Energy
         coef_ene => p_coef_ene(idim)%p
         do ipol=1,numpol_p**idim
            coef_ene(ipol,isurf) = coefs_all(1,ipol)
         enddo

         !Dipole
         icnt = 1
         do itype=2,min(ires_all_p(idim),4)
            if(icnt.eq.1) coef_dip => p_coef_dip(idim)%p
            do ipol=1,numpol_p**idim
               coef_dip(ipol,icnt,isurf) = coefs_all(itype,ipol)
            enddo
            icnt = icnt + 1
         enddo
         
         ! Polarisability
         icnt = 1
         do itype=5,min(ires_all_p(idim),10)
            if(icnt.eq.1) coef_pol => p_coef_pol(idim)%p
            do ipol=1,numpol_p**idim
               coef_pol(ipol,icnt,isurf) = coefs_all(itype,ipol)
            enddo
            icnt = icnt + 1
         enddo
      end subroutine
      
! --------------------------------------------------------------

      subroutine get_values_from_qvalp(qvalp,values,ires)
         use intc_info
         implicit none
         integer :: ires,ncrd_all,idim_ic,i,isrf_idx,kres,icnt,ispol,isdip
         integer, dimension(:),allocatable :: icrds_all
         integer, dimension(:),allocatable :: icrds,mpos,ipos
         double precision :: val
         double precision, dimension(ires) :: values
         double precision, dimension(ncrd_p) :: qvalp
         double precision, dimension(:), allocatable:: disfac_ic
         logical :: ilog
         values = 0.0d0
      
         isdip = 0
         ispol = 0
         if(ires.ge.2) isdip = 1
         if(ires.ge.5) ispol = 1
         ! Get the contributing coordinates
         allocate(icrds_all(ncrd_p))
         ncrd_all = 0
         call adapt_qvalp_elong(qvalp,ncrd_p,iout_p)
         do i=1,ncrd_p
            if(abs(qvalp(i)).gt.1.0d-10) then
               ncrd_all = ncrd_all + 1
               icrds_all(ncrd_all) = i
            endif
         enddo
         do idim_ic=1,min(ncrd_all,ncoup_p)
            !TODO:
            ! 1. Think about what happens, if a needed surface is not provided due
            !    to dimensionality or skipping
            ! 2. Also think about what happens if the desired value is outside the
            !    desired range
            allocate(icrds(idim_ic))
            allocate(mpos(idim_ic))
            allocate(ipos(idim_ic))
            allocate(disfac_ic(idim_ic))
            do i=1,idim_ic
               mpos(i) = ncrd_all-idim_ic+i
               ipos(i) = i
            enddo
            if(ires_all_p(idim_ic).ge.2) isdip = isdip * 1
            if(ires_all_p(idim_ic).ge.5) ispol = ispol * 1
            srfloop: do
               ! TODO: Make a check if this surface exists
               do i=1,idim_ic
                  icrds(i) = icrds_all(ipos(i))
                  disfac_ic(i) = qvalp(icrds(i))
               enddo
               call get_isurf(icrds,isrf_idx,idim_ic)
               ! Energy
               icnt = 1
               call fit_surface_ic_ene(idim_ic,disfac_ic,isrf_idx,val)
               values(icnt) = values(icnt) + val
               ! Dipole
               do i=1,3*isdip
                  icnt = icnt + 1
                  call fit_surface_ic_dip(idim_ic,i,disfac_ic,isrf_idx,val)
                  values(icnt) = values(icnt) + val
               enddo
               ! Polar
               do i=1,6*ispol
                  icnt = icnt + 1
                  call fit_surface_ic_pol(idim_ic,i,disfac_ic,isrf_idx,val)
                  values(icnt) = values(icnt) + val
               enddo
               call odo3(ilog,ipos,mpos,idim_ic)
               if(.not.ilog) exit srfloop
            enddo srfloop
            deallocate(icrds)
            deallocate(mpos)
            deallocate(ipos)
            deallocate(disfac_ic)
         enddo
         deallocate(icrds_all)
      end subroutine
      
! --------------------------------------------------------------

      subroutine get_isurf(icrds,idxsurf,idim)
         implicit none
         integer :: idim,isurf,idxsurf,i,icnt
         integer, dimension(idim) :: icrds
         integer, dimension(:,:), pointer, contiguous :: modes
         modes => p_modes(idim)%p

         idxsurf=0
         do isurf=1,ivnsurf_p(idim)
            icnt = 0
            do i=1,idim
               if(modes(i,isurf).ne.icrds(i)) exit
               icnt = icnt + 1
            enddo
            if(icnt.eq.idim) then
               idxsurf = isurf
               return
            endif
         enddo
         
         if(idxsurf.eq.0) then
            write(iout_p,*)'No surface selected'
            stop
         endif
         
      end subroutine

! --------------------------------------------------------------

      subroutine fit_surface_ic_ene(idim,disfac,isurf,eval)
         implicit none
         integer :: idim,isurf,i
         integer, dimension(idim) :: ivbas,ivgrid
         integer, dimension(idim) :: ifit_ene
         integer, dimension(:,:), pointer :: modes
         double precision :: eval
         double precision, dimension(idim) :: disfac
         double precision, dimension(:,:), pointer, contiguous :: coef_ene
         ivbas = numpol_p
         ivgrid = 1
         coef_ene => p_coef_ene(idim)%p
         modes => p_modes(idim)%p
         do i=1,idim
            ifit_ene(i) = ifit_fctp(modes(i,isurf),1)
         enddo
         call calc_V_kronecker(idim,ivgrid,ivbas,disfac,eval,ifit_ene&
                               ,coef_ene(:,isurf))
      
      end subroutine

! --------------------------------------------------------------

      subroutine fit_surface_ic_dip(idim,kres,disfac,isurf,eval)
         implicit none
         integer :: idim,isurf,kres,i
         integer, dimension(idim) :: ivbas,ivgrid
         integer, dimension(idim) :: ifit_dip
         integer, dimension(:,:), pointer :: modes
         double precision :: eval
         double precision, dimension(idim) :: disfac
         double precision, dimension(:,:,:), pointer, contiguous :: coef_dip
         ivbas = numpol_p
         ivgrid = 1
         coef_dip => p_coef_dip(idim)%p
         modes => p_modes(idim)%p
         do i=1,idim
            ifit_dip(i) = ifit_fctp(modes(i,isurf),2)
         enddo
         call calc_V_kronecker(idim,ivgrid,ivbas,disfac,eval,ifit_dip&
                               ,coef_dip(:,kres,isurf))
      
      end subroutine

! --------------------------------------------------------------

      subroutine fit_surface_ic_pol(idim,kres,disfac,isurf,eval)
         implicit none
         integer :: idim,isurf,kres,i
         integer, dimension(idim) :: ivbas,ivgrid
         integer, dimension(idim) :: ifit_pol
         integer, dimension(:,:), pointer :: modes
         double precision :: eval
         double precision, dimension(idim) :: disfac
         double precision, dimension(:,:,:), pointer, contiguous :: coef_pol
         ivbas = numpol_p
         ivgrid = 1
         coef_pol => p_coef_pol(idim)%p
         modes => p_modes(idim)%p
         do i=1,idim
            ifit_pol(i) = ifit_fctp(modes(i,isurf),2)
         enddo
         call calc_V_kronecker(idim,ivgrid,ivbas,disfac,eval,ifit_pol &
                               ,coef_pol(:,kres,isurf))
      
      end subroutine
end module
