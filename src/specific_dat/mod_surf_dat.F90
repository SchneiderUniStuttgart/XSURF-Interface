module surf_dat
   use gen_dat
   save
   ! objects
   type iptr_array_1
      integer, dimension(:), pointer, contiguous :: p
   end type iptr_array_1
   !type iptr_array_2
   !   integer, dimension(:,:), pointer, contiguous :: p
   !end type iptr_array_2
   type ptr_array_2
      double precision, dimension(:,:), pointer, contiguous :: p
   end type ptr_array_2
   type ptr_array_3
      double precision, dimension(:,:,:), pointer, contiguous :: p
   end type ptr_array_3

   ! General Surface informations
   integer :: ncrd_s                                                             ! Number of coordinates
   integer :: ncoup_s_pot                                                        ! Coupling order of energy in Pot-file
   integer :: ncoupdip_s_pot                                                     ! Coupling order of dipole in Pot-file
   integer :: ncouppol_s_pot                                                     ! Coupling order of polarisability in Pot-file
   integer,dimension(:), allocatable :: ires_all_s                               ! number of computed properties
   integer, dimension(:), allocatable :: ivnsurf_s                               ! number of surfaces per dimension
   double precision, dimension(:), allocatable, public :: xmin_s                 ! Min value of elongation
   double precision, dimension(:), allocatable, public :: xmax_s                 ! Max value of elongation
   type(iptr_array_1), dimension(:),allocatable :: i_labxd_s
                  !integer :: idirec_vmult                                       ! Is a Multi-Calc done?

   ! Grid points
   integer, dimension(:), allocatable :: igrid_all                               ! All Grid-Points depending on dimension
   type(iptr_array_1), dimension(:), allocatable:: p_igrid_total                 ! All Grid-Points in one dimension -> total number
   integer, dimension(:),pointer, contiguous :: igrid_total                      ! All Grid-Points in one dimension -> total number
   type(iptr_array_1), dimension(:), allocatable :: p_isurface_adr               ! Specifier for points
   integer, dimension(:), pointer, contiguous :: isurface_adr                    ! Specifier for points
   integer, dimension(:,:),allocatable :: igrid_nabi                             ! Number points
   type(ptr_array_2), dimension(:), allocatable :: p_grid_abi_points
   double precision, dimension(:,:), pointer, contiguous :: grid_abi_ptr         ! Complete grid points for one surface
   double precision, dimension(:,:), pointer, contiguous :: grid_s               ! Grid points
   double precision, dimension(:,:), pointer, contiguous :: abi_points           ! Values corresponding to grid points
   type(ptr_array_3), dimension(:), allocatable :: p_fine_grid_points
   double precision, dimension(:,:,:),pointer, contiguous :: fine_grid_points    !fine grid points

   ! Fitting of PES
   integer :: ngrid_max_s                                                        ! Maximum number of fine grid points
   integer, dimension(:),allocatable :: ivbas_f                                  ! Number of fit functions
   integer, dimension(:),allocatable :: ivbas_max                                ! Maximum number of basis functions
   integer, dimension(:,:),allocatable :: ivgrid_f                               ! Number of Grid Points
   integer, dimension(:,:),allocatable :: ifit_fcts                              ! Fit functions
   double precision, dimension(:,:), allocatable :: chivals_f                    ! Chi values
   contains
      subroutine release_surf_dat
         implicit none
         integer :: idim
         do idim=1,ncoup_s_pot
            call dealloc_i_labxd_s(idim)
            call dealloc_abi_grid(idim)
            call dealloc_fine_grid_points(idim)
            call dealloc_surfadr(idim)
         enddo
         call dealloc_fitting
         call dealloc_surf_dat_basic
      end subroutine
!========================================================================
!=== Initialisation and Allocation ======================================
!========================================================================
      subroutine print_basic_surface_information
         use intc_info, only: itypeint
         implicit none
         character(len=30) :: text30
         character(len=5) :: ct
         write(iout,'(/,1x,5x,a24)')'Information from potfile'
         write(text30,'(A)')'Order of energy coupling'
         write(iout,'(1x,7x,A30,a1,i2)')text30,'=',ncoup_s_pot
         write(text30,'(A)')'Order of dipole coupling'
         write(iout,'(1x,7x,A30,a1,i2)')text30,'=',ncoupdip_s_pot
         write(text30,'(A)')'Order of polar coupling'
         write(iout,'(1x,7x,A30,a1,i2)')text30,'=',ncouppol_s_pot
         write(text30,'(A)')'Type of int. Coord:'
         call get_typeint(itypeint,ct)
         write(iout,'(1x,7x,A30,a1,i2,1x,a1,a5,a1)')text30,'=',itypeint,'(',ct,')'
      end subroutine

!---------------------------------------------------------------------

      subroutine get_typeint(it,ct)
         implicit none
         integer :: it
         character(len=5) :: ct
         if(it.eq.1) then
            write(ct,'(a5)')'ZMAT'
         endif
      end subroutine

!---------------------------------------------------------------------

      subroutine print_dim_specific_headr(idim)
         implicit none
         integer :: idim
         character(len=30) :: text30
         write(iout,'(/,1x,7x,i1,a10)')idim,'D Surfaces'
         write(text30,'(A)')'Number of surfaces'
         write(iout,'(1x,7x,A30,a1,i10)')text30,'=',ncoup_s_pot
      end subroutine

!========================================================================
!=== Initialisation and Allocation ======================================
!========================================================================
      subroutine init_surf_dat_basic
         ! Initialisation and allocation of the basic surface informations
         implicit none
         allocate(ivnsurf_s(ncoup_s_pot))
         ivnsurf_s = 0
         allocate(igrid_all(ncoup_s_pot))
         igrid_all = 0
         allocate(xmin_s(ncrd_s))
         xmin_s = 0.0d0
         allocate(xmax_s(ncrd_s))
         xmax_s = 0.0d0
         allocate(i_labxd_s(ncoup_s_pot))
         allocate(ires_all_s(ncoup_s_pot))
         allocate(p_isurface_adr(ncoup_s_pot))
         allocate(p_igrid_total(ncoup_s_pot))
         allocate(p_grid_abi_points(ncoup_s_pot))
         allocate(p_fine_grid_points(ncoup_s_pot))
      end subroutine

!---------------------------------------------------------------------

      subroutine dealloc_surf_dat_basic
         ! Deallocation of the basic surface informations
         implicit none
         deallocate(ivnsurf_s)
         deallocate(ires_all_s)
         deallocate(igrid_all)
         deallocate(xmin_s)
         deallocate(xmax_s)
         deallocate(i_labxd_s)
         deallocate(p_isurface_adr)
         deallocate(p_igrid_total)
         deallocate(p_grid_abi_points)
         deallocate(p_fine_grid_points)
      end subroutine

!---- Fitting --------------------------------------------------------

      subroutine init_fitting
         implicit none
         ngrid_max_s = 16
         allocate(ivbas_f(ncoup_s_pot))
         ivbas_f = 9
         allocate(ivbas_max(ncoup_s_pot))
         ivbas_max = 10
         allocate(ivgrid_f(ncrd_s,ncoup_s_pot))
         ivgrid_f = 16
         allocate(ifit_fcts(ncrd_s,2))
      end subroutine

!---------------------------------------------------------------------

      subroutine dealloc_fitting
         implicit none
         integer :: idim
         deallocate(ivbas_f)
         deallocate(ivbas_max)
         deallocate(ivgrid_f)
         deallocate(ifit_fcts)
      end subroutine

!---- LabXD ----------------------------------------------------------

      subroutine alloc_i_labxd_s(idim)
         ! Allocation of the labeling for the surfaces
         implicit none
         integer :: idim
         allocate(i_labxd_s(idim)%p(ivnsurf_s(idim)*idim))
         i_labxd_s(idim)%p = 0
      end subroutine

!---------------------------------------------------------------------

      subroutine dealloc_i_labxd_s(idim)
         ! Deallocation of the labeling for the surfaces
         implicit none
         integer :: idim
         deallocate(i_labxd_s(idim)%p)
      end subroutine

!---- Dim specific ---------------------------------------------------

      subroutine init_dim_spec_stuff(idim,ires_read)
         implicit none
         integer :: idim,ires_read
         allocate(p_isurface_adr(idim)%p(ivnsurf_s(idim)+1))

         allocate(p_igrid_total(idim)%p(ivnsurf_s(idim)))

         allocate(p_grid_abi_points(idim)%p(idim+ires_read,igrid_all(idim)))
         call coarse_grid(idim,ires_read)
         isurface_adr = 0
         igrid_total = 0
         grid_abi_ptr = 0.0d0
         
         call alloc_igrid_nabi(idim)

         allocate(p_fine_grid_points(idim)%p(ngrid_max_s**idim,ires_read,ivnsurf_s(idim)))
         fine_grid_points => p_fine_grid_points(idim)%p
         fine_grid_points = 0.0d0
      end subroutine

!---------------------------------------------------------------------

      subroutine coarse_grid(idim,ires_read)
         implicit none
         integer :: idim,ires_read
         isurface_adr => p_isurface_adr(idim)%p

         igrid_total => p_igrid_total(idim)%p

         grid_abi_ptr => p_grid_abi_points(idim)%p
         grid_s => grid_abi_ptr(1:idim,:)
         abi_points => grid_abi_ptr(idim+1:idim+ires_read,:)
      end subroutine

!---------------------------------------------------------------------

      subroutine alloc_igrid_nabi(idim)
         implicit none
         integer :: idim
         allocate(igrid_nabi(idim,ivnsurf_s(idim)))
         igrid_nabi = 0
      end subroutine

!---------------------------------------------------------------------

      subroutine dealloc_fine_grid_points(idim)
         implicit none
         integer :: idim
         deallocate(p_fine_grid_points(idim)%p)
      end subroutine

!---------------------------------------------------------------------

      subroutine dealloc_surfadr(idim)
         implicit none
         integer :: idim
         deallocate(p_isurface_adr(idim)%p)
      end subroutine

!---------------------------------------------------------------------

      subroutine dealloc_igrid_nabi(idim)
         implicit none
         integer :: idim
         deallocate(igrid_nabi)
      end subroutine

!---------------------------------------------------------------------

      subroutine dealloc_abi_grid(idim)
         implicit none
         integer :: idim
         deallocate(p_igrid_total(idim)%p)
         deallocate(p_grid_abi_points(idim)%p)
      end subroutine
!========================================================================
!=== Set needed variables ======================================
!========================================================================
      subroutine iset_lab_grid_s(idim,isurf,idx)
         implicit none
         integer :: idim,isurf,i
         integer, dimension(idim) :: idx
         integer, dimension(:), pointer :: i_lab

         i_lab => i_labxd_s(idim)%p

         do i=1,idim
            i_lab(idim*(isurf-1)+i) = idx(i)
         enddo

      end subroutine

!---------------------------------------------------------------------

      subroutine iget_lab_grid_s(idim,isurf,idx)
         implicit none
         integer :: idim,isurf,i
         integer, dimension(idim) :: idx
         integer, dimension(:), pointer :: i_lab

         i_lab => i_labxd_s(idim)%p

         do i=1,idim
            idx(i)=i_lab(idim*(isurf-1)+i)
         enddo

      end subroutine

!---------------------------------------------------------------------

      subroutine select_fit_function
!         !TODO: Rewrite this to read from input
!         use intc_info, only: chose_fit_function_intc
!         use intc_info, only: itypeint
!         implicit none
!         integer :: i
!         integer, dimension(:,:), pointer :: ifit
!         ifit = 3
!         if(idim.eq.1.AND.itypeint.ge.1) then
!            do i=1,ivnsurf_s(idim)
!               call chose_fit_function_intc(i,ifit(1,i))
!            enddo
!         endif
      end subroutine

!---------------------------------------------------------------------

      subroutine get_right_surface_stuff(idim,ires_read)
         implicit none
         integer :: idim,ires_read
         igrid_total => p_igrid_total(idim)%p
         fine_grid_points => p_fine_grid_points(idim)%p
         grid_abi_ptr => p_grid_abi_points(idim)%p
         grid_s => grid_abi_ptr(1:idim,:)
         abi_points => grid_abi_ptr(idim+1:idim+ires_read,:)
         isurface_adr => p_isurface_adr(idim)%p
      end subroutine
end module
