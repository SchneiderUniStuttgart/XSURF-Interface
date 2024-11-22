module intc_info
   !use surf_dat
   implicit none
   save
   !Internal Coordinates generall
!   double precision, dimension(:),pointer,contiguous :: refcoord_ic     ! intern used reference geometry
   integer :: itypeint
   double precision, dimension(:),allocatable :: refqvals        ! reference qvals
!   double precision, dimension(:,:),pointer,contiguous :: chargemass    ! Charge and Mass of Atoms

   !Primitive Coordinates
   integer :: nstretch                                                  ! Number of pure Stretch
   integer, dimension(:,:),allocatable :: atstretch              ! Atoms connected to each other
   integer, dimension(:),allocatable :: iconnect                 ! Connectivity of Atoms
   integer :: nbending                                                  ! Number of pure Bond angle
   integer, dimension(:,:),allocatable :: atbending              ! Atoms for Bond angles
   integer :: ntorsion                                                  ! Number of pure Torsion angles
   integer, dimension(:,:),allocatable:: attorsion              ! Atoms for Torsion angles
   double precision, dimension(:),allocatable :: pf_intc       ! prefactor for the internal coordinates
   contains
!====================================================================================
!===============================General==============================================
!====================================================================================
      subroutine allocate_refqvals(qval_tmp,ncrd)
         implicit none
         integer :: icrd,ncrd
         double precision, dimension(ncrd) :: qval_tmp
         allocate(refqvals(ncrd))
         if(itypeint.eq.1) then
            do icrd=1,ncrd
               refqvals(icrd) = qval_tmp(icrd) * pf_intc(icrd)
            enddo
         endif
      end subroutine

!---------------------------------------------------------------------

      subroutine chose_fit_function_intc(icrd,icase)
         ! This Routine selects a fit function for fitting 1D-Surfaces
         ! Input:
         !     icrd: Number of Coordinate
         ! Output:
         !     icase: Type of fit function (compare to polytools.F90 -> calc_Ai)
         implicit none
         integer :: icrd,icase

         icase = -1
         if(itypeint.eq.1) then
            if(any(atstretch(3,:).eq.icrd)) then
               !Stretching (Morse)
               icase = 4
            elseif(any(atbending(4,:).eq.icrd)) then
               !Bending (Polynominals)
               icase = 3
            elseif(any(attorsion(5,:).eq.icrd)) then
               !Torsion (Trigonometric)
               icase = 5
            endif
         else
            !TODO: Extend to more coordinates
            icase = 3 !(Polynominals)
         endif

      end subroutine

!---------------------------------------------------------------------

      subroutine determine_pf_intcs(pi_i)
         ! Defining scaling factors to obtain right dimension of coordinates
         ! Important: The Internal coordinates need to be defined first
         ! TODO: Extend to more coordinate types
         implicit none
         integer :: icrd
         double precision :: pi_i

         if(itypeint.eq.1) then
            do icrd=1,nstretch
               pf_intc(atstretch(3,icrd)) = 1.0d0
            enddo
            do icrd=1,nbending
               pf_intc(atbending(4,icrd)) = 1.0d0/(2.0d0*pi_i)
            enddo
            do icrd=1,ntorsion
               pf_intc(attorsion(5,icrd)) = 1.0d0/(2.0d0*pi_i)
            enddo
         else
            pf_intc = 1.0d0
         endif

      end subroutine
!====================================================================================
!=============================Primitives=============================================
!====================================================================================
      subroutine init_itypeint_1_2(ncrd)
         ! Some additional initialisation for primitive and delocalised internal coordinates
         implicit none
         integer :: ncrd

         !allocate(iconnect(ncart))
         !iconnect = 0
         ncrd = nstretch + nbending + ntorsion
         allocate(atstretch(4,nstretch))
         allocate(atbending(4,nbending))
         allocate(attorsion(5,ntorsion))
         !if(itypeint.eq.2) kmat => memory_allocate(nicrd,nmodes_norm,description='kmat')
         allocate(pf_intc(ncrd))

      end subroutine

!---------------------------------------------------------------------

      subroutine release_itypeint_1_2
         implicit none
         !deallocate(iconnect)
         deallocate(atstretch)
         deallocate(atbending)
         deallocate(attorsion)
         deallocate(pf_intc)
         deallocate(refqvals)
      end subroutine

!---------------------------------------------------------------------

      subroutine adapt_qvalp_elong(qvalp,ncrd,iout)
         implicit none
         integer :: itor,icnt2,ncrd,iout
         double precision, dimension(ncrd) :: qvalp
         do itor=1,ntorsion
            icnt2 = 0
            correctloop: do
               if(icnt2.eq.20) then
                  write(iout,*)' The torsional angle is to far elongated'
                  write(iout,*)' mod_intc_head>adapt_torsion'
                  stop
               endif
               if(qvalp(attorsion(5,itor)).lt.-0.5d0) then
                  qvalp(attorsion(5,itor)) = qvalp(attorsion(5,itor)) + 1.0d0
                  icnt2 = icnt2 + 1
                  cycle correctloop
               elseif(qvalp(attorsion(5,itor)).gt.0.5d0) then
                  qvalp(attorsion(5,itor)) = qvalp(attorsion(5,itor)) - 1.0d0
                  icnt2 = icnt2 + 1
                  cycle correctloop
               endif
               exit correctloop
            enddo correctloop
         enddo
      end subroutine
end module
