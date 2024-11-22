subroutine iones(array,ndim)
   implicit none
   integer :: ndim,i
   integer, dimension(ndim) :: array
   do i=1,ndim
      array(i) = 1
   enddo
end subroutine

! ----------------------------------------------------------------------

subroutine odo(ilog,ix,m,n)
   !Routine for generalized loops, increase first index first
   implicit none
   integer :: n,i
   logical :: ilog
   integer, dimension(n) :: ix
   integer, dimension(n) :: m

   i=1
   do while (i<=n)
      ix(i)=ix(i)+1
      if(ix(i)<=m(i)) then
         ilog=.TRUE.
         return
      endif
      ix(i)=1
      i=i+1
   end do
   ilog=.FALSE.

end subroutine

! ----------------------------------------------------------------------

integer function iget_iad2(idim,idxj,nu_grid,ipos_zero,idx) result(iad2)
   !This Routine determines the position of a sub-surface
   implicit none
   integer :: iprod,idim,i
   integer, dimension(idim) :: idxj,nu_grid,ipos_zero,idx
   iprod=1
   iad2=1
   do i=1,idim
      if(idx(i)==1)then
         iad2=iad2+(idxj(i)-1)*iprod
      else
         iad2=iad2+(ipos_zero(i)-1)*iprod
      endif
      iprod=iprod*nu_grid(i)
   end do
end function

! ----------------------------------------------------------------------

subroutine print1Darray(array,idim1)
   use gen_dat
   implicit none
   integer, intent(in) :: idim1
   real(8) array(idim1)
   integer i
   character(len=1024) :: form

   write(form,'(a)')'(i5,f12.4)'
   do i=1,idim1
      write(iout,form) i,array(i)
   end do
   write(iout,'(1x)')
end subroutine

! ----------------------------------------------------------------------

subroutine print3Darray(array,idim1,idim2,idim3)
   use gen_dat
   implicit none
   integer, intent(in) :: idim1,idim2,idim3
   real(8) array(idim1,idim2,idim3)
   integer i,j,k
   character(len=1024) :: form

   write(form,'(a,i3,a)')'(i5,i5,',idim3,'f15.9)'
   do i=1,idim1
      do j=1,idim2
         write(iout,form) i,j,(array(i,j,k),k=1,idim3)
      end do
      write(iout,'(1x)')
   end do
   write(iout,'(1x)')
end subroutine

! ----------------------------------------------------------------------

subroutine index_gen_one_elem(n,idim,i,idx,isort)
   implicit none
   integer :: n,idim,i,isort,k
   integer, dimension(idim) :: idx

   if(isort==1)then
      do k=1,idim
         idx(idim-k+1)=mod(floor(dble(i-1)/dble(n**(idim-k))),n)+1
      enddo
   else
      do k=idim,1,-1
         idx(k)=mod(floor(dble(i-1)/dble(n**(idim-k))),n)+1
      enddo
   endif
end subroutine

! ----------------------------------------------------------------------

subroutine surf_to_poly
   use surf_dat
   use poly
   implicit none
   integer :: idim,isurf
   integer, dimension(:,:), pointer :: modes
   ncrd_p = ncrd_s
   ngrid_max_p = ngrid_max_s
   ivnsurf_p = ivnsurf_s
   ires_all_p = ires_all_s
   call init_poly_fitfct_and_modes
   do idim=1,ncrd_p
      xmin_p(idim) = xmin_s(idim)
      xmax_p(idim) = xmax_s(idim)
      ifit_fctp(idim,1) = ifit_fcts(idim,1)
      ifit_fctp(idim,2) = ifit_fcts(idim,2)
   enddo
   do idim=1,ncoup_p
      modes => p_modes(idim)%p
      do isurf=1,ivnsurf_p(idim)
         call iget_lab_grid_s(idim,isurf,modes(:,isurf))
      enddo
   enddo
end subroutine

! ----------------------------------------------------------------------

subroutine calc_coef_xpoly(ires_read)
   !Fitting of the fine grid
   use surf_dat, only: fine_grid_points,p_igrid_total,ifit_fcts
   use surf_dat, only: get_right_surface_stuff
   use gen_dat
   use poly
   implicit none
   integer :: ibas,ibas_work,idim,isurf,isum,i,iprod,ires_read,itype
   integer :: irespoly,ncoup_type,itl
   integer, dimension(ncoup_p) :: ivbas
   integer, dimension(ncoup_p) :: ivgrid ,ivfinegrid
   integer, dimension(:,:), pointer :: modes
   integer, dimension(:), pointer :: grid_total
   integer, dimension(:), allocatable :: ifit

   double precision :: chi_sum
   double precision, dimension(:,:), pointer :: fine_ene
   double precision, dimension(:,:), pointer :: coef_ene
   double precision, dimension(:,:,:), pointer :: fine_type
   double precision, dimension(:,:,:), pointer :: coef_type
   double precision, dimension(:), pointer :: chi_all
   double precision, dimension(:,:), pointer :: chi_type
   double precision, dimension(:), allocatable :: tmpgrid,tmpgrid2
   double precision, dimension(:), pointer :: Vpot

   do idim=1,ncoup_p
      ivfinegrid(idim)=ngrid_max_p
   enddo

   ! energies
   ! during the calculation of the energie coefficents, the number of basis
   ! function might be increased until the chis are smaller than the thresholds
   basloop : do ibas=ibas_min,ibas_max
      ibas_work=ibas
      do idim=1,ncoup_p
         ivbas(idim)=ibas
      enddo
      ! if the number of basis functions is increased reallocate the memory for the coefficients
      if(ibas/=ibas_min) then
         do idim=1,ncoup_p
            deallocate(p_coef_ene(idim)%p)
         enddo
      endif
      do idim=1,ncoup_p
         allocate(p_coef_ene(idim)%p(ivbas(1)**idim,ivnsurf_p(idim)))
      enddo

      ! coefficients for energy surfaces
      do idim=1,ncoup_p
         allocate(tmpgrid(ngrid_max_p*idim))
         allocate(ifit(idim))
         call get_right_surface_stuff(idim,ires_all_p(idim))
         modes => p_modes(idim)%p
         fine_ene => fine_grid_points(:,1,:)
         coef_ene => p_coef_ene(idim)%p
         grid_total => p_igrid_total(idim)%p
         chi_all => p_chi(idim)%p
         do isurf=1,ivnsurf_p(idim)
            ! tmp grid for solve_lin_eq_kronecker
            ! structure: grid( 1:16)=x_1^1 ... x_16^1
            !            grid(17:32)=x_1^2 ... x_16^2
            ! where _i denotes the grid point and ^1 the mode
            isum=1
            do i=1,idim
               call dcopy(ivfinegrid(i),fine_grid(1,modes(i,isurf)),1,tmpgrid(isum),1)
               isum=isum+ivfinegrid(i)
               ifit(i) = ifit_fctp(modes(i,isurf),1)
            enddo
            ! calculate energy coefficients on fine grid points
            call solve_lin_eq_kronecker(idim,ivfinegrid,ivbas,tmpgrid,&
               fine_ene(:,isurf),ifit,coef_ene(:,isurf))
            ! calculate chi on fine grid points an ab initio points
            isum=0
            iprod=1
            do i=1,idim
               !ivgrid(i)=nint(dble(grid_total(isurf))**(1.0d0/dble(idim)))
               call get_ivgrid_abi_pts(idim,i,isurf,ivgrid(i))
               isum=isum+ivgrid(i)
               iprod=iprod*ivgrid(i)
            enddo
            allocate(tmpgrid2(isum)) !abi initio grid
            allocate(Vpot(iprod)) !ab initio pot
            ! extract ab initio grid information for calc_chi_kronecker
            call extract_grid_information(tmpgrid2,isum,ivgrid,idim,isurf)
            ! calculate difference surfaces for ab initio points
            call calc_difference_surface_abi(Vpot,iprod,ivgrid,idim,isurf,0,1)
            ! calculate ab initio chi
            call calc_chi_kronecker(idim,ivgrid,ivbas,tmpgrid2,Vpot,&
               ifit,coef_ene(:,isurf),chi_all(isurf))
            ! check if the chi-value is to large
            if(ibas/=ibas_max.and.chi_all(isurf)>chi_max(idim)) then
               deallocate(tmpgrid)
               deallocate(tmpgrid2)
               deallocate(ifit)
               deallocate(Vpot)
               cycle basloop
            endif
            deallocate(Vpot)
            deallocate(tmpgrid2)
         enddo
         deallocate(ifit)
         deallocate(tmpgrid)
         chi_sum = 0.0d0
         do isurf=1,ivnsurf_p(idim)
            chi_sum = chi_sum + chi_all(isurf)
         enddo
         if(chi_sum/dble(ivnsurf_p(idim))>chi_av(idim).and.ibas_work<ibas_max) cycle basloop
      enddo
      exit
   enddo basloop
   numpol_p = ibas_work
   ! calculate coefficents for
   ! dipol, polarization, mu-tensor, G matrix
   do itype=1,2
      if(itype==1)then
         ! Dipole
         ncoup_type=ncoupdip_p
         irespoly=3
         do idim=1,ncoupdip_p
            allocate(p_coef_dip(idim)%p(ibas_work**idim,3,ivnsurf_p(idim)))
         enddo
      else if(itype==2)then
         ! Polarisability
         ncoup_type=ncouppol_p
         irespoly=6
         do idim=1,ncouppol_p
            allocate(p_coef_pol(idim)%p(ibas_work**idim,6,ivnsurf_p(idim)))
         enddo
      endif
      do idim=1,ncoup_type
         ivbas(idim)=ibas_work
      enddo
      ! coefficients
      do idim=1,ncoup_type
         call get_right_surface_stuff(idim,ires_all_p(idim))
         allocate(tmpgrid(ngrid_max_p*idim))
         allocate(ifit(idim))
         modes => p_modes(idim)%p
         if(itype==1)then
            fine_type => fine_grid_points(:,2:4,:)
            coef_type => p_coef_dip(idim)%p
            chi_type => p_chi_dip(idim)%p
         elseif(itype==2)then
            fine_type => fine_grid_points(:,5:10,:)
            coef_type => p_coef_pol(idim)%p
            chi_type => p_chi_pol(idim)%p
         endif
         do isurf=1,ivnsurf_p(idim)
            ! tmp grid
            isum=1
            do i=1,idim
               call dcopy(ivfinegrid(i),fine_grid(1,modes(i,isurf)),1,tmpgrid(isum),1)
               isum=isum+ivfinegrid(i)
            enddo
            isum = 0
            iprod=1
            do i=1,idim
               !ivgrid(i)=nint(dble(grid_total(isurf))**(1.0d0/dble(idim)))
               call get_ivgrid_abi_pts(idim,i,isurf,ivgrid(i))
               isum=isum+ivgrid(i)
               iprod=iprod*ivgrid(i)
               ifit(i) = ifit_fctp(modes(i,isurf),2)
            enddo
            allocate(tmpgrid2(isum))
            allocate(Vpot(iprod)) !ab initio pot
               ! extract ab initio grid information for calc_chi_kronecker
            call extract_grid_information(tmpgrid2,isum,ivgrid,idim,isurf)
            do itl=1,irespoly
               ! calc coef
               call solve_lin_eq_kronecker(idim,ivfinegrid,ivbas,tmpgrid,&
                  fine_type(:,itl,isurf),ifit,coef_type(:,itl,isurf))
               ! calculate difference surfaces for ab initio points
               call calc_difference_surface_abi(Vpot,iprod,ivgrid,idim,isurf,itype,itl)
               ! calculate ab initio chi
               call calc_chi_kronecker(idim,ivgrid,ivbas,tmpgrid2,Vpot,&
                  ifit,coef_type(:,itl,isurf),chi_type(itl,isurf))
            enddo
            deallocate(Vpot) !ab initio pot
            deallocate(tmpgrid2) !abi initio grid
         enddo
         deallocate(ifit)
         deallocate(tmpgrid)
      enddo
   enddo

end subroutine

! ----------------------------------------------------------------------

subroutine get_ivgrid_abi_pts(idim,i,isurf,ivgrid)
   ! Defines how many points are calculated for each dimension
   use poly
   use surf_dat
   implicit none
   integer :: idim,i,isurf,ivgrid,igrid
   integer, dimension(:), pointer :: surf_adr
   double precision :: thr_grid
   double precision, dimension(:,:), pointer :: abi_grid
   thr_grid = 1.0d-8
   ivgrid = 1
   gridloop: do igrid=2,igrid_total(isurf)
      if(abs(grid_s(i,isurface_adr(isurf)+igrid-1)-grid_s(i,isurface_adr(isurf)+igrid)).lt.thr_grid) cycle gridloop
      if(abs(grid_s(i,isurface_adr(isurf)+1)-grid_s(i,isurface_adr(isurf)+igrid)).lt.thr_grid) exit gridloop
      ivgrid = ivgrid + 1
   enddo gridloop
end subroutine

! ----------------------------------------------------------------------

! extract grid information, because the ab initio grid points are all stored
subroutine extract_grid_information(tmpgrid2,nn,ivgrid,idim,isurf)
   use surf_dat
   implicit none
   integer :: nn,idim,ijump,isum,j,imode,isurf,i
   integer, dimension(idim) :: ipos,imax,ivgrid
   integer, dimension(:), pointer :: surf_adr

   double precision, dimension(:,:), pointer :: abi_grid
   double precision, dimension(nn) :: tmpgrid2

   logical :: ilog

   do i=1,idim
      ipos(i)=1
      imax(i)=2
   enddo
   ipos(1)=2
   ijump=1
   isum=0
   loop1 :do
      if(sum(ipos)==idim+1)then
         do j=1,idim
           if(ipos(j)==2) imode=j
         enddo
         do j=1,ivgrid(imode)
            tmpgrid2(isum+j)=grid_s(imode,isurface_adr(isurf)+ijump*(j-1)+1)
         enddo
         isum=isum+ivgrid(imode)
         ijump=ijump*ivgrid(imode)
      endif
      call odo(ilog,ipos,imax,idim)
      if(.not.ilog) exit loop1
   enddo loop1
end subroutine

! ----------------------------------------------------------------------

! calculate difference surfaces to check whether the fine fits are reliable
subroutine calc_difference_surface_abi(Vpot,nn,ivgrid,idim,isurf,itype,itl)
   use surf_dat
   implicit none
   integer :: idim,isurf,nn,itl,iprod,iad,icntnottwo,i,j,itype,iad_type,iad2
   integer, dimension(idim) :: m,idx,idxj
   integer, dimension(idim) :: ipos_zero
   integer, dimension(idim) :: ivgrid
!
!   double precision, dimension(:,:), pointer :: abi_grid
!   integer, dimension(:), pointer :: surf_adr
!   double precision, dimension(:,:), pointer :: abi_type
   double precision, dimension(nn) :: Vpot
   logical :: ilog,ilogj

   if(itype==0)then
      iad_type = 0
   elseif(itype==1)then
      iad_type = 1
   else
      iad_type = 4
   endif

   iprod=1
   do i=1,idim
      do j=1,ivgrid(i)
         if(abs(grid_s(i,isurface_adr(isurf)+iprod*(j-1)+1))<1.0d-8)then
            ipos_zero(i)=j
            exit
         endif
      end do
      iprod=iprod*ivgrid(i)
   end do
   do i=1,idim
      idxj(i)=1
   end do
   jloop : do
      iprod=1
      iad=1
      do i=1,idim
         iad=iad+(idxj(i)-1)*iprod
         iprod=iprod*ivgrid(i)
      end do
      Vpot(iad)=abi_points(iad_type+itl,isurface_adr(isurf)+iad)
      do i=1,idim
         m(i)=2
         idx(i)=1
      end do
      idx(1)=2
      loop : do
         icntnottwo=0
         do i=1,idim
            if(idx(i)==1) icntnottwo=icntnottwo+1
         end do
         iprod=1
         iad2=1
         do i=1,idim
            if(idx(i)==1)then
               iad2=iad2+(idxj(i)-1)*iprod
            else
               iad2=iad2+(ipos_zero(i)-1)*iprod
            endif
            iprod=iprod*ivgrid(i)
         end do
         Vpot(iad)=Vpot(iad)+(-1.0d0)**(idim-icntnottwo)*abi_points(iad_type+itl,isurface_adr(isurf)+iad2)
         call odo(ilog,idx,m,idim)
         if(.not.ilog) exit loop
      end do loop
      call odo(ilogj,idxj,ivgrid,idim)
      if(.not.ilogj) exit jloop
   end do jloop
end subroutine

! ----------------------------------------------------------------------

subroutine fit_xdim_surface(idim,ires_read)
   use surf_dat
   implicit none
   integer :: idim,ibas_use,ivbas_max_use,ires_read,isurf,i,j
   integer :: iprod,kres,iad,icntnottwo,iad2,iget_iad2
   integer, dimension(idim) :: ipos_zero
   integer, dimension(idim) :: nu_grid
   integer, dimension(idim) :: m,idx,mj,idxj
   integer, dimension(idim) :: ibas,ngrid_v
   integer, dimension(idim) :: modes

   double precision :: thr_chi
   double precision, dimension(:),allocatable :: V
   double precision, dimension(:), pointer, contiguous :: V_fine

   logical :: ilog,ilogj
   thr_chi = 1.0d-9

   allocate(chivals_f(ires_read,ivnsurf_s(idim)))
   if(ifit_auto.eq.1) then
      ivbas_max_use = ivbas_max(idim)
   else
      ivbas_max_use = ivbas_f(idim)
   endif

   ibasloop: do ibas_use=ivbas_f(idim),ivbas_max_use
      chivals_f = 0.0d0
      if(idim>1)then
         do isurf=1,ivnsurf_s(idim)
            call iget_lab_grid_s(idim,isurf,modes)
            allocate(V(igrid_total(isurf)))
            do i=1,idim
               ibas(i) = ibas_use
               ngrid_v(i)=ivgrid_f(modes(i),idim)
            end do

            allocate(V_fine(max(product(ngrid_v),igrid_total(isurf))))

            do i=1,idim
               nu_grid(i)=igrid_nabi(i,isurf)
               mj(i)=igrid_nabi(i,isurf)
            end do
            iprod=1
            do i=1,idim
               do j=1,nu_grid(i)
                  if(abs(grid_s(i,isurface_adr(isurf)+iprod*(j-1)+1))<1.0d-8)then
                     ipos_zero(i)=j
                     exit
                  endif
               end do
               iprod=iprod*nu_grid(i)
            end do
            do kres=1,ires_read
               idxj = 1
               jloop : do
                  iprod=1
                  iad=1
                  do i=1,idim
                     iad=iad+(idxj(i)-1)*iprod
                     iprod=iprod*nu_grid(i)
                  end do
                  V(iad)=abi_points(kres,isurface_adr(isurf)+iad)
                  ! calc difference surface
                  do i=1,idim
                     m(i)=2
                     idx(i)=1
                  end do
                  idx(1)=2
                  loop : do
                     icntnottwo=0
                     do i=1,idim
                        if(idx(i)==1) icntnottwo=icntnottwo+1
                     end do
                     iad2 = iget_iad2(idim,idxj,nu_grid,ipos_zero,idx)
                     V(iad)=V(iad)+(-1.0d0)**(idim-icntnottwo)*abi_points(kres,isurface_adr(isurf)+iad2)
                     call odo(ilog,idx,m,idim)
                     if(.not.ilog) exit loop
                  end do loop
                  call odo(ilogj,idxj,mj,idim)
                  if(.not.ilogj) exit jloop
               end do jloop
               !if(iprint_fine_s.le.kres) call print_Vdiff(kres,idim,isurf,V)
               call fitting_multi_kron_one(idim,ngrid_v,ibas,V_fine,isurf,kres,V,1,modes)
               if(kres.eq.1.AND.ifit_auto.eq.1) then
                  if(ibas_use.ne.ivbas_max_use.AND.chivals_f(1,isurf).gt.thr_chi) then
                     deallocate(V)
                     cycle ibasloop
                  endif
               endif
               call dcopy(product(ngrid_v),V_fine,1,fine_grid_points(:,kres,isurf),1)
            enddo
            deallocate(V_fine)
            deallocate(V)
         enddo
         ivbas_f(idim) = ibas_use
         exit ibasloop
      else
         ibas(1)=ibas_use
         ngrid_v(1) = ngrid_max_s
         call fitting_multi_kron_all(1,ngrid_max_s,ibas,fine_grid_points,ires_read)
         if(maxval(chivals_f).gt.thr_chi.AND.ifit_auto.eq.1) then
            ivbas_f(idim) = ibas_use
            cycle ibasloop
         endif
         exit ibasloop
      endif
   enddo ibasloop
   deallocate(chivals_f)

end subroutine

! ----------------------------------------------------------------------

subroutine fitting_multi_kron_all(idim,ngrid,ibas_input,V,ires_read)
   use surf_dat
   implicit none
   integer :: idim,i,ngrid,isurf,isum,iprod,ires_read,ii
   integer :: iprod1,ip,j,kres,igrid,iloop
   integer, dimension(idim) :: ipos,max_pos,ipos_grid,max_pos_grid
   integer, dimension(idim) :: modes,itgrid,itbas,ibas_input
   integer, dimension(idim,2) :: ifit

   double precision :: ddot
   double precision, dimension(ngrid**idim,ires_read,ivnsurf_s(idim)) :: V
   double precision, dimension(idim) :: hgrid,xval
   double precision, dimension(:),allocatable :: tgrid,grid_tmp
   double precision, dimension(:),allocatable :: tmp_basis,tmp_array
   double precision, dimension(:,:),allocatable :: tcoef

   logical :: ilog
   logical :: ilog_grid

   do i=1,idim
      max_pos_grid(i)=ngrid
   end do
   do isurf=1,ivnsurf_s(idim)
      call iget_lab_grid_s(idim,isurf,modes)
      isum=0
      iprod=1
      do i=1,idim
         isum=isum+igrid_nabi(i,isurf)
         itgrid(i)=igrid_nabi(i,isurf)
         if(itgrid(i)<ibas_input(i))then
            itbas(i)=itgrid(i)
         else
            itbas(i)=ibas_input(i)
         endif
         max_pos(i)=itbas(i)
         iprod=iprod*itbas(i)
      end do
      allocate(tcoef(iprod,ires_read))
      allocate(tmp_basis(iprod))
      allocate(tgrid(isum))
      allocate(tmp_array(igrid_total(isurf)))
      do i=1,idim
         iprod1=1
         isum=0
         do ip=1,i-1
            iprod1=iprod1*igrid_nabi(ip,isurf)
            isum=isum+igrid_nabi(ip,isurf)
         end do
         do j=1,igrid_nabi(i,isurf)
            tgrid(isum+j)=grid_s(i,iprod1*(j-1)+1+isurface_adr(isurf))
         end do
         hgrid(i)=(grid_s(i,iprod1*(igrid_nabi(i,isurf)-1)+1+isurface_adr(isurf))&
                   -grid_s(i,1+isurface_adr(isurf)))/(dble(ngrid-1))
         ifit(i,1) = ifit_fcts(modes(i),1)
         ifit(i,2) = ifit_fcts(modes(i),2)
      end do
      do kres=1,ires_read
         do i=1,igrid_total(isurf)
            tmp_array(i)=abi_points(kres,isurface_adr(isurf)+i)
         end do
         call solve_lin_eq_kronecker(idim,itgrid,itbas,tgrid(:),tmp_array,ifit(:,min(kres,2)),tcoef(:,kres))
         call calc_chi_kronecker(idim,itgrid,itbas,tgrid(:),tmp_array,ifit(:,min(kres,2)),tcoef(:,kres),chivals_f(kres,isurf))
      end do
      deallocate(tmp_array)
      deallocate(tgrid)
      allocate(grid_tmp(ngrid))
      igrid = 0
      do i=1,idim
         do ii=1,max_pos_grid(i)
            igrid = igrid + 1
            grid_tmp(igrid) = grid_s(i,1+isurface_adr(isurf)) + (ii-1) * hgrid(i)
         enddo
      enddo
      do kres=1,ires_read
         call calc_V_kronecker(idim,max_pos_grid,itbas,grid_tmp,V(:,kres,isurf),&
                                             ifit(:,min(kres,2)),tcoef(:,kres))
      enddo
      deallocate(grid_tmp)
      deallocate(tmp_basis) 
      deallocate(tcoef)
   end do
end subroutine

! ----------------------------------------------------------------------

subroutine fitting_multi_kron_one(idim,ngrid_v,ibas_v,V,isurf,kres,V_inn,inn,modes)
   use surf_dat
   implicit none
   integer :: idim,kres,isurf,isum,iprod,i,inn,ip,iprod1,j,igrid,ii
   integer, dimension(idim) :: ibas_v,ngrid_v,itgrid,itbas,modes,ifit

   double precision, dimension(idim) :: hgrid
   double precision, dimension(product(ngrid_v)) :: V
   double precision, dimension(igrid_total(isurf)) :: V_inn
   double precision, dimension(:),allocatable :: tgrid
   double precision, dimension(:),allocatable :: finegrid
   double precision, dimension(:),allocatable :: tcoef
   double precision, dimension(:),allocatable :: V_tmp


   isum=0
   iprod=1
   do i=1,idim
      isum=isum+igrid_nabi(i,isurf)
      itgrid(i)=igrid_nabi(i,isurf)
      if(itgrid(i)<ibas_v(i))then
         itbas(i)=itgrid(i)
      else
         itbas(i)=ibas_v(i)
      endif
      iprod=iprod*itbas(i)
   end do

   allocate(V_tmp(igrid_total(isurf)))
   if(inn==1)then !ms Only first case needed?
      do i=1,igrid_total(isurf) 
         V_tmp(i)=V_inn(i)
      enddo
   else
      do i=1,igrid_total(isurf)
         V_tmp(i)=abi_points(kres,isurface_adr(isurf)+i)
      enddo
   endif


   allocate(tcoef(iprod))
   allocate(tgrid(isum))
   do i=1,idim
      iprod1=1
      isum=0
      do ip=1,i-1
         iprod1=iprod1*igrid_nabi(ip,isurf)
         isum=isum+igrid_nabi(ip,isurf)
      end do
      do j=1,igrid_nabi(i,isurf)
         tgrid(isum+j)=grid_s(i,iprod1*(j-1)+1+isurface_adr(isurf))
      end do
      hgrid(i)=(grid_s(i,iprod1*(igrid_nabi(i,isurf)-1)+1+isurface_adr(isurf))-&
                grid_s(i,1+isurface_adr(isurf)))/(dble(ngrid_v(i)-1))
      ifit(i) = ifit_fcts(modes(i),min(kres,2))
   enddo
   call solve_lin_eq_kronecker(idim,itgrid,itbas,tgrid(:),V_tmp,ifit,tcoef(:))
   call calc_chi_kronecker(idim,itgrid,itbas,tgrid(:),V_tmp,ifit,tcoef(:),chivals_f(kres,isurf))
   deallocate(V_tmp)
   deallocate(tgrid)
   allocate(finegrid(product(ngrid_v)))
   igrid=0
   do i=1,idim
      do ii=1,ngrid_v(i)
         igrid=igrid+1
         finegrid(igrid)=grid_s(i,1+isurface_adr(isurf))+(ii-1)*hgrid(i)
      enddo
   enddo
   call calc_V_kronecker(idim,ngrid_v,itbas,finegrid,V,ifit,tcoef(:))
   deallocate(tcoef)
   deallocate(finegrid)

end subroutine

! ----------------------------------------------------------------------

subroutine print_coeffs
   use intc_info
   use poly
   use gen_dat
   implicit none
   integer :: iprint,idim,isurf,i,ipol,icnt,itype
   integer, dimension(:,:), pointer :: coords
   character(len=30) :: text30
   character(len=16) :: text16
   integer, dimension(:), pointer :: ibasis_type
   double precision :: ddummy
   double precision, dimension(:,:), allocatable :: coefs_all
   double precision, dimension(:,:), pointer :: coef_ene
   double precision, dimension(:,:,:), pointer :: coef_dip
   double precision, dimension(:,:,:), pointer :: coef_pol
   call set_ci_from
   !TODO: Add Xmin/Xmax
   
   iprint=18
   open(iprint,file=outname,action='write')
   !TODO: Print out the following informations
   !  1. General data like dimensionality, number of surfaces, number of basis functions
   write(iprint,'(/,A)')'### PARAMETERS'
   write(text30,'(A)')'Number of cart. coordinates'
   write(iprint,ci_form)text30,'=',ncart_p
   write(text30,'(A)')'Order of energy coupling'
   write(iprint,ci_form)text30,'=',ncoup_p
   write(text30,'(A)')'Order of dipole coupling'
   write(iprint,ci_form)text30,'=',ncoupdip_p
   write(text30,'(A)')'Order of polar coupling'
   write(iprint,ci_form)text30,'=',ncouppol_p
   write(text30,'(A)')'Number of basis functions'
   write(iprint,ci_form)text30,'=',numpol_p
   write(text30,'(A)')'Number of coordinates'
   write(iprint,ci_form)text30,'=',ncrd_p
   write(text30,'(A)')'Type of int. coord.'
   write(iprint,ci_form)text30,'=',itypeint
   write(iprint,'(/,A)') '### PARAMETER INTERNAL COORDINATES'
   text16 = 'Stretchings'
   write(iprint,'(4x,2A,I6)') text16,'=',nstretch
   text16 = 'Bendings'
   write(iprint,'(4x,2A,I6)') text16,'=',nbending
   text16 = 'Torsions'
   write(iprint,'(4x,2A,I6)') text16,'=',ntorsion
   write(iprint,'(/,A)') '### DEFINITION INTERNAL COORDINATES'
   write(iprint,'(4x,A)')'Stretching'
   do i=1,nstretch
      ddummy = refqvals(atstretch(3,i))/pf_intc(atstretch(3,i))
      write(iprint,'(4x,I5,I8,I8,I8,F18.12)')atstretch(3,i),atstretch(1,i),atstretch(2,i),&
                                                            atstretch(4,i),ddummy
   enddo
   write(iprint,'(4x,A)')'Bending'
   do i=1,nbending
      ddummy = refqvals(atbending(4,i))/pf_intc(atbending(4,i))
      write(iprint,'(4x,I5,I8,I8,I8,F18.12)')atbending(4,i),atbending(1,i),atbending(2,i),&
                                 atbending(3,i),ddummy
   enddo
   write(iprint,'(4x,A)')'Torsion'
   do i=1,ntorsion
      ddummy = refqvals(attorsion(5,i))/pf_intc(attorsion(5,i))
      write(iprint,'(4x,I5,I8,I8,I8,I8,F18.12)')attorsion(5,i),attorsion(1,i),attorsion(2,i),&
                                                attorsion(3,i),attorsion(4,i),ddummy
   enddo
   write(iprint,'(/,A)') '### FIT FUNCTIONS'
   do i=1,ncrd_p
      write(iprint,'(4x,I5,I5,I5)')i,ifit_fctp(i,1),ifit_fctp(i,2)
   enddo
   write(iprint,'(/,A)')'### COEFFICIENTS'
   do idim=1,ncoup_p
      call set_dim_spec_forms(idim)
      allocate(coefs_all(ires_all_p(idim),numpol_p**idim))
      coefs_all = 0.0d0
      write(iprint,'(/,A,i1,A)')'##  ',idim,'D SURFACES'
      write(text30,'(A)')'Number of Surfaces'
      write(iprint,ci_form)text30,'=',ivnsurf_p(idim)
      write(text30,'(A)')'Number of Columns'
      write(iprint,ci_form)text30,'=',ires_all_p(idim)
      coef_ene => p_coef_ene(idim)%p
      if(ires_all_p(idim).ge.2) then
         coef_dip => p_coef_dip(idim)%p
      endif
      if(ires_all_p(idim).ge.5) then
         coef_pol => p_coef_pol(idim)%p
      endif
      do isurf=1,ivnsurf_p(idim)
         coords => p_modes(idim)%p
         write(iprint,form2) '#   ',idim,'D SURFACE:',' COORD:',(coords(i,isurf),i=1,idim),'(',isurf,')'
         ! Print out the fitting function
         do ipol=1,numpol_p**idim
            coefs_all(1,ipol) = coef_ene(ipol,isurf)
         enddo
         icnt = 1
         do itype=2,min(ires_all_p(idim),4)
            do ipol=1,numpol_p**idim
               coefs_all(itype,ipol) = coef_dip(ipol,icnt,isurf)
            enddo
            icnt = icnt + 1
         enddo
         icnt = 1
         do itype=5,min(ires_all_p(idim),10)
            do ipol=1,numpol_p**idim
               coefs_all(itype,ipol) = coef_pol(ipol,icnt,isurf)
            enddo
            icnt = icnt + 1
         enddo
         ! Print coefficients
         do ipol=1,numpol_p**idim
            write(iprint,formc)ipol,coefs_all(:,ipol)
         enddo
      enddo
      deallocate(coefs_all)
   enddo
   close(iprint)
end subroutine

! ----------------------------------------------------------------------

subroutine solve_lin_eq_kronecker(idim,ivgrid,ivbas,grid,V,info_basis,coef)
   !Main Routine for the kronecker fitting
   !Input:
   !  idim: dimension of the surface V
   !  ivgrid: vector of the number of gridpoints per direction
   !  ivbas: vector of the number of basis functions per direction
   !  grid: gridpoints for surface V
   !  V: surface to approximate
   !  info_basis=1 (Gaussians)
   !            =2 (B-splines)
   !            =3 (polynomials)
   !            =4 (Morse functions)
   !            =5 (trigonometric)
   !Output:
   !  coef: resulting coefficients
   !
   !Warning: Be aware of the fact, that the first coordinate belongs to the
   !         last matrix in the Kronecker product of the small matrices!
   use gen_dat, only:iout
   implicit none
   integer :: idim,i,isum_grid,info
   integer, dimension(idim) :: ivgrid,ivbas,info_basis
   integer, dimension(idim+1) :: ngridxnbas
   integer, dimension(:), allocatable :: ipiv,iwork
   double precision, dimension(sum(ivgrid)) :: grid
   double precision, dimension(product(ivgrid)) :: V
   double precision, dimension(product(ivbas)) :: coef
   double precision, dimension(:), allocatable :: iad_pinv_small_A,iad_small_A,iad_A,iad_ATA,iad_kron

   ! calculate the small matrices Ai and their pseudoinverses
   ngridxnbas(1)=0
   do i=2,idim+1
      ngridxnbas(i)=ngridxnbas(i-1)+ivgrid(idim-i+2)*ivbas(idim-i+2)
   enddo
   allocate(iad_pinv_small_A(ngridxnbas(idim+1)))
   allocate(iad_small_A(ngridxnbas(idim+1)))
   call fzero(iad_pinv_small_A,ngridxnbas(idim+1))
   call fzero(iad_small_A,ngridxnbas(idim+1))
   isum_grid=1+sum(ivgrid)
   do i=idim,1,-1
      isum_grid=isum_grid-ivgrid(i)
      call calc_Ai(ivbas(i),ivgrid(i),info_basis(i),grid(isum_grid),&
                   iad_small_A(1+ngridxnbas(idim-i+1)))

!     pseudo inverses via (A* A)^-1 A*
      allocate(iad_A(ivbas(i)*ivgrid(i)))
      call dcopy(ivgrid(i)*ivbas(i),iad_small_A(ngridxnbas(idim-i+1)+1),&
                   1,iad_A,1)
      allocate(iad_ATA(ivbas(i)*ivbas(i)))
      call dgemm('T','N',ivbas(i),ivbas(i),ivgrid(i),1.0d0,iad_A,ivgrid(i),iad_A,ivgrid(i),0.0d0,iad_ATA,ivbas(i))
      allocate(ipiv(ivbas(i)))
      call dgetrf(ivbas(i),ivbas(i),iad_ATA,ivbas(i),ipiv,info)
      if(info.ne.0) then
         write(iout,'(1x,a)')'Error in dgetrf -> solve_lin_eq_kronecker'
         stop
      endif
      allocate(iwork(ivbas(i)*ivbas(i)))
      call dgetri(ivbas(i),iad_ATA,ivbas(i),ipiv,iwork,ivbas(i)*ivbas(i),info)
      if(info.ne.0) then
         write(iout,'(1x,a)')'Error in dgetri -> solve_lin_eq_kronecker'
         stop
      endif
      call dgemm('N','T',ivbas(i),ivgrid(i),ivbas(i),1.0d0,iad_ATA,ivbas(i),iad_A,ivgrid(i),0.0d0,&
                  iad_pinv_small_A(ngridxnbas(idim-i+1)+1),&
                  ivbas(i))
      deallocate(iad_A)
      deallocate(ipiv)
      deallocate(iwork)
      deallocate(iad_ATA)
   enddo
   ! solve linear least square problem
   if(idim==1)then
      call dgemm('N','N',ivbas(1),1,ivgrid(1),1.0d0,iad_pinv_small_A,&
                   ivbas(1),V,ivgrid(1),0.0d0,coef,ivbas(1))
   elseif(idim==2)then
      allocate(iad_kron(product(ivgrid)*product(ivbas)))
      call fzero(iad_kron,product(ivgrid)*product(ivbas))
      call krongeneral(iad_pinv_small_A,iad_pinv_small_A(ngridxnbas(2)+1),&
                       ivbas(2),ivgrid(2),ivbas(1),ivgrid(1),iad_kron)
      call dgemm('N','N',product(ivbas),1,product(ivgrid),1.0d0,iad_kron,&
           product(ivbas),V,product(ivgrid),0.0d0,coef,product(ivbas))
      deallocate(iad_kron)
   else
      call kron_mult(idim,V,iad_pinv_small_A,ivgrid,ivbas,ngridxnbas,coef)
   endif
   deallocate(iad_pinv_small_A)
   deallocate(iad_small_A)
end subroutine

! ----------------------------------------------------------------------

subroutine calc_chi_kronecker(idim,ivgrid,ivbas,grid,V,info_basis,coef,chi)
   implicit none
   integer :: idim,i,isum_grid
   integer, dimension(idim) :: ivgrid,ivbas,info_basis
   integer, dimension(idim+1) :: ngridxnbas

   double precision :: chi
   double precision, dimension(sum(ivgrid)) :: grid
   double precision, dimension(product(ivgrid)) :: V
   double precision, dimension(product(ivbas)) :: coef
   double precision, dimension(:), allocatable :: iad_small_A
   double precision, dimension(:), allocatable :: iad_Vapprox
   double precision, dimension(:), allocatable :: iad_kron

   ! calculate the small matrices Ai
   ngridxnbas(1)=0
   do i=2,idim+1
      ngridxnbas(i)=ngridxnbas(i-1)+ivgrid(idim-i+2)*ivbas(idim-i+2)
   enddo
   allocate(iad_small_A(ngridxnbas(idim+1)))
   call fzero(iad_small_A,ngridxnbas(idim+1))
   isum_grid=1+sum(ivgrid)
   do i=idim,1,-1
      isum_grid=isum_grid-ivgrid(i)
      call calc_Ai(ivbas(i),ivgrid(i),info_basis(i),grid(isum_grid),&
                   iad_small_A(ngridxnbas(idim-i+1)+1))
   enddo
   ! calculate chi
   allocate(iad_Vapprox(product(ivgrid)))
   call fzero(iad_Vapprox,product(ivgrid))
   if(idim==1)then
      call dgemm('N','N',ivgrid(1),1,ivbas(1),1.0d0,iad_small_A,&
                   ivgrid(1),coef,ivbas(1),0.0d0,iad_Vapprox,ivgrid(1))
   elseif(idim==2)then
      allocate(iad_kron(product(ivgrid)*product(ivbas)))
      call fzero(iad_kron,product(ivgrid)*product(ivbas))
      call krongeneral(iad_small_A,iad_small_A(ngridxnbas(2)+1),&
                        ivgrid(2),ivbas(2),ivgrid(1),ivbas(1),iad_kron)
      call dgemm('N','N',product(ivgrid),1,product(ivbas),1.0d0,iad_kron,&
                   product(ivgrid),coef,product(ivbas),0.0d0,iad_Vapprox,product(ivgrid))
      deallocate(iad_kron)
   else
      call kron_mult(idim,coef,iad_small_A,ivbas,ivgrid,ngridxnbas,iad_Vapprox)
   endif
   deallocate(iad_small_A)
   chi=0.0d0
   do i=1,product(ivgrid)
      chi=chi+(abs(V(i)-iad_Vapprox(i)))**2
   enddo
   deallocate(iad_Vapprox)
   chi=chi/dble(product(ivgrid))
end subroutine
