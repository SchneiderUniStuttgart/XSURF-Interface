subroutine calc_Ai(nbas,ngrid,info_basis,grid,Ai)
   !Evalulates the basisfunctions on a defined grid
   !Input:
   !  nbas: Number of basis functions
   !  grid: Vector containing all ngrid grid points
   !  info_basis: Type of basis function
   !Output:
   !  Ai: Matrix containing all of the evaluated basis functions
   implicit none
   integer :: nbas,ngrid,info_basis,j,i,ii

   double precision :: jg,ju,pi,xmin,xmax,c0,c1,a_o,b_o,ct,at,bt,d,bspline,a,cx,xval
   double precision, dimension(ngrid,nbas) :: Ai
   double precision, dimension(ngrid) :: grid
   double precision, dimension(nbas+4) :: xi_BS

   jg = 0
   ju = 0
   pi=4.0d0*atan(1.0d0)
   do j=1,nbas
      do i=1,ngrid
         if(info_basis.eq.1)then
            !Gauss functions
            xmin = grid(1)
            xmax = grid(ngrid)
            c0 = -1.0258d0
            c1 = 0.4838d0
            a_o = (c0+c1*dble(nbas))**2
            b_o = -1.0d0 + (-2.0d0+j)*2.0d0/(dble(nbas)-3.0d0)
            ct = (2.0d0*a_o/pi)**0.25d0
            at = a_o*(2.0d0/(xmax-xmin))**2
            bt = (xmax+xmin)/2.0d0+(xmax-xmin)*b_o/2.0d0
            Ai(i,j)=ct*exp(-at*(grid(i)-bt)**2)
         elseif(info_basis.eq.2)then
            !bsplines
            xmin = grid(1)
            xmax = grid(ngrid)
            d = (xmax-xmin)/(dble(nbas)-3.0d0)
            do ii = 0,nbas+3
               xi_BS(ii+1) = xmin-3.0d0*d + dble(ii)*d
            end do
            Ai(i,j)=bspline(xi_BS(j),grid(i),5)
         elseif(info_basis.eq.3)then
            !Polynominals
            Ai(i,j)=grid(i)**(j-1)
         elseif(info_basis.eq.4) then
            !Morse Function
            a=1.0d0
            xval = grid(i)
            Ai(i,j) = (1-exp(-a*xval))**dble(j-1)
         elseif(info_basis.eq.5) then
            !Trigonometric function !ms
            xval = grid(i)
            cx = 2.0d0*pi
            if(j.le.floor(dble(nbas)/2.0d0)) then
               Ai(i,j) = sin(dble(j)*xval*cx)
            else
               Ai(i,j) = 1.0d0 - cos(dble(j-floor(dble(nbas)/2.0d0))*xval*cx)
            endif
         endif
      enddo
   enddo

end subroutine

! ----------------------------------------------------------------------

subroutine krongeneral(A,B,na,ma,nb,mb,Z)
   !Computes the kronecker product of two matricees
   !Input:
   !  A: Matrix with dimension (na,ma)
   !  B: Matrix with dimension (nb,mb)
   !Output:
   !  Z: Kronecker product of A and B
   implicit none
   integer :: i,j,k,l,na,ma,nb,mb

   double precision, dimension(na,ma) :: A
   double precision, dimension(nb,mb) :: B
   double precision, dimension(na*nb,ma*mb) :: Z

   do i = 1,na
      do j = 1,ma
         do k = 1,nb
            do l = 1,mb
               Z((i-1)*nb+k,(j-1)*mb+l) = A(i,j)*B(k,l)
            end do
         end do
      end do
   end do

end subroutine

! ----------------------------------------------------------------------

recursive subroutine kron_mult(idim,V,pinv_Ai,ivgrid,ivbas,ngridxnbas,coef)
   !Smallers the Matrix pinv_Ai down to 3 dimenstions and then makes the kronecker product
   !Input:
   !  idim: Dimension
   !  V: Potential values
   !  pinv_Ai: matrix for kronecker
   !  ivgrid(idim): Number of gridpoints per dimension
   !  ivbas(idim): Number of basis functions per dimension
   !  ngridxnbas(idim+1): product of ngrid and nabs per dimension
   !Output:
   !  coef: fitted coefficients
   implicit none
   integer :: idim,i,igrid_right,igrid_left,ibas_left,ibas_right
   integer, dimension(idim) :: ivgrid,ivbas
   integer, dimension(idim+1) :: ngridxnbas
   double precision, dimension(ngridxnbas(idim+1)) :: pinv_Ai
   double precision, dimension(idim-2) :: ngridxnbas_tmp
   double precision, dimension(product(ivgrid)) :: V
   double precision, dimension(product(ivbas)) :: coef
   double precision, dimension(:), allocatable :: iad_V1,iad_kron_left,iad_kron_right

   ibas_left=1
   ibas_right=1
   igrid_left=1
   igrid_right=1
   do i=1,2
      ibas_right=ibas_right*ivbas(idim-i+1)
      igrid_right=igrid_right*ivgrid(idim-i+1)
   enddo
   do i=3,idim
      ibas_left=ibas_left*ivbas(idim-i+1)
      igrid_left=igrid_left*ivgrid(idim-i+1)
   enddo
   allocate(iad_V1(ibas_left*igrid_right))
   call fzero(iad_V1,ibas_left*igrid_right)

   if(idim==3)then
      call dgemm('N','N',ibas_left,igrid_right,igrid_left,1.0d0,&
                   pinv_Ai(ngridxnbas(3)+1),ibas_left,V,igrid_left,0.0d0,&
                   iad_V1,ibas_left)
   elseif(idim==4)then
      allocate(iad_kron_left(ibas_left*igrid_left))
      call fzero(iad_kron_left,ibas_left*igrid_left)
      call krongeneral(pinv_Ai(ngridxnbas(3)+1),pinv_Ai(ngridxnbas(4)+1),&
                     ivbas(2),ivgrid(2),ivbas(1),ivgrid(1),iad_kron_left)
      call dgemm('N','N',ibas_left,igrid_right,igrid_left,1.0d0,&
                   iad_kron_left,ibas_left,V,igrid_left,0.0d0,&
                   iad_V1,ibas_left)
      deallocate(iad_kron_left)
   else
      do i=3,idim
         ngridxnbas_tmp(i-2)=ngridxnbas(i)-ngridxnbas(3)
      enddo
      do i=1,igrid_right
         call kron_mult(idim-2,V((i-1)*igrid_left+1),pinv_Ai(ngridxnbas(3)+1),&
                        ivgrid(1),ivbas(1),ngridxnbas,iad_V1((i-1)*ibas_left))
      enddo
   endif
   allocate(iad_kron_right(igrid_right*ibas_right))
   call fzero(iad_kron_right,igrid_right*ibas_right)
   call krongeneral(pinv_Ai(1),pinv_Ai(ngridxnbas(2)+1),ivbas(idim),ivgrid(idim),&
                     ivbas(idim-1),ivgrid(idim-1),iad_kron_right)
   call dgemm('N','T',ibas_left,ibas_right,igrid_right,1.0d0,iad_V1,&
                ibas_left,iad_kron_right,ibas_right,&
                0.0d0,coef,ibas_left)
   deallocate(iad_V1)
   deallocate(iad_kron_right)
end subroutine
!==============================================================================
!====Evaluation of Surfaces ======================================================
!==============================================================================
subroutine calc_V_kronecker(idim,ivgrid,ivbas,grid,V,info_basis,coef)
   !Evaluates a fitted surface at given points
   !Input:
   !  idim: Dimension
   !  ivgrid/ivbas: Number of grid points/basis functions per dimension
   !  grid: Grid points at which the surface needs to be evaluated
   !  info_basis: Type of fit function
   !  coef: Coefficients of the fit functions
   !Output:
   !  V: The energy points at the given grid points
   implicit none
   integer :: idim,i,isum_grid
   integer, dimension(idim) :: ivgrid,ivbas,info_basis
   integer, dimension(idim+1) :: ngridxnbas
   double precision, dimension(sum(ivgrid)) :: grid
   double precision, dimension(product(ivgrid)) :: V
   double precision, dimension(product(ivbas)) ::coef
   double precision, dimension(:), allocatable :: iad_small_A,iad_kron

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
   ! Evaluate Surface
   call fzero(V,product(ivgrid))
   if(idim==1)then
      call dgemm('N','N',ivgrid(1),1,ivbas(1),1.0d0,iad_small_A,&
                  ivgrid(1),coef,ivbas(1),0.0d0,V,ivgrid(1))
   elseif(idim==2)then
      allocate(iad_kron(product(ivgrid)*product(ivbas)))
      call fzero(iad_kron,product(ivgrid)*product(ivbas))
      call krongeneral(iad_small_A,iad_small_A(ngridxnbas(2)+1),&
                        ivgrid(2),ivbas(2),ivgrid(1),ivbas(1),iad_kron)
      call dgemm('N','N',product(ivgrid),1,product(ivbas),1.0d0,iad_kron,&
                  product(ivgrid),coef,product(ivbas),0.0d0,V,product(ivgrid))
      deallocate(iad_kron)
   else
      call kron_mult(idim,coef,iad_small_A,ivbas,ivgrid,ngridxnbas,V)
   endif

   deallocate(iad_small_A)

end subroutine

! ----------------------------------------------------------------------

subroutine odo3(ilog,ix,m,n)
   implicit none
   integer :: i,n,ii
   logical :: ilog
   integer, dimension(n) :: ix
   integer, dimension(n) :: m

   i=n
   do while (i>=1)
      ix(i)=ix(i)+1
      if(ix(i)<=m(i)) then
         ilog=.TRUE.
         return
      endif
      if(i.eq.1) exit
      ix(i) = ix(i-1) + 2
      i=i-1
      do ii=i+2,n
         ix(ii) = ix(ii-1) + 1
      enddo
   end do
   ilog=.FALSE.

end subroutine

! ----------------------------------------------------------------------

subroutine get_current_elongation_intc(geom,qvalp)
   use poly
   use intc_info
   implicit none
   integer :: icrd
   double precision, dimension(ncart_p) :: geom
   double precision, dimension(ncrd_p) :: qvalp
   qvalp = 0.0d0
   call update_primitives(geom,qvalp)
   do icrd=1,ncrd_p
      qvalp(icrd) = qvalp(icrd) * pf_intc(icrd) - refqvals(icrd)
   enddo
end subroutine

! ----------------------------------------------------------------------

subroutine update_primitives(geom,qvalp)
   !Calculates the values of the primitive internal Coordiantes
   use poly
   use intc_info
   implicit none
   integer :: icrd

   double precision, dimension(ncart_p) :: geom
   double precision, dimension(ncrd_p) :: qvalp

   qvalp = 0.0d0
   !Stretchings
   do icrd=1,nstretch
      call calc_r(ncart_p,geom,atstretch(1:2,icrd),qvalp(atstretch(3,icrd)))
   enddo
   !Bendings
   do icrd=1,nbending
      call calc_angle(ncart_p,geom,atbending(1:3,icrd),qvalp(atbending(4,icrd)))
   enddo
   !Torsions
   do icrd=1,ntorsion
      call calc_torsion(ncart_p,geom,attorsion(1:4,icrd),qvalp(attorsion(5,icrd)))
   enddo
!
end subroutine

! ----------------------------------------------------------------------

subroutine calc_r(ndg,geom,iatoms,r1)
   ! Computes the distance between to atoms
   ! Input:
   !     ndg: Number of Atoms
   !     geom(ndg): geometry
   !     iatoms(2): Considered atoms
   ! Output:
   !     r1: distance between atoms from iatoms
   implicit none
   integer :: ndg
   integer, dimension(2) :: iatoms

   double precision :: r1
   double precision, dimension(ndg) :: geom
   double precision, dimension(3) :: vec1

   vec1 = geom((iatoms(2)-1)*3+1:(iatoms(2)-1)*3+3) - geom((iatoms(1)-1)*3+1:(iatoms(1)-1)*3+3)
   r1 = sqrt(sum(vec1(:)**2))

end subroutine

! ----------------------------------------------------------------------

subroutine calc_angle(ndg,geom,iatoms,phi)
   ! Computes the bond angle between 3 atoms
   ! Input:
   !     ndg: Number of Atoms
   !     geom(ndg): geometry
   !     iatoms(3): Considered atoms (iatoms(2) is central atom)
   ! Output:
   !     phi: bond angle defined by atoms from iatoms
   implicit none
   integer :: ndg
   integer, dimension(3) :: iatoms

   double precision :: phi,arg
   double precision, dimension(ndg) :: geom
   double precision, dimension(3) :: vu,vv

   ! Generate vectors
   vu = geom((iatoms(1)-1)*3+1:(iatoms(1)-1)*3+3) - geom((iatoms(2)-1)*3+1:(iatoms(2)-1)*3+3)
   vu = vu/sqrt(sum(vu(:)**2))
   vv = geom((iatoms(3)-1)*3+1:(iatoms(3)-1)*3+3) - geom((iatoms(2)-1)*3+1:(iatoms(2)-1)*3+3)
   vv = vv/sqrt(sum(vv(:)**2))
   arg = dot_product(vu,vv)
   ! Some sepcial cases due to accuracy of cosine function in Fortran
   if (arg.gt.1.0d0.AND.arg.lt.(1.0d0+1.0d-6)) arg = 1.0d0
   if (arg.lt.-1.0d0.AND.arg.gt.(-1.0d0-1.0d-6)) arg = -1.0d0
   phi = acos(arg)

end subroutine

! ----------------------------------------------------------------------

subroutine calc_torsion(ndg,geom,iatoms,tau)
   !Computes the torsion angle between 4 atoms
   ! Input:
   !     ndg: Number of Atoms
   !     geom(ndg): geometry
   !     iatoms(3): Considered atoms stored as a chain (iatoms(1) is front atom)
   ! Output:
   !     tau: torsion angle defined by atoms from iatoms
   use poly, only: pi_p
   implicit none
   integer :: ndg
   integer, dimension(4) :: iatoms

   double precision :: arg,tau,thr
   double precision, dimension(ndg) :: geom
   double precision, dimension(3) :: vu,vw,vv
   thr = 1.0d-7

   ! Generate Vectors
   vu = geom((iatoms(1)-1)*3+1:(iatoms(1)-1)*3+3) - geom((iatoms(2)-1)*3+1:(iatoms(2)-1)*3+3)
   vu = vu/sqrt(sum(vu(:)**2))
   vw = geom((iatoms(3)-1)*3+1:(iatoms(3)-1)*3+3) - geom((iatoms(2)-1)*3+1:(iatoms(2)-1)*3+3)
   vw = vw/sqrt(sum(vw(:)**2))
   vv = geom((iatoms(4)-1)*3+1:(iatoms(4)-1)*3+3) - geom((iatoms(3)-1)*3+1:(iatoms(3)-1)*3+3)
   vv = vv/sqrt(sum(vv(:)**2))

   vu = vu-dot_product(vu,vw)*vw
   vu = vu/sqrt(sum(vu(:)**2))

   vv = vv-dot_product(vv,vw)*vw
   vv = vv/sqrt(sum(vv(:)**2))
   arg = dot_product(vu,vv)

   ! Special cases
   if (arg.gt.1.0d0.AND.arg.lt.(1.0d0+1.0d-6)) arg = 1.0d0
   if (arg.lt.-1.0d0.AND.arg.gt.(-1.0d0-1.0d-6)) arg = -1.0d0

   tau = acos(arg)

   if(abs(tau).le.thr) tau = 0.0d0
   if(abs(pi_p-tau).le.thr) tau = pi_p
   if(tau.gt.thr.AND.abs(tau-pi_p).gt.thr) then
      ! Consider the direction of the torsion
      call def_sign(tau,vu,vv,vw)
   endif

end subroutine

! ----------------------------------------------------------------------

subroutine def_sign(alpha,vu,vv,axis)
   ! Define the sign of the torsion angle.
   ! If the front atom needs to be rotated clockwise the angle is positive.
   ! If the rotation is counterclockwise it the torsion is negative.
   ! Input:
   !     vu,vv: bond vector of the first and last atom
   !     axis: rotation axis
   ! Input/Output:
   !     alpha: torsion angle
   implicit none
   integer :: iax,icntp

   double precision :: alpha,thr
   double precision, dimension(3) :: vu,vv,axis,vnew
   double precision, dimension(3,3) :: rotMat

   vnew = 0.0d0
   thr = 1.0d-7

   !if(dot_product(vu,axis).gt.1.0d-8.OR.dot_product(vv,axis).gt.1.0d-8) then
   !   write(ivout,*)'Something is wrong in the calculation'
   !   write(ivout,*)'of torsinal angle. The considered vector is not'
   !   write(ivout,*)'orthogonal to the axis.'
   !   write(ivout,*)'mod_intc_head>def_sign'
   !   call fehler()
   !endif

   call get_rot_ax(alpha,axis,rotmat)
   call dgemv('N',3,3,1.0d0,rotMat,3,vu,1,0.0d0,vnew,1)
   icntp = 0
   do iax=1,3
      if(abs(vv(iax)-vnew(iax)).lt.thr) icntp = icntp + 1
   enddo
   if(icntp.eq.3) return
   vnew = 0
   call dgemv('N',3,3,1.0d0,rotMat,3,vv,1,0.0d0,vnew,1)
   icntp = 0
   do iax=1,3
      if(abs(vu(iax)-vnew(iax)).lt.thr) icntp = icntp + 1
   enddo
   if(icntp.eq.3) then
      alpha = -alpha
      return
   endif
   !write(ivout,*)'Cannot define sign of torsional angle'
   !write(ivout,*)'mod_intc_head>def_sign'
   !call fehler()

end subroutine

! ----------------------------------------------------------------------

subroutine get_rot_ax(alpha,rotax,rotmat)
   ! Computing rotational matrix
   ! Input:
   !     alpha: rotation angle
   !     rotax: rotation axis
   ! Output:
   !     rotation matrix
   implicit none
   double precision :: alpha
   double precision, dimension(3) :: rotax
   double precision, dimension(3,3) :: rotmat

   rotMat=reshape( &
   (/rotax(1)*rotax(1)*(1-cos(alpha))+cos(alpha), &
     rotax(1)*rotax(2)*(1-cos(alpha))+rotax(3)*sin(alpha), &
     rotax(1)*rotax(3)*(1-cos(alpha))-rotax(2)*sin(alpha), &
     rotax(1)*rotax(2)*(1-cos(alpha))-rotax(3)*sin(alpha), &
     rotax(2)*rotax(2)*(1-cos(alpha))+cos(alpha), &
     rotax(2)*rotax(3)*(1-cos(alpha))+rotax(1)*sin(alpha), &
     rotax(1)*rotax(3)*(1-cos(alpha))+rotax(2)*sin(alpha), &
     rotax(2)*rotax(3)*(1-cos(alpha))-rotax(1)*sin(alpha), &
     rotax(3)*rotax(3)*(1-cos(alpha))+cos(alpha)/) &
     , (/3, 3/))

end subroutine

! ----------------------------------------------------------------------

double precision recursive function bspline(xi,x,nxi) result(value)
   ! xi = knotsequence
   ! nxi = starting index
   ! the spline degree n is nxi-2

   implicit none
   integer :: nxi
   double precision :: x,d
   double precision, dimension(nxi) :: xi

   value = 0.0d0
   if (nxi.eq.2) then
      ! constant B-spline
      if (x>=xi(1).and.x<xi(nxi)) then
         value = 1.0d0
      end if
   else
      ! recursion, discarding terms with zero denominator
      d = xi(nxi-1)-xi(1)
      if (d>0) then
         value = value + (x-xi(1))*bspline(xi(1),x,nxi-1)/d
      end if
      d = xi(nxi)-xi(2)
      if (d>0) then
         value = value + (xi(nxi)-x)*bspline(xi(2),x,nxi-1)/d
      end if
   end if

   return
end function bspline

! ----------------------------------------------------------------------

subroutine fzero(array,ndim)
   implicit none
   integer :: ndim,idim
   double precision, dimension(ndim) :: array
   do idim=1,ndim
      array(idim) = 0.0d0
   enddo
end subroutine
