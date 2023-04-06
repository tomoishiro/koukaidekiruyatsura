module mdl_001_setting
   implicit none
   integer, parameter :: nrmax =xxx
   double precision, parameter :: dr = 0.1d0
   double precision, parameter :: rmax = dr*(dble(nrmax) +1.d-11)

!--- Calculational setting.
   integer, parameter :: lmax = 5 !max # of spatial ang-momenta
   integer, parameter :: ndmax = 100 !max # of nodes for s.p.sts
   integer, parameter :: nspmax = 200 !max # of s.p.sts

!--- Define constants frequently used
   complex (kind=kind(0d0)), parameter :: ai = (0.d0, 1.d0)
   double precision, parameter :: hbarc = 197.329d0
   double precision, parameter :: pi = 3.141592653d0
   double precision, parameter :: alpha = 7.2973525698d-3
   double precision, parameter :: anmass = 938.767d0

   double precision, parameter :: nc = 8.d0
   double precision, parameter :: zc = 8.d0 !proton-number of core
   double precision, parameter :: ac = nc+zc !mass-number of core
   double precision, parameter :: anmass_p = 938.2720813d0 ! proton-mass
   double precision, parameter :: anmass_n = 939.5654133d0 !neutron-mass

   double precision, parameter :: BE_core = 7.976207d0*ac !MeV, BE for O-16 nucleus.
!!   double precision, parameter :: BE_core = 8.5513064d0*ac !MeV, BE for Ca-40 nucleus.
!!   double precision, parameter :: BE_core = 8.253d0*ac !MeV, BE for Sn-100 nucleus.

   double precision, parameter :: mass_c = zc*anmass_p + nc*anmass_n - BE_core !MeV

   double precision, parameter :: xmu_cp = anmass_p*mass_c/(anmass_p + mass_c)
   double precision, parameter :: xmu_cn = anmass_n*mass_c/(anmass_n + mass_c)
!!   double precision, parameter :: xmu = xmu_cp
   double precision, parameter :: xmu = xmu_cn

   double precision, parameter :: ecut = 20.d0
   double precision, parameter :: tol = 1.d-6

!--- parameters used in C-N WS-potential
   double precision, parameter :: r00 = 1.25d0
   double precision, parameter :: a00 = 0.65d0
end module mdl_001_setting
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

module mdl_002_routines!mod_01
contains
!**************************************************************
   subroutine deriv_2nd_5p(f, ddf, dr, nr)
! The second derivative of F (5 point formula). See Abramowicz, p. 914.
!**************************************************************
      implicit none
      integer, intent(in) :: nr
      double precision, intent(in) :: dr, f(nr)
      double precision, intent(out) :: ddf(nr)
      integer :: ir
         ddf(1) = (35.d0*f(1) - 104.d0*f(2) + 114.d0*f(3) - 56.d0*f(4) + 11.d0*f(5))/12.d0/dr/dr
         ddf(2) = (11.d0*f(1) -  20.d0*f(2) +   6.d0*f(3) +  4.d0*f(4) -  1.d0*f(5))/12.d0/dr/dr
      ddf(nr-1) = (-1.d0*f(nr-4) +  4.d0*f(nr-3) +   6.d0*f(nr-2) -  20.d0*f(nr-1) + 11.d0*f(nr))/12.d0/dr/dr
        ddf(nr) = (11.d0*f(nr-4) - 56.d0*f(nr-3) + 114.d0*f(nr-2) - 104.d0*f(nr-1) + 35.d0*f(nr))/12.d0/dr/dr
      do ir = 3, nr - 2
         ddf(ir) = - f(ir-2) + 16.d0*f(ir-1) - 30.d0*f(ir) + 16.d0*f(ir+1) - f(ir+2)
         ddf(ir) = ddf(ir)/12.d0/dr/dr
      end do
      return
   end subroutine deriv_2nd_5p
!**************************************************************
   subroutine deriv_5p(f, df, dr, nr)
! The first derivative of f(x) (5 point formula). See Abramowicz, p.914.
!**************************************************************
      implicit none
      integer, intent(in) :: nr
      double precision, intent(in) :: dr, f(nr)
      double precision, intent(out) :: df(nr)
      integer :: ir
         df(1) = (-50.d0*f(1) + 96.d0*f(2) - 72.d0*f(3) + 32.d0*f(4) - 6.d0*f(5))/24.d0/dr
         df(2) = ( -6.d0*f(1) - 20.d0*f(2) + 36.d0*f(3) - 12.d0*f(4) + 2.d0*f(5))/24.d0/dr
      df(nr-1) = ( -2.d0*f(nr-4) + 12.d0*f(nr-3) - 36.d0*f(nr-2) + 20.d0*f(nr-1) + 6.d0*f(nr) )/24.d0/dr
        df(nr) = (  6.d0*f(nr-4) - 32.d0*f(nr-3) + 72.d0*f(nr-2) - 96.d0*f(nr-1) + 50.d0*f(nr))/24.d0/dr
      do ir = 3, nr - 2
         df(ir) = f(ir-2) - 8.d0*f(ir-1) + 8.d0*f(ir+1) - f(ir+2)
         df(ir) = df(ir)/12.d0/dr
      end do
      return
   end subroutine deriv_5p

!**********************************************************
   function factlog(n) result (ff)
! factlog(N)=log(N!)
      implicit none
      integer, intent (in) :: n
      integer :: i
      double precision :: ff
      ff = 0.d0
      if (n<=0) go to 500
      do i = 1, n
         ff = ff + log(dble(i))
      end do
500   return
   end function factlog

!***************************************************************
   pure function skipfact(n) result (g)
! N!! calculater by Oishi
      implicit none
      integer, intent (in) :: n
      integer :: i
      double precision :: g

      g = 1.d0
      if (n<0) go to 100
      if (n>1 .and. mod(n,2)==0) then
         do i = 2, n, 2
            g = g*dble(i)
         end do
      else if (n>1 .and. mod(n,2)==1) then
         do i = 3, n, 2
            g = g*dble(i)
         end do
      end if
      return

100   g = 0.d0
      return
   end function skipfact

!***************************************************************
   function fact(n) result (g)
! N! calculater by Oishi
      implicit none
      integer, intent (in) :: n
      integer :: i
      double precision :: g

      g = 1.d0
      if (n<2) go to 100
      do i = 2, n
         g = g*dble(i)
      end do
      return

100   if (n<0) g = 0.d0
      return
   end function fact

!***************************************************************
   pure function f_npm(n, m) result (ff) != nPm
      implicit none
      integer, intent (in) :: n, m
      integer :: i
      double precision :: ff
      ff = 1.d0
      if (m==0) go to 500
      ff = dble(n-m+1)
      do i = n - m + 2, n
         ff = ff*dble(i)
      end do
500   return
   end function f_npm

!****************************************************************
   pure function dlt(i, j) result (ff)
! Kronecker's delta_ij
      implicit none
      integer, intent (in) :: i, j
      double precision :: ff
      ff = 0.d0
      if (i==j) ff = 1.d0
      return
   end function dlt

end module mdl_002_routines
!///////////////////////////////////////////////////////////////////////////////


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
module mdl_003_subroutines

   use mdl_001_setting, only: ac, r00, a00

!--- Minnesota parameters:
   double precision, dimension(3) :: vnn = (/200.d0, -91.85d0, -178.d0/) !MeV
   double precision, dimension(3) :: bv = (/1.487d0, 0.465d0, 0.639d0/) !fm^(-2)

!--- parameters of the core-n (-p) WS(+Coulomb) potential.
   double precision, parameter :: W_0  = -53.2d0 ![MeV]
   double precision, parameter :: W_ls =  22.1d0 ![MeV*fm^2]
!<40Ca>   double precision, parameter :: W_0  = -55.57d0 ![MeV]
!<40Ca>   double precision, parameter :: W_ls =  11.28d0 ![MeV*fm^2]
   double precision, parameter :: rc0 = r00*ac**(1.d0/3.d0)
   double precision, parameter :: ac0 = a00
   double precision, parameter :: g_coul = 1.0d0 !Factor for the core-p Coulomb potential. Default is "1.d0".

contains
   function vmin_u1_S0(r) result(ff)
      implicit none
      double precision, intent(IN) :: r
      double precision :: ff
      ff = vnn(1)*exp(-bv(1)*r*r) + vnn(2)*exp(-bv(2)*r*r)
      return
   end function
   function vmin_u1_S1(r) result(ff)
      implicit none
      double precision, intent(IN) :: r
      double precision :: ff
      ff = vnn(1)*exp(-bv(1)*r*r) + vnn(3)*exp(-bv(3)*r*r)
      return
   end function


!***********************************************************
   function vcn_ws_0(l, j, r) result (ff)
      implicit none
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: v0, v0_l0, v0_lq, r0, a0, ff
      r0 = rc0 ; a0 = ac0
      v0_l0 = W_0
      v0_lq = W_0

!--- Woods-Saxon
      v0 = v0_l0
      if (l.ne.0) v0 = v0_lq
      ff = v0/(1.d0+exp((r-r0)/a0))
      return
   end function vcn_ws_0
!***********************************************************
   function vcn_ws_ls(l, j, r) result (ff) ! potential between the Core and a nucleon
      implicit none
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: xls, yl, yj, vls, vls_l0, vls_lq, df, rd, bd, ff
      rd = rc0 ; bd = ac0

      vls_l0 = W_ls*0.d0
      vls_lq = W_ls

      vls = vls_l0
      if (l.ne.0) vls = vls_lq
      yl = dble(l) ; yj = dble(j)/2.d0
!--- Woods-Saxon
      df = -exp((r-rd)/bd)/(1.d0+exp((r-rd)/bd))**2/bd
      xls = 0.5d0*(yj*(yj+1.d0)-yl*(yl+1.d0)-0.75d0)
      ff = vls*xls*df/r
      return
   end function vcn_ws_ls


!***********************************************************
   function vcn_cent(l, r) result (ff) ! potential between the Core and a nucleon
      use mdl_001_setting, only: hbarc, xmu
      implicit none
      integer, intent (in) :: l
      double precision, intent (in) :: r
      double precision :: ff
!--- centrifugal
      ff = dble(l*(l+1))*0.5d0*hbarc*hbarc/xmu/r/r
      return
   end function vcn_cent

!***********************************************************
   function vcn_coul(r) result (ff) ! potential between the Core and a nucleon
      use mdl_001_setting, only: hbarc, alpha, zc
      implicit none
      double precision, intent (in) :: r
      double precision :: vcoul, rb, ff
      rb = rc0
!--- Coulomb (only for a proton)
      if (r<rb) then
         vcoul = zc*hbarc*alpha*(3.d0-(r/rb)**2)*0.5d0/rb
      else
         vcoul = zc*hbarc*alpha/r
      end if
      ff = vcoul*g_coul !for C-p
      return
   end function vcn_coul


!***********************************************************
   function vc_prot(l, j, r) result (ff) ! potential between the Core and a proton
      implicit none
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: ff
      ff = vcn_ws_0(l, j, r) + vcn_ws_ls(l, j, r) + vcn_cent(l, r) + vcn_coul(r) !WS
!!      ff = V2b(r) + vcn_cent(l,r) + vcn_coul(r) !V2b_any
      return
   end function
!***********************************************************
   function vc_neut(l, j, r) result (ff) ! potential between the Core and a neutron
      implicit none
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: ff
      ff = vcn_ws_0(l, j, r) + vcn_ws_ls(l, j, r) + vcn_cent(l, r) !WS
!!      ff = V2b(r) + vcn_cent(l,r) !V2b_any
   end function


!***********************************************************
   function V2b(r)
      implicit none
      double precision, intent (in) :: r
      double precision :: V2b
!!      V2b = vmin_u1_S1(r)
!!      V2b = vmin_u1_S0(r)
      V2b = vc_neut(2, 3, r) -vcn_cent(2, r)
   end function V2b


!***********************************************************
   function vcx(xind,l,j,r) result (ff) ! potential between the Core and X = p or n.
      implicit none
      character(4), intent(IN) :: xind
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: ff
      ff = vc_neut(l, j, r)
      IF(xind=='prot') ff = vc_prot(l,j,r)
      IF(r .Lt. 0.0000001) ff = 7.77d9
   end function


!***********************************************************
   function vcx_cnf(xind, l, j, r) result (ff) ! Coffin-potential between the Core and a Nucleon
      implicit none
      character(4), intent(INOUT) :: xind
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: ff, vb0, rb
      rb  = 10.9d0  !8.2d0
      ff  = vcx(xind, l, j, r)
      vb0 = vcx(xind, l, j, rb)
      if (r>rb) then
         ff = vcx(xind,l,j,r) + vb0 - vcx(xind, l,j, r)
      end if
      return
   end function vcx_cnf

!**************************************************************
   subroutine hsp_cnf(xind,nsp, ll, jj, esp, psi, hsp3)
      use mdl_001_setting
      implicit none
      character(4), intent(INOUT) :: xind
      integer, intent (in) :: nsp, ll(nspmax), jj(nspmax)
      double precision, intent (in) :: esp(nspmax), psi(nspmax, nrmax)
      double precision, intent (out) :: hsp3(nspmax, nspmax)

      integer :: ir, i1, i2, l1, l2, j1, j2
      double precision :: r, dw, sum, e1, e2

!--- s.p.elements of hsp_cnf = hsp + (V_mod - V_true)
      hsp3 = 0.d0
      do i1 = 1, nsp
         j1 = jj(i1)
         l1 = ll(i1)
         e1 = esp(i1)

         do i2 = 1, nsp
            j2 = jj(i2)
            l2 = ll(i2)
            e2 = esp(i2)

            if ((j1/=j2) .or. (l1/=l2)) go to 100

            if (i1==i2) then
               hsp3(i2, i1) = hsp3(i2, i1) + e1
            end if

            sum = 0.d0
            do ir = 1, nrmax
               r = dble(ir)*dr
               dw = vcx_cnf(xind,l1, j1, r) - vcx(xind,l1, j1, r)
               sum = sum + psi(i2, ir)*dw*psi(i1, ir)*dr
            end do

            hsp3(i2, i1) = hsp3(i2, i1) + sum
100      end do
      end do

   end subroutine hsp_cnf



!**************************************************************
   subroutine sort(nsp0, ll0, jj0, node0, esp0, psi0)
! A routine to sort the eigenvalues obtained in "spbsis" and
! the associated eigen functions.
!**************************************************************
      use mdl_001_setting
      implicit none
      integer, intent (inout) :: nsp0, ll0(nspmax), jj0(nspmax), node0(nspmax)
      double precision, intent (inout) :: esp0(nspmax), psi0(nspmax, nrmax)
      integer :: i, kk, j, ir, ip
      double precision :: p

      do i = 1, nsp0 - 1
         kk = i
         p = esp0(i)
         do j = i + 1, nsp0
            if (esp0(j)<=p) then
               kk = j
               p = esp0(j)
            end if
         end do

         if (kk/=i) then
            esp0(kk) = esp0(i)
            esp0(i) = p

            ip = jj0(i)
            jj0(i) = jj0(kk)
            jj0(kk) = ip

            ip = ll0(i)
            ll0(i) = ll0(kk)
            ll0(kk) = ip

            ip = node0(i)
            node0(i) = node0(kk)
            node0(kk) = ip

            do ir = 1, nrmax
               p = psi0(i, ir)
               psi0(i, ir) = psi0(kk, ir)
               psi0(kk, ir) = p
            end do
         end if

      end do
      return
   end subroutine sort



!*****************************************************************
   subroutine spbss(xind,nsp0, ll0, jj0, node0, esp0, psi0, lin, jin)
! single-particle wave functions
!*****************************************************************
      use mdl_001_setting
      implicit none

      character(4), intent(INOUT) :: xind
      integer, intent (in) :: lin, jin
      integer, intent (inout) :: nsp0, ll0(nspmax), jj0(nspmax), node0(nspmax)
      double precision, intent (inout) :: esp0(nspmax), psi0(nspmax, nrmax)

      integer :: j, k, l, m
      integer :: node00, nodemax, inode, ic, kmax, ir
      integer :: ncore, jjc(5), llc(5), nodec(5) !orbits occupied by the core

      double precision :: e, emax, emin, xkmax
      double precision, dimension (nrmax) :: psi00

      ncore = 1 !the number of the occupied state
!      llc(1) = 0 ; jjc(1) = 1 ; nodec(1) = 0
!      llc(2) = 1 ; jjc(2) = 3 ; nodec(2) = 0
!      llc(3) = 1 ; jjc(3) = 1 ; nodec(3) = 0

      m = 0
      do l = 0, lmax
         if (l/=lin) go to 102 !filter_1
         do j = 2*l - 1, 2*l + 1, 2
            if (j/=jin) go to 101 !filter_2
            if (j<0) go to 101
            emax = ecut
            emin = -abs(vcx(xind, 0, 1, tol))
            call numerov1(xind,inode, l, j, ecut, psi00) !determine nodemax from E_cut
            nodemax = inode - 1

            do node00 = 0, nodemax

               if (node00>ndmax) go to 100

!--- remove SP-states which have the same quantum-number of core
               do ic = 1, ncore
                  if (l==llc(ic) .and. j==jjc(ic) .and. node00==nodec(ic)) go to 100
               end do

!--- formula : log_(2)(x) = log_(10)(x) / log_(10)(2)
               xkmax = log10(abs(emax-emin)/tol)/log10(2.d0) + 0.5d0
               kmax = int(xkmax)

!--- iteration to approximate `E' by both-side-attack
               do k = 1, kmax
                  e = (emin+emax)/2.d0
                  call numerov1(xind,inode, l, j, e, psi00)
                  if (inode>node00) then
                     emax = e
                  else
                     emin = e
                  end if
               end do

               m = m + 1
               if (m>nspmax) then
                  write (6, *) 'Increase ``nspmax'' (~_~)'
                  stop
               end if

!--- matching at large distance for bound state
               if ((e-vcx(xind,l,j,rmax))<0.d0) call numerov2(xind,l, j, e, psi00)

               esp0(m) = e
               ll0(m) = l
               jj0(m) = j
               node0(m) = node00
               do ir = 1, nrmax
                  psi0(m, ir) = psi00(ir)
               end do

               emax = ecut
               emin = e

100         end do
101      end do
102   end do

      nsp0 = m
      call sort(nsp0, ll0, jj0, node0, esp0, psi0)

      return
   end subroutine spbss

!**********************************************************************
   subroutine numerov1(xind,inode, l, j, e, psi00)
! Subroutine for integration of the Schroedinger eq. by Numerov method
!**********************************************************************
      use mdl_001_setting, only : rmax,nrmax,dr,ecut,hbarc,xmu,tol,lmax,ndmax,ac
      implicit none
      character(4), intent (INOUT) :: xind
      integer, intent (inout) :: inode
      integer, intent (in) :: l, j
      double precision, intent (in) :: e
      double precision, intent (inout) :: psi00(nrmax)
      integer :: ir
      double precision :: r, r0, r1, fac, psi0, psi1, dd, cc, cin

      inode = 0
      fac = dr*dr*(2.d0*xmu/hbarc/hbarc)

      psi0 = dr**(l+1)
      psi1 = (2.d0+fac*(vcx(xind,l,j,dr)-e)*5.d0/6.d0)*psi0
      psi1 = psi1/(1.d0-fac*(vcx(xind,l,j,dr+dr)-e)/12.d0)
      psi00(1) = psi0
      psi00(2) = psi1
      do ir = 3, nrmax
         r = dr*dble(ir)
         r0 = dr*dble(ir-2)
         r1 = dr*dble(ir-1)
         dd = psi0 - 2.d0*psi1
         dd = dd - fac*(vcx(xind,l,j,r0)-e)*psi0/12.d0
         dd = dd - fac*(vcx(xind,l,j,r1)-e)*psi1*5.d0/6.d0

         cc = -fac*(vcx(xind,l,j,r)-e)/12.d0 + 1.d0
         cin = 1.d0/cc

         psi00(ir) = -cin*dd

         if (ir>5 .and. psi00(ir)*psi00(ir-1)<0.d0) inode = inode + 1

         psi0 = psi1
         psi1 = psi00(ir)
      end do

!--- normalization
      fac = 0.d0
      do ir = 1, nrmax
         fac = fac + psi00(ir)*psi00(ir)*dr
      end do
      psi00 = psi00/sqrt(fac)

      return
   end subroutine numerov1

!**************************************************************
   subroutine numerov2(xind,l, j, e, psi00)
! Solve the Schroedinger eq. backward in order to ensure the
! asymptotic form for E < 0.
!*************************************************************
      use mdl_001_setting, only : rmax,nrmax,dr,ecut,hbarc,xmu,tol,lmax,ndmax,ac
!<>      use functions, only : vcx
      implicit none
      character(4), intent (INOUT) :: xind
      integer, intent (in) :: l, j
      double precision, intent (in) :: e
      double precision, intent (inout) :: psi00(nrmax)

      integer :: ir, irmatch
      double precision :: r, r0, r1, rmatch
      double precision :: evl, fac, ak, psi0, psi1, dd, cc, cin, cmatch
      double precision, dimension (:), allocatable :: psid

      allocate (psid(nrmax))

      fac = dr*dr*(2.d0*xmu/hbarc/hbarc)
      evl = vcx(xind,l, j, rmax) - e
      ak = sqrt(abs(evl)*2.d0*xmu/hbarc/hbarc)

      rmatch = 1.5d0*ac**(1.d0/3.d0)
      irmatch = int(rmatch/dr+tol)

      psi0 = exp(-ak*dble(nrmax)*dr)
      psi1 = exp(-ak*dble(nrmax-1)*dr)

      psid(nrmax) = psi0
      psid(nrmax-1) = psi1

      do ir = nrmax - 2, irmatch, -1
         r = dr*dble(ir)
         r0 = dr*dble(ir+2)
         r1 = dr*dble(ir+1)

         dd = psi0 - 2.d0*psi1
         dd = dd - fac*(vcx(xind,l,j,r0)-e)*psi0/12.d0
         dd = dd - fac*(vcx(xind,l,j,r1)-e)*psi1*5.d0/6.d0

         cc = -fac*(vcx(xind,l,j,r)-e)/12.d0 + 1.d0
         cin = 1.d0/cc

         psid(ir) = -cin*dd

         psi0 = psi1
         psi1 = psid(ir)

      end do

!--- matching
      cmatch = psi00(irmatch)/psid(irmatch)
      do ir = irmatch, nrmax
         psi00(ir) = psid(ir)*cmatch
      end do

!--- normalization
      fac = 0.d0
      do ir = 1, nrmax
         fac = fac + psi00(ir)*psi00(ir)*dr
      end do
      psi00 = psi00/sqrt(fac)

      deallocate (psid)
      return
   end subroutine numerov2

!*****************************************************************
   subroutine spbss_cnf(xind,nsp0, ll0, jj0, node0, esp0, psi0, lin, jin)
! single-particle wave functions
!*****************************************************************
      use mdl_001_setting
      implicit none
      character(4), intent(INOUT) :: xind
      integer, intent (in) :: lin, jin
      integer, intent (inout) :: nsp0, ll0(nspmax), jj0(nspmax), node0(nspmax)
      double precision, intent (inout) :: esp0(nspmax), psi0(nspmax, nrmax)
      integer :: j, k, l, m
      integer :: node00, nodemax, inode, ic, kmax, ir
      integer :: ncore, jjc(5), llc(5), nodec(5) !orbits occupied by the core
      double precision :: e, emax, emin, xkmax
      double precision, dimension (nrmax) :: psi00
      ncore = 1 !the number of the occupied state
      llc(1) = 0
      jjc(1) = 1
      nodec(1) = 0
!      llc(2) = 1 ; jjc(2) = 3 ; nodec(1) = 0
!      llc(3) = 1 ; jjc(3) = 1 ; nodec(1) = 0

      m = 0
      do l = 0, lmax
         if (l/=lin) go to 102 !filter_1
         do j = 2*l - 1, 2*l + 1, 2
            if (j/=jin) go to 101 !filter_2
            if (j<0) go to 101
            emax = ecut
            emin = -abs(vcx_cnf(xind,0,1,tol))
            call numerov1_cnf(xind,inode, l, j, ecut, psi00) !determine nodemax from E_cut
            nodemax = inode - 1

            do node00 = 0, nodemax

               if (node00>ndmax) go to 100

!--- remove SP-states which have the same quantum-number of core
               do ic = 1, ncore
                  if (l==llc(ic) .and. j==jjc(ic) .and. node00==nodec(ic)) go to 100
               end do

!--- formula : log_(2)(x) = log_(10)(x) / log_(10)(2)
               xkmax = log10(abs(emax-emin)/tol)/log10(2.d0) + 0.5d0
               kmax = int(xkmax)

!--- iteration to approximate `E' by both-side-attack
               do k = 1, kmax
                  e = (emin+emax)/2.d0
                  call numerov1_cnf(xind,inode, l, j, e, psi00)

                  if (inode>node00) then
                     emax = e
                  else
                     emin = e
                  end if
               end do

               m = m + 1
               if (m>nspmax) then
                  write (6, *) 'Increase ``nspmax'' (~_~)'
                  stop
               end if

!--- matching at large distance for bound state
               if ((e-vcx_cnf(xind,l,j,rmax))<0.d0) call numerov2_cnf(xind,l, j, e, psi00)

               esp0(m) = e
               ll0(m) = l
               jj0(m) = j
               node0(m) = node00
               do ir = 1, nrmax
                  psi0(m, ir) = psi00(ir)
               end do

               emax = ecut
               emin = e

100         end do
101      end do
102   end do

      nsp0 = m
      call sort(nsp0, ll0, jj0, node0, esp0, psi0)
      return
   end subroutine spbss_cnf

!**********************************************************************
   subroutine numerov1_cnf(xind,inode, l, j, e, psi00)
! Subroutine for integration of the Schroedinger eq. by Numerov method
!**********************************************************************
      use mdl_001_setting
      implicit none
      character(4), intent(INOUT) :: xind
      integer, intent (inout) :: inode
      integer, intent (in) :: l, j
      double precision, intent (in) :: e
      double precision, intent (inout) :: psi00(nrmax)

      integer :: ir
      double precision :: r, r0, r1, fac, psi0, psi1, dd, cc, cin

      inode = 0
      fac = dr*dr*(2.d0*xmu/hbarc/hbarc)
      psi0 = dr**(l+1)
      psi1 = (2.d0+fac*(vcx_cnf(xind,l,j,dr)-e)*5.d0/6.d0)*psi0
      psi1 = psi1/(1.d0-fac*(vcx_cnf(xind,l,j,dr+dr)-e)/12.d0)

      psi00(1) = psi0
      psi00(2) = psi1

      do ir = 3, nrmax
         r = dr*dble(ir)
         r0 = dr*dble(ir-2)
         r1 = dr*dble(ir-1)

         dd = psi0 - 2.d0*psi1
         dd = dd - fac*(vcx_cnf(xind,l,j,r0)-e)*psi0/12.d0
         dd = dd - fac*(vcx_cnf(xind,l,j,r1)-e)*psi1*5.d0/6.d0

         cc = -fac*(vcx_cnf(xind,l,j,r)-e)/12.d0 + 1.d0
         cin = 1.d0/cc

         psi00(ir) = -cin*dd

         if (ir>5 .and. psi00(ir)*psi00(ir-1)<0.d0) inode = inode + 1

         psi0 = psi1
         psi1 = psi00(ir)
      end do

!--- normalization
      fac = 0.d0
      do ir = 1, nrmax
         fac = fac + psi00(ir)*psi00(ir)*dr
      end do
      psi00 = psi00/sqrt(fac)

      return
   end subroutine numerov1_cnf

!**************************************************************
   subroutine numerov2_cnf(xind,l, j, e, psi00)
! Solve the Schroedinger eq. backward in order to ensure the
! asymptotic form for E < 0.
!*************************************************************
      use mdl_001_setting
      implicit none
      character(4), intent(INOUT) :: xind
      integer, intent (in) :: l, j
      double precision, intent (in) :: e
      double precision, intent (inout) :: psi00(nrmax)

      integer :: ir, irmatch
      double precision :: r, r0, r1, rmatch
      double precision :: evl, fac, ak, psi0, psi1, dd, cc, cin, cmatch
      double precision, dimension (:), allocatable :: psid

      allocate (psid(nrmax))

      fac = dr*dr*(2.d0*xmu/hbarc/hbarc)
      evl = vcx_cnf(xind,l, j, rmax) - e
      ak = sqrt(abs(evl)*2.d0*xmu/hbarc/hbarc)

      rmatch = 1.5d0*ac**(1.d0/3.d0)
      irmatch = int(rmatch/dr+tol)

      psi0 = exp(-ak*dble(nrmax)*dr)
      psi1 = exp(-ak*dble(nrmax-1)*dr)

      psid(nrmax) = psi0
      psid(nrmax-1) = psi1

      do ir = nrmax - 2, irmatch, -1

         r = dr*dble(ir)
         r0 = dr*dble(ir+2)
         r1 = dr*dble(ir+1)

         dd = psi0 - 2.d0*psi1
         dd = dd - fac*(vcx_cnf(xind,l,j,r0)-e)*psi0/12.d0
         dd = dd - fac*(vcx_cnf(xind,l,j,r1)-e)*psi1*5.d0/6.d0

         cc = -fac*(vcx_cnf(xind,l,j,r)-e)/12.d0 + 1.d0
         cin = 1.d0/cc

         psid(ir) = -cin*dd

         psi0 = psi1
         psi1 = psid(ir)

      end do

!--- matching
      cmatch = psi00(irmatch)/psid(irmatch)
      do ir = irmatch, nrmax
         psi00(ir) = psid(ir)*cmatch
      end do

!--- normalization
      fac = 0.d0
      do ir = 1, nrmax
         fac = fac + psi00(ir)*psi00(ir)*dr
      end do
      psi00 = psi00/sqrt(fac)

      deallocate (psid)
      return
   end subroutine numerov2_cnf

end module mdl_003_subroutines
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


!***********************************************************
   subroutine CN_TD(sind, lsp,jsp)
   use mdl_002_routines
   use mdl_001_setting
   use mdl_003_subroutines
   implicit none
      character(4), intent(IN) :: sind
      integer, intent(IN) :: lsp, jsp
   double precision, parameter :: tmax = 10000.d0
   double precision, parameter :: dt = 1.d0
   integer, parameter :: ntmax = int(tmax/dt+1.d-6) !max # of t-grid
   !---
   integer :: i, j, k, l, m, n, ir, it, idate(8), nreso, ilab, IF_TEMPO
   integer :: nsp1, ll1(nspmax), jj1(nspmax), node1(nspmax)
   integer :: nsp2, ll2(nspmax), jj2(nspmax), node2(nspmax)
   double precision :: e, r, t, gamma
   double precision :: esp1(nspmax), psi1(nspmax, nrmax)
   double precision :: esp2(nspmax), psi2(nspmax, nrmax)
   double precision :: edy3(nspmax), q3(nspmax), psi3(nspmax, nrmax)
   double precision, dimension (nspmax, nspmax) :: hsp3, ceff3
   double precision, dimension (:), allocatable :: psdum, dpsdum
   double precision, dimension (:, :), allocatable :: psurv
   complex (kind=kind(0d0)) :: z1, zf3(nrmax), zd3(nrmax)
   complex (kind=kind(0d0)), dimension (:, :), allocatable :: cft3, cdt3, cfi3
   character(4) :: xind
   xind = sind
   call rectime(21)
   l = lsp ; j = jsp ; ilab = 1000 ; IF(xind == 'neut') ilab = 2000
   do ir = 0, nrmax
      r = dble(ir)*dr
      write (ilab+1, '(F7.3, 6(2X,G11.4))') r, vcx(xind,l,j,r), vcn_cent(l,r), vcn_coul(r), V2b(r)
   end do
   psi1 = 0.d0 ; esp1 = 0.d0 ; node1 = 0 ; ll1 = 0 ; jj1 = 0
   call spbss(xind,nsp1, ll1, jj1, node1, esp1, psi1, l, j)
   write (6, *) 'Single Particle States [1]'
   write (6, *) '   # Node    L    J*2      E(MeV)'
   write (6, *) ''
   write (8, *) 'Single Particle states [1]'
   write (8, *) '   # Node    L    J*2      E(MeV)'
   write (8, *) ''
   do i = 1, 8 !nsp1-1
      n = node1(i)
      e = esp1(i)
      write (6, '(4(2X,I3),F15.7)') i, n, l, j, e
      write (8, '(4(2X,I3),F15.7)') i, n, l, j, e
   end do
   i = nsp1
   write (6, '(4(2X,I3),F15.7)') i, node1(i), ll1(i), jj1(i), esp1(i)
   write (6, *) 'number of s.p.basis=', nsp1
   write (6, *) ''
   write (8, '(4(2X,I3),F15.7)') i, node1(i), ll1(i), jj1(i), esp1(i)
   write (8, *) 'number of s.p.basis=', nsp1
   write (8, *) ''
!<><><><><>
   open(8492, file="Rmax_and_Energies.txt", action="write", status="replace")
   write (8492, '( 7(X,F18.9) )') rmax, esp1(1),esp1(2),esp1(3),esp1(4),esp1(5),esp1(6)
   close(8492)
!<><><><><>

!--- TD calculation was eliminated.
   IF_TEMPO = 0
   IF(IF_TEMPO.eq.0) go to 87314
   87314 call rectime(29)
   write (6, *) 'The program terminated normally (^_^)'
   write (8, *) 'The program terminated normally (^_^)'
   end subroutine CN_TD

!***********************************************************
   subroutine rectime(id)
   implicit none
   integer, intent(in) :: id
   integer, dimension (8) :: md
   call date_and_time(values = md)
   write (6, *) '' ; write (8, *) ''
   write (6, '("*** time -",i2,2x,i4,"/",i2,"/",i2,", ",i2,":",i2,":",i2," ***")') id, md(1),md(2),md(3),md(5),md(6),md(7)
   write (8, '("*** time -",i2,2x,i4,"/",i2,"/",i2,", ",i2,":",i2,":",i2," ***")') id, md(1),md(2),md(3),md(5),md(6),md(7)
   write (6, *) '' ; write (8, *) ''
   end subroutine rectime


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ CN-MAIN
program MAIN_clock
   implicit none
   real :: c1, c2
   integer :: it1, it2, it_rate, it_max, idiff, l, jj
   open(8, file='Clock_Core_Nucleon.dat', action='write')
   call cpu_time( c1 )
   call system_clock(it1)

   !--- For "proton or neutron", change also the "xmu" at top.
   if(1.eq.0) then
   write(6,*) "**************************************************"
   write(6,*) "************* PROTON *****************************"
   write(6,*) "**************************************************"
      l = 2 ; jj = 3!---------------------------------filter
      call CN_TD('prot',l,jj)
   end if
   !----------------------------
   if(2.eq.2) then
   write(6,*) "**************************************************"
   write(6,*) "************* NEUTRON ****************************"
   write(6,*) "**************************************************"
      l = 2 ; jj = 3!---------------------------------filter
      call CN_TD('neut',l,jj)
   end if
   !----------------------------

   call system_clock(it2, it_rate, it_max)
   if ( it2 < it1 ) then
     idiff = (it_max - it1) + it2 + 1
   else
     idiff = it2 - it1
   endif
   call cpu_time( c2 )
   write(6, *) "system time :", idiff/it_rate, "seconds."
   write(6, *) "cpu time    :", c2-c1, "seconds."
   open( 8, file='Clock_Core_Nucleon.dat', action='write')
   write(8, *) "system time :", idiff/it_rate, "seconds."
   write(8, *) "cpu time    :", c2-c1, "seconds."
   close(8)
end program MAIN_clock


