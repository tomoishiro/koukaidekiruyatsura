!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
module mdl_001_setting
   implicit none

!--- Calculational setting.
   integer, parameter :: lmax = 5 !max # of spatial ang-momenta
   integer, parameter :: ndmax = 100 !max # of nodes for s.p.sts
   integer, parameter :: nspmax = 200 !max # of s.p.sts
   double precision, parameter :: rmax = 240.d0
   double precision, parameter :: dr = 0.1d0
   integer, parameter :: nrmax = int(rmax/dr+1.d-6) !max # of r-grid

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

module mdl_002_routines
contains
!**************************************************************
   subroutine deriv_2nd_5p(f, ddf, dr, nr)! The second derivative of F (5 point formula). See Abramowicz, p. 914.
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
   subroutine deriv_5p(f, df, dr, nr)!The first derivative of f(x) (5 point formula). See Abramowicz, p.914.
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
!*****************************************************************
   subroutine diajac(A, D, VV, np)
      implicit none
      integer, intent(IN) :: np
      double precision, intent(INOUT) :: A(np, np), D(np), VV(np, np)
      integer :: n, nrot, ip, iq, i, j
      double precision :: c, gx, hx, s, sm, t, tau, theta, tresh
      double precision, dimension (np) :: b, z
! The Jacobi method to compute all eigenvalues and eigenvectors of a real symmetric matrix A.
! On output, the elements of A above the diagonal line are destroyed.
! The I-th component of the eigenvector for the K-th eigenvalue is given by VV(I,K).
! D returns the eigenvalues of A in its first NP elements.
! NROT returns the number of Givens-Jacobi rotations required.
      n=np
      do ip = 1, n
         do iq = 1, n
            VV(ip, iq) = 0.d0
         end do
         VV(ip, ip) = 1.d0
      end do
      do ip = 1, n
         b(ip) = A(ip, ip)
         D(ip) = b(ip)
         z(ip) = 0.d0
      end do
      nrot = 0
      do 1000 i = 1, 50
         sm = 0.d0
         do ip = 1, n - 1
            do iq = ip + 1, n
               sm = sm + abs(A(ip,iq))
            end do
         end do
         if (sm==0.d0) then
            call eigsrt(D, VV, n)
            return
         end if
         if (i<4) then
            tresh = 0.2d0*sm/n**2
         else
            tresh = 0.d0
         end if
         do 286 ip = 1, n - 1
            do 287 iq = ip + 1, n
               gx = 100.d0*abs(A(ip,iq))
               !-----------------------------------------------------------------------<UTA starts.>
               if ((i>4) .and. (abs(D(ip))+gx==abs(D(ip))) .and. (abs(D(iq))+gx==abs(D(iq)))) then
                  A(ip, iq) = 0.d0
               else if (abs(A(ip,iq))>tresh) then
                  hx = D(iq) - D(ip)
                  if (abs(hx)+gx==abs(hx)) then
                     t = A(ip, iq)/hx
                  else
                     theta = 0.5d0*hx/A(ip, iq)
                     t = 1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                     if (theta<0.d0) t = -t
                  end if
                  c = 1.d0/sqrt(1.d0+t**2)
                  s = t*c
                  tau = s/(1.d0+c)
                  hx = t*A(ip, iq)
                  z(ip) = z(ip) -hx
                  z(iq) = z(iq) +hx
                  D(ip) = D(ip) -hx
                  D(iq) = D(iq) +hx
                  A(ip, iq) = 0.d0
                  do j = 1, ip - 1
                     gx = A(j, ip)
                     hx = A(j, iq)
                     A(j, ip) = gx -s*(hx +gx*tau)
                     A(j, iq) = hx +s*(gx -hx*tau)
                  end do
                  do j = ip + 1, iq - 1
                     gx = A(ip, j)
                     hx = A(j, iq)
                     A(ip, j) = gx - s*(hx +gx*tau)
                     A(j, iq) = hx + s*(gx -hx*tau)
                  end do
                  do j = iq + 1, n
                     gx = A(ip, j)
                     hx = A(iq, j)
                     A(ip, j) = gx - s*(hx +gx*tau)
                     A(iq, j) = hx + s*(gx -hx*tau)
                  end do
                  do j = 1, n
                     gx = VV(j, ip)
                     hx = VV(j, iq)
                     VV(j, ip) = gx - s*(hx +gx*tau)
                     VV(j, iq) = hx + s*(gx -hx*tau)
                  end do
                  nrot = nrot + 1
               end if
               !----------------------------------------------------------------<UTA ends.>
     287    end do
  286    end do
         !---
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            D(ip) = b(ip)
            z(ip) = 0.d0
         end do
1000  end do
!<>      pause '50 iterations should never happen'
      return
   end subroutine diajac
!**********************************************************
   subroutine eigsrt(d, v, n)
      implicit none
      integer :: n, i, j, k
      double precision :: p
      double precision, dimension (n) :: d
      double precision, dimension (n, n) :: v
      do i = 1, n - 1
         k = i
         p = d(i)
         do j = i + 1, n
            if (d(j)<=p) then
               k = j
               p = d(j)
            end if
         end do
         if (k/=i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j, i)
               v(j, i) = v(j, k)
               v(j, k) = p
            end do
         end if
      end do
      return
   end subroutine eigsrt
!**********************************************************
   function factlog(n) result (ff)! factlog(N)=log(N!)
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
   pure function skipfact(n) result (g)! N!! calculater by Oishi
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
   function fact(n) result (g)! N! calculater by Oishi
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
   pure function dlt(i, j) result (ff)! Kronecker's delta_ij
      implicit none
      integer, intent (in) :: i, j
      double precision :: ff
      ff = 0.d0
      if (i==j) ff = 1.d0
      return
   end function dlt
!****************************************************************
end module mdl_002_routines
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
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
      rb  = 11.2d0  !8.2d0 !for confining.
      ff  = vcx(xind, l, j, r)
      vb0 = vcx(xind, l, j, rb)
      if (r>rb) then
         ff = vcx(xind,l,j,r) + vb0 - vcx(xind, l,j, r)
      end if
      return
   end function vcx_cnf
!**************************************************************
   subroutine coulscat(xind, l, j, ein, zs, zf)
! subroutine for integration of the schroedinger eq. by
! modified numerov method. Note; s.p.W.F. is not saved.
!<>      use mdl_006_cwf
      use mdl_001_setting, only: hbarc, xmu, zc, tol, alpha, ai, dr, rmax, nrmax
      implicit none
      character (4), intent (inout) :: xind
      integer, intent (in) :: l, j
      double precision, intent (in) :: ein
      complex (kind=kind(0d0)), intent (out) :: zf, zs
      integer :: iterat, ir
      integer, dimension (0:l+1) :: iexp
      double precision :: rmind, rmaxd, r, ri, r1, fac, ak, eta, rho, vcb
      double precision :: psi, psi0, psi1, xi, xi0, xi1, dd, cc
      double precision, dimension (0:l+1) :: fcw, gcw, fpcw, gpcw, sigmad
      complex (kind=kind(0d0)) :: cwup0, cwdown0, cwup1, cwdown1, bb, bb2
      rmind = 0.d0
      rmaxd = rmax
      iterat = int((rmaxd-rmind)/dr+tol)
      !##########################################<2023.04.13> Eliminated.
      return
   end subroutine coulscat
!**************************************************************
   subroutine hspdiag(n, h, eall, veigs)! Diagonalize the s.p.Hamiltonian matrix
      use mdl_002_routines, only : diajac
      use mdl_001_setting
      implicit none
      integer :: i, j
      integer, intent (in) :: n
      double precision, intent (in) :: h(nspmax, nspmax)
      double precision, intent (out) :: eall(nspmax), veigs(nspmax, nspmax)
      double precision, dimension (:), allocatable :: d
      double precision, dimension (:, :), allocatable :: v, h0
      allocate (d(n), v(n,n), h0(n,n))
      do i = 1, n
         d(i) = 0.d0
         do j = 1, n
            v(j, i) = 0.d0
            h0(j, i) = h(j, i)
         end do
      end do
      call diajac(h0, d, v, n)
      do i = 1, n
         eall(i) = d(i)
         do j = 1, n
            veigs(j, i) = v(j, i) !c_{M-th}(k)=veigs(k,M)
         end do
      end do
      deallocate (d, v, h0)
      return
   end subroutine hspdiag
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
! A routine to sort the eigenvalues obtained in "spbsis" and the associated eigen functions.
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
   subroutine spbss(xind,nsp0, ll0, jj0, node0, esp0, psi0, lin, jin)! single-particle wave functions
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
! Solve the Schroedinger eq. backward in order to ensure the asymptotic form for E < 0.
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
   subroutine numerov1_cnf(xind,inode, l, j, e, psi00)! Subroutine for integration of the Schroedinger eq. by Numerov method
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

!--- Time-dependent calculation starts.
   IF_TEMPO = 1
   IF(IF_TEMPO.eq.0) go to 87314

   psi2 = 0.d0 ; esp2 = 0.d0 ; node2 = 0 ; ll2 = 0 ; jj2 = 0
   l = lsp ; j = jsp
   call spbss_cnf(xind,nsp2, ll2, jj2, node2, esp2, psi2, l, j)
   write (6, *) 'Single Particle States [2], Confined for TD Calc.'
   write (6, *) '   # Node    L    J*2      E(MeV)'
   write (6, *) ''
   write (8, *) 'Single Particle states [2], Confined for TD Calc.'
   write (8, *) '   # Node    L    J*2      E(MeV)'
   write (8, *) ''
   do i = 1, 8 !nsp2-1
      n = node2(i) !l = ll2(i) ; j = jj2(i)
      e = esp2(i)
      write (6, '(4(2X,I3),F15.7)') i, n, l, j, e
      write (8, '(4(2X,I3),F15.7)') i, n, l, j, e
   end do

   i = nsp2
   write (6, '(4(2X,I3),F15.7)') i, node2(i), ll2(i), jj2(i), esp2(i)
   write (6, *) 'number of s.p.basis=', nsp2
   write (6, *) ''
   write (8, '(4(2X,I3),F15.7)') i, node2(i), ll2(i), jj2(i), esp2(i)
   write (8, *) 'number of s.p.basis=', nsp2
   write (8, *) ''
   call rectime(22)

   call hsp_cnf(xind,nsp1, ll1, jj1, esp1, psi1, hsp3)
   call hspdiag(nsp1, hsp3, edy3, ceff3) !a_{A-th}(i)=c(i,A)
   psi3 = 0.d0
   q3 = 0.d0
   do i = 1, nsp1
      do j = 1, nsp1
         q3(i) = q3(i) + ceff3(j, i)*ceff3(j, i)*esp1(j)
         do ir = 1, nrmax
            psi3(i, ir) = psi3(i, ir) + ceff3(j, i)*psi1(j, ir)
         end do
      end do
   end do
   write (6, *) 'Expan. of h_2 via sp_1 sts.'
   write (6, *) '  #      Esp(num)    Esp(hdiag)   $2-$3    Q=<h_true>'
   write (8, *) 'Expan. of h_2 via sp_1 sts.'
   write (8, *) '  #      Esp(num)    Esp(hdiag)   $2-$3    Q=<h_true>'
   do i = 1, 1 !nsp1
      write (6, '(X,I3,4(X,F11.4))') i, edy3(i), esp2(i), edy3(i) - esp2(i), q3(i)
      write (8, '(X,I3,4(X,F11.4))') i, edy3(i), esp2(i), edy3(i) - esp2(i), q3(i)
      l = ll2(i) ; j = jj2(i) ; e = esp2(i)
      do ir = 1, nrmax
         r = dble(ir)*dr
         write (ilab+10+i, '(F6.2,X,F8.4,3(X,E12.5))') r, e, vcx_cnf(xind,l, j, r), psi2(i, ir)*psi2(i, ir), psi3(i, ir)*psi3(i, ir)
      end do
   end do
   write (6, *) ''
   write (8, *) ''
!--- Note: "psi2" is not used anymore.

   call rectime(23)
!--- preparation for time-evolution
   write (6, '(X,A6,X,F9.4)') 'tmax =', tmax
   write (6, '(X,A4,X,F5.2)') 'dt =', dt
   write (6, '(X,A7,X,I6)') 'ntmax =', ntmax
   write (6, *) ''
   write (8, '(X,A6,X,F9.4)') 'tmax =', tmax
   write (8, '(X,A4,X,F5.2)') 'dt =', dt
   write (8, '(X,A7,X,I6)') 'ntmax =', ntmax
   write (8, *) ''

   10009 allocate (psurv(0:ntmax,9), psdum(1:ntmax+1), dpsdum(1:ntmax+1), &
                   cfi3(nsp1,nsp1), cft3(nsp1,nsp1), cdt3(nsp1,nsp1))
   psurv = 1.d0
   psdum = 0.d0
   dpsdum = 0.d0

!--- initial wave functions (including resonance).
   do j = 1, nsp1 !cfi3_{reso'}(i) = veigs(i,N_reso') for t=0.
      do i = 1, nsp1 !N_reso' is in case by case.
         cfi3(i, j) = cmplx(ceff3(i,j), 0.d0, kind(0d0)) !via true basis
      end do
   end do

!### TE1 #######################################################
   do it = 1, ntmax
      t = dble(it)*dt
!--- {cft3} at t=it*dt
      do k = 1, 6 !nsp1
         do i = 1, nsp1
            cft3(i, k) = exp(-ai*t*esp1(i)/hbarc)*cfi3(i, k)
         end do
      end do
!--- suvival probability
      do k = 1, 6 !nsp1
         z1 = (0.d0, 0.d0)
         do i = 1, nsp1
            z1 = z1 + conjg(cft3(i,k))*cfi3(i, k)
         end do
         psurv(it, k) = abs(conjg(z1)*z1)
      end do
   end do
!### E. of TE1 ################################################
!### TE2 #######################################################
   m = -1 ; nreso = 1
   do it = 0, 8000, 2000
      t = dble(it)*dt ; m = m + 1
!--- {cft3} at t=it*dt, only for nreso
      do i = 1, nsp1
         cft3(i, nreso) = exp(-ai*t*esp1(i)/hbarc)*cfi3(i, nreso)
      end do
!--- W.F.(t) & z1 = <W.F.(0)|W.F.(t)>
      zf3 = 0.d0 ; z1 = (0.d0, 0.d0)
      do j = 1, nsp1
         do ir = 1, nrmax
            zf3(ir) = zf3(ir) + cft3(j, nreso)*psi1(j, ir)
         end do
         z1 = z1 + conjg(cfi3(j, nreso))*cft3(j, nreso)
      end do
!--- {cdt3} at t=it*dt for the decay-state, only for nreso
      if (it.ne.0) then
      do i = 1, nsp1
         cdt3(i, nreso) = cft3(i, nreso) - z1*cfi3(i, nreso)
      end do
!--- decayed W.F.(t) & z1 = Nd(t), only for nreso
      zd3 = 0.d0 ; z1 = (0.d0, 0.d0)
      do j = 1, nsp1
         do ir = 1, nrmax
            zd3(ir) = zd3(ir) + cdt3(j, nreso)*psi1(j, ir)
         end do
         z1 = z1 + conjg(cdt3(j, nreso))*cdt3(j, nreso)
      end do
!      zd3 = zd3/sqrt(abs(z1)) !normalization
      end if
      do ir = 1, nrmax
         write (ilab+500+m, '(F6.2,X,4(X,E12.5))') dble(ir)*dr, abs(conjg(zf3(ir))*zf3(ir)), abs(conjg(zd3(ir))*zd3(ir))
      end do
   end do
!### E. of TE2 ################################################
!--- write ``gamma''
   do k = 1, 1 !nsp1
      psdum = 0.d0 ; dpsdum = 0.d0
      do it = 0, ntmax
         psdum(it+1) = psurv(it, k)
      end do
      call deriv_5p(psdum, dpsdum, dt, ntmax+1)
      do it = 0, ntmax
         t = dble(it)*dt
         gamma = -hbarc*(dpsdum(it+1)/psdum(it+1))
         write (ilab+20+k, '(F12.5,2(X,G17.8))') t, psdum(it+1), gamma
      end do
   end do
   20009 deallocate (psurv, psdum, dpsdum, cfi3, cft3, cdt3)

   87314 call rectime(29)
   write (6, *) 'The program terminated normally (^_^)'
   write (8, *) 'The program terminated normally (^_^)'
   end subroutine CN_TD

!***********************************************************
   subroutine  CN_SC(sind)
   use mdl_002_routines
   use mdl_001_setting, only : ai, tol, pi
   use mdl_003_subroutines
   implicit none
   character(4), intent(IN) :: sind
   double precision, parameter :: de   = 0.001d0
   double precision, parameter :: emin = 0.3d0 !MeV
   double precision, parameter :: emax = 0.6d0 !MeV
   integer, parameter :: nemax = int((emax-emin)/de+tol)
   integer :: ie, idate(8), il
   double precision :: e
   double precision, dimension (0:nemax) :: delta_s1, dedel_s1
   double precision, dimension (0:nemax) :: delta_p3, dedel_p3, delta_p1, dedel_p1
   double precision, dimension (0:nemax) :: delta_d5, dedel_d5, delta_d3, dedel_d3
   double precision, dimension (0:nemax) :: delta_f7, dedel_f7, delta_f5, dedel_f5
   complex (kind=kind(0d0)) :: z1, facn, smat
   character(4) :: xind
   xind = sind
   call rectime(11)
   delta_s1 = 0.d0 ; dedel_s1 = 0.d0
   delta_p3 = 0.d0 ; dedel_p3 = 0.d0
   delta_p1 = 0.d0 ; dedel_p1 = 0.d0
   delta_d5 = 0.d0 ; dedel_d5 = 0.d0
   delta_d3 = 0.d0 ; dedel_d3 = 0.d0
   delta_f7 = 0.d0 ; dedel_f7 = 0.d0
   delta_f5 = 0.d0 ; dedel_f5 = 0.d0

!--- calculate phase-shift as a function of E
   do 1000 ie = 1, nemax
      e = dble(ie)*de +emin
!----------------------------------<2023.04.13> "COULSCAT" is stopped due to bugs.
!--- s(1/2)-channel:
!<>      call coulscat(xind, 0, 1, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_s1(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC, s(1/2)".'
      endif
!--- p(3/2)-channel:
!<>      call coulscat(xind,1, 3, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_p3(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC, p(3/2)".'
      endif
!--- d(5/2)-channel:
!<>      call coulscat(xind, 2, 5, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_d5(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC, d(5/2)".'
      endif
!--- d(3/2)-channel:
!<>      call coulscat(xind, 2, 3, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_d3(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC".'
      endif
!--- f(7/2)-channel:
!<>      call coulscat(xind, 3, 7, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_f7(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC,".'
      endif
!--- f(5/2)-channel:
!<>      call coulscat(xind, 3, 5, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_f5(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC,".'
      endif
   1000  end do

   call deriv_5p(delta_s1, dedel_s1, de, nemax)
   call deriv_5p(delta_p3, dedel_p3, de, nemax)
   call deriv_5p(delta_p1, dedel_p1, de, nemax)
   call deriv_5p(delta_d5, dedel_d5, de, nemax)
   call deriv_5p(delta_d3, dedel_d3, de, nemax)
   call deriv_5p(delta_f7, dedel_f7, de, nemax)
   call deriv_5p(delta_f5, dedel_f5, de, nemax)

   il = 190000; IF(xind == 'neut') il = 290000
   do ie = 1, nemax
      e = de*dble(ie) +emin
      write (il+ 1, '(3(1X,G14.7))') e, delta_s1(ie), dedel_s1(ie)
      write (il+13, '(3(1X,G14.7))') e, delta_p3(ie), dedel_p3(ie)
!!      write (il+25, '(3(1X,G14.7))') e, delta_d5(ie), dedel_d5(ie)
      write (il+23, '(3(1X,G14.7))') e, delta_d3(ie), dedel_d3(ie)
!!      write (il+37, '(3(1X,G14.7))') e, delta_f7(ie), dedel_f7(ie)
      write (il+35, '(3(1X,G14.7))') e, delta_f5(ie), dedel_f5(ie)
   end do
   call rectime(19)
   write (6, *) 'The program CN_SC terminated normally (^_^)'
   end subroutine cn_sc

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
!<>      call CN_SC('prot')!---<2023.04.13> Stopped with minor bugs.
                    l = 9999 ; jj = -9999!---------------------------------filter
      call CN_TD('prot',l,jj)
   end if
   !----------------------------
   if(2.eq.2) then
   write(6,*) "**************************************************"
   write(6,*) "************* NEUTRON ****************************"
   write(6,*) "**************************************************"
!<>      call CN_SC('neut')!---<2023.04.13> Stopped with minor bugs.
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

