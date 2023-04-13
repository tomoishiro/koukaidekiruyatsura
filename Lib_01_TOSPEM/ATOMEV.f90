!///////////////////////////////////////////////////////////////////////////////
!<2023/03/22> T. Oishi
! SR ang_mag is replaced to dbmmag_*_RS, which are copied from CPNM1.f90 code.
!///////////////////////////////////////////////////////////////////////////////
module mdl_001_setting
   implicit none
   !---
   double precision, parameter :: pi = 3.141592653d0
   double precision, parameter :: hbarc = 197.329d0 !MeV.fm = keV.pm = eV.nm
   double precision, parameter :: alpha = 7.2973525698d-3 !fine structure constant
   double precision, parameter :: tol = 1.d-8
   complex (kind=kind(0d0)), parameter :: ai = (0.d0, 1.d0)
   !---
   double precision, parameter :: zc =  20.d0 ! proton-number of core
   double precision, parameter :: nc =  20.d0 !neutron-number of core
   double precision, parameter :: ac = zc + nc
   double precision, parameter :: anmass_p =  510998.995d0 !938.2720813d0 ! proton-mass
   double precision, parameter :: anmass_n =  939.5654133d0 !neutron-mass
!!   double precision, parameter :: BE_core = 7.976207d0*ac !MeV, BE for O-16 nucleus.
   double precision, parameter :: BE_core = 8.551305d0*ac !MeV, BE for Ca-40 nucleus.
!!   double precision, parameter :: BE_core = 8.253d0*ac !MeV, BE for Sn-100 nucleus.
   double precision, parameter :: mass_c = zc*anmass_p + nc*anmass_n - BE_core !MeV, for the core.
!<>   double precision, parameter :: xmu_cp = anmass_p*mass_c/(anmass_p + mass_c)
   double precision, parameter :: xmu_cp = anmass_p
   double precision, parameter :: xmu_cn = anmass_n*mass_c/(anmass_n + mass_c)
   !---
   double precision, parameter :: ecut = tol !MeV
   double precision, parameter :: RMAX = 4.d0 !120.d0
   double precision, parameter :: dr = 0.001d0
   integer, parameter :: nrmax = int(rmax/dr + 1.d-6) !max # of r-grid
   integer, parameter :: lmax = 5 !15 !max # of spatial ang-momenta
   integer, parameter :: jmax = 2*lmax+1
   integer, parameter :: ndmax = 40 !max # of nodes for s.p.sts
   integer, parameter :: nspmax = 300 !max # of s.p.sts
end module mdl_001_setting
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
module mdl_002_potentials
   use mdl_001_setting, only: ac; implicit none
   !--- parameters of the core-n (-p) WS(+Coulomb) potentials.
   double precision, parameter :: r00 = 1.25d0
   double precision, parameter :: a00 = 0.65d0
   double precision, parameter :: rc0 = r00*ac**(1.d0/3.d0)
   double precision, parameter :: ac0 = a00
   double precision, parameter :: W_0  = -55.57d0 !-53.2d0 ![MeV]
   double precision, parameter :: W_ls =  11.28d0 ! 22.1d0 ![MeV*fm^2]
contains
   !***********************************************************
   function v_cn(l, j, r) result (ff); implicit none
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: ff
      ff = vsp_ws_0(l, j, r) + vsp_ws_ls(l, j, r) + vsp_cent_n(l, r)
   end function v_cn
   !***********************************************************
   function v_cp(l, j, r) result (ff); implicit none
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: ff
!<>      ff = vsp_ws_0(l, j, r) + vsp_ws_ls(l, j, r) + vsp_cent_p(l, r) + vsp_coul(r)
      ff = vatom(r) + vsp_cent_p(l, r)
   end function v_cp
   !***********************************************************
   function vsp_ws_0(l, j, r) result (ff) !Woods-Saxon potential between the Core and a nucleon
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: rc, ac, ff
      rc = rc0 ; ac = ac0
      ff = W_0/(1.d0+exp((r-rc)/ac))
   end function
   !***********************************************************
   function vsp_ws_ls(l, j, r) result (ff) !ls-WS potential between the Core and a nucleon
      integer, intent (in) :: l, j
      double precision, intent (in) :: r
      double precision :: xls, yl, yj, df, rc, ac, ff
      rc = rc0 ; ac = ac0 ; yl = dble(l) ; yj = dble(j)/2.d0
      df = -exp((r-rc)/ac)/(1.d0+exp((r-rc)/ac))**2/ac
      xls = 0.5d0*(yj*(yj+1.d0)-yl*(yl+1.d0)-0.75d0)
      ff = W_ls*xls*df/r
   end function
   !***********************************************************
   function vsp_cent_p(l, r) result (ff) !--- centrifugal
      use mdl_001_setting, only: hbarc, xmu_cp; implicit none
      integer, intent (in) :: l
      double precision, intent (in) :: r
      double precision :: ff
      ff = dble(l*(l+1))*0.5d0*hbarc*hbarc/xmu_cp/r/r
   end function
   !***********************************************************
   function vsp_cent_n(l, r) result (ff) !--- centrifugal
      use mdl_001_setting, only: hbarc, xmu_cn; implicit none
      integer, intent (in) :: l
      double precision, intent (in) :: r
      double precision :: ff
      ff = dble(l*(l+1))*0.5d0*hbarc*hbarc/xmu_cn/r/r
   end function
   !***********************************************************
   function vsp_coul(r) result (ff) !--- Coulomb (only for a proton)
      use mdl_001_setting, only: hbarc, alpha, ac, zc; implicit none
      double precision, intent (in) :: r
      double precision :: vcoul, rb, ff
      rb = r00*ac**(1.d0/3.d0)
      if (r<rb) then
         vcoul = zc*hbarc*alpha*(3.d0-(r/rb)**2)*0.5d0/rb
      else
         vcoul = zc*hbarc*alpha/r
      end if
      ff = vcoul !for C-p
   end function
   !***********************************************************
   function WSD(r); implicit none
      double precision, intent (in) :: r
      double precision :: rc, ac, WSD
      rc = rc0 ; ac = ac0
      WSD = 1.d0/(1.d0+exp((r-rc)/ac))
   end function WSD
   !***********************************************************
   function vatom(r) result(ff)!---Hydrogen-atom potential for single electron.
      use mdl_001_setting, only: hbarc, alpha; implicit none
      double precision, intent (in) :: r
      double precision :: ff, zz, rb, vcoul
      zz = 2.d0; rb = 1.d-9
      if (r<rb) then
         vcoul = -zz*hbarc*alpha*(3.d0-(r/rb)**2)*0.5d0/rb
      else
         vcoul = -zz*hbarc*alpha/r
      end if
      ff = vcoul
   end function vatom
   !***********************************************************
end module mdl_002_potentials
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
module mdl_007_solvers
   implicit none
   integer, parameter :: noc_p = 0 !# of Proton- orbits occupied by the core.
   integer, parameter :: noc_n = 0 !# of Neutron-orbits occupied by the core.
                                   !=  3 for O-16. != 11 for Sn-100.
contains
   !*****************************************************************
   subroutine spsts_np(nsp0, ll0, jj0, node0, esp0, psi0, xind) ! single-particle wave functions.
      use mdl_001_setting, only : ecut,tol,rmax,nrmax,nspmax,ndmax,lmax
      use mdl_002_potentials, only : v_cp, v_cn, vsp_ws_0
      implicit none
      integer, intent (inout) :: nsp0, ll0(nspmax), jj0(nspmax), node0(nspmax)
      double precision, intent (inout) :: esp0(nspmax), psi0(nspmax, 0:nrmax)
      character(4), intent (in) :: xind
      integer :: j, k, l, m
      integer :: node00, nodemax, inode, ic, nxkmax, ir
      integer :: ncore, jjc(0:15), llc(0:15), nodec(0:15) !orbits occupied by the core
      double precision :: e, emax, emin, xkmax, psi00(0:nrmax)
      double precision :: rmaxd
      rmaxd = rmax
      IF((xind.ne.'prot') .and. (xind.ne.'neut')) THEN
        WRITE(6,*) "Check 'xind' in sp-states." ; STOP
      END IF
      !---core-states to be excluded(?).
      llc(0) =-1 ; jjc(0) =-1 ; nodec(0) =-1 !for free
      llc(1) = 0 ; jjc(1) = 1 ; nodec(1) = 0 !0s_1/2 (2)
      llc(2) = 1 ; jjc(2) = 3 ; nodec(2) = 0 !0p_3/2 (4)
      llc(3) = 1 ; jjc(3) = 1 ; nodec(3) = 0 !0p_1/2 (2) --- (8)
      llc(5) = 0 ; jjc(5) = 1 ; nodec(5) = 1 !1s_1/2 (2)
      llc(4) = 2 ; jjc(4) = 5 ; nodec(4) = 0 !0d_5/2 (6)
      llc(6) = 2 ; jjc(6) = 3 ; nodec(6) = 0 !0d_3/2 (4) --- (20)
      llc(7) = 3 ; jjc(7) = 7 ; nodec(7) = 0 !0f_7/2 (8) --- (28)
      llc(8) = 3 ; jjc(8) = 5 ; nodec(8) = 0 !0f_5/2 (6)
      llc(9) = 1 ; jjc(9) = 3 ; nodec(9) = 1 !1p_3/2 (4) --- (38)
      llc(10)= 1 ; jjc(10)= 1 ; nodec(10)= 1 !1p_1/2 (2) --- (40)
      llc(11)= 4 ; jjc(11)= 9 ; nodec(11)= 0 !0g_9/2 (10)--- (50)
      ncore = noc_p
      if(xind=='neut') ncore = noc_n !the number of the occupied states
      !---
      m = 0
      do 102 l = 0, lmax
!!         if((xind.ne.'prot') .and. (l.ge.2)) goto 102 !filter_1
!!         if(mod(l,2).eq.1) goto 102 !filter_1
!!         if(l.ne.2) goto 102 !filter_1
         do 101 j = 2*l - 1, 2*l + 1, 2
!!            if(j/=3) goto 101 !filter_2
            if (j<0) go to 101
            emax = ecut ; emin = -abs(vsp_ws_0(0,1,tol))
            call numerov_for(inode, l, j, ecut, psi00, xind) !determine nodemax from E_cut
            nodemax = inode - 1
            do 100 node00 = 0, nodemax
               if (node00>ndmax) go to 100
               do ic = 0, ncore !--- Pauli rule to remove SP states with the same quantum-number in the core
                  if (l==llc(ic) .and. j==jjc(ic) .and. node00==nodec(ic)) go to 100
               end do
               xkmax = log10(abs(emax-emin)/tol)/log10(2.d0) + 0.5d0!--- formula : log_(2)(x) = log_(10)(x) / log_(10)(2)
               nxkmax = int(xkmax)
               do k = 1, nxkmax!--- iteration to approximate `E' by both-side-attack
                  e = (emin + emax)*0.5d0
                  call numerov_for(inode,l,j,e,psi00,xind)
                  if (inode>node00) then
                     emax = e
                  else
                     emin = e
                  end if
               end do
               m = m + 1
               if (m>nspmax) then
                  write (6, *) 'Increase ``nspmax'' (~_~)' ; stop
               end if
               !--- matching with U(r) from large distance for bound state.
               if((e-vsp_ws_0(l,j,rmaxd))<0.d0) call numerov_mat(l,j,e,psi00,xind)
               esp0(m) = e ; ll0(m) = l ; jj0(m) = j ; node0(m) = node00
               do ir = 0, nrmax
                  psi0(m, ir) = psi00(ir)
               end do
               emax = ecut ; emin = e
100         end do
101      end do
102   end do
      nsp0 = m
      call sort_np(nsp0, ll0, jj0, node0, esp0, psi0)
      return
   end subroutine spsts_np
!**********************************************************************
   subroutine numerov_for(inode, l, j, e, psi00, xind)! Subroutine for integration of the Schroedinger eq. by Numerov method
      use mdl_001_setting
      use mdl_002_potentials, only : v_cn, v_cp
      implicit none
      integer, intent (inout) :: inode
      integer, intent (in) :: l, j
      integer :: nrmaxd
      double precision, intent (in) :: e
      double precision, intent (inout) :: psi00(0:nrmax)
      character(4), intent(in) :: xind
      integer :: ir
      double precision :: r, r0, r1, fac, psi0, psi1, dd, cc, cin, rmaxd
      psi00 = 0.d0
      rmaxd = rmax ; nrmaxd = nrmax
      inode = 0
      fac = dr*dr*(2.d0*xmu_cn/hbarc/hbarc)
      if(xind=='prot') then
          fac = dr*dr*(2.d0*xmu_cp/hbarc/hbarc)
          goto 200
      end if
!--- for core-n
 100  psi0 = dr**(l+1)
      psi1 = (2.d0+fac*(v_cn(l,j,dr)-e)*5.d0/6.d0)*psi0
      psi1 = psi1/(1.d0-fac*(v_cn(l,j,dr+dr)-e)/12.d0)
      psi00(0) = 0.d0
      psi00(1) = psi0
      psi00(2) = psi1
      do ir = 3, nrmaxd
         r = dr*dble(ir)
         r0 = dr*dble(ir-2)
         r1 = dr*dble(ir-1)
         !---
         dd = psi0 - 2.d0*psi1
         dd = dd - fac*(v_cn(l,j,r0)-e)*psi0/12.d0
         dd = dd - fac*(v_cn(l,j,r1)-e)*psi1*5.d0/6.d0
         !---
         cc = -fac*(v_cn(l,j,r)-e)/12.d0 + 1.d0
         cin = 1.d0/cc
         psi00(ir) = -cin*dd
         if (ir>5 .and. psi00(ir)*psi00(ir-1)<0.d0) inode = inode + 1
         psi0 = psi1
         psi1 = psi00(ir)
      end do
      goto 300
!--- for core-p
 200  psi0 = dr**(l+1)
      psi1 = (2.d0+fac*(v_cp(l,j,dr)-e)*5.d0/6.d0)*psi0
      psi1 = psi1/(1.d0-fac*(v_cp(l,j,dr+dr)-e)/12.d0)
      psi00(0) = 0.d0
      psi00(1) = psi0
      psi00(2) = psi1
      do ir = 3, nrmaxd
         r = dr*dble(ir)
         r0 = dr*dble(ir-2)
         r1 = dr*dble(ir-1)
         dd = psi0 - 2.d0*psi1
         dd = dd - fac*(v_cp(l,j,r0)-e)*psi0/12.d0
         dd = dd - fac*(v_cp(l,j,r1)-e)*psi1*5.d0/6.d0
         cc = -fac*(v_cp(l,j,r)-e)/12.d0 + 1.d0
         cin = 1.d0/cc
         psi00(ir) = -cin*dd
         if (ir>5 .and. psi00(ir)*psi00(ir-1)<0.d0) inode = inode + 1
         psi0 = psi1
         psi1 = psi00(ir)
      end do
!--- normalization
 300  fac = 0.d0
      do ir = 0, nrmax
         fac = fac + psi00(ir)*psi00(ir)*dr
      end do
      psi00 = psi00/sqrt(fac)
      return
   end subroutine numerov_for
!**************************************************************
   subroutine numerov_mat(l, j, e, psi00, xind)! Solve the Schroedinger eq. backward in order to ensure the asymptotic form for E < 0.
      use mdl_001_setting, only : rmax,nrmax,dr,ecut,hbarc,xmu_cp,xmu_cn,tol,lmax,ndmax,ac
      use mdl_002_potentials
      implicit none
      integer, intent (in) :: l, j
      double precision, intent (in) :: e
      double precision, intent (inout) :: psi00(0:nrmax)
      character(4), intent(in) :: xind
      integer :: ir, irmatch, nrmaxd
      double precision :: r, r0, r1, rmatch, rmaxd
      double precision :: evl, fac, ak, psi0, psi1, dd, cc, cin, cmatch, psid(0:nrmax)
!!      psi00 = 0.d0
      rmaxd = rmax ; nrmaxd = nrmax
      rmatch = 1.5d0*ac**(1.d0/3.d0) !!0.5d0*rmaxd
      irmatch = int(rmatch/dr+tol)
      fac = dr*dr*(2.d0*xmu_cn/hbarc/hbarc)
      if(xind=='prot') then
          fac = dr*dr*(2.d0*xmu_cp/hbarc/hbarc)
          goto 200
      end if
!--- For core-neutron
 100  evl = v_cn(l, j, rmaxd) - e
      ak = sqrt(abs(evl)*2.d0*xmu_cn/hbarc/hbarc)
      psi0 = exp(-ak*dble(nrmaxd)*dr)
      psi1 = exp(-ak*dble(nrmaxd-1)*dr)
      !---
      psid(0) = 0.d0
      psid(nrmaxd) = psi0
      psid(nrmaxd-1) = psi1
      !---
      do ir = nrmaxd - 2, irmatch, -1
         r = dr*dble(ir)
         r0 = dr*dble(ir+2)
         r1 = dr*dble(ir+1)
         dd = psi0 - 2.d0*psi1
         dd = dd - fac*(v_cn(l,j,r0)-e)*psi0/12.d0
         dd = dd - fac*(v_cn(l,j,r1)-e)*psi1*5.d0/6.d0
         cc = -fac*(v_cn(l,j,r)-e)/12.d0 + 1.d0
         cin = 1.d0/cc
         psid(ir) = -cin*dd
         psi0 = psi1
         psi1 = psid(ir)
      end do
      !--- matching
      cmatch = psi00(irmatch)/psid(irmatch)
      do ir = irmatch, nrmaxd
         psi00(ir) = psid(ir)*cmatch
      end do
      goto 300
!--- For core-proton
 200  evl = v_cp(l, j, rmaxd) - e
      ak = sqrt(abs(evl)*2.d0*xmu_cp/hbarc/hbarc)
      psi0 = exp(-ak*dble(nrmaxd)*dr)
      psi1 = exp(-ak*dble(nrmaxd-1)*dr)
      !---
      psid(0) = 0.d0
      psid(nrmaxd) = psi0
      psid(nrmaxd-1) = psi1
      !---
      do ir = nrmaxd - 2, irmatch, -1
         r = dr*dble(ir)
         r0 = dr*dble(ir+2)
         r1 = dr*dble(ir+1)
         dd = psi0 - 2.d0*psi1
         dd = dd - fac*(v_cp(l,j,r0)-e)*psi0/12.d0
         dd = dd - fac*(v_cp(l,j,r1)-e)*psi1*5.d0/6.d0
         cc = -fac*(v_cp(l,j,r)-e)/12.d0 + 1.d0
         cin = 1.d0/cc
         psid(ir) = -cin*dd
         psi0 = psi1
         psi1 = psid(ir)
      end do
      !--- matching
      cmatch = psi00(irmatch)/psid(irmatch)
      do ir = irmatch, nrmaxd
         psi00(ir) = psid(ir)*cmatch
      end do
!--- normalization
 300  fac = 0.d0
      do ir = 0, nrmax
         fac = fac + psi00(ir)*psi00(ir)*dr
      end do
      psi00 = psi00/sqrt(fac)
      return
   end subroutine numerov_mat
!**************************************************************
   subroutine sort_np(nsp0, ll0, jj0, node0, esp0, psi0)! A routine to sort the eigenvalues obtained in "spbsis" and the associated eigen functions.
      use mdl_001_setting
      implicit none
      integer, intent (inout) :: nsp0, ll0(nspmax), jj0(nspmax), node0(nspmax)
      double precision, intent (inout) :: esp0(nspmax), psi0(nspmax, 0:nrmax)
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
         !---
         if (kk/=i) then
            esp0(kk) = esp0(i)
            esp0(i) = p
            !---
            ip = jj0(i)
            jj0(i) = jj0(kk)
            jj0(kk) = ip
            !---
            ip = ll0(i)
            ll0(i) = ll0(kk)
            ll0(kk) = ip
            !---
            ip = node0(i)
            node0(i) = node0(kk)
            node0(kk) = ip
            !---
            do ir = 0, nrmax
               p = psi0(i, ir)
               psi0(i, ir) = psi0(kk, ir)
               psi0(kk, ir) = p
            end do
         end if
      end do
      return
   end subroutine sort_np
!*****************************************************************
   subroutine spbasis_p(xind, nsp, ll, jj, node, esp, psi)
      use mdl_001_setting, only: nspmax, nrmax
      implicit none
      character (4), intent (in) :: xind
      integer, intent (inout) :: nsp, ll(nspmax), jj(nspmax), node(nspmax)
      double precision, intent (inout) :: esp(nspmax), psi(nspmax, 0:nrmax)
      integer, save :: nspd, lld(nspmax), jjd(nspmax), noded(nspmax)
      double precision, save :: espd(nspmax), psid(nspmax, 0:nrmax)
      if (xind=='save') then
         nspd = nsp
         lld = ll
         jjd = jj
         noded = node
         espd = esp
         psid = psi
      else if (xind=='load') then
         nsp = nspd
         ll = lld
         jj = jjd
         node = noded
         esp = espd
         psi = psid
      end if
   end subroutine spbasis_p
!*****************************************************************
   subroutine spbasis_n(xind, nsp, ll, jj, node, esp, psi)
      use mdl_001_setting, only: nspmax, nrmax
      implicit none
      character (4), intent (in) :: xind
      integer, intent (inout) :: nsp, ll(nspmax), jj(nspmax), node(nspmax)
      double precision, intent (inout) :: esp(nspmax), psi(nspmax, 0:nrmax)
      integer, save :: nspd, lld(nspmax), jjd(nspmax), noded(nspmax)
      double precision, save :: espd(nspmax), psid(nspmax, 0:nrmax)
      if (xind=='save') then
         nspd = nsp
         lld = ll
         jjd = jj
         noded = node
         espd = esp
         psid = psi
      else if (xind=='load') then
         nsp = nspd
         ll = lld
         jj = jjd
         node = noded
         esp = espd
         psi = psid
      end if
   end subroutine spbasis_n
!*****************************************************************
end module mdl_007_solvers
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
module mdl_008_elem
contains
!*************************************************************
   function c3j(j1, m1, j2, m2, j3, m3) result (ff)
   != (J1/2 J2/2 J3/2)
   !  (M1/2 M2/2 M3/2), the Racah formula for the 3j-symbol
      implicit none
      integer, intent (in) :: j1, m1, j2, m2, j3, m3
      integer :: iphase
      double precision :: ff
      iphase = int(0.5d0*dble(j1-j2-m3) + 0.000001d0)
      ff = phass(iphase)*cg(j1, m1, j2, m2, j3, -m3)/sqrt(dble(j3+1))
      return
   end function c3j
!**************************************************************
   function ang_ele(lf,jf,li,ji, ilam) result (ff)
      integer, intent(INOUT) :: lf,jf,li,ji, ilam
      integer :: iphase
      double precision :: ff, sr
      sr = 0.5d0*dble(1 + (-1)**(li+lf+ilam))
      iphase = int(0.5d0*dble(jf-1) + 0.000001d0)
      ff = c3j(jf, -1, ilam*2,0, ji, 1)*phass(iphase)*sr
   end function
!****************************************************************
   !Angular part of double-bar M.E. of magnetic transition,
   !<lf, jf/2 || M(lam) || li,ji>, given in Eq.(B.82) in R&S textbook.
   !This result is multiplied by sqrt(4\pi /(2*lam+1)).
   function dbmmag_l_RS(jf,lf,  lam,   ji,li) result (ff)
   implicit none
   integer, intent(IN) :: jf, lf, ji, li, lam
   double precision :: ff, gl
   integer :: k, njf, nji, mjf
   ff = 0.d0
   if(mod(lf+li+lam,2) .eq. 0) go to 500 !If "lf+li+lam"=even, it's zero.
   !---
   njf = (jf+1)/2 !="true j_f" + 1/2.
   mjf = (jf-1)/2 !="true j_f" - 1/2.
   nji = (ji+1)/2
   !mji = (ji-1)/2
   k = nji*(-1)**(li+nji) + njf*(-1)**(lf+njf)
   !---
   gl = 1.d0!dummy
   ff = -gl*(1.d0 + dble(k)/dble(lam+1))
   ff = ff * sqrt(dble((ji+1)*(jf+1))) * c3j(jf,-1, 2*lam,0, ji,1) * dble(lam-k)
   !---
   500 continue
   ff = ff * (-1.d0)**(mjf)
   end function dbmmag_l_RS
   !+++++++++++++++++++++++
   function dbmmag_s_RS(jf,lf,  lam,   ji,li) result (ff)
   implicit none
   integer, intent(IN) :: jf, lf, ji, li, lam
   double precision :: ff, gs
   integer :: k, njf, nji, mjf
   ff = 0.d0
   if(mod(lf+li+lam,2) .eq. 0) go to 500 !If "lf+li+lam"=even, it's zero.
   !---
   njf = (jf+1)/2 !="true j_f" + 1/2.
   mjf = (jf-1)/2 !="true j_f" - 1/2.
   nji = (ji+1)/2
   !mji = (ji-1)/2
   k = nji*(-1)**(li+nji) + njf*(-1)**(lf+njf)
   !---
   gs = 1.d0!dummy
   ff = 0.5d0*gs * sqrt(dble((ji+1)*(jf+1))) * c3j(jf,-1, 2*lam,0, ji,1) * dble(lam-k)
   !---
   500 continue
   ff = ff * (-1.d0)**(mjf)
   end function dbmmag_s_RS
!****************************************************************
!---<2023.03.22> This function does not coincide with dbmmag_*_RS, and thus,
!--- eliminated from the calculations.
   function ang_mag(lf,jf,li,ji, ilam) result (ff)
     integer, intent(INOUT) :: lf,jf,li,ji, ilam
     integer :: iphase
     double precision :: ff
                   sr = 0.5d0*dble(1 - (-1)**(li+lf+ilam))
     iphase = int(0.5d0*dble(jf-1) + 0.000001d0)
     ff = c3j(jf, -1, ilam*2,0, ji, 1)*phass(iphase)*sr
   end function ang_mag
!**************************************************************
   pure function dlt(i, j) result (ff)! Kronecker's delta_ij
      implicit none
      integer, intent (in) :: i, j
      double precision :: ff, sr
      ff = 0.d0
      if (i==j) ff = 1.d0
      return
   end function dlt
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
!****************************************************************
   function cg(j1, m1, j2, m2, j3, m3) result (ff)
   != Clebsh-Gordan Coefficient,
   !  C( j1/2 m1/2 ; j2/2 m2/2 | j3/2 m3/2 )
      implicit none
      integer, intent (in) :: j1, m1, j2, m2, j3, m3
      integer :: n, ka, kb, kc, kd, k1, k2, k3, k4, k5, k6, ka1, ka2, ka3, ka4, ka5
      double precision :: aj1, aj2, aj3, am1, am2, am3, an, del, ff
      ff = 0.d0
      if (m1+m2/=m3) go to 500
      if ((min(j1,j2,j3)<0) .or. (j3<abs(j1-j2)) .or. (j3>j1+j2)) go to 500
!      if ((abs(m1)>j1) .or. (abs(m2)>j2) .or. (abs(m3)>j3)) go to 500
      aj1 = dble(j1)*0.5d0
      aj2 = dble(j2)*0.5d0
      aj3 = dble(j3)*0.5d0
      am1 = dble(m1)*0.5d0
      am2 = dble(m2)*0.5d0
      am3 = dble(m3)*0.5d0
      ka = nint(aj1+aj2-aj3)
      kb = nint(aj3+aj1-aj2)
      kc = nint(aj2+aj3-aj1)
      kd = nint(aj1+aj2+aj3+1)
      k1 = nint(aj1+am1)
      k2 = nint(aj1-am1)
      k3 = nint(aj2+am2)
      k4 = nint(aj2-am2)
      k5 = nint(aj3+am3)
      k6 = nint(aj3-am3)
      del = ( factlog(ka) + factlog(kb) + factlog(kc) - factlog(kd) + factlog(k1) + factlog(k2) &
            + factlog(k3) + factlog(k4) + factlog(k5) + factlog(k6) )*0.5d0
      k1 = nint(aj1+aj2-aj3)
      k2 = nint(aj1-am1)
      k3 = nint(aj2+am2)
      do n = 0, max(k1, k2, k3)
         ka1 = nint(aj1+aj2-aj3-n)
         if (ka1<0.d0) go to 100
         ka2 = nint(aj3-aj2+am1+n)
         if (ka2<0.d0) go to 100
         ka3 = nint(aj3-aj1-am2+n)
         if (ka3<0.d0) go to 100
         ka4 = nint(aj1-am1-n)
         if (ka4<0.d0) go to 100
         ka5 = nint(aj2+am2-n)
         if (ka5<0.d0) go to 100
         an = n*1.d0
         ff = ff + (-1.d0)**n*exp(del-factlog(n)-factlog(ka1)-factlog(ka2)-factlog(ka3)-factlog(ka4)-factlog(ka5))
100   end do
      ff = ff*sqrt(2.d0*aj3+1.d0)
500   return
   end function cg
!****************************************************************
   pure function phass(n) result (ff) !=(-)**n
      implicit none
      integer, intent (in) :: n
      integer :: i
      double precision :: ff
      ff = -1.d0
      do i = 0, abs(n)
         ff = -1.d0*ff
      enddo
      return
   end function phass
!**************************************************************
  function ang_y2s(lf,jf,li,ji, ilam) result (ff)
    integer, intent(INOUT) :: lf,jf,li,ji, ilam
    integer :: iphase
    double precision :: ff
    ff = c9j(li*2,lf*2,ilam*2,  1,1,2,  ji,jf,2)
  end function
!*************************************************************
   function c9j(j11, j12, j13, j21, j22, j23, j31, j32, j33)
!  The Racah formula for the 9j-symbol
!  {J11/2  J12/2  J13/2}
!  {J21/2  J22/2  J23/2}
!  {J31/2  J32/2  J33/2}
      implicit none
      integer, intent (in) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
      integer :: kmin1, kmin2, kmin3, kmax1, kmax2, kmax3, kmin, kmax, k1, k
      double precision :: s1, s2, s3, p, c9j
      kmin1 = abs(j11-j33)
      kmin2 = abs(j32-j21)
      kmin3 = abs(j23-j12)
      kmax1 = j11 + j33
      kmax2 = j32 + j21
      kmax3 = j23 + j12
      kmin = max(kmin1, kmin2, kmin3) + 1
      kmax = min(kmax1, kmax2, kmax3) + 1
      c9j = 0.d0
      if (kmin>kmax) go to 100
      do k1 = kmin, kmax, 2
         k = k1 - 1
         s1 = c6j_v2(j11, j21, j31, j32, j33, k)
         s2 = c6j_v2(j12, j22, j32, j21, k, j23)
         s3 = c6j_v2(j13, j23, j33, k, j11, j12)
         p = (k*1.d0+1.d0)*(-1.d0)**k
         c9j = c9j + p*s1*s2*s3
      end do
100   return
   end function c9j
!****************************************************************
   function c6j_v2(j1, j2, j3, j4, jf, j6) result (ff)
!  The Racah formula for the 6j-symbol ; {j1/2  j2/2  j3/2}
!                                        {j4/2  jf/2  j6/2}
! = c6j(j1, j2, jf, j4, j3, j6)
      implicit none
      integer, intent (in) :: j1, j2, j3, j4, jf, j6
      integer :: numin, numax, nu, jy1, jy2, jy3, jy4, jy5, jy6, jy7
      double precision :: ff, fctot, aj1, aj2, aj3, aj4, aj5, aj6, z1, jt1, jt2, jt3, jt4, jz1, jz2, jz3
      ff = 0.0
      aj1 = dble(j1)*0.5d0
      aj2 = dble(j2)*0.5d0
      aj3 = dble(j3)*0.5d0
      aj4 = dble(j4)*0.5d0
      aj5 = dble(jf)*0.5d0
      aj6 = dble(j6)*0.5d0
      if (min(aj1,aj2,aj3,aj4,aj5,aj6)<0.d0) return
      z1 = delr(aj1, aj2, aj3)
      if (z1==0.0) go to 100
      z1 = delr(aj1, aj5, aj6)*z1
      if (z1==0.0) go to 100
      z1 = delr(aj4, aj2, aj6)*z1
      if (z1==0.0) go to 100
      z1 = delr(aj4, aj5, aj3)*z1
      if (z1==0.0) go to 100
      z1 = sqrt(z1)
      jt1 = aj1 + aj2 + aj3 + 1.d-6
      jt2 = aj1 + aj5 + aj6 + 1.d-6
      jt3 = aj4 + aj2 + aj6 + 1.d-6
      jt4 = aj4 + aj5 + aj3 + 1.d-6
      jz1 = aj1 + aj2 + aj4 + aj5 + 1.d-6
      jz2 = aj2 + aj3 + aj5 + aj6 + 1.d-6
      jz3 = aj3 + aj1 + aj6 + aj4 + 1.d-6
      numin = nint(max(jt1,jt2,jt3,jt4))
      numax = nint(min(jz1,jz2,jz3))
      if (numax<numin) go to 100
      do nu = numin, numax
         jy1 = nu - int(jt1)
         jy2 = nu - int(jt2)
         jy3 = nu - int(jt3)
         jy4 = nu - int(jt4)
         jy5 = int(jz1) - nu
         jy6 = int(jz2) - nu
         jy7 = int(jz3) - nu
         fctot = (-1.d0)**nu*fact(nu+1)/fact(jy1)/fact(jy2)/fact(jy3)/fact(jy4)/fact(jy5)/fact(jy6)/fact(jy7)
         ff = ff + fctot
      end do
      ff = ff*z1
100   return
   end function c6j_v2
   !--------------------------------------
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
   !--------------------------------------
   function delr(aj1, aj2, aj3)
      implicit none
      double precision, intent (in) :: aj1, aj2, aj3
      double precision :: delr, test
      integer :: jz1, jz2, jz3, jz4
      jz1 = int(aj1+aj2-aj3+1.d-6)
      if (jz1<0) go to 100
      jz2 = int(aj2+aj3-aj1+1.d-6)
      if (jz2<0) go to 100
      jz3 = int(aj3+aj1-aj2+1.d-6)
      if (jz3<0) go to 100
      jz4 = int(aj1+aj2+aj3+1+1.d-6)
      test = aj1 + aj2 + aj3 + 1.d0 - dble(jz4) + 1.d-6
      if (abs(test)>1.d-4) go to 100
      delr = fact(jz1)*fact(jz2)*fact(jz3)/fact(jz4)
      return
100   delr = 0.d0
      return
   end function delr
   !--------------------------------------
!****************************************************************
end module mdl_008_elem
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
module mdl_099_summary
   private :: rectime
contains
   !***********************************************************************
   subroutine SPEM
   use mdl_001_setting
   use mdl_002_potentials, only : v_cp, v_cn
   use mdl_007_solvers, only : spsts_np, spbasis_p, spbasis_n
   use mdl_008_elem, only : ang_ele, phass, dbmmag_l_RS, dbmmag_s_RS
   implicit none
   !---
   integer :: nsp1, ll1(nspmax), jj1(nspmax), node1(nspmax) !for `spbasis_p'
   double precision :: esp_p(nspmax), psi1(nspmax, 0:nrmax) !for `spbasis_p'
   !---
   integer :: nsp2, ll2(nspmax), jj2(nspmax), node2(nspmax) !for `spbasis_n'
   double precision :: esp_n(nspmax), psi2(nspmax, 0:nrmax) !for `spbasis_n'
   !---
   integer :: i,j1,j2,l1,l2,m,n1,n2,ir
   integer :: k1,k2,k3,k4!--- Labels of initial and final states of interest
   double precision :: x,y,p,r, qp,qn, vp,vn
   integer :: ni,li,ji,nf,lf,jf,nem,ilam
   double precision :: gl_p, gl_n, gs_p, gs_n, x1
   namelist/mine/ni,li,ji,nf,lf,jf,nem,ilam
   namelist/gfactors/gl_p, gl_n, gs_p, gs_n
   open( 711, file='PARAM.inp', status='old')
   read( 711, nml=mine)
   read( 711, nml=gfactors)
   close(711)
   open( 711, file='RESULT_RUN1.dat', status='replace')
   write(711, nml=mine)
   write(711, nml=gfactors)
   call rectime(0)
   !---Potentials
   open(8491, file='Yed_0013_SP_Potentials.dat', action='write')
   open(8492, file='Yed_0011_SP_Potentials.dat', action='write')
   open(8493, file='Yed_0001_SP_Potentials.dat', action='write')
   open(8495, file='Yed_0025_SP_Potentials.dat', action='write')
   open(8496, file='Yed_0023_SP_Potentials.dat', action='write')
      do ir = 0, nrmax
         x1 = dble(ir)*dr + tol
         write(8491,'(5(2X,F14.7))') x1, v_cp(1,3, x1), v_cp(1,3, x1), &
                                         v_cn(1,3, x1), v_cn(1,3, x1)
         write(8492,'(5(2X,F14.7))') x1, v_cp(1,1, x1), v_cp(1,1, x1), &
                                         v_cn(1,1, x1), v_cn(1,1, x1)
         write(8493,'(5(2X,F14.7))') x1, v_cp(0,1, x1), v_cp(0,1, x1), &
                                         v_cn(0,1, x1), v_cn(0,1, x1)
         write(8495,'(5(2X,F14.7))') x1, v_cp(2,5, x1), v_cp(2,5, x1), &
                                         v_cn(2,5, x1), v_cn(2,5, x1)
         write(8496,'(5(2X,F14.7))') x1, v_cp(2,3, x1), v_cp(2,3, x1), &
                                         v_cn(2,3, x1), v_cn(2,3, x1)
      end do
   close(8491)
   close(8492)
   close(8493)
   close(8495)
   close(8496)
   !--- s.p. states for core-proton
   call spsts_np(nsp1, ll1, jj1, node1, esp_p, psi1, 'prot')
   call spbasis_p('save', nsp1, ll1, jj1, node1, esp_p, psi1)
   write (711, *) '(proton-core) s.p. states'
   write (711, *) '   # Node    L    J*2      E(MeV)'
   do 11 i = 1, min(20,nsp1)
      n1 = node1(i) ; l1 = ll1(i) ; j1 = jj1(i) ; x = esp_p(i)
      write (711, '(4(2X,I3),F15.7)') i, n1, l1, j1, x
11 end do
   write (711, '(4(2X,I3),F15.7)') nsp1, node1(nsp1), ll1(nsp1), jj1(nsp1), esp_p(nsp1)
   write (711, *) 'number of s.p.basis=', nsp1
   write (711, *) ''

   !--- s.p. states for core-neutron
   call spsts_np(nsp2, ll2, jj2, node2, esp_n, psi2, 'neut')
   call spbasis_n('save', nsp2, ll2, jj2, node2, esp_n, psi2)
   write (711, *) '(neutron-core) s.p. states'
   write (711, *) '   # Node    L    J*2      E(MeV)'
   do 12 i = 1, min(20,nsp2)
      n2 = node2(i) ; l2 = ll2(i) ; j2 = jj2(i) ; x = esp_n(i)
      write (711, '(4(2X,I3),F15.7)') i, n2, l2, j2, x
12 end do
   write (711, '(4(2X,I3),F15.7)') nsp2, node2(nsp2), ll2(nsp2), jj2(nsp2), esp_n(nsp2)
   write (711, *) 'number of s.p.basis=', nsp2
   write (711, *) ''

   !--- find the initial & final states of interest
   k1 = -10; k2 = -10
   do i = 1, nsp1
      if (node1(i).eq.ni .and. ll1(i).eq.li .and. jj1(i).eq.ji) k1 = i!--- proton's initial
      if (node1(i).eq.nf .and. ll1(i).eq.lf .and. jj1(i).eq.jf) k2 = i!--- proton's final
   end do
   k3 = -10; k4 = -10
   do i = 1, nsp2
      if (node2(i).eq.ni .and. ll2(i).eq.li .and. jj2(i).eq.ji) k3 = i!--- neutron's initial
      if (node2(i).eq.nf .and. ll2(i).eq.lf .and. jj2(i).eq.jf) k4 = i!--- neutron's final
   end do

   !--- angular part of <f || Q || i>
   qp=1.d0 ; qn=1.d0 ; p=0.d0
   if (nem.eq.8) then
      p = sqrt( dble((2*ilam+1)*(ji+1)*(jf+1))*0.25d0/pi  )
      qp = p*ang_ele(lf,jf, li,ji, ilam)
      qn = qp!--- where the charge of neutron is assumed as e, too.
   else if (nem.eq.9) then
      p = sqrt( dble(2*ilam+1)*0.25d0/pi  )
      qp = p*( gl_p*dbmmag_l_RS(jf,lf, ilam, ji,li) +gs_p*dbmmag_s_RS(jf,lf, ilam, ji,li) )
      qn = p*( gl_n*dbmmag_l_RS(jf,lf, ilam, ji,li) +gs_n*dbmmag_s_RS(jf,lf, ilam, ji,li) )
!<>      i = li + (ji+1)/2
!<>      r =    0.5d0*dble(ji+1)*phass(i)
!<>      i = lf + (jf+1)/2
!<>      r = r +0.5d0*dble(jf+1)*phass(i)
!<>      m = int(r + 1.d-8)
!<>      qp = dble(ilam-m)*(  0.5d0*gs_p - gl_p*(1.d0 +dble(m)/dble(ilam+1))  )
!<>      qp = qp*ang_mag(lf,jf, li,ji, ilam)*p
!<>      qn = dble(ilam-m)*(  0.5d0*gs_n - gl_n*(1.d0 +dble(m)/dble(ilam+1))  )
!<>      qn = qn*ang_mag(lf,jf, li,ji, ilam)*p
      !!qp = qp*ang_y2s(lf,jf, li,ji, ilam)*p
      !!qn = qn*ang_y2s(lf,jf, li,ji, ilam)*p
   end if

   !--- radial integration
   open( 714, file='SPEMWF_PRO.dat', status='replace')
   open( 715, file='SPEMWF_NEU.dat', status='replace')
   vp = 0.d0
   vn = 0.d0
   if(nem.eq.8) m = ilam!--- for electric modes
   if(nem.eq.9) m = ilam-1!--- for magnetic modes
   do ir = 0, nrmax
      r = dr*dble(ir)
      vp = vp +psi1(k2,ir)*psi1(k1,ir)*(r**m +tol)!proton
      vn = vn +psi2(k4,ir)*psi2(k3,ir)*(r**m +tol)!neutron
      write (714, *) r, psi1(k2,ir), psi1(k1,ir), psi1(k2,ir)*psi1(k1,ir)*(r**m +tol)
      write (715, *) r, psi2(k4,ir), psi2(k3,ir), psi2(k4,ir)*psi2(k3,ir)*(r**m +tol)
   end do
   vp = vp*dr
   vn = vn*dr
   close(714)
   close(715)

   !--- results
   call units
   write (711, *) 'Mode of Q(J):',nem,'for (8:Elec, 9:Mag) with J =', ilam
   write (711, *) ""
   write (711, *) '|i>   with (n,l,j*2)=', node1(k1),ll1(k1),jj1(k1)
   write (711, *) '  <f| with (n,l,j*2)=', node1(k2),ll1(k2),jj1(k2), ' (proton)'
   write (711, *) '|i>   with (n,l,j*2)=   ', node2(k3),ll2(k3),jj2(k3)
   write (711, *) '  <f| with (n,l,j*2)=   ', node2(k4),ll2(k4),jj2(k4), ' (neutron)'
   write (711, *) ""
   write (711, *) "<><<><><><><<><><><><<><><><><<><><><><<><><><><<><><><><<><><>"
   write (711, *) " [A] Suhonen's formulas Eq. (6.23) and (6.24), where"
   write (711, *) "     both the radial and angular integrations are properly computed:"
   501 format("-->                       B_{QJ} =|<f|Q(J)|i>|**2 /(2j_i +1) = ",X,F13.6,X,"(proton)")
   502 format("-->                                                            ",X,F13.6,X,"(neutron)")
   x = vp*qp
   y = vn*qn
   write (711, 501)   (x**2)/dble(ji+1)
   write (711, 502)   (y**2)/dble(ji+1)
   write (711, *) " ---"
   write (711, *) " radial integration of <f|Q|i> =", vp, " (prot.)"
   write (711, *) "                               =", vn, " (neut.)"
   write (711, *) " ---"
   write (711, *) " angular part of <f|Q_J|i> =", qp, " (prot.)"
   write (711, *) "                         =", qn, " (neut.)"
      p = sqrt( dble(2*ilam+1)*0.25d0/pi  )
   write (711, *) " ---"
   write (711, *) " angular part of <f|Q_J|i>, but the factor,"
   write (711, *) " sqrt((2J+1)/4pi) is eliminated: =", qp/p, " (prot.)"
   write (711, *) "                                 =", qn/p, " (neut.)"

   765  continue
   close (711)
   end subroutine SPEM
   !***********************************************************************
   subroutine Weisskopf
   use mdl_001_setting
   use mdl_008_elem, only : ang_ele, dbmmag_l_RS, dbmmag_s_RS
   implicit none
   integer :: i,m
   double precision :: vp,vn,qp,qn,v,w,sr,p,r,ff1
   integer :: ni,li,ji,nf,lf,jf,nem,ilam
   double precision :: gl_p, gl_n, gs_p, gs_n
   namelist/mine/ni,li,ji,nf,lf,jf,nem,ilam
   namelist/gfactors/gl_p, gl_n, gs_p, gs_n
   open( 711, file='PARAM.inp', status='old')
   read( 711, nml=mine)
   read( 711, nml=gfactors)
   close(711)
   open( 711, file='RESULT_RUN2.dat', status='replace')
   !--- results-1 by Ring & Schuck
   write (711, *) ""
   write (711, *) "<><<><><><><<><><><><<><><><><<><><><><<><><><><<><><><><<><><>"
   write (711, *) " [B] VS Weisskopf estimate (1) [Phys. Rev. (1951) 83, 1073]:"
   r = 1.2d0*ac**(1.d0/3.d0)
   if (nem.eq.8) then!--- electric modes.
                             sr = 0.5d0*dble(1 + (-1)**(li+lf+ilam))
      v = r**(ilam+ilam)
      v = v*(3.d0/(3.d0 +dble(ilam)))**2
      v = v*0.25d0/pi
      w = v!--- where the charge of neutron is assumed as e, too.
   else if (nem.eq.9) then!--- magnetic modes.
                             sr = 0.5d0*dble(1 - (-1)**(li+lf+ilam))
      v = r**(ilam+ilam-2)
      v = v*(3.d0/(3.d0 +dble(ilam)))**2
      w = v*(0.d0 + 1.913d0*1.913d0)/pi!magnetic of N
      v = v*(1.d0 + 2.793d0*2.793d0)/pi!magnetic of P
!<>      v = v*(10.d0)/pi!magnetic of P
   end if
   write (711, *) ' Mode of Q(J):',nem,'for (8:Elec, 9:Mag) with J =', ilam
   501 format("-->                           B_{QJ} =|<f|Q(J)|i>|**2 /(2j_i +1) = ",X,F13.6,X,"(proton)")
   502 format("-->                                                                ",X,F13.6,X,"(neutron)")
   v=v*sr*sr
   w=w*sr*sr
   write (711, 501) v
   write (711, 502) w
   !---Weisskopf estimate without A-dependent factor
   ff1 = ac**(2.d0/3.d0)
   if (nem.eq.8) then!--- electric modes.
      ff1 = ff1**ilam
   else if (nem.eq.9) then!--- magnetic modes.
      ff1 = ff1**(ilam-1)
   end if
   write (711, *) " Same but without A-dep. factor f_J =x^{J}   for EJ mode,"
   write (711, *) "   where x=A^{2/3},        or   f_J =x^{J-1} for MJ mode:"
   503 format("-->                      B_{QJ} /f_J      = ",X,F13.6,X,"(proton)")
   504 format("-->                                         ",X,F13.6,X,"(neutron)")
   write (711, 503) v/ff1
   write (711, 504) w/ff1
   !--- result-2
   write (711, *) ""
   write (711, *) "<><<><><><><<><><><><<><><><><<><><><><<><><><><<><><><><<><><>"
   write (711, *) " [C] VS Weisskopf estimate (2), where"
   write (711, *) " radial integration is approximated as ~= (3*R**N)/(3 + ilam),"
   write (711, *) " where R=1.2*ac**(1/3) and N=ilam (ilam-1) for E (M),"
   write (711, *) " whereas the angular part is computed commonly as in [A]:"
   !--- angular part of <f || Q || i>
   qp=1.d0 ; qn=1.d0 ; p=0.d0
   if (nem.eq.8) then
      p = sqrt( dble((2*ilam+1)*(ji+1)*(jf+1))*0.25d0/pi  )
      qp = p*ang_ele(lf,jf, li,ji, ilam)
      qn = qp!--- where the charge of neutron is assumed as e, too.
   else if (nem.eq.9) then
      p = sqrt( dble(2*ilam+1)*0.25d0/pi  )
      qp = p*( gl_p*dbmmag_l_RS(jf,lf, ilam, ji,li) +gs_p*dbmmag_s_RS(jf,lf, ilam, ji,li) )
      qn = p*( gl_n*dbmmag_l_RS(jf,lf, ilam, ji,li) +gs_n*dbmmag_s_RS(jf,lf, ilam, ji,li) )
   end if
   !--- radial integration
   r = 1.2d0*ac**(1.d0/3.d0)
   if (nem.eq.8) then!--- electric modes.
      vp = r**(ilam+ilam)
   else if (nem.eq.9) then!--- magnetic modes.
      vp = r**(ilam+ilam-2)
   end if
   vp = vp*(3.d0/(3.d0 +dble(ilam)))**2
   vp = sqrt(vp)
   vn = vp
   !--- where the charge of neutron is assumed as e, too.
   v = vp*qp
   w = vn*qn
   write (711, 501) (v)**2/dble(ji+1)
   write (711, 502) (w)**2/dble(ji+1)
   write (711, *) " angular part of <f|Q|i> =", qp, " (prot.)"
   write (711, *) "                         =", qn, " (neut.)"
   write (711, *) " ---"
   write (711, *) " radial part of <f|Q|i>, approximated as"
   write (711, *) " ~= (3*R**N)/(3 + ilam) =", vp, " (prot.)"
   write (711, *) "                        =", vn, " (neut.)"
   write (711, *) ""
   765  continue
   close (711)
   end subroutine Weisskopf
   !***********************************************************************
   subroutine rectime(id)
   implicit none
   integer, intent(in) :: id
   integer, dimension (8) :: md
   call date_and_time(values = md)
   write (711, *) ''
   write (711, '("*** time -",i2,2x,i4,"/",i2,"/",i2,", ",i2,":",i2,":",i2," ***")') &
          id, md(1),md(2),md(3),md(5),md(6),md(7)
   write (711, *) ''
   end subroutine rectime
   !***********************************************************************
   subroutine units
     implicit none
     write( 711,*) "-----------------------------------------------------------------"
     write( 711,*) "Units for quantities in the CGS-Gauss system of units."
     write( 711,*) "- Electric mode with J: B_{EJ} in             [e^2.fm^(2J)],"
     write( 711,*) "                        angular part in       [e],"
     write( 711,*) "                        radial intergation in [fm^J]."
     write( 711,*) "- Magnetic mode with J: B_{MJ} in             [\mu^2.fm^(2J-2)],"
     write( 711,*) "                        angular part in       [\mu],"
     write( 711,*) "                        radial intergation in [fm^(J-1)],"
     write( 711,*) "               where \mu^2 ~= 1.102 x 10^{-2} [e^2.fm^2]."
     write( 711,*) "-----------------------------------------------------------------"
   end subroutine units
   !***********************************************************************
end module mdl_099_summary
!///////////////////////////////////////////////////////////////////////////////
         program HONTAI
         use mdl_099_summary
         implicit none
         real :: c1, c2
         integer :: it1, it2, it_rate, it_max, idiff
         call cpu_time( c1 )
         call system_clock(it1)
         call SPEM
!<>         call Weisskopf
         call system_clock(it2, it_rate, it_max)
         if ( it2 < it1 ) then
           idiff = (it_max - it1) + it2 + 1
         else
           idiff = it2 - it1
         endif
         call cpu_time( c2 )
         write(6,   *) "system time :", idiff/it_rate, "seconds."
         write(6,   *) "cpu time    :", c2-c1, "seconds."
!!!         open( 129, file='Clock.dat', action='write')
!!!            write(129, *) "system time :", idiff/it_rate, "seconds."
!!!            write(129, *) "cpu time    :", c2-c1, "seconds."
!!!         close(129)
         end program HONTAI
