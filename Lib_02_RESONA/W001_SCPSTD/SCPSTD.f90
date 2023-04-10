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
Module mdl_006_cwf!Mod_cwf
   Public :: Dfcoul
   Private :: Jflgam, Yfclen, Yfasym, Yfireg, Yfrica, Dfcz0, Jfdelg
!======================================================================!
!     * no implicit definition                                         !
!     * all internal functions are used as generic name                !
!     * modified constants(attached 'd0')                              !
!----------------------------------------------------------------------!
!     dfcoul = Subroutine for the Coulomb wave function                !
!                                                                      !
!     eta  --- Sommerfeld parameter. real(8), intent(in)               !
!     rho  --- Independent variable. real(8), intent(in)               !
!     fcw(0:L) --- Regular Coulomb wave function.                      !
!                  real(8), intent(out)                                !
!     fpcw(0:L) --- Derivative of regular Coulomb wave function.       !
!                  real(8), intent(out)                                !
!     gcw(0:L)  --- Irregular Coulomb wave function.                   !
!                   real(8), intent(out)                               !
!     gpcw(0:L) --- Derivative of irregular Coulomb wave function.     !
!                   real(8), intent(out)                               !
!     sigmad(0:L) --- Coulomb phase shift. real(8), intent(out)        !
!     L --- Angular momentum. We can compute above quantities up to    !
!           this angular momentum. integer, intent(in)                 !
!     iexp(0:L) --- ???. integer, intent(out)                          !
!----------------------------------------------------------------------!
!======================================================================!

Contains

!***********************************************************************
   Subroutine Dfcoul(Eta, Ro, F, Fp, G, Gp, Sigma, L, Iexp)
!***********************************************************************
      Implicit None
      Integer, Intent (In) :: L
      Integer, Intent (Out) :: Iexp(L+1)
      Double Precision, Intent (In) :: Eta, Ro
      Double Precision, Intent (Out), Dimension (L+1) :: F, Fp, G, Gp, Sigma
      Double Precision :: Depi = 6.283185307179586D0
      Double Precision :: Etac, F0, Fp0, G0, Gp0, Z, Roinf, Zig, Zag, Sig, Xm
      Double Precision :: Ftest, Fptest, Fi, Fpi, Fl, Fpl, Fact, Factp, A, B, S
      Integer :: I, J, L1, Ind, Linf, Lin, Ipi, Lmax, Indice, I1, I2

      Etac = Eta*Eta
      L1 = L + 1
      Call Dfcz0(Eta, Ro, F0, Fp0, G0, Gp0, S, I)
      F(1) = F0
      Fp(1) = Fp0
      G(1) = G0
      Gp(1) = Gp0
      Iexp(1) = I
      Sigma(1) = S
      If (L) 100, 100, 110
100   Return
110   Linf = 0
      Ind = 0
      If ((Eta>0) .And. (Ro<(Eta+Eta))) Go To 210
      Z = Eta + Sqrt(Etac+Dble(L*(L+1)))
      If (Ro>=Z) Go To 170
120   Roinf = Eta + Sqrt(Etac+Dble(Linf*(Linf+1)))
      If (Ro-Roinf) 150, 130, 130
130   If (Linf-L) 140, 160, 160
140   Linf = Linf + 1
      Go To 120
150   Ind = 1
160   Lin = Linf + 1
170   Xm = 1.D0
      If (Ind==0) Lin = L1
      Do J = 2, Lin
         Zig = (Sqrt(Etac+Xm*Xm))/Xm
         Zag = Eta/Xm + Xm/Ro
         F(J) = (Zag*F(J-1)-Fp(J-1))/Zig
         Fp(J) = Zig*F(J-1) - Zag*F(J)
         G(J) = (Zag*G(J-1)-Gp(J-1))/Zig
         Gp(J) = Zig*G(J-1) - Zag*G(J)
         Iexp(J) = I
         Sig = Sigma(J-1) + Atan(Eta/(J-1))
         Ipi = Sig/Depi
         Sig = Sig - Ipi*Depi
         If (Sig) 180, 200, 190
180      If (Sig<(-Depi/2.D0)) Sig = Sig + Depi
         Go To 200
190      If (Sig>(Depi/2.D0)) Sig = Sig - Depi
200      Sigma(J) = Sig
         Xm = Xm + 1.D0
      End Do
      If (Ind==0) Return
      Go To 220
210   Lin = 1
220   Ftest = F(Lin)
      Fptest = Fp(Lin)
      Lmax = Linf + 25 + Int(5.D0*Abs(Eta))
      If (Lmax-L) 230, 240, 240
230   Lmax = L
240   Fi = 1.D0
      Fpi = 1.D0
250   Xm = Lmax + 1
      Zig = (Sqrt(Etac+Xm*Xm))/Xm
      Zag = Eta/Xm + Xm/Ro
      Fl = (Zag*Fi+Fpi)/Zig
      Fpl = Zag*Fl - Zig*Fi
      If (Abs(Fl)-1.D15) 260, 270, 270
260   If (Abs(Fpl)-1.D15) 280, 270, 270
270   Fl = Fl*1.D-15
      Fpl = Fpl*1.D-15
280   Fi = Fl
      Fpi = Fpl
      If (Lmax-L) 300, 300, 290
290   Lmax = Lmax - 1
      Go To 250
300   F(Lmax+1) = Fl
      Fp(Lmax+1) = Fpl
      If (Lmax-Linf) 320, 320, 310
310   Go To 290
320   Fact = Ftest/F(Lin)
      Factp = Fptest/Fp(Lin)
      Indice = I/60
      Xm = Linf
      Do J = Lin, L1
         F(J) = F(J)*Fact
         Fp(J) = Fp(J)*Factp
         If (J==1) Go To 420
         Zig = (Sqrt(Etac+Xm*Xm))/Xm
         Zag = Eta/Xm + Xm/Ro
         G(J) = (Zag*G(J-1)-Gp(J-1))/Zig
         Gp(J) = Zig*G(J-1) - Zag*G(J)
         If (Abs(G(J))-1.D60) 330, 340, 340
330      If (Abs(Gp(J))-1.D60) 350, 340, 340
340      G(J) = G(J)/1.D60
         Gp(J) = Gp(J)/1.D60
         Indice = Indice + 1
350      Iexp(J) = Indice*60
         A = Fp(J)*G(J)
         B = -F(J)*Gp(J)
         If (A-1.D0) 360, 370, 370
360      I1 = Int(Log10(A))
         I2 = Int(Log10(B))
         Go To 380
370      I1 = Int(Log10(A)) + 1
         I2 = Int(Log10(B)) + 1
380      F(J) = F(J)*1.D1**(-I2)
         Fp(J) = Fp(J)*1.D1**(-I1)
         Sig = Sigma(J-1) + Atan(Eta/(J-1))
         Ipi = Sig/Depi
         Sig = Sig - Ipi*Depi
         If (Sig) 390, 410, 400
390      If (Sig<(-Depi/2.D0)) Sig = Sig + Depi
         Go To 410
400      If (Sig>(Depi/2.D0)) Sig = Sig - Depi
410      Sigma(J) = Sig
420      Xm = Xm + 1.D0
      End Do
      Return
   End Subroutine Dfcoul

!***********************************************************************
   Subroutine Dfcz0(Eta, Ro, F0, Fp0, G0, Gp0, Sigma0, Iexp)
!***********************************************************************
      Implicit None
      Integer, Intent (Out) :: Iexp
      Integer :: I, J, Ii, N, M, Ntruc
      Double Precision, Intent (In) :: Eta, Ro
      Double Precision, Intent (Out) :: F0, Fp0, G0, Gp0, Sigma0
      Double Precision, Dimension (110) :: A1, A2, B1, B2
      Double Precision :: Pi = 3.141592653589793D0
      Double Precision :: Borne, Tra, Ro2, Etap, Pieta, Pieta2, B, U0, U1, U2
      Double Precision :: Xn, U, Xn1, Up, H, Dh, Di, S, Q, D, C, X1, X2, T1, T2
      Double Precision :: Derive, Ti, Result, Z, Y, X

      If (Ro) 100, 100, 120
100   Continue
      Write (6, 110)
110   Format (' ro negatif ou nul **')
      F0 = 0.0D0
      Fp0 = 0.0D0
      G0 = 0.0D0
      Gp0 = 0.0D0
      Sigma0 = 0.0D0 !I added
      Iexp = 0
      Return
120   If (Eta-30.D0) 170, 170, 130
130   If (Abs(Eta)-5.D2) 140, 140, 150
140   Call Yfrica(Eta, Ro, F0, Fp0, G0, Gp0, Sigma0, Iexp)
      Return
150   Continue

      F0 = 0.0D0
      Fp0 = 0.0D0
      G0 = 0.0D0
      Gp0 = 0.0D0
      Sigma0 = 0.0D0 !I added
      Iexp = 0

      Write (6, 160)
160   Format (' valeur absolue de eta supe-&eu-e a 500 **')
      Return
170   If (Eta+8.D0) 130, 180, 180
180   If (Eta) 200, 190, 200
190   F0 = Sin(Ro)
      G0 = Cos(Ro)
      Fp0 = G0
      Gp0 = -F0
      Iexp = 0
      Sigma0 = 0.D0
      Return
200   Borne = 1.666666666666667D0*Abs(Eta) + 7.5D0
      If (Ro-Borne) 220, 210, 210
210   Call Yfasym(Eta, Ro, F0, Fp0, G0, Gp0, Sigma0, Iexp)
      Return
220   If (Eta-10.D0) 230, 250, 250
230   If (Eta) 260, 260, 240
240   If (Ro-2.D0) 310, 260, 260
250   If (Eta-(5.D0*Ro+6.D1)/7.D0) 260, 260, 310
260   Call Yfasym(Eta, Borne, F0, Fp0, G0, Gp0, Sigma0, Iexp)
      H = Borne
      Dh = F0/H
      If (Eta) 270, 280, 280
270   N = -0.5D0*Eta + 5.D0
      Go To 290
280   N = Eta/5.D0 + 5.D0
290   N = 5*(N+1)
      Z = 4.D0/H
      Y = 1.D0 - (Eta+Eta)*Z
      A1(N+2) = 1.D-55
      A1(N+3) = 0.D0
      A1(N+4) = 1.D-64
      B1(N+3) = 1.D-50
      B1(N+4) = 1.D-68
      A2(N+2) = 0.D0
      A2(N+3) = 1.D-74
      A2(N+4) = 1.D-53
      B2(N+3) = 0.D0
      B2(N+4) = 1.D-66
      M = N + 2
      Di = N
      Do Ii = 2, M
         I = M - Ii + 2
         B1(I) = B1(I+2) + Z*(Di+1.D0)*A1(I+1)
         S = A1(I+2) + Y*(A1(I+1)-A1(I))
         Q = (Di+2.D0)*B1(I) + (Di-1.D0)*B1(I+1)
         A1(I-1) = S - Z*Q
         B2(I) = B2(I+2) + Z*(Di+1.D0)*A2(I+1)
         S = A2(I+2) + Y*(A2(I+1)-A2(I))
         Q = (Di+2.D0)*B2(I) + (Di-1.D0)*B2(I+1)
         A2(I-1) = S - Z*Q
         If (I>=N) Go To 300
         D = -(B2(I+2)+B2(I))/(B1(I+2)+B1(I))
         Do J = I, M
            A2(J) = A2(J) + D*A1(J)
            B2(J) = B2(J) + D*B1(J)
         End Do
         A2(I-1) = A2(I-1) + D*A1(I-1)
300      Di = Di - 1.D0
      End Do
      Q = A1(3) - A1(1)
      C = A2(3) - A2(1)
      C = Q/C
      X1 = 0.5D0*(A1(2)-C*A2(2))
      Do I = 3, M
         X1 = X1 + A1(I) - C*A2(I)
      End Do
      X1 = Dh/X1
      X2 = -C*X1
      Do I = 2, M
         B1(I) = X1*B1(I) + X2*B2(I)
         A1(I) = X1*A1(I) + X2*A2(I)
      End Do
      A1(1) = X1*A1(1) + X2*A2(1)
      B1(1) = 0.D0
      X = Ro/H
      Y = 2.D0*X - 1.D0
      T1 = 1.D0
      T2 = Y
      Result = 0.5D0*A1(2) + Y*A1(3)
      Derive = 0.5D0*B1(2) + Y*B1(3)
      Do I = 2, N
         Ti = 2.D0*Y*T2 - T1
         T1 = T2
         T2 = Ti
         Result = Result + Ti*A1(I+2)
         Derive = Derive + Ti*B1(I+2)
      End Do
      F0 = Result*Ro
      Fp0 = Derive*Ro + Result
      Go To 330
!   serie origine reguliere
310   Pi = 3.141592653589793D0
      Call Jflgam(1.D0, Eta, Tra, Sigma0, Ntruc)
      Iexp = 0
      Ro2 = Ro*Ro
      Etap = Eta + Eta
      Pieta = Pi*Eta
      Pieta2 = 0.5D0*Pieta
      B = Exp(Pieta2)*Sqrt(Sinh(Pieta)/Pieta)
      U0 = 0.D0
      U1 = Ro
      U = U0 + U1
      Up = 1.D0
      Xn = 2.D0
      Do N = 2, 10000
         Xn1 = Xn*(Xn-1.D0)
         U2 = (Etap*Ro*U1-Ro2*U0)/Xn1
         U = U + U2
         Up = Up + Xn*U2/Ro
         If (Abs(U2/U)<1.D-10) Go To 320
         U0 = U1
         U1 = U2
         Xn = Xn + 1.D0
      End Do
320   F0 = U/B
      Fp0 = Up/B
330   Call Yfireg(Eta, Ro, G0, Gp0)
      Return
   End Subroutine Dfcz0

!***********************************************************************
   Subroutine Jflgam(Xd, Yd, Par, Pai, Nbchif)
!***********************************************************************
      Implicit None
      Integer :: I, K, Kr
      Integer, Intent (Out) :: Nbchif
      Double Precision, Intent (In) :: Xd, Yd
      Double Precision, Intent (Out) :: Par, Pai
      Double Precision, Parameter :: Hlo2pi = 0.918938533204672D0
      Double Precision, Parameter :: Pi = 3.141592653589793D0
      Double Precision, Parameter :: Pis2 = 1.570796326794897D0
      Double Precision, Parameter :: Pis4 = 0.785398163397448D0
      Double Precision, Parameter :: Alo2pi = 1.837877066409345D0
      Double Precision, Parameter :: Rac2 = 0.3465735902799726D0
      Double Precision, Parameter :: Depi = 6.283185307179586D0
      Double Precision, Parameter :: Alopi = 1.1447298858494002D0
      Double Precision, Parameter :: Supint = 2147483647.0D0
      Double Precision :: Test(7) = (/ 2.9152D+7, 2.2958D+3, 1.4124D+2, 3.9522D+1, 19.6611D0, 12.791D0, -10.D0 /)
      Double Precision :: C(6) = (/ 8.333333333333333D-2, -2.777777777777777D-3, 7.936507936507937D-4, &
         -5.952380952380952D-4, 8.417508417508418D-4, -1.917526917526918D-3 /)
      Double Precision :: X, Y, U, V, Tra, Tra1, Cosi, Sini, Cos2i, Sin2i
      Double Precision :: Zmod, Xx

      Nbchif = 15
      X = Abs(Xd)
      Xx = X
      If (Yd) 100, 470, 100
100   Y = Abs(Yd)
      Kr = 1
      I = Mod(10.99D0-X, Supint)
!     translation
      If (I) 120, 120, 110
110   Tra = I
      X = X + Tra
!     logarithme(x+iy) (x,y positifs)
120   If (X-Y) 130, 150, 190
130   Tra1 = X/Y
      If (Tra1) 140, 140, 160
140   U = Log(Y)
      V = Pis2
      Go To 180
150   U = Rac2 + Log(X)
      V = Pis4
      Go To 180
160   Tra = Y*Sqrt(1.D0+Tra1*Tra1)
      Tra1 = Y/X
170   U = Log(Tra)
      V = Atan(Tra1)
180   Go To (200, 250, 350), Kr
190   Tra1 = Y/X
      Tra = X*Sqrt(1.D0+Tra1*Tra1)
      Go To 170
200   Kr = 2
!     developpement asymptotique ( x superieur a 10 )
      Tra = X - 0.5D0
      Pai = V*Tra + Y*(U-1.D0)
      Par = -X + Hlo2pi + U*Tra - V*Y
      Zmod = X + Y
      If (Zmod-Test(1)) 210, 210, 260
210   Tra = X*X + Y*Y
      Cosi = X/Tra
      Sini = Y/Tra
      Sin2i = (Sini*Cosi) + (Sini*Cosi)
      Cos2i = (Cosi+Sini)*(Cosi-Sini)
      K = 1
      Go To 230
220   Tra = Cosi*Cos2i - Sini*Sin2i
      Sini = Sini*Cos2i + Cosi*Sin2i
      Cosi = Tra
230   Par = Par + C(K)*Cosi
      Pai = Pai - C(K)*Sini
      K = K + 1
      If (Zmod-Test(K)) 220, 220, 260
!     translation inverse
240   I = I - 1
      X = I
      X = Xx + X
      Go To 120
250   Par = Par - U
      Pai = Pai - V
260   If (I-1) 280, 270, 240
270   If (Xd) 240, 290, 240
!     controle du quadrant
280   If (Xd) 340, 290, 330
290   Tra = Pi*Y
      If (Tra-1.D-2) 300, 300, 310
300   Par = Tra*(2.D0+Tra*(-2.D0+Tra*(1.333333333333333D0+Tra*(-0.6666666666666666D0+Tra*(0.2666666666666666D0+Tra*(- &
         0.08888888888888888D0+Tra*0.02539682539682540D0))))))
      Tra1 = -Log(Y) - Log(Par)
      Go To 320
310   Par = 1.D0 - Exp(-Tra-Tra)
      Tra1 = -Log(Y*Par)
320   Par = 0.5D0*(Alo2pi-Tra+Tra1)
      Pai = Pai - Pis2
330   If (Yd) 400, 410, 410
!     x+iy change en -x-iy
340   Tra = Pi*Y
      Par = Alo2pi - U - Par - Tra
      Pai = Pi - V - Pai
      Tra = Exp(-Tra-Tra)
      X = Pi*Mod(X, 2.D0)
      Sini = (1.D0-Tra)*Cos(X)
      Cosi = (1.D0+Tra)*Sin(X)
      Kr = 3
      X = Abs(Cosi)
      Y = Abs(Sini)
      Go To 120
350   If (Cosi) 360, 370, 370
360   V = Pi - Sign(V, Sini)
      Go To 390
370   If (Sini) 380, 390, 390
380   V = -V
390   Par = Par - U
      Pai = Pai - V
      If (Yd) 410, 410, 400
400   Pai = -Pai
!     argument dans -pi,pi
410   Tra = Abs(Pai/Depi)
      If (Tra-1.D+15) 430, 420, 420
420   Nbchif = 0
      Pai = 0.D0
      Go To 730
430   If (Tra-1.D0) 450, 450, 440
440   Nbchif = Log10(Tra)
      Nbchif = 14 - Nbchif
      Tra = Mod(Tra, Supint)
      Pai = Mod(Tra, 1.D0)*Sign(Depi, Pai)
450   If (Abs(Pai)-Pi) 730, 730, 460
460   Pai = Pai - Sign(Depi, Pai)
      Go To 730
!     jflgam reel
470   Pai = 0.D0
      If (Xd) 500, 480, 590
!     conditions d existence
480   Write (6, 490)
490   Format (' jflgam(0) est infini')
      Go To 580
500   If (X-4503599627370496.D0) 530, 510, 510
510   Write (6, 520)
520   Format (' argument de jflgam trop grand')
      Go To 580
530   Y = Mod(X, Supint)
      If (Y) 540, 560, 540
540   If (Y-0.99D0) 590, 590, 550
550   Tra = Int(Y+0.1D0)
      If (Abs(Y-Tra)-5.D-15) 560, 560, 590
560   Write (6, 570)
570   Format (' jflgam (-entier) est infini')
580   Par = 1.D+74
      Nbchif = 0
      Go To 730
!     translation
590   I = Mod(10.99D0-X, Supint)
      If (I) 610, 610, 600
600   Tra = I
      X = X + Tra
!     developpement asymptotique
610   Y = Log(X)
      Par = -X + Hlo2pi + Y*(X-0.5D0)
      If (X-Test(1)) 620, 620, 690
620   Cosi = 1.D0/X
      Cos2i = Cosi*Cosi
      K = 1
      Go To 640
630   Cosi = Cosi*Cos2i
640   Par = Par + C(K)*Cosi
      K = K + 1
      If (X-Test(K)) 630, 630, 670
!     translation inverse
650   X = X - 1.D0
660   Y = Log(X)
      Par = Par - Y
      I = I - 1
670   If (I-1) 690, 680, 650
680   X = Abs(Xd)
      Go To 660
!     x negatif
690   If (Xd) 700, 730, 730
700   Par = Alopi - Par - Y
      Y = Pi*Mod(X, 2.D0)
      Y = -Sin(Y)
      If (Y) 710, 560, 720
710   Y = -Y
      Pai = Pi
720   Par = Par - Log(Y)
      Entry Jflgv1
730   Return
   End Subroutine Jflgam

!***********************************************************************
   Subroutine Yfclen(Eta, Ro, U, Up, V, Vp, Sigma0, Idiv, Nn)
!***********************************************************************
      Implicit None
      Integer, Intent (In) :: Nn
      Integer, Intent (Out) :: Idiv
      Integer :: Indg, Nb, N, M, K, J
      Double Precision, Intent (In) :: Eta, Ro, Sigma0
      Double Precision, Intent (Out) :: U, Up, V, Vp
      Double Precision :: Pi
      Double Precision :: Xa
      Double Precision :: Ro2, Etap, Pieta, Z, Pieta2, Par, Pai, R, X, Xx
      Double Precision :: Z1, U0, V0, V1, Xn, U2, V2, Pp, W, P, U1, Xn1, Wp
      Double Precision :: Eta2, T0, T1, Tj, Tm, Tl, Tk
      Complex (Kind=Kind(0D0)) :: Fa, Ak, Am, Al, Bj, Bn, Bm, Bl, Bk, F, G, Gp
      Complex (Kind=Kind(0D0)) :: C1, C2, C3, C4, C5, C6, D, Hk, Axpo
!
      If (Nn==1) Go To 180
!
      Eta2 = Eta*Eta
      Fa = Cmplx(1.D0, Eta, Kind(0D0))
      M = 0.25D0*Eta + 4.D1
!
!          polynomes de tchebichev jusqu'au rang m
!
      K = M + 2
      X = (Eta+Eta)/Ro
      Xx = X + X - 1.D0
      T0 = 1.D0
      T1 = Xx
      Xx = Xx + Xx
      Do J = 2, K
         Tj = Xx*T1 - T0
         T0 = T1
         T1 = Tj
      End Do
      Tm = T1
      Tl = T0
!
!          initialisation
!
      Am = (0.D0, 0.D0)
      Al = (1.D-40, 1.D-40)
      Bn = (0.D0, 0.D0)
      Bm = (1.D-40, 1.D-40)
      Bl = (0.D0, 0.D0)
      Bk = Cmplx(4.D0*Dble(M+1), 0.D0, Kind(0D0))*Al + Bm
      F = (0.D0, 0.D0)
      G = (0.D0, 0.D0)
      Gp = (0.D0, 0.D0)
!
!          recurrence descendante
!
      K = M
100   R = K
      Tk = Xx*Tl - Tm
      Tm = Tl
      Tl = Tk
      Hk = Cmplx(Tk, 0.D0, Kind(0D0))
      C1 = Cmplx(R*(R+1.D0)-Eta2, Eta*(R+R+1.D0), Kind(0D0))
      C2 = (4.D0, 0.D0)*Cmplx(R+1.D0, 0.D0, Kind(0D0))
      C2 = C2*Cmplx(-R-1.D0, Eta*3.D0, Kind(0D0))
      C3 = Fa*Cmplx(-R-R-4.D0, Eta, Kind(0D0))
      C4 = Cmplx((7.D0*R+5.D0)/4.D0, 0.D0, Kind(0D0))
      C5 = Cmplx(R+R+2.D0, 0.D0, Kind(0D0))
      C6 = Cmplx((R+3.D0)/4.D0, 0.D0, Kind(0D0))
      Ak = (C2*Al+C3*Am-C4*Bl-C5*Bm-C6*Bn)/C1
      J = K/2
      J = K - J - J
      If (J) 110, 120, 110
110   F = F - Ak
      Go To 130
120   F = F + Ak
130   Z = Abs(Ak)
      G = G + Hk*Ak
      Gp = Gp + Hk*Bk
!
!          f=a0/2-a1+a2-a3+a4-a5+...
!
!          congruence modulo 10**60
!
      If (Z-1.D60) 150, 140, 140
140   D = (1.D60, 0.D0)
      Ak = Ak/D
      Al = Al/D
      Am = Am/D
      Bk = Bk/D
      Bl = Bl/D
      Bm = Bm/D
      Bn = Bn/D
      F = F/D
      G = G/D
      Gp = Gp/D
150   If (K) 170, 170, 160
160   D = Cmplx(4.D0*R, 0.D0, Kind(0D0))
      Bj = D*Ak + Bl
      Am = Al
      Al = Ak
      Bn = Bm
      Bm = Bl
      Bl = Bk
      Bk = Bj
      K = K - 1
      Go To 100
!
!          normalisation et calcul de z(ro)
!
170   D = (0.5D0, 0.D0)*Ak
      F = F - D
      G = G - D
      Gp = Gp - (0.5D0, 0.D0)*Bk
      D = Cmplx(0.D0, -Eta*Log(2.D0)+Sigma0, Kind(0D0))
      Axpo = Exp(D)
      F = F/Axpo
      G = G/F
      Gp = Gp/F
!
!          calcul de f et g
!
      D = Cmplx(0.D0, Ro-Eta*Log(Ro), Kind(0D0))
      Axpo = Exp(D)
      D = G*Axpo
      Gp = Axpo*(Cmplx(0.D0,1.D0-Eta/Ro,Kind(0D0))*G-Cmplx(X/Ro,0.D0,Kind(0D0))*Gp)
      V = D
      D = (0.D0, -1.D0)*D
      U = D
      Vp = Gp
      Gp = (0.D0, -1.D0)*Gp
      Up = Gp
      Idiv = 0
      Return
!
!          serie origine
!
180   Pi = 3.141592653589793D0
      Xa = 0.577215664901533D0
      Ro2 = Ro*Ro
      Etap = Eta + Eta
      Pieta = Pi*Eta
      Z = 138.15510557964276D0
      Idiv = 0
      If (Abs(Pieta)-Z) 200, 200, 190
190   Indg = Int(Pieta/Z)
      Idiv = 60*Indg
      If (Eta<0) Idiv = 0
      Pieta = Pieta - Z*Dble(Indg)
200   Pieta2 = 0.5D0*Pieta
      P = Exp(Pieta2)*Sqrt(Sinh(Pieta)/Pieta)
      Call Jfdelg(1.D0, Eta, Par, Pai, Nb)
      Z1 = Etap*(Xa+Xa+Log(2.D0)-1.D0+Par)
      U0 = 0.D0
      U1 = Ro
      V0 = 1.D0
      V1 = Z1*Ro
      U = U0 + U1
      V = V0 + V1
      Up = 1.D0
      Vp = Z1
      Xn = 2.D0
      Do N = 2, 10000
         Xn1 = Xn*(Xn-1.D0)
         U2 = (Etap*Ro*U1-Ro2*U0)/Xn1
         U = U + U2
         V2 = (Etap*Ro*V1-Ro2*V0-Etap*(Xn+Xn-1.D0)*U2)/Xn1
         V = V + V2
         Up = Up + Xn*U2/Ro
         Vp = Vp + Xn*V2/Ro
         If (Abs(U2/U)>1.D-14) Go To 210
         If (Abs(V2/V)<=1.D-14) Go To 220
210      U0 = U1
         U1 = U2
         V0 = V1
         V1 = V2
         Xn = Xn + 1.D0
      End Do
220   Pp = V + Etap*U*Log(Ro)
      W = U/P
      Wp = Up/P
      V = P*Pp
      Vp = P*(Vp+Etap*(Up*Log(Ro)+U/Ro))
      U = W
      Up = Wp
      Return
   End Subroutine Yfclen

!***********************************************************************
   Subroutine Yfasym(Eta, Rau, Fo, Fpo, Go, Gpo, Sigo, Iexp)
!***********************************************************************
      Implicit None
      Integer, Intent (Out) :: Iexp
      Integer :: N, Ntruc
      Double Precision, Intent (In) :: Eta, Rau
      Double Precision, Intent (Out) :: Fo, Fpo, Go, Gpo, Sigo
      Double Precision :: Trb, Rau2, Etac, Tra, Ps, Pt, Gs, Gt, Sf, Sg, Spg
      Double Precision :: Spf, Denom, An, Bn, Ps1, Gs1, Pt1, Gt1, Test, Tetao

      Iexp = 0
      Trb = 0.D0
      Rau2 = Rau + Rau
      Etac = Eta*Eta
      Call Jflgam(1.D0, Eta, Tra, Sigo, Ntruc)
      N = 0
      Ps = 1.D0
      Gs = 0.D0
      Pt = 0.D0
      Gt = 1.D0 - Eta/Rau
      Sf = Ps
      Sg = Gs
      Spf = Pt
      Spg = Gt
100   Denom = Dble(N+1)*Rau2
      An = Dble(N+N+1)*Eta/Denom
      Bn = (Etac-Dble(N*(N+1)))/Denom
      Ps1 = An*Ps - Bn*Pt
      Gs1 = An*Gs - Bn*Gt - Ps1/Rau
      Pt1 = An*Pt + Bn*Ps
      Gt1 = An*Gt + Bn*Gs - Pt1/Rau
      Sf = Sf + Ps1
      Sg = Sg + Gs1
      Spf = Spf + Pt1
      Spg = Spg + Gt1
      N = N + 1
      If (N-17) 130, 110, 120
110   Tra = Ps*Ps + Pt*Pt
120   Trb = Ps1*Ps1 + Pt1*Pt1
      Test = Tra - Trb
      If (Test) 140, 140, 130
130   Ps = Ps1
      Gs = Gs1
      Pt = Pt1
      Gt = Gt1
      Tra = Trb
      Go To 100
140   Tetao = Rau - Eta*Log(Rau2) + Sigo
      Tra = Sin(Tetao)
      Trb = Cos(Tetao)
      Go = Sf*Trb - Spf*Tra
      Gpo = Sg*Trb - Spg*Tra
      Fo = Spf*Trb + Sf*Tra
      Fpo = Spg*Trb + Sg*Tra
      Return
   End Subroutine Yfasym

!***********************************************************************
   Subroutine Yfireg(Eta, Ro, G0, Gp0)
!***********************************************************************
      Implicit None
      Integer :: Iexp, N, Nb
      Double Precision, Intent (In) :: Eta, Ro
      Double Precision, Intent (Out) :: G0, Gp0
      Double Precision :: Rau0, F0, Fp0, Sigma0, X, X2, X3, Unr, Etr0, U0, U1
      Double Precision :: U2, S, V1, V2, T, Xn, Xn1, U3, V3, Pi, Ga, Eta2, Ro2
      Double Precision :: Etap, Pieta, Pieta2, B, Par, Pai, C1, V0, U, V, Up
      Double Precision :: Vp, Gp

      If (Eta<=0.D0) Go To 160
      If (Eta<=3.D0) Go To 170
      If (Eta<=1.D1) Go To 180
      If (Eta<=18.D0) Go To 190
      If (Eta<=22.D0) Go To 200
      If (Ro<=0.3D0+(3.D1-Eta)/8.D1) Go To 130
!   serie de taylor depart rau0
100   Continue
      Rau0 = 1.666666666666667D0*Abs(Eta) + 7.5D0
      Call Yfasym(Eta, Rau0, F0, Fp0, G0, Gp0, Sigma0, Iexp)
      X = Rau0 - Ro
      X2 = X*X
      X3 = X*X2
      Unr = 1.D0/Rau0
      Etr0 = 1.D0 - 2.D0*Eta*Unr
      U0 = G0
      U1 = -X*Gp0
      U2 = -0.5D0*Etr0*X2*U0
      S = U0 + U1 + U2
      V1 = U1/X
      V2 = 2.D0*U2/X
      T = V1 + V2
      Xn = 3.D0
      Do N = 3, 10000
!     n=n
         Xn1 = Xn - 1.D0
         Xn1 = Xn*Xn1
         U3 = X*U2*Unr*(1.D0-2.D0/Xn) - Etr0*U1*X2/Xn1 + X3*U0*Unr/Xn1
         S = S + U3
         V3 = Xn*U3/X
         T = T + V3
         If (Abs(U3/S)>1.D-11) Go To 110
         If (Abs(V3/T)<=1.D-11) Go To 120
110      U0 = U1
         U1 = U2
         U2 = U3
         Xn = Xn + 1.D0
      End Do
120   G0 = S
      Gp0 = -T
      Return
!   serie  origine
130   Continue
      Pi = 3.141592653589793D0
      Ga = 0.577215664901533D0
      Eta2 = Eta*Eta
      Ro2 = Ro*Ro
      Etap = Eta + Eta
      Pieta = Pi*Eta
      Pieta2 = 0.5D0*Pieta
      B = Exp(Pieta2)*Sqrt(Sinh(Pieta)/Pieta)
      Call Jfdelg(1.D0, Eta, Par, Pai, Nb)
      C1 = Etap*(Ga+Ga+Log(2.D0)-1.D0+Par)
      U0 = 0.D0
      U1 = Ro
      V0 = 1.D0
      V1 = C1*Ro
      U = U0 + U1
      V = V0 + V1
      Up = 1.D0
      Vp = C1
      Xn = 2.D0
      Do N = 2, 10000
         Xn1 = Xn*(Xn-1.D0)
         U2 = (Etap*Ro*U1-Ro2*U0)/Xn1
         U = U + U2
         V2 = (Etap*Ro*V1-Ro2*V0-Etap*(Xn+Xn-1.D0)*U2)/Xn1
         V = V + V2
         Up = Up + Xn*U2/Ro
         Vp = Vp + Xn*V2/Ro
         If (Abs(U2/U)>1.D-14) Go To 140
         If (Abs(V2/V)<=1.D-14) Go To 150
140      U0 = U1
         U1 = U2
         V0 = V1
         V1 = V2
         Xn = Xn + 1.D0
      End Do
150   Gp = V + Etap*U*Log(Ro)
      G0 = B*Gp
      Gp0 = B*(Vp+Etap*(Up*Log(Ro)+U/Ro))
      Return
160   If (Ro<=0.5D0*Eta+9.D0) Go To 130
      Go To 100
170   If (Ro<=2.25D0+7.35D0*(3.D0-Eta)) Go To 130
      Go To 100
180   If (Ro<=1.2D0+1.5D-1*(1.D1-Eta)) Go To 130
      Go To 100
190   If (Ro<=0.6D0+0.75D-1*(18.D0-Eta)) Go To 130
      Go To 100
200   If (Ro<=0.4D0+0.5D-1*(22.D0-Eta)) Go To 130
      Go To 100
   End Subroutine Yfireg

!***********************************************************************
   Subroutine Yfrica(Eta, Rau, Fo, Fpo, Go, Gpo, Sigma0, Idiv)
!***********************************************************************
      Implicit None
      Integer, Intent (Out) :: Idiv
      Integer :: N, Ind, Jnd, Ig, Nn, Indice, Indg
      Double Precision, Intent (In) :: Eta
      Double Precision :: Ro
      Double Precision, Intent (Out) :: Fo, Fpo, Go, Gpo, Sigma0
!
!        coefficients riccati
!
      Double Precision :: G61 = 0.1159057617187498D-1
      Double Precision :: G62 = 0.3863525390624998D-1
      Double Precision :: G63 = 0.4660034179687498D-1
      Double Precision :: G64 = 0.4858398437499998D-1
      Double Precision :: G65 = 0.1156514485677080D1
      Double Precision :: G66 = 0.5687475585937496D1
      Double Precision :: G67 = 0.1323888288225445D2
      Double Precision :: G68 = 0.1713083224826384D2
      Double Precision :: G69 = 0.1269003295898436D2
      Double Precision :: G610 = 0.5055236816406248D1
      Double Precision :: G611 = 0.8425394694010415D0
      Double Precision :: G81 = 0.1851092066083633D-01
      Double Precision :: G82 = 0.8638429641723630D-01
      Double Precision :: G83 = 0.1564757823944092D0
      Double Precision :: G84 = 0.1430139541625977D0
      Double Precision :: G85 = 0.1924622058868408D0
      Double Precision :: G86 = 0.8500803152720129D1
      Double Precision :: G87 = 0.7265429720878595D2
      Double Precision :: G88 = 0.3057942376817972D3
      Double Precision :: G89 = 0.7699689544836672D3
      Double Precision :: G810 = 0.1254157054424285D4
      Double Precision :: G811 = 0.1361719536066055D4
      Double Precision :: G812 = 0.9831831171035763D3
      Double Precision :: G813 = 0.4547869927883148D3
      Double Precision :: G814 = 0.1222640538215636D3
      Double Precision :: G815 = 0.1455524450256709D2
      Double Precision :: Gp61 = 0.2897644042968748D-01
      Double Precision :: Gp62 = 0.2318115234375000D0
      Double Precision :: Gp63 = 0.8056640625000000D0
      Double Precision :: Gp64 = 0.1601562499999998D1
      Double Precision :: Gp65 = 0.3046875000000000D0
      Double Precision :: Gp66 = 0.5624999999999998D1
      Double Precision :: Gp81 = 0.6478822231292720D-01
      Double Precision :: Gp82 = 0.6910743713378906D0
      Double Precision :: Gp83 = 0.3322952270507811D1
      Double Precision :: Gp84 = 0.9483032226562498D1
      Double Precision :: Gp85 = 0.1769653320312499D2
      Double Precision :: Gp86 = 0.3478710937499998D2
      Double Precision :: Gp87 = 0.5020312499999999D2
      Double Precision :: Gp88 = 0.7874999999999999D2
      Double Precision :: Q(5) = (/ 0.4959570165D-1, 0.8888888889D-2, 0.2455199181D-2, 0.9108958061D-3, &
         0.2534684115D-3 /)
      Double Precision :: Qp(5) = (/ 0.1728260369D0, 0.3174603174D-3, 0.3581214850D-2, 0.3117824680D-3, &
         0.9073966427D-3 /)
      Double Precision :: Tra, Rau2, Rauc, Etac, Eta2, Etaro, Etaro2, Pieta
      Double Precision :: Rau0, X, U, X2, Ru, Rx, Tre, Trb, Fi, Tr1, Tr2, Tr3, S
      Double Precision :: Tr4, Tr5, Tr6, Tr7, Tr8, Psip, Xxx, Psi, Fip, Et
      Double Precision :: Etad, Et0, Et1, Et2, Et3, Et4, Et5, X3, Unr, Etr0, U0
      Double Precision :: U1, U2, U3, V1, V2, V3, T, Xn, Xn1, Ho, Hpo, Trd, Trc
      Double Precision :: Rau

      Call Jflgam(1.D0, Eta, Tra, Sigma0, Ind)
      Rau2 = Rau + Rau
      Rauc = Rau*Rau
      Etac = Eta*Eta
      Eta2 = Eta + Eta
      Etaro = Eta*Rau
      Etaro2 = Etaro + Etaro
      Pieta = 3.141592653589793D0*Eta
      Ind = 0
      Jnd = 0
      Ig = 0
      If (Eta) 100, 100, 120
100   If (-Etaro-14.0625D0) 150, 110, 110
110   Indice = 1

!             riccati 3
      Idiv = 0
      Go To 290
120   If (Abs(Rau-Eta2)<=1.D-9) Go To 280
      If (Rau-Eta2) 180, 280, 130
130   If (Rau-Eta2-2.D1*(Eta**0.25D0)) 160, 140, 140
140   Indice = 0
!             riccati  2
      Idiv = 0
      Go To 290
150   Nn = 1
      Go To 170
160   Nn = 0
170   Call Yfclen(Eta, Rau, Fo, Fpo, Go, Gpo, Sigma0, Idiv, Nn)
      Return
180   If (Etaro-12.D0) 150, 150, 190
190   Tra = Eta2 - 6.75D0*(Eta**0.4D0)
      If (Rau-Tra) 210, 210, 200
200   Ind = 1
      Jnd = 1
      Ro = Rau
      Rau = Tra
      Rau0 = Tra
!             riccati  1

210   X = Rau/Eta2
      U = (1.D0-X)/X
      X2 = X*X
      Ru = Sqrt(U)
      Rx = Sqrt(X)
      Tre = 1.D0/(U*Ru*Eta2)
      Trb = Tre*Tre
      Fi = (Sqrt((1.D0-X)*X)+Asin(Rx)-1.570796326794897D0)*Eta2
      Tr1 = -0.25D0*Log(U)
      Tr2 = -((9.D0*U+6.D0)*U+5.D0)/48.D0
      Tr3 = ((((-3.D0*U-4.D0)*U+6.D0)*U+12.D0)*U+5.D0)/64.D0
      Tr4 = -((((((U+2.D0)*945.D0*U+1395.D0)*U+12300.D0)*U+25191.D0)*U+19890.D0)*U+5525.D0)/46080.D0
      Tr5 = ((((((((-27.D0*U-72.D0)*U-68.D0)*U+360.D0)*U+2190.D0)*U+4808.D0)*U+5148.D0)*U+2712.D0)*U+565.D0)/2048.D0
      Tr6 = -(((((((((G61*U+G62)*U+G63)*U+G64)*U+G65)*U+G66)*U+G67)*U+G68)*U+G69)*U+G610)*U + G611
      Tr7 = ((((((((((((-81.D0*U-324.D0)*U-486.D0)*U-404.D0)*U+4509.D0)*U+52344.D0)*U+233436.D0)*U+567864.D0)*U+ &
         838521.D0)*U+775884.D0)*U+441450.D0)*U+141660.D0)*U+19675.D0)/6144.D0
      Tr8 = (((((((((((((G81*U+G82)*U+G83)*U+G84)*U+G85)*U+G86)*U+G87)*U+G88)*U+G89)*U+G810)*U+G811)*U+G812)*U+G813)*U+ &
         G814)*U + G815
      Psip = Psip + Tra
      Xxx = 138.1551055796428D0
      Fi = Fi + Tre*(Tr2+Trb*(Tr4+Trb*(Tr6+Trb*Tr8)))
      Psi = -Fi
      Indg = Int(Psi/Xxx)
      Idiv = 60*Indg
      Tra = Tr1 + Trb*(Tr3+Trb*(Tr5+Trb*Tr7))
      Fi = Fi + Tra
      Psi = Psi + Tra

      Fip = Ru*Eta2
      Tra = 1.D0/(X2*U)
      Tr1 = 0.25D0
      Tre = Tre/(X2*X2*U)
      Trb = Trb/(X2*X2)
      Tr2 = -(8.D0*X-3.D0)/32.D0
      Tr3 = ((24.D0*X-12.D0)*X+3.D0)/64.D0
      Tr4 = (((-1536.D0*X+704.D0)*X-336.D0)*X+63.D0)/2048.D0
      Tr5 = ((((1920.D0*X-576.D0)*X+504.D0)*X-180.D0)*X+27.D0)/1024.D0
      Tr6 = ((((-Gp66*X+Gp65)*X-Gp64)*X+Gp63)*X-Gp62)*X + Gp61
      Tr7 = -((((((-40320.D0*X-10560.D0)*X-13248.D0)*X+7560.D0)*X-3132.D0)*X+756.D0)*X-81.D0)/2048.D0
      Tr8 = -((((((Gp88*X+Gp87)*X+Gp86)*X-Gp85)*X+Gp84)*X-Gp83)*X+Gp82)*X - Gp81
      Fip = Fip + Tre*(Tr2+Trb*(Tr4+Trb*(Tr6+Trb*Tr8)))
      Tra = Tra*(Tr1+Trb*(Tr3+Trb*(Tr5+Trb*Tr7)))
      Fip = Fip + Tra
      Psip = -Fip
      If (Indg==0) Go To 220
      Psi = Psi - Xxx*Dble(Indg)
      Fi = Fi + Xxx*Dble(Indg)
220   Fo = 0.5D0*Dexp(Fi)
      Go = Exp(Psi)
      Fpo = Fo*Fip/Eta2
      Gpo = Go*Psip/Eta2
      If (Jnd==0) Return
      Rau = Ro
      Go = Fo
      Gpo = Fpo
230   X = Rau0 - Ro
      X2 = X*X
      X3 = X*X2
      Unr = 1.D0/Rau0
      Etr0 = 1.D0 - 2.D0*Eta*Unr
      U0 = Go
      U1 = -X*Gpo
      U2 = -0.5D0*Etr0*X2*U0
      S = U0 + U1 + U2
      V1 = U1/X
      V2 = 2.D0*U2/X
      T = V1 + V2
      Xn = 3.D0

      Do N = 3, 10000
!     n=n
         Xn1 = Xn - 1.D0
         Xn1 = Xn*Xn1
         U3 = X*U2*Unr*(1.D0-2.D0/Xn) - Etr0*U1*X2/Xn1 + X3*U0*Unr/Xn1
         S = S + U3
         V3 = Xn*U3/X
         T = T + V3
         If (Abs(U3/S)>1.D-10) Go To 240
         If (Abs(V3/T)<=1.D-10) Go To 250
240      U0 = U1
         U1 = U2
         U2 = U3
         Xn = Xn + 1.D0
      End Do
250   If (Ig) 260, 270, 260
260   Go = S
      Gpo = -T
      Fo = Ho
      Fpo = Hpo
      Return
270   Ho = S
      Hpo = -T
280   Et0 = Eta**(0.166666666666667D0)
      Etad = Etac*Etac
      Et = Eta**(0.6666666666666667D0)
      Et1 = Et*Et
      Et2 = Et1*Et1
      Et3 = Et2*Et
      Et4 = Etad*Et
      Et5 = Et4*Et
      Fo = 1.D0 - Q(1)/Et1 - Q(2)/Etac - Q(3)/Et3 - Q(4)/Etad - Q(5)/Et5
      Go = 1.D0 + Q(1)/Et1 - Q(2)/Etac + Q(3)/Et3 - Q(4)/Etad + Q(5)/Et5
      Fpo = 1.D0 + Qp(1)/Et + Qp(2)/Etac + Qp(3)/Et2 + Qp(4)/Etad + Qp(5)/Et4
      Gpo = 1.D0 - Qp(1)/Et + Qp(2)/Etac - Qp(3)/Et2 + Qp(4)/Etad - Qp(5)/Et4
      Fo = 0.7063326373D0*Et0*Fo
      Go = 1.223404016D0*Et0*Go
      Fpo = 0.4086957323D0*Fpo/Et0
      Gpo = -0.7078817734D0*Gpo/Et0
      Idiv = 0
      If (Ind==0) Return
      Ig = 1
      Rau0 = Eta2
      Go To 230
290   X = Eta2/Rau
      X2 = X*X
      U = 1.D0 - X
      Ru = Sqrt(U)
      U3 = U*U*U
      Trd = 1.D0/(U3*Eta2*Eta2)
      Trc = X2*Trd
      Tre = 1.D0/(U*Ru*Eta2)
      Fi = -0.25D0*Log(U)
      Trb = Trd/64.D0
      Tr3 = (((3.D0*U-4.D0)*U-6.D0)*U+12.D0)*U - 5.D0
      Tr5 = ((((((((-27.D0*U+72.D0)*U-68.D0)*U-360.D0)*U+2190.D0)*U-4808.D0)*U+5148.D0)*U-2712.D0)*U+565.D0)/32.D0
      Tr7 = ((((((((((((81.D0*U-324.D0)*U+486.D0)*U-404.D0)*U-4509.D0)*U+52344.D0)*U-233436.D0)*U+567864.D0)*U &
        -838521.D0)*U+775884.D0)*U-441450.D0)*U+141660.D0)*U-19675.D0)/96.D0
      Fi = Fi + Trb*(Tr3+Trd*(Tr5+Trd*Tr7))

      Fip = 0.25D0/U
      Trb = 3.D0*Trc/(64.D0*U)
      Tr3 = (X-4.D0)*X + 8.D0
      Tr5 = ((((9.D0*X-60.D0)*X+168.D0)*X-192.D0)*X+640.D0)/16.D0
      Tr7 = ((((((-27.D0*X+252.D0)*X-1044.D0)*X+2520.D0)*X-4416.D0)*X-3520.D0)*X-13440.D0)/32.D0
      Fip = Fip + Trb*(Tr3+Trc*(Tr5+Trc*Tr7))
      Tra = Abs((Ru-1.D0)/(Ru+1.D0))
      Psi = (0.5D0*Log(Tra)+Ru/X)*Eta2 + 0.785398163397448D0
      Tr2 = -((9.D0*U-6.D0)*U+5.D0)/48.D0
      Tr4 = ((((((U-2.D0)*945.D0*U+1395.D0)*U-12300.D0)*U+25191.D0)*U-19890.D0)*U+5525.D0)/46080.D0
      Tr6 = (((((((((-G61*U+G62)*U-G63)*U+G64)*U-G65)*U+G66)*U-G67)*U+G68)*U-G69)*U+G610)*U - G611
      Tr8 = (((((((((((((G81*U-G82)*U+G83)*U-G84)*U+G85)*U-G86)*U+G87)*U-G88)*U+G89)*U-G810)*U+G811)*U-G812)*U+G813) &
        *U-G814)*U + G815
      Psi = Psi + Tre*(Tr2+Trd*(Tr4+Trd*(Tr6+Trd*Tr8)))
      Psip = -Ru*Eta2/X2
      Trb = Tre*X/U
      Tr2 = (3.D0*X-8.D0)/32.D0
      Tr4 = -(((63.D0*X-336.D0)*X+704.D0)*X-1536.D0)/2048.D0
      Tr6 = ((((Gp61*X-Gp62)*X+Gp63)*X-Gp64)*X+Gp65)*X - Gp66
      Tr8 = ((((((-Gp81*X+Gp82)*X-Gp83)*X+Gp84)*X-Gp85)*X+Gp86)*X+Gp87)*X + Gp88
      Psip = Psip + Trb*(Tr2+Trc*(Tr4+Trc*(Tr6+Trc*Tr8)))
      Tra = Exp(Fi)
      Fo = Tra*Sin(Psi)
      Go = Tra*Cos(Psi)
      If (Indice) 300, 310, 300
300   Tra = Fo
      Fo = Go
      Go = -Tra
310   Tra = -Eta2/Rauc
      Fpo = (Fip*Fo+Psip*Go)*Tra
      Gpo = (Fip*Go-Psip*Fo)*Tra
      Return
   End Subroutine Yfrica

!***********************************************************************
   Subroutine Jfdelg(Xd, Yd, Par, Pai, Nbchif)
!***********************************************************************
      Implicit None
      Integer, Intent (Out) :: Nbchif
      Integer :: Kr, I, K
      Double Precision, Intent (In) :: Xd, Yd
      Double Precision, Intent (Out) :: Par, Pai
      Double Precision, Parameter :: Rac2 = 0.3465735902799726D0
      Double Precision, Parameter :: Pis4 = 0.785398163397448D0
      Double Precision, Parameter :: Pi = 3.141592653589793D0
      Double Precision, Parameter :: Depi = 6.283185307179586D0
      Double Precision, Parameter :: Supint = 2147483647.0D0
      Double Precision :: Test(7) = (/ 2.9152D7, 2.2958D3, 1.4124D2, 3.9522D1, 19.6611D0, 12.791D0, -10.D0 /)
      Double Precision :: C(6) = (/ 8.333333333333333D-2, -8.33333333333333D-3, 3.968253968253968D-3, &
         -4.166666666666667D-3, 7.575757575757576D-3, -2.109279609279609D-2 /)
      Double Precision :: X, Y, U, V, Tra, Tra1, Trb, Cosi, Cos2i, Sini
      Double Precision :: Sin2i, Zmod, Xx

      X = Abs(Xd)
      Xx = X
      Nbchif = 15
      If (Yd) 100, 460, 100
100   Y = Abs(Yd)
      Kr = 1
      I = Mod(10.99D0-X, Supint)
!     translation
      If (I) 120, 120, 110
110   Tra = I
      X = X + Tra
!     logarithme(x+iy) (x,y positifs)
120   If (X-Y) 130, 140, 150
130   Tra1 = X/Y
      Trb = 1.D0 + Tra1*Tra1
      Tra = Y*Sqrt(Trb)
      Sini = 1./(Trb*Y)
      Cosi = Sini*Tra1
      Tra1 = Y/X
      Go To 160
140   U = Rac2 + Log(X)
      V = Pis4
      Sini = 0.5D0/X
      Cosi = Sini
      Go To 170
150   Tra1 = Y/X
      Trb = 1.D0 + Tra1*Tra1
      Tra = X*Sqrt(Trb)
      Cosi = 1./(Trb*X)
      Sini = Cosi*Tra1
160   U = Log(Tra)
      V = Atan(Tra1)
!     developpement asymptotique ( x superieur a 10 )
170   Par = U - 0.5*Cosi
      Pai = V + 0.5*Sini
      Zmod = X + Y
      If (Zmod-Test(1)) 180, 180, 250
180   Sin2i = (Sini*Cosi) + (Sini*Cosi)
      Cos2i = (Cosi+Sini)*(Cosi-Sini)
      Sini = Sin2i
      Cosi = Cos2i
      K = 1
      Go To 200
190   Tra = Cosi*Cos2i - Sini*Sin2i
      Sini = Sini*Cos2i + Cosi*Sin2i
      Cosi = Tra
200   Par = Par - C(K)*Cosi
      Pai = Pai + C(K)*Sini
      K = K + 1
      If (Zmod-Test(K)) 190, 190, 250
!     translation inverse
210   I = I - 1
      X = I
      X = Xx + X
      If (X-Y) 220, 220, 230
220   Tra1 = X/Y
      Trb = X*Tra1 + Y
      Sini = 1.D0/Trb
      Cosi = Tra1/Trb
      Go To 240
230   Tra1 = Y/X
      Trb = X + Y*Tra1
      Cosi = 1.D0/Trb
      Sini = Tra1/Trb
240   Par = Par - Cosi
      Pai = Pai + Sini
250   If (I) 260, 260, 210

!     controle du quadrant
260   If (Xd) 320, 270, 310
270   Tra = Pi*Y
      If (Tra-1.D-2) 280, 280, 290
280   Trb = Tra*(2.D0+Tra*(-2.D0+Tra*(1.333333333333333D0+Tra*(-0.6666666666666666D0+Tra*(0.2666666666666666D0 &
               +Tra*(-0.08888888888888888D0+Tra*0.02539682539682540D0))))))
      Trb = (2.D0-Trb)/Trb
      Go To 300
290   Trb = Exp(-Tra-Tra)
      Trb = (1.D0+Trb)/(1.D0-Trb)
300   Pai = 0.5D0*(1.D0/Y+Pi*Trb)
310   If (Yd) 330, 400, 400
!     x+iy change en -x-iy
320   Tra = Exp(-Depi*Y)
      Trb = Tra*Tra
      Cos2i = Depi*Mod(X, 1.D0)
      Sin2i = -2.D0*Tra*Cos(Cos2i) + 1.D0 + Trb
      Par = Par + Cosi + Depi*Tra*Sin(Cos2i)/Sin2i
      Pai = Pai - Sini + Pi*(Trb-1.D0)/Sin2i
      If (Yd) 400, 400, 330
330   Pai = -Pai
340   Write (6, 350)
350   Format (' argument de jfdelg trop grand')
      Go To 510
360   Y = Mod(X, Supint)
      If (Y) 370, 490, 370
370   If (Y-0.99D0) 520, 520, 380
380   Tra = Int(Y+0.1D0)
      If (Abs(Y-Tra)-5.D-15) 490, 490, 520
390   If (X-4503599627370496.D0) 360, 340, 340

!     argument dans -pi,pi
400   Tra = Abs(Pai/Depi)
      If (Tra-1.D+15) 420, 410, 410
410   Nbchif = 0
      Pai = 0.D0
      Go To 620
420   If (Tra-1.D0) 440, 440, 430
430   Nbchif = Log10(Tra)
      Nbchif = 14 - Nbchif
      Tra = Mod(Tra, Supint)
      Pai = Mod(Tra, 1.D0)*Sign(Depi, Pai)
440   If (Dabs(Pai)-Pi) 620, 620, 450
450   Pai = Pai - Sign(Depi, Pai)
      Go To 620
!        delgamma reel
460   Pai = 0.D0
      If (Xd) 390, 470, 520
!     conditions d existence
470   Write (6, 480)
480   Format (' jfdelg(0) est infini')
      Go To 510
490   Write (6, 500)
500   Format (' jfdelg (-entier) est infini')
510   Par = 1.D+74
      Nbchif = 0
      Go To 620
!     translation
520   I = Mod(10.99D0-X, Supint)
      If (I) 540, 540, 530
530   Tra = I
      X = X + Tra
!     developpement asymptotique
540   Y = Log(X)
      Par = Y - 0.5D0/X
      If (X-Test(1)) 550, 550, 600
550   Cos2i = 1.D0/(X*X)
      Cosi = Cos2i
      K = 1
      Go To 570
560   Cosi = Cosi*Cos2i
570   Par = Par - C(K)*Cosi
      K = K + 1
      If (X-Test(K)) 560, 560, 590
!     translation inverse
580   I = I - 1
      X = I
      X = Xx + X
      Par = Par - 1.D0/X
590   If (I) 600, 600, 580
!     x negatif
600   If (Xd) 610, 620, 620
610   Par = Par + 1.D0/X
      Y = Pi*Mod(X, 2.D0)
      Par = Par + Pi*Cos(Y)/Sin(Y)
!entry jfdev1
620   Return
   End Subroutine Jfdelg
!======================================================================!
End Module mdl_006_cwf
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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
! modified numerov method.
! Note; s.p.W.F. is not saved.
!**************************************************************
      use mdl_006_cwf
      use mdl_001_setting, only: hbarc, xmu, zc, tol, alpha, ai, dr, rmax, nrmax
!!      use mdl_003_subroutines, only: vcx
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

!--- for Coulomb WF's
      fcw = 0.d0
      gcw = 0.d0
      fpcw = 0.d0
      gpcw = 0.d0
      sigmad = 0.d0
      iexp = 0

!--- initial conditions
      psi0 = 0.d0 !at r=0
      psi1 = 1.d-2 !at r=1*dr
!      psi1 = dr**(l+1) !at r=1*dr

      fac = 2.d0*xmu*dr*dr/hbarc/hbarc
      xi0 = 0.d0
      xi1 = (1.d0-fac*(vcx(xind,l,j,rmind+dr)-ein)/12.d0)*psi1
!      xi1 = ( 1.d0 - fac*(v(rmind+dr) - ai*w(rmind+dr) - ein)/12.d0 )*psi1

!--- start of the iterations for s.p.W.F.
      do ir = 2, iterat + 1
         r = rmind + dr*dble(ir)
         ri = rmind + dr*dble(ir-2)
         r1 = rmind + dr*dble(ir-1)

         dd = fac*(vcx(xind,l,j,r1)-ein)/sqrt(12.d0) + sqrt(3.d0)
!         dd = fac*(v(r1) - ai*w(r1) - ein)/sqrt(12.d0) + sqrt(3.d0)

         xi = (dd*dd-1.d0)*xi1 - xi0

         if (ir==iterat+1) go to 100

         xi0 = xi1/xi
         xi1 = xi/xi
100   end do

      cc = -fac*(vcx(xind,l,j,rmaxd-dr)-ein)/12.d0 + 1.d0
!      cc = -fac*(v(rmaxd-dr) - ai*w(rmaxd-dr) - ein)/12.d0 + 1.d0
      psi0 = xi0/cc

      cc = -fac*(vcx(xind,l,j,rmaxd+dr)-ein)/12.d0 + 1.d0
!      cc = -fac/12.d0*(v(rmaxd+dr)-ai*w(rmaxd+dr)-ein) + 1.d0
      psi = xi/cc

!--- compute coulomb W.F.
      ak = sqrt(2.d0*xmu*ein)/hbarc
      eta = alpha*zc*sqrt(xmu/2.d0/ein) != (zp*zt/137.d0)*sqrt(xmu/2.d0/ein)
      IF(xind=='neut')  eta = 0.d0

      rho = (rmaxd-dr)*sqrt(2.d0*xmu*ein)/hbarc
      call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, l, iexp)
      cwup0 = gcw(l) + ai*fcw(l)
      cwdown0 = gcw(l) - ai*fcw(l)

      rho = (rmaxd+dr)*sqrt(2.d0*xmu*ein)/hbarc
      call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, l, iexp)
      cwup1 = gcw(l) + ai*fcw(l)
      cwdown1 = gcw(l) - ai*fcw(l)

!--- matching at rmaxd with Coulomb W.F.
      bb = (cwup0*psi-cwup1*psi0)/(cwup0*cwdown1-cwup1*cwdown0)
      bb2 = (cwdown0*psi-cwdown1*psi0)/(cwup0*cwdown1-cwup1*cwdown0)

!--- results
      zs = bb2/bb !S-matrix
      zf = dble(2*l+1)*exp(2.d0*ai*sigmad(l))*(zs-1.d0)/2.d0/ai/ak
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
   use mdl_001_setting
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

!--- s(1/2)-channel:
      call coulscat(xind, 0, 1, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_s1(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC, s(1/2)".'
      endif

!--- p(3/2)-channel:
      call coulscat(xind,1, 3, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_p3(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC, p(3/2)".'
      endif

!--- d(5/2)-channel:
      call coulscat(xind, 2, 5, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_d5(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC, d(5/2)".'
      endif

!--- d(3/2)-channel:
      call coulscat(xind, 2, 3, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_d3(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC".'
      endif

!--- f(7/2)-channel:
      call coulscat(xind, 3, 7, e, smat, facn)
      z1 = 0.5d0*log(smat)/ai !p.shift
      if (real(z1)<0.d0) z1 = z1 + pi
      delta_f7(ie) = real(z1)
      if( abs(1.d0-abs(smat)) .gt. tol) then
         write(6,*) 'abs(smat) is not 1, check "CN_SC,".'
      endif

!--- f(5/2)-channel:
      call coulscat(xind, 3, 5, e, smat, facn)
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
      call CN_SC('prot')
                    l = 9999 ; jj = -9999!---------------------------------filter
      call CN_TD('prot',l,jj)
   end if
   !----------------------------
   if(2.eq.2) then
   write(6,*) "**************************************************"
   write(6,*) "************* NEUTRON ****************************"
   write(6,*) "**************************************************"
      call CN_SC('neut')
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

