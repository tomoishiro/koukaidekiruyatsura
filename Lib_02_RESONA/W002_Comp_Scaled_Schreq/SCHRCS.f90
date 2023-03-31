!///////////////////////////////////////////////////////////////////////////////
module mdl_003_energies
   double precision, parameter :: pi = 3.141592653d0
   double precision, parameter :: theta = (pi*0.5d0)*0.123d0 !comple--scaling angle
   !---
   double precision, parameter :: recent =  0.4314d0!MeV
   double precision, parameter :: gapr   =  0.00018d0
   double precision, parameter ::   de   =  0.00001d0
   integer, parameter :: nemax = int(gapr/de +1.d-8)*2 +1
   double precision, parameter :: emin = recent -gapr
   double precision, parameter :: emax = recent +gapr
   !---
   double precision, parameter :: eicent  = -0.00764d0!MeV
   double precision, parameter :: gapi    =  0.00018d0
   double precision, parameter ::   dg    =  0.00001d0
   integer, parameter :: ngmax = int(gapi/dg +1.d-8)*2 +1
   double precision, parameter :: gmin = eicent -gapi
   double precision, parameter :: gmax = eicent +gapi
end module mdl_003_energies
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
module mdl_001_setting
   double precision, parameter :: pi = 3.141592653d0
   double precision, parameter :: hbarc = 197.329d0 !MeV.fm
   double precision, parameter :: alpha = 7.2973525698d-3 !fine structure constant
   double precision, parameter :: tol = 1.d-12 !tolerance
   complex(kind=kind(0d0)), parameter :: ai = (0.d0, 1.d0)
   !-------------------------------
   double precision, parameter :: nc = 8.d0
   double precision, parameter :: zc = 8.d0 !proton-number of core
   double precision, parameter :: ac = nc+zc !mass-number of core
   !-------------------------------
   double precision, parameter :: anmass_p =  938.2720813d0 ! proton-mass
   double precision, parameter :: anmass_n =  939.5654133d0 !neutron-mass
   double precision, parameter :: BE_core = 7.976207d0*ac !MeV, BE for O-16 nucleus.
!!   double precision, parameter :: BE_core = 8.5513064d0*ac !MeV, BE for Ca-40 nucleus.
!!   double precision, parameter :: BE_core = 8.253d0*ac !MeV, BE for Sn-100 nucleus.
   double precision, parameter :: mass_c = zc*anmass_p + nc*anmass_n - BE_core !MeV
   double precision, parameter :: xmu_cp = anmass_p*mass_c/(anmass_p + mass_c)
   double precision, parameter :: xmu_cn = anmass_n*mass_c/(anmass_n + mass_c)
!<>   double precision, parameter :: amcc = anmass_p !MeV
!!   double precision, parameter :: amcc = xmu_cp !MeV
   double precision, parameter :: amcc = xmu_cn !MeV
   !-------------------------------
!<>   double precision, parameter :: h2omu = 1.d0 !MeV.fm^2, =hbar*hbar/mu.
   double precision, parameter :: h2omu = hbarc*hbarc/amcc !MeV.fm^2, =hbar*hbar/mu.
   !-------------------------------
   double precision, parameter :: Rmax = 14.d0 !fm
   double precision, parameter :: dr = 0.05d0 !fm
   integer, parameter :: Nr = int(Rmax/dr + tol)
   integer, parameter :: nspmax = 100
   integer, parameter :: leo = 2
   integer, parameter :: jeo = 3
   !-------------------------------
end module mdl_001_setting
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
module mdl_004_potentials
   use mdl_001_setting, only: ac
   implicit none
   !--- parameters of the core-n (-p) WS(+Coulomb) potential.
   complex(kind=kind(0d0)), parameter :: rc0 = 1.25d0*ac**(1.d0/3.d0)
   complex(kind=kind(0d0)), parameter :: ac0 = 0.65d0
!<40Ca>   complex(kind=kind(0d0)), parameter :: W_0  = -55.57d0 ![MeV]
!<40Ca>   complex(kind=kind(0d0)), parameter :: W_ls =  11.28d0 ![MeV*fm^2]
   complex(kind=kind(0d0)), parameter :: W_0  = -53.2d0 ![MeV]
   complex(kind=kind(0d0)), parameter :: W_ls =  22.1d0 ![MeV*fm^2]
contains
!*****************************************************
function zwsp(zr) result(zf) !--- My Potential
  use mdl_001_setting, only : leo, jeo
  implicit none
  complex(kind=kind(0d0)), intent(IN) :: zr
  complex(kind=kind(0d0)) :: zf
  complex(kind=kind(0d0)) :: zc, z2, zls, zx
  double precision :: yl, yj, xls
  zf = 0.d0
  zx = exp((zr - rc0) / ac0)
  !--- central Woods-Saxon
  zc = W_0 / (1.d0 + zx)
  !--- LS Woods-Saxon
!!  z2  = -exp((zr - rc0)/ac0) / (1.d0 + exp((zr - rc0)/ac0))**2/ac0
  z2  = (1.d0 + zx)**2
  z2  = 1.d0/z2
  z2  = z2 * (-zx) / ac0
  yl = dble(leo) ; yj = dble(jeo)/2.d0
  xls = 0.5d0*(yj*(yj+1.d0)-yl*(yl+1.d0)-0.75d0)
  zls = W_ls * xls * z2 /(zr +1.d-12)
  !--- Note: centrigugal is added in the next routine "c21".
  zf = zc + zls
end function zwsp
!*****************************************************
function zvho(zr) result(zf) !--- HO potential for testing.
  use mdl_001_setting, only : h2omu
  implicit none
  complex(kind=kind(0d0)), intent(IN) :: zr
  complex(kind=kind(0d0)) :: zf
  double precision :: hbrom, V0, p
  hbrom = 0.5d0 !MeV
  p = 0.5d0*hbrom*hbrom/h2omu
  V0 = -4.d0 !MeV
  zf = cmplx(p)*zr*zr  +V0 !-dble(l*(l+1))/(zr*zr +1.d-9)
end function zvho
!*****************************************************
function zell(zr) result(zf) !--- Myo's Potential
  implicit none
  complex(kind=kind(0d0)), intent(IN) :: zr
  complex(kind=kind(0d0)) :: v1, v2, b1, b2
  complex(kind=kind(0d0)) :: zf
  v1 = -8.d0
  v2 =  4.d0
  b1 = -0.16d0
  b2 = -0.04d0
  zf = v1*exp(b1*zr*zr) +v2*exp(b2*zr*zr)
end function zell
!*****************************************************
end module mdl_004_potentials
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
module mdl_002_routines
contains
!*****************************************************
function zfre(zr) result(zf)
  use mdl_001_setting, only : h2omu
  use mdl_004_potentials
  implicit none
  complex(kind=kind(0d0)), intent(IN) :: zr
  complex(kind=kind(0d0)) :: zf, y
  y = zr
!!  zf = zvho(y)
!!  zf = zell(y)
  zf = zwsp(y)
end function zfre
!*****************************************************
function a11(zr) result(zf)
  implicit none
  complex(kind=kind(0d0)), intent(IN) :: zr
  complex(kind=kind(0d0)) :: zf
  zf = (0.d0, 0.d0)
end function a11
!*****************************************************
function b12(zr) result(zf)
  implicit none
  complex(kind=kind(0d0)), intent(IN) :: zr
  complex(kind=kind(0d0)) :: zf
  zf = (1.d0, 0.d0)
end function b12
!*****************************************************
function c21(zr, ein) result(zf)
  use mdl_001_setting, only : h2omu, leo
  implicit none
  complex(kind=kind(0d0)), intent(IN) :: zr, ein
  complex(kind=kind(0d0)) :: zf, z1
  z1 = zr
  zf = -2.d0*(ein -zfre(z1)) / h2omu  +dble(leo*(leo+1))/(zr*zr +1.d-9) !MeV
end function c21
!*****************************************************
function d22(zr) result(zf)
  implicit none
  complex(kind=kind(0d0)), intent(IN) :: zr
  complex(kind=kind(0d0)) :: zf
  zf = (0.d0, 0.d0)
  end function d22
!*****************************************************
!*****************************************************
!*****************************************************
! Subroutine to solve the 2*2 matrix 1st-derivative equation,
!   d/dr [A(r)] = {a11(r)   b12(r)} [A(r)]
!        [B(r)] = {c21(r)   d22(r)} [B(r)],
! by the Runge-Kutta-4 method.
! Note that A(r), B(r), a11(r), etc. are assumed as complex.
subroutine RK4_forward(N, A, B, ZR, D, E)
  implicit none
  integer, intent(IN) :: N !number of mesh points r=0~rmax
  complex(kind=kind(0d0)), intent(INOUT) :: A(0:N), B(0:N)
  complex(kind=kind(0d0)), intent(IN) :: ZR(0:N) !mesh points for A & B
  complex(kind=kind(0d0)), intent(IN) :: D !mesh gap of ZR(N)
  complex(kind=kind(0d0)), intent(IN) :: E !energy as complex variable
  complex(kind=kind(0d0)) :: za1, za2, za3, za4
  complex(kind=kind(0d0)) :: zb1, zb2, zb3, zb4
  complex(kind=kind(0d0)) :: z, wa, wb, zd
  integer :: i
  zd = D
  do i = 2, N-10!N
     !---RK1
     wa = A(i-1)
     wb = B(i-1)
     z = ZR(i-1)+            zd*0.d0
     za1 = a11(z)*wa   +b12(z)*wb
     zb1 = c21(z,E)*wa +d22(z)*wb
     !---RK2
     wa = A(i-1)  + za1*0.5d0*zd
     wb = B(i-1)  + zb1*0.5d0*zd
     z = ZR(i-1)+           zd*0.5d0
     za2 = a11(z)*wa   + b12(z)*wb
     zb2 = c21(z,E)*wa + d22(z)*wb
     !---RK3
     wa = A(i-1) + za2*0.5d0*zd
     wb = B(i-1) + zb2*0.5d0*zd
     z = ZR(i-1)+           zd*0.5d0
     za3 = a11(z)*wa   +b12(z)*wb
     zb3 = c21(z,E)*wa +d22(z)*wb
     !---RK3
     wa = A(i-1) + za3*zd
     wb = B(i-1) + zb3*zd
     z = ZR(i-1)+           zd*1.d0
     za4 = a11(z)*wa   + b12(z)*wb
     zb4 = c21(z,E)*wa + d22(z)*wb
     !---Total
     A(i) = A(i-1) +zd*(za1 +2.d0*za2 +2.d0*za3 +za4)/6.d0
     B(i) = B(i-1) +zd*(zb1 +2.d0*zb2 +2.d0*zb3 +zb4)/6.d0
  end do
!  z = 0.d0
!  do i = 0, N!--- rough normalization
!     z = z + conjg(A(i))*A(i)*D
!  end do
!  A = A/sqrt(abs(z))
!  B = B/sqrt(abs(z))
  return
end subroutine RK4_forward
!*****************************************************
subroutine RK4_backward(N, A, B, ZR, D, E)
  implicit none
  integer, intent(IN) :: N !number of mesh points r=0~rmax
  complex(kind=kind(0d0)), intent(INOUT) :: A(0:N), B(0:N)
  complex(kind=kind(0d0)), intent(IN) :: ZR(0:N) !mesh points for A & B
  complex(kind=kind(0d0)), intent(IN) :: D !mesh gap of ZR(N)
  complex(kind=kind(0d0)), intent(IN) :: E !energy as complex variable
  complex(kind=kind(0d0)) :: za1, za2, za3, za4
  complex(kind=kind(0d0)) :: zb1, zb2, zb3, zb4
  complex(kind=kind(0d0)) :: z, wa, wb, zd
  integer :: i
  zd = -D
  do i = N-1, 7, -1
     !---RK1
     z = ZR(i+1) +zd*0.d0
     wa = A(i+1)
     wb = B(i+1)
     za1 = a11(z)*wa  +b12(z)*wb
     zb1 = c21(z,E)*wa+d22(z)*wb
     !---RK2
     z = ZR(i+1) +zd*0.5d0
     wa = A(i+1) +za1*0.5d0*zd
     wb = B(i+1) +zb1*0.5d0*zd
     za2 = a11(z)*wa  +b12(z)*wb
     zb2 = c21(z,E)*wa+d22(z)*wb
     !---RK3
     z = ZR(i+1) +zd*0.5d0
     wa = A(i+1) +za2*0.5d0*zd
     wb = B(i+1) +zb2*0.5d0*zd
     za3 = a11(z)*wa  +b12(z)*wb
     zb3 = c21(z,E)*wa+d22(z)*wb
     !---RK3
     z = ZR(i+1) +zd
     wa = A(i+1) +za3*zd
     wb = B(i+1) +zb3*zd
     za4 = a11(z)*wa  +b12(z)*wb
     zb4 = c21(z,E)*wa+d22(z)*wb
     !---Total
     A(i) = A(i+1) +zd*(za1 +2.d0*za2 +2.d0*za3 +za4)/6.d0
     B(i) = B(i+1) +zd*(zb1 +2.d0*zb2 +2.d0*zb3 +zb4)/6.d0
  end do
!  z = 0.d0
!  do i = 0, N!--- rough normalization
!     z = z + conjg(A(i))*A(i)*D
!  end do
!  A = A/sqrt(abs(z))
!  B = B/sqrt(abs(z))
  return
end subroutine RK4_backward
!*****************************************************

end module mdl_002_routines
!///////////////////////////////////////////////////////////////////////////////


subroutine RK4_FG(zein, AF,BF,AG,BG, rm)
  use mdl_001_setting, only : Nr,dr,rmax, tol, pi,ai, h2omu, leo
  use mdl_002_routines, only : zfre, RK4_forward, RK4_backward
  use mdl_003_energies, only : theta
  implicit none
  double precision, intent(IN) :: rm !fm, matching point
  complex(kind=kind(0d0)), intent(INOUT) :: zein
  complex(kind=kind(0d0)), intent(OUT) :: AF,BF,AG,BG
  !---
!<>  double precision, parameter :: theta = 0.d0 !pi/16.d0 !comple--scaling angle
  double precision :: x1,x2
  complex(kind=kind(0d0)) :: zr,z1,z2,z3, csf, csg, zkap
  complex(kind=kind(0d0)) :: zrc(0:Nr), UR(0:Nr), US(0:Nr)
  complex(kind=kind(0d0)) :: VR(0:Nr), VS(0:Nr)
  integer :: i, j, match
  match = int(rm/dr +tol)
  !---Wave number, being possibly complex
  zkap = sqrt(2.d0*(-zein)/h2omu)
  if (real(zein) .ge. 0.d0) then
     zkap = (2.d0*(7.77777d0)/h2omu)
  end if
  !---Complex-scaling factor
  csf = exp( ai*theta)
  csg = exp(-ai*theta)
  !---Potential of interest
  open(11, file='W100_Poten.dat', action='write')
  do i=0,Nr
     x1 = dble(i)*dr
     z1 = x1*csf
     z2 = zfre(z1)
     write(11, *) x1, real(z2), aimag(z2)
  end do
  close(11)

  !---Coordinate points, being possibly complex
  do i = 0, Nr
     x1 = dble(i)*dr
     zrc(i) = x1*csf
  end do

  !---Forward Runge-Kutta from nearly r=0
  z1=cmplx(dr)*csf
  UR(0) = 0.d0!zr**(leo+1) !wave-function but multiplied by r, and
  US(0) = 0.d0!dble(leo+1)*zr**(leo) !its derivative in the asymptotic form
  UR(1) = z1**(leo+1) !wave-function but multiplied by r, and
  US(1) = dble(leo+1)*z1**(leo) !its derivative in the asymptotic form
  j=Nr
  call RK4_forward(j,UR,US,zrc,z1,zein)
  !---Convert from U(r) to the U(r)/r:
  do i=1,Nr
     z1 = zrc(i)
     z2 = UR(i)/z1
     z3 = US(i)/z1 -UR(i)/(z1*z1)
     UR(i) = z2
     US(i) = z3
  end do

  !---Backward Runge-Kutta
  z1=cmplx(dr)*csf
  z2=zrc(Nr)
  VR(Nr) = z2*exp(-z2*zkap) !wave-function but multiplied by r
  VS(Nr) = (1.d0 -z2*zkap)*exp(-z2*zkap) !derivative
  j=Nr
  call RK4_backward(j,VR,VS,zrc,z1,zein)
  !---Convert from U(r) to the U(r)/r:
  do i=Nr,0,-1
     z1 = zrc(i)
     z2 = VR(i)/z1
     z3 = VS(i)/z1 -VR(i)/(z1*z1)
     VR(i) = z2
     VS(i) = z3
  end do
  
  !---At matching point,
  z1=UR(match); UR=UR/z1; US=US/z1
  z2=VR(match); VR=VR/z2; VS=VS/z2
  AF = UR(match)!wave-function by RK4_FORWARD
  BF = US(match)!derivative by RK4_FORWARD
  AG = VR(match)!wave-function by RK4_BACKWARD
  BG = VS(match)!derivative by RK4_BACKWARD

  !---Plotting
  open(11, file='W201_WF.dat', action='write')
  open(12, file='W202_WFDERI.dat', action='write')
  do i=1,Nr
     z1 = zrc(i); z2 = UR(i); z3 = US(i)
     write(11, *) abs(z1), aimag(z1), real(z2), aimag(z2)
     write(12, *) abs(z1), aimag(z1), real(z3), aimag(z3)
  end do
  close(11)
  close(12)
  open(11, file='W301_WF.dat', action='write')
  open(12, file='W302_WFDERI.dat', action='write')
  do i=Nr,match,-1
     z1 = zrc(i); z2 = VR(i); z3 = VS(i)
     write(11, *) abs(z1), aimag(z1), real(z2), aimag(z2)
     write(12, *) abs(z1), aimag(z1), real(z3), aimag(z3)
  end do
  close(11)
  close(12)
  return
end subroutine RK4_FG

program main
  use mdl_001_setting, only : ai, tol
  use mdl_003_energies, only : emin,nemax,de, gmin,ngmax,dg
  implicit none
  complex(kind=kind(0d0)), dimension(:,:), allocatable :: zeins, gp, gp2
  complex(kind=kind(0d0)) :: AF,BF,AG,BG, ze, z1, werr
  complex(kind=kind(0d0)) :: AF2,BF2,AG2,BG2
  integer :: i,j,ne,ng
  double precision :: x1,x2,x3,x4
  double precision, parameter :: rmatch = 2.6d0 !fm
  double precision, parameter :: rmcorr = 5.2d0 !fm
  allocate(zeins(nemax,ngmax), gp(nemax,ngmax), gp2(nemax,ngmax))
  !---Loop for energies, being possibly complex
  do i = 1, nemax
     z1 = emin + de*dble(i-1)
     do j = 1, ngmax
        ze = z1 + ai*(gmin + dg*dble(j-1))
        zeins(i,j) = ze
        x1 = rmatch
        call RK4_FG(ze, AF,BF, AG,BG, x1)
        werr = BF*AG - AF*BG
        gp(i,j) = werr
        x1 = rmcorr
        call RK4_FG(ze, AF2,BF2, AG2,BG2, x1)
        werr = BF2*AG2 - AF2*BG2
        gp2(i,j) = werr
     end do
  end do
  !---Find the energy with minimum error "gp"
  x4 = 999999999.d0
  open(16, file='W621_ERROR.dat', action='write')
  do i = 1, nemax
     do j = 1, ngmax
        x1 = abs(gp(i,j))
        x2 = abs(gp2(i,j))
        x3 = sqrt(x1*x2)
        ze = zeins(i,j)
        write(16, *) real(ze), aimag(ze), x1,x2,x3
        if(x3 .lt. x4) then
           ne=i
           ng=j
           x4=x3
        end if
     end do
     write(16,*) ""!---for 3D plotting
  end do
  close(16)
  !---Repeat RK4 for the eigen-energy found
  ze=zeins(ne,ng)
  write(6,*) "Re(E) & Im(E) = ", real(ze)*(1.d0+tol), aimag(ze)*(1.d0+tol), " MeV"
  write(6,*) "ERROR = ", gp(ne,ng), gp2(ne,ng)
  x1 = rmatch
  call RK4_FG(ze, AF,BF,AG,BG, x1)
  deallocate(zeins,gp,gp2)
end program main


!//////////////////////////////////////////////////////////////
!NOTE:
! <2022/11/22>
! Forward & Backward RK4 are done.
! Energy plane is done.
! How to correctly evaluate the continuous condition?
!//////////////////////////////////////////////////////////////


