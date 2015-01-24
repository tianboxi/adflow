   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of invisciddissfluxscalarapprox in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: *w *fw
   !   with respect to varying inputs: gammainf rhoinf pinfcorr *p
   !                *w *radi *radj *radk
   !   Plus diff mem management of: p:in w:in fw:in radi:in radj:in
   !                radk:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          inviscidDissFluxScalar.f90                      *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 03-24-2003                                      *
   !      * Last modified: 10-29-2007                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE INVISCIDDISSFLUXSCALARAPPROX_D()
   !
   !      ******************************************************************
   !      *                                                                *
   !      * inviscidDissFluxScalar computes the scalar artificial          *
   !      * dissipation, see AIAA paper 81-1259, for a given block.        *
   !      * Therefore it is assumed that the pointers in  blockPointers    *
   !      * already point to the correct block.                            *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS
   USE CGNSGRID
   USE CONSTANTS
   USE FLOWVARREFSTATE
   USE INPUTDISCRETIZATION
   USE INPUTPHYSICS
   USE ITERATION
   IMPLICIT NONE
   !
   !      Local parameter.
   !
   REAL(kind=realtype), PARAMETER :: dssmax=0.25_realType
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, k, ind
   REAL(kind=realtype) :: sslim, rhoi
   REAL(kind=realtype) :: sslimd, rhoid
   REAL(kind=realtype) :: sfil, fis2, fis4
   REAL(kind=realtype) :: ppor, rrad, dis2
   REAL(kind=realtype) :: rradd, dis2d
   REAL(kind=realtype) :: dss1, dss2, ddw, fs
   REAL(kind=realtype) :: dss1d, dss2d, ddwd, fsd
   INTRINSIC ABS
   INTRINSIC MAX
   INTRINSIC MIN
   REAL(kind=realtype) :: pwr1
   REAL(kind=realtype) :: pwr1d
   REAL(kind=realtype) :: x6d
   REAL(kind=realtype) :: min3
   REAL(kind=realtype) :: min2
   REAL(kind=realtype) :: min1
   REAL(kind=realtype) :: x6
   REAL(kind=realtype) :: x5
   REAL(kind=realtype) :: min1d
   REAL(kind=realtype) :: x4
   REAL(kind=realtype) :: x3
   REAL(kind=realtype) :: x2
   REAL(kind=realtype) :: x2d
   REAL(kind=realtype) :: x1
   REAL(kind=realtype) :: x5d
   REAL(kind=realtype) :: y3d
   REAL(kind=realtype) :: x1d
   REAL(kind=realtype) :: min3d
   REAL(kind=realtype) :: x4d
   REAL(kind=realtype) :: y2d
   REAL(kind=realtype) :: abs0
   REAL(kind=realtype) :: min2d
   REAL(kind=realtype) :: y3
   REAL(kind=realtype) :: y2
   REAL(kind=realtype) :: x3d
   REAL(kind=realtype) :: y1
   REAL(kind=realtype) :: y1d
   IF (rfil .GE. 0.) THEN
   abs0 = rfil
   ELSE
   abs0 = -rfil
   END IF
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Check if rFil == 0. If so, the dissipative flux needs not to
   ! be computed.
   IF (abs0 .LT. thresholdreal) THEN
   fwd = 0.0_8
   RETURN
   ELSE
   ! Determine the variables used to compute the switch.
   ! For the inviscid case this is the pressure; for the viscous
   ! case it is the entropy.
   SELECT CASE  (equations) 
   CASE (eulerequations) 
   ! Inviscid case. Pressure switch is based on the pressure.
   ! Also set the value of sslim. To be fully consistent this
   ! must have the dimension of pressure and it is therefore
   ! set to a fraction of the free stream value.
   sslimd = 0.001_realType*pinfcorrd
   sslim = 0.001_realType*pinfcorr
   CASE (nsequations, ransequations) 
   !===============================================================
   ! Viscous case. Pressure switch is based on the entropy.
   ! Also set the value of sslim. To be fully consistent this
   ! must have the dimension of entropy and it is therefore
   ! set to a fraction of the free stream value.
   IF (rhoinf .GT. 0.0_8) THEN
   pwr1d = rhoinf**gammainf*(LOG(rhoinf)*gammainfd+gammainf*rhoinfd&
   &         /rhoinf)
   ELSE IF (rhoinf .EQ. 0.0_8) THEN
   IF (gammainf .EQ. 1.0) THEN
   pwr1d = rhoinfd
   ELSE
   pwr1d = 0.0_8
   END IF
   ELSE IF (gammainf .EQ. INT(gammainf)) THEN
   pwr1d = gammainf*rhoinf**(gammainf-1)*rhoinfd
   ELSE
   pwr1d = 0.0_8
   END IF
   pwr1 = rhoinf**gammainf
   sslimd = (0.001_realType*pinfcorrd*pwr1-0.001_realType*pinfcorr*&
   &       pwr1d)/pwr1**2
   sslim = 0.001_realType*pinfcorr/pwr1
   CASE DEFAULT
   sslimd = 0.0_8
   END SELECT
   ! Set a couple of constants for the scheme.
   fis2 = rfil*vis2
   fis4 = rfil*vis4
   sfil = one - rfil
   ! Replace the total energy by rho times the total enthalpy.
   ! In this way the numerical solution is total enthalpy preserving
   ! for the steady Euler equations. Also replace the velocities by
   ! the momentum. Only done for the entries used in the
   ! discretization, i.e. ignore the corner halo's.
   DO k=0,kb
   DO j=2,jl
   DO i=2,il
   wd(i, j, k, ivx) = wd(i, j, k, irho)*w(i, j, k, ivx) + w(i, j&
   &           , k, irho)*wd(i, j, k, ivx)
   w(i, j, k, ivx) = w(i, j, k, irho)*w(i, j, k, ivx)
   wd(i, j, k, ivy) = wd(i, j, k, irho)*w(i, j, k, ivy) + w(i, j&
   &           , k, irho)*wd(i, j, k, ivy)
   w(i, j, k, ivy) = w(i, j, k, irho)*w(i, j, k, ivy)
   wd(i, j, k, ivz) = wd(i, j, k, irho)*w(i, j, k, ivz) + w(i, j&
   &           , k, irho)*wd(i, j, k, ivz)
   w(i, j, k, ivz) = w(i, j, k, irho)*w(i, j, k, ivz)
   wd(i, j, k, irhoe) = wd(i, j, k, irhoe) + pd(i, j, k)
   w(i, j, k, irhoe) = w(i, j, k, irhoe) + p(i, j, k)
   END DO
   END DO
   END DO
   DO k=2,kl
   DO j=2,jl
   wd(0, j, k, ivx) = wd(0, j, k, irho)*w(0, j, k, ivx) + w(0, j, k&
   &         , irho)*wd(0, j, k, ivx)
   w(0, j, k, ivx) = w(0, j, k, irho)*w(0, j, k, ivx)
   wd(0, j, k, ivy) = wd(0, j, k, irho)*w(0, j, k, ivy) + w(0, j, k&
   &         , irho)*wd(0, j, k, ivy)
   w(0, j, k, ivy) = w(0, j, k, irho)*w(0, j, k, ivy)
   wd(0, j, k, ivz) = wd(0, j, k, irho)*w(0, j, k, ivz) + w(0, j, k&
   &         , irho)*wd(0, j, k, ivz)
   w(0, j, k, ivz) = w(0, j, k, irho)*w(0, j, k, ivz)
   wd(0, j, k, irhoe) = wd(0, j, k, irhoe) + pd(0, j, k)
   w(0, j, k, irhoe) = w(0, j, k, irhoe) + p(0, j, k)
   wd(1, j, k, ivx) = wd(1, j, k, irho)*w(1, j, k, ivx) + w(1, j, k&
   &         , irho)*wd(1, j, k, ivx)
   w(1, j, k, ivx) = w(1, j, k, irho)*w(1, j, k, ivx)
   wd(1, j, k, ivy) = wd(1, j, k, irho)*w(1, j, k, ivy) + w(1, j, k&
   &         , irho)*wd(1, j, k, ivy)
   w(1, j, k, ivy) = w(1, j, k, irho)*w(1, j, k, ivy)
   wd(1, j, k, ivz) = wd(1, j, k, irho)*w(1, j, k, ivz) + w(1, j, k&
   &         , irho)*wd(1, j, k, ivz)
   w(1, j, k, ivz) = w(1, j, k, irho)*w(1, j, k, ivz)
   wd(1, j, k, irhoe) = wd(1, j, k, irhoe) + pd(1, j, k)
   w(1, j, k, irhoe) = w(1, j, k, irhoe) + p(1, j, k)
   wd(ie, j, k, ivx) = wd(ie, j, k, irho)*w(ie, j, k, ivx) + w(ie, &
   &         j, k, irho)*wd(ie, j, k, ivx)
   w(ie, j, k, ivx) = w(ie, j, k, irho)*w(ie, j, k, ivx)
   wd(ie, j, k, ivy) = wd(ie, j, k, irho)*w(ie, j, k, ivy) + w(ie, &
   &         j, k, irho)*wd(ie, j, k, ivy)
   w(ie, j, k, ivy) = w(ie, j, k, irho)*w(ie, j, k, ivy)
   wd(ie, j, k, ivz) = wd(ie, j, k, irho)*w(ie, j, k, ivz) + w(ie, &
   &         j, k, irho)*wd(ie, j, k, ivz)
   w(ie, j, k, ivz) = w(ie, j, k, irho)*w(ie, j, k, ivz)
   wd(ie, j, k, irhoe) = wd(ie, j, k, irhoe) + pd(ie, j, k)
   w(ie, j, k, irhoe) = w(ie, j, k, irhoe) + p(ie, j, k)
   wd(ib, j, k, ivx) = wd(ib, j, k, irho)*w(ib, j, k, ivx) + w(ib, &
   &         j, k, irho)*wd(ib, j, k, ivx)
   w(ib, j, k, ivx) = w(ib, j, k, irho)*w(ib, j, k, ivx)
   wd(ib, j, k, ivy) = wd(ib, j, k, irho)*w(ib, j, k, ivy) + w(ib, &
   &         j, k, irho)*wd(ib, j, k, ivy)
   w(ib, j, k, ivy) = w(ib, j, k, irho)*w(ib, j, k, ivy)
   wd(ib, j, k, ivz) = wd(ib, j, k, irho)*w(ib, j, k, ivz) + w(ib, &
   &         j, k, irho)*wd(ib, j, k, ivz)
   w(ib, j, k, ivz) = w(ib, j, k, irho)*w(ib, j, k, ivz)
   wd(ib, j, k, irhoe) = wd(ib, j, k, irhoe) + pd(ib, j, k)
   w(ib, j, k, irhoe) = w(ib, j, k, irhoe) + p(ib, j, k)
   END DO
   END DO
   DO k=2,kl
   DO i=2,il
   wd(i, 0, k, ivx) = wd(i, 0, k, irho)*w(i, 0, k, ivx) + w(i, 0, k&
   &         , irho)*wd(i, 0, k, ivx)
   w(i, 0, k, ivx) = w(i, 0, k, irho)*w(i, 0, k, ivx)
   wd(i, 0, k, ivy) = wd(i, 0, k, irho)*w(i, 0, k, ivy) + w(i, 0, k&
   &         , irho)*wd(i, 0, k, ivy)
   w(i, 0, k, ivy) = w(i, 0, k, irho)*w(i, 0, k, ivy)
   wd(i, 0, k, ivz) = wd(i, 0, k, irho)*w(i, 0, k, ivz) + w(i, 0, k&
   &         , irho)*wd(i, 0, k, ivz)
   w(i, 0, k, ivz) = w(i, 0, k, irho)*w(i, 0, k, ivz)
   wd(i, 0, k, irhoe) = wd(i, 0, k, irhoe) + pd(i, 0, k)
   w(i, 0, k, irhoe) = w(i, 0, k, irhoe) + p(i, 0, k)
   wd(i, 1, k, ivx) = wd(i, 1, k, irho)*w(i, 1, k, ivx) + w(i, 1, k&
   &         , irho)*wd(i, 1, k, ivx)
   w(i, 1, k, ivx) = w(i, 1, k, irho)*w(i, 1, k, ivx)
   wd(i, 1, k, ivy) = wd(i, 1, k, irho)*w(i, 1, k, ivy) + w(i, 1, k&
   &         , irho)*wd(i, 1, k, ivy)
   w(i, 1, k, ivy) = w(i, 1, k, irho)*w(i, 1, k, ivy)
   wd(i, 1, k, ivz) = wd(i, 1, k, irho)*w(i, 1, k, ivz) + w(i, 1, k&
   &         , irho)*wd(i, 1, k, ivz)
   w(i, 1, k, ivz) = w(i, 1, k, irho)*w(i, 1, k, ivz)
   wd(i, 1, k, irhoe) = wd(i, 1, k, irhoe) + pd(i, 1, k)
   w(i, 1, k, irhoe) = w(i, 1, k, irhoe) + p(i, 1, k)
   wd(i, je, k, ivx) = wd(i, je, k, irho)*w(i, je, k, ivx) + w(i, &
   &         je, k, irho)*wd(i, je, k, ivx)
   w(i, je, k, ivx) = w(i, je, k, irho)*w(i, je, k, ivx)
   wd(i, je, k, ivy) = wd(i, je, k, irho)*w(i, je, k, ivy) + w(i, &
   &         je, k, irho)*wd(i, je, k, ivy)
   w(i, je, k, ivy) = w(i, je, k, irho)*w(i, je, k, ivy)
   wd(i, je, k, ivz) = wd(i, je, k, irho)*w(i, je, k, ivz) + w(i, &
   &         je, k, irho)*wd(i, je, k, ivz)
   w(i, je, k, ivz) = w(i, je, k, irho)*w(i, je, k, ivz)
   wd(i, je, k, irhoe) = wd(i, je, k, irhoe) + pd(i, je, k)
   w(i, je, k, irhoe) = w(i, je, k, irhoe) + p(i, je, k)
   wd(i, jb, k, ivx) = wd(i, jb, k, irho)*w(i, jb, k, ivx) + w(i, &
   &         jb, k, irho)*wd(i, jb, k, ivx)
   w(i, jb, k, ivx) = w(i, jb, k, irho)*w(i, jb, k, ivx)
   wd(i, jb, k, ivy) = wd(i, jb, k, irho)*w(i, jb, k, ivy) + w(i, &
   &         jb, k, irho)*wd(i, jb, k, ivy)
   w(i, jb, k, ivy) = w(i, jb, k, irho)*w(i, jb, k, ivy)
   wd(i, jb, k, ivz) = wd(i, jb, k, irho)*w(i, jb, k, ivz) + w(i, &
   &         jb, k, irho)*wd(i, jb, k, ivz)
   w(i, jb, k, ivz) = w(i, jb, k, irho)*w(i, jb, k, ivz)
   wd(i, jb, k, irhoe) = wd(i, jb, k, irhoe) + pd(i, jb, k)
   w(i, jb, k, irhoe) = w(i, jb, k, irhoe) + p(i, jb, k)
   END DO
   END DO
   ! Initialize the dissipative residual to a certain times,
   ! possibly zero, the previously stored value. Owned cells
   ! only, because the halo values do not matter.
   DO k=2,kl
   DO j=2,jl
   DO i=2,il
   fwd(i, j, k, irho) = 0.0_8
   fw(i, j, k, irho) = sfil*fw(i, j, k, irho)
   fwd(i, j, k, imx) = 0.0_8
   fw(i, j, k, imx) = sfil*fw(i, j, k, imx)
   fwd(i, j, k, imy) = 0.0_8
   fw(i, j, k, imy) = sfil*fw(i, j, k, imy)
   fwd(i, j, k, imz) = 0.0_8
   fw(i, j, k, imz) = sfil*fw(i, j, k, imz)
   fwd(i, j, k, irhoe) = 0.0_8
   fw(i, j, k, irhoe) = sfil*fw(i, j, k, irhoe)
   END DO
   END DO
   END DO
   fwd = 0.0_8
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Dissipative fluxes in the i-direction.                         *
   !      *                                                                *
   !      ******************************************************************
   !
   DO k=2,kl
   DO j=2,jl
   x1d = -((shocksensor(2, j, k)-two*shocksensor(1, j, k)+&
   &         shocksensor(0, j, k))*sslimd/(shocksensor(2, j, k)+two*&
   &         shocksensor(1, j, k)+shocksensor(0, j, k)+sslim)**2)
   x1 = (shocksensor(2, j, k)-two*shocksensor(1, j, k)+shocksensor(&
   &         0, j, k))/(shocksensor(2, j, k)+two*shocksensor(1, j, k)+&
   &         shocksensor(0, j, k)+sslim)
   IF (x1 .GE. 0.) THEN
   dss1d = x1d
   dss1 = x1
   ELSE
   dss1d = -x1d
   dss1 = -x1
   END IF
   ! Loop in i-direction.
   DO i=1,il
   x2d = -((shocksensor(i+2, j, k)-two*shocksensor(i+1, j, k)+&
   &           shocksensor(i, j, k))*sslimd/(shocksensor(i+2, j, k)+two*&
   &           shocksensor(i+1, j, k)+shocksensor(i, j, k)+sslim)**2)
   x2 = (shocksensor(i+2, j, k)-two*shocksensor(i+1, j, k)+&
   &           shocksensor(i, j, k))/(shocksensor(i+2, j, k)+two*&
   &           shocksensor(i+1, j, k)+shocksensor(i, j, k)+sslim)
   IF (x2 .GE. 0.) THEN
   dss2d = x2d
   dss2 = x2
   ELSE
   dss2d = -x2d
   dss2 = -x2
   END IF
   ! Compute the dissipation coefficients for this face.
   ppor = zero
   IF (pori(i, j, k) .EQ. normalflux) ppor = half
   rradd = ppor*(radid(i, j, k)+radid(i+1, j, k))
   rrad = ppor*(radi(i, j, k)+radi(i+1, j, k))
   IF (dss1 .LT. dss2) THEN
   y1d = dss2d
   y1 = dss2
   ELSE
   y1d = dss1d
   y1 = dss1
   END IF
   IF (dssmax .GT. y1) THEN
   min1d = y1d
   min1 = y1
   ELSE
   min1 = dssmax
   min1d = 0.0_8
   END IF
   ! Modification for FD Preconditioner Note: This lumping
   ! actually still results in a greater than 3 cell stencil
   ! in any direction. Since this seems to work slightly
   ! better than the dis2=sigma*fis4*rrad, we will just use
   ! a 5-cell stencil for doing the PC
   dis2d = fis2*(rradd*min1+rrad*min1d) + sigma*fis4*rradd
   dis2 = fis2*rrad*min1 + sigma*fis4*rrad
   ! Compute and scatter the dissipative flux.
   ! Density. Store it in the mass flow of the
   ! appropriate sliding mesh interface.
   ddwd = wd(i+1, j, k, irho) - wd(i, j, k, irho)
   ddw = w(i+1, j, k, irho) - w(i, j, k, irho)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i+1, j, k, irho) = fwd(i+1, j, k, irho) + fsd
   fw(i+1, j, k, irho) = fw(i+1, j, k, irho) + fs
   fwd(i, j, k, irho) = fwd(i, j, k, irho) - fsd
   fw(i, j, k, irho) = fw(i, j, k, irho) - fs
   ! X-momentum.
   ddwd = wd(i+1, j, k, ivx) - wd(i, j, k, ivx)
   ddw = w(i+1, j, k, ivx) - w(i, j, k, ivx)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i+1, j, k, imx) = fwd(i+1, j, k, imx) + fsd
   fw(i+1, j, k, imx) = fw(i+1, j, k, imx) + fs
   fwd(i, j, k, imx) = fwd(i, j, k, imx) - fsd
   fw(i, j, k, imx) = fw(i, j, k, imx) - fs
   ! Y-momentum.
   ddwd = wd(i+1, j, k, ivy) - wd(i, j, k, ivy)
   ddw = w(i+1, j, k, ivy) - w(i, j, k, ivy)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i+1, j, k, imy) = fwd(i+1, j, k, imy) + fsd
   fw(i+1, j, k, imy) = fw(i+1, j, k, imy) + fs
   fwd(i, j, k, imy) = fwd(i, j, k, imy) - fsd
   fw(i, j, k, imy) = fw(i, j, k, imy) - fs
   ! Z-momentum.
   ddwd = wd(i+1, j, k, ivz) - wd(i, j, k, ivz)
   ddw = w(i+1, j, k, ivz) - w(i, j, k, ivz)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i+1, j, k, imz) = fwd(i+1, j, k, imz) + fsd
   fw(i+1, j, k, imz) = fw(i+1, j, k, imz) + fs
   fwd(i, j, k, imz) = fwd(i, j, k, imz) - fsd
   fw(i, j, k, imz) = fw(i, j, k, imz) - fs
   ! Energy.
   ddwd = wd(i+1, j, k, irhoe) - wd(i, j, k, irhoe)
   ddw = w(i+1, j, k, irhoe) - w(i, j, k, irhoe)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i+1, j, k, irhoe) = fwd(i+1, j, k, irhoe) + fsd
   fw(i+1, j, k, irhoe) = fw(i+1, j, k, irhoe) + fs
   fwd(i, j, k, irhoe) = fwd(i, j, k, irhoe) - fsd
   fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
   ! Set dss1 to dss2 for the next face.
   dss1d = dss2d
   dss1 = dss2
   END DO
   END DO
   END DO
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Dissipative fluxes in the j-direction.                         *
   !      *                                                                *
   !      ******************************************************************
   !
   DO k=2,kl
   DO i=2,il
   x3d = -((shocksensor(i, 2, k)-two*shocksensor(i, 1, k)+&
   &         shocksensor(i, 0, k))*sslimd/(shocksensor(i, 2, k)+two*&
   &         shocksensor(i, 1, k)+shocksensor(i, 0, k)+sslim)**2)
   x3 = (shocksensor(i, 2, k)-two*shocksensor(i, 1, k)+shocksensor(&
   &         i, 0, k))/(shocksensor(i, 2, k)+two*shocksensor(i, 1, k)+&
   &         shocksensor(i, 0, k)+sslim)
   IF (x3 .GE. 0.) THEN
   dss1d = x3d
   dss1 = x3
   ELSE
   dss1d = -x3d
   dss1 = -x3
   END IF
   ! Loop in j-direction.
   DO j=1,jl
   x4d = -((shocksensor(i, j+2, k)-two*shocksensor(i, j+1, k)+&
   &           shocksensor(i, j, k))*sslimd/(shocksensor(i, j+2, k)+two*&
   &           shocksensor(i, j+1, k)+shocksensor(i, j, k)+sslim)**2)
   x4 = (shocksensor(i, j+2, k)-two*shocksensor(i, j+1, k)+&
   &           shocksensor(i, j, k))/(shocksensor(i, j+2, k)+two*&
   &           shocksensor(i, j+1, k)+shocksensor(i, j, k)+sslim)
   IF (x4 .GE. 0.) THEN
   dss2d = x4d
   dss2 = x4
   ELSE
   dss2d = -x4d
   dss2 = -x4
   END IF
   ! Compute the dissipation coefficients for this face.
   ppor = zero
   IF (porj(i, j, k) .EQ. normalflux) ppor = half
   rradd = ppor*(radjd(i, j, k)+radjd(i, j+1, k))
   rrad = ppor*(radj(i, j, k)+radj(i, j+1, k))
   IF (dss1 .LT. dss2) THEN
   y2d = dss2d
   y2 = dss2
   ELSE
   y2d = dss1d
   y2 = dss1
   END IF
   IF (dssmax .GT. y2) THEN
   min2d = y2d
   min2 = y2
   ELSE
   min2 = dssmax
   min2d = 0.0_8
   END IF
   ! Modification for FD Preconditioner
   dis2d = fis2*(rradd*min2+rrad*min2d) + sigma*fis4*rradd
   dis2 = fis2*rrad*min2 + sigma*fis4*rrad
   ! Compute and scatter the dissipative flux.
   ! Density. Store it in the mass flow of the
   ! appropriate sliding mesh interface.
   ddwd = wd(i, j+1, k, irho) - wd(i, j, k, irho)
   ddw = w(i, j+1, k, irho) - w(i, j, k, irho)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j+1, k, irho) = fwd(i, j+1, k, irho) + fsd
   fw(i, j+1, k, irho) = fw(i, j+1, k, irho) + fs
   fwd(i, j, k, irho) = fwd(i, j, k, irho) - fsd
   fw(i, j, k, irho) = fw(i, j, k, irho) - fs
   ! X-momentum.
   ddwd = wd(i, j+1, k, ivx) - wd(i, j, k, ivx)
   ddw = w(i, j+1, k, ivx) - w(i, j, k, ivx)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j+1, k, imx) = fwd(i, j+1, k, imx) + fsd
   fw(i, j+1, k, imx) = fw(i, j+1, k, imx) + fs
   fwd(i, j, k, imx) = fwd(i, j, k, imx) - fsd
   fw(i, j, k, imx) = fw(i, j, k, imx) - fs
   ! Y-momentum.
   ddwd = wd(i, j+1, k, ivy) - wd(i, j, k, ivy)
   ddw = w(i, j+1, k, ivy) - w(i, j, k, ivy)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j+1, k, imy) = fwd(i, j+1, k, imy) + fsd
   fw(i, j+1, k, imy) = fw(i, j+1, k, imy) + fs
   fwd(i, j, k, imy) = fwd(i, j, k, imy) - fsd
   fw(i, j, k, imy) = fw(i, j, k, imy) - fs
   ! Z-momentum.
   ddwd = wd(i, j+1, k, ivz) - wd(i, j, k, ivz)
   ddw = w(i, j+1, k, ivz) - w(i, j, k, ivz)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j+1, k, imz) = fwd(i, j+1, k, imz) + fsd
   fw(i, j+1, k, imz) = fw(i, j+1, k, imz) + fs
   fwd(i, j, k, imz) = fwd(i, j, k, imz) - fsd
   fw(i, j, k, imz) = fw(i, j, k, imz) - fs
   ! Energy.
   ddwd = wd(i, j+1, k, irhoe) - wd(i, j, k, irhoe)
   ddw = w(i, j+1, k, irhoe) - w(i, j, k, irhoe)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j+1, k, irhoe) = fwd(i, j+1, k, irhoe) + fsd
   fw(i, j+1, k, irhoe) = fw(i, j+1, k, irhoe) + fs
   fwd(i, j, k, irhoe) = fwd(i, j, k, irhoe) - fsd
   fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
   ! Set dss1 to dss2 for the next face.
   dss1d = dss2d
   dss1 = dss2
   END DO
   END DO
   END DO
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Dissipative fluxes in the k-direction.                         *
   !      *                                                                *
   !      ******************************************************************
   !
   DO j=2,jl
   DO i=2,il
   x5d = -((shocksensor(i, j, 2)-two*shocksensor(i, j, 1)+&
   &         shocksensor(i, j, 0))*sslimd/(shocksensor(i, j, 2)+two*&
   &         shocksensor(i, j, 1)+shocksensor(i, j, 0)+sslim)**2)
   x5 = (shocksensor(i, j, 2)-two*shocksensor(i, j, 1)+shocksensor(&
   &         i, j, 0))/(shocksensor(i, j, 2)+two*shocksensor(i, j, 1)+&
   &         shocksensor(i, j, 0)+sslim)
   IF (x5 .GE. 0.) THEN
   dss1d = x5d
   dss1 = x5
   ELSE
   dss1d = -x5d
   dss1 = -x5
   END IF
   ! Loop in k-direction.
   DO k=1,kl
   x6d = -((shocksensor(i, j, k+2)-two*shocksensor(i, j, k+1)+&
   &           shocksensor(i, j, k))*sslimd/(shocksensor(i, j, k+2)+two*&
   &           shocksensor(i, j, k+1)+shocksensor(i, j, k)+sslim)**2)
   x6 = (shocksensor(i, j, k+2)-two*shocksensor(i, j, k+1)+&
   &           shocksensor(i, j, k))/(shocksensor(i, j, k+2)+two*&
   &           shocksensor(i, j, k+1)+shocksensor(i, j, k)+sslim)
   IF (x6 .GE. 0.) THEN
   dss2d = x6d
   dss2 = x6
   ELSE
   dss2d = -x6d
   dss2 = -x6
   END IF
   ! Compute the dissipation coefficients for this face.
   ppor = zero
   IF (pork(i, j, k) .EQ. normalflux) ppor = half
   rradd = ppor*(radkd(i, j, k)+radkd(i, j, k+1))
   rrad = ppor*(radk(i, j, k)+radk(i, j, k+1))
   IF (dss1 .LT. dss2) THEN
   y3d = dss2d
   y3 = dss2
   ELSE
   y3d = dss1d
   y3 = dss1
   END IF
   IF (dssmax .GT. y3) THEN
   min3d = y3d
   min3 = y3
   ELSE
   min3 = dssmax
   min3d = 0.0_8
   END IF
   ! Modification for FD Preconditioner
   dis2d = fis2*(rradd*min3+rrad*min3d) + sigma*fis4*rradd
   dis2 = fis2*rrad*min3 + sigma*fis4*rrad
   ! Compute and scatter the dissipative flux.
   ! Density. Store it in the mass flow of the
   ! appropriate sliding mesh interface.
   ddwd = wd(i, j, k+1, irho) - wd(i, j, k, irho)
   ddw = w(i, j, k+1, irho) - w(i, j, k, irho)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j, k+1, irho) = fwd(i, j, k+1, irho) + fsd
   fw(i, j, k+1, irho) = fw(i, j, k+1, irho) + fs
   fwd(i, j, k, irho) = fwd(i, j, k, irho) - fsd
   fw(i, j, k, irho) = fw(i, j, k, irho) - fs
   ! X-momentum.
   ddwd = wd(i, j, k+1, ivx) - wd(i, j, k, ivx)
   ddw = w(i, j, k+1, ivx) - w(i, j, k, ivx)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j, k+1, imx) = fwd(i, j, k+1, imx) + fsd
   fw(i, j, k+1, imx) = fw(i, j, k+1, imx) + fs
   fwd(i, j, k, imx) = fwd(i, j, k, imx) - fsd
   fw(i, j, k, imx) = fw(i, j, k, imx) - fs
   ! Y-momentum.
   ddwd = wd(i, j, k+1, ivy) - wd(i, j, k, ivy)
   ddw = w(i, j, k+1, ivy) - w(i, j, k, ivy)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j, k+1, imy) = fwd(i, j, k+1, imy) + fsd
   fw(i, j, k+1, imy) = fw(i, j, k+1, imy) + fs
   fwd(i, j, k, imy) = fwd(i, j, k, imy) - fsd
   fw(i, j, k, imy) = fw(i, j, k, imy) - fs
   ! Z-momentum.
   ddwd = wd(i, j, k+1, ivz) - wd(i, j, k, ivz)
   ddw = w(i, j, k+1, ivz) - w(i, j, k, ivz)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j, k+1, imz) = fwd(i, j, k+1, imz) + fsd
   fw(i, j, k+1, imz) = fw(i, j, k+1, imz) + fs
   fwd(i, j, k, imz) = fwd(i, j, k, imz) - fsd
   fw(i, j, k, imz) = fw(i, j, k, imz) - fs
   ! Energy.
   ddwd = wd(i, j, k+1, irhoe) - wd(i, j, k, irhoe)
   ddw = w(i, j, k+1, irhoe) - w(i, j, k, irhoe)
   fsd = dis2d*ddw + dis2*ddwd
   fs = dis2*ddw
   fwd(i, j, k+1, irhoe) = fwd(i, j, k+1, irhoe) + fsd
   fw(i, j, k+1, irhoe) = fw(i, j, k+1, irhoe) + fs
   fwd(i, j, k, irhoe) = fwd(i, j, k, irhoe) - fsd
   fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
   ! Set dss1 to dss2 for the next face.
   dss1d = dss2d
   dss1 = dss2
   END DO
   END DO
   END DO
   ! Replace rho times the total enthalpy by the total energy and
   ! store the velocities again instead of the momentum. Only for
   ! those entries that have been altered, i.e. ignore the
   ! corner halo's.
   DO k=0,kb
   DO j=2,jl
   DO i=2,il
   rhoid = -(one*wd(i, j, k, irho)/w(i, j, k, irho)**2)
   rhoi = one/w(i, j, k, irho)
   wd(i, j, k, ivx) = wd(i, j, k, ivx)*rhoi + w(i, j, k, ivx)*&
   &           rhoid
   w(i, j, k, ivx) = w(i, j, k, ivx)*rhoi
   wd(i, j, k, ivy) = wd(i, j, k, ivy)*rhoi + w(i, j, k, ivy)*&
   &           rhoid
   w(i, j, k, ivy) = w(i, j, k, ivy)*rhoi
   wd(i, j, k, ivz) = wd(i, j, k, ivz)*rhoi + w(i, j, k, ivz)*&
   &           rhoid
   w(i, j, k, ivz) = w(i, j, k, ivz)*rhoi
   wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - pd(i, j, k)
   w(i, j, k, irhoe) = w(i, j, k, irhoe) - p(i, j, k)
   END DO
   END DO
   END DO
   DO k=2,kl
   DO j=2,jl
   rhoid = -(one*wd(0, j, k, irho)/w(0, j, k, irho)**2)
   rhoi = one/w(0, j, k, irho)
   wd(0, j, k, ivx) = wd(0, j, k, ivx)*rhoi + w(0, j, k, ivx)*rhoid
   w(0, j, k, ivx) = w(0, j, k, ivx)*rhoi
   wd(0, j, k, ivy) = wd(0, j, k, ivy)*rhoi + w(0, j, k, ivy)*rhoid
   w(0, j, k, ivy) = w(0, j, k, ivy)*rhoi
   wd(0, j, k, ivz) = wd(0, j, k, ivz)*rhoi + w(0, j, k, ivz)*rhoid
   w(0, j, k, ivz) = w(0, j, k, ivz)*rhoi
   wd(0, j, k, irhoe) = wd(0, j, k, irhoe) - pd(0, j, k)
   w(0, j, k, irhoe) = w(0, j, k, irhoe) - p(0, j, k)
   rhoid = -(one*wd(1, j, k, irho)/w(1, j, k, irho)**2)
   rhoi = one/w(1, j, k, irho)
   wd(1, j, k, ivx) = wd(1, j, k, ivx)*rhoi + w(1, j, k, ivx)*rhoid
   w(1, j, k, ivx) = w(1, j, k, ivx)*rhoi
   wd(1, j, k, ivy) = wd(1, j, k, ivy)*rhoi + w(1, j, k, ivy)*rhoid
   w(1, j, k, ivy) = w(1, j, k, ivy)*rhoi
   wd(1, j, k, ivz) = wd(1, j, k, ivz)*rhoi + w(1, j, k, ivz)*rhoid
   w(1, j, k, ivz) = w(1, j, k, ivz)*rhoi
   wd(1, j, k, irhoe) = wd(1, j, k, irhoe) - pd(1, j, k)
   w(1, j, k, irhoe) = w(1, j, k, irhoe) - p(1, j, k)
   rhoid = -(one*wd(ie, j, k, irho)/w(ie, j, k, irho)**2)
   rhoi = one/w(ie, j, k, irho)
   wd(ie, j, k, ivx) = wd(ie, j, k, ivx)*rhoi + w(ie, j, k, ivx)*&
   &         rhoid
   w(ie, j, k, ivx) = w(ie, j, k, ivx)*rhoi
   wd(ie, j, k, ivy) = wd(ie, j, k, ivy)*rhoi + w(ie, j, k, ivy)*&
   &         rhoid
   w(ie, j, k, ivy) = w(ie, j, k, ivy)*rhoi
   wd(ie, j, k, ivz) = wd(ie, j, k, ivz)*rhoi + w(ie, j, k, ivz)*&
   &         rhoid
   w(ie, j, k, ivz) = w(ie, j, k, ivz)*rhoi
   wd(ie, j, k, irhoe) = wd(ie, j, k, irhoe) - pd(ie, j, k)
   w(ie, j, k, irhoe) = w(ie, j, k, irhoe) - p(ie, j, k)
   rhoid = -(one*wd(ib, j, k, irho)/w(ib, j, k, irho)**2)
   rhoi = one/w(ib, j, k, irho)
   wd(ib, j, k, ivx) = wd(ib, j, k, ivx)*rhoi + w(ib, j, k, ivx)*&
   &         rhoid
   w(ib, j, k, ivx) = w(ib, j, k, ivx)*rhoi
   wd(ib, j, k, ivy) = wd(ib, j, k, ivy)*rhoi + w(ib, j, k, ivy)*&
   &         rhoid
   w(ib, j, k, ivy) = w(ib, j, k, ivy)*rhoi
   wd(ib, j, k, ivz) = wd(ib, j, k, ivz)*rhoi + w(ib, j, k, ivz)*&
   &         rhoid
   w(ib, j, k, ivz) = w(ib, j, k, ivz)*rhoi
   wd(ib, j, k, irhoe) = wd(ib, j, k, irhoe) - pd(ib, j, k)
   w(ib, j, k, irhoe) = w(ib, j, k, irhoe) - p(ib, j, k)
   END DO
   END DO
   DO k=2,kl
   DO i=2,il
   rhoid = -(one*wd(i, 0, k, irho)/w(i, 0, k, irho)**2)
   rhoi = one/w(i, 0, k, irho)
   wd(i, 0, k, ivx) = wd(i, 0, k, ivx)*rhoi + w(i, 0, k, ivx)*rhoid
   w(i, 0, k, ivx) = w(i, 0, k, ivx)*rhoi
   wd(i, 0, k, ivy) = wd(i, 0, k, ivy)*rhoi + w(i, 0, k, ivy)*rhoid
   w(i, 0, k, ivy) = w(i, 0, k, ivy)*rhoi
   wd(i, 0, k, ivz) = wd(i, 0, k, ivz)*rhoi + w(i, 0, k, ivz)*rhoid
   w(i, 0, k, ivz) = w(i, 0, k, ivz)*rhoi
   wd(i, 0, k, irhoe) = wd(i, 0, k, irhoe) - pd(i, 0, k)
   w(i, 0, k, irhoe) = w(i, 0, k, irhoe) - p(i, 0, k)
   rhoid = -(one*wd(i, 1, k, irho)/w(i, 1, k, irho)**2)
   rhoi = one/w(i, 1, k, irho)
   wd(i, 1, k, ivx) = wd(i, 1, k, ivx)*rhoi + w(i, 1, k, ivx)*rhoid
   w(i, 1, k, ivx) = w(i, 1, k, ivx)*rhoi
   wd(i, 1, k, ivy) = wd(i, 1, k, ivy)*rhoi + w(i, 1, k, ivy)*rhoid
   w(i, 1, k, ivy) = w(i, 1, k, ivy)*rhoi
   wd(i, 1, k, ivz) = wd(i, 1, k, ivz)*rhoi + w(i, 1, k, ivz)*rhoid
   w(i, 1, k, ivz) = w(i, 1, k, ivz)*rhoi
   wd(i, 1, k, irhoe) = wd(i, 1, k, irhoe) - pd(i, 1, k)
   w(i, 1, k, irhoe) = w(i, 1, k, irhoe) - p(i, 1, k)
   rhoid = -(one*wd(i, je, k, irho)/w(i, je, k, irho)**2)
   rhoi = one/w(i, je, k, irho)
   wd(i, je, k, ivx) = wd(i, je, k, ivx)*rhoi + w(i, je, k, ivx)*&
   &         rhoid
   w(i, je, k, ivx) = w(i, je, k, ivx)*rhoi
   wd(i, je, k, ivy) = wd(i, je, k, ivy)*rhoi + w(i, je, k, ivy)*&
   &         rhoid
   w(i, je, k, ivy) = w(i, je, k, ivy)*rhoi
   wd(i, je, k, ivz) = wd(i, je, k, ivz)*rhoi + w(i, je, k, ivz)*&
   &         rhoid
   w(i, je, k, ivz) = w(i, je, k, ivz)*rhoi
   wd(i, je, k, irhoe) = wd(i, je, k, irhoe) - pd(i, je, k)
   w(i, je, k, irhoe) = w(i, je, k, irhoe) - p(i, je, k)
   rhoid = -(one*wd(i, jb, k, irho)/w(i, jb, k, irho)**2)
   rhoi = one/w(i, jb, k, irho)
   wd(i, jb, k, ivx) = wd(i, jb, k, ivx)*rhoi + w(i, jb, k, ivx)*&
   &         rhoid
   w(i, jb, k, ivx) = w(i, jb, k, ivx)*rhoi
   wd(i, jb, k, ivy) = wd(i, jb, k, ivy)*rhoi + w(i, jb, k, ivy)*&
   &         rhoid
   w(i, jb, k, ivy) = w(i, jb, k, ivy)*rhoi
   wd(i, jb, k, ivz) = wd(i, jb, k, ivz)*rhoi + w(i, jb, k, ivz)*&
   &         rhoid
   w(i, jb, k, ivz) = w(i, jb, k, ivz)*rhoi
   wd(i, jb, k, irhoe) = wd(i, jb, k, irhoe) - pd(i, jb, k)
   w(i, jb, k, irhoe) = w(i, jb, k, irhoe) - p(i, jb, k)
   END DO
   END DO
   END IF
   END SUBROUTINE INVISCIDDISSFLUXSCALARAPPROX_D
