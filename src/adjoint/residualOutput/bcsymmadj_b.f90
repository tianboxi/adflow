!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 2.2.4 (r2308) - 03/04/2008 10:03
!  
!  Differentiation of bcsymmadj in reverse (adjoint) mode:
!   gradient, with respect to input variables: padj wadj normadj
!   of linear combination of output variables: padj wadj normadj
!
!      ******************************************************************
!      *                                                                *
!      * File:          bcSymmAdj.f90                                   *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 04-17-2008                                      *
!      * Last modified: 04-17-2008                                      *
!      *                                                                *
!      ******************************************************************
!
SUBROUTINE BCSYMMADJ_B(wadj, wadjb, padj, padjb, normadj, normadjb, &
&  icell, jcell, kcell, secondhalo)
  USE blockpointers, ONLY : ie, ib, il, je, jb, jl, ke, kb, kl, nbocos&
&  , gamma, bcfaceid, bctype, bcdata
  USE bctypes
  USE constants
  USE iteration
  USE flowvarrefstate
  IMPLICIT NONE
!  enddo nHalo
!
!      ******************************************************************
!      *                                                                *
!      * bcSymmAdj applies the symmetry boundary conditions to a single *
!      * cell stencil.
!      * It is assumed that the pointers in blockPointers are already   *
!      * set to the correct block on the correct grid level.            *
!      *                                                                *
!      * In case also the second halo must be set the loop over the     *
!      * boundary subfaces is executed twice. This is the only correct  *
!      * way in case the block contains only 1 cell between two         *
!      * symmetry planes, i.e. a 2D problem.                            *
!      *                                                                *
!      ******************************************************************
!
!nw
!nt1mg,nt2mg
!
!      Subroutine arguments.
!
  LOGICAL, INTENT(IN) :: secondhalo
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2, nw) :: wadj
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2, nw) :: wadjb
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2) :: padj
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2) :: padjb
  REAL(KIND=REALTYPE), DIMENSION(nbocos, -2:2, -2:2, 3) :: normadj
  REAL(KIND=REALTYPE), DIMENSION(nbocos, -2:2, -2:2, 3) :: normadjb
!
!      Local variables.
!
  INTEGER(KIND=INTTYPE) :: kk, mm, nn, i, j, l, ii, jj
  REAL(KIND=REALTYPE) :: vn, nnx, nny, nnz
  REAL(KIND=REALTYPE) :: vnb, nnxb, nnyb, nnzb
!real(kind=realType), dimension(:,:),   pointer :: gamma1, gamma2
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2) :: gamma1, gamma2
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, nw) :: wadj0, wadj1
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, nw) :: wadj0b, wadj1b
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, nw) :: wadj2, wadj3
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, nw) :: wadj2b, wadj3b
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2) :: padj0, padj1
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2) :: padj0b, padj1b
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2) :: padj2, padj3
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2) :: padj2b, padj3b
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2) :: rlvadj, revadj
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2) :: rlvadj1, rlvadj2
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2) :: revadj1, revadj2
  INTEGER(KIND=INTTYPE) :: icell, jcell, kcell
  INTEGER(KIND=INTTYPE) :: isbeg, jsbeg, ksbeg, isend, jsend, ksend
  INTEGER(KIND=INTTYPE) :: ibbeg, jbbeg, kbbeg, ibend, jbend, kbend
  INTEGER(KIND=INTTYPE) :: icbeg, jcbeg, kcbeg, icend, jcend, kcend
  INTEGER(KIND=INTTYPE) :: ioffset, joffset, koffset
  LOGICAL :: computebc
  INTEGER :: branch
  INTEGER :: ad_from
  INTEGER :: ad_to
  INTEGER :: ad_from0
  INTEGER :: ad_to0
  REAL(KIND=REALTYPE) :: tempb0
  REAL(KIND=REALTYPE) :: tempb
!!$!
!!$!      Interfaces
!!$!
!!$       interface
!!$         subroutine setBcPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
!!$                                  rev1, rev2, offset)
!!$           use blockPointers
!!$           implicit none
!!$
!!$           integer(kind=intType), intent(in) :: nn, offset
!!$           real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
!!$           real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
!!$           real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
!!$           real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
!!$         end subroutine setBcPointers
!!$       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!print *,'in bcSymmAdj'
!Accounted for inside bocos loop
!!$  ! Set the value of kk; kk == 0 means only single halo, kk == 1
!!$  ! double halo.
!!$
!!$  kk = 0
!!$  if( secondHalo ) kk = 1
!!$
!!$  ! Loop over the number of times the halo computation must be done.
!!$
!!$  nHalo: do mm=0,kk
! Loop over the boundary condition subfaces of this block.
bocos:DO nn=1,nbocos
    CALL PUSHINTEGER4(kbend)
    CALL PUSHINTEGER4(jbend)
    CALL PUSHINTEGER4(ibend)
    CALL PUSHINTEGER4(kbbeg)
    CALL PUSHINTEGER4(jbbeg)
    CALL PUSHINTEGER4(ibbeg)
    CALL PUSHINTEGER4(ksend)
    CALL PUSHINTEGER4(jsend)
    CALL PUSHINTEGER4(isend)
    CALL PUSHINTEGER4(ksbeg)
    CALL PUSHINTEGER4(jsbeg)
    CALL PUSHINTEGER4(isbeg)
    CALL CHECKOVERLAPADJ(nn, icell, jcell, kcell, isbeg, jsbeg, ksbeg, &
&                   isend, jsend, ksend, ibbeg, jbbeg, kbbeg, ibend, &
&                   jbend, kbend, computebc)
!(nn,icell,jcell,kcell,isbeg,jsbeg,&
!     ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,jbend,kbend,&
!     computeBC)
    IF (computebc) THEN
! Check for symmetry boundary condition.
      IF (bctype(nn) .EQ. symm) THEN
        CALL PUSHBOOLEAN(secondhalo)
        CALL PUSHREAL8ARRAY(wadj3, 5**2*nw)
        CALL PUSHREAL8ARRAY(wadj2, 5**2*nw)
!!$             ! Nullify the pointers, because some compilers require that.
!!$
!!$             nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
!!$              ! Set the pointers to the correct subface.
!!$              call setBcPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
!!$                                rev1, rev2, mm)
!Copy the states and other parameters to subfaces
        CALL EXTRACTBCSTATESADJ(nn, wadj, padj, wadj0, wadj1, wadj2, &
&                          wadj3, padj0, padj1, padj2, padj3, rlvadj, &
&                          revadj, rlvadj1, rlvadj2, revadj1, revadj2, &
&                          ioffset, joffset, koffset, icell, jcell, &
&                          kcell, isbeg, jsbeg, ksbeg, isend, jsend, &
&                          ksend, ibbeg, jbbeg, kbbeg, ibend, jbend, &
&                          kbend, icbeg, jcbeg, icend, jcend, secondhalo&
&                         )
!(nn,wAdj,pAdj, wAdj1, wAdj2, pAdj1, pAdj2,&
!     rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,iOffset,&
!     jOffset, kOffset,iCell, jCell,kCell,&
!     isbeg,jsbeg,ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,&
!     jbend,kbend,icbeg,jcbeg,icend,jcend)
! Set the additional pointers for gamma1 and gamma2.
!              select case (BCFaceID(nn))
!              case (iMin)
!                 gamma1 => gamma(1, 1:,1:); gamma2 => gamma(2, 1:,1:)
!              case (iMax)
!                 gamma1 => gamma(ie,1:,1:); gamma2 => gamma(il,1:,1:)
!              case (jMin)
!                 gamma1 => gamma(1:,1, 1:); gamma2 => gamma(1:,2, 1:)
!              case (jMax)
!                 gamma1 => gamma(1:,je,1:); gamma2 => gamma(1:,jl,1:)
!              case (kMin)
!                 gamma1 => gamma(1:,1:,1 ); gamma2 => gamma(1:,1:,2 )
!              case (kMax)
!                 gamma1 => gamma(1:,1:,ke); gamma2 => gamma(1:,1:,kl)
!              end select
        ad_from = jcbeg
! Loop over the generic subface to set the state in the
! halo cells.
!!$             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
!!$               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
        DO j=ad_from,jcend
          ad_from0 = icbeg
          DO i=ad_from0,icend
            CALL PUSHINTEGER4(ii)
            ii = i - ioffset
            CALL PUSHINTEGER4(jj)
            jj = j - joffset
            CALL PUSHREAL8(nnx)
!print *,'i,j',i,j,jcBeg,jcEnd,ii,jj
! Store the three components of the unit normal a
! bit easier.
!BCData(nn)%norm(i,j,1)
            nnx = normadj(nn, ii, jj, 1)
            CALL PUSHREAL8(nny)
!BCData(nn)%norm(i,j,2)
            nny = normadj(nn, ii, jj, 2)
            CALL PUSHREAL8(nnz)
!BCData(nn)%norm(i,j,3)
            nnz = normadj(nn, ii, jj, 3)
            CALL PUSHREAL8(vn)
! Determine twice the normal velocity component,
! which must be substracted from the donor velocity
! to obtain the halo velocity.
            vn = two*(wadj2(ii, jj, ivx)*nnx+wadj2(ii, jj, ivy)*nny+&
&              wadj2(ii, jj, ivz)*nnz)
! Determine the flow variables in the halo cell.
            wadj1(ii, jj, irho) = wadj2(ii, jj, irho)
            wadj1(ii, jj, ivx) = wadj2(ii, jj, ivx) - vn*nnx
            wadj1(ii, jj, ivy) = wadj2(ii, jj, ivy) - vn*nny
            wadj1(ii, jj, ivz) = wadj2(ii, jj, ivz) - vn*nnz
            wadj1(ii, jj, irhoe) = wadj2(ii, jj, irhoe)
! Simply copy the turbulent variables.
            DO l=nt1mg,nt2mg
              wadj1(ii, jj, l) = wadj2(ii, jj, l)
            END DO
! Set the pressure and gamma and possibly the
! laminar and eddy viscosity in the halo.
            padj1(ii, jj) = padj2(ii, jj)
            IF (viscous) rlvadj1(ii, jj) = rlvadj2(ii, jj)
            IF (eddymodel) revadj1(ii, jj) = revadj2(ii, jj)
            IF (secondhalo) THEN
              CALL PUSHREAL8(vn)
! Determine twice the normal velocity component,
! which must be substracted from the donor velocity
! to obtain the halo velocity.
              vn = two*(wadj3(ii, jj, ivx)*nnx+wadj3(ii, jj, ivy)*nny+&
&                wadj3(ii, jj, ivz)*nnz)
! Determine the flow variables in the halo cell.
              wadj0(ii, jj, irho) = wadj3(ii, jj, irho)
              wadj0(ii, jj, ivx) = wadj3(ii, jj, ivx) - vn*nnx
              wadj0(ii, jj, ivy) = wadj3(ii, jj, ivy) - vn*nny
              wadj0(ii, jj, ivz) = wadj3(ii, jj, ivz) - vn*nnz
              wadj0(ii, jj, irhoe) = wadj3(ii, jj, irhoe)
! Simply copy the turbulent variables.
              DO l=nt1mg,nt2mg
                wadj1(ii, jj, l) = wadj3(ii, jj, l)
              END DO
! Set the pressure and gamma and possibly the
! laminar and eddy viscosity in the halo.
              padj0(ii, jj) = padj3(ii, jj)
!!$                       if( viscous )   rlvAdj0(ii,jj) = rlvAdj3(ii,jj)
!!$                       if( eddyModel ) revAdj0(ii,jj) = revAdj3(ii,jj)
              CALL PUSHINTEGER4(2)
            ELSE
              CALL PUSHINTEGER4(1)
            END IF
          END DO
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from0)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from)
        CALL PUSHREAL8ARRAY(wadj, 5**3*nw)
        CALL REPLACEBCSTATESADJ(nn, wadj0, wadj1, wadj2, wadj3, padj0, &
&                          padj1, padj2, padj3, rlvadj1, rlvadj2, &
&                          revadj1, revadj2, icell, jcell, kcell, wadj, &
&                          padj, rlvadj, revadj, secondhalo)
        CALL PUSHINTEGER4(3)
      ELSE
        CALL PUSHINTEGER4(2)
      END IF
    ELSE
      CALL PUSHINTEGER4(1)
    END IF
  END DO bocos
  padj0b(-2:2, -2:2) = 0.0
  padj1b(-2:2, -2:2) = 0.0
  padj2b(-2:2, -2:2) = 0.0
  padj3b(-2:2, -2:2) = 0.0
  wadj0b(-2:2, -2:2, :) = 0.0
  wadj1b(-2:2, -2:2, :) = 0.0
  wadj2b(-2:2, -2:2, :) = 0.0
  wadj3b(-2:2, -2:2, :) = 0.0
  DO nn=nbocos,1,-1
    CALL POPINTEGER4(branch)
    IF (.NOT.branch .LT. 3) THEN
      CALL POPREAL8ARRAY(wadj, 5**3*nw)
      CALL REPLACEBCSTATESADJ_B(nn, wadj0, wadj0b, wadj1, wadj1b, wadj2&
&                          , wadj3, padj0, padj0b, padj1, padj1b, padj2&
&                          , padj3, rlvadj1, rlvadj2, revadj1, revadj2, &
&                          icell, jcell, kcell, wadj, wadjb, padj, padjb&
&                          , rlvadj, revadj, secondhalo)
      CALL POPINTEGER4(ad_from)
      CALL POPINTEGER4(ad_to)
      DO j=ad_to,ad_from,-1
        CALL POPINTEGER4(ad_from0)
        CALL POPINTEGER4(ad_to0)
        DO i=ad_to0,ad_from0,-1
          CALL POPINTEGER4(branch)
          IF (branch .LT. 2) THEN
            nnxb = 0.0
            nnyb = 0.0
            nnzb = 0.0
          ELSE
            padj3b(ii, jj) = padj3b(ii, jj) + padj0b(ii, jj)
            padj0b(ii, jj) = 0.0
            DO l=nt2mg,nt1mg,-1
              wadj3b(ii, jj, l) = wadj3b(ii, jj, l) + wadj1b(ii, jj, l)
              wadj1b(ii, jj, l) = 0.0
            END DO
            wadj3b(ii, jj, irhoe) = wadj3b(ii, jj, irhoe) + wadj0b(ii, &
&              jj, irhoe)
            wadj0b(ii, jj, irhoe) = 0.0
            wadj3b(ii, jj, ivz) = wadj3b(ii, jj, ivz) + wadj0b(ii, jj, &
&              ivz)
            vnb = -(nnz*wadj0b(ii, jj, ivz))
            nnzb = -(vn*wadj0b(ii, jj, ivz))
            wadj0b(ii, jj, ivz) = 0.0
            wadj3b(ii, jj, ivy) = wadj3b(ii, jj, ivy) + wadj0b(ii, jj, &
&              ivy)
            vnb = vnb - nny*wadj0b(ii, jj, ivy)
            nnyb = -(vn*wadj0b(ii, jj, ivy))
            wadj0b(ii, jj, ivy) = 0.0
            wadj3b(ii, jj, ivx) = wadj3b(ii, jj, ivx) + wadj0b(ii, jj, &
&              ivx)
            vnb = vnb - nnx*wadj0b(ii, jj, ivx)
            tempb0 = two*vnb
            nnxb = wadj3(ii, jj, ivx)*tempb0 - vn*wadj0b(ii, jj, ivx)
            wadj0b(ii, jj, ivx) = 0.0
            wadj3b(ii, jj, irho) = wadj3b(ii, jj, irho) + wadj0b(ii, jj&
&              , irho)
            wadj0b(ii, jj, irho) = 0.0
            CALL POPREAL8(vn)
            wadj3b(ii, jj, ivx) = wadj3b(ii, jj, ivx) + nnx*tempb0
            wadj3b(ii, jj, ivy) = wadj3b(ii, jj, ivy) + nny*tempb0
            nnyb = nnyb + wadj3(ii, jj, ivy)*tempb0
            wadj3b(ii, jj, ivz) = wadj3b(ii, jj, ivz) + nnz*tempb0
            nnzb = nnzb + wadj3(ii, jj, ivz)*tempb0
          END IF
          padj2b(ii, jj) = padj2b(ii, jj) + padj1b(ii, jj)
          padj1b(ii, jj) = 0.0
          DO l=nt2mg,nt1mg,-1
            wadj2b(ii, jj, l) = wadj2b(ii, jj, l) + wadj1b(ii, jj, l)
            wadj1b(ii, jj, l) = 0.0
          END DO
          wadj2b(ii, jj, irhoe) = wadj2b(ii, jj, irhoe) + wadj1b(ii, jj&
&            , irhoe)
          wadj1b(ii, jj, irhoe) = 0.0
          wadj2b(ii, jj, ivz) = wadj2b(ii, jj, ivz) + wadj1b(ii, jj, ivz&
&            )
          vnb = -(nnz*wadj1b(ii, jj, ivz))
          nnzb = nnzb - vn*wadj1b(ii, jj, ivz)
          wadj1b(ii, jj, ivz) = 0.0
          wadj2b(ii, jj, ivy) = wadj2b(ii, jj, ivy) + wadj1b(ii, jj, ivy&
&            )
          vnb = vnb - nny*wadj1b(ii, jj, ivy)
          nnyb = nnyb - vn*wadj1b(ii, jj, ivy)
          wadj1b(ii, jj, ivy) = 0.0
          wadj2b(ii, jj, ivx) = wadj2b(ii, jj, ivx) + wadj1b(ii, jj, ivx&
&            )
          vnb = vnb - nnx*wadj1b(ii, jj, ivx)
          tempb = two*vnb
          nnxb = nnxb + wadj2(ii, jj, ivx)*tempb - vn*wadj1b(ii, jj, ivx&
&            )
          wadj1b(ii, jj, ivx) = 0.0
          wadj2b(ii, jj, irho) = wadj2b(ii, jj, irho) + wadj1b(ii, jj, &
&            irho)
          wadj1b(ii, jj, irho) = 0.0
          CALL POPREAL8(vn)
          wadj2b(ii, jj, ivx) = wadj2b(ii, jj, ivx) + nnx*tempb
          wadj2b(ii, jj, ivy) = wadj2b(ii, jj, ivy) + nny*tempb
          nnyb = nnyb + wadj2(ii, jj, ivy)*tempb
          wadj2b(ii, jj, ivz) = wadj2b(ii, jj, ivz) + nnz*tempb
          nnzb = nnzb + wadj2(ii, jj, ivz)*tempb
          CALL POPREAL8(nnz)
          normadjb(nn, ii, jj, 3) = normadjb(nn, ii, jj, 3) + nnzb
          CALL POPREAL8(nny)
          normadjb(nn, ii, jj, 2) = normadjb(nn, ii, jj, 2) + nnyb
          CALL POPREAL8(nnx)
          normadjb(nn, ii, jj, 1) = normadjb(nn, ii, jj, 1) + nnxb
          CALL POPINTEGER4(jj)
          CALL POPINTEGER4(ii)
        END DO
      END DO
      CALL POPREAL8ARRAY(wadj2, 5**2*nw)
      CALL POPREAL8ARRAY(wadj3, 5**2*nw)
      CALL POPBOOLEAN(secondhalo)
      CALL EXTRACTBCSTATESADJ_B(nn, wadj, wadjb, padj, padjb, wadj0, &
&                          wadj0b, wadj1, wadj1b, wadj2, wadj2b, wadj3, &
&                          wadj3b, padj0, padj0b, padj1, padj1b, padj2, &
&                          padj2b, padj3, padj3b, rlvadj, revadj, &
&                          rlvadj1, rlvadj2, revadj1, revadj2, ioffset, &
&                          joffset, koffset, icell, jcell, kcell, isbeg&
&                          , jsbeg, ksbeg, isend, jsend, ksend, ibbeg, &
&                          jbbeg, kbbeg, ibend, jbend, kbend, icbeg, &
&                          jcbeg, icend, jcend, secondhalo)
    END IF
    CALL POPINTEGER4(isbeg)
    CALL POPINTEGER4(jsbeg)
    CALL POPINTEGER4(ksbeg)
    CALL POPINTEGER4(isend)
    CALL POPINTEGER4(jsend)
    CALL POPINTEGER4(ksend)
    CALL POPINTEGER4(ibbeg)
    CALL POPINTEGER4(jbbeg)
    CALL POPINTEGER4(kbbeg)
    CALL POPINTEGER4(ibend)
    CALL POPINTEGER4(jbend)
    CALL POPINTEGER4(kbend)
  END DO
END SUBROUTINE BCSYMMADJ_B
