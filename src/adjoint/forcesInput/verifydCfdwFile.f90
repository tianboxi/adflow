!
!     ******************************************************************
!     *                                                                *
!     * File:          verifydCfdx.f90                                 *
!     * Author:        Andre C. Marta, C.A.(Sandy) Mader               *
!     * Starting date: 01-15-2007                                      *
!     * Last modified: 05-13-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine verifydCfdwfile(level)
!
!     ******************************************************************
!     *                                                                *
!     * Compute all entries in dIdx (partial) using the automatically  *
!     * differentiated routines generated by Tapenade and compare      *
!     * them to the finite-difference results using the modified       *
!     * force routine. This is only executed in debug mode.            *
!     *                                                                *
!     ******************************************************************
!
      use adjointpetsc        !djdw
      use adjointVars         !nCellsGlobal
      use blockPointers
      use cgnsGrid            ! cgnsDoms
      use communication       ! procHalo(currentLevel)%nProcSend, myID
      use inputPhysics        ! equations
      use flowVarRefState     ! nw
      use inputDiscretization ! spaceDiscr, useCompactDiss
      use iteration           ! overset, currentLevel
      use inputTimeSpectral   ! nTimeIntervalsSpectral
      use section
      use monitor             ! monLoc, MonGlob, nMonSum
      use bcTypes             !imin,imax,jmin,jmax,kmin,kmax
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level
!
!     Local variables.
!
      integer(kind=intType) :: discr, nHalo, sps
      integer(kind=intType) :: icell, jcell, kcell, mm, nn, m, n
      integer(kind=intType) :: ii, jj, kk, i1, j1, k1, i2, j2, k2

      integer(kind=intType) ::  i2Beg,  i2End,  j2Beg,  j2End
      integer(kind=intType) :: iiBeg, iiEnd, jjBeg, jjEnd
      integer(kind=intType) :: i,j,k,l

      logical :: fineGrid,correctForK, exchangeTurb

      real(kind=realType) :: Cl,Cd,Cfx,Cfy,Cfz,Cmx,Cmy,Cmz
      real(kind=realType) :: ClAdj,CdAdj,CfxAdj,CfyAdj,CfzAdj,&
                             &CmxAdj,CmyAdj,CmzAdj  
      real(kind=realType) :: CLAdjP,CLAdjM,CDAdjP,CDAdjM,dCLdwFD,dCDdwFD,&
                             &dCmxdwFD,CmxAdjP,CmxAdjM
      real(kind=realType) :: ClAdjB,CdAdjB,CfxAdjB,CfyAdjB,CfzAdjB,&
                             &CmxAdjB,CmyAdjB,CmzAdjB 

      real(kind=realType) :: CLP,CLM,CDP,CDM,CmxP,CmxM,CmyP,CmyM,CmzP,CmzM,&
                             &CfxP,CfxM,CfyP,CfyM,CfzP,CfzM

      real(kind=realType), dimension(:,:,:,:), allocatable :: xAdj,xAdjB
      real(kind=realType), dimension(:,:,:,:), allocatable :: wAdj,wAdjB,wadjb2
      real(kind=realType), dimension(:,:,:), allocatable :: pAdj

      real(kind=realType), dimension(3) :: cFpAdj, cFvAdj
      real(kind=realType), dimension(3) :: cMpAdj, cMvAdj

      real(kind=realType), dimension(3) :: cFp, cFv
      real(kind=realType), dimension(3) :: cMp, cMv



      real(kind=realType) :: alphaAdj, betaAdj,MachAdj,machCoefAdj,MachGridAdj
      real(kind=realType) :: alphaAdjb, betaAdjb,MachAdjb,machCoefAdjb,MachGridAdjb
      REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
      REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
      REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
      REAL(KIND=REALTYPE) :: murefAdj, timerefAdj
      integer(kind=intType)::liftIndex
      real(kind=realType), dimension(3) ::rotRateAdj,rotCenterAdj,rotrateadjb
      
      real(kind=realType) :: factI, factJ, factK, tmp

      integer(kind=intType), dimension(0:nProc-1) :: offsetRecv

      real(kind=realType), dimension(4) :: time
      real(kind=realType)               :: timeAdj, timeFD

      ! > derivative output

      real(kind=realType), parameter :: deltaw = 1.e-8_realType

      real(kind=realType) :: wAdjRef, wref,test

      real(kind=realType), dimension(:,:,:), pointer :: norm
      real(kind=realType), dimension(:,:,:),allocatable:: normAdj
      real(kind=realType), dimension(3) :: refPoint
      real(kind=realType) :: yplusMax

      real(kind=realType),  dimension(:), allocatable :: monLoc1, monGlob1
      real(kind=realType),  dimension(:), allocatable :: monLoc2, monGlob2

      logical :: contributeToForce, viscousSubface,secondHalo,righthanded

      integer :: ierr,nmonsum1,nmonsum2,idxmgb

      character(len=2*maxStringLen) :: errorMessage

 

!File Parameters
      integer :: unitWcl = 8,unitWcd = 9,unitWcm = 10,ierror
      character(len = 10)::outfile
      
      outfile = "CSWcl.txt"
      
      open (UNIT=unitWcl,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydCfdwFile", &
           "Something wrong when &
           &calling open")

      outfile = "CSWcd.txt"
      
      open (UNIT=unitWcd,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydCfdwFile", &
           "Something wrong when &
           &calling open")

      outfile = "CSWcm.txt"
      
      open (UNIT=unitWcm,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydCfdwFile", &
           "Something wrong when &
           &calling open")
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      print *,'in verifydcfdw'
 !     if( myID==0 ) write(*,*) "Running verifydCfdx...",sps

      ! Set the grid level of the current MG cycle, the value of the
      ! discretization and the logical correctForK.

      currentLevel = level
      discr        = spaceDiscr
      correctForK  = .false.
      fineGrid     = .true.

      ! Determine whether or not the total energy must be corrected
      ! for the presence of the turbulent kinetic energy and whether
      ! or not turbulence variables should be exchanged.

      correctForK  = .false.
      exchangeTurb = .false.
      secondhalo = .true.

     
      nmonsum1 = 8
      nmonsum2 = 8

      allocate(monLoc1(nmonsum1), monGlob1(nmonsum1))
      allocate(monLoc2(nmonsum2), monGlob2(nmonsum2))

 

     ! Exchange the pressure if the pressure must be exchanged early.
      ! Only the first halo's are needed, thus whalo1 is called.
      ! Only on the fine grid.
      
      if(exchangePressureEarly .and. currentLevel <= groundLevel) &
           call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
           .false., .false.)
      
      ! Apply all boundary conditions to all blocks on this level.
      
      call applyAllBC(secondHalo)
      
      ! Exchange the solution. Either whalo1 or whalo2
      ! must be called.
      
      if( secondHalo ) then
         call whalo2(currentLevel, 1_intType, nMGVar, .true., &
              .true., .true.)
      else
         call whalo1(currentLevel, 1_intType, nMGVar, .true., &
              .true., .true.)
      endif

      call mpi_barrier(SUmb_comm_world, ierr)      


      print *,"halo's updated"
!
!     ******************************************************************
!     *                                                                *
!     * Compute the d(forces)/dx (partial) using the tapenade routines.*
!     *                                                                *
!     ******************************************************************
         
!*********************
      !from ForcesAndMoments.f90
       ! Determine the reference point for the moment computation in
       ! meters.
      print *,'setting refpoint'
       refPoint(1) = LRef*pointRef(1)
       refPoint(2) = LRef*pointRef(2)
       refPoint(3) = LRef*pointRef(3)

       ! Initialize the force and moment coefficients to 0 as well as
       ! yplusMax.

!!$       cFpAdj(1) = zero; cFpAdj(2) = zero; cFpAdj(3) = zero
!!$       cFvAdj(1) = zero; cFvAdj(2) = zero; cFvAdj(3) = zero
!!$       cMpAdj(1) = zero; cMpAdj(2) = zero; cMpAdj(3) = zero
!!$       cMvAdj(1) = zero; cMvAdj(2) = zero; cMvAdj(3) = zero
       print *,'zeroing output'
       ClAdj = zero
       CDAdj = zero
       CmxAdj = zero 
       CmyAdj = zero
       CmzAdj = zero
       CfxAdj = zero
       CfyAdj = zero
       CfzAdj = zero

       yplusMax = zero

   !    print *,'Adjoint forces initialized'
!
!***********************************
       print *,' computing adjoint derivatives'
       
   spectralLoopAdj: do sps=1,nTimeIntervalsSpectral
      
      print *,'zeroing force output'
      !zero the force components
      cFpAdj(1) = zero; cFpAdj(2) = zero; cFpAdj(3) = zero
      cFvAdj(1) = zero; cFvAdj(2) = zero; cFvAdj(3) = zero
      cMpAdj(1) = zero; cMpAdj(2) = zero; cMpAdj(3) = zero
      cMvAdj(1) = zero; cMvAdj(2) = zero; cMvAdj(3) = zero
      ! Loop over the number of local blocks.
      print *,'starting domain loop'
      domainLoopAD: do nn=1,nDom

        ! Set some pointers to make the code more readable.
         !print *,'setting pointers'
        call setPointersAdj(nn,level,sps)
        !print *,'allocating memory'
        !print *,'globalcell',globalcell
        !stop
        allocate(xAdj(0:ie,0:je,0:ke,3), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for xAdj.")

        allocate(xAdjB(0:ie,0:je,0:ke,3), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for xAdjB.")
        
        allocate(wAdj(0:ib,0:jb,0:kb,nw), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for wAdj.")

        allocate(wAdjB(0:ib,0:jb,0:kb,nw), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for wAdjB.")
        allocate(wAdjB2(0:ib,0:jb,0:kb,nw), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for wAdjB2.")
        
        allocate(pAdj(0:ib,0:jb,0:kb), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for pAdj.")
        
        !print *,'finished allocating',nn,level,sps
        righthanded = flowDoms(nn,level,sps)%righthanded
        
        !print *,'copying stencil'
        ! Copy the coordinates into xAdj and
        ! Compute the face normals on the subfaces
        call copyADjointForcesStencil(wAdj,xAdj,alphaAdj,betaAdj,&
           MachAdj,machCoefAdj,machGridAdj,prefAdj,rhorefAdj, pinfdimAdj,&
           rhoinfdimAdj,rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,murefAdj,&
           timerefAdj,pInfCorrAdj,nn,level,sps,liftIndex)
        wAdjB(:,:,:,:) = zero 
        wAdjB2(:,:,:,:) = zero 
        bocoLoop: do mm=1,nBocos

           ! Determine the range of cell indices of the owned cells
           ! Notice these are not the node indices
           iiBeg = BCData(mm)%icBeg
           iiEnd = BCData(mm)%icEnd
           jjBeg = BCData(mm)%jcBeg
           jjEnd = BCData(mm)%jcEnd
           
           i2Beg= BCData(mm)%inBeg+1; i2End = BCData(mm)%inEnd
           j2Beg= BCData(mm)%jnBeg+1; j2End = BCData(mm)%jnEnd

          
           ! Initialize the seed for reverse mode. Cl is the first one
           ClAdjB = 1
           CDAdjB = 0
           CmxAdjB = 0
           CmyAdjB = 0
           CmzAdjB = 0
!!$        CfxAdjB = 0
!!$        CfyAdjB = 0
!!$        CfzAdjB = 0
           
           !xAdjB(:,:,:,:) = zero ! > return dCf/dx
           !wAdjB(:,:,:,:) = zero ! > return dCf/dW
           !print *,'calling adjoint forces_b'
           !===============================================================
           !           
           !print *,'Initial Parameters Calculated,Computing Lift Partials...'
           
           !===============================================================
           ! Compute the force derivatives
           call  COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
&  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
&  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
&  cmpadj, righthanded, secondhalo, alphaadj, alphaadjb, betaadj, &
&  betaadjb, machadj, machadjb, machcoefadj, machcoefadjb, machgridadj, &
&  machgridadjb, prefadj, rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj&
&  , pinfadj, murefadj, timerefadj, pinfcorradj, rotcenteradj, &
&  rotrateadj, rotrateadjb, liftindex)

           wadjb2 = wadjb2+wadjb


        enddo bocoLoop
        do k = 0,kb
           do j = 0,jb
              do i = 0,ib
                 do l = 1,nw
                    !idxmgb = globalCell(i,j,k)
                    idxmgb   = globalCell(i,j,k)*nw+l
                    !test = wadjb(i,j,k,l)
                    
                    if( idxmgb>=0) then
                       !if ( test.ne.0 .and. idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<nCellsGlobal) then
                       !print *,'globalcell',idxmgb,wadjb(i,j,k,l),wadjb2(i,j,k,l)
                       !idxmgb = globalCell(i,j,k)*nw+l-1  !L minus 1 not 1 minus 1
                       write(unitwcl,10) wadjb2(i,j,k,l),nn,i,j,k,l,idxmgb
10                     format(1x,'wcl ',f18.10,7I8)
                    endif
                 enddo
              enddo
           enddo
        enddo

        xAdjB(:,:,:,:) = zero ! > return dCf/dx
        wAdjB(:,:,:,:) = zero ! > return dCf/dW
        wAdjB2(:,:,:,:) = zero
        bocoLoop2: do mm=1,nBocos
           
           ! Determine the range of cell indices of the owned cells
           ! Notice these are not the node indices
           iiBeg = BCData(mm)%icBeg
           iiEnd = BCData(mm)%icEnd
           jjBeg = BCData(mm)%jcBeg
           jjEnd = BCData(mm)%jcEnd
           
           i2Beg= BCData(mm)%inBeg+1; i2End = BCData(mm)%inEnd
           j2Beg= BCData(mm)%jnBeg+1; j2End = BCData(mm)%jnEnd
           ! Initialize the seed for reverse mode. Cd second
           ClAdjB = 0
           CDAdjB = 1
           CmxAdjB = 0
           CmyAdjB = 0
           CmzAdjB = 0
!!$        CfxAdjB = 0
!!$        CfyAdjB = 0
!!$        CfzAdjB = 0
       
           !xAdjB(:,:,:,:) = zero ! > return dCf/dx
           !wAdjB(:,:,:,:) = zero ! > return dCf/dW

           !===============================================================
           !           
!           print *,'Calculating Drag Partials...'
           
           !===============================================================
           ! Compute the force derivatives
            
           call  COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
&  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
&  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
&  cmpadj, righthanded, secondhalo, alphaadj, alphaadjb, betaadj, &
&  betaadjb, machadj, machadjb, machcoefadj, machcoefadjb, machgridadj, &
&  machgridadjb, prefadj, rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj&
&  , pinfadj, murefadj, timerefadj, pinfcorradj, rotcenteradj, &
&  rotrateadj, rotrateadjb, liftindex)
           wadjb2 = wadjb2+wadjb  
        enddo bocoLoop2
        do k = 0,kb
           do j = 0,jb
              do i = 0,ib
                 do l = 1,nw
                    !idxmgb = globalCell(i,j,k)
                    idxmgb   = globalCell(i,j,k)*nw+l
                    !test = wadjb(i,j,k,l)
                    
                    if( idxmgb>=0) then
                       !if ( test.ne.0 .and. idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<nCellsGlobal) then
                       
                       !idxmgb = globalCell(i,j,k)*nw+l-1  !L minus 1 not 1 minus 1
                       write(unitwcd,11) wadjb2(i,j,k,l),nn,i,j,k,l,idxmgb
11                     format(1x,'wcd ',f18.10,7I8)
                    endif
                 enddo
              enddo
           enddo
        enddo
        xAdjB(:,:,:,:) = zero ! > return dCf/dx
        wAdjB(:,:,:,:) = zero ! > return dCf/dW
        wAdjB2(:,:,:,:) = zero ! > return dCf/dW
        bocoLoop3: do mm=1,nBocos

           ! Determine the range of cell indices of the owned cells
           ! Notice these are not the node indices
           iiBeg = BCData(mm)%icBeg
           iiEnd = BCData(mm)%icEnd
           jjBeg = BCData(mm)%jcBeg
           jjEnd = BCData(mm)%jcEnd
           
           i2Beg= BCData(mm)%inBeg+1; i2End = BCData(mm)%inEnd
           j2Beg= BCData(mm)%jnBeg+1; j2End = BCData(mm)%jnEnd
           ! Initialize the seed for reverse mode. Cmx third
           ClAdjB = 0
           CDAdjB = 0
           CmxAdjB = 1
           CmyAdjB = 0
           CmzAdjB = 0
 !       CfxAdjB = 0
 !       CfyAdjB = 0
 !       CfzAdjB = 0
       
           !xAdjB(:,:,:,:) = zero ! > return dCf/dx
           !wAdjB(:,:,:,:) = zero ! > return dCf/dW


           !===============================================================
           !           
!           print *,'Calculating MomX Partials...'
           
           !===============================================================
           ! Compute the force derivatives
           call  COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
&  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
&  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
&  cmpadj, righthanded, secondhalo, alphaadj, alphaadjb, betaadj, &
&  betaadjb, machadj, machadjb, machcoefadj, machcoefadjb, machgridadj, &
&  machgridadjb, prefadj, rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj&
&  , pinfadj, murefadj, timerefadj, pinfcorradj, rotcenteradj, &
&  rotrateadj, rotrateadjb, liftindex)
           wadjb2 = wadjb2+wadjb
           
        enddo bocoLoop3
                do k = 0,kb
           do j = 0,jb
              do i = 0,ib
                 do l = 1,nw
                    !idxmgb = globalCell(i,j,k)
                    idxmgb   = globalCell(i,j,k)*nw+l
                    !test = wadjb(i,j,k,l)
                    
                    if( idxmgb>=0) then
                       !if ( test.ne.0 .and. idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<nCellsGlobal) then
                       
                       !idxmgb = globalCell(i,j,k)*nw+l-1  !L minus 1 not 1 minus 1
                       write(unitwcm,12) wadjb2(i,j,k,l),nn,i,j,k,l,idxmgb
12                     format(1x,'wcm ',f18.10,7I8)
                    endif
                 enddo
              enddo
           enddo
        enddo
        xAdjB(:,:,:,:) = zero ! > return dCf/dx
        wAdjB(:,:,:,:) = zero ! > return dCf/dW
        !===============================================================
        
        !print *,' deallocating'
        ! Deallocate the xAdj.
        deallocate(pAdj, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for xAdj.")
             
        ! Deallocate the xAdj.
        deallocate(wAdj, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for wAdj.") 

         deallocate(wAdjB, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for wAdjb.") 
        deallocate(wAdjB2, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for wAdjb.") 
        ! Deallocate the xAdj.
        deallocate(xAdj, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for xAdj.") 

         deallocate(xAdjB, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for xAdjb.") 
        !print *,'finishhed deallocating'
       
      enddo domainLoopAD

   enddo spectralLoopAdj
      
   print *,'AD loop finished'
   !stop
   ! Get new time and compute the elapsed AD time.
   
   call mpi_barrier(SUmb_comm_world, ierr)
   if(myID == 0) then
      call cpu_time(time(2))
      timeAdj = time(2)-time(1)
   endif
   
   
   ! Flush the output buffer and synchronize the processors.
   
   call f77flush()
   call mpi_barrier(SUmb_comm_world, ierr)
!
!     ******************************************************************
!

      
      ! Deallocate memory for the temporary arrays.


  
      ! Output formats.

!  10  format(1x,a,1x,i3,1x,i3,1x,a,1x,i3,1x,i3,1x,i3)           
  20  format(1x,(e18.6),2x,(e18.6),2x,(e18.6))
  30  format(1x,a,1x,i3,2x,e13.6,1x,5(i2,1x),3x,e13.6,1x,5(i2,1x))
  99  format(a,1x,i6)
    end subroutine verifydCfdwfile
    
