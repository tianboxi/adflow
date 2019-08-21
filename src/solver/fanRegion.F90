module fanRegion

  use constants
  use communication, only : commType, internalCommType
  use fanRegionData
  implicit none

contains
  subroutine addFanRegion(pts, conn, donorData, omega, B, famName, famID, &
      relaxStart, relaxEnd, nPts, nConn)

    use communication, only : myID, adflow_comm_world
    use constants
    use adtBuild, only : buildSerialHex, destroySerialHex
    use adtLocalSearch, only : containmentTreeSearchSinglePoint
    use ADTUtils, only : stack
    use ADTData
    use blockPointers, only : x, il, jl, kl, nDom, iBlank, vol, volRef, si, sj, sk
    use adjointVars, only : nCellsLocal
    use utils, only : setPointers, EChk
    implicit none

    ! Input variables
    real(kind=realType), dimension(3, nPts), intent(in), target :: pts
    integer(kind=intType), dimension(8, nConn), intent(in), target :: conn
    real(kind=realType), dimension(4, nPts), intent(in), target :: donorData 
    integer(kind=intType), intent(in) :: nPts, nConn, famID
    real(kind=realType), intent(in) :: omega, B, relaxStart, relaxEnd
    character(len=*) :: famName

    ! Working variables
    integer(kind=intType) :: i, j, k, nn, iDim, cellID, intInfo(3), sps, level, nData, ierr, iii
    real(kind=realType) :: volLocal
    type(fanRegionType), pointer :: region
    type(adtType) :: ADT
    integer(kind=intType), dimension(:,:), pointer :: tmp
    real(kind=realType), dimension(:,:), pointer :: tmpr
    real(kind=realType), dimension(:), pointer :: tmp2
    logical :: failed

    ! ADT Type required data
    integer(kind=intType), dimension(:), pointer :: frontLeaves, frontLeavesNew
    integer(kind=intType), dimension(:), pointer :: BB
    real(kind=realType) :: coor(3), uvw(7)

    nFanRegions = nFanRegions + 1
    
    ! Set pointers and save fan region info
    region => fanRegions(nFanRegions)
    region%famName = famName
    region%famID = famID
    region%Omega = omega
    region%B = B
    region%relaxStart = relaxStart
    region%relaxEnd = relaxEnd

    allocate(region%blkPtr(0:nDom))
    region%blkPtr(0) = 0
    region%nCellIDs = 0

    ! Allocate the extra data the tree search requires.
    allocate(stack(100), BB(100), frontLeaves(100), frontLeavesNew(100))

    ! Allocate sufficient space for the maximum possible number of cellIDs
    allocate(region%cellIDs(3, nCellsLocal(1)))
    allocate(region%normals(3, nCellsLocal(1)))
    allocate(region%blockage(nCellsLocal(1)))

    ! Build the ADT tree.
    call buildSerialHex(nConn, nPts, pts, conn, ADT)

    ! Only work for single grid now
    sps = 1
    level = 1
    ! Three components of normal vector and the blockage factor need to be interpolated
    nData = 4 
    failed = .True.
    ! Loop for searching
    do nn=1, nDom
       call setPointers(nn, level, sps)
       do k=2, kl
          do j=2, jl
             do i=2, il
                ! Only check real cells
                if (iblank(i,j,k) == 1) then
                   ! Compute the cell center
                   coor = eighth*(x(i-1,j-1,k-1,:) + x(i,j-1,k-1,:)  &
                        +         x(i-1,j,  k-1,:) + x(i,j,  k-1,:)  &
                        +         x(i-1,j-1,k,  :) + x(i,j-1,k,  :)  &
                        +         x(i-1,j,  k,  :) + x(i,j,  k,  :))
                   intInfo(3) = 0
                   uvw(1) = -1.0_realType
                   uvw(2) = -1.0_realType
                   uvw(3) = -1.0_realType
                   call containmentTreeSearchSinglePoint(ADT, coor, intInfo, &
                        uvw, donorData, nData, BB, frontLeaves, frontLeavesNew, failed)
                   if(uvw(1) > zero .and. uvw(1) < one .and. &
                      uvw(2) > zero .and. uvw(2) < one .and. &
                      uvw(3) > zero .and. uvw(3) < one .and. .not. failed) then
                      ! Add cell to the list
                      region%nCellIDs = region%nCellIDs + 1
                      region%cellIDs(:, region%nCellIDs) = (/i, j, k/)
                      region%normals(:, region%nCellIDs) = uvw(4:6)
                      region%blockage(region%nCellIDs) = uvw(7)
                   end if 
                end if
             end do
          end do
       end do
       region%blkPtr(nn) = region%nCellIDs
    end do

    
    ! Resize the cellIDs and normals to the correct size now that we know the
    ! correct exact number.
    tmp => region%cellIDs
    allocate(region%cellIDs(3, region%nCellIDs))
    region%cellIDs = tmp(:, 1:region%nCellIDs)
    deallocate(tmp)
    tmpr => region%normals
    allocate(region%normals(3, region%nCellIDs))
    region%normals = tmpr(:, 1:region%nCellIDs)
    deallocate(tmpr)
    tmp2 => region%blockage
    allocate(region%blockage(region%nCellIDs))
    region%blockage = tmp2(1:region%nCellIDs)
    deallocate(tmp2)

    allocate(region%F(3, region%nCellIDs))
    allocate(region%W(3, region%nCellIDs))
    allocate(region%Wt(3, region%nCellIDs))
 
    ! Perform corrections for the cell vols and face area projections
    ! in i and j direction (due to blockage)
    volLocal = zero

    do nn=1, nDom
       call setPointers(nn, level, sps)

       ! Loop over the region for this block
       do iii=region%blkPtr(nn-1) + 1, region%blkPtr(nn)
          i = region%cellIDs(1, iii)
          j = region%cellIDs(2, iii)
          k = region%cellIDs(3, iii)
          vol(i,j,k) = vol(i,j,k) * region%blockage(iii)
          volRef(i,j,k) = volRef(i,j,k) * region%blockage(iii)
          si(i,j,k,:) = si(i,j,k,:) * region%blockage(iii)
          sj(i,j,k,:) = sj(i,j,k,:) * region%blockage(iii)
          volLocal = volLocal + vol(i, j, k)
       end do
    end do

    call mpi_allreduce(volLocal, region%volume, 1, adflow_real, &
         MPI_SUM, adflow_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)


    ! Final memory cleanup
    deallocate(frontLeaves, frontLeavesNew, BB)
    call destroySerialHex(ADT)

  end subroutine addFanRegion

  subroutine writeFanRegions(fileName)

    ! This a (mostly) debug routine that is used to verify to the user
    ! the that the cells that the user thinks should be specified as
    ! being inside the actuator region actually are. We will dump a
    ! hex unstructured ascii tecplot file with all the zones we
    ! found. We won't be super concerned about efficiency here.

    use constants
    use utils, only : EChk, pointReduce, setPointers
    use communication, only : myID, adflow_comm_world, nProc
    use blockPointers, only : x, nDom
    implicit none

    ! Input
    character(len=*) :: fileName

    ! Working
    integer(kind=intType) :: iRegion, nn, i, j, k, ii, jj, kk, iii, kkk, iDim
    integer(kind=intType) :: level, sps, iProc, ierr, totalCount, offset, nUnique
    integer(kind=intType), dimension(:), allocatable :: sizesProc, cumSizesProc
    real(kind=realType) , dimension(:), allocatable :: pts, allPts, blockage
    real(kind=realType) , dimension(:,:), allocatable :: cellNormals
    real(kind=realType) , dimension(:,:), allocatable :: tmp, uniquePts
    real(kind=realType), parameter :: tol=1e-8
    integer(kind=intType), dimension(:), allocatable :: conn, allConn, link
    character(80) :: zoneName
    type(fanRegionType), pointer :: region

    ! Before we start the main region loop the root procesoor has to
    ! open up the tecplot file and write the header

    if (myid == 0) then
       open(unit=101, file=trim(fileName), form='formatted')
       write(101,*)'# vtk DataFile Version 3.0'
       write(101,*)'Fan regions'
       write(101,*)'ASCII'
       write(101,*)'DATASET UNSTRUCTURED_GRID'
    end if

    ! Region Loop
    regionLoop: do iRegion=1, nFanRegions

       ! Only for the finest grid level.
       level =1
       sps = 1

       ! Do an allgather with the number of actuator cells on each
       ! processor so that everyone knows the sizes and can compute the offsets,
       region => fanRegions(iRegion)

       allocate(sizesProc(nProc), cumSizesProc(0:nProc))

       call mpi_allgather(region%nCellIDs, 1, adflow_integer, sizesProc, 1, &
            adflow_integer, adflow_comm_world, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       cumSizesProc(0) = 0
       do iProc=1, nProc
          cumSizesProc(iProc) = cumSizesProc(iProc-1) + sizesProc(iProc)
       end do

       ! Fill up our own nodes/conn with the nodes we have here.
       allocate(conn(8*region%nCellIDs), pts(24*region%nCellIDs), cellNormals(3,region%nCellIDs))
       allocate(blockage(region%nCellIDs))

       kkk = 0
       do nn=1, nDom
          call setPointers(nn, level, sps)
          ! Loop over the ranges for this block
          do iii=region%blkPtr(nn-1) + 1, region%blkPtr(nn)

             ! Carful with the conn values! They need to be in counter clock wise ordering!
             offset = (iii-1)*8 + cumSizesProc(myID)*8
             conn((iii-1)*8 + 1) = 1 + offset
             conn((iii-1)*8 + 2) = 2 + offset
             conn((iii-1)*8 + 3) = 4 + offset
             conn((iii-1)*8 + 4) = 3 + offset
             conn((iii-1)*8 + 5) = 5 + offset
             conn((iii-1)*8 + 6) = 6 + offset
             conn((iii-1)*8 + 7) = 8 + offset
             conn((iii-1)*8 + 8) = 7 + offset

             ! Add in the 24 values for the nodal coordinates in coordinate
             ! ordering. Do all the coordinates interlaced
             do kk=-1, 0
                do jj=-1, 0
                   do ii=-1, 0
                      do iDim=1,3
                         i = region%cellIDs(1, iii)
                         j = region%cellIDs(2, iii)
                         k = region%cellIDs(3, iii)
                         kkk = kkk + 1
                         pts(kkk) = x(i+ii, j+jj, k+kk, iDim)
                      end do
                   end do
                end do
             end do
          end do
       end do
       

       ! Now that we've filled up our array, we can allocate the total
       ! space we need on the root proc and it
       if (myid == 0) then
          totalCount = sum(sizesProc)
          allocate(allConn(8*totalCount), allPts(24*totalCount))
       end if

       ! Perform the two gatherV's
       call mpi_gatherV(pts, region%nCellIDs*24, adflow_real, &
            allPts, 24*sizesProc, 24*cumSizesProc, adflow_real, &
            0, adflow_comm_world, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       call mpi_gatherV(conn, region%nCellIDs*8, adflow_integer, &
            allConn, 8*sizesProc, 8*cumSizesProc, adflow_integer, &
            0, adflow_comm_world, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       ! We can deallocate all the per-proc memory now
       deallocate(sizesProc, cumSizesProc, pts, conn)

       if (myid == 0) then

          ! Now the poor root processor dumps everything out to a
          ! file. To help cut down on the already bloated file size,
          ! we'll point reduce it which will help tecplot display them
          ! better as well.

          allocate(tmp(3, totalCount*8))
          do i=1, totalCount*8
             do iDim=1, 3
                tmp(iDim, i) = allPts((i-1)*3 + iDim)
             end do
          end do
          deallocate(allPts)
          allocate(uniquePts(3, totalCount*8), link(totalCount*8))

          ! Get unique set of nodes.
          call pointReduce(tmp, totalCount*8, tol, uniquePts, link, nUnique)
          write(101,*) "POINTS", nUnique, "double"

13        format (F20.12)

          ! Write all the coordinates...this is horrendously slow...

          do i=1, nUnique
             do iDim=1, 3
                write(101,13,advance='no') uniquePts(iDim, i)
             end do
             write(101,"(1x)")
          end do

          ! Write out the connectivity
          write(101, *) "CELLS", totalcount, totalcount*9
15        format(I8)
          do i=1, totalCount
             write(101, 15, advance='no') 8
             do j=1,8
                write(101, 15, advance='no') link(allConn((i-1)*8 + j))-1
             end do
             write(101,"(1x)")
          end do
          write(101, *) "CELL_TYPES", totalcount

          do i=1, totalcount
             write(101, *) 12
          end do

          write(101,*) "CELL_DATA", totalcount
          write(101,*) "VECTORS normals double"
          do i=1, totalcount
             do iDim =1,3
                write(101,13,advance='no') region%normals(iDim, i)
             end do
             write(101,"(1x)")
          end do
          
          write(101,*) "VECTORS forces double"
          do i=1, totalcount
             do iDim =1,3
                write(101,13,advance='no') region%F(iDim, i)
             end do
             write(101,"(1x)")
          end do

          write(101,*) "VECTORS W double"
          do i=1, totalcount
             do iDim =1,3
                write(101,13,advance='no') region%W(iDim, i)
             end do
             write(101,"(1x)")
          end do
          write(101,*) "VECTORS Wt double"
          do i=1, totalcount
             do iDim =1,3
                write(101,13,advance='no') region%Wt(iDim, i)
             end do
             write(101,"(1x)")
          end do

          write(101,*) "SCALARS blockage double 1"
          write(101,*) "LOOKUP_TABLE default"
          do i=1, totalcount
              write(101,13) region%blockage(i)
          end do

          ! Ditch the memory only allocated on this proc
          deallocate(allConn, link, tmp, uniquePts)
       end if
    end do regionLoop

    ! Close the output file on the root proc
    if (myid == 0) then
       close(101)
    end if
  end subroutine writeFanRegions

end module fanRegion

