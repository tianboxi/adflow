module fanRegionData
  use constants

  type fanRegionType

     character(len=maxStringLen) :: famName
     integer(kind=intType) :: famID
     ! The block indexes of the cells included in this region
     integer(kind=intType), dimension(:, :), pointer :: cellIDs

     ! Switches for three main conponents of the model
     logical :: flowTurning
     logical :: flowBlockage
     logical :: loss

     ! The total number of cells included this proc has
     integer(kind=intType) :: nCellIDs

     ! Blade normal vectors for each cell
     real(kind=realType), dimension(:,:), pointer :: normals

     ! Omega is the fan rotational speed
     real(kind=realType) :: Omega 

     ! B is the number of blades in fan 
     real(kind=realType) :: B

     ! blokage factor
     real(kind=realType), dimension(:), pointer :: blockage 

     ! F is the cell bodyforce
     real(kind=realType), dimension(:,:), pointer :: F

     ! W is the relative flow direction
     real(kind=realType), dimension(:,:), pointer :: W
     real(kind=realType), dimension(:,:), pointer :: Wt

     ! Total volume of all cells in the fan region 
     real(kind=realType) :: volume

     integer(kind=intType), dimension(:), allocatable :: blkPtr
     
     ! Set the defaults for solution relaxation 
     real(kind=realType) :: relaxStart = -one
     real(kind=realType) :: relaxEnd = -one
  end type fanRegionType

  integer(kind=intType), parameter :: nFanRegionsMax=10
  type(fanRegionType), dimension(nFanRegionsMax), target :: fanRegions
  integer(kind=intTYpe) :: nFanRegions=0

end module fanRegionData
