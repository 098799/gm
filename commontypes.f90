module commontypes
  use precision
  implicit none
  type SystemData
     real(prec) :: g
     integer    :: nbas
     integer    :: nocc
     integer    :: nprim
     real(prec) :: alpha
  end type SystemData
  type SCFData
     real(prec), allocatable :: eval(:)
     real(prec), allocatable :: evec(:,:)
     real(prec)              :: ener
     real(prec), allocatable :: pair_energy(:,:,:)
  end type SCFData
end module commontypes
