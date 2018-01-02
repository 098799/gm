module commontypes
  use precision
  implicit none
  type SystemData
     real(prec) :: g
     integer    :: nbas
     integer    :: nocc
     integer    :: nprim
     real(prec) :: alpha
     integer    :: how_many
  end type SystemData
  type SCFData
     real(prec), allocatable :: eval(:)
     real(prec), allocatable :: evec(:,:)
     real(prec)              :: ener
     real(prec), allocatable :: pair_energy(:,:,:)
  end type SCFData

contains

  integer function how_many_pairs(nprim)
    implicit none
    integer             :: i,j,counter
    integer, intent(in) :: nprim
    i = -1
    j = -1
    how_many_pairs = 0
    do while (next_pair(i,j,nprim))
       how_many_pairs = how_many_pairs + 1
    end do
  end function how_many_pairs

  logical function next_pair(i,j,nprim)
    implicit none
    integer :: i,j,nprim
    next_pair = .TRUE.
    if ((i.LT.0).OR.(j.LT.0)) then
       i = 0
       j = 0
       return
    end if
    i = i+1
    if ((i.GT.j).OR.(i+j.GT.nprim)) then
       i = 0
       j = j+1
       if (j.GT.nprim) then
          next_pair = .FALSE.
          i = -1
          j = -1
          return
       end if
    end if
  end function next_pair
end module commontypes
