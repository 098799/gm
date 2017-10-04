module f12
  use precision
  use commontypes
  implicit none

contains

  subroutine do_F12(System)
    implicit none
    type(SystemData) :: System
    integer          :: i,j
    i=-1
    j=-1
    print*, System%nprim
    ! do while (next_pair(i,j,System%nprim))
    !    print*, i,j
    ! end do
  end subroutine do_F12

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
end module f12
