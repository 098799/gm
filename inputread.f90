module inputread
  use precision
  use commontypes
  implicit none

contains

  subroutine read_Input(System)
    implicit none
    type(SystemData) :: System
    open(10,file="input.in",status="old")
    read(10,*) System%g
    read(10,*) System%nbas
    read(10,*) System%norb
    read(10,*) System%nprim
    read(10,*) System%alpha
    close(10)
  end subroutine read_Input
end module inputread
