module scf
  use precision
  use commontypes
  use integrals
  implicit none

contains

  subroutine do_SCF(System)
    implicit none
    type(SystemData)        :: System
    real(prec), allocatable :: matH(:,:)
    allocate(matH(System%nbas,System%nbas))

    call SCF_matH(System,matH)
    print*, intgh(1,2,3,4)


    deallocate(matH)
  end subroutine do_SCF

end module scf
