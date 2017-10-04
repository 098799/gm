program gm
  use precision
  use commontypes
  use inputread
  use scf
  use f12
  implicit none
  type(SystemData) :: System

  call read_Input(System)
  call import_gh()
  call gener_hermiteh_gh()
  call make_norm()
  call do_SCF(System)
  call do_F12(System)

end program gm
