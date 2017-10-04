program gm
  use precision
  use commontypes
  use inputread
  use diis
  use scf_driver
  use f12
  implicit none
  type(SystemData) :: System
  type(SCFData)    :: SCF

  call read_Input(System)
  call make_norm()
  call import_gh()
  call gener_hermiteh_gh()

  call do_SCF(System,SCF)
  call do_F12(System)

end program gm
