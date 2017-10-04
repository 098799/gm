module commontypes
  use precision
  implicit none
  type SystemData
     real(prec) :: g
     integer    :: nbas
     integer    :: norb
     integer    :: nprim
     real(prec) :: alpha
  end type SystemData
end module commontypes
