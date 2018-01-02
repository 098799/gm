module f12
  use precision
  use commontypes
  use integrals
  implicit none

contains

  subroutine do_F12(System,SCF)
    implicit none
    type(SystemData) :: System
    type(SCFData)    :: SCF
    integer          :: lwork,liwork,info
    integer          :: i,j,k,l
    real(prec), allocatable :: matS_S(:,:),matS_T(:,:),matH_S(:,:),matH_T(:,:)
    real(prec), allocatable :: intf122table(:,:,:,:),intdf122table(:,:,:,:)
    real(prec), allocatable :: eigen(:)
    real(prec), allocatable :: work(:),total
    integer, allocatable    :: iwork(:)

    System%how_many = how_many_pairs(System%nprim)

    allocate(intf122table(0:50,0:50,0:50,0:50))
    allocate(intdf122table(0:50,0:50,0:50,0:50))
    call read_intf122(intf122table)
    call read_intdf122(intdf122table)

    allocate(matS_S(System%how_many,System%how_many))
    allocate(matS_T(System%how_many,System%how_many))
    allocate(matH_S(System%how_many,System%how_many))
    allocate(matH_T(System%how_many,System%how_many))

    call make_matS(System,matS_S,matS_T,intf122table)
    call make_matH(System,matH_S,matH_T,intf122table,intdf122table)

    print*,    norm2(matS_S-transpose(matS_S))
    print*,    norm2(matS_T-transpose(matS_T))
    print*,    norm2(matH_S-transpose(matH_S))
    print*,    norm2(matH_T-transpose(matH_T))
    print*, "---"

    do i = 1, 5
       do j = 1, 5
          print*, i,j,matH_S(i,j)
       end do
    end do

    allocate(eigen(System%how_many))
    allocate(work(2),iwork(2))
    call dsyevd('V','U',System%how_many,matS_S,System%how_many,eigen,work,-1,iwork,-1,info)
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),iwork(liwork))
    call dsyevd('V','U',System%how_many,matS_S,System%how_many,eigen,work,lwork,iwork,liwork,info)
    deallocate(work,iwork)

    do i = 1, System%how_many
       print*, eigen(i)
    end do
    print*, "---"

    call make_matS(System,matS_S,matS_T,intf122table)
    eigen = 0._prec

    allocate(work(2),iwork(2))
    call dsygvd(1,'V','U',System%how_many,matH_S,System%how_many,matS_S,System%how_many,eigen,work,-1,iwork,-1,info)
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),iwork(liwork))
    call dsygvd(1,'V','U',System%how_many,matH_S,System%how_many,matS_S,System%how_many,eigen,work,lwork,iwork,iwork,info)
    deallocate(work,iwork)

    do i = 1, System%how_many
       print*, eigen(i)
    end do
    print*, "---"

  end subroutine do_F12
end module f12
