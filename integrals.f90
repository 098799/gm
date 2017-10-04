module integrals
  use precision
  use commontypes
  implicit none
  integer, parameter, private :: tho=400,points=400,order=400,RECS=12
  real(prec), parameter       :: sqrtpi=1.77245385090551602729816748334114518_prec
  real(prec), private         :: norm(0:tho),ghabsciss(order),ghweights(order),hermiteh_gh(0:order,points)

contains

  subroutine SCF_matH(System,matH)
    implicit none
    type(SystemData) :: System
    integer          :: i
    real(prec)       :: matH(:,:)
    matH = 0._prec
    do i = 1, System%nbas
       matH(i,i) = (i-1)+0.5_prec
    end do
  end subroutine SCF_matH

  subroutine read_SCF_file(System,twoel)
    implicit none
    type(SystemData) :: System
    integer :: length
    integer :: number,i,j,k,l,ii
    real(prec) :: val
    real(prec) :: twoel(:,:,:,:)
    open(10,file='bas/file2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    inquire(10,size=length)
    twoel = 0._prec
    do ii = 1, length/RECS
       read(10,rec=ii) number,val
       i=iand(number,255)+1
       j=iand(ishft(number,-8),255)+1
       k=iand(ishft(number,-16),255)+1
       l=ishft(number,-24)+1
       val = val * System%g
       twoel(i,j,k,l) = val
       twoel(i,j,l,k) = val
       twoel(i,k,j,l) = val
       twoel(i,k,l,j) = val
       twoel(i,l,j,k) = val
       twoel(i,l,k,j) = val
       twoel(j,i,k,l) = val
       twoel(j,i,l,k) = val
       twoel(j,k,i,l) = val
       twoel(j,k,l,i) = val
       twoel(j,l,i,k) = val
       twoel(j,l,k,i) = val
       twoel(k,i,j,l) = val
       twoel(k,i,l,j) = val
       twoel(k,j,i,l) = val
       twoel(k,j,l,i) = val
       twoel(k,l,i,j) = val
       twoel(k,l,j,i) = val
       twoel(l,i,j,k) = val
       twoel(l,i,k,j) = val
       twoel(l,j,i,k) = val
       twoel(l,j,k,i) = val
       twoel(l,k,i,j) = val
       twoel(l,k,j,i) = val
    end do
    close(10)
  end subroutine read_SCF_file

  subroutine create_SCF_file(System)
    implicit none
    integer          :: i,j,k,l,counter
    type(SystemData) :: System
    open(10,file='bas/file2.F',status='replace',form='unformatted',access='direct',RECL=RECS)
    counter = 0
    do l = 0, System%nbas-1
       do k = 0, l
          do j = 0, k
             do i = 0, j
                if ((modulo(i+j+k+l,2).EQ.0)) then
                   counter = counter + 1
                   write(10,rec=counter) i+ishft(j,8)+ishft(k,16)+ishft(l,24),intgh(i,j,k,l)
                end if
             end do
          end do
       end do
    end do
    close(10)
  end subroutine create_SCF_file

  real(prec) function intgh(mm,oo,pp,rr)
    implicit none
    integer    :: m,o,p,r,i
    integer    :: mm,oo,pp,rr,vec(4),vect(4)
    real(prec) :: beta
    beta  = 2._prec
    intgh = 0._prec
    if (modulo(mm+oo+pp+rr,2).EQ.0) then
       vect(1) = mm
       vect(2) = oo
       vect(3) = pp
       vect(4) = rr
       vec = sortt(vect)
       m = vec(1)
       o = vec(2)
       p = vec(3)
       r = vec(4)
       do i = 1, points
          intgh = intgh + hermiteh_gh(m,i)*hermiteh_gh(o,i)*hermiteh_gh(p,i)*hermiteh_gh(r,i)*ghweights(i)
       end do
       intgh = intgh/sqrt(beta)
    end if
  end function intgh

  subroutine gener_hermiteh_gh()
    implicit none
    integer         :: i,j
    real(prec)      :: beta
    hermiteh_gh=0._prec
    do i = 1, points
       hermiteh_gh(0,i)=1._prec
       hermiteh_gh(1,i)=2._prec*ghabsciss(i)
    end do
    do j = 2, order
       do i = 1, points
          hermiteh_gh(j,i) = 2._prec*ghabsciss(i)*hermiteh_gh(j-1,i)-2._prec*(j-1)*hermiteh_gh(j-2,i)
       end do
    end do
    do j = 0, order
       do i = 1, points
          hermiteh_gh(j,i) = hermiteh_gh(j,i)*norm(j)
       end do
    end do
  end subroutine gener_hermiteh_gh

  subroutine import_gh
    implicit none
    integer            :: i
    open(10,file='dat/rtssq.dat',status='unknown',form='formatted')
    ghabsciss = 0._prec
    do i = 1, order
       read(10,*) ghabsciss(i)
    end do
    close(10)

    open(11,file='dat/weights.dat',status='unknown',form='formatted')
    ghweights = 0._prec
    do i = 1, order
       read(11,*) ghweights(i)
    end do
    close(11)
  end subroutine import_gh

  function sortt(vec)
    implicit none
    integer               :: n,j,temp,bubble
    integer, dimension(:) :: vec
    integer, dimension(4) :: sortt
    n = 4
    do while (n.GT.1)
       bubble = 0
       do j = 1, (n-1)
          if (vec(j)>vec(j+1)) then
             temp = vec(j)
             vec(j) = vec(j+1)
             vec(j+1) = temp
             bubble = j
          end if
       end do
       n = bubble
    end do
    sortt = vec
  end function sortt

  subroutine make_norm()
    implicit none
    integer :: m
    norm(0)=1._prec/sqrt(sqrtpi)
    do m = 1, tho
       norm(m) = norm(m-1)/sqrt(2._prec*m)
    end do
  end subroutine make_norm
end module integrals
