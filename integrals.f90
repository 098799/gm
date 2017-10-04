module integrals
  use precision
  use commontypes
  implicit none
  integer, parameter, private :: tho=300,ghn=100,hn=200
  real(prec), parameter       :: sqrtpi=1.77245385090551602729816748334114518_prec
  real(prec), private         :: norm(0:tho),ghabsciss(ghn),ghweights(ghn),hermiteh_gh(0:hn,ghn)


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
       do i = 1, ghn
          intgh = intgh + hermiteh_gh(m,i)*hermiteh_gh(o,i)*hermiteh_gh(p,i)*hermiteh_gh(r,i)*ghweights(i)
       end do
       intgh = intgh*norm(m)*norm(o)*norm(p)*norm(r)/sqrt(beta)
    end if
  end function intgh

  subroutine gener_hermiteh_gh()
    implicit none
    integer         :: i,j
    real(prec)      :: beta
    beta = 2._prec
    hermiteh_gh=0._prec
    do i = 1, ghn
       hermiteh_gh(0,i)=1._prec
       hermiteh_gh(1,i)=2._prec*ghabsciss(i)/sqrt(beta)
    end do
    do j=2,hn
       do i = 1, ghn
          hermiteh_gh(j,i) = 2._prec*ghabsciss(i)/sqrt(beta)*hermiteh_gh(j-1,i)-2._prec*(j-1)*hermiteh_gh(j-2,i)
       end do
    end do
  end subroutine gener_hermiteh_gh

  subroutine import_gh
    implicit none
    integer            :: i
    open(10,file='dat/ghabsciss.dat')
    ghabsciss=0._prec
    do i = 1, ghn
       read(10,*) ghabsciss(i)
    end do
    close(10)

    open(10,file='dat/ghweights.dat')
    ghweights=0._prec
    do i = 1, ghn
       read(10,*) ghweights(i)
    end do
    close(10)
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
