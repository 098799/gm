module integrals
  use precision
  use commontypes
  implicit none
  integer, parameter, private :: tho=400,points=400,order=400,RECS=12
  real(prec), parameter       :: sqrtpi=1.77245385090551602729816748334114518_prec
  real(prec), private         :: norm(0:tho),ghabsciss(order),ghweights(order),hermiteh_gh(0:order,points)
  real(prec), private         :: factorials(0:tho),binomials(0:tho,0:tho),dzejmu(4,-5:tho)

contains

  subroutine make_matH(System,matH_S,matH_T,intf122table,intdf122table)
    implicit none
    integer    :: i,j,k,l,iprim,kprim
    type(SystemData) :: System
    real(prec) :: matH_S(:,:),matH_T(:,:),val1,val2
    real(prec) :: intdf122table(0:,0:,0:,0:),intf122table(0:,0:,0:,0:)
    matH_S = 0._prec
    matH_T = 0._prec
    i = -1
    j = -1
    k = -1
    l = -1
    kprim = 0
    do while (next_pair(k,l,System%nprim))
       kprim = kprim + 1
       iprim = 0
       do while (next_pair(i,j,System%nprim))
          iprim = iprim + 1
          val1=(intdf122table(i,k,j,l)+0.5_prec*(i+j+k+l+2)*intf122table(i,k,j,l))
          val2=(intdf122table(i,l,j,k)+0.5_prec*(i+j+k+l+2)*intf122table(i,l,j,k))
          matH_S(iprim,kprim) = val1+val2
          matH_T(iprim,kprim) = val1-val2
       end do
    end do
  end subroutine make_matH

  subroutine make_matS(System,matS_S,matS_T,intf122table)
    implicit none
    integer    :: i,j,k,l,iprim,kprim
    type(SystemData) :: System
    real(prec) :: matS_S(:,:),matS_T(:,:),val1,val2
    real(prec) :: intf122table(0:,0:,0:,0:)
    matS_S = 0._prec
    matS_T = 0._prec
    i = -1
    j = -1
    k = -1
    l = -1
    kprim = 0
    do while (next_pair(k,l,System%nprim))
       kprim = kprim + 1
       iprim = 0
       do while (next_pair(i,j,System%nprim))
          iprim = iprim + 1
          val1 = intf122table(i,k,j,l)
          val2 = intf122table(i,l,j,k)
          matS_S(iprim,kprim) = val1+val2
          matS_T(iprim,kprim) = val1-val2
       end do
    end do
  end subroutine make_matS

  subroutine read_intf122(intf122table)
    implicit none
    integer    :: i,j,k,l,ii,length,number
    real(prec) :: val
    real(prec) :: intf122table(0:,0:,0:,0:)
    open(10,file='dat/file_f122_50.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    inquire(10,size=length)
    intf122table = 0._prec
    do ii = 1, length/RECS
       read(10,rec=ii) number,val
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       intf122table(i,j,k,l) = val
       intf122table(j,i,k,l) = val
       intf122table(i,j,l,k) = val
       intf122table(j,i,l,k) = val
    end do
    close(10)
  end subroutine read_intf122

  subroutine read_intdf122(intdf122table)
    implicit none
    integer    :: i,j,k,l,ii,length,number
    real(prec) :: val
    real(prec) :: intdf122table(0:,0:,0:,0:)
    open(10,file='dat/file_df122_50.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    inquire(10,size=length)
    intdf122table = 0._prec
    do ii = 1, length/RECS
       read(10,rec=ii) number,val
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       intdf122table(i,j,k,l) = val
       intdf122table(j,i,k,l) = val
       intdf122table(i,j,l,k) = val
       intdf122table(j,i,l,k) = val
    end do
    close(10)
  end subroutine read_intdf122

  real(prec) function intdf122(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: tempf12,beta
    intdf122 = 0._prec
    if (modulo(mm+oo+pp+rr,2).EQ.0) then
       if (pp.lt.mm) then
          m=pp
          p=mm
       else
          m=mm
          p=pp
       end if
       if (rr.lt.oo) then
          r=oo
          o=rr
       else
          r=rr
          o=oo
       end if
       do s = 0, m
          tempf12 = 0._prec
          do t = 0, o
             tempf12 = tempf12 + coeff(t,o,r) * imunudf122norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          intdf122 = intdf122 + coeff(s,m,p) * tempf12/norm(m+p-2*s)
       end do
       intdf122 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intdf122
    end if
  end function intdf122

  real(prec) function imunudf122norm(m,n,beta)
    integer    :: m,n,b
    real(prec) :: beta,mn
    !some super bogus code...
    if (m.EQ.0.and.n.eq.0) then
       imunudf122norm = 0.5390125124744174978149239120563448922870409268698824918081870965_prec
    else if (m.EQ.2.and.n.eq.0.or.m.EQ.0.and.n.eq.2) then
       imunudf122norm = -0.1345197891935502887834707784128512214525037706050587983590004814_prec
    else if (m.EQ.1.and.n.eq.1) then
       imunudf122norm = 0.1902397102850885286405613807257687855130732683070173500499483870_prec
    else
       !/bogus
       mn = real(m+n,prec)
       imunudf122norm = 0._prec
       if (mod(m+n,2).NE.1) then
          b = bet(beta)
          imunudf122norm = norm(n)*norm(m)/norm(m+n)*&
               & 2._prec*sqrtpi*(-1._prec)**n*2._prec**(-(m+n)/2._prec)*&
               &(dzejmu(b,m+n)-&
               & 2._prec*(beta-1._prec)/4._prec*(dzejmu(b,m+n+2)*abn(m+n,2,.TRUE.)+&
               & (4._prec*(m+n)+2._prec)*dzejmu(b,m+n)+&
               & (4._prec*mn**2-4._prec*mn)*dzejmu(b,m+n-2)*abn(m+n,2,.FALSE.))+&
               & ((beta-1._prec)/4._prec)**2*(dzejmu(b,m+n+4)*abn(m+n,4,.TRUE.)+&
               & (8._prec*mn+12._prec)*dzejmu(b,m+n+2)*abn(m+n,2,.TRUE.)+&
               & (24._prec*mn**2+24._prec*mn+12._prec)*dzejmu(b,m+n)+&
               & (32._prec*mn**3-48._prec*mn**2+16._prec*mn)*dzejmu(b,m+n-2)*abn(m+n,2,.FALSE.)+&
               & (16._prec*mn*(mn-1._prec)*(mn-2._prec)*(mn-3._prec))*dzejmu(b,m+n-4)*abn(m+n,4,.FALSE.)))
       end if
    end if
  end function imunudf122norm

  real(prec) function abn(m,n,pm)
    implicit none
    integer    :: m,n,i
    ! real(prec) :: abn
    logical    :: pm
    abn = 1._prec
    if (n.LE.m) then
       if (pm) then
          do i = 1, n
             abn = abn * real(m+i,prec)
          end do
          abn = sqrt(abn) * sqrt(2._prec**real(n,prec))
       else
          do i = 1, n
             abn = abn / real(m-i+1,prec)
          end do
          abn = sqrt(abn) / sqrt(2._prec**real(n,prec))
       end if
    end if
  end function abn

  real(prec) function intf122(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: tempf12,beta
    intf122 = 0._prec
    if (modulo(mm+oo+pp+rr,2).EQ.0) then
       if (pp.lt.mm) then
          m=pp
          p=mm
       else
          m=mm
          p=pp
       end if
       if (rr.lt.oo) then
          r=oo
          o=rr
       else
          r=rr
          o=oo
       end if
       do s = 0, m
          tempf12 = 0._prec
          do t = 0, o
             tempf12 = tempf12 + coeff(t,o,r) * imunuf122norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          intf122 = intf122 + coeff(s,m,p) * tempf12/norm(m+p-2*s)
       end do
       intf122 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intf122
    end if
  end function intf122

  real(prec) function imunuf122norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: beta
    imunuf122norm = 0._prec
    if (mod(mu+nu,2).NE.1) then
       b = bet(beta)
       imunuf122norm = norm(nu)*norm(mu)/norm(mu+nu)*4._prec*sqrtpi*(-1._prec)**nu*2._prec**(-(mu+nu)/2._prec)&
            &*(0.25_prec*dzejmu(b,mu+nu+2)*sqrt(4._prec*real(mu+nu+2)*real(mu+nu+1))&
            &+(real(mu+nu,prec)+0.5_prec)*dzejmu(b,mu+nu)&
            &+sqrt(real(mu+nu,prec)*real(mu+nu-1,prec))*dzejmu(b,mu+nu-2)/2._prec)
    end if
  end function imunuf122norm

  real(prec) function coeff(s,i,k)
    implicit none
    integer :: i,k,s
    coeff = factorials(s)*2._prec**real(s,prec)*binomials(i,s)*binomials(k,s)
  end function coeff

  subroutine import_factorials
    implicit none
    integer            :: i
    open(10,file='dat/factorials.dat')
    factorials=0._prec
    do i = 0, tho
       read(10,*) factorials(i)
    end do
    close(10)
  end subroutine import_factorials

  subroutine import_binomials
    implicit none
    integer :: i,j
    open(10,file='dat/binomials.dat')
    binomials=0._prec
    do i = 0, tho
       do j = 0, i
          read(10,*) binomials(i,j)
       end do
    end do
    close(10)
  end subroutine import_binomials


  subroutine import_dzejmu()
    implicit none
    integer    :: mu,b
    open(10,file='dat/dz_2.5.dat')
    b = bet(2.5_prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    open(10,file='dat/dz_3.dat')
    b = bet(3._prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    open(10,file='dat/dz_4.dat')
    b = bet(4._prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    open(10,file='dat/dz_5.dat')
    b = bet(5._prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
  end subroutine import_dzejmu

  integer function bet(beta)
    real(prec) :: beta
    if (abs(beta-3._prec).LT.1e-8) then
       bet = 1
    else if (abs(beta-5._prec).LT.1e-8) then
       bet = 2
    else if (abs(beta-2.5_prec).LT.1e-8) then
       bet = 3
    else if (abs(beta-4._prec).LT.1e-8) then
       bet = 4
    else
       print*, "Unsupported beta. Bombing out."
       call exit()
    end if
  end function bet

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
