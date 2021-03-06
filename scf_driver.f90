module scf_driver
  use precision
  use commontypes
  use integrals
  use diis
  implicit none

contains

  subroutine do_SCF(System,SCF)
    implicit none
    type(SystemData)        :: System
    type(SCFData)           :: SCF
    type(DIISData)          :: DIIS
    integer                 :: lwork,liwork,info,iter
    integer                 :: i
    integer, allocatable    :: iwork(:)
    real(prec), allocatable :: matH(:,:),matC(:,:),matD(:,:),matF(:,:),matP(:,:)
    real(prec), allocatable :: twoel(:,:,:,:),moint(:,:,:,:),eval(:),work(:)
    real(prec)              :: ener,err
    allocate(matH(System%nbas,System%nbas),matC(System%nbas,System%nbas))
    allocate(matD(System%nbas,System%nbas),matF(System%nbas,System%nbas))
    allocate(matP(System%nbas,System%nbas))
    allocate(twoel(System%nbas,System%nbas,System%nbas,System%nbas))
    allocate(moint(System%nbas,System%nbas,System%nbas,System%nbas))
    allocate(eval(System%nbas))

    allocate(work(2),iwork(2))
    call dsyevd('V','U',System%nbas,matC,System%nbas,eval,work,-1,iwork,-1,info)
    lwork=int(work(1))
    liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),iwork(liwork))

    call SCF_matH(System,matH)
    call create_SCF_file(System)
    call read_SCF_file(System,twoel)

    matC(:,:) = matH
    call dsyevd('V','U',System%nbas,matC,System%nbas,eval,work,lwork,iwork,liwork,info)

    call init_DIIS(DIIS,System%nbas**2,system%nbas**2,6)
    write(*,*) "*** SCF iterations                                    ***"
    write(*,*) "       Iter   Energy                   Error"
    iter = 0
    do
       iter = iter + 1
       call dgemm('N','T',System%nbas,System%nbas,System%nocc,&
            1._prec,matC,System%nbas,matC,System%nbas,0._prec,matD,System%nbas)
       matF(:,:) = matH
       call dgemv('N',System%nbas**2,System%nbas**2,&
            1._prec,twoel,System%nbas**2,matD,1,1._prec,matF,1)
       call dgemm('N','N',System%nbas,System%nbas,System%nbas,&
            1._prec,matF,System%nbas,matD,System%nbas,0._prec,matP,System%nbas)
       call dgemm('N','N',System%nbas,System%nbas,System%nbas,&
            1._prec,matD,System%nbas,matF,System%nbas,-1._prec,matP,System%nbas)
       ener=sum((matH+matF)*matD)
       err=norm2(matP)
       write(*,*) iter,ener,err
       if (err.LT.1e-12) exit
       call use_DIIS(DIIS,matF,matP)
       matC(:,:) = matF
       call dsyevd('V','U',System%nbas,matC,System%nbas,eval,work,lwork,iwork,liwork,info)
    enddo
    ! write(*,*) "*** End of SCF iterations                             ***"
    call free_DIIS(DIIS)

    allocate(SCF%eval(System%nbas),SCF%evec(System%nbas,System%nbas))
    SCF%eval(:)   = eval
    SCF%evec(:,:) = matC
    SCF%ener      = ener
    write(*,*) "*** SCF Energy                                        ***"
    write(*,*) SCF%ener

    call four_index_transform(System,SCF,twoel,moint)

    call mp2(System,SCF,moint)

  end subroutine do_SCF

  subroutine mp2(System,SCF,moint)
    implicit none
    type(SystemData)        :: System
    type(SCFData)           :: SCF
    real(prec)              :: moint(:,:,:,:)
    real(prec)              :: emp2
    integer                 :: a,b,i,j
    emp2 = 0._prec
    do a = System%nocc+1, System%nbas
       do b = System%nocc+1, System%nbas
          do i = 1, System%nocc
             do j = 1, System%nocc
                emp2 = emp2 + (2._prec*moint(a,i,b,j)**2-moint(a,i,b,j)*moint(b,i,a,j))&
                     /(SCF%eval(i)+SCF%eval(j)-SCF%eval(a)-SCF%eval(b))
             end do
          end do
       end do
    end do
    write(*,*) "*** MP2 Energy,                MP2+SCF                ***"
    print*, "  ",emp2,emp2+SCF%ener
  end subroutine mp2

  subroutine four_index_transform(System,SCF,twoel,moint)
    implicit none
    type(SystemData)        :: System
    type(SCFData)           :: SCF
    integer                 :: k,l
    real(prec)              :: twoel(:,:,:,:),moint(:,:,:,:)
    real(prec),allocatable  :: tmp1(:,:),tmp2(:,:)
    write(*,*) "*** 4-index-transform commencing!!!                   ***"
    allocate(tmp1(System%nbas,System%nbas),tmp2(System%nbas,System%nbas))
    do l=1,System%nbas
       do k=1,System%nbas
          call dgemm('T','N',System%nbas,System%nbas,System%nbas,&
               1._prec,SCF%evec,System%nbas,twoel(:,:,k,l),System%nbas,0._prec,tmp1,System%nbas)
          call dgemm('N','N',System%nbas,System%nbas,System%nbas,&
               1._prec,tmp1,System%nbas,SCF%evec,System%nbas,0._prec,tmp2,System%nbas)
          moint(k,l,:,:) = tmp2
       enddo
    enddo
    do l=1,System%nbas
       do k=1,System%nbas
          call dgemm('T','N',System%nbas,System%nbas,System%nbas,&
               1._prec,SCF%evec,System%nbas,moint(:,:,k,l),System%nbas,0._prec,tmp1,System%nbas)
          call dgemm('N','N',System%nbas,System%nbas,System%nbas,&
               1._prec,tmp1,System%nbas,SCF%evec,System%nbas,0._prec,moint(:,:,k,l),System%nbas)
       enddo
    enddo

    ! block
    !   integer :: i,j,k,l,a,b,c,d
    !   real(prec),allocatable :: in1(:,:,:,:)
    !   real(prec) :: temp
    !   allocate(in1(System%nbas,System%nbas,System%nbas,System%nbas))
    !   do l = 1, System%nbas
    !      do k = 1, System%nbas
    !         do j = 1, System%nbas
    !            do a = 1, System%nbas
    !               temp = 0._prec
    !               do i = 1, System%nbas
    !                  temp = temp + SCF%evec(i,a)*twoel(i,j,k,l)
    !               end do
    !               in1(a,j,k,l) = temp
    !            end do
    !         end do
    !      end do
    !   end do
    !   do l = 1, System%nbas
    !      do k = 1, System%nbas
    !         do a = 1, System%nbas
    !            do b = 1, System%nbas
    !               temp = 0._prec
    !               do j = 1, System%nbas
    !                  temp = temp + SCF%evec(j,b)*in1(a,j,k,l)
    !               end do
    !               moint(a,b,k,l) = temp
    !            end do
    !         end do
    !      end do
    !   end do
    !   do l = 1, System%nbas
    !      do a = 1, System%nbas
    !         do b = 1, System%nbas
    !            do c = 1, System%nbas
    !               temp = 0._prec
    !               do k = 1, System%nbas
    !                  temp = temp + SCF%evec(k,c)*moint(a,b,k,l)
    !               end do
    !               in1(a,b,c,l) = temp
    !            end do
    !         end do
    !      end do
    !   end do
    !   do a = 1, System%nbas
    !      do b = 1, System%nbas
    !         do c = 1, System%nbas
    !            do d = 1, System%nbas
    !               temp = 0._prec
    !               do l = 1, System%nbas
    !                  temp = temp + SCF%evec(l,d)*in1(a,b,c,l)
    !               end do
    !               moint(a,b,c,d) = temp
    !            end do
    !         end do
    !      end do
    !   end do
    ! end block
  end subroutine four_index_transform

end module scf_driver
