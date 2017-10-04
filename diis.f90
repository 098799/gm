module diis
use precision
implicit none

private
public DIISData
public init_DIIS,restart_DIIS,free_DIIS,use_DIIS

type DIISData
integer :: LDA,LDE,LDN
integer :: nA,nE
integer :: max_n,act_n,act_i
integer :: lwork
real(prec),allocatable :: A(:,:),E(:,:)
real(prec),allocatable :: mat(:,:),mat_solve(:,:),vec_solve(:)
real(prec),allocatable :: work(:)
integer,allocatable    :: ipiv(:)
end type DIISData

contains

subroutine init_DIIS(DIIS,nA,nE,max_n)
implicit none
type(DIISData) :: DIIS
integer,intent(in) :: nA,nE
integer,intent(in) :: max_n
real(prec) :: dum
integer :: i,idum,info

DIIS%LDA = nA
DIIS%LDE = nE
DIIS%LDN = max_n

DIIS%nA = nA
DIIS%nE = nE

DIIS%max_n = max_n
DIIS%act_n = 0
DIIS%act_i = 0

DIIS%lwork = 0
do i=2,DIIS%LDN+1
   call dsysv('U',i,1,dum,DIIS%LDN+1,idum,dum,DIIS%LDN+1,dum,-1,info)
   DIIS%lwork = max(DIIS%lwork,int(dum))
enddo

if(DIIS%LDE<0) return

allocate(DIIS%A(DIIS%LDA,DIIS%LDN))
allocate(DIIS%E(DIIS%LDE,DIIS%LDN))

allocate(DIIS%mat(DIIS%LDN,DIIS%LDN))
allocate(DIIS%mat_solve(DIIS%LDN+1,DIIS%LDN+1))
allocate(DIIS%vec_solve(DIIS%LDN+1))

allocate(DIIS%work(DIIS%lwork))
allocate(DIIS%ipiv(DIIS%LDN+1))

! call mem_alloc(DIIS%A,DIIS%LDA,DIIS%LDN)
! call mem_alloc(DIIS%E,DIIS%LDE,DIIS%LDN)

! call mem_alloc(DIIS%mat,DIIS%LDN,DIIS%LDN)
! call mem_alloc(DIIS%mat_solve,DIIS%LDN+1,DIIS%LDN+1)
! call mem_alloc(DIIS%vec_solve,DIIS%LDN+1)

! call mem_alloc(DIIS%work,DIIS%lwork)
! call mem_alloc(DIIS%ipiv,DIIS%LDN+1)

end subroutine init_DIIS

subroutine restart_DIIS(DIIS,nA,nE,max_n)
implicit none
type(DIISData) :: DIIS
integer,intent(in) :: nA,nE
integer,intent(in) :: max_n
logical :: valid

valid = .true.

if(nA>DIIS%LDA) then
   write(*,'(a)') 'Size of the A matrix is too large in restart_DIIS!'
   valid = .false.
endif
if(nE>DIIS%LDE) then
   write(*,'(a)') 'Size of the E matrix is too large in restart_DIIS!'
   valid = .false.
endif
if(max_n>DIIS%LDN) then
   write(*,'(a)') 'Maximal number of vectors is too large in restart_DIIS!'
   valid = .false.
endif

if(.not.valid) stop

DIIS%nA = nA
DIIS%nE = nE

DIIS%max_n = max_n
DIIS%act_n = 0
DIIS%act_i = 0

end subroutine restart_DIIS

subroutine free_DIIS(DIIS)
implicit none
type(DIISData) :: DIIS

if(DIIS%LDE<0) return

deallocate(DIIS%ipiv,DIIS%work,DIIS%vec_solve,DIIS%mat_solve,DIIS%mat,DIIS%E,DIIS%A)

! call mem_dealloc(DIIS%ipiv)
! call mem_dealloc(DIIS%work)

! call mem_dealloc(DIIS%vec_solve)
! call mem_dealloc(DIIS%mat_solve)
! call mem_dealloc(DIIS%mat)

! call mem_dealloc(DIIS%E)
! call mem_dealloc(DIIS%A)

end subroutine free_DIIS

subroutine use_DIIS(DIIS,A,E)
implicit none
type(DIISData) :: DIIS
real(prec),intent(inout) :: A(*)
real(prec),intent(in) :: E(*)
integer :: i,info

if(DIIS%LDE<0) return

if(DIIS%act_n<DIIS%max_n) DIIS%act_n = DIIS%act_n + 1

DIIS%act_i = DIIS%act_i + 1
if(DIIS%act_i>DIIS%max_n) DIIS%act_i = 1

DIIS%A(1:DIIS%nA,DIIS%act_i) = A(1:DIIS%nA)
DIIS%E(1:DIIS%nE,DIIS%act_i) = E(1:DIIS%nE)

call dgemv('T',DIIS%nE,DIIS%act_n,&
     1._prec,DIIS%E,DIIS%LDE,DIIS%E(:,DIIS%act_i),1,&
     0._prec,DIIS%mat(:,DIIS%act_i),1)
do i=1,DIIS%act_n
   if(i/=DIIS%act_i) DIIS%mat(DIIS%act_i,i) = DIIS%mat(i,DIIS%act_i)
enddo

DIIS%mat_solve(1:DIIS%act_n,1:DIIS%act_n) = DIIS%mat(1:DIIS%act_n,1:DIIS%act_n)
DIIS%mat_solve(1:DIIS%act_n,DIIS%act_n+1) = -1
DIIS%mat_solve(DIIS%act_n+1,1:DIIS%act_n) = -1
DIIS%mat_solve(DIIS%act_n+1,DIIS%act_n+1) = 0

DIIS%vec_solve(1:DIIS%act_n) = 0
DIIS%vec_solve(DIIS%act_n+1) = -1

call dsysv('U',DIIS%act_n+1,1,&
     DIIS%mat_solve,DIIS%LDN+1,DIIS%ipiv,DIIS%vec_solve,DIIS%LDN+1,&
     DIIS%work,DIIS%lwork,info)
if(info/=0) then
   write(*,'(a)') 'Linear problem not solved in use_DIIS!'
   stop
endif

call dgemv('N',DIIS%nA,DIIS%act_n,&
     1._prec,DIIS%A,DIIS%LDA,DIIS%vec_solve,1,0._prec,A,1)

end subroutine use_DIIS

end module diis
