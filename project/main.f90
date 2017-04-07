! program to solve the heat equation (Dirichlet boundary conditions only)
PROGRAM main

! turn off implicit typing
implicit none

! define double precision with 15 digits of accuracy
! and an exponent range of +- 307
integer, parameter :: rk = selected_real_kind(15, 307)

integer             :: AllocateStatus ! variable to hold memory allocation success
integer             :: i, j, q        ! loop iteration variables

integer  :: n_el               ! number of elements
integer  :: n_en               ! number of nodes per element
integer  :: n_nodes            ! total number of nodes
integer  :: order              ! polynomial order
integer  :: n_qp               ! number of quadrature points
real(rk) :: length             ! length of the domain (1-D)
real(rk) :: h                  ! length of one element
real(rk) :: k                  ! thermal conductivity
real(rk) :: source             ! uniform heat source
real(rk) :: theta              ! conjugate gradient coefficient
real(rk) :: lambda             ! conjugate gradient coefficient
real(rk) :: convergence        ! difference between successive iterations
real(rk) :: tol                ! conjugate gradient convergence tolerance

integer,  dimension(:, :), allocatable :: LM     ! location matrix
real(rk), dimension(:, :), allocatable :: qp     ! quadrature points and weights
real(rk), dimension(:),    allocatable :: x      ! coordinates of the nodes
real(rk), dimension(:, :), allocatable :: kel    ! elemental stiffness matrix
real(rk), dimension(:),    allocatable :: rel    ! elemental load vector
real(rk), dimension(:, :), allocatable :: phi    ! shape functions
real(rk), dimension(:, :), allocatable :: dphi   ! shape function derivatives
real(rk), dimension(:, :), allocatable :: kglob  ! global stiffness matrix
real(rk), dimension(:),    allocatable :: rglob  ! global load vector
real(rk), dimension(:),    allocatable :: a      ! conjugate gradient solution iterates
real(rk), dimension(:),    allocatable :: aprev  ! conjugate gradient solution iterates
real(rk), dimension(:),    allocatable :: z      ! conjugate gradient update iterates
real(rk), dimension(:),    allocatable :: zprev  ! conjugate gradient update iterates
real(rk), dimension(:),    allocatable :: res    ! solution residual

! initialize the thermal conductivity and heat source
k = 1.0
source = 0.0
tol = 0.001

call commandline(length, n_el, order)              ! parse command line arguments
call initialize(h, x, n_en, n_el, order, n_nodes)  ! initialize problem variables
call quadrature(order, n_qp, qp)                   ! initialize quadrature rule
call phi_val(order, qp, phi, dphi)                 ! initialize shape functions

! form the elemental stiffness matrix and load vector
allocate(kel(n_en, n_en), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of kel array failed."
allocate(rel(n_en), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of rel array failed."

kel = 0.0
rel = 0.0
! only works for 1-D elements
do q = 1, n_qp
  do i = 1, n_en
    rel(i) = rel(i) + qp(q, 2) * source * phi(q, i) * h / 2.0
    do j = 1, n_en
      kel(i, j) = kel(i, j) + qp(q, 2) * dphi(q, i) * k * dphi(q, j) * 2.0 / h
    end do
  end do
end do

! form the location matrix
call locationmatrix(LM, n_el, n_en)

! form the global stiffness matrix and global load vector
allocate(kglob(n_nodes, n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of kglob array failed."
allocate(rglob(n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of rglob array failed."

kglob = 0.0
rglob = 0.0
do q = 1, n_el
  do i = 1, n_en
    rglob(LM(q, i)) = rglob(LM(q, i)) + rel(i)
    do j = 1, n_en
      kglob(LM(q, i), LM(q, j)) = kglob(LM(q, i), LM(q, j)) + kel(i, j)
    end do
  end do
end do


allocate(a(n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of a array failed."
allocate(aprev(n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of aprev array failed."
allocate(z(n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of z array failed."
allocate(zprev(n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of zprev array failed."

allocate(res(n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of res array failed."


kglob(1, :) = 0.0
kglob(1, 1) = 1.0                ! set up left BC              

kglob(n_nodes, :) = 0.0
kglob(n_nodes, n_nodes) = 1.0    ! set up right BC

rglob(1) = 2.0                   ! left BC value
rglob(n_nodes) = 1.0             ! right BC value     


!print *, 'solving matrix system: kglob: ', kglob(:, :)
!print *, 'with load vector: ', rglob(:)

! conjugate gradient solver for solving kglob * a = rglob
aprev  = 0.0
print *, 'first guess: ', aprev(:)
res    = rglob - matmul(kglob, aprev)
print *, 'matrix-vector product: ', matmul(kglob, aprev)
zprev  = res

! if we hit it on the head, then don't do any more iterations
! check to make sure the residual is not already zero



lambda = dot_product(zprev, res) / dot_product(zprev, matmul(kglob, zprev))
a      = aprev + lambda * zprev
print *, 'first residual: ', res, '\n'
print *, 'lambda: ', lambda
print *, 'a: ', a

convergence = 0.0
do i = 1, n_nodes
  convergence = convergence + abs(a(i) - aprev(i))
end do

do while (convergence > tol)
  aprev  = a
  res    = rglob - matmul(kglob, a)
  theta  = - dot_product(res, matmul(kglob, zprev)) / dot_product(zprev, matmul(kglob, zprev))
  z      = res + theta * zprev
  lambda = dot_product(z, res) / dot_product(z, matmul(kglob, z))
  a      = a + lambda * z
  zprev  = z
  
  print *, 'residual: ', res(:)
  !print *, 'theta: ', theta
  !print *, 'z: ', z(:)
  !print *, 'lambda: ', lambda
  !print *, 'a: ', a

  convergence = 0.0
  do i = 1, n_nodes
    convergence = convergence + abs(a(i) - aprev(i))
  end do
  
  print *, 'error: ', convergence
end do

! write to an output file for plotting. If this file exists, it will be re-written.
open(1, file='output.txt', iostat=AllocateStatus, status="replace")
if (AllocateStatus /= 0) STOP "output.txt file opening failed."

write(1, *) a(:)



! write program variables to user
!print *, 'Quadrature points:'
!print *, qp
!print *, 'Phi values:', phi(:, :)
!print *, 'Derivative values:', dphi(:, :)
!print *, 'Elemental stiffness matrix:', kel(:, :)
!print *, 'Elemental load vector:', rel(:)
!print *, 'Global vector: ', rglob(:) 

!print *, 'Global matrix:', kglob(:, :)
!do i = 1, n_nodes
!  print *, kglob(:, i)
!end do

! deallocate memory 
deallocate(qp, x, kel, rel, phi, dphi, kglob, rglob, a, aprev, z, zprev, res)



CONTAINS

subroutine locationmatrix(LM, n_el, n_en)
  implicit none
  integer, dimension(:, :), allocatable, intent(out) :: LM
  integer, intent(in) :: n_el
  integer, intent(in) :: n_en
  integer :: i, j

  allocate(LM(n_el, n_en), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of LM array failed."
  
  do i = 1, n_el
    do j = 1, n_en
LM(i, j) = 1 + (i - 1) * (n_en - 1) + (j - 1)
    end do
  end do
end subroutine locationmatrix

subroutine phi_val(order, qp, phi, dphi)
  implicit none
  integer, intent(in)  :: order
  real(rk), dimension(:, :), intent(in)  :: qp
  real(rk), dimension(:, :), allocatable, intent(inout) :: phi
  real(rk), dimension(:, :), allocatable, intent(inout) :: dphi

  ! allocate memory for the shape functions - quadrature points
  ! do not change throughout the simulation
  allocate(phi(n_qp, order + 1), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of phi array failed."

  ! allocate memory for the shape functions - quadrature points
  ! do not change throughout the simulation
  allocate(dphi(n_qp, order + 1), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of dphi array failed."

  select case(order)
    case(1)
      phi(:, 1) = (1.0 - qp(:, 1)) / 2.0
      phi(:, 2) = (1.0 + qp(:, 1)) / 2.0
      dphi(:, 1) = -1.0 / 2.0
      dphi(:, 2) =  1.0 / 2.0
    case(2)
      phi(:, 1) = qp(:, 1) * (qp(:, 1) - 1.0) / 2.0
      phi(:, 2) = (1.0 - qp(:, 1)) * (1.0 + qp(:, 1))
      phi(:, 3) = qp(:, 1) * (qp(:, 1) + 1.0) / 2.0
      dphi(:, 1) = (2.0 * qp(:, 1) - 1.0) / 2.0
      dphi(:, 2) = 1.0 - qp(:, 1) * qp(:, 1)
      dphi(:, 3) = (2.0 * qp(:, 1) + 1.0) / 2.0
    case default
      write(*,*) "phi and dphi not initialized due to polynomial order not being supported."
  end select
end subroutine phi_val


subroutine quadrature(order, n_qp, qp)
  implicit none

  integer, intent(in)  :: order
  integer, intent(out) :: n_qp
  real (rk), dimension(:, :), allocatable :: qp
  
  n_qp = 2
  !n_qp = ceiling((real(order) + 1.0) / 2.0)

  allocate(qp(n_qp, 2), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of qp array failed."
  
  select case(n_qp)
    case(1)
      qp(:, 1) = (/ 0.0 /)
      qp(:, 2) = (/ 2.0 /)
    case(2)
      qp(:, 1) = (/ -1.0/sqrt(3.0), 1.0/sqrt(3.0) /)
      qp(:, 2) = (/ 1.0, 1.0 /)
    case(3)
      qp(:, 1) = (/ -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0) /)
      qp(:, 2) = (/ 5.0/9.0, 8.0/9.0, 5.0/9.0 /)
    case default
      write(*,*) "Error in selecting quadrature rule."
  end select
end subroutine quadrature


subroutine commandline(length, n_el, order)
  implicit none
  ! define subroutine parameters
  real(rk), intent(out) :: length
  integer, intent(out)  :: n_el
  integer, intent(out)  :: order
 
  ! define local variables 
  integer :: nargs            ! number of command line arguments
  integer :: i                ! looping variable
  character(len = 12) :: args ! command line argument

  nargs = command_argument_count()

  do i = 1, nargs
    call get_command_argument(i, args)
  
    ! use internal reads to convert from character to numeric types
    ! (read from the character buffer into the numeric variable)
    select case (i)
      case(1)
        read(args, *) length  
      case(2)
        read(args, *) n_el
      case(3)
        read(args, *) order
      case default
        write(*,*) "Too many command line parameters."
    end select  
  enddo
end subroutine commandline


subroutine initialize(h, x, n_en, n_el, order, n_nodes)
  implicit none
  real(rk), intent(out) :: h 
  integer, intent(out) :: n_en    
  integer, intent(in)  :: n_el    
  integer, intent(in)  :: order   
  integer, intent(out) :: n_nodes 
  real (rk), dimension(:), allocatable, intent(out) :: x 

  integer :: i ! looping variable

  h = length / real(n_el)
  n_en = order + 1
  n_nodes = (order + 1) * n_el - (n_el - 1)
   
  ! allocate memory for the vector of node coordinates
  allocate(x(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of x array failed."

  do i = 1, size(x)
    x(i) = real((i - 1)) * h / real((n_en - 1))
  end do
end subroutine initialize


END PROGRAM main
