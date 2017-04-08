! program to solve the heat equation (Dirichlet boundary conditions only)

PROGRAM main

implicit none

! define double precision with 15 digits of accuracy
! and an exponent range of +- 307
integer, parameter :: rk = selected_real_kind(15, 307)

integer  :: AllocateStatus ! variable to hold memory allocation success
integer  :: i, j, q        ! loop iteration variables

integer  :: n_el               ! number of elements
integer  :: n_en               ! number of nodes per element
integer  :: n_nodes            ! total number of nodes
integer  :: order              ! polynomial order
integer  :: n_qp               ! number of quadrature points
real(rk) :: length             ! length of the domain (1-D)
real(rk) :: h                  ! length of one element
real(rk) :: k                  ! thermal conductivity
real(rk) :: source             ! uniform heat source
real(rk) :: leftBC             ! left Dirichlet boundary condition value
real(rk) :: rightBC            ! right Dirichlet boundary condition value
real(rk) :: theta              ! conjugate gradient coefficient
real(rk) :: lambda             ! conjugate gradient coefficient
real(rk) :: convergence        ! difference between successive iterations
real(rk) :: tol                ! conjugate gradient convergence tolerance

integer,  dimension(:, :), allocatable :: LM     ! location matrix
real(rk), dimension(:),    allocatable :: qp     ! quadrature points
real(rk), dimension(:),    allocatable :: wt     ! quadrature weights
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
real(rk) :: val
!integer, dimension(:), allocatable :: LMcount ! number of times each node number appears
!integer, dimension(:), allocatable  :: sparse_mult_result
integer, dimension(5) :: p
integer, dimension(5, 5) :: mat
! initialize the thermal conductivity and heat source
k = 1.0
source = 1.0
tol = 0.001

! test out array modifications manually

p = 0.0
mat(1, :) = (/ 1, 2, 3, 4, 5/)
mat(2, :) = (/6, 7, 8, 9, 1/)
mat(3, :) = (/2, 3, 4, 5, 6/)
mat(4, :) = (/ 7, 8, 9, 1, 2/)
mat(5, :) = (/3, 4, 5, 6, 7/)

do i = 1, 5 ! why won't this add the value of the matrix to the vector?
  p(i) = p(i) + mat(2, i)
  print *, 'p: ', p
end do

print *, 'final p: ', p










call commandline(length, n_el, order, leftBC, rightBC) ! parse command line arguments
call initialize(h, x, n_en, n_el, order, n_nodes)      ! initialize problem variables
call quadrature(order, n_qp, qp, wt)                   ! initialize quadrature rule

! allocate memory for the shape functions - quadrature points
! do not change throughout the simulation
allocate(phi(order + 1, n_qp), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of phi array failed."

! allocate memory for the shape functions - quadrature points
! do not change throughout the simulation
allocate(dphi(order + 1, n_qp), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of dphi array failed."





call phi_val(order, qp)                     ! initialize shape functions

! form the elemental stiffness matrix and load vector
allocate(kel(n_en, n_en), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of kel array failed."
allocate(rel(n_en), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of rel array failed."
!allocate(sparse_mult_result(n_nodes), stat = AllocateStatus)
!if (AllocateStatus /= 0) STOP "Allocation of sparse_mult_result array failed."

kel = 0.0
rel = 0.0
do q = 1, n_qp ! quadrature point is slowest varying
  do i = 1, n_en
    rel(i) = rel(i) + wt(q) * source * phi(i, q) * h *h / 2.0
    do j = 1, n_en
      kel(i, j) = kel(i, j) + wt(q) * dphi(i, q) * k * dphi(j, q) * 2.0
    end do
  end do
end do

! form the location matrix
allocate(LM(n_en, n_el), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of LM array failed."
!allocate(LMcount(n_nodes), stat = AllocateStatus)
!if (AllocateStatus /= 0) STOP "Allocation of LMcount array failed."
call locationmatrix()

! form the global stiffness matrix and global load vector
allocate(kglob(n_nodes, n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of kglob array failed."
allocate(rglob(n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of rglob array failed."


! populate global matrices element by element so that
! elements in LM are accessed by column (fastest)
kglob = 0.0
rglob = 0.0
do q = 1, n_el
  do i = 1, n_en
    rglob(LM(i, q)) = rglob(LM(i, q)) + rel(i)
    do j = 1, n_en
      kglob(LM(i, q), LM(j, q)) = kglob(LM(i, q), LM(j, q)) + kel(i, j)
    end do
  end do
end do


do i = 1, n_nodes
  print *, kglob(i,:)
end do

! for now, this only works for 1-D FEM problems where we know the 
! global array structure

! loop over the LM columns and count number of times each node appears
! (for each element)
!LMcount = 0
!do i = 1, n_en
!  do j = 1, n_el
!    LMcount(LM(i, j)) = LMcount(LM(i, j)) + 1
!  end do
!end do

!j = 1
!do i = 1, n_nodes + 1
!  pt(i) = j
!  j = j + n_en + LMcount(i) - 1
!end do

!print *, 'starting numbers of each row: ', pt


! compute by looping over the elements
!re = 0.0

do q = 1, n_el ! loop over the elements
  do i = 1, n_en ! loop over all entries in kel
    do j = 1, n_en
      !re(LM(i, q)) = re(LM(i, q)) + kel(i, j) * D(LM(j, q))
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

! apply boundary conditions
kglob(1, :) = 0.0
kglob(n_nodes, :) = 0.0
kglob(1, 1) = 1.0                ! set up left BC              
kglob(n_nodes, n_nodes) = 1.0    ! set up right BC
rglob(1) = leftBC                ! left BC value
rglob(n_nodes) = rightBC         ! right BC value     




! conjugate gradient solver for solving kglob * a = rglob
aprev       = 1.0




print *, 'using Fortran intrinsic: ', matmul(kglob, aprev)
print *, 'my function: ', sparse_mult(kel, LM, aprev)


res         = rglob - matmul(kglob, aprev)
zprev       = res
lambda      = dot_product(zprev, res) / dot_product(zprev, matmul(kglob, zprev))
a           = aprev + lambda * zprev
convergence = 0.0

do i = 1, n_nodes
  convergence = convergence + abs(a(i) - aprev(i))
end do

do while (convergence > tol)
  aprev  = a
  res    = rglob - matmul(kglob, aprev)
  theta  = - dot_product(res, matmul(kglob, zprev)) / dot_product(zprev, matmul(kglob, zprev))
  z      = res + theta * zprev
  lambda = dot_product(z, res) / dot_product(z, matmul(kglob, z))
  a      = aprev + lambda * z
  zprev  = z

  convergence = 0.0
  do i = 1, n_nodes
    convergence = convergence + abs(a(i) - aprev(i))
  end do
end do

! write to an output file for plotting. If this file exists, it will be re-written.
open(1, file='output.txt', iostat=AllocateStatus, status="replace")
if (AllocateStatus /= 0) STOP "output.txt file opening failed."

write(1, *) a(:)
! deallocate memory 
deallocate(qp, wt, x, kel, rel, phi, dphi, kglob, rglob, a, aprev, z, zprev, res, LM)



CONTAINS ! define all internal procedures

function sparse_mult(matrix, LM, vector)
 ! pass in an elementary matrix and the full vector
 real(rk) :: matrix(:, :) ! define the function parameters as assumed-shape arrays  
 real(rk) :: vector(:)
 integer  :: LM(:, :)
 !integer  :: n_el, n_en

 ! return value of function, as an automatic array
 real(rk) :: sparse_mult(size(vector))
 
 integer :: i, j, q ! looping variables
 sparse_mult = 0.0
  
  do q = 1, n_el ! loop over the elements
    do i = 1, n_en ! loop over all entries in kel
      do j = 1, n_en
        if ((LM(i, q) == 1) .and. (LM(j, q) /= 1)) then ! left boundary condition
        else if ((LM(i, q) == n_nodes) .and. (LM(j, q) /= n_nodes)) then ! right boundary condition
        else
          sparse_mult(LM(i, q)) = sparse_mult(LM(i, q)) + matrix(i, j) * vector(LM(j, q))
        end if
      end do
    end do
  end do
end function sparse_mult


subroutine locationmatrix()
  ! forms the location matrix, which is global in the calling program
  ! fills column-by-column (each column pertains to an element)
  implicit none
  integer :: i, j       ! looping variables
  
  do j = 1, n_el
    do i = 1, n_en
      LM(i, j) = (j - 1) * (n_en - 1) + i
    end do
  end do
end subroutine locationmatrix


subroutine phi_val(order, qp)
! populate phi and dphi, which are global to the calling program
  implicit none
  integer,  intent(in)  :: order
  real(rk), intent(in)  :: qp(:)

  select case(order)
    case(1)
      phi(1, :)  = (1.0 - qp(:)) / 2.0
      phi(2, :)  = (1.0 + qp(:)) / 2.0
      dphi(1, :) = -1.0 / 2.0
      dphi(2, :) =  1.0 / 2.0
    case(2)
      phi(1, :)  = qp(:) * (qp(:) - 1.0) / 2.0
      phi(2, :)  = (1.0 - qp(:)) * (1.0 + qp(:))
      phi(3, :)  = qp(:) * (qp(:) + 1.0) / 2.0
      dphi(1, :) = (2.0 * qp(:) - 1.0) / 2.0
      dphi(2, :) = 1.0 - qp(:) * qp(:)
      dphi(3, :) = (2.0 * qp(:) + 1.0) / 2.0
    case default
      write(*,*) "phi and dphi not initialized due to polynomial order not being supported."
  end select
end subroutine phi_val


subroutine quadrature(order, n_qp, qp, wt)
  implicit none

  integer, intent(in)  :: order
  integer, intent(out) :: n_qp
  real (rk), dimension(:), allocatable :: qp
  real (rk), dimension(:), allocatable :: wt
  
  n_qp = 2
  !n_qp = ceiling((real(order) + 1.0) / 2.0)

  allocate(qp(n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of qp array failed."
  allocate(wt(n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of wt array failed."
  
  select case(n_qp)
    case(1)
      qp = (/ 0.0 /)
      wt = (/ 2.0 /)
    case(2)
      qp = (/ -1.0/sqrt(3.0), 1.0/sqrt(3.0) /)
      wt = (/ 1.0, 1.0 /)
    case(3)
      qp = (/ -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0) /)
      wt = (/ 5.0/9.0, 8.0/9.0, 5.0/9.0 /)
    case default
      write(*,*) "Error in selecting quadrature rule."
  end select
end subroutine quadrature


subroutine commandline(length, n_el, order, leftBC, rightBC)
  implicit none
  real(rk), intent(out) :: length
  integer, intent(out)  :: n_el
  integer, intent(out)  :: order
  real(rk), intent(out) :: leftBC
  real(rk), intent(out) :: rightBC
 
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
      case(4)
        read(args, *) leftBC
      case(5)
        read(args, *) rightBC
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
