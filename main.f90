! program to solve the heat equation (Dirichlet boundary conditions only)


PROGRAM main

! turn off implicit typing
implicit none

! define double precision with 15 digits of accuracy
! and an exponent range of +- 307
integer, parameter :: rk = selected_real_kind(15, 307)

! program variables related to execution
integer             :: AllocateStatus ! variable to hold memory allocation success
integer             :: i, j, q        ! loop iteration variables

! program variables related to FEM solve
integer  :: n_el               ! number of elements
integer  :: n_en               ! number of nodes per element
integer  :: n_nodes            ! total number of nodes
integer  :: order              ! polynomial order
integer  :: n_qp               ! number of quadrature points
real(rk) :: length             ! length of the domain (1-D)
real(rk) :: h                  ! length of one element
real(rk) :: k                  ! thermal conductivity

integer, dimension(:, :), allocatable  :: LM     ! location matrix

real(rk), dimension(:, :), allocatable :: qp     ! quadrature points and weights
real(rk), dimension(:), allocatable :: x         ! coordinates of the nodes
real(rk), dimension(:, :), allocatable :: kel    ! elemental stiffness matrix
real(rk), dimension(:, :), allocatable :: phi    ! shape functions
real(rk), dimension(:, :), allocatable :: dphi   ! shape function derivatives


! initialize the thermal conductivity
k = 1.0

! parse command line arguments
call commandline(length, n_el, order)

! initialize problem variables
call initialize(h, x, n_en, n_el, order, n_nodes)

! initialize the quadrature rule
call quadrature(order, n_qp, qp)

! initialize shape function values
call phi_val(order, qp, phi, dphi)

! form the elemental stiffness matrix
allocate(kel(n_en, n_en), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of kel array failed."

kel = 0.0

do q = 1, n_qp
  do i = 1, n_en
    do j = 1, n_en
      ! only works for 1-D elements
      kel(i, j) = kel(i, j) + qp(q, 2) * dphi(q, i) * k * dphi(q, j) * h / 2.0
    end do
  end do
end do

! form the location matrix (each column belongs to an element)
allocate(LM(n_el, n_en), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of LM array failed."

do i = 1, n_el
  do j = 1, n_en
    LM(i, j) = 1 + (i - 1) * (n_en - 1) + (j - 1)
  end do
end do








! write program variables to user
print *, 'Beginning FEM solve of the heat equation with:'
print *, 'length: ', length
print *, 'number of elements: ', n_el
print *, 'polynomial order: ', order
print *, 'Total nodes: ', n_nodes
print *, 'Number of nodes per element: ', n_en
!print *, 'Domain node locations: ', x
print *, 'Quadrature points:'
print *, qp
print *, 'Phi values:'
print *, phi(:, 1)
print *, phi(:, 2)
print *, 'Derivative values:'
print *, dphi(:, 1)
print *, dphi(:, 2)
print *, 'Elemental stiffness matrix:'
print *, kel
print *, 'Location matrix:'
print *, LM(1, :)
print *, LM(2, :)



! deallocate memory 
deallocate(x)
deallocate(qp)
deallocate(phi)
deallocate(dphi)
deallocate(LM)
CONTAINS

subroutine phi_val(order, qp, phi, dphi)
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
  
  n_qp = ceiling((real(order) + 1.0) / 2.0)

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
