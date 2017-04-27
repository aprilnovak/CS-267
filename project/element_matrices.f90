module element_matrices

implicit none

! shape function type
type shapefunction
  real(8), allocatable :: phi(:, :)
  real(8), allocatable :: dphi(:, :)
end type shapefunction

! elemental load vector, global variable
real(8), save :: rel(2)

! elemental stiffness matrix, global variable
real(8), save :: kel(2, 2)

! shape functions, private variable (only needed within this module).
! You can declare a variable as private, but still use the module
! procedures that directly reference that variable, you just can't 
! _directly_ modify that private variable
type(shapefunction), private :: func
integer, private             :: AllocateStatus

contains

subroutine define_shapefunctions()
! assign the values for the shape functions and their derivatives
  use quadrature, only: set ! quadrature set

  allocate(func%phi(2, set%n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of phi array failed."
  allocate(func%dphi(2, set%n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of dphi array failed."

  func%phi(1, :)  = (1.0 - set%qp(:)) / 2.0
  func%phi(2, :)  = (1.0 + set%qp(:)) / 2.0
  func%dphi(1, :) = -1.0 / 2.0
  func%dphi(2, :) =  1.0 / 2.0
end subroutine define_shapefunctions


subroutine elementalload()
! assemble the elemental load vector
  use quadrature, only: set   ! quadrature set
  use read_data, only: source ! heat source
  use mesh, only: global      ! contains global element width

  integer             :: q, i
  
  rel = 0.0
  do q = 1, set%n_qp
    do i = 1, 2
      rel(i) = rel(i) + set%wt(q) * source * func%phi(i, q) * global%h * global%h / 2.0
    end do
  end do
end subroutine elementalload


subroutine elementalstiffness()
! assemble the elemental stiffness matrix  
  use quadrature, only: set ! quadrature set
  use read_data, only: k
  
  integer             :: q, i

  kel = 0.0
  do q = 1, set%n_qp
    do i = 1, 2
      kel(i, 1) = kel(i, 1) + set%wt(q) * func%dphi(i, q) * k * func%dphi(1, q) * 2.0
      kel(i, 2) = kel(i, 2) + set%wt(q) * func%dphi(i, q) * k * func%dphi(2, q) * 2.0
    end do
  end do
end subroutine elementalstiffness


function globalvector(lmatrix, elemvec, n_nodes, factor) result(vector)
! assemble global vector from location matrix and elemental vector
  use mesh, only: LM, global

  type(LM), intent(in) :: lmatrix
  real(8), intent(in)  :: elemvec(:)
  real(8), intent(in), optional  :: factor(:) ! multiplicative factor (lengths)
  integer, intent(in)  :: n_nodes

  real(8)                :: vector(n_nodes)
  integer                :: q, n

  n = size(lmatrix%matrix(1, :))

  vector = 0.0

  if (present(factor)) then
    do q = 1, n
      vector(lmatrix%matrix(:, q)) = vector(lmatrix%matrix(:, q)) + factor(q) * factor(q) * elemvec(:) / (global%h * global%h)
    end do
  else
    do q = 1, n
       vector(lmatrix%matrix(:, q)) = vector(lmatrix%matrix(:, q)) + elemvec(:)
    end do
  end if
end function globalvector


subroutine dealloc_shapefunctions()
! deallocate the shape functions
  deallocate(func%phi, func%dphi)
end subroutine dealloc_shapefunctions

end module element_matrices
