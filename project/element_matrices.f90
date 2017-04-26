module element_matrices

use quadrature, only: quadset ! included for type definition
implicit none

integer, private     :: AllocateStatus

type shapefunction
  real(8), allocatable :: phi(:, :)
  real(8), allocatable :: dphi(:, :)
end type shapefunction

contains

function define_shapefunctions(set) result(func)
  type(quadset), intent(in)  :: set
  type(shapefunction)        :: func

  allocate(func%phi(2, set%n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of phi array failed."
  allocate(func%dphi(2, set%n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of dphi array failed."

  func%phi(1, :)  = (1.0 - set%qp(:)) / 2.0
  func%phi(2, :)  = (1.0 + set%qp(:)) / 2.0
  func%dphi(1, :) = -1.0 / 2.0
  func%dphi(2, :) =  1.0 / 2.0
end function define_shapefunctions

end module element_matrices
