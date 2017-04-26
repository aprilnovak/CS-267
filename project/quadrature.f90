module quadrature

! specification statements
! a module may use another module, but not itself
implicit none

integer, private :: AllocateStatus

! save statement means that these are visible to all program units that
! use this module
type quadset
  real(8), allocatable :: wt(:)
  real(8), allocatable :: qp(:)
  integer              :: n_qp
end type quadset

contains ! define module procedures (methods that operate on the quadset class)


function definequadset(n) result(set)
  type(quadset)       :: set
  integer, intent(in) :: n

  allocate(set%wt(n), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of set%wt array failed."
  allocate(set%qp(n), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of set%qp array failed."

  set%n_qp = n
  
  if (n == 2) then
    set%qp   = (/ -1.0/sqrt(3.0), 1.0/sqrt(3.0) /)
    set%wt   = (/ 1.0, 1.0 /)
  else
    print *, 'error in quadrature rule selection.'
  end if
end function definequadset


subroutine dealloc_quad(set)
  type(quadset) :: set
  deallocate(set%wt, set%qp)
end subroutine dealloc_quad

end module quadrature
