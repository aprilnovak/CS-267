module quadrature

implicit none

type quadset
  real(8), allocatable :: wt(:)
  real(8), allocatable :: qp(:)
  integer              :: n_qp
end type quadset

type(quadset), save :: set         ! quadrature set
integer, private :: AllocateStatus ! allocation check

contains

subroutine define_quadset()
  use read_data, only: n=>n_qp     ! number of quadrature points

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
end subroutine define_quadset


subroutine dealloc_quadset()
  deallocate(set%wt, set%qp)
end subroutine dealloc_quadset

end module quadrature
