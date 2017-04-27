module mesh

implicit none

real(8), save :: h                 ! element length
integer, save :: n_nodes_global    ! number of global nodes
real(8), save, allocatable :: x(:) ! global mesh coordinates

type LM
  integer, allocatable :: matrix(:, :) ! connectivity matrix
  integer, allocatable :: cnt(:)       ! node counts
end type LM

integer, private :: AllocateStatus

contains

subroutine initialize_global_mesh()
  use read_data, only: n_el=>n_el_global, length

  integer              :: i
 
  h = length / real(n_el)
  n_nodes_global = n_el + 1
 
  allocate(x(n_nodes_global), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of x array failed."
 
  do i = 1, n_nodes_global
    x(i) = real(i - 1) * h
  end do
end subroutine initialize_global_mesh


function locationmatrix(n_el, n_nodes) result(lmatrix)
  type(LM) :: lmatrix
  integer, intent(in) :: n_el, n_nodes
  integer :: j
 
  allocate(lmatrix%matrix(2, n_el), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of lmatrix%matrix array failed."
  allocate(lmatrix%cnt(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of lmatrix%cnt array failed."

  lmatrix%cnt = 1

  do j = 1, n_el
    lmatrix%matrix(:, j) = (/ j, j + 1 /)
    lmatrix%cnt(lmatrix%matrix(1, j)) = lmatrix%cnt(lmatrix%matrix(1, j)) + 1
    lmatrix%cnt(lmatrix%matrix(2, j)) = lmatrix%cnt(lmatrix%matrix(2, j)) + 1
  end do  
end function locationmatrix


subroutine dealloc_x()
  deallocate(x)
end subroutine dealloc_x

end module mesh
