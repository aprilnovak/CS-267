module mesh

implicit none

real(8), save :: h                 ! element length
integer, save :: n_nodes_global    ! number of global nodes
real(8), save, allocatable :: x(:) ! global mesh coordinates

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

end module mesh
