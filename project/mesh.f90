module mesh

implicit none

real(8), save :: h                 ! element length
integer, save :: n_nodes_global    ! number of global nodes
real(8), save, allocatable :: x(:) ! global mesh coordinates

type LM
  integer, allocatable :: matrix(:, :) ! connectivity matrix
  integer, allocatable :: cnt(:)       ! node counts
end type LM

type decomp
  integer, allocatable :: elems(:)          ! elements per domain
  integer, allocatable :: numnodes(:)       ! nodes per domain
  integer, allocatable :: edges(:, :)       ! nodes on edges of each domain
  integer, allocatable :: recv_displs(:)    ! displacements from start
end type decomp

type(decomp), save :: domains ! holds domain decomposition information

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


subroutine initialize_domain_decomposition(numprocs)
! elems: elements per parallel MPI process
! numnodes: nodes per parallel MPI process
! edges: nodes on edge of each parallel MPI process 
  use read_data, only: n_el=>n_el_global

  integer, intent(in) :: numprocs
  integer :: mx, i, j
 
  allocate(domains%elems(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocate of elems array failed."
  allocate(domains%numnodes(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocate of numnodes array failed."
  allocate(domains%edges(2, numprocs + 1), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocate of edges array failed."
  allocate(domains%recv_displs(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocate of recv_displs array failed."

  mx = (n_el + numprocs - 1) / numprocs
  domains%elems = mx
 
  i = 1
  do j = mx * numprocs - n_el, 1, -1
    domains%elems(i) = domains%elems(i) - 1
    i = i + 1
    if (i == numprocs + 1) i = 1
  end do
 
  do j = 1, numprocs
    domains%numnodes(j) = domains%elems(j) * 2 - (domains%elems(j) - 1)
  end do
 
  domains%edges(:, 1) = (/1, domains%elems(1) * 2 - (domains%elems(1) - 1)/)
  do i = 2, numprocs
    domains%edges(:, i) = (/domains%edges(2, i - 1), domains%edges(2, i - 1) + domains%elems(i) * 2 - domains%elems(i) /)
  end do
 
  domains%recv_displs = 0
  do i = 2, numprocs
    domains%recv_displs(i) = domains%recv_displs(i - 1) + domains%elems(i - 1)
  end do
end subroutine initialize_domain_decomposition


subroutine dealloc_x()
  deallocate(x)
end subroutine dealloc_x

subroutine dealloc_domains()
  deallocate(domains%elems, domains%numnodes, domains%edges, domains%recv_displs)
end subroutine dealloc_domains

end module mesh
