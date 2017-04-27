module mesh

implicit none

type LM ! location matrix
  integer, allocatable :: matrix(:, :) ! connectivity matrix
  integer, allocatable :: cnt(:)       ! node counts
end type LM

type decomp ! decomposed domain
  integer, allocatable :: elems(:)          ! elements per domain
  integer, allocatable :: edges(:, :)       ! nodes on edges of each domain
  integer, allocatable :: recv_displs(:)    ! displacements from start
end type decomp

type geom
  integer :: n_nodes
  integer :: n_el
  real(8) :: length
  real(8) :: h
  real, allocatable :: x(:) 
end type geom

type(decomp), save :: domains   ! holds domain decomposition information
type(geom),   save :: global    ! holds global mesh information
type(geom),   save :: coarse    ! holds coarse mesh information
type(geom), allocatable,  save :: dd(:)        ! holds domain decomposition information

integer, private :: AllocateStatus

contains

subroutine initialize_global_mesh()
  integer :: i
  global%h = global%length / real(global%n_el)
  global%n_nodes = global%n_el + 1
 
  allocate(global%x(global%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of global%x array failed."
 
  do i = 1, global%n_nodes
    global%x(i) = real(i - 1) * global%h
  end do
end subroutine initialize_global_mesh


subroutine initialize_coarse_mesh()
  coarse%n_el = size(dd)
  coarse%n_nodes = coarse%n_el + 1
end subroutine initialize_coarse_mesh


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

  integer, intent(in) :: numprocs
  integer :: mx, i, j
 
  allocate(domains%elems(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocate of elems array failed."
  allocate(domains%edges(2, numprocs + 1), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocate of edges array failed."
  allocate(domains%recv_displs(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocate of recv_displs array failed."

  allocate(dd(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocate of dd array failed."
 

  mx = (global%n_el + numprocs - 1) / numprocs
  
  do i = 1, numprocs
    dd(i)%n_el = mx
  end do

  domains%elems = mx
 
  i = 1
  do j = mx * numprocs - global%n_el, 1, -1
    domains%elems(i) = domains%elems(i) - 1
    dd(i)%n_el = dd(i)%n_el - 1
    i = i + 1
    if (i == numprocs + 1) i = 1
  end do
 
  do i = 1, numprocs
    dd(i)%n_nodes = dd(i)%n_el + 1
  end do 

  domains%edges(:, 1) = (/1, domains%elems(1) + 1/)
  do i = 2, numprocs
    domains%edges(:, i) = (/domains%edges(2, i - 1), &
                            domains%edges(2, i - 1) + domains%elems(i)/)
  end do
  
  do i = 1, numprocs
    allocate(dd(i)%x(dd(i)%n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocate of dd(i)%x array failed."
    dd(i)%x = global%x(domains%edges(1, i):domains%edges(2, i))
  end do 
 
  domains%recv_displs = 0
  do i = 2, numprocs
    domains%recv_displs(i) = domains%recv_displs(i - 1) + domains%elems(i - 1)
  end do
end subroutine initialize_domain_decomposition


subroutine dealloc_global_mesh()
  deallocate(global%x)
end subroutine dealloc_global_mesh

subroutine dealloc_domains()
  deallocate(domains%elems, domains%edges, domains%recv_displs)
  deallocate(dd)
end subroutine dealloc_domains

end module mesh
