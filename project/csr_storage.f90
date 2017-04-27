module csr_storage

use element_matrices, only: kel
implicit none

! data type that holds the values and columns for each value for a row
! in a matrix
type row
  real(8), allocatable, dimension(:) :: values     ! values in each row
  integer, allocatable, dimension(:) :: columns    ! nonzero col nums per row
  integer                            :: entri = 1  ! entry next to be filled
end type row

integer, private :: AllocateStatus

contains

function form_csr(l, n_nodes) result(rows)
  use mesh, only: LM
  type(LM), intent(in)     :: l
  integer, intent(in)      :: n_nodes
 
  type(row), allocatable :: rows(:)
  integer   :: i, j, q, n_el
 
  allocate(rows(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rows array failed."

  n_el    = size(l%matrix(1, :))
 
  ! allocate space for the elements of the rows data structure
  ! The formula used to determine how many contributions are made in a row
  ! would need to be redetermined for higher-dimensional meshes.
  do i = 1, n_nodes
    allocate(rows(i)%values(2 * l%cnt(i) - 2), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of values array failed."
    allocate(rows(i)%columns(2 * l%cnt(i) - 2), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of columns array failed."
  end do
 
  ! populate vectors of values and the columns they belong in
  do q = 1, n_el
    do i = 1, 2
      do j = 1, 2
        rows(l%matrix(i, q))%columns(rows(l%matrix(i, q))%entri) = l%matrix(j, q)
        rows(l%matrix(i, q))%values(rows(l%matrix(i, q))%entri)  = kel(i, j)
        rows(l%matrix(i, q))%entri = rows(l%matrix(i, q))%entri + 1
      end do
    end do
  end do
end function form_csr


end module csr_storage
