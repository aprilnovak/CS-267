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

function form_csr(LM, LMcount, n_nodes) result(rows)
  integer, intent(in)      :: LM(:, :)
  integer, intent(in)      :: LMcount(:)
  integer, intent(in)      :: n_nodes
 
  type(row) :: rows(n_nodes)
  integer   :: i, j, q, n_el
 
  n_el    = size(LM(1, :))
 
  ! allocate space for the elements of the rows data structure
  ! The formula used to determine how many contributions are made in a row
  ! would need to be redetermined for higher-dimensional meshes.
  do i = 1, n_nodes
    allocate(rows(i)%values(2 * LMcount(i) - 2), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of values array failed."
    allocate(rows(i)%columns(2 * LMcount(i) - 2), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of columns array failed."
  end do
 
  ! populate vectors of values and the columns they belong in
  do q = 1, n_el
    do i = 1, 2
      do j = 1, 2
        rows(LM(i, q))%columns(rows(LM(i, q))%entri) = LM(j, q)
        rows(LM(i, q))%values(rows(LM(i, q))%entri)  = kel(i, j)
        rows(LM(i, q))%entri = rows(LM(i, q))%entri + 1
      end do
    end do
  end do
end function form_csr

end module csr_storage
