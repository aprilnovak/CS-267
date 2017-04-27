module solvers

implicit none

contains

real function dotprod(vec1, vec2)
  real(8)  :: vec1(:), vec2(:)
  integer  :: i, n

  n = size(vec1)
  dotprod = 0.0

  do i = 1, n
    dotprod = dotprod + vec1(i) * vec2(i)
  end do
end function dotprod


function csr_mult(rows, a, BCs)
  use csr_storage, only: row ! need to know data type row

  type(row), intent(in) :: rows(:)
  real(8), intent(in)   :: a(:)
  integer, intent(in)   :: BCs(:)

  real(8)   :: csr_mult(size(a))
  integer   :: n, i, j
  real(8)   :: accum

  n = size(a)

  do i = 1, n
    accum = 0.0
    do j = 1, size(rows(i)%columns(:))
      accum = accum + rows(i)%values(j) * a(rows(i)%columns(j))
    end do
    csr_mult(i) = accum
  end do

  ! apply boundary conditions
  do i = 1, size(BCs)
    csr_mult(BCs(i)) = a(BCs(i))
  end do
end function csr_mult


end module solvers
