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


function straightline(geometry)
  use mesh, only: geom
  type(geom), intent(in) :: geometry

  real(8) :: m
  real(8) :: straightline(geometry%n_nodes)
  m = (geometry%BCvals(2) - geometry%BCvals(1)) / &
      (geometry%x(geometry%n_nodes) - geometry%x(1))
  straightline = m * (geometry%x - geometry%x(1)) + geometry%BCvals(1)
end function straightline


function conjugategradient(rows, guess, rglob, BCs, reltol) result(a)
  use csr_storage, only: row
  real(8), intent(in)    :: guess(:)      ! starting guess
  real(8), intent(in)    :: rglob(:)      ! rhs vector
  integer, intent(in)    :: BCs(:)        ! BC nodes
  type(row), intent(in)  :: rows(:)       ! kglob in CSR form
  real(8), intent(in), optional :: reltol ! relative CG tolerance

  real(8) :: a(size(guess))

  ! local variables
  real(8) :: z(size(a)), res(size(a))
  real(8) :: lambda, theta, internal, tol, conv
  integer :: cnt, n

  n = size(a)

  res         = rglob - csr_mult(rows, a, BCs)
  internal    = sum(abs(res))

  ! set relative tolerance for convergence, using 0.001 as default
  if (present(reltol)) then
    tol = reltol * internal
  else
    tol = 0.001 * internal
  end if

  z      = res
  lambda = dotprod(z, res) / dotprod(csr_mult(rows, z, BCs), z)

  a      = guess + lambda * z
  res    = rglob - csr_mult(rows, a, BCs)
  conv   = sum(abs(res))

  cnt = 0
  do while (conv > tol)
    theta  = dotprod(csr_mult(rows, z, BCs), res) / dotprod(csr_mult(rows, z, BCs), z)
    z      = res - theta * z
    lambda = dotprod(z, res) / dotprod(csr_mult(rows, z, BCs), z)
    a      = a + lambda * z
    res    = rglob - csr_mult(rows, a, BCs)
    conv   = sum(abs(res))

  cnt = cnt + 1
  end do
end function conjugategradient


end module solvers
