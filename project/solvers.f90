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

end module solvers
