! program to solve the heat equation (Dirichlet boundary conditions only)
! using domain decomposition. This version forms the foundation from which
! I created the actual parallel code - this code is serial, and simulates
! the presence of multiple processes by looping over ranks.

! serial interface problem version

PROGRAM main

implicit none

! define double precision with 8 digits of accuracy
! and an exponent range of +- 30
integer, parameter :: rk = selected_real_kind(8, 30)

! variables for overall execution
integer  :: AllocateStatus     ! variable to hold memory allocation success
integer  :: i, j, q            ! loop iteration variables
real(rk) :: start              ! holds start run time
real(rk) :: finish             ! holds end run time

! variables to define the global problem
integer  :: n_el               ! number of (global) elements
integer  :: n_en               ! number of nodes per element
integer  :: n_nodes            ! number of (global) nodes
integer  :: order              ! polynomial order
integer  :: n_qp               ! number of quadrature points
real(rk) :: length             ! length of the domain (1-D)
real(rk) :: h                  ! length of one element
real(rk) :: k                  ! thermal conductivity
real(rk) :: source             ! uniform heat source
real(rk) :: leftBC             ! left Dirichlet boundary condition value
real(rk) :: rightBC            ! right Dirichlet boundary condition value

! variables to define the CG solver
integer  :: cnt                ! number of CG iterations
real(rk) :: theta              ! CG coefficient
real(rk) :: lambda             ! CG coefficient
real(rk) :: convergence        ! difference between CG iterations
real(rk) :: tol                ! CG convergence tolerance
real(rk) :: startCG            ! start, CG
real(rk) :: endCG              ! end, CG
real(rk) :: m                  ! slope of line

! parallel variables
integer  :: numprocs           ! number of processors
integer  :: maxperproc         ! maximum number of elements per processor
integer  :: rank               ! processor rank
integer  :: face               ! ID number of interface problem
integer  :: iter               ! domain decomposition solution iteration counter
integer  :: sidenum            ! half-thickness of interface layer
integer  :: leftnode           ! node number at left of interface layer
integer  :: rightnode          ! node number at right of interface layer
integer  :: ddcnt              ! domain decomposition counter
real(rk) :: itererror          ! whole-loop iteration error
real(rk) :: ddtol              ! domain decomposition loop tolerance

real(rk), dimension(:),    allocatable :: prev     ! previous interface values
real(rk), dimension(:),    allocatable :: soln     ! global solution vector
real(rk), dimension(:, :), allocatable :: BCvals   ! values of BCs for each domain
real(rk), dimension(:),    allocatable :: xel      ! coordinates in each domain
integer,  dimension(:),    allocatable :: numnodes ! number of nodes in each domain
integer,  dimension(:),    allocatable :: elems    ! n_el in each domain
integer,  dimension(:, :), allocatable :: edges    ! nodes on edge of each domain
integer,  dimension(:, :), allocatable :: LM       ! location matrix
integer,  dimension(2)                 :: BCs      ! boundary condition nodes
real(rk), dimension(:),    allocatable :: qp       ! quadrature points
real(rk), dimension(:),    allocatable :: wt       ! quadrature weights
real(rk), dimension(:),    allocatable :: x        ! coordinates of the nodes
real(rk), dimension(:, :), allocatable :: kel      ! elemental stiffness matrix
real(rk), dimension(:),    allocatable :: rel      ! elemental load vector
real(rk), dimension(:, :), allocatable :: phi      ! shape functions
real(rk), dimension(:, :), allocatable :: dphi     ! shape function derivatives
real(rk), dimension(:),    allocatable :: rglob    ! global load vector
real(rk), dimension(:),    allocatable :: a        ! CG solution iterates
real(rk), dimension(:),    allocatable :: z        ! CG update iterates
real(rk), dimension(:),    allocatable :: res      ! solution residual

k = 1.0        ! thermal conductivity
source = 10.0  ! heat source
tol = 0.0001   ! CG convergence tolerance
ddtol = 0.0005 ! domain decomposition loop tolerance
sidenum = 1    ! elements on each side of the layers

call cpu_time(start)
order = 1 ! only works for linear elements

call commandline(n_el, length, leftBC, rightBC)   ! parse command line args
call initialize(h, x, n_en, n_el, order, n_nodes) ! initialize problem vars
call quadrature(order, n_qp)                      ! initialize quadrature
call phi_val(order, qp)                           ! initialize shape functions
call elementalmatrices()                          ! form elemental matrices and vectors

! allocate the final solution vector before we change the value of n_nodes
allocate(soln(n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of soln array failed."
soln = 0.0

numprocs = 4
call initializedecomp()                   ! initialize domain decomposition

ddcnt = 0
do while (itererror > ddtol)
  ! save the previous values of the interface BCs
  prev = BCvals(1, 2:numprocs)

  do rank = 1, numprocs ! solve over each domain
    n_el = elems(rank)
    n_nodes = numnodes(rank)
    
    allocate(xel(numnodes(rank)), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of xel array failed."
    allocate(LM(n_en, n_el), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of LM array failed."
    allocate(rglob(n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of rglob array failed."
    allocate(a(n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of a array failed."
    allocate(z(n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of z array failed."
    allocate(res(n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of res array failed."
    
    xel = x(edges(1, rank):edges(2, rank))
    
    call locationmatrix()                       ! form the location matrix
    call globalload()                           ! form the global load vector
    
    BCs = (/ 1, n_nodes /)                      ! BCs are always applied on end nodes
    
    rglob(BCs(1)) = BCvals(1, rank)
    rglob(BCs(2)) = BCvals(2, rank)   
   
    ! perform CG solve 
    call conjugategradient(BCvals(2, rank), BCvals(1, rank))
    
    ! save results to the global solution vector
    soln(edges(1, rank):edges(2, rank)) = a
  
    ! deallocate memory before next processor begins its solve
    deallocate(xel, LM, rglob, a, z, res)
  end do ! ends solution of all decomposed domains

  ! solve the numprocs - 1 interface problems
  do face = 1, numprocs - 1
    n_el = sidenum * 2
    n_nodes = n_el * n_en - (n_el - 1)
   
    allocate(xel(n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of xel array failed."
    allocate(LM(n_en, n_el), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of LM array failed."
    allocate(rglob(n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of rglob array failed."
    allocate(a(n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of a array failed."
    allocate(z(n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of z array failed."
    allocate(res(n_nodes), stat = AllocateStatus)
    if (AllocateStatus /= 0) STOP "Allocation of res array failed."
   
    leftnode  = edges(2, face) - sidenum 
    rightnode = edges(2, face) + sidenum
    xel = x(leftnode:rightnode)
    
    call locationmatrix()                       ! form the location matrix
    call globalload()                           ! form the global load vector
    
    BCs = (/1, n_nodes/)                        ! assign BCs on ends of domain 
    rglob(BCs(1)) = soln(leftnode)
    rglob(BCs(2)) = soln(rightnode)
    
    call conjugategradient(rglob(1), rglob(n_nodes))
    
    ! update the BCvals matrix
    BCvals(2, face) = a(n_en * sidenum)
    BCvals(1, face + 1) = a(n_en * sidenum)
  
    deallocate(xel, LM, rglob, a, z, res)
  end do ! ends solution of all interface problems

  ! compute itererror to evaluate whether to continue the iterations
  itererror = 0.0
  do i = 1, numprocs - 1
    itererror = itererror + abs(BCvals(1, i + 1) - prev(i))
  end do
  
  if (ddcnt == 0) then
    ! write to an output file. If this file exists, it will be re-written.
    open(1, file='output.txt', iostat=AllocateStatus, status="replace")
    if (AllocateStatus /= 0) STOP "output.txt file opening failed."
  end if

  write(1, *) soln(:)
 
  open(2, file='timing.txt', status='old', action='write', &
    form='formatted', position='append')
  write(2, *), n_el, finish - start, endCG - startCG, cnt

ddcnt = ddcnt + 1
end do ! ends outermost domain decomposition loop



call cpu_time(finish)
print *, 'runtime: ', finish - start
print *, 'DD iterations: ', ddcnt

! deallocate variables shared by all processors
deallocate(qp, wt, x, kel, rel, phi, dphi)
deallocate(elems, edges, BCvals, prev)





CONTAINS ! define all internal procedures

subroutine initializedecomp()
  implicit none

  allocate(numnodes(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of numnodes array failed."
  allocate(elems(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of elems array failed."
  allocate(edges(2, numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of edges array failed."
  allocate(BCvals(2, numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of BCvals array failed."
  
  ! distribute the elements among the processors 
  maxperproc = (n_el + numprocs - 1) / numprocs
  elems = maxperproc
  j = maxperproc * numprocs - n_el
  
  i = 1
  do while (j > 0)
    elems(i) = elems(i) - 1
    i = i + 1
    j = j - 1
    if (i == numprocs + 1) i = 1
  end do
 
  ! assign the numbers of nodes in each domain 
  do j = 1, numprocs
    numnodes(j) = elems(j) * n_en - (elems(j) - 1)
  end do
  
  ! assign the global node numbers that correspond to the edges of each domain
  edges(:, i) = (/1, elems(i) * n_en - (elems(i) - 1)/)
  do i = 2, numprocs
    edges(:, i) = (/edges(2, i - 1), edges(2, i - 1) + elems(i) * n_en - elems(i) /)
  end do
  
  ! assign an initial guess to all boundaries, and then assign
  ! user-defined BCs to the very edge nodes. This initial guess
  ! is a straight line between the two endpoints
  m = (rightBC - leftBC) / length
  BCvals(1, 1) = leftBC
  BCvals(2, numprocs) = rightBC

  do i = 1, numprocs - 1
    BCvals(2, i) = m * x(edges(2, i)) + leftBC
    BCvals(1, i + 1) = BCvals(2, i)
  end do

  allocate(prev(numprocs - 1), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of prev array failed."
  
  ! assign an initial itererror (dummy value to enter the loop)
  itererror = 1
end subroutine initializedecomp



subroutine globalload()
  implicit none
  rglob = 0.0
  do q = 1, n_el
    do i = 1, n_en
      rglob(LM(i, q)) = rglob(LM(i, q)) + rel(i)
    end do
  end do
end subroutine globalload


subroutine elementalmatrices()
  implicit none

  allocate(kel(n_en, n_en), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of kel array failed."
  allocate(rel(n_en), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rel array failed."

  kel = 0.0
  rel = 0.0
  do q = 1, n_qp
    do i = 1, n_en
      rel(i) = rel(i) + wt(q) * source * phi(i, q) * h *h / 2.0
      do j = 1, n_en
        kel(i, j) = kel(i, j) + wt(q) * dphi(i, q) * k * dphi(j, q) * 2.0
      end do
    end do
  end do
end subroutine elementalmatrices


subroutine conjugategradient(rightBC, leftBC)
  implicit none
  real(rk), intent(in) :: rightBC, leftBC

  ! initial guess is a straight line between the two endpoints
  m           = (rightBC - leftBC) / (xel(n_nodes) - xel(1))
  a           = m * (xel - xel(1)) + leftBC
  res         = rglob - sparse_mult(kel, LM, a)
  z           = res
  lambda      = dotprod(z, res)/dotprod(z, sparse_mult(kel, LM, z))
  a           = a + lambda * z
  res         = rglob - sparse_mult(kel, LM, a)
  convergence = 0.0
 
  ! convergence is assessed by the magnitude of the residual
  do i = 1, n_nodes
    convergence = convergence + abs(res(i))
  end do
  
  cnt = 0
  do while (convergence > tol)
    theta    = sparse_mult_dot(kel, LM, z, res) / sparse_mult_dot(kel, LM, z, z)
    z        = res - theta * z
    lambda   = dotprod(z, res) / sparse_mult_dot(kel, LM, z, z)
    a        = a + lambda * z
    res      = rglob - sparse_mult(kel, LM, a)
    
    convergence = 0.0
    do i = 1, n_nodes
      convergence = convergence + abs(res(i))
    end do
    
  cnt = cnt + 1
  end do
end subroutine conjugategradient


integer function kronecker(i, j)
  implicit none
  integer :: i, j
  kronecker = int((float((i + j) - abs(i - j))) / (float((i + j) + abs(i - j))))
end function kronecker


real function dotprod(vec1, vec2)
  implicit none
  real(rk) :: vec1(:)
  real(rk) :: vec2(:)

  integer  :: i ! looping variable

  dotprod = 0.0  
  do i = 1, size(vec1)
    dotprod = dotprod + vec1(i) * vec2(i)
  end do
end function dotprod

function sparse_mult_dot(matrix, LM, vector, vecdot)
  implicit none
  real(rk) :: matrix(:, :) ! elementary matrix (assumed-shape array)
  real(rk) :: vector(:)    ! full vector (assumed-shape array)
  real(rk) :: vecdot(:)
  integer  :: LM(:, :)     ! location matrix
 
  ! return value of function
  real(rk) :: sparse_mult_dot
  integer  :: i, j, q ! looping variables
 
  sparse_mult_dot = 0.0
  do q = 1, n_el ! loop over the elements
    do i = 1, n_en ! loop over all entries in kel
      if (any(BCs == LM(i, q))) then 
        do j = 1, n_en
          sparse_mult_dot = sparse_mult_dot + vecdot(LM(i, q)) * kronecker(LM(i, q), LM(j, q)) * vector(LM(j, q))
        end do
      else
        ! implicitly assumes that the matrix is symmetric (ok for this application)
        do j = 1, n_en
          sparse_mult_dot = sparse_mult_dot + vecdot(LM(i, q)) * matrix(j, i) * vector(LM(j, q))
        end do
      end if
    end do
  end do
end function sparse_mult_dot


function sparse_mult(matrix, LM, vector)
  implicit none
  real(rk) :: matrix(:, :) ! elementary matrix (assumed-shape array)
  real(rk) :: vector(:)    ! full vector (assumed-shape array)
  integer  :: LM(:, :)     ! location matrix
 
  ! return value of function, as an automatic array
  real(rk) :: sparse_mult(size(vector))
  
  integer :: i, j, q ! looping variables
  sparse_mult = 0.0
   
  do q = 1, n_el ! loop over the elements
    do i = 1, n_en ! loop over all entries in kel
      if (any(BCs == LM(i, q))) then 
        do j = 1, n_en
          ! diagonal terms set to 1.0, off-diagonal set to 0.0
          sparse_mult(LM(i, q)) = sparse_mult(LM(i, q)) + &
                       kronecker(LM(i, q), LM(j, q)) * vector(LM(j, q))
        end do
      else
        do j = 1, n_en
          sparse_mult(LM(i, q)) = sparse_mult(LM(i, q)) + matrix(i, j) * vector(LM(j, q))
        end do
      end if
    end do
  end do
end function sparse_mult


subroutine locationmatrix()
  ! forms the location matrix, which is global in the calling program
  ! fills column-by-column (each column pertains to an element)
  implicit none
  integer :: i, j       ! looping variables
  
  do j = 1, n_el
    do i = 1, n_en
      LM(i, j) = (j - 1) * (n_en - 1) + i
    end do
  end do
end subroutine locationmatrix


subroutine phi_val(order, qp)
! populate phi and dphi, which are global to the calling program
  implicit none
  integer,  intent(in)  :: order
  real(rk), intent(in)  :: qp(:)

  allocate(phi(order + 1, n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of phi array failed."
  allocate(dphi(order + 1, n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of dphi array failed."

  select case(order)
    case(1)
      phi(1, :)  = (1.0 - qp(:)) / 2.0
      phi(2, :)  = (1.0 + qp(:)) / 2.0
      dphi(1, :) = -1.0 / 2.0
      dphi(2, :) =  1.0 / 2.0
    case(2)
      phi(1, :)  = qp(:) * (qp(:) - 1.0) / 2.0
      phi(2, :)  = (1.0 - qp(:)) * (1.0 + qp(:))
      phi(3, :)  = qp(:) * (qp(:) + 1.0) / 2.0
      dphi(1, :) = (2.0 * qp(:) - 1.0) / 2.0
      dphi(2, :) = 1.0 - qp(:) * qp(:)
      dphi(3, :) = (2.0 * qp(:) + 1.0) / 2.0
    case default
      write(*,*) "polynomial order not supported."
  end select
end subroutine phi_val


subroutine quadrature(order, n_qp)
  implicit none

  integer, intent(in)  :: order
  integer, intent(out) :: n_qp
  
  n_qp = 2
  !n_qp = ceiling((real(order) + 1.0) / 2.0)

  allocate(qp(n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of qp array failed."
  allocate(wt(n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of wt array failed."
  
  select case(n_qp)
    case(1)
      qp = (/ 0.0 /)
      wt = (/ 2.0 /)
    case(2)
      qp = (/ -1.0/sqrt(3.0), 1.0/sqrt(3.0) /)
      wt = (/ 1.0, 1.0 /)
    case(3)
      qp = (/ -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0) /)
      wt = (/ 5.0/9.0, 8.0/9.0, 5.0/9.0 /)
    case default
      write(*,*) "Error in selecting quadrature rule."
  end select
end subroutine quadrature


subroutine commandline(n_el, length, leftBC, rightBC)
  implicit none
  integer, intent(out)  :: n_el
  real(rk), intent(out) :: length
  real(rk), intent(out) :: leftBC
  real(rk), intent(out) :: rightBC
 
  integer :: nargs            ! number of command line arguments
  integer :: i                ! looping variable
  character(len = 12) :: args ! command line argument

  nargs = command_argument_count()

  do i = 1, nargs
    call get_command_argument(i, args)
  
    ! use internal reads to convert from character to numeric types
    ! (read from the character buffer into the numeric variable)
    select case (i)
      case(1)
        read(args, *), length
      case(2)
        read(args, *) n_el
      case(3)
        read(args, *) leftBC
      case(4)
        read(args, *) rightBC
      case default
        write(*,*) "Too many command line parameters."
    end select  
  enddo
end subroutine commandline


subroutine initialize(h, x, n_en, n_el, order, n_nodes)
  implicit none
  real(rk), intent(out) :: h 
  integer, intent(out) :: n_en    
  integer, intent(in)  :: n_el    
  integer, intent(in)  :: order   
  integer, intent(out) :: n_nodes 
  real (rk), dimension(:), allocatable, intent(out) :: x 

  integer :: i ! looping variable

  h = length / real(n_el)
  n_en = order + 1
  n_nodes = (order + 1) * n_el - (n_el - 1)

  ! allocate memory for the vector of node coordinates
  allocate(x(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of x array failed."

  do i = 1, size(x)
    x(i) = real(i - 1) * h / real(n_en - 1)
  end do
end subroutine initialize


END PROGRAM main
