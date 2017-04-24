! non-overlapping interface
! each processor computes one interface problem
! builds from bmain by adding openMP parallelism

! program to solve the heat equation (Dirichlet boundary conditions only)
! using domain decomposition. This version has the domain to the left of
! each interface send a single interface value to the processor to the
! right of it. This process then computes the interface problem, and 
! sends a single udpated value to the processor to the left.

PROGRAM main

implicit none
include 'mpif.h'

! variables for overall execution
integer  :: AllocateStatus     ! variable to hold memory allocation success
integer  :: i, j, q            ! loop iteration variables
real(8)  :: start              ! holds start run time
real(8)  :: finish             ! holds end run time

! variables to define the global problem
integer                               :: n_el_global    ! global elements
integer                               :: n_nodes_global ! global nodes
integer                               :: n_qp           ! number of quadrature points
real(8)                               :: length         ! length of the domain (1-D)
real(8)                               :: h              ! length of one element
real(8)                               :: k              ! thermal conductivity
real(8)                               :: source         ! uniform heat source
real(8)                               :: leftBC         ! left Dirichlet BC value
real(8)                               :: rightBC        ! right Dirichlet BC value
integer, dimension(2)                 :: BCs            ! boundary condition nodes
real(8), dimension(2, 2)              :: kel            ! elemental stiffness matrix
real(8), dimension(2)                 :: rel            ! elemental load vector
real(8), dimension(:),    allocatable :: soln           ! global solution vector
real(8), dimension(:),    allocatable :: qp             ! quadrature points
real(8), dimension(:),    allocatable :: wt             ! quadrature weights
real(8), dimension(:),    allocatable :: x              ! coordinates of the nodes
real(8), dimension(:, :), allocatable :: phi            ! shape functions
real(8), dimension(:, :), allocatable :: dphi           ! shape function derivatives
integer, dimension(:, :), allocatable :: LM             ! location matrix

! variables to define the CG solver
integer                            :: cnt         ! number of CG iterations
real(8)                            :: convergence ! difference between CG iterations
real(8)                            :: tol         ! CG convergence tolerance
real(8)                            :: m           ! slope of line
real(8), dimension(:), allocatable :: z           ! CG update iterates
real(8), dimension(:), allocatable :: res         ! solution residual

! variables to define the local problem
integer                               :: n_el        ! number of (local) elements
integer                               :: n_nodes     ! number of (local) nodes
integer                               :: numprocs    ! number of processors
integer                               :: maxperproc  ! maximum number of elements per processor
integer                               :: rank        ! processor rank
integer                               :: ddcnt       ! domain decomposition counter
real(8)                               :: itererror   ! whole-loop iteration error
real(8)                               :: ddtol       ! domain decomposition loop tolerance
integer                               :: ierr        ! holds error state for MPI calls
integer, dimension(:, :), allocatable :: edges       ! nodes on edge of each domain
integer, dimension(:),    allocatable :: recv_displs ! displacement of each domain
real(8), dimension(:),    allocatable :: xel         ! coordinates in each domain
integer, dimension(:),    allocatable :: numnodes    ! number of nodes in each domain
integer, dimension(:),    allocatable :: elems       ! n_el in each domain
real(8), dimension(:),    allocatable :: rglob       ! global load vector
real(8), dimension(:),    allocatable :: a           ! CG solution iterates
integer, dimension(mpi_status_size)   :: stat        ! MPI send/receive status
real(8), dimension(2)                 :: prev        ! previous interface values
real(8), dimension(2)                 :: BClocals    ! values of BCs for each interface
real(8), dimension(2)                 :: BCvals      ! values of BCs for each domain

! variables to define the coarse-mesh solution
real(8), dimension(:),    allocatable :: hlocal      ! coarse element lengths
real(8), dimension(:),    allocatable :: zcoarse     ! coarse CG vector
real(8), dimension(:),    allocatable :: rescoarse   ! coarse CG residual vector 
real(8), dimension(:),    allocatable :: rglobcoarse ! global load vector, coarse mesh
real(8), dimension(:),    allocatable :: xcoarse     ! coordinates of the shared nodes
real(8), dimension(:),    allocatable :: acoarse     ! coarse-mesh solution
real(8), dimension(:, :), allocatable :: BCcoarse    ! coarse solution BCs
integer, dimension(:, :), allocatable :: LMcoarse    ! location matrix of coarse problem

k = 1.0        ! thermal conductivity
source = 10.0  ! heat source
tol = 0.0001   ! CG convergence tolerance
ddtol = 0.0001 ! domain decomposition loop tolerance

call cpu_time(start)

! initialize variables that are known to all MPI processes
call commandline(n_el, length, leftBC, rightBC)   ! parse command line args
call initialize(h, x, n_el, n_nodes)              ! initialize problem vars
call quadrature(n_qp)                             ! initialize quadrature
call phi_val(qp)                                  ! initialize shape functions
call elementalmatrices()                          ! form elemental matrices and vectors

n_nodes_global = n_nodes
n_el_global = n_el

! initialize the parallel MPI environment
call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world, numprocs, ierr)
call mpi_comm_rank(mpi_comm_world, rank, ierr)

! only the rank 0 process knows about the whole solution vector
if (rank == 0) then
  allocate(soln(n_nodes_global), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of soln array failed."
  soln = 0.0
end if

call allocatedecomp()                     ! allocate space for decompsition data
call initializedecomp()                   ! initialize domain decomposition


! perform a coarse solution to get initial guesses for the interface values
! this only needs to be performed by the rank 0 process
if (rank == 0) then
  n_el    = numprocs
  n_nodes = numprocs + 1
  
  ! allocate storage for the coarse-mesh CG solve
  allocate(xcoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of xcoarse array failed."
  allocate(hlocal(n_el), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of hlocal array failed."
  allocate(LMcoarse(2, n_el), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of LMcoarse array failed."
  allocate(rglobcoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rglobcoarse array failed."
  allocate(acoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of acoarse array failed."
  allocate(zcoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of zcoarse array failed."
  allocate(rescoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rescoarse array failed."

  xcoarse(1) = x(edges(1, 1))
  do i = 2, n_nodes
    xcoarse(i) = x(edges(2, i - 1))
    hlocal(i - 1) = xcoarse(i) - xcoarse(i - 1)
  end do

  BCs    = (/ 1, n_nodes /) 
  BCvals = (/ leftBC, rightBC /)

  call locationmatrix(LMcoarse, n_el)      ! form the location matrix
  
  ! form the global load vector
  rglobcoarse = 0.0
  do q = 1, n_el
    do i = 1, 2
      rglobcoarse(LMcoarse(i, q)) = rglobcoarse(LMcoarse(i, q)) &
                                    + hlocal(q) * hlocal(q) * rel(i) / (h * h)
    end do
  end do
  rglobcoarse(BCs(1)) = BCvals(1)
  rglobcoarse(BCs(2)) = BCvals(2)   

  ! initial guess is a straight line between the two endpoints
  m = (rightBC - leftBC) / length
  acoarse = m * (xcoarse - xcoarse(1)) + BCvals(1)

  call conjugategradient(kel, acoarse, LMcoarse, rglobcoarse, zcoarse, rescoarse, BCs)
  
  ! insert first-guess BCs into BCcoarse array, then broadcast to all processes
  BCcoarse(1, 1) = acoarse(1)
  BCcoarse(2, numprocs) = acoarse(n_nodes)

  do i = 1, numprocs - 1
    BCcoarse(2, i) = acoarse(i + 1)
    BCcoarse(1, i + 1) = acoarse(i + 1)
  end do

  deallocate(xcoarse, LMcoarse, rglobcoarse, acoarse, hlocal)  
end if

call mpi_bcast(BCcoarse, 2 * numprocs, mpi_real8, 0, mpi_comm_world, ierr)

! Specify the domain decomposition parameters for each separate domain
n_el    = elems(rank + 1)
n_nodes = numnodes(rank + 1)

call allocateDDdata()

xel       = x(edges(1, rank + 1):edges(2, rank + 1))
BCs       = (/1, n_nodes/)
BCvals(1) = BCcoarse(1, rank + 1)
BCvals(2) = BCcoarse(2, rank + 1)

call locationmatrix(LM, n_el)            ! form the location matrix
call globalload()                              ! form the global load vector

! initial guess is a straight line between the two endpoints
m = (BCvals(2) - BCvals(1)) / (xel(n_nodes) - xel(1))
a = m * (xel - xel(1)) + BCvals(1)

ddcnt = 0
do while (itererror > ddtol)
  ! save the previous values of the interface BCs
  prev = BCvals

  ! each processor solves for its domain ------------------------------------------
  rglob(BCs(1)) = BCvals(1)
  rglob(BCs(2)) = BCvals(2)   
 
  call conjugategradient(kel, a, LM, rglob, z, res, BCs)
  
  ! each processor sends a boundary value to the processor to the right -----------
  if (rank /= numprocs - 1) then
    call mpi_send(a(n_nodes - 1), 1, mpi_real8, rank + 1, rank, mpi_comm_world, ierr)
  end if

  ! processor to the right receives the message -----------------------------------
  if (rank /= 0) then
    call mpi_recv(BClocals(1), 1, mpi_real8, rank - 1, rank - 1, mpi_comm_world, stat, ierr)
    ! assign other local boundary condition
    BClocals(2) = a(2)
  end if

  ! each processor solves its interface problem -----------------------------------
  if (rank /= 0) then
    BCvals(1) = (rel(2) + rel(1) - kel(2, 1) * BClocals(2) &
                      - kel(1, 2) * BClocals(1)) / (kel(2, 2) + kel(1, 1))
    ! send new interface result to rank - 1 process (to BCvals(2) of rank - 1)
    call mpi_send(BCvals(1), 1, mpi_real8, rank - 1, rank, mpi_comm_world, ierr)
  end if

  ! rank - 1 process receives from the process to the right -----------------------
  if (rank /= numprocs - 1) then
    call mpi_recv(BCvals(2), 1, mpi_real8, rank + 1, rank + 1, mpi_comm_world, stat, ierr)
  end if

  ! compute iteration error to determine whether to continue looping -------------- 
  call mpi_barrier(mpi_comm_world, ierr)

  call mpi_allreduce(abs(BCvals(2) - prev(2) + BCvals(1) - prev(1)), itererror, 1, & 
                     mpi_real8, mpi_sum, mpi_comm_world, ierr)  
 
  call mpi_barrier(mpi_comm_world, ierr)
  ddcnt = ddcnt + 1
end do ! ends outermost domain decomposition loop

! each processor broadcasts its final solution to the rank 0 process --------------
call mpi_gatherv(a(1:(n_nodes - 1)), n_nodes - 1, mpi_real8, &
                    soln(1:(n_nodes_global - 1)), elems, recv_displs, mpi_real8, &
                    0, mpi_comm_world, ierr)

! write results to output file ----------------------------------------------------
if (rank == 0) then
  soln(n_nodes_global) = rightBC 
  ! write to an output file. If this file exists, it will be re-written.
  open(1, file='output.txt', iostat=AllocateStatus, status="replace")
  if (AllocateStatus /= 0) STOP "output.txt file opening failed."

  write(1, *) numprocs, n_el_global, ddcnt, soln
end if

! final timing results ----------------------------------------------------------
if (rank == 0) then
  call cpu_time(finish)
  print *, 'P: ', numprocs, 'n_el: ', n_el_global, 'runtime: ', finish - start
end if

! deallocate memory -------------------------------------------------------------
deallocate(xel, LM, rglob, a, z, res)
deallocate(numnodes, elems, edges, recv_displs, BCcoarse)

if (rank == 0) deallocate(soln)

call mpi_finalize(ierr)

deallocate(x, qp, wt, phi, dphi)

! ------------------------------------------------------------------------------


CONTAINS ! define all internal procedures

subroutine allocateDDdata()
  implicit none
  allocate(xel(numnodes(rank + 1)), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of xel array failed."
  allocate(LM(2, n_el), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of LM array failed."
  allocate(rglob(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rglob array failed."
  allocate(a(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of a array failed."
  allocate(z(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of z array failed."
  allocate(res(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of res array failed."
end subroutine allocateDDdata


subroutine allocatedecomp()
  implicit none
  ! each process has their own copy of this DD information
  allocate(numnodes(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of numnodes array failed."
  allocate(elems(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of elems array failed."
  allocate(edges(2, numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of edges array failed."
  allocate(recv_displs(numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of recv_displs array failed."
  allocate(BCcoarse(2, numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of BCcoarse array failed."
end subroutine allocatedecomp


subroutine initializedecomp()
  implicit none

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
    numnodes(j) = elems(j) * 2 - (elems(j) - 1)
  end do
  
  ! assign the global node numbers that correspond to the edges of each domain
  edges(:, 1) = (/1, elems(1) * 2 - (elems(1) - 1)/)
  do i = 2, numprocs
    edges(:, i) = (/edges(2, i - 1), edges(2, i - 1) + elems(i) * 2 - elems(i) /)
  end do
  
  ! assign an initial itererror (dummy value to enter the loop)
  itererror = 1

  recv_displs = 0
  do i = 2, numprocs
    recv_displs(i) = recv_displs(i - 1) + elems(i - 1)
  end do
end subroutine initializedecomp



subroutine globalload()
  implicit none
  rglob = 0.0
  do q = 1, n_el
      rglob(LM(:, q)) = rglob(LM(:, q)) + rel(:)
  end do
end subroutine globalload


subroutine elementalmatrices()
  implicit none

  kel = 0.0
  rel = 0.0
  do q = 1, n_qp
    do i = 1, 2
      rel(i) = rel(i) + wt(q) * source * phi(i, q) * h * h / 2.0
      
      kel(i, 1) = kel(i, 1) + wt(q) * dphi(i, q) * k * dphi(1, q) * 2.0
      kel(i, 2) = kel(i, 2) + wt(q) * dphi(i, q) * k * dphi(2, q) * 2.0
    end do
  end do
end subroutine elementalmatrices


subroutine conjugategradient(kel, a, LM, rglob, z, res, BCs)
! solves K * a = rglob using the conjugate gradient method

  implicit none
  real(8), intent(inout) :: a(:)         ! resultant vector
  real(8), intent(inout) :: z(:), res(:) ! CG vectors
  real(8), intent(in)    :: rglob(:)     ! rhs vector
  real(8), intent(in)    :: kel(:, :)    ! elementary stiffness matrix
  integer, intent(in)    :: LM(:, :)     ! location matrix for multiplication
  integer, intent(in)    :: BCs(:)       ! BC nodes 
  
  ! local variables
  real(8) :: lambda, theta
  integer :: cnt

  res         = rglob - sparse_mult(kel, LM, a)
  z           = res
  lambda      = dotprod(z, res)/dotprod(z, sparse_mult(kel, LM, z))
  a           = a + lambda * z
  res         = rglob - sparse_mult(kel, LM, a)
  convergence = 0.0
 
  do i = 1, size(res)
    convergence = convergence + abs(res(i))
  end do
  
  cnt = 0
  do while (convergence > tol)
    theta       = sparse_mult_dot(kel, LM, z, res, BCs) &
                  / sparse_mult_dot(kel, LM, z, z, BCs)
    z           = res - theta * z
    lambda      = dotprod(z, res) / sparse_mult_dot(kel, LM, z, z, BCs)
    a           = a + lambda * z
    res         = rglob - sparse_mult(kel, LM, a)
    convergence = 0.0
    
    do i = 1, size(res)
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
  real(8)  :: vec1(:), vec2(:)
  integer  :: i

  dotprod = 0.0  
  do i = 1, size(vec1)
    dotprod = dotprod + vec1(i) * vec2(i)
  end do
end function dotprod


function sparse_mult_dot(matrix, LM, vector, vecdot, BCs)
  implicit none
  real(8) :: matrix(:, :) ! elementary matrix
  real(8) :: vector(:)    ! full vector
  real(8) :: vecdot(:)    ! vector to take dot product with
  integer :: LM(:, :)     ! location matrix
  integer :: BCs(:)       ! list of BC nodes
 
  ! return value of function
  real(8) :: sparse_mult_dot

  integer  :: i, j, q ! looping variables
 
  sparse_mult_dot = 0.0
  do q = 1, n_el ! loop over the elements
    do i = 1, 2 ! loop over all entries in kel
      if (any(BCs == LM(i, q))) then ! apply boundary conditions
          sparse_mult_dot = sparse_mult_dot + & 
                            vecdot(LM(i, q)) * kronecker(LM(i, q), LM(1, q)) * vector(LM(1, q)) + &
                            vecdot(LM(i, q)) * kronecker(LM(i, q), LM(2, q)) * vector(LM(2, q))
      else
        ! implicitly assumes that the matrix is symmetric (ok for this application)
          sparse_mult_dot = sparse_mult_dot + &
                            vecdot(LM(i, q)) * matrix(1, i) * vector(LM(1, q)) + &
                            vecdot(LM(i, q)) * matrix(2, i) * vector(LM(2, q))
      end if
    end do
  end do
end function sparse_mult_dot


function sparse_mult(matrix, LM, vector)
  implicit none
  real(8) :: matrix(:, :) ! elementary matrix (assumed-shape array)
  real(8) :: vector(:)    ! full vector (assumed-shape array)
  integer  :: LM(:, :)     ! location matrix
 
  ! return value of function, as an automatic array
  real(8) :: sparse_mult(size(vector))
  
  integer :: i, j, q ! looping variables
  sparse_mult = 0.0
   
  do q = 1, n_el ! loop over the elements
    do i = 1, 2 ! loop over all entries in kel
      if (any(BCs == LM(i, q))) then 
          ! diagonal terms set to 1.0, off-diagonal set to 0.0
          sparse_mult(LM(i, q)) = sparse_mult(LM(i, q)) + &
                                  kronecker(LM(i, q), LM(1, q)) * vector(LM(1, q)) + &
                                  kronecker(LM(i, q), LM(2, q)) * vector(LM(2, q))
      else
          sparse_mult(LM(i, q)) = sparse_mult(LM(i, q)) + &
                                  matrix(i, 1) * vector(LM(1, q)) + &
                                  matrix(i, 2) * vector(LM(2, q))
      end if
    end do
  end do
end function sparse_mult


subroutine locationmatrix(LM, n_el)
  implicit none
  integer, intent(inout) :: LM(:, :)
  integer, intent(in)    :: n_el
  integer                :: j
  
  do j = 1, n_el
    LM(:, j) = (/ j, j + 1 /)
  end do
end subroutine locationmatrix

subroutine phi_val(qp)
! populate phi and dphi, which are global to the calling program
  implicit none
  real(8), intent(in)  :: qp(:)

  allocate(phi(2, n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of phi array failed."
  allocate(dphi(2, n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of dphi array failed."

  phi(1, :)  = (1.0 - qp(:)) / 2.0
  phi(2, :)  = (1.0 + qp(:)) / 2.0
  dphi(1, :) = -1.0 / 2.0
  dphi(2, :) =  1.0 / 2.0
end subroutine phi_val


subroutine quadrature(n_qp)
  implicit none
  integer, intent(out) :: n_qp
  
  n_qp = 2

  allocate(qp(n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of qp array failed."
  allocate(wt(n_qp), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of wt array failed."
  
  qp = (/ -1.0/sqrt(3.0), 1.0/sqrt(3.0) /)
  wt = (/ 1.0, 1.0 /)
end subroutine quadrature


subroutine commandline(n_el, length, leftBC, rightBC)
  implicit none
  integer, intent(out)  :: n_el
  real(8), intent(out) :: length
  real(8), intent(out) :: leftBC
  real(8), intent(out) :: rightBC
 
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


subroutine initialize(h, x, n_el, n_nodes)
  implicit none
  real(8), intent(out) :: h 
  integer, intent(in)  :: n_el    
  integer, intent(out) :: n_nodes 
  real(8), dimension(:), allocatable, intent(out) :: x 

  integer :: i ! looping variable

  h = length / real(n_el)
  n_nodes = n_el + 1 

  ! allocate memory for the vector of node coordinates
  allocate(x(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of x array failed."

  do i = 1, size(x)
    x(i) = real(i - 1) * h
  end do
end subroutine initialize


END PROGRAM main
