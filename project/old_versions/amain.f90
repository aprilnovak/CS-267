! alternative Schwarz method, where the domains initially overlap

PROGRAM main

implicit none
include 'mpif.h'

! variables for overall execution
integer  :: AllocateStatus     ! variable to hold memory allocation success
integer  :: i, j, q            ! loop iteration variables
real(8) :: start               ! holds start run time
real(8) :: finish              ! holds end run time

! variables to define the global problem
integer  :: n_el               ! number of (local) elements
integer  :: n_el_global        ! number of (global) elements
integer  :: n_en               ! number of nodes per element
integer  :: n_nodes            ! number of (local) nodes
integer  :: n_nodes_global     ! number of (global) nodes
integer  :: order              ! polynomial order
integer  :: n_qp               ! number of quadrature points
real(8) :: length              ! length of the domain (1-D)
real(8) :: h                   ! length of one element
real(8) :: k                   ! thermal conductivity
real(8) :: source              ! uniform heat source
real(8) :: leftBC              ! left Dirichlet boundary condition value
real(8) :: rightBC             ! right Dirichlet boundary condition value

! variables to define the CG solver
integer  :: cnt                ! number of CG iterations
real(8) :: theta               ! CG coefficient
real(8) :: lambda              ! CG coefficient
real(8) :: convergence         ! difference between CG iterations
real(8) :: tol                 ! CG convergence tolerance
real(8) :: startCG             ! start, CG
real(8) :: endCG               ! end, CG
real(8) :: m                   ! slope of line

! parallel variables
integer  :: numprocs           ! number of processors
integer  :: maxperproc         ! maximum number of elements per processor
integer  :: rank               ! processor rank
integer  :: face               ! ID number of interface problem
integer  :: iter               ! domain decomposition solution iteration counter
integer  :: leftnode           ! node number at left of interface layer
integer  :: rightnode          ! node number at right of interface layer
integer  :: ddcnt              ! domain decomposition counter
real(8) :: itererror           ! whole-loop iteration error
real(8) :: ddtol               ! domain decomposition loop tolerance
integer :: ierr                ! holds error state for MPI calls
integer :: overlap             ! number of elements to overlap each domain
integer :: sendlength          ! rank's length of solution to plot

integer, dimension(mpi_status_size) :: stat ! MPI send/receive status

integer, dimension(:), allocatable :: recv_displs ! displacement of each domain
real(8), dimension(:), allocatable :: prev        ! previous interface values
real(8), dimension(:), allocatable :: soln        ! global solution vector
real(8), dimension(:), allocatable :: xel         ! coordinates in each domain
integer, dimension(:), allocatable :: numnodes    ! number of nodes in each domain
integer, dimension(:), allocatable :: elemspre    ! n_el in each domain before overlap
integer, dimension(:), allocatable :: elems       ! n_el in each domain
integer, dimension(2)              :: BCs         ! boundary condition nodes
real(8), dimension(:), allocatable :: qp          ! quadrature points
real(8), dimension(:), allocatable :: wt          ! quadrature weights
real(8), dimension(:), allocatable :: x           ! coordinates of the nodes
real(8), dimension(:), allocatable :: rel         ! elemental load vector
real(8), dimension(:), allocatable :: rglob       ! global load vector
real(8), dimension(:), allocatable :: a           ! CG solution iterates
real(8), dimension(:), allocatable :: z           ! CG update iterates
real(8), dimension(:), allocatable :: res         ! solution residual
real(8), dimension(:), allocatable :: send        ! final solution values to send

real(8), dimension(2)                 :: BCvals   ! values of BCs for each domain
real(8), dimension(:, :), allocatable :: kel      ! elemental stiffness matrix
real(8), dimension(:, :), allocatable :: phi      ! shape functions
real(8), dimension(:, :), allocatable :: dphi     ! shape function derivatives
integer, dimension(:, :), allocatable :: edges    ! nodes on edge of each domain
integer, dimension(:, :), allocatable :: LM       ! location matrix

k = 1.0        ! thermal conductivity
source = 10.0  ! heat source
tol = 0.0001   ! CG convergence tolerance
ddtol = 0.0005 ! domain decomposition loop tolerance
order = 1      ! only works for linear elements

call cpu_time(start)

! initialize variables that are known to all MPI processes
call commandline(n_el, length, leftBC, rightBC)   ! parse command line args
call initialize(h, x, n_en, n_el, order, n_nodes) ! initialize problem vars
call quadrature(order, n_qp)                      ! initialize quadrature
call phi_val(order, qp)                           ! initialize shape functions
call elementalmatrices()                          ! form elemental matrices and vectors

n_nodes_global = n_nodes
n_el_global    = n_el

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

! each process has their own copy of this DD information
allocate(numnodes(numprocs), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of numnodes array failed."
allocate(elems(numprocs), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of elems array failed."
allocate(elemspre(numprocs), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of elemspre array failed."
allocate(edges(2, numprocs), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of edges array failed."
allocate(recv_displs(numprocs), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of recv_displs array failed."
allocate(prev(numprocs), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of prev array failed."
  
overlap = 1 ! number of elements to overlap each domain with the next (half-overlap)

call initializedecomp()                   ! initialize domain decomposition

if (rank == 0) then
  print *, 'numnodes: ', numnodes
  print *, 'elems: ', elems
  print *, 'elemspre: ', elemspre
  print *, 'edges: ', edges
  print *, 'recv_displs: ', recv_displs
end if

allocate(xel(numnodes(rank + 1)), stat = AllocateStatus)
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

n_el    = elems(rank + 1)
n_nodes = numnodes(rank + 1)

! determine the length of my solution that is saved for plotting
allocate(send(sendlength), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of send array failed."

! assign the initial boundary conditions given the rank of the calling process
! only the left processors do anything the first iteration
!if (mod(rank + 1, 2) /= 0) then ! odd processes
  m = (rightBC - leftBC) / length
  BCvals(1) = m * x(edges(1, rank + 1)) + leftBC
  BCvals(2) = m * x(edges(2, rank + 1)) + leftBC
!end if

print *, 'rank: ', rank, 'initial BCs: ', BCvals

! each processor prepares its solution needs  
xel = x(edges(1, rank + 1):edges(2, rank + 1)) ! domain sizes don't change
BCs = (/1, n_nodes/)                           ! BCs are always applied on end nodes
call locationmatrix()                          ! form the location matrix
call globalload()                              ! form the global load vector

ddcnt = 0
do while (itererror > ddtol)
!do while (ddcnt < 10)
  ! save the previous values of the interface BCs
  prev = BCvals

  ! each odd processor solves for its domain ------------------------------------------

  if (mod(rank + 1, 2) /= 0) then
    rglob(BCs(1)) = BCvals(1)
    rglob(BCs(2)) = BCvals(2)   
  
    call conjugategradient(BCvals(2), BCvals(1))
    print *, 'rank: ', rank, 'odd solution: ', a

    ! each odd processor sends a value to the processor to the right ------------------

    if (rank + 1 /= numprocs) then ! no one to send to if last
      ! tag of the message is _rank_
      call mpi_send(a(n_nodes - 2 * overlap), 1, mpi_real8, rank + 1, rank, mpi_comm_world, ierr)
    end if
 
    ! each odd processor sends a value to the processor to the left -------------------

    if (rank + 1 /= 1) then
      ! tag of the message is _rank_ + 1000
      call mpi_send(a(2 * overlap), 1, mpi_real8, rank - 1, rank + 1000, mpi_comm_world, ierr)
    end if

  else ! even processes receive from the left -----------------------------------------
    ! tag of the message is rank - 1, since it is the rank of the sending process
    call mpi_recv(BCvals(1), 1, mpi_real8, rank - 1, rank - 1, mpi_comm_world, stat, ierr)

    ! even processes receive from the right (if last process is even, there is nothing to receive)
    if (rank + 1 /= numprocs) then
      call mpi_recv(BCvals(2), 1, mpi_real8, rank + 1, rank + 1 + 1000, mpi_comm_world, stat, ierr)
    end if
    
   print *, 'rank: ', rank, 'I just received: ', BCvals

  end if

  call mpi_barrier(mpi_comm_world, ierr)

  ! each even processor solves its domain ----------------------------------------------
  if (mod(rank + 1, 2) == 0) then
    rglob(BCs(1)) = BCvals(1)
    rglob(BCs(2)) = BCvals(2)   
  
    call conjugategradient(BCvals(2), BCvals(1))
  
    print *, 'rank: ', rank, 'first even solution: ', a

    ! each even processor sends a value to the right -----------------------------------
    if (rank + 1 /= numprocs) then
      ! tag is _rank_
      call mpi_send(a(n_nodes - 2 * overlap), 1, mpi_real8, rank + 1, rank, mpi_comm_world, ierr)
    end if 

    ! each even processor sends a value to the left ------------------------------------
    call mpi_send(a(2 * overlap), 1, mpi_real8, rank - 1, rank + 2000, mpi_comm_world, ierr)

  else
    ! each odd processor receives a value from the left --------------------------------
    if (rank /= 0) then ! first process has no one to the left
      call mpi_recv(BCvals(1), 1, mpi_real8, rank - 1, rank - 1, mpi_comm_world, stat, ierr)
    end if
   
    ! each odd processor receives a value from the right -------------------------------
    if (rank + 1 /= numprocs) then
      call mpi_recv(BCvals(2), 1, mpi_real8, rank + 1, rank + 1 + 2000, mpi_comm_world, stat, ierr)
    end if
  end if

  ! compute iteration error to determine whether to continue looping -------------- 
  call mpi_barrier(mpi_comm_world, ierr)

  call mpi_allreduce(abs(BCvals(2) - prev(2) + BCvals(1) - prev(1)), itererror, 1, mpi_real8, mpi_sum, &
             mpi_comm_world, ierr)  
 
  call mpi_barrier(mpi_comm_world, ierr)
  ddcnt = ddcnt + 1
end do ! ends outermost domain decomposition loop

call mpi_barrier(mpi_comm_world, ierr)

if (rank == 0) then
  send = a(1:(n_nodes - overlap - 1))
else if (rank == numprocs - 1) then
  send = a((overlap + 1):(n_nodes - 1))
else
  send = a((overlap + 1):(n_nodes - overlap - 1))
end if


if (rank == 0) print *, 'count: ', ddcnt

! each processor broadcasts its final solution to the rank 0 process --------------
call mpi_gatherv(send, sendlength, mpi_real8, &
                    soln(1:(n_nodes_global - 1)), elemspre, recv_displs, mpi_real8, &
                    0, mpi_comm_world, ierr)

! write results to output file ----------------------------------------------------

if (rank == 0) then
  soln(n_nodes_global) = rightBC 
  ! write to an output file. If this file exists, it will be re-written.
  open(1, file='output.txt', iostat=AllocateStatus, status="replace")
  if (AllocateStatus /= 0) STOP "output.txt file opening failed."

  write(1, *) numprocs, n_el_global, ddcnt, soln
  write(1, *) soln(:)
end if

! final timing results ----------------------------------------------------------
if (rank == 0) then
  call cpu_time(finish)
  print *, 'P: ', numprocs, 'n_el: ', n_el, 'runtime: ', finish - start
end if

! deallocate memory -------------------------------------------------------------
deallocate(xel, LM, rglob, a, z, res)
deallocate(numnodes, elems, edges, recv_displs, prev)

if (rank == 0) deallocate(soln)

call mpi_finalize(ierr)

deallocate(x, qp, wt, phi, dphi, kel, rel)

! ------------------------------------------------------------------------------


CONTAINS ! define all internal procedures

subroutine initializedecomp()
  implicit none
  ! distribute the elements among the processors, initially assuming no overlap
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

  if (rank == 0) print *, 'elems: ', elems
  elemspre = elems

  ! number of pieces to send (last entry will need to be manually input as rightBC)
  sendlength = elems(rank + 1)
  
  ! set the displacements according to non-overlapping domains
  recv_displs = 0
  do i = 2, numprocs
    recv_displs(i) = recv_displs(i - 1) + elems(i - 1)
  end do
   
  ! assign the global node numbers that correspond to the edges of each domain
  ! this must take into account that the domains overlap, so this uses the 
  ! element counts before manually adding in the overlap
  edges(:, 1) = (/1, elems(1) * n_en - (elems(1) - 1)/)
  do i = 2, numprocs
    edges(:, i) = (/edges(2, i - 1), edges(2, i - 1) + elems(i) * n_en - elems(i) /)
  end do
  
  edges(2, 1) = edges(2, 1) + overlap;
  edges(1, numprocs) = edges(1, numprocs) - overlap;
  
  do i = 2, numprocs - 1
    edges(1, i) = edges(1, i) - overlap;
    edges(2, i) = edges(2, i) + overlap;
  end do
 
  ! add overlap to the elems() array 
  do i = 1, numprocs
    if ((i == 1) .or. (i == numprocs)) then
      elems(i) = elems(i) + overlap
    else
      elems(i) = elems(i) + 2 * overlap
    end if
  end do
  
  ! assign the numbers of nodes in each domain 
  do j = 1, numprocs
    numnodes(j) = elems(j) * n_en - (elems(j) - 1)
  end do
  
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
  real(8), intent(in) :: rightBC, leftBC

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
    
    if (rank == 0) print *, convergence
 
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
  real(8) :: vec1(:)
  real(8) :: vec2(:)

  integer  :: i ! looping variable

  dotprod = 0.0  
  do i = 1, size(vec1)
    dotprod = dotprod + vec1(i) * vec2(i)
  end do
end function dotprod

function sparse_mult_dot(matrix, LM, vector, vecdot)
  implicit none
  real(8) :: matrix(:, :) ! elementary matrix (assumed-shape array)
  real(8) :: vector(:)    ! full vector (assumed-shape array)
  real(8) :: vecdot(:)
  integer  :: LM(:, :)     ! location matrix
 
  ! return value of function
  real(8) :: sparse_mult_dot
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
  real(8) :: matrix(:, :) ! elementary matrix (assumed-shape array)
  real(8) :: vector(:)    ! full vector (assumed-shape array)
  integer  :: LM(:, :)     ! location matrix
 
  ! return value of function, as an automatic array
  real(8) :: sparse_mult(size(vector))
  
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
  real(8), intent(in)  :: qp(:)

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


subroutine initialize(h, x, n_en, n_el, order, n_nodes)
  implicit none
  real(8), intent(out) :: h 
  integer, intent(out) :: n_en    
  integer, intent(in)  :: n_el    
  integer, intent(in)  :: order   
  integer, intent(out) :: n_nodes 
  real(8), dimension(:), allocatable, intent(out) :: x 

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
