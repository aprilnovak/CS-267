! non-overlapping interface
! each processor computes one interface problem
! builds from bmain by adding CSR matrix storage

! program to solve the heat equation (Dirichlet boundary conditions only)
! using domain decomposition. This version has the domain to the left of
! each interface send a single interface value to the processor to the
! right of it. This process then computes the interface problem, and 
! sends a single udpated value to the processor to the left.

PROGRAM main

! read input information
use read_data

! define quantities related to the mesh
use mesh

! define the quadrature rule
use quadrature

! define the elemental matrices
use element_matrices

! define the CSR matrix storage
use csr_storage

! provide interface to solvers
use solvers

implicit none

include 'mpif.h'


! variables for overall execution
integer  :: AllocateStatus     ! variable to hold memory allocation success
integer  :: i, j, q            ! loop iteration variables
real(8)  :: start              ! start run time
real(8)  :: finish             ! end run time
real(8)  :: startCG            ! start CG time
real(8)  :: endCG              ! end CG time
real(8)  :: startCSR           ! start CSR time
real(8)  :: endCSR             ! end CSR time

! variables to define the global problem
integer, dimension(2)                 :: BCs            ! BC nodes
real(8), dimension(:),    allocatable :: soln           ! global soln vector

! variables to define the CG solver
integer                            :: cnt         ! number of CG iterations
real(8)                            :: convergence ! difference b/w CG iterations
real(8)                            :: m           ! slope of line
real(8), dimension(:), allocatable :: z           ! CG update iterates
real(8), dimension(:), allocatable :: res         ! solution residual

! variables to define the local problem
integer                               :: n_el        ! number of local elements
integer                               :: n_nodes     ! number of (local) nodes
integer                               :: numprocs    ! number of processors
integer                               :: rank        ! processor rank
integer                               :: ddcnt       ! DD counter
real(8)                               :: itererror   ! whole-loop iter error
integer                               :: ierr        ! error for MPI calls
real(8), dimension(:),    allocatable :: xel         ! coordinates in domain
real(8), dimension(:),    allocatable :: rglob       ! global load vector
real(8), dimension(:),    allocatable :: a           ! CG solution iterates
integer, dimension(mpi_status_size)   :: stat        ! MPI send/receive status
real(8), dimension(2)                 :: prev        ! prev interface values
real(8), dimension(2)                 :: BClocals    ! BCs for each interface
real(8), dimension(2)                 :: BCvals      ! BCs for each domain

! variables to define the coarse-mesh solution
real(8), dimension(:),    allocatable :: hlocal        ! coarse element lengths
real(8), dimension(:),    allocatable :: zcoarse       ! CG vector
real(8), dimension(:),    allocatable :: rescoarse     ! CG residual vector
real(8), dimension(:),    allocatable :: rglobcoarse   ! global load vector
real(8), dimension(:),    allocatable :: xcoarse       ! coords of shared nodes
real(8), dimension(:),    allocatable :: acoarse       ! coarse-mesh solution
real(8), dimension(:, :), allocatable :: BCcoarse      ! coarse solution BCs

! variables to define OpenMP thread parallelism
integer :: numthreads ! number of OpenMP threads
integer :: mythread   ! current thread number
integer :: provided   ! holds provided level of thread support
integer :: omp_get_thread_num, omp_get_num_threads ! OpenMP routines

type(row), allocatable, dimension(:) :: rows       ! rows in global K
type(row), allocatable, dimension(:) :: rowscoarse ! rows in coarse global K

type(LM) :: LMfine
type(LM) :: LMcoarse

call cpu_time(start)

! initialize variables that are known to all MPI processes
! instead of having one process compute these and then broadcast

! read information from namelist
call read_namelist()

! read information from the command line
call read_commandline()

! determine global mesh quantities
call initialize_global_mesh()

! define the quadrature rule
call define_quadset()

! define the shape functions
call define_shapefunctions()

! determine the elemental load vector
call elementalload()

! determine the elemental stiffness matrix
call elementalstiffness()

! initialize the parallel MPI environment with mpi_thread_single thread support
call mpi_init_thread(0, provided, ierr)
call mpi_comm_size(mpi_comm_world, numprocs, ierr)
call mpi_comm_rank(mpi_comm_world, rank, ierr)

! only the rank 0 process knows about the whole solution vector
if (rank == 0) then
  allocate(soln(global%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of soln array failed."
  soln = 0.0
end if

call allocatedecomp()           ! allocate space for decompsition data
call initialize_domain_decomposition(numprocs)

! perform a coarse solution to get initial guesses for the interface values
! this only needs to be performed by the rank 0 process
if (rank == 0) then
  n_el    = numprocs
  n_nodes = numprocs + 1
  
  call allocatecoarse()

  xcoarse(1) = global%x(domains%edges(1, 1))
  do i = 2, n_nodes
    xcoarse(i) = global%x(domains%edges(2, i - 1))
    hlocal(i - 1) = xcoarse(i) - xcoarse(i - 1)
  end do

  BCs    = (/ 1, n_nodes /) 
  BCvals = (/ leftBC, rightBC /)

  LMcoarse = locationmatrix(n_el, n_nodes)

  ! form the global load vector
  rglobcoarse = 0.0
  do q = 1, n_el
    do i = 1, 2
      rglobcoarse(LMcoarse%matrix(i, q)) = rglobcoarse(LMcoarse%matrix(i, q)) &
                                    + hlocal(q) * hlocal(q) * rel(i) / (global%h * global%h)
    end do
  end do
  rglobcoarse(BCs(1)) = BCvals(1)
  rglobcoarse(BCs(2)) = BCvals(2)   

  rowscoarse = form_csr(LMcoarse%matrix, LMcoarse%cnt, n_nodes)

  ! initial guess is a straight line between the two endpoints
  m = (rightBC - leftBC) / global%length
  acoarse = m * (xcoarse - xcoarse(1)) + BCvals(1)

  call conjugategradient(rowscoarse, acoarse, rglobcoarse, zcoarse, rescoarse, BCs, reltol = reltol)
  
  ! insert first-guess BCs into BCcoarse array, then broadcast to all processes
  BCcoarse(1, 1) = acoarse(1)
  BCcoarse(2, numprocs) = acoarse(n_nodes)

  do i = 1, numprocs - 1
    BCcoarse(2, i) = acoarse(i + 1)
    BCcoarse(1, i + 1) = acoarse(i + 1)
  end do

  deallocate(xcoarse, LMcoarse%matrix, LMcoarse%cnt, rglobcoarse, acoarse, hlocal)  
  deallocate(rowscoarse)
end if

call mpi_bcast(BCcoarse, 2 * numprocs, mpi_real8, 0, mpi_comm_world, ierr)

! Specify the domain decomposition parameters for each domain -----------------
!n_el    = domains%elems(rank + 1)
!n_nodes = domains%numnodes(rank + 1)

call allocateDDdata()

xel       = global%x(domains%edges(1, rank + 1):domains%edges(2, rank + 1))
BCs       = (/1, dd(rank + 1)%n_nodes/)
BCvals(1) = BCcoarse(1, rank + 1)
BCvals(2) = BCcoarse(2, rank + 1)

LMfine = locationmatrix(dd(rank + 1)%n_el, dd(rank + 1)%n_nodes)
rglob = globalload(LMfine%matrix, rel, dd(rank + 1)%n_el, dd(rank + 1)%n_nodes)      
rows = form_csr(LMfine%matrix, LMfine%cnt, dd(rank + 1)%n_nodes)

! initial guess is a straight line between the two endpoints
m = (BCvals(2) - BCvals(1)) / (xel(n_nodes) - xel(1))
a = m * (xel - xel(1)) + BCvals(1)

itererror = 1
ddcnt = 0
do while (itererror > ddtol)
  ! save the previous values of the interface BCs
  prev = BCvals

  ! each processor solves for its domain --------------------------------------
  rglob(BCs(1)) = BCvals(1)
  rglob(BCs(2)) = BCvals(2)   
  
  call conjugategradient(rows, a, rglob, z, res, BCs, reltol = reltol)
  
  ! each processor sends a boundary value to the processor to the right -------
  if (rank /= numprocs - 1) then
    call mpi_send(a(dd(rank + 1)%n_nodes - 1), 1, mpi_real8, rank + 1, rank, mpi_comm_world, ierr)
  end if

  ! processor to the right receives the message -------------------------------
  if (rank /= 0) then
    call mpi_recv(BClocals(1), 1, mpi_real8, rank - 1, rank - 1, mpi_comm_world, stat, ierr)
    ! assign other local boundary condition
    BClocals(2) = a(2)
  end if

  call mpi_barrier(mpi_comm_world, ierr)

  ! each processor solves its interface problem -------------------------------
  if (rank /= 0) then
    BCvals(1) = (rel(2) + rel(1) - kel(2, 1) * BClocals(2) &
                      - kel(1, 2) * BClocals(1)) / (kel(2, 2) + kel(1, 1))
    ! send new interface result to rank - 1 process (to BCvals(2) of rank - 1)
    call mpi_send(BCvals(1), 1, mpi_real8, rank - 1, rank, mpi_comm_world, ierr)
  end if

  ! rank - 1 process receives from the process to the right -------------------
  if (rank /= numprocs - 1) then
    call mpi_recv(BCvals(2), 1, mpi_real8, rank + 1, rank + 1, mpi_comm_world, stat, ierr)
  end if

  ! compute iteration error to determine whether to continue looping ---------- 
  call mpi_barrier(mpi_comm_world, ierr)

  call mpi_allreduce(abs(BCvals(2) - prev(2) + BCvals(1) - prev(1)), itererror,&
                     1, mpi_real8, mpi_sum, mpi_comm_world, ierr)  
 
  call mpi_barrier(mpi_comm_world, ierr)
  ddcnt = ddcnt + 1
end do ! ends outermost domain decomposition loop

! each processor broadcasts its final solution to the rank 0 process ----------
call mpi_gatherv(a(1:(dd(rank + 1)%n_nodes - 1)), dd(rank + 1)%n_nodes - 1, mpi_real8, &
                 soln(1:(global%n_nodes - 1)), domains%elems, domains%recv_displs, mpi_real8, &
                 0, mpi_comm_world, ierr)

! write results to output file ------------------------------------------------
if (rank == 0) then
  soln(global%n_nodes) = rightBC 
  ! write to an output file. If this file exists, it will be re-written.
  open(1, file='output.txt', iostat=AllocateStatus, status="replace")
  if (AllocateStatus /= 0) STOP "output.txt file opening failed."

  write(1, *) numprocs, global%n_el, ddcnt, soln
end if

! final timing results --------------------------------------------------------
if (rank == 0) then
  call cpu_time(finish)
  print *, 'P: ', numprocs, 'n_el: ', global%n_el, 'runtime: ', finish - start
end if

! deallocate memory -----------------------------------------------------------
deallocate(xel, LMfine%matrix, LMfine%cnt, rglob, a, z, res)
deallocate(BCcoarse)
!deallocate(numnodes, elems, edges, recv_displs, BCcoarse)
deallocate(rows)
if (rank == 0) deallocate(soln)

call mpi_finalize(ierr)

call dealloc_domains()
call dealloc_global_mesh()
call dealloc_shapefunctions()
call dealloc_quadset()
! ------------------------------------------------------------------------------


CONTAINS ! define all internal procedures

subroutine allocatecoarse()
  implicit none
  allocate(xcoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of xcoarse array failed."
  allocate(hlocal(n_el), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of hlocal array failed."
  allocate(rglobcoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rglobcoarse array failed."
  allocate(acoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of acoarse array failed."
  allocate(zcoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of zcoarse array failed."
  allocate(rescoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rescoarse array failed."
end subroutine allocatecoarse


subroutine allocateDDdata()
  implicit none
  allocate(xel(dd(rank + 1)%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of xel array failed."
  allocate(rglob(dd(rank + 1)%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rglob array failed."
  allocate(a(dd(rank + 1)%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of a array failed."
  allocate(z(dd(rank + 1)%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of z array failed."
  allocate(res(dd(rank + 1)%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of res array failed."
  allocate(rows(dd(rank + 1)%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rows array failed."
end subroutine allocateDDdata


subroutine allocatedecomp()
  implicit none
  allocate(BCcoarse(2, numprocs), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of BCcoarse array failed."
end subroutine allocatedecomp


function globalload(LM, rel, n_el, n_nodes)
  implicit none
  integer, intent(in)    :: LM(:, :)
  real(8), intent(in)    :: rel(:)
  integer, intent(in)    :: n_el, n_nodes
  real(8)                :: globalload(n_nodes)
  integer                :: q

  globalload = 0.0
  do q = 1, n_el
      globalload(LM(:, q)) = globalload(LM(:, q)) + rel(:)
  end do
end function globalload


subroutine conjugategradient(rows, a, rglob, z, res, BCs, reltol)
  implicit none
  real(8), intent(inout) :: a(:)          ! resultant vector
  real(8), intent(inout) :: z(:), res(:)  ! CG vectors
  real(8), intent(in)    :: rglob(:)      ! rhs vector
  integer, intent(in)    :: BCs(:)        ! BC nodes 
  type(row), intent(in)  :: rows(:)       ! kglob in CSR form
  real(8), intent(in), optional :: reltol ! relative CG tolerance

  ! local variables
  real(8) :: lambda, theta, internal, tol
  integer :: cnt, n
  
  n = size(a)

  res         = rglob - csr_mult(rows, a, BCs)

  internal = 0.0
  do i = 1, n
    internal = internal + abs(res(i))
  end do
  
  ! set relative tolerance for convergence, using 0.001 as default
  if (present(reltol)) then
    tol = reltol * internal  
  else
    tol = 0.001 * internal
  end if

  z           = res
  lambda      = dotprod(z, res) / csr_mult_dot(rows, z, BCs, z)
  a           = a + lambda * z
  res         = rglob - csr_mult(rows, a, BCs)
  convergence = 0.0
 
  do i = 1, n
    convergence = convergence + abs(res(i))
  end do

  cnt = 0
  do while (convergence > tol)
    theta       = csr_mult_dot(rows, z, BCs, res) / csr_mult_dot(rows, z, BCs, z)
    z           = res - theta * z
    lambda      = dotprod(z, res) / csr_mult_dot(rows, z, BCs, z)
    a           = a + lambda * z
    res         = rglob - csr_mult(rows, a, BCs)
    convergence = 0.0
    
    do i = 1, n
      convergence = convergence + abs(res(i))
    end do
    
  cnt = cnt + 1
  end do
end subroutine conjugategradient


!real function dotprod(vec1, vec2)
!  implicit none
!  real(8)  :: vec1(:), vec2(:)
!  integer  :: i, n
!  
!  n = size(vec1)
!  dotprod = 0.0  
!
!  do i = 1, n
!    dotprod = dotprod + vec1(i) * vec2(i)
!  end do
!end function dotprod


function csr_mult_dot(rows, a, BCs, vector)
  implicit none
  type(row) :: rows(:)
  real(8)   :: a(:)
  integer   :: BCs(:)
  real(8)   :: vector(:)
  integer   :: n, i, j
  real(8)   :: accum
  real(8)   :: temp(size(a))

  ! return value of function
  real(8)   :: csr_mult_dot 

  n = size(a)
 
  do i = 1, n
    accum = 0.0
    do j = 1, size(rows(i)%columns(:))
      accum = accum + rows(i)%values(j) * a(rows(i)%columns(j))
    end do  
    temp(i) = accum * vector(i)
  end do
 
  ! apply boundary conditions
  do i = 1, size(BCs)
    temp(BCs(i)) = a(BCs(i)) * vector(BCs(i))
  end do

  csr_mult_dot = 0.0
  do i = 1, n
    csr_mult_dot = csr_mult_dot + temp(i)
  end do
end function csr_mult_dot


function csr_mult(rows, a, BCs)
  implicit none
  type(row) :: rows(:)
  real(8)   :: a(:)
  integer   :: BCs(:)
  integer   :: n, i, j
  real(8)   :: accum
  
  ! return value of function, as an automatic array
  real(8)   :: csr_mult(size(a))  

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


END PROGRAM main
