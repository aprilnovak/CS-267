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
real(8), dimension(:),    allocatable :: soln           ! global soln vector

! variables to define the CG solver
real(8)                            :: m           ! slope of line

! variables to define the local problem
integer                               :: numprocs    ! number of processors
integer                               :: rank        ! processor rank
integer                               :: ddcnt       ! DD counter
real(8)                               :: itererror   ! whole-loop iter error
integer                               :: ierr        ! error for MPI calls
real(8), dimension(:),    allocatable :: rglob       ! global load vector
real(8), dimension(:),    allocatable :: a           ! CG solution iterates
integer, dimension(mpi_status_size)   :: stat        ! MPI send/receive status
real(8), dimension(2)                 :: prev        ! prev interface values
real(8), dimension(2)                 :: BClocals    ! BCs for each interface

! variables to define the coarse-mesh solution
real(8), dimension(:),    allocatable :: hlocal        ! coarse element lengths
real(8), dimension(:),    allocatable :: rglobcoarse   ! global load vector
real(8), dimension(:),    allocatable :: acoarse       ! coarse-mesh solution
real(8), dimension(:, :), allocatable :: BCcoarse      ! coarse solution BCs

! variables to define OpenMP thread parallelism
integer :: provided   ! holds provided level of thread support
integer :: omp_get_thread_num, omp_get_num_threads ! OpenMP routines

type(row), allocatable, dimension(:) :: rows       ! rows in global K
type(row), allocatable, dimension(:) :: rowscoarse ! rows in coarse global K
type(LM)                             :: LMfine, LMcoarse

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

allocate(BCcoarse(2, numprocs), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of BCcoarse array failed."

call initialize_domain_decomposition(numprocs)

! perform a coarse solution to get initial guesses for the interface values
! this only needs to be performed by the rank 0 process
if (rank == 0) then
  call initialize_coarse_mesh()
  
  allocate(hlocal(coarse%n_el), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of hlocal array failed."
  allocate(rglobcoarse(coarse%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rglobcoarse array failed."
  allocate(acoarse(coarse%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of acoarse array failed."

  do i = 2, coarse%n_nodes
    hlocal(i - 1) = coarse%x(i) - coarse%x(i - 1)
  end do

  LMcoarse = locationmatrix(coarse%n_el, coarse%n_nodes)
  rglobcoarse = globalvector(LMcoarse, rel, coarse%n_nodes, hlocal)
  
  rglobcoarse(coarse%BCs(1)) = coarse%BCvals(1)
  rglobcoarse(coarse%BCs(2)) = coarse%BCvals(2)   

  rowscoarse = form_csr(LMcoarse, coarse%n_nodes)

  ! initial guess is a straight line between the two endpoints
  acoarse = straightline(coarse)

  call conjugategradient(rowscoarse, acoarse, rglobcoarse, coarse%BCs, reltol = reltol)
  
  ! insert first-guess BCs into BCcoarse array, then broadcast to all processes
  do i = 1, numprocs
    BCcoarse(:, i) = (/acoarse(i), acoarse(i + 1)/)
  end do

  deallocate(LMcoarse%matrix, LMcoarse%cnt, rglobcoarse, acoarse, hlocal)  
  deallocate(rowscoarse)
end if

call mpi_bcast(BCcoarse, 2 * numprocs, mpi_real8, 0, mpi_comm_world, ierr)

! Specify the domain decomposition parameters for each domain -----------------
allocate(rglob(dd(rank + 1)%n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of rglob array failed."
allocate(a(dd(rank + 1)%n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of a array failed."
allocate(rows(dd(rank + 1)%n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of rows array failed."

dd(rank + 1)%BCvals(1) = BCcoarse(1, rank + 1)
dd(rank + 1)%BCvals(2) = BCcoarse(2, rank + 1)

LMfine = locationmatrix(dd(rank + 1)%n_el, dd(rank + 1)%n_nodes)
rglob  = globalvector(LMfine, rel, dd(rank + 1)%n_nodes)
rows   = form_csr(LMfine, dd(rank + 1)%n_nodes)

! initial guess is a straight line between the two endpoints
m = (dd(rank + 1)%BCvals(2) - dd(rank + 1)%BCvals(1)) / &
     (dd(rank + 1)%x(dd(rank + 1)%n_nodes) - dd(rank + 1)%x(1))
a = m * (dd(rank + 1)%x - dd(rank + 1)%x(1)) + dd(rank + 1)%BCvals(1)

itererror = 1
ddcnt     = 0
do while (itererror > ddtol)
  ! save the previous values of the interface BCs
  prev = dd(rank + 1)%BCvals

  ! each processor solves for its domain --------------------------------------
  rglob(dd(rank + 1)%BCs(1)) = dd(rank + 1)%BCvals(1)
  rglob(dd(rank + 1)%BCs(2)) = dd(rank + 1)%BCvals(2)   
  
  call conjugategradient(rows, a, rglob, dd(rank + 1)%BCs, reltol = reltol)
  
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
    dd(rank + 1)%BCvals(1) = (rel(2) + rel(1) - kel(2, 1) * BClocals(2) &
                      - kel(1, 2) * BClocals(1)) / (kel(2, 2) + kel(1, 1))
    ! send new interface result to rank - 1 process (to BCvals(2) of rank - 1)
    call mpi_send(dd(rank + 1)%BCvals(1), 1, mpi_real8, rank - 1, rank, mpi_comm_world, ierr)
  end if

  ! rank - 1 process receives from the process to the right -------------------
  if (rank /= numprocs - 1) then
    call mpi_recv(dd(rank + 1)%BCvals(2), 1, mpi_real8, rank + 1, rank + 1, mpi_comm_world, stat, ierr)
  end if

  ! compute iteration error to determine whether to continue looping ---------- 
  call mpi_barrier(mpi_comm_world, ierr)

  call mpi_allreduce(abs(dd(rank + 1)%BCvals(2) - prev(2) + dd(rank + 1)%BCvals(1) - prev(1)), itererror,&
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
deallocate(LMfine%matrix, LMfine%cnt, rglob, a)
deallocate(BCcoarse)
deallocate(rows)
if (rank == 0) deallocate(soln)

call mpi_finalize(ierr)

call dealloc_domains()
call dealloc_global_mesh()
call dealloc_shapefunctions()
call dealloc_quadset()
! ------------------------------------------------------------------------------


CONTAINS ! define all internal procedures

subroutine conjugategradient(rows, a, rglob, BCs, reltol)
  real(8), intent(inout) :: a(:)          ! resultant vector
  real(8), intent(in)    :: rglob(:)      ! rhs vector
  integer, intent(in)    :: BCs(:)        ! BC nodes 
  type(row), intent(in)  :: rows(:)       ! kglob in CSR form
  real(8), intent(in), optional :: reltol ! relative CG tolerance

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

  a      = a + lambda * z
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
end subroutine conjugategradient

END PROGRAM main
