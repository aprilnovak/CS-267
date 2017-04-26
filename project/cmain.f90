! non-overlapping interface
! each processor computes one interface problem
! builds from bmain by adding CSR matrix storage

! program to solve the heat equation (Dirichlet boundary conditions only)
! using domain decomposition. This version has the domain to the left of
! each interface send a single interface value to the processor to the
! right of it. This process then computes the interface problem, and 
! sends a single udpated value to the processor to the left.

PROGRAM main

! must precede any implicit statements
!use matrices, only : hello

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
integer                               :: n_el_global    ! global elements
integer                               :: n_nodes_global ! global nodes
integer                               :: n_qp           ! number of quad points
real(8)                               :: length         ! length of the domain
real(8)                               :: h              ! length of one element
real(8)                               :: k              ! thermal conductivity
real(8)                               :: source         ! uniform heat source
real(8)                               :: leftBC         ! left Dirichlet BC
real(8)                               :: rightBC        ! right Dirichlet BC
integer, dimension(2)                 :: BCs            ! BC nodes
real(8), dimension(2, 2)              :: kel            ! element stiffness mat
real(8), dimension(2)                 :: rel            ! elemental load vector
real(8), dimension(:),    allocatable :: soln           ! global soln vector
real(8), dimension(2)                 :: qp             ! quadrature points
real(8), dimension(2)                 :: wt             ! quadrature weights
real(8), dimension(:),    allocatable :: x              ! node coordinates
real(8), dimension(:, :), allocatable :: phi            ! shape functions
real(8), dimension(:, :), allocatable :: dphi           ! shape function deriv
integer, dimension(:, :), allocatable :: LM             ! location matrix

! variables to define the CG solver
integer                            :: cnt         ! number of CG iterations
real(8)                            :: convergence ! difference b/w CG iterations
real(8)                            :: reltol      ! CG relative tolerance
real(8)                            :: m           ! slope of line
real(8), dimension(:), allocatable :: z           ! CG update iterates
real(8), dimension(:), allocatable :: res         ! solution residual

! variables to define the local problem
integer                               :: n_el        ! number of local elements
integer                               :: n_nodes     ! number of (local) nodes
integer                               :: numprocs    ! number of processors
integer                               :: maxperproc  ! max elems per processor
integer                               :: rank        ! processor rank
integer                               :: ddcnt       ! DD counter
real(8)                               :: itererror   ! whole-loop iter error
real(8)                               :: ddtol       ! DD loop tolerance
integer                               :: ierr        ! error for MPI calls
integer, dimension(:, :), allocatable :: edges       ! nodes on edge of domain
integer, dimension(:),    allocatable :: recv_displs ! displacement of domain
real(8), dimension(:),    allocatable :: xel         ! coordinates in domain
integer, dimension(:),    allocatable :: numnodes    ! nodes in domain
integer, dimension(:),    allocatable :: LMcount     ! elements in row i of K
integer, dimension(:),    allocatable :: elems       ! n_el in each domain
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
integer, dimension(:),    allocatable :: LMcountcoarse ! elements in row i of K
integer, dimension(:, :), allocatable :: LMcoarse      ! LM of coarse problem

! variables to define OpenMP thread parallelism
integer :: numthreads ! number of OpenMP threads
integer :: mythread   ! current thread number
integer :: provided   ! holds provided level of thread support
integer :: omp_get_thread_num, omp_get_num_threads ! OpenMP routines

! variables to implement CSR matrix storage
type row
  real(8), allocatable, dimension(:) :: values     ! values in each row
  integer, allocatable, dimension(:) :: columns    ! nonzero col nums per row
  integer                            :: entri = 1  ! entry next to be filled    
end type row

type(row), allocatable, dimension(:) :: rows       ! rows in global K
type(row), allocatable, dimension(:) :: rowscoarse ! rows in coarse global K

! read in variable values for simulation from namelist
namelist /FEM/ k, source, n_qp, reltol, ddtol
open(20, file='setup.nml')
read(20, FEM)
close(20)


!call hello()
call cpu_time(start)

! initialize variables that are known to all MPI processes
! instead of having one process compute these and then broadcast
call commandline(n_el, length, leftBC, rightBC)   ! parse command line args
call initialize(h, x, n_el, n_nodes)              ! initialize problem vars
call quadrature(qp, wt, n_qp)                     ! initialize quadrature
call phi_val(qp)                                  ! initialize shape funcs
rel = elementalload(wt, phi, source, h) 
kel = elementalstiffness(wt, dphi, k)

n_nodes_global = n_nodes
n_el_global    = n_el

! initialize the parallel MPI environment with mpi_thread_single thread support
call mpi_init_thread(0, provided, ierr)
call mpi_comm_size(mpi_comm_world, numprocs, ierr)
call mpi_comm_rank(mpi_comm_world, rank, ierr)

! only the rank 0 process knows about the whole solution vector
if (rank == 0) then
  allocate(soln(n_nodes_global), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of soln array failed."
  soln = 0.0
end if

call allocatedecomp()           ! allocate space for decompsition data
call initializedecomp()         ! initialize domain decomposition

! perform a coarse solution to get initial guesses for the interface values
! this only needs to be performed by the rank 0 process
if (rank == 0) then
  n_el    = numprocs
  n_nodes = numprocs + 1
  
  call allocatecoarse()

  xcoarse(1) = x(edges(1, 1))
  do i = 2, n_nodes
    xcoarse(i) = x(edges(2, i - 1))
    hlocal(i - 1) = xcoarse(i) - xcoarse(i - 1)
  end do

  BCs    = (/ 1, n_nodes /) 
  BCvals = (/ leftBC, rightBC /)

  call locationmatrix(LMcoarse, LMcountcoarse, n_el)      ! form the location matrix

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

  call csr(rowscoarse, kel, LMcoarse, LMcountcoarse)

  ! initial guess is a straight line between the two endpoints
  m = (rightBC - leftBC) / length
  acoarse = m * (xcoarse - xcoarse(1)) + BCvals(1)

  call conjugategradient(rowscoarse, acoarse, rglobcoarse, zcoarse, rescoarse, BCs, reltol = reltol)
  
  ! insert first-guess BCs into BCcoarse array, then broadcast to all processes
  BCcoarse(1, 1) = acoarse(1)
  BCcoarse(2, numprocs) = acoarse(n_nodes)

  do i = 1, numprocs - 1
    BCcoarse(2, i) = acoarse(i + 1)
    BCcoarse(1, i + 1) = acoarse(i + 1)
  end do

  deallocate(xcoarse, LMcoarse, rglobcoarse, acoarse, hlocal)  
  deallocate(LMcountcoarse, rowscoarse)
end if

call mpi_bcast(BCcoarse, 2 * numprocs, mpi_real8, 0, mpi_comm_world, ierr)

! Specify the domain decomposition parameters for each domain -----------------
n_el    = elems(rank + 1)
n_nodes = numnodes(rank + 1)

call allocateDDdata()

xel       = x(edges(1, rank + 1):edges(2, rank + 1))
BCs       = (/1, n_nodes/)
BCvals(1) = BCcoarse(1, rank + 1)
BCvals(2) = BCcoarse(2, rank + 1)

call locationmatrix(LM, LMcount, n_el)          ! form LM and count entries
rglob = globalload(LM, rel, n_el, n_nodes)      ! form global load vector
call csr(rows, kel, LM, LMcount)                ! form CSR storage

! initial guess is a straight line between the two endpoints
m = (BCvals(2) - BCvals(1)) / (xel(n_nodes) - xel(1))
a = m * (xel - xel(1)) + BCvals(1)

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
    call mpi_send(a(n_nodes - 1), 1, mpi_real8, rank + 1, rank, mpi_comm_world, ierr)
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
call mpi_gatherv(a(1:(n_nodes - 1)), n_nodes - 1, mpi_real8, &
                 soln(1:(n_nodes_global - 1)), elems, recv_displs, mpi_real8, &
                 0, mpi_comm_world, ierr)

! write results to output file ------------------------------------------------
if (rank == 0) then
  soln(n_nodes_global) = rightBC 
  ! write to an output file. If this file exists, it will be re-written.
  open(1, file='output.txt', iostat=AllocateStatus, status="replace")
  if (AllocateStatus /= 0) STOP "output.txt file opening failed."

  write(1, *) numprocs, n_el_global, ddcnt, soln
end if

! final timing results --------------------------------------------------------
if (rank == 0) then
  call cpu_time(finish)
  print *, 'P: ', numprocs, 'n_el: ', n_el_global, 'runtime: ', finish - start
end if

! deallocate memory -----------------------------------------------------------
deallocate(xel, LM, rglob, a, z, res)
deallocate(numnodes, elems, edges, recv_displs, BCcoarse, LMcount)
deallocate(rows)
if (rank == 0) deallocate(soln)

call mpi_finalize(ierr)

deallocate(x, phi, dphi)

! ------------------------------------------------------------------------------


CONTAINS ! define all internal procedures

subroutine allocatecoarse()
  implicit none
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
  allocate(LMcountcoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of LMcountcoarse array failed."
  allocate(rowscoarse(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rowscoarse array failed."
end subroutine allocatecoarse


subroutine allocateDDdata()
  implicit none
  allocate(xel(n_nodes), stat = AllocateStatus)
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
  allocate(LMcount(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of LMcount array failed."
  allocate(rows(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rows array failed."
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


subroutine csr(rows, kel, LM, LMcount)
  implicit none
  type(row), intent(inout) :: rows(:)
  real(8), intent(in)      :: kel(:, :)
  integer, intent(in)      :: LM(:, :)
  integer, intent(in)      :: LMcount(:)

  integer :: i, j, q, n_el, n_nodes

  n_el    = size(LM(1, :))
  n_nodes = size(rows)

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
end subroutine csr


subroutine initializedecomp()
  implicit none

  ! distribute the elements among the processors 
  maxperproc = (n_el + numprocs - 1) / numprocs
  elems = maxperproc
 
  i = 1
  do j = maxperproc * numprocs - n_el, 1, -1
    elems(i) = elems(i) - 1
    i = i + 1
    if (i == numprocs + 1) i = 1
  end do
 
  ! assign the numbers of nodes in each domain 
  !!$omp parallel do default(private) shared(numprocs, elems, numnodes) private(j)
  do j = 1, numprocs
    numnodes(j) = elems(j) * 2 - (elems(j) - 1)
  end do
  !!$omp end parallel do
  
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


function elementalstiffness(wt, dphi, k) result(kelem)
  implicit none
  real(8), intent(in) :: k
  real(8), intent(in) :: wt(:)
  real(8), intent(in) :: dphi(:, :)
  real(8)             :: kelem(2, 2)
  integer             :: q, i, n

  n     = size(wt)
  kelem = 0.0

  do q = 1, n
    do i = 1, 2
      kelem(i, 1) = kelem(i, 1) + wt(q) * dphi(i, q) * k * dphi(1, q) * 2.0
      kelem(i, 2) = kelem(i, 2) + wt(q) * dphi(i, q) * k * dphi(2, q) * 2.0
    end do
  end do
end function elementalstiffness


function elementalload(wt, phi, source, h) result(relem)
  implicit none
  real(8), intent(in) :: source, h
  real(8), intent(in) :: wt(:)
  real(8), intent(in) :: phi(:, :)
  real(8)             :: relem(2)
  integer             :: q, i, n

  n     = size(wt)
  relem = 0.0

  do q = 1, n
    do i = 1, 2
      relem(i) = relem(i) + wt(q) * source * phi(i, q) * h * h / 2.0
    end do
  end do
end function elementalload


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


real function dotprod(vec1, vec2)
  implicit none
  real(8)  :: vec1(:), vec2(:)
  integer  :: i, n
  
  n = size(vec1)
  dotprod = 0.0  

  do i = 1, n
    dotprod = dotprod + vec1(i) * vec2(i)
  end do
end function dotprod


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


subroutine locationmatrix(LM, LMcount, n_el)
  implicit none
  integer, intent(inout) :: LM(:, :)
  integer, intent(inout) :: LMcount(:)
  integer, intent(in)    :: n_el
  integer                :: j
  
  ! Determine the number of elements that contain each node (given as a 
  ! vector of length n_nodes). Then, the formula to obtain the number of 
  ! nonzero entries per row is (LMcount - 1 + n_en), so compute the number
  ! of nonzero entries per row in the global stiffness matrix.
 
  LMcount = 1
 
  do j = 1, n_el
    LM(:, j) = (/ j, j + 1 /)
    LMcount(LM(1, j)) = LMcount(LM(1, j)) + 1
    LMcount(LM(2, j)) = LMcount(LM(2, j)) + 1
  end do
end subroutine locationmatrix

subroutine phi_val(qp)
  implicit none
  real(8), intent(in)  :: qp(:)

  allocate(phi(2, size(qp)), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of phi array failed."
  allocate(dphi(2, size(qp)), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of dphi array failed."

  phi(1, :)  = (1.0 - qp(:)) / 2.0
  phi(2, :)  = (1.0 + qp(:)) / 2.0
  dphi(1, :) = -1.0 / 2.0
  dphi(2, :) =  1.0 / 2.0
end subroutine phi_val


subroutine quadrature(qp, wt, n_qp)
  implicit none
  real(8), intent(inout) :: qp(:), wt(:)
  integer, intent(in)    :: n_qp
 
  if (n_qp == 2) then 
    qp   = (/ -1.0/sqrt(3.0), 1.0/sqrt(3.0) /)
    wt   = (/ 1.0, 1.0 /)
  else
    print *, 'error in quadrature rule selection.'
  end if
end subroutine quadrature


subroutine commandline(n_el, length, leftBC, rightBC)
  implicit none
  integer, intent(out) :: n_el
  real(8), intent(out) :: length
  real(8), intent(out) :: leftBC
  real(8), intent(out) :: rightBC
  integer              :: nargs
  integer              :: i
  character(len = 12)  :: args

  nargs = command_argument_count()

  do i = 1, nargs
    call get_command_argument(i, args)
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
  integer              :: i

  h = length / real(n_el)
  n_nodes = n_el + 1 

  allocate(x(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of x array failed."

  !!$omp parallel do default(private) shared(x, h) private(i)
  do i = 1, size(x)
    x(i) = real(i - 1) * h
  end do
  !!$omp end parallel do
end subroutine initialize


END PROGRAM main
