! non-overlapping interface
! each processor computes one interface problem
! builds from bmain by adding CSR matrix storage

! program to solve the heat equation (Dirichlet boundary conditions only)
! using domain decomposition. This version has the domain to the left of
! each interface send a single interface value to the processor to the
! right of it. This process then computes the interface problem, and 
! sends a single udpated value to the processor to the left.

PROGRAM main

use read_data
use mesh
use quadrature
use element_matrices
use csr_storage
use solvers

implicit none

! variables for overall execution
integer  :: AllocateStatus     ! variable to hold memory allocation success
integer  :: i, j, q            ! loop iteration variables
real(8)  :: start              ! start run time
real(8)  :: finish             ! end run time
real(8)  :: startCG            ! start CG time
real(8)  :: endCG              ! end CG time
real(8)  :: startCSR           ! start CSR time
real(8)  :: endCSR             ! end CSR time
real(8) :: mm

! variables to define the global problem
real(8), dimension(2)                 :: prev        ! prev interface values
real(8), dimension(2)                 :: BClocals    ! BCs for each interface
real(8), dimension(:),    allocatable :: hlocal      ! coarse element lengths
real(8), dimension(:, :), allocatable :: BCcoarse    ! coarse solution BCs

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

allocate(global%a(global%n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of global%a array failed."

allocate(BCcoarse(2, pretend_procs), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of BCcoarse array failed."

call initialize_domain_decomposition(pretend_procs)

! perform a coarse solution to get initial guesses for the interface values
! this only needs to be performed by the rank 0 process
if (.true.) then
  call initialize_coarse_mesh()
  
  allocate(rowscoarse(coarse%n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of rows array failed."
  
  allocate(hlocal(coarse%n_el), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of hlocal array failed."

  do i = 2, coarse%n_nodes
    hlocal(i - 1) = coarse%x(i) - coarse%x(i - 1)
  end do

  LMcoarse = locationmatrix(coarse%n_el, coarse%n_nodes)
  coarse%rglob = globalvector(LMcoarse, rel, coarse%n_nodes, hlocal)
  
  coarse%rglob(coarse%BCs(1)) = coarse%BCvals(1)
  coarse%rglob(coarse%BCs(2)) = coarse%BCvals(2)   

  call form_csr(LMcoarse, coarse%n_nodes, rowscoarse)

  ! initial guess for CG is a straight line between the two endpoints
  coarse%a = straightline(coarse)
  coarse%a = conjugategradient(rowscoarse, coarse%a, coarse%rglob, coarse%BCs, reltol)

  ! insert first-guess BCs into initial guess for solution
  do i = 1, pretend_procs
    BCcoarse(:, i) = (/coarse%a(i), coarse%a(i + 1)/)
    mm = (BCcoarse(2, i) - BCcoarse(1, i)) / (coarse%x(i + 1) - coarse%x(i)) 
    global%a(domains%edges(1, i):(domains%edges(2, i) - 1)) = &
        mm * (global%x(domains%edges(1, i):(domains%edges(2, i) - 1)) &
              - global%x(domains%edges(1, i))) + BCcoarse(1, i)
  end do
  
  global%a(global%n_nodes) = rightBC 
end if


! Specify the domain decomposition parameters for each domain -----------------
allocate(rows(global%n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of rows array failed."

LMfine = locationmatrix(global%n_el, global%n_nodes)
global%rglob  = globalvector(LMfine, rel, global%n_nodes)
global%rglob(global%BCs(1)) = global%BCvals(1)
global%rglob(global%BCs(2)) = global%BCvals(2)

call form_csr(LMfine, global%n_nodes, rows)

global%a = conjugategradient(rows, global%a, global%rglob, &
                                     global%BCs, reltol=reltol)
  
! write results to output file ------------------------------------------------
if (.false.) then
  global%a(global%n_nodes) = rightBC 
  ! write to an output file. If this file exists, it will be re-written.
  open(1, file='output.txt', iostat=AllocateStatus, status="replace")
  if (AllocateStatus /= 0) STOP "output.txt file opening failed."

  write(1, *) pretend_procs, global%n_el, 0, global%a
end if

! final timing results --------------------------------------------------------
call cpu_time(finish)
print *, 'SERAIL P: ', pretend_procs, 'n_el: ', global%n_el, 'runtime: ', finish - start

! deallocate memory -----------------------------------------------------------
deallocate(LMfine%matrix, LMfine%cnt, BCcoarse, rows)
deallocate(LMcoarse%matrix, LMcoarse%cnt, hlocal, rowscoarse)  

deallocate(global%a, global%rglob)
call dealloc_coarse_mesh()
  
call dealloc_domains()
call dealloc_global_mesh()
call dealloc_shapefunctions()
call dealloc_quadset()

END PROGRAM
