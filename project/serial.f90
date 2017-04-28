! non-overlapping interface
! each processor computes one interface problem
! builds from bmain by adding CSR matrix storage
! SERIAL

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
real(8) :: m

! variables to define the global problem
real(8), dimension(:),    allocatable :: soln        ! global soln vector
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

allocate(soln(global%n_nodes), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of soln array failed."
soln = 0.0

allocate(BCcoarse(2, pretend_procs), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of BCcoarse array failed."

call initialize_domain_decomposition(pretend_procs)
call initialize_coarse_mesh()

allocate(hlocal(coarse%n_el), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of hlocal array failed."

do i = 2, coarse%n_nodes
  hlocal(i - 1) = coarse%x(i) - coarse%x(i - 1)
end do

LMcoarse = locationmatrix(coarse%n_el, coarse%n_nodes)
coarse%rglob = globalvector(LMcoarse, rel, coarse%n_nodes, hlocal)

coarse%rglob(coarse%BCs(1)) = coarse%BCvals(1)
coarse%rglob(coarse%BCs(2)) = coarse%BCvals(2)   

rowscoarse = form_csr(LMcoarse, coarse%n_nodes)

! initial guess for CG is a straight line between the two endpoints
coarse%a = straightline(coarse)
coarse%a = conjugategradient(rowscoarse, coarse%a, coarse%rglob, coarse%BCs, reltol)

! insert first-guess BCs into BCcoarse array, then broadcast to all processes
do i = 1, pretend_procs
  BCcoarse(:, i) = (/coarse%a(i), coarse%a(i + 1)/)
end do

print *, BCcoarse

! give initial guess for the solution
do i = 1, pretend_procs
  m = (BCcoarse(2, i) - BCcoarse(1, i)) / (coarse%x(i + 1) - coarse%x(i)) 
!  print *, 'slope: ', m
!  print *, 'edges: ', domains%edges(1, i), 'right: ', (domains%edges(2, i))
  soln(domains%edges(1, i):(domains%edges(2, i) - 1)) = m * (global%x(domains%edges(1, i):(domains%edges(2, i) - 1)) - global%x(domains%edges(1, i))) + BCcoarse(1, i)
!  print *, 'solution: ', soln
end do


!allocate(rows(dd(rank + 1)%n_nodes), stat = AllocateStatus)
!if (AllocateStatus /= 0) STOP "Allocation of rows array failed."
!
!dd(rank + 1)%BCvals(1) = BCcoarse(1, rank + 1)
!dd(rank + 1)%BCvals(2) = BCcoarse(2, rank + 1)
!
!LMfine = locationmatrix(dd(rank + 1)%n_el, dd(rank + 1)%n_nodes)
!dd(rank + 1)%rglob  = globalvector(LMfine, rel, dd(rank + 1)%n_nodes)
!rows   = form_csr(LMfine, dd(rank + 1)%n_nodes)
!
!! initial guess is a straight line between the two endpoints
!dd(rank + 1)%a = straightline(dd(rank + 1))

!do while (itererror > ddtol)
!  ! save the previous values of the interface BCs
!  prev = dd(rank + 1)%BCvals
!
!  ! each processor solves for its domain --------------------------------------
!  dd(rank + 1)%rglob(dd(rank + 1)%BCs(1)) = dd(rank + 1)%BCvals(1)
!  dd(rank + 1)%rglob(dd(rank + 1)%BCs(2)) = dd(rank + 1)%BCvals(2)   
!
!  dd(rank + 1)%a = conjugategradient(rows, dd(rank + 1)%a, dd(rank + 1)%rglob, &
!                                     dd(rank + 1)%BCs, reltol=reltol)
!  
!  ! each processor sends a boundary value to the processor to the right -------
!  if (rank /= pretend_procs - 1) then
!    call mpi_send(dd(rank + 1)%a(dd(rank + 1)%n_nodes - 1), 1, mpi_real8, &
!                  rank + 1, rank, mpi_comm_world, ierr)
!  end if
!
!  ! processor to the right receives the message -------------------------------
!  if (rank /= 0) then
!    call mpi_recv(BClocals(1), 1, mpi_real8, rank - 1, rank - 1, mpi_comm_world, stat, ierr)
!    BClocals(2) = dd(rank + 1)%a(2)
!  end if
!
!  call mpi_barrier(mpi_comm_world, ierr)
!
!  ! each processor solves its interface problem -------------------------------
!  if (rank /= 0) then
!    dd(rank + 1)%BCvals(1) = (rel(2) + rel(1) - kel(2, 1) * BClocals(2) &
!                      - kel(1, 2) * BClocals(1)) / (kel(2, 2) + kel(1, 1))
!    ! send new interface result to rank - 1 process (to BCvals(2) of rank - 1)
!    call mpi_send(dd(rank + 1)%BCvals(1), 1, mpi_real8, rank - 1, rank, mpi_comm_world, ierr)
!  end if
!
!  ! rank - 1 process receives from the process to the right -------------------
!  if (rank /= pretend_procs - 1) then
!    call mpi_recv(dd(rank + 1)%BCvals(2), 1, mpi_real8, rank + 1, rank + 1, &
!                  mpi_comm_world, stat, ierr)
!  end if
!
!  ! compute iteration error to determine whether to continue looping ---------- 
!  call mpi_barrier(mpi_comm_world, ierr)
!
!  call mpi_allreduce(abs(dd(rank + 1)%BCvals(2) - prev(2) + dd(rank + 1)%BCvals(1) - prev(1)),&
!                     itererror, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)  
! 
!  call mpi_barrier(mpi_comm_world, ierr)
!  ddcnt = ddcnt + 1
!end do ! ends outermost domain decomposition loop
!
!! each processor broadcasts its final solution to the rank 0 process ----------
!call mpi_gatherv(dd(rank + 1)%a(1:(dd(rank + 1)%n_nodes - 1)), dd(rank + 1)%n_nodes - 1, mpi_real8, &
!                 soln(1:(global%n_nodes - 1)), domains%elems, domains%recv_displs, mpi_real8, &
!                 0, mpi_comm_world, ierr)
!
! write results to output file ------------------------------------------------
!if (.false.) then
soln(global%n_nodes) = rightBC 
! write to an output file. If this file exists, it will be re-written.
open(1, file='output.txt', iostat=AllocateStatus, status="replace")
if (AllocateStatus /= 0) STOP "output.txt file opening failed."

write(1, *) pretend_procs, global%n_el, 0, soln
!end if

! final timing results --------------------------------------------------------
  call cpu_time(finish)
  print *, 'serial procs: ', pretend_procs, 'n_el: ', global%n_el, 'runtime: ', finish - start

! deallocate memory -----------------------------------------------------------
!deallocate(LMfine%matrix, LMfine%cnt, BCcoarse, rows)
deallocate(LMfine%cnt, BCcoarse)
deallocate(soln)
deallocate(LMcoarse%matrix, LMcoarse%cnt, hlocal, rowscoarse)  
call dealloc_coarse_mesh()
call dealloc_domains()
call dealloc_global_mesh()
call dealloc_shapefunctions()
call dealloc_quadset()

END PROGRAM
