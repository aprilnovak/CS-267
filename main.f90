PROGRAM main

! turn off implicit typing
implicit none

! define double precision with 15 digits of accuracy
! and an exponent range of +- 307
integer, parameter :: rk = selected_real_kind(15, 307)

! program variables related to execution
integer             :: AllocateStatus ! variable to hold memory allocation success

! program variables related to FEM solve
real(rk) :: length             ! length of the domain (1-D)
integer  :: n_el               ! number of elements
integer  :: n_en               ! number of nodes per element
integer  :: n_nodes            ! total number of nodes
integer  :: order              ! polynomial order
real (rk), dimension(:), allocatable :: x  ! coordinates of the nodes
real(rk) :: h                  ! length of one element
real (rk), dimension(:, :), allocatable :: qp ! quadrature points and weights

! parse command line arguments
call commandline(length, n_el, order)

! initialize problem variables
call initialize(h, x, n_en, n_el, order, n_nodes)





allocate(qp(ceiling((real(order) + 1.0) / 2.0), 2), stat = AllocateStatus)
if (AllocateStatus /= 0) STOP "Allocation of qp array failed."


! initialize the quadrature rule
call quadrature(order)







! deallocate memory 
deallocate(x)

CONTAINS

subroutine quadrature(order)
  integer, intent(in) :: order
end subroutine quadrature


subroutine commandline(length, n_el, order)
  ! define subroutine parameters
  real(rk), intent(out) :: length
  integer, intent(out)  :: n_el
  integer, intent(out)  :: order
 
  ! define local variables 
  integer :: nargs            ! number of command line arguments
  integer :: i                ! looping variable
  character(len = 12) :: args ! command line argument

  nargs = command_argument_count()

  do i = 1, nargs
    call get_command_argument(i, args)
  
    ! use internal reads to convert from character to numeric types
    ! (read from the character buffer into the numeric variable)
    select case (i)
      case(1) ! length of domain
        read(args, *) length  
      case(2) ! number of elements
        read(args, *) n_el
      case(3) ! polynomial order
        read(args, *) order
      case default
        write(*,*) "Too many command line parameters."
    end select  
  enddo
end subroutine commandline


subroutine initialize(h, x, n_en, n_el, order, n_nodes)
  real(rk), intent(out) :: h ! element length
  real (rk), dimension(:), allocatable, intent(out) :: x ! node coordinates
  integer, intent(out) :: n_en    ! number of nodes per element
  integer, intent(in)  :: n_el    ! number of elements
  integer, intent(in)  :: order   ! polynomial order
  integer, intent(out) :: n_nodes ! number of total nodes

  h = length / real(n_el)
  n_en = order + 1
  n_nodes = (order + 1) * n_el - (n_el - 1)
   
  ! allocate memory for the vector of node coordinates
  allocate(x(n_nodes), stat = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Allocation of x array failed."

end subroutine initialize



END PROGRAM main
