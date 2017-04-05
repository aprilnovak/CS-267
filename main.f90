PROGRAM main

! turn off implicit typing
implicit none

! define double precision with 15 digits of accuracy
! and an exponent range of +- 307
integer, parameter :: rk = selected_real_kind(15, 307)

! program variables related to execution
integer             :: nargs   ! number of command line parameters
integer             :: i       ! looping variable for reading command line params
character(len = 12) :: args    ! command line parameters

! program variables related to FEM solve
real(rk) :: length             ! length of the domain (1-D)
integer  :: n_el               ! number of elements
integer  :: order              ! polynomial order
real (rk), dimension(13) :: x  ! coordinates of the nodes

! obtain number of command line arguments
nargs = command_argument_count()

! parse command line arguments
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

write(*,*) length


END PROGRAM main
