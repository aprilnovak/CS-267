module read_data

implicit none

integer, save :: n_el_global ! number of global elements
integer, save :: n_qp        ! number of quadrature points
real(8), save :: length      ! domain length
real(8), save :: leftBC      ! left Dirichlet boundary condition
real(8), save :: rightBC     ! right Dirichlet boundary condition
real(8), save :: k           ! thermal conductivity
real(8), save :: source      ! heat source
real(8), save :: reltol      ! CG relative convergence tolerance
real(8), save :: ddtol       ! DD convergence tolerance

contains

subroutine read_commandline()
  integer              :: nargs
  integer              :: i
  character(len = 12)  :: args

  nargs = command_argument_count()

  do i = 1, nargs
    call get_command_argument(i, args)
    select case (i)
      case(1)
        read(args, *) length
      case(2)
        read(args, *) n_el_global
      case(3)
        read(args, *) leftBC
      case(4)
        read(args, *) rightBC
      case default
        write(*,*) "Too many command line parameters."
    end select
  enddo
end subroutine read_commandline


subroutine read_namelist()
  namelist /FEM/ k, source, n_qp, reltol, ddtol
  open(20, file='setup.nml')
  read(20, FEM)
  close(20)
end subroutine read_namelist


end module read_data
