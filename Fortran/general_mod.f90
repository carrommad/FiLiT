!------------------------------------------------------------------------------
! MODULE - general_mod.f90
!------------------------------------------------------------------------------
!> 
!> DESCRIPTION:
!> Module containing general shared constants and simple functions.
!> 
!> REVISION HISTORY:
!> 22 02 2021 - Carlos Romero
!> 
!> TO DO:
!> - Add subroutines
!> 
!------------------------------------------------------------------------------
module general_mod
    
    implicit none
    
    public

    integer, parameter :: dp = selected_real_kind(15, 307)
  
contains
    !------------------------------------------------------------------------------
    ! SUBROUTINE - linspace
    !------------------------------------------------------------------------------
    !> 
    !> DESCRIPTION:
    !> Create a vector of N points from values start to end
    !> 
    !> 
    !> INPUT/OUTPUT:
    !> - [IN]  start - initial point of array
    !> - [IN]  end   - final point of array
    !> - [IN]  N     - number of equispaced points of array
    !> - [OUT] array - equispaced array of length N from start to end
    !> 
    !> EXAMPLE:
    !> - Example:
    !>      start = 0, end = 1; N = 5
    !>      call linspace(start, end, N, array)
    !>      $ array = [0. 0.25 0.5 0.75 1.]
    !> 
    !> REVISION HISTORY:
    !> 14 02 2021 - Carlos Romero
    !> 
    !> TO DO:
    !> - 
    !> 
    !------------------------------------------------------------------------------
    subroutine linspace(start, end, N, array)

        implicit none 
        integer, intent(in) :: N
        real(dp), intent(in) :: start, end
        real(dp), intent(out), allocatable :: array(:)
        integer :: i

        allocate(array(N))

        if (N == 0) return

        if (N == 1) then
            array(1) = start
            return
        end if

        do i=1, N
            array(i) = start + (end - start) * (i - 1) / (N - 1)
        end do

    end subroutine linspace

    !------------------------------------------------------------------------------
    ! SUBROUTINE - meshgrid_param
    !------------------------------------------------------------------------------
    !> 
    !> DESCRIPTION:
    !> Create 2-dimensional uniform rectangular grid provided grid size parameters.
    !> 
    !> INPUT/OUTPUT:
    !> - [IN]  N_x - grid size in the x-direction
    !> - [IN]  N_y - grid size in the y-direction
    !> - [IN]  L_x - size of the rectangular grid in the x-direction
    !> - [IN]  L_y - size of the rectangular grid in the y-direction
    !> - [OUT] G_x - grid of size (N_x, N_y) containing x-direction equispacing
    !> - [OUT] G_y - grid of size (N_x, N_y) containing y-direction equispacing 
    !> 
    !> EXAMPLE:
    !> - Example with small grid sizes (remember to declare real(dp) variables):
    !>      N_x = 6                  
    !>      N_y = 5
    !>      L_x = 1
    !>      L_y = 1
    !>      G = grid_2D(N_x, N_y, L_x, L_y)
    !>      $ G_x = [0.   0.2  0.4  0.6  0.8  1.;
    !>               0.   0.2  0.4  0.6  0.8  1.;
    !>               0.   0.2  0.4  0.6  0.8  1.;
    !>               0.   0.2  0.4  0.6  0.8  1.;
    !>               0.   0.2  0.4  0.6  0.8  1.]
    !>      $ G_y = [0.   0.   0.   0.   0.   0.  ;
    !>               0.25 0.25 0.25 0.25 0.25 0.25;
    !>               0.5  0.5  0.5  0.5  0.5  0.5 ;
    !>               0.75 0.75 0.75 0.75 0.75 0.75;
    !>               1.   1.   1.   1.   1.   1.  ]
    !> 
    !> REVISION HISTORY:
    !> 22 02 2021 - Carlos Romero
    !> 
    !> TO DO:
    !> -
    !> 
    !------------------------------------------------------------------------------
    subroutine meshgrid_param(N_x, N_y, L_x, L_y, G_x, G_y)

        implicit none
        integer, intent(in) :: N_x, N_y
        real(dp), intent(in) :: L_x, L_y
        real(dp), intent(out), allocatable :: G_x(:,:), G_y(:,:)
        integer :: i, j
        real(dp), allocatable :: array_x(:), array_y(:)

        allocate(G_x(N_x,N_y))
        allocate(G_y(N_x,N_y))

        call linspace(0._dp, L_x, N_x, array_x)
        do j=1, N_y
            G_x(j,:) = array_x
        end do
        
        call linspace(0._dp, L_y, N_y, array_y)
        do i=1, N_x
            G_y(:,i) = array_y
        end do
    
    end subroutine meshgrid_param
    
    !------------------------------------------------------------------------------
    ! SUBROUTINE - meshgrid_array
    !------------------------------------------------------------------------------
    !> 
    !> DESCRIPTION:
    !> Create 2-dimensional uniform rectangular grid provided two arrays.
    !> 
    !> INPUT/OUTPUT:
    !> - [IN]  array_x - array in the x-direction to be meshed
    !> - [IN]  array_y - array in the y-direction to be meshed
    !> - [OUT] G_x     - grid of size (size(array_x), size(array_y)) containing 
    !>                 x-direction equispacing
    !> - [OUT] G_y     - grid of size (size(array_x), size(array_y)) containing 
    !>                 y-direction equispacing
    !> 
    !> EXAMPLE:
    !> - Example with small grid sizes (remember to declare real(dp) variables):
    !>      array_x = [0.   0.2  0.4  0.6  0.8  1.]
    !>  	array_y = [0.   0.25 0.5  0.75 1.]
    !>      G = grid_2D(array_x, array_y, L_x, L_y)
    !>      $ G_x = [0.   0.2  0.4  0.6  0.8  1.;
    !>               0.   0.2  0.4  0.6  0.8  1.;
    !>               0.   0.2  0.4  0.6  0.8  1.;
    !>               0.   0.2  0.4  0.6  0.8  1.;
    !>               0.   0.2  0.4  0.6  0.8  1.]
    !>      $ G_y = [0.   0.   0.   0.   0.   0.  ;
    !>               0.25 0.25 0.25 0.25 0.25 0.25;
    !>               0.5  0.5  0.5  0.5  0.5  0.5 ;
    !>               0.75 0.75 0.75 0.75 0.75 0.75;
    !>               1.   1.   1.   1.   1.   1.  ]
    !> 
    !> REVISION HISTORY:
    !> 22 02 2021 - Carlos Romero
    !> 
    !> TO DO:
    !> -
    !> 
    !------------------------------------------------------------------------------
    subroutine meshgrid_array(array_x, array_y, G_x, G_y)

        implicit none

		real(dp), allocatable, intent(in)  :: array_x(:), array_y(:)
        real(dp), allocatable, intent(out) :: G_x(:,:), G_y(:,:)
        integer :: N_x, N_y, i, j
        
        N_x = size(array_x)
        N_y = size(array_y)

        allocate(G_x(N_x,N_y))
        allocate(G_y(N_x,N_y))

        do j=1, N_y
            G_x(j,:) = array_x
        end do
        
        do i=1, N_x
            G_y(:,i) = array_y
        end do
    
    end subroutine meshgrid_array
    
end module general_mod
