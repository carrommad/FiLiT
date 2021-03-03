!------------------------------------------------------------------------------
! MODULE - field_line_mod.f90
!------------------------------------------------------------------------------
!> 
!> DESCRIPTION:
!> Contains auxiliary functions to trace magnetic field lines
!>
!> REVISION HISTORY:
!> 22 02 2021 - Carlos Romero
!> 
!------------------------------------------------------------------------------
module field_line_mod

	use netcdf
	use general_mod

	integer :: nr,   &
			   nphi, &
			   nz
			 
	real(8), parameter :: pi  = 3.14159265358979323846d0, & 
						  pi2 = 2.0d0 * pi

contains 
	!------------------------------------------------------------------------------
	! SUBROUTINE - check
	!------------------------------------------------------------------------------
	!> 
	!> DESCRIPTION:
	!> Checks if a netCDF function outputs an error and prints on screen
	!> consequently.
	!> 
	!> INPUT/OUTPUT:
	!> - [IN] status - output status of nf90 functions 
	!> 
	!> EXAMPLE:
	!> - Example:
	!>      call check(nc
	!>
	!> REVISION HISTORY:
	!> 22 02 2021 - Carlos Romero
	!> 
	!> TO DO:
	!> -
	!> 
	!------------------------------------------------------------------------------
	subroutine check(status)
		
		integer, intent(in) :: status
		
		if (status /= nf90_noerr) then
			print *, trim(nf90_strerror(status))
			stop 'Stopped due to netCDF error!'
		end if

	end subroutine check

	!------------------------------------------------------------------------------
	! SUBROUTINE - read_data_hint
	!------------------------------------------------------------------------------
	!> 
	!> DESCRIPTION:
	!> Reads netcdf output data from HINT snapfiles.
	!> Adding or reducing returned variables is possible.
	!> 
	!> INPUT/OUTPUT:
	!> - [IN]  filepath - path of netcdf data structure of an output (snapfile)
	!>                  from HINT
	!> - [OUT] rminb    - minimum value of R for the grid
	!> - [OUT] rmaxb    - maximum value of R for the grid
	!> - [OUT] zminb    - minimum value of Z for the grid
	!> - [OUT] zmaxb    - maximum value of Z for the grid
	!> - [OUT] B_R      - 3-tensor of radial magnetic field in the reactor
	!> - [OUT] B_phi    - 3-tensor of toroidal magnetic field in the reactor
	!> - [OUT] B_Z      - 3-tensor of vertical magnetic field in the reactor
	!> 
	!> EXAMPLE:
	!> - Example:
	!>      filepath = 'folder/filename.nc'
	!>      call read_data_hint(filepath, rminb, rmaxb, zminb, zmaxb, B_R, B_phi, B_Z)
	!>      $ rminb = 2.52
	!>
	!> REVISION HISTORY:
	!> 22 02 2021 - Carlos Romero
	!> 
	!> TO DO:
	!> - Tensors v_R, v_phi, v_Z and P can also be extracted and declared as 
	!> outputs (for field line tracing is not necessary). 
	!> - Test
	!> -
	!>
	!------------------------------------------------------------------------------
	subroutine read_data_hint(filepath, rminb, rmaxb, zminb, zmaxb, B_R, B_phi, B_Z)

		character(len=50), intent(in) :: filepath

		! - Boundaries
		real(dp), intent(out) :: rminb, &
				   			 	rmaxb, &
				  		 		zminb, &
				   				zmaxb
				
		! - Magnetic (velocity and pressure) field tensors
		real(dp), allocatable, intent(out) :: B_R(:,:,:),   &
								 			  B_phi(:,:,:), &
								 			  B_Z(:,:,:)
	
		! - Dummy netCDF indices	   
		integer :: ncid,       &
				   r_dimid,    & 
				   p_dimid,    &
				   z_dimid,    &
				   rmin_varid, &
				   rmax_varid, & 
				   zmin_varid, &
				   zmax_varid, & 
				   br_varid,   &
				   bp_varid,   &
				   bz_varid

		call check( nf90_open(filepath, nf90_nowrite, ncid) )

		call check( nf90_inq_dimid(ncid, 'R',   r_dimid) )
		call check( nf90_inq_dimid(ncid, 'phi', p_dimid) )
		call check( nf90_inq_dimid(ncid, 'Z',   z_dimid) )

		call check( nf90_inquire_dimension(ncid, r_dimid, len=nr) )
		call check( nf90_inquire_dimension(ncid, p_dimid, len=np) )
		call check( nf90_inquire_dimension(ncid, z_dimid, len=nz) )

		call check( nf90_inq_varid(ncid, "rminb", rmin_varid) )
		call check( nf90_inq_varid(ncid, "rmaxb", rmax_varid) )
		call check( nf90_inq_varid(ncid, "zminb", zmin_varid) )
		call check( nf90_inq_varid(ncid, "zmaxb", zmax_varid) )


		call check( nf90_get_var(ncid, rmin_varid, rminb) )
		call check( nf90_get_var(ncid, rmax_varid, rmaxb) )
		call check( nf90_get_var(ncid, zmin_varid, zminb) )
		call check( nf90_get_var(ncid, zmax_varid, zmaxb) )

		call check( nf90_inq_varid(ncid, 'B_R',   br_varid) )
		call check( nf90_inq_varid(ncid, 'B_phi', bp_varid) )
		call check( nf90_inq_varid(ncid, 'B_Z',   bz_varid) )

		allocate(B_R(nr,np,nz))
		allocate(B_phi(nr,np,nz))
		allocate(B_Z(nr,np,nz))

		call check( nf90_get_var(ncid, br_varid, B_R) )
		call check( nf90_get_var(ncid, br_varid, B_phi) )
		call check( nf90_get_var(ncid, br_varid, B_Z) )
		
	end subroutine read_data_hint

	!------------------------------------------------------------------------------
	! SUBROUTINE - create_targets
	!------------------------------------------------------------------------------
	!> 
	!> DESCRIPTION:
	!> Creates starting target points provided boundaries a number of points in 
	!> one dimension. If 'matrix' mode selected, total targets equals ntargets^2.
	!> 
	!> INPUT/OUTPUT:
	!> - [IN]  mode       - distribution of target particles: 'line' or 'matrix'
	!> - [IN]  ntargets1D - number of target particles of one dimension
	!> - [IN]  rmin       - minimum radial coordinate
	!> - [IN]  rmax       - maximum radial coordinate
	!> - [IN]  zmin       - minimum vertical coordinate
	!> - [IN]  zmax       - maximum vertical coordinate
	!> - [OUT] targets_R  - R-coordinates of target particles
	!> - [OUT] targets_Z  - Z-coordinates of target particles
	!> 
	!> EXAMPLE:
	!> - Example:
	!>
	!> REVISION HISTORY:
	!> 22 02 2021 - Carlos Romero
	!> 
	!> TO DO:
	!> - Create example
	!> - Test
	!> -
	!> 
	!------------------------------------------------------------------------------
	subroutine create_targets(mode, ntargets1D, rmin, rmax, zmin, zmax, targets_R, targets_Z)

		character(len=6), intent(in) :: mode
		
		integer, intent(in) :: ntargets1D
		
		real(dp), intent(in) :: rmin, &
								rmax, &
								zmin, &
								zmax
		
		real(dp), allocatable, intent(out) :: targets_R(:), &
											  targets_Z(:)
		
		real(dp), allocatable :: array_R(:),    &
								 array_Z(:),    &
								 matrix_R(:,:), &
								 matrix_Z(:,:)
		
		logical, allocatable :: mask(:,:)
		
		integer :: m_shape(2), &
				   ntargets
		
		if (mode == 'line') then
			allocate(array_R(ntargets1D))
			allocate(array_Z(1))
		    call linspace(rmin, rmax, ntargets1D, array_R)
		    array_Z = (zmax - zmin) /2
		    
		    print *
		    print *, 'Default line choice: Z = (zmax - zmin) / 2'
		    print *

		else if (mode == 'matrix') then
			allocate(array_R(ntargets1D))
			allocate(array_Z(ntargets1D))
		    call linspace(rmin, rmax, ntargets1D, array_R)
		    call linspace(zmin, zmax, ntargets1D, array_Z)
		    
		end if
		
		
		call meshgrid_array(array_R, array_Z, matrix_R, matrix_Z)

		ntargets = size(array_R) * size(array_Z)
		allocate(targets_R(ntargets))
		allocate(targets_Z(ntargets))
		
		m_shape = shape(matrix_R)
		allocate(mask(m_shape(1),m_shape(2)))
		mask = matrix_R == matrix_R
		
		targets_R = pack(matrix_R, mask)
		targets_Z = pack(matrix_Z, mask)
		
		deallocate(array_R,array_Z,matrix_R,matrix_Z)
		
	end subroutine create_targets

	!------------------------------------------------------------------------------
	! SUBROUTINE - grid_neighbours
	!------------------------------------------------------------------------------
	!> 
	!> DESCRIPTION:
	!> Given an array of coordinates (R,Z) and a grid of R-Z space, computes the
	!> indices of grid for which the coordinates are closest and the coordinate 
	!> value is just above the grid's value  
	!> 
	!> INPUT/OUTPUT:
	!> - [IN]  targets_R - R-coordinates of target points
	!> - [IN]  targets_Z - Z-coordinates of target points
	!> - [IN]  grid_R   - R-coordinates of grid points
	!> - [IN]  grid_Z   - Z-coordinates of grid points
	!> - [IN]  myrank   - particular rank of a process
	!> - [IN]  nprocs   - number of MPI processes
	!> - [OUT] index_R  - R-indices of grid points closest but below target point 
	!> - [OUT] index_Z  - Z-indices of grid points closest but below target point 
	!> - [OUT] dist_R   - distance from target point to closest but smallergrid
	!>                  point in R
	!> - [OUT] dist_Z   - distance from target point to closest but smallergrid
	!>                  point in Z
	!> 
	!> EXAMPLE:
	!> - Example:
	!>      grid_R   = [0.  0.25 0.5  0.75 1.]
	!>      grid_Z   = [0.  0.5  1.   1.5  2.]
	!>      targets_R = [0.4 0.8]
	!>      targets_Z = [0.  0.4]
	!>      index_R, index_Z, dist_R, dist_Z = &
	!>                          grid_neighbours(targets_R, targets_Z, grid_R, grid_Z)
	!>      $ index_R = [1 3]
	!>      $ index_Z = [1 0]
	!>      $ dist_R  = [0.15 0.05]
	!>      $ dist_Z  = [0.   0.4 ]
	!>      $ grid_R[index_R] = [0.25 0.75]
	!>      $ grid_Z[index_Z] = [0.   0.  ]
	!> 
	!> REVISION HISTORY:
	!> 22 02 2021 - Carlos Romero
	!> 
	!> TO DO:
	!> - Why the int() on the index? Because it is giving me floats, Idk why!!!
	!> - Test
	!> - 
	!> 
	!------------------------------------------------------------------------------
	subroutine grid_neighbours(targets_R, targets_Z, grid_R, grid_Z, myrank, nprocs, index_R, index_Z, dist_R, dist_Z)
		
		real(dp), intent(in) :: targets_R(:), &
								targets_Z(:), &
								grid_R(:),    &
								grid_Z(:)
											 
	    integer, intent(in) :: myrank, &
	    					   nprocs
		
		integer, allocatable, intent(out)  :: index_R(:), &
							   				  index_Z(:)
							   	
		real(dp), allocatable, intent(out) :: dist_R(:),  &
											  dist_Z(:)

		integer :: ntargets, &
				   m,        &
				   r,        &
				   z
		
		ntargets = size(targets_R)
		allocate(index_R(ntargets))
		allocate(index_R(ntargets))
		allocate(dist_R(ntargets))
		allocate(dist_R(ntargets))
        
		do m = myrank+1, ntargets, nprocs
		
			r = minloc(abs(targets_R(m) - grid_R), 1)
			if (targets_R(m) - grid_R(r) < 0) then
				r = r - 1
			end if
				
			!abs_diff_z = 
			z = minloc(abs(targets_Z(m) - grid_Z), 1)
			if (targets_Z(m) - grid_Z(z) < 0) then
				z = z - 1
			end if
			
			index_R(m) = r
			index_Z(m) = z
			dist_R(m)  = targets_R(m) - grid_R(r)
			dist_Z(m)  = targets_Z(m) - grid_Z(z)
		
		end do
		
			
	end subroutine grid_neighbours

	!------------------------------------------------------------------------------
	! SUBROUTINE - interpolate_field
	!------------------------------------------------------------------------------
	!> 
	!> DESCRIPTION:
	!> Provided indices of a 2-dimensional grid in space, previously computed
	!> distances of target points to points in that grid, and a tangent field,
	!> computes for an interpolation of the field in the target point, which is not
	!> necessarily in the grid.
	!> 
	!> INPUT/OUTPUT:
	!> - [IN]  index_R - R-indices of grid points closest but below target
	!> - [IN]  index_Z - Z-indices of grid points closest but below target
	!> - [IN]  dist_R  - distance from target point to closest but smaller grid
	!>                 point in R
	!> - [IN]  dist_Z  - distance from target point to closest but smaller grid
	!>                 point in Z
	!> - [IN]  F_grid  - original scalar field matrix
	!> - [IN]  step_R  - distance between two grid points in the R-coordinate
	!> - [IN]  step_Z  - distance between two grid points in the Z-coordinate
	!> - [IN]  myrank  - particular rank of a process
	!> - [IN]  nprocs  - number of MPI processes
	!> - [OUT] F_inter - interpolated scalar field for targets
	!> 
	!> EXAMPLE:
	!> - Example:
	!>
	!> REVISION HISTORY:
	!> 25 02 2021 - Carlos Romero
	!> 
	!> TO DO:
	!> - Change to Fortran!!!
	!> - Create example
	!> -
	!> 
	!------------------------------------------------------------------------------
	subroutine interpolate_field(index_R, index_Z, dist_R, dist_Z, F_grid, step_R, step_Z, myrank, nprocs, F_inter)
	
		integer, intent(in) :: index_R(:), &
							   index_Z(:)
							   
		real(dp), intent(in) :: dist_R(:),   &
								dist_Z(:),   &
								F_grid(:,:), &
								step_R,   &
								step_Z
								
	    integer, intent(in) :: myrank, &
	    					   nprocs
	    					   
		real(dp), allocatable, intent(out) :: F_inter(:)
		
		integer :: ntargets, &
				   m,        &
				   r,        &
				   z
		
		ntargets = size(index_R)
		allocate(F_inter(ntargets))
        
		do m = myrank+1, ntargets, nprocs
		
			r = index_R(m)
			z = index_Z(m)
			F_inter(m) = (                                                  &
		                  F_grid(r,z)     * (1-dist_R(m)) * (1-dist_Z(m)) + &
		                  F_grid(r+1,z)   * (1-dist_R(m)) * dist_Z(m)     + &
		                  F_grid(r,z+1)   * dist_R(m)     * (1-dist_Z(m)) + &
		                  F_grid(r+1,z+1) * dist_R(m)     * dist_Z(m)       &
		                 )                                                  &
		                  / (step_R * step_Z)    
		                 
		end do
		                    
	end subroutine interpolate_field

	!------------------------------------------------------------------------------
	! SUBROUTINE - step_phi
	!------------------------------------------------------------------------------
	!> 
	!> DESCRIPTION:
	!> Traces discretized field lines in cylindrical geometry: field line equation.
	!> 
	!> INPUT/OUTPUT:
	!> - [IN]     B_R_inter   - interpolated radial magnetic field for targets
	!> - [IN]     B_phi_inter - interpolated toroidal magnetic field for targets
	!> - [IN]     B_Z_inter   - interpolated vertical magnetic field for targets
	!> - [IN]     nphi        - number of points of the grid in the phi direction
	!> - [INOUT]  targets_R   - R-coordinates of target points
	!> - [INOUT]  targets_Z   - Z-coordinates of target points
	!> 
	!> EXAMPLE:
	!> - Example: unnecessary
	!>
	!> REVISION HISTORY:
	!> 24 02 2021 - Carlos Romero
	!> 
	!> TO DO:
	!> - What if B_phi_inter[k] = 0?
	!> -
	!> 
	!------------------------------------------------------------------------------
	subroutine step_phi(B_R_inter, B_phi_inter, B_Z_inter, nphi, targets_R, targets_Z)
	

							       
		real(dp), intent(in) :: B_R_inter(:),   &
								B_phi_inter(:), &
								B_Z_inter(:)
								
		integer, intent(in) :: nphi
		
		real(dp), intent(inout) :: targets_R(:), &
							       targets_Z(:)
		
		real(dp) :: delta_phi
		
		delta_phi = pi2 / (nphi - 1)
		
		targets_R = targets_R + delta_phi * ( B_R_inter / B_phi_inter ) * targets_R
		targets_Z = targets_Z + delta_phi * ( B_Z_inter / B_phi_inter ) * targets_R
		
	end subroutine step_phi

	!------------------------------------------------------------------------------
	! SUBROUTINE - delete_outer_particles
	!------------------------------------------------------------------------------
	!> 
	!> DESCRIPTION:
	!> Deletes particles that are outside simulation boundaries.
	!> 
	!> INPUT/OUTPUT:
	!> - [IN]    rmin      - minimum radial coordinate
	!> - [IN]    rmax      - maximum radial coordinate
	!> - [IN]    zmin      - minimum vertical coordinate
	!> - [IN]    zmax      - maximum vertical coordinate
	!> - [IN]    step_R    - distance between two grid points in the R-coordinate
	!> - [IN]    step_Z    - distance between two grid points in the Z-coordinate
	!> - [INOUT] targets_R - R-coordinates of target points
	!> - [INOUT] targets_Z - Z-coordinates of target points
	!> 
	!> EXAMPLE:
	!> - Example: omitting unnecessary variables for illustration
	!>      rmin = 0
	!>      targets_R = [2. 8.3 -1.3 2.7]
	!>      targets_R = delete_outer_particles(targets_R, ~, rmin, ~, ~, ~)
	!>      $ targets_R = [2. 8.3 2.7]
	!>
	!> REVISION HISTORY:
	!> 24 02 2021 - Carlos Romero
	!> 
	!> TO DO:
	!> - Reallocate only if bigger array size
	!> -
	!> 
	!------------------------------------------------------------------------------
	subroutine delete_outer_particles(rmin, rmax, zmin, zmax, step_R, step_Z, targets_R, targets_Z)
	
		real(dp), intent(in) :: rmin,   &
								rmax,   &
								zmin,   &
								zmax,   &
								step_R, &
								step_Z
		
		real(dp), intent(inout) :: targets_R(:), &
							       targets_Z(:)
		
		integer :: ntargets
		
		logical, allocatable :: condition_R(:), &
							    condition_Z(:), &
							    condition(:)
							    
		ntargets = size(targets_R)
		allocate(condition_R(ntargets))
		allocate(condition_Z(ntargets))
		allocate(condition(ntargets))
		
		condition_R = (rmin < targets_R) .and. (targets_R < rmax)
		condition_Z = (zmin < targets_Z) .and. (targets_Z < zmax)
		condition   = condition_R .and. condition_Z
		
		targets_R   = pack(targets_R, condition)
		targets_Z   = pack(targets_Z, condition)
		
	end subroutine delete_outer_particles
	
	!------------------------------------------------------------------------------
	! SUBROUTINE - save_intersection
	!------------------------------------------------------------------------------
	!> 
	!> DESCRIPTION:
	!> Appends target vectors only when a turn has been completed (see conditional)
	!> in loop. 
	!> 
	!> INPUT/OUTPUT:
	!> - [IN]    targets_R   - R-coordinates of target points
	!> - [IN]    targets_Z   - Z-coordinates of target points
	!> - [INOUT] intersect_R - appended R-coordinates of target points at every turn
	!> - [INOUT] intersect_Z - appended Z-coordinates of target points at every turn
	!> 
	!> EXAMPLE:
	!> - Example: unnecessary
	!>
	!> REVISION HISTORY:
	!> 24 02 2021 - Carlos Romero
	!> 
	!> TO DO:
	!> - 
	!> 
	!------------------------------------------------------------------------------
	subroutine save_intersection(targets_R, targets_Z, intersect_R, intersect_Z)
	
		real(dp), intent(in) :: targets_R(:), &
							    targets_Z(:)
							    
		real(dp), allocatable, intent(inout) :: intersect_R(:), &
					       					    intersect_Z(:)
	
		integer :: length_old,    &
				   length_append, &
				   length_new
		
		length_old    = size(intersect_R)
		length_append = size(targets_R)
		length_new    = length_old + length_append
		
		allocate(intersect_R(length_new))
		allocate(intersect_R(length_new))

		intersect_R(length_old+1:length_new) = targets_R
		intersect_Z(length_old+1:length_new) = targets_Z

	end subroutine save_intersection

	!------------------------------------------------------------------------------
	
end module field_line_mod
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
