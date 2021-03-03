!------------------------------------------------------------------------------
! PROGRAM - field_line_tracer.f90
!------------------------------------------------------------------------------
!> 
!> DESCRIPTION:
!> Contains necessary constants and subroutines to trace field lines in 2- or 
!> 3-dimensional systems to visualize flow structures.
!>
!> Tracing field lines requires:
!> - starting point of the curve: consider a line or a matrix of points
!> - vector field considered pointwise tangent to the curve
!> 
!> As we are dealing with a discretized space, both the space and the vector 
!> field are given in a multidimensional grid. To compute a realistic field
!> line, we will need to do:
!> - multidimensional interpolation to obtain precise field values
!> - coordinate transformations to account for specific geometries
!> 
!> REVISION HISTORY:
!> 24 02 2021 - Carlos Romero
!> 
!> TO DO:
!> - What if R and Z dimensions are different?
!> -
!> 
!------------------------------------------------------------------------------
program field_line_tracer

!------------------------------------------------------------------------------
! Modules
use netcdf
use mpi
use field_line_mod
use general_mod

implicit none

!------------------------------------------------------------------------------
! Declaration
! - File path
character(len=50) :: input_file  = 'example.nc', &
! - Targets' mode
					 mode        = 'matrix',     &
! - Output file name
					 output_file = 'test.nc'
					 
! - Number of targets, turns and steps
integer :: ntargets1D    = 1000, &
		   nturns        = 100,  &
		   size_new_grid = 1000, &
		   nsteps,   &
		   ntargets, &
		   angle,	 &
		   turn, 	 &
		   i,        &
		   j,        &
		   step,     &
		   m,        &
! - netCDF
		   ncid,       &
		   r_dimid,    &
		   z_dimid,    &
		   d_dimid(2), &
		   varid1,     &
		   varid2,     &
		   varid3,     &
! - MPI
		   nprocs,     &
		   myrank,     &
		   ierr
		   
! - Measuring time
real(dp) :: tic,      &
			toc,      &
			tic_do,   &
			toc_do,   &
			time_exp, & 
! - Boundaries
			rminb,  &
		    rmaxb,  &
		    zminb,  &
		    zmaxb,  &
! - Spatial step of the grid
			step_R, &
			step_Z
		
! -> Are these declarations necessary if they are already done inside the 
! correspondent subroutines?
! - Magnetic, velocity and pressure field tensors

real(dp), allocatable :: B_R(:,:,:),     &
						 B_phi(:,:,:),   &
						 B_Z(:,:,:),     &
! - Grid in the correspondent dimension
						 grid_R(:),      &
						 grid_Z(:),      &					 
! - Targets whose magnetic field lines will be traced 
						 targets_R(:),   &
						 targets_Z(:),   &
! - Distances target to lower grid point
						 dist_R(:),      &
						 dist_Z(:),      &
! - Interpolated field for targets
						 B_R_inter(:),   &
						 B_phi_inter(:), &
						 B_Z_inter(:),   &
! - Trajectories intersections with toroidal phi = 0
						 intersect_R(:), &
						 intersect_Z(:), &
! - Plotting intersection density in thinner grid
						 grid_R_update(:), &
						 grid_Z_update(:), &
						 density(:,:) 
						 
! - Index of the lower grid point to target
integer, allocatable ::  index_R(:),     &
						 index_Z(:)
				
		
!------------------------------------------------------------------------------
! Tic
call cpu_time(tic)

!------------------------------------------------------------------------------
! Welcome and terminal output
print *
print *,'==============================='
print *,' Welcome to field line tracer! '
print *,'==============================='
print *

print *
print *,'==================================='
print *,' You have selected:                '
print *
print *,'  	- input file  -> ', input_file 
print *,'   - mode        -> ', mode
print *,'   - output file -> ', output_file
print *,'   - ntargets1D  -> ', ntargets1D
print *,'   - nturns      -> ', nturns
print *
print *,'==================================='
print *

!------------------------------------------------------------------------------
! Read data
call read_data_hint(input_file, rminb, rmaxb, zminb, zmaxb, B_R, B_phi, B_Z)

print *
print *,'========================='
print *,' Data successfully read! '
print *,'========================='
print *

!------------------------------------------------------------------------------
! Extracting dimensions from data and creating original grid
call linspace(rminb, rmaxb, nr, grid_R)
call linspace(zminb, zmaxb, nz, grid_Z)

step_R = (rmaxb - rminb ) / (nr - 1)
step_Z = (zmaxb - zminb ) / (nz - 1)

!------------------------------------------------------------------------------
! Creating initital target points
call create_targets(mode, ntargets1D, rminb + 2*step_R, &
									  rmaxb - 2*step_R, &
									  zminb + 2*step_Z, &
									  zmaxb - 2*step_Z, targets_R, targets_Z)

!------------------------------------------------------------------------------
! Iterate phi a certain number of toroidal turns and save at every turn
ntargets = size(targets_R)
allocate(intersect_R(ntargets))
allocate(intersect_Z(ntargets))

intersect_R = targets_R
intersect_Z = targets_Z

nsteps = nturns * nphi

time_exp = 0.
turn     = 0

!------------------------------------------------------------------------------		 
! MPI Initialization   
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

! - This loop cannot be parallelized
do step = 1, nsteps

	!------------------------------------------------------------------------------
	! Barrier parallel processes so that this non-parallelized loop is fine
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)

	! Tic
	call cpu_time(tic_do)
	
	angle = mod(step, nphi)
	
	!------------------------------------------------------------------------------
	! Parallel
	! - Compute the grid neighbours to our targets
	call grid_neighbours(targets_R, targets_Z, grid_R, grid_Z, myrank, nprocs, index_R, index_Z, dist_R, dist_Z)
	
	! - Interpolate radial magnetic field
	call interpolate_field(index_R, index_Z, dist_R, dist_Z, B_R(angle,:,:),   step_R, step_Z, myrank, nprocs, B_R_inter)
	
	! - Interpolate toroidal magnetic field
	call interpolate_field(index_R, index_Z, dist_R, dist_Z, B_phi(angle,:,:), step_R, step_Z, myrank, nprocs, B_phi_inter)
	
	! - Interpolate vertical magnetic field
	call interpolate_field(index_R, index_Z, dist_R, dist_Z, B_Z(angle,:,:),   step_R, step_Z, myrank, nprocs, B_Z_inter)
	
	
	!!!!!!!!!!!!!!!!!!!!!
	if (myrank == 0) then 
		! - Solve finite differences
		call step_phi(B_R_inter, B_phi_inter, B_Z_inter, nphi, targets_R, targets_Z)
		
		! - Delete particles outside simulation box
		call delete_outer_particles(rminb, rmaxb, zminb, zmaxb, step_R, step_Z, targets_R, targets_Z)
		
		!------------------------------------------------------------------------------
		! Save intersections at every turn (intersection with phi = 0 - not every step)
		if (angle == 0) then
		
			call save_intersection(targets_R, targets_Z, intersect_R, intersect_Z)
			
			! Update turn
			turn = turn + 1
			
			print *
			print *,'==========================================='
			print *,' Finished turn number ', turn, '/', nturns
			print *,'==========================================='
			print *
			
			! Toc:
			call cpu_time(toc_do)
			
			! Mean of expected time (using recursive expected time to be more precise)
			time_exp = (nturns - turn) * ( time_exp * (turn - 1) + (toc_do - tic_do) ) / turn
			
			print *
			print *,'==========================================='
			print *,' Expected remaining time: ', time_exp, ' s '
			print *,'==========================================='
			print *

		end if
	end if
	
end do

!------------------------------------------------------------------------------
! Barrier parallel processes so that this non-parallelized loop is fine
call MPI_BARRIER(MPI_COMM_WORLD, ierr)


!!!!!!!!!!!!!!!!!!!!!
if (myrank == 0) then
	!------------------------------------------------------------------------------
	! Compute and save intersection density with a thinner grid
	allocate(grid_R_update(size_new_grid))
	allocate(grid_Z_update(size_new_grid))
	allocate(density(size_new_grid,size_new_grid))

	call linspace(rminb, rmaxb, size_new_grid, grid_R_update)
	call linspace(zminb, zmaxb, size_new_grid, grid_Z_update)
end if


! Parallel loop?
do m = myrank+1, size(intersect_R), nprocs

	print *
	print *,'====================================================='
	print *,' Computing density in process ', myrank, '/', nprocs 	
	print *,'====================================================='
	print *
	
	i = minloc(abs(grid_R_update - intersect_R(m)), 1)
	j = minloc(abs(grid_Z_update - intersect_Z(m)), 1)
	
	density(i,j) = density(i,j) + 1
	
end do

!!!!!!!!!!!!!!!!!!!!!
if (myrank == 0) then
	!------------------------------------------------------------------------------
	! Write to variable

	print *
	print *,'======================================'
	print *,' Writing intersection density to file '
	print *,'======================================'
	print *
		
	call check( nf90_create(output_file, NF90_CLOBBER, ncid) )

	call check( nf90_def_dim(ncid, "R", size_new_grid, r_dimid) )
	call check( nf90_def_dim(ncid, "Z", size_new_grid, z_dimid) )

	call check( nf90_def_var(ncid, "grid_R", NF90_REAL8, r_dimid, varid1) )
	call check( nf90_def_var(ncid, "grid_Z", NF90_REAL8, r_dimid, varid2) )

	d_dimid =  (/r_dimid, z_dimid/)

	call check( nf90_def_var(ncid, "density", NF90_REAL8, d_dimid, varid3) )

	call check( nf90_enddef(ncid) )

	call check( nf90_put_var(ncid, varid1, grid_R_update) )
	call check( nf90_put_var(ncid, varid2, grid_R_update) )
	call check( nf90_put_var(ncid, varid3, density) )

	call check( nf90_close(ncid) )



	!------------------------------------------------------------------------------
	! Ending program
	call cpu_time(toc)
	print *
	print *,'============================================='
	print *,' Total elapsed time: ', toc - tic, ' seconds '
	print *,'============================================='
	print *
end if


!------------------------------------------------------------------------------
! MPI Finalize
call MPI_FINALIZE(ierr)

!------------------------------------------------------------------------------
! End program
end program field_line_tracer

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
