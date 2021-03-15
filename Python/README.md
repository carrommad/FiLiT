# Available versions
Two versions are provided, differing in structure but with very similar content (functions and plotting).

## 1. Turns -> Particles
The functions named ***field_line_tracer*** are structured in this way: given a step size for a given major radius, the toroidal angle increment is computed so that the length increment for a particle in a major radius positions coincides with a desired value (in mm).

Pros:
- Except for two *for* loops in the *particles* or *targets* (see ***field_line_tracer.mod***: *grid_neighbours* and *interpolate_field*), the rest of the operations is element-wise

Cons:
- Improvements need to be made when saving the intersections so that we can distinguish **which** particle is intersecting **where** (easy improvement)
- Step size not the same for every particle as it depends on its radius (cylindrical geometry)
- Difficult to parallelize (I've been told)

## 2. Particles -> Turns
The functions name ***filit*** are structured the opposite way: we loop over the particles and compute for each step the toroidal angle that we need to advance in order to fulfill the step size condition.

Pros:
- Seems to be a more intuitive approach
- Easy to follow specific targets
- Easy to parallelize (I've been told)

Cons:
- Slower than the other version, I believe because there are no element-wise operations
