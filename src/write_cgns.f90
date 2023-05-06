subroutine write_cgns(counter)
   use, intrinsic :: iso_fortran_env
   use cgns
   use sizes
   use types
   implicit none
   integer, intent(in) :: counter
   integer(int64), allocatable :: isize(:,:)
   integer :: ier, cgfile, cgbase, cgzone,cgcoord, cgsol, cgsola
   character*80 :: outfile, solname
   write(outfile,'(a,I10.10,a)') 'out',counter,'.cgns'
   call cg_open_f(trim(adjustl(outfile)), CG_MODE_WRITE, cgfile, ier)
   call cg_base_write_f(cgfile, 'Base', 2, 2, cgbase, ier)
   allocate(isize(1,6))
   isize(1,1) = Nxnod
   isize(1,2) = Nynod
   isize(1,3) = Nx
   isize(1,4) = Ny
   isize(1,5:6) = 0
   call cg_zone_write_f(cgfile, cgbase, 'Zone', isize, Structured, cgzone, ier)
   call cg_simulation_type_write_f(cgfile, cgbase, TimeAccurate, ier)	
   call cg_coord_write_f(cgfile, cgbase, cgzone, RealDouble, 'CoordinateX', xnod, cgcoord, ier)
   call cg_coord_write_f(cgfile, cgbase, cgzone, RealDouble, 'CoordinateY', ynod, cgcoord, ier)

   write(solname,'(a)') 'FlowSolution'
   call cg_sol_write_f(cgfile, cgbase, cgzone, trim(solname), Vertex, cgsol, ier)
   call cg_field_write_f(cgfile, cgbase, cgzone, cgsol, RealDouble, "Elevation", u, cgsola, ier)
   call cg_field_write_f(cgfile, cgbase, cgzone, cgsol, RealDouble, "Damping", s12, cgsola, ier)
   call cg_field_write_f(cgfile, cgbase, cgzone, cgsol, RealDouble, "Bathymetry", h, cgsola, ier)
   call cg_field_write_f(cgfile, cgbase, cgzone, cgsol, Integer, "PML", PML, cgsola, ier)
   call cg_field_write_f(cgfile, cgbase, cgzone, cgsol, RealDouble, "phi", phi, cgsola, ier)
   call cg_field_write_f(cgfile, cgbase, cgzone, cgsol, RealDouble, "Hamiltonian", Ham, cgsola, ier)

   call cg_close_f(cgfile, ier)
endsubroutine write_cgns






