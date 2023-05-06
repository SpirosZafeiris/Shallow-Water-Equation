subroutine write_outfile(counter)
   use types
   use time 
   use space
   use sizes
   implicit none
   integer,      intent(in)  :: counter
   character(:), allocatable :: varhead
   character*80 :: outfile
   integer      :: i,j
   write(outfile,'(a,I10.10,a)') 'out',counter,'.dat'
   open (1, file = outfile)
      varhead = 'Variables = "X","Y"'
      varhead = varhead//',"Elevation"'
      varhead = varhead//',"Damping"'
      !varhead = varhead//',"Damping y"'
      varhead = varhead//',"bathymetry"'
      !varhead = varhead//',"Elevation temporal velocity"'
      varhead = varhead//',"PML"'
      varhead = varhead//',"phi"'
      varhead = varhead//',"Hamiltonian"'
      write(1,'(a)') varhead
      write(1, '(a,I6.6,a,I6.6,a)') 'ZONE T="BIG ZONE", I=',Nxnod,', J=',Nynod,', DATAPACKING=POINT'
      write(1,*) 'Solutiontime=', t
      do i=1,Nx+1
         do j=1,Ny+1
            write(1, *) xnod(i,j),ynod(i,j),u(i,j),' ', &
                        s12(i,j),' ',h(i,j),' ',PML(i,j),' ',phi(i,j),' ',Ham(i,j)
            !write(1, "(8(F64.32,a1))") xnod(i,j),' ',ynod(i,j),' ',u(i,j),' ', &
            !            s12(i,j),' ',h(i,j),' ',PML(i,j),' ',phi(i,j),' ',Ham(i,j)
         enddo !'(6(F8.4))'
      enddo
   close(1)
endsubroutine write_outfile

