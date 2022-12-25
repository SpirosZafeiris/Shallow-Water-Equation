Program Swe
   use, intrinsic :: iso_fortran_env, only : real64, int32
   use types
   use time 
   use space
   use sizes
   use ellipsis
   implicit none
   external linspace
   external write_outfile
   external explicit
   external calc_Ham
   external integrate_ham
   integer :: i,j,coun,counter
   real(real64), allocatable :: x(:),y(:)
   real(real64) :: temp, len_PML, t1, t2, full_time
   !open(4,file='input')
   !read(4,*) a
   !read(4,*) b
   !close(4)
   !print *, a, b
      !! determine total number of segments per direction
   Nx = 300
   Ny = 300
   Nxnod = Nx+1
   Nynod = Ny+1
   dof = Nxnod*Nynod
      !! determine boundary positions
   xedge(1) = -300._real64; xedge(2)=300._real64
   yedge(1) = -250._real64; yedge(2)=250._real64
   len_PML = 60.d0

      ! determine PML coefficient
   const = -2.5e-4 ! for quadratic sigma distribution
   const = -8._real64 ! for inverse linear sigma distribution
   !const = -8.5e-6 ! for qubic sigma distribution
   !const = -2.5e-3 ! for qubic sigma distribution
      !! value allocation
   allocate(x(Nx+1),y(Ny+1),xnod(Nx+1,Ny+1), ynod(Nx+1,Ny+1))
   call linspace(Nx+1,xedge(1),xedge(2), x)
   call linspace(Ny+1,yedge(1),yedge(2), y)
   ! ---------------Construct grid ---------------------
   do j=1,Ny+1
      xnod(:,j) = x
   enddo

   do i=1,Nx+1
      ynod(i,:) = y
   enddo

   dx = x(2)-x(1)
      !! grid x-spacing 
   dy = y(2)-y(1)
      !! grid y-spacing 
   deallocate(x, y)

   Npml = len_PML/dx+1


   print *, 'number of PML layers per region', Npml-1
   dt = 1e-2
   dt22 = dt**2 / 2.d0
   dx2 = 2*dx; dx22 = dx**2; dy2 = 2*dy; dy22 = dy**2

   lenx(1) = abs(abs(xedge(1))-len_PML)
   lenx(2) = abs(abs(xedge(2))-len_PML)
   leny(1) = abs(abs(yedge(1))-len_PML)
   leny(2) = abs(abs(yedge(2))-len_PML)


   ! ------------initial values--------------------------
   allocate(u(Nx+1,Ny+1) , du(Nx+1,Ny+1) , ddu(Nx+1,Ny+1),  &
            u1(Nx+1,Ny+1), du1(Nx+1,Ny+1), ddu1(Nx+1,Ny+1), &
            u2(Nx+1,Ny+1), du2(Nx+1,Ny+1), ddu2(Nx+1,Ny+1), &
            u3(Nx+1,Ny+1), du3(Nx+1,Ny+1), ddu3(Nx+1,Ny+1), &
            a2(Nx+1,Ny+1), da2(Nx+1,Ny+1), dda2(Nx+1,Ny+1), &
            a3(Nx+1,Ny+1), da3(Nx+1,Ny+1), dda3(Nx+1,Ny+1), &
            h(Nx+1,Ny+1) , h_x(Nx+1,Ny+1), h_y(Nx+1,Ny+1) , &
            PML(Nx+1,Ny+1),phi(Nxnod,Nynod),Ham(Nx+1,Ny+1)  )
   PML = 0
   do j=1,Ny+1
      do i=1,Nx+1
            !! initial conditions
         call phi0(xnod(i,j), temp)
         phi(i,j) = temp
         call u0(xnod(i,j),ynod(i,j), temp)!initiate elevation & velocity
         u(i,j)  = temp!; du(i,j)  = 0.d0
         u1(i,j) = temp/3!; du1(i,j) = 0.d0
         u2(i,j) = temp/3!; du2(i,j) = 0.d0
         u3(i,j) = temp/3!; du3(i,j) = 0.d0
         call du0(xnod(i,j),ynod(i,j), temp)
         du(i,j) = temp
         du1(i,j) = temp/3
         du2(i,j) = temp/3
         du3(i,j) = temp/3
         if (i.lt.Npml.or.i.gt.(Nx+1)-Npml) then
            PML(i,j) = 1
            cycle
         endif
         if (j.lt.Npml.or.j.gt.(Ny+1)-Npml) then
            PML(i,j) = 1
         endif
      enddo
   enddo
   do j=1,Nynod
      do i=1,Nxnod
            !! assign depth formation on the grid
         call fh(xnod(i,j),ynod(i,j), temp)!assign depth
         h(i,j) = temp
         if (PML(i,j).eq.0) then
            call fhx(xnod(i,j),ynod(i,j), temp)
            h_x(i,j) = temp
            call fhy(xnod(i,j),ynod(i,j), temp)
            h_y(i,j) = temp
         else
            h_x(i,j) = 0._real64
            h_y(i,j) = 0._real64
         endif
      enddo
   enddo
   allocate(s1(Nx+1,Ny+1)  , s2(Nx+1,Ny+1),  &
            s1_x(Nx+1,Ny+1), s2_y(Nx+1,Ny+1) )
   allocate(s12(Nxnod,Nynod))

   !do i=Npml,Nxnod-Npml
   !   do j=Npml, Nynod-Npml
   !      if (PML(i,j).eq.1) print *, i, j, PML(i,j)
   !   enddo
   !enddo

   s1 = 0.d0; s2 = 0.d0; s1_x = 0.d0; s2_y = 0.d0
   do j=1,Ny+1
      do i=1,Nx+1
            !! calculate PML sinks
         if (PML(i,j).eq.1) then
            if (xnod(i,j).lt.0.d0.and.abs(xnod(i,j)).ge.lenx(1).or. &
                xnod(i,j).gt.0.d0.and.abs(xnod(i,j)).ge.lenx(2)) then
               call fs1(xnod(i,j), temp)
               s1(i,j) = temp
               call fs1_x(xnod(i,j), temp)
               s1_x(i,j) = temp
            endif
            !!if (xnod(i,j).gt.0.d0.and.abs(xnod(i,j)).ge.lenx(2)) then
            !!   call fs1(xnod(i,j), temp)
            !!   s1(i,j) = temp
            !!   call fs1_x(xnod(i,j), temp)
            !!   s1_x(i,j) = temp
            !!   print *, 'ok', s1(i,j)
            !!endif

            if (ynod(i,j).lt.0.d0.and.abs(ynod(i,j)).ge.leny(1)) then
               call fs2(ynod(i,j), temp)
               s2(i,j) = temp
               call fs2_y(ynod(i,j), temp)
               s2_y(i,j) = temp
            endif
            if (ynod(i,j).gt.0.d0.and.abs(ynod(i,j)).ge.leny(2)) then
               call fs2(ynod(i,j), temp)
               s2(i,j) = temp
               call fs2_y(ynod(i,j), temp)
               s2_y(i,j) = temp
            endif
         endif
         s12(i,j) = s1(i,j)+s2(i,j)
      enddo
   enddo
   !allocate(ux(Nx+1,Ny+1), uy(Nx+1,Ny+1),uxx(Nx+1,Ny+1), uyy(Nx+1,Ny+1))
   !ux = 0; uy = 0; uxx = 0; uyy = 0;
   coun = 0
   t = -dt
   counter = -1
   t1=0.d0; t2=0.d0
   full_time = 0.d0
   write(*,'(a)') 'Begin loop'

   !call map_global

   do while (t.lt.89)!1205) 
      counter = counter + 1
      t = t + dt
      ! ---------write output file -----------------------
       !! configure to specify at how many steps an output is created
      if (mod(counter,100).eq.0) then 
         call write_outfile(counter)
         ! call write_cgns(counter) !! uncomment for fast CGNS I/O (HDF5 required)
      endif
      full_time = full_time + (t2-t1)
      write(*,'(a,i10.10,a,f12.7,a,f8.6)')'iteration: ',counter,', time of simulation: ',t,' sec, elapsed time: ',t2-t1
     ! write(75,'(i10.10,1x,f12.7,1x,f8.6,1x,f16.5)')counter,t,t2-t1,Total_energy
      !--------begin main calculations-----------------
      call cpu_time(t1)
      call explicit
      !!!call calc_Ham
      !!!call integrate_ham

      !call general_implicit
      call cpu_time(t2)
   enddo
   write(*,'(a,f22.7,a)') 'Simulation Ended, total elapsed time: ', full_time, ' sec'  


   ! ------functions-------------------------
   contains
      !! velocity initial condition, user defined
      subroutine u0(x,y, yy)
         use Param
         implicit none
         real(real64), intent(in)  :: x, y
         real(real64), intent(out) :: yy
          real(real64) :: A, k, om
         A = 0.1
         k = 0.03*pi
         om = 0.7*pi
         if (x.gt.0.d0.and.x.lt.2.d0) then
         yy = A*cos(k*x)
         else
            yy = 0.d0
         endif
         !yy = 0.1*exp(-((x+100)**2+y**2)/100)
         !yy = 0.1*exp(-(x**2+y**2)/100)
         !yy = 0.d0*x*y
      endsubroutine u0

      !! velocity time derivative initial condition, user defined
      subroutine du0(x,y, yy)
         use Param
         implicit none
         real(real64), intent(in)  :: x, y
         real(real64), intent(out) :: yy
         !!!!double precision              :: A, k, om
         !!!!A = 0.1
         !!!!k = 0.03*pi
         !!!!om = 0.7*pi
         !if (x.gt.-0.d0.and.x.lt.2.d0) then
         !yy = -om*A*sin(k*x)
         !else
         !   yy = 0.d0
         !endif
         yy = 0.d0*x*y
      endsubroutine du0
      
      !! depth equation, user defined
      subroutine fh(x,y, yy)
         use Param
         use ellipsis
         implicit none
         real(real64), intent(in)  :: x,y
         real(real64), intent(out) :: yy
         real(real64) :: x_dash,y_dash,Aw,k
         !a = 10
         !b = 20
         !c = 1
         !yy = h0 ! ellispis
         !x_dash = x - x0
         !y_dash = y - y0
         !if ((x_dash/a)**2+(y_dash/b)**2.lt.1) then
         !   yy = sqrt( 1 - (x_dash/a)**2 - (y_dash/b)**2 )
         !   yy = yy * c ! relative depth
         !   yy = h0 + yy
         !endif! ellipsis
         !Aw = 2 !ripple bed
         !k = pi/10
         !if (x.gt.(-20).and.x.lt.100) then
         !   yy = Aw*sin(k*x)
         !else
         !   yy = 0.d0
         !endif
         !yy = yy + h0
         call slit(x,y, yy)
         yy = yy + h0
         !yy = h0
      endsubroutine fh
      
      !! depth equation spacial derivative-x, fh manually differentiated
      subroutine fhx(x,y, yy)
         use ellipsis
         use space
         real(real64), intent(in)  :: x,y
         real(real64), intent(out) :: yy
         real(real64) :: tmp1,tmp2
         tmp1 = 0
         tmp2 = 0
         !a = 5
         !b = 15
         !c = 1
         !!!x0=20; x_dash = x - x0
         !!!y0=0; y_dash = y - y0
         !if (abs(x).lt.a.and.abs(y).lt.b) then
         yy = 0
         !if ((x_dash/a)**2+(y_dash/b)**2.lt.1) then
  !      !    call fh(x,y, temp)
 !       !    temp = h0 - temp
         !   yy = c**2/sqrt( 1 - (x_dash/a)**2 + (y_dash/b)**2 ) * (-x0)/a
         !endif
         !!!!!Aw = 2
         !!!!!k = pi / 10
         !!!!!if (x.gt.(-20).and.x.lt.100) then
         !!!!!   yy = Aw*k*cos(k*x)
         !!!!!else
         !!!!!   yy = 0.d0
         !!!!!endif
         call fh(x+dx,y, tmp1)
         call fh(x-dx,y, tmp2)
         yy = (tmp1-tmp2)/dx2
         !yy = 0
      endsubroutine fhx
      
      !! depth equation spacial derivative-y, fh manually differentiated
      subroutine fhy(x,y, yy)
         use ellipsis
         use space
         implicit none
         real(real64), intent(in)  :: x,y
         real(real64), intent(out) :: yy
         real(real64) :: x_dash,y_dash,tmp1,tmp2
         !a = 2
         !b = 10
         !c = 2
         x_dash = x - x0
         y_dash = y - y0
         !if (abs(x).lt.a.and.abs(y).lt.b) then
         yy = 0
         !if ((x0/a)**2+(y0/b)**2.lt.1) then
         !   !call fh(x,y, temp)
         !   !temp = h0 - temp
         !   yy = c**2/sqrt( 1 - (x_dash/a)**2 + (y_dash/b)**2 ) * (-y0)/b
         !endif
         call fh(x,y+dy, tmp1)
         call fh(x,y-dy, tmp2)
         yy = (tmp1-tmp2)/dy2
         !yy = 0
      endsubroutine fhy

      !! equation of x-PML sink function, user defined 
      subroutine fs1(x, yy)
         use types
         implicit none
         real(real64), intent(in)  :: x
         real(real64), intent(out) :: yy
         !yy = const*abs(lenx - abs(x)) ! linear
         !yy = const*(lenx-abs(x))**2
         !yy = const*abs((lenx-abs(x))**3)
         yy = const / (abs(xedge(2))+dx - abs(x))
      endsubroutine fs1

      !! equation of y-PML sink function, user defined 
      subroutine fs2(y, yy)
         use types
         implicit none
         real(real64), intent(in) :: y
         real(real64), intent(out) :: yy
         !yy = const*abs(leny - abs(y)) ! linear
         !yy = const*(leny-abs(y))**2  ! quadratic
         !yy = const*abs((leny-abs(y))**3)  ! qubic
         yy = const / (abs(yedge(2))+dy - abs(y)) ! inverse linear
      endsubroutine fs2

      !! equation of x-PML sink function x-derivative, depends on fs1
      subroutine fs1_x(x, yy)
         use types
         implicit none
         real(real64), intent(in) :: x
         real(real64), intent(out) :: yy
         !yy = const * abs(x)
         !yy = 2*const*(abs(x)-lenx)
         !yy = 3*const*(lenx-abs(x))**2
         yy = -const / (abs(xedge(2))+dx - abs(x))**2
      endsubroutine fs1_x

      !! equation of y-PML sink function y-derivative, depends on fs2
      subroutine fs2_y(y, yy)
         use types
         implicit none
         real(real64), intent(in) :: y
         real(real64), intent(out) :: yy
         !yy = const * abs(y)
         !yy = 2*const*(abs(y)-leny)
         !yy = 3*const*(leny-abs(y))**2
         yy = -const / (abs(yedge(2))+dy - abs(y))**2
      endsubroutine fs2_y

      subroutine phi0(x, yy)
         use Param
         implicit none
         real(real64), intent(in)  :: x
         real(real64), intent(out) :: yy
         real(real64)              :: A, k, om
         A = 0.1
         k = 0.03*pi
         om = 0.7*pi
         yy = 0.d0
         yy = A*9.81d0/om *sin(k*x)
      endsubroutine phi0
Endprogram Swe
