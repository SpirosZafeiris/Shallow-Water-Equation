subroutine map_global
   use sizes
   use space
   implicit none
   integer :: i,j,var
   allocate(gl(Nxnod,Nynod,2))
   do i=1,Nxnod
      do j=1,Nynod
         do var=1,2
            gl(i,j,var) = dof*(var-1) + i + (j-1)*Nxnod
         enddo
      enddo
   enddo
endsubroutine map_global


subroutine explicit
   use types
   use sizes
   use time
   use space
   implicit none
   external wave_source
   integer :: i,j
   double precision :: temp,rhs0,rhs
   do i=1,Nx+1 ! calculate first & second derivatives
      if (t.eq.0) exit
      do j=1,Ny+1
         if (PML(i,j).eq.0) then
            u(i,j) = u(i,j) + dt*du(i,j)! + dt22 * ddu(i,j)
            if (abs(u(i,j)).gt.100) then
               print *, 'diverges'
               stop
            endif
            !du(i,j) = du(i,j) + dt * ddu(i,j)
         elseif (PML(i,j).eq.1) then
            u1(i,j) = u1(i,j) + dt*du1(i,j)! + dt22 * ddu1(i,j)
            a2(i,j) = a2(i,j) + dt*da2(i,j)! + dt22 * dda2(i,j)
            a3(i,j) = a3(i,j) + dt*da3(i,j)! + dt22 * dda3(i,j)
            u2(i,j) = u2(i,j) + dt*(a2(i,j) + s1(i,j)*u2(i,j))
            u3(i,j) = u3(i,j) + dt*(a3(i,j) + s2(i,j)*u3(i,j))
            u(i,j) = u1(i,j) + u2(i,j) + u3(i,j)
            if (i.eq.(Nx+1).or.j.eq.(Ny+1).or.i.eq.1.or.j.eq.1) then
               u(i,j) = 0.d0 ! Dirichlet
            endif
            !if (i.eq.Nx+1) u(i,j) = u(Nx,j) !Neumann
            !if (j.eq.Ny+1) u(i,j) = u(i,Ny)
            !if (i.eq.1)    u(i,j) = u(2,j)
            !if (j.eq.1)    u(i,j) = u(i,2)

            !du1(i,j) = du1(i,j) + dt * ddu1(i,j)
            !da2(i,j) = da2(i,j) + dt * dda2(i,j)
            !da3(i,j) = da3(i,j) + dt * dda3(i,j)
            !du2(i,j) = du2(i,j) + dt * ddu2(i,j)
            !du3(i,j) = du3(i,j) + dt * ddu3(i,j)
         endif
         !phi(i,j) = phi(i,j) + dt*(-9.81d0*u(i,j))
      enddo
   enddo
   do i=1,Nx+1
      do j=1,Ny+1
         if (i.gt.1.and.i.lt.(Nx+1)) then
            ux  = (u(i+1,j) - u(i-1,j)) / (dx2)
            uxx = (u(i-1,j) - 2*u(i,j) + u(i+1,j)) / (dx22)
         elseif (i.eq.1) then
            ux  = (u(2,j) - u(1,j)) / dx
            uxx = (u(1,j) - 2*u(2,j) + u(3,j)) / (dx22)
         elseif (i.eq.(Nx+1)) then
            ux  = (u(Nx+1,j) - u(Nx,j)) / dx
            uxx = (u(Nx-1,j) - 2*u(Nx,j) + u(Nx+1,j)) / (dx22)
         endif
         if (j.gt.1.and.j.lt.(Ny+1)) then
            uy  = (u(i,j+1) - u(i,j-1)) / (dy2)
            uyy = (u(i,j-1) - 2*u(i,j) + u(i,j+1)) / (dy22)
         elseif (j.eq.1) then
            uy  = (u(i,2) - u(i,1)) / dy
            uyy = (u(i,1) - 2*u(i,2) + u(i,3)) / (dy22)
         elseif (j.eq.(Ny+1)) then
            uy  = (u(i,Ny+1) - u(i,Ny)) / dy
            uyy = (u(i,Ny-1) - 2*u(i,Ny) + u(i,Ny+1)) / (dy22)
         endif
         ! ---------calculate second derivatives
         if (PML(i,j).eq.0) then
            call wave_source(xnod(i,j), ynod(i,j),t, temp)
            rhs = h_x(i,j)*ux+h_y(i,j)*uy + h(i,j)*(uxx+uyy) + temp
            rhs = rhs*9.81d0
            du(i,j) = du(i,j) + dt*rhs
            !ddu(i,j) = ddu(i,j)*9.81d0
         elseif (PML(i,j).eq.1) then
            rhs0 = ( h_x(i,j)*ux+h_y(i,j)*uy + h(i,j)*(uxx+uyy) ) * 9.81d0
            !rhs0 =  h(i,j)*(uxx+uyy) * 9.81d0
            rhs0 = rhs0 + 2*(s1(i,j) + s2(i,j))*du1(i,j) &
                      - (s1(i,j)**2 + s2(i,j)**2) * u1(i,j)
            du1(i,j) = du1(i,j) + dt * rhs0
            rhs0 = 9.81d0*s1_x(i,j)*h(i,j)*ux
            rhs0 = rhs0 + 2*s1(i,j)*da2(i,j) - (s1(i,j)**2)*a2(i,j)
            da2(i,j) = da2(i,j) + dt*rhs0
            !dda2(i,j) = rhs0 + 2*s1(i,j)*da2(i,j) - (s1(i,j)**2)*a2(i,j)

            rhs0 = 9.81d0*s2_y(i,j)*h(i,j)*uy
            rhs0 = rhs0 + 2*s2(i,j)*da3(i,j) - (s2(i,j)**2)*a3(i,j)
            da3(i,j) = da3(i,j) + dt*rhs0
         endif
      enddo
   enddo
endsubroutine explicit

subroutine calc_ham
   use Param
   use types
   use sizes
   use time
   use space
   implicit none
   integer :: i,j
   double precision :: phix, phiy
   do i=1,Nxnod
      do j=1,Nynod
         if (i.gt.1.and.i.lt.(Nxnod)) then
            phix = (phi(i+1,j) - phi(i-1,j)) / (dx2)
         elseif (i.eq.1) then
            phix = (phi(2,j) - phi(1,j)) / dx
         elseif (i.eq.(Nxnod)) then
            phix = (u(Nx+1,j) - u(Nx,j)) / dx
         endif
         if (j.gt.1.and.j.lt.(Nynod)) then
            phiy  = (phi(i,j+1) - phi(i,j-1)) / (dy2)
         elseif (j.eq.1) then
            phiy  = (phi(i,2) - phi(i,1)) / dy
         elseif (j.eq.(Nynod)) then
            phiy  = (phi(i,Ny+1) - phi(i,Ny)) / dy
         endif
         
         Ham(i,j) = 0.5*rho * ( h(i,j)*(phix**2+phiy**2)+9.81d0*u(i,j)**2)
      enddo
   enddo
endsubroutine calc_ham

subroutine integrate_ham
   use Param
   use types
   use sizes
   use time
   use space
   implicit none
   integer :: i,j
   double precision :: vol1, vol2, vol, sm
   vol = 0.d0
   sm  = 0.d0
   do i=1,Nx
      do j=1,Ny
         if (PML(i,j).eq.1) cycle
         vol1 = (Ham(i,j)+Ham(i+1,j))*0.5*dx
         vol2 = (Ham(i,j+1)+Ham(i+1,j+1))*0.5*dx
         vol  = (vol1 + vol2)*0.5*dy
         sm = sm + vol
      enddo
   enddo
   total_energy = sm
endsubroutine integrate_ham
