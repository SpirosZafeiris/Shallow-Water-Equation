module Param
   double precision, parameter :: pi = acos(-1.d0), limit = 1e-20
   double precision, parameter :: rho = 1024
end module Param

module time
   double precision, save :: dt,t,beta,gamma
   !double precision, allocatable, save :: rhs(:)
   !integer, allocatable :: line(:),column(:), num_col(:), num_col_acc(:)
   double precision, allocatable :: val(:)
   integer :: dof2,dof6
endmodule time

module sizes
   integer, save :: Nx,Ny,Nxnod,Nynod,dof,Npml
endmodule sizes

module space
   double precision, save :: dx,dy,dx2,dx22,dy2,dy22,dt2,dt22,ux,uxx,uy,uyy
   integer, allocatable, save :: gl(:,:,:)
   integer, allocatable, save :: lci(:),lcj(:),lcvar(:)
endmodule space

module types
   double precision, save :: lenx(2),leny(2),const,xedge(2),yedge(2), total_energy
   double precision, allocatable, save :: xnod(:,:),ynod(:,:)
   double precision, allocatable, save :: u(:,:),du(:,:),ddu(:,:)
   double precision, allocatable, save :: u1(:,:), du1(:,:), ddu1(:,:), &
                                          u2(:,:), du2(:,:), ddu2(:,:), &
                                          u3(:,:), du3(:,:), ddu3(:,:), &  
                                          a2(:,:), da2(:,:), dda2(:,:), &  
                                          a3(:,:), da3(:,:), dda3(:,:)
   double precision, allocatable, save :: h(:,:) , h_x(:,:), h_y(:,:)
   double precision, allocatable, save :: s1(:,:),s2(:,:),s1_x(:,:),s2_y(:,:), s12(:,:)
   double precision, allocatable, save :: phi(:,:)
   double precision, allocatable, save :: Ham(:,:)
   integer,          allocatable, save :: PML(:,:)
endmodule types

module ellipsis
   double precision, save :: a  !!= 50.d0   !10!
   double precision, save :: b  !!=  5.d0   !20!
   double precision, save :: c  = 3.d0   !0.5
   double precision, save :: x0 = -20.d0
   double precision, save :: y0 = 0.d0 
   double precision, save :: h0 = 5.d0
endmodule ellipsis
