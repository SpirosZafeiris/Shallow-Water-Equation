module Param
   use, intrinsic :: iso_fortran_env, only : real64
   real(real64), parameter :: pi = acos(-1._real64), limit = 1e-20_real64
   real(real64), parameter :: rho = 1024._real64
end module Param

module time
   use, intrinsic :: iso_fortran_env, only : real64
   real(real64), save :: dt,t,beta,gamma
   !real(real64), allocatable, save :: rhs(:)
   !integer, allocatable :: line(:),column(:), num_col(:), num_col_acc(:)
   real(real64), allocatable :: val(:)
   integer :: dof2,dof6
endmodule time

module sizes
   integer, save :: Nx,Ny,Nxnod,Nynod,dof,Npml
endmodule sizes

module space
   use, intrinsic :: iso_fortran_env, only : real64
   real(real64), save :: dx,dy,dx2,dx22,dy2,dy22,dt2,dt22,ux,uxx,uy,uyy
   integer, allocatable, save :: gl(:,:,:)
   integer, allocatable, save :: lci(:),lcj(:),lcvar(:)
endmodule space

module types
   use, intrinsic :: iso_fortran_env, only : real64
   real(real64), save :: lenx(2),leny(2),const,xedge(2),yedge(2), total_energy
   real(real64), allocatable, save :: xnod(:,:),ynod(:,:)
   real(real64), allocatable, save :: u(:,:),du(:,:),ddu(:,:)
   real(real64), allocatable, save :: u1(:,:), du1(:,:), ddu1(:,:), &
                                          u2(:,:), du2(:,:), ddu2(:,:), &
                                          u3(:,:), du3(:,:), ddu3(:,:), &  
                                          a2(:,:), da2(:,:), dda2(:,:), &  
                                          a3(:,:), da3(:,:), dda3(:,:)
   real(real64), allocatable, save :: h(:,:) , h_x(:,:), h_y(:,:)
   real(real64), allocatable, save :: s1(:,:),s2(:,:),s1_x(:,:),s2_y(:,:), s12(:,:)
   real(real64), allocatable, save :: phi(:,:)
   real(real64), allocatable, save :: Ham(:,:)
   integer,          allocatable, save :: PML(:,:)
endmodule types

module ellipsis
   use, intrinsic :: iso_fortran_env, only : real64
   real(real64), save :: a  !!= 50._real64   !10!
   real(real64), save :: b  !!=  5._real64   !20!
   real(real64), save :: c  = 3._real64   !0.5
   real(real64), save :: x0 = -20._real64
   real(real64), save :: y0 = 0._real64 
   real(real64), save :: h0 = 5._real64
endmodule ellipsis
