!subroutine wave_source(x,y,t, yy)
!   use Param
!   implicit none
!   double precision, intent(in)  :: x,y,t
!   double precision, intent(out) :: yy
!   double precision              :: A0,om
!   A0 = 0.01
!   om = 0.4*pi
!   !yy = exp(-50*tempx**2)*(0.5+0.5*tanh(50*(tempy+0.5)))*(0.5-0.5*tanh(50*(tempy-0.5)))
!   !yy = A0*yy*cos(om*t)
!   yy = exp(-0.01* (x+100)**2)
!   !yy = 0.05*(1-(x+100)**2)*exp(-0.01* (x+100)**2)
!   yy = yy*(0.5+0.5*tanh(0.50*(y+120)))*(0.5-0.5*tanh(0.50*(y-120)))
!   yy = A0*yy*cos(om*t)
!   !yy = 0.d0
!endsubroutine wave_source


subroutine slit(x,y, yy)
   use Param
   implicit none
   double precision, intent(in) :: x,y
   double precision, intent(out) :: yy
   yy  = 4*exp(-0.02*x **2)*sin(2* pi/60*y )*(0.5+0.5*tanh(0.50*(y +275)))*(0.5-0.5* tanh(0.50*(y -275)))
endsubroutine slit

