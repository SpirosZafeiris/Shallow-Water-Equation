subroutine linspace(N, a, b, x) ! checked
    implicit none
    integer,          intent(in)  :: N
    double precision, intent(in)  :: a, b
    double precision, intent(out) :: x(N)
    integer                       :: i
    double precision              :: step, junk
    step = (b-a)/(N-1)
    junk = a-step
    do i=1,N
        junk = junk + step
        x(i) = junk
    enddo
endsubroutine linspace

subroutine hstack(N1,N2,N3,x1,x2,x3,y) ! checked
    implicit none
    integer,          intent(in)  :: N1,N2,N3
    double precision, intent(in)  ::x1(N1), x2(N2), x3(N3)
    double precision, intent(out) :: y(N1+N2+N3)

    y(1:N1)                 = x1
    y(N1+1:(N1+N2))         = x2
    y((N1+N2+1):(N1+N2+N3)) = x3
endsubroutine hstack

subroutine chint(a, b, x, mask)
   use Param
   implicit none
   double precision, intent(in) :: a, b, x
   integer, intent(out) :: mask
   if ((x-a).lt.0.d0.or.(x-b).gt.0.d0) then
      mask = 0
   else
      mask = 1
   endif
   if (abs(x-b).lt.limit) then
      mask = 0
   endif
endsubroutine chint

subroutine unique(Nc,p,kv, y)
    implicit none
    integer,          intent(in)  :: Nc, p
    double precision, intent(in)  :: kv(Nc+p+1)
    integer                       :: i, u
    double precision, intent(out) :: y(Nc-p+1)
    u = Nc - p + 1
    do i=1,Nc-p+1
       y(i) = kv(i+p)
    enddo
endsubroutine


!subroutine belongsort(Nc,degree,Ncalc,xcalc,x, y)
!    use Param
!    integer,          intent(in)  :: Nc,degree,Ncalc
!    double precision, intent(in)  :: xcalc(Ncalc), x(Nc+degree+1)
!    integer,          intent(out) :: y(Ncalc)
!    integer                       :: i, N
!    prev = p
!    do i=1,Ncalc
!       do while (mask.ne.1)
!          prev = prev + 1
!          call chint(kv(prev), kv(prev+1), x(i), mask)
!          if (mask.eq.1) then
!             y(i) = prev
!          endif
!          if (prev.eq.(Nc+1)) then
!             exit
!          endif
!endsubroutine belongsort


subroutine belongkv(Nc,degree,Ncalc,xcalc,x, y)
    use Param
    implicit none
    external chint
    external unique
    integer,          intent(in)  :: Nc, degree, Ncalc
    double precision, intent(in)  :: xcalc(Ncalc), x(Nc+degree+1)
    integer,          intent(out) :: y(Ncalc)
    integer                       :: i,N,a,b,c,mask
    double precision              :: uniq(Nc-degree+1)
    N = Nc - degree + 1
    call unique(Nc,degree,x, uniq)
    do i=1,Ncalc
       a = 1
       b = N
       do while ((b-a).gt.1)
          c = int((a+b)/2)
          call chint(uniq(a),uniq(c),xcalc(i), mask)
          if (mask.eq.1) then
             b = c
          else
             a = c
          endif
       enddo
       y(i) = a
       if (.not.xcalc(i).lt.uniq(N)) then
          print *, "OK"
          y(i) = N-1
       endif
       y(i) = y(i) + degree
    enddo
endsubroutine belongkv   

          

!subroutine belongkv(Na, Nc, p, a, x, y)!x knot vector 
!    use Param
!    integer,          intent(in)  :: Na, Nc, p
!    double precision, intent(in)  :: a(Na), x(Nc+p+1)
!    integer,          intent(out) :: y(Na)
!    integer                       :: i, j
!    !print *, Na, Nc, p
!    !print *, a
!    !print *, Nc
!    !print *, x
!    do i=1,Na
!       do j=p+1,Nc
!          !print *, j
!          !print *, x(j)
!          if (a(i).ge.x(j).and.a(i).lt.x(j+1)) then
!             y(i) = j
!          endif
!          if (abs(a(i)-x(Nc+1)).le.limit) then
!             y(i) = Nc
!          endif
!       enddo
!       !print *, y(i)
!    enddo
!end subroutine belongkv

!subroutine unique(N, x, y, yindex, length)
!    use Param
!    integer,          intent(in)  :: N
!    double precision, intent(in)  :: x(N)
!    double precision, intent(out) :: y(N)
!    integer,          intent(out) :: length, yindex(N)
!    integer                       :: i, coun
!
!    call qsort(x, N)
!    coun = 0
!    call zeros(N, y)
!    call zerosi(N, yindex)
!    do i=2,N
!       if (abs(x(i) - x(i-1)).le.limit) then
!          coun = coun + 1
!       else
!          y(coun) = x(i)
!          yindex(coun) = i
!       endif
!    enddo
!    length = N - coun
!end subroutine unique

! -----------------------------------------------
! ----------------USEFUL FUNCTIONS---------------
! minloc(array, dim) (location of minimum)
! maxloc(array, dim) (location of maximum)
! findloc(array, value, dim) (location of value)
! minval(array, dim) (minimum)
! maxval(array, dim) (maximum)
!
!
!
! -----------------------------------------------
