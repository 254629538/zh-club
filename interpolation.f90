module interpolation
use blas95
use f95_precision
use lapack95
contains
subroutine solve(n,x,y,da,db,t,ty)
!--------------------------------------
!input parameters:
!n:  Number of interpolation nodes -1
!x:  argument   (:n)dimensional vector
!y:
!nt: dimensional
!t:  vectors to be calculated  (1:nt)
!da:  f'(x(0))
!db:  f'(x(n))
!output parameter:
!ty: vectors to be calculated  (:nt)
!--------------------------------------
implicit none
integer::n,nt
integer::i,j,k
real*8::da,db
real*8::x(0:n),y(0:n),t,ty,h(0:n-1),f1(0:n-1),f2(1:n-1)
real*8::u(1:n-1),namda(1:n-1),d(0:n),M(0:n),A(0:n,0:n)
do i=0,n-1
  h(i)=x(i+1)-x(i)
  f1(i)=(y(i+1)-y(i))/h(i)
end do
!set the boundary conditions
d(0)=6d0/h(0)*(f1(0)-da)
d(n)=6d0/h(n-1)*(db-f1(n-1))
!get u, namda, d
do i=1,n-1
  u(i)=h(i-1)/(h(i-1)+h(i))
  namda(i)=1-u(i)
  f2(i)=(f1(i-1)-f1(i))/(x(i-1)-x(i+1))
  d(i)=6d0*f2(i)
end do
!set matrix_A
A=0
do i=1,n-1
  A(i,i)=2D0
end do
do i=2,n-1
 A(i,i-1)=u(i)
end do
do i=1,n-2
 A(i,i+1)=namda(i)
end do
A(0,0)=2D0;A(0,1)=1D0
A(n,n-1)=1D0;A(n,n)=2D0
!set right vector
d(1)=d(1)-u(1)*M(0)
d(n-1)=d(n-1)-namda(n-1)*M(n)
!get M(n) and complete the establishment of interpolation polynomial
call gesvx(A,d,M)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do k=1,nt
! find the location in the data for each component of the interpolation vector
 do i=1,n-1
  if (t<x(i+1)) exit
 end do
  ty=M(i)*(x(i+1)-t)**3D0/(6D0*h(i))+M(i+1)*(t-x(i))**3D0/(6D0*h(i))+(y(i)-M(i)*h(i)**2D0/6D0)*(x(i+1)-t)/h(i)+(y(i+1)-M(i+1)*h(i)**2d0/6D0)*(t-x(i))/h(i)
!end do
end subroutine solve
end module interpolation


! program main
! use interpolation
! implicit none
! integer::i,j
! real*8::x(0:8),y(0:8),t(5),ty(5),simty(5),da,db
! x=(/-2d0,-1.5D0,-1D0,-0.5d0,0d0,0.5d0,1d0,1.5d0,2.0d0/)
! y=(/ -0.9093d0,-0.9975d0,-0.8415d0,-0.4794d0,0d0,0.4794d0,0.8415d0,0.9975d0,0.9093d0/)
! da=0d0 !-0.4161d0
! db=0d0 !-0.4161d0
! t=(/0.4,-0.6,1.7,0.8,1.8/)
! !x=(/1d0,2d0,5d0,6d0,7d0,8d0,10d0,13d0,17d0/)
! !y=(/3d0,3.7d0,3.9d0,4.2d0,5.7d0,6.6d0,7.1d0,6.7d0,4.5d0/)
! !da=0d0;db=0d0
! !t=(/3.5d0,1.5d0,5.5d0,6.5d0,7.5d0,9.0d0,11.5d0,15d0/)
! call solve(8,x,y,da,db,5,t,ty)
! simty=dsin(t)
! write(*,103)((i,t(i),ty(i),simty(i)),i=1,5)
! 103 format (I3,3F16.8)
! end program main




!ifort Source1.f90 -I/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/include/intel64/lp64 -L/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -liomp5
