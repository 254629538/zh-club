module vec
  integer,parameter::n_=300
  real*8 k(n_),tau(n_),fome(n_),bome(n_)
  save k,tau,fome,bome
  parameter (pi=3.141592653589793238462643383279502884197D0)
end module vec

module formidconstant
real*8 t,eta,omega0,lambda,ntot
save t,eta,omega0,lambda,ntot
end module formidconstant

subroutine vector(xx,m,beta)
!give 4 different vector for SFT, used as x_i in f(x_i)
!n give the dimention of the vector. 40 <= n <= 300
!m is the kind of vector
!1 named 0~beta,2 named  sum (2n+1)Pi/beta, 3 named sum (2n)Pi/beta
!4 named int infinity
use vec
implicit none
real*8,intent(out):: xx(n_)
real*8,intent(in):: beta
integer,intent(in):: m
integer j,n
select case(m)
  !   tau from 0-beta in logarithmic scales
case(1)
  n=n_/2
do j = 1,n
   xx(j) =sinh(0.02D0*dble(j))/sinh(0.02D0*dble(n))*(beta/2.0D0-0.03D0)+1.0D0/beta/150.0D0!dble(j)/dble(n)*exp(-dble(n-j)/dble(n))*(beta/2.0D0-0.02D0)
  !xx(n+1-j)=(1.0D0-exp(-dble(n-j)/dble(n)*5.0D0))*beta/2.0D0+0.01D0
enddo
do j=n+1,n_
  xx(j) = beta - xx(n_+1-j)
end do

! do j=1,n_
!   xx(j)=(beta-0.01D0)*dble(j)/dble(n_)
! end do
! do j=1,n_
!   xx(j)=(beta-0.02D0)*exp(-dble(n_-j)/dble(n_)*5.5D0)
! end do
  !   omegaF  (2n+1)pi/beta  n: e^0~e^8
case(2)
  n=20!n_/2
  do j = 1,n
    xx(j)=dble(2*j-1)*pi/beta
  enddo
  ! !xx(n+1)=dble(2*(n+3)-1)*pi/beta
  ! do j = n+1,n_
  !   xx(j) = (2.0D0*dble(dint(exp(dble(j-n)/dble(n_-n)*5.0D0+4.5D0)-exp(4.5D0))+j)-1.0D0)*pi/beta
  ! enddo
  do j = n+1,n_
    xx(j)=(2.0D0*dble(dint(sinh(0.020D0*dble(j))/0.020D0)+j-n)-1.0D0)*pi/beta
  end do
! do j = n+1,n_
!   xx(j) = (2.0D0*(exp(dble(j-n)/dble(n_-n)*8.0D0)+dble(j)-1.0D0)-1.0D0)*pi/beta
! enddo
! do j = 1,n_
!   xx(j) = (2.0D0*(exp(dble(j-1)/dble(n_-1)*12.0D0)+dble(j)-1.0D0)-1.0D0)*pi/beta
! enddo
  !   omegaB  (2n)pi/beta  n: e^0~e^8
case(3)
  n=20
  do j = 1,n
    xx(j)=dble(2*j)*pi/beta
  enddo
  ! !xx(n+1)=dble(2*(n+3))*pi/beta
  ! do j = n+1,n_-1
  !   xx(j) = 2.0D0*dble(dint(exp(dble(j-n)/dble(n_-n)*5.0D0+4.5D0)-exp(4.5D0))+j)*pi/beta
  ! enddo
  ! xx(n_)=0D0
  do j = n+1,n_-1
    xx(j)=2.0D0*dble(dint(sinh(0.020D0*dble(j))/0.020D0)+j-n)*pi/beta
  end do
  xx(n_)=0D0
! do j = n+1,n_
!   xx(j) = 2.0D0*(exp(dble(j-n)/dble(n_-n)*8.0D0)+dble(j)-1.0D0)*pi/beta
! enddo

  !   k  0-e^6 in logarithmic scales
case(4)
do j = 1,n_
  !xx(j) = dble(j)/dble(n_)*exp(dble(j)/dble(n_)*7.0D0)!-1.0D0
   xx(j) = exp(dble(j)/dble(n_)*6.5D0)-1.02D0                           !         !!! n=500::  -1.01D0 ； n=300:  -1.02D0
enddo

end select
end subroutine vector

SUBROUTINE spline(x,y,n,y2)
  !give an array y2(1:n) of length n which contains the second derivatives of the interpolating function at the tabulated points x
INTEGER n,NMAX
REAL*8 x(n),y(n),y2(n)
PARAMETER (NMAX=2000)
INTEGER i,k
REAL*8 p,qn,sig,un,u(NMAX)

y2(1)=0.0D0
u(1)=0.0D0

do i=2,n-1
sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
p=sig*y2(i-1)+2.0D0
y2(i)=(sig-1.0D0)/p
u(i)=(6.0D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p

enddo

qn=0.0D0
un=0.0D0

y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0D0)
do  k=n-1,1,-1
y2(k)=y2(k)*y2(k+1)+u(k)
enddo
return
end subroutine


subroutine SFT_ktox(inmatrix,outmatrix)
!SFT_ktox transfer a function from k-space to r-space
!E^(ikr)
!inmatrix is i*j function value,i->k
!outmatrix is m*n function value,m->r
use vec
use omp_lib
implicit none
real*8,intent(in):: inmatrix(n_,n_,2)
real*8,intent(out):: outmatrix(n_,n_,2)
real*8::x(n_),y1(2*n_),y2(2*n_),y11(2*n_),y22(2*n_),xx(2*n_)
real*8::mid(4,2),kkk
integer i1,i2,i3,ncho,ncho1
complex*16::coe1,coe2,In0(n_),In1(n_),In2(n_),In3(n_),a(n_),b(n_),c(n_),d(n_),f1(n_)
complex*16::co1,co2,co3,co4,step
! i1 index k, i2 index tau ,i3 index x
 x=k
! ! do i3=1,n_
! !   ncho1=i3-1
! !   if (x(i3)>1.0D-1) exit
! ! end do
do i3=1,n_
  ncho=i3-1
  if (k(i3)>1.0D0/k(1)) exit
end do
! ncho=n_

xx(1:n_)=-k(n_:1:-1)
xx(n_+1:2*n_)=k(1:n_:1)
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,coe1,coe2,In0,In1,In2,In3,a,b,c,d,f1,co1,co2,co3,co4,y1,y2,y11,y22,mid,kkk,step) SCHEDULE(DYNAMIC,1)
do i2=1,n_
  a=0;b=0;c=0;d=0;In0=0;In1=0;In2=0;In3=0
  do i1=1,n_
    y1(n_-i1+1)=inmatrix(i1,i2,1)*k(i1)     !  y1: real part
    y1(n_+i1)=inmatrix(i1,i2,1)*k(i1)
    y2(n_-i1+1)=inmatrix(i1,i2,2)*k(i1)     !  y2: imaginary part
    y2(n_+i1)=inmatrix(i1,i2,2)*k(i1)
  end do

call spline(xx,y1,2*n_,y11)            ! y11: y1''
call spline(xx,y2,2*n_,y22)            ! y22: y2''

do i1=n_+1,2*n_-1
  mid(1,1)=y1(i1)
  mid(1,2)=y2(i1)
  mid(2,1)=(y1(i1+1)-y1(i1))/(xx(i1+1)-xx(i1))-(xx(i1+1)-xx(i1))*(2.0D0*y11(i1)+y11(i1+1))/6.0D0
  mid(2,2)=(y2(i1+1)-y2(i1))/(xx(i1+1)-xx(i1))-(xx(i1+1)-xx(i1))*(2.0D0*y22(i1)+y22(i1+1))/6.0D0
  mid(3,1)=y11(i1)/2.0D0
  mid(3,2)=y22(i1)/2.0D0
  mid(4,1)=(y11(i1+1)-y11(i1))/(xx(i1+1)-xx(i1))/6.0D0
  mid(4,2)=(y22(i1+1)-y22(i1))/(xx(i1+1)-xx(i1))/6.0D0

  a(i1-n_)=cmplx(mid(1,1),0D0)!,mid(1,2))                  !!!!!!!!!!
  b(i1-n_)=cmplx(mid(2,1),0D0)!mid(2,2))                  !!!!!!!!!!
  c(i1-n_)=cmplx(mid(3,1),0D0)!mid(3,2))                  !!!!!!!!!!
  d(i1-n_)=cmplx(mid(4,1),0D0)!mid(4,2))                  !!!!!!!!!!
end do

! do i3=1,ncho1
!   do i1=1,n_-1
!     coe1=exp(cmplx(0.0D0,x(i3)*k(i1)))
!     coe2=exp(cmplx(0.0D0,x(i3)*(k(i1+1)-k(i1))))
!       In0(i1)=coe1*coe2*cmplx(k(i1+1)-k(i1),0.0D0)
!       In1(i1)=coe1*coe2*cmplx((k(i1+1)-k(i1))**2,0.0D0)
!       In2(i1)=coe1*coe2*cmplx((k(i1+1)-k(i1))**3,0.0D0)
!       In3(i1)=coe1*coe2*cmplx((k(i1+1)-k(i1))**4,0.0D0)
!   end do
!   f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)
!
!   outmatrix(i3,i2,1)=2.0D0*aimag(f1(i3))/x(i3)
!   outmatrix(i3,i2,2)=0.0D0
! end do

do i3=1,ncho
  f1(i3)=0
  do i1=1,n_-1
    coe1=exp(cmplx(0.0D0,x(i3)*k(i1)))
    coe2=exp(cmplx(0.0D0,x(i3)*(k(i1+1)-k(i1))))
   if ((abs(x(i3)*(k(i1+1)-k(i1))).gt.1.0D-1)) then
      In0(i1)=coe1*(coe2-cmplx(1.0D0,0.0D0))*cmplx(0.0D0,-1.0D0/x(i3))
      In1(i1)=coe1*((coe2-cmplx(1.0D0,0.0D0))*cmplx(1.0D0/x(i3)**2,0.0D0)-cmplx(0.0D0,1.0D0)*coe2*cmplx((k(i1+1)-k(i1))/x(i3),0.0D0))
      In2(i1)=coe1*(cmplx(2.0D0*(k(i1+1)-k(i1))/x(i3)**2,0.0D0)*coe2+cmplx(0.0D0,1.0D0)*(cmplx(2.0D0/x(i3)**3,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))-coe2*cmplx((k(i1+1)-k(i1))**2/x(i3),0.0D0)))
      In3(i1)=coe1*(-cmplx(6.0D0/x(i3)**4,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))+cmplx(3.0D0*(k(i1+1)-k(i1))**2/x(i3)**2,0.0D0)*coe2+cmplx(0.0D0,1.0D0)*coe2*(cmplx(6.0D0*(k(i1+1)-k(i1))/x(i3)**3-(k(i1+1)-k(i1))**3/x(i3),0.0D0)))
    else
      kkk=k(i1+1)-k(i1)
      In0(i1)=coe1*cmplx(kkk-kkk**3*x(i3)**2/6.0D0+kkk**5*x(i3)**4/120.0D0-kkk**7*x(i3)**6/5040.0D0,kkk**2*x(i3)/2.0D0-kkk**4*x(i3)**3/24.0D0+kkk**6*x(i3)**5/720.0D0)
      In1(i1)=coe1*cmplx(kkk**2/2.0D0-kkk**4*x(i3)**2/8.0D0+kkk**6*x(i3)**4/144.0D0-kkk**8*x(i3)**6/5760.0D0,kkk**3*x(i3)/3.0D0-kkk**5*x(i3)**3/30.0D0+kkk**7*x(i3)**5/840.0D0)
      In2(i1)=coe1*cmplx(kkk**3/3.0D0-kkk**5*x(i3)**2/10.0D0+kkk**7*x(i3)**4/168.0D0-kkk**9*x(i3)**6/6480.0D0,kkk**4*x(i3)/4.0D0-kkk**6*x(i3)**3/36.0D0+kkk**8*x(i3)**5/960.0D0)
      In3(i1)=coe1*cmplx(kkk**4/4.0D0-kkk**6*x(i3)**2/12.0D0+kkk**8*x(i3)**4/192.0D0-kkk**10*x(i3)**6/7200.0D0+kkk**12*x(i3)**8/483840.0D0,kkk**5*x(i3)/5.0D0-kkk**7*x(i3)**3/42.0D0+kkk**9*x(i3)**5/1080.0D0-kkk**11*x(i3)**7/55440.0D0)
    end if
    step=In0(i1)*a(i1)+In1(i1)*b(i1)+In2(i1)*c(i1)+In3(i1)*d(i1)
    f1(i3)=f1(i3)+step
  end do
  ! f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)

  outmatrix(i3,i2,1)=2.0D0*aimag(f1(i3))/x(i3)
  outmatrix(i3,i2,2)=0.0D0
end do

do i3=ncho+1,n_
  co1 = cmplx(0.0D0,-1.0D0/x(i3))
  co2 = cmplx(1.0D0/x(i3)**2,0.0D0)
  co3 = cmplx(0.0D0,2.0D0/x(i3)**3)
  co4 = cmplx(-6.0D0/x(i3)**4,0.0D0)
  step=cmplx(0.0D0,0.0D0)
  do i1=1,n_-1
    step = step+d(i1)*(exp(cmplx(0.0D0,x(i3)*k(i1+1)))-exp(cmplx(0.0D0,x(i3)*k(i1))))
  end do
  coe1=exp(cmplx(0.0D0,x(i3)*k(n_)))
  coe2=exp(cmplx(0.0D0,x(i3)*k(1)))
  f1(i3)=co1*(coe1*a(n_-1)-coe2*a(1))+co2*(coe1*b(n_-1)-coe2*b(1))+co3*(coe1*c(n_-1)-coe2*c(1))+co4*step

  outmatrix(i3,i2,1)=2.0D0*aimag(f1(i3))/x(i3)
  outmatrix(i3,i2,2)=0.0D0
end do
end do
!$OMP END PARALLEL DO
outmatrix=outmatrix*0.025330295911D0    !  1/(2pi)^2
end subroutine

subroutine SFT_xtok(inmatrix,outmatrix)
use vec
implicit none
real*8,intent(in):: inmatrix(n_,n_,2)
real*8,intent(out):: outmatrix(n_,n_,2)
call SFT_ktox(inmatrix,outmatrix)
outmatrix=outmatrix*248.0502134D0    ! (2pi)^3
end subroutine SFT_xtok


subroutine SFT_tautoomegaF(inmatrix,outmatrix,beta)
!SFT_tautoomegaF transfer a function from tau-space to omega_n-space, the second F means fermion
!E^(i*omega_n*tau)
!inmatrix is i*j function value,j->tau
!outmatrix is m*n function value,n->omega
use vec
use omp_lib
implicit none
real*8,intent(in):: inmatrix(n_,n_,2),beta
real*8,intent(out):: outmatrix(n_,n_,2)
real*8::x(n_),y1(n_),y2(n_),y11(n_),y22(n_)
real*8::mid(4,2),kkk
integer i1,i2,i3,ncho
complex*16::coe1,coe2,In0(n_),In1(n_),In2(n_),In3(n_),a(n_),b(n_),c(n_),d(n_),f1(n_),f2(n_)
complex*16::co1,co2,co3,co4,step
! i1 index k, i2 index tau ,i3 index fomega

do i3=1,n_
  ncho=i3-1
  if (fome(i3)>1.0D0/tau(1)) exit
end do
ncho=n_
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,coe1,coe2,In0,In1,In2,In3,a,b,c,d,f1,co1,co2,co3,co4,y1,y2,y11,y22,mid,kkk,step) SCHEDULE(DYNAMIC,1)
do i1=1,n_
  a=0;b=0;c=0;d=0;In0=0;In1=0;In2=0;In3=0
  do i2=1,n_
    y1(i2)=inmatrix(i1,i2,1)     !  y1: real part
    y2(i2)=inmatrix(i1,i2,2)     !  y2: imaginary part
  end do

call spline(tau,y1,n_,y11)            ! y11: y1''
call spline(tau,y2,n_,y22)            ! y22: y2''

do i2=1,n_-1
  mid(1,1)=y1(i2)
  mid(1,2)=y2(i2)
  mid(2,1)=(y1(i2+1)-y1(i2))/(tau(i2+1)-tau(i2))-(tau(i2+1)-tau(i2))*(2.0D0*y11(i2)+y11(i2+1))/6.0D0
  mid(2,2)=(y2(i2+1)-y2(i2))/(tau(i2+1)-tau(i2))-(tau(i2+1)-tau(i2))*(2.0D0*y22(i2)+y22(i2+1))/6.0D0
  mid(3,1)=y11(i2)/2.0D0
  mid(3,2)=y22(i2)/2.0D0
  mid(4,1)=(y11(i2+1)-y11(i2))/(tau(i2+1)-tau(i2))/6.0D0
  mid(4,2)=(y22(i2+1)-y22(i2))/(tau(i2+1)-tau(i2))/6.0D0

  a(i2)=cmplx(mid(1,1),mid(1,2))
  b(i2)=cmplx(mid(2,1),mid(2,2))
  c(i2)=cmplx(mid(3,1),mid(3,2))
  d(i2)=cmplx(mid(4,1),mid(4,2))
end do

do i3=1,ncho
  f1(i3)=0
  do i2=1,n_-1
    coe1=exp(cmplx(0.0D0,fome(i3)*tau(i2)))
    coe2=exp(cmplx(0.0D0,fome(i3)*(tau(i2+1)-tau(i2))))
    if ((abs(fome(i3)*(tau(i2+1)-tau(i2))).gt.1.0D-1)) then
      In0(i2)=coe1*(coe2-cmplx(1.0D0,0.0D0))*cmplx(0.0D0,-1.0D0/fome(i3))
      In1(i2)=coe1*((coe2-cmplx(1.0D0,0.0D0))*cmplx(1.0D0/fome(i3)**2,0.0D0)-cmplx(0.0D0,1.0D0)*coe2*cmplx((tau(i2+1)-tau(i2))/fome(i3),0.0D0))
      In2(i2)=coe1*(cmplx(2.0D0*(tau(i2+1)-tau(i2))/fome(i3)**2,0.0D0)*coe2+cmplx(0.0D0,1.0D0)*(cmplx(2.0D0/fome(i3)**3,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))-coe2*cmplx((tau(i2+1)-tau(i2))**2/fome(i3),0.0D0)))
      In3(i2)=coe1*(-cmplx(6.0D0/fome(i3)**4,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))+cmplx(3.0D0*(tau(i2+1)-tau(i2))**2/fome(i3)**2,0.0D0)*coe2+cmplx(0.0D0,1.0D0)*coe2*(cmplx(6.0D0*(tau(i2+1)-tau(i2))/fome(i3)**3-(tau(i2+1)-tau(i2))**3/fome(i3),0.0D0)))
    else
      kkk=tau(i2+1)-tau(i2)
      In0(i2)=coe1*cmplx(kkk-kkk**3*fome(i3)**2/6.0D0+kkk**5*fome(i3)**4/120.0D0-kkk**7*fome(i3)**6/5040.0D0,kkk**2*fome(i3)/2.0D0-kkk**4*fome(i3)**3/24.0D0+kkk**6*fome(i3)**5/720.0D0-fome(i3)**7*kkk**8/40320.0D0)
      In1(i2)=coe1*cmplx(kkk**2/2.0D0-kkk**4*fome(i3)**2/8.0D0+kkk**6*fome(i3)**4/144.0D0-kkk**8*fome(i3)**6/5760.0D0,kkk**3*fome(i3)/3.0D0-kkk**5*fome(i3)**3/30.0D0+kkk**7*fome(i3)**5/840.0D0)
      In2(i2)=coe1*cmplx(kkk**3/3.0D0-kkk**5*fome(i3)**2/10.0D0+kkk**7*fome(i3)**4/168.0D0-kkk**9*fome(i3)**6/6480.0D0,kkk**4*fome(i3)/4.0D0-kkk**6*fome(i3)**3/36.0D0+kkk**8*fome(i3)**5/960.0D0)
      In3(i2)=coe1*cmplx(kkk**4/4.0D0-kkk**6*fome(i3)**2/12.0D0+kkk**8*fome(i3)**4/192.0D0-kkk**10*fome(i3)**6/7200.0D0+kkk**12*fome(i3)**8/483840.0D0,kkk**5*fome(i3)/5.0D0-kkk**7*fome(i3)**3/42.0D0+kkk**9*fome(i3)**5/1080.0D0-kkk**11*fome(i3)**7/55440.0D0)
    end if
    step=In0(i2)*a(i2)+In1(i2)*b(i2)+In2(i2)*c(i2)+In3(i2)*d(i2)
    f1(i3)=f1(i3)+step
  end do
  ! f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)

  outmatrix(i1,i3,1)=real(f1(i3))
  outmatrix(i1,i3,2)=aimag(f1(i3))
end do

do i3=ncho+1,n_
  co1 = cmplx(0.0D0,-1.0D0/fome(i3))
  co2 = cmplx(1.0D0/fome(i3)**2,0.0D0)
  co3 = cmplx(0.0D0,2.0D0/fome(i3)**3)
  co4 = cmplx(-6.0D0/fome(i3)**4,0.0D0)
  step=cmplx(0.0D0,0.0D0)
  do i2=1,n_-1
    step = step+d(i2)*(exp(cmplx(0.0D0,fome(i3)*tau(i2+1)))-exp(cmplx(0.0D0,fome(i3)*tau(i2))))
  end do
  coe1=exp(cmplx(0.0D0,fome(i3)*tau(n_-1)))
  coe2=exp(cmplx(0.0D0,fome(i3)*tau(1)))
  f1(i3)=co1*(coe1*a(n_-1)-coe2*a(1))+co2*(coe1*b(n_-1)-coe2*b(1))+co3*(coe1*c(n_-1)-coe2*c(1))+co4*step

  outmatrix(i1,i3,1)=real(f1(i3))
  outmatrix(i1,i3,2)=aimag(f1(i3))
end do
end do
!$OMP END PARALLEL DO
end subroutine SFT_tautoomegaF


subroutine SFT_tautoomegaB(inmatrix,outmatrix,beta)

use vec
implicit none
real*8,intent(in):: inmatrix(n_,n_,2),beta
real*8,intent(out):: outmatrix(n_,n_,2)
real*8::x(n_),y1(n_),y2(n_),y11(n_),y22(n_)
real*8::mid(4,2),kkk
integer i1,i2,i3,ncho
complex*16::coe1,coe2,In0(n_),In1(n_),In2(n_),In3(n_),a(n_),b(n_),c(n_),d(n_),f1(n_),f2(n_)
complex*16::co1,co2,co3,co4,step
! i1 index k, i2 index tau ,i3 index bomega

do i3=1,n_
  ncho=i3-1
  if (bome(i3)>1.0D0/tau(1)) exit
end do
ncho=n_-1
do i1=1,n_
  a=0;b=0;c=0;d=0;In0=0;In1=0;In2=0;In3=0
  do i2=1,n_
    y1(i2)=inmatrix(i1,i2,1)     !  y1: real part
    y2(i2)=inmatrix(i1,i2,2)     !  y2: imaginary part
  end do

call spline(tau,y1,n_,y11)            ! y11: y1''
call spline(tau,y2,n_,y22)            ! y22: y2''

do i2=1,n_-1
  mid(1,1)=y1(i2)
  mid(1,2)=y2(i2)
  mid(2,1)=(y1(i2+1)-y1(i2))/(tau(i2+1)-tau(i2))-(tau(i2+1)-tau(i2))*(2.0D0*y11(i2)+y11(i2+1))/6.0D0
  mid(2,2)=(y2(i2+1)-y2(i2))/(tau(i2+1)-tau(i2))-(tau(i2+1)-tau(i2))*(2.0D0*y22(i2)+y22(i2+1))/6.0D0
  mid(3,1)=y11(i2)/2.0D0
  mid(3,2)=y22(i2)/2.0D0
  mid(4,1)=(y11(i2+1)-y11(i2))/(tau(i2+1)-tau(i2))/6.0D0
  mid(4,2)=(y22(i2+1)-y22(i2))/(tau(i2+1)-tau(i2))/6.0D0

  a(i2)=cmplx(mid(1,1),mid(1,2))
  b(i2)=cmplx(mid(2,1),mid(2,2))
  c(i2)=cmplx(mid(3,1),mid(3,2))
  d(i2)=cmplx(mid(4,1),mid(4,2))
end do

do i3=1,ncho
  f1(i3)=0
  do i2=1,n_-1
    coe1=exp(cmplx(0.0D0,bome(i3)*tau(i2)))
    coe2=exp(cmplx(0.0D0,bome(i3)*(tau(i2+1)-tau(i2))))
    if (abs(bome(i3)*(tau(i2+1)-tau(i2))) .gt. 1.0D-1) then
      In0(i2)=coe1*(coe2-cmplx(1.0D0,0.0D0))*cmplx(0.0D0,-1.0D0/bome(i3))
      In1(i2)=coe1*((coe2-cmplx(1.0D0,0.0D0))*cmplx(1.0D0/bome(i3)**2,0.0D0)-cmplx(0.0D0,1.0D0)*coe2*cmplx((tau(i2+1)-tau(i2))/bome(i3),0.0D0))
      In2(i2)=coe1*(cmplx(2.0D0*(tau(i2+1)-tau(i2))/bome(i3)**2,0.0D0)*coe2+cmplx(0.0D0,1.0D0)*(cmplx(2.0D0/bome(i3)**3,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))-coe2*cmplx((tau(i2+1)-tau(i2))**2/bome(i3),0.0D0)))
      In3(i2)=coe1*(-cmplx(6.0D0/bome(i3)**4,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))+cmplx(3.0D0*(tau(i2+1)-tau(i2))**2/bome(i3)**2,0.0D0)*coe2+cmplx(0.0D0,1.0D0)*coe2*(cmplx(6.0D0*(tau(i2+1)-tau(i2))/bome(i3)**3-(tau(i2+1)-tau(i2))**3/bome(i3),0.0D0)))
    else
      kkk=tau(i2+1)-tau(i2)
      In0(i2)=coe1*cmplx(kkk-kkk**3*bome(i3)**2/6.0D0+kkk**5*bome(i3)**4/120.0D0-kkk**7*bome(i3)**6/5040.0D0,kkk**2*bome(i3)/2.0D0-kkk**4*bome(i3)**3/24.0D0+kkk**6*bome(i3)**5/720.0D0-bome(i3)**7*kkk**8/40320.0D0)
      In1(i2)=coe1*cmplx(kkk**2/2.0D0-kkk**4*bome(i3)**2/8.0D0+kkk**6*bome(i3)**4/144.0D0-kkk**8*bome(i3)**6/5760.0D0,kkk**3*bome(i3)/3.0D0-kkk**5*bome(i3)**3/30.0D0+kkk**7*bome(i3)**5/840.0D0)
      In2(i2)=coe1*cmplx(kkk**3/3.0D0-kkk**5*bome(i3)**2/10.0D0+kkk**7*bome(i3)**4/168.0D0-kkk**9*bome(i3)**6/6480.0D0,kkk**4*bome(i3)/4.0D0-kkk**6*bome(i3)**3/36.0D0+kkk**8*bome(i3)**5/960.0D0)
      In3(i2)=coe1*cmplx(kkk**4/4.0D0-kkk**6*bome(i3)**2/12.0D0+kkk**8*bome(i3)**4/192.0D0-kkk**10*bome(i3)**6/7200.0D0+kkk**12*bome(i3)**8/483840.0D0,kkk**5*bome(i3)/5.0D0-kkk**7*bome(i3)**3/42.0D0+kkk**9*bome(i3)**5/1080.0D0-kkk**11*bome(i3)**7/55440.0D0)
    end if
    step=In0(i2)*a(i2)+In1(i2)*b(i2)+In2(i2)*c(i2)+In3(i2)*d(i2)
    f1(i3)=f1(i3)+step
  end do
  ! f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)

  outmatrix(i1,i3,1)=real(f1(i3))
  outmatrix(i1,i3,2)=aimag(f1(i3))
end do

do i3=ncho+1,n_-1
  co1 = cmplx(0.0D0,-1.0D0/bome(i3))
  co2 = cmplx(1.0D0/bome(i3)**2,0.0D0)
  co3 = cmplx(0.0D0,2.0D0/bome(i3)**3)
  co4 = cmplx(-6.0D0/bome(i3)**4,0.0D0)
  step=cmplx(0.0D0,0.0D0)
  do i2=1,n_-1
    step = step+d(i2)*(exp(cmplx(0.0D0,bome(i3)*tau(i2+1)))-exp(cmplx(0.0D0,bome(i3)*tau(i2))))
  end do
  coe1=exp(cmplx(0.0D0,bome(i3)*tau(n_)))
  coe2=exp(cmplx(0.0D0,bome(i3)*tau(1)))
  f1(i3)=co1*(coe1*a(n_-1)-coe2*a(1))+co2*(coe1*b(n_-1)-coe2*b(1))+co3*(coe1*c(n_-1)-coe2*c(1))+co4*step

  outmatrix(i1,i3,1)=real(f1(i3))
  outmatrix(i1,i3,2)=aimag(f1(i3))
end do
!--------------   bomega(n_)=0.0D0  -------------
step=cmplx(0.0D0,0.0D0)
do i2=1,n_-1
   coe2=cmplx(tau(i2+1)-tau(i2),0.0D0)
   step=step+a(i2)*coe2+b(i2)*coe2**2/cmplx(2.0D0,0.0D0)+c(i2)*coe2**3/cmplx(3.0D0,0.0D0)+d(i2)*coe2**4/cmplx(4.0D0,0.0D0)
 enddo
 outmatrix(i1,n_,1)=real(step)
 outmatrix(i1,n_,2)=aimag(step)
!--------------
end do
end subroutine SFT_tautoomegaB


subroutine SFT_omegaFtotau(inmatrix,outmatrix,beta)
!SFT_omegaFtotau transfer a function from omega_n-space to tau-space, the second F means Fermion
use vec
implicit none
real*8,intent(in):: inmatrix(n_,n_,2),beta
real*8,intent(out):: outmatrix(n_,n_,2)
real*8::y1(2*n_),y2(2*n_),y11(2*n_),y22(2*n_),xx(2*n_)
real*8::mid(4,2),dw,ttt,ooo
integer i1,i2,i3,ncho,ncho1,ncho0
complex*16::In0(n_),In1(n_),In2(n_),In3(n_),a(n_),b(n_),c(n_),d(n_),f1(n_)
complex*16::co1,co2,co3,co4,step,coe1,coe2,coe3,tan3
! i1 index k, i2 index fomega ,i3 index tau
! do i3=1,n_
!   ncho0=i3-1
!   if (tau(i3)>1.0D-1) exit
! end do
outmatrix = 0.0D0

do i3=1,n_
  ncho=i3-1
  if (tau(i3)>beta/pi) exit
end do

do i3=1,n_
  ncho1=i3-1
  if ((beta-tau(i3))*pi/beta<=1.0D0) exit
end do
 ncho=ncho1

 xx(1:n_)=-fome(n_:1:-1)
 xx((n_+1):2*n_)=fome(1:n_)
do i1=1,n_
  a=0;b=0;c=0;d=0;In0=0;In1=0;In2=0;In3=0
  do i2=1,n_
    y1(n_-i2+1)=inmatrix(i1,i2,1)     !  y1: real part
    y1(n_+i2)=inmatrix(i1,i2,1)
    y2(n_-i2+1)=-inmatrix(i1,i2,2)     !  y2: imaginary part
    y2(i2+n_)=inmatrix(i1,i2,2)
  end do

call spline(xx,y1,2*n_,y11)            ! y11: y1''
call spline(xx,y2,2*n_,y22)            ! y22: y2''

dw=2.0D0*pi/beta

do i2=n_+1,2*n_-1
  mid(1,1)=y1(i2)
  mid(1,2)=y2(i2)
  mid(2,1)=(y1(i2+1)-y1(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y11(i2)+y11(i2+1))/6.0D0
  mid(2,2)=(y2(i2+1)-y2(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y22(i2)+y22(i2+1))/6.0D0
  mid(3,1)=y11(i2)/2.0D0
  mid(3,2)=y22(i2)/2.0D0
  mid(4,1)=(y11(i2+1)-y11(i2))/(xx(i2+1)-xx(i2))/6.0D0
  mid(4,2)=(y22(i2+1)-y22(i2))/(xx(i2+1)-xx(i2))/6.0D0

  a(i2-n_)=cmplx(mid(1,1),mid(1,2))
  b(i2-n_)=cmplx(mid(2,1),mid(2,2))!*dw
  c(i2-n_)=cmplx(mid(3,1),mid(3,2))!*dw**2
  d(i2-n_)=cmplx(mid(4,1),mid(4,2))!*dw**3
end do

do i3=1,n_!ncho1
  f1(i3)=0
  do i2=1,19!(n_/2)
  outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+2.0D0*real(cmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(cmplx(0.0D0,-fome(i2)*tau(i3))))/beta
  end do
  do i2=20,n_-1!(n_/2+1),n_-1
    coe1=exp(cmplx(0.0D0,-fome(i2)*tau(i3)))
    coe2=exp(cmplx(0.0D0,-tau(i3)*(fome(i2+1)-fome(i2))))
    coe3=cmplx(dw*tau(i3)/2.0D0,0.0D0)
    tan3=sin(coe3)/cos(coe3)
    ooo=fome(i2+1)-fome(i2)
    if ( abs(tau(i3))-beta/2.0D0 .le. 0.0D0) then
      if (abs(tau(i3)*(fome(i2+1)-fome(i2))) .gt. 1.0D-1) then
        In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2+cmplx(-1.0D0,0.0D0)))
        In1(i2)=coe1*cmplx(0.0D0,1.0D0)*(cmplx(0.5D0*dw*ooo,0.0D0)*coe2/tan3-cmplx(0.0D0,0.25D0*dw**2)*(coe2+cmplx(-1.0D0,0.0D0))/sin(coe3)**2)
        In2(i2)=coe1*cmplx(0.0D0,-0.5D0*dw)*(-coe2*cmplx(ooo,0.0D0)**2/tan3+cmplx(0.0D0,dw*ooo)*coe2/sin(coe3)**2+cmplx(0.5D0*dw**2,0.0D0)*(coe2+cmplx(-1.0D0,0.0D0))/tan3/sin(coe3)**2)
        In3(i2)=exp(cmplx(0.0D0,-fome(i2+1)*tau(i3)))*dw/8D0*(cmplx(0D0,4D0)*ooo**3/tan3+dw/sin(coe3)**2*(-2D0*dw**2*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))+6D0*ooo**2+3D0*dw/sin(coe3)**2*(dw*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))-cmplx(0D0,1D0)*ooo*sin(dw*tau(i3)))))
        ! In3(i2)=coe1*cmplx(0.5D0*dw,0.0D0)*(cmplx(0.0D0,ooo**3)/tan3*coe2+cmplx(1.5D0*dw*ooo**2,0.0D0)*coe2/sin(coe3)**2&
        !    &-cmplx(0.0D0,1.5D0*dw**2*ooo)*coe2/tan3/sin(coe3)**2+(coe2+cmplx(-1.0D0,0.0D0))*(-cmplx(0.5D0*dw**3,0.0D0)/tan3**2/sin(coe3)**2-cmplx(0.25D0*dw**3,0.0D0)/sin(coe3)**4))
      else
        In0(i2)=coe1*cmplx(ooo-ooo*tau(i3)**2*(2.0D0*ooo**2+dw**2)/12.0D0+tau(i3)**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*tau(i3)**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,&
              &-ooo**2*tau(i3)/2.0D0+ooo**2*tau(i3)**3*(ooo**2+dw**2)/24.0D0+tau(i3)**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
        In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*tau(i3)**2*(ooo**2+dw**2)/8.0D0+tau(i3)**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,&
              &-tau(i3)*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+tau(i3)**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*tau(i3)**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
        In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-tau(i3)**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*tau(i3)**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008.0D0,&
              &-ooo**2*tau(i3)*(ooo**2+dw**2)/4.0D0+tau(i3)**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-tau(i3)**5*(3.0D0*ooo**8+14.0D0*ooo**6*dw**2-7.0D0*ooo**4*dw**4+2.0D0*ooo**2*dw**6)/2880.0D0)
        In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-tau(i3)**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+tau(i3)**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,&
              &-tau(i3)*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*tau(i3)**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252.0D0)
      end if
    else
      if (abs((beta-tau(i3))*(fome(i2+1)-fome(i2))) .gt. 1.0D-1) then
        In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2+cmplx(-1.0D0,0.0D0)))
        In1(i2)=exp(cmplx(0.0D0,-fome(i2+1)*tau(i3)))*dw/4D0/sin(coe3)**2*(dw-dw*exp(cmplx(0.0D0,tau(i3)*ooo))+cmplx(0D0,ooo)*sin(dw*tau(i3)))
        In2(i2)=exp(cmplx(0.0D0,-fome(i2+1)*tau(i3)))*dw/4D0*(cmplx(0D0,2D0)*ooo**2/tan3+dw*(2D0*ooo+cmplx(0D0,dw)*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))/tan3)/sin(coe3)**2)
        In3(i2)=exp(cmplx(0.0D0,-fome(i2+1)*tau(i3)))*dw/8D0*(cmplx(0D0,4D0)*ooo**3/tan3+dw/sin(coe3)**2*(-2D0*dw**2*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))+6D0*ooo**2+3D0*dw/sin(coe3)**2*(dw*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))-cmplx(0D0,1D0)*ooo*sin(dw*tau(i3)))))
      else
        ttt=beta-tau(i3)
        In0(i2)=coe1*cmplx(ooo-ooo*ttt**2*(2.0D0*ooo**2+dw**2)/12.0D0+ttt**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*ttt**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,-ooo**2*ttt/2.0D0+ooo**2*ttt**3*(ooo**2+dw**2)/24.0D0+ttt**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
        In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*ttt**2*(ooo**2+dw**2)/8.0D0+ttt**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,-ttt*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+ttt**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*ttt**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
        In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-ttt**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*ttt**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008D0,-ooo**2*ttt*(ooo**2+dw**2)/4.0D0+ttt**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-ttt**5*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/2880D0)
        In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-ttt**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+ttt**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,-ttt*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*ttt**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252D0)
      end if
    end if

      step=In0(i2)*a(i2)+In1(i2)*b(i2)+In2(i2)*c(i2)+In3(i2)*d(i2)
      f1(i3)=f1(i3)+step
   end do
  !  f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)

   outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+(2.0D0*real(f1(i3)))/(2.0D0*pi)                    !n from -inf to inf
   outmatrix(i1,i3,2)=0.0D0!aimag(f1(i3))
end do

! do i3=ncho+1,ncho1
!   coe3=cmplx(dw*tau(i3)/2.0D0,0.0D0)
!   co1 = cmplx(0.0D0,0.5D0*dw)*cos(coe3)/sin(coe3)
!   co2 = cmplx(dw**2/4.0D0,0.0D0)/sin(coe3)**2
!   co3 = cmplx(0.0D0,-dw**3/4.0D0)*cos(coe3)/sin(coe3)**3
!   co4 = cmplx(-dw**4/4.0D0,0.0D0)*(cos(coe3)**2+cmplx(0.5D0,0.0D0))/sin(coe3)**4
!   step=cmplx(0.0D0,0.0D0)
!   do i2=1,n_-1
!     step = step+d(i2)*(exp(cmplx(0.0D0,-tau(i3)*fome(i2+1)))-exp(cmplx(0.0D0,-tau(i3)*fome(i2))))
!   end do
!   coe1=exp(cmplx(0.0D0,-fome(n_)*tau(i3)))
!   coe2=exp(cmplx(0.0D0,-fome(1)*tau(i3)))
!   f1(i3)=co1*(coe1*a(n_-1)-coe2*a(1))+co2*(coe1*b(n_-1)-coe2*b(1))+co3*(coe1*c(n_-1)-coe2*c(1))+co4*step
!
!   outmatrix(i1,i3,1)=2.0D0*real(f1(i3))
!   outmatrix(i1,i3,2)=0.0D0!aimag(f1(i3))
! end do
!-----------------------------tau->beta,  sin(dw*tau/2)->0  ------------------
! do i3=ncho1+1,n_
!   In0=0;In1=0;In2=0;In3=0
!   do i2=1,(n_/2)
!   outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+2.0D0*real(cmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(cmplx(0.0D0,-fome(i2)*tau(i3))))/beta
!   end do
!   do i2=(n_/2+1),n_-1
!     ttt=beta-tau(i3)
!     coe1=exp(cmplx(0.0D0,-fome(i2)*ttt))
!     coe2=exp(cmplx(0.0D0,-ttt*(fome(i2+1)-fome(i2))))
!     coe3=cmplx(dw*ttt/2.0D0,0.0D0)
!     tan3=sin(coe3)/cos(coe3)
!     ooo=fome(i2+1)-fome(i2)
!     if (abs(ttt*(fome(i2+1)-fome(i2))) .gt. 1.0D-1) then
!       In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2-cmplx(1.0D0,0.0D0)))
!       In1(i2)=coe1*(cmplx(0.5D0*dw*(fome(i2+1)-fome(i2)),0.0D0)*coe2/tan3*cmplx(0.0D0,1.0D0)+cmplx(dw**2/4.0D0,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))/sin(coe3)**2)
!       In2(i2)=coe1*cmplx(0.0D0,-0.5D0*dw)*(coe2*cmplx(-(fome(i2+1)-fome(i2))**2,0.0D0)/tan3+cmplx(0.0D0,dw*(fome(i2+1)-fome(i2)))*coe2/sin(coe3)**2+cmplx(0.5D0*dw**2,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))/tan3/sin(coe3)**2)
!       In3(i2)=coe1*cmplx(0.5*dw,0.0D0)*(cmplx(0.0D0,(fome(i2+1)-fome(i2))**3)/tan3*coe2+cmplx(1.5D0*dw*(fome(i2+1)-fome(i2))**2,0.0D0)*coe2/sin(coe3)**2&
!          &-cmplx(0.0D0,1.5D0*dw**2*(fome(i2+1)-fome(i2)))*coe2/tan3/sin(coe3)**2-cmplx(dw**3/4.0D0,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))*(cmplx(2.0D0,0.0D0)/tan3**2/sin(coe3)**2+cmplx(1.0D0,0.0D0)/sin(coe3)**4))
!     else
!       In0(i2)=coe1*cmplx(ooo-ooo*ttt**2*(2.0D0*ooo**2+dw**2)/12.0D0+ttt**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*ttt**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,-ooo**2*ttt/2.0D0+ooo**2*ttt**3*(ooo**2+dw**2)/24.0D0+ttt**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
!       In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*ttt**2*(ooo**2+dw**2)/8.0D0+ttt**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,-ttt*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+ttt**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*ttt**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
!       In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-ttt**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*ttt**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008D0,-ooo**2*ttt*(ooo**2+dw**2)/4.0D0+ttt**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-ttt**5*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/2880D0)
!       In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-ttt**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+ttt**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,-ttt*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*ttt**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252D0)
!     end if
!    end do
!    f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)
!
!    outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+(2.0D0*real(f1(i3)))/(2.0D0*pi)
!    outmatrix(i1,i3,2)=0.0D0
! end do
!---------------------------------------------------
end do
!outmatrix=outmatrix/(2.0D0*pi)
end subroutine SFT_omegaFtotau


subroutine SFT_omegaBtotau(inmatrix,outmatrix,beta)
!SFT_omegaBtotau transfer a function from omega_n-space to tau-space
use vec
use omp_lib
implicit none
real*8,intent(in):: inmatrix(n_,n_,2),beta
real*8,intent(out):: outmatrix(n_,n_,2)
real*8::xx(2*n_-2),y1(2*n_-2),y2(2*n_-2),y11(2*n_-2),y22(2*n_-2)
! real*8::y1(2*n_),y2(2*n_),y11(2*n_),y22(2*n_),xx(2*n_)
real*8::mid(4,2),dw,ttt,ooo,bo0
integer i1,i2,i3,ncho,ncho1,ncho0
complex*16::In0(n_),In1(n_),In2(n_),In3(n_),a(n_),b(n_),c(n_),d(n_),f1(n_),f0(n_)
complex*16::co1,co2,co3,co4,step,coe1,coe2,coe3,tan3
! i1 index k, i2 index bomega ,i3 index tau
outmatrix=0.0D0

! do i3=1,n_
!   ncho0=i3-1
!   if (tau(i3)>1.0D-1) exit
! end do
! ncho0=0
! do i3=1,n_
!   ncho=i3-1
!   if (tau(i3)>beta/(2.0D0*pi)) exit
! end do
!
! do i3=1,n_
!   ncho1=i3-1
!   if ((beta-tau(i3))*2.0D0*pi/beta<=1.0D0) exit
! end do
! ncho=ncho1
!!!$OMP PARALLEL DO PRIVATE(i1,i2,i3,coe1,coe2,In0,In1,In2,In3,a,b,c,d,f1,co1,co2,co3,co4,y1,y2,y11,y22,mid,step,dw,tan3,ttt,ooo) SCHEDULE(DYNAMIC,1) REDUCTION(+:outmatrix)
dw=2.0D0*pi/beta

xx(1:n_-1)=-bome((n_-1):1:-1)
xx(n_:2*n_-2)=bome(1:n_-1)
do i1=1,n_
 a=0;b=0;c=0;d=0;In0=0;In1=0;In2=0;In3=0
 do i2=1,n_-1
   y1(n_-i2)=inmatrix(i1,i2,1)     !  y1: real part
   y1(n_+i2-1)=inmatrix(i1,i2,1)
   y2(n_-i2)=-inmatrix(i1,i2,2)     !  y2: imaginary part
   y2(i2+n_-1)=inmatrix(i1,i2,2)
 end do
 call spline(xx,y1,2*n_-2,y11)            ! y11: y1''
 call spline(xx,y2,2*n_-2,y22)            ! y22: y2''

 do i2=n_,2*n_-3
   mid(1,1)=y1(i2)
   mid(1,2)=y2(i2)
   mid(2,1)=(y1(i2+1)-y1(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y11(i2)+y11(i2+1))/6.0D0
   mid(2,2)=(y2(i2+1)-y2(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y22(i2)+y22(i2+1))/6.0D0
   mid(3,1)=y11(i2)/2.0D0
   mid(3,2)=y22(i2)/2.0D0
   mid(4,1)=(y11(i2+1)-y11(i2))/(xx(i2+1)-xx(i2))/6.0D0
   mid(4,2)=(y22(i2+1)-y22(i2))/(xx(i2+1)-xx(i2))/6.0D0

   a(i2-n_+1)=cmplx(mid(1,1),mid(1,2))
   b(i2-n_+1)=cmplx(mid(2,1),mid(2,2))!*dw
   c(i2-n_+1)=cmplx(mid(3,1),mid(3,2))!*dw**2
   d(i2-n_+1)=cmplx(mid(4,1),mid(4,2))!*dw**3
 end do

 mid(1,1)=y1(n_-1)
 mid(2,1)=(y1(n_)-y1(n_-1))/(xx(n_)-xx(n_-1))-(xx(n_)-xx(n_-1))*(2.0D0*y11(n_-1)+y11(n_))/6.0D0
 mid(3,1)=y11(n_-1)/2.0D0
 mid(4,1)=(y11(n_)-y11(n_-1))/(xx(n_)-xx(n_-1))/6.0D0
 bo0=mid(1,1)-mid(2,1)*xx(n_-1)+mid(3,1)*xx(n_-1)**2-mid(4,1)*xx(n_-1)**3
 ! bo0=(2D0*inmatrix(i1,n_,1))/(1D0-sqrt(1D0-4D0*pi**2*inmatrix(i1,n_,1)**2*bome(n_)**2))

do i3=1,n_!ncho1
  In0=0;In1=0;In2=0;In3=0
  f1(i3)=0
  ! outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+real(cmplx(inmatrix(i1,1,1),inmatrix(i1,1,2))*exp(cmplx(0.0D0,-bome(1)*tau(i3))))/beta        ! bome=0
  ! outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+inmatrix(i1,1,1)/beta
  do i2=1,19!(n_/2)
  outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+2.0D0*real(cmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(cmplx(0.0D0,-bome(i2)*tau(i3))))/beta
  end do
  do i2=20,n_-1!(n_/2+1),n_-2
    coe1=exp(cmplx(0.0D0,-bome(i2)*tau(i3)))
    coe2=exp(cmplx(0.0D0,-tau(i3)*(bome(i2+1)-bome(i2))))
    coe3=cmplx(dw*tau(i3)/2.0D0,0.0D0)
    tan3=sin(coe3)/cos(coe3)
    ooo=bome(i2+1)-bome(i2)
    if ( abs(tau(i3))-beta/2.0D0 .le. 0.0D0) then
      if (abs(tau(i3)*(bome(i2+1)-bome(i2))) .gt. 1.0D-1) then
        In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2+cmplx(-1.0D0,0.0D0)))
        In1(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/4D0/sin(coe3)**2*(dw-dw*exp(cmplx(0.0D0,tau(i3)*ooo))+cmplx(0D0,ooo)*sin(dw*tau(i3)))
        In2(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/4D0*(cmplx(0D0,2D0)*ooo**2/tan3+dw*(2D0*ooo+cmplx(0D0,dw)*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))/tan3)/sin(coe3)**2)
        In3(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/8D0*(cmplx(0D0,4D0)*ooo**3/tan3+dw/sin(coe3)**2*(-2D0*dw**2*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))+6D0*ooo**2+3D0*dw/sin(coe3)**2*(dw*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))-cmplx(0D0,1D0)*ooo*sin(dw*tau(i3)))))
        ! In1(i2)=coe1*cmplx(0.0D0,1.0D0)*(cmplx(0.5D0*dw*ooo,0.0D0)*coe2/tan3-cmplx(0.0D0,0.25D0*dw**2)*(coe2+cmplx(-1.0D0,0.0D0))/sin(coe3)**2)
        ! In2(i2)=coe1*cmplx(0.0D0,-0.5D0*dw)*(-coe2*cmplx(ooo,0.0D0)**2/tan3+cmplx(0.0D0,dw*ooo)*coe2/sin(coe3)**2+cmplx(0.5D0*dw**2,0.0D0)*(coe2+cmplx(-1.0D0,0.0D0))/tan3/sin(coe3)**2)
        ! In3(i2)=coe1*cmplx(0.5D0*dw,0.0D0)*(cmplx(0.0D0,ooo**3)/tan3*coe2+cmplx(1.5D0*dw*ooo**2,0.0D0)*coe2/sin(coe3)**2&
        !    &-cmplx(0.0D0,1.5D0*dw**2*ooo)*coe2/tan3/sin(coe3)**2+(coe2+cmplx(-1.0D0,0.0D0))*(-cmplx(0.5D0*dw**3,0.0D0)/tan3**2/sin(coe3)**2-cmplx(0.25D0*dw**3,0.0D0)/sin(coe3)**4))
      else
        ttt=tau(i3)
        In0(i2)=coe1*cmplx(ooo-ooo*ttt**2*(2.0D0*ooo**2+dw**2)/12.0D0+ttt**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*ttt**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,-ooo**2*ttt/2.0D0+ooo**2*ttt**3*(ooo**2+dw**2)/24.0D0+ttt**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
        In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*ttt**2*(ooo**2+dw**2)/8.0D0+ttt**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,-ttt*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+ttt**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*ttt**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
        In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-ttt**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*ttt**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008D0,-ooo**2*ttt*(ooo**2+dw**2)/4.0D0+ttt**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-ttt**5*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/2880D0)
        In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-ttt**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+ttt**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,-ttt*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*ttt**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252D0)
      end if
    else
      if (abs((beta-tau(i3))*(bome(i2+1)-bome(i2))) .gt. 1.0D-1) then
        In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2+cmplx(-1.0D0,0.0D0)))
        In1(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/4D0/sin(coe3)**2*(dw-dw*exp(cmplx(0.0D0,tau(i3)*ooo))+cmplx(0D0,ooo)*sin(dw*tau(i3)))
        In2(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/4D0*(cmplx(0D0,2D0)*ooo**2/tan3+dw*(2D0*ooo+cmplx(0D0,dw)*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))/tan3)/sin(coe3)**2)
        In3(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/8D0*(cmplx(0D0,4D0)*ooo**3/tan3+dw/sin(coe3)**2*(-2D0*dw**2*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))+6D0*ooo**2+3D0*dw/sin(coe3)**2*(dw*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))-cmplx(0D0,1D0)*ooo*sin(dw*tau(i3)))))
      else
        ttt=beta-tau(i3)
        In0(i2)=coe1*cmplx(ooo-ooo*ttt**2*(2.0D0*ooo**2+dw**2)/12.0D0+ttt**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*ttt**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,-ooo**2*ttt/2.0D0+ooo**2*ttt**3*(ooo**2+dw**2)/24.0D0+ttt**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
        In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*ttt**2*(ooo**2+dw**2)/8.0D0+ttt**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,-ttt*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+ttt**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*ttt**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
        In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-ttt**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*ttt**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008D0,-ooo**2*ttt*(ooo**2+dw**2)/4.0D0+ttt**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-ttt**5*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/2880D0)
        In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-ttt**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+ttt**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,-ttt*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*ttt**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252D0)
      end if
    end if
      step=In0(i2)*a(i2)+In1(i2)*b(i2)+In2(i2)*c(i2)+In3(i2)*d(i2)
      f1(i3)=f1(i3)+step
   end do
  !  f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)
   !f0(i3)=a(1)*In0(1)+b(1)*In1(1)+c(1)*In2(1)+d(1)*In3(1)
   !outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+(2.0D0*real(f1(i3))-real(f0(i3)))/(2.0D0*pi)                    !n from -inf to inf  except 0!!
   outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+(2.0D0*real(f1(i3)))/(2.0D0*pi)+inmatrix(i1,n_,1)/beta!bo0/beta
   outmatrix(i1,i3,2)=0.0D0!aimag(cmplx(inmatrix(i1,1,1),inmatrix(i1,1,2))*exp(cmplx(0.0D0,-bome(1)*tau(i3))))/beta
end do

! do i3=ncho+1,ncho1
!   coe3=cmplx(dw*tau(i3)/2.0D0,0.0D0)
!   co1 = cmplx(0.0D0,0.5D0*dw)*cos(coe3)/sin(coe3)
!   co2 = cmplx(dw**2/4.0D0,0.0D0)/sin(coe3)**2
!   co3 = cmplx(0.0D0,-dw**3/4.0D0)*cos(coe3)/sin(coe3)**3
!   co4 = cmplx(-dw**4/4.0D0,0.0D0)*(cos(coe3)**2+cmplx(0.5D0,0.0D0))/sin(coe3)**4
!   step=cmplx(0.0D0,0.0D0)
!   do i2=2,n_-1
!     step = step+d(i2)*(exp(cmplx(0.0D0,-tau(i3)*bome(i2+1)))-exp(cmplx(0.0D0,-tau(i3)*bome(i2))))
!   end do
!   coe1=exp(cmplx(0.0D0,-bome(n_)*tau(i3)))
!   coe2=exp(cmplx(0.0D0,-bome(2)*tau(i3)))
!   f1(i3)=co1*(coe1*a(n_-1)-coe2*a(2))+co2*(coe1*b(n_-1)-coe2*b(2))+co3*(coe1*c(n_-1)-coe2*c(2))+co4*step
!
!   outmatrix(i1,i3,1)=2.0D0*real(f1(i3))
!   outmatrix(i1,i3,2)=0.0D0!aimag(f1(i3))
!
!   !-----------------bomega(1)=0-----
!   coe1=exp(cmplx(0.0D0,-bome(1)*tau(i3)))
!   coe2=exp(cmplx(0.0D0,-tau(i3)*(bome(2)-bome(1))))
!   coe3=cmplx(dw*tau(i3)/2.0D0,0.0D0)
!   In0(1)=coe1*(cmplx(0.0D0,0.5D0*dw)*cos(coe3)/sin(coe3)*(coe2-cmplx(1.0D0,0.0D0)))
!   In1(1)=coe1*(cmplx(0.5D0*dw*(bome(2)-bome(1)),0.0D0)*coe2*cos(coe3)/sin(coe3)*cmplx(0.0D0,1.0D0)+cmplx(dw**2/4.0D0,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))/sin(coe3)**2)
!   In2(1)=coe1*cmplx(0.0D0,-0.5D0*dw)*(coe2*cmplx(-(bome(2)-bome(1))**2,0.0D0)*cos(coe3)/sin(coe3)+cmplx(0.0D0,dw*(bome(2)-bome(1)))*coe2/sin(coe3)**2+cmplx(0.5D0*dw**2,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))*cos(coe3)/sin(coe3)**3)
!   In3(1)=coe1*cmplx(0.5*dw,0.0D0)*(cmplx(0.0D0,(bome(2)-bome(1))**3)*cos(coe3)/sin(coe3)*coe2+cmplx(1.5D0*dw*(bome(2)-bome(1))**2,0.0D0)*coe2/sin(coe3)**2&
!          &-cmplx(0.0D0,1.5D0*dw**2*(bome(2)-bome(1)))*coe2*cos(coe3)/sin(coe3)**3-cmplx(dw**3/4.0D0,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))*(cmplx(2.0D0,0.0D0)*cos(coe3)**2+cmplx(1.0D0,0.0D0))/sin(coe3)**4)
!   f0(i3)=a(1)*In0(1)+b(1)*In1(1)+c(1)*In2(1)+d(1)*In3(1)
!
!   outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+real(f0(i3))                   !  n=0
!   outmatrix(i1,i3,2)=0.0D0!outmatrix(i1,i3,2)+aimag(f0(i3))
!   !-------------------------------------
! end do
! !-----------------------------tau->beta,  sin(dw*tau/2)->0  ------------------
! do i3=ncho1+1,n_
!   In0=0;In1=0;In2=0;In3=0
!   f1(i3)=0
!   ! outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+real(cmplx(inmatrix(i1,1,1),inmatrix(i1,1,2))*exp(cmplx(0.0D0,-bome(1)*tau(i3))))/beta        ! bome=0
!   do i2=1,(n_/2)
!   outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+2.0D0*real(cmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(cmplx(0.0D0,-bome(i2)*tau(i3))))/beta
!   end do
!   ttt=beta-tau(i3)
!   do i2=(n_/2+1),n_-1
!     coe1=exp(cmplx(0.0D0,-bome(i2)*ttt))
!     coe2=exp(cmplx(0.0D0,-ttt*(bome(i2+1)-bome(i2))))
!     coe3=cmplx(dw*ttt/2.0D0,0.0D0)
!     tan3=sin(coe3)/cos(coe3)
!     ooo=bome(i2+1)-bome(i2)
!     if (abs(ttt*(bome(i2+1)-bome(i2))) .gt. 1.0D-1) then
!       In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2-cmplx(1.0D0,0.0D0)))
!       In1(i2)=coe1*(cmplx(0.5D0*dw*(bome(i2+1)-bome(i2)),0.0D0)*coe2/tan3*cmplx(0.0D0,1.0D0)+cmplx(dw**2/4.0D0,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))/sin(coe3)**2)
!       In2(i2)=coe1*cmplx(0.0D0,-0.5D0*dw)*(coe2*cmplx(-(bome(i2+1)-bome(i2))**2,0.0D0)/tan3+cmplx(0.0D0,dw*(bome(i2+1)-bome(i2)))*coe2/sin(coe3)**2+cmplx(0.5D0*dw**2,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))/tan3/sin(coe3)**2)
!       In3(i2)=exp(cmplx(0.0D0,-bome(i2+1)*ttt))*dw/8D0*(cmplx(0D0,4D0)*ooo**3/tan3+dw/sin(coe3)**2*(-2D0*dw**2*(exp(cmplx(0D0,ttt*ooo))-cmplx(1D0,0D0))+6D0*ooo**2+3D0*dw/sin(coe3)**2*(dw*(exp(cmplx(0D0,ttt*ooo))-cmplx(1D0,0D0))-cmplx(0D0,1D0)*ooo*sin(dw*ttt))))
!       ! In3(i2)=coe1*cmplx(0.5*dw,0.0D0)*(cmplx(0.0D0,(bome(i2+1)-bome(i2))**3)/tan3*coe2+cmplx(1.5D0*dw*(bome(i2+1)-bome(i2))**2,0.0D0)*coe2/sin(coe3)**2&
!       !         &-cmplx(0.0D0,1.5D0*dw**2*(bome(i2+1)-bome(i2)))*coe2/tan3/sin(coe3)**2-cmplx(dw**3/4.0D0,0.0D0)*(coe2-cmplx(1.0D0,0.0D0))*(cmplx(2.0D0,0.0D0)/tan3**2/sin(coe3)**2+cmplx(1.0D0,0.0D0)/sin(coe3)**4))
!     else
!       In0(i2)=coe1*cmplx(ooo-ooo*ttt**2*(2.0D0*ooo**2+dw**2)/12.0D0+ttt**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*ttt**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,-ooo**2*ttt/2.0D0+ooo**2*ttt**3*(ooo**2+dw**2)/24.0D0+ttt**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
!       In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*ttt**2*(ooo**2+dw**2)/8.0D0+ttt**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,-ttt*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+ttt**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*ttt**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
!       In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-ttt**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*ttt**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008D0,-ooo**2*ttt*(ooo**2+dw**2)/4.0D0+ttt**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-ttt**5*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/2880D0)
!       In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-ttt**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+ttt**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,-ttt*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*ttt**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252D0)
!     end if
!     step=In0(i2)*a(i2)+In1(i2)*b(i2)+In2(i2)*c(i2)+In3(i2)*d(i2)
!     f1(i3)=f1(i3)+step
!    end do
!   !  f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)
!    !f0(i3)=a(1)*In0(1)+b(1)*In1(1)+c(1)*In2(1)+d(1)*In3(1)
!
!    !outmatrix(i1,i3,1)=2.0D0*real(f1(i3))-real(f0(i3))                    !n from -inf to inf  except 0!!
!    !outmatrix(i1,i3,2)=aimag(f0(i3))
!    outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+(2.0D0*real(f1(i3)))/(2.0D0*pi)
!    outmatrix(i1,i3,2)=0.0D0!aimag(cmplx(inmatrix(i1,1,1),inmatrix(i1,1,2))*exp(cmplx(0.0D0,-bome(1)*tau(i3))))/beta
! end do
!---------------------------------------------------
end do
!!$OMP END PARALLEL DO
!outmatrix=outmatrix/(2.0D0*pi)
end subroutine SFT_omegaBtotau

subroutine SFT_omegaBtotau1(inmatrix,outmatrix,beta)
!SFT_omegaBtotau transfer a function from omega_n-space to tau-space
use vec
use omp_lib
implicit none
real*8,intent(in):: inmatrix(n_,n_,2),beta
real*8,intent(out):: outmatrix(n_,n_,2)
real*8::xx(2*n_-2),y1(2*n_-2),y2(2*n_-2),y11(2*n_-2),y22(2*n_-2)
! real*8::y1(2*n_),y2(2*n_),y11(2*n_),y22(2*n_),xx(2*n_)
real*8::mid(4,2),dw,ttt,ooo
integer i1,i2,i3,ncho,ncho1,ncho0
complex*16::In0(n_),In1(n_),In2(n_),In3(n_),a(n_),b(n_),c(n_),d(n_),f1(n_),f0(n_)
complex*16::co1,co2,co3,co4,step,coe1,coe2,coe3,tan3
! i1 index k, i2 index bomega ,i3 index tau
outmatrix=0.0D0

! do i3=1,n_
!   ncho0=i3-1
!   if (tau(i3)>1.0D-1) exit
! end do
! ncho0=0
! do i3=1,n_
!   ncho=i3-1
!   if (tau(i3)>beta/(2.0D0*pi)) exit
! end do
!
! do i3=1,n_
!   ncho1=i3-1
!   if ((beta-tau(i3))*2.0D0*pi/beta<=1.0D0) exit
! end do
! ncho=ncho1
!!!$OMP PARALLEL DO PRIVATE(i1,i2,i3,coe1,coe2,In0,In1,In2,In3,a,b,c,d,f1,co1,co2,co3,co4,y1,y2,y11,y22,mid,step,dw,tan3,ttt,ooo) SCHEDULE(DYNAMIC,1) REDUCTION(+:outmatrix)
dw=2.0D0*pi/beta

xx(1:n_-1)=-bome((n_-1):1:-1)
xx(n_:2*n_-2)=bome(1:n_-1)
do i1=1,n_
 a=0;b=0;c=0;d=0;In0=0;In1=0;In2=0;In3=0
 do i2=1,n_-1
   y1(n_-i2)=inmatrix(i1,i2,1)     !  y1: real part
   y1(n_+i2-1)=inmatrix(i1,i2,1)
   y2(n_-i2)=-inmatrix(i1,i2,2)     !  y2: imaginary part
   y2(i2+n_-1)=inmatrix(i1,i2,2)
 end do
 call spline(xx,y1,2*n_-2,y11)            ! y11: y1''
 call spline(xx,y2,2*n_-2,y22)            ! y22: y2''

 do i2=n_,2*n_-3
   mid(1,1)=y1(i2)
   mid(1,2)=y2(i2)
   mid(2,1)=(y1(i2+1)-y1(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y11(i2)+y11(i2+1))/6.0D0
   mid(2,2)=(y2(i2+1)-y2(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y22(i2)+y22(i2+1))/6.0D0
   mid(3,1)=y11(i2)/2.0D0
   mid(3,2)=y22(i2)/2.0D0
   mid(4,1)=(y11(i2+1)-y11(i2))/(xx(i2+1)-xx(i2))/6.0D0
   mid(4,2)=(y22(i2+1)-y22(i2))/(xx(i2+1)-xx(i2))/6.0D0

   a(i2-n_+1)=cmplx(mid(1,1),mid(1,2))
   b(i2-n_+1)=cmplx(mid(2,1),mid(2,2))!*dw
   c(i2-n_+1)=cmplx(mid(3,1),mid(3,2))!*dw**2
   d(i2-n_+1)=cmplx(mid(4,1),mid(4,2))!*dw**3
 end do


do i3=1,n_!ncho1
  In0=0;In1=0;In2=0;In3=0
  f1(i3)=0
  ! outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+real(cmplx(inmatrix(i1,1,1),inmatrix(i1,1,2))*exp(cmplx(0.0D0,-bome(1)*tau(i3))))/beta        ! bome=0
  ! outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+inmatrix(i1,1,1)/beta
  do i2=1,19!(n_/2)
  outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+2.0D0*real(cmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(cmplx(0.0D0,-bome(i2)*tau(i3))))/beta
  end do
  do i2=20,n_-1!(n_/2+1),n_-2
    coe1=exp(cmplx(0.0D0,-bome(i2)*tau(i3)))
    coe2=exp(cmplx(0.0D0,-tau(i3)*(bome(i2+1)-bome(i2))))
    coe3=cmplx(dw*tau(i3)/2.0D0,0.0D0)
    tan3=sin(coe3)/cos(coe3)
    ooo=bome(i2+1)-bome(i2)
    if ( abs(tau(i3))-beta/2.0D0 .le. 0.0D0) then
      if (abs(tau(i3)*(bome(i2+1)-bome(i2))) .gt. 1.0D-1) then
        In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2+cmplx(-1.0D0,0.0D0)))
        In1(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/4D0/sin(coe3)**2*(dw-dw*exp(cmplx(0.0D0,tau(i3)*ooo))+cmplx(0D0,ooo)*sin(dw*tau(i3)))
        In2(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/4D0*(cmplx(0D0,2D0)*ooo**2/tan3+dw*(2D0*ooo+cmplx(0D0,dw)*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))/tan3)/sin(coe3)**2)
        In3(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/8D0*(cmplx(0D0,4D0)*ooo**3/tan3+dw/sin(coe3)**2*(-2D0*dw**2*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))+6D0*ooo**2+3D0*dw/sin(coe3)**2*(dw*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))-cmplx(0D0,1D0)*ooo*sin(dw*tau(i3)))))
        ! In1(i2)=coe1*cmplx(0.0D0,1.0D0)*(cmplx(0.5D0*dw*ooo,0.0D0)*coe2/tan3-cmplx(0.0D0,0.25D0*dw**2)*(coe2+cmplx(-1.0D0,0.0D0))/sin(coe3)**2)
        ! In2(i2)=coe1*cmplx(0.0D0,-0.5D0*dw)*(-coe2*cmplx(ooo,0.0D0)**2/tan3+cmplx(0.0D0,dw*ooo)*coe2/sin(coe3)**2+cmplx(0.5D0*dw**2,0.0D0)*(coe2+cmplx(-1.0D0,0.0D0))/tan3/sin(coe3)**2)
        ! In3(i2)=coe1*cmplx(0.5D0*dw,0.0D0)*(cmplx(0.0D0,ooo**3)/tan3*coe2+cmplx(1.5D0*dw*ooo**2,0.0D0)*coe2/sin(coe3)**2&
        !    &-cmplx(0.0D0,1.5D0*dw**2*ooo)*coe2/tan3/sin(coe3)**2+(coe2+cmplx(-1.0D0,0.0D0))*(-cmplx(0.5D0*dw**3,0.0D0)/tan3**2/sin(coe3)**2-cmplx(0.25D0*dw**3,0.0D0)/sin(coe3)**4))
      else
        ttt=tau(i3)
        In0(i2)=coe1*cmplx(ooo-ooo*ttt**2*(2.0D0*ooo**2+dw**2)/12.0D0+ttt**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*ttt**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,-ooo**2*ttt/2.0D0+ooo**2*ttt**3*(ooo**2+dw**2)/24.0D0+ttt**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
        In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*ttt**2*(ooo**2+dw**2)/8.0D0+ttt**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,-ttt*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+ttt**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*ttt**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
        In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-ttt**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*ttt**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008D0,-ooo**2*ttt*(ooo**2+dw**2)/4.0D0+ttt**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-ttt**5*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/2880D0)
        In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-ttt**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+ttt**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,-ttt*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*ttt**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252D0)
      end if
    else
      if (abs((beta-tau(i3))*(bome(i2+1)-bome(i2))) .gt. 1.0D-1) then
        In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2+cmplx(-1.0D0,0.0D0)))
        In1(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/4D0/sin(coe3)**2*(dw-dw*exp(cmplx(0.0D0,tau(i3)*ooo))+cmplx(0D0,ooo)*sin(dw*tau(i3)))
        In2(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/4D0*(cmplx(0D0,2D0)*ooo**2/tan3+dw*(2D0*ooo+cmplx(0D0,dw)*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))/tan3)/sin(coe3)**2)
        In3(i2)=exp(cmplx(0.0D0,-bome(i2+1)*tau(i3)))*dw/8D0*(cmplx(0D0,4D0)*ooo**3/tan3+dw/sin(coe3)**2*(-2D0*dw**2*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))+6D0*ooo**2+3D0*dw/sin(coe3)**2*(dw*(exp(cmplx(0D0,tau(i3)*ooo))-cmplx(1D0,0D0))-cmplx(0D0,1D0)*ooo*sin(dw*tau(i3)))))
      else
        ttt=beta-tau(i3)
        In0(i2)=coe1*cmplx(ooo-ooo*ttt**2*(2.0D0*ooo**2+dw**2)/12.0D0+ttt**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*ttt**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,-ooo**2*ttt/2.0D0+ooo**2*ttt**3*(ooo**2+dw**2)/24.0D0+ttt**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
        In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*ttt**2*(ooo**2+dw**2)/8.0D0+ttt**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,-ttt*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+ttt**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*ttt**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
        In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-ttt**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*ttt**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008D0,-ooo**2*ttt*(ooo**2+dw**2)/4.0D0+ttt**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-ttt**5*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/2880D0)
        In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-ttt**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+ttt**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,-ttt*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*ttt**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252D0)
      end if
    end if
      step=In0(i2)*a(i2)+In1(i2)*b(i2)+In2(i2)*c(i2)+In3(i2)*d(i2)
      f1(i3)=f1(i3)+step
   end do
  !  f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)
   !f0(i3)=a(1)*In0(1)+b(1)*In1(1)+c(1)*In2(1)+d(1)*In3(1)
   !outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+(2.0D0*real(f1(i3))-real(f0(i3)))/(2.0D0*pi)                    !n from -inf to inf  except 0!!
   outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+(2.0D0*real(f1(i3)))/(2.0D0*pi)+inmatrix(i1,n_,1)/beta
   outmatrix(i1,i3,2)=0.0D0!aimag(cmplx(inmatrix(i1,1,1),inmatrix(i1,1,2))*exp(cmplx(0.0D0,-bome(1)*tau(i3))))/beta
end do

end do
!!$OMP END PARALLEL DO
!outmatrix=outmatrix/(2.0D0*pi)
end subroutine SFT_omegaBtotau1
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------

program main
  use Solve_NonLin
  use formidconstant
  implicit none
  real*8 beta,mu,num
  real*8 x(1),fvec(1),diag(1)
  integer info,i

  t=0.2D0
  ! beta=1.0D0/t
  eta=-0.1D0!-0.07D0
  ntot=2.0D5
  lambda= 0.045D0
  omega0=1.0D0/(3.0D0*lambda*ntot)**(1.0D0/3.0D0)

mu=0.5D0
do i=1,15
  write(14,*)i
  write(14,*)'t=',t
  write(13,*)i
  write(13,*)'t=',t
x(1)=mu
! do i=1,5
! call partinum(mu,num)
! fvec(1)=(num-ntot)/ntot
call hbrd(formu,1,x,fvec,1.0D-2,1.0D-3,info,diag)
write(14,*)'mu=',x(1)
mu=x(1)-0.05D0
write(*,*)'fvec=',fvec(1)
t=t+0.05D0
end do
contains
  subroutine formu(n0,x,fvec,iflag)
    use formidconstant
  implicit none
  integer,intent(in):: n0
  real*8,intent(in):: x(n0)
  real*8,intent(out):: fvec(n0)
  integer,intent(in out):: iflag
  real*8:: num

  call partinum(x(1),num)
  ! fvec(1)=(num-ntot)/ntot
  fvec(1)=1D0-num

  end subroutine formu
end program main

! subroutine partinum(mu,num)
!   use vec
!   use formidconstant
!   use gsld
!   use omp_lib
!   implicit none
!   real*8,intent(in):: mu
!   real*8,intent(out)::num
!   real*8 beta,vext,num1,num2
!   real*8::x(30),rho(30),u1(30),s1(30),mid
!   Real*8::a,b,answer(3),minl,maxl,r0(3),step
!   integer::n0
!   integer i,j,i1
!
!   beta=1.0D0/t
!   !$OMP PARALLEL SECTIONS
!   !$OMP SECTION
!   call vector(tau,1,beta)
!   !$OMP SECTION
!   call vector(fome,2,beta)    !(2n+1)pi/beta
!   !$OMP SECTION
!   call vector(bome,3,beta)    !(2n)pi/beta
!   !$OMP SECTION
!   call vector(k,4,beta)
!   !$OMP END PARALLEL SECTIONS
!   write(11,*)'tau=',tau
!   write(11,*)'fome=',fome
!   write(11,*)'bome=',bome
!   write(11,*)'k=',k
!   Call fn0(fn1, ak1)
!   write(*,*)'mu=',mu
!   ! do j = 1,30
!   !   x(j) = exp(dble(j)/dble(30)*5.0D0)-1.0D0
!   ! enddo
!   !!!!!!$OMP DO
!   ! do i=1,30
!   ! write(12,*)'i=',i
!   ! write (*,*) "i=",i
!   ! vext = 0.25D0*omega0**2*x(i)**2
!   ! call modeleq(mu,beta,vext,rho(i),u1(i),s1(i))
!   ! write (*,*) "x=",x(i),"rho=",rho(i),"u1=",u1(i),"s1=",s1(i)
!   ! end do
!   ! !!!!$OMP END DO
!   ! num=0.0D0
!   ! do i = 1,29
!   ! mid = (rho(i)*x(i)**2+rho(i+1)*x(i+1)**2)*(x(i+1)-x(i))/2.0D0
!   ! num = num+mid
!   ! enddo
!   !  num = -2.0D0*12.56637061D0*num/lambda
!   !  write(*,*)'num',num
!   !
!   !  num1=0.0D0
!   !  do i = 1,29
!   !  mid = (u1(i)*x(i)**2+u1(i+1)*x(i+1)**2)*(x(i+1)-x(i))/2.0D0
!   !  num1 = num1+mid
!   !  enddo
!   !   num1 = 12.56637061D0*num1/lambda/ntot!+mu
!   !   write(14,*)'U=',num1
!   !   write(*,*)'U',num1
!   !
!   !   num2=0.0D0
!   !   do i = 1,29
!   !   mid = (s1(i)*x(i)**2+s1(i+1)*x(i+1)**2)*(x(i+1)-x(i))/2.0D0
!   !   num2 = num2+mid
!   !   enddo
!   !    num2 = 12.56637061D0*num2/lambda/ntot
!   !    write(14,*)'s=',num2
!   !    write(*,*)'s',num2
!
!     r0 = 0.0D0
!     minl=0.0D0
!     maxl=3D0*sqrt(2D0*t/omega0)
!     n0=4
!     step=(maxl-minl)/dble(n0)
!     a=minl
!     b=minl+step
!     i1=0
!     do j=1,n0
!       answer=0
!       Do i=1,nk
!         i1=i1+1
!         write(12,*)'i=',i1
!         write (*,*) "i=",i1
!         answer=answer+ak1(i)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),mu,beta)
!       End Do
!       answer=answer*(b-a)/2.0D0
!       r0=r0+answer
!       a=b
!       b=b+step
!     end do
!     r0=r0*8.0D0*pi*sqrt(6.0D0*ntot/lambda)
!     num=r0(1)
!     write(*,*)'num',num
!     num1=r0(2)+mu
!     write(14,*)'U=',num1
!     write(*,*)'U',num1
!     num2=r0(3)
!     write(14,*)'S=',num2
!     write(*,*)'S',num2
!   contains
!     function y(r,mu,beta)
!       Implicit none
!       Real*8:: y(3)
!       real*8::r,mu,beta,vext,rho,u1,s1
!       vext = 0.5D0*omega0*r**2
!       call modeleq(mu,beta,vext,rho,u1,s1)
!       write (*,*) "r=",r,"rho=",rho,"u1=",u1,"s1=",s1
!       y(1)=-2.0D0*rho*r**2
!       y(2)=u1*r**2/ntot
!       y(3)=s1*r**2/ntot
!     end function y
! end subroutine partinum


subroutine partinum(mu,num)
  use vec
  use formidconstant
  use gsld
  use omp_lib
  implicit none
  real*8,intent(in):: mu
  real*8,intent(out)::num
  real*8 beta,vext,num1,num2
  real*8::x,rho,u1,s1,mid
  Real*8::a,b,answer(3),minl,maxl,r0(3),step
  integer::n0
  integer i,j,i1

  beta=1.0D0/t
  !$OMP PARALLEL SECTIONS
  !$OMP SECTION
  call vector(tau,1,beta)
  !$OMP SECTION
  call vector(fome,2,beta)    !(2n+1)pi/beta
  !$OMP SECTION
  call vector(bome,3,beta)    !(2n)pi/beta
  !$OMP SECTION
  call vector(k,4,beta)
  !$OMP END PARALLEL SECTIONS
  write(11,*)'tau=',tau
  write(11,*)'fome=',fome
  write(11,*)'bome=',bome
  write(11,*)'k=',k
  Call fn0(fn1, ak1)
  write(*,*)'mu=',mu
  call modeleq(mu,beta,0D0,rho,u1,s1)
  num = -2.0D0*rho*3D0*pi**2
  num1=u1*3D0*pi**2+mu
  num2=s1*3D0*pi**2
  write (*,*) "rho=",num,"u1=",num1,"s1=",num2
  write(14,*)'U=',num1
  write(14,*)'S=',num2
end subroutine partinum
