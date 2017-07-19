module vec
  integer,parameter::n_=500
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
   xx(j) =exp(-dble(n-j)/dble(n)*6.0D0)*(beta/2.0D0-0.02D0)
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
  n=n_/2
  do j = 1,n
    xx(j)=dble(2*j-1)*pi/beta
  enddo
  xx(n+1)=dble(2*(n+2)-1)*pi/beta
  do j = n+2,n_
    xx(j) = (2.0D0*dble(dint(exp(dble(j-n)/dble(n_-n)*6.0D0+3.5D0)-exp(3.5D0))+j+1)-1.0D0)*pi/beta
  enddo
! do j = n+1,n_
!   xx(j) = (2.0D0*(exp(dble(j-n)/dble(n_-n)*8.0D0)+dble(j)-1.0D0)-1.0D0)*pi/beta
! enddo
! do j = 1,n_
!   xx(j) = (2.0D0*(exp(dble(j-1)/dble(n_-1)*12.0D0)+dble(j)-1.0D0)-1.0D0)*pi/beta
! enddo
  !   omegaB  (2n)pi/beta  n: e^0~e^8
case(3)
  n=n_/2
  do j = 1,n
    xx(j)=dble(2*j)*pi/beta
  enddo
  xx(n+1)=dble(2*(n+2))*pi/beta
  do j = n+2,n_-1
    xx(j) = 2.0D0*dble(dint(exp(dble(j-n)/dble(n_-n)*6.0D0+3.5D0)-exp(3.5D0))+j+1)*pi/beta
  enddo
  xx(n_)=0D0
! do j = n+1,n_
!   xx(j) = 2.0D0*(exp(dble(j-n)/dble(n_-n)*8.0D0)+dble(j)-1.0D0)*pi/beta
! enddo

  !   k  0-e^6 in logarithmic scales
case(4)
do j = 1,n_
  !xx(j) = dble(j)/dble(n_)*exp(dble(j)/dble(n_)*7.0D0)!-1.0D0
   xx(j) = exp(dble(j)/dble(n_)*6.5D0)-1.01D0
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
  do i2=1,(n_/2)
  outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+2.0D0*real(cmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(cmplx(0.0D0,-fome(i2)*tau(i3))))/beta
  end do
  do i2=(n_/2+1),n_-1
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

! xx(1:n_-1)=-bome(n_:2:-1)
! xx(n_:2*n_-1)=bome(1:n_)
! do i1=1,n_
!   a=0;b=0;c=0;d=0;
!   do i2=2,n_
!     y1(n_-i2+1)=inmatrix(i1,i2,1)     !  y1: real part
!     y2(n_-i2+1)=-inmatrix(i1,i2,2)     !  y2: imaginary part
!   end do
!   do i2=1,n_
!     y1(n_+i2-1)=inmatrix(i1,i2,1)
!     y2(i2+n_-1)=inmatrix(i1,i2,2)
!   end do
!   y2(n_)=-y2(n_)
  ! call spline(xx,y1,2*n_-1,y11)            ! y11: y1''
  ! call spline(xx,y2,2*n_-1,y22)            ! y22: y2''
! write(13,*)'y11',y11,'y22',y22


! do i2=n_,2*n_-2
!   mid(1,1)=y1(i2)
!   mid(1,2)=y2(i2)
!   mid(2,1)=(y1(i2+1)-y1(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y11(i2)+y11(i2+1))/6.0D0
!   mid(2,2)=(y2(i2+1)-y2(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y22(i2)+y22(i2+1))/6.0D0
!   mid(3,1)=y11(i2)/2.0D0
!   mid(3,2)=y22(i2)/2.0D0
!   mid(4,1)=(y11(i2+1)-y11(i2))/(xx(i2+1)-xx(i2))/6.0D0
!   mid(4,2)=(y22(i2+1)-y22(i2))/(xx(i2+1)-xx(i2))/6.0D0
!
!   a(i2-n_+1)=cmplx(mid(1,1),mid(1,2))
!   b(i2-n_+1)=cmplx(mid(2,1),mid(2,2))!*dw
!   c(i2-n_+1)=cmplx(mid(3,1),mid(3,2))!*dw**2
!   d(i2-n_+1)=cmplx(mid(4,1),mid(4,2))!*dw**3
! end do

do i3=1,n_!ncho1
  In0=0;In1=0;In2=0;In3=0
  f1(i3)=0
  ! outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+real(cmplx(inmatrix(i1,1,1),inmatrix(i1,1,2))*exp(cmplx(0.0D0,-bome(1)*tau(i3))))/beta        ! bome=0
  ! outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+inmatrix(i1,1,1)/beta
  do i2=1,(n_/2)
  outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+2.0D0*real(cmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(cmplx(0.0D0,-bome(i2)*tau(i3))))/beta
  end do
  do i2=(n_/2+1),n_-2
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
   outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+(2.0D0*real(f1(i3)))/(2.0D0*pi)!+inmatrix(i1,n_,1)/beta
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

!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------

! subroutine SFT_omegaFtotau(inmatrix,outmatrix,beta)
! !SFT_omegaFtotau transfer a function from omega_n-space to tau-space, the second F means Fermion
! use vec
! implicit none
! real*8,intent(in):: inmatrix(n_,n_,2),beta
! real*8,intent(out):: outmatrix(n_,n_,2)
! real*8::y1(2*n_),y2(2*n_),y11(2*n_),y22(2*n_),xx(2*n_)
! real*8::mid(4,2),dw,ttt,ooo
! integer i1,i2,i3,ncho,ncho1,ncho0
! complex*16::In0(n_),In1(n_),In2(n_),In3(n_),a(n_),b(n_),c(n_),d(n_),f1(n_)
! complex*16::co1,co2,co3,co4,step,coe1,coe2,coe3,tan3
! ! i1 index k, i2 index fomega ,i3 index tau
!
! outmatrix = 0.0D0
!
!
!  xx(1:n_)=-fome(n_:1:-1)
!  xx((n_+1):2*n_)=fome(1:n_)
! do i1=1,n_
!   a=0;b=0;c=0;d=0;In0=0;In1=0;In2=0;In3=0
!   do i2=1,n_
!     y1(n_-i2+1)=inmatrix(i1,i2,1)     !  y1: real part
!     y1(n_+i2)=inmatrix(i1,i2,1)
!     y2(n_-i2+1)=-inmatrix(i1,i2,2)     !  y2: imaginary part
!     y2(i2+n_)=inmatrix(i1,i2,2)
!   end do
!
! call spline(xx,y1,2*n_,y11)            ! y11: y1''
! call spline(xx,y2,2*n_,y22)            ! y22: y2''
!
! dw=2.0D0*pi/beta
!
! do i2=n_+1,2*n_-1
!   mid(1,1)=y1(i2)
!   mid(1,2)=y2(i2)
!   mid(2,1)=(y1(i2+1)-y1(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y11(i2)+y11(i2+1))/6.0D0
!   mid(2,2)=(y2(i2+1)-y2(i2))/(xx(i2+1)-xx(i2))-(xx(i2+1)-xx(i2))*(2.0D0*y22(i2)+y22(i2+1))/6.0D0
!   mid(3,1)=y11(i2)/2.0D0
!   mid(3,2)=y22(i2)/2.0D0
!   mid(4,1)=(y11(i2+1)-y11(i2))/(xx(i2+1)-xx(i2))/6.0D0
!   mid(4,2)=(y22(i2+1)-y22(i2))/(xx(i2+1)-xx(i2))/6.0D0
!
!   a(i2-n_)=cmplx(mid(1,1),mid(1,2))
!   b(i2-n_)=cmplx(mid(2,1),mid(2,2))!*dw
!   c(i2-n_)=cmplx(mid(3,1),mid(3,2))!*dw**2
!   d(i2-n_)=cmplx(mid(4,1),mid(4,2))!*dw**3
! end do
!
! do i3=1,n_
!   do i2=1,(n_/2)
!   outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+2.0D0*real(cmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(cmplx(0.0D0,-fome(i2)*tau(i3))))/beta
!   end do
!   do i2=(n_/2+1),n_-1
!     coe1=exp(cmplx(0.0D0,-fome(i2)*tau(i3)))
!     coe2=exp(cmplx(0.0D0,-tau(i3)*(fome(i2+1)-fome(i2))))
!     coe3=cmplx(dw*tau(i3)/2.0D0,0.0D0)
!     tan3=sin(coe3)/cos(coe3)
!     ooo=fome(i2+1)-fome(i2)
!   if ( abs(tau(i3))-beta/2.0D0 .lt. 0.0D0) then
!       if (abs(tau(i3)*(fome(i2+1)-fome(i2))) .gt. 1.0D-1) then
!         In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2+cmplx(-1.0D0,0.0D0)))
!         In1(i2)=coe1*cmplx(0.0D0,1.0D0)*(cmplx(0.5D0*dw*ooo,0.0D0)*coe2/tan3-cmplx(0.0D0,0.25D0*dw**2)*(coe2+cmplx(-1.0D0,0.0D0))/sin(coe3)**2)
!         In2(i2)=coe1*cmplx(0.0D0,-0.5D0*dw)*(-coe2*cmplx(ooo,0.0D0)**2/tan3+cmplx(0.0D0,dw*ooo)*coe2/sin(coe3)**2+cmplx(0.5D0*dw**2,0.0D0)*(coe2+cmplx(-1.0D0,0.0D0))/tan3/sin(coe3)**2)
!         In3(i2)=coe1*cmplx(0.5D0*dw,0.0D0)*(cmplx(0.0D0,ooo**3)/tan3*coe2+cmplx(1.5D0*dw*ooo**2,0.0D0)*coe2/sin(coe3)**2&
!            &-cmplx(0.0D0,1.5D0*dw**2*ooo)*coe2/tan3/sin(coe3)**2+(coe2+cmplx(-1.0D0,0.0D0))*(-cmplx(0.5D0*dw**3,0.0D0)/tan3**2/sin(coe3)**2-cmplx(0.25D0*dw**3,0.0D0)/sin(coe3)**4))
!       else
!         In0(i2)=coe1*cmplx(ooo-ooo*tau(i3)**2*(2.0D0*ooo**2+dw**2)/12.0D0+tau(i3)**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*tau(i3)**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,&
!               &-ooo**2*tau(i3)/2.0D0+ooo**2*tau(i3)**3*(ooo**2+dw**2)/24.0D0+tau(i3)**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
!         In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*tau(i3)**2*(ooo**2+dw**2)/8.0D0+tau(i3)**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,&
!               &-tau(i3)*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+tau(i3)**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*tau(i3)**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
!         In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-tau(i3)**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*tau(i3)**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008.0D0,&
!               &-ooo**2*tau(i3)*(ooo**2+dw**2)/4.0D0+tau(i3)**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-tau(i3)**5*(3.0D0*ooo**8+14.0D0*ooo**6*dw**2-7.0D0*ooo**4*dw**4+2.0D0*ooo**2*dw**6)/2880.0D0)
!         In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-tau(i3)**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+tau(i3)**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,&
!               &-tau(i3)*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*tau(i3)**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252.0D0)
!       end if
!     else
!       ttt=beta-tau(i3)
!       coe1=exp(cmplx(0.0D0,-fome(i2)*ttt))
!       coe2=exp(cmplx(0.0D0,-ttt*(fome(i2+1)-fome(i2))))
!       coe3=cmplx(dw*ttt/2.0D0,0.0D0)
!       tan3=sin(coe3)/cos(coe3)
!       ooo=fome(i2+1)-fome(i2)
!       if (abs((beta-tau(i3))*(fome(i2+1)-fome(i2))) .gt. 1.0D-1) then
!         In0(i2)=coe1*(cmplx(0.0D0,0.5D0*dw)/tan3*(coe2+cmplx(-1.0D0,0.0D0)))
!         In1(i2)=coe1*cmplx(0.0D0,1.0D0)*(cmplx(0.5D0*dw*ooo,0.0D0)*coe2/tan3-cmplx(0.0D0,0.25D0*dw**2)*(coe2+cmplx(-1.0D0,0.0D0))/sin(coe3)**2)
!         In2(i2)=coe1*cmplx(0.0D0,-0.5D0*dw)*(-coe2*cmplx(ooo,0.0D0)**2/tan3+cmplx(0.0D0,dw*ooo)*coe2/sin(coe3)**2+cmplx(0.5D0*dw**2,0.0D0)*(coe2+cmplx(-1.0D0,0.0D0))/tan3/sin(coe3)**2)
!         In3(i2)=coe1*cmplx(0.5D0*dw,0.0D0)*(cmplx(0.0D0,ooo**3)/tan3*coe2+cmplx(1.5D0*dw*ooo**2,0.0D0)*coe2/sin(coe3)**2&
!            &-cmplx(0.0D0,1.5D0*dw**2*ooo)*coe2/tan3/sin(coe3)**2+(coe2+cmplx(-1.0D0,0.0D0))*(-cmplx(0.5D0*dw**3,0.0D0)/tan3**2/sin(coe3)**2-cmplx(0.25D0*dw**3,0.0D0)/sin(coe3)**4))
!       else
!
!         In0(i2)=coe1*cmplx(ooo-ooo*ttt**2*(2.0D0*ooo**2+dw**2)/12.0D0+ttt**4*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/720.0D0-ooo*ttt**6*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/30240.0D0,-ooo**2*ttt/2.0D0+ooo**2*ttt**3*(ooo**2+dw**2)/24.0D0+ttt**5*(ooo**2*dw**4-2.0D0*ooo**6-5.0D0*ooo**4*dw**2)/1440.0D0)
!         In1(i2)=coe1*cmplx(ooo**2/2.0D0-ooo**2*ttt**2*(ooo**2+dw**2)/8.0D0+ttt**4*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/288.0D0,-ttt*ooo*(2.0D0*ooo**2+dw**2)/6.0D0+ttt**3*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/180.0D0-ooo*ttt**5*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/5040.0D0)
!         In2(i2)=coe1*cmplx(ooo*(2.0D0*ooo**2+dw**2)/6.0D0-ttt**2*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/60.0D0+ooo*ttt**4*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/1008D0,-ooo**2*ttt*(ooo**2+dw**2)/4.0D0+ttt**3*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/72.0D0-ttt**5*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/2880D0)
!         In3(i2)=coe1*cmplx(ooo**2*(ooo**2+dw**2)/4.0D0-ttt**2*(2.0D0*ooo**6+5.0D0*ooo**4*dw**2-ooo**2*dw**4)/24.0D0+ttt**4*(3D0*ooo**8+14D0*ooo**6*dw**2-7D0*ooo**4*dw**4+2D0*ooo**2*dw**6)/576D0,-ttt*(6.0D0*ooo**5+10.0D0*ooo**3*dw**2-ooo*dw**4)/30.0D0+ooo*ttt**3*(6D0*ooo**6+21D0*ooo**4*dw**2-7D0*ooo**2*dw**4+dw**6)/252D0)
!        endif
!       endif
!    end do
!    f1(i3)=dot_product(a,In0)+dot_product(b,In1)+dot_product(c,In2)+dot_product(d,In3)
!
!    outmatrix(i1,i3,1)=outmatrix(i1,i3,1)+(2.0D0*real(f1(i3)))/(2.0D0*pi)                    !n from -inf to inf
!    outmatrix(i1,i3,2)=0.0D0!aimag(f1(i3))
! end do
!
! end do
! !---------------------------------------------------
! end subroutine SFT_omegaFtotau
!

! subroutine SFT_omegaFtotau(inmatrix,outmatrix,beta)
! use vec
! implicit none
! real*8,intent(in):: inmatrix(n_,n_,2)
! real*8,intent(out):: outmatrix(n_,n_,2)
! real*8 xx(2*n_),yy1(2*n_),yy2(2*n_),yy21(2*n_),yy22(2*n_)
! real*8 a,b,c,d,e,f,midmatrix(n_,n_,8),c1,midmatrix1(n_,n_,2),aa,bb,cc,dd,ee,ff
! complex*16 c2,step,c3,c4,ii,fm
! integer i1,i2,i3,nu,nnn(n_),nm
! integer*8 nn
! real*8 mc0,mc1,ms0,ms1,rr,tf0,tf1,tf2,tf3,tf4,beta
! real*8 dx,kkk,xi,xj,xxx
! complex*16 co1,co2,co3,co4,ccc1,ccc2,ccc3,tan3
!
! outmatrix = 0.0D0
!
! ! FT at omega_n
!
!    xx(1:n_) = -fome(n_:1:-1)
!
!    xx((n_+1):2*n_) = fome(1:n_)
!
!    do i1 = 1,n_
!    !give the spline number
!    yy1(1:n_) = inmatrix(i1,n_:1:-1,1)
!    !yy1(42) = inmatrix(i1,42,1)
!    yy1(n_+1:2*n_) = inmatrix(i1,1:n_,1)
!    yy2(1:n_) = -inmatrix(i1,n_:1:-1,2)
!   ! yy2(42) = 0.0D0
!    yy2(n_+1:2*n_) = inmatrix(i1,1:n_,2)
!
!   ! yy2(1:301) = inmatrix(i1,1:301,2)
!    call spline(xx,yy1,2*n_,yy21)
!    call spline(xx,yy2,2*n_,yy22)
!    !give the interpolation function
!
!      do i2 = n_+1,2*n_-1
!    !here i2 = j-1
!        midmatrix(i1,i2-n_,1) = yy1(i2)
!        midmatrix(i1,i2-n_,2) = xx(i2)*yy21(i2)/3.0D0-xx(i2+1)*yy21(i2)/3.0D0+xx(i2)*yy21(i2+1)/6.0D0-xx(i2+1)*yy21(i2+1)/6.0D0+yy1(i2)/(xx(i2)-xx(i2+1))+yy1(i2+1)/(-xx(i2)+xx(i2+1))
!        midmatrix(i1,i2-n_,3) = yy21(i2)/2.0D0
!        midmatrix(i1,i2-n_,4) = (yy21(i2)-yy21(i2+1))/(6.0D0*(xx(i2)-xx(i2+1)))
!        midmatrix(i1,i2-n_,5) = yy2(i2)
!        midmatrix(i1,i2-n_,6) = xx(i2)*yy22(i2)/3.0D0-xx(i2+1)*yy22(i2)/3.0D0+xx(i2)*yy22(i2+1)/6.0D0-xx(i2+1)*yy22(i2+1)/6.0D0+yy2(i2)/(xx(i2)-xx(i2+1))+yy2(i2+1)/(-xx(i2)+xx(i2+1))
!        midmatrix(i1,i2-n_,7) = yy22(i2)/2.0D0
!        midmatrix(i1,i2-n_,8) = (yy22(i2)-yy22(i2+1))/(6.0D0*(xx(i2)-xx(i2+1)))
!      enddo
!    enddo
!
!   ! do i1 = 1,301
!   !write (*,*) midmatrix(1,i1,1),midmatrix(1,i1,5),i1
!   ! enddo
!    c1 = 3.141592654D0/beta
!    dx = 2.0D0*c1
!
!    do i1 = 1,n_
!    !integral for every k
!      do i3 = 1,n_
!
!      !integral for every tau
!
!        do i2 = 1,n_/2                            !!!!!!! 18->20  19->21    !2017/4/13
!        outmatrix(i1,i3,1) = outmatrix(i1,i3,1)+real(dcmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(dcmplx(0.0D0,-fome(i2)*tau(i3))))
!        enddo
!       !  outmatrix(i1,i3,1) = outmatrix(i1,i3,1)+0.5D0*real(dcmplx(inmatrix(i1,21,1),inmatrix(i1,21,2))*exp(dcmplx(0.0D0,-fome(21)*tau(i3))))
!        do i2 = n_/2+1,n_-1
!          a = midmatrix(i1,i2,1)!+midmatrix(i1,i2,2)*c1+midmatrix(i1,i2,3)*c1**2+midmatrix(i1,i2,4)*c1**3
!          b = 2.0D0*midmatrix(i1,i2,2)*c1!+4.0D0*midmatrix(i1,i2,3)*c1**2+6.0D0*midmatrix(i1,i2,4)*c1**3
!          c = 4.0D0*midmatrix(i1,i2,3)*c1**2!+12.0D0*midmatrix(i1,i2,4)*c1**3
!          d = 8.0D0*midmatrix(i1,i2,4)*c1**3
!
!          aa = midmatrix(i1,i2,5)!+midmatrix(i1,i2,6)*c1+midmatrix(i1,i2,7)*c1**2+midmatrix(i1,i2,8)*c1**3
!          bb = 2.0D0*midmatrix(i1,i2,6)*c1!+4.0D0*midmatrix(i1,i2,7)*c1**2+6.0D0*midmatrix(i1,i2,8)*c1**3
!          cc = 4.0D0*midmatrix(i1,i2,7)*c1**2!+12.0D0*midmatrix(i1,i2,8)*c1**3
!          dd = 8.0D0*midmatrix(i1,i2,8)*c1**3
!
!          ccc1 = exp(dcmplx(0.0D0,-tau(i3)*fome(i2)))
!          ccc2 = exp(dcmplx(0.0D0,-tau(i3)*(fome(i2+1)-fome(i2))))
!          ccc3 = dcmplx(dx*tau(i3)/2.0D0,0.0D0)
!          tan3 = sin(ccc3)/cos(ccc3)
!          kkk = tau(i3)
!          if ( abs(kkk) - beta/2.0D0 .lt. 0.0D0) then
!            if (abs(kkk*(fome(i2+1)-fome(i2))) .gt. 1.0D-1) then
!               co1 = dcmplx(0.0D0,0.5D0*dx)*ccc1*(dcmplx(-1.0D0,0.0D0)+ccc2)/tan3
!               co2 = dcmplx(0.0D0,1.0D0)*ccc1*(dcmplx(0.5D0*dx*(fome(i2+1)-fome(i2)),0.0D0)*ccc2/tan3-dcmplx(0.0D0,0.25D0*dx**2)*(dcmplx(-1.0D0,0.0D0)+ccc2)/sin(ccc3)**2)
!               co3 = dcmplx(0.0D0,-0.5D0*dx)*ccc1*(-dcmplx(fome(i2+1)-fome(i2),0.0D0)**2*ccc2/tan3+dcmplx(0.0D0,dx*(fome(i2+1)-fome(i2)))*ccc2/sin(ccc3)**2+dcmplx(0.5D0*dx**2,0.0D0)*(-dcmplx(1.0D0,0.0D0)+ccc2)/tan3/sin(ccc3)**2)
!               co4 = dcmplx(0.5D0*dx,0.0D0)*ccc1*(dcmplx(0.0D0,(fome(i2+1)-fome(i2))**3)*ccc2/tan3+dcmplx(1.5D0*dx*(fome(i2+1)-fome(i2))**2,0.0D0)*ccc2/sin(ccc3)**2-dcmplx(0.0D0,1.5D0*dx**2*(fome(i2+1)-fome(i2)))*ccc2/tan3/sin(ccc3)**2+(dcmplx(-1.0D0,0.0D0)+ccc2)*(-dcmplx(0.5D0*dx**3,0.0D0)/tan3**2/sin(ccc3)**2-dcmplx(0.25D0,0.0D0)*dx**3/sin(ccc3)**4))
!             else
!               xi = fome(i2)
!               xj = fome(i2+1)
!               xxx = xi-xj
!               co1 = dcmplx(0.0D0,0.5D0)*ccc1*dcmplx(-xxx**2*kkk+kkk**3*(dx**2*xxx**2/12.0D0+xxx**4/12.0D0),2.0D0*xxx+(-dx**2*xxx/6.0D0-xxx**3/3.0D0)*kkk**2+kkk**4*(-dx**4*xxx/360.0D0+dx**2*xxx**3/36.0D0+xxx**5/60.0D0))
!               co2 = dcmplx(0.0D0,1.0D0)*ccc1*dcmplx(kkk*(dx**2*xi+2.0D0*xi**3-dx**2*xj-6.0D0*xi**2*xj+6.0D0*xi*xj**2-2.0D0*xj**3)/6.0D0,-0.5D0*xxx**2+kkk**2*(dx**2*xi**2+xi**4-2.0D0*dx**2*xi*xj-4.0D0*xi**3*xj+dx**2*xj**2+6.0D0*xi**2*xj**2-4.0D0*xi*xj**3+xj**4)/8.0D0)
!               co3 = dcmplx(0.0D0,-0.5D0)*ccc1*dcmplx((dx**2*xi**2+xi**4-2.0D0*dx**2*xi*xj-4.0D0*xi**3*xj+dx**2*xj**2+6.0D0*xi**2*xj**2-4.0D0*xi*xj**3+xj**4)*kkk/2.0D0,-(dx**2*xi+2.0D0*xi**3-dx**2*xj-6.0D0*xi**2*xj+6.0D0*xi*xj**2-2.0D0*xj**3)/3.0D0-(dx**4*xi-10.0D0*dx**2*xi**3-6.0D0*xi**5-dx**4*xj+30.0D0*dx**2*xi**2*xj+30.0D0*xi**4*xj-30.0D0*dx**2*xi*xj**2-60.0D0*xi**3*xj**2+10.0D0*dx**2*xj**3+60.0D0*xi**2*xj**3-30.0D0*xi*xj**4+6.0D0*xj**5)*kkk**2/30.0D0)
!               co4 = dcmplx(0.5D0,0.0D0)*ccc1*dcmplx((xxx**2*(dx**2+xi**2-2.0D0*xi*xj+xj**2))/2.0D0,-(dx**4*xi-10.0D0*dx**2*xi**3-6.0D0*xi**5-dx**4*xj+30.0D0*dx**2*xi**2*xj+30.0D0*xi**4*xj-30.0D0*dx**2*xi*xj**2-60.0D0*xi**3*xj**2+10.0D0*dx**2*xj**3+60.0D0*xi**2*xj**3-30.0D0*xi*xj**4+6.0D0*xj**5)*kkk/15.0D0)
!            endif
!           else
!            if (abs((abs(kkk)-beta)*(fome(i2+1)-fome(i2))) .gt. 1.0D-1) then
!               co1 = dcmplx(0.0D0,0.5D0*dx)*ccc1*(dcmplx(-1.0D0,0.0D0)+ccc2)/tan3
!               co2 = dcmplx(0.0D0,1.0D0)*ccc1*(dcmplx(0.5D0*dx*(fome(i2+1)-fome(i2)),0.0D0)*ccc2/tan3-dcmplx(0.0D0,0.25D0*dx**2)*(dcmplx(-1.0D0,0.0D0)+ccc2)/sin(ccc3)**2)
!               co3 = dcmplx(0.0D0,-0.5D0*dx)*ccc1*(-dcmplx(fome(i2+1)-fome(i2),0.0D0)**2*ccc2/tan3+dcmplx(0.0D0,dx*(fome(i2+1)-fome(i2)))*ccc2/sin(ccc3)**2+dcmplx(0.5D0*dx**2,0.0D0)*(-dcmplx(1.0D0,0.0D0)+ccc2)/tan3/sin(ccc3)**2)
!               co4 = dcmplx(0.5D0*dx,0.0D0)*ccc1*(dcmplx(0.0D0,(fome(i2+1)-fome(i2))**3)*ccc2/tan3+dcmplx(1.5D0*dx*(fome(i2+1)-fome(i2))**2,0.0D0)*ccc2/sin(ccc3)**2-dcmplx(0.0D0,1.5D0*dx**2*(fome(i2+1)-fome(i2)))*ccc2/tan3/sin(ccc3)**2+(dcmplx(-1.0D0,0.0D0)+ccc2)*(-dcmplx(0.5D0*dx**3,0.0D0)/tan3**2/sin(ccc3)**2-dcmplx(0.25D0,0.0D0)*dx**3/sin(ccc3)**4))
!             else
!               xi = fome(i2)
!               xj = fome(i2+1)
!               xxx = xi-xj
!               if (kkk .gt. 0.0D0) then
!                 kkk = beta - kkk
!               else
!                 kkk = (-beta-kkk)
!               endif
!               co1 = dcmplx(0.0D0,0.5D0)*ccc1*dcmplx(-(-xxx**2*kkk+kkk**3*(dx**2*xxx**2/12.0D0+xxx**4/12.0D0)),2.0D0*xxx+(-dx**2*xxx/6.0D0-xxx**3/3.0D0)*kkk**2+kkk**4*(-dx**4*xxx/360.0D0+dx**2*xxx**3/36.0D0+xxx**5/60.0D0))
!               co2 = dcmplx(0.0D0,1.0D0)*ccc1*dcmplx(-(kkk*(dx**2*xi+2.0D0*xi**3-dx**2*xj-6.0D0*xi**2*xj+6.0D0*xi*xj**2-2.0D0*xj**3)/6.0D0),-0.5D0*xxx**2+kkk**2*(dx**2*xi**2+xi**4-2.0D0*dx**2*xi*xj-4.0D0*xi**3*xj+dx**2*xj**2+6.0D0*xi**2*xj**2-4.0D0*xi*xj**3+xj**4)/8.0D0)
!               co3 = dcmplx(0.0D0,-0.5D0)*ccc1*dcmplx(-((dx**2*xi**2+xi**4-2.0D0*dx**2*xi*xj-4.0D0*xi**3*xj+dx**2*xj**2+6.0D0*xi**2*xj**2-4.0D0*xi*xj**3+xj**4)*kkk/2.0D0),-(dx**2*xi+2.0D0*xi**3-dx**2*xj-6.0D0*xi**2*xj+6.0D0*xi*xj**2-2.0D0*xj**3)/3.0D0-(dx**4*xi-10.0D0*dx**2*xi**3-6.0D0*xi**5-dx**4*xj+30.0D0*dx**2*xi**2*xj+30.0D0*xi**4*xj-30.0D0*dx**2*xi*xj**2-60.0D0*xi**3*xj**2+10.0D0*dx**2*xj**3+60.0D0*xi**2*xj**3-30.0D0*xi*xj**4+6.0D0*xj**5)*kkk**2/30.0D0)
!               co4 = dcmplx(0.5D0,0.0D0)*ccc1*dcmplx((xxx**2*(dx**2+xi**2-2.0D0*xi*xj+xj**2))/2.0D0,-(-(dx**4*xi-10.0D0*dx**2*xi**3-6.0D0*xi**5-dx**4*xj+30.0D0*dx**2*xi**2*xj+30.0D0*xi**4*xj-30.0D0*dx**2*xi*xj**2-60.0D0*xi**3*xj**2+10.0D0*dx**2*xj**3+60.0D0*xi**2*xj**3-30.0D0*xi*xj**4+6.0D0*xj**5)*kkk/15.0D0))
!            endif
!           endif
!
!          step = dcmplx(a,aa)*co1+dcmplx(b,bb)*co2+dcmplx(c,cc)*co3+dcmplx(d,dd)*co4
!          outmatrix(i1,i3,1) = outmatrix(i1,i3,1)+real(step)
!          outmatrix(i1,i3,2) = 0.0D0
!
!        enddo
!          outmatrix(i1,i3,1) = 2.0D0*outmatrix(i1,i3,1) !+ real(dcmplx(inmatrix(i1,301,1),inmatrix(i1,301,2))*exp(dcmplx(0.0D0,-fome(301)*tau(i3))))
!          outmatrix(i1,i3,2) = 0.0D0
!
!      enddo
!    enddo
!   !  outmatrix(42,:,:) = 0.0D0
!   !  outmatrix(:,42,:) = 0.0D0
!    outmatrix = outmatrix/beta
! end subroutine SFT_omegaFtotau


! subroutine SFT_tautoomegaF(inmatrix,outmatrix,beta)
! !SFT_tautoomegaB transfer a function from tau-space to omega_n-space, the second B means bosion
! !E^(i omega_n tau)
! !inmatrix is i*j function value,j->tau
! !outmatrix is m*n function value,n->omega
! use vec
! implicit none
! real*8,intent(in):: inmatrix(n_,n_,2)
! real*8,intent(out):: outmatrix(n_,n_,2)
! real*8 xx(n_),yy1(n_),yy2(n_),yy21(n_),yy22(n_)
! real*8 a,b,c,d,e,f,midmatrix(n_,n_,8),c1,midmatrix1(n_,n_,2),aa,bb,cc,dd,ee,ff
! complex*16 c2,step,c3,c4,ii,fm
! integer nn,i1,i2,i3
! real*8 mc0,mc1,ms0,ms1,rr,tf0,tf1,tf2,tf3,tf4,beta
! real*8 dx
! complex*16 co1,co2,co3,co4,ccc1,ccc2,ccc3
!
! outmatrix = 0.0D0
!
! ! FT at tau
!
!    xx(1:n_) = tau(1:n_)
!    do i1 = 1,n_
!    !give the spline number
!    yy1(1:n_) = inmatrix(i1,1:n_,1)
!    yy2(1:n_) = inmatrix(i1,1:n_,2)
!    call spline(xx,yy1,n_,yy21)
!    call spline(xx,yy2,n_,yy22)
!    !give the interpolation function
!
!      do i2 = 1,n_-1
!    !here i2 = j-1
!        midmatrix(i1,i2,1) = yy1(i2)
!        midmatrix(i1,i2,2) = xx(i2)*yy21(i2)/3.0D0-xx(i2+1)*yy21(i2)/3.0D0+xx(i2)*yy21(i2+1)/6.0D0-xx(i2+1)*yy21(i2+1)/6.0D0+yy1(i2)/(xx(i2)-xx(i2+1))+yy1(i2+1)/(-xx(i2)+xx(i2+1))
!        midmatrix(i1,i2,3) = yy21(i2)/2.0D0
!        midmatrix(i1,i2,4) = (yy21(i2)-yy21(i2+1))/(6.0D0*(xx(i2)-xx(i2+1)))
!        midmatrix(i1,i2,5) = yy2(i2)
!        midmatrix(i1,i2,6) = xx(i2)*yy22(i2)/3.0D0-xx(i2+1)*yy22(i2)/3.0D0+xx(i2)*yy22(i2+1)/6.0D0-xx(i2+1)*yy22(i2+1)/6.0D0+yy2(i2)/(xx(i2)-xx(i2+1))+yy2(i2+1)/(-xx(i2)+xx(i2+1))
!        midmatrix(i1,i2,7) = yy22(i2)/2.0D0
!        midmatrix(i1,i2,8) = (yy22(i2)-yy22(i2+1))/(6.0D0*(xx(i2)-xx(i2+1)))
!      enddo
!    enddo
!
!
!    c1 = 3.141592654D0/beta
!    !integral
!    do i1 = 1,n_
!    !integral for every k
!      do i3 = 1,n_
!      !integral for every omega_n
!          do i2 = 1,n_-1
!            c2 =  dcmplx(0.0D0,fome(i3))
!            c4 =  dcmplx(0.0D0,fome(i3)*tau(i2))
! !        nn = nint((fome(i2+1)-fome(i2))/c1)/2-1
!            ii = dcmplx(tau(i2+1)-tau(i2),0.0D0)
!            step = dcmplx(0.0D0,0.0D0)
!            fm = c2*ii
!
!            a = midmatrix(i1,i2,1)!+midmatrix(i1,i2,2)*c1+midmatrix(i1,i2,3)*c1**2+midmatrix(i1,i2,4)*c1**3
!            b = midmatrix(i1,i2,2)!*c1!+4.0D0*midmatrix(i1,i2,3)*c1**2+6.0D0*midmatrix(i1,i2,4)*c1**3
!            c = midmatrix(i1,i2,3)!*c1**2!+12.0D0*midmatrix(i1,i2,4)*c1**3
!            d = midmatrix(i1,i2,4)!*c1**3
!
!            aa = midmatrix(i1,i2,5)!+midmatrix(i1,i2,6)*c1+midmatrix(i1,i2,7)*c1**2+midmatrix(i1,i2,8)*c1**3
!            bb = midmatrix(i1,i2,6)!*c1!+4.0D0*midmatrix(i1,i2,7)*c1**2+6.0D0*midmatrix(i1,i2,8)*c1**3
!            cc = midmatrix(i1,i2,7)!*c1**2!+12.0D0*midmatrix(i1,i2,8)*c1**3
!            dd = midmatrix(i1,i2,8)!*c1**3
!
!            ccc1 = exp(dcmplx(0.0D0,-fome(i3)*tau(i2)))
!            ccc2 = exp(dcmplx(0.0D0,-fome(i3)*(tau(i2+1)-tau(i2))))
!            ccc3 = dcmplx(dx*tau(i3)/2.0D0,0.0D0)
!            co1 = dcmplx(0.0D0,1.0D0/fome(i3))*ccc1*(dcmplx(-1.0D0,0.0D0)+ccc2)
!            co2 = dcmplx(0.0D0,1.0D0)*ccc1*(dcmplx(0.0D0,-1.0D0/fome(i3)**2)*(dcmplx(-1.0D0,0.0D0)+ccc2)+dcmplx((tau(i2+1)-tau(i2))/fome(i3),0.0D0)*ccc2)
!            co3 = dcmplx(-1.0D0,0.0D0)*ccc1*(dcmplx(0.0D0,2.0D0/fome(i3)**3)*(dcmplx(-1.0D0,0.0D0)+ccc2)-dcmplx(2.0D0*(tau(i2+1)-tau(i2))/fome(i3)**2,0.0D0)*ccc2-dcmplx(0.0D0,(tau(i2+1)-tau(i2))**2/fome(i3))*ccc2)
!            co4 = dcmplx(0.0D0,1.0D0)*ccc1*(dcmplx(0.0D0,-6.0D0/fome(i3)**4)*(dcmplx(-1.0D0,0.0D0)+ccc2)+dcmplx(6.0D0*(tau(i2+1)-tau(i2))/fome(i3)**3,0.0D0)*ccc2+dcmplx(0.0D0,3.0D0*(tau(i2+1)-tau(i2))**2/fome(i3)**2)*ccc2-dcmplx((tau(i2+1)-tau(i2))**3/fome(i3),0.0D0)*ccc2)
!            step = dcmplx(a,aa)*co1+dcmplx(b,bb)*co2+dcmplx(c,cc)*co3+dcmplx(d,dd)*co4
!            outmatrix(i1,i3,1) = outmatrix(i1,i3,1)+real(step)
!            outmatrix(i1,i3,2) = outmatrix(i1,i3,2)+aimag(step)
!          enddo
!      enddo
!    enddo
!    outmatrix = outmatrix
! end subroutine SFT_tautoomegaF


! subroutine SFT_tautoomegaB(inmatrix,outmatrix,beta)
! !SFT_tautoomegaB transfer a function from tau-space to omega_n-space, the second B means bosion
! !E^(i omega_n tau)
! !inmatrix is i*j function value,j->tau
! !outmatrix is m*n function value,n->omega
! use vec
! implicit none
! real*8,intent(in):: inmatrix(n_,n_,2)
! real*8,intent(out):: outmatrix(n_,n_,2)
! real*8 xx(n_),yy1(n_),yy2(n_),yy21(n_),yy22(n_)
! real*8 a,b,c,d,e,f,midmatrix(n_,n_,8),c1,midmatrix1(n_,n_,2),aa,bb,cc,dd,ee,ff
! complex*16 c2,step,c3,c4,ii,fm
! integer nn,i1,i2,i3
! real*8 mc0,mc1,ms0,ms1,rr,tf0,tf1,tf2,tf3,tf4,beta
! real*8 dx,kkk,xi,xj,xxx
! complex*16 co1,co2,co3,co4,ccc1,ccc2,ccc3
!
! outmatrix = 0.0D0
!
! ! FT at tau
!
!    xx(1:n_) = tau(1:n_)
!    do i1 = 1,n_
!    !give the spline number
!    yy1(1:n_) = inmatrix(i1,1:n_,1)
!    yy2(1:n_) = inmatrix(i1,1:n_,2)
!    call spline(xx,yy1,n_,yy21)
!    call spline(xx,yy2,n_,yy22)
!    !give the interpolation function
!
!      do i2 = 1,n_-1
!    !here i2 = j-1
!        midmatrix(i1,i2,1) = yy1(i2)
!        midmatrix(i1,i2,2) = xx(i2)*yy21(i2)/3.0D0-xx(i2+1)*yy21(i2)/3.0D0+xx(i2)*yy21(i2+1)/6.0D0-xx(i2+1)*yy21(i2+1)/6.0D0+yy1(i2)/(xx(i2)-xx(i2+1))+yy1(i2+1)/(-xx(i2)+xx(i2+1))
!        midmatrix(i1,i2,3) = yy21(i2)/2.0D0
!        midmatrix(i1,i2,4) = (yy21(i2)-yy21(i2+1))/(6.0D0*(xx(i2)-xx(i2+1)))
!        midmatrix(i1,i2,5) = yy2(i2)
!        midmatrix(i1,i2,6) = xx(i2)*yy22(i2)/3.0D0-xx(i2+1)*yy22(i2)/3.0D0+xx(i2)*yy22(i2+1)/6.0D0-xx(i2+1)*yy22(i2+1)/6.0D0+yy2(i2)/(xx(i2)-xx(i2+1))+yy2(i2+1)/(-xx(i2)+xx(i2+1))
!        midmatrix(i1,i2,7) = yy22(i2)/2.0D0
!        midmatrix(i1,i2,8) = (yy22(i2)-yy22(i2+1))/(6.0D0*(xx(i2)-xx(i2+1)))
!      enddo
!    enddo
!
!
!    c1 = 3.141592654D0/beta
!    !integral
!    do i1 = 1,n_
!    !integral for every k
!      do i3 = 2,n_
!      !integral for every omega_n
!          do i2 = 1,n_-1
!            c2 =  dcmplx(0.0D0,bome(i3))
!            c4 =  dcmplx(0.0D0,bome(i3)*tau(i2))
! !        nn = nint((bome(i2+1)-bome(i2))/c1)/2-1
!            ii = dcmplx(tau(i2+1)-tau(i2),0.0D0)
!            step = dcmplx(0.0D0,0.0D0)
!            fm = c2*ii
!
!            a = midmatrix(i1,i2,1)!+midmatrix(i1,i2,2)*c1+midmatrix(i1,i2,3)*c1**2+midmatrix(i1,i2,4)*c1**3
!            b = midmatrix(i1,i2,2)!*c1!+4.0D0*midmatrix(i1,i2,3)*c1**2+6.0D0*midmatrix(i1,i2,4)*c1**3
!            c = midmatrix(i1,i2,3)!*c1**2!+12.0D0*midmatrix(i1,i2,4)*c1**3
!            d = midmatrix(i1,i2,4)!*c1**3
!
!            aa = midmatrix(i1,i2,5)!+midmatrix(i1,i2,6)*c1+midmatrix(i1,i2,7)*c1**2+midmatrix(i1,i2,8)*c1**3
!            bb = midmatrix(i1,i2,6)!*c1!+4.0D0*midmatrix(i1,i2,7)*c1**2+6.0D0*midmatrix(i1,i2,8)*c1**3
!            cc = midmatrix(i1,i2,7)!*c1**2!+12.0D0*midmatrix(i1,i2,8)*c1**3
!            dd = midmatrix(i1,i2,8)!*c1**3
!
!            ccc1 = exp(dcmplx(0.0D0,bome(i3)*tau(i2)))
!            ccc2 = exp(dcmplx(0.0D0,bome(i3)*(tau(i2+1)-tau(i2))))
!            ccc3 = dcmplx(dx*tau(i3)/2.0D0,0.0D0)
!            kkk = bome(i3)
!            if (abs(bome(i3)*(tau(i2+1)-tau(i2))) .gt. 1.0D-3) then
!              co1 = dcmplx(0.0D0,-1.0D0/bome(i3))*ccc1*(dcmplx(-1.0D0,0.0D0)+ccc2)
!              co2 = dcmplx(0.0D0,1.0D0)*ccc1*(dcmplx(0.0D0,1.0D0/bome(i3)**2)*(dcmplx(-1.0D0,0.0D0)+ccc2)+dcmplx((tau(i2+1)-tau(i2))/bome(i3),0.0D0)*ccc2)
!              co3 = dcmplx(-1.0D0,0.0D0)*ccc1*(dcmplx(0.0D0,-2.0D0/bome(i3)**3)*(dcmplx(-1.0D0,0.0D0)+ccc2)-dcmplx(2.0D0*(tau(i2+1)-tau(i2))/bome(i3)**2,0.0D0)*ccc2+dcmplx(0.0D0,(tau(i2+1)-tau(i2))**2/bome(i3))*ccc2)
!              co4 = dcmplx(0.0D0,-1.0D0)*ccc1*(dcmplx(0.0D0,6.0D0/bome(i3)**4)*(dcmplx(-1.0D0,0.0D0)+ccc2)+dcmplx(6.0D0*(tau(i2+1)-tau(i2))/bome(i3)**3,0.0D0)*ccc2-dcmplx(0.0D0,3.0D0*(tau(i2+1)-tau(i2))**2/bome(i3)**2)*ccc2-dcmplx((tau(i2+1)-tau(i2))**3/bome(i3),0.0D0)*ccc2)
!            else
!              xi = tau(i2)
!              xj = tau(i2+1)
!              xxx = xj-xi
!              co1 = dcmplx(0.0D0,-1.0D0)*ccc1*dcmplx(-xxx**2*kkk/2.0D0+kkk**3*xxx**4/24.0D0,xxx-kkk**2*xxx**3/6.0D0)
!              co2 = dcmplx(0.0D0,1.0D0)*ccc1*dcmplx(-kkk*xxx**3/3.0D0+kkk**3*xxx**5/30.0D0,xxx**2/2.0D0-kkk**2*xxx**4/8.0D0)
!              co3 = dcmplx(-1.0D0,0.0D0)*ccc1*dcmplx(-xxx**3/3.0D0+kkk**2*xxx**5/10.0D0,-kkk*xxx**4/4.0D0+kkk**3*xxx**6/36.0D0)
!              co4 = dcmplx(0.0D0,-1.0D0)*ccc1*dcmplx(kkk*xxx**5/5.0D0-xxx**7*kkk**3/42.0D0,-xxx**4/4.0D0+kkk**2*xxx**6/12.0D0)
!            endif
!             step = dcmplx(a,aa)*co1+dcmplx(b,bb)*co2+dcmplx(c,cc)*co3+dcmplx(d,dd)*co4
!            outmatrix(i1,i3,1) = outmatrix(i1,i3,1)+real(step)
!            outmatrix(i1,i3,2) = outmatrix(i1,i3,2)+aimag(step)
!          enddo
!      enddo
!
!
!     do i2 = 1,n_-1
!       ! c2 =  dcmplx(0.0D0,bome(i3))
!       ! c4 =  dcmplx(0.0D0,bome(i3)*tau(i2))
! !      nn = nint((bome(i2+1)-bome(i2))/c1)/2-1
!        ii = dcmplx(tau(i2+1)-tau(i2),0.0D0)
!        step = dcmplx(0.0D0,0.0D0)
!
!        a = midmatrix(i1,i2,1)!+midmatrix(i1,i2,2)*c1+midmatrix(i1,i2,3)*c1**2+midmatrix(i1,i2,4)*c1**3
!        b = midmatrix(i1,i2,2)!*c1!+4.0D0*midmatrix(i1,i2,3)*c1**2+6.0D0*midmatrix(i1,i2,4)*c1**3
!        c = midmatrix(i1,i2,3)!*c1**2!+12.0D0*midmatrix(i1,i2,4)*c1**3
!        d = midmatrix(i1,i2,4)!*c1**3
!
!        aa = midmatrix(i1,i2,5)!+midmatrix(i1,i2,6)*c1+midmatrix(i1,i2,7)*c1**2+midmatrix(i1,i2,8)*c1**3
!        bb = midmatrix(i1,i2,6)!*c1!+4.0D0*midmatrix(i1,i2,7)*c1**2+6.0D0*midmatrix(i1,i2,8)*c1**3
!        cc = midmatrix(i1,i2,7)!*c1**2!+12.0D0*midmatrix(i1,i2,8)*c1**3
!        dd = midmatrix(i1,i2,8)!*c1**3
!
!        step = step+dcmplx(a,aa)*(ii)
!        step = step+dcmplx(b,bb)*((ii**2)/dcmplx(2.0D0,0.0D0))
!        step = step+dcmplx(c,cc)*((ii**3)/dcmplx(3.0D0,0.0D0))
!        step = step+dcmplx(d,dd)*((ii**4)/dcmplx(4.0D0,0.0D0))
!
!
!        outmatrix(i1,1,1) = outmatrix(i1,1,1) + real(step)
!        outmatrix(i1,1,2) = outmatrix(i1,1,2) + aimag(step)
!      enddo
!
!    enddo
!    outmatrix = outmatrix
!
! end subroutine SFT_tautoomegaB

! subroutine SFT_omegaBtotau(inmatrix,outmatrix,beta)
! use vec
! implicit none
! real*8,intent(in):: inmatrix(n_,n_,2)
! real*8,intent(out):: outmatrix(n_,n_,2)
! real*8 xx(2*n_),yy1(2*n_),yy2(2*n_),yy21(2*n_),yy22(2*n_)
! real*8 a,b,c,d,e,f,midmatrix(n_,n_,8),c1,midmatrix1(n_,n_,2),aa,bb,cc,dd,ee,ff
! complex*16 c2,step,c3,c4,ii,fm
! integer i1,i2,i3
! real*8 rr,beta
! real*8 dx,xj,xi,xxx,kkk
! complex*16 co1,co2,co3,co4,ccc1,ccc2,ccc3,tan3
!
! outmatrix = 0.0D0
!
! ! FT at omega_n
!
!    xx(1:n_) = -bome(n_:1:-1)
!  !  xx(42) = bome(42)
!    xx(n_+1:2*n_) = bome(1:n_)
!
!    do i1 = 1,n_
!    !give the spline number
!    yy1(1:n_) = inmatrix(i1,n_:1:-1,1)
!   ! yy1(42) = inmatrix(i1,42,1)
!    yy1(n_+1:2*n_) = inmatrix(i1,1:n_,1)
!    yy2(1:n_) = -inmatrix(i1,n_:1:-1,2)
!  !  yy2(42) = 0.0D0
!    yy2(n_+1:2*n_) = inmatrix(i1,1:n_,2)
!   ! yy2(1:301) = inmatrix(i1,1:301,2)
!    yy21 = 0.0D0
!    yy22 = 0.0D0
!    call spline(xx,yy1,2*n_,yy21)
!    call spline(xx,yy2,2*n_,yy22)
!    !give the interpolation function
!
!      do i2 = n_+1,2*n_-1
!    !here i2 = j-1
!        midmatrix(i1,i2-n_,1) = yy1(i2)
!        midmatrix(i1,i2-n_,2) = xx(i2)*yy21(i2)/3.0D0-xx(i2+1)*yy21(i2)/3.0D0+xx(i2)*yy21(i2+1)/6.0D0-xx(i2+1)*yy21(i2+1)/6.0D0+yy1(i2)/(xx(i2)-xx(i2+1))+yy1(i2+1)/(-xx(i2)+xx(i2+1))
!        midmatrix(i1,i2-n_,3) = yy21(i2)/2.0D0
!        midmatrix(i1,i2-n_,4) = (yy21(i2)-yy21(i2+1))/(6.0D0*(xx(i2)-xx(i2+1)))
!        midmatrix(i1,i2-n_,5) = yy2(i2)
!        midmatrix(i1,i2-n_,6) = xx(i2)*yy22(i2)/3.0D0-xx(i2+1)*yy22(i2)/3.0D0+xx(i2)*yy22(i2+1)/6.0D0-xx(i2+1)*yy22(i2+1)/6.0D0+yy2(i2)/(xx(i2)-xx(i2+1))+yy2(i2+1)/(-xx(i2)+xx(i2+1))
!        midmatrix(i1,i2-n_,7) = yy22(i2)/2.0D0
!        midmatrix(i1,i2-n_,8) = (yy22(i2)-yy22(i2+1))/(6.0D0*(xx(i2)-xx(i2+1)))
!      enddo
!    enddo
!
!   ! do i1 = 1,301
!   !write (*,*) midmatrix(1,i1,1),midmatrix(1,i1,5),i1
!   ! enddo
!    c1 = 3.141592654D0/beta
!    dx = 2.0D0*c1
!
!    do i1 = 1,n_
!    !integral for every k
!      do i3 = 1,n_
!      !integral for every tau
!           do i2 = 1,(n_/2)                               !!!!!!! 18->20  19->21    !2017/4/13
!             outmatrix(i1,i3,1) = outmatrix(i1,i3,1)+real(dcmplx(inmatrix(i1,i2,1),inmatrix(i1,i2,2))*exp(dcmplx(0.0D0,-bome(i2)*tau(i3))))
!           enddo
!             ! outmatrix(i1,i3,1) = outmatrix(i1,i3,1)+0.5D0*real(dcmplx(inmatrix(i1,21,1),inmatrix(i1,21,2))*exp(dcmplx(0.0D0,-bome(21)*tau(i3))))
!          do i2 = (n_/2+1),n_-1
!          a = midmatrix(i1,i2,1)!+midmatrix(i1,i2,2)*c1+midmatrix(i1,i2,3)*c1**2+midmatrix(i1,i2,4)*c1**3
!          b = 2.0D0*midmatrix(i1,i2,2)*c1!+4.0D0*midmatrix(i1,i2,3)*c1**2+6.0D0*midmatrix(i1,i2,4)*c1**3
!          c = 4.0D0*midmatrix(i1,i2,3)*c1**2!+12.0D0*midmatrix(i1,i2,4)*c1**3
!          d = 8.0D0*midmatrix(i1,i2,4)*c1**3
!
!          aa = midmatrix(i1,i2,5)!+midmatrix(i1,i2,6)*c1+midmatrix(i1,i2,7)*c1**2+midmatrix(i1,i2,8)*c1**3
!          bb = 2.0D0*midmatrix(i1,i2,6)*c1!+4.0D0*midmatrix(i1,i2,7)*c1**2+6.0D0*midmatrix(i1,i2,8)*c1**3
!          cc = 4.0D0*midmatrix(i1,i2,7)*c1**2!+12.0D0*midmatrix(i1,i2,8)*c1**3
!          dd = 8.0D0*midmatrix(i1,i2,8)*c1**3
!       !   write (*,*) a,b,c,d,aa,bb,cc,dd
!       !   pause
!          ccc1 = exp(dcmplx(0.0D0,-tau(i3)*bome(i2)))
!          ccc2 = exp(dcmplx(0.0D0,-tau(i3)*(bome(i2+1)-bome(i2))))
!          ccc3 = dcmplx(dx*tau(i3)/2.0D0,0.0D0)
!          tan3 = sin(ccc3)/cos(ccc3)
!          kkk = tau(i3)
!          co1 = dcmplx(0.0D0,0.0D0)
!          co2 = dcmplx(0.0D0,0.0D0)
!          co3 = dcmplx(0.0D0,0.0D0)
!          co4 = dcmplx(0.0D0,0.0D0)
!          if ( (kkk - beta/2.0D0) .le. 0.0D0) then
!            if (abs(kkk*(bome(i2+1)-bome(i2))) .gt. 1.0D-1) then
!               co1 = dcmplx(0.0D0,0.5D0*dx)*ccc1*(dcmplx(-1.0D0,0.0D0)+ccc2)/tan3
!               co2 = dcmplx(0.0D0,1.0D0)*ccc1*(dcmplx(0.5D0*dx*(bome(i2+1)-bome(i2)),0.0D0)*ccc2/tan3-dcmplx(0.0D0,0.25D0*dx**2)*(dcmplx(-1.0D0,0.0D0)+ccc2)/sin(ccc3)**2)
!               co3 = dcmplx(0.0D0,-0.5D0*dx)*ccc1*(-dcmplx(bome(i2+1)-bome(i2),0.0D0)**2*ccc2/tan3+dcmplx(0.0D0,dx*(bome(i2+1)-bome(i2)))*ccc2/sin(ccc3)**2+dcmplx(0.5D0*dx**2,0.0D0)*(-dcmplx(1.0D0,0.0D0)+ccc2)/tan3/sin(ccc3)**2)
!               co4 = dcmplx(0.5D0*dx,0.0D0)*ccc1*(dcmplx(0.0D0,(bome(i2+1)-bome(i2))**3)*ccc2/tan3+dcmplx(1.5D0*dx*(bome(i2+1)-bome(i2))**2,0.0D0)*ccc2/sin(ccc3)**2-dcmplx(0.0D0,1.5D0*dx**2*(bome(i2+1)-bome(i2)))*ccc2/tan3/sin(ccc3)**2+(dcmplx(-1.0D0,0.0D0)+ccc2)*(-dcmplx(0.5D0*dx**3,0.0D0)/tan3**2/sin(ccc3)**2-dcmplx(0.25D0,0.0D0)*dx**3/sin(ccc3)**4))
!             else
!               xi = bome(i2)
!               xj = bome(i2+1)
!               xxx = xi-xj
!               co1 = dcmplx(0.0D0,0.5D0)*ccc1*dcmplx(-xxx**2*kkk+kkk**3*(dx**2*xxx**2/12.0D0+xxx**4/12.0D0),2.0D0*xxx+(-dx**2*xxx/6.0D0-xxx**3/3.0D0)*kkk**2+kkk**4*(-dx**4*xxx/360.0D0+dx**2*xxx**3/36.0D0+xxx**5/60.0D0))
!               co2 = dcmplx(0.0D0,1.0D0)*ccc1*dcmplx(kkk*(dx**2*xi+2.0D0*xi**3-dx**2*xj-6.0D0*xi**2*xj+6.0D0*xi*xj**2-2.0D0*xj**3)/6.0D0,-0.5D0*xxx**2+kkk**2*(dx**2*xi**2+xi**4-2.0D0*dx**2*xi*xj-4.0D0*xi**3*xj+dx**2*xj**2+6.0D0*xi**2*xj**2-4.0D0*xi*xj**3+xj**4)/8.0D0)
!               co3 = dcmplx(0.0D0,-0.5D0)*ccc1*dcmplx((dx**2*xi**2+xi**4-2.0D0*dx**2*xi*xj-4.0D0*xi**3*xj+dx**2*xj**2+6.0D0*xi**2*xj**2-4.0D0*xi*xj**3+xj**4)*kkk/2.0D0,-(dx**2*xi+2.0D0*xi**3-dx**2*xj-6.0D0*xi**2*xj+6.0D0*xi*xj**2-2.0D0*xj**3)/3.0D0-(dx**4*xi-10.0D0*dx**2*xi**3-6.0D0*xi**5-dx**4*xj+30.0D0*dx**2*xi**2*xj+30.0D0*xi**4*xj-30.0D0*dx**2*xi*xj**2-60.0D0*xi**3*xj**2+10.0D0*dx**2*xj**3+60.0D0*xi**2*xj**3-30.0D0*xi*xj**4+6.0D0*xj**5)*kkk**2/30.0D0)
!               co4 = dcmplx(0.5D0,0.0D0)*ccc1*dcmplx((xxx**2*(dx**2+xi**2-2.0D0*xi*xj+xj**2))/2.0D0,-(dx**4*xi-10.0D0*dx**2*xi**3-6.0D0*xi**5-dx**4*xj+30.0D0*dx**2*xi**2*xj+30.0D0*xi**4*xj-30.0D0*dx**2*xi*xj**2-60.0D0*xi**3*xj**2+10.0D0*dx**2*xj**3+60.0D0*xi**2*xj**3-30.0D0*xi*xj**4+6.0D0*xj**5)*kkk/15.0D0)
!            endif
!           else
!            if (abs((kkk-beta)*(bome(i2+1)-bome(i2))) .gt. 1.0D-1) then
!               co1 = dcmplx(0.0D0,0.5D0*dx)*ccc1*(dcmplx(-1.0D0,0.0D0)+ccc2)/tan3
!               co2 = dcmplx(0.0D0,1.0D0)*ccc1*(dcmplx(0.5D0*dx*(bome(i2+1)-bome(i2)),0.0D0)*ccc2/tan3-dcmplx(0.0D0,0.25D0*dx**2)*(dcmplx(-1.0D0,0.0D0)+ccc2)/sin(ccc3)**2)
!               co3 = dcmplx(0.0D0,-0.5D0*dx)*ccc1*(-dcmplx(bome(i2+1)-bome(i2),0.0D0)**2*ccc2/tan3+dcmplx(0.0D0,dx*(bome(i2+1)-bome(i2)))*ccc2/sin(ccc3)**2+dcmplx(0.5D0*dx**2,0.0D0)*(-dcmplx(1.0D0,0.0D0)+ccc2)/tan3/sin(ccc3)**2)
!               co4 = dcmplx(0.5D0*dx,0.0D0)*ccc1*(dcmplx(0.0D0,(bome(i2+1)-bome(i2))**3)*ccc2/tan3+dcmplx(1.5D0*dx*(bome(i2+1)-bome(i2))**2,0.0D0)*ccc2/sin(ccc3)**2-dcmplx(0.0D0,1.5D0*dx**2*(bome(i2+1)-bome(i2)))*ccc2/tan3/sin(ccc3)**2+(dcmplx(-1.0D0,0.0D0)+ccc2)*(-dcmplx(0.5D0*dx**3,0.0D0)/tan3**2/sin(ccc3)**2-dcmplx(0.25D0,0.0D0)*dx**3/sin(ccc3)**4))
!             else
!               xi = bome(i2)
!               xj = bome(i2+1)
!               xxx = xi-xj
!               kkk = beta - kkk
!               co1 = dcmplx(0.0D0,0.5D0)*ccc1*dcmplx(-(-xxx**2*kkk+kkk**3*(dx**2*xxx**2/12.0D0+xxx**4/12.0D0)),2.0D0*xxx+(-dx**2*xxx/6.0D0-xxx**3/3.0D0)*kkk**2+kkk**4*(-dx**4*xxx/360.0D0+dx**2*xxx**3/36.0D0+xxx**5/60.0D0))
!               co2 = dcmplx(0.0D0,1.0D0)*ccc1*dcmplx(-(kkk*(dx**2*xi+2.0D0*xi**3-dx**2*xj-6.0D0*xi**2*xj+6.0D0*xi*xj**2-2.0D0*xj**3)/6.0D0),-0.5D0*xxx**2+kkk**2*(dx**2*xi**2+xi**4-2.0D0*dx**2*xi*xj-4.0D0*xi**3*xj+dx**2*xj**2+6.0D0*xi**2*xj**2-4.0D0*xi*xj**3+xj**4)/8.0D0)
!               co3 = dcmplx(0.0D0,-0.5D0)*ccc1*dcmplx(-((dx**2*xi**2+xi**4-2.0D0*dx**2*xi*xj-4.0D0*xi**3*xj+dx**2*xj**2+6.0D0*xi**2*xj**2-4.0D0*xi*xj**3+xj**4)*kkk/2.0D0),-(dx**2*xi+2.0D0*xi**3-dx**2*xj-6.0D0*xi**2*xj+6.0D0*xi*xj**2-2.0D0*xj**3)/3.0D0-(dx**4*xi-10.0D0*dx**2*xi**3-6.0D0*xi**5-dx**4*xj+30.0D0*dx**2*xi**2*xj+30.0D0*xi**4*xj-30.0D0*dx**2*xi*xj**2-60.0D0*xi**3*xj**2+10.0D0*dx**2*xj**3+60.0D0*xi**2*xj**3-30.0D0*xi*xj**4+6.0D0*xj**5)*kkk**2/30.0D0)
!               co4 = dcmplx(0.5D0,0.0D0)*ccc1*dcmplx((xxx**2*(dx**2+xi**2-2.0D0*xi*xj+xj**2))/2.0D0,-(-(dx**4*xi-10.0D0*dx**2*xi**3-6.0D0*xi**5-dx**4*xj+30.0D0*dx**2*xi**2*xj+30.0D0*xi**4*xj-30.0D0*dx**2*xi*xj**2-60.0D0*xi**3*xj**2+10.0D0*dx**2*xj**3+60.0D0*xi**2*xj**3-30.0D0*xi*xj**4+6.0D0*xj**5)*kkk/15.0D0))
!            endif
!           endif
!
!          step = dcmplx(a,aa)*co1+dcmplx(b,bb)*co2+dcmplx(c,cc)*co3+dcmplx(d,dd)*co4
!          outmatrix(i1,i3,1) = outmatrix(i1,i3,1)+real(step)
!          ! write (*,*) step,i2,'1',bome(i2),co1,co2,co3,co4,a,b,c,d,ccc1,ccc2,ccc3,tan3,outmatrix(i1,i3,1)
!          ! pause
!          outmatrix(i1,i3,2) = 0.0D0
!        enddo
!          outmatrix(i1,i3,1) = 2.0D0*outmatrix(i1,i3,1) !+ real(dcmplx(inmatrix(i1,42,1),inmatrix(i1,42,2)))
!          outmatrix(i1,i3,2) = 0.0D0
!     !      write (*,*) outmatrix(1,:,1)
!      enddo
!    enddo
!
!    outmatrix = outmatrix/beta
!
! end subroutine SFT_omegaBtotau


! subroutine SFT_ktox(inmatrix,outmatrix)
! !SFT_ktox transfer a function from k-space to r-space, the second F means Fermion
! !E^(i k r)
! !inmatrix is i*j function value,i->k
! !outmatrix is m*n function value,m->r
! use vec
! implicit none
! real*8,intent(in):: inmatrix(n_,n_,2)
! real*8,intent(out):: outmatrix(n_,n_,2)
! real*8 xx(2*n_),yy1(2*n_),yy2(2*n_),yy21(2*n_),yy22(2*n_)
! real*8 a,b,c,d,e,f,midmatrix(n_,n_,8),c1,midmatrix1(n_,n_,2),aa,bb,cc,dd,ee,ff
! complex*16 c2,step,c3,c4
! integer nn,i1,i2,i3
! real*8 mc0,mc1,ms0,ms1,rr,tf0,tf1,tf2,tf3,tf4,ii
! real*8 dx,kkk,xi,xj,xxx
! complex*16 co1,co2,co3,co4,co5,ccc1,ccc2,ccc3,mmmm
! outmatrix = 0.0D0
!
! midmatrix1 = inmatrix
!    midmatrix = 0.0D0
! ! FT at k
!    xx(1:n_) = -k(n_:1:-1)
!    xx(n_+1:2*n_) = k(1:n_:1)
!
!    do i2 = 1,n_
!    !give the spline number
!    yy1(n_:1:-1) = midmatrix1(1:n_,i2,1)*k(1:n_)
!    yy1(n_+1:2*n_:1) = midmatrix1(1:n_,i2,1)*k(1:n_)
!    yy2(n_+1:2*n_:1) = midmatrix1(1:n_,i2,2)*k(1:n_)
!    yy2(n_:1:-1) = midmatrix1(1:n_,i2,2)*k(1:n_)
!
!    call spline(xx,yy1,2*n_,yy21)
!    call spline(xx,yy2,2*n_,yy22)
!
!
!    !give the interpolation function
!      do i1 = n_+1,2*n_-1
!    !here i1 = i-1
!        midmatrix(i1-n_,i2,1) = yy1(i1)
!        midmatrix(i1-n_,i2,2) = xx(i1)*yy21(i1)/3.0D0-xx(i1+1)*yy21(i1)/3.0D0+xx(i1)*yy21(i1+1)/6.0D0-xx(i1+1)*yy21(i1+1)/6.0D0+yy1(i1)/(xx(i1)-xx(i1+1))+yy1(i1+1)/(-xx(i1)+xx(i1+1))
!        midmatrix(i1-n_,i2,3) = yy21(i1)/2.0D0
!        midmatrix(i1-n_,i2,4) = (yy21(i1)-yy21(i1+1))/(6.0D0*(xx(i1)-xx(i1+1)))
!      !  write (*,*) yy21(i1),yy21(i1+1),xx(i1),xx(i1+1)
!       ! pause
!        midmatrix(i1-n_,i2,5) = yy2(i1)
!        midmatrix(i1-n_,i2,6) = xx(i1)*yy22(i1)/3.0D0-xx(i1+1)*yy22(i1)/3.0D0+xx(i1)*yy22(i1+1)/6.0D0-xx(i1+1)*yy22(i1+1)/6.0D0+yy2(i1)/(xx(i1)-xx(i1+1))+yy2(i1+1)/(-xx(i1)+xx(i1+1))
!        midmatrix(i1-n_,i2,7) = yy22(i1)/2.0D0
!        midmatrix(i1-n_,i2,8) = (yy22(i1)-yy22(i1+1))/(6.0D0*(xx(i1)-xx(i1+1)))
!      enddo
!    enddo
!
!    !integral
!    do i2 = 1,n_
!      !integral for every tau/omega_n
!      do i3 = 1,n_
!      !integral for every r
!      c3 = dcmplx(k(i3),0.0D0)
!
!          c2 = dcmplx(0.0D0,k(i3))
!
!          !integral k
!          do i1 = 1,n_-1
!
!          a = midmatrix(i1,i2,1)!*k(i1)
!          b = midmatrix(i1,i2,2)!+midmatrix(i1,i2,2)*k(i1)
!          c = midmatrix(i1,i2,3)!+midmatrix(i1,i2,3)*k(i1)
!          d = midmatrix(i1,i2,4)!+midmatrix(i1,i2,4)*k(i1)
!          e = midmatrix(i1,i2,4)
!
!          aa = midmatrix(i1,i2,5)*k(i1)
!          bb = midmatrix(i1,i2,5)+midmatrix(i1,i2,6)*k(i1)
!          cc = midmatrix(i1,i2,6)+midmatrix(i1,i2,7)*k(i1)
!          dd = midmatrix(i1,i2,7)+midmatrix(i1,i2,8)*k(i1)
!          ee = midmatrix(i1,i2,8)
!
! !          write (*,*) a,b,c,d,e,aa,bb,cc,dd,ee,i1,k(i1)
! !          pause
!
!          step = dcmplx(0.0D0,0.0D0)
!          ccc1 = exp(dcmplx(0.0D0,k(i3)*k(i1)))
!          ccc2 = exp(dcmplx(0.0D0,k(i3)*(k(i1+1)-k(i1))))
!         ! ccc3 = dcmplx(dx*k(i3)/2.0D0,0.0D0)
!          kkk = k(i3)
!          xi = k(i1)
!          xj = k(i1+1)
!          xxx = xj-xi
!          if ((abs(kkk*xxx) .gt. 1.0D-1)) then
!            co1 = dcmplx(0.0D0,-1.0D0/kkk)*ccc1*(dcmplx(-1.0D0,0.0D0)+ccc2)
!            co2 = dcmplx(0.0D0,-1.0D0)*ccc1*(dcmplx(0.0D0,1.0D0/kkk**2)*(dcmplx(-1.0D0,0.0D0)+ccc2)+dcmplx(xxx/kkk,0.0D0)*ccc2)
!            co3 = dcmplx(-1.0D0,0.0D0)*ccc1*(dcmplx(0.0D0,-2.0D0/kkk**3)*(dcmplx(-1.0D0,0.0D0)+ccc2)-dcmplx(2.0D0*xxx/kkk**2,0.0D0)*ccc2+dcmplx(0.0D0,xxx**2/kkk)*ccc2)
!            co4 = dcmplx(0.0D0,1.0D0)*ccc1*(dcmplx(0.0D0,6.0D0/kkk**4)*(dcmplx(-1.0D0,0.0D0)+ccc2)+dcmplx(6.0D0*xxx/kkk**3,0.0D0)*ccc2+dcmplx(0.0D0,-3.0D0*xxx**2/kkk**2)*ccc2-dcmplx(xxx**3/kkk,0.0D0)*ccc2)
!          !  co5 = dcmplx(1.0D0,0.0D0)*ccc1*(dcmplx(0.0D0,-24.0D0/kkk**5)*(dcmplx(-1.0D0,0.0D0)+ccc2)+dcmplx(-24.0D0*xxx/kkk**4,0.0D0)*ccc2+dcmplx(0.0D0,12.0D0*xxx**2/kkk**3)*ccc2+dcmplx(4.0D0*xxx**3/kkk**2,0.0D0)*ccc2-dcmplx(0.0D0,1.0D0*xxx**4/kkk)*ccc2)
!          else
!            co1 = dcmplx(0.0D0,-1.0D0)*ccc1*dcmplx(-kkk*xxx**2/2.0D0+xxx**4*kkk**3/24.0D0-xxx**6*kkk**5/720.0D0,xxx-xxx**3*kkk**2/6.0D0+xxx**5*kkk**4/120.0D0-xxx**7*kkk**6/5040.0D0)
!            co2 = dcmplx(0.0D0,-1.0D0)*ccc1*dcmplx(-xxx**3*kkk/3.0D0+xxx**5*kkk**3/30.0D0-xxx**7*kkk**5/840.0D0,xxx**2/2.0D0-xxx**4*kkk**2/8.0D0+xxx**6*kkk**4/144.0D0-xxx**8*kkk**6/5760.0D0)
!            co3 = dcmplx(-1.0D0,0.0D0)*ccc1*dcmplx(-xxx**3/3.0D0+xxx**5*kkk**2/10.0D0-xxx**7*kkk**4/168.0D0+xxx**9*kkk**6/6480.0D0,-xxx**4*kkk/4.0D0+xxx**6*kkk**3/36.0D0-xxx**8*kkk**5/960.0D0)
!            co4 = dcmplx(0.0D0,1.0D0)*ccc1*dcmplx(xxx**5*kkk/5.0D0-xxx**7*kkk**3/42.0D0+xxx**9*kkk**5/1080.0D0,-xxx**4/4.0D0+xxx**6*kkk**2/12.0D0-xxx**8*kkk**4/192.0D0+xxx**10*kkk**6/7200.0D0)
!         !   co5 = dcmplx(1.0D0,0.0D0)*ccc1*dcmplx(xxx**5/5.0D0-xxx**7*kkk**2/14.0D0+xxx**9*kkk**4/216.0D0-xxx**11*kkk**6/7920.0D0,xxx**6*kkk/6.0D0-xxx**8*kkk**3/48.0D0+xxx**10*kkk**5/1200.0D0)
!          endif
!
!
!          step = dcmplx(a,0.0D0)*co1+dcmplx(b,0.0D0)*co2+dcmplx(c,0.0D0)*co3+dcmplx(d,0.0D0)*co4!+dcmplx(e,ee)*co5
!          step = step/c3
!
!          outmatrix(i3,i2,1) = outmatrix(i3,i2,1) + aimag(step)
!          outmatrix(i3,i2,2) = 0.0D0
!          enddo
!
!       enddo
!
!     enddo
!
!     outmatrix = outmatrix*0.05066059182D0
!     outmatrix(:,:,2) = 0.0D0
!     !0.05066059182D0= 4Pi/(2Pi)^3
!     ! cut tails
!
! end subroutine SFT_ktox


program main
  use Solve_NonLin
  use formidconstant
  implicit none
  real*8 beta,mu,num
  real*8 x(1),fvec(1),diag(1)
  integer info,i

  t=0.3D0
  ! beta=1.0D0/t
  eta=0D0!-0.07D0
  ntot=2.0D5
  lambda= 0.045D0
  omega0=1.0D0/(3.0D0*lambda*ntot)**(1.0D0/3.0D0)

mu=0.5D0
do i=1,5
  write(14,*)i
  write(14,*)'t=',t
  write(13,*)i
  write(13,*)'t=',t
x(1)=mu
! do i=1,5
! call partinum(mu,num)
! fvec(1)=(num-ntot)/ntot
call hbrd(formu,1,x,fvec,1.0D-3,1.0D-3,info,diag)
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
!   do j = 1,30
!     x(j) = exp(dble(j)/dble(30)*5.0D0)-1.0D0
!   enddo
!   !!!!!!$OMP DO
!   do i=1,30
!   write(12,*)'i=',i
!   write (*,*) "i=",i
!   vext = 0.25D0*omega0**2*x(i)**2
!   call modeleq(mu,beta,vext,rho(i),u1(i),s1(i))
!   write (*,*) "x=",x(i),"rho=",rho(i),"u1=",u1(i),"s1=",s1(i)
!   end do
!   !!!!$OMP END DO
!   num=0.0D0
!   do i = 1,29
!   mid = (rho(i)*x(i)**2+rho(i+1)*x(i+1)**2)*(x(i+1)-x(i))/2.0D0
!   num = num+mid
!   enddo
!    num = -2.0D0*12.56637061D0*num/lambda
!    write(*,*)'num',num
!
!    num1=0.0D0
!    do i = 1,29
!    mid = (u1(i)*x(i)**2+u1(i+1)*x(i+1)**2)*(x(i+1)-x(i))/2.0D0
!    num1 = num1+mid
!    enddo
!     num1 = 12.56637061D0*num1/lambda/ntot!+mu
!     write(14,*)'U=',num1
!     write(*,*)'U',num1
!
!     num2=0.0D0
!     do i = 1,29
!     mid = (s1(i)*x(i)**2+s1(i+1)*x(i+1)**2)*(x(i+1)-x(i))/2.0D0
!     num2 = num2+mid
!     enddo
!      num2 = 12.56637061D0*num2/lambda/ntot
!      write(14,*)'s=',num2
!      write(*,*)'s',num2
!
!     ! r0 = 0.0D0
!     ! minl=0.0D0
!     ! maxl=3D0*sqrt(2D0*t/omega0)
!     ! n0=4
!     ! step=(maxl-minl)/dble(n0)
!     ! a=minl
!     ! b=minl+step
!     ! i1=0
!     ! do j=1,n0
!     !   answer=0
!     !   Do i=1,nk
!     !     i1=i1+1
!     !     write(12,*)'i=',i1
!     !     write (*,*) "i=",i1
!     !     answer=answer+ak1(i)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),mu,beta)
!     !   End Do
!     !   answer=answer*(b-a)/2.0D0
!     !   r0=r0+answer
!     !   a=b
!     !   b=b+step
!     ! end do
!     ! r0=r0*8.0D0*pi*sqrt(6.0D0*ntot/lambda)
!     ! num=r0(1)
!     ! write(*,*)'num',num
!     ! num1=r0(2)
!     ! write(14,*)'U=',num1
!     ! write(*,*)'U',num1
!     ! num2=r0(3)
!     ! write(14,*)'S=',num2
!     ! write(*,*)'S',num2
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
