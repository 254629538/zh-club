subroutine modeleq(mu0,beta,vext,rho,u1,s1)
use vec
use formidconstant
use ieee_arithmetic
implicit none
real*8,intent(in):: mu0,beta,vext
real*8,intent(out)::rho,u1,s1
real*8:: mu
real*8:: g0k_t(n_,n_,2),g0_x_t(n_,n_,2),Mpair(n_,n_,2),Mpair0(n_,n_,2),Mpair_1(n_,n_,2),pi_q(n_,n_,2),gamma_q(n_,n_,2),gamma0(n_,n_,2),gamma0_xt(n_,n_,2),gamma1_xt(n_,n_,2)
real*8:: g0k(n_,n_,2),g0kt(n_,n_,2),self_xt(n_,n_,2),self0(n_,n_,2),self1(n_,n_,2),g_k(n_,n_,2),g1_t(n_,n_,2),g1k_t(n_,n_,2),g1k(n_,n_,2)
real*8:: midmatrix1(n_,n_,2),midmatrix2(n_,n_,2),midmatrix3(n_,n_,2),midmatrix4(n_,n_,2),iny(n_),iw_G0(n_),lnG0(n_),g_old(n_,n_),g_bijiao(n_,n_),dU1(n_,n_,2),gamma_g(n_)
real*8:: ep,u0,u2,u3,num1,num2,s0,s2,s3,s4,s5
real*8::time1,time2
integer::i1,i2,loop,I_compare
complex*16::mid,mid1

mu=mu0-vext

!---------defination of by hand functions --------------------------------
!-----------------------------------------
!define G0(k,-tau)=g0k_t  G0(k,fome)=g0k   G0(k,tau)=g0kt
do i1=1,n_
  do i2=1,n_
    ep = k(i1)**2-mu!+vext
    mid= cmplx(1.0D0,0.0D0)/cmplx(ep,-fome(i2))
    g0k(i1,i2,1)=real(mid)
    g0k(i1,i2,2)=aimag(mid)
    g0k_t(i1,i2,1)=-1.0D0/(exp((beta-tau(i2))*ep)+exp(-ep*tau(i2)))
    g0k_t(i1,i2,2)=0.0D0
    g0kt(i1,i2,1)=exp(-ep*tau(i2))/(exp(-beta*ep)+1.0D0)
    g0kt(i1,i2,2)=0.0D0
  enddo
enddo
!---------------------------------------
! define  M(q,omrgaB)=M0+M1
 call intMpair(mu,beta,Mpair)
 write(12,*)'Mpair=',Mpair(:,1,1)
!define M0(q,omegaB) = -1/16pi*sqrt(ep-2i*omegaB)
ep = 0.0D0
mid = cmplx(0.0D0,0.0D0)
do i1 = 1,n_
  ep  = k(i1)**2-4.0D0*mu
  do i2 = 1,n_
    mid = cmplx(0.01989436789D0,0.0D0)*sqrt(cmplx(ep,-2.0D0*bome(i2)))
    Mpair0(i1,i2,1) = -real(mid)
    Mpair0(i1,i2,2) = -aimag(mid)
  enddo
enddo
write(12,*)'Mpair0=',Mpair0(:,1,1)
!-------------------------------------------
call intgamma(mu,gamma0)
write(12,*)'gamma0',gamma0(1,:,1)
call SFT_ktox(gamma0,gamma0_xt)
!---------------------------------------
! do i1 = 1,n_
!   do i2 = 1,n_
!     if (ieee_is_nan(gamma0_xt(i1,i2,1))) then
!       write(10,*)'gamma0',gamma0(:,:,1)
!       pause
!     end if
!   enddo
! enddo
!------------------------------------------------------------
! write(13,*)'gamma0_xt',gamma0_xt(1,:,1)

!define G0(-x,-tau)=RE[G0(x,-tau)]-i*IM[G0(x,-tau)]      now  IM=0
call SFT_ktox(g0k_t,g0_x_t)         ! FT G0(k,-tau) from k to -x

call intself(mu,self0)
!---------------------------------------

!----------------------------------------
!define initial green function
g_k=g0k
midmatrix1=0.0D0;midmatrix2=0.0D0;midmatrix3=0.0D0;midmatrix4=0.0D0
!--------------main loop--------------------------------------
!---------------------------------------------------------------
!--------------------------------------------------------------
g_old=0D0
do loop=1,20
  write(*,*)'loop',loop
  write(12,*)'loop',loop

! define g1k(k,fome)
do i1=1,n_
  do i2=1,n_
    midmatrix1(i1,i2,1)=g_k(i1,i2,1)-g0k(i1,i2,1)
    midmatrix1(i1,i2,2)=g_k(i1,i2,2)-g0k(i1,i2,2)
  enddo
enddo
g1k=midmatrix1

call SFT_omegaFtotau(midmatrix1,midmatrix2,beta)    ! FT G1(k,fome) from omega to tau
call SFT_ktox(midmatrix2,midmatrix3)    ! FT G(k,tau) from k to x      !!!!!! midmatrix3 = G1(x,tau) in x-space
! write(13,*)'g1(x,tau)',midmatrix3(1,:,1)

midmatrix1 = 0.0D0
midmatrix2 = 0.0D0

midmatrix1=g0kt
call SFT_ktox(midmatrix1,midmatrix4)            !!!!!! midmatrix4 = G0(x,tau) in x-space
! write(13,*)'g0(x,tau)',midmatrix4(1,:,1)
!-----------------------------------------------------------
!define Mpair_1(x,tau)
do i1 = 1,n_
  do i2 = 1,n_
    Mpair_1(i1,i2,1)=2.0D0*midmatrix4(i1,i2,1)*midmatrix3(i1,i2,1)+midmatrix3(i1,i2,1)**2
    Mpair_1(i1,i2,2)=0D0!midmatrix3(i1,i2,2)**2    ! 0.0D0
  enddo
enddo

midmatrix1 = 0.0D0
midmatrix2 = 0.0D0
midmatrix3 = 0.0D0
midmatrix4 = 0.0D0

call SFT_xtok(Mpair_1,midmatrix1)      !  FT Mpair_1(x,tau) from x to q
call SFT_tautoomegaB(midmatrix1,midmatrix2,beta)          ! FT Mpair_1(q,tau) from tau to omegaB
!------------------------------------------------------------------

!define gamma^-1(q,omegaB) = eta/(8pi)+M
do i1 = 1,n_
  do i2 = 1,n_
    pi_q(i1,i2,1) = eta+8D0*pi*(Mpair(i1,i2,1)+Mpair0(i1,i2,1)+midmatrix2(i1,i2,1))
    pi_q(i1,i2,2) = 8D0*pi*(Mpair(i1,i2,2)+Mpair0(i1,i2,2)+midmatrix2(i1,i2,2))
  enddo
enddo

!-----------------------------------------------------------------
!define gamma(q,omegaB) = gamma_q = 1/(eta/(8pi)+Mpair)
do i1 = 1,n_
  do i2 = 1,n_
  mid = cmplx(8.0D0*pi,0.0D0)/cmplx(pi_q(i1,i2,1),pi_q(i1,i2,2))
  gamma_q(i1,i2,1) = real(mid)
  gamma_q(i1,i2,2) = aimag(mid)
  mid = cmplx(8.0D0*pi,0.0D0)/cmplx(eta+8D0*pi*Mpair0(i1,i2,1),8D0*pi*Mpair0(i1,i2,2))
  midmatrix3(i1,i2,1) = gamma_q(i1,i2,1) - real(mid)
  midmatrix3(i1,i2,2) = gamma_q(i1,i2,2) - aimag(mid)
  enddo
enddo

write(12,*)'gamma_q',gamma_q(1,:,1)
write(12,*)'midmatrix3',midmatrix3(:,n_,1)
!define gamma(q,tau) = gamma1+gamma0
call SFT_omegaBtotau(midmatrix3,midmatrix4,beta)
write(12,*)'gamma1',midmatrix4(1,:,1)
!-----------

midmatrix1 = 0.0D0
midmatrix2 = 0.0D0
! FT   define gamma(x,tau) = gamma_xt

call SFT_ktox(midmatrix4,midmatrix1)
gamma1_xt=midmatrix1
write(12,*)'gamma1_xt',gamma1_xt(1,:,1)

midmatrix1 = 0.0D0
midmatrix2 = 0.0D0
midmatrix3 = 0.0D0
midmatrix4 = 0.0D0

midmatrix1 = g1k
  tau = -tau
call SFT_omegaFtotau(midmatrix1,midmatrix2,beta)       ! FT G1(k,fome) from omega to -tau
  tau = -tau
midmatrix1 = 0.0D0
call SFT_ktox(midmatrix2,midmatrix1)         ! FT G1(k,-tau) from k to -x       !!!!!  midmatrix1 = G1(-x,-tau) in x-space
midmatrix2 = 0.0D0



!define Self_xt(x,tau)=G(-x,-tau)*gamma_xt(x,tau)
do i1=1,n_
  do i2=1,n_
    self_xt(i1,i2,1)=g0_x_t(i1,i2,1)*gamma1_xt(i1,i2,1)+midmatrix1(i1,i2,1)*(gamma1_xt(i1,i2,1)+gamma0_xt(i1,i2,1))
    self_xt(i1,i2,2)=0D0!midmatrix1(i1,i2,2)*(gamma1_xt(i1,i2,2)+gamma0_xt(i1,i2,2))    ! 0.0D0
  enddo
enddo

write(12,*)'g0_x_t',g0_x_t(1,:,1)
write(12,*)'self_xt',self_xt(1,:,1)

 call SFT_xtok(self_xt,midmatrix3)
 write(12,*)'midmatrix3',midmatrix3(1,:,1)
 call SFT_tautoomegaF(midmatrix3,midmatrix4,beta)
 self1=midmatrix4
 write(12,*)'self1',self1(:,1,1)
 write(12,*)'self0',self0(:,1,1)
 write(12,*)'self0+self1',self0(:,1,1)+self1(:,1,1)
!----------------------------
! call CPU_TIME(time2)
! write(*,*)'time',time2
!----------------------------
midmatrix1 = 0.0D0
midmatrix2 = 0.0D0
midmatrix3 = 0.0D0
midmatrix4 = 0.0D0
do i1 = 1,n_
  do i2 = 1,n_
    mid = cmplx(1.0D0,0.0D0)/cmplx(g0k(i1,i2,1),g0k(i1,i2,2))-cmplx(self0(i1,i2,1),self0(i1,i2,2))-cmplx(self1(i1,i2,1),self1(i1,i2,2))
    mid = cmplx(1.0D0,0.0D0)/mid
    g_k(i1,i2,1) = real(mid)
    g_k(i1,i2,2) = aimag(mid)
  enddo
enddo
write(12,*)'g_k',g_k(:,1,1)
!-----------------------------------------------
!----------------------------------------------------------------

    I_compare=0
   do i1=1,n_
     do i2=1,n_
       g_bijiao(i1,i2)=g_k(i1,i2,1)-g_old(i1,i2)
       if(abs(g_bijiao(i1,i2))>1d-4) then
         I_compare=I_compare+1
       end if
     end do
   end do
  if (I_compare==0) exit
    g_old=g_k(:,:,1)
  end do
!------------------------------
!------------end of main loop -----------------
!-----------------------------
g1k=0
do i1 = 1,n_
  do i2 = 1,n_
    mid1=cmplx(self0(i1,i2,1)+self1(i1,i2,1),self0(i1,i2,2)+self1(i1,i2,2))*cmplx(g_k(i1,i2,1),g_k(i1,i2,2))
    ! midmatrix1(i1,i2,1) = g_k(i1,i2,1) - g0k(i1,i2,1)
    ! midmatrix1(i1,i2,2) = g_k(i1,i2,2) - g0k(i1,i2,2)
    mid=mid1*cmplx(g0k(i1,i2,1),g0k(i1,i2,2))
    g1k(i1,i2,1)=real(mid)
    g1k(i1,i2,2)=aimag(mid)


    ep=k(i1)**2-mu
    mid=cmplx(ep,fome(i2))*mid1*cmplx(g0k(i1,i2,1),g0k(i1,i2,2))
    !mid=cmplx(ep,0.0D0)*cmplx(midmatrix1(i1,i2,1),midmatrix1(i1,i2,2))
    midmatrix3(i1,i2,1)=-real(mid)
    midmatrix3(i1,i2,2)=-aimag(mid)

    ! if (sqrt(real(mid1)**2+aimag(mid1)**2)>1D-2) then
      mid=-2D0*log(cmplx(1D0,0D0)+mid1)+2D0*mid1
    ! else
      ! mid=2D0*(0.5D0*mid1**2-mid1**3/3D0+mid1**4/4D0-mid1**5/5D0+mid1**6/6D0-mid1**7/7D0+mid1**8/8D0)
    ! end if
    midmatrix2(i1,i2,1)=real(mid)
    midmatrix2(i1,i2,2)=aimag(mid)

    ! mid=2D0*cmplx(g_k(i1,i2,1),g_k(i1,i2,2))*cmplx(self0(i1,i2,1)+self1(i1,i2,1),self0(i1,i2,2)+self1(i1,i2,2))
    ! midmatrix1(i1,i2,1)=real(mid)
    ! midmatrix1(i1,i2,2)=aimag(mid)

    ! mid=log(cmplx(pi_q(i1,i2,1),pi_q(i1,i2,2))/cmplx(eta+8D0*pi*Mpair0(i1,i2,1),8D0*pi*Mpair0(i1,i2,2)))
    ! midmatrix4(i1,i2,1)=-real(mid)
    ! midmatrix4(i1,i2,2)=-aimag(mid)

    ! iw_G0(i1)=2.0D0*ep/(exp(ep/t)+1.0D0)
    ! lnG0(i1)=2D0*t*log(exp(-ep/t)+1.0D0)
  enddo
enddo
! write(13,*)pi_q(:,1,1)
!------------------------------
tau=-tau
call SFT_omegaFtotau(g1k,g1k_t,beta)
write(12,*)'g1k_t',g1k_t(:,1,1)
dU1=0
call SFT_omegaFtotau(midmatrix3,dU1,beta)
! write(13,*)'u1=',dU1(:,1,1)
midmatrix3=0.0D0
call SFT_omegaFtotau(midmatrix2,midmatrix3,beta)
! write(13,*)'s1=',midmatrix3(:,1,1)
! midmatrix2=0.0D0
! call SFT_omegaFtotau(midmatrix1,midmatrix2,beta)
! ! write(13,*)'s2=',midmatrix2(:,1,1)
! midmatrix1=0.0D0
! call SFT_omegaBtotau(midmatrix4,midmatrix1,beta)
! write(13,*)'s3=',midmatrix1(:,1,1)
tau=-tau

! iny=0
! iny=g0k_t(:,1,1)
! call intk(n_,k,iny,num1)
! iny=0
! iny=g1k_t(:,1,1)
! call intk(n_,k,iny,num2)
! rho=num1+num2
call SFT_ktox(g1k_t,g1_t)
! write(13,*)'g1_t',g1_t(:,1,1)
rho=g1_t(1,1,1)+g0_x_t(1,1,1)
write(12,*)'rho',rho

!-------------------------------------

! call SFT_ktox(midmatrix3,midmatrix4)
! midmatrix3=0
! call SFT_ktox(midmatrix2,midmatrix3)
! midmatrix2=0
! call SFT_ktox(midmatrix1,midmatrix2)

! iny=0
! iny=midmatrix4(:,1,1)
! call intk(n_,k,iny,u2)
! call SFT_ktox(midmatrix1,midmatrix2)
iny=0

call fun_eng0(mu,u0)
! iny=dU1(:,1,1)
! call intk(n_,k,iny,u2)
midmatrix4=0.0D0
call SFT_ktox(dU1,midmatrix4)
u2=midmatrix4(1,1,1)
write(13,*)'u0',u0,'u2',u2
u1=u0+u2

call fun_entropy0(mu,s0)
midmatrix4=0
call SFT_ktox(midmatrix3,midmatrix4)
s3=midmatrix4(1,1,1)

! iny=0
! iny=midmatrix3(:,1,1)
! call intk(n_,k,iny,s3)
! midmatrix4=0
! call SFT_ktox(midmatrix1,midmatrix4)
! s5=midmatrix4(1,1,1)
call gammg_z(mu,s5)

write(13,*)'s0',s0
write(13,*)'s2',s2,'s3',s3,'s4',s4,'s5',s5
s1=s0+u2+s5+s3!+s4
s1=s1/t
end subroutine modeleq

subroutine fun_eng0(mu,eng)
use formidconstant
implicit none
real*8,intent(in):: mu
real*8,intent(out):: eng
real*8 c1,c2,c3,cp,cm,minl,maxl,step,a,b,sun,r0,f1,f2,f3,k1,k2,k3
integer n,i
real*8 ep,xi,nf
parameter  (c1=0.5555555556D0)
parameter  (c2=0.8888888889D0)
parameter  (c3=0.7745966692D0)
!vext = 0.25D0*omega**2*r**2
minl = 0.0D0
maxl = 1.0D3
n = 100000
step = (maxl-minl)/dble(n)
a = minl
b = minl + step

sun = 0.0D0
r0 = 0.0D0

do i = 1,n
  cp = (a+b)/2.0D0
  cm = (a-b)/2.0D0
  k1 = cp-cm*c3
  k2 = cp
  k3 = cp + cm*c3

  ep = k1**2-mu
  !xi = sqrt(ep**2+phi**2)
  nf = 1D0/(1.0D0+exp(ep/t))
  f1 = k1**2*ep*nf

  ep = k2**2-mu
  nf = 1.0D0/(1.0D0+exp(ep/t))
  f2 = k2**2*ep*nf

  ep = k3**2-mu
  nf = 1.0D0/(1.0D0+exp(ep/t))
  f3 = k3**2*ep*nf

  sun = -cm*(c1*f1+c2*f2+c1*f3)
  r0 = r0 + sun
  a = b
  b = b + step
enddo
   eng = r0*2.0D0*0.05066059182D0
end subroutine fun_eng0

subroutine fun_entropy0(mu,s)
  use formidconstant
  implicit none
  real*8,intent(in):: mu
  real*8,intent(out):: s
  real*8 c1,c2,c3,cp,cm,minl,maxl,step,a,b,sun,r0,f1,f2,f3,k1,k2,k3
  integer n,i
  real*8 ep,xi,nf
  parameter  (c1=0.5555555556D0)
  parameter  (c2=0.8888888889D0)
  parameter  (c3=0.7745966692D0)
  !vext = 0.25D0*omega**2*r**2
  minl = 0.0D0
  maxl = 1.0D3
  n = 100000
  step = (maxl-minl)/dble(n)
  a = minl
  b = minl + step

  sun = 0.0D0
  r0 = 0.0D0

  do i = 1,n
    cp = (a+b)/2.0D0
    cm = (a-b)/2.0D0
    k1 = cp-cm*c3
    k2 = cp
    k3 = cp + cm*c3

    ep = k1**2-mu
    nf = 1.0D0/(1.0D0+exp(ep/t))
    f1 = k1**2*(ep*nf+log(1.0D0+exp(-ep/t))*t)

    ep = k2**2-mu
    nf = 1.0D0/(1.0D0+exp(ep/t))
    f2 = k2**2*(ep*nf+log(1.0D0+exp(-ep/t))*t)

    ep = k3**2-mu
    nf = 1.0D0/(1.0D0+exp(ep/t))
    f3 = k3**2*(ep*nf+log(1.0D0+exp(-ep/t))*t)

    sun = -cm*(c1*f1+c2*f2+c1*f3)
    r0 = r0 + sun
    a = b
    b = b + step
  enddo
     s = 2.0D0*r0*0.05066059182D0
  end subroutine fun_entropy0

subroutine gammg_z(mu,gamma_g)
  !for t*Sum[ln(gamma0/g)]
  use omp_lib
  Use gsld
  use vec
  use formidconstant
  implicit none
  real*8,intent(in):: mu
  real*8,intent(out):: gamma_g
  real*8::minl,maxl,minl1,maxl1,step,a,b,a1,b1,answer,r0
  integer n0,i,i1,i3,j

  minl = 0.0D0
  maxl = 1.0D2
  minl1=0.0D0
  maxl1=1.0D2
  n0 = 1000
  step = (maxl-minl)/dble(n0)


if (abs(eta)<=0.05D0) then
  a = minl
  b = minl + step
  r0 = 0.0D0
  do i3 = 1,n0
    answer = 0
      Do i = 1,nk
        answer = answer + ak1(i)*y1((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),mu)
      end do
    answer = answer*(b-a)/2.0D0
    r0 = r0+answer
    a = b
    b = b+step
  end do
  gamma_g=-0.5D0*t*r0*0.05066059182D0
else
  a1= minl1
  b1= minl1+step
  r0 = 0.0D0
do i1 = 1,n0
  a = minl
  b = minl + step
  do i3 = 1,n0
    answer = 0
    Do j=1,nk
      Do i = 1,nk
        answer = answer + ak1(j)*ak1(i)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu)
      end do
    End Do
    answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
    r0 = r0+answer
    a = b
    b = b+step
  end do
  a1=b1
  b1=b1+step
end do
  gamma_g=r0
  gamma_g=eta*t/pi*gamma_g*0.05066059182D0
end if
contains
  Function y(z,q,mu)
     Implicit none
     Real*8:: y
     Real*8:: z,mu,beta,ep,q

     beta=1.0D0/t
     ep=q**2-4.0D0*mu
     y=q**2*log(dabs(1D0-exp(-beta/2D0*(z+ep))))/(4D0*eta**2+z)/sqrt(z)
   End Function y
   Function y1(q,mu)
      Implicit none
      Real*8:: y1
      Real*8:: mu,beta,ep,q

      beta=1.0D0/t
      ep=q**2-4.0D0*mu
      y1=q**2*log(dabs(1D0-exp(-beta/2D0*ep)))
    End Function y1
end subroutine gammg_z

subroutine intk(n,x,y,r)
  implicit none
integer,intent(in)::n
real*8,intent(in):: x(n),y(n)
real*8,intent(out):: r
integer n0
real*8 mid
r = 0.0D0
do n0 = 1,n-1
mid = (y(n0)*x(n0)**2+y(n0+1)*x(n0+1)**2)*(x(n0+1)-x(n0))/2.0D0
r = r+mid
enddo
 r = r*0.05066059182D0
end subroutine intk

Subroutine intMpair(mu,beta,Mpair)
  Use gsld
  use vec
  use formidconstant
  use omp_lib
  Implicit None
  Real*8:: a, b, a1, b1, answer(2)
  real*8:: minl,maxl,r0(2),step,minl1,maxl1
  real*8::mu,beta,Mpair(n_,n_,2),ep
  Integer :: i, j,n0,i1,i2,i3,ncho
  integer ndim,ncomp,flags,mineval,maxeval,key,nregions,neval,fail
  real*8 userdata(3),epsrel,epsabs,integral,error,prob
  integer seed,nnew,nvec
  ! real*8 flatness
  character*(*) statefile
  integer*8 spin
  parameter (spin = -1)
  parameter (statefile = "")
  parameter (nvec = 1)
  parameter (ndim = 2)
  parameter (ncomp = 1)
  parameter (epsrel = 1.0D-5)
  parameter (epsabs = 1.0D-5)
  parameter (flags = 4)
  parameter (mineval = 1000)
  parameter (maxeval = 80000)
  parameter (key = 13)
  parameter (seed = 0)
  ! parameter (nnew = 100)
  ! parameter (flatness = 1.0D2)
  real*8::time0,time1
!-----------------------------------------------
! do i1=1,n_
!   ncho=i1-1
!   if (4.0D0*mu <= k(i1)**2) exit
! end do
!-----------------------------------------------
  !mu = mu-vext

  minl = 0.0D0
  maxl = 1.0D2
  a1 = -1.0D0
  b1 = 1.0D0

  n0 = 1000
  step = (maxl-minl)/dble(n0)

 ! call OMP_set_dynamic(.true.)
 !call OMP_set_num_threads(8)
 time0=omp_get_wtime()
 !$OMP PARALLEL
 !$OMP DO PRIVATE(i1,i2,i3,i,j,r0,a,b,answer)  SCHEDULE(DYNAMIC,1)
  do i1=1,n_
    do i2=1,n_-1
      a = minl
      b = minl + step
      r0 = 0.0D0

      do i3 = 1,n0
        answer = 0
        Do j = 1, nk
          Do i = 1, nk
    	      answer = answer + ak1(i)*ak1(j)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,beta,k(i1),bome(i2))
          End Do
        End Do
        answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
        r0 = r0+answer
        a = b
        b = b+step
      end do
    Mpair(i1,i2,1)=r0(1)
    Mpair(i1,i2,2)=r0(2)
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  ! time1=omp_get_wtime()
  ! write(*,*) "Work took", time1-time0, "seconds"
! ----------------- bome(n_)=0-----------
  ! do i1=1,ncho
      ! userdata(1) = mu
      ! userdata(2) = beta
      ! userdata(3) = k(i1)
      ! call cuhre(ndim,ncomp,fpi,userdata,nvec,epsrel,epsabs,flags,mineval,maxeval,key,statefile,spin,nregions,neval,fail,integral,error,prob)
      ! Mpair(i1,n_,1)=integral
      ! Mpair(i1,n_,2)=0.0D0

      ! r0 = 0.0D0
      ! minl1=4.0D0
      ! maxl1=1.0D3
      ! n0=1000
      ! step=(maxl1-minl1)/dble(n0)
      ! a = minl1
      ! b = minl1 + step
      ! do i3 = 1,n0
      !   answer = 0
      !   Do j = 1, nk
      !     Do i = 1, nk
      !       answer = answer + ak1(i)*ak1(j)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,beta,k(i1),0.0D0)
      !     End Do
      !   End Do
      !   answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
      !   r0 = r0+answer
      !   a = b
      !   b = b+step
      ! end do
      ! Mpair(i1,n_,1)=Mpair(i1,n_,1)+r0(1)
      ! Mpair(i1,n_,2)=0.0D0!Mpair(i1,1,2)+r0(2)
      ! Mpair(i1,n_,1)=0D0
      ! Mpair(i1,n_,2)=0D0
  ! end do
  !!$OMP DO PRIVATE(i1,i3,i,j,ep,minl,maxl,n0,step,r0,a,b,answer) REDUCTION(+:Mpair) SCHEDULE(DYNAMIC,1)
do i1=1,n_
  ep=k(i1)**2-4.0D0*mu
  if (ep<0) then
  minl = 0.0D0
  maxl = sqrt(mu-k(i1)**2/4D0)-1D-3
  n0 = 10
  step = (maxl-minl)/dble(n0)
  a = minl
  b = minl + step
  r0 = 0.0D0
  do i3 = 1,n0
    answer = 0
    Do j=1,nk
      Do i=1,nk
        answer(1)=answer(1)+ak1(i)*ak1(j)*y1((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,beta,k(i1))
      End Do
    End Do
    answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
    r0 = r0+answer
    a = b
    b = b+step
  end do
  Mpair(i1,n_,1)=r0(1)
  Mpair(i1,n_,2)=0.0D0

  minl = sqrt(mu-k(i1)**2/4D0)+1D-3
  maxl = 2.0D0*sqrt(mu-k(i1)**2/4D0)
  n0 = 10
  step = (maxl-minl)/dble(n0)
  a = minl
  b = minl + step
  r0 = 0.0D0
  do i3 = 1,n0
    answer = 0
    Do j=1,nk
      Do i=1,nk
        answer(1) = answer(1) + ak1(i)*ak1(j)*y1((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,beta,k(i1))
      End Do
    End Do
    answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
    r0 = r0+answer
    a = b
    b = b+step
  end do
  Mpair(i1,n_,1)=Mpair(i1,n_,1)+r0(1)
  Mpair(i1,n_,2)=0.0D0

  minl = 2.0D0*sqrt(mu-k(i1)**2/4D0)
  maxl = 1.0D2
  n0 = 1000
  step = (maxl-minl)/dble(n0)
  a = minl
  b = minl + step
  r0 = 0.0D0
  do i3 = 1,n0
    answer = 0
    Do j=1,nk
      Do i=1,nk
        answer(1) = answer(1) + ak1(i)*ak1(j)*y1((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,beta,k(i1))
      End Do
    End Do
    answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
    r0 = r0+answer
    a = b
    b = b+step
  end do
  Mpair(i1,n_,1)=Mpair(i1,n_,1)+r0(1)
  Mpair(i1,n_,2)=0.0D0
else
  minl = 0.0D0
  maxl = 1.0D2
  n0=1000
  step=(maxl-minl)/dble(n0)

  ! do i1=ncho+1,n_
      a = minl
      b = minl + step
      r0 = 0.0D0
      do i3 = 1,n0
        answer = 0
        Do j = 1, nk
          Do i = 1, nk
    	      answer = answer + ak1(i)*ak1(j)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,beta,k(i1),0.0D0)
          End Do
        End Do
        answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
        r0 = r0+answer
        a = b
        b = b+step
      end do
        Mpair(i1,n_,1)=r0(1)
        Mpair(i1,n_,2)=0D0!r0(2)
  end if
end do
   !!$OMP END DO
   !!$OMP END PARALLEL
   time1=omp_get_wtime()
   write(*,*) "Work took", time1-time0, "seconds"
!----------------------------------------
  Mpair=Mpair*0.025330296D0
contains
 Function y(k0,the,mu,beta,q,omega) !被积函数
    Implicit none
    Real*8:: y(2)
    Real*8:: k0,the,mu,ep1,ep2,ferm1,ferm2,beta,q,omega
    complex*16::mid
    ep1=k0**2-mu+q**2/4D0-k0*q*the
    ep2=k0**2-mu+q**2/4D0+k0*q*the
    ferm1=1.0D0/(exp(ep1*beta)+1.0D0)
    ferm2=1.0D0/(exp(ep2*beta)+1.0D0)
    mid=cmplx((-ferm1-ferm2)*k0**2,0.0D0)/cmplx(ep1+ep2,-omega)

    y(1)=real(mid)
    y(2)=aimag(mid)
  End Function y

  Function y1(k0,the,mu,beta,q)
     Implicit none
     Real*8:: y1
     Real*8:: k0,the,mu,ep1,ep2,ferm1,ferm2,beta,q
     ep1=k0**2-mu+q**2/4D0-k0*q*the
     ep2=k0**2-mu+q**2/4D0+k0*q*the
     ferm1=1.0D0/(exp(ep1*beta)+1.0D0)
     ferm2=1.0D0/(exp(ep2*beta)+1.0D0)
     y1=(-ferm1-ferm2)*k0**2/(2D0*k0**2-2D0*mu+q**2/2D0)
   End Function y1

  function fpi(ndim,x,ncomp,rr,userdata)
  integer fpi
  real*8,intent(in):: x(2)
  integer,intent(in):: ndim,ncomp
  real*8,intent(out):: rr
  real*8 userdata(3),mu,beta,q
  real*8 ep,ep1,ferm,ferm1,k0,the
  mu=userdata(1)
  beta=userdata(2)
  q=userdata(3)

  k0=4.0D0*x(1)
  the=x(2)*2.0D0-1.0D0

  ep=k0**2-mu
  ep1=k0**2-mu+q**2-2.0D0*k0*q*the
  ferm=1.0D0/(exp(ep*beta)+1.0D0)
  ferm1=1.0D0/(exp(ep1*beta)+1.0D0)
  if (abs(ep+ep1) .gt. 1.0D-3) then
    rr=(-ferm-ferm1)*k0**2/(ep+ep1)*2.0D0*4.0D0
  else
    rr=0.0D0
  endif
  fpi = 0
  end function fpi
End subroutine intMpair

subroutine intgamma(mu,gamma0)
  Use gsld
  use formidconstant
  use vec
  use omp_lib
  Implicit None
  Real*8 :: a, b, answer
  real*8 minl,maxl,r0,step,ep
  real*8::mu,vext,gamma0(n_,n_,2)
  Integer :: i,j,n0,i1,i2,i3
  real*8::time0,time1
  !-------------------------
  !mu=mu-vext
  gamma0=0.0D0

  time0=omp_get_wtime()
    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i,j,ep,minl,maxl,r0,a,b,answer) SCHEDULE(DYNAMIC,1)
  do i1=1,n_
    ep=k(i1)**2-4.0D0*mu
    if (ep<0) then
      do i2=1,n_
        minl = 0.0D0
        maxl = -ep-1D-6
        n0 = 10
        step = (maxl-minl)/dble(n0)
        a = minl
        b = minl + step
        r0 = 0.0D0
        do i3 = 1,n0
          answer = 0
          Do i = 1, nk
            answer = answer + ak1(i)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),mu,k(i1),tau(i2))
          End Do
          answer = answer*(b-a)/2.0D0
          r0 = r0+answer
          a = b
          b = b+step
        end do
        gamma0(i1,i2,1)=r0

        minl = -ep+1D-6
        maxl = -2.0D0*ep
        n0 = 10
        step = (maxl-minl)/dble(n0)
        a = minl
        b = minl + step
        r0 = 0.0D0
        do i3 = 1,n0
          answer = 0
          Do i = 1, nk
            answer = answer + ak1(i)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),mu,k(i1),tau(i2))
          End Do
          answer = answer*(b-a)/2.0D0
          r0 = r0+answer
          a = b
          b = b+step
        end do
        gamma0(i1,i2,1)=gamma0(i1,i2,1)+r0

        minl = -2.0D0*ep
        maxl = 1.0D4
        n0 = 8000
        step = (maxl-minl)/dble(n0)
        a = minl
        b = minl + step
        r0 = 0.0D0
        do i3 = 1,n0
          answer = 0
          Do i = 1, nk
            answer = answer + ak1(i)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),mu,k(i1),tau(i2))
          End Do
          answer = answer*(b-a)/2.0D0
          r0 = r0+answer
          a = b
          b = b+step
        end do
        gamma0(i1,i2,1)=gamma0(i1,i2,1)+r0
      end do
    else
      do i2=1,n_
        minl = 0.0D0
        maxl = 1.0D4
        n0 = 8000
        step = (maxl-minl)/dble(n0)
        a = minl
        b = minl + step
        r0 = 0.0D0
        do i3 = 1,n0
          answer = 0
          Do i = 1, nk
            answer = answer + ak1(i)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),mu,k(i1),tau(i2))
          End Do
          answer = answer*(b-a)/2.0D0
          r0 = r0+answer
          a = b
          b = b+step
        end do
        gamma0(i1,i2,1)=r0
      end do
    end if
  end do
  !$OMP END PARALLEL DO
  time1=omp_get_wtime()
  write(*,*) "Work took", time1-time0, "seconds"

contains
  Function y(z,mu,q,tau0) !被积函数
     Implicit none
     Real*8:: y
     Real*8:: z,mu,beta,ep,bose,q,tau0
     beta=1.0D0/t
     ep=q**2-4.0D0*mu
     bose=exp(-0.5D0*(ep+z)*tau0)/(exp(-0.5D0*(ep+z)*beta)-1.0D0)
    !  bose=1.0D0/(exp(0.5D0*(ep+z)*(beta-tau0))-exp(-0.5D0*(ep+z)*tau0))
     y=8.0D0*sqrt(abs(z))/(4.0D0*eta**2+z)*bose
   End Function y
 end subroutine intgamma

 subroutine intself(mu,self0)
   Use gsld
   use formidconstant
   use vec
   use omp_lib
   Implicit None
   Real*8:: a, b, answer(2),a1,b1,sun(2)
   real*8 minl,maxl,r0(2),r1(2),step,minl1,maxl1,ep
   real*8::mu,vext,self0(n_,n_,2),self1(n_,n_,2),self2(n_,n_,2),self_temp(n_,2*n_)
   Integer :: i,j,n0,i1,i2,i3,l
   !-------------------------
   integer ndim,ncomp,flags,mineval,maxeval,key,nregions,neval,fail,core
  !  parameter (ncomp=2)
  !  real*8 userdata(3),epsrel,epsabs,integral(2),error(2),prob(2)
   parameter (ncomp = 2*n_)
   real*8 userdata(2),epsrel,epsabs,integral(ncomp),error(ncomp),prob(ncomp)
   integer seed,nnew,nvec
   real*8 flatness
   character*(*) statefile
   integer*8 spin
   parameter (core=-2)
   parameter (spin=-1)
   parameter (statefile = "")
   parameter (nvec = 1)
   parameter (ndim = 3)
   parameter (epsrel = 1.0D-3)
   parameter (epsabs = 1.0D-3)
   parameter (flags = 4)
   parameter (mineval = 5000)
   parameter (maxeval = 50000)    !!!!!!!!!!
   parameter (key = 11)
   parameter (seed = 0)
  !  parameter (nnew = 1000)
  !  parameter (flatness = 25D0)
   real*8::time0,time1,time2
   !!!!!----------------------
   !mu=mu-vext
   self0=0;self1=0;self2=0

   minl = 0.0D0
   maxl = 1.0D3
   a1 = -1.0D0
   b1 = 1.0D0

   n0 = 1000
   step = (maxl-minl)/dble(n0)

   time0=omp_get_wtime()
  !$OMP PARALLEL
  !$OMP DO PRIVATE(i1,i2,i3,i,j,r0,a,b,answer) SCHEDULE(DYNAMIC,1)
   do i1=1,n_
     do i2=1,n_
       a = minl
       b = minl + step
       r0=0.0D0
       !r1=0.0D0
       do i3 = 1,n0
         answer=0
         !sun=0
           Do j = 1, nk
             Do i = 1, nk
       	      answer=answer+ak1(i)*ak1(j)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,k(i1),fome(i2))
              !sun=sun+ak1(i)*ak1(j)*y2((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,k(i1),fome(i2))
             End Do
           End Do
           answer=answer*(b-a)/2.0D0*(b1-a1)/2.0D0
           !sun=sun*(b-a)/2.0D0*(b1-a1)/2.0D0
         r0=r0+answer
         !r1=r1+sun
         a = b
         b = b+step
       end do
       self1(i1,i2,1)=r0(1)!+r1(1)
       self1(i1,i2,2)=r0(2)!+r1(2)
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  time1=omp_get_wtime()
  write(*,*) "Work took", time1-time0, "seconds"

!--------------------fast but not sure if converged------------------------------------------
  do i1 = 1,n_
    userdata(1) = mu
    userdata(2) = k(i1)
    ! call cubafork(spin)
    ! spin=0
    call cuhre(ndim,ncomp,y1,userdata,nvec,epsrel,epsabs,flags,mineval,maxeval,key,statefile,spin,nregions,neval,fail,integral,error,prob)
    self_temp(i1,:) = integral
    ! call cubawait(spin)
  enddo

  ! l=0
  do i1=1,n_
    do i2=1,n_
      ! l=l+1
      self2(i1,i2,1)=self_temp(i1,i2)
      self2(i1,i2,2)=self_temp(i1,n_+i2)
    end do
  end do
!------------------------------------------------------------------
! do i1=1,n_
!   do i2=1,n_
!     userdata(1) = mu
!     userdata(2) = k(i1)
!     userdata(3) = fome(i2)
!     call cuhre(ndim,ncomp,y3,userdata,nvec,epsrel,epsabs,flags,mineval,maxeval,key,statefile,spin,nregions,neval,fail,integral,error,prob)
!     self2(i1,i2,1)=integral(1)
!     self2(i1,i2,2)=integral(2)
!   end do
! enddo



  time2=omp_get_wtime()
  write(*,*) "Work took", time2-time1, "seconds"
  ! write(10,*)'self1',self1(1,:,1)
  ! write(10,*)'self2',self2
  self0=self1+self2
  self0=self0*0.025330296D0
contains
 Function y(q,the,mu,k0,omega) !被积函数
    Implicit none
    Real*8:: y(2)
    Real*8:: k0,the,mu,ep,ep1,ferm,q,omega
    complex*16::mid
    ! ep=q**2+k0**2+2.0D0*k0*q*the-4.0D0*mu
    ! ep1=q**2-mu
     ep=q**2-4.0D0*mu
     ep1=k0**2-mu+q**2-2.0D0*k0*q*the
    ferm=1.0D0/(exp(ep1/t)+1.0D0)
    mid=cmplx(-16.0D0*pi*ferm*q**2,0.0D0)/(cmplx(2.0D0*eta,0.0D0)-sqrt(cmplx(-2.0D0*ep1+ep,-2.0D0*omega)))!+cmplx(0D0,1D0)*sqrt(cmplx(2.0D0*ep1-ep,2.0D0*omega)))

    y(1)=real(mid)
    y(2)=aimag(mid)
  End Function y

  function y1(ndim,x,ncomp,rr,userdata,nvec,core)
    implicit none
    integer y1
    real*8,intent(in):: x(ndim,nvec)
    integer,intent(in):: ndim,ncomp,nvec,core
    real*8,intent(out):: rr(ncomp,nvec)
    real*8:: userdata(2),mu,q,k0,the
    Real*8:: z,ep,ep1,bose,omega
    complex*16::mid
    integer::i1,i2,l

    mu=userdata(1)
    k0=userdata(2)

    q=1.0D3*x(1,nvec)
    the=x(2,nvec)*2.0D0-1.0D0
    z=1.0D4*x(3,nvec)

    ep=q**2-4.0D0*mu
    ! l=0
    ! do i1=1,n_
    ep1=k0**2-mu+q**2-2.0D0*k0*q*the
    ! ep=q**2+k0**2+2.0D0*k0*q*the-4.0D0*mu
    ! ep1=q**2-mu
    bose=1.0D0/(exp(0.5D0*(ep+z)/t)-1.0D0)
      do i2=1,n_
        ! l=l+1
        if (abs(ep) .gt. 1.0D-5) then
          mid=cmplx(16.0D0*sqrt(z)/(4.0D0*eta**2+z)*bose*q**2,0.0D0)/cmplx(-z+2.0D0*ep1-ep,2.0D0*fome(i2))
        else
          mid=cmplx(0.0D0,0.0D0)
        end if
        rr(i2,nvec)=real(mid)
        rr(n_+i2,nvec)=aimag(mid)
      end do
    ! end do
      rr = rr*1.0D3*2.0D0*1.0D4

    y1=0
  End Function y1

! function y2(q,the,mu,k0,omega)
!   implicit none
!   real*8::y2(2)
!   real*8::ep,ep1,z,q,mu,k0,the,omega,bose
!   Integer:: i,n0,i3
!   complex*16::mid
!
!   ep=q**2-4.0D0*mu
!   ep1=k0**2-mu+q**2-2.0D0*k0*q*the
!   ! bose=1.0D0/(exp(0.5D0*(ep+z)/t)-1.0D0)
!   ! mid=cmplx(16.0D0*sqrt(z)/(4.0D0*eta**2+z)*bose,0.0D0)/cmplx(-z+2.0D0*ep1-ep,2.0D0*omega)
!   ! y2(1)=real(mid)
!   ! y2(2)=aimag(mid)
!   ! y2=y2*q**2
!     if (ep<0.0D0) then
!         a = 0.0D0
!         b = -ep-1D-3
!         r0 = 0.0D0
!           answer = 0
!           Do i = 1, nk
!             answer = answer + ak1(i)*y3((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),ep,ep1,omega)
!           End Do
!           answer = answer*(b-a)/2.0D0
!           r0 = r0+answer
!           a = b
!           b = b+step
!         y2=r0
!
!         a = -ep+1D-3
!         b = -2.0D0*ep
!         r0 = 0.0D0
!         do i3 = 1,n0
!           answer = 0
!           Do i = 1, nk
!             answer = answer + ak1(i)*y3((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),ep,ep1,omega)
!           End Do
!           answer = answer*(b-a)/2.0D0
!           r0 = r0+answer
!           a = b
!           b = b+step
!         end do
!         y2=y2+r0
!
!         minl = -2.0D0*ep
!         maxl = 1.0D2
!         n0 = 100
!         step = (maxl-minl)/dble(n0)
!         a = minl
!         b = minl + step
!         r0 = 0.0D0
!         do i3 = 1,n0
!           answer = 0
!           Do i = 1, nk
!             answer = answer + ak1(i)*y3((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),ep,ep1,omega)
!           End Do
!           answer = answer*(b-a)/2.0D0
!           r0 = r0+answer
!           a = b
!           b = b+step
!         end do
!         y2=y2+r0
!     else
!         minl = 0.0D0
!         maxl = 1.0D2
!         n0 = 100
!         step = (maxl-minl)/dble(n0)
!         a = minl
!         b = minl + step
!         r0 = 0.0D0
!         do i3 = 1,n0
!           answer = 0
!           Do i = 1, nk
!             answer = answer + ak1(i)*y3((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),ep,ep1,omega)
!           End Do
!           answer = answer*(b-a)/2.0D0
!           r0 = r0+answer
!           a = b
!           b = b+step
!         end do
!         y2=r0
!     end if
!     y2=y2*q**2
! end function y2

 ! function y3(z,ep,ep1,omega)
 !   implicit none
 !   Real*8:: y3(2)
 !   Real*8:: z,ep,ep1,bose,omega
 !   complex*16::mid
 !   bose=1.0D0/(exp(0.5D0*(ep+z)/t)-1.0D0)
 !   mid=cmplx(16.0D0*sqrt(z)/(4.0D0*eta**2+z)*bose,0.0D0)/cmplx(-z+2.0D0*ep1-ep,2.0D0*omega)
 !   y3(1)=real(mid)
 !   y3(2)=aimag(mid)
 ! End Function y3
 function y3(ndim,x,ncomp,rr,userdata,nvec,core)
   implicit none
   integer y3
   real*8,intent(in):: x(ndim,nvec)
   integer,intent(in):: ndim,ncomp,nvec,core
   real*8,intent(out):: rr(2,nvec)
   real*8:: userdata(3),mu,q,k0,the
   Real*8:: z,ep,ep1,bose,omega
   complex*16::mid
   integer::i1,i2,l

   mu=userdata(1)
   k0=userdata(2)
   omega=userdata(3)

   q=1.0D3*x(1,nvec)
   the=x(2,nvec)*2.0D0-1.0D0
   z=1.0D4*x(3,nvec)

  !  ep=q**2-4.0D0*mu
  !  ep1=k0**2-mu+q**2-2.0D0*k0*q*the
   ep=q**2+k0**2+2.0D0*k0*q*the-4.0D0*mu
   ep1=q**2-mu
   bose=1.0D0/(exp(0.5D0*(ep+z)/t)-1.0D0)
       if (abs(ep) .gt. 1.0D-5) then
         mid=cmplx(16.0D0*sqrt(z)/(4.0D0*eta**2+z)*bose*q**2,0.0D0)/cmplx(-z+2.0D0*ep1-ep,2.0D0*omega)
       else
         mid=cmplx(0.0D0,0.0D0)
       end if
       rr(1,nvec)=real(mid)
       rr(2,nvec)=aimag(mid)
     rr = rr*1.0D3*2.0D0*1.0D4

   y3=0
 End Function y3
end subroutine intself
