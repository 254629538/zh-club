
  subroutine intd(mu,f)
    use formidconstant
    use gauss_int
    implicit none
    real*8,intent(in):: mu
    real*8,intent(out):: f
    integer n,i
    real*8 minl,maxl,r0,userdata(10)
    call gausspoint80(gp0,gr0)
    userdata(1)=mu
    minl = 0.0D0
    maxl = 5D0*sqrt(2D0*t/omega0)
    call gauss_int_1d80(minl,maxl,func1,gp0,gr0,r0,userdata)
    f=0.5D0*r0
  contains
  FUNCTION func(minl,maxl,x,f,userdata)
    implicit none
    integer func
    real*8,intent(in)::minl,maxl,x,userdata(10)
    real*8,intent(out)::f
    real*8::a,b,y,k,mu,r,vext

    r=userdata(1)
    mu=userdata(2)

    a=minl
    b=maxl
    k=((b-a)*x+(b+a))/2.0D0

    vext=0.5D0*omega0*r**2!mu*r**2
    f=k**1.5D0/(exp(k-(mu-vext)/t)+1.0D0)

    f=(b-a)*f/2.0D0
    func=0
  end FUNCTION func
  FUNCTION func1(minl,maxl,x,f,userdata)
    implicit none
    integer func1
    real*8,intent(in)::minl,maxl,x,userdata(10)
    real*8,intent(out)::f
    real*8::a,b,y,k,r,userdata1(10),mu
    real*8::minl1,maxl1,step,a1,b1,sun,r0
    integer n,i

    mu=userdata(1)

    a=minl
    b=maxl
    r=((b-a)*x+(b+a))/2.0D0

    userdata1(1)=r
    userdata1(2)=mu

    minl1 = 0.0D0
    maxl1 = 1.0D3
    n = 1000
    step = (maxl1-minl1)/dble(n)
    a1 = minl1
    b1 = minl1 + step
    sun = 0.0D0
    r0 = 0.0D0
    do i = 1,n
    call gauss_int_1d80(a1,b1,func,gp0,gr0,sun,userdata1)
    r0 = r0+sun
    a1 = b1
    b1 = b1+step
    enddo

    f=r**2*r0
    f=(b-a)*f/2.0D0
    func1=0
  end FUNCTION func1
  end subroutine intd

!omega define the potential field
!mu is the chemical potential
!r give the damping rate 1/\tau_T
subroutine rate(mu,tau22)
  use gsld
  use formidconstant
  implicit none
real*8,intent(in):: mu
real*8,intent(out):: tau22
integer ndim,ncomp,flags,mineval,maxeval,key,nregions,neval,fail
real*8 userdata(4),epsrel,epsabs,integral(1),error(1),prob(1),r,ero,sanse,c
integer seed,nnew,nvec
real*8 flatness,pi
character*(*) statefile
integer*8 spin
parameter (pi=3.141592653589793238462643383279502884197D0)
parameter (spin = -1)
parameter (statefile = "")
parameter (nvec = 1)
parameter (ndim = 5)
parameter (ncomp = 1)
parameter (epsrel = 1.0D-5)
parameter (epsabs = 1.0D-5)
parameter (flags = 6)
parameter (mineval = 5000)
parameter (maxeval = 500000)
parameter (key = 13)
parameter (seed = 0)
parameter (nnew = 1400)              !nnew = 1000          for vacuum tmatrix
parameter (flatness = 80.0D0)      !flatness = 50.0D0

userdata(1) = mu
 Call fn0(fn1, ak1)
 call intd(mu,r)
 call suave(ndim,ncomp,integrated,userdata,nvec,epsrel,epsabs,flags,seed,mineval,maxeval,nnew,flatness,statefile,spin,nregions,neval,fail,integral,error,prob)
 sanse = integral(1)
 ero = error(1)

 tau22 = sqrt(2.0D0)/(5.0D0*pi)*(t/eta)**2*sanse/r*(3.0D0*lambda*ntot)**(1.0D0/3.0D0)/2.0D0!/(4.0D0*mu)
 tau22=1.0D0/tau22
 ! write(11,*) 'D=',r
 ! write(11,*) 'Ic=',sanse
 ! write(11,*)'tau22=',tau22
 ! write(11,*)'nregions,neval,error,prob',nregions,neval,ero,prob
 !*************************************
 !
 !*************************************
 contains
 function integrated(ndim,x,ncomp,f,userdata)
   use formidconstant
   implicit none
 integer integrated
 integer,intent(in):: ndim,ncomp
 real*8,intent(in):: x(ndim),userdata(4)
 real*8,intent(out):: f(1)
 real*8 mu,z,vext,a,b
 real*8 x0,x1,r,y,y1,fu,p0,p1
 real*8 x12,x22,x32,x42
 real*8 re,im,tmatrix
 real*8 dd

 mu = userdata(1)

 a = 1.0D2
 b = 2.0D0                          !beware of this change!!!!!!!
 c = 5D0*sqrt(2D0*t/omega0)
 x0 = x(1)*a
 x1 = x(2)*a

 y = x(3)*b-1.0D0
 y1 = x(4)*b-1.0D0
 r = x(5)*c

 !x0 = sqrt(0.5D0*p0**2/t)
 !x1 = sqrt(2.0D0*p1**2/t)

 vext = 0.5D0*omega0*r**2!mu*r**2!
 z = exp((mu-vext)/t)
 x12 = 0.5D0*(x0**2+x1**2+2.0D0*x0*x1*y)
 x22 = 0.5D0*(x0**2+x1**2-2.0D0*x0*x1*y)
 x32 = 0.5D0*(x0**2+x1**2+2.0D0*x0*x1*y1)
 x42 = 0.5D0*(x0**2+x1**2-2.0D0*x0*x1*y1)
 fu = z**2*exp(-x0**2-x1**2)/((1.0D0+z*exp(-x12))*(1.0D0+z*exp(-x22))*(1.0D0+z*exp(-x32))*(1.0D0+z*exp(-x42)))
 !0.07957747155 = 1/(4Pi)

  p1=x1*sqrt(t/2.0D0)
  p0=x0*sqrt(t*2.0D0)
  call pair0(mu,vext,p0,p1,re,im)
  tmatrix = 1.0D0/((1.0D0+8.0D0*pi/eta*re)**2+(8.0D0*pi/eta*im-p1/eta)**2)
  ! call pair_gai(p1,p0,vext,mu,re,im)
  ! tmatrix = 1.0D0/((1.0D0+8.0D0*pi/eta*re)**2+(8.0D0*pi/eta*im)**2)


 ! tmatrix = 1.0D0/(1.0D0+(x1/eta)**2*t/2.0D0)
 f(1) = a**2*b**2*c*x0**2*x1**7*tmatrix*fu*(1.0D0+y**2+y1**2-3.0D0*y**2*y1**2)*r**2
! write (11,*) x0,x1,r,y,y1
 !write (11,*) f(1)
 integrated = 0
 end function integrated
end subroutine rate

Subroutine pair0(mu0,vext,k0,kr,re,im)
  Use gsld
  use formidconstant
  Implicit None
  real*8,intent(in):: mu0,vext,k0,kr
  real*8,intent(out):: re,im
  Real*8:: a, b, a1, b1, answer
  real*8:: minl,maxl,r0,step,minl1,maxl1,mu
  Integer :: i, j,n0,i3

  mu=mu0-vext

  a1 = -1.0D0
  b1 = 1.0D0

  if (kr<1D2) then
      minl = 0.0D0
      maxl = kr-1D-4
      n0 = 40
      step = (maxl-minl)/dble(n0)
      a = minl
      b = minl + step
      r0 = 0.0D0
      do i3 = 1,n0
        answer = 0
        Do j = 1, nk
          Do i = 1, nk
           answer = answer + ak1(i)*ak1(j)*y1((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,k0,kr)
          End Do
        End Do
        answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
        r0 = r0+answer
        a = b
        b = b+step
      end do
      re=r0

      minl = kr+1D-4
      maxl = 2.0D0*kr
      n0 = 40
      step = (maxl-minl)/dble(n0)
      a = minl
      b = minl + step
      r0 = 0.0D0
      do i3 = 1,n0
        answer = 0
        Do j = 1, nk
          Do i = 1, nk
           answer = answer + ak1(i)*ak1(j)*y1((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,k0,kr)
          End Do
        End Do
        answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
        r0 = r0+answer
        a = b
        b = b+step
      end do
      re=re+r0

      minl = 2.0D0*kr
      maxl = 1.0D2
      n0 = 100
      step = (maxl-minl)/dble(n0)
      a = minl
      b = minl + step
      r0 = 0.0D0
      do i3 = 1,n0
        answer = 0
        Do j = 1, nk
          Do i = 1, nk
           answer = answer + ak1(i)*ak1(j)*y1((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,k0,kr)
          End Do
        End Do
        answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
        r0 = r0+answer
        a = b
        b = b+step
      end do
      re=re+r0
  else
    minl = 0.0D0
    maxl = 1.0D2
    n0 = 100
    step = (maxl-minl)/dble(n0)

        a = minl
        b = minl + step
        r0 = 0.0D0
      do i3 = 1,n0
        answer = 0
        Do j = 1, nk
          Do i = 1, nk
           answer = answer + ak1(i)*ak1(j)*y1((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),(a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(j),mu,k0,kr)
          End Do
        End Do
        answer = answer*(b-a)/2.0D0*(b1-a1)/2.0D0
        r0 = r0+answer
        a = b
        b = b+step
      end do
      re=r0
  end if
!------ imaginary part --------
    answer = 0
      Do i = 1, nk
       answer = answer + ak1(i)*y2((a1+b1)/2.0D0+(b1-a1)/2.0D0*fn1(i),mu,k0,kr)
      End Do
    answer = answer*(b1-a1)/2.0D0
    im=answer
!------------------------------
  re=re*0.025330296D0
  im=im*0.025330296D0
contains
  Function y1(q,the,mu,k0,kr)
    Implicit none
    Real*8:: y1
    Real*8:: k0,kr,the,mu,ep1,ep2,ferm1,ferm2,q,eps
    complex*16::mid
    eps = 1.0D-6
    ep1 = k0**2/4.0D0+q**2+q*k0*the-mu
    ep2 = k0**2/4.0D0+q**2-q*k0*the-mu
    ferm1=1.0D0/(exp(ep1/t)+1.0D0)
    ferm2=1.0D0/(exp(ep2/t)+1.0D0)
    mid=cmplx(q**2*(-ferm1-ferm2),0.0D0)/cmplx(2.0D0*kr**2-2.0D0*q**2,eps)!+cmplx(q**2,0.0D0)/cmplx(2.0D0*q**2,-eps)
    y1 = real(mid)
  End Function y1
  Function y2(the,mu,k0,kr)
    Implicit none
    Real*8:: y2
    Real*8:: k0,kr,the,mu,ferm1,ferm2,en1,en2
    en1 = k0**2/4.0D0+kr**2+kr*k0*the-mu
    en2 = k0**2/4.0D0+kr**2-kr*k0*the-mu
    ferm1=1.0D0/(exp(en1/t)+1.0D0)
    ferm2=1.0D0/(exp(en2/t)+1.0D0)
    y2 = 0.7853981634D0*kr*(ferm1+ferm2)
  End Function y2
end subroutine pair0


!give a pair popergater for cooper pair
!q is the momentent for center of mass
!k is the momentent for two partile
!vext is potential, mu is chemical potentail, t is temperture
!re, im is the real part and imagine part of pair popergater
subroutine pair(k,q,vext,mu,re,im)
  use formidconstant
  implicit none
real*8,intent(in):: k,q,vext,mu
real*8,intent(out):: re,im
integer ndim,ncomp,flags,mineval,maxeval,key,nregions,neval,fail
real*8 userdata(6),epsrel,epsabs,integral(2),error(2),prob(2)
integer seed,nnew,nvec
real*8 flatness
integer*8 spin
character*(*) statefile
parameter (statefile = "")
parameter (spin = -1)
parameter (nvec = 1)
parameter (ndim = 2)
parameter (ncomp = 2)
parameter (epsrel = 1.0D-4)
parameter (epsabs = 1.0D-4)
parameter (flags = 0)
parameter (mineval = 10000)
parameter (maxeval = 40000)
parameter (key = 13)
parameter (seed = 0)
parameter (nnew = 50)
parameter (flatness = 1.0D2)

userdata(1) = k
userdata(2) = q
userdata(3) = vext
userdata(4) = mu
! userdata(5) = t
! userdata(6) = eta

 call cuhre(ndim,ncomp,reandim,userdata,nvec,epsrel,epsabs,flags,mineval,maxeval,key,statefile,spin,nregions,neval,fail,integral,error,prob)
 re = integral(1)
 im = integral(2)

 contains
 function reandim(ndim,x,ncomp,f,userdata)
   use formidconstant
   implicit none
 integer,intent(in):: ndim,ncomp
 real*8,intent(in):: x(ndim),userdata(6)
 real*8,intent(out):: f(2)
 real*8 en1,en2,a,b,q1,the,eps,phi
 real*8 ferm1,ferm2
 real*8 k,q,vext,mu
 integer reandim
 a = 1.0D3
 b = 2.0D0
 q1 = x(1)*a
 the = x(2)*b-1.0D0
 k = userdata(1)
 q = userdata(2)
 vext = userdata(3)
 mu = userdata(4)

 ! k=k*sqrt(t/2.0D0)
 ! q=q*sqrt(t/2.0D0)

 eps = 1.0D-6

 !real part

 en1 = q**2/4.0D0+q1**2+q*q1*the
 en2 = q**2/4.0D0+q1**2-q*q1*the

 ferm1 = fermion(en1,vext,mu)
 ferm2 = fermion(en2,vext,mu)
 f(1) = real(cmplx(q1**2*(1.0D0-ferm1-ferm2),0.0D0)/cmplx(2.0D0*k**2-2.0D0*q1**2,eps)+cmplx(q1**2,0.0D0)/cmplx(2.0D0*q1**2,-eps))
 f(1) = a*b*f(1)*0.02533029591D0

 !imagine part

 en1 = q**2/4.0D0+k**2+q*k*the
 en2 = q**2/4.0D0+k**2-q*k*the

 ferm1 = fermion(en1,vext,mu)
 ferm2 = fermion(en2,vext,mu)
 f(2) = -0.7853981634D0*k*(1.0D0-ferm1-ferm2)!-Pi/4
 f(2) = b*f(2)*0.02533029591D0!2Pi/(2Pi)^3


 reandim = 0
 end function reandim

 !en is the momentent energy
 function fermion(en,vext,mu)
   use formidconstant
   implicit none
 real*8,intent(in):: en,vext,mu
 real*8 fermion
 fermion = 1.0D0/(exp((en-mu+vext)/t)+1.0D0)
 end function fermion

end subroutine pair

subroutine pair_gai(k,q,vext,mu,re,im)
  use formidconstant
  implicit none
real*8,intent(in):: k,q,vext,mu
real*8,intent(out):: re,im
real*8 userdata(6)
real*8 c1,c2,c3,cp,cm,minl,maxl,step,a,b,sun(2),r0(2)
real*8 k1,k2,k3,f1(2),f2(2),f3(2)
integer n,i
 parameter  (c1=0.5555555556D0)
 parameter  (c2=0.8888888889D0)
 parameter  (c3=0.7745966692D0)

 userdata(1) = k
 userdata(2) = q
 userdata(3) = vext
 userdata(4) = mu

 sun = 0.0D0
 r0 = 0.0D0

 minl = 0.0D0
 maxl = 1.0D-3
 n = 80
 step = (maxl-minl)/dble(n)
 a = minl
 b = minl + step
 do i = 1,n
   cp = (a+b)/2.0D0
   cm = (a-b)/2.0D0
   k1 = cp - cm*c3
   k2 = cp
   k3 = cp + cm*c3

   call intx2(k1,f1,userdata)
   call intx2(k2,f2,userdata)
   call intx2(k3,f3,userdata)
   sun = -cm*(c1*f1+c2*f2+c1*f3)
   r0 = r0+sun
   a = b
   b = b+step
 enddo

 minl = 1.0D-3
 maxl = 10.0D-3
 n = 25
 step = log(maxl/minl)/dble(n)
 a = minl
 b = minl*exp(step)
 do i = 1,n
   cp = (a+b)/2.0D0
   cm = (a-b)/2.0D0
   k1 = cp - cm*c3
   k2 = cp
   k3 = cp + cm*c3

   call intx2(k1,f1,userdata)
   call intx2(k2,f2,userdata)
   call intx2(k3,f3,userdata)
   sun = -cm*(c1*f1+c2*f2+c1*f3)
   r0 = r0+sun
   a = b
   b = b*exp(step)
 enddo

  minl = 10.0D-3
  maxl = 1.0D0
  n = 35
  step = log(maxl/minl)/dble(n)
  a = minl
  b = minl*exp(step)
  do i = 1,n
    cp = (a+b)/2.0D0
    cm = (a-b)/2.0D0
    k1 = cp - cm*c3
    k2 = cp
    k3 = cp + cm*c3

    call intx2(k1,f1,userdata)
    call intx2(k2,f2,userdata)
    call intx2(k3,f3,userdata)
    sun = -cm*(c1*f1+c2*f2+c1*f3)
    r0 = r0+sun
    a = b
    b = b*exp(step)
  enddo
  re = r0(1)
  im = r0(2)
 contains

 subroutine intx2(x,f,userdata)
 real*8,intent(in):: x,userdata(6)
 real*8,intent(out):: f(2)
 real*8 c1,c2,c3,cp,cm,minl,maxl,step,a,b,sun(2),r0(2)
 real*8 k1(2),k2(2),k3(2),f1(2),f2(2),f3(2)
 integer n,i
 parameter  (c1=0.5555555556D0)
 parameter  (c2=0.8888888889D0)
 parameter  (c3=0.7745966692D0)

 minl = 0.0D0
 maxl = 1.0D0
 n = 50
 step = (maxl-minl)/dble(n)
 a = minl
 b = minl + step

 sun = 0.0D0
 r0 = 0.0D0
 k1(1) = x
 k2(1) = x
 k3(1) = x
 do i = 1,n
 cp = (a+b)/2.0D0
 cm = (a-b)/2.0D0
 k1(2) = cp - cm*c3
 k2(2) = cp
 k3(2) = cp + cm*c3

 call reandim(k1,f1,userdata)
 call reandim(k2,f2,userdata)
 call reandim(k3,f3,userdata)
 sun = -cm*(c1*f1+c2*f2+c1*f3)
 r0 = r0+sun
 a = b
 b = b+step
 enddo
 f = r0

 end subroutine intx2

 subroutine reandim(x,f,userdata)
 real*8,intent(in):: x(2),userdata(6)
 real*8,intent(out):: f(2)
 real*8 en1,en2,a,b,q1,the,eps
 real*8 ferm1,ferm2
 real*8 k,q,vext,mu
 a = 1.0D3
 b = 2.0D0
 q1 = x(1)*a
 the = x(2)*b-1.0D0
 eps = 1.0D-6

 k = userdata(1)
 q = userdata(2)
 vext = userdata(3)
 mu =  userdata(4)

 ! k=k*sqrt(t/2.0D0)
 ! q=q*sqrt(t/2.0D0)

 !real part

 en1 = q**2/4.0D0+q1**2+q*q1*the
 en2 = q**2/4.0D0+q1**2-q*q1*the

 ferm1 = fermion(en1,vext,mu)
 ferm2 = fermion(en2,vext,mu)
 f(1) = real(cmplx(q1**2*(1.0D0-ferm1-ferm2),0.0D0)/cmplx(2.0D0*k**2-2.0D0*q1**2,eps)+cmplx(q1**2,0.0D0)/cmplx(2.0D0*q1**2,-eps))
 f(1) = a*b*f(1)*0.02533029591D0

 !imagine part

 en1 = q**2/4.0D0+k**2+q*k*the
 en2 = q**2/4.0D0+k**2-q*k*the

 ferm1 = fermion(en1,vext,mu)
 ferm2 = fermion(en2,vext,mu)
 f(2) = -0.7853981634D0*k*(1.0D0-ferm1-ferm2)!-Pi/4
 f(2) = b*f(2)*0.02533029591D0!2Pi/(2Pi)^3

 end subroutine reandim

 !en is the momentent energy
 function fermion(en,vext,mu)
 real*8,intent(in):: en,vext,mu
 real*8 fermion
 fermion = 1.0D0/(exp((en+vext-mu)/t)+1.0D0)
 end function fermion
end subroutine pair_gai

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine intd1(mu,r)
  use formidconstant
  implicit none
real*8,intent(in):: mu
real*8,intent(out):: r
integer ndim,ncomp,flags,mineval,maxeval,key,nregions,neval,fail
real*8 userdata(3),epsrel,epsabs,integral(1),error(1),prob(1)
integer seed,nnew,nvec
integer*8 spin
real*8 flatness
character*(*) statefile
parameter (statefile = "")
parameter (nvec = 1)
parameter (ndim = 2)
parameter (ncomp = 1)
parameter (epsrel = 1.0D-5)
parameter (epsabs = 1.0D-5)
parameter (flags = 6)
parameter (mineval = 5000)
parameter (maxeval = 500000)
parameter (key = 9)
parameter (seed = 3)
parameter (nnew = 500)
parameter (flatness = 5.0D2)
parameter (spin = -1)

userdata(2) = mu

call cuhre(ndim,ncomp,f_rr,userdata,nvec,epsrel,epsabs,flags,mineval,maxeval,key,statefile,spin,nregions,neval,fail,integral,error,prob)

r = integral(1)
contains
  function f_rr(ndim,x,ncomp,f,userdata)
    use formidconstant
    implicit none
  integer f_rr
  integer,intent(in):: ndim,ncomp
  real*8,intent(in):: x(ndim),userdata(3)
  real*8,intent(out):: f(ncomp)
  real*8 mu,vext,a,b,r,k,pf

  mu = userdata(2)

  a = sqrt(40.0D0*t/omega0**2)
  b = 1.0D1
  r = a*x(1)
  k = b*x(2)

  vext = 0.25D0*omega0**2*r**2
  pf = 1.0D0/(exp((k-mu+vext)/t)+1.0D0)
  f(1) = pf*k**1.5D0*r**2*a*b

  f_rr = 0
  end function f_rr
end subroutine intd1

subroutine def_tau(mu,ub,gfb,tau12)
!define tau22 and tau12, and give value of uf,ub,ufb,gff,gfb,gbb,v,t,u
use formidconstant
implicit none
real*8,intent(in):: mu,ub,gfb
real*8,intent(out):: tau12
real*8 eb,r,c
integer ndim,ncomp,flags,mineval,maxeval,key,nregions,neval,fail
real*8 userdata(7),epsrel,epsabs,integral(1),error(1),prob(1)
integer seed,nnew,nvec
real*8 flatness
character*(*) statefile
integer*8 spin
parameter (spin = -1)
parameter (statefile = "")
parameter (nvec = 1)
parameter (ndim = 6)
parameter (ncomp = 1)
parameter (epsrel = 1.0D-5)
parameter (epsabs = 1.0D-5)
parameter (flags = 6)
parameter (mineval = 5000000)
parameter (maxeval = 10000000)
parameter (key = 13)
parameter (seed = 3)
parameter (nnew = 500)
parameter (flatness = 5.0D2)

call intd1(mu,r)
c = 1.0D0/(60.0D0*t*9.869604401D0*r)                 !pi^2=9.869604401
c = c*16.0D0
!write(11,*) "c=",c


  eb = 2.0D0*mu-ub
  userdata(1) = mu
  userdata(6) = gfb
  userdata(7) = eb
  if (t<0.3D0) then
    call suave(ndim,ncomp,value_tau12,userdata,nvec,epsrel,epsabs,flags,seed,mineval,maxeval,nnew,flatness,statefile,spin,nregions,neval,fail,integral,error,prob)
    integral(1) = integral(1)*c/(320.0D0*2.0D0)
    !write (14,*) integral(1)
    integral(1) = max(abs(integral(1)),1.0D-7)
    tau12 = 1.0D0*omega0/integral(1)
  else
    tau12 = 1.0D7
  endif

 contains
  function value_tau12(ndim,x,ncomp,f,userdata)
    use formidconstant
    implicit none
  integer value_tau12
  integer,intent(in):: ndim,ncomp
  real*8,intent(in):: x(ndim),userdata(7)
  real*8,intent(out):: f(ncomp)
  real*8 mu,beta,gfb
  real*8 a,b,c,r,po,pr,y1,y2,vext,pd,pr0,pmid,mid,midpd
  real*8 p1,p2,p3,p4
  real*8 aa,bb,cc,k1,k2,k3,k4,k5,k6,k7,up,down,mide,gap1,gap2,gap3,gap4,gap5
   mu = userdata(1)
   gfb = userdata(6)
   eb = userdata(7)

  beta = 1.0D0/t
  a = 0.4D1
  b = 2.0D0
  c = sqrt(40.0D0/(omega0**2*beta))

  po = x(1)*a
  r = x(6)*c

  if (po**2+eb .lt. po**2/3.0D0) then
    f(1) = 0.0D0
    return
  else


  mide = po**2+eb
  up = (po+sqrt(6.0D0*eb+4.0D0*po**2))/3.0D0
  down = (po-sqrt(6.0D0*eb+4.0D0*po**2))/3.0D0
  gap1 = up - down
  k1 = down + x(2)*gap1
   if (mide .lt. 0.0D0) then
   write (*,*) mide,'11'
   pause
   else
   endif
   if (6.0D0*eb+4.0D0*po**2 .lt. 0.0D0) then
   write (*,*) mide,'12'
   pause
   else
   endif

  mide = po**2+eb-k1**2
  aa = 2.0D0
  bb = 2.0D0*k1-2.0D0*po
  cc = po**2+k1**2-2.0D0*po*k1-mide
  midpd = max(bb**2-4.0D0*aa*cc,0.0D0)
  up = (-bb+sqrt(midpd))/(2.0D0*aa)
  down = (-bb-sqrt(midpd))/(2.0D0*aa)
  gap2 = up - down
  k2 = down +x(3)*gap2
   if (mide .lt. 0.0D0) then
   write (*,*) mide,'21'
   pause
   else
   endif

  mide = max(po**2+eb-k1**2-k2**2-(po-k1-k2)**2,0.0D0)
  up = sqrt(mide/2.0D0)
  down = 0.0D0
  gap3 = up - down
  k3 = down +x(4)*gap3
  if (mide .lt. 0.0D0) then
   write (*,*) mide,'31'
   pause
   else
   endif

  aa = 2.0D0
  bb = -2.0D0*k3
  cc = 2.0D0*k3**2-mide

  up = (-bb+sqrt(bb**2-4.0D0*aa*cc))/(2.0D0*aa)
  down =(-bb-sqrt(bb**2-4.0D0*aa*cc))/(2.0D0*aa)
  gap4 = up - down
  k4 = down+ x(5)*gap4

   if (bb**2-4.0D0*aa*cc .lt. 0.0D0) then
   write (*,*) mide,'42'
   pause
   else
   endif

  k5 = abs(k3 - k4)

  mide = max(mide - k3**2 - k4**2 - k5**2,0.0D0)
  k6 = sqrt(mide/2.0D0)
  k7 = sqrt(mide/2.0D0)
  if (mide .lt. 0.0D0) then
   write (*,*) mide,'51'
   pause
   else
   endif

  vext = 0.25D0*omega0**2*r**2


    p1 = beta*(po**2+vext-mu)
    p2 = beta*((k1**2+k3**2)+vext-mu)
    p3 = beta*((k2**2+k4**2+k6**2)+vext-mu)
    p4 = beta*(((po-k1-k2)**2+k5**2+k7**2)+vext-mu)
    mid = exp(-(p2+p3+p4))/((1.0D0+exp(-p1))*(1.0D0+exp(-p2))*(1.0D0+exp(-p3))*(1.0D0+exp(-p4)))
    mid = mid*a*c*gap1*gap2*gap3*gap4*0.0000206934D0                                !!!!!!!!!!!!!0.0000206934=16/(2pi)^7/2
    f(1) = mid*gfb**2*po**2*k3*r**2*2.0D0*3.141592654D0
    f(1) = f(1)*(3.0D0*po**2-8.0D0*(po**2*eb)/9.0D0)!*4.0D0*3.141592654D0

  endif

  value_tau12 = 0
  end function value_tau12
end subroutine  def_tau



subroutine rate0(mu,tau22)
     use gsld
  use formidconstant
  implicit none
real*8,intent(in):: mu
real*8,intent(out):: tau22
integer ndim,ncomp,flags,mineval,maxeval,key,nregions,neval,fail
real*8 userdata(4),epsrel,epsabs,integral(1),error(1),prob(1),r,ero,sanse,c
integer seed,nnew,nvec
real*8 flatness,pi
character*(*) statefile
integer*8 spin
parameter (pi=3.141592653589793238462643383279502884197D0)
parameter (spin = -1)
parameter (statefile = "")
parameter (nvec = 1)
parameter (ndim = 4)
parameter (ncomp = 1)
parameter (epsrel = 1.0D-5)
parameter (epsabs = 1.0D-5)
parameter (flags = 6)
parameter (mineval = 5000)
parameter (maxeval = 500000)
parameter (key = 13)
parameter (seed = 0)
parameter (nnew = 1400)              !nnew = 1000          for vacuum tmatrix
parameter (flatness = 80.0D0)      !flatness = 50.0D0

userdata(1) = mu
 Call fn0(fn1, ak1)
 call DDD(mu,r)
 call suave(ndim,ncomp,integrated,userdata,nvec,epsrel,epsabs,flags,seed,mineval,maxeval,nnew,flatness,statefile,spin,nregions,neval,fail,integral,error,prob)
 sanse = integral(1)
 ero = error(1)

 tau22 = sqrt(2.0D0)*(5.0D0*pi)*(eta/t)**2*r/sanse

 !*************************************
 !
 !*************************************
 contains
 function integrated(ndim,x,ncomp,f,userdata)
   use formidconstant
   implicit none
 integer integrated
 integer,intent(in):: ndim,ncomp
 real*8,intent(in):: x(ndim),userdata(4)
 real*8,intent(out):: f(1)
 real*8 mu,z,vext,a,b
 real*8 x0,x1,r,y,y1,fu,p0,p1
 real*8 x12,x22,x32,x42
 real*8 re,im,tmatrix
 real*8 dd

 mu = userdata(1)

 a = 1.0D2
 b = 2.0D0                          !beware of this change!!!!!!!
 ! c = 5D0*sqrt(2D0*t/omega0)
 x0 = x(1)*a
 x1 = x(2)*a

 y = x(3)*b-1.0D0
 y1 = x(4)*b-1.0D0
 ! r = x(5)*c

 !x0 = sqrt(0.5D0*p0**2/t)
 !x1 = sqrt(2.0D0*p1**2/t)

 ! vext = 0.5D0*omega0*r**2!mu*r**2!
 z = exp(mu/t)
 x12 = 0.5D0*(x0**2+x1**2+2.0D0*x0*x1*y)
 x22 = 0.5D0*(x0**2+x1**2-2.0D0*x0*x1*y)
 x32 = 0.5D0*(x0**2+x1**2+2.0D0*x0*x1*y1)
 x42 = 0.5D0*(x0**2+x1**2-2.0D0*x0*x1*y1)
 fu = z**2*exp(-x0**2-x1**2)/((1.0D0+z*exp(-x12))*(1.0D0+z*exp(-x22))*(1.0D0+z*exp(-x32))*(1.0D0+z*exp(-x42)))
 !0.07957747155 = 1/(4Pi)

  p1=x1*sqrt(t/2.0D0)
  p0=x0*sqrt(t*2.0D0)
 call pair0(mu,0D0,p0,p1,re,im)
 tmatrix = 1.0D0/((1.0D0+8.0D0*pi/eta*re)**2+(8.0D0*pi/eta*im-p1/eta)**2)

 ! tmatrix = 1.0D0/(1.0D0+(x1/eta)**2*t/2.0D0)
 f(1) = a**2*b**2*x0**2*x1**7*tmatrix*fu*(1.0D0+y**2+y1**2-3.0D0*y**2*y1**2)
! write (11,*) x0,x1,r,y,y1
 !write (11,*) f(1)
 integrated = 0
 end function integrated

 subroutine DDD(mu,fun)
   use formidconstant
   use gsld
   implicit none
   real*8,intent(in):: mu
   real*8,intent(out):: fun
   Real*8:: a,b,answer
   real*8:: minl,maxl,r0,step
   Integer :: i, j,n0,i3

   minl = 0.0D0
   maxl = 1.0D2

   n0 = 1000
   step = (maxl-minl)/dble(n0)

       a = minl
       b = minl + step
       r0 = 0.0D0
       do i3 = 1,n0
         answer = 0
           Do i = 1, nk
     	      answer = answer + ak1(i)*y((a+b)/2.0D0+(b-a)/2.0D0*fn1(i),mu)
           End Do
         answer = answer*(b-a)/2.0D0
         r0 = r0+answer
         a = b
         b = b+step
       end do
   fun=0.5D0*r0
 end subroutine DDD

 Function y(x,mu)
    Implicit none
    Real*8:: y
    Real*8:: x,mu
    y=x**1.5D0/(exp(x-mu/t)+1.0D0)
  End Function y

end subroutine rate0
