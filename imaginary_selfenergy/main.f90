! Code made by Shaozhi Li, 01/13/2023
! Maximum Entropy method to calculate the Eliashberg function from the imaginary part of the self-energy
module cluster
implicit none
integer, parameter:: runid=1
integer,parameter::L=120
integer,parameter:: momentum=0
integer momenta
integer band
double precision,parameter::beta=1.0d0/(6.0d0*0.08617)
double precision,dimension(0:L)::tau,green,green_err

integer,parameter:: num_w=160
double precision,parameter:: wmax=80.0d0
double precision,parameter::wmin=0.01d0
double precision,dimension(1:num_w)::freq
double precision,dimension(1:num_w):: dw
double precision,dimension(1:num_w)::spec_default
double precision,dimension(0:L,1:num_w)::ker

double precision:: debaifreq=20.0d0
double precision:: m0=0.5d0

integer,parameter::num_alpha=50
double precision,parameter::alphamin=-1.0d0
double precision,parameter::alphamax=2.0d0
double precision, dimension(1:num_alpha)::alphas

double precision,parameter::tol=1.0d-4
double precision:: mu=5000.0

double precision,dimension(1:num_w,1:num_alpha)::allspec
double precision, dimension(1:num_alpha)::allprob

double precision, dimension(1:num_w)::aver_spec
double precision ybase
end module cluster
!----------------------------------------------------------------
module Kernel
use cluster
implicit none
double precision,allocatable::xi(:,:),U(:,:),VT(:,:),UT(:,:),V(:,:),M(:,:)
double precision,dimension(1:L+1,1:L+1):: stdG
! matrix XI,U,VT,UT,V,M size "N"
integer N

contains
subroutine createkernel()
use cluster
implicit none
integer, parameter:: Nt=40000
integer i,j,LDA,lwork,ii,jj, k
double precision,dimension(1:num_w,0:L)::kerT
double precision,allocatable::xi1(:),work(:)
double precision,dimension(1:num_w,1:num_w)::U1
double precision,dimension(1:L+1,1:L+1)::VT1
double precision,dimension(1:num_w)::dw1
double precision deltaw
double precision y, yp, x, fy
double precision pi, w1, w2, fw1, fw2, nw1, sigma
double precision, dimension(1:Nt):: func
double complex eye
integer is,info

deltaw=(wmax-wmin)/dble(num_w-1)
do i=1,num_w
   freq(i)=wmin+deltaw*(i-1)
end do

dw=0.0d0
do i=1,num_w-1
  dw(i)=freq(i+1)-freq(i)
end do
dw=dw/2.0d0
dw1=0.0d0
do i=2,num_w
  dw1(i)=dw(i-1)
end do
dw(:)=dw(:)+dw1(:)

LDA=3*L
lwork=10*num_w
allocate(work(lwork))
is=min(L+1,num_w)
allocate(xi1(1:is))

!do i=0,L
!  do j=1,num_w
!    ker(i,j)=dexp(-tau(i)*freq(j))/(1+dexp(-beta*freq(j)))
!  end do
!end do
! kernel for real self energy
goto 4311
do i=0,L
 do j=1,num_w
    y=tau(i)*beta
    yp=freq(j)*beta
    do k=1,Nt
       x=0.02d0*k
       func(k)= 1.0d0/(2.0d0+dexp(-yp+y+x)+dexp(yp-y-x))+1.0d0/(2.0d0+dexp(-yp+y-x)+dexp( yp-y+x))-&
                1.0d0/(2.0d0+dexp(yp+y+x)+dexp(-yp-y-x))-1.0d0/(2.0d0+dexp( yp+y-x)+dexp(-yp-y+x))
       func(k)=func(k)*dlog(x)
    enddo
    fy=0.0d0
    do k=2,Nt-2,2
       fy=fy+0.02d0*4.0d0*func(k)+2.0d0*0.02*func(k+1)
    enddo
    fy=fy+0.02*(func(1)+func(Nt))
    ker(i,j)=fy/3.0d0
 enddo
enddo

! kernel for imaginary self energy
4311 continue
pi=2.0d0*asin(1.0d0)
sigma=0.0d0
do i=0, L
 w1=tau(i)*beta
 do j=1, num_w
    w2=freq(j)*beta
    fw1=1.0d0/(1.0d0+dexp(w1-w2))
    fw2=1.0d0/(1.0d0+dexp(w1+w2))
    nw1=1.0d0/(dexp(w2)-1.0d0)
    ker(i,j)=pi*(1.0d0-fw1+fw2+2.0d0*nw1)!*freq(j)
 enddo
enddo


do i=0,L
  do j=1,num_w
      kerT(j,i)=ker(i,j)
  end do
end do

call dgesvd('A','A',num_w,L+1,kerT,num_w,xi1,U1,num_w,VT1,L+1,work,lwork,info)

ii=0
do i=1,is
  if(xi1(i).gt.(xi1(1)*max(L+1,num_w)*2.22d-18)) ii=ii+1
end do

allocate(xi(1:ii,1:ii))
allocate(U(1:num_w,1:ii))
allocate(UT(1:ii,1:num_w))
allocate(V(1:L+1,1:ii))
allocate(VT(1:ii,1:L+1))
xi=0.0d0
U=0.0d0
VT=0.0d0
UT=0.0d0
V=0.0d0
do i=1,ii
  xi(i,i)=xi1(i)
  U(:,i)=U1(:,i)
  vt(I,:)=VT1(i,:)
end do

do i=1,ii
  do j=1,num_w
     UT(i,j)=U(j,i)
  end do
end do

do i=1,L+1
  do j=1,ii
     V(i,j)=VT(j,i)
  end do
end do         

stdG=0.0d0
do i=1,L+1
  stdG(i,i)=1.0d0/(green_err(i-1)**2)
end do

ALLOCATE(M(ii,ii))
M=matmul(matmul(matmul(matmul(xi,Vt),stdG),V),xi)
N=ii
deallocate(work)
deallocate(xi1)

return
end subroutine createkernel

end module Kernel
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
program main
use cluster
use Kernel
implicit none
integer i,j
character filename1*100,filename2*100
double precision, dimension(0:L):: fitgreen
double precision nw, xp
mu=5000.0d0

print*, "=============================================="
print*, "Maximum Entropy method to calculate the Eliashberg function from the imaginary part of the self-energy"
print*, "Code can only by used with the permission of Shaozhi Li"
print*, "Contact lishaozhiphys@gmail.com"
print*, "=============================================="


call init()
call readfile()
call createkernel()
print*, "finish kernel"
!call gaussian()
call set_default_model()
call get_all_spec()

xp=aver_spec(1)
aver_spec(1)=0.0d0

write(unit=filename1,fmt="('Eliashberg_',i6.6,'.dat')") runid
open(file=filename1,unit=60,action="write")
do i=1,num_w
  write(60,*) freq(i),aver_spec(i)
end do
close(60)

write(unit=filename2,fmt="(i2.2,'alpha_Prob.dat')") runid
open(file=filename2,unit=60,action="write")
do i=1,num_alpha
  write(60,*) alphas(i),allprob(i)
end do
close(60)

aver_spec(1)=xp
call cal_green(aver_spec,fitgreen)
write(unit=filename1,fmt="('back_selfenergy_',i6.6,'.dat')") runid
open(file=filename1,unit=60,action="write")
do i=0,L
  write(60,*) tau(i), fitgreen(i)+ybase
end do
close(60)


call dealloca()


end program main
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine init()
use cluster
implicit none
freq=0.0d0
spec_default=0.0d0
dw=0.0d0
allspec=0.0d0
allprob=0.0d0
tau=0.0d0
green=0.0d0
green_err=0.0d0
ker=0.0d0
aver_spec=0.0d0


return
end subroutine
!-----------------------------------------------------------------------------
subroutine readfile()
use cluster
implicit none
integer i,j
character filename*100
character line1*100
double precision, dimension(0:201):: x, y
!write(unit=filename,fmt="('alpha1.dat')")
if(runid.eq.1) filename="alpha_band_input.dat"
if(runid.eq.2) filename="beta_band_input.dat"
open(file=filename,unit=60,action='read')
!read(60,*) line1
do i=190,0,-1
  read(60,*) y(i),x(i)
end do
close(60)

if(runid.eq.2) then
   y(:)=y(:)*(-1.0d0)
endif

ybase=0.0d0
if(runid.eq.1) then
  do i=0, 9
     ybase=ybase+y(i)
  enddo
  ybase=ybase/10.0d0
else
   ybase=y(0)
endif

do i=0,L
   tau(i)=x(i+0)
   green(i)=abs(y(i+0)-ybase)*1000.0d0
   green_err(i)=2.0d0
enddo

ybase=ybase*1000.0d0

return 
end subroutine readfile
!------------------------------------------------------------------------------
! from spectral function calculate green function via
!   G(tau)=\int K(tau,w)A(w) dw
subroutine cal_green(specf,greena)
use cluster,only: num_w,freq,L,dw,ker
implicit none
double precision,dimension(0:L)::greena
double precision,dimension(1:num_w)::specf
integer i,j

greena=0.0d0
do i=0,L
 do j=1,num_w
    greena(i)=greena(i)+dw(j)*ker(i,j)*specf(j)
 end do
end do
return
end subroutine cal_green
!------------------------------------------------------------------------------
! calculate kai^2
double precision function cal_kai(greena)
use cluster, only: L, green, green_err
implicit none
double precision,dimension(0:L)::greena
integer i

cal_kai=0.0d0
do i=0,L
  cal_kai=cal_kai+(green(i)-greena(i))**2/(green_err(i)**2)
end do
!cal_kai=cal_kai/dble(L+1)
return
end function cal_kai
!------------------------------------------------------------------------------
subroutine gaussian()
use cluster,only:num_w,freq,spec_default
implicit none
integer i,j
double precision,parameter::mu=0.0d0
double precision,parameter::sig=4.0d0
double precision pi

pi=2.0d0*asin(1.0d0)
spec_default(:)=1.0d0/(dsqrt(2.0d0*pi)*sig)*dexp(-(freq(:)-mu)**2/2.0d0/(sig**2))
!spec_default(:)=1.0d0/40.0d0
return
end subroutine gaussian
!----------------------------------------------------------------------------------
subroutine set_default_model()
use cluster, only: runid, num_w, freq, spec_default, m0, debaifreq, beta
implicit none
integer i, j
double precision nw
character file1*100
double precision x

!do i=1,num_w
!  if(freq(i)<debaifreq) then
!     spec_default(i)=m0*(freq(i)/debaifreq)**2
!  else
!     spec_default(i)=m0
!  endif

!  if(i.eq.1) then
!     spec_default(i)=0.001*(freq(i)/debaifreq)**2/nw
!  endif
!enddo
if(runid.eq.1) then
   file1="real_alpha_default_spec.dat"
else if(runid.eq.2) then
   file1="real_beta_default_spec.dat"
endif

open(file=file1,unit=43,action='read')
do i=1, num_w
   read(43,*) x, spec_default(i)
enddo
close(43)
end subroutine
!----------------------------------------------------------------------------------
! calculate entropy S
! S=\int dw [A(w)-m(w)-A(w)ln(A(w)/m(w))]
double precision function cal_entropy(specf)
use cluster,only: num_w, freq, spec_default,dw
implicit none
double precision, dimension(1:num_w)::specf
integer i

cal_entropy=0.0d0
do i=1,num_w
 cal_entropy=cal_entropy+(specf(i)-spec_default(i)-specf(i)*dlog(specf(i)/spec_default(i)))*dw(i)
end do


return
end function cal_entropy
!------------------------------------------------------------------------------------
subroutine getspecf(alpha,specf)
use cluster,only: L, num_w, freq, dw,mu, spec_default,green,green_err,tol
use Kernel
implicit none
double precision, dimension(1:N,1:N):: R1,identity
double precision, dimension(1:L+1)::R2
double precision, dimension(1:N)::btemp,deltab,g
double precision, dimension(1:num_w)::specf
double precision, dimension(1:num_w,1:num_w)::diag_specf
double precision, dimension(0:L)::greena
double precision, dimension(1:N,1:N)::LHS
double precision, dimension(1:N)::RHS
integer, dimension(1:N)::IPIV
double precision, dimension(1:N):: work
double precision, dimension(1:N)::crit1
double precision, dimension(1:num_w)::al
double precision cal_kai, cal_entropy
integer INFO
double precision Qold, Qnew, criteria
double precision alpha
integer i,j,i1,j1
logical judgement

! inital values
identity=0.0d0
do i=1,N
 identity(i,i)=1.0d0
end do

btemp=0.0d0
deltab=0.0d0
 criteria=0.0d0

specf(:)=spec_default(:)
call cal_green(specf,greena)
Qold=0.5d0*cal_kai(greena)-alpha*cal_entropy(specf)
judgement=.true.


! use newton method find btemp
do while(judgement)
   diag_specf=0.0d0
   do i=1,num_w
     diag_specf(i,i)=specf(i)
   end do
   call cal_green(specf,greena)
   LHS=0.0d0
   RHS=0.0d0

   R1=matmul(UT, matmul(diag_specf, U))
   do i=1,L+1
    R2(i)=-1.0d0*(green(i-1)-greena(i-1))/(green_err(i-1)**2)
   end do
   g=matmul(Xi, matmul(VT,R2))

   RHS(:)=-alpha*btemp(:)-g(:)

   LHS=matmul(M, R1)
   LHS(:,:)=LHS(:,:)+(alpha+mu)*identity(:,:)
   CALL DGETRF(N, N, LHS, N, IPIV, INFO)
   call dgetri(N, LHS, N, IPIV, work, N, INFO)
   deltab=matmul(LHS,RHS)
   
   crit1=0.0d0
   crit1=matmul(deltab,R1)
   criteria=0.0d0
   do i=1,N
     criteria=criteria+crit1(i)*deltab(i)
   end do
   if(criteria<0.1*sum(spec_default(:))) then
      btemp(:)=btemp(:)+deltab(:)
      al=matmul(U,btemp)
      specf(:)=spec_default(:)*dexp(al(:))
      call cal_green(specf,greena)
      Qnew=0.5*cal_kai(greena)-alpha*cal_entropy(specf)
      if(abs(Qnew-Qold)/Qold.lt.tol) then
            judgement=.false.
            write(*,*) "alpha = ", alpha, "Q = ", Qnew
      end if
      Qold=Qnew
   else
      mu=mu*2.0d0
      specf(:)=spec_default(:)
      call cal_green(specf, greena)
      Qold=0.5d0*cal_kai(greena)-alpha*cal_entropy(specf)
      btemp=0.0d0
      write(*,*) "mu is ", mu
   end if
end do

!do i=1,num_w
!  print*, i, specf(i)
!enddo
!stop
return
end subroutine getspecf
!------------------------------------------------------------------
subroutine cal_cova(covariance)
use cluster, only: L,green_err
implicit none
double precision, dimension(0:L,0:L):: covariance
integer i,j

covariance=0.0d0
do i=0,L
     covariance(i,i)=green_err(i)**2
end do

return
end subroutine cal_cova
!-------------------------------------------------------------------
double precision function cal_prob(alpha,specf)
use cluster, only: num_w, L, ker
implicit none
double precision, dimension(1:num_w)::specf
double precision, dimension(0:L):: greena
double precision, dimension(1:num_w,0:L)::kerT
double precision, dimension(1:num_w,1:num_w)::mat_a,mat_b
double precision, dimension(0:L,0:L)::covariance
double precision, dimension(1:num_w)::eigen
double precision, dimension(1:3*num_w)::work
double precision cal_kai,cal_entropy,s1,s2,s3
double precision alpha
integer info
integer i,j

do i=0,L
  do j=1,num_w
      kerT(j,i)=ker(i,j)
  end do
end do
call cal_cova(covariance)
do i=0,L
 covariance(i,i)=1.0/covariance(i,i)
end do
mat_a=matmul(kert,matmul(covariance,ker))
mat_b=0.0d0
do i=1,num_w
  do j=1,num_w
     mat_b(i,j)=dsqrt(specf(i))*mat_a(i,j)*dsqrt(specf(j))
  end do
end do

call dsyev('N','U', num_w, mat_b, num_w, eigen, work, 3*num_w, INFO)
s1=1.0d0
do i=1,num_w
  s1=s1*alpha/(alpha+eigen(i))
end do
s1=dsqrt(s1)

call cal_green(specf,greena)
S2=0.5d0*cal_kai(greena)-alpha*cal_entropy(specf)
S2=dexp(-S2)
cal_prob=s1*s2/alpha
return
end function cal_prob
!-----------------------------------------------------------------------------------
subroutine get_all_spec()
use cluster, only: alphamin, alphamax, num_alpha, num_w,allprob,allspec,aver_spec,&
                   alphas
implicit none
double precision, dimension(1:num_w)::specf
double precision, dimension(1:num_alpha)::dalpha,dalpha_1
double precision sum_prob
double precision da,maxa,mina
double precision cal_prob
integer i,j

maxa=10.0**alphamax
mina=10.0**alphamin

da=(alphamax-alphamin)/dble(num_alpha-1)
do i=1,num_alpha
  alphas(i)=10**(alphamin+da*(i-1))
end do
do i=1,num_alpha-1
  dalpha(i)=alphas(i+1)-alphas(i)
end do
dalpha=dalpha/2.0d0

dalpha_1=0.0d0
do i=2,num_alpha
  dalpha_1(i)=dalpha(i-1)
end do
do i=1,num_alpha
  dalpha(i)=dalpha(i)+dalpha_1(i)
end do

do i=1,num_alpha
   call getspecf(alphas(i),specf)
   allprob(i)=cal_prob(alphas(i),specf)
   allspec(:,i)=specf(:)
   write(*,*) allprob(i),sum(specf)
end do !i

sum_prob=0.0d0
do i=1,num_alpha
  sum_prob=sum_prob+dalpha(i)*allprob(i)
end do
allprob(:)=allprob(:)/sum_prob

do i=1,num_w
  aver_spec(i)=0.0d0
  do j=1,num_alpha
     aver_spec(i)=aver_spec(i)+dalpha(j)*allspec(i,j)*allprob(j)
  end do
end do

return
end subroutine get_all_spec

subroutine dealloca()
use Kernel
implicit none

deallocate(xi)
deallocate(U)
deallocate(UT)
deallocate(V)
deallocate(VT)
deallocate(M)

end subroutine dealloca
