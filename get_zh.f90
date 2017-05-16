subroutine getzh(sig,fre,nsamp,dt,period,zh,dphi,phi_h,phi_z)
!     get Z/H ratio for three component microtremor at period
parameter (nmax=4000000)
integer i,j,ne,nb,k,marki,markii,nf_id
integer:: nsamp,nerr,shift,nn,npow,ndelta
real:: sig(nsamp,3),period,convdeg,x1,y1
real:: dt,begin_t,df,per(100),fre(4),ddf,dphi
real:: depmax1,depmax2,depmax3,depmax4,y(nmax)
real:: phi,pi,ii,imax,phi_h,phi_z,zh,f1,f2,zh1
complex:: iphi,rr,zero
complex:: sigc_z(nmax),sigc_n(nmax),sigc_e(nmax)
!ddf=1.0/period
Zero=cmplx(0.0,0.0)
pi=4.0*atan(1.0)
convdeg=pi/180.0
sigc_z=Zero
sigc_n=Zero
sigc_e=Zero
nn=2;npow=1
do while (nn<nsamp)
   nn=nn*2
   npow=npow+1
enddo
nk=nn/2+1
df=1.0/(nn*dt)
sigc_z(1:nsamp)=cmplx(sig(1:nsamp,1),0)
sigc_n(1:nsamp)=cmplx(sig(1:nsamp,2),0)
sigc_e(1:nsamp)=cmplx(sig(1:nsamp,3),0)
!      write(*,*)'hello'
! do fast fourier transform to these three component data
call clogc(npow,sigc_z,1,dt)
call clogc(npow,sigc_n,1,dt)
call clogc(npow,sigc_e,1,dt)
!call zfour(sigc_z,nn,1,dt,df)
!call zfour(sigc_n,nn,1,dt,df)
!call zfour(sigc_e,nn,1,dt,df)
!call flt4(fre(1),fre(2),fre(3),fre(4),df,nk,npow,sigc_z)
!call flt4(fre(1),fre(2),fre(3),fre(4),df,nk,npow,sigc_n)
!call flt4(fre(1),fre(2),fre(3),fre(4),df,nk,npow,sigc_e)

!call smooth(fre(1),fre(2),fre(3),fre(4),df,nk,sigc_z,20)
!call smooth(fre(1),fre(2),fre(3),fre(4),df,nk,sigc_n,20)
!call smooth(fre(1),fre(2),fre(3),fre(4),df,nk,sigc_e,20)
nf_id=int(1.0/period/df)+1         !     the id of frequency 1.0/period is nf
!     frequency band between f1 and f2
ndelta=int(0.01/df)
nb=nf_id-ndelta
ne=nf_id+ndelta
!write(*,*)'nb=',nb,'ne=',ne
!      write(*,*)'f1=',(nb-1)*df,'f2=',(ne-1)*df,'fc=',1/period,'ndelta=',ndelta
! considering the finite frequency effect, in the frequency band between f1 and f2
! write(*,*)"f1=",f1,"f2=",f2
! search for the horizontal angel.
!open(100,file='test.ii')
imax=0 
do k=1,1800  ! loop over all angle
   ii=0;phi=k*convdeg/10.0
   do j=nb,ne
      iphi=sigc_n(j)*cos(phi)+sigc_e(j)*sin(phi)
      ii=ii+cabs(iphi)*cabs(iphi)
   enddo
!      write(100,*)phi/convdeg,ii
   if(imax<ii)then     ! search for the maximum of integration
      marki=k
      imax=ii
   endif
enddo
!      close(100)
phi=marki*convdeg
!write(*,*)'phi0=',phi/convdeg
!nf=int(ddf/df)+1
y1=aimag(sigc_n(nf_id)*cos(phi)+sigc_e(nf_id)*sin(phi))
x1=real(sigc_n(nf_id)*cos(phi)+sigc_e(nf_id)*sin(phi))
!phi_h=(aimag(sigc_n(nf))*cos(phi)+aimag(sigc_e(nf))*sin(phi))/(real(sigc_n(nf))*cos(phi)+real(sigc_e(nf))*sin(phi))
phi_h=atan(y1/x1)/convdeg
if(x1==0.and.y1.ge.0)phi_h=0
if(x1==0.and.y1.lt.0)phi_h=180
if(x1<0)phi_h=270-phi_h
if(x1>0)phi_h=90-phi_h
y1=aimag(sigc_z(nf_id))
x1=real(sigc_z(nf_id))
phi_z=atan(y1/x1)/convdeg
if(x1==0.and.y1.ge.0)phi_z=0
if(x1==0.and.y1.lt.0)phi_z=180
if(x1<0)phi_z=270-phi_z
if(x1>0)phi_z=90-phi_z
rr=sigc_n(nf_id)*cos(phi)+sigc_e(nf_id)*sin(phi)
zh=abs(cabs(sigc_z(nf_id))/cabs(rr))
zh1=abs(real(sigc_z(nf_id))/real(rr))
dphi=phi_h-phi_z
if (dphi.lt.-180)then
   dphi=dphi+360
elseif (dphi.gt.180)then
   dphi=dphi-360
endif
return
end subroutine
