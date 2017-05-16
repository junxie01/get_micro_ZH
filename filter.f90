! ==========================================================
! Function filter4. Broadband filreting.
! ==========================================================
! Parameters for filter4 function:
! Input parameters:
! f1,f2   - low corner frequences, f2 > f1, Hz, (double)
! f3,f4   - high corner frequences, f4 > f3, Hz, (double)
! npow    - power of cosine tapering,  (int)
! dt      - sampling rate of seismogram in seconds, (double)
! n       - number of input samples, (int)
! seis_in - input array length of n, (float)
! Output parameters:
! seis_out - output array length of n, (float)
! ==========================================================
subroutine filter(sig,sigout,fre,dt,npts)
implicit none
integer*4 nmax
parameter(nmax=4000000)
integer*4 k,nk,npts,ns,npow
real*4    fre(4),dt,sig(nmax),sigout(nmax)
real*8    dom
complex   czero,s(nmax)
! ---
czero = (0.0,0.0)
ns=1
npow=0
do while(ns.lt.npts)
   ns=ns*2
   npow=npow+1
enddo
dom = dble(1.0/dt/ns)
s=czero
s(1:npts) = cmplx(sig(1:npts),0)
!write(*,*)' do fft'
call clogc(npow,s,1,dt)
nk = ns/2+1
s(nk+1:ns) = czero
!sf(1) = sf(1)/2.0d0
!sf(nk) = cmplx(real(sf(nk)),0.0d0)
!write(*,*)'taper the waveform'
call taperf(s,fre,dom,nk,npow)
!write(*,*)'taper the waveform done'
! make forward FFT for seismogram: sf ==> s
call clogc(npow,s,-1,dt)
sigout(1:npts)=real(s(1:npts))
return
end subroutine
