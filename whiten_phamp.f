c ==========================================================
c Function filter4. Broadband filreting.
c ==========================================================
c Parameters for filter4 function:
c Input parameters:
c f1,f2   - low corner frequences, f2 > f1, Hz, (double)
c f3,f4   - high corner frequences, f4 > f3, Hz, (double)
c npow    - power of cosine tapering,  (int)
c dt      - sampling rate of seismogram in seconds, (double)
c n       - number of input samples, (int)
c seis_in - input array length of n, (float)
c Output parameters:
c seis_out - output array length of n, (float)
c ==========================================================
      subroutine filter4(dt,n,fre,seis1,seis2,seis3,
     &   seis_out1,seis_out2,seis_out3)
      parameter(nmax=10000000)
      real*4    fre(4),dt,f1,f2,f3,f4
      real*4    seis1(nmax),seis2(nmax),seis3(nmax),a
      real*4    seis_out1(nmax),seis_out2(nmax),seis_out3(nmax)
      real*4    dom
      integer*4   k,ns,nk,i,npow,n
      complex czero,s1(nmax),s2(nmax),s3(nmax),sf(nmax)
      czero = (0.0,0.0)
      f1=fre(1)
      f2=fre(2)
      f3=fre(3)
      f4=fre(4)
      npow=0
      ns=1
      do while(ns.lt.n)
         ns=ns*2
         npow=npow+1
      enddo
      dom = 1.0/dt/ns
c      write(*,*)'ns=',ns,'dom=',dom,'npow=',npow,'2**npow=',2**npow
      do k=1,ns
         s1(k)=czero
         s2(k)=czero
         s3(k)=czero
      enddo
      do k=1,n
         s1(k)=seis1(k)
         s2(k)=seis2(k)
         s3(k)=seis3(k)
      enddo
c    do fft
      call clogc(npow,s1,1,dt)
      call clogc(npow,s2,1,dt)
      call clogc(npow,s3,1,dt)
      nk = ns/2+1
      do k = nk+1,ns
        s1(k) = czero
        s2(k) = czero
        s3(k) = czero
      enddo
      s1(1) = s1(1)/2.0
      s2(1) = s2(1)/2.0
      s3(1) = s3(1)/2.0
      s1(nk) = cmplx(real(s1(n)),0.0)
      s2(nk) = cmplx(real(s1(n)),0.0)
      s3(nk) = cmplx(real(s1(n)),0.0)
c     do smoothing on sf equivalent to do " smooth mean h 20" in SAC
      call smooth(f1,f2,f3,f4,dom,nk,s1,20)
      call smooth(f1,f2,f3,f4,dom,nk,s2,20)
      call smooth(f1,f2,f3,f4,dom,nk,s3,20)
C=============================================================
c===============================================================
c   make tapering
      call flt4(f1,f2,f3,f4,dom,nk,npow,s1)
      call flt4(f1,f2,f3,f4,dom,nk,npow,s2)
      call flt4(f1,f2,f3,f4,dom,nk,npow,s3)
!      do i=1,nk
!         a=(real(s1(i))+real(s2(i))+real(s3(i)))/3.0
!         if (a.ne.0)then
!         s1(i)=cmplx(real(s1(i))/a,aimag(s1(i)))
!         s2(i)=cmplx(real(s2(i))/a,aimag(s2(i)))
!         s3(i)=cmplx(real(s3(i))/a,aimag(s3(i)))
!         endif
!      enddo
c       do i = 1,nk
c        seis_outamp(i)= real(dsqrt(dreal(sf(i))**2 +
c     1                        dimag(sf(i))**2))
c        seis_outph(i) = real(datan2(dimag(sf(i)),dreal(sf(i))))
c       enddo
c make forward FFT for seismogram: sf ==> s
!      call dfftw_plan_dft_1d(plan2,ns,sf,s,
!     *                         FFTW_FORWARD, FFTW_ESTIMATE)
!      call dfftw_execute(plan2)
!      call dfftw_destroy_plan(plan2)
c forming final result
!      call clogc(npow,s1,-1,dt)
!      call clogc(npow,s2,-1,dt)
!      call clogc(npow,s3,-1,dt)
!      do k = 1,n
!         seis_out1(k) = 2.0*real(real(s1(k)))/!ns
!         seis_out2(k) = 2.0*real(real(s2(k)))/ns
!         seis_out3(k) = 2.0*real(real(s3(k)))/ns
!      enddo
      return
      end
c===============================================================
c Tapering subroutine itself
c===============================================================
      subroutine flt4(f1,f2,f3,f4,dom,nk,npow,sf)
      parameter(nmax=1000000)
      real*4    f1,f2,f3,f4
      real*4    dom
      real*4    d1,d2,f,dpi,ss,s(nmax)
      integer   i,j
      integer   nk,npow
      complex   sf(nmax)
c ---
      dpi = atan(1.0)*4.0
      do i = 1,nk
         s(i) = 0.0
      enddo
      do i = 1,nk
        f = (i-1)*dom
        if(f.le.f1) then
          goto 1
        else if(f.le.f2) then
          d1 = dpi/(f2-f1)
          ss = 1.0
          do j = 1,npow
            ss = ss*(1-cos(d1*(f1-f)))/2.0
          enddo
          s(i) = ss
        else if(f.le.f3) then
           s(i) = 1.0
        else if(f.le.f4) then
          d2 = dpi/(f4-f3)
          ss = 1.0
          do j = 1,npow
            ss = ss*(1+cos(d2*(f3-f)))/2.0
          enddo
          s(i) = ss
        endif
  1     continue
      enddo
      do i = 1,nk
        sf(i) = sf(i)*s(i)
      enddo
      return
      end


c===============================================================
c  smoothing routine      call smooth(f1,f2,f3,f4,dom,nk,sf,20)
c=s==============================================================
      subroutine smooth(f1,f2,f3,f4,dom,nk,sf,number)
      parameter(nmax=1000000)
      integer  number,nk
      complex  sf(nmax)
      real*4   sorig(nmax), sout(nmax)
      real*4   f1,f2,f3,f4
      real*4 dom
      real*4   f,sum, avg
c ---
      do i = 1,nk
         sorig(i) = sqrt(real(sf(i))**2+aimag(sf(i))**2)
      enddo
      do i = 1,nk
         f = (i-1)*dom
        if( f .ge. f1 .and. f .le. f4 ) then
            sum = 0. 
          do jk = -number,number
             ijk = i+jk
             sum = sum + sorig(ijk)
          enddo
            sout(i) = sum/(2.*number+1.)
        else
            sout(i) = sorig(i)
        endif
      enddo
       do i = 1,nk
          f = (i-1)*dom
          if( f .ge. f1 .and. f .le. f4 ) then
              sout(i) = 1.0/sout(i)
          else
              sout(i) = 0.0
          endif
       enddo
       do i = 1,nk
           sf(i) = sf(i)*sout(i)
       enddo
       return
       end
