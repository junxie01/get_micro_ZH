      subroutine clogc(n,x,sigm,dt)
!subroutine clogc(n,x,sigm,dt)
!--- performs fft on signals with length 2**n and sampling interval
!--- of dt seconds (if in the time domain; notice that dt*df=1/2**n).
!--- the signal is stored in x. it may be complex.
!--- the spectrum is returned in x. it is almost always complex.
!--- a time-to-frequency transform is done with sign=+1. (conform
!--- the convention adopted in aki and richards - the alternative
!--- convention may be obtained by taking complex conjugates after
!--- the call to clogc).
!--- the normalization factor 1./twopi occurs in the frequency-to
!--- time transform (again aki&richards).
!--- normalization is such that physical dimensions are respected.
!--- thus, if the time signal is dimensioned in meters, the
!--- resulting spectral density in x is in meters/hz. for example,
!--- if the time signal is the unit sinc function of width dt, centered
!--- at t=0, the spectral density is dt for all values of the frequency.c
!--- array locations: if x contains the spectrum, it has the spectrum
!--- for positive frequencies in the first 2**n/2+1 elements, such that
!--- x(1) is at 0 hz, x(2) at df hertz, and x(2**n/2+1) at the nyquist,
!--- where df=1./(2**n*dt) and the nyquist is 1./(2*dt) hz.
!--- the second half of x contains the spectrum for negative frequencies
!--- such that x(2**n) is at -df, x(2**n-1) at -2*df hz etcetera.
!--- if x contains the time signal, x(1) is at time 0, x(2)
!--- at time dt etc.
!
      dimension x(1),m(25)
      complex x,wk,hold,q
      if(sigm.ge.0.) then
        sign=1.
      else
        sign=-1.
      endif
      lx=2**n
      do 1 i=1,n
    1    m(i)=2**(n-i)
      do 4 l=1,n
      nblock=2**(l-1)
      lblock=lx/nblock
      lbhalf=lblock/2
      k=0
      do 4 iblock=1,nblock
      fk=k
      flx=lx
      v=sign*6.283185308*fk/flx
      wk=cmplx(cos(v),sin(v))
      istart=lblock*(iblock-1)
      do 2 i=1,lbhalf
      j=istart+i
      jh=j+lbhalf
      q=x(jh)*wk
      x(jh)=x(j)-q
      x(j)=x(j)+q
    2 continue
      do 3 i=2,n
      ii=i
      if(k.lt.m(i)) go to 4
    3 k=k-m(i)
    4 k=k+m(ii)
      k=0
      do 7 j=1,lx
      if(k.lt.j) go to 5
      hold=x(j)
      x(j)=x(k+1)
      x(k+1)=hold
    5 do 6 i=1,n
      ii=i
      if(k.lt.m(i)) go to 7
    6 k=k-m(i)
    7 k=k+m(ii)
      if(sign.gt.0.) go to 9
      if(sign.lt.0)then
      flx=flx*dt
      do 8 i=1,lx
    8 x(i)=x(i)/flx
      return
    9 do 10 i=1,lx
   10 x(i)=x(i)*dt
      return
      endif
      end subroutine
