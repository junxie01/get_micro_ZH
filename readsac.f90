subroutine sactraces(name,sig,nlen,dt,beg,stla,stlo,nzhour,nzmin,nzsec,nzmsec,nerr)
parameter (nn=5000000)
real    sig(nn),dt,beg,stla,stlo
integer nzhour,nzmin,nzsec,nzmsec
integer nlen,nerr
character*180 name
!write(*,*)'name is ',name
call rsac1(name,sig,nlen,beg,dt,nn,nerr)
if (nerr.ne.0) then
       write(*,*) 'nerr=',nerr
       write(*,*) 'Error in reading sac file ',name
       return
endif
call getfhv('stla',  stla,  nerr)
call getfhv('stlo',  stlo,  nerr)
call getnhv('nzhour',nzhour,nerr)
call getnhv('nzmin', nzmin, nerr)
call getnhv('nzsec', nzsec, nerr)
call getnhv('nzmsec',nzmsec,nerr)
!write(*,*)'nerr=',nerr
return
end subroutine
