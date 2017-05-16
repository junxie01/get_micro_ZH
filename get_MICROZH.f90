! get microtremor Z/H ratio directly
! 2016/07/31 introduce sacio
! make sure the sampling rate of the three components are the same
program get_MICROZH
use sacio
implicit none
type(sac_head) :: sachead
integer    nmax,nstmax
parameter (nmax=4000000,nstmax=1000)
real sigo(nmax,3)
real sig(nmax,3),sig_s(nmax,3)
real zh,dphi,pmin,pmax,p1,p2
real periods(100),period,dhour,dmin
real stla,stlo,evla,evlo,fre(4),f1,f2
real phi_h,phi_z
character (180)dirinn,dirout,output1,input,list,output2
character (180)name1,name2,dir_day,dir,command
character (20)year_day,nd,pd
character (10)sta(nstmax)
character (2) net(nstmax)
character (180)sac(3)
character (3)com(3),co
integer year_b,year_e,day_b,day_e,ic,idd,ib,ie
integer iy,id,ih,ist,nsamp,is,np,ip
integer nst,i,nh,j,nnh,nlen,nmin
integer begday,endday,multpt,im
integer jday,npts,npt2,nerr,ns
logical ext
if (iargc().ne.1)then
   write(*,*)'Usage: get_micro_zh param.dat'
   write(*,*)'param.dat is like:'
   write(*,*)'station_list '
   write(*,*)'year_from day_from year_to day_to'
   write(*,*)'dhour,dmin,multpt (length of SAC file in hour. e.g. 2)'
   write(*,*)'p1 p2             (period band'
   write(*,*)'component         (e.g. BH)'
   write(*,*)'data_dir          (SAC file directory'
   write(*,*)'dir_out           (output directory'
   call exit(-1)
endif
call getarg(1,input)
open(10,file=input)                        ! read in control parameter
read(10,*)list
read(10,*)year_b,day_b,year_e,day_e
read(10,*)dhour,dmin,multpt
read(10,*)p1,p2
read(10,*)co
read(10,'(a180)')dirinn
read(10,'(a180)')dirout
close(10)
com(1)=trim(co)//'Z'
com(2)=trim(co)//'N'
com(3)=trim(co)//'E'
open(11,file=list)       ! read in station list
do i=1,nstmax
12 read(11,*,err=13,end=13) net(i),sta(i)
enddo
13 close(11)
nst=i-1
write(*,'("There are ",i0," stations")')nst

fre(2)=1.0/p2
fre(1)=fre(2)*0.98
fre(3)=1.0/p1
fre(4)=fre(3)*1.02
pmin=1.0/fre(4)
pmax=1.0/fre(1)

write(command,'("mkdir -p",1x,1a)')trim(dirout)
call system(command)
call getper(periods,100,pmin,pmax,np,p1,p2) ! get peirod
inquire(file='pfile',exist=ext)
if(.not.ext)then
   open(11,file='pfile')
   do i=1,np
      write(11,*)periods(i)
!      write(*,*)1.0/periods(i)
   enddo
   close(11)
endif
nh=int(24/dhour)                            ! number of segments per day
!nmin=int(dhour*60/dmin)                     ! number segments per segments
do ip=1,np                                  ! loop over period
   period=periods(ip)
   write(pd,'(i3.3)')ip
   write(pd,'(f7.1)')period
   pd=trim(adjustl(pd))//'_list'
   do is=1,nst                         ! loop over station
      output1=trim(dirout)//'/'//trim(sta(is))//'_'//trim(pd)
      open(20,file=output1)
      do iy=year_b,year_e         ! loop over year
         jday=365
         if(mod(iy,4).eq.0.and.mod(iy,100).ne.0.or.mod(iy,400).eq.0)jday=366
         endday=day_e
         if(iy.ne.year_e)endday=jday
         begday=day_b
         if(iy.ne.year_b)begday=1
         do id=begday,endday                      ! loop over day
            write(year_day,'(i0,"_",i3.3)')iy,id
            do ih=1,nh                       ! loop over daily segment
               nnh=int((ih-1)*dhour)
               write(nd,'("_",i2.2)')nnh
               do ic=1,3
                  !sac(ic)=trim(dirinn)//'/'//trim(year_day)//'/'//trim(year_day)//trim(nd)//'_'//&
!trim(sta(is))//'_'//com(ic)//'.SAC'
                  write(sac(ic),'(1a,"/",i4.4,"_",i3.3,"/",i4.4,"_",i3.3,"_",i2.2,"_",1a,"_",1a,"_",1a,".SAC")')&
trim(dirinn),iy,id,iy,id,nnh,net(is),trim(sta(is)),com(ic)
                  inquire(file=sac(ic),exist=ext)
                  write(*,*)period,trim(sac(ic))
                  if(.not.ext)then
                     idd=1000
                     exit
                  endif
                  call read_sachead(sac(ic),sachead,nerr)
                  if(nerr.eq.-1)then
                     idd=1000
                     exit
                  endif
                  call read_sac(sac(ic),sig(:,ic),sachead,nerr)
                  if(nerr.eq.-1)then
                     idd=1000
                     exit
                  endif
                  idd=ic
               enddo
               if(idd.ne.3)cycle          ! if all three component exist
               npts=int(dmin*60/sachead%delta)+1
               nmin=int((sachead%npts-npts))/int((1-multpt/100.0)/npts) +1
               do im=1,nmin      ! loop over segment in minutes 
                  ns=npts*(im-1)+1
                  ib=int((1-multpt/100.0)*npts)*(im-1)+1
                  ie=ib+npts
                  sig_s(1:npts,:)=sig(ib:ie,:)
                  call getzh(sig_s(1:npts,1:3),fre,npts,sachead%delta,period,zh,dphi,phi_h,phi_z)
                  write(20,'(2i4,5f9.2)')iy,id,nnh+(1-multpt/100.0)*dmin*(im-1)/60,zh,phi_h,phi_z,dphi
               enddo
            enddo       ! end loop daily segment
         enddo          ! end loop day
      enddo             ! end loop year
      close(20)
   enddo                ! end loop station
enddo                   ! loop over period
end program
