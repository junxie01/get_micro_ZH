module sacio
type :: sac_head
real(4):: delta
real(4):: depmin
real(4):: depmax
real(4):: scale
real(4):: odelta
real(4):: b
real(4):: e
real(4):: o
real(4):: a
real(4):: internal1
real(4):: t0
real(4):: t1
real(4):: t2
real(4):: t3
real(4):: t4
real(4):: t5
real(4):: t6
real(4):: t7
real(4):: t8
real(4):: t9
real(4):: f
real(4):: resp0
real(4):: resp1
real(4):: resp2
real(4):: resp3
real(4):: resp4
real(4):: resp5
real(4):: resp6
real(4):: resp7
real(4):: resp8
real(4):: resp9
real(4):: stla
real(4):: stlo
real(4):: stel
real(4):: stdp
real(4):: evla
real(4):: evlo
real(4):: evel
real(4):: evdp
real(4):: mag
real(4):: user0
real(4):: user1
real(4):: user2
real(4):: user3
real(4):: user4
real(4):: user5
real(4):: user6
real(4):: user7
real(4):: user8
real(4):: user9
real(4):: dist
real(4):: az
real(4):: baz
real(4):: gcarc
real(4):: internal2
real(4):: internal3
real(4):: depmen
real(4):: cmpaz
real(4):: cmpinc
real(4):: unused2
real(4):: unused3
real(4):: unused4
real(4):: unused5
real(4):: unused6
real(4):: unused7
real(4):: unused8
real(4):: unused9
real(4):: unused10
real(4):: unused11
real(4):: unused12
integer:: nzyear
integer:: nzjday
integer:: nzhour
integer:: nzmin
integer:: nzsec
integer:: nzmsec
integer:: internal4
integer:: internal5
integer:: internal6
integer:: npts
integer:: internal7
integer:: internal8
integer:: unused13
integer:: unused14
integer:: unused15
integer:: iftype
integer:: idep
integer:: iztype
integer:: unused16
integer:: iinst
integer:: istreg
integer:: ievreg
integer:: ievtyp
integer:: iqual
integer:: isynth
integer:: unused17
integer:: unused18
integer:: unused19
integer:: unused20
integer:: unused21
integer:: unused22
integer:: unused23
integer:: unused24
integer:: unused25
integer:: unused26
integer:: leven
integer:: lpspol
integer:: lovrok
integer:: lcalda
integer:: unused27
character(8)::  kstnm
character(16):: kevnm
character(8)::  khole
character(8)::  ko
character(8)::  ka
character(8)::  kt0
character(8)::  kt1
character(8)::  kt2
character(8)::  kt3
character(8)::  kt4
character(8)::  kt5
character(8)::  kt6
character(8)::  kt7
character(8)::  kt8
character(8)::  kt9
character(8)::  kf
character(8)::  kuser0
character(8)::  kuser1
character(8)::  kuser2
character(8)::  kcmpnm
character(8)::  knetwk
character(8)::  kdatrd
character(8)::  kinst
end type sac_head

end module sacio

subroutine read_sachead(fname,sachead,nerr)
use sacio
character*(*):: fname
type(sac_head):: sachead
integer:: funit,iso,nerr
funit=10
nerr=0
open(funit,file=fname,form='unformatted',access='direct',recl=158*4,iostat=iso)
if(iso/=0)then
   write(*,*)"Error in openning file ",trim(fname)
   write(*,*)"Stop the subroutine read_sachead"
   nerr=-1
   return
endif
read(funit,rec=1,iostat=iso)sachead
if(iso/=0)then
   write(*,*)"Error in read head value of file ",trim(fname)
   write(*,*)"Stop the subroutine read_sachead"
   nerr=-1
   return
endif
close(funit)
return
end subroutine read_sachead

subroutine read_sac(fname,array,sachead,nerr)
use sacio
integer:: nmax,funit,iso,nerr,npts,i
parameter(nmax=4000000)
character*(*):: fname
type(sac_head):: sachead
real,dimension(nmax):: array
funit=10
nerr=0
npts=sachead%npts
!write(*,*)'npts=',npts
open(funit,file=fname,form='unformatted',access='direct',recl=158*4+npts*4,iostat=iso)
if(iso/=0)then
   write(*,*)"Error in openning file ",trim(fname)
   write(*,*)"Stop the subroutine read_sac"
   nerr=-1
   return
endif
!do i=1,npts
!   read(funit,rec=158+i,iostat=iso)array(i)
read(funit,rec=1,iostat=iso)sachead,(array(i),i=1,npts)
if(iso/=0)then
   write(*,*)"Error in read file ",trim(fname)
   write(*,*)"Stop the subroutine read_sac"
   nerr=-1
   write(*,*)i,array(i)
   return
endif
!enddo
close(funit)
return
end subroutine read_sac

subroutine write_sac(fname,array,sachead,nerr) ! problem with nerr=-1
use sacio
parameter(nmax=4000000)
character*(*):: fname
type(sac_head):: sachead
real,dimension(nmax):: array
integer:: funit,iso,nerr,npts,i
nerr=0
funit=10
npts=sachead%npts
open(funit,file=fname,form='unformatted',action='write',recl=158*4+npts*4,access='direct',iostat=iso)
if(iso/=0)then
   write(*,*)"Error in openning file ",trim(fname)
   write(*,*)"Stop the subroutine write_sac"
   nerr=-1
   return
endif
write(funit,rec=1,iostat=iso)sachead,(array(i),i=1,npts)
if(iso/=0)then
   write(*,*)"Error in write file ",trim(fname)
   write(*,*)"Stop the subroutine write_sac"
   nerr=-1
   return
endif
close(funit)
return
end subroutine write_sac

subroutine initial_sachead(sachead)
use sacio
type(sac_head):: sachead
sachead=sac_head(-12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345., -12345., -12345., -12345., -12345.,&
                 -12345,  -12345,  -12345,  -12345,  -12345,&
                 -12345,  -12345,  -12345,  -12345,  -12345,&
                 -12345,  -12345,  -12345,  -12345,  -12345,&
                 -12345,  -12345,  -12345,  -12345,  -12345,&
                 -12345,  -12345,  -12345,  -12345,  -12345,&
                 -12345,  -12345,  -12345,  -12345,  -12345,&
                 -12345,  -12345,  -12345,  -12345,  -12345,&
                 -12345,  -12345,  -12345,  -12345,  -12345,&
                 "-12345  ","-12345          ",   &
                 "-12345  ","-12345  ","-12345  ",&
                 "-12345  ","-12345  ","-12345  ","-12345  ","-12345  ",&
                 "-12345  ","-12345  ","-12345  ","-12345  ","-12345  ",&
                 "-12345  ","-12345  ","-12345  ","-12345  ","-12345  ",&
                 "-12345  ","-12345  ","-12345  ")
return
end subroutine initial_sachead


subroutine write_ncf_sac(fname,array,sachead1,sachead2,npts,num,nerr)
use sacio
implicit none
character*(*) fname
integer:: nmax,npts,nerr,num,i,funit,iso
type(sac_head):: sachead1,sachead2,sachead
parameter(nmax=4000000)
real:: array(nmax)
funit=10
!call initial_sachead(sachead)
sachead=sachead1
sachead%npts =npts
sachead%delta=sachead1%delta
sachead%evla =sachead1%stla
sachead%evlo =sachead1%stlo
sachead%stla =sachead2%stla
sachead%stlo =sachead2%stlo
sachead%unused13=num
sachead%b=-int(npts/2)*sachead%delta
sachead%e= int(npts/2)*sachead%delta
sachead%kevnm=trim(sachead1%kstnm)//"_"//trim(sachead2%kstnm)
sachead%o=0
!write(*,*)sachead%delta,sachead%evla,sachead%evlo
!write(*,*)sachead%b,sachead%o,sachead%unused13
!do i=1,npts
!   write(*,*)array(i)
!enddo
call write_sac(fname,array,sachead,nerr)
!open(funit,file=fname,form='unformatted',action='write',recl=158*4+npts*4,access='direct',iostat=iso)
!if(iso/=0)then
!   write(*,*)"Error in openning file ",trim(fname)
!   write(*,*)"Stop the subroutine write_sac"
!   nerr=-1
!   return
!endif
!write(funit,rec=1,iostat=iso)sachead,(array(i),i=1,npts)
!if(iso/=0)then
!   write(*,*)"Error in write file ",trim(fname)
!   write(*,*)"Stop the subroutine write_sac"
!   nerr=-1
!   return
!endif
!close(funit)
!return
if(nerr.eq.-1)then
   write(*,*) 'Something wrong with ncf sac writing for file: ',trim(fname)
endif
end subroutine write_ncf_sac
