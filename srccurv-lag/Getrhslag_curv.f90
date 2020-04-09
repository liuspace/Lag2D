!
!...subroutine: Calculate the F^* N ds for interior faces using SMS for cubic cell
!...with different skills
!
subroutine getfnds_lagsmsif_hybrid_cubic_l(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac), intent(inout)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(inout)::gesgq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,isv,ishp
integer::iv
integer::nintf
integer::ipoin
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac)      :: ipf
real*8,dimension(1:ndimn,  1:2)::xpf
real*8,dimension(1:nvqua)::dsprq,dspsq,shpq
real*8,dimension(nvtri):: shpt
real*8, dimension(nvfac,2):: shpf,dsprf
!...local real array
real*8::posif(4)
!...local real number
real*8::xsc,ysc
real*8::anx, any,dxdr,dydr
real*8::dr, ds,r,s,djaco
real*8::c16,c10
!...Allocatble
integer,allocatable::ipift(:,:),ipifq(:,:)
real*8, allocatable::xpq(:,:),xpt(:,:),posiq(:,:),xpq_if(:,:),posqc(:,:),xpqv(:,:)
!
data c16   /0.1666666666666666d0 /
data c10   /1.d0 /

!
!...Part I: Preliminary setup
!
!...Specify nintf
if(ncurv.eq.1)then
nintf = 4

allocate(ipifq(nvfac,nintf))

allocate(xpt(1:ndimn, 1:nvtri+1))
allocate(xpq(1:ndimn, 1:(nvqua)))
elseif(ncurv.eq.2)then
nintf = 8

allocate(ipifq(2,nintf))

allocate(xpq(1:ndimn, 1:(nvqua)))
allocate(xpq_if(1:ndimn, 1:(nvqua+4)))
endif
!
ipifq(1, 1)= 5; ipifq(2, 1)=11;
ipifq(1, 2)= 9; ipifq(2, 2)= 7;
ipifq(1, 3)= 6; ipifq(2, 3)=12;
ipifq(1, 4)=10; ipifq(2, 4)= 8;
ipifq(1, 5)= 7; ipifq(2, 5)= 9;
ipifq(1, 6)=11; ipifq(2, 6)= 5;
ipifq(1, 7)= 8; ipifq(2, 7)=10;
ipifq(1, 8)=12; ipifq(2, 8)= 6;

!
!...Part II: Area vector for linear boundary face...
!
do ifa =1 , nbfac
!
iel = intfac(1, ifa)
!
ipf(1:2) = intfac(3:4, ifa)
!
xpf(1, 1:2) = coord(1, ipf(1:2))
xpf(2, 1:2) = coord(2, ipf(1:2))
!
anx =   xpf(2, 2) - xpf(2, 1)
any = -(xpf(1, 2) - xpf(1, 1))
!
gflag(1, ifa) = anx/sqrt(anx**2 + any**2)
gflag(2, ifa) = any/sqrt(anx**2 + any**2)
gflag(3, ifa) = sqrt(anx**2 + any**2)
!
enddo
!
!...Part III: Triangular edge
!
do ie=1, ntria
!
ipt(1:nvtri) = iptri(1:nvtri,ie)

!...coordinates
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
print*,'Triangle will be implmented in future in getfnds_lagsmsif_hybrid_cubic!'
stop

enddo !enddo for tria
!
!...Part IIII: Quad edge
!
do 55 ie=1, nquad

ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq = 0.d0
xpq_if = 0.d0

xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpq_if(1, 1:nvqua) = xpq(1, 1:nvqua)
xpq_if(2, 1:nvqua) = xpq(2, 1:nvqua)

!...The normal vector for interior faces
xpf = 0.d0

do ifa = 1,nintf
!...Interior face
xpf(1:2, 1) =  xpq_if(1:2, ipifq(1, ifa))
xpf(1:2, 2) =  xpq_if(1:2, ipifq(2, ifa))

!
r=-1.d0
call getshapfct_edge(0,2,shpf(1:2,1), dsprf(1:2,1), r)

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0

!...Get the determinat Jacobian with linear meshes
do ishp = 1, 2
dxdr = dxdr + dsprf(ishp, 1)*xpf(1, ishp)
dydr = dydr + dsprf(ishp, 1)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
gesgq(1, 8+ifa , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 8+ifa , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 8+ifa , ie) = djaco*2.d0*1.5d0!2.d0/(1.d0-sqrt(5.d0)/5.d0)

!
enddo !do ifa = 1,4

!
!if(ie.eq.19.or.ie.eq.380)then
!print*,'face',xpf(1:2, 1:4)
!print*,'ges',ie,gesgq(1:3,10 , ie),gesgq(1:3, 16 , ie)
!endif
!
55 enddo  !...(1)ie=1,nelem
!
end subroutine getfnds_lagsmsif_hybrid_cubic_l
!
!...subroutine: Calculate the F^* N ds for interior faces using SMS
!for cubic cell...
!
subroutine getfnds_lagsmsif_hybrid_cubic4(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag
!...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac), intent(inout)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(inout)::gesgq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,isv,ishp
integer::iv
integer::nintf
integer::ipoin
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq,idmpp
integer,dimension(1:nvfac)      :: ipf
integer,dimension(1:nvfac, 1:3) :: ipft
integer,dimension(1:nvfac, 1:4) :: ipfq
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8,dimension(1:nvqua)::dshprq,dshpsq,shpq
real*8,dimension(nvtri):: shpt
real*8, dimension(nvfac,2):: shpf,dshprf
!...local real array
real*8,dimension(2,2)::amtsp,rhssp
real*8::posif(4)
!...local real number
real*8::detmt,rf
real*8::anx, any,dxdr,dydr
real*8::dr, ds,r,s,djaco
real*8::c16,c10
!...Allocatble
integer,allocatable::ipift(:,:),ipifq(:,:)
real*8, allocatable::xpq(:,:),xpt(:,:),posiq(:,:),xpq_if(:,:)
!
data c16   /0.1666666666666666d0 /
data c10   /1.d0 /

!
!...Part I: Preliminary setup

!...Specify nintf
if(ncurv.eq.1)then
nintf = 4

allocate(ipifq(nvfac,nintf))

allocate(xpt(1:ndimn, 1:nvtri+1))
allocate(xpq(1:ndimn, 1:(nvqua)))
elseif(ncurv.eq.2)then
nintf = 8

allocate(ipifq(nvfac,nintf))

allocate(xpq(1:ndimn, 1:(nvqua)))
allocate(xpq_if(1:ndimn, 1:(nvqua+4)))
allocate(posiq(1:ndimn, 1:(nvqua+4)))
endif
!
ipifq(1, 1)= 5; ipifq(2, 1)=13; ipifq(3, 1)=13; ipifq(4, 1)=16;
ipifq(1, 2)= 9; ipifq(2, 2)=14; ipifq(3, 2)=14; ipifq(4, 2)=15;
ipifq(1, 3)= 6; ipifq(2, 3)=14; ipifq(3, 3)=14; ipifq(4, 3)=13;
ipifq(1, 4)=10; ipifq(2, 4)=15; ipifq(3, 4)=15; ipifq(4, 4)=16;
ipifq(1, 5)= 7; ipifq(2, 5)=15; ipifq(3, 5)=15; ipifq(4, 5)=14;
ipifq(1, 6)=11; ipifq(2, 6)=16; ipifq(3, 6)=16; ipifq(4, 6)=13;
ipifq(1, 7)= 8; ipifq(2, 7)=16; ipifq(3, 7)=16; ipifq(4, 7)=15;
ipifq(1, 8)=12; ipifq(2, 8)=13; ipifq(3, 8)=13; ipifq(4, 8)=14;
!
ipfq(1, 1)= 1; ipfq(2, 1)= 2; ipfq(3, 1)= 5; ipfq(4, 1)= 9;
ipfq(1, 2)= 2; ipfq(2, 2)= 3; ipfq(3, 2)= 6; ipfq(4, 2)=10;
ipfq(1, 3)= 3; ipfq(2, 3)= 4; ipfq(3, 3)= 7; ipfq(4, 3)=11;
ipfq(1, 4)= 4; ipfq(2, 4)= 1; ipfq(3, 4)= 8; ipfq(4, 4)=12;
!
posif(1)=-1.d0; posif(2)= 1.d0; posif(3)=-sqrt(5.d0)/5.d0;
posif(4)=sqrt(5.d0)/5.d0;
posiq(1, 1) = -1.d0; posiq(2, 1) = -1.d0;
posiq(1, 2) =  1.d0; posiq(2, 2) = -1.d0;
posiq(1, 3) =  1.d0; posiq(2, 3) =  1.d0;
posiq(1, 4) = -1.d0; posiq(2, 4) =  1.d0;
posiq(1, 5) = posif(3); posiq(2, 5) = -1.d0;
posiq(1, 6) =  1.d0;      posiq(2, 6) = posif(3);
posiq(1, 7) = posif(4); posiq(2, 7) =  1.d0;
posiq(1, 8) = -1.d0;      posiq(2, 8) = posif(4);
posiq(1, 9) = posif(4); posiq(2, 9) = -1.d0;
posiq(1,10) =  1.d0;      posiq(2,10) =  posif(4);
posiq(1,11) = posif(3); posiq(2,11) =  1.d0;
posiq(1,12) = -1.d0;      posiq(2,12) = posif(3);
posiq(1,13) = posif(3); posiq(2,13) = posif(3);
posiq(1,14) = posif(4); posiq(2,14) = posif(3);
posiq(1,15) = posif(4); posiq(2,15) = posif(4);
posiq(1,16) = posif(3); posiq(2,16) =  posif(4);
!
idmpp(1:4) =  (/1,2,3,4/)
idmpp(5:8) =  (/5,9,7,11/)
idmpp(9:12) = (/6,10,8,12/)

!
!...Part II: Area vector for linear boundary face...
!
do ifa =1 , nbfac
!
iel = intfac(1, ifa)
!
ipf(1:2) = intfac(3:4, ifa)
!
xpf(1, 1:2) = coord(1, ipf(1:2))
xpf(2, 1:2) = coord(2, ipf(1:2))
!
anx =   xpf(2, 2) - xpf(2, 1)
any = -(xpf(1, 2) - xpf(1, 1))
!
gflag(1, ifa) = anx/sqrt(anx**2 + any**2)
gflag(2, ifa) = any/sqrt(anx**2 + any**2)
gflag(3, ifa) = sqrt(anx**2 + any**2)
!
enddo
!
!...Part III: Triangular edge
!
do ie=1, ntria
!
ipt(1:nvtri) = iptri(1:nvtri,ie)

!...coordinates
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
print*,'Triangle will be implmented in future in getfnds_lagsmsif_hybrid_cubic!'
stop

enddo !enddo for tria
!
!...Part IIII: Quad edge
!
do 55 ie=1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

xpq = 0.d0
xpq_if = 0.d0

!...coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpq_if(1, 1:nvqua) = xpq(1, 1:nvqua)
xpq_if(2, 1:nvqua) = xpq(2, 1:nvqua)

!...Coordinates based on shape function
if(ncurv.eq.2)then

!...Get the physical coordinates for FE shape functions
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!...Get the interior points for internal face
!
xpf = 0.d0
!...Get the normal vector and Jacobian for curved nodes...
do ifa = 1,8

!...Edge 9 (p95)
xpf(1, 1:2) =  xpq_if(1, ipifq(1:2, ifa))
xpf(2, 1:2) =  xpq_if(2, ipifq(1:2, ifa))
!
r  = 0.d0
call getshapfct_edge(0,2,shpf(1:2,1), dshprf(1:2,1), r)
!
r = 0.d0
s = 0.d0
do ishp = 1, 2
r = r + shpf(ishp, 1)*posiq(1, ipifq(ishp, ifa))
s = s + shpf(ishp, 1)*posiq(2, ipifq(ishp, ifa))
enddo
!
call getshapfct_quad(ncurv,nvqua,shpq, dshprq, dshpsq, r, s)
!
xpf(1:2, 3) = 0.d0
do ishp = 1, nvqua
xpf(1:2, 3) = xpf(1:2, 3) + shpq(ishp)*xpq(1:2, ishp)
enddo


!if(ie.eq.19.and.ifa.eq.8)then
!print*,'gesfacew2',ie,ifa,r,s,shpq(:)*xpq(1,:)
!endif
!
r=-1.d0

call getshapfct_edge(1,3,shpf(:,1), dshprf(:,1), r)
!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dshprf(ishp, 1)*xpf(1, ishp)
dydr = dydr + dshprf(ishp, 1)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
gesgq(1, 8+ifa , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 8+ifa , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 8+ifa , ie) = djaco*2.d0*1.d0
!if(ie.eq.1)then
!print*,'face',xpq_if(1:2, ipifq(1:4, 1))
!print*,'gesfacew',ie,ifa,gesgq(1:3, 8+ifa , ie)
!endif
!
enddo !do ifa = 1,4
!
endif
!
!if(ie.eq.19.or.ie.eq.380)then
!print*,'face',xpf(1:2, 1:4)
!print*,'ges',ie,gesgq(1:3,10 , ie),gesgq(1:3, 16 , ie)
!endif
!
55 enddo  !...(1)ie=1,nelem
!
end subroutine getfnds_lagsmsif_hybrid_cubic4
!
!...subroutine: Calculate the F^* N ds for interior faces using SMS for cubic cell
!...with different skills
!
subroutine getfnds_lagsmsif_hybrid_cubic2(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac), intent(inout)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(inout)::gesgq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,isv,ishp
integer::iv
integer::nintf
integer::ipoin
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac)      :: ipf
integer,dimension(1:nvfac, 1:3) :: ipft
integer,dimension(1:nvfac, 1:4) :: ipfq
real*8,dimension(1:ndimn,  1:3)::xpf
real*8,dimension(1:nvqua)::dshprq,dshpsq,shpq
real*8,dimension(nvtri):: shpt
real*8, dimension(nvfac,2):: shpf,dshprf
!...local real array
real*8::posif(4)
!...local real number
real*8::xsc,ysc
real*8::anx, any,dxdr,dydr
real*8::dr, ds,r,s,djaco
real*8::c16,c10
!...Allocatble
integer,allocatable::ipift(:,:),ipifq(:,:)
real*8, allocatable::xpq(:,:),xpt(:,:),posiq(:,:),xpq_if(:,:),posqc(:,:),xpqv(:,:)
!
data c16   /0.1666666666666666d0 /
data c10   /1.d0 /

!
!...Part I: Preliminary setup
!
!...Specify nintf
if(ncurv.eq.1)then
nintf = 4

allocate(ipifq(nvfac,nintf))

allocate(xpt(1:ndimn, 1:nvtri+1))
allocate(xpq(1:ndimn, 1:(nvqua)))
elseif(ncurv.eq.2)then
nintf = 8

allocate(ipifq(3,nintf))

allocate(posqc(2, 8),xpqv(2, 8))
allocate(xpq(1:ndimn, 1:(nvqua)))
allocate(xpq_if(1:ndimn, 1:(nvqua+4)))
allocate(posiq(1:ndimn, 1:(nvqua+4)))
endif
!
ipifq(1, 1)= 5; ipifq(2, 1)=13; ipifq(3, 1)=1;
ipifq(1, 2)= 9; ipifq(2, 2)=14; ipifq(3, 2)=2;
ipifq(1, 3)= 6; ipifq(2, 3)=14; ipifq(3, 3)=3;
ipifq(1, 4)=10; ipifq(2, 4)=15; ipifq(3, 4)=4;
ipifq(1, 5)= 7; ipifq(2, 5)=15; ipifq(3, 5)=5;
ipifq(1, 6)=11; ipifq(2, 6)=16; ipifq(3, 6)=6;
ipifq(1, 7)= 8; ipifq(2, 7)=16; ipifq(3, 7)=7;
ipifq(1, 8)=12; ipifq(2, 8)=13; ipifq(3, 8)=8;
!
ipifq(1, 1)= 5; ipifq(2, 1)=11; ipifq(3, 1)=1;
ipifq(1, 2)= 9; ipifq(2, 2)= 7; ipifq(3, 2)=2;
ipifq(1, 3)= 6; ipifq(2, 3)=12; ipifq(3, 3)=3;
ipifq(1, 4)=10; ipifq(2, 4)= 8; ipifq(3, 4)=4;
ipifq(1, 5)= 7; ipifq(2, 5)= 9; ipifq(3, 5)=5;
ipifq(1, 6)=11; ipifq(2, 6)= 5; ipifq(3, 6)=6;
ipifq(1, 7)= 8; ipifq(2, 7)=10; ipifq(3, 7)=7;
ipifq(1, 8)=12; ipifq(2, 8)= 6; ipifq(3, 8)=8;
!
ipfq(1, 1)= 1; ipfq(2, 1)= 2; ipfq(3, 1)= 5; ipfq(4, 1)= 9;
ipfq(1, 2)= 2; ipfq(2, 2)= 3; ipfq(3, 2)= 6; ipfq(4, 2)=10;
ipfq(1, 3)= 3; ipfq(2, 3)= 4; ipfq(3, 3)= 7; ipfq(4, 3)=11;
ipfq(1, 4)= 4; ipfq(2, 4)= 1; ipfq(3, 4)= 8; ipfq(4, 4)=12;
!
posif(1)=-1.d0; posif(2)= 1.d0; posif(3)=-sqrt(5.d0)/5.d0; posif(4)=sqrt(5.d0)/5.d0;
!posif(3)=-1.d0/3.d0; posif(4)=1.d0/3.d0;

posiq(1, 1) = -1.d0; posiq(2, 1) = -1.d0;
posiq(1, 2) =  1.d0; posiq(2, 2) = -1.d0;
posiq(1, 3) =  1.d0; posiq(2, 3) =  1.d0;
posiq(1, 4) = -1.d0; posiq(2, 4) =  1.d0;
posiq(1, 5) = posif(3);   posiq(2, 5) = -1.d0;
posiq(1, 6) =  1.d0;      posiq(2, 6) = posif(3);
posiq(1, 7) = posif(4);   posiq(2, 7) =  1.d0;
posiq(1, 8) = -1.d0;      posiq(2, 8) = posif(4);
posiq(1, 9) = posif(4);   posiq(2, 9) = -1.d0;
posiq(1,10) =  1.d0;      posiq(2,10) =  posif(4);
posiq(1,11) = posif(3);   posiq(2,11) =  1.d0;
posiq(1,12) = -1.d0;      posiq(2,12) = posif(3);
posiq(1,13) = posif(3); posiq(2,13) = posif(3);
posiq(1,14) = posif(4); posiq(2,14) = posif(3);
posiq(1,15) = posif(4); posiq(2,15) = posif(4);
posiq(1,16) = posif(3); posiq(2,16) =  posif(4);

!...Additional points defining the interior faces
posqc(1, 1) =-sqrt(5.d0)/5.d0;             posqc(2, 1) = -(sqrt(5.d0)/5.d0+1.d0)/2.d0;
posqc(1, 2) = sqrt(5.d0)/5.d0;             posqc(2, 2) = -(sqrt(5.d0)/5.d0+1.d0)/2.d0;
posqc(1, 3) = (sqrt(5.d0)/5.d0+1.d0)/2.d0; posqc(2, 3) = -sqrt(5.d0)/5.d0;
posqc(1, 4) = (sqrt(5.d0)/5.d0+1.d0)/2.d0; posqc(2, 4) =  sqrt(5.d0)/5.d0;
posqc(1, 5) = sqrt(5.d0)/5.d0;             posqc(2, 5) =  (sqrt(5.d0)/5.d0+1.d0)/2.d0;
posqc(1, 6) =-sqrt(5.d0)/5.d0;             posqc(2, 6) =  (sqrt(5.d0)/5.d0+1.d0)/2.d0;
posqc(1, 7) =-(sqrt(5.d0)/5.d0+1.d0)/2.d0; posqc(2, 7) =  sqrt(5.d0)/5.d0;
posqc(1, 8) =-(sqrt(5.d0)/5.d0+1.d0)/2.d0; posqc(2, 8) = -sqrt(5.d0)/5.d0;

!
!...Part II: Area vector for linear boundary face...
!
do ifa =1 , nbfac
!
iel = intfac(1, ifa)
!
ipf(1:2) = intfac(3:4, ifa)
!
xpf(1, 1:2) = coord(1, ipf(1:2))
xpf(2, 1:2) = coord(2, ipf(1:2))
!
anx =   xpf(2, 2) - xpf(2, 1)
any = -(xpf(1, 2) - xpf(1, 1))
!
gflag(1, ifa) = anx/sqrt(anx**2 + any**2)
gflag(2, ifa) = any/sqrt(anx**2 + any**2)
gflag(3, ifa) = sqrt(anx**2 + any**2)
!
enddo
!
!...Part III: Triangular edge
!
do ie=1, ntria
!
ipt(1:nvtri) = iptri(1:nvtri,ie)

!...coordinates
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
print*,'Triangle will be implmented in future in getfnds_lagsmsif_hybrid_cubic!'
stop

enddo !enddo for tria
!
!...Part IIII: Quad edge
!
do 55 ie=1, nquad

ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq = 0.d0
xpq_if = 0.d0

xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpq_if(1, 1:nvqua) = xpq(1, 1:nvqua)
xpq_if(2, 1:nvqua) = xpq(2, 1:nvqua)

!...Get the physical coordinates for FE shape functions
if(ncurv.eq.2)then
call getcoord_fe(ncurv, nvfac, nvqua, xpq)
endif

!...The interior nodes defining Lagrangian cell
do isv = 13, 16
r = posiq(1,isv)
s = posiq(2,isv)

!...  shape function
call getshapfct_quad(ncurv,nvqua,shpq, dshprq, dshpsq, r, s)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, nvqua
xsc = xsc + shpq(ishp)*xpq(1,ishp)
ysc = ysc + shpq(ishp)*xpq(2,ishp)
enddo
!
xpq_if(1, isv) = xsc
xpq_if(2, isv) = ysc
!
enddo

!...The nodes defining interior curved faces
do isv = 1, 8
r = posqc(1,isv)
s = posqc(2,isv)

!...  shape function
call getshapfct_quad(ncurv,nvqua,shpq, dshprq, dshpsq, r, s)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, nvqua
xsc = xsc + shpq(ishp)*xpq(1,ishp)
ysc = ysc + shpq(ishp)*xpq(2,ishp)
enddo
!
xpqv(1, isv) = xsc
xpqv(2, isv) = ysc
!
enddo

!...The normal vector for interior faces
xpf = 0.d0
do ifa = 1,8

!...Interior face
xpf(1:2, 1) =  xpq_if(1:2, ipifq(1, ifa))
xpf(1:2, 2) =  xpq_if(1:2, ipifq(2, ifa))
xpf(1:2, 3) =  xpqv(1:2, ipifq(3, ifa))
!
r=-1.d0

call getshapfct_edge(0,2,shpf(1:2,1), dshprf(1:2,1), r)
!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 2
dxdr = dxdr + dshprf(ishp, 1)*xpf(1, ishp)
dydr = dydr + dshprf(ishp, 1)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
gesgq(1, 8+ifa , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 8+ifa , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 8+ifa , ie) = djaco*2.d0*1.d0!2.d0/(1.d0-sqrt(5.d0)/5.d0)

!if(ie.eq.1)then
!print*,'gesfacew',ie,ifa,gesgq(1:3, 8+ifa , ie)
!endif
!
enddo !do ifa = 1,4

!
!if(ie.eq.19.or.ie.eq.380)then
!print*,'face',xpf(1:2, 1:4)
!print*,'ges',ie,gesgq(1:3,10 , ie),gesgq(1:3, 16 , ie)
!endif
!
55 enddo  !...(1)ie=1,nelem
!
end subroutine getfnds_lagsmsif_hybrid_cubic2

!
!...subroutine: Calculate the F^* N ds for interior faces using SMS for cubic cell...
!
subroutine getfnds_lagsmsif_hybrid_cubic(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac), intent(inout)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(inout)::gesgq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,isv,ishp
integer::iv
integer::nintf
integer::ipoin
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq,idmpp
integer,dimension(1:nvfac)      :: ipf
integer,dimension(1:nvfac, 1:3) :: ipft
integer,dimension(1:nvfac, 1:4) :: ipfq
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8,dimension(1:nvqua)::dshprq,dshpsq,shpq
real*8,dimension(nvtri):: shpt
real*8, dimension(nvfac,2):: shpf,dshprf
!...local real array
real*8,dimension(2,2)::amtsp,rhssp
real*8::posif(4)
!...local real number
real*8::detmt,rf
real*8::anx, any,dxdr,dydr
real*8::dr, ds,r,s,djaco
real*8::c16,c10
!...Allocatble
integer,allocatable::ipift(:,:),ipifq(:,:)
real*8, allocatable::xpq(:,:),xpt(:,:),posiq(:,:),xpq_if(:,:)
!
data c16   /0.1666666666666666d0 /
data c10   /1.d0 /

!
!...Part I: Preliminary setup

!...Specify nintf
if(ncurv.eq.1)then
nintf = 4

allocate(ipifq(nvfac,nintf))

allocate(xpt(1:ndimn, 1:nvtri+1))
allocate(xpq(1:ndimn, 1:(nvqua)))
elseif(ncurv.eq.2)then
nintf = 8

allocate(ipifq(nvfac,nintf))

allocate(xpq(1:ndimn, 1:(nvqua)))
allocate(xpq_if(1:ndimn, 1:(nvqua+4)))
allocate(posiq(1:ndimn, 1:(nvqua+4)))
endif
!
ipifq(1, 1)= 5; ipifq(2, 1)=11; ipifq(3, 1)=13; ipifq(4, 1)=16;
ipifq(1, 2)= 9; ipifq(2, 2)= 7; ipifq(3, 2)=14; ipifq(4, 2)=15;
ipifq(1, 3)= 6; ipifq(2, 3)=12; ipifq(3, 3)=14; ipifq(4, 3)=13;
ipifq(1, 4)=10; ipifq(2, 4)= 8; ipifq(3, 4)=15; ipifq(4, 4)=16;
ipifq(1, 5)= 7; ipifq(2, 5)= 9; ipifq(3, 5)=15; ipifq(4, 5)=14;
ipifq(1, 6)=11; ipifq(2, 6)= 5; ipifq(3, 6)=16; ipifq(4, 6)=13;
ipifq(1, 7)= 8; ipifq(2, 7)=10; ipifq(3, 7)=16; ipifq(4, 7)=15;
ipifq(1, 8)=12; ipifq(2, 8)= 6; ipifq(3, 8)=13; ipifq(4, 8)=14;
!
ipfq(1, 1)= 1; ipfq(2, 1)= 2; ipfq(3, 1)= 5; ipfq(4, 1)= 9;
ipfq(1, 2)= 2; ipfq(2, 2)= 3; ipfq(3, 2)= 6; ipfq(4, 2)=10;
ipfq(1, 3)= 3; ipfq(2, 3)= 4; ipfq(3, 3)= 7; ipfq(4, 3)=11;
ipfq(1, 4)= 4; ipfq(2, 4)= 1; ipfq(3, 4)= 8; ipfq(4, 4)=12;
!
posif(1)=-1.d0; posif(2)= 1.d0; posif(3)=-sqrt(5.d0)/5.d0; posif(4)=sqrt(5.d0)/5.d0;
!posif(3)=-1.d0/3.d0; posif(4)=1.d0/3.d0;

posiq(1, 1) = -1.d0; posiq(2, 1) = -1.d0;
posiq(1, 2) =  1.d0; posiq(2, 2) = -1.d0;
posiq(1, 3) =  1.d0; posiq(2, 3) =  1.d0;
posiq(1, 4) = -1.d0; posiq(2, 4) =  1.d0;
posiq(1, 5) = posif(3); posiq(2, 5) = -1.d0;
posiq(1, 6) =  1.d0;      posiq(2, 6) = posif(3);
posiq(1, 7) = posif(4); posiq(2, 7) =  1.d0;
posiq(1, 8) = -1.d0;      posiq(2, 8) = posif(4);
posiq(1, 9) = posif(4); posiq(2, 9) = -1.d0;
posiq(1,10) =  1.d0;      posiq(2,10) =  posif(4);
posiq(1,11) = posif(3); posiq(2,11) =  1.d0;
posiq(1,12) = -1.d0;      posiq(2,12) = posif(3);
posiq(1,13) = posif(3); posiq(2,13) = posif(3);
posiq(1,14) = posif(4); posiq(2,14) = posif(3);
posiq(1,15) = posif(4); posiq(2,15) = posif(4);
posiq(1,16) = posif(3); posiq(2,16) =  posif(4);
!
idmpp(1:4) =  (/1,2,3,4/)
idmpp(5:8) =  (/5,9,7,11/)
idmpp(9:12) = (/6,10,8,12/)

!
!...Part II: Area vector for linear boundary face...
!
do ifa =1 , nbfac
!
iel = intfac(1, ifa)
!
ipf(1:2) = intfac(3:4, ifa)
!
xpf(1, 1:2) = coord(1, ipf(1:2))
xpf(2, 1:2) = coord(2, ipf(1:2))
!
anx =   xpf(2, 2) - xpf(2, 1)
any = -(xpf(1, 2) - xpf(1, 1))
!
gflag(1, ifa) = anx/sqrt(anx**2 + any**2)
gflag(2, ifa) = any/sqrt(anx**2 + any**2)
gflag(3, ifa) = sqrt(anx**2 + any**2)
!
enddo
!
!...Part III: Triangular edge
!
do ie=1, ntria
!
ipt(1:nvtri) = iptri(1:nvtri,ie)

!...coordinates
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
print*,'Triangle will be implmented in future in getfnds_lagsmsif_hybrid_cubic!'
stop

enddo !enddo for tria
!
!...Part IIII: Quad edge
!
do 55 ie=1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

 xpq = 0.d0
 xpq_if = 0.d0

!...coordinates
 xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
 xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
 xpq_if(1, 1:nvqua) = xpq(1, 1:nvqua)
 xpq_if(2, 1:nvqua) = xpq(2, 1:nvqua)

!...Coordinates based on shape function
if(ncurv.eq.2)then

!...Get the physical coordinates for FE shape functions
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!...Get the interior points for internal face
!
xpf = 0.d0
!...Get the normal vector and Jacobian for curved nodes...
do ifa = 1,8

!...Edge 9 (p95)
xpf(1, 1:2) =  xpq_if(1, ipifq(1:2, ifa))
xpf(2, 1:2) =  xpq_if(2, ipifq(1:2, ifa))
!
r  = -1.d0/3.d0
call getshapfct_edge(0,2,shpf(1:2,1), dshprf(1:2,1), r)
!
r = 0.d0
s = 0.d0
do ishp = 1, 2
 r = r + shpf(ishp, 1)*posiq(1, ipifq(ishp, ifa))
 s = s + shpf(ishp, 1)*posiq(2, ipifq(ishp, ifa))
enddo
!
call getshapfct_quad(ncurv,nvqua,shpq, dshprq, dshpsq, r, s)
!
xpf(1:2, 3) = 0.d0
do ishp = 1, nvqua
 xpf(1:2, 3) = xpf(1:2, 3) + shpq(ishp)*xpq(1:2, ishp)
enddo

r  = 1.d0/3.d0
call getshapfct_edge(0,2,shpf(1:2,1), dshprf(1:2,1), r)
!
r = 0.d0
s = 0.d0
do ishp = 1, 2
 r = r + shpf(ishp, 1)*posiq(1, ipifq(ishp, ifa))
 s = s + shpf(ishp, 1)*posiq(2, ipifq(ishp, ifa))
enddo
!
!if(ie.eq.19.and.ifa.eq.8)then
!print*,'gesfacew',ie,ifa,r,s,xpf(1:2, 1:2)
!endif
!
call getshapfct_quad(ncurv,nvqua,shpq, dshprq, dshpsq, r, s)
!
xpf(1:2, 4) = 0.d0
do ishp = 1, nvqua
 xpf(1:2, 4) = xpf(1:2, 4) + shpq(ishp)*xpq(1:2, ishp)
enddo

!if(ie.eq.19.and.ifa.eq.8)then
!print*,'gesfacew2',ie,ifa,r,s,shpq(:)*xpq(1,:)
!endif
!
r=-1.d0

call getshapfct_edge(ncurv,nvfac,shpf(:,1), dshprf(:,1), r)
!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 4
dxdr = dxdr + dshprf(ishp, 1)*xpf(1, ishp)
dydr = dydr + dshprf(ishp, 1)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
gesgq(1, 8+ifa , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 8+ifa , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 8+ifa , ie) = djaco*2.d0

!if(ie.eq.380.and.ifa.eq.2)then
!print*,'face',xpq_if(1:2, ipifq(1:4, 1))
!print*,'gesfacew',ie,ifa,xpf(1:2, 1:4)
!endif
!
enddo !do ifa = 1,4
!
endif
!
!if(ie.eq.19.or.ie.eq.380)then
!print*,'face',xpf(1:2, 1:4)
!print*,'ges',ie,gesgq(1:3,10 , ie),gesgq(1:3, 16 , ie)
!endif
!
55 enddo  !...(1)ie=1,nelem
!
end subroutine getfnds_lagsmsif_hybrid_cubic
!
!...subroutine: Calculate the F^* N dsfor all faces for hybrid grids...
!
subroutine getfnds_lagsmsef_cubic(gflag,gesgq,gesgt,intfac,inpoel,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac),intent(inout)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(out)::gesgq
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv,ishp
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:2, 1:2)    ::comatr !...cofactor matrix...
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8,dimension(1:ndimn, 1:nvtri)::xpt
real*8,dimension(1:ndimn, 1:nvqua)::xpq
real*8::vnorm(1:2)
real*8,dimension(1:nvtri):: xv, yv
real*8::posit(2, 9)
real*8::posiq(2, 16)
real*8::posif(4)
real*8::shp(nptri),dspr(nptri),dsps(nptri)
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::dxdr,dxds,dydr,dyds
!...local real number
real*8::dwav1,dwav2,larea,farea
real*8::anx, any
real*8::r,s
real*8::dr, ds, rc, sc
real*8::c16, c10
!
data c16   /0.1666666666666666d0 /
data c10   /1.0d0 /
!
!...Specify the posi for 9 simpson points...
!
!...Triangle

!
!...8-nodes quad
!
posif(1)=-1.d0; posif(2)= 1.d0; posif(3)=-sqrt(5.d0)/5.d0; posif(4)=sqrt(5.d0)/5.d0;

!...Uniform space
!posif(3)=-1.d0/3.d0; posif(4)=1.d0/3.d0;
!
posiq(1, 1) =-1.d0;   posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0;   posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0;   posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0;   posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0;   posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0;   posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0;   posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0;   posiq(2, 8) =-1.d0;
posiq(1, 9) =posif(3); posiq(2, 9) =-1.d0;
posiq(1,10) =posif(4); posiq(2,10) =-1.d0;
posiq(1,11) = 1.d0;   posiq(2,11) =posif(3);
posiq(1,12) = 1.d0;   posiq(2,12) =posif(4);
posiq(1,13) =posif(4); posiq(2, 13) = 1.d0;
posiq(1,14) =posif(3); posiq(2, 14) = 1.d0;
posiq(1,15) =-1.d0;   posiq(2,15) = posif(4);
posiq(1,16) =-1.d0;   posiq(2,16) = posif(3);
!
do 100 ie=1, ntria!...(1)ifa=1,nafac
!...First step: calcualte the F^* NdS for every face of all the cells...
ipt(1:nvtri) = inpoel(1:nvtri,ie)
print*,'Tria will be implemented in future for getfnds_lagsmsef_cubic!'
!
100 enddo  !...(1)ie=1, ntria
!
!...Part III: Quad
!
do 200 ie=1, nquad !...(1)ie=1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...coordinates
xpq(1, 1:nvqua) = coord(1, ipqua(1:nvqua,ie))
xpq(2, 1:nvqua) = coord(2, ipqua(1:nvqua,ie))

!...Get coordinates for FEM
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!
do ig = 1, 16!...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)

!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!...Jacobian determinate
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo

!...Cofactor matrix for left cell
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa

!...Identify the local No. of one internal face for left cell...
if(ig.eq.1.or.ig.eq.2.or.ig.eq.9.or.ig.eq.10)then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 2.d0;
!
idfal = ig;
!
elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.11.or.ig.eq.12)then
!
vnorm(1) = 1.d0; vnorm(2) = 0.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.13.or.ig.eq.14)then
!
vnorm(1) = 0.d0;           vnorm(2) = 1.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.7.or.ig.eq.8.or.ig.eq.15.or.ig.eq.16)then
!
vnorm(1) =-1.d0;           vnorm(2) = 0.d0;
larea    = 2.d0
!
idfal = ig;
endif
!
anx = comatr(1, 1)*vnorm(1) + comatr(1, 2)*vnorm(2)
any = comatr(2, 1)*vnorm(1) + comatr(2, 2)*vnorm(2)
!
vnorm(1) = anx*larea
vnorm(2) = any*larea
!
!...Unit vector...
farea    = sqrt(vnorm(1)**2 + vnorm(2)**2) !...farea: face area...
!
gelagq(1, idfal, ie) = anx*larea/farea
gelagq(2, idfal, ie) = any*larea/farea
gelagq(3, idfal, ie) = farea
!
!if(ie.eq.1)print*,'gelag',xpq(1,:)
!
enddo !...ig from 1...9

!
gesgq(1, 1:8, ie) = gelagq(1, 1:8, ie)
gesgq(2, 1:8, ie) = gelagq(2, 1:8, ie)
gesgq(3, 1:8, ie) = gelagq(3, 1:8, ie)
!
gesgq(1, 17:24, ie) = gelagq(1, 9:16, ie)
gesgq(2, 17:24, ie) = gelagq(2, 9:16, ie)
gesgq(3, 17:24, ie) = gelagq(3, 9:16, ie)
!
200 enddo  !...(1)ie=1,nquad
!
!do ig = 1,16
!print*,'Inside getfnds_lag',ig,gesgq(3,1:24,1)
!enddo
!
end subroutine getfnds_lagsmsef_cubic
!
!...subroutine: Calculate the Riemann solutions for curved meshes with SMS...
!
subroutine getRiemann_lag_sms(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, fstart, fstarq, aflim, afvec, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac), intent(in)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:3,1:nsubt, 1:ntria),  intent(out)::fstart !...Riemann forces
real*8,dimension(1:ndimn,1:2,1:4,1:nsubq, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),             intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec

integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop, isg, ivsg
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(4, nsubq)::ipqsg
integer,dimension(3, nsubt)::iptsg
integer::indnd(npoin)

!...local real array
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8::munaci(2, 2)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: munacn(:,:,:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: munacu(:,:), snsigm(:,:)
real*8,allocatable:: munaclq(:,:,:,:,:,:),munaclt(:,:,:,:,:,:)
real*8,allocatable:: snsigmlq(:,:,:,:,:), munaulq(:,:,:,:,:)
real*8,allocatable:: snsigmlt(:,:,:,:,:), munault(:,:,:,:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:2, 1:2, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munaclq(1:2, 1:2, 1:2, 1:4, 1:nsubq, 1:nquad), munaulq(1:ndimn, 1:2, 1:4, 1:nsubq,  1:nquad),&
snsigmlq(1:ndimn, 1:2,  1:4, 1:nsubq,  1:nquad))
allocate (munaclt(1:2, 1:2, 1:2, 1:3, 1:nsubt, 1:ntria), munault(1:ndimn, 1:2, 1:3, 1:nsubt,  1:ntria),&
snsigmlt(1:ndimn, 1:2,  1:3, 1:nsubt,  1:ntria))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
!
!...Part I: Specify some gauss points...
!
if(ncurv.eq.1)then
!...Local vertex No. of gauss points in one unit
!...0 means the subcell vertex is void and inside the triangle
iptsg(1, 1) = 1; iptsg(2, 1) = 4; iptsg(3, 1) = 6;
iptsg(1, 2) = 4; iptsg(2, 2) = 5; iptsg(3, 2) = 6;
iptsg(1, 3) = 4; iptsg(2, 3) = 2; iptsg(3, 3) = 5;
iptsg(1, 4) = 6; iptsg(2, 4) = 5; iptsg(3, 4) = 3;

!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
elseif(ncurv.eq.2)then
!...Local vertex No. of gauss points in one unit
!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8

endif
!...Zero out vlave
vlave = 0.d0
indnd = 0

!...Mark the boundary nodes...
if(ncase.eq.2)then
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
indnd(ipf(1:nvfac)) = 1
enddo
endif
!
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
indnd(ipf(1:nvfac)) = 1
endif
enddo

!...Get averaged velocity at the node
do 950 ifa = 1 , nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
if(bface(3,ifa).eq.22)then
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(4,ifa).eq.221)then
vlave(2,ipf(1:nvfac)) = 0.d0
elseif(bface(4,ifa).eq.222)then
vlave(1,ipf(1:nvfac)) = 0.d0
endif
elseif(bface(3,ifa).eq.25)then
vlave(2,ipf(1:nvfac)) = 0.d0
vlave(1,ipf(1:nvfac)) = 0.d0
endif
950 enddo
!
!...Part II: Loop to get the info from Riemann solver
!
do iloop= 1, 1

!...Give vlave
vlave= ustar

!...Zero out munacn
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0

!...Tria will be implemented in future!

!...Quad
if(nquad.gt.0) call getriem_quad_sms2(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)

!...Boundary condition
!call getbc_lagc_general_curv(bface, intfac, fpres, coord, munacn, munacu, snsigm, itime)

!...Periodic boundary condition for 1D isentropic Sin problem...
!if(ncase.eq.12) call getbc_prdic(bface, munacn, munacu, snsigm)


!...Update the velocity at the vertex...
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
detma = munacn(1, 1, ipoin)*munacn(2, 2, ipoin) - munacn(2, 1, ipoin)*munacn(1, 2, ipoin)
munaci(1, 1) = munacn(2, 2, ipoin)/detma
munaci(1, 2) =-munacn(1, 2, ipoin)/detma
munaci(2, 1) =-munacn(2, 1, ipoin)/detma
munaci(2, 2) = munacn(1, 1, ipoin)/detma
!
rhsu1 = munacu(1, ipoin) - snsigm(1, ipoin) !- fpres(1, ipoin)
rhsu2 = munacu(2, ipoin) - snsigm(2, ipoin) !- fpres(2, ipoin)
!
ustar(1, ipoin) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ustar(2, ipoin) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2

!if(ipoin.eq.33)print*,itime,ipoin,ustar(1:2,ipoin),detma,munacn(1,1,ipoin),snsigm(1:2, ipoin),&
!                 munacu(1:2, ipoin)
endif
enddo
!
!print*,'node',ustar(1,301)
!...Get the vertex velocity at the boundary
do 900 ifa = 1 , nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
if(bface(3,ifa).eq.22)then
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(4,ifa).eq.221)then
ustar(2,ipf(1:nvfac)) = 0.d0
elseif(bface(4,ifa).eq.222)then
ustar(1,ipf(1:nvfac)) = 0.d0
endif
elseif(bface(3,ifa).eq.25)then
ustar(2,ipf(1:nvfac)) = 0.d0
ustar(1,ipf(1:nvfac)) = 0.d0

!...Slatzman
if(ncase.eq.13)then
ustar(2,ipf(1:nvfac)) = 0.d0
ustar(1,ipf(1:nvfac)) = 1.d0
elseif(ncase.eq.15)then
ustar(2,ipf(1:nvfac)) = 0.d0
ustar(1,ipf(1:nvfac)) = 2.629369d0
endif
endif

!...Specify the exact velocity for Taylor-Green vortex
if(ncase.eq.-1)then
ustar(1, ipf(1:nvfac)) = sin(pi*coord(1,ipf(1:nvfac)))*cos(pi*coord(2,ipf(1:nvfac)))
ustar(2, ipf(1:nvfac)) =-cos(pi*coord(1,ipf(1:nvfac)))*sin(pi*coord(2,ipf(1:nvfac)))
endif
!
900 enddo

enddo !iloop
!
!...Part III: Get the Riemann forces for face integral
!
!...Tria
do ie = 1, ntria
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
do isg = 1, nsubt
do ivsg = 1, 3
!
iv =  iptsg(ivsg, isg)
!
do ifa =1, 2
fstart(1, ifa, ivsg, isg, ie) = snsigmlt(1, ifa, ivsg, isg, ie) + &
munaclt(1, 1, ifa, ivsg, isg,ie)*ustar(1, ipt(iv))+&
munaclt(2, 1, ifa, ivsg, isg,ie)*ustar(2, ipt(iv)) - munault(1, ifa, ivsg, isg, ie)
fstart(2, ifa, ivsg, isg,ie) = snsigmlt(2, ifa, ivsg, isg, ie) + &
munaclt(1, 2, ifa, ivsg, isg, ie)*ustar(1, ipt(iv))+&
munaclt(2, 2, ifa, ivsg, isg, ie)*ustar(2, ipt(iv))- munault(2, ifa, ivsg, isg, ie)
!
!if(ie.eq.1.and.isg.eq.1.and.ivsg.eq.1.and.ifa.eq.2)then
!print*,'fstar',snsigmlq(2, ifa, ivsg, isg, ie),munaclq(1, 2, ifa, ivsg, isg, ie),ustar(1, ipq(iv)),&
!munaclq(2, 2, ifa, ivsg, isg, ie),ustar(2, ipq(iv)), munaulq(2, ifa, ivsg, isg, ie)
!endif
!
enddo
enddo
enddo
!
enddo


!...Quad
do ie = 1, nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
do isg = 1, nsubq
do ivsg = 1, 4
!
iv =  ipqsg(ivsg, isg)

if(ipqsg(ivsg, isg).gt.12) cycle
!
do ifa =1, 2
fstarq(1, ifa, ivsg, isg, ie) = snsigmlq(1, ifa, ivsg, isg, ie) + &
munaclq(1, 1, ifa, ivsg, isg,ie)*ustar(1, ipq(iv))+&
munaclq(2, 1, ifa, ivsg, isg,ie)*ustar(2, ipq(iv)) - munaulq(1, ifa, ivsg, isg, ie)
fstarq(2, ifa, ivsg, isg,ie) = snsigmlq(2, ifa, ivsg, isg, ie) + &
munaclq(1, 2, ifa, ivsg, isg, ie)*ustar(1, ipq(iv))+&
munaclq(2, 2, ifa, ivsg, isg, ie)*ustar(2, ipq(iv))- munaulq(2, ifa, ivsg, isg, ie)
!
!if(ie.eq.2100.and.isg.eq.1.and.ivsg.eq.1)then
!print*,'fstar',ifa,fstarq(1, ifa, ivsg, isg, ie),snsigmlq(1, ifa, ivsg, isg, ie),munaclq(1, 1, ifa, ivsg, isg, ie),&
!ustar(1, ipq(iv)),&
!munaclq(2, 1, ifa, ivsg, isg, ie),ustar(2, ipq(iv)), munaulq(1, ifa, ivsg, isg, ie)
!endif
!
enddo
enddo
enddo
!
enddo
!
!
deallocate (munacn, bpres, fpres)
deallocate (munacu, snsigm, bnorm)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getRiemann_lag_sms
!
!...subroutine: Riemann input for hybrid curved quads using SMS....
!
subroutine getriem_quad_sms(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:ngesgq,1:nquad),       intent(in)::gesgq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),             intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8, dimension(1:2, 1:2, 1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
real*8, dimension(1:2, 1:2, 1:2, 1:4, 1:nsubq, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:nsubq, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:nsubq, 1:nquad), intent(out)::snsigmlq

!...Local integer
integer::ie, ideg, ielem, ifa, iv, isg, ivsg, ifsg

!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(8, nsubq)::fnqsg
integer,dimension(4, nsubq)::ipqsg
!...local real array
real*8,dimension(1:ndegr, 1:nvqua)::bq
real*8,dimension(1:ndegr, 1:4)::bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nq,1:4)::unsgq
real*8::aujmp(1:3, 1:4)
real*8::vnorm(1:3, 1:2, 1:6)
real*8::sigma(1:2, 1:2, 1:4)
real*8,dimension(1:2, 1:4)::murie
real*8,dimension(1:2, 1:nvqua):: xrq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndegr, 1:nsubq)::unksgq
real*8,dimension(1:nq+1, 1:nsubq)::uqsgc
real*8,dimension(1:4, 1:4)::prsgq
real*8,dimension(1:4, 1:4)::bqvp
real*8,dimension(1:4)::prqz
real*8,dimension(1:7, 1:4)::geoq_sub
real*8,dimension(2, 4, nsubq)::wfgsq
real*8,dimension(nvfac)::wegtf,posif
real*8::bqsg(ndegr)

!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::xcrho, ycrho
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx,rhovsg,dpvt
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::rcsg,scsg
real*8::dxp,dyp,xc,yc
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon
!
data eps   / 1.0d-6/
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Specify some Gauss points...
!
if(ncurv.eq.1)then
!...Local vertex No. of gauss points in one unit
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
fnqsg(1, 1) =  8; fnqsg(2, 1) = 1;  fnqsg(3, 1) = 13;
fnqsg(4, 1) =  -9;
fnqsg(5, 1) = -9; fnqsg(6, 1) =12;  fnqsg(7, 1) = 12; fnqsg(8, 1) =  16;

fnqsg(1, 2) =  9; fnqsg(2, 2) =13;
fnqsg(3, 2) =  2; fnqsg(4, 2) =  3;
fnqsg(5, 2) =  14;
fnqsg(6, 2) =-10;  fnqsg(7, 2) =-10; fnqsg(8, 2) = 9;

fnqsg(1, 3) =-11; fnqsg(2, 3) = 10; fnqsg(3, 3) = 10; fnqsg(4, 3) =14;
fnqsg(5, 3) =  4; fnqsg(6, 3) =  5; fnqsg(7, 3) = 15;
fnqsg(8, 3) =-11;

fnqsg(1, 4) = 16;
fnqsg(2, 4) =-12; fnqsg(3, 4) =-12; fnqsg(4, 4) = 11;
fnqsg(5, 4) = 11; fnqsg(6, 4) = 15;
fnqsg(7, 4) = 6;  fnqsg(8, 4) = 7;

!...4/6 for internal face...This is used for our work
wfgsq(1, 1, 1) = 4.d0/6.d0; wfgsq(2, 1, 1) = 4.d0/6.d0;
wfgsq(1, 2, 1) = 4.d0/3.d0; wfgsq(2, 2, 1) = 8.d0/3.d0;
wfgsq(1, 3, 1) = 1.d0;      wfgsq(2, 3, 1) = 1.d0;
wfgsq(1, 4, 1) = 8.d0/3.d0;      wfgsq(2, 4, 1) = 4.d0/3.d0;

wfgsq(1, 1, 2) = 8.d0/3.d0;      wfgsq(2, 1, 2) = 4.d0/3.d0;
wfgsq(1, 2, 2) = 4.d0/6.d0; wfgsq(2, 2, 2) = 4.d0/6.d0;
wfgsq(1, 3, 2) = 4.d0/3.d0; wfgsq(2, 3, 2) = 8.d0/3.d0;
wfgsq(1, 4, 2) = 1.d0;      wfgsq(2, 4, 2) = 1.d0;

wfgsq(1, 1, 3) = 1.d0;      wfgsq(2, 1, 3) = 1.d0;
wfgsq(1, 2, 3) = 8.d0/3.d0;      wfgsq(2, 2, 3) = 4.d0/3.d0;
wfgsq(1, 3, 3) = 4.d0/6.d0; wfgsq(2, 3, 3) = 4.d0/6.d0;
wfgsq(1, 4, 3) = 4.d0/3.d0; wfgsq(2, 4, 3) = 8.d0/3.d0

wfgsq(1, 1, 4) = 4.d0/3.d0; wfgsq(2, 1, 4) = 8.d0/3.d0;
wfgsq(1, 2, 4) = 1.d0;      wfgsq(2, 2, 4) = 1.d0;
wfgsq(1, 3, 4) = 8.d0/3.d0;      wfgsq(2, 3, 4) = 4.d0/3.d0;
wfgsq(1, 4, 4) = 4.d0/6.d0; wfgsq(2, 4, 4) = 4.d0/6.d0;

!...no internal contribution
wfgsq = 0.25d0*wfgsq

!...Vertex coordinate
xrq(1, 1) = -1.d0; xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0; xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0; xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0; xrq(2, 4) =  1.d0;
xrq(1, 5) =  0.d0; xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0; xrq(2, 6) =  0.d0;
xrq(1, 7) =  0.d0; xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0; xrq(2, 8) =  0.d0;
xrq(1, 9) =  0.d0; xrq(2, 9) =  0.d0;

elseif(ncurv.eq.2)then
!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8
!
fnqsg(1, 1) =  8; fnqsg(2, 1) =  1;  fnqsg(3, 1) = 17;  fnqsg(4, 1) =  -9;
fnqsg(5, 1) = -9; fnqsg(6, 1) = 16;  fnqsg(7, 1) = 16;  fnqsg(8, 1) =  24;

fnqsg(1, 2) =  9; fnqsg(2, 2) = 17;  fnqsg(3, 2) = 18;  fnqsg(4, 2) = -10;
fnqsg(5, 2) = -10;fnqsg(6, 2) =-10;  fnqsg(7, 2) =  9;  fnqsg(8, 2) = 9;

fnqsg(1, 3) = 10; fnqsg(2, 3) = 18;  fnqsg(3, 3) =  2;  fnqsg(4, 3) = 3;
fnqsg(5, 3) = 19; fnqsg(6, 3) =-11;  fnqsg(7, 3) =-11;  fnqsg(8, 3) =10;

fnqsg(1, 4) = 11; fnqsg(2, 4) = 11;  fnqsg(3, 4) = 11;  fnqsg(4, 4) = 19;
fnqsg(5, 4) = 20; fnqsg(6, 4) =-12;  fnqsg(7, 4) =-12;  fnqsg(8, 4)=-12;

fnqsg(1, 5) =-13; fnqsg(2, 5) = 12;  fnqsg(3, 5) = 12;  fnqsg(4, 5) = 20;
fnqsg(5, 5) =  4; fnqsg(6, 5) =  5;  fnqsg(7, 5) = 21;  fnqsg(8, 5)=-13;

fnqsg(1, 6) =-14; fnqsg(2, 6) =-14;  fnqsg(3, 6) = 13;  fnqsg(4, 6) = 13;
fnqsg(5, 6) = 13; fnqsg(6, 6) = 21;  fnqsg(7, 6) = 22;  fnqsg(8, 6) =-14;

fnqsg(1, 7) = 23; fnqsg(2, 7) =-15;  fnqsg(3, 7) =-15;  fnqsg(4, 7) = 14;
fnqsg(5, 7) = 14; fnqsg(6, 7) = 22;  fnqsg(7, 7) =  6;  fnqsg(8, 7) =  7;

fnqsg(1, 8) = 24; fnqsg(2, 8) =-16;  fnqsg(3, 8) =-16;  fnqsg(4, 8) =-16;
fnqsg(5, 8) =-16; fnqsg(6, 8) = 15;  fnqsg(7, 8) = 15;  fnqsg(8, 8) =23;

!...
wegtf(1) = 1.d0/12.d0;
wegtf(2) = 5.d0/12.d0;
wegtf(3) = 5.d0/12.d0;
wegtf(4) = 1.d0/12.d0;

!...Uniform 
wegtf = 1.d0

!...4/6 for internal face...This is used for our work
wfgsq(1, 1, 1) = wegtf(1);      wfgsq(2, 1, 1) = wegtf(1);
wfgsq(1, 2, 1) = wegtf(2)/1.d0; wfgsq(2, 2, 1) = wegtf(2)/1.d0
wfgsq(1, 3, 1) = 1.d0;          wfgsq(2, 3, 1) = 1.d0;
wfgsq(1, 4, 1) = wegtf(2)/1.d0; wfgsq(2, 4, 1) = wegtf(2)/1.d0;

wfgsq(1, 1, 2) = wegtf(2)/1.d0; wfgsq(2, 1, 2) = wegtf(2)/1.d0;
wfgsq(1, 2, 2) = wegtf(2)/1.d0; wfgsq(2, 2, 2) = wegtf(2)/1.d0;
wfgsq(1, 3, 2) = 1.d0;          wfgsq(2, 3, 2) = 1.d0;
wfgsq(1, 4, 2) = 1.d0;          wfgsq(2, 4, 2) = 1.d0;

wfgsq(1, 1, 3) = wegtf(2)/1.d0; wfgsq(2, 1, 3) = wegtf(2)/1.d0;
wfgsq(1, 2, 3) = wegtf(1);      wfgsq(2, 2, 3) = wegtf(1);
wfgsq(1, 3, 3) = wegtf(2)/1.d0; wfgsq(2, 3, 3) = wegtf(2)/1.d0;
wfgsq(1, 4, 3) = 1.d0;          wfgsq(2, 4, 3) = 1.d0

wfgsq(1, 1, 4) = 1.d0;          wfgsq(2, 1, 4) = 1.d0;
wfgsq(1, 2, 4) = wegtf(2)/1.d0; wfgsq(2, 2, 4) = wegtf(2)/1.d0;
wfgsq(1, 3, 4) = wegtf(2)/1.d0; wfgsq(2, 3, 4) = wegtf(2)/1.d0;
wfgsq(1, 4, 4) = 1.d0;          wfgsq(2, 4, 4) = 1.d0;

wfgsq(1, 1, 5) = 1.d0;          wfgsq(2, 1, 5) = 1.d0;
wfgsq(1, 2, 5) = wegtf(2)/1.d0; wfgsq(2, 2, 5) = wegtf(2)/1.d0
wfgsq(1, 3, 5) = wegtf(1);      wfgsq(2, 3, 5) = wegtf(1);
wfgsq(1, 4, 5) = wegtf(2)/1.d0; wfgsq(2, 4, 5) = wegtf(2)/1.d0;

wfgsq(1, 1, 6) = 1.d0;          wfgsq(2, 1, 6) = 1.d0;
wfgsq(1, 2, 6) = 1.d0;          wfgsq(2, 2, 6) = 1.d0;
wfgsq(1, 3, 6) = wegtf(2)/1.d0; wfgsq(2, 3, 6) = wegtf(2)/1.d0;
wfgsq(1, 4, 6) = wegtf(2)/1.d0; wfgsq(2, 4, 6) = wegtf(2)/1.d0;

wfgsq(1, 1, 7) = wegtf(2)/1.d0; wfgsq(2, 1, 7) = wegtf(2)/1.d0;
wfgsq(1, 2, 7) = 1.d0;          wfgsq(2, 2, 7) = 1.d0;
wfgsq(1, 3, 7) = wegtf(2)/1.d0; wfgsq(2, 3, 7) = wegtf(2)/1.d0;
wfgsq(1, 4, 7) = wegtf(1);      wfgsq(2, 4, 7) = wegtf(1)

wfgsq(1, 1, 8) = wegtf(2)/1.d0; wfgsq(2, 1, 8) = wegtf(2)/1.d0;
wfgsq(1, 2, 8) = 1.d0;          wfgsq(2, 2, 8) = 1.d0;
wfgsq(1, 3, 8) = 1.d0;          wfgsq(2, 3, 8) = 1.d0;
wfgsq(1, 4, 8) = wegtf(2)/1.d0; wfgsq(2, 4, 8) = wegtf(2)/1.d0;

!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;

!...Uniform space
!posif(2)=-1.d0/3.d0; posif(3)=1.d0/3.d0;
!...Vertex coordinate
xrq(1, 1) = -1.d0; xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0; xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0; xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0; xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;      xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;      xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;      xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;      xrq(2,12) = posif(2);

endif
!
!...Part II: Loop over every quad...
!
do 350 ie = 1,nquad !...(1)ie = 1,nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria

!...shape functions
dr = 1.0d0
ds = 1.0d0

!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...Basis function
do iv =1 ,nvqua
bq(1, iv) = 1.d0

if(npoly.ge.1)then
bq(2, iv) = (xrq(1, iv)-rc)/dr
bq(3, iv) = (xrq(2, iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
endif
endif
enddo

!...cell averaged value...
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

!... u, v, e at every vertex....
unknvq = 0.d0
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
enddo

!...Get density for sub-cells...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!...Update the coordinates for FE
!call getcoord_fe(ncurv, nvfac, nvqua, xpq)
!call getcoord_fe(ncurv, nvfac, nvqua, xpqi)
!
call getcoord_fe(ncurv, nvfac, nvqua, xpq)
call getcoord_fe(ncurv, nvfac, nvqua, xpqi)
!
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

!...Get density correction
call getrhosubcell_sms_general5(xcrho, ycrho, xpq,xpqi,unksgq, ielem)
!
!call getcoord_fe(ncurv, nvfac, nvqua, xpq)
!call getcoord_fe(ncurv, nvfac, nvqua, xpqi)

call getrhosubcell_daverg(rc,sc, geoel(19:21, ielem), xpq, unksgq,  unkno(:,:,ielem), aflim(:,ielem), uqsgc,ielem)

!...Output for debugging
!if(ielem.eq.2100)then
!print*,'Variabe',ielem,ivsg,unksgq(1,1:8)
!print*,'Variab2',ielem,ivsg,1.d0/uqsgc(1, 1:8)
!endif


!...II.1: Loop over sub-cells....
do isg = 1, nsubq

!...Normal vector for every quadrature point...
vnorm(1:2, 1, 1) = sign(1, fnqsg(1, isg))*gesgq(1:2,abs(fnqsg(1, isg)), ie); vnorm(  3, 1, 1) = gesgq(3,abs(fnqsg(1, isg)), ie);
vnorm(1:2, 2, 1) = sign(1, fnqsg(2, isg))*gesgq(1:2,abs(fnqsg(2, isg)), ie); vnorm(  3, 2, 1) = gesgq(3,abs(fnqsg(2, isg)), ie);

vnorm(1:2, 1, 2) = sign(1, fnqsg(3, isg))*gesgq(1:2,abs(fnqsg(3, isg)), ie); vnorm(  3, 1, 2) = gesgq(3,abs(fnqsg(3, isg)), ie);
vnorm(1:2, 2, 2) = sign(1, fnqsg(4, isg))*gesgq(1:2,abs(fnqsg(4, isg)), ie); vnorm(  3, 2, 2) = gesgq(3,abs(fnqsg(4, isg)), ie);

vnorm(1:2, 1, 3) = sign(1, fnqsg(5, isg))*gesgq(1:2,abs(fnqsg(5, isg)), ie); vnorm(  3, 1, 3) = gesgq(3,abs(fnqsg(5, isg)), ie);
vnorm(1:2, 2, 3) = sign(1, fnqsg(6, isg))*gesgq(1:2,abs(fnqsg(6, isg)), ie); vnorm(  3, 2, 3) = gesgq(3,abs(fnqsg(6, isg)), ie);

vnorm(1:2, 1, 4) = sign(1, fnqsg(7, isg))*gesgq(1:2,abs(fnqsg(7, isg)), ie); vnorm(  3, 1, 4) = gesgq(3,abs(fnqsg(7, isg)), ie);
vnorm(1:2, 2, 4) = sign(1, fnqsg(8, isg))*gesgq(1:2,abs(fnqsg(8, isg)), ie); vnorm(  3, 2, 4) = gesgq(3,abs(fnqsg(8, isg)), ie);

!...Get weighted area normal vector
do ifsg =1,2
vnorm(3, ifsg, 1) =  vnorm(  3, ifsg, 1)*wfgsq(ifsg, 1, isg);
vnorm(3, ifsg, 2) =  vnorm(  3, ifsg, 2)*wfgsq(ifsg, 2, isg);
vnorm(3, ifsg, 3) =  vnorm(  3, ifsg, 3)*wfgsq(ifsg, 3, isg);
vnorm(3, ifsg, 4) =  vnorm(  3, ifsg, 4)*wfgsq(ifsg, 4, isg);
enddo

!...Correct the density for Reimann input
do ivsg = 1,4

if(ipqsg(ivsg, isg).gt.12) cycle

if(ndens.eq.1) then
!print*,'rhovt',ie,ipqsg(ivsg, isg)
rhovt = 1.d0/unknvq(1, ipqsg(ivsg, isg))
rhovsg = rhovt +cdrho*(unksgq(1,isg)-1.d0/uqsgc(1, isg))
!rhovsg = 0.5d0*(rhovt+unksgq(1,isg))


!if(rhovsg.le.1d-6)then
!rhovsg = rhovt
!endif
!rhovsg = 1.d0/uqsgc(1, isg)
elseif(ndens.eq.3)then
print*,'ndens.eq.3 will be implemented in future in getriem_quad_sms!'
endif

uvtx = unknvq(2, ipqsg(ivsg, isg))
vvtx = unknvq(3, ipqsg(ivsg, isg))
evtx = unknvq(4, ipqsg(ivsg, isg))
!
pvtx = max(eps, (gamlg-1.d0)*rhovsg*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx

!...Output for debugging
!if(ielem.eq.2100.and.isg.eq.1.and.ivsg.eq.1)then
!print*,'Variabe',ielem,isg,ivsg,rhovsg,evtx,uvtx, vvtx,pvtx,rhovt,unksgq(1,isg)
!print*,'limit',afvec(:,:, ielem)
!endif
!
!if(ipq(ipqsg(ivsg, isg)).eq.31.or.ipq(ipqsg(ivsg, isg)).eq.107.or.ipq(ipqsg(ivsg, isg)).eq.242)then
!print*,'vaue',ipq(ipqsg(ivsg, isg)),ielem,isg,ivsg
!print*,'uvelo',unsgq(2, ivsg),unkno(:,2,ielem)
!print*,'vvelo',unsgq(3, ivsg),unkno(:,3,ielem)
!endif

!...Limiter
!...For high-order method(>P1), if the cell includes shock, then
!...the high-order terms are chopped, and only limit the linear DG(P1) part.

if(nlimi.eq.6.and.geoel(10, ielem).gt.10.d0)then
if(ndens.eq.1)then
 rhovsg = 1.d0/(rhomc + aflim(1, ielem)*(1.d0/rhovsg - rhomc))
!rhovsg = 1.d0/(rhomc + aflim(1, ielem)*(1.d0/rhovt - rhomc)) +cdrho*(unksgq(1,isg)-1.d0/uqsgc(1, isg))
elseif(ndens.eq.3)then
rhovsg = unksgq(1,isg) + aflim(1, ielem)*(rhovsg - unksgq(1,isg))
endif

!...Output for debugging
!print*,'Shock identification',ielem,unkno(:,:,ielem)
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bq(2, ipqsg(ivsg, isg)) + duds*bq(3, ipqsg(ivsg, isg))
vvtx = unkno(1,3,ielem)  + dvdr*bq(2, ipqsg(ivsg, isg)) + dvds*bq(3, ipqsg(ivsg, isg))
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr) !+ dpvt

!...Fitted pressure...
!pvtx = prsgq(ivsg, isg)

!...Updtae velocity
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx
endif

!...Get stress tensor at one vertex
sigma(1, 1, ivsg) = -pvtx
sigma(1, 2, ivsg) = 0.d0
sigma(2, 1, ivsg) = 0.d0
sigma(2, 2, ivsg) = -pvtx
!
!...Output for debugging
!if(ielem.eq.2)then
!print*,'smooth2',ielem,unkno(1,2:3,ielem)
!print*,'Variabe2',ielem,isg,ivsg,rhovsg,pvtx,unsgq(2:3, ivsg),sigma(1, 1, ivsg)
!endif

!...Get the a_c (unit vector)
aujmp(1:2, ivsg) = vlave(1:2, ipq(ipqsg(ivsg, isg))) - unsgq(2:3, ivsg)

acnx = aujmp(1, ivsg)
acny = aujmp(2, ivsg)
!
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ivsg) = 1.e-11!0.d0;
else
aujmp(1:2, ivsg) = aujmp(1:2, ivsg)/sqrt(acnx**2 + acny**2)
endif
!
aujmp(3, ivsg) = sqrt(acnx**2 + acny**2)
enddo

!...Get the variables at the center...
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedence coefficient...
do ivsg   = 1, 4
!
if(ipqsg(ivsg, isg).gt.12) cycle
!
dux= vlave(1, ipq(ipqsg(ivsg, isg)))-unsgq(2, ivsg)
duy= vlave(2, ipq(ipqsg(ivsg, isg)))-unsgq(3, ivsg)
!deltu = sqrt(dux**2 + duy**2)
do ifa = 1, 2
deltu = abs(dux*vnorm(1, ifa, ivsg) + duy*vnorm(2, ifa, ivsg))
murie(ifa, ivsg) = rhoct*sdctr + cimpd*rhoct*slpdu*deltu
!murie(ifa, ivsg) = uqsgc(1, isg)*sdctr + uqsgc(1, isg)*slpdu*deltu
!...The exact shock impedance of 2 shock model
!murie(ifa, ivsg) = cimpd*rhoct*slpdu*deltu/2.d0+&
!     rhoct*sqrt((slpdu*deltu/2.d0)**2 + gamlg*pctr/rhoct)
enddo
enddo

!...Feed the input into Riemann solver
do ivsg  = 1, 4

!...Local vertex No. of gauss points in one unit
iv = ipqsg(ivsg, isg)

if(ipqsg(ivsg, isg).gt.12) cycle
!
do ifa = 1, 2 !...Every corner consists of 2 faces...

!...Call Riemann solver...
call getriecoef_matrixnew(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:3, ivsg), &
unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:2, ivsg), &
!unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipq(iv)) = munacn(1:2, 1, ipq(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipq(iv)) = munacn(1:2, 2, ipq(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipq(iv)) = munacu(1:2, ipq(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipq(iv)) = snsigm(1:2, ipq(iv)) + snsigm_rie(1:2)

!...Output for debugging
if(ipq(iv).eq.2) print*,'p36 muacn(vv) post',ipq(iv),ie,ifa,isg,ivsg,snsigm_rie(1:2),&
vnorm(1:3, ifa, ivsg),&
sigma(1, 1, ivsg),&
murie(ifa, ivsg),unsgq(2:3, ivsg),ipq(iv)

!...Local variable...
munaclq(1:2, 1, ifa, ivsg, isg, ie) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, ifa, ivsg, isg, ie) =  munacn_rie(1:2, 2)
!
munaulq(1:2,    ifa, ivsg, isg, ie) =  munacu_rie(1:2)
!
snsigmlq(1:2,   ifa, ivsg, isg, ie)=  snsigm_rie(1:2)
!
enddo
enddo
!
enddo
!
350 enddo  !...(1)ie = 1,nquad

end subroutine getriem_quad_sms
!
!...Face integral for hybrid quad using SMS...
!
subroutine rhsifacedg_lag_sms(ipqua, unkno, ustar,fstarq, gesgq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:4, 1:nsubq, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem,isg,ivsg,ifsg
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua)  :: ipq
integer,dimension(8, nsubq)::fnqsg
integer,dimension(4, nsubq)::ipqsg
real*8,dimension(2,4,nsubq)::wfgsq,wfgsqm,wfgsqm1,wfgsqm2
real*8,dimension(1:3,1:2,1:4)::vnorm
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8, dimension(1:ndimn, 1:ndegr, 1:2, 1:nvqua)::lpnpq
real*8, dimension(1:ndimn,nvqua)::xrq
real*8::bq(1:ndegr,1:nvqua)
real*8,dimension(nvfac)::wegtf,posif

!...local real number
real*8::eps,c00,c05,c10,c20,c13,c16
real*8::dr,ds,rc,sc
real*8::dwav1,dwav2
real*8::anx, any
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
data c16   / 0.1666666666666666d0 /
!
!...Part I: preliminary setup
!
if(ncurv.eq.1)then
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
!
fnqsg(1, 1) =  8; fnqsg(2, 1) = 1;  fnqsg(3, 1) = 13;
fnqsg(4, 1) =  -9;
fnqsg(5, 1) = -9; fnqsg(6, 1) =12;  fnqsg(7, 1) = 12; fnqsg(8, 1) =  16;

fnqsg(1, 2) =  9; fnqsg(2, 2) =13;
fnqsg(3, 2) =  2; fnqsg(4, 2) =  3;
fnqsg(5, 2) =  14;
fnqsg(6, 2) =-10;  fnqsg(7, 2) =-10; fnqsg(8, 2) = 9;

fnqsg(1, 3) =-11; fnqsg(2, 3) = 10; fnqsg(3, 3) = 10; fnqsg(4, 3) =14;
fnqsg(5, 3) =  4; fnqsg(6, 3) =  5; fnqsg(7, 3) = 15;
fnqsg(8, 3) =-11;

fnqsg(1, 4) = 16;
fnqsg(2, 4) =-12; fnqsg(3, 4) =-12; fnqsg(4, 4) = 11;
fnqsg(5, 4) = 11; fnqsg(6, 4) = 15;
fnqsg(7, 4) = 6;  fnqsg(8, 4) = 7;
!
wfgsq(1, 1, 1) = 1.d0; wfgsq(2, 1, 1) = 1.d0;
wfgsq(1, 2, 1) = 1.d0; wfgsq(2, 2, 1) = 1.d0;
wfgsq(1, 3, 1) = 0.d0; wfgsq(2, 3, 1) = 0.d0;
wfgsq(1, 4, 1) = 1.d0; wfgsq(2, 4, 1) = 1.d0;

wfgsq(1, 1, 2) = 1.d0; wfgsq(2, 1, 2) =1.d0;
wfgsq(1, 2, 2) = 1.d0; wfgsq(2, 2, 2) =1.d0;
wfgsq(1, 3, 2) = 1.d0; wfgsq(2, 3, 2) =1.d0;
wfgsq(1, 4, 2) = 0.d0; wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 1, 3) = 0.d0; wfgsq(2, 1, 3) = 0.d0;
wfgsq(1, 2, 3) = 1.d0; wfgsq(2, 2, 3) = 1.d0;
wfgsq(1, 3, 3) = 1.d0; wfgsq(2, 3, 3) = 1.d0;
wfgsq(1, 4, 3) = 1.d0; wfgsq(2, 4, 3) = 1.d0;

wfgsq(1, 1, 4) = 1.d0; wfgsq(2, 1, 4) = 1.d0;
wfgsq(1, 2, 4) = 0.d0; wfgsq(2, 2, 4) = 0.d0;
wfgsq(1, 3, 4) = 1.d0; wfgsq(2, 3, 4) = 1.d0;
wfgsq(1, 4, 4) = 1.d0; wfgsq(2, 4, 4) = 1.d0;!
!...New
!
!
!...Taking the weight into consideration when calculate the ustar...
!
!wfgsq = 1.d0!0.25d0*wfgsq1
!
!...Conitnuity equation...
!
wfgsqm(1, 1, 1) = 4.d0/6.d0; wfgsqm(2, 1, 1) = 4.d0/6.d0;
wfgsqm(1, 2, 1) = 4.d0/3.d0; wfgsqm(2, 2, 1) = 0.d0;
wfgsqm(1, 3, 1) = 0.d0;      wfgsqm(2, 3, 1) = 0.d0;
wfgsqm(1, 4, 1) = 0.d0;      wfgsqm(2, 4, 1) = 4.d0/3.d0;

wfgsqm(1, 1, 2) = 0.d0/3.d0; wfgsqm(2, 1, 2) = 4.d0/3.d0;
wfgsqm(1, 2, 2) = 4.d0/6.d0; wfgsqm(2, 2, 2) = 4.d0/6.d0;
wfgsqm(1, 3, 2) = 4.d0/3.d0; wfgsqm(2, 3, 2) = 0.d0/3.d0;
wfgsqm(1, 4, 2) = 0.d0;      wfgsqm(2, 4, 2) = 0.d0;

wfgsqm(1, 1, 3) = 0.d0;      wfgsqm(2, 1, 3) = 0.d0;
wfgsqm(1, 2, 3) = 0.d0/3.d0; wfgsqm(2, 2, 3) = 4.d0/3.d0;
wfgsqm(1, 3, 3) = 4.d0/6.d0; wfgsqm(2, 3, 3) = 4.d0/6.d0;
wfgsqm(1, 4, 3) = 4.d0/3.d0; wfgsqm(2, 4, 3) = 0.d0/3.d0;

wfgsqm(1, 1, 4) = 4.d0/3.d0; wfgsqm(2, 1, 4) = 0.d0/3.d0;
wfgsqm(1, 2, 4) = 0.d0;      wfgsqm(2, 2, 4) = 0.d0;
wfgsqm(1, 3, 4) = 0.d0/3.d0; wfgsqm(2, 3, 4) = 4.d0/3.d0;
wfgsqm(1, 4, 4) = 4.d0/6.d0; wfgsqm(2, 4, 4) = 4.d0/6.d0;
!
wfgsqm = 0.25d0*wfgsqm

!...Vertex coordinate
xrq(1, 1) = -1.d0; xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0; xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0; xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0; xrq(2, 4) =  1.d0;
xrq(1, 5) =  0.d0; xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0; xrq(2, 6) =  0.d0;
xrq(1, 7) =  0.d0; xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0; xrq(2, 8) =  0.d0;
xrq(1, 9) =  0.d0; xrq(2, 9) =  0.d0;

elseif(ncurv.eq.2)then

ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8
!
fnqsg(1, 1) =  8; fnqsg(2, 1) =  1;  fnqsg(3, 1) = 17;  fnqsg(4, 1) =  -9;
fnqsg(5, 1) = -9; fnqsg(6, 1) = 16;  fnqsg(7, 1) = 16;  fnqsg(8, 1) =  24;

fnqsg(1, 2) =  9; fnqsg(2, 2) = 17;  fnqsg(3, 2) = 18;  fnqsg(4, 2) = -10;
fnqsg(5, 2) = -10;fnqsg(6, 2) =-10;  fnqsg(7, 2) =  9;  fnqsg(8, 2) = 9;

fnqsg(1, 3) = 10; fnqsg(2, 3) = 18;  fnqsg(3, 3) =  2;  fnqsg(4, 3) = 3;
fnqsg(5, 3) = 19; fnqsg(6, 3) =-11;  fnqsg(7, 3) =-11;  fnqsg(8, 3) =10;

fnqsg(1, 4) = 11; fnqsg(2, 4) = 11;  fnqsg(3, 4) = 11;  fnqsg(4, 4) = 19;
fnqsg(5, 4) = 20; fnqsg(6, 4) =-12;  fnqsg(7, 4) =-12;  fnqsg(8, 4)=-12;

fnqsg(1, 5) =-13; fnqsg(2, 5) = 12;  fnqsg(3, 5) = 12;  fnqsg(4, 5) = 20;
fnqsg(5, 5) =  4; fnqsg(6, 5) =  5;  fnqsg(7, 5) = 21;  fnqsg(8, 5)=-13;

fnqsg(1, 6) =-14; fnqsg(2, 6) =-14;  fnqsg(3, 6) = 13;  fnqsg(4, 6) = 13;
fnqsg(5, 6) = 13; fnqsg(6, 6) = 21;  fnqsg(7, 6) = 22;  fnqsg(8, 6) =-14;

fnqsg(1, 7) = 23; fnqsg(2, 7) =-15;  fnqsg(3, 7) =-15;  fnqsg(4, 7) = 14;
fnqsg(5, 7) = 14; fnqsg(6, 7) = 22;  fnqsg(7, 7) =  6;  fnqsg(8, 7) =  7;

fnqsg(1, 8) = 24; fnqsg(2, 8) =-16;  fnqsg(3, 8) =-16;  fnqsg(4, 8) =-16;
fnqsg(5, 8) =-16; fnqsg(6, 8) = 15;  fnqsg(7, 8) = 15;  fnqsg(8, 8) =23;
!
wegtf(1) = 1.d0/12.d0;
wegtf(2) = 5.d0/12.d0;
wegtf(3) = 5.d0/12.d0;
wegtf(4) = 1.d0/12.d0;

!...Newton-cotes
!wegtf(1) = 1.d0/8.d0;
!wegtf(2) = 3.d0/8.d0;
!wegtf(3) = 3.d0/8.d0;
!wegtf(4) = 1.d0/8.d0;
!
wfgsq = 1.d0
!...4/6 for internal face...This is used for our work
wfgsq(1, 3, 1) = 0.d0;          wfgsq(2, 3, 1) = 0.d0;

wfgsq(1, 3, 2) = 0.d0;          wfgsq(2, 3, 2) = 0.d0;
wfgsq(1, 4, 2) = 0.d0;          wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 4, 3) = 0.d0;          wfgsq(2, 4, 3) = 0.d0

wfgsq(1, 1, 4) = 0.d0;          wfgsq(2, 1, 4) = 0.d0;
wfgsq(1, 4, 4) = 0.d0;          wfgsq(2, 4, 4) = 0.d0;

wfgsq(1, 1, 5) = 0.d0;          wfgsq(2, 1, 5) = 0.d0;

wfgsq(1, 1, 6) = 0.d0;          wfgsq(2, 1, 6) = 0.d0;
wfgsq(1, 2, 6) = 0.d0;          wfgsq(2, 2, 6) = 0.d0;

wfgsq(1, 2, 7) = 0.d0;          wfgsq(2, 2, 7) = 0.d0;

wfgsq(1, 2, 8) = 0.d0;          wfgsq(2, 2, 8) = 0.d0;
wfgsq(1, 3, 8) = 0.d0;          wfgsq(2, 3, 8) = 0.d0;

!
!...4/6 for internal face...This is used for our work
wfgsqm(1, 1, 1) = wegtf(1);      wfgsqm(2, 1, 1) = wegtf(1);
wfgsqm(1, 2, 1) = wegtf(2)/2.d0; wfgsqm(2, 2, 1) = wegtf(2)/2.d0
wfgsqm(1, 3, 1) = 1.d0;          wfgsqm(2, 3, 1) = 1.d0;
wfgsqm(1, 4, 1) = wegtf(2)/2.d0; wfgsqm(2, 4, 1) = wegtf(2)/2.d0;

wfgsqm(1, 1, 2) = wegtf(2)/2.d0; wfgsqm(2, 1, 2) = wegtf(2)/2.d0;
wfgsqm(1, 2, 2) = wegtf(2)/2.d0; wfgsqm(2, 2, 2) = wegtf(2)/2.d0;
wfgsqm(1, 3, 2) = 1.d0;          wfgsqm(2, 3, 2) = 1.d0;
wfgsqm(1, 4, 2) = 1.d0;          wfgsqm(2, 4, 2) = 1.d0;

wfgsqm(1, 1, 3) = wegtf(2)/2.d0; wfgsqm(2, 1, 3) = wegtf(2)/2.d0;
wfgsqm(1, 2, 3) = wegtf(1);      wfgsqm(2, 2, 3) = wegtf(1);
wfgsqm(1, 3, 3) = wegtf(2)/2.d0; wfgsqm(2, 3, 3) = wegtf(2)/2.d0;
wfgsqm(1, 4, 3) = 1.d0;          wfgsqm(2, 4, 3) = 1.d0

wfgsqm(1, 1, 4) = 1.d0;          wfgsqm(2, 1, 4) = 1.d0;
wfgsqm(1, 2, 4) = wegtf(2)/2.d0; wfgsqm(2, 2, 4) = wegtf(2)/2.d0;
wfgsqm(1, 3, 4) = wegtf(2)/2.d0; wfgsqm(2, 3, 4) = wegtf(2)/2.d0;
wfgsqm(1, 4, 4) = 1.d0;          wfgsqm(2, 4, 4) = 1.d0;

wfgsqm(1, 1, 5) = 1.d0;          wfgsqm(2, 1, 5) = 1.d0;
wfgsqm(1, 2, 5) = wegtf(2)/2.d0; wfgsqm(2, 2, 5) = wegtf(2)/2.d0
wfgsqm(1, 3, 5) = wegtf(1);      wfgsqm(2, 3, 5) = wegtf(1);
wfgsqm(1, 4, 5) = wegtf(2)/2.d0; wfgsqm(2, 4, 5) = wegtf(2)/2.d0;

wfgsqm(1, 1, 6) = 1.d0;          wfgsqm(2, 1, 6) = 1.d0;
wfgsqm(1, 2, 6) = 1.d0;          wfgsqm(2, 2, 6) = 1.d0;
wfgsqm(1, 3, 6) = wegtf(2)/2.d0; wfgsqm(2, 3, 6) = wegtf(2)/2.d0;
wfgsqm(1, 4, 6) = wegtf(2)/2.d0; wfgsqm(2, 4, 6) = wegtf(2)/2.d0;

wfgsqm(1, 1, 7) = wegtf(2)/2.d0; wfgsqm(2, 1, 7) = wegtf(2)/2.d0;
wfgsqm(1, 2, 7) = 1.d0;          wfgsqm(2, 2, 7) = 1.d0;
wfgsqm(1, 3, 7) = wegtf(2)/2.d0; wfgsqm(2, 3, 7) = wegtf(2)/2.d0;
wfgsqm(1, 4, 7) = wegtf(1);      wfgsqm(2, 4, 7) = wegtf(1)

wfgsqm(1, 1, 8) = wegtf(2)/2.d0; wfgsqm(2, 1, 8) = wegtf(2)/2.d0;
wfgsqm(1, 2, 8) = 1.d0;          wfgsqm(2, 2, 8) = 1.d0;
wfgsqm(1, 3, 8) = 1.d0;          wfgsqm(2, 3, 8) = 1.d0;
wfgsqm(1, 4, 8) = wegtf(2)/2.d0; wfgsqm(2, 4, 8) = wegtf(2)/2.d0;

!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;

!...Uniform space
!posif(2)=-1.d0/3.d0; posif(3)=1.d0/3.d0;

!...Vertex coordinate
xrq(1, 1) = -1.d0; xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0; xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0; xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0; xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;      xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;      xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;      xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;      xrq(2,12) = posif(2);
endif
!
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria

!...The vertex constituting one cell...
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvqua
!
bq(1, iv) = 1.d0
!
if(npoly.ge.1)then
bq(2, iv) = (xrq(1, iv)-rc)/dr
bq(3, iv) = (xrq(2, iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
endif
endif
enddo

!...Initialize ulnpn, plnpn, elnpn
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
do isg = 1, nsubq

!...Normal vector for every face...
vnorm(1:2, 1, 1) = sign(1, fnqsg(1, isg))*gesgq(1:2,abs(fnqsg(1, isg)), ie); vnorm(  3, 1, 1) = gesgq(3,abs(fnqsg(1, isg)), ie);
vnorm(1:2, 2, 1) = sign(1, fnqsg(2, isg))*gesgq(1:2,abs(fnqsg(2, isg)), ie); vnorm(  3, 2, 1) = gesgq(3,abs(fnqsg(2, isg)), ie);

vnorm(1:2, 1, 2) = sign(1, fnqsg(3, isg))*gesgq(1:2,abs(fnqsg(3, isg)), ie); vnorm(  3, 1, 2) = gesgq(3,abs(fnqsg(3, isg)), ie);
vnorm(1:2, 2, 2) = sign(1, fnqsg(4, isg))*gesgq(1:2,abs(fnqsg(4, isg)), ie); vnorm(  3, 2, 2) = gesgq(3,abs(fnqsg(4, isg)), ie);

vnorm(1:2, 1, 3) = sign(1, fnqsg(5, isg))*gesgq(1:2,abs(fnqsg(5, isg)), ie); vnorm(  3, 1, 3) = gesgq(3,abs(fnqsg(5, isg)), ie);
vnorm(1:2, 2, 3) = sign(1, fnqsg(6, isg))*gesgq(1:2,abs(fnqsg(6, isg)), ie); vnorm(  3, 2, 3) = gesgq(3,abs(fnqsg(6, isg)), ie);

vnorm(1:2, 1, 4) = sign(1, fnqsg(7, isg))*gesgq(1:2,abs(fnqsg(7, isg)), ie); vnorm(  3, 1, 4) = gesgq(3,abs(fnqsg(7, isg)), ie);
vnorm(1:2, 2, 4) = sign(1, fnqsg(8, isg))*gesgq(1:2,abs(fnqsg(8, isg)), ie); vnorm(  3, 2, 4) = gesgq(3,abs(fnqsg(8, isg)), ie);
!
do ifsg =1,2
vnorm(3, ifsg, 1) =  vnorm(  3, ifsg, 1)*wfgsqm(ifsg, 1, isg);
vnorm(3, ifsg, 2) =  vnorm(  3, ifsg, 2)*wfgsqm(ifsg, 2, isg);
vnorm(3, ifsg, 3) =  vnorm(  3, ifsg, 3)*wfgsqm(ifsg, 3, isg);
vnorm(3, ifsg, 4) =  vnorm(  3, ifsg, 4)*wfgsqm(ifsg, 4, isg);
enddo
!
!...loop over every subgrid vertex...
!
do ivsg = 1, 4
!
if(ipqsg(ivsg, isg).gt.12) cycle
!
ulnpn(1:ndegr) = ulnpn(1:ndegr)  +&
ustar(1, ipq(ipqsg(ivsg, isg)))*vnorm(1, 1, ivsg)*vnorm(  3, 1, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg)) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*vnorm(2, 1, ivsg)*vnorm(  3, 1, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg)) +&
ustar(1, ipq(ipqsg(ivsg, isg)))*vnorm(1, 2, ivsg)*vnorm(  3, 2, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg)) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*vnorm(2, 2, ivsg)*vnorm(  3, 2, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg))
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, 1, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsqm(1, ivsg, isg) +&
fstarq(1, 2, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsqm(2, ivsg, isg)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, 1, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsqm(1, ivsg, isg) +&
fstarq(2, 2, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsqm(2, ivsg, isg)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(ipqsg(ivsg, isg)))*fstarq(1, 1, ivsg, isg,  ie)*&
bq(1:ndegr, ipqsg(ivsg, isg))*wfgsqm(1, ivsg, isg) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*fstarq(2, 1, ivsg, isg,  ie)*&
bq(1:ndegr, ipqsg(ivsg, isg))*wfgsqm(1, ivsg, isg) +&
ustar(1, ipq(ipqsg(ivsg, isg)))*fstarq(1, 2, ivsg, isg, ie)*&
bq(1:ndegr, ipqsg(ivsg, isg))*wfgsqm(2, ivsg, isg)  +&
ustar(2, ipq(ipqsg(ivsg, isg)))*fstarq(2, 2, ivsg, isg,  ie)*&
bq(1:ndegr, ipqsg(ivsg, isg))*wfgsqm(2, ivsg, isg)

!
!...Output foe debugging
!if(ie==2100)  print*,'rhs iface2',isg,ivsg,ipq(ipqsg(ivsg, isg)),ie,fstarq(1, 1:2, ivsg, isg, ie)
enddo
!
enddo
!
!...Distribute to every corner...
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)

!...Output foe debugging
!if(ie==1)  print*,'rhs iface',ielem, ie, plnpn(1, 1:ndegr)
650 enddo
!
end subroutine rhsifacedg_lag_sms
!
!....domain integral for arbitrary hybrid curved quad cells
!
subroutine rhsdomndg_lag_quadc(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8, dimension(1:nq+1, 1:nsize),      intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),     intent(in)::afvec
real*8,dimension(1:ndegr,1:nq,1:ncell),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem
!
!...local integer array
!
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds, bv
real*8:: unknod(1:nq)
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::rhoad,uadv,vadv,eadv,rhoma
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr
real*8:: rhoi, rhon
!
!
data eps   / 1.0d-6/
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
call ruqope_lobatto(2, ngausdq, posiq, weighq)
!
!...Loop over quad
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

!...Points consitituting one element...
xpq(1, 1:nvqua) = coord(1,ipqua(1:nvqua, ie))
xpq(2, 1:nvqua) = coord(2,ipqua(1:nvqua, ie))

!...Get the points for FEM
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!...Geometry parameters for reference cell...
dr = 1.d0
ds = 1.d0
!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...The derivatives of basis function...
!...Here dbdx means dbd(xsi), dbdy means dbd(eta)
dbdr(1)= 0.d0
dbdr(2)= 1.d0/dr
dbdr(3)= 0.d0

dbds(1)= 0.d0
dbds(2)= 0.d0
dbds(3)= 1.0/ds

!...Gauss loop
do ig = 1,ngausdq !...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
wi  = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi

!...The high-order derivative of basis function
if(npoly.eq.2)then
dbdr(4)= (r-rc)/dr**2
dbdr(5)= 0.d0
dbdr(6)= (s-sc)/dr/ds

dbds(4)= 0.d0
dbds(5)= (s-sc)/ds**2
dbds(6)= (r-rc)/dr/ds
endif

!...Jacobian transformation matrix
jacbf(1, 1) = dxdr; jacbf(1, 2) = dxds
jacbf(2, 1) = dydr; jacbf(2, 2) = dyds

!...Cofactor matrix of Jacobian transformation matrix
jacbg(1, 1) = dyds; jacbg(1, 2) =-dydr
jacbg(2, 1) =-dxds; jacbg(2, 2) = dxdr

!...Calculate G dot dbdx or dbdy
do ideg = 1, ndegr
gdshp(1, ideg) = jacbg(1, 1)*dbdr(ideg) + jacbg(1, 2)*dbds(ideg)
gdshp(2, ideg) = jacbg(2, 1)*dbdr(ideg) + jacbg(2, 2)*dbds(ideg)
enddo

!...Gauss points...
xg = r
yg = s

!...Basis function for solutions...
b(1) = 1.d0
b(2) = (xg-rc)/dr
b(3) = (yg-sc)/ds
!DGP2
if(npoly.eq.2)then
b(4) = 0.5d0*b(2)*b(2) - geoel(19, ielem)
b(5) = 0.5d0*b(3)*b(3) - geoel(20, ielem)
b(6) =       b(2)*b(3) - geoel(21, ielem)
endif

!...Solution at the Gauss points...
unknod = 0.d0
do ideg =1,ndegr
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo

!...Primitive variables...
if(ndens.eq.1)then
rhoma = unknod(1)
rhoad  = 1.d0/rhoma
elseif(ndens.eq.2)then
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quadcurv(rhoi, xpqi)!
call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)

!call getrhoig_quad(rhoi,r,s, xpqi(1:2,1:4))!
!call getdensity_quadllnl(r, s, xpq(1:2,1:4), xpqi(1:2,1:4), rhoi, rhon)
rhoma = 1.d0/rhon
rhoad = rhon
elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
bv(1) = 1.d0
bv(2) = (xg-rcv)/dr
bv(3) = (yg-scv)/ds
!
unknod(1) =0.d0
do ideg = 1,mdegr
unknod(1) = unknod(1) + unkno(ideg,1,ielem)*bv(ideg)
enddo
!
rhoma = 1.d0/unknod(1)
rhoad  = unknod(1)
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = max(eps,(gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2)))

!if(ie==2100) print*,'rhs iface idegr1',ie,ig,uadv,vadv,unkno(2,2,ielem)


!...Limiting
if(nlimi.eq.6.and.geoel(10,ielem).gt.10.d0)then
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

!
if(ndens.eq.1)then
rhoma = rhomc + aflim(1 ,ielem)*(rhoma - rhomc)
rhoad = 1.d0/rhoma
elseif(ndens.eq.2)then
rhoad = rhoct + aflim(1 ,ielem)*(rhoad - rhoct)
elseif(ndens.eq.3)then
rhoad = rhoct + aflim(1 ,ielem)*(rhoad - rhoct)
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)

uadv = unkno(1,2,ielem)  + dudr*b(2) + duds*b(3)
vadv = unkno(1,3,ielem)  + dvdr*b(2) + dvds*b(3)
!
!eadv = ectr + aflim(5, ielem)*(eadv- ectr)
!...Null mode at gauss points....
uadv = uadv !- vgnul(1)
vadv = vadv !- vgnul(2)
pres = pctr + aflim(4, ielem)*(pres- pctr)
!pres = max(eps, (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2)))
!
endif

!...High-order
do ideg = 1, ndegr
fluxd(ideg,1) = gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv
fluxd(ideg,2) = gdshp(1, ideg)*(-pres)
fluxd(ideg,3) = gdshp(2, ideg)*(-pres)
fluxd(ideg,4) = (gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv)*(-pres)
enddo

!finally, scatter the contribution to the RHS
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) - fluxd(ideg,1:nq)*djaco
enddo

!...Output for debugging
!if(ie==2100) print*,'rhs iface idegr',ie,ig,gdshp(1:2,4:5),uadv,vadv,djaco
!
enddo !...(2)ig = 1,ngausd
!
650 enddo

end subroutine rhsdomndg_lag_quadc
!
!...Get the node coordinates based on the shape function
!
subroutine getcoord_fe(ncurv, nvfac, nvqua, xpq)
implicit none
integer, intent(in):: ncurv, nvfac, nvqua
real*8, dimension(1:2, nvqua), intent(inout):: xpq

!...Local
integer,dimension(1:nvfac, 1:4) :: ipfq
real*8, dimension(2, nvfac)  :: xpf
real*8, dimension(1:2, nvqua):: xpq_fe
!
integer::ifa
real*8::rf
real*8::detmt
real*8::rhssp(2, 2),amtsp(2,2)
real*8, dimension(nvfac, 2)::shpf, dshprf
!
!...Part I: Preliminarty setup
!
ipfq(1, 1)= 1; ipfq(2, 1)= 2; ipfq(3, 1)= 5; ipfq(4, 1)= 9;
ipfq(1, 2)= 2; ipfq(2, 2)= 3; ipfq(3, 2)= 6; ipfq(4, 2)=10;
ipfq(1, 3)= 3; ipfq(2, 3)= 4; ipfq(3, 3)= 7; ipfq(4, 3)=11;
ipfq(1, 4)= 4; ipfq(2, 4)= 1; ipfq(3, 4)= 8; ipfq(4, 4)=12;

!...Give xpq to xpq_fe
xpq_fe = xpq

!
!...Part II: Update the edge points
!
do ifa = 1, 4
!
xpf(1, 1:nvfac) = xpq_fe(1, ipfq(1:nvfac,ifa))
xpf(2, 1:nvfac) = xpq_fe(2, ipfq(1:nvfac,ifa))

!...Call shape functions
rf =-sqrt(5.d0)/5.d0

!rf =-1.d0/3.d0
call getshapfct_edge(ncurv,nvfac,shpf(:,1), dshprf(:,1), rf)

rf = sqrt(5.d0)/5.d0

!rf = 1.d0/3.d0
call getshapfct_edge(ncurv,nvfac,shpf(:,2), dshprf(:,2), rf)
!
rhssp(1, 1) = xpf(1, 3) - shpf(1, 1)*xpf(1, 1) - shpf(2, 1)*xpf(1, 2)
rhssp(2, 1) = xpf(1, 4) - shpf(1, 2)*xpf(1, 1) - shpf(2, 2)*xpf(1, 2)
!
rhssp(1, 2) = xpf(2, 3) - shpf(1, 1)*xpf(2, 1) - shpf(2, 1)*xpf(2, 2)
rhssp(2, 2) = xpf(2, 4) - shpf(1, 2)*xpf(2, 1) - shpf(2, 2)*xpf(2, 2)
!
detmt = shpf(3 ,1)*shpf(4 ,2) - shpf(3 ,2)*shpf(4 ,1)
!
amtsp(1 ,1)= shpf(4 ,2); amtsp(1 ,2)=-shpf(4 ,1)
amtsp(2 ,1)=-shpf(3 ,2); amtsp(2 ,2)= shpf(3 ,1)
!
amtsp = amtsp/detmt

!...New coordinates
xpf(1, 3) = amtsp(1, 1)*rhssp(1, 1) + amtsp(1, 2)*rhssp(2, 1)
xpf(1, 4) = amtsp(2, 1)*rhssp(1, 1) + amtsp(2, 2)*rhssp(2, 1)

xpf(2, 3) = amtsp(1, 1)*rhssp(1, 2) + amtsp(1, 2)*rhssp(2, 2)
xpf(2, 4) = amtsp(2, 1)*rhssp(1, 2) + amtsp(2, 2)*rhssp(2, 2)

!...Update the physical coordinates for FE shape functions...
xpq_fe(1:2, ipfq(3,ifa)) = xpf(1:2, 3)
xpq_fe(1:2, ipfq(4,ifa)) = xpf(1:2, 4)
!
enddo

!
!...Part III: Update xpq
!
xpq = xpq_fe
end subroutine getcoord_fe
!
!...Calculate geoel initiallyfor general cuvrecd mesh...
!
subroutine getgeoel_initial_lag_curv(iptri, ipqua, geoel, coord)
use constant
implicit none
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:ngeel, 1:nsize), intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord

!...local array...
real*8,dimension(1:2, 1:nvtri)::xp
real*8,dimension(1:2, 1:nvqua)::xpq
real*8,dimension(1:nvtri)::shp, dspr, dsps
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8:: weight(ngausd_geo), posit(2, ngausd_geo)
real*8:: weighq(ngausd_geoq), posiq(2, ngausd_geoq)
!...local real number
real*8:: dxdr,dxds,dydr,dyds
real*8:: xc, yc, xcel, ycel
real*8:: eps,c00,c10,c05,c20
real*8:: r, s, djaco, volel, masel
real*8:: xcerz, ycerz, maselrz, volelrz
real*8:: wi, xcent, ycent,xgaus, ygaus, rhog
real*8::radie, radii,radie2,radii2,radic2, radig2,sentr,rhoin,rhoex,rho0ba
real*8::preex,prein
real*8::rcoef !...Coefficient R of RZ or XY coordinate system...
real*8::areel !...2D area of the cell
real*8::xceli, yceli !...The initial cell center(volume)
!
integer::ielem, igaus, ishp
integer::ie
!
data eps / 1.0d-06 /
data c00 / 0.0d0 /
data c10 / 1.0d0 /
data c05 / 0.5d0 /
data c20 / 2.0d0 /
!
!...Find weight and position for gauss points...
call rutope(2, ngausd_geo, posit, weight)
call ruqope(2, ngausd_geoq, posiq, weighq)
!
!...1st loop to find center and volulme...
!
do ie  = 1, ntria !...(1)ie = 1,nelem
!
ielem = ie
!
xp(1, 1:nvtri) = coord(1,iptri(1:nvtri, ie))
xp(2, 1:nvtri) = coord(2,iptri(1:nvtri, ie))

!...Get the coordinates for FE


!...Physical coordinate
r = 1.d0/3.d0; s= 1.d0/3.d0
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nptri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo

!...
xcent = 0.d0
ycent = 0.d0
xcel  = 0.d0
ycel  = 0.d0
volel = 0.d0
masel = 0.d0
areel = 0.d0
!
xcerz = 0.d0
ycerz = 0.d0
!
xceli = 0.d0
yceli = 0.d0
!
maselrz = 0.d0
volelrz = 0.d0
!
do igaus =1,ngausd_geo
!
r  = posit(1,igaus)
s  = posit(2,igaus)
wi  = weight(igaus)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
!
dspr(1) = c10-4.d0*(c10-r-s)
dspr(2) = -1.d0 + 4.d0*r
dspr(3) = 0.d0
dspr(4) = 4.d0*(1.d0-2.d0*r-s)
dspr(5) = 4.d0*s
dspr(6) =-4.d0*s

!
dsps(1) = c10-4.d0*(c10-r-s)
dsps(2) = 0.d0
dsps(3) = -1.d0 + 4.d0*s
dsps(4) = -4.d0*r
dsps(5) =  4.d0*r
dsps(6) =  4.d0*(1.d0-r-2.d0*s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*xp(1,ishp)
dxds = dxds + dsps(ishp)*xp(1,ishp)

dydr = dydr + dspr(ishp)*xp(2,ishp)
dyds = dyds + dsps(ishp)*xp(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!volel = volel + djaco
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhog,  xgaus, ygaus, xc, yc)
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
!...Restore r, s
!
r  = posit(1,igaus)
s  = posit(2,igaus)
!
!...Mass...
!
masel = masel + rhog*djaco*rcoef
volel = volel + djaco*rcoef
!
xcent = xcent + rhog*r*djaco*rcoef !...mass center at the reference element...
ycent = ycent + rhog*s*djaco*rcoef
!
xcel = xcel + djaco*r*rcoef!xgaus
ycel = ycel + djaco*s*rcoef!ygaus
!
if(nrz.eq.2)then
!
maselrz = maselrz + rhog*djaco
volelrz = volelrz + djaco
!
xcerz = xcerz + rhog*r*djaco !...mass center at the reference element...
ycerz = ycerz + rhog*s*djaco
!
xceli = xceli + djaco*r
yceli = yceli + djaco*s

areel = areel + djaco !...The area in R-Z coordinates
endif
!
!if(ielem.eq.1882) print*,'geoel',xcent,ycent,xcel,ycel,r,s,djaco
!
enddo
!
geoel(1, ielem) = xcent/masel
geoel(2, ielem) = ycent/masel
geoel(3, ielem) = volel
geoel(4, ielem) = masel
!
geoel(7, ielem) = xcel/volel
geoel(8, ielem) = ycel/volel
!
if(nrz.eq.2)then
geoel(1, ielem) = xcerz/maselrz
geoel(2, ielem) = ycerz/maselrz
geoel(3, ielem) = volel!rz
geoel(4, ielem) = maselrz
geoel(7, ielem) = xceli/areel
geoel(8, ielem) = yceli/areel
geoel(11, ielem) = volel/areel !...Cell averaged radius
geoel(12, ielem) = areel
endif
!
!...Other components are same as those of Eulerian framework...
!
geoel(5, ielem) = maxval(xp(2, 1:3))-minval(xp(2, 1:3))
!
!
!if(ielem.eq.1884) print*,'geoel',ielem,geoel(1:2, ielem),xcent,ycent,geoel(7:8, ielem)
!    print*,'old geoel',ielem, geoel(1:5, ielem)
!    print*,'new geoel',ielem, xcent/volel, ycent/volel, volel, maxval(xp(1, 1:3))-minval(xp(1, 1:3)), &
!                         maxval(xp(2, 1:3))- minval(xp(2, 1:3))
!
enddo !...(1)ie = 1,nelem
!
!...2nd loop to find  mass center and volulme for quads...
!
do ie  = 1, nquad !...(1)ie = 1,nelem
!
!...global numbering
!
ielem = ie + ntria
!
xpq(1, 1:nvqua) = coord(1,ipqua(1:nvqua, ie))
xpq(2, 1:nvqua) = coord(2,ipqua(1:nvqua, ie))

!...Get coordinates for FE
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!...Physical coordinates
r=0.d0;s=0.d0
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

xc = 0.d0
yc = 0.d0
!
do ishp = 1, nvqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo

!...
xcent = 0.d0
ycent = 0.d0
xcel = 0.d0
ycel = 0.d0
volel = 0.d0
masel = 0.d0
areel = 0.d0
!
xcerz = 0.d0
ycerz = 0.d0
!
xceli = 0.d0
yceli = 0.d0
maselrz = 0.d0
volelrz = 0.d0
!
do igaus =1,ngausd_geoq
!
r  = posiq(1,igaus)
s  = posiq(2,igaus)
wi  = weighq(igaus)
!
!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!volel = volel + djaco
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhog,  xgaus, ygaus, xc, yc)

!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!...Restore r, s
r  = posiq(1,igaus)
s  = posiq(2,igaus)
!
!...Mass...
!
masel = masel + rhog*djaco*rcoef
volel = volel + djaco*rcoef
!
!if(ielem.eq.1) print*,'bad',rhog,djaco,wi
!
xcent = xcent + rhog*r*djaco*rcoef !...mass center at the reference element...
ycent = ycent + rhog*s*djaco*rcoef
!
xcel = xcel + djaco*r*rcoef!xgaus
ycel = ycel + djaco*s*rcoef!ygaus
!
if(nrz.eq.2)then
!
maselrz = maselrz + rhog*djaco
volelrz = volelrz + djaco
!
xcerz = xcerz + rhog*r*djaco !...mass center at the reference element...
ycerz = ycerz + rhog*s*djaco
!
xceli = xceli + djaco*r
yceli = yceli + djaco*s
areel = areel + djaco !...The area in R-Z coordinates
endif
!
!if(ie.eq.1) print*,'ieeq1',djaco,r,s
!
enddo
!
geoel(1, ielem) = xcent/masel
geoel(2, ielem) = ycent/masel
geoel(3, ielem) = volel
!
geoel(4, ielem) = masel
geoel(5, ielem) = maxval(xpq(2, 1:4))-minval(xpq(2, 1:4))
!
geoel(7, ielem) = xcel/volel
geoel(8, ielem) = ycel/volel
!
if(nrz.eq.2)then
geoel(1, ielem) = xcerz/maselrz
geoel(2, ielem) = ycerz/maselrz
geoel(3, ielem) = volel!rz
geoel(4, ielem) = maselrz
geoel(7, ielem) = xceli/areel
geoel(8, ielem) = yceli/areel
geoel(11, ielem) = volel/areel !...Cell averaged radius
geoel(12, ielem) = areel
endif
!
enddo !...(1)ie = 1,nelem
!
print*,'geoel',geoel(3:4,1)
!
end subroutine getgeoel_initial_lag_curv
!
!...Find terms in geoel for high-order for the gradient deformation for general curved cells...
!
subroutine getgeoelho_initial_lag_curv(iptri, ipqua, geoel, coord)
use constant
implicit none
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:ngeel, 1:nsize), intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
!
!...local array...
real*8,dimension(1:2, 1:nvtri)::xp
real*8,dimension(1:2, 1:nvqua)::xpq
real*8,dimension(1:nvtri)::shp, dspr, dsps
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8:: weight(ngausd_geo), posit(2, ngausd_geo)
real*8:: weighq(ngausd_geoq), posiq(2, ngausd_geoq)
real*8:: bt(ndegr), bq(ndegr)
real*8:: bti(ndegr), bqi(ndegr)
!...local real number
real*8:: dxdr,dxds,dydr,dyds
real*8:: xc, yc, xcel, ycel
real*8:: eps,c00,c10,c05,c20
real*8:: r, s, djaco, volel, masel
real*8:: dr, ds, rc,sc, rci, sci
real*8:: xcerz, ycerz, maselrz, volelrz
real*8:: wi, xcent, ycent,xgaus, ygaus, rhog
real*8::radie, radii,radie2,radii2,radic2, radig2,sentr,rhoin,rhoex,rho0ba
real*8::preex,prein
real*8::rcoef !...Coefficient R of RZ or XY coordinate system...
real*8::areel !...2D area of the cell
real*8::xceli, yceli !...The initial cell center(volume)
!
integer::ielem, igaus, ishp
integer::ie
!
data eps / 1.0d-06 /
data c00 / 0.0d0 /
data c10 / 1.0d0 /
data c05 / 0.5d0 /
data c20 / 2.0d0 /
!
!...Find weight and position for gauss points...
call rutope(2, ngausd_geo, posit, weight)
call ruqope(2, ngausd_geoq, posiq, weighq)
!
!...Part I: Get the basis function terms for DGP2
!
!...I.1: Triangle
do ie  = 1, ntria !...(1)ie = 1,nelem
!
ielem = ie
!
xp(1, 1:nvtri) = coord(1,iptri(1:nvtri, ie))
xp(2, 1:nvtri) = coord(2,iptri(1:nvtri, ie))

!...Mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...Initial Volume center
rci = geoel(7, ielem)
sci = geoel(8, ielem)
!
dr = 0.5d0
ds = 0.5d0

!...Cell volume(initial) and mass
volel = geoel(3, ielem)
masel = geoel(4, ielem)

!...Physical coordinate for the mass center (rc, sc)
r = rc; s= sc
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nvtri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo

!...Zero out bt, bti
bt  = 0.d0
bti = 0.d0
!
do igaus =1,ngausd_geo
!
r  = posit(1,igaus)
s  = posit(2,igaus)
wi  = weight(igaus)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
!
dspr(1) = c10-4.d0*(c10-r-s)
dspr(2) = -1.d0 + 4.d0*r
dspr(3) = 0.d0
dspr(4) = 4.d0*(1.d0-2.d0*r-s)
dspr(5) = 4.d0*s
dspr(6) =-4.d0*s

!
dsps(1) = c10-4.d0*(c10-r-s)
dsps(2) = 0.d0
dsps(3) = -1.d0 + 4.d0*s
dsps(4) = -4.d0*r
dsps(5) =  4.d0*r
dsps(6) =  4.d0*(1.d0-r-2.d0*s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*xp(1,ishp)
dxds = dxds + dsps(ishp)*xp(1,ishp)

dydr = dydr + dspr(ishp)*xp(2,ishp)
dyds = dyds + dsps(ishp)*xp(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!volel = volel + djaco
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo
!
!...Get initial density
call getrhog_initial(rhog,  xgaus, ygaus, xc, yc)
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus

!...Restore r, s
r  = posit(1,igaus)
s  = posit(2,igaus)

!...Basis function
bt(1) = 1.d0
bt(2) = (r-rc)/dr
bt(3) = (s-sc)/ds
!
bt(4) = bt(4) + 0.5d0*rhog*bt(2)*bt(2)*djaco*rcoef
bt(5) = bt(5) + 0.5d0*rhog*bt(3)*bt(3)*djaco*rcoef
bt(6) = bt(6) + rhog*bt(2)*bt(3)*djaco*rcoef

!...Basis function at the initial domn
bti(1) = 1.d0
bti(2) = (r-rci)/dr
bti(3) = (s-sci)/ds
!
bti(4) = bti(4) + 0.5d0*bti(2)*bti(2)*djaco*rcoef
bti(5) = bti(5) + 0.5d0*bti(3)*bti(3)*djaco*rcoef
bti(6) = bti(6) +       bti(2)*bti(3)*djaco*rcoef

!
if(nrz.eq.2)then
print*,'DG(P2) for AW-RZ will be implemented in the future!'
stop
endif
!
enddo
!
geoel(19, ielem) = bt(4)/masel
geoel(20, ielem) = bt(5)/masel
geoel(21, ielem) = bt(6)/masel
!
geoel(22, ielem) = bti(4)/volel
geoel(23, ielem) = bti(5)/volel
geoel(24, ielem) = bti(6)/volel
!
if(nrz.eq.2)then
print*,'DG(P2) for AW-RZ will be implemented in the future!'
stop
endif
!
enddo !...(1)ie = 1,ntria
!
!...I.2:Quads...
!
do ie  = 1, nquad !...(1)ie = 1,nelem
!
!...global numbering
!
ielem = ie + ntria
!
xpq(1, 1:nvqua) = coord(1,ipqua(1:nvqua, ie))
xpq(2, 1:nvqua) = coord(2,ipqua(1:nvqua, ie))

!...Get coordinates for FE
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!...Mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...Initial Volume center
rci = geoel(7, ielem)
sci = geoel(8, ielem)
!
dr = 1.d0
ds = 1.d0

!...Cell volume(initial) and mass
volel = geoel(3, ielem)
masel = geoel(4, ielem)

!...Physical coordinates
r=rc;s=sc
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nvqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo

!...Zero out bq, bqi
bq  = 0.d0
bqi = 0.d0
!
do igaus =1,ngausd_geoq
!
r  = posiq(1,igaus)
s  = posiq(2,igaus)
wi  = weighq(igaus)
!
!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!volel = volel + djaco
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
!...Get initial density
call getrhog_initial(rhog,  xgaus, ygaus, xc, yc)
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus

!...Restore r, s
r  = posiq(1,igaus)
s  = posiq(2,igaus)

!...Basis function
bq(1) = 1.d0
bq(2) = (r-rc)/dr
bq(3) = (s-sc)/ds
!
bq(4) = bq(4) + 0.5d0*rhog*bq(2)*bq(2)*djaco*rcoef
bq(5) = bq(5) + 0.5d0*rhog*bq(3)*bq(3)*djaco*rcoef
bq(6) = bq(6) + rhog*bq(2)*bq(3)*djaco*rcoef

!...Basis function for initial domain
bqi(1) = 1.d0
bqi(2) = (r-rci)/dr
bqi(3) = (s-sci)/ds
!
bqi(4) = bqi(4) + 0.5d0*bqi(2)*bqi(2)*djaco*rcoef
bqi(5) = bqi(5) + 0.5d0*bqi(3)*bqi(3)*djaco*rcoef
bqi(6) = bqi(6) +       bqi(2)*bqi(3)*djaco*rcoef
!
if(nrz.eq.2)then
print*,'DG(P2) for AW-RZ will be implemented in the future!'
stop
endif
!
!if(ie.eq.1) print*,'ieeq1',djaco,r,s
!
enddo
!
geoel(19, ielem) = bq(4)/masel
geoel(20, ielem) = bq(5)/masel
geoel(21, ielem) = bq(6)/masel
!
geoel(22, ielem) = bqi(4)/volel
geoel(23, ielem) = bqi(5)/volel
geoel(24, ielem) = bqi(6)/volel
!
if(nrz.eq.2)then
print*,'DG(P2) for AW-RZ will be implemented in the future!'
stop
endif
!
enddo !...(1)ie = 1,nquad
!
end subroutine getgeoelho_initial_lag_curv
!
!...Subroutine: Set the initial flow field for general curved mesh...
!
subroutine inifield_curv(unkno,uchar,geoel,coord,inpoel,iptri,ipqua)
use constant
implicit none
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(out)::unkno
real*8,dimension(1:nq)::uchar
real*8,dimension(1:ngeel,1:nsize),      intent(in)::geoel
real*8,dimension(1:ndimn, 1:npoin),intent(in)::coord
integer,  dimension(1:nvtri, 1:ntria), intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!...Local
!...local real array
real*8,dimension(1:2, 1:nvtri)::xp
real*8,dimension(1:nvtri):: shp, dspr, dsps
real*8,dimension(1:2, 1:nvqua)::xpq
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
!
real*8::unkct(1:ndegr, 1:nq), unkin(1:nq)
!...local real

real*8:: dr, ds
real*8:: dxdr,dxds,dydr,dyds
real*8:: volel,wi,xc,yc,xgaus,ygaus
real*8:: djaco
real*8:: rhomvel, rhouvel, rhovvel, rhoevel
real*8:: du0dx,du0dy,dv0dx,dv0dy
real*8:: de0dx,de0dy
real*8:: drhomdx, drhomdy
real*8:: drhomdx2,drhomdy2,drhomdxy,du0dx2,du0dy2,du0dxy,dv0dx2,dv0dy2,dv0dxy,de0dx2,de0dy2,de0dxy
real*8:: masel,denom, volsed
real*8:: rho, uini, vini ,pini

real*8:: c00,c10,c05,c20,eps
real*8::r, s, rc, sc
real*8:: xcel, ycel
real*8:: rcoef
!...local integer
integer*4::ie,id, ig, ishp
integer::ielem, igaus
!
real*8:: weight(ngausd), posit(2, ngausd)
real*8:: weighq(ngausdqi), posiq(2, ngausdqi)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
call rutope(2, ngausd, posit, weight)
call ruqope(2, ngausdqi, posiq, weighq)
!
!...Zero out everything...
!
unkno = 0.d0
!
volsed = 0.d0
!
!...initilization of flow field...
!
if(nmeth.eq.1)then !...nmeth=1 for Euleriam frame...
!
print*,'...................Notice................... !'
print*,'Initial condition for Eulerian frame for mass center is not allowed!'
stop
!
elseif(nmeth.eq.2)then !...nmeth=2 for Lagrangian frame...
!
if(nint==1)then
!
!...Triangles...
!
do ie = 1, ntria
!
ielem = ie
!
xp(1, 1:nvtri) = coord(1,iptri(1:nvtri, ie))
xp(2, 1:nvtri) = coord(2,iptri(1:nvtri, ie))
!
!...Reference coordinates
!
rc = geoel(1,ielem)
sc = geoel(2,ielem)
volel = geoel(3, ielem)
masel = geoel(4, ielem)
!
!xcel = geoel(7, ielem)
!ycel = geoel(8, ielem)
!
!if(ielem.ge.1.and.ielem.le.20)then
!volsed=volsed + volel
!print*,'ielem',volsed,volel
!endif
!
!...Avoid the centroid out of one curved cell...
!
xcel = (-xp(1, 1)-xp(1, 2)-xp(1, 3)+4.d0*xp(1, 4)+4.d0*xp(1, 5)+4.d0*xp(1, 6))/9.d0
ycel = (-xp(2, 1)-xp(2, 2)-xp(2, 3)+4.d0*xp(2, 4)+4.d0*xp(2, 5)+4.d0*xp(2, 6))/9.d0

!
dr = 0.5d0
ds = 0.5d0
!
rhomvel = 0.d0
rhouvel = 0.d0
rhovvel = 0.d0
rhoevel = 0.d0
!
do igaus =1, ngausd
!
r  = posit(1,igaus)
s  = posit(2,igaus)
wi  = weight(igaus)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
dspr(1) = c10-4.d0*(c10-r-s)
dspr(2) = -1.d0 + 4.d0*r
dspr(3) = 0.d0
dspr(4) = 4.d0*(1.d0-2.d0*r-s)
dspr(5) = 4.d0*s
dspr(6) =-4.d0*s

!
dsps(1) = c10-4.d0*(c10-r-s)
dsps(2) = 0.d0
dsps(3) = -1.d0 + 4.d0*s
dsps(4) = -4.d0*r
dsps(5) =  4.d0*r
dsps(6) =  4.d0*(1.d0-r-2.d0*s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nptri
dxdr = dxdr + dspr(ishp)*xp(1,ishp)
dxds = dxds + dsps(ishp)*xp(1,ishp)

dydr = dydr + dspr(ishp)*xp(2,ishp)
dyds = dyds + dsps(ishp)*xp(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo
!
call getunkin(unkin, xgaus, ygaus, ielem, volel,xcel,ycel)
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
!print*,'bad'
!
if(nrz.eq.0.or.nrz.eq.1)then
rhomvel = rhomvel + unkin(1)*djaco*rcoef
rhouvel = rhouvel + unkin(2)*djaco*rcoef
rhovvel = rhovvel + unkin(3)*djaco*rcoef
rhoevel = rhoevel + unkin(4)*djaco*rcoef
elseif(nrz.eq.2)then
rhomvel = rhomvel + unkin(1)*djaco
rhouvel = rhouvel + unkin(2)*djaco
rhovvel = rhovvel + unkin(3)*djaco
rhoevel = rhoevel + unkin(4)*djaco
endif
!
! if(ie==21.or.ie==5.or.ie==4)print*,'initial ig',ie,ig,uini,-cos(pi*xg)*sin(pi*yg),xg,yg,xc,yc,nvtri,ig,ngausd
!
enddo
!
unkno(1, 1, ielem) = rhomvel/masel;
unkno(1, 2, ielem) = rhouvel/masel;
unkno(1, 3, ielem) = rhovvel/masel;
unkno(1, 4, ielem) = rhoevel/masel
!
enddo
!
!...Second loop for quads...
!
do ie = 1, nquad
!
!...global numbering
!
ielem = ie + ntria
!
xpq(1, 1:nvqua) = coord(1,ipqua(1:nvqua, ie))
xpq(2, 1:nvqua) = coord(2,ipqua(1:nvqua, ie))

!...Get coordinates for FE
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!...Reference coordinates
rc = geoel(1,ielem)
sc = geoel(2,ielem)
volel = geoel(3, ielem)
masel = geoel(4, ielem)
!xcel = geoel(7, ielem)
!ycel = geoel(8, ielem)
!
!if(ielem.ge.1.and.ielem.le.20)then
!  volsed=volsed + volel
!  print*,'ielem',volsed,volel
!endif

!...Physical coordinates
r=rc;s=sc
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xcel = 0.d0
ycel = 0.d0
!
do ishp = 1, nvqua
xcel = xcel + shpq(ishp)*xpq(1,ishp)
ycel = ycel + shpq(ishp)*xpq(2,ishp)
enddo

!
dr = 1.d0  !...Refrence cell (Delt x)
ds = 1.d0
!
rhomvel = 0.d0
rhouvel = 0.d0
rhovvel = 0.d0
rhoevel = 0.d0
!
do igaus =1, ngausdqi
!
r  = posiq(1,igaus)
s  = posiq(2,igaus)
wi  = weighq(igaus)
!
!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
call getunkin(unkin, xgaus, ygaus, ielem, volel,xcel,ycel)
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
!djaco = volel !...assuming degenerate...
!
if(nrz.eq.0.or.nrz.eq.1)then
rhomvel = rhomvel + unkin(1)*djaco*rcoef
rhouvel = rhouvel + unkin(2)*djaco*rcoef
rhovvel = rhovvel + unkin(3)*djaco*rcoef
rhoevel = rhoevel + unkin(4)*djaco*rcoef
elseif(nrz.eq.2)then
rhomvel = rhomvel + unkin(1)*djaco
rhouvel = rhouvel + unkin(2)*djaco
rhovvel = rhovvel + unkin(3)*djaco
rhoevel = rhoevel + unkin(4)*djaco
endif
!
!
! if(ie==21.or.ie==5.or.ie==4)print*,'initial ig',ie,ig,uini,-cos(pi*xg)*sin(pi*yg),xg,yg,xc,yc,nvtri,ig,ngausd
!
enddo
!
unkno(1, 1, ielem) = rhomvel/masel;
unkno(1, 2, ielem) = rhouvel/masel;
unkno(1, 3, ielem) = rhovvel/masel;
unkno(1, 4, ielem) = rhoevel/masel

if(nmatel.eq.1)then
if(ncase.eq.3.and.ielem.gt.1)then
!unkno(1, 2, ielem) =-xcel/sqrt(xcel**2 + ycel**2);
!unkno(1, 3, ielem) =-ycel/sqrt(xcel**2 + ycel**2);
endif
endif
!
!...Record the local basis
!
!geoel(5:6, ielem) = unkno(1, 2:3, ielem)
!if(ielem.eq.100) print*,'variable',ielem,unkno(1, 1:4, ielem),masel,0.5d0*(unkno(1, 2, ielem)**2+unkno(1, 3, ielem)**2),&
!(gamlg-1.d0)*(unkno(1, 4, ielem)-0.5d0*(unkno(1, 2, ielem)**2+unkno(1, 3, ielem)**2))
!
enddo
!
elseif(nint==2)then
!
!Input restarted file...
!
print*,'...................Notice................... !'
print*,'Restart function is not allowed by Lagrangian framework right now !'
!
open(7,file='restartdg.dat')
do ie=1,ncell
read(7,*)unkno(1:ndegr,1,ie)
read(7,*)unkno(1:ndegr,2,ie)
read(7,*)unkno(1:ndegr,3,ie)
read(7,*)unkno(1:ndegr,4,ie)
enddo
close(7)

endif
!
endif
end subroutine inifield_curv
!
!...Get the initial mass matrix for lagrangian curved meshes...
!
subroutine  getamatr_initial_lag_curv(unkno,amatr,geoel,coord,iptri, ipqua)
use constant
implicit none
!...Input
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize)::geoel
real*8,dimension(1:ndimn,1:npoin),intent(in)::coord
real*8,dimension(1:nmatr,1:ncell),intent(out)::amatr
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, ielem
!...Local real array
real*8::xp(1:2, 1:nvtri)
real*8,dimension(1:2, 1:nvqua)::xpq
real*8,dimension(1:nvtri)::shp, dspr, dsps
real*8:: weight(ngausd), posit(2, ngausd)
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
!...Local real
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8::xc, yc
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8::c10
real*8::f0,f1,f2,f3,f4
real*8::f5,f6,f7,f8
real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::b2,b3,b4,b5,b6
real*8::radie, radii,radie2,radii2,radic2, radig2,sentr,rhoin,rhoex,rho0ba
real*8::masel,xgaus,ygaus
real*8::rcoef
real*8::preex,prein
!
real*8,allocatable::x5(:,:),b55(:),mmatr(:,:)
!
data c10 / 1.0d0 /
!
if(npoly==2) allocate(x5(5,5), mmatr(5,5), b55(5))
!
call rutope(ndimn, ngausd, posit, weight)
call ruqope(2, ngausdq, posiq, weighq)
!
!...get amatr...
!...Note: The first term of mass matrix, mass in one cell,
!...is stored in the last term of amatr for convenience...
!
do ie = 1,ntria !...(2)ie = 1,nelem
!
ielem = ie
!
xp(1, 1:nvtri) = coord(1,iptri(1:nvtri, ie))
xp(2, 1:nvtri) = coord(2,iptri(1:nvtri, ie))

!...Mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)
!
dr = 0.5d0
ds = 0.5d0
!
masel = geoel(4, ielem)

!...Physical coordinate for the mass center (rc, sc)
r = rc; s= sc
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nvtri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo

!
f0 = 0.d0
f1 = 0.d0
f2 = 0.d0
f3 = 0.d0
f4 = 0.d0

f22= 0.d0
f23= 0.d0
f24= 0.d0
f25= 0.d0
f26= 0.d0


f33= 0.d0
f34= 0.d0
f35= 0.d0
f36= 0.d0

f44= 0.d0
f45= 0.d0
f46= 0.d0

f55= 0.d0
f56= 0.d0

f66 = 0.d0

do ig =1,ngausd
!
r = posit(1,ig)
s = posit(2,ig)
wi = weight(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
!
dspr(1) = c10-4.d0*(c10-r-s)
dspr(2) = -1.d0 + 4.d0*r
dspr(3) = 0.d0
dspr(4) = 4.d0*(1.d0-2.d0*r-s)
dspr(5) = 4.d0*s
dspr(6) =-4.d0*s

!
dsps(1) = c10-4.d0*(c10-r-s)
dsps(2) = 0.d0
dsps(3) = -1.d0 + 4.d0*s
dsps(4) = -4.d0*r
dsps(5) =  4.d0*r
dsps(6) =  4.d0*(1.d0-r-2.d0*s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*xp(1,ishp)
dxds = dxds + dsps(ishp)*xp(1,ishp)

dydr = dydr + dspr(ishp)*xp(2,ishp)
dyds = dyds + dsps(ishp)*xp(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo
!
xg = r
yg = s

!...Basis function
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
b4 = 0.5d0*b2*b2 - geoel(19, ielem)
b5 = 0.5d0*b3*b3 - geoel(20, ielem)
b6 =       b2*b3 - geoel(21, ielem)
!
call getrhog_initial(rho0, xgaus, ygaus, xc, yc)
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
if(nrz.eq.0.or.nrz.eq.1)then
f0 = f0 + rho0*djaco*rcoef
f1 = f1 + rho0*(xg-rc)/dr*(xg-rc)/dr*djaco*rcoef
f2 = f2 + rho0*(xg-rc)/dr*(yg-sc)/ds*djaco*rcoef
f3 = f3 + rho0*(yg-sc)/ds*(yg-sc)/ds*djaco*rcoef
elseif(nrz.eq.2)then
f0 = f0 + rho0*djaco
f1 = f1 + rho0*(xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + rho0*(xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + rho0*(yg-sc)/ds*(yg-sc)/ds*djaco
endif

if(npoly==2)then
f22 = f22 + rho0*b2*b2*djaco*rcoef
f23 = f23 + rho0*b2*b3*djaco*rcoef
f24 = f24 + rho0*b2*b4*djaco*rcoef
f25 = f25 + rho0*b2*b5*djaco*rcoef
f26 = f26 + rho0*b2*b6*djaco*rcoef

f33 = f33 + rho0*b3*b3*djaco*rcoef
f34 = f34 + rho0*b3*b4*djaco*rcoef
f35 = f35 + rho0*b3*b5*djaco*rcoef
f36 = f36 + rho0*b3*b6*djaco*rcoef

f44 = f44 + rho0*b4*b4*djaco*rcoef
f45 = f45 + rho0*b4*b5*djaco*rcoef
f46 = f46 + rho0*b4*b6*djaco*rcoef

f55 = f55 + rho0*b5*b5*djaco*rcoef
f56 = f56 + rho0*b5*b6*djaco*rcoef

f66 = f66 + rho0*b6*b6*djaco*rcoef
endif
enddo
!
!if(ie.eq.1) print*,'ielem', ie, f22,f23,f24,f25,f26
if(npoly==0)then
amatr(1,ielem) = 1.d0/f0
elseif(npoly==1)then
det = f1*f3-f2**2

amatr(1,ielem) = f3/det
amatr(2,ielem) = -f2/det
amatr(3,ielem) = f1/det
amatr(4,ielem) = 1.d0/f0
elseif(npoly==2)then
mmatr(1,1) = f22
mmatr(1,2) = f23
mmatr(1,3) = f24
mmatr(1,4) = f25
mmatr(1,5) = f26

mmatr(2,1) = mmatr(1,2)
mmatr(2,2) = f33
mmatr(2,3) = f34
mmatr(2,4) = f35
mmatr(2,5) = f36

mmatr(3,1) = mmatr(1,3)
mmatr(3,2) = mmatr(2,3)
mmatr(3,3) = f44
mmatr(3,4) = f45
mmatr(3,5) = f46

mmatr(4,1) = mmatr(1,4)
mmatr(4,2) = mmatr(2,4)
mmatr(4,3) = mmatr(3,4)
mmatr(4,4) = f55
mmatr(4,5) = f56

mmatr(5,1) = mmatr(1,5)
mmatr(5,2) = mmatr(2,5)
mmatr(5,3) = mmatr(3,5)
mmatr(5,4) = mmatr(4,5)
mmatr(5,5) = f66

!...Invert matrix
x5 = 0.d0
b55 = 0.d0
call getinvmat(5, mmatr, x5, b55)
!
amatr(1,ielem) = x5(1,1)
amatr(2,ielem) = x5(1,2)
amatr(3,ielem) = x5(1,3)
amatr(4,ielem) = x5(1,4)
amatr(5,ielem) = x5(1,5)

amatr(6,ielem) = x5(2,2)
amatr(7,ielem) = x5(2,3)
amatr(8,ielem) = x5(2,4)
amatr(9,ielem) = x5(2,5)

amatr(10,ielem) = x5(3,3)
amatr(11,ielem) = x5(3,4)
amatr(12,ielem) = x5(3,5)

amatr(13,ielem) = x5(4,4)
amatr(14,ielem) = x5(4,5)

amatr(15,ielem) = x5(5,5)
!
amatr(16,ielem) = 1.d0/f0
!
endif

!...Treatment for AW RZ
if(nrz.eq.2) amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)
!
enddo !...(2)ie = 1,nelem

!
!...Second loop for quads
!
do ie = 1,nquad !...(2)ie = 1,nelem
!
ielem = ie + ntria
!
xpq(1, 1:nvqua) = coord(1,ipqua(1:nvqua, ie))
xpq(2, 1:nvqua) = coord(2,ipqua(1:nvqua, ie))

!...Get coordinates for FE
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!...Mass center...
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 1.d0
ds = 1.d0!
!
masel = geoel(4, ielem)

!..Physical coordinate
r=rc;s=sc
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nvqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo


f0 = 0.d0
f1 = 0.d0
f2 = 0.d0
f3 = 0.d0
f4 = 0.d0

f22= 0.d0
f23= 0.d0
f24= 0.d0
f25= 0.d0
f26= 0.d0


f33= 0.d0
f34= 0.d0
f35= 0.d0
f36= 0.d0

f44= 0.d0
f45= 0.d0
f46= 0.d0

f55= 0.d0
f56= 0.d0

f66 = 0.d0

do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
b4 = 0.5d0*b2*b2 - geoel(19, ielem)
b5 = 0.5d0*b3*b3 - geoel(20, ielem)
b6 =       b2*b3 - geoel(21, ielem)
!
!...Get the initial density
call getrhog_initial(rho0, xgaus, ygaus, xc, yc)
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
if(nrz.eq.0.or.nrz.eq.1)then
f0 = f0 + rho0*djaco*rcoef
f1 = f1 + rho0*(xg-rc)/dr*(xg-rc)/dr*djaco*rcoef
f2 = f2 + rho0*(xg-rc)/dr*(yg-sc)/ds*djaco*rcoef
f3 = f3 + rho0*(yg-sc)/ds*(yg-sc)/ds*djaco*rcoef
elseif(nrz.eq.2)then
f0 = f0 + rho0*djaco
f1 = f1 + rho0*(xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + rho0*(xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + rho0*(yg-sc)/ds*(yg-sc)/ds*djaco
endif

if(npoly==2)then
f22 = f22 + rho0*b2*b2*djaco*rcoef
f23 = f23 + rho0*b2*b3*djaco*rcoef
f24 = f24 + rho0*b2*b4*djaco*rcoef
f25 = f25 + rho0*b2*b5*djaco*rcoef
f26 = f26 + rho0*b2*b6*djaco*rcoef

f33 = f33 + rho0*b3*b3*djaco*rcoef
f34 = f34 + rho0*b3*b4*djaco*rcoef
f35 = f35 + rho0*b3*b5*djaco*rcoef
f36 = f36 + rho0*b3*b6*djaco*rcoef

f44 = f44 + rho0*b4*b4*djaco*rcoef
f45 = f45 + rho0*b4*b5*djaco*rcoef
f46 = f46 + rho0*b4*b6*djaco*rcoef

f55 = f55 + rho0*b5*b5*djaco*rcoef
f56 = f56 + rho0*b5*b6*djaco*rcoef

f66 = f66 + rho0*b6*b6*djaco*rcoef
endif
!
!if(ie.eq.1) print*,'ielem', ie, f25,rho0*b2*b5*djaco,rho0,b2,b5,djaco,xg,yg
enddo
!
if(npoly==0)then
amatr(1,ielem) = 1.d0/f0
elseif(npoly==1)then
det = f1*f3-f2**2

amatr(1,ielem) = f3/det
amatr(2,ielem) = -f2/det
amatr(3,ielem) = f1/det
amatr(4,ielem) = 1.d0/f0
elseif(npoly==2)then
mmatr(1,1) = f22
mmatr(1,2) = f23
mmatr(1,3) = f24
mmatr(1,4) = f25
mmatr(1,5) = f26

mmatr(2,1) = mmatr(1,2)
mmatr(2,2) = f33
mmatr(2,3) = f34
mmatr(2,4) = f35
mmatr(2,5) = f36

mmatr(3,1) = mmatr(1,3)
mmatr(3,2) = mmatr(2,3)
mmatr(3,3) = f44
mmatr(3,4) = f45
mmatr(3,5) = f46

mmatr(4,1) = mmatr(1,4)
mmatr(4,2) = mmatr(2,4)
mmatr(4,3) = mmatr(3,4)
mmatr(4,4) = f55
mmatr(4,5) = f56

mmatr(5,1) = mmatr(1,5)
mmatr(5,2) = mmatr(2,5)
mmatr(5,3) = mmatr(3,5)
mmatr(5,4) = mmatr(4,5)
mmatr(5,5) = f66
!
!if(ielem.eq.1)print*,'initialmat',mmatr,rc,sc

!...Invert matrix
x5 = 0.d0
b55 = 0.d0
call getinvmat(5, mmatr, x5, b55)
!
amatr(1,ielem) = x5(1,1)
amatr(2,ielem) = x5(1,2)
amatr(3,ielem) = x5(1,3)
amatr(4,ielem) = x5(1,4)
amatr(5,ielem) = x5(1,5)

amatr(6,ielem) = x5(2,2)
amatr(7,ielem) = x5(2,3)
amatr(8,ielem) = x5(2,4)
amatr(9,ielem) = x5(2,5)

amatr(10,ielem) = x5(3,3)
amatr(11,ielem) = x5(3,4)
amatr(12,ielem) = x5(3,5)

amatr(13,ielem) = x5(4,4)
amatr(14,ielem) = x5(4,5)

amatr(15,ielem) = x5(5,5)
!
amatr(16,ielem) = 1.d0/f0
!
endif

!...Treatment for AW RZ
if(nrz.eq.2) amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)
!
enddo !...(2)ie = 1,nelem

end subroutine  getamatr_initial_lag_curv
!
!...Get the initial condition for arbitrary curved mesh
!
subroutine  getunkno_initial_lag_curv(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
use constant
implicit none
!...Input
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
real*8,dimension(1:ngeel,1:nsize)::geoel
real*8,dimension(1:ndimn,1:npoin),intent(in)::coord
real*8,dimension(1:nmatr,1:ncell),intent(out)::amatr
integer,  dimension(1:nvtri,1:ntria), intent(in):: inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, ielem, iq, ideg,idegx,idegy, imx, imy
!...Local real array
real*8::xp(1:2, 1:nvtri)
real*8,dimension(1:2, 1:nvqua)::xpq
real*8,dimension(1:nvtri)::shp, dspr, dsps
real*8:: weight(ngausd), posit(2, ngausd)
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:nq)::unkin
real*8,dimension(1:ndegr, 1:nq)::rhsini
real*8,dimension(1:ndegr)::bt, bq
real*8,dimension(ndegr-1, ndegr-1)::matri, mmatr
real*8,dimension(ndegr-1)::rhsmat

!...Local real
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8::xc, yc
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8::c10
real*8::f0,f1,f2,f3,f4
real*8::f5,f6,f7,f8
real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::b2,b3,b4,b5,b6
real*8::radie, radii,radie2,radii2,radic2, radig2,sentr,rhoin,rhoex,rho0ba
real*8::masel,xgaus,ygaus
real*8::rcoef
real*8::preex,prein
!
data c10 / 1.0d0 /
!
call rutope(ndimn, ngausd, posit, weight)
call ruqope(2, ngausdq, posiq, weighq)
!
!...get amatr...
!...Note: The first term of mass matrix, mass in one cell,
!...is stored in the last term of amatr for convenience...
!
do ie = 1,ntria !...(2)ie = 1,nelem
!
ielem = ie
!
xp(1, 1:nvtri) = coord(1,iptri(1:nvtri, ie))
xp(2, 1:nvtri) = coord(2,iptri(1:nvtri, ie))

!...Mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)
!
dr = 0.5d0
ds = 0.5d0
!
volel = geoel(3, ielem)
masel = geoel(4, ielem)

!...Physical coordinate for the mass center (rc, sc)
r = rc; s= sc
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nvtri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo

!...Zero out the initial rhs
rhsini = 0.d0

!...Zero out local mass matrix
mmatr =.0d0

!...Zero out f0

f0 = 0.d0

do ig =1,ngausd
!
r = posit(1,ig)
s = posit(2,ig)
wi = weight(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
!
dspr(1) = c10-4.d0*(c10-r-s)
dspr(2) = -1.d0 + 4.d0*r
dspr(3) = 0.d0
dspr(4) = 4.d0*(1.d0-2.d0*r-s)
dspr(5) = 4.d0*s
dspr(6) =-4.d0*s

!
dsps(1) = c10-4.d0*(c10-r-s)
dsps(2) = 0.d0
dsps(3) = -1.d0 + 4.d0*s
dsps(4) = -4.d0*r
dsps(5) =  4.d0*r
dsps(6) =  4.d0*(1.d0-r-2.d0*s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*xp(1,ishp)
dxds = dxds + dsps(ishp)*xp(1,ishp)

dydr = dydr + dspr(ishp)*xp(2,ishp)
dyds = dyds + dsps(ishp)*xp(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)

!...Physcial gauss points...
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo

!...Basis function
bt(1) = 1.d0
bt(2) = (r - rc)/dr
bt(3) = (s - sc)/ds

if(npoly.eq.2)then
bt(4) = 0.5d0*bt(2)*bt(2) - geoel(19, ielem)
bt(5) = 0.5d0*bt(3)*bt(3) - geoel(20, ielem)
bt(6) =       bt(2)*bt(3) - geoel(21, ielem)
endif

!...Get the density weighted variables
call getunkin(unkin, xgaus, ygaus, ielem, volel,xc,yc)

!...Get the initial RHS
do iq = 1, nq
do ideg = 1, ndegr
rhsini(ideg, iq) = rhsini(ideg, iq) + unkin(iq)*bt(ideg)*djaco
enddo
enddo

!...Get the initial density
call getrhog_initial(rho0, xgaus, ygaus, xc, yc)

!...Coefficient R of RZ or XY system...
rcoef = 1.d0 - alfrz + alfrz*ygaus

!....P1
if(npoly.eq.1)then

if(nrz.eq.0.or.nrz.eq.1)then
f0 = f0 + rho0*djaco*rcoef
mmatr(1, 1) = mmatr(1, 1) + rho0*bt(2)*bt(2)*djaco*rcoef
mmatr(1, 2) = mmatr(1, 2) + rho0*bt(2)*bt(3)*djaco*rcoef

mmatr(2, 1) = mmatr(2, 1) + rho0*bt(2)*bt(3)*djaco*rcoef
mmatr(2, 2) = mmatr(2, 2) + rho0*bt(3)*bt(3)*djaco*rcoef
!
elseif(nrz.eq.2)then
f0 = f0 + rho0*djaco

mmatr(1, 1) = mmatr(1, 1) + rho0*bt(2)*bt(2)*djaco
mmatr(1, 2) = mmatr(1, 2) + rho0*bt(2)*bt(3)*djaco

mmatr(2, 1) = mmatr(2, 1) + rho0*bt(2)*bt(3)*djaco
mmatr(2, 2) = mmatr(2, 2) + rho0*bt(3)*bt(3)*djaco
endif

!...P2
elseif(npoly==2)then
f0 = f0 + rho0*djaco*rcoef

do imx = 1, ndegr-1
do imy = 1, ndegr-1
mmatr(imx, imy) = mmatr(imx, imy) + rho0*bt(imx+1)*bt(imy+1)*djaco*rcoef
enddo
enddo

endif
enddo
!
!if(ie.eq.1) print*,'ielem', ie, f22,f23,f24,f25,f26

if(npoly==1)then
det = mmatr(1, 1)*mmatr(2, 2) - mmatr(1, 2)**2

amatr(1,ielem) =  mmatr(2, 2)/det
amatr(2,ielem) = -mmatr(1, 2)/det
amatr(3,ielem) =  mmatr(1, 1)/det
amatr(4,ielem) = 1.d0/f0
!
matri(1, 1) = amatr(1,ielem)
matri(1, 2) = amatr(2,ielem)

matri(2, 1) = amatr(2,ielem)
matri(2, 2) = amatr(3,ielem)
elseif(npoly==2)then

!...Invert matrix
matri = 0.d0
rhsmat = 0.d0

call getinvmat(5, mmatr, matri, rhsmat)
!
amatr(1,ielem) = matri(1,1)
amatr(2,ielem) = matri(1,2)
amatr(3,ielem) = matri(1,3)
amatr(4,ielem) = matri(1,4)
amatr(5,ielem) = matri(1,5)

amatr(6,ielem) = matri(2,2)
amatr(7,ielem) = matri(2,3)
amatr(8,ielem) = matri(2,4)
amatr(9,ielem) = matri(2,5)

amatr(10,ielem) = matri(3,3)
amatr(11,ielem) = matri(3,4)
amatr(12,ielem) = matri(3,5)

amatr(13,ielem) = matri(4,4)
amatr(14,ielem) = matri(4,5)

amatr(15,ielem) = matri(5,5)
!
amatr(16,ielem) = 1.d0/f0
!
endif

!...Treatment for AW RZ
if(nrz.eq.2) amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)

!...Zero out the high-order terms
unkno(2:ndegr, :, ielem) = 0.d0

!...Get the initial unkno
do iq = 1, nq
do idegx = 1, ndegr-1
do idegy = 1, ndegr-1
unkno(idegx+1, iq, ielem) = unkno(idegx+1, iq, ielem) + matri(idegx, idegy)*rhsini(idegy+1,iq)
enddo
enddo
enddo
!
enddo !...(2)ie = 1,nelem

!
!...Second loop for quads
!
do ie = 1,nquad !...(2)ie = 1,nelem
!
ielem = ie + ntria
!
xpq(1, 1:nvqua) = coord(1,ipqua(1:nvqua, ie))
xpq(2, 1:nvqua) = coord(2,ipqua(1:nvqua, ie))

!...Get coordinates for FE
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!...Mass center...
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 1.d0
ds = 1.d0!

!..Physical coordinate
r=rc;s=sc
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nvqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo

!
volel = geoel(3, ielem)
masel = geoel(4, ielem)

!...Zero out the initial rhs
rhsini = 0.d0

!...Zero out local mass matrix
mmatr =.0d0

!...Zero out f0

f0 = 0.d0

do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
!...Get the density weighted variables
call getunkin(unkin, xgaus, ygaus, ielem, volel,xc,yc)
!
!...Basis function
bq(1) = 1.d0
bq(2) = (r - rc)/dr
bq(3) = (s - sc)/ds

if(npoly.eq.2)then
bq(4) = 0.5d0*bq(2)*bq(2) - geoel(19, ielem)
bq(5) = 0.5d0*bq(3)*bq(3) - geoel(20, ielem)
bq(6) =       bq(2)*bq(3) - geoel(21, ielem)
endif
!
do iq = 1, nq
do ideg = 1, ndegr
rhsini(ideg, iq) = rhsini(ideg, iq) + unkin(iq)*bq(ideg)*djaco
enddo
enddo

!...Get the initial density
call getrhog_initial(rho0, xgaus, ygaus, xc, yc)


!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
!....P1
if(npoly.eq.1)then

if(nrz.eq.0.or.nrz.eq.1)then
f0 = f0 + rho0*djaco*rcoef
mmatr(1, 1) = mmatr(1, 1) + rho0*bq(2)*bq(2)*djaco*rcoef
mmatr(1, 2) = mmatr(1, 2) + rho0*bq(2)*bq(3)*djaco*rcoef

mmatr(2, 1) = mmatr(2, 1) + rho0*bq(2)*bq(3)*djaco*rcoef
mmatr(2, 2) = mmatr(2, 2) + rho0*bq(3)*bq(3)*djaco*rcoef
!
elseif(nrz.eq.2)then
f0 = f0 + rho0*djaco

mmatr(1, 1) = mmatr(1, 1) + rho0*bq(2)*bq(2)*djaco
mmatr(1, 2) = mmatr(1, 2) + rho0*bq(2)*bq(3)*djaco

mmatr(2, 1) = mmatr(2, 1) + rho0*bq(2)*bq(3)*djaco
mmatr(2, 2) = mmatr(2, 2) + rho0*bq(3)*bq(3)*djaco
endif

!...P2
elseif(npoly==2)then
f0 = f0 + rho0*djaco*rcoef

do imx = 1, ndegr-1
do imy = 1, ndegr-1
mmatr(imx, imy) = mmatr(imx, imy) + rho0*bq(imx+1)*bq(imy+1)*djaco*rcoef
enddo
enddo

endif
enddo
!
!if(ie.eq.1) print*,'ielem', ie, f22,f23,f24,f25,f26

if(npoly==1)then
det = mmatr(1, 1)*mmatr(2, 2) - mmatr(1, 2)**2

amatr(1,ielem) =  mmatr(2, 2)/det
amatr(2,ielem) = -mmatr(1, 2)/det
amatr(3,ielem) =  mmatr(1, 1)/det
amatr(4,ielem) = 1.d0/f0
!
matri(1, 1) = amatr(1,ielem)
matri(1, 2) = amatr(2,ielem)

matri(2, 1) = amatr(2,ielem)
matri(2, 2) = amatr(3,ielem)
elseif(npoly==2)then

!...Invert matrix
matri = 0.d0
rhsmat = 0.d0

call getinvmat(5, mmatr, matri, rhsmat)
!
amatr(1,ielem) = matri(1,1)
amatr(2,ielem) = matri(1,2)
amatr(3,ielem) = matri(1,3)
amatr(4,ielem) = matri(1,4)
amatr(5,ielem) = matri(1,5)

amatr(6,ielem) = matri(2,2)
amatr(7,ielem) = matri(2,3)
amatr(8,ielem) = matri(2,4)
amatr(9,ielem) = matri(2,5)

amatr(10,ielem) = matri(3,3)
amatr(11,ielem) = matri(3,4)
amatr(12,ielem) = matri(3,5)

amatr(13,ielem) = matri(4,4)
amatr(14,ielem) = matri(4,5)

amatr(15,ielem) = matri(5,5)
!
amatr(16,ielem) = 1.d0/f0
!
endif

!...Treatment for AW RZ
if(nrz.eq.2) amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)

!...Zero out the high-order terms
unkno(2:ndegr, :, ielem) = 0.d0
!
if(ielem.eq.1)print*,'ielem',matri

!...Get the initial unkno
do iq = 1, nq
do idegx = 1, ndegr-1
do idegy = 1, ndegr-1
unkno(idegx+1, iq, ielem) = unkno(idegx+1, iq, ielem) + matri(idegx, idegy)*rhsini(idegy+1,iq)
enddo
enddo
enddo
!
enddo !...(2)ie = 1,nelem

end subroutine  getunkno_initial_lag_curv
!
!
!...Update the parameters for volumes for Lag...
!
subroutine getcell_evol_lag(iptri, ipqua, geoel, coord)
use constant
implicit none
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:ngeel, 1:nsize), intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
!
!...local array...
real*8,dimension(1:2, 1:nvtri)::xp
real*8,dimension(1:2, 1:nvqua)::xpq
real*8,dimension(1:nvtri)::shp, dspr, dsps
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8:: weight(ngausd_geo), posit(2, ngausd_geo)
real*8:: weighq(ngausd_geoq), posiq(2, ngausd_geoq)
!...local real number
real*8:: dxdr,dxds,dydr,dyds
real*8:: eps,c00,c10,c05,c20
real*8:: r, s, djaco, volel, masel,wi
real*8:: rcv, scv
real*8:: rcvrz,scvrz
real*8:: xgaus, ygaus
real*8:: rcoef
real*8:: areel
!
integer::ielem, igaus, ishp
integer::ie
!
data eps / 1.0d-06 /
data c00 / 0.0d0 /
data c10 / 1.0d0 /
data c05 / 0.5d0 /
data c20 / 2.0d0 /
!
!...Find weight and position for gauss points...
call rutope(2, ngausd_geo, posit, weight)
call ruqope(2, ngausd_geoq, posiq, weighq)
!
!...1st loop to find center and volulme...
!
do ie  = 1, ntria !...(1)ie = 1,nelem
!
ielem = ie
!
xp(1, 1:nvtri) = coord(1,iptri(1:nvtri, ie))
xp(2, 1:nvtri) = coord(2,iptri(1:nvtri, ie))

!...Zero out some variables

rcv = 0.d0
scv = 0.d0
volel = 0.d0
areel = 0.d0
rcvrz = 0.d0
scvrz = 0.d0
!
do igaus =1,ngausd_geo
!
r  = posit(1,igaus)
s  = posit(2,igaus)
wi  = weight(igaus)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
!
dspr(1) = c10-4.d0*(c10-r-s)
dspr(2) = -1.d0 + 4.d0*r
dspr(3) = 0.d0
dspr(4) = 4.d0*(1.d0-2.d0*r-s)
dspr(5) = 4.d0*s
dspr(6) =-4.d0*s

!
dsps(1) = c10-4.d0*(c10-r-s)
dsps(2) = 0.d0
dsps(3) = -1.d0 + 4.d0*s
dsps(4) = -4.d0*r
dsps(5) =  4.d0*r
dsps(6) =  4.d0*(1.d0-r-2.d0*s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*xp(1,ishp)
dxds = dxds + dsps(ishp)*xp(1,ishp)

dydr = dydr + dspr(ishp)*xp(2,ishp)
dyds = dyds + dsps(ishp)*xp(2,ishp)
enddo
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
volel = volel + djaco*rcoef
rcv = rcv + r*djaco*rcoef
scv = scv + s*djaco*rcoef
!
if(nrz.eq.2)then
areel = areel + djaco !...The area in R-Z coordinates
rcvrz = rcvrz + r*djaco
scvrz = scvrz + s*djaco
endif
!
enddo
!
geoel(3, ielem) = volel
geoel(5, ielem) = rcv/volel
geoel(6, ielem) = scv/volel
!geoel(5, ielem) = geoel(1, ielem)
!geoel(6, ielem) = geoel(2, ielem)
!
if(nrz.eq.2)then
geoel(5, ielem) = rcvrz/areel
geoel(6, ielem) = scvrz/areel
geoel(11, ielem) = volel/areel !...Cell averaged radius
geoel(12, ielem) = areel
geoel(15, ielem) = rcvrz/areel
geoel(16, ielem) = scvrz/areel
endif
!
enddo !...(1)ie = 1,nelem
!
!...2nd loop to find  mass center and volulme for quads...
!
do ie  = 1, nquad !...(1)ie = 1,nelem
!
!...global numbering
!
ielem = ie + ntria
!
xpq(1, 1:nvqua) = coord(1,ipqua(1:nvqua, ie))
xpq(2, 1:nvqua) = coord(2,ipqua(1:nvqua, ie))

!...Get coordinates for FE
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!...Zero out some variables...
rcv = 0.d0
scv = 0.d0
volel = 0.d0
areel = 0.d0
rcvrz = 0.d0
scvrz = 0.d0
!
do igaus =1,ngausd_geoq
!
r  = posiq(1,igaus)
s  = posiq(2,igaus)
wi  = weighq(igaus)

!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
djaco = wi*(dxdr*dyds - dydr*dxds)
volel = volel + djaco*rcoef
rcv = rcv + r*djaco*rcoef
scv = scv + s*djaco*rcoef
!
if(nrz.eq.2)then
areel = areel + djaco !...The area in R-Z coordinates
rcvrz = rcvrz + r*djaco
scvrz = scvrz + s*djaco
endif
!
enddo
!
geoel(3, ielem) = volel
geoel(5, ielem) = rcv/volel
geoel(6, ielem) = scv/volel
!geoel(5, ielem) = geoel(1, ielem)
!geoel(6, ielem) = geoel(2, ielem)
!
if(nrz.eq.2)then
geoel(5, ielem) = rcvrz/areel
geoel(6, ielem) = scvrz/areel
geoel(11, ielem) = volel/areel !...Cell averaged radius
geoel(12, ielem) = areel
geoel(15, ielem) = rcvrz/areel !...Mass center with the averaged density <===> volume center
geoel(16, ielem) = scvrz/areel
endif
!
enddo !...(1)ie = 1,nelem

end subroutine getcell_evol_lag
!
!...Source term integration for curv quad grids...
!
subroutine rhsdomnsrcdg_lag_quadc(intfac, ipqua, coord, geoel,rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:ncell),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem
!
!...local integer array
!
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)

!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8::rm,sm,rp,sp
real*8::dr,ds,rc,sc
real*8::xg, yg
real*8::djaco, wi,src

!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
call ruqope(2, ngausdq, posiq, weighq)
!
!...Loop over quads
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
xpq(1, 1:nvqua) = coord(1,ipqua(1:nvqua, ie))
xpq(2, 1:nvqua) = coord(2,ipqua(1:nvqua, ie))

!...Get the points for FEM
call getcoord_fe(ncurv, nvfac, nvqua, xpq)

!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!rxx=
!rxy=
!ryy=
!
!...Gauss loop
!
do ig = 1,ngausdq !...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
wi  = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
!...Gauss points...
!
xg =0.d0; yg= 0.d0
!
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpq(1,ishp)
yg = yg + shpq(ishp)*xpq(2,ishp)
enddo
!
!
!...Basis function for solutions...
!
b(1) = 1.d0

if(npoly.ge.1)then
 b(2) = (r-rc)/dr
 b(3) = (s-sc)/ds

!DGP2
 if(npoly.eq.2)then
 b(4) = 0.5d0*b(2)*b(2) - geoel(19, ielem)
 b(5) = 0.5d0*b(3)*b(3) - geoel(20, ielem)
 b(6) =       b(2)*b(3) - geoel(21, ielem)
 endif
endif

!...source term
src = 0.25d0*pi*(cos(3.d0*pi*xg)*cos(pi*yg) - cos(3.d0*pi*yg)*cos(pi*xg))/(gamlg-1.d0)
!
do ideg = 1,ndegr
rhsel(ideg,4,ielem)=rhsel(ideg,4,ielem) + src*b(ideg)*djaco
enddo
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
end subroutine rhsdomnsrcdg_lag_quadc
!
!...Get the averaged-density within one sub-cell using SMS ...
!
subroutine  getrhosubcell_sms_general(xcrho, ycrho, xpq,xpqi,unksgq, ielem)
use constant
implicit none
!...Input
integer, intent(in)::ielem
real*8,  intent(in)::xcrho,ycrho
real*8,dimension(1:2, 1:nvqua), intent(in)::xpq
real*8,dimension(1:2, 1:nvqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:nsubq), intent(out)::unksgq
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc

integer,dimension(4, nsubq)::ipqsg
!...Local real array
real*8,dimension(1:2, 1:nvqua)::xpq_sg,xpqi_sg
real*8,dimension(1:ndegr,1:nsubq)::rhsel
real*8,dimension(1:nmatr, 1:nsubq)::matin
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:4)    ::shpql, dsprql, dspsql
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: xrqc(2, 16),xrq(2, nvqua)
real*8:: posif(nvfac)
real*8::bq(ndegr)
!...Local real
real*8::r, s, rc,sc,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8::c10
real*8::f0, wt
real*8::masel,xgaus,ygaus

!
real*8,allocatable::mminv(:,:),binv(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
!...Part 0: Preliminary setup...
!
if(npoly==2) allocate(mminv(5,5), mmatr(5,5), binv(5))
!
if(ncurv.eq.1)then
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4

elseif(ncurv.eq.2)then
!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8
!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;

!posif(2) = -1.d0/3.d0; posif(3) = 1.d0/3.d0;
!...Vertex coordinate
xrq(1, 1) = -1.d0;    xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0;    xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0;    xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0;    xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;    xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;    xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;    xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;    xrq(2,12) = posif(2);
!
xrqc(:, 1:12) = xrq(:, 1:12)
xrqc(1, 13) = posif(2); xrqc(2, 13) = posif(2);
xrqc(1, 14) = posif(3); xrqc(2, 14) = posif(2);
xrqc(1, 15) = posif(3); xrqc(2, 15) = posif(3);
xrqc(1, 16) = posif(2); xrqc(2, 16) = posif(3);

endif

!...Get the coordinates for FE
xpq_sg  = xpq
xpqi_sg = xpqi

!...Get the Gauss
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0
!
!...Part I: Get cell averaged r and s for one subgrid...
!
do isg = 1, nsubq

!...Initialze parameters...
if(npoly.eq.2)  mmatr = 0.d0
f0 = 0.d0
wt = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!... update the reference coordinates (r, s)
call getshapfct_quad(0, 4, shpql, dsprql, dspsql, r, s)
!
r = 0.d0
s = 0.d0
do ishp = 1, 4
r =  r + shpql(ishp)*xrqc(1, ipqsg(ishp, isg))
s =  s + shpql(ishp)*xrqc(2, ipqsg(ishp, isg))
enddo

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq_sg(1,ishp)
dxds = dxds + dspsq(ishp)*xpq_sg(1,ishp)

dydr = dydr + dsprq(ishp)*xpq_sg(2,ishp)
dyds = dyds + dspsq(ishp)*xpq_sg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco
wt = wt +wi
!
!if(ielem.eq.1) print*,'djaco',isg,ig, wi,dxdr*dyds - dydr*dxds,xpq_sg
enddo
!if(ielem.eq.1) print*,'volume',f0,isg,wt
!
!stop
!
if(npoly==1)then
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then
matin(16,isg) = 1.d0/f0
endif

enddo
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, nsubq
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!... update the reference coordinates (r, s)
call getshapfct_quad(0, 4, shpql, dsprql, dspsql, r, s)
!
r = 0.d0
s = 0.d0
do ishp = 1, 4
r =  r + shpql(ishp)*xrqc(1, ipqsg(ishp, isg))
s =  s + shpql(ishp)*xrqc(2, ipqsg(ishp, isg))
enddo

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpqi_sg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqi_sg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqi_sg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqi_sg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
bq(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpqi_sg(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqi_sg(2,ishp)
enddo

!...Get the initial density
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*bq(1)*djaco
enddo
!
if(npoly==1)then
m(1,1) = matin(4, isg)
elseif(npoly==2)then
m(1,1) = matin(16, isg)
endif

!...Update the subcell density ditribution
unint = 0.d0
unint(1) = unint(1) + m(1, 1)*rhsel(1,isg)

!if(ielem.eq.1) print*,'isg4',isg,m(1, 1),rhsel(1,isg),m(1, 1)*rhsel(1,isg)
unksgq(1, isg)= unint(1)
enddo
!
end subroutine  getrhosubcell_sms_general
!
!...Get the averaged density at the sub cell using the whole cell distribution...
!
subroutine  getrhosubcell_daverg(rc,sc, bqho, xpq, unksgq, unkno, aflim, uqsgc,ielem)
use constant
implicit none
!...Input
real*8, intent(in)::rc, sc
integer, intent(in)::ielem
real*8,dimension(3),            intent(in)::bqho
real*8,dimension(1:2, 1:nvqua), intent(in)::xpq
real*8,dimension(1:ndegr, 1:nsubq), intent(in)::unksgq
real*8,dimension(1:ndegr,1:nq),      intent(in)::unkno
real*8,dimension(1:nq+1),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:nq+1, 1:nsubq),       intent(out)::uqsgc
!
!...Local integer
integer :: ie, ig, ishp,id, isg, ideg
integer,dimension(4, nsubq)::ipqsg
!...Local real array
real*8,dimension(1:2, 1:nvqua)::xpq_sg
real*8,dimension(1:ndegr)::bqsg
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:nvqua)::shpql, dsprql, dspsql
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: xrqc(2, 16),xrq(2, nvqua)
real*8:: posif(nvfac)
!...Local real
real*8::r, s, dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhomsg,rhomc
real*8::wi,djaco,volel, masel
real*8::c10
!
!-xxx-real contant
data c10 / 1.0d0 /
!
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0
!
!...Now configuration...
xpq_sg = xpq
!
if(ncurv.eq.1)then
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4

elseif(ncurv.eq.2)then
!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8

!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;

!posif(2) = -1.d0/3.d0; posif(3) = 1.d0/3.d0;
!...Vertex coordinate
xrq(1, 1) = -1.d0;    xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0;    xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0;    xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0;    xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;    xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;    xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;    xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;    xrq(2,12) = posif(2);
!
xrqc(:, 1:12) = xrq(:, 1:12)
xrqc(1, 13) = posif(2); xrqc(2, 13) = posif(2);
xrqc(1, 14) = posif(3); xrqc(2, 14) = posif(2);
xrqc(1, 15) = posif(3); xrqc(2, 15) = posif(3);
xrqc(1, 16) = posif(2); xrqc(2, 16) = posif(3);
endif
!
rhomc = unkno(1,1)
!
!...Part 0: Get cell averaged r and s...
!
do isg = 1, nsubq

!...Initialze parameters...
masel = 0.d0
volel = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1, ig)
s = posiq(2, ig)
wi = weighq(ig)

!... update the reference coordinates (r, s)
call getshapfct_quad(0, 4, shpql, dsprql, dspsql, r, s)
!
r = 0.d0
s = 0.d0
do ishp = 1, 4
r =  r + shpql(ishp)*xrqc(1, ipqsg(ishp, isg))
s =  s + shpql(ishp)*xrqc(2, ipqsg(ishp, isg))
enddo

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!Basis function
bqsg(1) = 1.d0

if(npoly.ge.1)then
 bqsg(2) = (r-rc)/dr
 bqsg(3) = (s-sc)/ds

!...DGP2
 if(npoly.eq.2)then
  bqsg(4) = 0.5d0*bqsg(2)*bqsg(2) - bqho(1)
  bqsg(5) = 0.5d0*bqsg(3)*bqsg(3) - bqho(2)
  bqsg(6) =       bqsg(2)*bqsg(3) - bqho(3)
 endif
endif
!
rhomsg = 0.d0
do ideg = 1, mdegr
rhomsg = rhomsg + unkno(ideg,1)*bqsg(ideg)
enddo
!
!if(ielem.eq.1)print*,'drho',ielem,isg,ig,rhomsg
!
if(nlimi.eq.6)then
rhomsg = rhomc + aflim(1)*(rhomsg - rhomc)
endif

!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq_sg(1,ishp)
dxds = dxds + dspsq(ishp)*xpq_sg(1,ishp)

dydr = dydr + dsprq(ishp)*xpq_sg(2,ishp)
dyds = dyds + dspsq(ishp)*xpq_sg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)*0.25d0
!
!...Parameters for now configuration...
!
masel = masel + 1.d0/rhomsg*djaco
volel = volel + djaco
!
enddo

!...The avearged specific volume at the sub-cell...
uqsgc(1, isg) = volel/masel
!
!if(ielem.eq.1)print*,'drho',ielem,isg,uqsgc(1, isg),unkno(:,1)
!
enddo
!
end subroutine getrhosubcell_daverg
!
!...Get the averaged-density within one sub-cell using SMS ...
!
subroutine  getrhosubcell_sms_general2(xcrho, ycrho, xpq,xpqi,unksgq, ielem)
use constant
implicit none
!...Input
integer, intent(in)::ielem
real*8,  intent(in)::xcrho,ycrho
real*8,dimension(1:2, 1:nvqua), intent(in)::xpq
real*8,dimension(1:2, 1:nvqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:nsubq), intent(out)::unksgq
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc,ip,ifa

integer,dimension(4, nsubq)::ipqsg
integer, dimension(1:4,1:4)::ipfsg
!...Local real array
real*8,dimension(1:2, 1:nvqua)::xpq_sg,xpqi_sg,xpq_fe,xpqi_fe
real*8,dimension(1:ndegr,1:nsubq)::rhsel
real*8,dimension(1:nmatr, 1:nsubq)::matin
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:4)    ::shpql, dsprql, dspsql
real*8,dimension(1:2)    ::shpfl, dsprfl
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: xrqc(2, 16),xrq(2, nvqua),xrq_sg(2, nvqua)
real*8:: posif(nvfac)
real*8::bq(ndegr)
!...Local real
real*8::r, s, rc,sc,dr,ds,xg,yg
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8::c10
real*8::f0, wt
real*8::masel,xgaus,ygaus

!
real*8,allocatable::mminv(:,:),binv(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
!...Part 0: Preliminary setup...
!
if(npoly==2) allocate(mminv(5,5), mmatr(5,5), binv(5))
!
if(ncurv.eq.1)then
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4

elseif(ncurv.eq.2)then
!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8

!...face for subcell
ipfsg(1, 1) = 1; ipfsg(2, 1) = 2; ipfsg(3, 1) = 5; ipfsg(4, 1) = 9; 
ipfsg(1, 2) = 2; ipfsg(2, 2) = 3; ipfsg(3, 2) = 6; ipfsg(4, 2) =10;  
ipfsg(1, 3) = 3; ipfsg(2, 3) = 4; ipfsg(3, 3) = 7; ipfsg(4, 3) =11; 
ipfsg(1, 4) = 4; ipfsg(2, 4) = 1; ipfsg(3, 4) = 8; ipfsg(4, 4) =12;  

!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;

!posif(2) = -1.d0/3.d0; posif(3) = 1.d0/3.d0;
!...Vertex coordinate
xrq(1, 1) = -1.d0;    xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0;    xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0;    xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0;    xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;    xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;    xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;    xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;    xrq(2,12) = posif(2);
!
xrqc(:, 1:12) = xrq(:, 1:12)
xrqc(1, 13) = posif(2); xrqc(2, 13) = posif(2);
xrqc(1, 14) = posif(3); xrqc(2, 14) = posif(2);
xrqc(1, 15) = posif(3); xrqc(2, 15) = posif(3);
xrqc(1, 16) = posif(2); xrqc(2, 16) = posif(3);

endif

!...Get the coordinates for FE
xpq_fe  = xpq
xpqi_fe = xpqi

!
call getcoord_fe(ncurv, nvfac, nvqua, xpq_fe)
call getcoord_fe(ncurv, nvfac, nvqua, xpqi_fe)

!...Get the Gauss
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0
!
!...Part I: Get cell averaged r and s for one subgrid...
!
do isg = 1, nsubq

!...Reference coordinates
!
do ifa = 1, 4
   r = -1.d0/3.d0
   call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
   xrq_sg(1:2, ipfsg(3, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1, ifa), isg)) + &
                                shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
!
   r =  1.d0/3.d0
   call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
   xrq_sg(1:2, ipfsg(4, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1, ifa), isg)) + &
                                shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
enddo

!...Physical coordinates
!
do ip = 1,4
 r = xrqc(1, ipqsg(ip, isg))
 s = xrqc(2, ipqsg(ip, isg))
!
 call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
 xg = 0.d0
 yg = 0.d0
do ishp = 1, nvqua
 xg = xg + shpq(ishp)*xpq_fe(1, ishp)
 yg = yg + shpq(ishp)*xpq_fe(2, ishp)
enddo
!
 xpq_sg(1, ip) = xg
 xpq_sg(2, ip) = yg
enddo
!
do ip = 5, 12
 r = xrq_sg(1, ip)
 s = xrq_sg(2, ip)
!
 call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
 xg = 0.d0
 yg = 0.d0
do ishp = 1, nvqua
 xg = xg + shpq(ishp)*xpq_fe(1, ishp)
 yg = yg + shpq(ishp)*xpq_fe(2, ishp)
enddo
!
 xpq_sg(1, ip) = xg
 xpq_sg(2, ip) = yg
!
enddo

!...Initialze parameters...
if(npoly.eq.2)  mmatr = 0.d0
f0 = 0.d0
wt = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq_sg(1,ishp)
dxds = dxds + dspsq(ishp)*xpq_sg(1,ishp)

dydr = dydr + dsprq(ishp)*xpq_sg(2,ishp)
dyds = dyds + dspsq(ishp)*xpq_sg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco
wt = wt +wi
!
!if(ielem.eq.1) print*,'djaco',isg,ig, wi,dxdr*dyds - dydr*dxds,xpq_sg
enddo
!if(ielem.eq.1) print*,'volume',f0,isg,wt
!
!stop
!
if(npoly==1)then
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then
matin(16,isg) = 1.d0/f0
endif

enddo
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, nsubq
!...Reference coordinates
!
do ifa = 1, 4
   r = -1.d0/3.d0
   call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
   xrq_sg(1:2, ipfsg(3, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1, ifa), isg)) + &
                                shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
!
   r =  1.d0/3.d0
   call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
   xrq_sg(1:2, ipfsg(4, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1, ifa), isg)) + &
                                shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
enddo

!...Physical coordinates
do ip = 1,4
 r = xrqc(1, ipqsg(ip, isg))
 s = xrqc(2, ipqsg(ip, isg))
!
 call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
 xg = 0.d0
 yg = 0.d0
do ishp = 1, nvqua
 xg = xg + shpq(ishp)*xpqi_fe(1, ishp)
 yg = yg + shpq(ishp)*xpqi_fe(2, ishp)
enddo
!
 xpqi_sg(1, ip) = xg
 xpqi_sg(2, ip) = yg
enddo

do ip = 5, 12
 r = xrq_sg(1, ip)
 s = xrq_sg(2, ip)
!
 call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
 xg = 0.d0
 yg = 0.d0
do ishp = 1, nvqua
 xg = xg + shpq(ishp)*xpqi_fe(1, ishp)
 yg = yg + shpq(ishp)*xpqi_fe(2, ishp)
enddo
!
 xpqi_sg(1, ip) = xg
 xpqi_sg(2, ip) = yg
!
enddo
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpqi_sg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqi_sg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqi_sg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqi_sg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
bq(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpqi_sg(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqi_sg(2,ishp)
enddo

!...Get the initial density
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*bq(1)*djaco
enddo
!
if(npoly==1)then
m(1,1) = matin(4, isg)
elseif(npoly==2)then
m(1,1) = matin(16, isg)
endif

!...Update the subcell density ditribution
unint = 0.d0
unint(1) = unint(1) + m(1, 1)*rhsel(1,isg)

!if(ielem.eq.1) print*,'isg4',isg,m(1, 1),rhsel(1,isg),m(1, 1)*rhsel(1,isg)
unksgq(1, isg)= unint(1)
enddo
!
end subroutine  getrhosubcell_sms_general2
!
!...Get the averaged-density within one sub-cell using SMS ...
!
subroutine  getrhosubcell_sms_general3(xcrho, ycrho, xpq_fe,xpqi_fe,unksgq, ielem)
use constant
implicit none
!...Input
integer, intent(in)::ielem
real*8,  intent(in)::xcrho,ycrho
real*8,dimension(1:2, 1:nvqua), intent(in)::xpq_fe
real*8,dimension(1:2, 1:nvqua), intent(in)::xpqi_fe
real*8,dimension(1:ndegr, 1:nsubq), intent(out)::unksgq
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc,ip,ifa

integer,dimension(4, nsubq)::ipqsg
integer, dimension(1:4,1:4)::ipfsg
!...Local real array
real*8,dimension(1:2, 1:nvqua)::xpq_sg,xpqi_sg
real*8,dimension(1:ndegr,1:nsubq)::rhsel
real*8,dimension(1:nmatr, 1:nsubq)::matin
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:4)    ::shpql, dsprql, dspsql
real*8,dimension(1:2)    ::shpfl, dsprfl
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: xrqc(2, 16),xrq(2, nvqua),xrq_sg(2, nvqua)
real*8:: posif(nvfac)
real*8::bq(ndegr)
!...Local real
real*8::r, s, rc,sc,dr,ds,xg,yg
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8::c10
real*8::f0, wt
real*8::masel,xgaus,ygaus

!
real*8,allocatable::mminv(:,:),binv(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
!...Part 0: Preliminary setup...
!
if(npoly==2) allocate(mminv(5,5), mmatr(5,5), binv(5))
!
if(ncurv.eq.1)then
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4

elseif(ncurv.eq.2)then
!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8

!...face for subcell
ipfsg(1, 1) = 1; ipfsg(2, 1) = 2; ipfsg(3, 1) = 5; ipfsg(4, 1) = 9; 
ipfsg(1, 2) = 2; ipfsg(2, 2) = 3; ipfsg(3, 2) = 6; ipfsg(4, 2) =10;  
ipfsg(1, 3) = 3; ipfsg(2, 3) = 4; ipfsg(3, 3) = 7; ipfsg(4, 3) =11; 
ipfsg(1, 4) = 4; ipfsg(2, 4) = 1; ipfsg(3, 4) = 8; ipfsg(4, 4) =12;  

!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;

!posif(2) = -1.d0/3.d0; posif(3) = 1.d0/3.d0;
!...Vertex coordinate
xrq(1, 1) = -1.d0;    xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0;    xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0;    xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0;    xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;    xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;    xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;    xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;    xrq(2,12) = posif(2);
!
xrqc(:, 1:12) = xrq(:, 1:12)
xrqc(1, 13) = posif(2); xrqc(2, 13) = posif(2);
xrqc(1, 14) = posif(3); xrqc(2, 14) = posif(2);
xrqc(1, 15) = posif(3); xrqc(2, 15) = posif(3);
xrqc(1, 16) = posif(2); xrqc(2, 16) = posif(3);

endif

!...Get the Gauss
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0
!
!...Part I: Get cell averaged r and s for one subgrid...
!
do isg = 1, nsubq

!...Reference coordinates
!
do ifa = 1, 4
   r = -1.d0/3.d0
   call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
   xrq_sg(1:2, ipfsg(3, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1, ifa), isg)) + &
                                shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
!
   r =  1.d0/3.d0
   call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
   xrq_sg(1:2, ipfsg(4, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1, ifa), isg)) + &
                                shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
enddo

!...Physical coordinates
!
do ip = 1,4
 r = xrqc(1, ipqsg(ip, isg))
 s = xrqc(2, ipqsg(ip, isg))
!
 call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
 xg = 0.d0
 yg = 0.d0
do ishp = 1, nvqua
 xg = xg + shpq(ishp)*xpq_fe(1, ishp)
 yg = yg + shpq(ishp)*xpq_fe(2, ishp)
enddo
!
 xpq_sg(1, ip) = xg
 xpq_sg(2, ip) = yg
enddo
!
do ip = 5, 12
 r = xrq_sg(1, ip)
 s = xrq_sg(2, ip)
!
 call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
 xg = 0.d0
 yg = 0.d0
do ishp = 1, nvqua
 xg = xg + shpq(ishp)*xpq_fe(1, ishp)
 yg = yg + shpq(ishp)*xpq_fe(2, ishp)
enddo
!
 xpq_sg(1, ip) = xg
 xpq_sg(2, ip) = yg
!
enddo

!...Initialze parameters...
if(npoly.eq.2)  mmatr = 0.d0
f0 = 0.d0
wt = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq_sg(1,ishp)
dxds = dxds + dspsq(ishp)*xpq_sg(1,ishp)

dydr = dydr + dsprq(ishp)*xpq_sg(2,ishp)
dyds = dyds + dspsq(ishp)*xpq_sg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco
wt = wt +wi
!
!if(ielem.eq.1) print*,'djaco',isg,ig, wi,dxdr*dyds - dydr*dxds,xpq_sg
enddo
!if(ielem.eq.1) print*,'volume',f0,isg,wt
!
!stop
!
if(npoly==1)then
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then
matin(16,isg) = 1.d0/f0
endif

enddo
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, nsubq
!...Reference coordinates
!
do ifa = 1, 4
   r = -1.d0/3.d0
   call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
   xrq_sg(1:2, ipfsg(3, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1, ifa), isg)) + &
                                shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
!
   r =  1.d0/3.d0
   call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
   xrq_sg(1:2, ipfsg(4, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1, ifa), isg)) + &
                                shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
enddo

!...Physical coordinates
do ip = 1,4
 r = xrqc(1, ipqsg(ip, isg))
 s = xrqc(2, ipqsg(ip, isg))
!
 call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
 xg = 0.d0
 yg = 0.d0
do ishp = 1, nvqua
 xg = xg + shpq(ishp)*xpqi_fe(1, ishp)
 yg = yg + shpq(ishp)*xpqi_fe(2, ishp)
enddo
!
 xpqi_sg(1, ip) = xg
 xpqi_sg(2, ip) = yg
enddo

do ip = 5, 12
 r = xrq_sg(1, ip)
 s = xrq_sg(2, ip)
!
 call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
 xg = 0.d0
 yg = 0.d0
do ishp = 1, nvqua
 xg = xg + shpq(ishp)*xpqi_fe(1, ishp)
 yg = yg + shpq(ishp)*xpqi_fe(2, ishp)
enddo
!
 xpqi_sg(1, ip) = xg
 xpqi_sg(2, ip) = yg
!
enddo
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpqi_sg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqi_sg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqi_sg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqi_sg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
bq(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpqi_sg(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqi_sg(2,ishp)
enddo

!...Get the initial density
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*bq(1)*djaco
enddo
!
if(npoly==1)then
m(1,1) = matin(4, isg)
elseif(npoly==2)then
m(1,1) = matin(16, isg)
endif

!...Update the subcell density ditribution
unint = 0.d0
unint(1) = unint(1) + m(1, 1)*rhsel(1,isg)

!if(ielem.eq.1) print*,'isg4',isg,m(1, 1),rhsel(1,isg),m(1, 1)*rhsel(1,isg)
unksgq(1, isg)= unint(1)
enddo
!
end subroutine  getrhosubcell_sms_general3
!
!...Get the averaged-density within one sub-cell using SMS ...
!
subroutine  getrhosubcell_sms_general4(xcrho, ycrho,xpq_fe,xpqi_fe,unksgq, ielem)
use constant
implicit none
!...Input
integer, intent(in)::ielem
real*8,  intent(in)::xcrho,ycrho
real*8,dimension(1:2, 1:nvqua), intent(in)::xpq_fe
real*8,dimension(1:2, 1:nvqua), intent(in)::xpqi_fe
real*8,dimension(1:ndegr, 1:nsubq), intent(out)::unksgq
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc,ip,ifa

integer,dimension(4, nsubq)::ipqsg
integer, dimension(1:4,1:4)::ipfsg
!...Local real array
real*8,dimension(1:2, 1:nvqua)::xpq_sg,xpqi_sg
real*8,dimension(1:ndegr,1:nsubq)::rhsel
real*8,dimension(1:nmatr, 1:nsubq)::matin
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:4)    ::shpql, dsprql, dspsql
real*8,dimension(1:2)    ::shpfl, dsprfl
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: xrqc(2, 16),xrq(2, nvqua),xrq_sg(2, nvqua)
real*8:: posif(nvfac)
real*8::bq(ndegr)
!...Local real
real*8::r, s, rc,sc,dr,ds,xg,yg
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8::c10
real*8::f0, wt
real*8::masel,xgaus,ygaus

!
real*8,allocatable::mminv(:,:),binv(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
!...Part 0: Preliminary setup...
!
if(npoly==2) allocate(mminv(5,5), mmatr(5,5), binv(5))
!
if(ncurv.eq.1)then
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4

elseif(ncurv.eq.2)then
!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8

!...face for subcell
ipfsg(1, 1) = 1; ipfsg(2, 1) = 2; ipfsg(3, 1) = 5; ipfsg(4, 1) = 9;
ipfsg(1, 2) = 2; ipfsg(2, 2) = 3; ipfsg(3, 2) = 6; ipfsg(4, 2) =10;
ipfsg(1, 3) = 3; ipfsg(2, 3) = 4; ipfsg(3, 3) = 7; ipfsg(4, 3) =11;
ipfsg(1, 4) = 4; ipfsg(2, 4) = 1; ipfsg(3, 4) = 8; ipfsg(4, 4) =12;

!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;
!...Vertex coordinate
xrq(1, 1) = -1.d0;    xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0;    xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0;    xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0;    xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;    xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;    xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;    xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;    xrq(2,12) = posif(2);
!
xrqc(:, 1:12) = xrq(:, 1:12)
xrqc(1, 13) = posif(2); xrqc(2, 13) = posif(2);
xrqc(1, 14) = posif(3); xrqc(2, 14) = posif(2);
xrqc(1, 15) = posif(3); xrqc(2, 15) = posif(3);
xrqc(1, 16) = posif(2); xrqc(2, 16) = posif(3);

endif

!...Get the Gauss
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0
!
!...Part I: Get cell averaged r and s for one subgrid...
!
do isg = 1, nsubq

!...Reference coordinates
!
do ifa = 1, 4
r =  0.d0
call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
xrq_sg(1:2, ipfsg(3, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1,ifa), isg)) + &
shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
!
r =  0.d0
call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
xrq_sg(1:2, ipfsg(4, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1,ifa), isg)) + &
shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
enddo

!...Physical coordinates
!
do ip = 1,4
r = xrqc(1, ipqsg(ip, isg))
s = xrqc(2, ipqsg(ip, isg))
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xg = 0.d0
yg = 0.d0
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpq_fe(1, ishp)
yg = yg + shpq(ishp)*xpq_fe(2, ishp)
enddo
!
xpq_sg(1, ip) = xg
xpq_sg(2, ip) = yg
enddo
!
do ip = 5, 12
r = xrq_sg(1, ip)
s = xrq_sg(2, ip)
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xg = 0.d0
yg = 0.d0
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpq_fe(1, ishp)
yg = yg + shpq(ishp)*xpq_fe(2, ishp)
enddo
!
xpq_sg(1, ip) = xg
xpq_sg(2, ip) = yg
!
enddo

!...the 9th nodes
xpq_sg(1:2, 9) = -0.25d0*(xpq_sg(1:2, 1) + xpq_sg(1:2, 2) +xpq_sg(1:2, 3) + xpq_sg(1:2, 4)) +&
0.5d0*(xpq_sg(1:2, 5) + xpq_sg(1:2, 6) +xpq_sg(1:2, 7) + xpq_sg(1:2, 8))

!...Initialze parameters...
if(npoly.eq.2)  mmatr = 0.d0
f0 = 0.d0
wt = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(1,9,shpq(1:9), dsprq(1:9), dspsq(1:9), r, s)

dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 9
dxdr = dxdr + dsprq(ishp)*xpq_sg(1,ishp)
dxds = dxds + dspsq(ishp)*xpq_sg(1,ishp)

dydr = dydr + dsprq(ishp)*xpq_sg(2,ishp)
dyds = dyds + dspsq(ishp)*xpq_sg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco
wt = wt +wi
!
!if(ielem.eq.1) print*,'djaco',isg,ig, wi,dxdr*dyds - dydr*dxds,xpq_sg
enddo
!if(ielem.eq.1) print*,'volume',f0,isg,wt
!
!stop
!
if(npoly==1)then
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then
matin(16,isg) = 1.d0/f0
endif

enddo
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, nsubq
!...Reference coordinates
!
do ifa = 1, 4
r =  0.d0
call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
xrq_sg(1:2, ipfsg(3, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1,ifa), isg)) + &
shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
!
r =  0.d0
call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
xrq_sg(1:2, ipfsg(4, ifa)) = shpfl(1)*xrqc(1:2, ipqsg(ipfsg(1,ifa), isg)) + &
shpfl(2)*xrqc(1:2, ipqsg(ipfsg(2, ifa), isg))
enddo

!...Physical coordinates
do ip = 1,4
r = xrqc(1, ipqsg(ip, isg))
s = xrqc(2, ipqsg(ip, isg))
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xg = 0.d0
yg = 0.d0
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpqi_fe(1, ishp)
yg = yg + shpq(ishp)*xpqi_fe(2, ishp)
enddo
!
xpqi_sg(1, ip) = xg
xpqi_sg(2, ip) = yg
enddo

do ip = 5, 12
r = xrq_sg(1, ip)
s = xrq_sg(2, ip)
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xg = 0.d0
yg = 0.d0
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpqi_fe(1, ishp)
yg = yg + shpq(ishp)*xpqi_fe(2, ishp)
enddo
!
xpqi_sg(1, ip) = xg
xpqi_sg(2, ip) = yg
!
enddo
!
!...the 9th nodes
xpqi_sg(1:2, 9) = -0.25d0*(xpqi_sg(1:2, 1) + xpqi_sg(1:2, 2) +xpqi_sg(1:2, 3) + xpqi_sg(1:2, 4)) +&
0.5d0*(xpqi_sg(1:2, 5) + xpqi_sg(1:2, 6) +xpqi_sg(1:2, 7) + xpqi_sg(1:2, 8))
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(1,9,shpq(1:9), dsprq(1:9), dspsq(1:9), r, s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 9
dxdr = dxdr + dsprq(ishp)*xpqi_sg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqi_sg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqi_sg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqi_sg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
bq(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, 9
xgaus = xgaus + shpq(ishp)*xpqi_sg(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqi_sg(2,ishp)
enddo

!...Get the initial density
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*bq(1)*djaco
!
!if(ielem.eq.2100.and.isg.eq.1) print*,'isg4',isg,ig,m(1, 1),rhoi,bq(1),djaco,xcrho,ycrho

enddo
!
if(npoly==1)then
m(1,1) = matin(4, isg)
elseif(npoly==2)then
m(1,1) = matin(16, isg)
endif

!...Update the subcell density ditribution
unint = 0.d0
unint(1) = unint(1) + m(1, 1)*rhsel(1,isg)

unksgq(1, isg)= unint(1)
enddo
!
end subroutine  getrhosubcell_sms_general4
!
!...subroutine: Riemann input for hybrid curved quads using SMS....
!
subroutine getriem_quad_sms2(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:ngesgq,1:nquad),       intent(in)::gesgq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),             intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8, dimension(1:2, 1:2, 1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
real*8, dimension(1:2, 1:2, 1:2, 1:4, 1:nsubq, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:nsubq, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:nsubq, 1:nquad), intent(out)::snsigmlq

!...Local integer
integer::ie, ideg, ielem, ifa, iv, isg, ivsg, ifsg

!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(8, nsubq)::fnqsg
integer,dimension(4, nsubq)::ipqsg
!...local real array
real*8,dimension(1:ndegr, 1:nvqua)::bq
real*8,dimension(1:ndegr, 1:4)::bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nq,1:4)::unsgq
real*8::aujmp(1:3, 1:4)
real*8::vnorm(1:3, 1:2, 1:6)
real*8::sigma(1:2, 1:2, 1:4)
real*8,dimension(1:2, 1:4)::murie
real*8,dimension(1:2, 1:nvqua):: xrq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndegr, 1:nsubq)::unksgq
real*8,dimension(1:nq+1, 1:nsubq)::uqsgc
real*8,dimension(1:4, 1:4)::prsgq
real*8,dimension(1:4, 1:4)::bqvp
real*8,dimension(1:4)::prqz
real*8,dimension(1:7, 1:4)::geoq_sub
real*8,dimension(2, 4, nsubq)::wfgsq
real*8,dimension(nvfac)::wegtf,posif
real*8::bqsg(ndegr)

!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::xcrho, ycrho
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx,rhovsg
real*8::dpvt, drhov
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::rcsg,scsg
real*8::dxp,dyp,xc,yc
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon
!
data eps   / 1.0d-6/
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Specify some Gauss points...
!
if(ncurv.eq.1)then
!...Local vertex No. of gauss points in one unit
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
fnqsg(1, 1) =  8; fnqsg(2, 1) = 1;  fnqsg(3, 1) = 13;
fnqsg(4, 1) =  -9;
fnqsg(5, 1) = -9; fnqsg(6, 1) =12;  fnqsg(7, 1) = 12; fnqsg(8, 1) =  16;

fnqsg(1, 2) =  9; fnqsg(2, 2) =13;
fnqsg(3, 2) =  2; fnqsg(4, 2) =  3;
fnqsg(5, 2) =  14;
fnqsg(6, 2) =-10;  fnqsg(7, 2) =-10; fnqsg(8, 2) = 9;

fnqsg(1, 3) =-11; fnqsg(2, 3) = 10; fnqsg(3, 3) = 10; fnqsg(4, 3) =14;
fnqsg(5, 3) =  4; fnqsg(6, 3) =  5; fnqsg(7, 3) = 15;
fnqsg(8, 3) =-11;

fnqsg(1, 4) = 16;
fnqsg(2, 4) =-12; fnqsg(3, 4) =-12; fnqsg(4, 4) = 11;
fnqsg(5, 4) = 11; fnqsg(6, 4) = 15;
fnqsg(7, 4) = 6;  fnqsg(8, 4) = 7;

!...4/6 for internal face...This is used for our work
wfgsq(1, 1, 1) = 4.d0/6.d0; wfgsq(2, 1, 1) = 4.d0/6.d0;
wfgsq(1, 2, 1) = 4.d0/3.d0; wfgsq(2, 2, 1) = 8.d0/3.d0;
wfgsq(1, 3, 1) = 1.d0;      wfgsq(2, 3, 1) = 1.d0;
wfgsq(1, 4, 1) = 8.d0/3.d0;      wfgsq(2, 4, 1) = 4.d0/3.d0;

wfgsq(1, 1, 2) = 8.d0/3.d0;      wfgsq(2, 1, 2) = 4.d0/3.d0;
wfgsq(1, 2, 2) = 4.d0/6.d0; wfgsq(2, 2, 2) = 4.d0/6.d0;
wfgsq(1, 3, 2) = 4.d0/3.d0; wfgsq(2, 3, 2) = 8.d0/3.d0;
wfgsq(1, 4, 2) = 1.d0;      wfgsq(2, 4, 2) = 1.d0;

wfgsq(1, 1, 3) = 1.d0;      wfgsq(2, 1, 3) = 1.d0;
wfgsq(1, 2, 3) = 8.d0/3.d0;      wfgsq(2, 2, 3) = 4.d0/3.d0;
wfgsq(1, 3, 3) = 4.d0/6.d0; wfgsq(2, 3, 3) = 4.d0/6.d0;
wfgsq(1, 4, 3) = 4.d0/3.d0; wfgsq(2, 4, 3) = 8.d0/3.d0

wfgsq(1, 1, 4) = 4.d0/3.d0; wfgsq(2, 1, 4) = 8.d0/3.d0;
wfgsq(1, 2, 4) = 1.d0;      wfgsq(2, 2, 4) = 1.d0;
wfgsq(1, 3, 4) = 8.d0/3.d0;      wfgsq(2, 3, 4) = 4.d0/3.d0;
wfgsq(1, 4, 4) = 4.d0/6.d0; wfgsq(2, 4, 4) = 4.d0/6.d0;

!...no internal contribution
wfgsq = 0.25d0*wfgsq

!...Vertex coordinate
xrq(1, 1) = -1.d0; xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0; xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0; xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0; xrq(2, 4) =  1.d0;
xrq(1, 5) =  0.d0; xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0; xrq(2, 6) =  0.d0;
xrq(1, 7) =  0.d0; xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0; xrq(2, 8) =  0.d0;
xrq(1, 9) =  0.d0; xrq(2, 9) =  0.d0;

elseif(ncurv.eq.2)then
!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8
!
fnqsg(1, 1) =  8; fnqsg(2, 1) =  1;  fnqsg(3, 1) = 17;  fnqsg(4, 1) =  -9;
fnqsg(5, 1) = -9; fnqsg(6, 1) = 16;  fnqsg(7, 1) = 16;  fnqsg(8, 1) =  24;

fnqsg(1, 2) =  9; fnqsg(2, 2) = 17;  fnqsg(3, 2) = 18;  fnqsg(4, 2) = -10;
fnqsg(5, 2) = -10;fnqsg(6, 2) =-10;  fnqsg(7, 2) =  9;  fnqsg(8, 2) = 9;

fnqsg(1, 3) = 10; fnqsg(2, 3) = 18;  fnqsg(3, 3) =  2;  fnqsg(4, 3) = 3;
fnqsg(5, 3) = 19; fnqsg(6, 3) =-11;  fnqsg(7, 3) =-11;  fnqsg(8, 3) =10;

fnqsg(1, 4) = 11; fnqsg(2, 4) = 11;  fnqsg(3, 4) = 11;  fnqsg(4, 4) = 19;
fnqsg(5, 4) = 20; fnqsg(6, 4) =-12;  fnqsg(7, 4) =-12;  fnqsg(8, 4)=-12;

fnqsg(1, 5) =-13; fnqsg(2, 5) = 12;  fnqsg(3, 5) = 12;  fnqsg(4, 5) = 20;
fnqsg(5, 5) =  4; fnqsg(6, 5) =  5;  fnqsg(7, 5) = 21;  fnqsg(8, 5)=-13;

fnqsg(1, 6) =-14; fnqsg(2, 6) =-14;  fnqsg(3, 6) = 13;  fnqsg(4, 6) = 13;
fnqsg(5, 6) = 13; fnqsg(6, 6) = 21;  fnqsg(7, 6) = 22;  fnqsg(8, 6) =-14;

fnqsg(1, 7) = 23; fnqsg(2, 7) =-15;  fnqsg(3, 7) =-15;  fnqsg(4, 7) = 14;
fnqsg(5, 7) = 14; fnqsg(6, 7) = 22;  fnqsg(7, 7) =  6;  fnqsg(8, 7) =  7;

fnqsg(1, 8) = 24; fnqsg(2, 8) =-16;  fnqsg(3, 8) =-16;  fnqsg(4, 8) =-16;
fnqsg(5, 8) =-16; fnqsg(6, 8) = 15;  fnqsg(7, 8) = 15;  fnqsg(8, 8) =23;

!...
wegtf(1) = 1.d0/12.d0;
wegtf(2) = 5.d0/12.d0;
wegtf(3) = 5.d0/12.d0;
wegtf(4) = 1.d0/12.d0;
!...4/6 for internal face...This is used for our work
wfgsq(1, 1, 1) = wegtf(1);      wfgsq(2, 1, 1) = wegtf(1);
wfgsq(1, 2, 1) = wegtf(2)/2.d0; wfgsq(2, 2, 1) = wegtf(2)/2.d0
wfgsq(1, 3, 1) = 1.d0;          wfgsq(2, 3, 1) = 1.d0;
wfgsq(1, 4, 1) = wegtf(2)/2.d0; wfgsq(2, 4, 1) = wegtf(2)/2.d0;

wfgsq(1, 1, 2) = wegtf(2)/2.d0; wfgsq(2, 1, 2) = wegtf(2)/2.d0;
wfgsq(1, 2, 2) = wegtf(2)/2.d0; wfgsq(2, 2, 2) = wegtf(2)/2.d0;
wfgsq(1, 3, 2) = 1.d0;          wfgsq(2, 3, 2) = 1.d0;
wfgsq(1, 4, 2) = 1.d0;          wfgsq(2, 4, 2) = 1.d0;

wfgsq(1, 1, 3) = wegtf(2)/2.d0; wfgsq(2, 1, 3) = wegtf(2)/2.d0;
wfgsq(1, 2, 3) = wegtf(1);      wfgsq(2, 2, 3) = wegtf(1);
wfgsq(1, 3, 3) = wegtf(2)/2.d0; wfgsq(2, 3, 3) = wegtf(2)/2.d0;
wfgsq(1, 4, 3) = 1.d0;          wfgsq(2, 4, 3) = 1.d0

wfgsq(1, 1, 4) = 1.d0;          wfgsq(2, 1, 4) = 1.d0;
wfgsq(1, 2, 4) = wegtf(2)/2.d0; wfgsq(2, 2, 4) = wegtf(2)/2.d0;
wfgsq(1, 3, 4) = wegtf(2)/2.d0; wfgsq(2, 3, 4) = wegtf(2)/2.d0;
wfgsq(1, 4, 4) = 1.d0;          wfgsq(2, 4, 4) = 1.d0;

wfgsq(1, 1, 5) = 1.d0;          wfgsq(2, 1, 5) = 1.d0;
wfgsq(1, 2, 5) = wegtf(2)/2.d0; wfgsq(2, 2, 5) = wegtf(2)/2.d0
wfgsq(1, 3, 5) = wegtf(1);      wfgsq(2, 3, 5) = wegtf(1);
wfgsq(1, 4, 5) = wegtf(2)/2.d0; wfgsq(2, 4, 5) = wegtf(2)/2.d0;

wfgsq(1, 1, 6) = 1.d0;          wfgsq(2, 1, 6) = 1.d0;
wfgsq(1, 2, 6) = 1.d0;          wfgsq(2, 2, 6) = 1.d0;
wfgsq(1, 3, 6) = wegtf(2)/2.d0; wfgsq(2, 3, 6) = wegtf(2)/2.d0;
wfgsq(1, 4, 6) = wegtf(2)/2.d0; wfgsq(2, 4, 6) = wegtf(2)/2.d0;

wfgsq(1, 1, 7) = wegtf(2)/2.d0; wfgsq(2, 1, 7) = wegtf(2)/2.d0;
wfgsq(1, 2, 7) = 1.d0;          wfgsq(2, 2, 7) = 1.d0;
wfgsq(1, 3, 7) = wegtf(2)/2.d0; wfgsq(2, 3, 7) = wegtf(2)/2.d0;
wfgsq(1, 4, 7) = wegtf(1);      wfgsq(2, 4, 7) = wegtf(1)

wfgsq(1, 1, 8) = wegtf(2)/2.d0; wfgsq(2, 1, 8) = wegtf(2)/2.d0;
wfgsq(1, 2, 8) = 1.d0;          wfgsq(2, 2, 8) = 1.d0;
wfgsq(1, 3, 8) = 1.d0;          wfgsq(2, 3, 8) = 1.d0;
wfgsq(1, 4, 8) = wegtf(2)/2.d0; wfgsq(2, 4, 8) = wegtf(2)/2.d0;

!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;
!...Vertex coordinate
xrq(1, 1) = -1.d0; xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0; xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0; xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0; xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;      xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;      xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;      xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;      xrq(2,12) = posif(2);

endif
!
!...Part II: Loop over every quad...
!
do 350 ie = 1,nquad !...(1)ie = 1,nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria

!...shape functions
dr = 1.0d0
ds = 1.0d0

!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...Basis function
do iv =1 ,nvqua
bq(1, iv) = 1.d0

if(npoly.ge.1)then
bq(2, iv) = (xrq(1, iv)-rc)/dr
bq(3, iv) = (xrq(2, iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
endif
endif
enddo

!...cell averaged value...
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

!... u, v, e at every vertex....
unknvq = 0.d0
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
enddo

!...Get density for sub-cells...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!...Update the coordinates for FE
!call getcoord_fe(ncurv, nvfac, nvqua, xpq)
!call getcoord_fe(ncurv, nvfac, nvqua, xpqi)
!
call getcoord_fe(ncurv, nvfac, nvqua, xpq)
call getcoord_fe(ncurv, nvfac, nvqua, xpqi)
!
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

!...Get density correction
call getrhosubcell_sms_general5(xcrho, ycrho, xpq,xpqi,unksgq, ielem)
!
!call getcoord_fe(ncurv, nvfac, nvqua, xpq)
!call getcoord_fe(ncurv, nvfac, nvqua, xpqi)

call getrhosubcell_daverg(rc,sc, geoel(19:21, ielem), xpq, unksgq,  unkno(:,:,ielem), aflim(:,ielem), uqsgc,ielem)

!...Output for debugging
!if(ielem.eq.1)then
!print*,'Variabe',ielem,ivsg,unksgq(1,1:8)
!print*,'Variab2',ielem,ivsg,1.d0/uqsgc(1, 1:8)
!endif


!...II.1: Loop over sub-cells....
do isg = 1, nsubq

!...Normal vector for every quadrature point...
vnorm(1:2, 1, 1) = sign(1, fnqsg(1, isg))*gesgq(1:2,abs(fnqsg(1, isg)), ie); vnorm(  3, 1, 1) = gesgq(3,abs(fnqsg(1, isg)), ie);
vnorm(1:2, 2, 1) = sign(1, fnqsg(2, isg))*gesgq(1:2,abs(fnqsg(2, isg)), ie); vnorm(  3, 2, 1) = gesgq(3,abs(fnqsg(2, isg)), ie);

vnorm(1:2, 1, 2) = sign(1, fnqsg(3, isg))*gesgq(1:2,abs(fnqsg(3, isg)), ie); vnorm(  3, 1, 2) = gesgq(3,abs(fnqsg(3, isg)), ie);
vnorm(1:2, 2, 2) = sign(1, fnqsg(4, isg))*gesgq(1:2,abs(fnqsg(4, isg)), ie); vnorm(  3, 2, 2) = gesgq(3,abs(fnqsg(4, isg)), ie);

vnorm(1:2, 1, 3) = sign(1, fnqsg(5, isg))*gesgq(1:2,abs(fnqsg(5, isg)), ie); vnorm(  3, 1, 3) = gesgq(3,abs(fnqsg(5, isg)), ie);
vnorm(1:2, 2, 3) = sign(1, fnqsg(6, isg))*gesgq(1:2,abs(fnqsg(6, isg)), ie); vnorm(  3, 2, 3) = gesgq(3,abs(fnqsg(6, isg)), ie);

vnorm(1:2, 1, 4) = sign(1, fnqsg(7, isg))*gesgq(1:2,abs(fnqsg(7, isg)), ie); vnorm(  3, 1, 4) = gesgq(3,abs(fnqsg(7, isg)), ie);
vnorm(1:2, 2, 4) = sign(1, fnqsg(8, isg))*gesgq(1:2,abs(fnqsg(8, isg)), ie); vnorm(  3, 2, 4) = gesgq(3,abs(fnqsg(8, isg)), ie);

!...Get weighted area normal vector
do ifsg =1,2
vnorm(3, ifsg, 1) =  vnorm(  3, ifsg, 1)*wfgsq(ifsg, 1, isg);
vnorm(3, ifsg, 2) =  vnorm(  3, ifsg, 2)*wfgsq(ifsg, 2, isg);
vnorm(3, ifsg, 3) =  vnorm(  3, ifsg, 3)*wfgsq(ifsg, 3, isg);
vnorm(3, ifsg, 4) =  vnorm(  3, ifsg, 4)*wfgsq(ifsg, 4, isg);
enddo

!...Correct the density for Reimann input
do ivsg = 1,4

if(ipqsg(ivsg, isg).gt.12) cycle

if(ndens.eq.1) then
!print*,'rhovt',ie,ipqsg(ivsg, isg)
rhovt = 1.d0/unknvq(1, ipqsg(ivsg, isg))
rhovsg = rhovt +cdrho*(unksgq(1,isg)-1.d0/uqsgc(1, isg))
!
!print*,'bad'
!rhovsg = unksgq(1,isg)!-1.d0/uqsgc(1, isg))
!if(rhovsg.le.1d-6)then
!rhovsg = rhovt
!endif
!rhovsg = 1.d0/uqsgc(1, isg)
elseif(ndens.eq.3)then
print*,'ndens.eq.3 will be implemented in future in getriem_quad_sms!'
endif

uvtx = unknvq(2, ipqsg(ivsg, isg))
vvtx = unknvq(3, ipqsg(ivsg, isg))
evtx = unknvq(4, ipqsg(ivsg, isg))
!
pvtx = max(eps, (gamlg-1.d0)*rhovsg*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx

!...Correction density
!cdrho = 0.3d0
!drhov = cdrho*(unksgq(1,isg)-1.d0/uqsgc(1, isg))
!dpvt = (gamlg-1.d0)*drhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2))
!
!cimpd = 1.d0

!...Output for debugging
!if(ielem.eq.1)then
!print*,'smooth',ielem,unkno(:,1,ielem)
!print*,'Variabe',ielem,ivsg,rhovsg,unksgq(1,isg),1.d0/uqsgc(1, isg),pvtx,evtx,uvtx, vvtx
!endif
!
!...Output for debugging
!if(ipq(ipqsg(ivsg, isg)).eq.33)then
!print*,'vaue',ipq(ipqsg(ivsg, isg)),ielem,isg,ivsg,pctr
!print*,'uvelo',rhovsg,uvtx,vvtx,evtx,pvtx,rhovt,unksgq(1,isg),1.d0/uqsgc(1, isg)
!print*,'pres',pctr + aflim(4, ielem)*(pvtx - pctr),&
!(gamlg-1.d0)*cdrho*(unksgq(1,isg)-1.d0/uqsgc(1, isg))*(evtx - 0.5d0*(uvtx**2 + vvtx**2))
!endif
!...Limiter
!...For high-order method(>P1), if the cell includes shock, then
!...the high-order terms are chopped, and only limit the linear DG(P1) part.


if(nlimi.eq.6.and.geoel(10, ielem).gt.10.d0)then
if(ndens.eq.1)then
rhovsg = 1.d0/(rhomc + aflim(1, ielem)*(1.d0/rhovsg - rhomc))
elseif(ndens.eq.3)then
rhovsg = unksgq(1,isg) + aflim(1, ielem)*(rhovsg - unksgq(1,isg))
endif

!...Output for debugging
!print*,'Shock identification',ielem,unkno(:,:,ielem)
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bq(2, ipqsg(ivsg, isg)) + duds*bq(3, ipqsg(ivsg, isg))
vvtx = unkno(1,3,ielem)  + dvdr*bq(2, ipqsg(ivsg, isg)) + dvds*bq(3, ipqsg(ivsg, isg))
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
evtx = ectr + aflim(5, ielem)*(evtx - ectr)

!...Correction density
!cdrho = 0.8d0
!drhov = cdrho*(unksgq(1,isg)-1.d0/uqsgc(1, isg))
!cimpd = 2.d0
!pvtx = 0.d0
!drhov = unksgq(1,isg)


!...Correction pressure
!dpvt = (gamlg-1.d0)*drhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2))

!...Updtae velocity
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx
endif

!...Final pressure
!pvtx = pvtx + dpvt
!...Get stress tensor at one vertex
sigma(1, 1, ivsg) = -pvtx
sigma(1, 2, ivsg) = 0.d0
sigma(2, 1, ivsg) = 0.d0
sigma(2, 2, ivsg) = -pvtx
!

!...Get the a_c (unit vector)
aujmp(1:2, ivsg) = vlave(1:2, ipq(ipqsg(ivsg, isg))) - unsgq(2:3, ivsg)

acnx = aujmp(1, ivsg)
acny = aujmp(2, ivsg)
!
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ivsg) = 1.e-11!0.d0;
else
aujmp(1:2, ivsg) = aujmp(1:2, ivsg)/sqrt(acnx**2 + acny**2)
endif
!
aujmp(3, ivsg) = sqrt(acnx**2 + acny**2)
enddo

!...Get the variables at the center...
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedence coefficient...
do ivsg   = 1, 4
!
if(ipqsg(ivsg, isg).gt.12) cycle
!
dux= vlave(1, ipq(ipqsg(ivsg, isg)))-unsgq(2, ivsg)
duy= vlave(2, ipq(ipqsg(ivsg, isg)))-unsgq(3, ivsg)
!deltu = sqrt(dux**2 + duy**2)
do ifa = 1, 2
deltu = abs(dux*vnorm(1, ifa, ivsg) + duy*vnorm(2, ifa, ivsg))
murie(ifa, ivsg) = rhoct*sdctr + cimpd*rhoct*slpdu*deltu
!murie(ifa, ivsg) = uqsgc(1, isg)*sdctr + uqsgc(1, isg)*slpdu*deltu
!...The exact shock impedance of 2 shock model
!murie(ifa, ivsg) = rhoct*slpdu*deltu/2.d0+&
!     rhoct*sqrt((slpdu*deltu/2.d0)**2 + gamlg*pctr/rhoct)
enddo
enddo

!...Feed the input into Riemann solver
do ivsg  = 1, 4

!...Local vertex No. of gauss points in one unit
iv = ipqsg(ivsg, isg)

if(ipqsg(ivsg, isg).gt.12) cycle
!
do ifa = 1, 2 !...Every corner consists of 2 faces...

!...Call Riemann solver...
call getriecoef_matrixnew(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:3, ivsg), &
unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:2, ivsg), &
!unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipq(iv)) = munacn(1:2, 1, ipq(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipq(iv)) = munacn(1:2, 2, ipq(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipq(iv)) = munacu(1:2, ipq(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipq(iv)) = snsigm(1:2, ipq(iv)) + snsigm_rie(1:2)

!...Output for debugging
!if(ipq(iv).eq.33) print*,'p36 muacn(vv) post',ipq(iv),ie,ifa,isg,ivsg,snsigm_rie(1:2),&
!vnorm(1:3, ifa, ivsg),&
!sigma(1, 1, ivsg),snsigm(1, ipq(iv)),&
!murie(ifa, ivsg),unsgq(2:3, ivsg)

!...Local variable...
munaclq(1:2, 1, ifa, ivsg, isg, ie) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, ifa, ivsg, isg, ie) =  munacn_rie(1:2, 2)
!
munaulq(1:2,    ifa, ivsg, isg, ie) =  munacu_rie(1:2)
!
snsigmlq(1:2,   ifa, ivsg, isg, ie)=  snsigm_rie(1:2)
!
enddo
enddo
!
enddo
!
350 enddo  !...(1)ie = 1,nquad

end subroutine getriem_quad_sms2
!
!...Face integral for hybrid quad using SMS...
!
subroutine rhsifacedg_lag_sms2(ipqua, unkno, ustar,fstarq, gesgq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:4, 1:nsubq, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem,isg,ivsg,ifsg
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua)  :: ipq
integer,dimension(8, nsubq)::fnqsg
integer,dimension(4, nsubq)::ipqsg
real*8,dimension(2,4,nsubq)::wfgsq,wfgsqm
real*8,dimension(1:3,1:2,1:4)::vnorm
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8, dimension(1:ndimn, 1:ndegr, 1:2, 1:nvqua)::lpnpq
real*8, dimension(1:ndimn,nvqua)::xrq
real*8::bq(1:ndegr,1:nvqua)
real*8,dimension(nvfac)::wegtf,posif

!...local real number
real*8::eps,c00,c05,c10,c20,c13,c16
real*8::dr,ds,rc,sc
real*8::dwav1,dwav2
real*8::anx, any
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
data c16   / 0.1666666666666666d0 /
!
!...Part I: preliminary setup
!
if(ncurv.eq.1)then
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
!
fnqsg(1, 1) =  8; fnqsg(2, 1) = 1;  fnqsg(3, 1) = 13;
fnqsg(4, 1) =  -9;
fnqsg(5, 1) = -9; fnqsg(6, 1) =12;  fnqsg(7, 1) = 12; fnqsg(8, 1) =  16;

fnqsg(1, 2) =  9; fnqsg(2, 2) =13;
fnqsg(3, 2) =  2; fnqsg(4, 2) =  3;
fnqsg(5, 2) =  14;
fnqsg(6, 2) =-10;  fnqsg(7, 2) =-10; fnqsg(8, 2) = 9;

fnqsg(1, 3) =-11; fnqsg(2, 3) = 10; fnqsg(3, 3) = 10; fnqsg(4, 3) =14;
fnqsg(5, 3) =  4; fnqsg(6, 3) =  5; fnqsg(7, 3) = 15;
fnqsg(8, 3) =-11;

fnqsg(1, 4) = 16;
fnqsg(2, 4) =-12; fnqsg(3, 4) =-12; fnqsg(4, 4) = 11;
fnqsg(5, 4) = 11; fnqsg(6, 4) = 15;
fnqsg(7, 4) = 6;  fnqsg(8, 4) = 7;
!
wfgsq(1, 1, 1) = 1.d0; wfgsq(2, 1, 1) = 1.d0;
wfgsq(1, 2, 1) = 1.d0; wfgsq(2, 2, 1) = 1.d0;
wfgsq(1, 3, 1) = 0.d0; wfgsq(2, 3, 1) = 0.d0;
wfgsq(1, 4, 1) = 1.d0; wfgsq(2, 4, 1) = 1.d0;

wfgsq(1, 1, 2) = 1.d0; wfgsq(2, 1, 2) =1.d0;
wfgsq(1, 2, 2) = 1.d0; wfgsq(2, 2, 2) =1.d0;
wfgsq(1, 3, 2) = 1.d0; wfgsq(2, 3, 2) =1.d0;
wfgsq(1, 4, 2) = 0.d0; wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 1, 3) = 0.d0; wfgsq(2, 1, 3) = 0.d0;
wfgsq(1, 2, 3) = 1.d0; wfgsq(2, 2, 3) = 1.d0;
wfgsq(1, 3, 3) = 1.d0; wfgsq(2, 3, 3) = 1.d0;
wfgsq(1, 4, 3) = 1.d0; wfgsq(2, 4, 3) = 1.d0;

wfgsq(1, 1, 4) = 1.d0; wfgsq(2, 1, 4) = 1.d0;
wfgsq(1, 2, 4) = 0.d0; wfgsq(2, 2, 4) = 0.d0;
wfgsq(1, 3, 4) = 1.d0; wfgsq(2, 3, 4) = 1.d0;
wfgsq(1, 4, 4) = 1.d0; wfgsq(2, 4, 4) = 1.d0;!
!...New
!
!
!...Taking the weight into consideration when calculate the ustar...
!
!wfgsq = 1.d0!0.25d0*wfgsq1
!
!...Conitnuity equation...
!
wfgsqm(1, 1, 1) = 4.d0/6.d0; wfgsqm(2, 1, 1) = 4.d0/6.d0;
wfgsqm(1, 2, 1) = 4.d0/3.d0; wfgsqm(2, 2, 1) = 0.d0;
wfgsqm(1, 3, 1) = 0.d0;      wfgsqm(2, 3, 1) = 0.d0;
wfgsqm(1, 4, 1) = 0.d0;      wfgsqm(2, 4, 1) = 4.d0/3.d0;

wfgsqm(1, 1, 2) = 0.d0/3.d0; wfgsqm(2, 1, 2) = 4.d0/3.d0;
wfgsqm(1, 2, 2) = 4.d0/6.d0; wfgsqm(2, 2, 2) = 4.d0/6.d0;
wfgsqm(1, 3, 2) = 4.d0/3.d0; wfgsqm(2, 3, 2) = 0.d0/3.d0;
wfgsqm(1, 4, 2) = 0.d0;      wfgsqm(2, 4, 2) = 0.d0;

wfgsqm(1, 1, 3) = 0.d0;      wfgsqm(2, 1, 3) = 0.d0;
wfgsqm(1, 2, 3) = 0.d0/3.d0; wfgsqm(2, 2, 3) = 4.d0/3.d0;
wfgsqm(1, 3, 3) = 4.d0/6.d0; wfgsqm(2, 3, 3) = 4.d0/6.d0;
wfgsqm(1, 4, 3) = 4.d0/3.d0; wfgsqm(2, 4, 3) = 0.d0/3.d0;

wfgsqm(1, 1, 4) = 4.d0/3.d0; wfgsqm(2, 1, 4) = 0.d0/3.d0;
wfgsqm(1, 2, 4) = 0.d0;      wfgsqm(2, 2, 4) = 0.d0;
wfgsqm(1, 3, 4) = 0.d0/3.d0; wfgsqm(2, 3, 4) = 4.d0/3.d0;
wfgsqm(1, 4, 4) = 4.d0/6.d0; wfgsqm(2, 4, 4) = 4.d0/6.d0;
!
wfgsqm = 0.25d0*wfgsqm

!...Vertex coordinate
xrq(1, 1) = -1.d0; xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0; xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0; xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0; xrq(2, 4) =  1.d0;
xrq(1, 5) =  0.d0; xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0; xrq(2, 6) =  0.d0;
xrq(1, 7) =  0.d0; xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0; xrq(2, 8) =  0.d0;
xrq(1, 9) =  0.d0; xrq(2, 9) =  0.d0;

elseif(ncurv.eq.2)then

ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8
!
fnqsg(1, 1) =  8; fnqsg(2, 1) =  1;  fnqsg(3, 1) = 17;  fnqsg(4, 1) =  -9;
fnqsg(5, 1) = -9; fnqsg(6, 1) = 16;  fnqsg(7, 1) = 16;  fnqsg(8, 1) =  24;

fnqsg(1, 2) =  9; fnqsg(2, 2) = 17;  fnqsg(3, 2) = 18;  fnqsg(4, 2) = -10;
fnqsg(5, 2) = -10;fnqsg(6, 2) =-10;  fnqsg(7, 2) =  9;  fnqsg(8, 2) = 9;

fnqsg(1, 3) = 10; fnqsg(2, 3) = 18;  fnqsg(3, 3) =  2;  fnqsg(4, 3) = 3;
fnqsg(5, 3) = 19; fnqsg(6, 3) =-11;  fnqsg(7, 3) =-11;  fnqsg(8, 3) =10;

fnqsg(1, 4) = 11; fnqsg(2, 4) = 11;  fnqsg(3, 4) = 11;  fnqsg(4, 4) = 19;
fnqsg(5, 4) = 20; fnqsg(6, 4) =-12;  fnqsg(7, 4) =-12;  fnqsg(8, 4)=-12;

fnqsg(1, 5) =-13; fnqsg(2, 5) = 12;  fnqsg(3, 5) = 12;  fnqsg(4, 5) = 20;
fnqsg(5, 5) =  4; fnqsg(6, 5) =  5;  fnqsg(7, 5) = 21;  fnqsg(8, 5)=-13;

fnqsg(1, 6) =-14; fnqsg(2, 6) =-14;  fnqsg(3, 6) = 13;  fnqsg(4, 6) = 13;
fnqsg(5, 6) = 13; fnqsg(6, 6) = 21;  fnqsg(7, 6) = 22;  fnqsg(8, 6) =-14;

fnqsg(1, 7) = 23; fnqsg(2, 7) =-15;  fnqsg(3, 7) =-15;  fnqsg(4, 7) = 14;
fnqsg(5, 7) = 14; fnqsg(6, 7) = 22;  fnqsg(7, 7) =  6;  fnqsg(8, 7) =  7;

fnqsg(1, 8) = 24; fnqsg(2, 8) =-16;  fnqsg(3, 8) =-16;  fnqsg(4, 8) =-16;
fnqsg(5, 8) =-16; fnqsg(6, 8) = 15;  fnqsg(7, 8) = 15;  fnqsg(8, 8) =23;
!
wegtf(1) = 1.d0/12.d0;
wegtf(2) = 5.d0/12.d0;
wegtf(3) = 5.d0/12.d0;
wegtf(4) = 1.d0/12.d0;
!
wfgsq = 1.d0
!...4/6 for internal face...This is used for our work
wfgsq(1, 3, 1) = 0.d0;          wfgsq(2, 3, 1) = 0.d0;

wfgsq(1, 3, 2) = 0.d0;          wfgsq(2, 3, 2) = 0.d0;
wfgsq(1, 4, 2) = 0.d0;          wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 4, 3) = 0.d0;          wfgsq(2, 4, 3) = 0.d0

wfgsq(1, 1, 4) = 0.d0;          wfgsq(2, 1, 4) = 0.d0;
wfgsq(1, 4, 4) = 0.d0;          wfgsq(2, 4, 4) = 0.d0;

wfgsq(1, 1, 5) = 0.d0;          wfgsq(2, 1, 5) = 0.d0;

wfgsq(1, 1, 6) = 0.d0;          wfgsq(2, 1, 6) = 0.d0;
wfgsq(1, 2, 6) = 0.d0;          wfgsq(2, 2, 6) = 0.d0;

wfgsq(1, 2, 7) = 0.d0;          wfgsq(2, 2, 7) = 0.d0;

wfgsq(1, 2, 8) = 0.d0;          wfgsq(2, 2, 8) = 0.d0;
wfgsq(1, 3, 8) = 0.d0;          wfgsq(2, 3, 8) = 0.d0;

!
!...4/6 for internal face...This is used for our work
wfgsqm(1, 1, 1) = wegtf(1);      wfgsqm(2, 1, 1) = wegtf(1);
wfgsqm(1, 2, 1) = wegtf(2)/2.d0; wfgsqm(2, 2, 1) = wegtf(2)/2.d0
wfgsqm(1, 3, 1) = 1.d0;          wfgsqm(2, 3, 1) = 1.d0;
wfgsqm(1, 4, 1) = wegtf(2)/2.d0; wfgsqm(2, 4, 1) = wegtf(2)/2.d0;

wfgsqm(1, 1, 2) = wegtf(2)/2.d0; wfgsqm(2, 1, 2) = wegtf(2)/2.d0;
wfgsqm(1, 2, 2) = wegtf(2)/2.d0; wfgsqm(2, 2, 2) = wegtf(2)/2.d0;
wfgsqm(1, 3, 2) = 1.d0;          wfgsqm(2, 3, 2) = 1.d0;
wfgsqm(1, 4, 2) = 1.d0;          wfgsqm(2, 4, 2) = 1.d0;

wfgsqm(1, 1, 3) = wegtf(2)/2.d0; wfgsqm(2, 1, 3) = wegtf(2)/2.d0;
wfgsqm(1, 2, 3) = wegtf(1);      wfgsqm(2, 2, 3) = wegtf(1);
wfgsqm(1, 3, 3) = wegtf(2)/2.d0; wfgsqm(2, 3, 3) = wegtf(2)/2.d0;
wfgsqm(1, 4, 3) = 1.d0;          wfgsqm(2, 4, 3) = 1.d0

wfgsqm(1, 1, 4) = 1.d0;          wfgsqm(2, 1, 4) = 1.d0;
wfgsqm(1, 2, 4) = wegtf(2)/2.d0; wfgsqm(2, 2, 4) = wegtf(2)/2.d0;
wfgsqm(1, 3, 4) = wegtf(2)/2.d0; wfgsqm(2, 3, 4) = wegtf(2)/2.d0;
wfgsqm(1, 4, 4) = 1.d0;          wfgsqm(2, 4, 4) = 1.d0;

wfgsqm(1, 1, 5) = 1.d0;          wfgsqm(2, 1, 5) = 1.d0;
wfgsqm(1, 2, 5) = wegtf(2)/2.d0; wfgsqm(2, 2, 5) = wegtf(2)/2.d0
wfgsqm(1, 3, 5) = wegtf(1);      wfgsqm(2, 3, 5) = wegtf(1);
wfgsqm(1, 4, 5) = wegtf(2)/2.d0; wfgsqm(2, 4, 5) = wegtf(2)/2.d0;

wfgsqm(1, 1, 6) = 1.d0;          wfgsqm(2, 1, 6) = 1.d0;
wfgsqm(1, 2, 6) = 1.d0;          wfgsqm(2, 2, 6) = 1.d0;
wfgsqm(1, 3, 6) = wegtf(2)/2.d0; wfgsqm(2, 3, 6) = wegtf(2)/2.d0;
wfgsqm(1, 4, 6) = wegtf(2)/2.d0; wfgsqm(2, 4, 6) = wegtf(2)/2.d0;

wfgsqm(1, 1, 7) = wegtf(2)/2.d0; wfgsqm(2, 1, 7) = wegtf(2)/2.d0;
wfgsqm(1, 2, 7) = 1.d0;          wfgsqm(2, 2, 7) = 1.d0;
wfgsqm(1, 3, 7) = wegtf(2)/2.d0; wfgsqm(2, 3, 7) = wegtf(2)/2.d0;
wfgsqm(1, 4, 7) = wegtf(1);      wfgsqm(2, 4, 7) = wegtf(1)

wfgsqm(1, 1, 8) = wegtf(2)/2.d0; wfgsqm(2, 1, 8) = wegtf(2)/2.d0;
wfgsqm(1, 2, 8) = 1.d0;          wfgsqm(2, 2, 8) = 1.d0;
wfgsqm(1, 3, 8) = 1.d0;          wfgsqm(2, 3, 8) = 1.d0;
wfgsqm(1, 4, 8) = wegtf(2)/2.d0; wfgsqm(2, 4, 8) = wegtf(2)/2.d0;

!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;

!...Vertex coordinate
xrq(1, 1) = -1.d0; xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0; xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0; xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0; xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;      xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;      xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;      xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;      xrq(2,12) = posif(2);
endif
!
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...The vertex constituting one cell...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvqua
!
bq(1, iv) = 1.d0
!
if(npoly.ge.1)then
bq(2, iv) = (xrq(1, iv)-rc)/dr
bq(3, iv) = (xrq(2, iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
endif
endif
enddo

!...Initialize ulnpn, plnpn, elnpn
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
do isg = 1, nsubq

!...Normal vector for every face...
vnorm(1:2, 1, 1) = sign(1, fnqsg(1, isg))*gesgq(1:2,abs(fnqsg(1, isg)), ie); vnorm(  3, 1, 1) = gesgq(3,abs(fnqsg(1, isg)), ie);
vnorm(1:2, 2, 1) = sign(1, fnqsg(2, isg))*gesgq(1:2,abs(fnqsg(2, isg)), ie); vnorm(  3, 2, 1) = gesgq(3,abs(fnqsg(2, isg)), ie);

vnorm(1:2, 1, 2) = sign(1, fnqsg(3, isg))*gesgq(1:2,abs(fnqsg(3, isg)), ie); vnorm(  3, 1, 2) = gesgq(3,abs(fnqsg(3, isg)), ie);
vnorm(1:2, 2, 2) = sign(1, fnqsg(4, isg))*gesgq(1:2,abs(fnqsg(4, isg)), ie); vnorm(  3, 2, 2) = gesgq(3,abs(fnqsg(4, isg)), ie);

vnorm(1:2, 1, 3) = sign(1, fnqsg(5, isg))*gesgq(1:2,abs(fnqsg(5, isg)), ie); vnorm(  3, 1, 3) = gesgq(3,abs(fnqsg(5, isg)), ie);
vnorm(1:2, 2, 3) = sign(1, fnqsg(6, isg))*gesgq(1:2,abs(fnqsg(6, isg)), ie); vnorm(  3, 2, 3) = gesgq(3,abs(fnqsg(6, isg)), ie);

vnorm(1:2, 1, 4) = sign(1, fnqsg(7, isg))*gesgq(1:2,abs(fnqsg(7, isg)), ie); vnorm(  3, 1, 4) = gesgq(3,abs(fnqsg(7, isg)), ie);
vnorm(1:2, 2, 4) = sign(1, fnqsg(8, isg))*gesgq(1:2,abs(fnqsg(8, isg)), ie); vnorm(  3, 2, 4) = gesgq(3,abs(fnqsg(8, isg)), ie);
!
do ifsg =1,2
vnorm(3, ifsg, 1) =  vnorm(  3, ifsg, 1)*wfgsqm(ifsg, 1, isg);
vnorm(3, ifsg, 2) =  vnorm(  3, ifsg, 2)*wfgsqm(ifsg, 2, isg);
vnorm(3, ifsg, 3) =  vnorm(  3, ifsg, 3)*wfgsqm(ifsg, 3, isg);
vnorm(3, ifsg, 4) =  vnorm(  3, ifsg, 4)*wfgsqm(ifsg, 4, isg);
enddo
!
!...loop over every subgrid vertex...
!
do ivsg = 1, 4
!
if(ipqsg(ivsg, isg).gt.12) cycle
!
ulnpn(1:ndegr) = ulnpn(1:ndegr)  +&
ustar(1, ipq(ipqsg(ivsg, isg)))*vnorm(1, 1, ivsg)*vnorm(  3, 1, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg)) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*vnorm(2, 1, ivsg)*vnorm(  3, 1, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg)) +&
ustar(1, ipq(ipqsg(ivsg, isg)))*vnorm(1, 2, ivsg)*vnorm(  3, 2, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg)) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*vnorm(2, 2, ivsg)*vnorm(  3, 2, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg))
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, 1, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(1, ivsg, isg) +&
fstarq(1, 2, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(2, ivsg, isg)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, 1, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(1, ivsg, isg) +&
fstarq(2, 2, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(2, ivsg, isg)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(ipqsg(ivsg, isg)))*fstarq(1, 1, ivsg, isg,  ie)*&
bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(1, ivsg, isg) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*fstarq(2, 1, ivsg, isg,  ie)*&
bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(1, ivsg, isg) +&
ustar(1, ipq(ipqsg(ivsg, isg)))*fstarq(1, 2, ivsg, isg, ie)*&
bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(2, ivsg, isg)  +&
ustar(2, ipq(ipqsg(ivsg, isg)))*fstarq(2, 2, ivsg, isg,  ie)*&
bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(2, ivsg, isg)

!
!...Output foe debugging
!if(ie==1)  print*,'rhs iface2',isg,ivsg,ipq(ipqsg(ivsg, isg)),ie,fstarq(1:2, 1:2, ivsg, isg, ie),&
!ustar(1:2, ipq(ipqsg(ivsg, isg)))

enddo
!
enddo
!
!...Distribute to every corner...
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)

!...Output foe debugging
!if(ie==1)  print*,'rhs iface',ielem, ie, plnpn(1, 1:ndegr)
650 enddo
!
end subroutine rhsifacedg_lag_sms2
!
!...Get the averaged-density within one sub-cell using SMS ...
!
subroutine  getrhosubcell_sms_general5(xcrho, ycrho, xpq_fe,xpqi_fe,unksgq, ielem)
use constant
implicit none
!...Input
integer, intent(in)::ielem
real*8,  intent(in)::xcrho,ycrho
real*8,dimension(1:2, 1:nvqua), intent(in)::xpq_fe
real*8,dimension(1:2, 1:nvqua), intent(in)::xpqi_fe
real*8,dimension(1:ndegr, 1:nsubq), intent(out)::unksgq
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc,ip,ifa,ifv
integer :: nvfsg, nvqsg

integer,dimension(4, nsubq)::ipqsg
!...Local real array
real*8,dimension(1:ndegr,1:nsubq)::rhsel
real*8,dimension(1:nmatr, 1:nsubq)::matin
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:4)    ::shpql, dsprql, dspsql
real*8,dimension(1:2)    ::shpfl, dsprfl
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: xrqc(2, 16),xrq(2, nvqua)
real*8:: posif(nvfac)
real*8::bq(ndegr)
!...Local real
real*8::r, s, rc,sc,dr,ds,xg,yg
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8::c10
real*8::f0, wt
real*8::masel,xgaus,ygaus

!
integer, allocatable::ipfsg(:, :)
real*8,allocatable::pofsg(:)
real*8,allocatable::mminv(:,:),binv(:),mmatr(:,:)
real*8,allocatable::xrq_sg(:, :, :),xpq_sg(:,:,:),xpqi_sg(:,:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
!...Part 0: Preliminary setup...
!
if(npoly==2) allocate(mminv(5,5), mmatr(5,5), binv(5))
!
if(ncurv.eq.1)then
!

elseif(ncurv.eq.2)then
!...Parmeters
nvfsg = 3
nvqsg = 9
!
allocate(pofsg(1:nvfsg))
allocate(xrq_sg(1:2, 1:nvqsg, 1:nsubq))
allocate(xpq_sg(1:2, 1:nvqsg, 1:nsubq))
allocate(xpqi_sg(1:2, 1:nvqsg, 1:nsubq))
allocate(ipfsg(nvfsg, 4))
!...quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) =13; ipqsg(4, 1) =12
ipqsg(1, 2) = 5; ipqsg(2, 2) = 9; ipqsg(3, 2) =14; ipqsg(4, 2) =13
ipqsg(1, 3) = 9; ipqsg(2, 3) = 2; ipqsg(3, 3) = 6; ipqsg(4, 3) =14
ipqsg(1, 4) =14; ipqsg(2, 4) = 6; ipqsg(3, 4) =10; ipqsg(4, 4) =15
ipqsg(1, 5) =15; ipqsg(2, 5) =10; ipqsg(3, 5) = 3; ipqsg(4, 5) = 7
ipqsg(1, 6) =16; ipqsg(2, 6) =15; ipqsg(3, 6) = 7; ipqsg(4, 6) =11
ipqsg(1, 7) = 8; ipqsg(2, 7) =16; ipqsg(3, 7) =11; ipqsg(4, 7) = 4
ipqsg(1, 8) =12; ipqsg(2, 8) =13; ipqsg(3, 8) =16; ipqsg(4, 8) = 8

!
if(nvfsg==3)then
pofsg(1) = -1.d0
pofsg(2) = 1.d0
pofsg(3) = 0.d0
!...face for subcell
ipfsg(1, 1) = 1; ipfsg(2, 1) = 2; ipfsg(3, 1) = 5;
ipfsg(1, 2) = 2; ipfsg(2, 2) = 3; ipfsg(3, 2) = 6;
ipfsg(1, 3) = 3; ipfsg(2, 3) = 4; ipfsg(3, 3) = 7;
ipfsg(1, 4) = 4; ipfsg(2, 4) = 1; ipfsg(3, 4) = 8;
elseif(nvfsg==4)then
pofsg(1) = -1.d0
pofsg(2) = 1.d0
pofsg(3) =-1.d0/3.d0
pofsg(4) = 1.d0/3.d0
!...face for subcell
ipfsg(1, 1) = 1; ipfsg(2, 1) = 2; ipfsg(3, 1) = 5; ipfsg(4, 1) = 9;
ipfsg(1, 2) = 2; ipfsg(2, 2) = 3; ipfsg(3, 2) = 6; ipfsg(4, 2) =10;
ipfsg(1, 3) = 3; ipfsg(2, 3) = 4; ipfsg(3, 3) = 7; ipfsg(4, 3) =11;
ipfsg(1, 4) = 4; ipfsg(2, 4) = 1; ipfsg(3, 4) = 8; ipfsg(4, 4) =12;
endif

!
posif(1) = -1.d0;
posif(2) = -sqrt(5.d0)/5.d0; posif(3) = sqrt(5.d0)/5.d0;
posif(4) = 1.d0;
!...Vertex coordinate
xrq(1, 1) = -1.d0;    xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0;    xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0;    xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0;    xrq(2, 4) =  1.d0;
xrq(1, 5) = posif(2); xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;    xrq(2, 6) = posif(2);
xrq(1, 7) = posif(3); xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;    xrq(2, 8) = posif(3);
xrq(1, 9) = posif(3); xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;    xrq(2,10) = posif(3);
xrq(1,11) = posif(2); xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;    xrq(2,12) = posif(2);
!
xrqc(:, 1:12) = xrq(:, 1:12)
xrqc(1, 13) = posif(2); xrqc(2, 13) = posif(2);
xrqc(1, 14) = posif(3); xrqc(2, 14) = posif(2);
xrqc(1, 15) = posif(3); xrqc(2, 15) = posif(3);
xrqc(1, 16) = posif(2); xrqc(2, 16) = posif(3);

endif

!...Get the Gauss
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0
!
!...Part 0: The reference coordinates
!
do isg = 1, nsubq
!...4 quad vertices
xrq_sg(1:2, 1:4 ,isg) = xrqc(1:2, ipqsg(1:4, isg))

do ifa = 1, 4
do ifv = 1, nvfsg
r =  pofsg(ifv)
call getshapfct_edge(0, 2 ,shpfl, dsprfl, r)
xrq_sg(1:2, ipfsg(ifv, ifa), isg) = shpfl(1)*xrq_sg(1:2, ipfsg(1, ifa), isg) +&
shpfl(2)*xrq_sg(1:2, ipfsg(2, ifa), isg)
enddo
enddo

enddo
!
!...Part I: Get cell averaged r and s for one subgrid...
!
do isg = 1, nsubq

!...Physical coordinates
do ip = 1, 8
!
r = xrq_sg(1, ip, isg)
s = xrq_sg(2, ip, isg)
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xg = 0.d0
yg = 0.d0
!
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpq_fe(1, ishp)
yg = yg + shpq(ishp)*xpq_fe(2, ishp)
enddo
!
xpq_sg(1, ip, isg) = xg
xpq_sg(2, ip, isg) = yg
enddo

!...the 9th nodes
xpq_sg(1:2, 9, isg) = -0.25d0*(xpq_sg(1:2, 1, isg) + xpq_sg(1:2, 2, isg) + xpq_sg(1:2, 3, isg) + xpq_sg(1:2, 4, isg)) +&
0.5d0*(xpq_sg(1:2, 5, isg) + xpq_sg(1:2, 6, isg) + xpq_sg(1:2, 7, isg) + xpq_sg(1:2, 8, isg))

!...Initialze parameters...
if(npoly.eq.2)  mmatr = 0.d0
f0 = 0.d0
wt = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(1,9,shpq, dsprq, dspsq, r, s)

dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 9
dxdr = dxdr + dsprq(ishp)*xpq_sg(1,ishp,isg)
dxds = dxds + dspsq(ishp)*xpq_sg(1,ishp,isg)

dydr = dydr + dsprq(ishp)*xpq_sg(2,ishp,isg)
dyds = dyds + dspsq(ishp)*xpq_sg(2,ishp,isg)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco
wt = wt +wi
!
!if(ielem.eq.1) print*,'djaco',isg,ig, wi,dxdr*dyds - dydr*dxds,xpq_sg
enddo
!if(ielem.eq.1) print*,'volume',f0,isg,wt
!
!stop
!
if(npoly==1)then
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then
matin(16,isg) = 1.d0/f0
endif

enddo
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, nsubq

!...Physical coordinates
do ip = 1, 8
!
r = xrq_sg(1, ip, isg)
s = xrq_sg(2, ip, isg)
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xg = 0.d0
yg = 0.d0
!
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpqi_fe(1, ishp)
yg = yg + shpq(ishp)*xpqi_fe(2, ishp)
enddo
!
xpqi_sg(1, ip, isg) = xg
xpqi_sg(2, ip, isg) = yg
enddo
!...the 9th nodes
xpqi_sg(1:2, 9, isg) = -0.25d0*(xpqi_sg(1:2, 1, isg) + xpqi_sg(1:2, 2, isg) + xpqi_sg(1:2, 3, isg) + xpqi_sg(1:2, 4, isg)) +&
0.5d0*(xpqi_sg(1:2, 5, isg) + xpqi_sg(1:2, 6, isg) + xpqi_sg(1:2, 7, isg) + xpqi_sg(1:2, 8, isg))
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinate
call getshapfct_quad(1,9,shpq, dsprq, dspsq, r, s)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 9
dxdr = dxdr + dsprq(ishp)*xpqi_sg(1,ishp,isg)
dxds = dxds + dspsq(ishp)*xpqi_sg(1,ishp,isg)

dydr = dydr + dsprq(ishp)*xpqi_sg(2,ishp,isg)
dyds = dyds + dspsq(ishp)*xpqi_sg(2,ishp,isg)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
bq(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, 9
xgaus = xgaus + shpq(ishp)*xpqi_sg(1,ishp,isg)
ygaus = ygaus + shpq(ishp)*xpqi_sg(2,ishp,isg)
enddo

!...Get the initial density
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*bq(1)*djaco
enddo
!
if(npoly==1)then
m(1,1) = matin(4, isg)
elseif(npoly==2)then
m(1,1) = matin(16, isg)
endif

!...Update the subcell density ditribution
unint = 0.d0
unint(1) = unint(1) + m(1, 1)*rhsel(1,isg)

!if(ielem.eq.1) print*,'isg4',isg,m(1, 1),rhsel(1,isg),m(1, 1)*rhsel(1,isg)
unksgq(1, isg)= unint(1)
enddo
!
end subroutine  getrhosubcell_sms_general5
!
!...Gte pressure distribution...
!
subroutine getunpre(bface, intfac, ipqua, esqua, unkno, unpre, geoel, coord, cooro)
!
use constant
implicit none
!
integer,dimension(1:nbfai,nbfac), intent(in)::bface
integer,dimension(1:nifai,1:nafac), intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndegr,1:nsize),      intent(inout)::unpre
real*8,dimension(1:ndimn,1:npoin),             intent(in) :: coord, cooro
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
!
!...  local arrays
!
integer::iel, ier
integer:: ipq(nvqua)
integer:: mapfe(1:2,1:nfqua)
!
real*8,dimension(1:2, 1:nvqua)::xpq, xpqi

real*8,dimension(3, nsize)::dmets
real*8,dimension(1:2,1:nsize) ::geoph
real*8,dimension(1:4) ::rhsmt, lhsmt, mtinv
real*8,dimension(1:2)::rhs
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:ndegr)::bq,bqp
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: weigh, wi

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dr, ds, rc, sc, r, s
real*8:: xmc, ymc
real*8:: masel
real*8:: xcrho, ycrho
real*8:: xgaus, ygaus, xgausi, ygausi
real*8:: dx, dy
real*8::dxwei,dywei
real*8::du1we
real*8::rxx,ryy,rxy
real*8::rco1,rco2,rco3,rco4
real*8::ru1x,ru1y
real*8::rjac1
real*8:: detma
real*8:: djaco,dxdr,dxds,dydr,dyds
real*8:: djacoi,dxdri,dxdsi,dydri,dydsi
real*8:: rhogi
real*8:: rhoc, rhom, uctr, vctr, ectr, pctr

!
integer:: ie, ies, iq, is, isten, ntrsf, ielem, jelem, ifa, igaus, ishp,ifprj,neaj,ieaj,istor,iv
integer:: nweno
real*8::  rpowe
!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /
!
!...Part 0: Basis parameters setup
!
!...Initial pressure
unpre = 0.d0

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...Scaling parameters
dr   = 1.d0
ds   = 1.d0

!...Mapping array
mapfe(1, 1) = 1; mapfe(2, 1) = 2;
mapfe(1, 2) = 2; mapfe(2, 2) = 3;
mapfe(1, 3) = 3; mapfe(2, 3) = 4;
mapfe(1, 4) = 4; mapfe(2, 4) = 1;
!
!...Part I: Get the L2 projection matrix
!

!...I.1:Get the physical mass center xmc, ymc
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

xpqi(1, 1:nvqua) = cooro(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = cooro(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

!...mass
masel = geoel(4, ielem)

!...physical mass center
xmc = 0.d0
ymc = 0.d0
!
do igaus =1,ngausdq
!
r  = posiq(1,igaus)
s  = posiq(2,igaus)
wi  = weighq(igaus)

!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!...Current domain
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)

!...Initial domain
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, nvqua
dxdri = dxdri + dsprq(ishp)*xpqi(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqi(1,ishp)

dydri = dydri + dsprq(ishp)*xpqi(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqi(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqi(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqi(2,ishp)

xgaus  = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus  = ygaus + shpq(ishp)*xpq(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...physical mass center
xmc = xmc + rhogi*djacoi*xgaus
ymc = ymc + rhogi*djacoi*ygaus
enddo

!...Physical mass center
geoph(1, ielem) = xmc/masel
geoph(2, ielem) = ymc/masel
enddo !do ie = 1, nquad

!
!...LS reconstruction
!
dmets = 0.d0
unpre = 0.d0

!...Give the cell-averaged density, veloccity, pressure....
do ie = 1, nquad
!
ielem=ie+ntria
!
rhom = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoc = 1.d0/rhom
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
unpre(1, ielem) = pctr;
enddo
!
do ifa = 1, nbfac

!... end-elements of the faces
iel   = intfac(1,ifa)
ier   = intfac(2,ifa)

if(bface(3, ifa).eq.21)then
cycle
elseif(bface(3, ifa).eq.22)then
if(bface(4, ifa).eq.222)then
geoph(1, ier) =-geoph(1, iel)
geoph(2, ier) = geoph(2, iel)
elseif(bface(4, ifa).eq.221)then
geoph(1, ier) = geoph(1, iel)
geoph(2, ier) =-geoph(2, iel)
endif
endif
!
unpre(1, ier) = unpre(1, iel)
enddo

!...Start LS
do 1400 ifa =1, nafac
!
if(ifa.le.nbfac)then
if(bface(3, ifa).eq.21) cycle
endif

!... end-elements of the faces
iel   = intfac(1,ifa)
ier   = intfac(2,ifa)
!
dx  = geoph(1,ier) - geoph(1,iel)
dy  = geoph(2,ier) - geoph(2,iel)

!...  weighting for this edge
weigh = 1.0
!
dxwei = weigh*dx
dywei = weigh*dy
!
dmets(1,iel) = dmets(1,iel) + dxwei*dxwei
dmets(2,iel) = dmets(2,iel) + dywei*dywei
dmets(3,iel) = dmets(3,iel) + dxwei*dywei
!
dmets(1,ier) = dmets(1,ier) + dxwei*dxwei
dmets(2,ier) = dmets(2,ier) + dywei*dywei
dmets(3,ier) = dmets(3,ier) + dxwei*dywei
!
du1we = weigh*(unpre(1, ier) - unpre(1, iel))
!
unpre(2,iel) = unpre(2,iel) + dxwei*du1we
unpre(3,iel) = unpre(3,iel) + dywei*du1we

!  excluding the ghost element gradient
if(ifa.gt.nbfac) then
unpre(2,ier) = unpre(2,ier) + dxwei*du1we
unpre(3,ier) = unpre(3,ier) + dywei*du1we
!
endif
!
1400 continue

!...
do 2000 ie  = 1, nquad
!
rxx         = dmets(1,ie)
ryy         = dmets(2,ie)
rxy         = dmets(3,ie)
!
rco1  = ryy
rco2  = -rxy
!
rco3  = -rxy
rco4  = rxx
!
rjac1       = 1.0/(rxx*ryy - rxy*rxy)
!
ru1x        = unpre(2,ie)
ru1y        = unpre(3,ie)
!
unpre(2,ie) = (ru1x*rco1 + ru1y*rco2 )*rjac1
unpre(3,ie) = (ru1x*rco3 + ru1y*rco4 )*rjac1
!
2000 continue

!...I.2: L2 projection matrix

rhsmt= 0.d0
lhsmt=0.d0
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

xpqi(1, 1:nvqua) = cooro(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = cooro(2, ipq(1:nvqua))
!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...physical mass center
xmc =geoph(1, ielem)
ymc =geoph(2, ielem)

!...scale
dx = maxval(xpq(1, 1:nvqua)) - minval(xpq(1, 1:nvqua))
dy = maxval(xpq(2, 1:nvqua)) - minval(xpq(2, 1:nvqua))

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

!...mass
masel = geoel(4, ielem)

!
do igaus =1,ngausdq
!
r  = posiq(1,igaus)
s  = posiq(2,igaus)
wi  = weighq(igaus)

!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!...Initial domain
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, nvqua
dxdri = dxdri + dsprq(ishp)*xpqi(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqi(1,ishp)

dydri = dydri + dsprq(ishp)*xpqi(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqi(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqi(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqi(2,ishp)

xgaus  = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus  = ygaus + shpq(ishp)*xpq(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bq(1) = 1.d0
bq(2) = (r-rc)/dr
bq(3) = (s-sc)/ds

bqp(1) = 1.d0
bqp(2) = (xgaus-xmc)/dx
bqp(3) = (ygaus-ymc)/dy

!...Matrix
rhsmt(1) = rhsmt(1) + rhogi*djacoi*bq(2)*bqp(2)
rhsmt(2) = rhsmt(2) + rhogi*djacoi*bq(3)*bqp(2)
rhsmt(3) = rhsmt(3) + rhogi*djacoi*bq(2)*bqp(3)
rhsmt(4) = rhsmt(4) + rhogi*djacoi*bq(3)*bqp(3)
!
lhsmt(1) = lhsmt(1)  + rhogi*djacoi*bq(2)*bq(2)
lhsmt(2) = lhsmt(2)  + rhogi*djacoi*bq(2)*bq(3)
lhsmt(3) = lhsmt(3)  + rhogi*djacoi*bq(2)*bq(3)
lhsmt(4) = lhsmt(4)  + rhogi*djacoi*bq(3)*bq(3)
!
enddo

!...Physical unkno
detma = lhsmt(1)*lhsmt(4)-lhsmt(2)*lhsmt(3)
mtinv(1) = lhsmt(4)
mtinv(2) =-lhsmt(2)
mtinv(3) =-lhsmt(3)
mtinv(4) = lhsmt(1)
!
mtinv = mtinv/detma
!
rhs(1) = unpre(2,ielem)*rhsmt(1) + unpre(3,ielem)*rhsmt(2)
rhs(2) = unpre(2,ielem)*rhsmt(3) + unpre(3,ielem)*rhsmt(4)

!...Unkowns on physical domain
unpre(2, ielem) = mtinv(1)*rhs(1) + mtinv(2)*rhs(2)
unpre(3, ielem) = mtinv(3)*rhs(1) + mtinv(4)*rhs(2)

enddo !do ie = 1, nquad

return
end subroutine getunpre
!
!...high order subroutine for barth limiter on curved quads with symmetry preserving....
!
subroutine barthlimitpres_lag_curvquadcb_HO(geoel, coord, coold, ustar, unpre, ipqua, &
bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nsize),             intent(inout)::unpre
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord, coold
real*8,dimension(1:nq+1, 1:nsize),             intent(inout)::aflim
integer,  dimension(1:nvqua,1:nquad),          intent(in):: ipqua
integer*4,dimension(1:nbfai,nbfac),           intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
real*8, dimension(1:nq+2, 1:npoin),           intent(in):: unmax, unmin
real*8,dimension(1:2, 1:2, 1:nsize),          intent(inout)::afvec
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
!
!...Local
!
integer:: ipq(nvqua)
integer:: indpt(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa
integer:: ielem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(1:ndegr, 1:nvqua),bqv(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xpqi
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy,rcv,scv
real*8:: rhov, rhoct, rhom
real*8:: afbar
real*8:: dudr, duds, dvdr, dvds
real*8:: r, s, rhoi, rhon
!
!...Part I: Some preparing work
!

!...Coloring the boundary node
indpt = 0  !...indpt represents index of boundary node
do ifa =1 ,nbfac
indpt(intfac(3:(2+nvfac), ifa)) = 1
enddo
!
eps = 1.e-6

!...shape functions
dr = 1.d0
ds = 1.d0
!
if(ncurv.eq.1)then
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
!
xvq(9) =  0.d0; yvq(9) =  0.d0
!
elseif(ncurv.eq.2)then
!
xvq(1) = -1.d0;            yvq(1) = -1.d0
xvq(2) =  1.d0;            yvq(2) = -1.d0
xvq(3) =  1.d0;            yvq(3) =  1.d0
xvq(4) = -1.d0;            yvq(4) =  1.d0
xvq(5) = -sqrt(5.d0)/5.d0; yvq(5) = -1.d0
xvq(6) =  1.d0;            yvq(6) =  -sqrt(5.d0)/5.d0
xvq(7) =  sqrt(5.d0)/5.d0; yvq(7) =  1.d0
xvq(8) = -1.d0;            yvq(8) =  sqrt(5.d0)/5.d0
!
xvq(9) =  sqrt(5.d0)/5.d0;  yvq(9)  = -1.d0
xvq(10) =  1.d0;            yvq(10) =  sqrt(5.d0)/5.d0
xvq(11) = -sqrt(5.d0)/5.d0; yvq(11) =  1.d0
xvq(12) = -1.d0;            yvq(12) = -sqrt(5.d0)/5.d0

endif

!
!...Part II: Impose limiter
!
do ie = 1, nquad

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Skip the smooth cell
if(geoel(10, ielem).lt.10.d0) cycle
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo

!...zero out unknv
unknvq = 0.d0
do iv   = 1,nvqua
do ideg = 1,3 !mdegr
unknvq(4, iv) = unknvq(4, iv) + unpre(ideg,ielem)*bq(ideg, iv)
enddo
enddo !do iv   = 1,nvqua

unctr(nq)  = unpre(1,ielem)

do iv = 1, nvqua
do iq =nq,nq
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)

alfa(iq, iv) = afbar
enddo
!
enddo ! do iv = 1, nvqua

!...Get the minimum value for one cell
do iq = nq,nq
if(geoel(10, ielem).gt.10.d0)then
aflim(iq, ielem) = minval(alfa(iq, 1:4))
else
aflim(iq, ielem) = 1.d0
endif
enddo
!
enddo


!
!...Part 4: Annihilate the HO moments of DGp2
!
do ie = 1,nquad

ielem = ie + ntria

!...Bad cell
if(geoel(10, ielem).gt.10.d0)then
unpre(4:6, ielem) = 0.d0
endif

enddo

end subroutine barthlimitpres_lag_curvquadcb_HO
