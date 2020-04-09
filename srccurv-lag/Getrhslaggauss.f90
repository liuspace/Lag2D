!
!...subroutine: Calculate the F^* N ds for all face gauss points for hybrid grids...
!
subroutine getfnds_lag_gaussh(gflag,gelag,gelagq,intfac,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(inout)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad),      intent(inout)::gelagq
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
real*8::posit(2, ngelg)
real*8::posiq(2, ngelgq)
real*8::shp(nptri),dspr(nptri),dsps(nptri)
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
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
if(ngausf.eq.3)then
!...Triangle
posit(1, 1) = 0.d0; posit(2, 1) = 0.d0;
posit(1, 2) = 1.d0; posit(2, 2) = 0.d0;
posit(1, 3) = 1.d0; posit(2, 3) = 0.d0;
posit(1, 4) = 0.d0; posit(2, 4) = 1.d0;
posit(1, 5) = 0.d0; posit(2, 5) = 1.d0;
posit(1, 6) = 0.d0; posit(2, 6) = 0.d0;
posit(1, 7) = .5d0; posit(2, 7) = 0.d0;
posit(1, 8) = .5d0; posit(2, 8) = .5d0;
posit(1, 9) = 0.d0; posit(2, 9) = .5d0;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;
!
elseif(ngausf.eq.4)then
!
!...Add triangle in the future...
!
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.4472135954999579392818d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.4472135954999579392818d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.4472135954999579392818d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.4472135954999579392818d0;

posiq(1,13) = 0.4472135954999579392818d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.4472135954999579392818d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.4472135954999579392818d0;
posiq(1,16) =-1.d0; posiq(2,16) =-0.4472135954999579392818d0;

endif
!
do 100 ie=1, ntria!...(1)ifa=1,nafac
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
!...coordinates
!
if(ncurv==0)then
xpt(1, 1:3) = coord(1, iptri(1:3,ie))
xpt(2, 1:3) = coord(2, iptri(1:3,ie))
!
xpt(1:2,4) = 0.5d0*(xpt(1:2,1) + xpt(1:2,2))
xpt(1:2,5) = 0.5d0*(xpt(1:2,2) + xpt(1:2,3))
xpt(1:2,6) = 0.5d0*(xpt(1:2,1) + xpt(1:2,3))
elseif(ncurv==1)then
xpt(1, 1:nptri) = coord(1,iptri(1:nptri, ie))
xpt(2, 1:nptri) = coord(2,iptri(1:nptri, ie))
endif
!
do ig = 1, ngelg!...(2)ig = 1,ngausd
!
r  = posit(1,ig)
s  = posit(2,ig)
!
!...Shape functions and their derivatives...
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
do ishp = 1, nptri
dxdr = dxdr + dspr(ishp)*xpt(1,ishp)
dxds = dxds + dsps(ishp)*xpt(1,ishp)

dydr = dydr + dspr(ishp)*xpt(2,ishp)
dyds = dyds + dsps(ishp)*xpt(2,ishp)
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
enddo
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa
!
!...Identify the local No. of one internal face for left cell...
!
if(ig.eq.1.or.ig.eq.2.or.ig.eq.7)then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 1.d0;
!
idfal = ig;
!
elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.8)then
!
vnorm(1) = sqrt(2.d0)*0.5d0; vnorm(2) = sqrt(2.d0)*0.5d0;
larea    = sqrt(2.d0)
!
idfal = ig;
elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.9)then
!
vnorm(1) = -1.d0;           vnorm(2) = 0.d0;
larea    = 1.d0
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
gelag(1, idfal, ie) = anx*larea/farea
gelag(2, idfal, ie) = any*larea/farea
gelag(3, idfal, ie) = farea
!
enddo !...ig from 1...9
!
100 enddo  !...(1)ifa=1,nafac
!
! print*,'vnotmfn',gelag(1:3, 1:6, 170)
!
!
do 200 ie=1, nquad !...(1)ifa=1,nquad
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...coordinates
!
if(ncurv==0)then
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif
!
do ig = 1, ngelgq!...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa
!
!...Identify the local No. of one internal face for left cell...
!
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
enddo !...ig from 1...9
!
200 enddo  !...(1)ifa=1,nquad
!
!print*,'Inside getfnds_lag',gelagq(1:3,2,2434)
!
end subroutine getfnds_lag_gaussh
!
!...Get the nodal velocity based on Lobatto Gauss integration...
!
subroutine getndvelo_lag_gaussh(gflag,gelag,gelagq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, ufgaus, fstart, fstarq, aflim, afvec, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3, 1:ngelg, 1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:ngelg, 1:ntria),   intent(out)::fstart !...Riemann forces
real*8,dimension(1:ndimn,1:ngelgq, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngelgq,1:nquad),       intent(out)::ufgaus

integer:: itime,ip
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)
integer, dimension(3, 4)::fglvq
integer, dimension(3, 3)::fglvt
!
integer, dimension(ngausf, 4):: fglgq
integer, dimension(ngausf, 3):: fglgt

!...local real array
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:ndimn,1:nvqua,1:nquad)::vnulq
real*8,  dimension(1:nquad)::gqdmp
real*8::munaci(2, 2)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2
real*8::r,shp1,shp2,shp3
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: munacn(:,:,:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: munacu(:,:), snsigm(:,:)
real*8,allocatable:: munaclt(:,:,:,:),munaclq(:,:,:,:)
real*8,allocatable:: snsigmlt(:,:,:), munault(:,:,:)
real*8,allocatable:: snsigmlq(:,:,:), munaulq(:,:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:2, 1:2, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munaclt(1:2, 1:2, 1:ngelg, 1:ntria), munault(1:ndimn, 1:ngelg,  1:ntria),&
snsigmlt(1:ndimn, 1:ngelg,  1:ntria))
allocate (munaclq(1:2, 1:2, 1:ngelgq, 1:nquad), munaulq(1:ndimn, 1:ngelgq,  1:nquad),&
snsigmlq(1:ndimn, 1:ngelgq,  1:nquad))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
!
!...Part I: Specify some gauss points...
!

!...Local vertex No. of gauss points in one unit quad...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;

!...Local vertex No. of gauss points in one unit tria...
fglvt(1, 1) = 1;  fglvt(2, 1) = 2; fglvt(3, 1) = 4;
fglvt(1, 2) = 2;  fglvt(2, 2) = 3; fglvt(3, 2) = 5;
fglvt(1, 3) = 3;  fglvt(2, 3) = 1; fglvt(3, 3) = 6;

!...Local gauss point No. of any gauss point in one face...
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;
!
fglgt(1, 1) = 1;  fglgt(2, 1) = 2; fglgt(3, 1) = 7;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4; fglgt(3, 2) = 8;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6; fglgt(3, 3) = 9;
!
elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
!
endif
!
!...Zero out vlave
!
vlave = 0.d0
indnd = 0
!
!...Mark the boundary nodes...
!
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

!...Get averaged velocity at nodes...
!call getvlavenew(iptri, ipqua, geoel, vlave, unkno, aflim, afvec)
!
do 950 ifa = 1 , nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
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
endif!
950 enddo
!
!...Part II: Loop to get information from Riemann solver
!
do iloop= 1, 1!5

!...Give value
vlave= ustar
vnulq = 0.d0

!...Zero out munacn
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0

!...Triangle in the future...

!...Quad
if(nquad.gt.0) call getriem_quad_gauss(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)

!....Update the velocity at the end points...
do ifa = 1, nafac
do ip = 1, 2

 ipoin = intfac(2+ip, ifa)
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

!if(ipoin.eq.1921)print*,'ipnode',ustar(1:2,ipoin)!,detma,munacn(:,:,ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin)
endif
enddo
enddo

!...Get the velocity at the middle point...
call getvelo_mpt_vtx(ustar,intfac,ipqua, geoel, vlave,coord,unkno,indnd, afvec, aflim)

!...Get the vertex velocity at the boundary....
do 900 ifa = 1 , nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
!if(bface(3,ifa).eq.22)then
!ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!if(bface(4,ifa).eq.221)then
!ustar(2,ipf(1:nvfac)) = 0.d0
!elseif(bface(4,ifa).eq.222)then
!ustar(1,ipf(1:nvfac)) = 0.d0
!endif
!elseif(bface(3,ifa).eq.25)then
!ustar(2,ipf(1:nvfac)) = 0.d0
!ustar(1,ipf(1:nvfac)) = 0.d0
!endif
!
!ustar(1, ipf(3)) = 0.5d0*(ustar(1, ipf(1))+ustar(1, ipf(2)))
!ustar(2, ipf(3)) = 0.5d0*(ustar(2, ipf(1))+ustar(2, ipf(2)))
!
if(ncase.eq.1)then
ustar(1, ipf(1:nvfac)) = sin(pi*coord(1,ipf(1:nvfac)))*cos(pi*coord(2,ipf(1:nvfac)))
ustar(2, ipf(1:nvfac)) =-cos(pi*coord(1,ipf(1:nvfac)))*sin(pi*coord(2,ipf(1:nvfac)))
endif
900 enddo

!...Give ufgaus initially...
if(iloop.eq.1)then
!
do 300 ie=1, nquad !...(1)ifa=1,nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
do ifa = 1,4
!
ufgaus(1:ndimn,fglgq(1, ifa),  ie) = ustar(1:2, ipq(fglvq(1, ifa)))
ufgaus(1:ndimn,fglgq(2, ifa),  ie) = ustar(1:2, ipq(fglvq(2, ifa)))
!
r = -0.4472135954999579392818d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
!
ufgaus(1:ndimn,fglgq(3, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
shp3*ustar(1:2, ipq(fglvq(3, ifa)))
!
r =  0.4472135954999579392818d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
!
ufgaus(1:ndimn,fglgq(4, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
shp3*ustar(1:2, ipq(fglvq(3, ifa)))
!
enddo
300 enddo

endif
!
call getvelo_mpt_gauss2(ufgaus,gelagq,intfac,ipqua,coord,unkno,indnd, &
munaclq, munaulq, snsigmlq,afvec, aflim)

enddo
!
!...Update ufgaus using array ustar ...
!
do 100 ie=1, nquad !...(1)ifa=1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
do ifa = 1,4
!
!ufgaus(1:ndimn,fglgq(1, ifa),  ie) = ustar(1:2, ipq(fglvq(1, ifa)))
!ufgaus(1:ndimn,fglgq(2, ifa),  ie) = ustar(1:2, ipq(fglvq(2, ifa)))
!
r = -0.4472135954999579392818d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
!
if(ufgaus(1,fglgq(3, ifa),  ie).gt.1.d5)then
  ufgaus(1:ndimn,fglgq(3, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
                                        shp3*ustar(1:2, ipq(fglvq(3, ifa)))
endif
r =  0.4472135954999579392818d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
!
if(ufgaus(1,fglgq(4, ifa),  ie).gt.1.d5)then
  ufgaus(1:ndimn,fglgq(4, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
                                       shp3*ustar(1:2, ipq(fglvq(3, ifa)))
endif
!
!ustar(1:2, ipq(fglvq(3, ifa))) = 0.625d0*(ufgaus(1:ndimn,fglgq(3, ifa),  ie)+ufgaus(1:ndimn,fglgq(4, ifa),  ie))-&
!                                 0.125d0*(ufgaus(1:ndimn,fglgq(1, ifa),  ie)+ufgaus(1:ndimn,fglgq(2, ifa),  ie))
enddo
!!
!if(ie.eq.1830) print*,'ie1830',ufgaus(1:2,7,ie),fglgq(1, 4),ipq(fglvq(1, 4))
!
100 enddo

!
!...4.2: Update the Riemann forces at every node...
!...Tria
!
!...Quad
!
do ie = 1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!!
do ifa = 1, 4
do ig =1, ngausf
!
!...Basis function!
!
fstarq(1, fglgq(ig, ifa), ie) = snsigmlq(1, fglgq(ig, ifa), ie) +&
munaclq(1,1, fglgq(ig, ifa), ie)*ufgaus(1, fglgq(ig, ifa), ie)+&
munaclq(2,1, fglgq(ig, ifa), ie)*ufgaus(2, fglgq(ig, ifa), ie)-&
munaulq(1, fglgq(ig, ifa), ie)
!
!if(ie.eq.40) print*,'fstarq',ie,ifa,ig,snsigmlq(1, fglgq(ig, ifa), ie),munaclq(1,1, fglgq(ig, ifa), ie),&
!munaulq(1, fglgq(ig, ifa), ie),fstarq(1, fglgq(ig, ifa), ie),ufgaus(1:2, fglgq(ig, ifa), ie)
!
fstarq(2, fglgq(ig, ifa), ie) = snsigmlq(2, fglgq(ig, ifa), ie) +&
munaclq(2,2,fglgq(ig, ifa), ie)*ufgaus(2, fglgq(ig, ifa), ie)+&
munaclq(1,2,fglgq(ig, ifa), ie)*ufgaus(1, fglgq(ig, ifa), ie)-&
munaulq(2, fglgq(ig, ifa), ie)
!
enddo
enddo
!
!if(ie.eq.379) print*,'ie379',fstarq(1:2, fglgq(1, 3), ie)
enddo
!
deallocate (munacn, bpres, fpres)
deallocate (munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lag_gaussh
!
!...Riemann input for curved quad using gauss rule...
!
subroutine getriem_quad_gauss(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),           intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),         intent(in)::afvec
!
real*8, dimension(1:2, 1:2, 1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
!
real*8, dimension(1:2, 1:2, 1:ngelgq, 1:nquad),       intent(out)::munaclq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(out)::snsigmlq
!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer, dimension(nvfac, 4)::fglvq
integer, dimension(ngausf, 4):: fglgq
!...local real array
real*8,dimension(1:ndegr)::bg, bgv
real*8,dimension(1:nq,1:ngausf)::unkng
real*8::aujmp(1:3, 1:ngausf)
real*8::vnorm(1:3, 1:ngelgq)
real*8::posiq(1:2, 1:ngelgq)
real*8::sigmg(1:2, 1:2, 1:ngausf)
real*8,dimension(1:ngausf)::murie
real*8,dimension(1:nvqua):: xvq,  yvq
!...arraies for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhog,rhomg,ug,vg,eg, pg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, rg, sg,rcv,scv
real*8:: dudr, duds, dvdr, dvds
real*8::acnx, acny
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Local vertex No. of gauss points in one unit ...
!
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
!...Local gauss point No. of any gauss point in one face...
!
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;
!
elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.4472135954999579392818d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.4472135954999579392818d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.4472135954999579392818d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.4472135954999579392818d0;

posiq(1,13) = 0.4472135954999579392818d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.4472135954999579392818d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.4472135954999579392818d0;
posiq(1,16) =-1.d0; posiq(2,16) =-0.4472135954999579392818d0;

endif
!
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
dr = 1.d0
ds = 1.d0
!
!...Triangle...
!
do 250 ie = 1,nquad !...(1)ie = 1,nelem
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
!...shape functions
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
!...Give the normal vector of every face...
!
vnorm(1:3,  1:ngelgq) = gelagq(1:3, 1:ngelgq, ie);
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...Give the normal vector of every face...
!
do ifa =1, 4
!
!...zero out unkno
!
unkng = 0.d0
!
do ig   = 1,  2!ngausf
!
rg = posiq(1, fglgq(ig, ifa))
sg = posiq(2, fglgq(ig, ifa))
!
bg(1) = 1.d0
bg(2) = (rg-rc)/dr
bg(3) = (sg-sc)/ds
!
do ideg = 1,mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq,ie)*bg(ideg)
enddo
!
!if(ip(fagsp(ig, ifa)).eq.155) print*,'ndegr',unkno(1:3,2,ie),bg(1:3),ie,ifa,rg,sg
!
if(ndens.eq.1)then
rhog  = 1.d0/unkng(1, ig)
elseif(ndens.eq.3)then
!
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
bgv(1) = 1.d0
bgv(2) = (rg-rcv)/dr
bgv(3) = (sg-scv)/ds
!
unkng(1, ig) =0.d0
!
do ideg = 1,mdegr
unkng(1, ig) = unkng(1, ig) + unkno(ideg,1,ielem)*bgv(ideg)
enddo
!
rhog = unkng(1, ig)
endif
!
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
pg = max(eps, (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2)))
!
!if(ie.eq.1829)print*,ie,ig,ifa,pg,rhog,ug,vg,unkno(1:3,2,ie)
!
!...Limiter
!
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomg = rhomc + aflim(1, ielem)*(unkng(1, ig) - rhomc)
rhog = 1.d0/rhomg
elseif(ndens.eq.3)then
rhog = 1.d0/rhomc + aflim(1, ielem)*(rhog - 1.d0/rhomc)
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
ug = unkno(1,2,ielem)  + dudr*bg(2) + duds*bg(3)
vg = unkno(1,3,ielem)  + dvdr*bg(2) + dvds*bg(3)
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pg = pctr + aflim(4, ielem)*(pg - pctr)
!!
!...updtae unknv(2:3,:)
unkng(2, ig) = ug
unkng(3 ,ig) = vg
!
endif
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig) = -pg
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pg
!
!...Get the a_c (unit vector)
!
aujmp(1:2, ig) = vlave(1:2, ipq(fglvq(ig, ifa))) - unkng(2:3, ig)
!
acnx = aujmp(1, ig)
acny = aujmp(2, ig)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig) = 1.e-11!0.d0;
else
aujmp(1:2, ig) = aujmp(1:2, ig)/sqrt(acnx**2 + acny**2)
endif
!
aujmp(3, ig) = sqrt(acnx**2 + acny**2)
!
!if(ie.eq.1829)print*,ie,ig,ifa,pg,pctr,aflim(4, ielem)
!
enddo !....ig
!
!...Get the variables at the center...
!
!rhoct = 1.d0/unkno(1, 1, ielem)         !...ct denots center of one cell; cn denotes corner of one cell.
!uctr  = unkno(1, 2, ielem)
!vctr  = unkno(1, 3, ielem)
!ectr  = unkno(1, 4, ielem)
!pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
aujmp(3,:)=aujmp(3,:)/sdctr
!
!...Get impedence coefficient...
!
do ig   = 1, 2
dux= vlave(1, ipq(fglvq(ig, ifa)))-unkng(2, ig)
duy= vlave(2, ipq(fglvq(ig, ifa)))-unkng(3, ig)
!deltu = sqrt(dux**2 + duy**2)
deltu = abs(dux*vnorm(1,fglgq(ig, ifa)) + duy*vnorm(2,fglgq(ig, ifa)))
murie(ig) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!if(ipq(fglvq(ig, ifa)).eq.62) print*,'murie22', sdctr,rhoct,vlave(1:2, ip(fagsp(ig, ifa))),unkng(2:3, ig),unkno(1,2,ie),ie
enddo
!
!if(ie==1831) print*,'vnotm',unkno(1,1:4,ie)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do ig  = 1,  2!ngausf
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
call getriecoef_matrixnew(murie(ig), vnorm(3, fglgq(ig, ifa)), vnorm(1:2, fglgq(ig, ifa)), aujmp(1:3, ig), &
unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ig), vnorm(3, fglgq(ig, ifa)), vnorm(1:2, fglgq(ig, ifa)), aujmp(1:2, ig), &
!unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipq(fglvq(ig, ifa))) = munacn(1:2, 1, ipq(fglvq(ig, ifa))) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipq(fglvq(ig, ifa))) = munacn(1:2, 2, ipq(fglvq(ig, ifa))) + munacn_rie(1:2, 2)
!
munacu(1:2, ipq(fglvq(ig, ifa))) = munacu(1:2, ipq(fglvq(ig, ifa))) + munacu_rie(1:2)
!
snsigm(1:2,ipq(fglvq(ig, ifa))) = snsigm(1:2, ipq(fglvq(ig, ifa))) + snsigm_rie(1:2)!
!
munaclq(1:2, 1, fglgq(ig, ifa), ie) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fglgq(ig, ifa), ie) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fglgq(ig, ifa), ie)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fglgq(ig, ifa), ie)    =  snsigm_rie(1:2)
!
!if(ipq(fglvq(ig, ifa)).eq.1921) print*,'ep', murie(ig),ie,ifa,ig,fglvq(ig, ifa),sigmg(1:2, 1:2, ig),vnorm(1:3, fglgq(ig, ifa))
!
enddo
!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,nelem!
!
!print*,'ep1862',munacn(1:2, 1:2, 1862),munacu(1:2, 1862), snsigm(1:2, 1862)
!
end subroutine getriem_quad_gauss
!
!...Idensity the smae gaus point
!
subroutine getfgauss(ipqua, ipf, fgaus)
use constant
implicit none
integer, dimension(nvqua), intent(in):: ipqua
integer, dimension(nvfac), intent(in):: ipf
integer:: ig, nsum
integer, dimension(1:ngausf), intent(out)::fgaus
!
integer::fglgq(ngausf,4)
!
!...Specify fglgq
!
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

elseif(ngausf.eq.4)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
endif
!
!...Specify the No. of face gauss points in fglgq...
!
nsum = ngausf +3
!
do ig = 3, ngausf
!
if(ipqua(5).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(1))then
fgaus(ig) = fglgq(ig,1)
else
fgaus(ig) = fglgq(nsum-ig,1)
endif
!
elseif(ipqua(6).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(2))then
fgaus(ig) = fglgq(ig,2)
else
fgaus(ig) = fglgq(nsum-ig,2)
endif
!
elseif(ipqua(7).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(3))then
fgaus(ig) = fglgq(ig,3)
else
fgaus(ig) = fglgq(nsum-ig,3)
endif
!
elseif(ipqua(8).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(4))then
fgaus(ig) = fglgq(ig,4)
else
fgaus(ig) = fglgq(nsum-ig,4)
endif
!
endif

enddo
end subroutine getfgauss
!
!...Idensity the vertex point
!
subroutine getfvtx(ipqua, ipf, fvtx)
use constant
implicit none
integer, dimension(nvqua), intent(in):: ipqua
integer, dimension(nvfac), intent(in):: ipf
integer:: ig, nsum
integer, dimension(1:nvfac), intent(out)::fvtx
!
integer::fglvq(nvfac,4)
!
!...Specify fglgq
!
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
!...Specify the No. of face gauss points in fglgq...
!
nsum = ngausf +3
!
do ig = 3, nvfac
!
if(ipqua(5).eq.ipf(3))then
fvtx(ig) = fglvq(ig,1)
elseif(ipqua(6).eq.ipf(3))then
fvtx(ig) = fglvq(ig,2)
elseif(ipqua(7).eq.ipf(3))then
fvtx(ig) = fglvq(ig,3)
elseif(ipqua(8).eq.ipf(3))then
fvtx(ig) = fglvq(ig,4)
endif

enddo
end subroutine getfvtx
!
!...Calculate the velocity at the middle point for Marie method...
!
subroutine getvelo_mpt_vtx(ustar,intfac,ipqua,geoel, vlave,coord,unkno,indnd, afvec, aflim)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:2, 1:npoin),              intent(inout)::ustar
integer*4,dimension(1:npoin),                intent(in)::indnd
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic,ishp
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(nvfac, 4)::fglvq
integer,dimension(nvfac)::fgausl, fgausr
!
real*8::eps
real*8::comatr(2,2)
real*8,dimension(1:3,1:nvqua,1:nquad)::vnormq
real*8::unkng(1:nq+3, 1:nvfac, 1:2)
real*8::sigmg(1:2, 1:2, 1:nvfac, 1:2)
real*8::vngs(3, nvfac, 2)
real*8::aujmp(1:3,1:nvfac,1:2)
real*8::vnorm(1:2)
real*8::xvq(nvqua), yvq(nvqua), b(ndegr, nvfac)
real*8::xpq(2, npqua)
real*8::xpf(1:2, 1:nvfac)
real*8,dimension(npqua)::shpq,dsprq,dspsq
!
!...Riemann variables...
!
real*8::munacn_rie(2,2), munacu_rie(2), snsigm_rie(2)
real*8::munacn(1:2, 1:2, 1:nvfac), munacu(1:2, 1:nvfac), snsigm(1:2, 1:nvfac)
real*8::munaci(2, 2)
!
!...Local real number
!
real*8::dxdr,dxds,dydr,dyds
real*8::detma,dudr,duds,dvdr,dvds
real*8::pvtx,rhovt, rhomc, rhomv, rhovl, rhovr,rhsu1,rhsu2
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
real*8::othog
real*8::acnx,acny,shp1,shp2,shp3
real*8::c10,r,s
real*8::anx,any,farea,larea,delu,muasp

!...Constant definition...
 eps = 1.d-6
 c10 = 1.d0

!...dr,ds,rc,zc
dr = 1.d0
ds = 1.d0
rc = 0.0
sc = rc
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
xvq(9) =  0.d0; yvq(9) =  0.d0
!
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
!...Part I:...Get the normal vector at the gauss...
!
do 200 ie=1, nquad !...(1)ifa=1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...coordinates
if(ncurv==0)then
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif

!
do ig = 5, 8!...(2)ig = 1,ngausd
!
r  = xvq(ig)
s  = yvq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa
!
!...Identify the local No. of one internal face for left cell...
!
if(ig.eq.5)then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 2.d0;
!
idfal = ig;
!
elseif(ig.eq.6)then
!
vnorm(1) = 1.d0; vnorm(2) = 0.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.7)then
!
vnorm(1) = 0.d0;           vnorm(2) = 1.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.8)then
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
!
farea    = sqrt(vnorm(1)**2 + vnorm(2)**2) !...farea: face area...
!
vnormq(1, ig, ie) = anx*larea/farea
vnormq(2, ig, ie) = any*larea/farea
vnormq(3, ig, ie) = farea
!
!if(ipq(ig).eq.7473) print*,'node',vnormq(1:3, ig, ie),ie
!
enddo
200 enddo !...ig from 1...9
!
!...Part II: Get the middle point velocity...
!
fgausl = 0.d0
fgausr = 0.d0
!
do 450 ifa = 1, nafac !...(1)ie = 1,nelem
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))

!...For the linear PP+
ftx = xpf(1 ,2)- xpf(1, 1)
fty = xpf(2, 2)- xpf(2, 1)
!
fnx = -fty/sqrt(ftx**2 + fty**2)
fny =  ftx/sqrt(ftx**2 + fty**2)

!...For linear PM
ftx = xpf(1 ,3)- xpf(1, 1)
fty = xpf(2, 3)- xpf(2, 1)
!
lenmc = sqrt(ftx**2 + fty**2)
!
ftx = ftx/lenmc
fty = fty/lenmc
!
!...The orthogonality criterion....
!
othog = abs(fnx*ftx + fny*fty)
!
if(ifa.le.nbfac)then
!
!...Calculate the velocity of boundary face point...
!
elseif(ifa.gt.nbfac)then

!...xxx: Left cell
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...

!...Identify the fgaus for left cell
call getfvtx(ipqua(1:nvqua, iel), ipf, fgausl)
!
unkng = 0.d0

!...Basis fctn...
do ig =3 ,nvfac
b(1, ig) = 1.d0
b(2, ig) = (xvq(fgausl(ig))-rc)/dr
b(3, ig) = (yvq(fgausl(ig))-sc)/ds
enddo

!...Unkno at gauss points...
do ig = 3, nvfac
do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*b(ideg, ig)
enddo
vngs(1:3, ig, 1) = vnormq(1:3, fgausl(ig), iel) !...face normal vector
enddo

!...Limiter...
do ig = 3, nvfac
!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, iel)*(unkng(1, ig, 1) - rhomc)
unkng(1, ig, 1) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, iel)*unkno(2,2,iel) +  afvec(1, 2, iel)*unkno(2,3,iel)
duds = afvec(1, 1, iel)*unkno(3,2,iel) +  afvec(1, 2, iel)*unkno(3,3,iel)
dvdr = afvec(2, 1, iel)*unkno(2,2,iel) +  afvec(2, 2, iel)*unkno(2,3,iel)
dvds = afvec(2, 1, iel)*unkno(3,2,iel) +  afvec(2, 2, iel)*unkno(3,3,iel)
!
uvtx = unkno(1,2,iel)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,iel)  + dvdr*b(2, ig) + dvds*b(3, ig)
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pvtx = pctr + aflim(4, iel)*(pvtx - pctr)
!...updtae unknv(2:3,:)
unkng(2, ig, 1) = uvtx
unkng(3 ,ig, 1) = vvtx
!
endif
!
unkng(5, ig, 1) = pvtx
unkng(6, ig, 1) = rhoct*sdctr
unkng(7, ig, 1) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig, 1) = -pvtx
sigmg(1, 2, ig, 1) = 0.d0
sigmg(2, 1, ig, 1) = 0.d0
sigmg(2, 2, ig, 1) = -pvtx!
!
enddo

!...xxx: Right cell
rhoct = 1.d0/unkno(1, 1, ier)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ier)
vctr  = unkno(1, 3, ier)
ectr  = unkno(1, 4, ier)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center..

!...Identify the fgaus for right cell
call getfvtx(ipqua(1:nvqua, ier), ipf, fgausr)
!
do ig =3 ,nvfac
b(1, ig) = 1.d0
b(2, ig) = (xvq(fgausr(ig))-rc)/dr
b(3, ig) = (yvq(fgausr(ig))-sc)/ds
enddo
!
do ig = 3, nvfac
do ideg =1, mdegr
unkng(1:nq, ig, 2) = unkng(1:nq, ig, 2) + unkno(ideg,1:nq,ier)*b(ideg, ig)
enddo
vngs(1:3, ig, 2) = vnormq(1:3, fgausr(ig), ier) !...face normal vector
enddo

!...Limiter...
do ig = 3, nvfac
!
rhovt = 1.d0/unkng(1, ig, 2)
uvtx = unkng(2, ig, 2)
vvtx = unkng(3, ig, 2)
evtx = unkng(4, ig, 2)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, ier)*(unkng(1, ig, 2) - rhomc)
unkng(1, ig, 2) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, ier)*unkno(2,2,ier) +  afvec(1, 2, ier)*unkno(2,3,ier)
duds = afvec(1, 1, ier)*unkno(3,2,ier) +  afvec(1, 2, ier)*unkno(3,3,ier)
dvdr = afvec(2, 1, ier)*unkno(2,2,ier) +  afvec(2, 2, ier)*unkno(2,3,ier)
dvds = afvec(2, 1, ier)*unkno(3,2,ier) +  afvec(2, 2, ier)*unkno(3,3,ier)
!
uvtx = unkno(1,2,ier)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,ier)  + dvdr*b(2, ig) + dvds*b(3, ig)

!if(iel.eq.2581) print*,'ie26252628',ier,pctr,aflim(4, ier),pvtx
pvtx = pctr + aflim(4, ier)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig, 2) = uvtx
unkng(3, ig, 2) = vvtx
!
endif
!
unkng(5, ig, 2) = pvtx
unkng(6, ig, 2) = rhoct*sdctr
unkng(7, ig, 2) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig, 2) = -pvtx
sigmg(1, 2, ig, 2) = 0.d0
sigmg(2, 1, ig, 2) = 0.d0
sigmg(2, 2, ig, 2) = -pvtx!

!if(iel.eq.469) print*,'NAn',vlave(1:2, fgausr(ig), ier) , unkng(2:3, ig, 2)
enddo

!...Get the averaged velocity at the gauss quadrature points...
do ig = 3, nvfac
do ic = 1, 2
!
aujmp(1:2, ig, ic) = vlave(1:2, ipf(3)) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic) 
!
!...Impedence
!
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(1, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) !+ unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = 3, nvfac !...ig =3, nvfac
do ic = 1, 2
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
!f(ipf(3).eq.7594)print*,ic,murie(ic), ic,iel,ier,unkno(1, 1:4,ier)
!
enddo
!
!...Get the middle point velocity...
!
if(indnd(ipf(3)).eq.0)then
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig) !- fpres(1, ipoin)
rhsu2 = munacu(2, ig) - snsigm(2, ig) !- fpres(2, ipoin)
!
ustar(1, ipf(3) ) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ustar(2, ipf(3) ) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!
!if(ier.eq.2433) print*,'iel469umpt',fgausl(ig),ig,munacu(1, ig), snsigm(1, ig),sigmg(1:2, 1:2, ig, 2)
!
endif
enddo !...ig =3, nvfac
!
!...The special treatment for munacn = 0.0
!
do ig = 3, nvfac
!
muasp = unkng(6, ig, 1)*unkng(7, ig, 1) + unkng(6, ig, 2)*unkng(7, ig, 2)
!
!if(othog.lt.1.d-4)then
!print*,ifa,munacn(1, 1, ig),unkng(7, ig, 1),unkng(7, ig, 2)
if(munacn(1, 1, ig)/muasp.lt.1.d-6)then
!
rhovl = 1.d0/unkng(1, ig, 1)
uvtxl = unkng(2, ig, 1)
vvtxl = unkng(3, ig, 1)
evtxl = unkng(4, ig, 1)
pvtxl = unkng(5, ig, 1)
mufal = unkng(6, ig, 1)
!
rhovr = 1.d0/unkng(1, ig, 2)
uvtxr = unkng(2, ig, 2)
vvtxr = unkng(3, ig, 2)
evtxr = unkng(4, ig, 2)
pvtxr = unkng(5, ig, 2)
mufar = unkng(6, ig, 2)
!
fnx = vngs(1, ig, 1) !...face normal vector for iel...
fny = vngs(2, ig, 1)
!
! print*,'fnx',vngs(1:2, ig, 1)
!
ftx = -fny
fty =  fnx
!...
if(indnd(ipf(3)).eq.0)then
!
ustar(1, ipf(3) )  = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (pvtxr- pvtxl)/(mufal+mufar)*fnx
ustar(2, ipf(3) )  = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (pvtxr- pvtxl)/(mufal+mufar)*fny
!
!if(ipf(3).eq.7353) print*,'mpt1862',ustar(1:2, ipf(3) ) ,mufal,uvtxl , mufar,uvtxr , (pvtxr- pvtxl),fnx, fny
!
endif
!
endif
enddo
!
endif
!
450 enddo
!
end subroutine getvelo_mpt_vtx
!
!...Calculate the velocity at the middle point for Marie method...
!
subroutine getvelo_mpt_gauss(ustar,ufgaus,gelagq,intfac,ipqua,coord,unkno,indnd, &
munaclq, munaulq, snsigmlq, afvec, aflim)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:2, 1:npoin),              intent(inout)::ustar
integer*4,dimension(1:npoin),                intent(in)::indnd
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngelgq,nquad),        intent(out)::ufgaus
real*8, dimension(1:2, 1:2, 1:ngelgq, 1:nquad), intent(inout)::munaclq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(inout)::munaulq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(inout)::snsigmlq
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(ngausf)::fgausl, fgausr
!
real*8::eps
real*8,dimension(1:2,1:ngelgq,nquad)::vlave
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::sigmg(1:2, 1:2, 1:ngausf, 1:2)
real*8::vngs(3, ngausf,2)
real*8::aujmp(1:3,1:ngausf,1:2)
real*8::murie(2)
real*8::vnorm(1:3,1:nvqua)
real*8::xvq(nvqua), yvq(nvqua), b(ndegr, ngausf)
real*8::posiq(1:2, 1:ngelgq)
integer,dimension(nvfac, 4)::fglvq
integer,dimension(ngausf,4)::fglgq
real*8::xpf(1:2, 1:nvfac)
!
!...Riemann parameters...
!
real*8::munacn_rie(2,2), munacu_rie(2), snsigm_rie(2)
real*8::munacn(1:2, 1:2, 1:ngausf), munacu(1:2, 1:ngausf), snsigm(1:2, 1:ngausf)
real*8::munaci(2, 2)
!
!...Local real number
!
real*8::detma,dudr,duds,dvdr,dvds
real*8::pvtx,rhovt, rhomc, rhomv, rhovl, rhovr,rhsu1,rhsu2
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
real*8::othog
real*8::acnx,acny,shp1,shp2,shp3,r,delu
!
!...Initial velocity at gauss point
!
ufgaus = 0.d0
!
eps = 1.d-6
!
!...Zero out plnpn, ulnpn
!
dr = 1.d0
ds = 1.d0
rc = 0.0
sc = rc
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
xvq(9) =  0.d0; yvq(9) =  0.d0
!
!...Local vertex No. of gauss points in one unit ...
!
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
if(ngausf.eq.3)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.4472135954999579392818d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.4472135954999579392818d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.4472135954999579392818d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.4472135954999579392818d0;

posiq(1,13) = 0.4472135954999579392818d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.4472135954999579392818d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.4472135954999579392818d0;
posiq(1,16) =-1.d0; posiq(2,16) =-0.4472135954999579392818d0;
endif
!
!...Part I: Get ustar for aujmp...
!
do 100 ie=1, -nquad !...(1)ifa=1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
do ifa = 1,4
!
vlave(1:ndimn,fglgq(1, ifa),  ie) = ustar(1:2, ipq(fglvq(1, ifa)))
vlave(1:ndimn,fglgq(2, ifa),  ie) = ustar(1:2, ipq(fglvq(2, ifa)))
!
r = -0.4472135954999579392818d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
!
vlave(1:ndimn,fglgq(3, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
shp3*ustar(1:2, ipq(fglvq(3, ifa)))
!
!if(ie.eq.469)print*,'vlave',ustar(1:2, ipq(fglvq(3, ifa))),fglgq(3, ifa)
!
r =  0.4472135954999579392818d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
vlave(1:ndimn,fglgq(4, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
shp3*ustar(1:2, ipq(fglvq(3, ifa)))
enddo
!
100 enddo
!
!...Part II: Get the middle point velocity...
!
fgausl = 0.d0
fgausr = 0.d0
!
do 450 ifa = 1, nafac !...(1)ie = 1,nelem
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))
!
!...For the linear PP+
ftx = xpf(1 ,2)- xpf(1, 1)
fty = xpf(2, 2)- xpf(2, 1)
!
fnx = -fty/sqrt(ftx**2 + fty**2)
fny =  ftx/sqrt(ftx**2 + fty**2)
!
!...For linear PM
ftx = xpf(1 ,3)- xpf(1, 1)
fty = xpf(2, 3)- xpf(2, 1)
!
lenmc = sqrt(ftx**2 + fty**2)
!
ftx = ftx/lenmc
fty = fty/lenmc
!
othog = abs(fnx*ftx + fny*fty)
!
if(ifa.le.nbfac)then
!
!...(ifa.gt.nbfac) For interior cells...
!

elseif(ifa.gt.nbfac)then
!
!...Parameters for Left cell
!
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...Identify the fgaus for left cell
!
call getfgauss(ipqua(1:nvqua, iel), ipf, fgausl)
!
unkng = 0.d0
!
!...Basis fctn...
!
do ig =3 ,ngausf
b(1, ig) = 1.d0
b(2, ig) = (posiq(1, fgausl(ig))-rc)/dr
b(3, ig) = (posiq(2, fgausl(ig))-sc)/ds
enddo
!
!...Unkno at gauss points...
!
do ig = 3, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*b(ideg, ig)
enddo
vngs(1:3, ig, 1) = gelagq(1:3, fgausl(ig), iel) !...face normal vector
enddo
!
!...Limiter...
!
do ig = 3, ngausf
!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, iel)*(unkng(1, ig, 1) - rhomc)
unkng(1, ig, 1) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, iel)*unkno(2,2,iel) +  afvec(1, 2, iel)*unkno(2,3,iel)
duds = afvec(1, 1, iel)*unkno(3,2,iel) +  afvec(1, 2, iel)*unkno(3,3,iel)
dvdr = afvec(2, 1, iel)*unkno(2,2,iel) +  afvec(2, 2, iel)*unkno(2,3,iel)
dvds = afvec(2, 1, iel)*unkno(3,2,iel) +  afvec(2, 2, iel)*unkno(3,3,iel)
!
uvtx = unkno(1,2,iel)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,iel)  + dvdr*b(2, ig) + dvds*b(3, ig)
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pvtx = pctr + aflim(4, iel)*(pvtx - pctr)
!!
!...updtae unknv(2:3,:)
unkng(2, ig, 1) = uvtx
unkng(3 ,ig, 1) = vvtx
!
endif
!
unkng(5, ig, 1) = pvtx
unkng(6, ig, 1) = rhoct*sdctr
unkng(7, ig, 1) = sdctr
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig, 1) = -pvtx
sigmg(1, 2, ig, 1) = 0.d0
sigmg(2, 1, ig, 1) = 0.d0
sigmg(2, 2, ig, 1) = -pvtx!
!
enddo
!
!...Right cell
!
rhoct = 1.d0/unkno(1, 1, ier)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ier)
vctr  = unkno(1, 3, ier)
ectr  = unkno(1, 4, ier)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center..
!
!...Identify the fgaus for right cell
!
call getfgauss(ipqua(1:nvqua, ier), ipf, fgausr)
!
!if(iel.eq.2433) print*,'fgaussr',ipqua(1:nvqua, iel),ipf(1:3),fgausl
!
do ig =3 ,ngausf
b(1, ig) = 1.d0
b(2, ig) = (posiq(1, fgausr(ig))-rc)/dr
b(3, ig) = (posiq(2, fgausr(ig))-sc)/ds
enddo
!
do ig = 3, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 2) = unkng(1:nq, ig, 2) + unkno(ideg,1:nq,ier)*b(ideg, ig)
enddo
vngs(1:3, ig, 2) = gelagq(1:3, fgausr(ig), ier) !...face normal vector
enddo
!
!...Limiter...
!
do ig = 3, ngausf
!
rhovt = 1.d0/unkng(1, ig, 2)
uvtx = unkng(2, ig, 2)
vvtx = unkng(3, ig, 2)
evtx = unkng(4, ig, 2)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, ier)*(unkng(1, ig, 2) - rhomc)
unkng(1, ig, 2) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, ier)*unkno(2,2,ier) +  afvec(1, 2, ier)*unkno(2,3,ier)
duds = afvec(1, 1, ier)*unkno(3,2,ier) +  afvec(1, 2, ier)*unkno(3,3,ier)
dvdr = afvec(2, 1, ier)*unkno(2,2,ier) +  afvec(2, 2, ier)*unkno(2,3,ier)
dvds = afvec(2, 1, ier)*unkno(3,2,ier) +  afvec(2, 2, ier)*unkno(3,3,ier)
!
uvtx = unkno(1,2,ier)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,ier)  + dvdr*b(2, ig) + dvds*b(3, ig)
!
!if(iel.eq.2581) print*,'ie26252628',ier,pctr,aflim(4, ier),pvtx
!
pvtx = pctr + aflim(4, ier)*(pvtx - pctr)
!!
!...updtae unknv(2:3,:)
unkng(2, ig, 2) = uvtx
unkng(3, ig, 2) = vvtx
!
endif
!
unkng(5, ig, 2) = pvtx
unkng(6, ig, 2) = rhoct*sdctr
unkng(7, ig, 2) = sdctr
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig, 2) = -pvtx
sigmg(1, 2, ig, 2) = 0.d0
sigmg(2, 1, ig, 2) = 0.d0
sigmg(2, 2, ig, 2) = -pvtx!
!
enddo
!
!...Get the averaged velocity at the gauss quadrature points...
!
do ig = 3, ngausf
!!
vlave(1, fgausl(ig), iel) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
vlave(2, fgausl(ig), iel) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
!!
vlave(1, fgausr(ig), ier) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
vlave(2, fgausr(ig), ier) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
enddo
!
!...Impedence...
!
do ig = 3, ngausf
do ic = 1, 2
!
aujmp(1:2, ig, ic) = vlave(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic) 
!
!...Impedence
!
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) !+ unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo
!
!...Summation over corners
!
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = 3, ngausf
do ic = 1, 2
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
!if(iel.eq.1892) print*,'mptl', iel,fgausl(ig), murie(ic), vngs(1:3, ig, ic),aujmp(1:3, ig, ic)
!if(ier.eq.1892) print*,'mptr', ic, ig,ier,fgausr(ig), murie(ic), vngs(1:3, ig, ic),aujmp(1:3, ig, ic)
!
if(ic.eq.1)then
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
!
!if(iel.eq.1892) print*,'mptr', ic, ig,iel, fgausl(ig),snsigmlq(1:2, fgausl(ig), iel)
! 
else
!
munaclq(1:2, 1, fgausr(ig), ier) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausr(ig), ier) =  munacn_rie(1:2, 2)
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaulq(1:2, fgausr(ig), ier)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausr(ig), ier)    = snsigm_rie(1:2)
!
!if(ier.eq.1892) print*,'mptr', ic, ig,ier, fgausr(ig),snsigmlq(1:2,  fgausr(ig), ier) 
endif
!
enddo
!
!
if(indnd(ipf(3)).eq.0)then
!
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig) !- fpres(1, ipoin)
rhsu2 = munacu(2, ig) - snsigm(2, ig) !- fpres(2, ipoin)
!
ufgaus(1, fgausl(ig), iel) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ufgaus(2, fgausl(ig), iel) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!
ufgaus(1, fgausr(ig), ier) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ufgaus(2, fgausr(ig), ier) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!
!if(ier.eq.2433) print*,'iel469umpt',fgausl(ig),ig,munacu(1, ig), snsigm(1, ig),sigmg(1:2, 1:2, ig, 2)
!
endif
!
enddo
!
!...The special treatment for munacn = 0.0
!
do ig = 3, ngausf
!
!if(othog.lt.1.d-1)then
!if(munacn(1, 1, ig).lt.1.d-3)then
!
rhovl = 1.d0/unkng(1, ig, 1)
uvtxl = unkng(2, ig, 1)
vvtxl = unkng(3, ig, 1)
evtxl = unkng(4, ig, 1)
pvtxl = unkng(5, ig, 1)
mufal = unkng(6, ig, 1)
!
rhovr = 1.d0/unkng(1, ig, 2)
uvtxr = unkng(2, ig, 2)
vvtxr = unkng(3, ig, 2)
evtxr = unkng(4, ig, 2)
pvtxr = unkng(5, ig, 2)
mufar = unkng(6, ig, 2)
!
fnx = vngs(1, ig, 1) !...face normal vector for iel...
fny = vngs(2, ig, 1)
!
ftx = -fny
fty =  fnx
!
!...Mar
!
if(indnd(ipf(3)).eq.0)then
ufgaus(1, fgausl(ig), iel) = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (pvtxr- pvtxl)/(mufal+mufar)*fnx
ufgaus(2, fgausl(ig), iel) = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (pvtxr- pvtxl)/(mufal+mufar)*fny
!
ufgaus(1, fgausr(ig), ier) = ufgaus(1, fgausl(ig), iel)
ufgaus(2, fgausr(ig), ier) = ufgaus(2, fgausl(ig), iel)
!
endif
!
!endif
enddo
endif
450 enddo
!
end subroutine getvelo_mpt_gauss
!
!...Calculate the velocity at the Gauss point for Marie method...
!
subroutine getvelo_mpt_gauss2(ufgaus,gelagq,intfac,ipqua,coord,unkno,indnd, &
munaclq, munaulq, snsigmlq, afvec, aflim)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
integer*4,dimension(1:npoin),                intent(in)::indnd
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngelgq,nquad),      intent(inout)::ufgaus
real*8, dimension(1:2, 1:2, 1:ngelgq, 1:nquad), intent(inout)::munaclq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(inout)::munaulq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(inout)::snsigmlq
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(ngausf)::fgausl, fgausr
!
real*8::eps
real*8,dimension(1:2,1:ngelgq,nquad)::vlave
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::sigmg(1:2, 1:2, 1:ngausf, 1:2)
real*8::vngs(3, ngausf,2)
real*8::aujmp(1:3,1:ngausf,1:2)
real*8::murie(2)
real*8::vnorm(1:3,1:nvqua)
real*8::xvq(nvqua), yvq(nvqua), b(ndegr, ngausf)
real*8::posiq(1:2, 1:ngelgq)
integer,dimension(nvfac, 4)::fglvq
integer,dimension(ngausf,4)::fglgq
real*8::xpf(1:2, 1:nvfac)

!...Riemann parameters...
real*8::munacn_rie(2,2), munacu_rie(2), snsigm_rie(2)
real*8::munacn(1:2, 1:2, 1:ngausf), munacu(1:2, 1:ngausf), snsigm(1:2, 1:ngausf)
real*8::munaci(2, 2)

!...Local real number
real*8::detma,dudr,duds,dvdr,dvds
real*8::pvtx,rhovt, rhomc, rhomv, rhovl, rhovr,rhsu1,rhsu2
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
real*8::othog
real*8::acnx,acny,shp1,shp2,shp3,r,delu

!...Initial velocity at gauss point
!..Initial Gauss point velocity from the interpolation of the 3-node velocity...
!
eps = 1.d-6
!
dr = 1.d0
ds = 1.d0
rc = 0.0
sc = rc
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
xvq(9) =  0.d0; yvq(9) =  0.d0

!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

!...8-nodes quad
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.4472135954999579392818d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.4472135954999579392818d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.4472135954999579392818d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.4472135954999579392818d0;

posiq(1,13) = 0.4472135954999579392818d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.4472135954999579392818d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.4472135954999579392818d0;
posiq(1,16) =-1.d0; posiq(2,16) =-0.4472135954999579392818d0;
endif
!
!...Part II: Get the middle point velocity...
!
fgausl = 0.d0
fgausr = 0.d0
!
do 450 ifa = 1, nafac !...(1)ie = 1,nelem
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))

!...For the linear PP+
ftx = xpf(1 ,2)- xpf(1, 1)
fty = xpf(2, 2)- xpf(2, 1)
!
fnx = -fty/sqrt(ftx**2 + fty**2)
fny =  ftx/sqrt(ftx**2 + fty**2)

!...For linear PM
ftx = xpf(1 ,3)- xpf(1, 1)
fty = xpf(2, 3)- xpf(2, 1)
!
lenmc = sqrt(ftx**2 + fty**2)
!
ftx = ftx/lenmc
fty = fty/lenmc
!
othog = abs(fnx*ftx + fny*fty)

!...Boundary face
if(ifa.le.nbfac)then

!...Parameters for Left cell
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...

!...Identify the fgaus for left cell
call getfgauss(ipqua(1:nvqua, iel), ipf, fgausl)
!
unkng = 0.d0

!...Basis fctn...
do ig =3 ,ngausf
b(1, ig) = 1.d0
b(2, ig) = (posiq(1, fgausl(ig))-rc)/dr
b(3, ig) = (posiq(2, fgausl(ig))-sc)/ds
enddo

!...Unkno at gauss points...
do ig = 3, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*b(ideg, ig)
enddo
vngs(1:3, ig, 1) = gelagq(1:3, fgausl(ig), iel) !...face normal vector
enddo

!...Limiter...
do ig = 3, ngausf
!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, iel)*(unkng(1, ig, 1) - rhomc)
unkng(1, ig, 1) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, iel)*unkno(2,2,iel) +  afvec(1, 2, iel)*unkno(2,3,iel)
duds = afvec(1, 1, iel)*unkno(3,2,iel) +  afvec(1, 2, iel)*unkno(3,3,iel)
dvdr = afvec(2, 1, iel)*unkno(2,2,iel) +  afvec(2, 2, iel)*unkno(2,3,iel)
dvds = afvec(2, 1, iel)*unkno(3,2,iel) +  afvec(2, 2, iel)*unkno(3,3,iel)
!
uvtx = unkno(1,2,iel)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,iel)  + dvdr*b(2, ig) + dvds*b(3, ig)

!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx

pvtx = pctr + aflim(4, iel)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig, 1) = uvtx
unkng(3 ,ig, 1) = vvtx
endif

unkng(5, ig, 1) = pvtx
unkng(6, ig, 1) = rhoct*sdctr
unkng(7, ig, 1) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig, 1) = -pvtx
sigmg(1, 2, ig, 1) = 0.d0
sigmg(2, 1, ig, 1) = 0.d0
sigmg(2, 2, ig, 1) = -pvtx!
enddo

!...Get the averaged velocity at the gauss quadrature points...
do ig = 3, ngausf
vlave(1, fgausl(ig), iel) = unkng(2, ig, 1)
vlave(2, fgausl(ig), iel) = unkng(3, ig, 1)
enddo

!...Impedence...
do ig = 3, ngausf
do ic = 1, 1
!
!aujmp(1:2, ig, ic) = vlave(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgaus(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)
!
!...Impedence
!
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) !+ unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = 3, ngausf
do ic = 1, 1
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
!
enddo
!

if(indnd(ipf(3)).eq.0)then
!
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig) !- fpres(1, ipoin)
rhsu2 = munacu(2, ig) - snsigm(2, ig) !- fpres(2, ipoin)
!
ufgaus(1, fgausl(ig), iel) = 1.d10!munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ufgaus(2, fgausl(ig), iel) = 1.d10!munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!
!if(ier.eq.2433) print*,'iel469umpt',fgausl(ig),ig,munacu(1, ig), snsigm(1, ig),sigmg(1:2, 1:2, ig, 2)
!
endif
enddo

!...Interior face
elseif(ifa.gt.nbfac)then

!...Parameters for Left cell
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...

!...Identify the fgaus for left cell
call getfgauss(ipqua(1:nvqua, iel), ipf, fgausl)
!
unkng = 0.d0

!...Basis fctn...
do ig =3 ,ngausf
b(1, ig) = 1.d0
b(2, ig) = (posiq(1, fgausl(ig))-rc)/dr
b(3, ig) = (posiq(2, fgausl(ig))-sc)/ds
enddo

!...Unkno at gauss points...
do ig = 3, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*b(ideg, ig)
enddo
vngs(1:3, ig, 1) = gelagq(1:3, fgausl(ig), iel) !...face normal vector
enddo

!...Limiter...
do ig = 3, ngausf
!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, iel)*(unkng(1, ig, 1) - rhomc)
unkng(1, ig, 1) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, iel)*unkno(2,2,iel) +  afvec(1, 2, iel)*unkno(2,3,iel)
duds = afvec(1, 1, iel)*unkno(3,2,iel) +  afvec(1, 2, iel)*unkno(3,3,iel)
dvdr = afvec(2, 1, iel)*unkno(2,2,iel) +  afvec(2, 2, iel)*unkno(2,3,iel)
dvds = afvec(2, 1, iel)*unkno(3,2,iel) +  afvec(2, 2, iel)*unkno(3,3,iel)
!
uvtx = unkno(1,2,iel)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,iel)  + dvdr*b(2, ig) + dvds*b(3, ig)
!
pvtx = pctr + aflim(4, iel)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig, 1) = uvtx
unkng(3 ,ig, 1) = vvtx
endif
!
unkng(5, ig, 1) = pvtx
unkng(6, ig, 1) = rhoct*sdctr
unkng(7, ig, 1) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig, 1) = -pvtx
sigmg(1, 2, ig, 1) = 0.d0
sigmg(2, 1, ig, 1) = 0.d0
sigmg(2, 2, ig, 1) = -pvtx!
enddo

!...Right cell
rhoct = 1.d0/unkno(1, 1, ier)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ier)
vctr  = unkno(1, 3, ier)
ectr  = unkno(1, 4, ier)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center..

!...Identify the fgaus for right cell
call getfgauss(ipqua(1:nvqua, ier), ipf, fgausr)
!
do ig =3 ,ngausf
b(1, ig) = 1.d0
b(2, ig) = (posiq(1, fgausr(ig))-rc)/dr
b(3, ig) = (posiq(2, fgausr(ig))-sc)/ds
enddo
!
do ig = 3, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 2) = unkng(1:nq, ig, 2) + unkno(ideg,1:nq,ier)*b(ideg, ig)
enddo
vngs(1:3, ig, 2) = gelagq(1:3, fgausr(ig), ier) !...face normal vector
enddo

!...Limiter...
do ig = 3, ngausf
!
rhovt = 1.d0/unkng(1, ig, 2)
uvtx = unkng(2, ig, 2)
vvtx = unkng(3, ig, 2)
evtx = unkng(4, ig, 2)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, ier)*(unkng(1, ig, 2) - rhomc)
unkng(1, ig, 2) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, ier)*unkno(2,2,ier) +  afvec(1, 2, ier)*unkno(2,3,ier)
duds = afvec(1, 1, ier)*unkno(3,2,ier) +  afvec(1, 2, ier)*unkno(3,3,ier)
dvdr = afvec(2, 1, ier)*unkno(2,2,ier) +  afvec(2, 2, ier)*unkno(2,3,ier)
dvds = afvec(2, 1, ier)*unkno(3,2,ier) +  afvec(2, 2, ier)*unkno(3,3,ier)
!
uvtx = unkno(1,2,ier)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,ier)  + dvdr*b(2, ig) + dvds*b(3, ig)
!
!if(iel.eq.2581) print*,'ie26252628',ier,pctr,aflim(4, ier),pvtx
!
pvtx = pctr + aflim(4, ier)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig, 2) = uvtx
unkng(3, ig, 2) = vvtx
!
endif
!
unkng(5, ig, 2) = pvtx
unkng(6, ig, 2) = rhoct*sdctr
unkng(7, ig, 2) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig, 2) = -pvtx
sigmg(1, 2, ig, 2) = 0.d0
sigmg(2, 1, ig, 2) = 0.d0
sigmg(2, 2, ig, 2) = -pvtx!
enddo

!...Get the averaged velocity at the gauss quadrature points...
do ig = 3, ngausf
vlave(1, fgausl(ig), iel) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
vlave(2, fgausl(ig), iel) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
!!
vlave(1, fgausr(ig), ier) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
vlave(2, fgausr(ig), ier) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
enddo
!
!...Impedence...
!
do ig = 3, ngausf
do ic = 1, 2
!
!aujmp(1:2, ig, ic) = vlave(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgaus(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)

!...Impedence
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) !+ unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = 3, ngausf
do ic = 1, 2
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)

!
if(ic.eq.1)then
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)

!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
else
!
munaclq(1:2, 1, fgausr(ig), ier) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausr(ig), ier) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausr(ig), ier)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausr(ig), ier)    = snsigm_rie(1:2)
!
!if(ier.eq.1892) print*,'mptr', ic, ig,ier, fgausr(ig),snsigmlq(1:2,  fgausr(ig), ier)
endif
!
enddo
!
if(indnd(ipf(3)).eq.0)then
!
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig) !- fpres(1, ipoin)
rhsu2 = munacu(2, ig) - snsigm(2, ig) !- fpres(2, ipoin)
!
ufgaus(1, fgausl(ig), iel) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ufgaus(2, fgausl(ig), iel) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!
ufgaus(1, fgausr(ig), ier) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ufgaus(2, fgausr(ig), ier) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!
endif
!
enddo
!
!...The special treatment for munacn = 0.0
!
do ig = 3, ngausf
!
!if(othog.lt.1.d-1)then
!if(munacn(1, 1, ig).lt.1.d-3)then
!
rhovl = 1.d0/unkng(1, ig, 1)
uvtxl = unkng(2, ig, 1)
vvtxl = unkng(3, ig, 1)
evtxl = unkng(4, ig, 1)
pvtxl = unkng(5, ig, 1)
mufal = unkng(6, ig, 1)
!
rhovr = 1.d0/unkng(1, ig, 2)
uvtxr = unkng(2, ig, 2)
vvtxr = unkng(3, ig, 2)
evtxr = unkng(4, ig, 2)
pvtxr = unkng(5, ig, 2)
mufar = unkng(6, ig, 2)
!
fnx = vngs(1, ig, 1) !...face normal vector for iel...
fny = vngs(2, ig, 1)
!
ftx = -fny
fty =  fnx
!
!...Mar
!
if(indnd(ipf(3)).eq.0)then
ufgaus(1, fgausl(ig), iel) = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (pvtxr- pvtxl)/(mufal+mufar)*fnx
ufgaus(2, fgausl(ig), iel) = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (pvtxr- pvtxl)/(mufal+mufar)*fny
!
ufgaus(1, fgausr(ig), ier) = ufgaus(1, fgausl(ig), iel)
ufgaus(2, fgausr(ig), ier) = ufgaus(2, fgausl(ig), iel)
!
endif
!
!endif
enddo
endif
450 enddo
!
end subroutine getvelo_mpt_gauss2
!
!...Rhsiface for gauss
!
subroutine rhsifacedg_lagquadgauss(ipqua,  unkno, ustar, ufgaus, geoel, gelagq, fstarq, coord, coold,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar, coord, coold!...nodal velocity
real*8,dimension(1:ndimn,1:ngelgq, 1:nquad), intent(in)  ::ufgaus
real*8,dimension(1:ndimn,1:ngelgq,1:nquad),  intent(in) ::fstarq !...Riemann forces


real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::rhsel
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ielem
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer, dimension(nvfac, 4)::fglvq
integer,dimension(ngausf, 4)::fglgq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::posiq(1:2, 1:ngelgq)
real*8::xvq(nvqua), yvq(nvqua)
real*8::bg(3, ngausf)
!
real*8::weigh(ngausf), posi(1,ngausf)
!...local real number
real*8::eps,c00,c05,c10,c20,c13
real*8::wi
real*8:: nx, ny ,sa
real*8:: rg, sg
real*8::dr,ds,rc,sc
real*8::dwav1,dwav2
real*8::anx, any
real*8::gpnx, gpny, gpsa
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
!
!...Gauss quadrature...
!
!print*,'Inside getfnds_lag',gelagq(1:3,2,2434)
!
!...Local vertex No. of gauss points in one unit ...
!
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
if(ngausf.eq.3)then
!
posi(1, 1) = -1.d0; posi(1 ,2)= 1.d0;
posi(1 ,3)= 0.d0;

!
weigh(  1) = 1.d0/6.d0
weigh(  3) = 4.d0/6.d0
weigh(  2) = 1.d0/6.d0
!
!...Local gauss point No. of any gauss point in one face...
!
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;
elseif(ngausf.eq.4)then

!
posi(1, 1) = -1.d0; posi(1 ,2)= 1.d0;
posi(1 ,3)= -0.4472135954999579392818d0;
posi(1 ,4)=  0.4472135954999579392818d0;
!
weigh(  1) = 1.d0/12.d0
weigh(  2) = 1.d0/12.d0
weigh(  3) = 5.d0/12.d0
weigh(  4) = 5.d0/12.d0
!
!...Local gauss point No. of any gauss point in one face...
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.4472135954999579392818d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.4472135954999579392818d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.4472135954999579392818d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.4472135954999579392818d0;

posiq(1,13) = 0.4472135954999579392818d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.4472135954999579392818d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.4472135954999579392818d0;
posiq(1,16) =-1.d0; posiq(2,16) =-0.4472135954999579392818d0;
endif
!
!...Zero out plnpn, ulnpn
!
dr = 1.d0
ds = 1.d0
!
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
!print*,'Inside getfnds_lag',gelagq(1:3,2,2434)
!
do 550 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
ipq(1:nvqua) = ipqua(1:nvqua ,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do ifa = 1, 4
!
do ig   = 1, ngausf
!
wi  = weigh(ig)
!
rg = posiq(1,fglgq(ig, ifa))
sg = posiq(2,fglgq(ig, ifa))
!
bg(1,  ig) = 1.d0
bg(2,  ig) = (rg-rc)/dr
bg(3,  ig) = (sg-sc)/ds
!
!if(ie.eq.1) print*,'bg',ig,ifa,infag(ig, ifa), rg, sg,bg(1:3,3-ig, infag(ig, ifa))
!
!
!...The vertex constituting one cell...
!d0
gpnx = gelagq(1, fglgq(ig, ifa),ie)
gpny = gelagq(2, fglgq(ig, ifa),ie)
gpsa = gelagq(3, fglgq(ig, ifa),ie)
!
!...Distribute to every corner...
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ufgaus(1, fglgq(ig, ifa),ie)*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig) +&
ufgaus(2, fglgq(ig, ifa),ie)*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)
!
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
elnpn(1:ndegr)=elnpn(1:ndegr) +&
ufgaus(1, fglgq(ig, ifa), ie)*&
fstarq(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig) +&
ufgaus(2, fglgq(ig, ifa), ie)*&
fstarq(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
!if(ie.eq.22.or.ie==379) print*,'rhs iface idegr',ie,ifa,ig,fglgq(ig, ifa),elnpn(1),&
!ufgaus(1, fglgq(ig, ifa), ie)*fstarq(1, fglgq(ig, ifa), ie),&
!ufgaus(2, fglgq(ig, ifa), ie)*fstarq(2, fglgq(ig, ifa), ie),&
!ufgaus(1, fglgq(ig, ifa), ie)*fstarq(1, fglgq(ig, ifa), ie)+&
!ufgaus(2, fglgq(ig, ifa), ie)*fstarq(2, fglgq(ig, ifa), ie)
!ufgaus(2, fglgq(ig, ifa),ie)*gpny*gpsa*bg(1, ig)*weigh(ig),ufgaus(1:2, fglgq(ig, ifa),ie)
!if(ie==1831) print*,'rhs ifaceidegr2',ifa,ig,fglgq(ig, ifa),ulnpn(1:ndegr),ustar(1, ipq(fglvq(1, 2))),ipq(fglvq(1, 2))
!
enddo
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1829) print*,'rhs iface idegr 1',ie,rhsel(1:3, 2, ielem)
!
550 enddo
!
!   open(8,file='rhstxt.dat')
!    do ie = 1, nquad
!   write(8,*) ie, rhsel(1, 1:4, ie)
!    enddo
!   close(8)
!
end subroutine rhsifacedg_lagquadgauss
!
!...Get FN dS for sub grid method...
!
subroutine getfnds_lagsubg(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
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
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
real*8,dimension(1:ndimn, 1:nvfac)::xpf
!...local real array
real*8,dimension(1:ndimn, 1:nvtri)::xp
real*8,dimension(1:ndimn, 1:nvqua)::xpq
!...local real number
real*8::anx, any
real*8::dr, ds, rc, sc
real*8::c16
!
data c16   /0.1666666666666666d0 /
!
!...Boundary face...
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
!print*,ifa,iel,gflag(1:3, ifa)
!
enddo
!
!...Quads
!
do 55 ie=1, nquad!...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...coordinates
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Edge 1 (p15)
!
anx =   xpq(2, 5) - xpq(2, 1)
any = -(xpq(1, 5) - xpq(1, 1))
!
gesgq(1, 1 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 1 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 1 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 2 (p52)
!
anx =   xpq(2, 2) - xpq(2, 5)
any = -(xpq(1, 2) - xpq(1, 5))
!
gesgq(1, 2 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 2 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 2 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 3 (p26)
!
anx =   xpq(2, 6) - xpq(2, 2)
any = -(xpq(1, 6) - xpq(1, 2))
!
gesgq(1, 3 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 3 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 3 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 4 (p63)
!
anx =   xpq(2, 3) - xpq(2, 6)
any = -(xpq(1, 3) - xpq(1, 6))
!
gesgq(1, 4 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 4 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 4 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 5 (p37)
!
anx =   xpq(2, 7) - xpq(2, 3)
any = -(xpq(1, 7) - xpq(1, 3))
!
gesgq(1, 5 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 5 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 5 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 6 (p74)
!
anx =   xpq(2, 4) - xpq(2, 7)
any = -(xpq(1, 4) - xpq(1, 7))
!
gesgq(1, 6 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 6 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 6 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 7 (p48)
!
anx =   xpq(2, 8) - xpq(2, 4)
any = -(xpq(1, 8) - xpq(1, 4))
!
gesgq(1, 7 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 7 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 7 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 8 (p81)
!
anx =   xpq(2, 1) - xpq(2, 8)
any = -(xpq(1, 1) - xpq(1, 8))
!
gesgq(1, 8 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 8 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 8 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 9 (p95)
!
anx =   xpq(2, 5) - xpq(2, 9)
any = -(xpq(1, 5) - xpq(1, 9))
!
gesgq(1, 9 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 9 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 9 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 10 (p96)
!
anx =   xpq(2, 6) - xpq(2, 9)
any = -(xpq(1, 6) - xpq(1, 9))
!
gesgq(1, 10 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 10 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 10 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 11 (p97)
!
anx =   xpq(2, 7) - xpq(2, 9)
any = -(xpq(1, 7) - xpq(1, 9))
!
gesgq(1, 11 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 11 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 11 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 12 (p98)
!
anx =   xpq(2, 8) - xpq(2, 9)
any = -(xpq(1, 8) - xpq(1, 9))
!
gesgq(1, 12 , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 12 , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 12 , ie) = sqrt(anx**2 + any**2)
!
55 enddo  !...(1)ifa=1,nelem
!
! print*,'vnotmfn',gelag(1:3, 1:9, 50)
!
end subroutine getfnds_lagsubg
!
!...subroutine: Calculate the F^* N ds for interior faces for hybrid grids...
!
subroutine getfnds_lagsubgc(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
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
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(1:nvfac,1:4)::ipfsc
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8::dshpr(3),shpq(9)
!...local real array
real*8,dimension(1:ndimn, 1:nvtri)::xp
real*8,dimension(1:ndimn, 1:nvqua)::xpq
real*8::posqc(1:2,1:4),xpqv(1:2,1:4)
!...local real number
real*8::anx, any,xsc,ysc,dxdr,dydr
real*8::dr, ds, rc, sc,r,s,djaco
real*8::c16,c10
!
data c16   /0.1666666666666666d0 /
data c10   /1.d0 /
!
posqc(1, 1)=  0.0d0; posqc(2, 1)= -.5d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)=  0.d0;
posqc(1, 3)=  0.0d0; posqc(2, 3)=  .5d0;
posqc(1, 4)= -0.5d0; posqc(2, 4)=  0.d0;
!
ipfsc(1, 1)=5; ipfsc(2, 1)=9; ipfsc(3, 1)=1;
ipfsc(1, 2)=6; ipfsc(2, 2)=9; ipfsc(3, 2)=2;
ipfsc(1, 3)=7; ipfsc(2, 3)=9; ipfsc(3, 3)=3;
ipfsc(1, 4)=8; ipfsc(2, 4)=9; ipfsc(3, 4)=4;
!
!...Boundary face...
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
!print*,ifa,iel,gflag(1:3, ifa)
!
enddo
!
!...Quads
!
do 55 ie=1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
do isv = 1, 4
r = posqc(1,isv)
s = posqc(2,isv)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpq(1,ishp)
ysc = ysc + shpq(ishp)*xpq(2,ishp)
enddo
!
xpqv(1, isv) = xsc
xpqv(2, isv) = ysc
!
enddo

!...Get the normal vector and Jacobian for curved nodes...
do ifa = 1,4

!...Edge 9 (p95)
xpf(1:2, 1) =  xpq(1:2, ipfsc(1, ifa))
xpf(1:2, 2) =  xpq(1:2, ipfsc(2, ifa))
xpf(1:2, 3) = xpqv(1:2, ipfsc(3, ifa))
!
r  = -1.d0
!
dshpr(1) = -0.5d0 + r
dshpr(2) =  0.5d0 + r
dshpr(3) = -2.d0*r

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dshpr(ishp)*xpf(1, ishp)
dydr = dydr + dshpr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
gesgq(1, 8+ifa , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 8+ifa , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 8+ifa , ie) = djaco*2.d0

enddo !do ifa = 1,4
!
55 enddo  !...(1)ifa=1,nelem
!
end subroutine getfnds_lagsubgc
!
!...subroutine: Calculate the F^* N ds for interior faces using SMS...
!
subroutine getfnds_lagsmsif_hybrid(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
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
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(1:nvfac,1:3)::ipift
integer,dimension(1:nvfac,1:4)::ipifq
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8::dshpr(3),shpq(9)
real*8, dimension(nvtri):: shpt
!...local real array
real*8,dimension(1:ndimn, 1:nvtri+1)::xpt
real*8,dimension(1:ndimn, 1:nvqua)::xpq
real*8::posqc(1:2,1:4),xpqv(1:2,1:4)
real*8::postc(1:2,1:3),xptv(1:2,1:3)
!...local real number
real*8::anx, any,xsc,ysc,dxdr,dydr
real*8::dr, ds, rc, sc,r,s,djaco
real*8::c16,c10
!
data c16   /0.1666666666666666d0 /
data c10   /1.d0 /
!
postc(1, 1)=  0.25d0; postc(2, 1)= 0.25d0;
postc(1, 2)=  0.5d0;  postc(2, 2)= 0.25d0;
postc(1, 3)=  0.25d0; postc(2, 3)= 0.5d0;
!
ipift(1, 1)=4; ipift(2, 1)=6; ipift(3, 1)=1;
ipift(1, 2)=5; ipift(2, 2)=4; ipift(3, 2)=2;
ipift(1, 3)=6; ipift(2, 3)=5; ipift(3, 3)=3;
!
posqc(1, 1)=  0.0d0; posqc(2, 1)= -.5d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)=  0.d0;
posqc(1, 3)=  0.0d0; posqc(2, 3)=  .5d0;
posqc(1, 4)= -0.5d0; posqc(2, 4)=  0.d0;
!
ipifq(1, 1)=5; ipifq(2, 1)=9; ipifq(3, 1)=1;
ipifq(1, 2)=6; ipifq(2, 2)=9; ipifq(3, 2)=2;
ipifq(1, 3)=7; ipifq(2, 3)=9; ipifq(3, 3)=3;
ipifq(1, 4)=8; ipifq(2, 4)=9; ipifq(3, 4)=4;
!
!...Boundary face...
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
!print*,ifa,iel,gflag(1:3, ifa)
!
enddo
!
!...Trias
!
do ie=1, ntria
!
ipt(1:nvtri) = iptri(1:nvtri,ie)

!...coordinates
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
do isv = 1, 3
r = postc(1,isv)
s = postc(2,isv)

!...  shape function
shpt(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shpt(2) = -r*(c10-2.d0*r)
shpt(3) = -s*(c10-2.d0*s)
shpt(4) = 4.d0*r*(c10-r-s)
shpt(5) = 4.d0*r*s
shpt(6) = 4.d0*s*(c10-r-s)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, nptri
xsc = xsc + shpt(ishp)*xpt(1,ishp)
ysc = ysc + shpt(ishp)*xpt(2,ishp)
enddo
!
xptv(1, isv) = xsc
xptv(2, isv) = ysc
!
enddo

!...Get the normal vector and Jacobian for curved nodes...
do ifa = 1,3

!...Edge 9 (p95)
xpf(1:2, 1) =  xpt(1:2, ipift(1, ifa))
xpf(1:2, 2) =  xpt(1:2, ipift(2, ifa))
xpf(1:2, 3) = xptv(1:2, ipift(3, ifa))
!
r  = -1.d0
!
dshpr(1) = -0.5d0 + r
dshpr(2) =  0.5d0 + r
dshpr(3) = -2.d0*r

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dshpr(ishp)*xpf(1, ishp)
dydr = dydr + dshpr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
gesgt(1, 6+ifa , ie) = anx/sqrt(anx**2 + any**2)
gesgt(2, 6+ifa , ie) = any/sqrt(anx**2 + any**2)
gesgt(3, 6+ifa , ie) = djaco*2.d0
!
r  = 1.d0
!
dshpr(1) = -0.5d0 + r
dshpr(2) =  0.5d0 + r
dshpr(3) = -2.d0*r

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dshpr(ishp)*xpf(1, ishp)
dydr = dydr + dshpr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
gesgt(1, 12+ifa , ie) = anx/sqrt(anx**2 + any**2)
gesgt(2, 12+ifa , ie) = any/sqrt(anx**2 + any**2)
gesgt(3, 12+ifa , ie) = djaco*2.d0

enddo

enddo !enddo for tria
!
!...Quads
!
do 55 ie=1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
do isv = 1, 4
r = posqc(1,isv)
s = posqc(2,isv)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpq(1,ishp)
ysc = ysc + shpq(ishp)*xpq(2,ishp)
enddo
!
xpqv(1, isv) = xsc
xpqv(2, isv) = ysc
!
enddo

!...Get the normal vector and Jacobian for curved nodes...
do ifa = 1,4

!...Edge 9 (p95)
xpf(1:2, 1) =  xpq(1:2, ipifq(1, ifa))
xpf(1:2, 2) =  xpq(1:2, ipifq(2, ifa))
xpf(1:2, 3) = xpqv(1:2, ipifq(3, ifa))
!
r  = -1.d0
!
dshpr(1) = -0.5d0 + r
dshpr(2) =  0.5d0 + r
dshpr(3) = -2.d0*r

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dshpr(ishp)*xpf(1, ishp)
dydr = dydr + dshpr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
gesgq(1, 8+ifa , ie) = anx/sqrt(anx**2 + any**2)
gesgq(2, 8+ifa , ie) = any/sqrt(anx**2 + any**2)
gesgq(3, 8+ifa , ie) = djaco*2.d0

enddo !do ifa = 1,4
!
55 enddo  !...(1)ifa=1,nelem
!
end subroutine getfnds_lagsmsif_hybrid

!
!...subroutine: Calculate the F^* N dsfor all faces for hybrid grids...
!
subroutine getfnds_lag_simpsubg(gflag,gesgq, gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(inout)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad),      intent(inout)::gelagq
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(inout)::gesgq
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
real*8::posiq(2, 12)
real*8::shp(nptri),dspr(nptri),dsps(nptri)
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
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
posit(1, 1) = 0.d0; posit(2, 1) = 0.d0;
posit(1, 2) = 1.d0; posit(2, 2) = 0.d0;
posit(1, 3) = 1.d0; posit(2, 3) = 0.d0;
posit(1, 4) = 0.d0; posit(2, 4) = 1.d0;
posit(1, 5) = 0.d0; posit(2, 5) = 1.d0;
posit(1, 6) = 0.d0; posit(2, 6) = 0.d0;
posit(1, 7) = .5d0; posit(2, 7) = 0.d0;
posit(1, 8) = .5d0; posit(2, 8) = .5d0;
posit(1, 9) = 0.d0; posit(2, 9) = .5d0;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;
!
do 100 ie=1, ntria!...(1)ifa=1,nafac
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipt(1:nvtri) = inpoel(1:nvtri,ie)
!
!...coordinates
!
if(ncurv==0)then
xpt(1, 1:3) = coord(1, iptri(1:3,ie))
xpt(2, 1:3) = coord(2, iptri(1:3,ie))
!
xpt(1:2,4) = 0.5d0*(xpt(1:2,1) + xpt(1:2,2))
xpt(1:2,5) = 0.5d0*(xpt(1:2,2) + xpt(1:2,3))
xpt(1:2,6) = 0.5d0*(xpt(1:2,1) + xpt(1:2,3))
elseif(ncurv==1)then
xpt(1, 1:nptri) = coord(1,iptri(1:nptri, ie))
xpt(2, 1:nptri) = coord(2,iptri(1:nptri, ie))
endif
!
do ig = 1, 9!...(2)ig = 1,ngausd
!
r  = posit(1,ig)
s  = posit(2,ig)
!
!...Shape functions and their derivatives...
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
do ishp = 1, nptri
dxdr = dxdr + dspr(ishp)*xpt(1,ishp)
dxds = dxds + dsps(ishp)*xpt(1,ishp)

dydr = dydr + dspr(ishp)*xpt(2,ishp)
dyds = dyds + dsps(ishp)*xpt(2,ishp)
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
enddo
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa
!
!...Identify the local No. of one internal face for left cell...
!
if(ig.eq.1.or.ig.eq.2.or.ig.eq.7)then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 1.d0;
!
idfal = ig;
!
elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.8)then
!
vnorm(1) = sqrt(2.d0)*0.5d0; vnorm(2) = sqrt(2.d0)*0.5d0;
larea    = sqrt(2.d0)
!
idfal = ig;
elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.9)then
!
vnorm(1) = -1.d0;           vnorm(2) = 0.d0;
larea    = 1.d0
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
gelag(1, idfal, ie) = anx*larea/farea
gelag(2, idfal, ie) = any*larea/farea
gelag(3, idfal, ie) = farea
!
enddo !...ig from 1...9
!
100 enddo  !...(1)ifa=1,nafac
!
! print*,'vnotmfn',gelag(1:3, 1:6, 170)
!
!
do 200 ie=1, nquad !...(1)ifa=1,nquad
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...coordinates
!
if(ncurv==0)then
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif
!
do ig = 1, 12!...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa
!
!...Identify the local No. of one internal face for left cell...
!
if(ig.eq.1.or.ig.eq.2.or.ig.eq.9)then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 2.d0;
!
idfal = ig;
!
elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.10)then
!
vnorm(1) = 1.d0; vnorm(2) = 0.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.11)then
!
vnorm(1) = 0.d0;           vnorm(2) = 1.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.7.or.ig.eq.8.or.ig.eq.12)then
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
enddo !...ig from 1...9
!
!gesgq(1:2, 1 , ie) = gelagq(1:2, 1, ie); gesgq(3, 1 , ie) = 0.5d0*gelagq(3, 1, ie)
!gesgq(1:2, 2 , ie) = gelagq(1:2, 2, ie); gesgq(3, 2 , ie) = 0.5d0*gelagq(3, 2, ie)

!gesgq(1:2, 3 , ie) = gelagq(1:2, 3, ie); gesgq(3, 3 , ie) = 0.5d0*gelagq(3, 3, ie)
!gesgq(1:2, 4 , ie) = gelagq(1:2, 4, ie); gesgq(3, 4 , ie) = 0.5d0*gelagq(3, 4, ie)

!gesgq(1:2, 5 , ie) = gelagq(1:2, 5, ie); gesgq(3, 5 , ie) = 0.5d0*gelagq(3, 5, ie)
!gesgq(1:2, 6 , ie) = gelagq(1:2, 6, ie); gesgq(3, 6 , ie) = 0.5d0*gelagq(3, 6, ie)

!gesgq(1:2, 7 , ie) = gelagq(1:2, 7, ie); gesgq(3, 7 , ie) = 0.5d0*gelagq(3, 7, ie)
!gesgq(1:2, 8 , ie) = gelagq(1:2, 8, ie); gesgq(3, 8 , ie) = 0.5d0*gelagq(3, 8, ie)
!
gesgq(1:2, 1 , ie) = gelagq(1:2, 1, ie); gesgq(3, 1 , ie) = 0.5d0*gelagq(3, 1, ie)
gesgq(1:2, 2 , ie) = gelagq(1:2, 2, ie); gesgq(3, 2 , ie) = 0.5d0*gelagq(3, 2, ie)

gesgq(1:2, 3 , ie) = gelagq(1:2, 3, ie); gesgq(3, 3 , ie) = 0.5d0*gelagq(3, 3, ie)
gesgq(1:2, 4 , ie) = gelagq(1:2, 4, ie); gesgq(3, 4 , ie) = 0.5d0*gelagq(3, 4, ie)

gesgq(1:2, 5 , ie) = gelagq(1:2, 5, ie); gesgq(3, 5 , ie) = 0.5d0*gelagq(3, 5, ie)
gesgq(1:2, 6 , ie) = gelagq(1:2, 6, ie); gesgq(3, 6 , ie) = 0.5d0*gelagq(3, 6, ie)

gesgq(1:2, 7 , ie) = gelagq(1:2, 7, ie); gesgq(3, 7 , ie) = 0.5d0*gelagq(3, 7, ie)
gesgq(1:2, 8 , ie) = gelagq(1:2, 8, ie); gesgq(3, 8 , ie) = 0.5d0*gelagq(3, 8, ie)
!
!gesgq(1:3, 9:12 , ie) = 0.d0
!
200 enddo  !...(1)ifa=1,nquad
!
!  print*,'Inside getfnds_lag'
!
end subroutine getfnds_lag_simpsubg
!
!...subroutine: Calculate the Riemann input for hybrid curved quad grids using subgrid....
!
subroutine getriem_quadsubg(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
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
!
real*8, dimension(1:2, 1:2, 1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
!
real*8, dimension(1:2, 1:2, 1:2, 1:4, 1:4, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:4, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:4, 1:nquad), intent(out)::snsigmlq
!...Local integer
integer::ie, ideg, ielem, ifa, iv, isg, ivsg, ifsg
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(4, 4)::fnqsg
integer,dimension(4, 4)::ipqsg
!...local real array
real*8,dimension(1:3, 1:nvqua)::bq
real*8,dimension(1:3, 1:4)::bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nq,1:4)::unsgq
real*8::aujmp(1:3, 1:4)
real*8::vnorm(1:3, 1:2, 1:4)
real*8::sigma(1:2, 1:2, 1:4)
real*8,dimension(1:4)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndegr, 1:4)::unksgq
real*8,dimension(1:4, 1:4)::geoq_sub
real*8,dimension(1:2,1:4, 1:4)::wfgsq
!...arraies for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx,rhovsg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Geometry information for subgrid quad...
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
fnqsg(1, 1) =  1; fnqsg(2, 1) = -9; fnqsg(3, 1) = 12; fnqsg(4, 1) =  8;
fnqsg(1, 2) =  2; fnqsg(2, 2) =  3; fnqsg(3, 2) =-10; fnqsg(4, 2) =  9;
fnqsg(1, 3) = 10; fnqsg(2, 3) =  4; fnqsg(3, 3) =  5; fnqsg(4, 3) =-11;
fnqsg(1, 4) =-12; fnqsg(2, 4) = 11; fnqsg(3, 4) =  6; fnqsg(4, 4) =  7;
!
wfgsq(1, 1, 1) = 1.d0; wfgsq(2, 1, 1) = 1.d0;
wfgsq(1, 2, 1) = 1.d0; wfgsq(2, 2, 1) = 0.d0;
wfgsq(1, 3, 1) = 0.d0; wfgsq(2, 3, 1) = 0.d0;
wfgsq(1, 4, 1) = 0.d0; wfgsq(2, 4, 1) = 1.d0;

wfgsq(1, 1, 2) = 0.d0; wfgsq(2, 1, 2) = 1.d0;
wfgsq(1, 2, 2) = 1.d0; wfgsq(2, 2, 2) = 1.d0;
wfgsq(1, 3, 2) = 1.d0; wfgsq(2, 3, 2) = 0.d0;
wfgsq(1, 4, 2) = 0.d0; wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 1, 3) = 0.d0; wfgsq(2, 1, 3) = 0.d0;
wfgsq(1, 2, 3) = 0.d0; wfgsq(2, 2, 3) = 1.d0;
wfgsq(1, 3, 3) = 1.d0; wfgsq(2, 3, 3) = 1.d0;
wfgsq(1, 4, 3) = 1.d0; wfgsq(2, 4, 3) = 0.d0;

wfgsq(1, 1, 4) = 1.d0; wfgsq(2, 1, 4) = 0.d0;
wfgsq(1, 2, 4) = 0.d0; wfgsq(2, 2, 4) = 0.d0;
wfgsq(1, 3, 4) = 0.d0; wfgsq(2, 3, 4) = 1.d0;
wfgsq(1, 4, 4) = 1.d0; wfgsq(2, 4, 4) = 1.d0;
!
wfgsq(1, 1, 1) = 4.d0/6.d0; wfgsq(2, 1, 1) = 4.d0/6.d0;
wfgsq(1, 2, 1) = 4.d0/3.d0; wfgsq(2, 2, 1) = 0.d0/3.d0;
wfgsq(1, 3, 1) = 1.d0;      wfgsq(2, 3, 1) = 1.d0;
wfgsq(1, 4, 1) = 0.d0/3.d0; wfgsq(2, 4, 1) = 4.d0/3.d0;

wfgsq(1, 1, 2) = 0.d0/3.d0; wfgsq(2, 1, 2) = 4.d0/3.d0;
wfgsq(1, 2, 2) = 4.d0/6.d0; wfgsq(2, 2, 2) = 4.d0/6.d0;
wfgsq(1, 3, 2) = 4.d0/3.d0; wfgsq(2, 3, 2) = 0.d0/3.d0;
wfgsq(1, 4, 2) = 1.d0;      wfgsq(2, 4, 2) = 1.d0;

wfgsq(1, 1, 3) = 1.d0;      wfgsq(2, 1, 3) = 1.d0;
wfgsq(1, 2, 3) = 0.d0/3.d0; wfgsq(2, 2, 3) = 4.d0/3.d0;
wfgsq(1, 3, 3) = 4.d0/6.d0; wfgsq(2, 3, 3) = 4.d0/6.d0;
wfgsq(1, 4, 3) = 4.d0/3.d0; wfgsq(2, 4, 3) = 0.d0/3.d0;

wfgsq(1, 1, 4) = 4.d0/3.d0; wfgsq(2, 1, 4) = 0.d0/3.d0;
wfgsq(1, 2, 4) = 1.d0;      wfgsq(2, 2, 4) = 1.d0;
wfgsq(1, 3, 4) = 0.d0/3.d0; wfgsq(2, 3, 4) = 4.d0/3.d0;
wfgsq(1, 4, 4) = 4.d0/6.d0; wfgsq(2, 4, 4) = 4.d0/6.d0;!
!
wfgsq(1, 1, 1) = 1.d0; wfgsq(2, 1, 1) = 1.d0;
wfgsq(1, 2, 1) = 1.d0; wfgsq(2, 2, 1) = 0.d0;
wfgsq(1, 3, 1) = 0.d0; wfgsq(2, 3, 1) = 0.d0;
wfgsq(1, 4, 1) = 0.d0; wfgsq(2, 4, 1) = 1.d0;

wfgsq(1, 1, 2) = 0.d0; wfgsq(2, 1, 2) = 1.d0;
wfgsq(1, 2, 2) = 1.d0; wfgsq(2, 2, 2) = 1.d0;
wfgsq(1, 3, 2) = 1.d0; wfgsq(2, 3, 2) = 0.d0;
wfgsq(1, 4, 2) = 0.d0; wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 1, 3) = 0.d0; wfgsq(2, 1, 3) = 0.d0;
wfgsq(1, 2, 3) = 0.d0; wfgsq(2, 2, 3) = 1.d0;
wfgsq(1, 3, 3) = 1.d0; wfgsq(2, 3, 3) = 1.d0;
wfgsq(1, 4, 3) = 1.d0; wfgsq(2, 4, 3) = 0.d0;

wfgsq(1, 1, 4) = 1.d0; wfgsq(2, 1, 4) = 0.d0;
wfgsq(1, 2, 4) = 0.d0; wfgsq(2, 2, 4) = 0.d0;
wfgsq(1, 3, 4) = 0.d0; wfgsq(2, 3, 4) = 1.d0;
wfgsq(1, 4, 4) = 1.d0; wfgsq(2, 4, 4) = 1.d0;

!
!...Quads...
!
do 350 ie = 1,nquad !...(1)ie = 1,nelem
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
!...shape functions
!
dr = 1.0d0
ds = 1.0d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
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
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...cell averaged value...
!
if(ndens.eq.3)then
!
!...now configuration
!
rhomc = 1.d0/unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc !...Cell center density...
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!if(ie.eq.1)print*,'sound speed',ie,rhoct, pctr,ectr
!
!...1. Get u, v, e at vertices....
!
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
enddo
!
!...Get density for subgrids...
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhosubg2(xpq,xpqi,unksgq, geoq_sub)
!
!...2. Loop over subgrid....
!
do isg = 1, 4
!
!...Normal vector for every face...
!
vnorm(1:2, 1, 1) = sign(1, fnqsg(4, isg))*gesgq(1:2,abs(fnqsg(4, isg)), ie); vnorm(  3, 1, 1) = gesgq(3,abs(fnqsg(4, isg)), ie);
vnorm(1:2, 2, 1) = sign(1, fnqsg(1, isg))*gesgq(1:2,abs(fnqsg(1, isg)), ie); vnorm(  3, 2, 1) = gesgq(3,abs(fnqsg(1, isg)), ie);

vnorm(1:2, 1, 2) = sign(1, fnqsg(1, isg))*gesgq(1:2,abs(fnqsg(1, isg)), ie); vnorm(  3, 1, 2) = gesgq(3,abs(fnqsg(1, isg)), ie);
vnorm(1:2, 2, 2) = sign(1, fnqsg(2, isg))*gesgq(1:2,abs(fnqsg(2, isg)), ie); vnorm(  3, 2, 2) = gesgq(3,abs(fnqsg(2, isg)), ie);

vnorm(1:2, 1, 3) = sign(1, fnqsg(2, isg))*gesgq(1:2,abs(fnqsg(2, isg)), ie); vnorm(  3, 1, 3) = gesgq(3,abs(fnqsg(2, isg)), ie);
vnorm(1:2, 2, 3) = sign(1, fnqsg(3, isg))*gesgq(1:2,abs(fnqsg(3, isg)), ie); vnorm(  3, 2, 3) = gesgq(3,abs(fnqsg(3, isg)), ie);

vnorm(1:2, 1, 4) = sign(1, fnqsg(3, isg))*gesgq(1:2,abs(fnqsg(3, isg)), ie); vnorm(  3, 1, 4) = gesgq(3,abs(fnqsg(3, isg)), ie);
vnorm(1:2, 2, 4) = sign(1, fnqsg(4, isg))*gesgq(1:2,abs(fnqsg(4, isg)), ie); vnorm(  3, 2, 4) = gesgq(3,abs(fnqsg(4, isg)), ie);
!
do ifsg =1 ,2

vnorm(  3, ifsg, 1) = vnorm(  3, ifsg, 1)*wfgsq(ifsg, 1, isg)
vnorm(  3, ifsg, 2) = vnorm(  3, ifsg, 2)*wfgsq(ifsg, 2, isg)
vnorm(  3, ifsg, 3) = vnorm(  3, ifsg, 3)*wfgsq(ifsg, 3, isg)
vnorm(  3, ifsg, 4) = vnorm(  3, ifsg, 4)*wfgsq(ifsg, 4, isg)
enddo
!
do ivsg = 1,4
!
if(ndens.eq.3)then
!
rcv = geoq_sub(1, isg); scv = geoq_sub(2, isg)
!
bqv(1, ivsg) = 1.d0
bqv(2, ivsg) = (xvq(ivsg)-rcv)/dr  !...bqv ....
bqv(3, ivsg) = (yvq(ivsg)-scv)/ds
!
rhovsg =0.d0
!
do ideg = 1,mdegr
rhovsg =  unksgq(1,isg)!!rhovsg + unksgq(ideg,isg)*bqv(ideg, ivsg)
enddo
!
endif
rhovsg = unknvq(1, ipqsg(ivsg, isg))
uvtx = unknvq(2, ipqsg(ivsg, isg))
vvtx = unknvq(3, ipqsg(ivsg, isg))
evtx = unknvq(4, ipqsg(ivsg, isg))
!
pvtx = max(eps, (gamlg-1.d0)*rhovsg*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!pvtx = 0.9d0*pvtx + 0.1d0*rhovsg*gamlg*pctr/rhoct

!pvtx = pctr/rhoct*rhovsg
!
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx
!
!...Limiter
!
if(nlimi.eq.6)then
!
!if(ndens.eq.3)then
!rhovsg = unksgq(1,isg) + aflim(1, ielem)*(rhovsg - unksgq(1,isg))
!endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bq(2, ipqsg(ivsg, isg)) + duds*bq(3, ipqsg(ivsg, isg))
vvtx = unkno(1,3,ielem)  + dvdr*bq(2, ipqsg(ivsg, isg)) + dvds*bq(3, ipqsg(ivsg, isg))
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!!
!...updtae unknv(2:3,:)
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx
!

endif
!
!...Get stress tensor at nodes
!
sigma(1, 1, ivsg) = -pvtx
sigma(1, 2, ivsg) = 0.d0
sigma(2, 1, ivsg) = 0.d0
sigma(2, 2, ivsg) = -pvtx!
!
!if(ipq(iv).eq.2) print*,'velocity 8',ie, afvec(1:2,1:2,ie)
!
!...Get the a_c (unit vector)
!
aujmp(1:2, ivsg) = vlave(1:2, ipq(ipqsg(ivsg, isg))) - unsgq(2:3, ivsg)
!if(ie==22) print*,'adjumpxxx22', vlave(1:2, ip(iv)) , unknv(2:3, iv), ip(iv)
acnx = aujmp(1, ivsg)
acny = aujmp(2, ivsg)
!
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ivsg) = 1.e-11!0.d0;
!print*,'point are reset', ip(iv)
else
aujmp(1:2, ivsg) = aujmp(1:2, ivsg)/sqrt(acnx**2 + acny**2)
endif
!
aujmp(3, ivsg) = sqrt(acnx**2 + acny**2)
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the variables at the center...
!
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
! print*,'sound speed2',sdctr,ie ,rhoct, uctr,vctr,pctr
!
aujmp(3,:)=aujmp(3,:)/sdctr
!
!...Get impedence coefficient...
!
do ivsg   = 1, 4
dux= vlave(1, ipq(ipqsg(ivsg, isg)))-unsgq(2, ivsg)
duy= vlave(2, ipq(ipqsg(ivsg, isg)))-unsgq(3, ivsg)
deltu = sqrt(dux**2 + duy**2)
murie(ivsg) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
! if(ipq(ipqsg(ivsg, isg)).eq.53) then
!print*,'murie22ang',ie,vlave(1:2, ipq(ipqsg(ivsg, isg))),unsgq(2:3, ivsg)
! endif
enddo
!
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do ivsg  = 1, 4
!
!...Local No. of subgrid No.
!
iv = ipqsg(ivsg, isg)
!
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!...Call Riemann solver...
!
call getriecoef_matrixnew(murie(ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:3, ivsg), &
unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:2, ivsg), &
!unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipq(iv)) = munacn(1:2, 1, ipq(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipq(iv)) = munacn(1:2, 2, ipq(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipq(iv)) = munacu(1:2, ipq(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipq(iv)) = snsigm(1:2, ipq(iv)) + snsigm_rie(1:2)
!
!if(ipq(iv).eq.23) print*,'p36 muacn(vv) post',murie(ivsg),ie,isg,ivsg,ipq(iv),sigma(1:2, 1:2, ivsg),vnorm(1:3, ifa, ivsg)
!
!
!
!...Local variable...
!
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
350 enddo  !...(1)ie = 1,nelem!

end subroutine getriem_quadsubg
!
!...subroutine: Calculate the nodal velocity U_p^* (mass center) for quadratic meshes with one sub-cell scheme...
!
subroutine getndvelo_lagsubg(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
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
real*8,dimension(1:ndimn,1:2,1:3,1:4, 1:ntria),  intent(out)::fstart !...Riemann forces
real*8,dimension(1:ndimn,1:2,1:4,1:4, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),             intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec

integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop, isg, ivsg
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(4, 4)::ipqsg
integer,dimension(3, 4)::iptsg
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
allocate (munaclq(1:2, 1:2, 1:2, 1:4, 1:4, 1:nquad), munaulq(1:ndimn, 1:2, 1:4, 1:4,  1:nquad),&
snsigmlq(1:ndimn, 1:2,  1:4, 1:4,  1:nquad))
!allocate (munaclt(1:2, 1:2, 1:2, 1:4, 1:3, 1:ntria), munault(1:ndimn, 1:2, 1:4, 1:3,  1:ntria),&
!snsigmlt(1:ndimn, 1:2,  1:4, 1:3,  1:ntria))
allocate (munaclt(1:2, 1:2, 1:2, 1:3, 1:4, 1:ntria), munault(1:ndimn, 1:2, 1:3, 1:4,  1:ntria),&
snsigmlt(1:ndimn, 1:2,  1:3, 1:4,  1:ntria))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
!
!...Part I: Specify some gauss points...
!
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

!...Tria
if(ntria.gt.0) call getriem_triasubg3(iptri, geoel, gesgt, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt,coord, coold, aflim, afvec)

!...Quad
if(nquad.gt.0) call getriem_quadsubg2(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec, itime)

!...Boundary condition
!'In the future'
!!call getbccurv_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)
!!call getbc_lagc(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)
call getbc_lagc_general(bface, intfac, fpres, coord, munacn, munacu, snsigm, itime)

!...Periodic boundary condition for 1D isentropic Sin problem...
if(ncase.eq.12) call getbc_prdic(bface, munacn, munacu, snsigm)


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

!if(itime.ge.60.and.ipoin.eq.90)print*,itime,ipoin,ustar(1:2,ipoin),detma,munacn(1,1,ipoin),snsigm(1:2, ipoin),&
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
do isg = 1, 4
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
do isg = 1, 4
do ivsg = 1, 4
!
iv =  ipqsg(ivsg, isg)
!
do ifa =1, 2
fstarq(1, ifa, ivsg, isg, ie) = snsigmlq(1, ifa, ivsg, isg, ie) + &
munaclq(1, 1, ifa, ivsg, isg,ie)*ustar(1, ipq(iv))+&
munaclq(2, 1, ifa, ivsg, isg,ie)*ustar(2, ipq(iv)) - munaulq(1, ifa, ivsg, isg, ie)
fstarq(2, ifa, ivsg, isg,ie) = snsigmlq(2, ifa, ivsg, isg, ie) + &
munaclq(1, 2, ifa, ivsg, isg, ie)*ustar(1, ipq(iv))+&
munaclq(2, 2, ifa, ivsg, isg, ie)*ustar(2, ipq(iv))- munaulq(2, ifa, ivsg, isg, ie)
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
!
!
deallocate (munacn, bpres, fpres)
deallocate (munacu, snsigm, bnorm)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lagsubg
!
!...Get the mass matrix for lagrangian based on mass center...
!
subroutine  getrhosubg(xpq,xpqi,unksgq, geoq_sub)
use constant
implicit none
!...Input
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:2, 1:npqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:4), intent(out)::unksgq
real*8,dimension(1:4, 1:4), intent(out)::geoq_sub
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, ielem, id, isg, ideg, iunk

integer,dimension(4, 4)::ipqsg
!...Local real array
real*8,dimension(1:ndegr,1:4)::rhsel
real*8,dimension(1:nmati, 1:4)::matin
real*8,dimension(1:2, 1:npqua)::xpqsg,xpqisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(3, 3)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
!...Local real
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::b2,b3,b4,b5,b6
real*8::masel,xgaus,ygaus
!
real*8,allocatable::x5(:,:),b55(:),mmatr(:,:)
!
data c10 / 1.0d0 /
!
if(npoly==1) allocate(x5(3,3), mmatr(3,3), b55(3))
if(npoly==2) allocate(x5(5,5), mmatr(5,5), b55(5))
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
call ruqope(2, ngausdq, posiq, weighq)
!
!...dr and ds...
!
dr = 1.d0
ds = 1.d0
!
!...Part I: Get mass matrix for 4 subgrids...
!
do isg = 1, 4
!
!...Now configuration...
!
xpqsg(1,1:4) =  xpq(1, ipqsg(1:4, isg))
xpqsg(2,1:4) =  xpq(2, ipqsg(1:4, isg))
!
xpqsg(1:2,5) = 0.5d0*(xpqsg(1:2,1)+xpqsg(1:2,2))
xpqsg(1:2,6) = 0.5d0*(xpqsg(1:2,2)+xpqsg(1:2,3))
xpqsg(1:2,7) = 0.5d0*(xpqsg(1:2,3)+xpqsg(1:2,4))
xpqsg(1:2,8) = 0.5d0*(xpqsg(1:2,4)+xpqsg(1:2,1))
xpqsg(1:2,9) = 0.5d0*(xpqsg(1:2,5)+xpqsg(1:2,7))
!
!...Initial configuration...
!
xpqisg(1,1:4) =  xpqi(1, ipqsg(1:4, isg))
xpqisg(2,1:4) =  xpqi(2, ipqsg(1:4, isg))
!
xpqisg(1:2,5) = 0.5d0*(xpqisg(1:2,1)+xpqisg(1:2,2))
xpqisg(1:2,6) = 0.5d0*(xpqisg(1:2,2)+xpqisg(1:2,3))
xpqisg(1:2,7) = 0.5d0*(xpqisg(1:2,3)+xpqisg(1:2,4))
xpqisg(1:2,8) = 0.5d0*(xpqisg(1:2,4)+xpqisg(1:2,1))
xpqisg(1:2,9) = 0.5d0*(xpqisg(1:2,5)+xpqisg(1:2,7))
!
!...Initialze parameters...
!
mmatr = 0.d0
!
volel = 0.d0
xcel = 0.d0
ycel = 0.d0
!
voleli = 0.d0
xceli = 0.d0
yceli = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, npqua
dxdri = dxdri + dsprq(ishp)*xpqisg(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqisg(1,ishp)

dydri = dydri + dsprq(ishp)*xpqisg(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqisg(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)
!
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
mmatr(1,1) = mmatr(1,1) + djaco
mmatr(2,1) = mmatr(2,1) + b2*djaco
mmatr(3,1) = mmatr(3,1) + b3*djaco

mmatr(2,2) = mmatr(2,2) + b2*b2*djaco
mmatr(3,2) = mmatr(3,2) + b2*b3*djaco

mmatr(3,3) = mmatr(3,3) + b3*b3*djaco
!
!...Parameters for now configuration...
!
volel = volel + djaco
xcel = xcel + djaco*r!xgaus
ycel = ycel + djaco*s!ygaus
!
!...Parameters for initial configuration...
!
voleli = voleli + djacoi
xceli = xceli + djacoi*r!xgaus
yceli = yceli + djacoi*s!ygaus
!
enddo
!
!...Volume center for subgrids...
!
geoq_sub(1, isg) = xcel/volel
geoq_sub(2, isg) = ycel/volel
geoq_sub(3, isg) = xceli/voleli
geoq_sub(4, isg) = yceli/voleli
!
!...Mass matrix
!
mmatr(1, 2) = mmatr(2, 1)
mmatr(1, 3) = mmatr(3, 1)
mmatr(2, 3) = mmatr(3, 2)
!
!print*,'ielem', ie, f0, volel

if(npoly==1)then
!   if(ie==2) print*,'ie',mmatr(:,:)
x5 = 0.d0
b55 = 0.d0
call getinvmat(3, mmatr, x5, b55)
!
matin(1,isg) = x5(1,1)
matin(2,isg) = x5(1,2)
matin(3,isg) = x5(1,3)
matin(4,isg) = x5(2,2)
matin(5,isg) = x5(2,3)
matin(6,isg) = x5(3,3)
!
endif
enddo
!
!...Part II: Get unksg...
!
rhsel = 0.d0
!
do isg = 1, 4
!
xpqsg(1,1:4) =  xpq(1, ipqsg(1:4, isg))
xpqsg(2,1:4) =  xpq(2, ipqsg(1:4, isg))
!
xpqsg(1:2,5) = 0.5d0*(xpqsg(1:2,1)+xpqsg(1:2,2))
xpqsg(1:2,6) = 0.5d0*(xpqsg(1:2,2)+xpqsg(1:2,3))
xpqsg(1:2,7) = 0.5d0*(xpqsg(1:2,3)+xpqsg(1:2,4))
xpqsg(1:2,8) = 0.5d0*(xpqsg(1:2,4)+xpqsg(1:2,1))
xpqsg(1:2,9) = 0.5d0*(xpqsg(1:2,5)+xpqsg(1:2,7))
!
rc = geoq_sub(3, isg)
sc = geoq_sub(4, isg)
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
xg = r
yg = s
!
b(1) = 1.d0
b(2) = (xg-rc)/dr
b(3) = (yg-sc)/ds
!
call getrhoig_quadcurv(rhoi, xpqi)
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg, isg)=rhsel(ideg, isg) + rhoi*b(ideg)*djaco
enddo
enddo
!
!
if(npoly==1)then
m(1,1) = matin(1, isg)
m(1,2) = 0.d0!matin(2, ie)
m(1,3) = 0.d0!matin(3, ie)

m(2,1) = m(1,2)
m(2,2) = matin(4, isg)
m(2,3) = matin(5, isg)

m(3,1) = m(1,3)
m(3,2) = m(2,3)
m(3,3) = matin(6, isg)
endif
!
do id =1,ndegr !...(4)id =1,ndegr
!
unint = 0.d0
do iunk = 1,ndegr
unint(1) = unint(1) + m(id, iunk)*rhsel(iunk,isg)
enddo
!
!...step 2
!
unksgq(id, isg)= unint(1)
!
enddo !...(4)id =1,ndegr
!
enddo
!
end subroutine  getrhosubg
!
!...Get the mass matrix for lagrangian based on mass center...
!
subroutine  getrhosubg2(xpq,xpqi,unksgq, geoq_sub)
use constant
implicit none
!...Input
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:2, 1:npqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:nq), intent(out)::unksgq
real*8,dimension(1:4, 1:4), intent(out)::geoq_sub
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, ielem, id, isg, ideg, iunk

integer,dimension(4, 4)::ipqsg
!...Local real array
real*8,dimension(1:ndegr,1:4)::rhsel
real*8,dimension(1:nmati, 1:4)::matin
real*8,dimension(1:2, 1:npqua)::xpqsg,xpqisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(3, 3)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
!...Local real
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::b2,b3,b4,b5,b6
real*8::masel,xgaus,ygaus
!
real*8,allocatable::x5(:,:),b55(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
if(npoly==1) allocate(x5(3,3), mmatr(3,3), b55(3))
if(npoly==2) allocate(x5(5,5), mmatr(5,5), b55(5))
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
call ruqope(2, ngausdq, posiq, weighq)
!
!...dr and ds...
!
dr = 1.d0
ds = 1.d0
!
!...Part 0: Get cell averaged r and s...
!
do isg = 1, 4
!
!...Now configuration...
!
xpqsg(1,1:4) =  xpq(1, ipqsg(1:4, isg))
xpqsg(2,1:4) =  xpq(2, ipqsg(1:4, isg))
!
xpqsg(1:2,5) = 0.5d0*(xpqsg(1:2,1)+xpqsg(1:2,2))
xpqsg(1:2,6) = 0.5d0*(xpqsg(1:2,2)+xpqsg(1:2,3))
xpqsg(1:2,7) = 0.5d0*(xpqsg(1:2,3)+xpqsg(1:2,4))
xpqsg(1:2,8) = 0.5d0*(xpqsg(1:2,4)+xpqsg(1:2,1))
xpqsg(1:2,9) = 0.5d0*(xpqsg(1:2,5)+xpqsg(1:2,7))

xpqsg(1:2,9) = -0.25d0*(xpqsg(1:2,1) + xpqsg(1:2,2) + xpqsg(1:2,3) + xpqsg(1:2,4)) +&
0.5d0*(xpqsg(1:2,5) + xpqsg(1:2,6) + xpqsg(1:2,7) + xpqsg(1:2,8))
!
!...Initial configuration...
!
xpqisg(1,1:4) =  xpqi(1, ipqsg(1:4, isg))
xpqisg(2,1:4) =  xpqi(2, ipqsg(1:4, isg))
!
xpqisg(1:2,5) = 0.5d0*(xpqisg(1:2,1)+xpqisg(1:2,2))
xpqisg(1:2,6) = 0.5d0*(xpqisg(1:2,2)+xpqisg(1:2,3))
xpqisg(1:2,7) = 0.5d0*(xpqisg(1:2,3)+xpqisg(1:2,4))
xpqisg(1:2,8) = 0.5d0*(xpqisg(1:2,4)+xpqisg(1:2,1))
xpqisg(1:2,9) = 0.5d0*(xpqisg(1:2,5)+xpqisg(1:2,7))

xpqisg(1:2,9) = -0.25d0*(xpqisg(1:2,1) + xpqisg(1:2,2) + xpqisg(1:2,3) + xpqisg(1:2,4)) +&
0.5d0*(xpqisg(1:2,5) + xpqisg(1:2,6) + xpqisg(1:2,7) + xpqisg(1:2,8))
!
!...Initialze parameters...
!
mmatr = 0.d0
!
volel = 0.d0
xcel = 0.d0
ycel = 0.d0
!
voleli = 0.d0
xceli = 0.d0
yceli = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, npqua
dxdri = dxdri + dsprq(ishp)*xpqisg(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqisg(1,ishp)

dydri = dydri + dsprq(ishp)*xpqisg(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqisg(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)
!
!...Parameters for now configuration...
!
volel = volel + djaco
xcel = xcel + djaco*r
ycel = ycel + djaco*s
!
!...Parameters for initial configuration...
!
voleli = voleli + djacoi
xceli = xceli + djacoi*r
yceli = yceli + djacoi*s
!
enddo
!
!...Volume center for subgrids...
!
geoq_sub(1, isg) = xcel/volel
geoq_sub(2, isg) = ycel/volel
geoq_sub(3, isg) = xceli/voleli
geoq_sub(4, isg) = yceli/voleli
!
enddo
!
!...Part I: Get mass matrix for 4 subgrids...
!
do isg = 1, 4
!
!...Now configuration...
!
xpqsg(1,1:4) =  xpq(1, ipqsg(1:4, isg))
xpqsg(2,1:4) =  xpq(2, ipqsg(1:4, isg))
!
xpqsg(1:2,5) = 0.5d0*(xpqsg(1:2,1)+xpqsg(1:2,2))
xpqsg(1:2,6) = 0.5d0*(xpqsg(1:2,2)+xpqsg(1:2,3))
xpqsg(1:2,7) = 0.5d0*(xpqsg(1:2,3)+xpqsg(1:2,4))
xpqsg(1:2,8) = 0.5d0*(xpqsg(1:2,4)+xpqsg(1:2,1))
xpqsg(1:2,9) = 0.5d0*(xpqsg(1:2,5)+xpqsg(1:2,7))
!
xpqsg(1:2,9) = -0.25d0*(xpqsg(1:2,1) + xpqsg(1:2,2) + xpqsg(1:2,3) + xpqsg(1:2,4)) +&
0.5d0*(xpqsg(1:2,5) + xpqsg(1:2,6) + xpqsg(1:2,7) + xpqsg(1:2,8))
!
!...Initialze parameters...
!
mmatr = 0.d0
!
rc = geoq_sub(1, isg)
sc = geoq_sub(2, isg)
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
mmatr(1,1) = mmatr(1,1) + djaco

mmatr(2,2) = mmatr(2,2) + b2*b2*djaco
mmatr(3,2) = mmatr(3,2) + b2*b3*djaco

mmatr(3,3) = mmatr(3,3) + b3*b3*djaco
!
enddo
!
!...Mass matrix
!
mmatr(2, 1) = 0.d0
mmatr(3, 1) = 0.d0
mmatr(1, 2) = mmatr(2, 1)
mmatr(1, 3) = mmatr(3, 1)
mmatr(2, 3) = mmatr(3, 2)
!
!print*,'ielem', ie, f0, volel

if(npoly==1)then
!   if(ie==2) print*,'ie',mmatr(:,:)
x5 = 0.d0
b55 = 0.d0
call getinvmat(3, mmatr, x5, b55)
!
matin(1,isg) = x5(1,1)
matin(2,isg) = x5(1,2)
matin(3,isg) = x5(1,3)
matin(4,isg) = x5(2,2)
matin(5,isg) = x5(2,3)
matin(6,isg) = x5(3,3)
!
endif
enddo
!
!...Part II: Get unksg...
!
rhsel = 0.d0
!
do isg = 1, 4
!
!...Initial configuration...
!
xpqisg(1,1:4) =  xpqi(1, ipqsg(1:4, isg))
xpqisg(2,1:4) =  xpqi(2, ipqsg(1:4, isg))
!
xpqisg(1:2,5) = 0.5d0*(xpqisg(1:2,1)+xpqisg(1:2,2))
xpqisg(1:2,6) = 0.5d0*(xpqisg(1:2,2)+xpqisg(1:2,3))
xpqisg(1:2,7) = 0.5d0*(xpqisg(1:2,3)+xpqisg(1:2,4))
xpqisg(1:2,8) = 0.5d0*(xpqisg(1:2,4)+xpqisg(1:2,1))
xpqisg(1:2,9) = 0.5d0*(xpqisg(1:2,5)+xpqisg(1:2,7))
!
xpqisg(1:2,9) = -0.25d0*(xpqisg(1:2,1) + xpqisg(1:2,2) + xpqisg(1:2,3) + xpqisg(1:2,4)) +&
0.5d0*(xpqisg(1:2,5) + xpqisg(1:2,6) + xpqisg(1:2,7) + xpqisg(1:2,8))
!
!
rc = geoq_sub(3, isg)
sc = geoq_sub(4, isg)
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqisg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqisg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqisg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqisg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
xg = r
yg = s
!
b(1) = 1.d0
b(2) = (xg-rc)/dr
b(3) = (yg-sc)/ds
!
call getrhoig_quadcurv(rhoi, xpqi)
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg, isg)=rhsel(ideg, isg) + rhoi*b(ideg)*djaco
enddo
enddo
!
!
if(npoly==1)then
m(1,1) = matin(1, isg)
m(1,2) = 0.d0!matin(2, ie)
m(1,3) = 0.d0!matin(3, ie)

m(2,1) = m(1,2)
m(2,2) = matin(4, isg)
m(2,3) = matin(5, isg)

m(3,1) = m(1,3)
m(3,2) = m(2,3)
m(3,3) = matin(6, isg)
endif
!
do id =1,ndegr !...(4)id =1,ndegr
!
unint = 0.d0
do iunk = 1,ndegr
unint(1) = unint(1) + m(id, iunk)*rhsel(iunk,isg)
enddo
!
!...step 2
!
unksgq(id, isg)= unint(1)
!
enddo !...(4)id =1,ndegr
!
enddo
!
end subroutine  getrhosubg2
!
!...Get the mass matrix for lagrangian based on mass center using exact curvec corner...
!
subroutine  getrhosubg3(xpq,xpqi,unksgq, geoq_sub)
use constant
implicit none
!...Input
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:2, 1:npqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:nq), intent(out)::unksgq
real*8,dimension(1:4, 1:4), intent(out)::geoq_sub
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, ielem, id, isg, ideg, iunk

integer,dimension(4, 4)::ipqsg
!...Local real array
real*8,dimension(1:ndegr,1:4)::rhsel
real*8,dimension(1:nmati, 1:4)::matin
real*8,dimension(1:2, 1:npqua)::xpqsg,xpqisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(3, 3)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
!...Local real
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::b2,b3,b4,b5,b6
real*8::masel,xgaus,ygaus
!
real*8,allocatable::x5(:,:),b55(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
if(npoly==1) allocate(x5(3,3), mmatr(3,3), b55(3))
if(npoly==2) allocate(x5(5,5), mmatr(5,5), b55(5))
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
call ruqope(2, ngausdq, posiq, weighq)
!
!...dr and ds...
!
dr = 1.d0
ds = 1.d0
!
!...Now configuration...
!
xpqsg(1,1:9) =  xpq(1, 1:9)
xpqsg(2,1:9) =  xpq(2, 1:9)
!
!...Initial configuration...
!
xpqisg(1,1:9) =  xpqi(1, 1:9)
xpqisg(2,1:9) =  xpqi(2, 1:9)
!
!...Part 0: Get cell averaged r and s...
!
do isg = 1, 4
!
!...Initialze parameters...
!
mmatr = 0.d0
!
volel = 0.d0
xcel = 0.d0
ycel = 0.d0
!
voleli = 0.d0
xceli = 0.d0
yceli = 0.d0
!
do ig =1,ngausdq
!
if(isg.eq.1)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.2)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.3)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
elseif(isg.eq.4)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
endif
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)*0.25d0
!
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, npqua
dxdri = dxdri + dsprq(ishp)*xpqisg(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqisg(1,ishp)

dydri = dydri + dsprq(ishp)*xpqisg(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqisg(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)*0.25d0
!
!...Parameters for now configuration...
!
volel = volel + djaco
xcel = xcel + djaco*r
ycel = ycel + djaco*s
!
!...Parameters for initial configuration...
!
voleli = voleli + djacoi
xceli = xceli + djacoi*r
yceli = yceli + djacoi*s
!
enddo
!
!...Volume center for subgrids...
!
geoq_sub(1, isg) = xcel/volel
geoq_sub(2, isg) = ycel/volel
geoq_sub(3, isg) = xceli/voleli
geoq_sub(4, isg) = yceli/voleli
!
enddo
!
!...Part I: Get mass matrix for 4 subgrids...
!
do isg = 1, 4
!
!...Initialze parameters...
!
mmatr = 0.d0
!
rc = geoq_sub(1, isg)
sc = geoq_sub(2, isg)
!
do ig =1,ngausdq
!
!
if(isg.eq.1)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.2)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.3)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
elseif(isg.eq.4)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
endif

wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)*0.25d0
!
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
mmatr(1,1) = mmatr(1,1) + djaco

mmatr(2,2) = mmatr(2,2) + b2*b2*djaco
mmatr(3,2) = mmatr(3,2) + b2*b3*djaco

mmatr(3,3) = mmatr(3,3) + b3*b3*djaco
!
enddo
!
!...Mass matrix
!
mmatr(2, 1) = 0.d0
mmatr(3, 1) = 0.d0
mmatr(1, 2) = mmatr(2, 1)
mmatr(1, 3) = mmatr(3, 1)
mmatr(2, 3) = mmatr(3, 2)
!
!print*,'ielem', ie, f0, volel

if(npoly==1)then
!   if(ie==2) print*,'ie',mmatr(:,:)
x5 = 0.d0
b55 = 0.d0
call getinvmat(3, mmatr, x5, b55)
!
matin(1,isg) = x5(1,1)
matin(2,isg) = x5(1,2)
matin(3,isg) = x5(1,3)
matin(4,isg) = x5(2,2)
matin(5,isg) = x5(2,3)
matin(6,isg) = x5(3,3)
!
endif
enddo
!
!...Part II: Get unksg...
!
rhsel = 0.d0
!
do isg = 1, 4
!
!
rc = geoq_sub(3, isg)
sc = geoq_sub(4, isg)
!
do ig =1,ngausdq
!
!
if(isg.eq.1)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.2)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.3)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
elseif(isg.eq.4)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
endif

wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqisg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqisg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqisg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqisg(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)*0.25d0
!
xg = r
yg = s
!
b(1) = 1.d0
b(2) = (xg-rc)/dr
b(3) = (yg-sc)/ds
!
call getrhoig_quadcurv(rhoi, xpqi)
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg, isg)=rhsel(ideg, isg) + rhoi*b(ideg)*djaco
enddo
enddo
!
!
if(npoly==1)then
m(1,1) = matin(1, isg)
m(1,2) = 0.d0!matin(2, ie)
m(1,3) = 0.d0!matin(3, ie)

m(2,1) = m(1,2)
m(2,2) = matin(4, isg)
m(2,3) = matin(5, isg)

m(3,1) = m(1,3)
m(3,2) = m(2,3)
m(3,3) = matin(6, isg)
endif
!
do id =1,ndegr !...(4)id =1,ndegr
!
unint = 0.d0
do iunk = 1,ndegr
unint(1) = unint(1) + m(id, iunk)*rhsel(iunk,isg)
enddo
!
!...step 2
!
unksgq(id, isg)= unint(1)
!
enddo !...(4)id =1,ndegr
!
enddo
!
end subroutine  getrhosubg3
!
!...Get the density distribution with curved subgrid...
!
subroutine  getrhosubg4(xpq,xpqi,unksgq, geoq_sub)
use constant
implicit none
!...Input
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:2, 1:npqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:nq), intent(out)::unksgq
real*8,dimension(1:7, 1:4), intent(out)::geoq_sub
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, ielem, id, isg, ideg, iunk,isc

integer,dimension(4, 4)::ipqsg, ipqsc
!...Local real array
real*8,dimension(1:ndegr,1:4)::rhsel
real*8,dimension(1:nmati, 1:4)::matin
real*8,dimension(1:2, 1:npqua)::xpqsg,xpqisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(3, 3)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: posqc(2, 12)
real*8:: xpqc(2, 12),xpqic(2, 12)
real*8:: xpqsc(2, npqua), xpqisc(2, npqua)
!...Local real
real*8:: xsc, ysc, xsci, ysci
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::b2,b3,b4,b5,b6
real*8::masel,xgaus,ygaus
!
real*8,allocatable::x5(:,:),b55(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
if(npoly==1) allocate(x5(3,3), mmatr(3,3), b55(3))
if(npoly==2) allocate(x5(5,5), mmatr(5,5), b55(5))
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
posqc(1, 1)= -0.5d0; posqc(2, 1)= -1.d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)= -1.d0;
posqc(1, 3)=  1.0d0; posqc(2, 3)= -.5d0;
posqc(1, 4)=  1.0d0; posqc(2, 4)=  .5d0;
posqc(1, 5)=  0.5d0; posqc(2, 5)=  1.d0;
posqc(1, 6)= -0.5d0; posqc(2, 6)=  1.d0;
posqc(1, 7)= -1.0d0; posqc(2, 7)=  .5d0;
posqc(1, 8)= -1.0d0; posqc(2, 8)= -.5d0;
posqc(1, 9)=  0.0d0; posqc(2, 9)= -.5d0;
posqc(1,10)=  0.5d0; posqc(2,10)=  0.d0;
posqc(1,11)=  0.0d0; posqc(2,11)=  .5d0;
posqc(1,12)= -0.5d0; posqc(2,12)=  0.d0;
!
ipqsc(1, 1)= 1; ipqsc(2, 1)= 9; ipqsc(3, 1)= 12; ipqsc(4, 1)= 8;
ipqsc(1, 2)= 2; ipqsc(2, 2)= 3; ipqsc(3, 2)= 10; ipqsc(4, 2)= 9;
ipqsc(1, 3)=10; ipqsc(2, 3)= 4; ipqsc(3, 3)=  5; ipqsc(4, 3)=11;
ipqsc(1, 4)=12; ipqsc(2, 4)=11; ipqsc(3, 4)=  6; ipqsc(4, 4)= 7;
!
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0

!...Now configuration...
xpqsg(1,1:9) =  xpq(1, 1:9)
xpqsg(2,1:9) =  xpq(2, 1:9)

!...Initial configuration...
xpqisg(1,1:9) =  xpqi(1, 1:9)
xpqisg(2,1:9) =  xpqi(2, 1:9)

!...Get the 12 high-order curved nodes
do isc = 1, 12
r = posqc(1,isc)
s = posqc(2,isc)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpqsg(1,ishp)
ysc = ysc + shpq(ishp)*xpqsg(2,ishp)
enddo
!
xsci = 0.d0
ysci = 0.d0
!
do ishp = 1, npqua
xsci = xsci + shpq(ishp)*xpqisg(1,ishp)
ysci = ysci + shpq(ishp)*xpqisg(2,ishp)
enddo
!
xpqc(1, isc) = xsc
xpqc(2, isc) = ysc
!
xpqic(1, isc) = xsci
xpqic(2, isc) = ysci
enddo
!
!print*,'test getrhosubg4'
!
!...Part I: Get cell averaged r and s for one subgrid...
!
do isg = 1, 4

!...Coordinates
xpqsc(1, 1:4) = xpqsg(1, ipqsg(1:4, isg))
xpqsc(2, 1:4) = xpqsg(2, ipqsg(1:4, isg))

xpqsc(1, 5:8) = xpqc(1, ipqsc(1:4, isg))
xpqsc(2, 5:8) = xpqc(2, ipqsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initial coordinates
xpqisc(1, 1:4) = xpqisg(1, ipqsg(1:4, isg))
xpqisc(2, 1:4) = xpqisg(2, ipqsg(1:4, isg))

xpqisc(1, 5:8) = xpqic(1, ipqsc(1:4, isg))
xpqisc(2, 5:8) = xpqic(2, ipqsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))

!...Initialze parameters...
mmatr = 0.d0
!
volel = 0.d0
xcel = 0.d0
ycel = 0.d0
!
voleli = 0.d0
xceli = 0.d0
yceli = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, npqua
dxdri = dxdri + dsprq(ishp)*xpqisc(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqisc(1,ishp)

dydri = dydri + dsprq(ishp)*xpqisc(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqisc(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)
!
!...Parameters for now configuration...
!
volel = volel + djaco
xcel = xcel + djaco*r
ycel = ycel + djaco*s
!
!...Parameters for initial configuration...
!
voleli = voleli + djacoi
xceli = xceli + djacoi*r
yceli = yceli + djacoi*s
!
enddo
!
!...Volume center for subgrids...
!
geoq_sub(1, isg) = xcel/volel
geoq_sub(2, isg) = ycel/volel
geoq_sub(3, isg) = xceli/voleli
geoq_sub(4, isg) = yceli/voleli
!
enddo
!
!...Part II: Get mass matrix for four subgrids...
!
do isg = 1, 4
!...Coordinates
xpqsc(1, 1:4) = xpqsg(1, ipqsg(1:4, isg))
xpqsc(2, 1:4) = xpqsg(2, ipqsg(1:4, isg))

xpqsc(1, 5:8) = xpqc(1, ipqsc(1:4, isg))
xpqsc(2, 5:8) = xpqc(2, ipqsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initialze parameters...
mmatr = 0.d0
!
rc = geoq_sub(1, isg)
sc = geoq_sub(2, isg)
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
mmatr(1,1) = mmatr(1,1) + djaco

mmatr(2,2) = mmatr(2,2) + b2*b2*djaco
mmatr(3,2) = mmatr(3,2) + b2*b3*djaco

mmatr(3,3) = mmatr(3,3) + b3*b3*djaco
!
enddo

!...Mass matrix
mmatr(2, 1) = 0.d0
mmatr(3, 1) = 0.d0
mmatr(1, 2) = mmatr(2, 1)
mmatr(1, 3) = mmatr(3, 1)
mmatr(2, 3) = mmatr(3, 2)
!

if(npoly==1)then
x5 = 0.d0
b55 = 0.d0
call getinvmat(3, mmatr, x5, b55)
!
matin(1,isg) = x5(1,1)
matin(2,isg) = x5(1,2)
matin(3,isg) = x5(1,3)
matin(4,isg) = x5(2,2)
matin(5,isg) = x5(2,3)
matin(6,isg) = x5(3,3)
!
endif
enddo
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, 4

!...Initial coordinates
xpqisc(1, 1:4) = xpqisg(1, ipqsg(1:4, isg))
xpqisc(2, 1:4) = xpqisg(2, ipqsg(1:4, isg))

xpqisc(1, 5:8) = xpqic(1, ipqsc(1:4, isg))
xpqisc(2, 5:8) = xpqic(2, ipqsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))

!
rc = geoq_sub(3, isg)
sc = geoq_sub(4, isg)
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqisc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqisc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqisc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqisc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
xg = r
yg = s
!
b(1) = 1.d0
b(2) = (xg-rc)/dr
b(3) = (yg-sc)/ds
!
call getrhoig_quadcurv(rhoi, xpqi)
!
do ideg = 1,ndegr
rhsel(ideg, isg)=rhsel(ideg, isg) + rhoi*b(ideg)*djaco
enddo
enddo
!
if(npoly==1)then
m(1,1) = matin(1, isg)
m(1,2) = 0.d0!matin(2, ie)
m(1,3) = 0.d0!matin(3, ie)

m(2,1) = m(1,2)
m(2,2) = matin(4, isg)
m(2,3) = matin(5, isg)

m(3,1) = m(1,3)
m(3,2) = m(2,3)
m(3,3) = matin(6, isg)
endif

!
do id =1,ndegr
unint = 0.d0
do iunk = 1,ndegr
unint(1) = unint(1) + m(id, iunk)*rhsel(iunk,isg)
enddo

!...Get the density for every subgrid
unksgq(id, isg)= unint(1)
enddo
enddo
!
end subroutine  getrhosubg4
!
!...Get the density distribution with curved subgrid...
!
subroutine  getrhosubg5(xpq,xpqi,unksgq, geoq_sub)
use constant
implicit none
!...Input
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:2, 1:npqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:nq), intent(out)::unksgq
real*8,dimension(1:7, 1:4), intent(out)::geoq_sub
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, ielem, id, isg, ideg, iunk,isc

integer,dimension(4, 4)::ipqsg, ipqsc
!...Local real array
real*8,dimension(1:ndegr,1:4)::rhsel
real*8,dimension(1:nmatr, 1:4)::matin
real*8,dimension(1:2, 1:npqua)::xpqsg,xpqisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: posqc(2, 12)
real*8:: xpqc(2, 12),xpqic(2, 12)
real*8:: xpqsc(2, npqua), xpqisc(2, npqua)
!...Local real
real*8:: xsc, ysc, xsci, ysci
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::b2,b3,b4,b5,b6
real*8::bq22,bq33,bq23
real*8::f0,f1,f2,f3
real*8:: f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::masel,xgaus,ygaus
real*8::rcsg,scsg
real*8::bq(ndegr)
!
real*8,allocatable::mminv(:,:),binv(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
if(npoly==2) allocate(mminv(5,5), mmatr(5,5), binv(5))
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
posqc(1, 1)= -0.5d0; posqc(2, 1)= -1.d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)= -1.d0;
posqc(1, 3)=  1.0d0; posqc(2, 3)= -.5d0;
posqc(1, 4)=  1.0d0; posqc(2, 4)=  .5d0;
posqc(1, 5)=  0.5d0; posqc(2, 5)=  1.d0;
posqc(1, 6)= -0.5d0; posqc(2, 6)=  1.d0;
posqc(1, 7)= -1.0d0; posqc(2, 7)=  .5d0;
posqc(1, 8)= -1.0d0; posqc(2, 8)= -.5d0;
posqc(1, 9)=  0.0d0; posqc(2, 9)= -.5d0;
posqc(1,10)=  0.5d0; posqc(2,10)=  0.d0;
posqc(1,11)=  0.0d0; posqc(2,11)=  .5d0;
posqc(1,12)= -0.5d0; posqc(2,12)=  0.d0;
!
ipqsc(1, 1)= 1; ipqsc(2, 1)= 9; ipqsc(3, 1)= 12; ipqsc(4, 1)= 8;
ipqsc(1, 2)= 2; ipqsc(2, 2)= 3; ipqsc(3, 2)= 10; ipqsc(4, 2)= 9;
ipqsc(1, 3)=10; ipqsc(2, 3)= 4; ipqsc(3, 3)=  5; ipqsc(4, 3)=11;
ipqsc(1, 4)=12; ipqsc(2, 4)=11; ipqsc(3, 4)=  6; ipqsc(4, 4)= 7;
!
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0

!...Now configuration...
xpqsg(1,1:9) =  xpq(1, 1:9)
xpqsg(2,1:9) =  xpq(2, 1:9)

!...Initial configuration...
xpqisg(1,1:9) =  xpqi(1, 1:9)
xpqisg(2,1:9) =  xpqi(2, 1:9)

!...Get the 12 high-order curved nodes
do isc = 1, 12
r = posqc(1,isc)
s = posqc(2,isc)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpqsg(1,ishp)
ysc = ysc + shpq(ishp)*xpqsg(2,ishp)
enddo
!
xsci = 0.d0
ysci = 0.d0
!
do ishp = 1, npqua
xsci = xsci + shpq(ishp)*xpqisg(1,ishp)
ysci = ysci + shpq(ishp)*xpqisg(2,ishp)
enddo
!
xpqc(1, isc) = xsc
xpqc(2, isc) = ysc
!
xpqic(1, isc) = xsci
xpqic(2, isc) = ysci
enddo
!
!print*,'test getrhosubg4'
!
!...Part I: Get cell averaged r and s for one subgrid...
!
do isg = 1, 4

!...Coordinates
xpqsc(1, 1:4) = xpqsg(1, ipqsg(1:4, isg))
xpqsc(2, 1:4) = xpqsg(2, ipqsg(1:4, isg))

xpqsc(1, 5:8) = xpqc(1, ipqsc(1:4, isg))
xpqsc(2, 5:8) = xpqc(2, ipqsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initial coordinates
xpqisc(1, 1:4) = xpqisg(1, ipqsg(1:4, isg))
xpqisc(2, 1:4) = xpqisg(2, ipqsg(1:4, isg))

xpqisc(1, 5:8) = xpqic(1, ipqsc(1:4, isg))
xpqisc(2, 5:8) = xpqic(2, ipqsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))

!...Initialze parameters...
if(npoly.eq.2) mmatr = 0.d0
!
volel = 0.d0
xcel = 0.d0
ycel = 0.d0
!
voleli = 0.d0
xceli = 0.d0
yceli = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, npqua
dxdri = dxdri + dsprq(ishp)*xpqisc(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqisc(1,ishp)

dydri = dydri + dsprq(ishp)*xpqisc(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqisc(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)
!
!...Parameters for now configuration...
!
volel = volel + djaco
xcel = xcel + djaco*r
ycel = ycel + djaco*s
!
!...Parameters for initial configuration...
!
voleli = voleli + djacoi
xceli = xceli + djacoi*r
yceli = yceli + djacoi*s
!
enddo
!
!...Volume center for subgrids...
!
geoq_sub(1, isg) = xcel/volel
geoq_sub(2, isg) = ycel/volel
geoq_sub(3, isg) = xceli/voleli
geoq_sub(4, isg) = yceli/voleli
!
enddo
!
!...Part II: Get high-order terms for basis function...
!
if(npoly.eq.2)then

do isg = 1, 4

!...Coordinates
xpqsc(1, 1:4) = xpqsg(1, ipqsg(1:4, isg))
xpqsc(2, 1:4) = xpqsg(2, ipqsg(1:4, isg))

xpqsc(1, 5:8) = xpqc(1, ipqsc(1:4, isg))
xpqsc(2, 5:8) = xpqc(2, ipqsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initial coordinates
xpqisc(1, 1:4) = xpqisg(1, ipqsg(1:4, isg))
xpqisc(2, 1:4) = xpqisg(2, ipqsg(1:4, isg))

xpqisc(1, 5:8) = xpqic(1, ipqsc(1:4, isg))
xpqisc(2, 5:8) = xpqic(2, ipqsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))
!
rcsg = geoq_sub(1, isg)
scsg = geoq_sub(2, isg)
!...Initialze parameters...
mmatr = 0.d0
!...Zeo our volel
volel = 0.d0

!...Zero out bq22,bq33,bq23
bq22 = 0.d0
bq33 = 0.d0
bq23 = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
volel = volel + djaco

!...Basis function
bq(1) = 1.d0
bq(2) = (r-rcsg)/dr
bq(3) = (s-scsg)/ds
!
bq22 = bq22 + 0.5d0*bq(2)*bq(2)*djaco
bq33 = bq33 + 0.5d0*bq(3)*bq(3)*djaco
bq23 = bq23 + bq(2)*bq(3)*djaco
!
enddo
!
!...Volume center for subgrids...
!
geoq_sub(5, isg) = bq22/volel
geoq_sub(6, isg) = bq33/volel
geoq_sub(7, isg) = bq23/volel
!
enddo

endif !...if (npoly.eq.2)
!
!...Part II: Get mass matrix for four subgrids...
!
do isg = 1, 4
!...Coordinates
xpqsc(1, 1:4) = xpqsg(1, ipqsg(1:4, isg))
xpqsc(2, 1:4) = xpqsg(2, ipqsg(1:4, isg))

xpqsc(1, 5:8) = xpqc(1, ipqsc(1:4, isg))
xpqsc(2, 5:8) = xpqc(2, ipqsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initialze parameters...
if(npoly.eq.2)  mmatr = 0.d0
f0 = 0.d0
f1 = 0.d0
f2 = 0.d0
f3 = 0.d0

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
!
rc = geoq_sub(1, isg)
sc = geoq_sub(2, isg)
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
print*,'djaco',djaco,wi,(dxdr*dyds - dydr*dxds)
!
xg = r
yg = s

!Basis function
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds

!...DGP2
if(npoly.eq.2)then
b4 = 0.5d0*b2*b2 - geoq_sub(5, isg)
b5 = 0.5d0*b3*b3 - geoq_sub(6, isg)
b6 =       b2*b3 - geoq_sub(7, isg)
endif
!
f0 = f0 + djaco
f1 = f1 + b2*b2*djaco
f2 = f2 + b2*b3*djaco
f3 = f3 + b3*b3*djaco
!
if(npoly==2)then
f22 = f22 + b2*b2*djaco
f23 = f23 + b2*b3*djaco
f24 = f24 + b2*b4*djaco
f25 = f25 + b2*b5*djaco
f26 = f26 + b2*b6*djaco

f33 = f33 + b3*b3*djaco
f34 = f34 + b3*b4*djaco
f35 = f35 + b3*b5*djaco
f36 = f36 + b3*b6*djaco

f44 = f44 + b4*b4*djaco
f45 = f45 + b4*b5*djaco
f46 = f46 + b4*b6*djaco

f55 = f55 + b5*b5*djaco
f56 = f56 + b5*b6*djaco

f66 = f66 + b6*b6*djaco
endif
!
enddo

!
if(npoly==1)then
det = f1*f3 - f2**2

matin(1,isg) = f3/det
matin(2,isg) =-f2/det
matin(3,isg) = f1/det
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then

!...Mass matrix
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
binv = 0.d0
mminv = 0.d0
call getinvmat(5, mmatr, mminv, binv)
!
matin(1,isg) = mminv(1,1)
matin(2,isg) = mminv(1,2)
matin(3,isg) = mminv(1,3)
matin(4,isg) = mminv(1,4)
matin(5,isg) = mminv(1,5)

matin(6,isg) = mminv(2,2)
matin(7,isg) = mminv(2,3)
matin(8,isg) = mminv(2,4)
matin(9,isg) = mminv(2,5)

matin(10,isg) = mminv(3,3)
matin(11,isg) = mminv(3,4)
matin(12,isg) = mminv(3,5)

matin(13,isg) = mminv(4,4)
matin(14,isg) = mminv(4,5)

matin(15,isg) = mminv(5,5)
!
matin(16,isg) = 1.d0/f0
!
endif

enddo

!
print*,'f0',f0
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, 4

!...Initial coordinates
xpqisc(1, 1:4) = xpqisg(1, ipqsg(1:4, isg))
xpqisc(2, 1:4) = xpqisg(2, ipqsg(1:4, isg))

xpqisc(1, 5:8) = xpqic(1, ipqsc(1:4, isg))
xpqisc(2, 5:8) = xpqic(2, ipqsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))

!
rcsg = geoq_sub(3, isg)
scsg = geoq_sub(4, isg)
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqisc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqisc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqisc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqisc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
xg = r
yg = s
!
b(1) = 1.d0
b(2) = (xg-rcsg)/dr
b(3) = (yg-scsg)/ds

!...DGP2
if(npoly.eq.2)then
b(4) = 0.5d0*b(2)*b(2) - geoq_sub(5, isg)
b(5) = 0.5d0*b(3)*b(3) - geoq_sub(6, isg)
b(6) =       b(2)*b(3) - geoq_sub(7, isg)
endif

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpqisc(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqisc(2,ishp)
enddo
!
call getrhoig_quadcurv2(rhoi, xpqi, xgaus, ygaus)
!
do ideg = 1,ndegr
rhsel(ideg, isg)=rhsel(ideg, isg) + rhoi*b(ideg)*djaco
enddo
!if(isg.eq.4)print*,'isgideg',rhsel(6, isg),b(6),b(2),b(3), geoq_sub(7, isg)
enddo
!
if(npoly==1)then
m(1,1) = matin(4, isg)
m(1,2) = 0.d0
m(1,3) = 0.d0

m(2,1) = m(1,2)
m(2,2) = matin(1, isg)
m(2,3) = matin(2, isg)

m(3,1) = m(1,3)
m(3,2) = matin(2, isg)
m(3,3) = matin(3, isg)

elseif(npoly==2)then
m(1,1) = matin(16, isg)
m(1,2) = 0.d0
m(1,3) = 0.d0
m(1,4) = 0.d0
m(1,5) = 0.d0
m(1,6) = 0.d0

m(2,1) = m(1,2)
m(2,2) = matin(1, isg)
m(2,3) = matin(2, isg)
m(2,4) = matin(3, isg)
m(2,5) = matin(4, isg)
m(2,6) = matin(5, isg)

m(3,1) = m(1,3)
m(3,2) = m(2,3)
m(3,3) = matin(6, isg)
m(3,4) = matin(7, isg)
m(3,5) = matin(8, isg)
m(3,6) = matin(9, isg)

m(4,1) = m(1,4)
m(4,2) = m(2,4)
m(4,3) = m(3,4)
m(4,4) = matin(10,isg)
m(4,5) = matin(11,isg)
m(4,6) = matin(12,isg)

m(5,1) = m(1,5)
m(5,2) = m(2,5)
m(5,3) = m(3,5)
m(5,4) = m(4,5)
m(5,5) = matin(13,isg)
m(5,6) = matin(14,isg)

m(6,1) = m(1,6)
m(6,2) = m(2,6)
m(6,3) = m(3,6)
m(6,4) = m(4,6)
m(6,5) = m(5,6)
m(6,6) = matin(15,isg)

endif

!...Update the subcell density ditribution
do id =1,ndegr
unint = 0.d0
do iunk = 1,ndegr
unint(1) = unint(1) + m(id, iunk)*rhsel(iunk,isg)
enddo

!if(isg.eq.4.and.id.eq.1) print*,'isg4',m(1, :),rhsel(:,4),unint(1)
unksgq(id, isg)= unint(1)
enddo
enddo
!
!print*,'rhosubg',m(1, 1),rhsel(1,1),m(1, 1)*rhsel(1,1)
!
end subroutine  getrhosubg5
!
!...Get the averaged-density within one sub-cell...
!
subroutine  getrhosubg_averg(xpq,xpqi,unksgq, geoq_sub, ielem)
use constant
implicit none
!...Input
integer, intent(in)::ielem
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:2, 1:npqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:nq), intent(out)::unksgq
real*8,dimension(1:7, 1:4), intent(out)::geoq_sub
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc

integer,dimension(4, 4)::ipqsg, ipqsc
!...Local real array
real*8,dimension(1:ndegr,1:4)::rhsel
real*8,dimension(1:nmatr, 1:4)::matin
real*8,dimension(1:2, 1:npqua)::xpqsg,xpqisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: posqc(2, 12)
real*8:: xpqc(2, 12),xpqic(2, 12)
real*8:: xpqsc(2, npqua), xpqisc(2, npqua)
!...Local real
real*8:: xsc, ysc, xsci, ysci
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::b2,b3,b4,b5,b6
real*8::bq22,bq33,bq23
real*8::f0,f1,f2,f3
real*8:: f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::masel,xgaus,ygaus
real*8::rcsg,scsg
real*8::bq(ndegr)
!
real*8,allocatable::mminv(:,:),binv(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
if(npoly==2) allocate(mminv(5,5), mmatr(5,5), binv(5))
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
posqc(1, 1)= -0.5d0; posqc(2, 1)= -1.d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)= -1.d0;
posqc(1, 3)=  1.0d0; posqc(2, 3)= -.5d0;
posqc(1, 4)=  1.0d0; posqc(2, 4)=  .5d0;
posqc(1, 5)=  0.5d0; posqc(2, 5)=  1.d0;
posqc(1, 6)= -0.5d0; posqc(2, 6)=  1.d0;
posqc(1, 7)= -1.0d0; posqc(2, 7)=  .5d0;
posqc(1, 8)= -1.0d0; posqc(2, 8)= -.5d0;
posqc(1, 9)=  0.0d0; posqc(2, 9)= -.5d0;
posqc(1,10)=  0.5d0; posqc(2,10)=  0.d0;
posqc(1,11)=  0.0d0; posqc(2,11)=  .5d0;
posqc(1,12)= -0.5d0; posqc(2,12)=  0.d0;
!
ipqsc(1, 1)= 1; ipqsc(2, 1)= 9; ipqsc(3, 1)= 12; ipqsc(4, 1)= 8;
ipqsc(1, 2)= 2; ipqsc(2, 2)= 3; ipqsc(3, 2)= 10; ipqsc(4, 2)= 9;
ipqsc(1, 3)=10; ipqsc(2, 3)= 4; ipqsc(3, 3)=  5; ipqsc(4, 3)=11;
ipqsc(1, 4)=12; ipqsc(2, 4)=11; ipqsc(3, 4)=  6; ipqsc(4, 4)= 7;
!
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0

!...Now configuration...
xpqsg(1,1:9) =  xpq(1, 1:9)
xpqsg(2,1:9) =  xpq(2, 1:9)

!...Initial configuration...
xpqisg(1,1:9) =  xpqi(1, 1:9)
xpqisg(2,1:9) =  xpqi(2, 1:9)

!...Get the 12 high-order curved nodes
do isc = 1, 12
r = posqc(1,isc)
s = posqc(2,isc)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpqsg(1,ishp)
ysc = ysc + shpq(ishp)*xpqsg(2,ishp)
enddo
!
xsci = 0.d0
ysci = 0.d0
!
do ishp = 1, npqua
xsci = xsci + shpq(ishp)*xpqisg(1,ishp)
ysci = ysci + shpq(ishp)*xpqisg(2,ishp)
enddo
!
xpqc(1, isc) = xsc
xpqc(2, isc) = ysc
!
xpqic(1, isc) = xsci
xpqic(2, isc) = ysci
enddo
!
!print*,'test getrhosubg4'
!
!...Part I: Get cell averaged r and s for one subgrid...
!

!
!...Part II: Get high-order terms for basis function...
!

!
!...Part II: Get mass matrix for four subgrids...
!
do isg = 1, 4
!...Coordinates
xpqsc(1, 1:4) = xpqsg(1, ipqsg(1:4, isg))
xpqsc(2, 1:4) = xpqsg(2, ipqsg(1:4, isg))

xpqsc(1, 5:8) = xpqc(1, ipqsc(1:4, isg))
xpqsc(2, 5:8) = xpqc(2, ipqsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initialze parameters...
if(npoly.eq.2)  mmatr = 0.d0
f0 = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco
!
!if(ielem.eq.1) print*,'djaco',f0,wi,dxdr,dyds,dydr,dxds,xpqsc(1:2,:)
enddo

!
if(npoly==0)then
matin(1,isg) = 1.d0/f0
elseif(npoly==1)then
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then
matin(16,isg) = 1.d0/f0
endif

enddo
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, 4

!...Initial coordinates
xpqisc(1, 1:4) = xpqisg(1, ipqsg(1:4, isg))
xpqisc(2, 1:4) = xpqisg(2, ipqsg(1:4, isg))

xpqisc(1, 5:8) = xpqic(1, ipqsc(1:4, isg))
xpqisc(2, 5:8) = xpqic(2, ipqsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqisc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqisc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqisc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqisc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
b(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpqisc(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqisc(2,ishp)
enddo
!
call getrhoig_quadcurv2(rhoi, xpqi, xgaus, ygaus)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*b(1)*djaco
!if(isg.eq.4)print*,'isgideg',rhsel(6, isg),b(6),b(2),b(3), geoq_sub(7, isg)
enddo
!
if(npoly==0)then
m(1,1) = matin(1, isg)
elseif(npoly==1)then
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
end subroutine  getrhosubg_averg
!
!...Get the averaged-density within one sub-cell using SMS with MEM...
!
subroutine  getrhosmsgd_averg(rci,sci,bqhoi,xpq,xpqi,unksgq, unkgd, geoq_sub, ielem)
use constant
implicit none
!...Input
real*8, intent(in)::rci, sci
real*8,dimension(3),            intent(in):: bqhoi
integer, intent(in)::ielem
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:2, 1:npqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:nq), intent(out)::unksgq
real*8,dimension(1:ndegr,1:4),  intent(in)::unkgd
real*8,dimension(1:7, 1:4), intent(out)::geoq_sub
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc

integer,dimension(4, 4)::ipqsg, ipqsc
!...Local real array
real*8,dimension(1:ndegr,1:4)::rhsel
real*8,dimension(1:nmatr, 1:4)::matin
real*8,dimension(1:2, 1:npqua)::xpqsg,xpqisg
real*8::b(ndegr),bi(ndegr)
real*8::jacbf(2,2)
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: posqc(2, 12)
real*8:: xpqc(2, 12),xpqic(2, 12)
real*8:: xpqsc(2, npqua), xpqisc(2, npqua)
!...Local real
real*8:: xsc, ysc, xsci, ysci
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::b2,b3,b4,b5,b6
real*8::bq22,bq33,bq23
real*8::f0,f1,f2,f3
real*8:: f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::a11,a22,a12,a21
real*8::masel,xgaus,ygaus
real*8::rcsg,scsg
real*8::bq(ndegr)
!
real*8,allocatable::mminv(:,:),binv(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
if(npoly==2) allocate(mminv(5,5), mmatr(5,5), binv(5))
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
posqc(1, 1)= -0.5d0; posqc(2, 1)= -1.d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)= -1.d0;
posqc(1, 3)=  1.0d0; posqc(2, 3)= -.5d0;
posqc(1, 4)=  1.0d0; posqc(2, 4)=  .5d0;
posqc(1, 5)=  0.5d0; posqc(2, 5)=  1.d0;
posqc(1, 6)= -0.5d0; posqc(2, 6)=  1.d0;
posqc(1, 7)= -1.0d0; posqc(2, 7)=  .5d0;
posqc(1, 8)= -1.0d0; posqc(2, 8)= -.5d0;
posqc(1, 9)=  0.0d0; posqc(2, 9)= -.5d0;
posqc(1,10)=  0.5d0; posqc(2,10)=  0.d0;
posqc(1,11)=  0.0d0; posqc(2,11)=  .5d0;
posqc(1,12)= -0.5d0; posqc(2,12)=  0.d0;
!
ipqsc(1, 1)= 1; ipqsc(2, 1)= 9; ipqsc(3, 1)= 12; ipqsc(4, 1)= 8;
ipqsc(1, 2)= 2; ipqsc(2, 2)= 3; ipqsc(3, 2)= 10; ipqsc(4, 2)= 9;
ipqsc(1, 3)=10; ipqsc(2, 3)= 4; ipqsc(3, 3)=  5; ipqsc(4, 3)=11;
ipqsc(1, 4)=12; ipqsc(2, 4)=11; ipqsc(3, 4)=  6; ipqsc(4, 4)= 7;
!
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0

!...Now configuration...
xpqsg(1,1:9) =  xpq(1, 1:9)
xpqsg(2,1:9) =  xpq(2, 1:9)

!...Initial configuration...
xpqisg(1,1:9) =  xpqi(1, 1:9)
xpqisg(2,1:9) =  xpqi(2, 1:9)

!...Get the 12 high-order curved nodes
do isc = 1, 12
r = posqc(1,isc)
s = posqc(2,isc)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpqsg(1,ishp)
ysc = ysc + shpq(ishp)*xpqsg(2,ishp)
enddo
!
xsci = 0.d0
ysci = 0.d0
!
do ishp = 1, npqua
xsci = xsci + shpq(ishp)*xpqisg(1,ishp)
ysci = ysci + shpq(ishp)*xpqisg(2,ishp)
enddo
!
xpqc(1, isc) = xsc
xpqc(2, isc) = ysc
!
xpqic(1, isc) = xsci
xpqic(2, isc) = ysci
enddo
!
!print*,'test getrhosubg4'
!
!...Part I: Get cell averaged r and s for one subgrid...
!

!
!...Part II: Get high-order terms for basis function...
!

!
!...Part II: Get mass matrix for four subgrids...
!
do isg = 1, 4
!...Coordinates
xpqsc(1, 1:4) = xpqisg(1, ipqsg(1:4, isg))
xpqsc(2, 1:4) = xpqisg(2, ipqsg(1:4, isg))

xpqsc(1, 5:8) = xpqic(1, ipqsc(1:4, isg))
xpqsc(2, 5:8) = xpqic(2, ipqsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initialze parameters...
if(npoly.eq.2)  mmatr = 0.d0
f0 = 0.d0
!
do ig =1,ngausdq
!
if(isg.eq.1)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.2)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.3)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
elseif(isg.eq.4)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
endif
!
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!...Jacobian transformation matrix
jacbf = 0.d0
!
bi(1) = 1.d0
bi(2) = (r - rci)/dr
bi(3) = (s - sci)/ds
!DGP2
if(npoly.eq.2)then
bi(4) = 0.5d0*bi(2)*bi(2)  - bqhoi(1)
bi(5) = 0.5d0*bi(3)*bi(3)  - bqhoi(2)
bi(6) =       bi(2)*bi(3)  - bqhoi(3)
endif
!
do ideg = 1, ndegr
jacbf(1,1) = jacbf(1, 1) + unkgd(ideg, 1)*bi(ideg)
jacbf(1,2) = jacbf(1, 2) + unkgd(ideg, 2)*bi(ideg)
jacbf(2,1) = jacbf(2, 1) + unkgd(ideg, 3)*bi(ideg)
jacbf(2,2) = jacbf(2, 2) + unkgd(ideg, 4)*bi(ideg)
enddo
!
a11 = jacbf(1, 1)*dxdr + jacbf(1, 2)*dydr
a12 = jacbf(1, 1)*dxds + jacbf(1, 2)*dyds

a21 = jacbf(2, 1)*dxdr + jacbf(2, 2)*dydr
a22 = jacbf(2, 1)*dxds + jacbf(2, 2)*dyds
!
djaco = wi*(a11*a22 - a12*a21)!*0.25d0
!
f0 = f0 + djaco
!
!if(ielem.eq.1) print*,'djaco',f0,wi,dxdr,dyds,dydr,dxds,xpqsc(1:2,:)
enddo

!
if(npoly==1)then
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then
matin(16,isg) = 1.d0/f0
endif

enddo
!
!if(ielem.eq.1) print*,'djaco',matin(16,:)
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, 4

!...Initial coordinates
xpqisc(1, 1:4) = xpqisg(1, ipqsg(1:4, isg))
xpqisc(2, 1:4) = xpqisg(2, ipqsg(1:4, isg))

xpqisc(1, 5:8) = xpqic(1, ipqsc(1:4, isg))
xpqisc(2, 5:8) = xpqic(2, ipqsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqisc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqisc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqisc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqisc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
b(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpqisc(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqisc(2,ishp)
enddo
!
call getrhoig_quadcurv2(rhoi, xpqi, xgaus, ygaus)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*b(1)*djaco
!if(isg.eq.4)print*,'isgideg',rhsel(6, isg),b(6),b(2),b(3), geoq_sub(7, isg)
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
end subroutine  getrhosmsgd_averg

!
!...Get the averaged-density within one sub-cell for triangles...
!
subroutine  getrhosubgt_averg(xpt,xpti,unksgt, ielem)
use constant
implicit none
!...Input
integer, intent(in)::ielem
real*8,dimension(1:2, 1:nptri), intent(in)::xpt
real*8,dimension(1:2, 1:nptri), intent(in)::xpti
real*8,dimension(1:ndegr, 1:3), intent(out)::unksgt
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc

integer,dimension(4, 3)::iptsg, iptsc
!...Local real array
real*8,dimension(1:ndegr,1:3)::rhsel
real*8,dimension(1:nmatr, 1:3)::matin
real*8,dimension(1:2, 1:nptri+1)::xptsg,xptisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:nptri)::shpt, dsprt, dspst
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: postc(2, 9)
real*8:: xptc(2, 9),xptic(2, 9)
real*8:: xpqsc(2, npqua), xpqisc(2, npqua)
!...Local real
real*8:: xsc, ysc, xsci, ysci
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::f0,f1,f2,f3
real*8::masel,xgaus,ygaus
!
!-xxx-real contant
data c10 / 1.0d0 /
!
iptsg(1, 1) = 1; iptsg(2, 1) = 4; iptsg(3, 1) = 7; iptsg(4, 1) = 6
iptsg(1, 2) = 4; iptsg(2, 2) = 2; iptsg(3, 2) = 5; iptsg(4, 2) = 7
iptsg(1, 3) = 7; iptsg(2, 3) = 5; iptsg(3, 3) = 3; iptsg(4, 3) = 6
!
postc(1, 1)=  0.25d0; postc(2, 1)=   0.d0;
postc(1, 2)=  0.75d0; postc(2, 2)=   0.d0;
postc(1, 3)=  0.75d0; postc(2, 3)= 0.25d0;
postc(1, 4)=  0.25d0; postc(2, 4)= 0.75d0;
postc(1, 5)=  0.d0;   postc(2, 5)= 0.75d0;
postc(1, 6)=  0.d0;   postc(2, 6)= 0.25d0;
postc(1, 7)=  5.d0/12.d0;   postc(2, 7)= 1.d0/6.d0;
postc(1, 8)=  5.d0/12.d0;   postc(2, 8)= 5.d0/12.d0;
postc(1, 9)=  1.d0/6.d0;    postc(2, 9)= 5.d0/12.d0;
!
iptsc(1, 1)= 1; iptsc(2, 1)= 7; iptsc(3, 1)= 9; iptsc(4, 1)= 6;
iptsc(1, 2)= 2; iptsc(2, 2)= 3; iptsc(3, 2)= 8; iptsc(4, 2)= 7;
iptsc(1, 3)= 8; iptsc(2, 3)= 4; iptsc(3, 3)= 5; iptsc(4, 3)= 9;
!
call ruqope(2, ngausdq, posiq, weighq)

!...Now configuration...
xptsg(1,1:6) =  xpt(1, 1:6)
xptsg(2,1:6) =  xpt(2, 1:6)

!...Initial configuration...
xptisg(1,1:6) =  xpti(1, 1:6)
xptisg(2,1:6) =  xpti(2, 1:6)

!...Get the 12 high-order curved nodes
do isc = 1, 9
r = postc(1,isc)
s = postc(2,isc)

!...  shape function
shpt(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shpt(2) = -r*(c10-2.d0*r)
shpt(3) = -s*(c10-2.d0*s)
shpt(4) = 4.d0*r*(c10-r-s)
shpt(5) = 4.d0*r*s
shpt(6) = 4.d0*s*(c10-r-s)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, nptri
xsc = xsc + shpt(ishp)*xptsg(1,ishp)
ysc = ysc + shpt(ishp)*xptsg(2,ishp)
enddo
!
xsci = 0.d0
ysci = 0.d0
!
do ishp = 1, nptri
xsci = xsci + shpt(ishp)*xptisg(1,ishp)
ysci = ysci + shpt(ishp)*xptisg(2,ishp)
enddo
!
xptc(1, isc) = xsc
xptc(2, isc) = ysc
!
xptic(1, isc) = xsci
xptic(2, isc) = ysci
enddo
!
!...Part I: Construct the nodal coordinates
!
xptsg(1, 1:6) = xpt(1, 1:6)
xptsg(2, 1:6) = xpt(2, 1:6)
!
xptsg(1, 7) = -1.d0/9.d0*xpt(1, 1)-1.d0/9.d0*xpt(1, 2)-1.d0/9.d0*xpt(1, 3) +&
               4.d0/9.d0*xpt(1, 4)+4.d0/9.d0*xpt(1, 5)+4.d0/9.d0*xpt(1, 6)

xptsg(2, 7) = -1.d0/9.d0*xpt(2, 1)-1.d0/9.d0*xpt(2, 2)-1.d0/9.d0*xpt(2, 3) +&
               4.d0/9.d0*xpt(2, 4)+4.d0/9.d0*xpt(2, 5)+4.d0/9.d0*xpt(2, 6)
!
xptisg(1, 1:6) = xpti(1, 1:6)
xptisg(2, 1:6) = xpti(2, 1:6)
!
xptisg(1, 7) = -1.d0/9.d0*xpti(1, 1)-1.d0/9.d0*xpti(1, 2)-1.d0/9.d0*xpti(1, 3) +&
                4.d0/9.d0*xpti(1, 4)+4.d0/9.d0*xpti(1, 5)+4.d0/9.d0*xpti(1, 6)

xptisg(2, 7) = -1.d0/9.d0*xpti(2, 1)-1.d0/9.d0*xpti(2, 2)-1.d0/9.d0*xpti(2, 3) +&
                4.d0/9.d0*xpti(2, 4)+4.d0/9.d0*xpti(2, 5)+4.d0/9.d0*xpti(2, 6)

!
!...Part II: Get mass matrix for four subgrids...
!
do isg = 1, 3
!...Coordinates
xpqsc(1, 1:4) = xptsg(1, iptsg(1:4, isg))
xpqsc(2, 1:4) = xptsg(2, iptsg(1:4, isg))

xpqsc(1, 5:8) = xptc(1, iptsc(1:4, isg))
xpqsc(2, 5:8) = xptc(2, iptsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initialze parameters...
f0 = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco
!
!if(ielem.eq.1) print*,'djaco',f0,wi,dxdr,dyds,dydr,dxds,xpqsc(1:2,:)
enddo

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

do isg = 1, 3

!...Initial coordinates
xpqisc(1, 1:4) = xptisg(1, iptsg(1:4, isg))
xpqisc(2, 1:4) = xptisg(2, iptsg(1:4, isg))

xpqisc(1, 5:8) = xptic(1, iptsc(1:4, isg))
xpqisc(2, 5:8) = xptic(2, iptsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqisc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqisc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqisc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqisc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
b(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpqisc(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqisc(2,ishp)
enddo
!
call getrhoig_triacurv2(rhoi, xpti, xgaus, ygaus)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*b(1)*djaco
!if(isg.eq.4)print*,'isgideg',rhsel(6, isg),b(6),b(2),b(3), geoq_sub(7, isg)
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

!if(isg.eq.4.and.id.eq.1) print*,'isg4',m(1, :),rhsel(:,4),unint(1)
unksgt(1, isg)= unint(1)
enddo
!
end subroutine  getrhosubgt_averg
!
!...Face integral (mass center) for hybrid quad...
!
subroutine rhsifacedg_lagsubgt2(iptri, unkno, ustar,fstart, gesgt, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:4, 1:3, 1:ntria),  intent(in)::fstart !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac), intent(in)::gesgt
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem,isg,ivsg,ifsg
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(8, 3)::fntsg
integer,dimension(4, 3)::iptsg
real*8,dimension(2,4,3)::wfgst,wfgstm
real*8,dimension(1:3,1:2,1:4)::vnorm
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::xvt(nvtri+1), yvt(nvtri+1),bt(1:ndegr,1:nvtri)

!...local real number
real*8::eps,c00,c05,c10,c20,c13,c16
real*8::dr,ds,rc,sc
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
data c16   / 0.1666666666666666d0 /
!
!...Zero out plnpn, ulnpn
!
iptsg(1, 1) = 1; iptsg(2, 1) = 4; iptsg(3, 1) = 0; iptsg(4, 1) = 6
iptsg(1, 2) = 4; iptsg(2, 2) = 2; iptsg(3, 2) = 5; iptsg(4, 2) = 0
iptsg(1, 3) = 0; iptsg(2, 3) = 5; iptsg(3, 3) = 3; iptsg(4, 3) = 6
!
fntsg(1, 1) =  6; fntsg(2, 1) = 1;
fntsg(3, 1) = 10; fntsg(4, 1) =-7;
fntsg(5, 1) = -7; fntsg(6, 1) = 9;
fntsg(7, 1) = 9;  fntsg(8, 1) = 12;

fntsg(1, 2) =  7; fntsg(2, 2) =10;
fntsg(3, 2) = 2;  fntsg(4, 2) =  3;
fntsg(5, 2) = 11; fntsg(6, 2) =-8;
fntsg(7, 2) = -8; fntsg(8, 2) = 7;

fntsg(1, 3) =-9; fntsg(2, 3) = 8;
fntsg(3, 3) = 8; fntsg(4, 3) =11;
fntsg(5, 3) = 4; fntsg(6, 3) = 5;
fntsg(7, 3) =12; fntsg(8, 3) =-9;

!...For momentum and energy eqaution
wfgst(1, 1, 1) = 1.d0; wfgst(2, 1, 1) = 1.d0;
wfgst(1, 2, 1) = 1.d0; wfgst(2, 2, 1) = 1.d0;
wfgst(1, 3, 1) = 0.d0; wfgst(2, 3, 1) = 0.d0;
wfgst(1, 4, 1) = 1.d0; wfgst(2, 4, 1) = 1.d0;

wfgst(1, 1, 2) = 1.d0; wfgst(2, 1, 2) =1.d0;
wfgst(1, 2, 2) = 1.d0; wfgst(2, 2, 2) =1.d0;
wfgst(1, 3, 2) = 1.d0; wfgst(2, 3, 2) =1.d0;
wfgst(1, 4, 2) = 0.d0; wfgst(2, 4, 2) = 0.d0;

wfgst(1, 1, 3) = 0.d0; wfgst(2, 1, 3) = 0.d0;
wfgst(1, 2, 3) = 1.d0; wfgst(2, 2, 3) = 1.d0;
wfgst(1, 3, 3) = 1.d0; wfgst(2, 3, 3) = 1.d0;
wfgst(1, 4, 3) = 1.d0; wfgst(2, 4, 3) = 1.d0;

!...Conitnuity equation...
wfgstm(1, 1, 1) = 4.d0/6.d0; wfgstm(2, 1, 1) = 4.d0/6.d0;
wfgstm(1, 2, 1) = 4.d0/3.d0; wfgstm(2, 2, 1) = 0.d0;
wfgstm(1, 3, 1) = 0.d0;      wfgstm(2, 3, 1) = 0.d0;
wfgstm(1, 4, 1) = 0.d0;      wfgstm(2, 4, 1) = 4.d0/3.d0;

wfgstm(1, 1, 2) = 0.d0;      wfgstm(2, 1, 2) = 4.d0/3.d0;
wfgstm(1, 2, 2) = 4.d0/6.d0; wfgstm(2, 2, 2) = 4.d0/6.d0;
wfgstm(1, 3, 2) = 4.d0/3.d0; wfgstm(2, 3, 2) = 0.d0;
wfgstm(1, 4, 2) = 0.d0;      wfgstm(2, 4, 2) = 0.d0;

wfgstm(1, 1, 3) = 0.d0;      wfgstm(2, 1, 3) = 0.d0;
wfgstm(1, 2, 3) = 0.d0/3.d0; wfgstm(2, 2, 3) = 4.d0/3.d0;
wfgstm(1, 3, 3) = 4.d0/6.d0; wfgstm(2, 3, 3) = 4.d0/6.d0;
wfgstm(1, 4, 3) = 4.d0/3.d0; wfgstm(2, 4, 3) = 0.d0/3.d0;
!
wfgstm = 0.25d0*wfgstm
!
!...Quads...
!
do 650 ie = 1,ntria !...(1)ie = 1,nelem
!
ielem = ie
!
!...The vertex constituting one cell...
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
xvt(4) = 0.5d0; yvt(4) = 0.d0
xvt(5) = 0.5d0; yvt(5) = 0.5d0
xvt(6) = 0.d0;  yvt(6) = 0.5d0

!...The pseudo-point 7
xvt(7) = 1.d0/3.d0; yvt(7) = 1.d0/3.d0
!
!
do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bt(4, iv) = 0.5d0*bt(2, iv)*bt(2, iv) - geoel(19, ielem)
bt(5, iv) = 0.5d0*bt(3, iv)*bt(3, iv) - geoel(20, ielem)
bt(6, iv) =       bt(2, iv)*bt(3, iv) - geoel(21, ielem)
endif
enddo
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
do isg = 1, 3
!
!...Normal vector for every face...
!
vnorm(1:2, 1, 1) = sign(1, fntsg(1, isg))*gesgt(1:2,abs(fntsg(1, isg)), ie); vnorm(  3, 1, 1) = gesgt(3,abs(fntsg(1, isg)), ie);
vnorm(1:2, 2, 1) = sign(1, fntsg(2, isg))*gesgt(1:2,abs(fntsg(2, isg)), ie); vnorm(  3, 2, 1) = gesgt(3,abs(fntsg(2, isg)), ie);

vnorm(1:2, 1, 2) = sign(1, fntsg(3, isg))*gesgt(1:2,abs(fntsg(3, isg)), ie); vnorm(  3, 1, 2) = gesgt(3,abs(fntsg(3, isg)), ie);
vnorm(1:2, 2, 2) = sign(1, fntsg(4, isg))*gesgt(1:2,abs(fntsg(4, isg)), ie); vnorm(  3, 2, 2) = gesgt(3,abs(fntsg(4, isg)), ie);

vnorm(1:2, 1, 3) = sign(1, fntsg(5, isg))*gesgt(1:2,abs(fntsg(5, isg)), ie); vnorm(  3, 1, 3) = gesgt(3,abs(fntsg(5, isg)), ie);
vnorm(1:2, 2, 3) = sign(1, fntsg(6, isg))*gesgt(1:2,abs(fntsg(6, isg)), ie); vnorm(  3, 2, 3) = gesgt(3,abs(fntsg(6, isg)), ie);

vnorm(1:2, 1, 4) = sign(1, fntsg(7, isg))*gesgt(1:2,abs(fntsg(7, isg)), ie); vnorm(  3, 1, 4) = gesgt(3,abs(fntsg(7, isg)), ie);
vnorm(1:2, 2, 4) = sign(1, fntsg(8, isg))*gesgt(1:2,abs(fntsg(8, isg)), ie); vnorm(  3, 2, 4) = gesgt(3,abs(fntsg(8, isg)), ie);
!
do ifsg =1,2
vnorm(3, ifsg, 1) =  vnorm(  3, ifsg, 1)*wfgstm(ifsg, 1, isg);
vnorm(3, ifsg, 2) =  vnorm(  3, ifsg, 2)*wfgstm(ifsg, 2, isg);
vnorm(3, ifsg, 3) =  vnorm(  3, ifsg, 3)*wfgstm(ifsg, 3, isg);
vnorm(3, ifsg, 4) =  vnorm(  3, ifsg, 4)*wfgstm(ifsg, 4, isg);
enddo
!
!...loop over every subgrid vertex...
!
do ivsg = 1, 4
!...Skip the void center node
if(iptsg(ivsg, isg).eq.0) cycle
!
ulnpn(1:ndegr) = ulnpn(1:ndegr)  +&
ustar(1, ipt(iptsg(ivsg, isg)))*vnorm(1, 1, ivsg)*vnorm(  3, 1, ivsg)*bt(1:ndegr,  iptsg(ivsg, isg)) +&
ustar(2, ipt(iptsg(ivsg, isg)))*vnorm(2, 1, ivsg)*vnorm(  3, 1, ivsg)*bt(1:ndegr,  iptsg(ivsg, isg)) +&
ustar(1, ipt(iptsg(ivsg, isg)))*vnorm(1, 2, ivsg)*vnorm(  3, 2, ivsg)*bt(1:ndegr,  iptsg(ivsg, isg)) +&
ustar(2, ipt(iptsg(ivsg, isg)))*vnorm(2, 2, ivsg)*vnorm(  3, 2, ivsg)*bt(1:ndegr,  iptsg(ivsg, isg))
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstart(1, 1, ivsg, isg, ie)*bt(1:ndegr, iptsg(ivsg, isg))*wfgst(1, ivsg, isg) +&
fstart(1, 2, ivsg, isg, ie)*bt(1:ndegr, iptsg(ivsg, isg))*wfgst(2, ivsg, isg)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstart(2, 1, ivsg, isg, ie)*bt(1:ndegr, iptsg(ivsg, isg))*wfgst(1, ivsg, isg) +&
fstart(2, 2, ivsg, isg, ie)*bt(1:ndegr, iptsg(ivsg, isg))*wfgst(2, ivsg, isg)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipt(iptsg(ivsg, isg)))*fstart(1, 1, ivsg, isg,  ie)*&
bt(1:ndegr, iptsg(ivsg, isg))*wfgst(1, ivsg, isg) +&
ustar(2, ipt(iptsg(ivsg, isg)))*fstart(2, 1, ivsg, isg,  ie)*&
bt(1:ndegr, iptsg(ivsg, isg))*wfgst(1, ivsg, isg) +&
ustar(1, ipt(iptsg(ivsg, isg)))*fstart(1, 2, ivsg, isg, ie)*&
bt(1:ndegr, iptsg(ivsg, isg))*wfgst(2, ivsg, isg)  +&
ustar(2, ipt(iptsg(ivsg, isg)))*fstart(2, 2, ivsg, isg,  ie)*&
bt(1:ndegr, iptsg(ivsg, isg))*wfgst(2, ivsg, isg)

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
!if(ie==1)  print*,'rhs iface',ielem, ie, plnpn(1, 1:ndegr),fstarq(1, 1, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(1, ivsg, isg)
650 enddo
!
end subroutine rhsifacedg_lagsubgt2
!
!...Face integral (mass center) for hybrid quad...
!
subroutine rhsifacedg_lagsubgt3(iptri, unkno, ustar,fstart, gesgt, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:3, 1:4, 1:ntria),  intent(in)::fstart !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac), intent(in)::gesgt
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem,isg,ivsg,ifsg
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(6, 4)::fntsg
integer,dimension(3, 4)::iptsg
real*8,dimension(2,3,4)::wfgst,wfgstm
real*8,dimension(1:3,1:2,1:4)::vnorm
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::xvt(nvtri+1), yvt(nvtri+1),bt(1:ndegr,1:nvtri)

!...local real number
real*8::eps,c00,c05,c10,c20,c13,c16
real*8::dr,ds,rc,sc
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
data c16   / 0.1666666666666666d0 /
!
!...Zero out plnpn, ulnpn
!
iptsg(1, 1) = 1; iptsg(2, 1) = 4; iptsg(3, 1) = 6;
iptsg(1, 2) = 4; iptsg(2, 2) = 5; iptsg(3, 2) = 6;
iptsg(1, 3) = 4; iptsg(2, 3) = 2; iptsg(3, 3) = 5;
iptsg(1, 4) = 6; iptsg(2, 4) = 5; iptsg(3, 4) = 3;

!
!
fntsg(1, 1) =  6;  fntsg(2, 1) = 1;
fntsg(3, 1) = 10;  fntsg(4, 1) =-7;
fntsg(5, 1) = -13; fntsg(6, 1) =12;

fntsg(1, 2) = 7;  fntsg(2, 2) = 14;
fntsg(3, 2) = 8;  fntsg(4, 2) = 15;
fntsg(5, 2) = 9;  fntsg(6, 2) = 13;

fntsg(1, 3) =-14; fntsg(2, 3) =10;
fntsg(3, 3) = 2; fntsg(4, 3) = 3;
fntsg(5, 3) =11; fntsg(6, 3) =-8;

fntsg(1, 4) =12; fntsg(2, 4) =-9;
fntsg(3, 4) =-15; fntsg(4, 4) =11;
fntsg(5, 4) = 4; fntsg(6, 4) =5;

!...For momentum and energy eqaution
wfgst= 1.d0;

!...Conitnuity equation...
wfgstm(1, 1, 1) = 4.d0/6.d0; wfgstm(2, 1, 1) = 4.d0/6.d0;
wfgstm(1, 2, 1) = 4.d0/3.d0; wfgstm(2, 2, 1) = 8.d0/3.d0;
wfgstm(1, 3, 1) = 8.d0/3.d0; wfgstm(2, 3, 1) = 4.d0/3.d0;

wfgstm(1, 1, 2) = 8.d0/3.d0; wfgstm(2, 1, 2) = 8.d0/3.d0;
wfgstm(1, 2, 2) = 8.d0/3.d0; wfgstm(2, 2, 2) = 8.d0/3.d0;
wfgstm(1, 3, 2) = 8.d0/3.d0; wfgstm(2, 3, 2) = 8.d0/3.d0;

wfgstm(1, 1, 3) = 8.d0/3.d0; wfgstm(2, 1, 3) = 4.d0/3.d0;
wfgstm(1, 2, 3) = 4.d0/6.d0; wfgstm(2, 2, 3) = 4.d0/6.d0;
wfgstm(1, 3, 3) = 4.d0/3.d0; wfgstm(2, 3, 3) = 8.d0/3.d0;

wfgstm(1, 1, 4) = 4.d0/3.d0; wfgstm(2, 1, 4) = 8.d0/3.d0;
wfgstm(1, 2, 4) = 8.d0/3.d0; wfgstm(2, 2, 4) = 4.d0/3.d0;
wfgstm(1, 3, 4) = 4.d0/6.d0; wfgstm(2, 3, 4) = 4.d0/6.d0;
!
wfgstm = 0.25d0*wfgstm
!
!...Quads...
!
do 650 ie = 1,ntria !...(1)ie = 1,nelem
!
ielem = ie
!
!...The vertex constituting one cell...
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
xvt(4) = 0.5d0; yvt(4) = 0.d0
xvt(5) = 0.5d0; yvt(5) = 0.5d0
xvt(6) = 0.d0;  yvt(6) = 0.5d0

!...The pseudo-point 7
xvt(7) = 1.d0/3.d0; yvt(7) = 1.d0/3.d0
!
!
do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bt(4, iv) = 0.5d0*bt(2, iv)*bt(2, iv) - geoel(19, ielem)
bt(5, iv) = 0.5d0*bt(3, iv)*bt(3, iv) - geoel(20, ielem)
bt(6, iv) =       bt(2, iv)*bt(3, iv) - geoel(21, ielem)
endif
enddo
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
do isg = 1, 4
!
!...Normal vector for every face...
!
vnorm(1:2, 1, 1) = sign(1, fntsg(1, isg))*gesgt(1:2,abs(fntsg(1, isg)), ie); vnorm(  3, 1, 1) = gesgt(3,abs(fntsg(1, isg)), ie);
vnorm(1:2, 2, 1) = sign(1, fntsg(2, isg))*gesgt(1:2,abs(fntsg(2, isg)), ie); vnorm(  3, 2, 1) = gesgt(3,abs(fntsg(2, isg)), ie);

vnorm(1:2, 1, 2) = sign(1, fntsg(3, isg))*gesgt(1:2,abs(fntsg(3, isg)), ie); vnorm(  3, 1, 2) = gesgt(3,abs(fntsg(3, isg)), ie);
vnorm(1:2, 2, 2) = sign(1, fntsg(4, isg))*gesgt(1:2,abs(fntsg(4, isg)), ie); vnorm(  3, 2, 2) = gesgt(3,abs(fntsg(4, isg)), ie);

vnorm(1:2, 1, 3) = sign(1, fntsg(5, isg))*gesgt(1:2,abs(fntsg(5, isg)), ie); vnorm(  3, 1, 3) = gesgt(3,abs(fntsg(5, isg)), ie);
vnorm(1:2, 2, 3) = sign(1, fntsg(6, isg))*gesgt(1:2,abs(fntsg(6, isg)), ie); vnorm(  3, 2, 3) = gesgt(3,abs(fntsg(6, isg)), ie);

!
do ifsg =1,2
vnorm(3, ifsg, 1) =  vnorm(  3, ifsg, 1)*wfgstm(ifsg, 1, isg);
vnorm(3, ifsg, 2) =  vnorm(  3, ifsg, 2)*wfgstm(ifsg, 2, isg);
vnorm(3, ifsg, 3) =  vnorm(  3, ifsg, 3)*wfgstm(ifsg, 3, isg);
enddo
!
!...loop over every subgrid vertex...
!
do ivsg = 1, 3
!
ulnpn(1:ndegr) = ulnpn(1:ndegr)  +&
ustar(1, ipt(iptsg(ivsg, isg)))*vnorm(1, 1, ivsg)*vnorm(  3, 1, ivsg)*bt(1:ndegr,  iptsg(ivsg, isg)) +&
ustar(2, ipt(iptsg(ivsg, isg)))*vnorm(2, 1, ivsg)*vnorm(  3, 1, ivsg)*bt(1:ndegr,  iptsg(ivsg, isg)) +&
ustar(1, ipt(iptsg(ivsg, isg)))*vnorm(1, 2, ivsg)*vnorm(  3, 2, ivsg)*bt(1:ndegr,  iptsg(ivsg, isg)) +&
ustar(2, ipt(iptsg(ivsg, isg)))*vnorm(2, 2, ivsg)*vnorm(  3, 2, ivsg)*bt(1:ndegr,  iptsg(ivsg, isg))
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstart(1, 1, ivsg, isg, ie)*bt(1:ndegr, iptsg(ivsg, isg))*wfgst(1, ivsg, isg) +&
fstart(1, 2, ivsg, isg, ie)*bt(1:ndegr, iptsg(ivsg, isg))*wfgst(2, ivsg, isg)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstart(2, 1, ivsg, isg, ie)*bt(1:ndegr, iptsg(ivsg, isg))*wfgst(1, ivsg, isg) +&
fstart(2, 2, ivsg, isg, ie)*bt(1:ndegr, iptsg(ivsg, isg))*wfgst(2, ivsg, isg)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipt(iptsg(ivsg, isg)))*fstart(1, 1, ivsg, isg,  ie)*&
bt(1:ndegr, iptsg(ivsg, isg))*wfgst(1, ivsg, isg) +&
ustar(2, ipt(iptsg(ivsg, isg)))*fstart(2, 1, ivsg, isg,  ie)*&
bt(1:ndegr, iptsg(ivsg, isg))*wfgst(1, ivsg, isg) +&
ustar(1, ipt(iptsg(ivsg, isg)))*fstart(1, 2, ivsg, isg, ie)*&
bt(1:ndegr, iptsg(ivsg, isg))*wfgst(2, ivsg, isg)  +&
ustar(2, ipt(iptsg(ivsg, isg)))*fstart(2, 2, ivsg, isg,  ie)*&
bt(1:ndegr, iptsg(ivsg, isg))*wfgst(2, ivsg, isg)

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
!if(ie==1)  print*,'rhs iface',ielem, ie, plnpn(1, 1:ndegr),fstarq(1, 1, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(1, ivsg, isg)
650 enddo
!
end subroutine rhsifacedg_lagsubgt3

!
!...Face integral (mass center) for hybrid quad...
!
subroutine rhsifacedg_lagsubgq2(ipqua, unkno, ustar,fstarq, gesgq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:4, 1:4, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem,isg,ivsg,ifsg
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:2, 1:nvqua) :: ipfq
integer,dimension(8, 4)::fnqsg
integer,dimension(4, 4)::ipqsg
real*8,dimension(2,4,4)::wfgsq,wfgsqm
real*8,dimension(1:3,1:2,1:4)::vnorm
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvqua)::lpnpq
real*8::xvq(nvqua), yvq(nvqua),bq(1:ndegr,1:nvqua)

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
!...Zero out plnpn, ulnpn
!
ipfq(1, 1) = 4; ipfq(2, 1) = 2
ipfq(1, 2) = 1; ipfq(2, 2) = 3
ipfq(1, 3) = 2; ipfq(2, 3) = 4
ipfq(1, 4) = 3; ipfq(2, 4) = 1
!
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
!
do iv =1 ,nvqua
!
bq(1, iv) = 1.d0
!
 if(npoly.ge.1)then
  bq(2, iv) = (xvq(iv)-rc)/dr
  bq(3, iv) = (yvq(iv)-sc)/ds

!DGP2
  if(npoly.eq.2)then
   bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
   bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
   bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
  endif
 endif
enddo
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
do isg = 1, 4
!
!...Normal vector for every face...
!
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
end subroutine rhsifacedg_lagsubgq2
!
!...Face surrounding point...
!
subroutine getfasup(intfac, fasup)
use constant
implicit none
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:4,1:npoin),            intent(out)::fasup
!
!...Local integer
!
integer::ifa
integer::ipf(2)
integer, allocatable::indpt(:)
!
allocate(indpt(npoin))
!
indpt = 0
fasup = 0
!
do ifa = 1, nafac
!
ipf(1:2) = intfac(3:4, ifa)
!
indpt(ipf(1)) = indpt(ipf(1)) + 1
indpt(ipf(2)) = indpt(ipf(2)) + 1
!
fasup(indpt(ipf(1)), ipf(1)) = ifa
fasup(indpt(ipf(2)), ipf(2)) = ifa
enddo
!
deallocate(indpt)
end subroutine getfasup
!
!....domain integral for hybrid curv quad cells
!
subroutine rhsdomndg_lagmc_qsubg(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
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
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds, bv
real*8:: unknod(1:nq)
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
real*8, dimension(1: ndimn):: vgnul
real*8,dimension(1:ndegr,1:4)::unksgq
real*8,dimension(1:4, 1:4)::geoq_sub
real*8::xgp(2)
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
!...Quads
!
!...Loop over elements
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua, ie)
!
!...Points consitituting one element...
!
if(ncurv==0)then
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
endif
!
call getrhosubg(xpq,xpqi,unksgq, geoq_sub)
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
!...The derivatives of basis function...
!...Here dbdx means dbd(xsi), dbdy means dbd(eta)
!
dbdr(1)= 0.d0
dbdr(2)= 1.d0/dr
dbdr(3)= 0.d0

dbds(1)= 0.d0
dbds(2)= 0.d0
dbds(3)= 1.0/ds
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
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi
!
!...Jacobian transformation matrix
!
jacbf(1, 1) = dxdr; jacbf(1, 2) = dxds
jacbf(2, 1) = dydr; jacbf(2, 2) = dyds
!
!...Cofactor matrix of Jacobian transformation matrix
!
jacbg(1, 1) = dyds; jacbg(1, 2) =-dydr
jacbg(2, 1) =-dxds; jacbg(2, 2) = dxdr
!
!...Calculate G dot dbdx or dbdy
!
do ideg = 1, ndegr
gdshp(1, ideg) = jacbg(1, 1)*dbdr(ideg) + jacbg(1, 2)*dbds(ideg)
gdshp(2, ideg) = jacbg(2, 1)*dbdr(ideg) + jacbg(2, 2)*dbds(ideg)
enddo
!
!...Gauss points...
!
xg = r
yg = s
!
!...Basis function for solutions...
!
b(1) = 1.d0
b(2) = (xg-rc)/dr
b(3) = (yg-sc)/ds
!
!...Solution at the Gauss points...
!
unknod = 0.d0
!
do ideg =1,mdegr
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo
!
!...Primitive variables...
!
if(ndens.eq.3)then
!
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
bv(1) = 1.d0
bv(2) = (xg-rcv)/dr
bv(3) = (yg-scv)/ds
!
unknod(1) =0.d0
!
do ideg = 1,mdegr
unknod(1) = unknod(1) + unkno(ideg,1,ielem)*bv(ideg)
enddo
!
rhoma = 1.d0/unknod(1)
rhoad  = unknod(1)
!
xgp(:) = 0.d0
!
do ishp = 1, npqua
xgp(1) = xgp(1) + shpq(ishp)*xpq(1,ishp)
xgp(2) = xgp(2) + shpq(ishp)*xpq(2,ishp)
enddo
!
call  getrhog(xgp, xpq, rhoad, unksgq)
!
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = max(eps,(gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2)))
!
if(nlimi.eq.6)then
!
! if(ie.ge.990.and.ie.le.1000) print*,'domn unk',ie,unkno(1, 1, ie)
!
if(ndens.eq.3)then
!
rhomc = 1.d0/unkno(1, 1, ielem)
!
endif
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
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
!
uadv = unkno(1,2,ielem)  + dudr*b(2) + duds*b(3)
vadv = unkno(1,3,ielem)  + dvdr*b(2) + dvds*b(3)
!
pres = pctr + aflim(4, ielem)*(pres- pctr)
!
endif
!
fluxd(1,1) = gdshp(1, 1)*uadv + gdshp(2, 1)*vadv
fluxd(2,1) = gdshp(1, 2)*uadv + gdshp(2, 2)*vadv
fluxd(3,1) = gdshp(1, 3)*uadv + gdshp(2, 3)*vadv
!
fluxd(1,2) = gdshp(1, 1)*(-pres)
fluxd(2,2) = gdshp(1, 2)*(-pres)
fluxd(3,2) = gdshp(1, 3)*(-pres)
!
fluxd(1,3) = gdshp(2, 1)*(-pres)
fluxd(2,3) = gdshp(2, 2)*(-pres)
fluxd(3,3) = gdshp(2, 3)*(-pres)
!
fluxd(1,4) = (gdshp(1, 1)*uadv + gdshp(2, 1)*vadv)*(-pres)
fluxd(2,4) = (gdshp(1, 2)*uadv + gdshp(2, 2)*vadv)*(-pres)
fluxd(3,4) = (gdshp(1, 3)*uadv + gdshp(2, 3)*vadv)*(-pres)
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) - fluxd(ideg,1:nq)*djaco
enddo
!
!if(ie==1772) print*,'rhs iface idegr',ie,ig,fluxd(1:3,4),pres,djaco
!
enddo !...(2)ig = 1,ngausd
!
650 enddo

end subroutine rhsdomndg_lagmc_qsubg
!
!...Get density at the gauss point...
!
subroutine getrhog(xgp, xpq, rhog, unksgq)
use constant
real*8, dimension(1:2),                  intent(in)::xgp
real*8,dimension(1:ndimn, 1:npqua),      intent(in) :: xpq
real*8,dimension(1:ndegr,1:4),          intent(in)::unksgq
real*8,           intent(out)::rhog

integer,dimension(4, 4)::ipqsg
integer,dimension(2, 4)::itqsg

real*8,dimension(1:ndimn, 1:4) :: xqs
real*8,dimension(1:ndimn, 1:3) :: xt
!
integer::ivs, itri
real*8:: dx13, dy13, dx24, dy24
real*8:: dx12, dy12
real*8:: areaq, areat
!
!...Geometry information for subgrid quad...
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
itqsg(1, 1) = 1; itqsg(2, 1) = 2;
itqsg(1, 2) = 2; itqsg(2, 2) = 3;
itqsg(1, 3) = 3; itqsg(2, 3) = 4;
itqsg(1, 4) = 4; itqsg(2, 4) = 1;
!
do ivs =1, 4
!
xqs(1, 1:4) = xpq(1, ipqsg(1:4, ivs))
xqs(2, 1:4) = xpq(2, ipqsg(1:4, ivs))
!
dx13 = xqs(1, 3)-xqs(1, 1)
dy13 = xqs(2, 3)-xqs(2, 1)

dx24 = xqs(1, 4)-xqs(1, 2)
dy24 = xqs(2, 4)-xqs(2, 2)

areaq = 0.5d0*abs(dx13*dy24 - dx24*dy13)
!
areat = 0.d0
!
xt(1:2, 1) = xgp(1:2)
!
do itri = 1, 4
!
xt(1:2, 2) = xpq(1:2, ipqsg(itqsg(1, itri), ivs))
xt(1:2, 3) = xpq(1:2, ipqsg(itqsg(2, itri), ivs))
!
dx12 = xt(1, 2) - xt(1, 1)
dy12 = xt(2, 2) - xt(2, 1)

dx13 = xt(1, 3) - xt(1, 1)
dy13 = xt(2, 3) - xt(2, 1)
!
areat = areat + 0.5d0*abs(dx12*dy13 - dy12*dx13)
enddo
!
if(abs(areat-areaq)/areaq.le.1.d-12)then
rhog = unksgq(1, ivs)
exit
endif

enddo
end subroutine getrhog
!
!...subroutine: Calculate the F^* N dsfor all faces for hybrid grids...
!
subroutine getfnds_lag_simpsubg3(gflag,gesgq, gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(inout)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad),      intent(inout)::gelagq
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(inout)::gesgq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ifg,ifq,ihf,nhf
integer::iv,ishp
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:ndimn, 1:nvqua)::xpq
real*8,dimension(1:ndimn, 1:nvfac, 1:4)::xpfq

real*8::vnorm(1:2)
real*8,dimension(1:nvtri):: xv, yv
real*8::weigh(ngausf), posi(1,ngausf)
real*8::shp(nvfac),dshpr(nvfac)
real*8::dxdr,dxds,dydr,dyds
!...local real number
real*8::dwav1,dwav2,larea
real*8::anx, any, farea(3, 8)
real*8::r,s,djaco,wi
real*8::dr, ds, rc, sc
real*8::c16, c10
!
data c16   /0.1666666666666666d0 /
data c10   /1.0d0 /
!
call rutope(1, ngausf, posi, weigh)
!
do 200 ie=1, nquad !...(1)ifa=1,nquad
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...coordinates
!
if(ncurv==0)then
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif
!-xxx-face 1
xpfq(1:2, 1, 1) = xpq(1:2, 1)
xpfq(1:2, 2, 1) = xpq(1:2, 2)
xpfq(1:2, 3, 1) = xpq(1:2, 5)
!-xxx-face 2
xpfq(1:2, 1, 2) = xpq(1:2, 2)
xpfq(1:2, 2, 2) = xpq(1:2, 3)
xpfq(1:2, 3, 2) = xpq(1:2, 6)
!-xxx-face 3
xpfq(1:2, 1, 3) = xpq(1:2, 3)
xpfq(1:2, 2, 3) = xpq(1:2, 4)
xpfq(1:2, 3, 3) = xpq(1:2, 7)
!-xxx-face 4
xpfq(1:2, 1, 4) = xpq(1:2, 4)
xpfq(1:2, 2, 4) = xpq(1:2, 1)
xpfq(1:2, 3, 4) = xpq(1:2, 8)

farea = 0.d0
!
nhf = 0
!
do ifq = 1, 4!...(2)ig = 1,ngausd
do ihf = 1, 2

nhf = ihf + (ifq-1)*2
!
anx = 0.d0
any = 0.d0
!
do ifg = 1, ngausf
!
if(ihf.eq.1)then
r  = (posi(1,ifg)-1.d0)/2.d0
else
r  = (posi(1,ifg)+1.d0)/2.d0
endif
!
wi = weigh(ifg)
!
shp(1) =  -0.5d0*(1.d0-r)*r
shp(2) =   0.5d0*(1.d0+r)*r
shp(3) =         (1.d0+r)*(1.d0-r)
!
dshpr(1) = -0.5d0 + r
dshpr(2) =  0.5d0 + r
dshpr(3) = -2.d0*r
!
!...Jacobian determinant...
!
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, nptfa
dxdr = dxdr + dshpr(ishp)*xpfq(1, ishp, ifq)
dydr = dydr + dshpr(ishp)*xpfq(2, ishp, ifq)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
!
dwav1 = dydr/djaco
dwav2 =-dxdr/djaco
!
anx = anx + dwav1*0.5d0*djaco*wi
any = any + dwav2*0.5d0*djaco*wi
!
farea(3, nhf) = farea(3, nhf) + 0.5d0*djaco*wi !...farea: face area...
!
enddo
!
farea(1, nhf) = anx/farea(3, nhf)
farea(2, nhf) = any/farea(3, nhf)!sqrt(anx**2 + any**2)
!
enddo
!
enddo !...ifq from 1...4
!
gesgq(1, 1:8, ie) = farea(1, 1:8)
gesgq(2, 1:8, ie) = farea(2, 1:8)
gesgq(3, 1:8, ie) = farea(3, 1:8)
!
!gesgq(1:3, 9:12 , ie) = 0.d0
!
200 enddo  !...(1)ifa=1,nquad
!
!print*,'Inside getfnds_lag',gesgq(1:3,1,1)
!
end subroutine getfnds_lag_simpsubg3
!
!...subroutine: Calculate the F^* N dsfor all faces for hybrid grids...
!
subroutine getfnds_lag_simpsonh3(gflag,gesgq,intfac,inpoel,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
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
real*8::posiq(2, 12)
real*8::shp(nptri),dspr(nptri),dsps(nptri)
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
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
posit(1, 1) = 0.d0; posit(2, 1) = 0.d0;
posit(1, 2) = 1.d0; posit(2, 2) = 0.d0;
posit(1, 3) = 1.d0; posit(2, 3) = 0.d0;
posit(1, 4) = 0.d0; posit(2, 4) = 1.d0;
posit(1, 5) = 0.d0; posit(2, 5) = 1.d0;
posit(1, 6) = 0.d0; posit(2, 6) = 0.d0;
posit(1, 7) = .5d0; posit(2, 7) = 0.d0;
posit(1, 8) = .5d0; posit(2, 8) = .5d0;
posit(1, 9) = 0.d0; posit(2, 9) = .5d0;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;
!
do 100 ie=1, ntria!...(1)ifa=1,nafac
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipt(1:nvtri) = inpoel(1:nvtri,ie)
!
!...coordinates
!
if(ncurv==0)then
xpt(1, 1:3) = coord(1, iptri(1:3,ie))
xpt(2, 1:3) = coord(2, iptri(1:3,ie))
!
xpt(1:2,4) = 0.5d0*(xpt(1:2,1) + xpt(1:2,2))
xpt(1:2,5) = 0.5d0*(xpt(1:2,2) + xpt(1:2,3))
xpt(1:2,6) = 0.5d0*(xpt(1:2,1) + xpt(1:2,3))
elseif(ncurv==1)then
xpt(1, 1:nptri) = coord(1,iptri(1:nptri, ie))
xpt(2, 1:nptri) = coord(2,iptri(1:nptri, ie))
endif
!
do ig = 1, 9!...(2)ig = 1,ngausd
!
r  = posit(1,ig)
s  = posit(2,ig)
!
!...Shape functions and their derivatives...
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
do ishp = 1, nptri
dxdr = dxdr + dspr(ishp)*xpt(1,ishp)
dxds = dxds + dsps(ishp)*xpt(1,ishp)

dydr = dydr + dspr(ishp)*xpt(2,ishp)
dyds = dyds + dsps(ishp)*xpt(2,ishp)
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
enddo
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa
!
!...Identify the local No. of one internal face for left cell...
!
if(ig.eq.1.or.ig.eq.2.or.ig.eq.7)then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 1.d0;
!
idfal = ig;
!
elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.8)then
!
vnorm(1) = sqrt(2.d0)*0.5d0; vnorm(2) = sqrt(2.d0)*0.5d0;
larea    = sqrt(2.d0)
!
idfal = ig;
elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.9)then
!
vnorm(1) = -1.d0;           vnorm(2) = 0.d0;
larea    = 1.d0
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
gelag(1, idfal, ie) = anx*larea/farea
gelag(2, idfal, ie) = any*larea/farea
gelag(3, idfal, ie) = farea
!
enddo !...ig from 1...9
!
100 enddo  !...(1)ifa=1,nafac
!
! print*,'vnotmfn',gelag(1:3, 1:6, 170)
!
!
do 200 ie=1, nquad !...(1)ifa=1,nquad
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...coordinates
!
if(ncurv==0)then
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif
!
do ig = 1, 12!...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa
!
!...Identify the local No. of one internal face for left cell...
!
if(ig.eq.1.or.ig.eq.2.or.ig.eq.9)then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 2.d0;
!
idfal = ig;
!
elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.10)then
!
vnorm(1) = 1.d0; vnorm(2) = 0.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.11)then
!
vnorm(1) = 0.d0;           vnorm(2) = 1.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.7.or.ig.eq.8.or.ig.eq.12)then
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
gesgq(1:2, 1, ie) = gelagq(1:2, 1, ie);gesgq(3, 1, ie) = gelagq(3, 1, ie);
gesgq(1:2, 2, ie) = gelagq(1:2, 2, ie);gesgq(3, 2, ie) = gelagq(3, 2, ie);
gesgq(1:2, 3, ie) = gelagq(1:2, 3, ie);gesgq(3, 3, ie) = gelagq(3, 3, ie);
gesgq(1:2, 4, ie) = gelagq(1:2, 4, ie);gesgq(3, 4, ie) = gelagq(3, 4, ie);
gesgq(1:2, 5, ie) = gelagq(1:2, 5, ie);gesgq(3, 5, ie) = gelagq(3, 5, ie);
gesgq(1:2, 6, ie) = gelagq(1:2, 6, ie);gesgq(3, 6, ie) = gelagq(3, 6, ie);
gesgq(1:2, 7, ie) = gelagq(1:2, 7, ie);gesgq(3, 7, ie) = gelagq(3, 7, ie);
gesgq(1:2, 8, ie) = gelagq(1:2, 8, ie);gesgq(3, 8, ie) = gelagq(3, 8, ie);
!
gesgq(1:2, 13, ie) = gelagq(1:2, 9, ie);gesgq(3, 13, ie) = gelagq(3, 9, ie);
gesgq(1:2, 14, ie) = gelagq(1:2,10, ie);gesgq(3, 14, ie) = gelagq(3,10, ie);
gesgq(1:2, 15, ie) = gelagq(1:2,11, ie);gesgq(3, 15, ie) = gelagq(3,11, ie);
gesgq(1:2, 16, ie) = gelagq(1:2,12, ie);gesgq(3, 16, ie) = gelagq(3,12, ie);
!
!gesgq(1:2, 13, ie) = gelagq(1:2, 1, ie);gesgq(3, 13, ie) = 0.5d0*gelagq(3,1, ie);
!gesgq(1:2, 14, ie) = gelagq(1:2, 3, ie);gesgq(3, 14, ie) = 0.5d0*gelagq(3,3, ie);
!gesgq(1:2, 15, ie) = gelagq(1:2, 5, ie);gesgq(3, 15, ie) = 0.5d0*gelagq(3,5, ie);
!gesgq(1:2, 16, ie) = gelagq(1:2, 7, ie);gesgq(3, 16, ie) = 0.5d0*gelagq(3,7, ie);
!
!gesgq(:, 9:12, ie) = 0.d0
!
200 enddo  !...(1)ifa=1,nquad
!
!do ig = 1,16
!print*,'Inside getfnds_lag',ig,gesgq(1:3,ig,1435)
!enddo
!
end subroutine getfnds_lag_simpsonh3
!
!...subroutine: Calculate the F^* N dsfor all faces for hybrid grids...
!
subroutine getfnds_lagsms_simpsonhybrid(gflag,gesgq,gesgt,intfac,inpoel,iptri,ipqua,coord)
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
real*8::posiq(2, 12)
real*8::shp(nptri),dspr(nptri),dsps(nptri)
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
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
posit(1, 1) = 0.d0; posit(2, 1) = 0.d0;
posit(1, 2) = 1.d0; posit(2, 2) = 0.d0;
posit(1, 3) = 1.d0; posit(2, 3) = 0.d0;
posit(1, 4) = 0.d0; posit(2, 4) = 1.d0;
posit(1, 5) = 0.d0; posit(2, 5) = 1.d0;
posit(1, 6) = 0.d0; posit(2, 6) = 0.d0;
posit(1, 7) = .5d0; posit(2, 7) = 0.d0;
posit(1, 8) = .5d0; posit(2, 8) = .5d0;
posit(1, 9) = 0.d0; posit(2, 9) = .5d0;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;
!
do 100 ie=1, ntria!...(1)ifa=1,nafac
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipt(1:nvtri) = inpoel(1:nvtri,ie)
!
!...coordinates
!
if(ncurv==0)then
xpt(1, 1:3) = coord(1, iptri(1:3,ie))
xpt(2, 1:3) = coord(2, iptri(1:3,ie))
!
xpt(1:2,4) = 0.5d0*(xpt(1:2,1) + xpt(1:2,2))
xpt(1:2,5) = 0.5d0*(xpt(1:2,2) + xpt(1:2,3))
xpt(1:2,6) = 0.5d0*(xpt(1:2,1) + xpt(1:2,3))
elseif(ncurv==1)then
xpt(1, 1:nptri) = coord(1,iptri(1:nptri, ie))
xpt(2, 1:nptri) = coord(2,iptri(1:nptri, ie))
endif
!
do ig = 1, 9!...(2)ig = 1,ngausd
!
r  = posit(1,ig)
s  = posit(2,ig)
!
!...Shape functions and their derivatives...
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
do ishp = 1, nptri
dxdr = dxdr + dspr(ishp)*xpt(1,ishp)
dxds = dxds + dsps(ishp)*xpt(1,ishp)

dydr = dydr + dspr(ishp)*xpt(2,ishp)
dyds = dyds + dsps(ishp)*xpt(2,ishp)
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
enddo
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa
!
!...Identify the local No. of one internal face for left cell...
!
if(ig.eq.1.or.ig.eq.2.or.ig.eq.7)then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 1.d0;
!
idfal = ig;
!
elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.8)then
!
vnorm(1) = sqrt(2.d0)*0.5d0; vnorm(2) = sqrt(2.d0)*0.5d0;
larea    = sqrt(2.d0)
!
idfal = ig;
elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.9)then
!
vnorm(1) = -1.d0;           vnorm(2) = 0.d0;
larea    = 1.d0
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
gelag(1, idfal, ie) = anx*larea/farea
gelag(2, idfal, ie) = any*larea/farea
gelag(3, idfal, ie) = farea
enddo !...ig from 1...9
!
gesgt(1:2, 1, ie) = gelag(1:2, 1, ie);gesgt(3, 1, ie) = gelag(3, 1, ie);
gesgt(1:2, 2, ie) = gelag(1:2, 2, ie);gesgt(3, 2, ie) = gelag(3, 2, ie);
gesgt(1:2, 3, ie) = gelag(1:2, 3, ie);gesgt(3, 3, ie) = gelag(3, 3, ie);
gesgt(1:2, 4, ie) = gelag(1:2, 4, ie);gesgt(3, 4, ie) = gelag(3, 4, ie);
gesgt(1:2, 5, ie) = gelag(1:2, 5, ie);gesgt(3, 5, ie) = gelag(3, 5, ie);
gesgt(1:2, 6, ie) = gelag(1:2, 6, ie);gesgt(3, 6, ie) = gelag(3, 6, ie);

!
gesgt(1:2, 10, ie) = gelag(1:2, 7, ie);gesgt(3, 10, ie) = gelag(3, 7, ie);
gesgt(1:2, 11, ie) = gelag(1:2, 8, ie);gesgt(3, 11, ie) = gelag(3, 8, ie);
gesgt(1:2, 12, ie) = gelag(1:2, 9, ie);gesgt(3, 12, ie) = gelag(3, 9, ie);
!
100 enddo  !...(1)ifa=1,nafac
!
! print*,'vnotmfn',gelag(1:3, 1:6, 170)
!
!
do 200 ie=1, nquad !...(1)ifa=1,nquad
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...coordinates
!
if(ncurv==0)then
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif
!
do ig = 1, 12!...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
!if(ie==170) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa
!
!...Identify the local No. of one internal face for left cell...
!
if(ig.eq.1.or.ig.eq.2.or.ig.eq.9)then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 2.d0;
!
idfal = ig;
!
elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.10)then
!
vnorm(1) = 1.d0; vnorm(2) = 0.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.11)then
!
vnorm(1) = 0.d0;           vnorm(2) = 1.d0;
larea    = 2.d0
!
idfal = ig;
elseif(ig.eq.7.or.ig.eq.8.or.ig.eq.12)then
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
gesgq(1:2, 1, ie) = gelagq(1:2, 1, ie);gesgq(3, 1, ie) = gelagq(3, 1, ie);
gesgq(1:2, 2, ie) = gelagq(1:2, 2, ie);gesgq(3, 2, ie) = gelagq(3, 2, ie);
gesgq(1:2, 3, ie) = gelagq(1:2, 3, ie);gesgq(3, 3, ie) = gelagq(3, 3, ie);
gesgq(1:2, 4, ie) = gelagq(1:2, 4, ie);gesgq(3, 4, ie) = gelagq(3, 4, ie);
gesgq(1:2, 5, ie) = gelagq(1:2, 5, ie);gesgq(3, 5, ie) = gelagq(3, 5, ie);
gesgq(1:2, 6, ie) = gelagq(1:2, 6, ie);gesgq(3, 6, ie) = gelagq(3, 6, ie);
gesgq(1:2, 7, ie) = gelagq(1:2, 7, ie);gesgq(3, 7, ie) = gelagq(3, 7, ie);
gesgq(1:2, 8, ie) = gelagq(1:2, 8, ie);gesgq(3, 8, ie) = gelagq(3, 8, ie);
!
gesgq(1:2, 13, ie) = gelagq(1:2, 9, ie);gesgq(3, 13, ie) = gelagq(3, 9, ie);
gesgq(1:2, 14, ie) = gelagq(1:2,10, ie);gesgq(3, 14, ie) = gelagq(3,10, ie);
gesgq(1:2, 15, ie) = gelagq(1:2,11, ie);gesgq(3, 15, ie) = gelagq(3,11, ie);
gesgq(1:2, 16, ie) = gelagq(1:2,12, ie);gesgq(3, 16, ie) = gelagq(3,12, ie);
!
!gesgq(1:2, 13, ie) = gelagq(1:2, 1, ie);gesgq(3, 13, ie) = 0.5d0*gelagq(3,1, ie);
!gesgq(1:2, 14, ie) = gelagq(1:2, 3, ie);gesgq(3, 14, ie) = 0.5d0*gelagq(3,3, ie);
!gesgq(1:2, 15, ie) = gelagq(1:2, 5, ie);gesgq(3, 15, ie) = 0.5d0*gelagq(3,5, ie);
!gesgq(1:2, 16, ie) = gelagq(1:2, 7, ie);gesgq(3, 16, ie) = 0.5d0*gelagq(3,7, ie);
!
!gesgq(:, 9:12, ie) = 0.d0
!
200 enddo  !...(1)ifa=1,nquad
!
!do ig = 1,16
!print*,'Inside getfnds_lag',ig,gesgq(1:3,ig,1435)
!enddo
!
end subroutine getfnds_lagsms_simpsonhybrid

!
!...subroutine: Riemann input for hybrid curved tria using sub-cell scheme....
!
subroutine getriem_triasubg2(iptri, geoel, gesgt, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac),intent(in)::gesgt
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),         intent(in)::afvec
real*8, dimension(1:2, 1:2, 1:npoin),        intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin),         intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin),         intent(inout)::snsigm
real*8, dimension(1:2, 1:2, 1:2, 1:4, 1:3, 1:ntria), intent(out)::munaclt
real*8, dimension(1:ndimn, 1:2,  1:4, 1:3, 1:ntria), intent(out)::munault
real*8, dimension(1:ndimn, 1:2,  1:4, 1:3, 1:ntria), intent(out)::snsigmlt

!...Local integer
integer::ie, ideg, ielem, ifa, iv, isg, ivsg, ifsg

!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvfac) :: ipf
integer,dimension(8, 3)::fntsg
integer,dimension(4, 3)::iptsg
!...local real array
real*8,dimension(1:ndegr, 1:nvtri)::bt
real*8,dimension(1:ndegr, 1:4)::btv
real*8,dimension(1:nq,1:nvtri)::unknvt
real*8,dimension(1:nq,1:4)::unsgt
real*8::aujmp(1:3, 1:4)
real*8::vnorm(1:3, 1:2, 1:6)
real*8::sigma(1:2, 1:2, 1:4)
real*8,dimension(1:2, 1:4)::murie
real*8,dimension(1:nvtri):: xvt, yvt
real*8,dimension(1:ndimn, 1:nvtri) :: xpt
real*8,dimension(1:ndimn, 1:nvtri) :: xpti
real*8,dimension(1:ndegr, 1:3)::unksgt
real*8,dimension(1:nq+1, 1:3)::utsgc
real*8,dimension(2, 4, 3)::wfgst

!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx,rhovsg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::rcsg,scsg
real*8::dxp,dyp,xc,yc
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon
!
data eps   / 1.0d-15/
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Specify some Gauss points...
!

!...Local vertex No. of gauss points in one unit
!...0 means the subcell vertex is void and inside the triangle
iptsg(1, 1) = 1; iptsg(2, 1) = 4; iptsg(3, 1) = 0; iptsg(4, 1) = 6
iptsg(1, 2) = 4; iptsg(2, 2) = 2; iptsg(3, 2) = 5; iptsg(4, 2) = 0
iptsg(1, 3) = 0; iptsg(2, 3) = 5; iptsg(3, 3) = 3; iptsg(4, 3) = 6
!
fntsg(1, 1) =  6; fntsg(2, 1) = 1;
fntsg(3, 1) = 10; fntsg(4, 1) =-7;
fntsg(5, 1) = -7; fntsg(6, 1) = 9;
fntsg(7, 1) = 9;  fntsg(8, 1) = 12;

fntsg(1, 2) =  7; fntsg(2, 2) =10;
fntsg(3, 2) = 2;  fntsg(4, 2) =  3;
fntsg(5, 2) = 11; fntsg(6, 2) =-8;
fntsg(7, 2) = -8; fntsg(8, 2) = 7;

fntsg(1, 3) =-9; fntsg(2, 3) = 8;
fntsg(3, 3) = 8; fntsg(4, 3) =11;
fntsg(5, 3) = 4; fntsg(6, 3) = 5;
fntsg(7, 3) =12; fntsg(8, 3) =-9;

!...4/6 for internal face...This is used for our work
wfgst(1, 1, 1) = 4.d0/6.d0; wfgst(2, 1, 1) = 4.d0/6.d0;
wfgst(1, 2, 1) = 4.d0/3.d0; wfgst(2, 2, 1) = 8.d0/3.d0;
wfgst(1, 3, 1) = 1.d0;      wfgst(2, 3, 1) = 1.d0;
wfgst(1, 4, 1) = 8.d0/3.d0; wfgst(2, 4, 1) = 4.d0/3.d0;

wfgst(1, 1, 2) = 8.d0/3.d0; wfgst(2, 1, 2) = 4.d0/3.d0;
wfgst(1, 2, 2) = 4.d0/6.d0; wfgst(2, 2, 2) = 4.d0/6.d0;
wfgst(1, 3, 2) = 4.d0/3.d0; wfgst(2, 3, 2) = 8.d0/3.d0;
wfgst(1, 4, 2) = 1.d0;      wfgst(2, 4, 2) = 1.d0;

wfgst(1, 1, 3) = 1.d0;      wfgst(2, 1, 3) = 1.d0;
wfgst(1, 2, 3) = 8.d0/3.d0; wfgst(2, 2, 3) = 4.d0/3.d0;
wfgst(1, 3, 3) = 4.d0/6.d0; wfgst(2, 3, 3) = 4.d0/6.d0;
wfgst(1, 4, 3) = 4.d0/3.d0; wfgst(2, 4, 3) = 8.d0/3.d0
!
wfgst = 0.25d0*wfgst
!
!...Part II: Loop over every quad...
!
do 350 ie = 1,ntria !...(1)ie = 1,nquad
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie

!...shape functions
dr = 0.5d0
ds = 0.5d0

!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...Vertex coordinate
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
xvt(4) = .5d0; yvt(4) = 0.d0
xvt(5) = .5d0; yvt(5) = .5d0
xvt(6) = 0.d0; yvt(6) = .5d0

!...Basis function
do iv =1 ,nvtri
bt(1, iv) = 1.d0
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bt(4, iv) = 0.5d0*bt(2, iv)*bt(2, iv) - geoel(19, ielem)
bt(5, iv) = 0.5d0*bt(3, iv)*bt(3, iv) - geoel(20, ielem)
bt(6, iv) =       bt(2, iv)*bt(3, iv) - geoel(21, ielem)
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
unknvt = 0.d0
do iv   = 1,nvtri
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
enddo
enddo

!...Get density for sub-cells...
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
!call getrhosubgt_averg(xpt,xpti,unksgt,ielem)

!...Get density correction
!call getrhosubgt_daverg(rc,sc,geoel(19:21, ielem), xpt, unksgt, unkno(:,:,ielem), aflim(:,ielem), utsgc ,ielem)

call getrhosubgt_average(rc, sc, geoel(19:21, ielem), xpt,xpti,unksgt, unkno(:,:,ielem), aflim(:,ielem), utsgc, ielem)

!...Output for debugging
!if(ielem.eq.12)then
!print*,'Variabe',ielem,unksgt(1,:),1.d0/utsgc(1, :),1.d0/unkno(1,1,1)
!endif


!...II.1: Loop over sub-cells....
do isg = 1, 3

!...Normal vector for every quadrature point...
vnorm(1:2, 1, 1) = sign(1, fntsg(1, isg))*gesgt(1:2,abs(fntsg(1, isg)), ie); vnorm(  3, 1, 1) = gesgt(3,abs(fntsg(1, isg)), ie);
vnorm(1:2, 2, 1) = sign(1, fntsg(2, isg))*gesgt(1:2,abs(fntsg(2, isg)), ie); vnorm(  3, 2, 1) = gesgt(3,abs(fntsg(2, isg)), ie);

vnorm(1:2, 1, 2) = sign(1, fntsg(3, isg))*gesgt(1:2,abs(fntsg(3, isg)), ie); vnorm(  3, 1, 2) = gesgt(3,abs(fntsg(3, isg)), ie);
vnorm(1:2, 2, 2) = sign(1, fntsg(4, isg))*gesgt(1:2,abs(fntsg(4, isg)), ie); vnorm(  3, 2, 2) = gesgt(3,abs(fntsg(4, isg)), ie);

vnorm(1:2, 1, 3) = sign(1, fntsg(5, isg))*gesgt(1:2,abs(fntsg(5, isg)), ie); vnorm(  3, 1, 3) = gesgt(3,abs(fntsg(5, isg)), ie);
vnorm(1:2, 2, 3) = sign(1, fntsg(6, isg))*gesgt(1:2,abs(fntsg(6, isg)), ie); vnorm(  3, 2, 3) = gesgt(3,abs(fntsg(6, isg)), ie);

vnorm(1:2, 1, 4) = sign(1, fntsg(7, isg))*gesgt(1:2,abs(fntsg(7, isg)), ie); vnorm(  3, 1, 4) = gesgt(3,abs(fntsg(7, isg)), ie);
vnorm(1:2, 2, 4) = sign(1, fntsg(8, isg))*gesgt(1:2,abs(fntsg(8, isg)), ie); vnorm(  3, 2, 4) = gesgt(3,abs(fntsg(8, isg)), ie);

!...Get weighted area normal vector
do ifsg =1,2
vnorm(3, ifsg, 1) =  vnorm(  3, ifsg, 1)*wfgst(ifsg, 1, isg);
vnorm(3, ifsg, 2) =  vnorm(  3, ifsg, 2)*wfgst(ifsg, 2, isg);
vnorm(3, ifsg, 3) =  vnorm(  3, ifsg, 3)*wfgst(ifsg, 3, isg);
vnorm(3, ifsg, 4) =  vnorm(  3, ifsg, 4)*wfgst(ifsg, 4, isg);
enddo

!...Correct the density for Reimann input
do ivsg = 1,4
!...Skip the void center node
if(iptsg(ivsg, isg).eq.0) cycle
!
if(ndens.eq.1)then
rhovt = 1.d0/unknvt(1, iptsg(ivsg, isg))
rhovsg = rhovt +1.0d0*(unksgt(1,isg)-1.d0/utsgc(1, isg))
!rhovsg = unksgt(1,isg)

elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)

btv(1, ivsg) = 1.d0
btv(2, ivsg) = (xvt(iptsg(ivsg, isg))-rcv)/dr  !...bqv ....
btv(3, ivsg) = (yvt(iptsg(ivsg, isg))-scv)/ds
!
rhovt = 0.d0
do ideg = 1, mdegr
rhovt = rhovt + unkno(ideg,1,ielem)*btv(ideg, ivsg)
enddo
rhovsg = rhovt + (unksgt(1,isg)-rhoct)

!...Based on the sub-cell
!rcv = geoq_sub(1, isg); scv = geoq_sub(2, isg)
!bqv(1, ivsg) = 1.d0
!bqv(2, ivsg) = (xvq(ivsg)-rcv)/dr  !...bqv ....
!bqv(3, ivsg) = (yvq(ivsg)-scv)/ds

!rhovsg =0.d0
!do ideg = 1,ndegr
!rhovsg =  rhovsg + unksgq(ideg,isg)*bqv(ideg, ivsg)
!enddo

!rhovsg = unksgq(1,isg)! unkno(1,1,ielem)!
endif

uvtx = unknvt(2, iptsg(ivsg, isg))
vvtx = unknvt(3, iptsg(ivsg, isg))
evtx = unknvt(4, iptsg(ivsg, isg))
!
pvtx = max(eps, (gamlg-1.d0)*rhovsg*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unsgt(2, ivsg) = uvtx
unsgt(3 ,ivsg) = vvtx

!...Output for debugging
!if(ielem.eq.1)then
!print*,'smooth',ielem,unkno(:,1,ielem)
!print*,'Variabe',ielem,ivsg,rhovsg,unksgq(1,isg),1.d0/uqsgc(1, isg),pvtx,evtx,uvtx, vvtx
!endif

!...Limiter
!...For high-order method(>P1), if the cell includes shock, then
!...the high-order terms are chopped, and only limit the linear DG(P1) part.

if(nlimi.eq.6.and.geoel(10, ielem).gt.10.d0)then
if(ndens.eq.1)then
rhovsg = 1.d0/(rhomc + aflim(1, ielem)*(1.d0/rhovsg - rhomc))
elseif(ndens.eq.3)then
rhovsg = unksgt(1,isg) + aflim(1, ielem)*(rhovsg - unksgt(1,isg))
endif

!...Output for debugging
!print*,'Shock identification',ielem,unkno(:,:,ielem)
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bt(2, iptsg(ivsg, isg)) + duds*bt(3, iptsg(ivsg, isg))
vvtx = unkno(1,3,ielem)  + dvdr*bt(2, iptsg(ivsg, isg)) + dvds*bt(3, iptsg(ivsg, isg))
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)

!...Updtae velocity
unsgt(2, ivsg) = uvtx
unsgt(3 ,ivsg) = vvtx
endif

!...Get stress tensor at one vertex
sigma(1, 1, ivsg) = -pvtx
sigma(1, 2, ivsg) = 0.d0
sigma(2, 1, ivsg) = 0.d0
sigma(2, 2, ivsg) = -pvtx
!
!...Output for debugging
!if(ielem.eq.1)then
!print*,'smooth2',ielem,unkno(:,:,ielem)
!print*,'Variabe2',ielem,ivsg,rhovsg,pvtx,unsgq(2:3, ivsg),sigma(1, 1, ivsg)
!endif

!...Get the a_c (unit vector)
aujmp(1:2, ivsg) = vlave(1:2, ipt(iptsg(ivsg, isg))) - unsgt(2:3, ivsg)
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

!...Skip the void center node
if(iptsg(ivsg, isg).eq.0) cycle
!
dux= vlave(1, ipt(iptsg(ivsg, isg)))-unsgt(2, ivsg)
duy= vlave(2, ipt(iptsg(ivsg, isg)))-unsgt(3, ivsg)
deltu = sqrt(dux**2 + duy**2)
do ifa = 1, 2
deltu = 0.d0*abs(dux*vnorm(1, ifa, ivsg) + duy*vnorm(2, ifa, ivsg))
murie(ifa, ivsg) = rhoct*sdctr + rhoct*slpdu*deltu
!murie(ifa, ivsg) = uqsgc(1, isg)*sdctr + uqsgc(1, isg)*slpdu*deltu
!...The exact shock impedance of 2 shock model
!murie(ifa, ivsg) = rhoct*slpdu*deltu/2.d0+&
!     rhoct*sqrt((slpdu*deltu/2.d0)**2 + gamlg*pctr/rhoct)
enddo
!
enddo

!...Feed the input into Riemann solver
do ivsg  = 1, 4

!...Skip the void center node
if(iptsg(ivsg, isg).eq.0) cycle

!...Local vertex No. of gauss points in one unit
iv = iptsg(ivsg, isg)
!
do ifa = 1, 2 !...Every corner consists of 2 faces...

!...Call Riemann solver...
call getriecoef_matrixnew(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:3, ivsg), &
unsgt(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:2, ivsg), &
!unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipt(iv)) = munacn(1:2, 1, ipt(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipt(iv)) = munacn(1:2, 2, ipt(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipt(iv)) = munacu(1:2, ipt(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipt(iv)) = snsigm(1:2, ipt(iv)) + snsigm_rie(1:2)

!...Output for debugging
!if(ipq(iv).eq.1) print*,'p36 muacn(vv) post',ie,ifa,isg,ivsg,snsigm_rie(1:2),vnorm(1:3, ifa, ivsg),sigma(1:2, 1:2, ivsg)

!...Local variable...
munaclt(1:2, 1, ifa, ivsg, isg, ie) =  munacn_rie(1:2, 1)
munaclt(1:2, 2, ifa, ivsg, isg, ie) =  munacn_rie(1:2, 2)
!
munault(1:2,    ifa, ivsg, isg, ie) =  munacu_rie(1:2)
!
snsigmlt(1:2,   ifa, ivsg, isg, ie)=  snsigm_rie(1:2)
!
enddo
enddo
!
enddo
!
350 enddo  !...(1)ie = 1,nquad

end subroutine getriem_triasubg2
!
!...subroutine: Riemann input for hybrid curved quads using sub-cell scheme....
!
subroutine getriem_quadsubg2(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq,coord, coold, aflim, afvec,itime)
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
real*8, dimension(1:2, 1:2, 1:2, 1:4, 1:4, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:4, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:4, 1:nquad), intent(out)::snsigmlq
integer,      intent(in):: itime


!...Local integer
integer::ie, ideg, ielem, ifa, iv, isg, ivsg, ifsg

!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(8, 4)::fnqsg
integer,dimension(4, 4)::ipqsg
!...local real array
real*8,dimension(1:ndegr, 1:nvqua)::bq
real*8,dimension(1:ndegr, 1:4)::bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nq,1:4)::unsgq
real*8::aujmp(1:3, 1:4)
real*8::vnorm(1:3, 1:2, 1:6)
real*8::sigma(1:2, 1:2, 1:4)
real*8,dimension(1:2, 1:4)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndegr, 1:4)::unksgq
!real*8,dimension(1:ndegr, 1:4)::bqsg
real*8,dimension(1:nq+1, 1:4)::uqsgc
real*8,dimension(1:4, 1:4)::prsgq
real*8,dimension(1:4, 1:4)::bqvp
real*8,dimension(1:4)::prqz
real*8,dimension(1:7, 1:4)::geoq_sub
real*8,dimension(2, 4, 4)::wfgsq
real*8::bqsg(ndegr)

!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx,rhovsg
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

!...1/2 for internal face
wfgsq(1, 1, 1) = 4.d0/6.d0; wfgsq(2, 1, 1) = 4.d0/6.d0;
wfgsq(1, 2, 1) = 4.d0/3.d0; wfgsq(2, 2, 1) = 6.d0/3.d0;
wfgsq(1, 3, 1) = 2.d0;      wfgsq(2, 3, 1) = 2.d0;
wfgsq(1, 4, 1) = 6.d0/3.d0;      wfgsq(2, 4, 1) = 4.d0/3.d0;

wfgsq(1, 1, 2) = 6.d0/3.d0;      wfgsq(2, 1, 2) = 4.d0/3.d0;
wfgsq(1, 2, 2) = 4.d0/6.d0; wfgsq(2, 2, 2) = 4.d0/6.d0;
wfgsq(1, 3, 2) = 4.d0/3.d0; wfgsq(2, 3, 2) = 6.d0/3.d0;
wfgsq(1, 4, 2) = 1.d0;      wfgsq(2, 4, 2) = 1.d0;

wfgsq(1, 1, 3) = 1.d0;      wfgsq(2, 1, 3) = 1.d0;
wfgsq(1, 2, 3) = 6.d0/3.d0;      wfgsq(2, 2, 3) = 4.d0/3.d0;
wfgsq(1, 3, 3) = 4.d0/6.d0; wfgsq(2, 3, 3) = 4.d0/6.d0;
wfgsq(1, 4, 3) = 4.d0/3.d0; wfgsq(2, 4, 3) = 6.d0/3.d0

wfgsq(1, 1, 4) = 4.d0/3.d0; wfgsq(2, 1, 4) = 6.d0/3.d0;
wfgsq(1, 2, 4) = 1.d0;      wfgsq(2, 2, 4) = 1.d0;
wfgsq(1, 3, 4) = 6.d0/3.d0;      wfgsq(2, 3, 4) = 4.d0/3.d0;
wfgsq(1, 4, 4) = 4.d0/6.d0; wfgsq(2, 4, 4) = 4.d0/6.d0;

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
!
wfgsq = 0.25d0*wfgsq
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

!...Vertex coordinate
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

!...Basis function
do iv =1 ,nvqua
bq(1, iv) = 1.d0

 if(npoly.ge.1)then
   bq(2, iv) = (xvq(iv)-rc)/dr
   bq(3, iv) = (yvq(iv)-sc)/ds

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
!
call getrhosubg_averg(xpq,xpqi,unksgq, geoq_sub,ielem)

!...Get density correction
!!call getdens_quadsubg(rc, sc, geoel(19:21, ielem), unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc)
call getrhosubg_daverg(rc,sc, geoel(19:21, ielem), xpq, unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc, ielem)

!...Output for debugging
!if(ielem.eq.1)then
!print*,'Variabe',ielem,ivsg,unksgq(1,1),1.d0/uqsgc(1, 1)
!endif


!...II.1: Loop over sub-cells....
do isg = 1, 4

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
if(ndens.eq.1)then
 rhovt = 1.d0/unknvq(1, ipqsg(ivsg, isg))
! rhovsg = unksgq(1,isg)
!if(rkstg.eq.1)then
 rhovsg = rhovt +cdrho*(unksgq(1,isg)-1.d0/uqsgc(1, isg))
!else
! rhovsg = rhovt
!endif
!
! rcsg = geoq_sub(1, isg)
! scsg = geoq_sub(2, isg)
!Basis function
! bqsg(2) = (xvq(ivsg)-rcsg)/dr
! bqsg(3) = (yvq(ivsg)-scsg)/ds

!...DGP2
!if(npoly.eq.2)then
!bqsg(4) = 0.5d0* bqsg(2)* bqsg(2) - geoq_sub(5, isg)
!bqsg(5) = 0.5d0* bqsg(3)* bqsg(3) - geoq_sub(6, isg)
!bqsg(6) =        bqsg(2)* bqsg(3) - geoq_sub(7, isg)
!endif
!
! rhovsg = 0.d0
!do ideg = 1, mdegr
! rhovsg = rhovsg + unksgq(ideg,isg)*bqsg(ideg)
!enddo

elseif(ndens.eq.3)then
 rcv = geoel(5, ielem); scv = geoel(6, ielem)

 bqv(1, ivsg) = 1.d0

 if(npoly.ge.1)then
  bqv(2, ivsg) = (xvq(ipqsg(ivsg, isg))-rcv)/dr  !...bqv ....
  bqv(3, ivsg) = (yvq(ipqsg(ivsg, isg))-scv)/ds
 endif
!
 rhovt = 0.d0
 do ideg = 1, mdegr
 rhovt = rhovt + unkno(ideg,1,ielem)*bqv(ideg, ivsg)
 enddo
 rhovsg = rhovt + (unksgq(1,isg)-rhoct)

!...Based on the sub-cell
!rcv = geoq_sub(1, isg); scv = geoq_sub(2, isg)
!bqv(1, ivsg) = 1.d0
!bqv(2, ivsg) = (xvq(ivsg)-rcv)/dr  !...bqv ....
!bqv(3, ivsg) = (yvq(ivsg)-scv)/ds

!rhovsg =0.d0
!do ideg = 1,ndegr
!rhovsg =  rhovsg + unksgq(ideg,isg)*bqv(ideg, ivsg)
!enddo

!rhovsg = unksgq(1,isg)! unkno(1,1,ielem)!
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
!if(ielem.eq.1)then
!print*,'smooth',ielem,unkno(:,1,ielem)
!print*,'Variabe',ielem,ivsg,rhovsg,unksgq(1,isg),1.d0/uqsgc(1, isg),pvtx,evtx,uvtx, vvtx
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
!if(ielem.eq.1)then
!print*,'smooth2',ielem,unkno(:,:,ielem)
!print*,'Variabe2',ielem,ivsg,rhovsg,pvtx,unsgq(2:3, ivsg),sigma(1, 1, ivsg)
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
dux= vlave(1, ipq(ipqsg(ivsg, isg)))-unsgq(2, ivsg)
duy= vlave(2, ipq(ipqsg(ivsg, isg)))-unsgq(3, ivsg)
deltu = sqrt(dux**2 + duy**2)
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
!if(ipq(iv).eq.90.and.itime.ge.60) print*,'p36 muacn(vv) post',ie,ifa,isg,ivsg,snsigm_rie(1:2),&
!vnorm(1:3, ifa, ivsg),&
!sigma(1, 1, ivsg),&
!murie(ifa, ivsg),unsgq(2:3, ivsg),ipq(iv),rhoct,sdctr

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

end subroutine getriem_quadsubg2
!
!...Get the desnity at subgrid center...
!
subroutine getdens_quadsubg(rc, sc, bqho, unksgq, unkno, aflim, uqsgc)
use constant
implicit none
real*8, intent(in)::rc, sc
real*8,dimension(3),            intent(in)::bqho
real*8,dimension(1:ndegr, 1:4), intent(in)::unksgq
real*8,dimension(1:ndegr,1:nq),      intent(in)::unkno
real*8,dimension(1:nq+1),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:nq+1, 1:4),       intent(out)::uqsgc
!
!...Local arrays and real number
!
integer::isgc,isg,ideg
real*8,dimension(1:2,1:4)::rsgc
real*8,dimension(1:ndegr,1:4)::bqsgc
!
real*8::eps,c10
real*8::rhosgc
real*8::dr,ds
real*8::ectr,rhoct,uctr,vctr,rhomc
!
eps = 1.d-6
c10 = 1.d0
!
dr = 1.0d0
ds = 1.0d0
!
rsgc(1, 1) = -0.5d0; rsgc(2, 1) = -0.5d0;
rsgc(1, 2) =  0.5d0; rsgc(2, 2) = -0.5d0;
rsgc(1, 3) =  0.5d0; rsgc(2, 3) =  0.5d0;
rsgc(1, 4) = -0.5d0; rsgc(2, 4) =  0.5d0;
!
!...Cell averaged variables...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1)
elseif(ndens.eq.3)then
!-xxx-now configuration
rhomc = 1.d0/unkno(1, 1)
endif
!
uctr = unkno(1, 2)
vctr = unkno(1, 3)
ectr = unkno(1, 4)
!
rhoct  = 1.d0/rhomc !...Cell center density...
!
!...subgrid cell center basis function...
do isgc = 1, 4
bqsgc(1, isgc) = 1.d0;
bqsgc(2, isgc) = (rsgc(1, isgc)-rc)/dr;
bqsgc(3, isgc) = (rsgc(2, isgc)-sc)/ds;
!...DGP2
if(npoly.eq.2)then
bqsgc(4, isgc) = 0.5d0*bqsgc(2, isgc)*bqsgc(2, isgc) - bqho(1)
bqsgc(5, isgc) = 0.5d0*bqsgc(3, isgc)*bqsgc(3, isgc) - bqho(2)
bqsgc(6, isgc) =       bqsgc(2, isgc)*bqsgc(3, isgc) - bqho(3)
endif
enddo
!
!...The unknown at the subgrid cell center...
!
uqsgc = 0.d0
!
!...Density, velocity and total energy at subgrid center...
!
do isg = 1,4
do ideg = 1, mdegr
uqsgc(1:nq, isg) = uqsgc(1:nq, isg) + unkno(ideg,1:nq)*bqsgc(ideg, isg)
enddo
!
if(nlimi.eq.6)then
uqsgc(1, isg) = rhomc + aflim(1)*(uqsgc(1, isg) - rhomc)
endif
!
!print*,'pctr2'!,isg,rhosgc,usgc,vsgc,esgc,psgc,aflim(4),uqsgc(nq+1, isg)
enddo
end subroutine getdens_quadsubg
!
!...Get the averaged density at the sub cell using the whole cell distribution...
!
subroutine  getrhosubg_daverg(rc,sc, bqho, xpq, unksgq, unkno, aflim, uqsgc,ielem)
use constant
implicit none
!...Input
real*8, intent(in)::rc, sc
integer, intent(in)::ielem
real*8,dimension(3),            intent(in)::bqho
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:ndegr, 1:nq), intent(in)::unksgq
real*8,dimension(1:ndegr,1:nq),      intent(in)::unkno
real*8,dimension(1:nq+1),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:nq+1, 1:4),       intent(out)::uqsgc
!
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp,id, isg, ideg

!...Local real array
real*8,dimension(1:2, 1:npqua)::xpqsg
real*8,dimension(1:ndegr)::bqsg
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
!...Local real
real*8::r, s, xg,yg,dr,ds
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
!
xpqsg(1,1:9) =  xpq(1, 1:9)
xpqsg(2,1:9) =  xpq(2, 1:9)
!
rhomc = unkno(1,1)
!
!...Part 0: Get cell averaged r and s...
!
do isg = 1, 4

!...Initialze parameters...
masel = 0.d0
volel = 0.d0
!
do ig =1,ngausdq
!
if(isg.eq.1)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.2)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.3)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
elseif(isg.eq.4)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
endif
wi = weighq(ig)
!
xg = r
yg = s

!Basis function
bqsg(1) = 1.d0

if(npoly.ge.1)then
 bqsg(2) = (xg-rc)/dr
 bqsg(3) = (yg-sc)/ds

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

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsg(2,ishp)
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
end subroutine getrhosubg_daverg
!
!...Get the averaged density at the sub cell using the whole cell distribution with MEM...
!
subroutine  getrhosms_daverg(rc,sc, rci, sci, bqho, bqhoi, xpq, unksgq, unkno, unkgd,aflim, uqsgc,ielem)
use constant
implicit none
!...Input
real*8, intent(in)::rc, sc, rci, sci
integer, intent(in)::ielem
real*8,dimension(3),            intent(in)::bqho, bqhoi
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:ndegr, 1:nq), intent(in)::unksgq
real*8,dimension(1:ndegr,1:nq),         intent(in)::unkno
real*8,dimension(1:ndegr,1:4),  intent(in)::unkgd
real*8,dimension(1:nq+1),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:nq+1, 1:4),       intent(out)::uqsgc
!
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp,id, isg, ideg

!...Local real array
real*8,dimension(1:2, 1:npqua)::xpqsg
real*8,dimension(1:ndegr)::bqsg,bi
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8, dimension(1: ndimn, 1:ndimn)::jacbf
real*8:: weighq(ngausdq), posiq(2, ngausdq)
!...Local real
real*8::a11,a22,a12,a21
real*8::r, s, xg,yg,dr,ds
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
!
xpqsg(1,1:9) =  xpq(1, 1:9)
xpqsg(2,1:9) =  xpq(2, 1:9)
!
rhomc = unkno(1,1)
!
!...Part 0: Get cell averaged r and s...
!
do isg = 1, 4

!...Initialze parameters...
masel = 0.d0
volel = 0.d0
!
do ig =1,ngausdq
!
if(isg.eq.1)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.2)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)-1.d0)/2.d0
elseif(isg.eq.3)then
r = (posiq(1,ig)+1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
elseif(isg.eq.4)then
r = (posiq(1,ig)-1.d0)/2.d0
s = (posiq(2,ig)+1.d0)/2.d0
endif
wi = weighq(ig)
!
xg = r
yg = s

!Basis function
bqsg(1) = 1.d0
bqsg(2) = (xg-rc)/dr
bqsg(3) = (yg-sc)/ds

!...DGP2
if(npoly.eq.2)then
bqsg(4) = 0.5d0*bqsg(2)*bqsg(2) - bqho(1)
bqsg(5) = 0.5d0*bqsg(3)*bqsg(3) - bqho(2)
bqsg(6) =       bqsg(2)*bqsg(3) - bqho(3)
endif
!
rhomsg = 0.d0
do ideg = 1, mdegr
rhomsg = rhomsg + unkno(ideg,1)*bqsg(ideg)
enddo
!
if(nlimi.eq.6)then
rhomsg = rhomc + aflim(1)*(rhomsg - rhomc)
endif

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsg(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsg(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsg(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsg(2,ishp)
enddo

!...Jacobian transformation matrix
jacbf = 0.d0
!
bi(1) = 1.d0
bi(2) = (xg - rci)/dr
bi(3) = (yg - sci)/ds
!DGP2
if(npoly.eq.2)then
bi(4) = 0.5d0*bi(2)*bi(2)  - bqhoi(1)
bi(5) = 0.5d0*bi(3)*bi(3)  - bqhoi(2)
bi(6) =       bi(2)*bi(3)  - bqhoi(3)
endif
!
do ideg = 1, 1!ndegr
jacbf(1,1) = jacbf(1, 1) + unkgd(ideg, 1)*bi(ideg)
jacbf(1,2) = jacbf(1, 2) + unkgd(ideg, 2)*bi(ideg)
jacbf(2,1) = jacbf(2, 1) + unkgd(ideg, 3)*bi(ideg)
jacbf(2,2) = jacbf(2, 2) + unkgd(ideg, 4)*bi(ideg)
enddo
!
a11 = jacbf(1, 1)*dxdr + jacbf(1, 2)*dydr
a12 = jacbf(1, 1)*dxds + jacbf(1, 2)*dyds

a21 = jacbf(2, 1)*dxdr + jacbf(2, 2)*dydr
a22 = jacbf(2, 1)*dxds + jacbf(2, 2)*dyds
!
djaco = wi*(a11*a22 - a12*a21)*0.25d0
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
!if(ielem.eq.32)print*,'drho',isg,volel,masel
!
enddo
!
end subroutine getrhosms_daverg
!
!...Get the averaged density at the sub cell using the whole cell distribution...
!
subroutine  getrhosubg_daverg2(rc,sc, bqho, xpq, unksgq, unkno, aflim, uqsgc, ielem)
use constant
implicit none
!...Input
real*8, intent(in)::rc, sc
integer, intent(in)::ielem
real*8,dimension(3),            intent(in)::bqho
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:ndegr, 1:4), intent(in)::unksgq
real*8,dimension(1:ndegr,1:nq),      intent(in)::unkno
real*8,dimension(1:nq+1),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:nq+1, 1:4),       intent(out)::uqsgc

!...Local integer
integer :: ie, ig, ishp,  id, isg, ideg, isc
integer:: ipqsc(4, 4), ipqsg(4 ,4)

!...Local real array
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:ndegr)::bqsg
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8,dimension(1:2, 1:npqua)::xpqsg
real*8:: posqc(2, 12)
real*8:: xpqc(2, 12)
real*8:: xpqsc(2, npqua)
!
real*8,dimension(1:4)::shpqs
real*8:: weighq(ngausdq), posiq(2, ngausdq)
!...Local real
real*8::xsc,ysc
real*8::r, s, xg,yg,dr,ds
real*8::rg, sg
real*8::rp, rm ,sp, sm
real*8:: dxdr,dxds,dydr,dyds
real*8::rhomsg,rhomc
real*8::wi,djaco,volel, masel
real*8::c10
!
!-xxx-real contant
data c10 / 1.0d0 /

!...Get Guass point
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0

!...The vertex No. or the local quad inside one quad
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
posqc(1, 1)= -0.5d0; posqc(2, 1)= -1.d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)= -1.d0;
posqc(1, 3)=  1.0d0; posqc(2, 3)= -.5d0;
posqc(1, 4)=  1.0d0; posqc(2, 4)=  .5d0;
posqc(1, 5)=  0.5d0; posqc(2, 5)=  1.d0;
posqc(1, 6)= -0.5d0; posqc(2, 6)=  1.d0;
posqc(1, 7)= -1.0d0; posqc(2, 7)=  .5d0;
posqc(1, 8)= -1.0d0; posqc(2, 8)= -.5d0;
posqc(1, 9)=  0.0d0; posqc(2, 9)= -.5d0;
posqc(1,10)=  0.5d0; posqc(2,10)=  0.d0;
posqc(1,11)=  0.0d0; posqc(2,11)=  .5d0;
posqc(1,12)= -0.5d0; posqc(2,12)=  0.d0;
!
ipqsc(1, 1)= 1; ipqsc(2, 1)= 9; ipqsc(3, 1)= 12; ipqsc(4, 1)= 8;
ipqsc(1, 2)= 2; ipqsc(2, 2)= 3; ipqsc(3, 2)= 10; ipqsc(4, 2)= 9;
ipqsc(1, 3)=10; ipqsc(2, 3)= 4; ipqsc(3, 3)=  5; ipqsc(4, 3)=11;
ipqsc(1, 4)=12; ipqsc(2, 4)=11; ipqsc(3, 4)=  6; ipqsc(4, 4)= 7;
!
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

!...Now configuration...
xpqsg(1,1:9) =  xpq(1, 1:9)
xpqsg(2,1:9) =  xpq(2, 1:9)

!...Get the 12 high-order curved nodes
do isc = 1, 12
r = posqc(1,isc)
s = posqc(2,isc)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)

!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpqsg(1,ishp)
ysc = ysc + shpq(ishp)*xpqsg(2,ishp)
enddo
!
xpqc(1, isc) = xsc
xpqc(2, isc) = ysc
!
enddo

!...Mass average specific volume
rhomc = unkno(1,1)
!
!...Part 0: Get cell averaged r and s...
!
do isg = 1, 4

!...Initialze parameters...
masel = 0.d0
volel = 0.d0
!
!...Coordinates
xpqsc(1, 1:4) = xpqsg(1, ipqsg(1:4, isg))
xpqsc(2, 1:4) = xpqsg(2, ipqsg(1:4, isg))

xpqsc(1, 5:8) = xpqc(1, ipqsc(1:4, isg))
xpqsc(2, 5:8) = xpqc(2, ipqsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
!
rp = c10 + posiq(1,ig)
rm = c10 - posiq(1,ig)
sp = c10 + posiq(2,ig)
sm = c10 - posiq(2,ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpqs(1) = 0.25d0*rm*sm
shpqs(2) = 0.25d0*rp*sm
shpqs(3) = 0.25d0*rp*sp
shpqs(4) = 0.25d0*rm*sp
!
wi = weighq(ig)

!...Get the mapping gauss points in the standard triangle...
rg = 0.d0
sg = 0.d0

do ishp = 1, 4
rg = rg + shpqs(ishp)*xvq(ipqsg(ishp, isg))
sg = sg + shpqs(ishp)*yvq(ipqsg(ishp, isg))
enddo
!
!if(ielem.eq.1.and.isg.eq.1)print*,'ielem1',ig, r, s,posiq(1:2,ig)
!
xg = rg
yg = sg

!Basis function
bqsg(1) = 1.d0
bqsg(2) = (xg-rc)/dr
bqsg(3) = (yg-sc)/ds

!...DGP2
if(npoly.eq.2)then
bqsg(4) = 0.5d0*bqsg(2)*bqsg(2) - bqho(1)
bqsg(5) = 0.5d0*bqsg(3)*bqsg(3) - bqho(2)
bqsg(6) =       bqsg(2)*bqsg(3) - bqho(3)
endif
!
rhomsg = 0.d0
do ideg = 1, mdegr
rhomsg = rhomsg + unkno(ideg,1)*bqsg(ideg)
enddo
!
if(nlimi.eq.6)then
rhomsg = rhomc + aflim(1)*(rhomsg - rhomc)
endif

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)!*0.25d0
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
!if(ielem.eq.32)print*,'drho',isg,volel,masel
!
enddo
!
end subroutine getrhosubg_daverg2

!
!...Get the averaged density at the sub cell using the whole cell distribution...
!
subroutine  getrhosubgt_daverg(rc,sc, btho, xpt, unksgt, unkno, aflim, utsgc ,ielem)
use constant
implicit none
!...Input
real*8, intent(in)::rc, sc
integer, intent(in)::ielem
real*8,dimension(3),            intent(in)::btho
real*8,dimension(1:2, 1:nptri), intent(in)::xpt
real*8,dimension(1:ndegr, 1:3), intent(in)::unksgt
real*8,dimension(1:ndegr,1:nq),      intent(in)::unkno
real*8,dimension(1:nq+1),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:nq+1, 1:3),       intent(out)::utsgc

!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, isc
integer:: iptsc(4, 3), iptsg(4 ,3)

!...Local real array
real*8:: xvt(nvtri+1), yvt(nvtri+1)
real*8,dimension(1:ndegr)::btsg
real*8,dimension(1:nptri)::shpt
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8,dimension(1:2, 1:nptri+1)::xptsg
real*8:: postc(2, 9)
real*8:: xptc(2, 9)
real*8:: xpqsc(2, npqua)
!
real*8,dimension(1:4)::shpqs
real*8:: weighq(ngausdq), posiq(2, ngausdq)
!...Local real
real*8::xsc,ysc
real*8::r, s, xg,yg,dr,ds
real*8::rg, sg
real*8::rp, rm ,sp, sm
real*8:: dxdr,dxds,dydr,dyds
real*8::rhomsg,rhomc
real*8::wi,djaco,volel, masel
real*8::c10
!
!-xxx-real contant
data c10 / 1.0d0 /

!...Get Guass point
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 0.5d0
ds = 0.5d0

!...The vertex No. or the local quad inside one triangle
iptsg(1, 1) = 1; iptsg(2, 1) = 4; iptsg(3, 1) = 7; iptsg(4, 1) = 6
iptsg(1, 2) = 4; iptsg(2, 2) = 2; iptsg(3, 2) = 5; iptsg(4, 2) = 7
iptsg(1, 3) = 7; iptsg(2, 3) = 5; iptsg(3, 3) = 3; iptsg(4, 3) = 6
!
postc(1, 1)=  0.25d0; postc(2, 1)=   0.d0;
postc(1, 2)=  0.75d0; postc(2, 2)=   0.d0;
postc(1, 3)=  0.75d0; postc(2, 3)= 0.25d0;
postc(1, 4)=  0.25d0; postc(2, 4)= 0.75d0;
postc(1, 5)=  0.d0;   postc(2, 5)= 0.75d0;
postc(1, 6)=  0.d0;   postc(2, 6)= 0.25d0;
postc(1, 7)=  5.d0/12.d0;   postc(2, 7)= 1.d0/6.d0;
postc(1, 8)=  5.d0/12.d0;   postc(2, 8)= 5.d0/12.d0;
postc(1, 9)=  1.d0/6.d0;    postc(2, 9)= 5.d0/12.d0;
!
iptsc(1, 1)= 1; iptsc(2, 1)= 7; iptsc(3, 1)= 9; iptsc(4, 1)= 6;
iptsc(1, 2)= 2; iptsc(2, 2)= 3; iptsc(3, 2)= 8; iptsc(4, 2)= 7;
iptsc(1, 3)= 8; iptsc(2, 3)= 4; iptsc(3, 3)= 5; iptsc(4, 3)= 9;
!
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
xvt(4) = 0.5d0; yvt(4) = 0.d0
xvt(5) = 0.5d0; yvt(5) = 0.5d0
xvt(6) = 0.d0;  yvt(6) = 0.5d0

!...The pseudo-point 7
xvt(7) = 1.d0/3.d0; yvt(7) = 1.d0/3.d0

!...Now configuration...
xptsg(1,1:6) =  xpt(1, 1:6)
xptsg(2,1:6) =  xpt(2, 1:6)

!...Build the 7th node
xptsg(1, 7) = -1.d0/9.d0*(xpt(1, 1)+xpt(1, 2)+xpt(1, 3)) +&
               4.d0/9.d0*(xpt(1, 4)+xpt(1, 5)+xpt(1, 6))

xptsg(2, 7) = -1.d0/9.d0*(xpt(2, 1)+xpt(2, 2)+xpt(2, 3)) +&
               4.d0/9.d0*(xpt(2, 4)+xpt(2, 5)+xpt(2, 6))

!...Get the 12 high-order curved nodes
do isc = 1, 9
r = postc(1,isc)
s = postc(2,isc)

!...  shape function
shpt(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shpt(2) = -r*(c10-2.d0*r)
shpt(3) = -s*(c10-2.d0*s)
shpt(4) = 4.d0*r*(c10-r-s)
shpt(5) = 4.d0*r*s
shpt(6) = 4.d0*s*(c10-r-s)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, nptri
xsc = xsc + shpt(ishp)*xptsg(1,ishp)
ysc = ysc + shpt(ishp)*xptsg(2,ishp)
enddo
!
xptc(1, isc) = xsc
xptc(2, isc) = ysc
!
enddo

!...Mass average specific volume
rhomc = unkno(1,1)
!
!...Part 0: Get cell averaged r and s...
!
do isg = 1, 3

!...Initialze parameters...
masel = 0.d0
volel = 0.d0
!
!...Coordinates
xpqsc(1, 1:4) = xptsg(1, iptsg(1:4, isg))
xpqsc(2, 1:4) = xptsg(2, iptsg(1:4, isg))

xpqsc(1, 5:8) = xptc(1, iptsc(1:4, isg))
xpqsc(2, 5:8) = xptc(2, iptsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
!
rp = c10 + posiq(1,ig)
rm = c10 - posiq(1,ig)
sp = c10 + posiq(2,ig)
sm = c10 - posiq(2,ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpqs(1) = 0.25d0*rm*sm
shpqs(2) = 0.25d0*rp*sm
shpqs(3) = 0.25d0*rp*sp
shpqs(4) = 0.25d0*rm*sp
!
wi = weighq(ig)

!...Get the mapping gauss points in the standard triangle...
rg = 0.d0
sg = 0.d0

do ishp = 1, 4
rg = rg + shpqs(ishp)*xvt(iptsg(ishp, isg))
sg = sg + shpqs(ishp)*yvt(iptsg(ishp, isg))
enddo
!
xg = rg
yg = sg

!Basis function
btsg(1) = 1.d0
btsg(2) = (xg-rc)/dr
btsg(3) = (yg-sc)/ds

!...DGP2
if(npoly.eq.2)then
btsg(4) = 0.5d0*btsg(2)*btsg(2) - btho(1)
btsg(5) = 0.5d0*btsg(3)*btsg(3) - btho(2)
btsg(6) =       btsg(2)*btsg(3) - btho(3)
endif
!
rhomsg = 0.d0
do ideg = 1, mdegr
rhomsg = rhomsg + unkno(ideg,1)*btsg(ideg)
enddo
!
if(nlimi.eq.6)then
rhomsg = rhomc + aflim(1)*(rhomsg - rhomc)
endif

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)!*0.25d0
!
!...Parameters for now configuration...
!
masel = masel + 1.d0/rhomsg*djaco
volel = volel + djaco
!
enddo
!
!if(ielem.eq.12)print*,'ielem',volel,masel
!...The avearged specific volume at the sub-cell...
utsgc(1, isg) = volel/masel
!
enddo
!
end subroutine getrhosubgt_daverg
!
!...Get the averaged-density within one sub-cell for triangles...
!
subroutine  getrhosubgt_average(rc, sc,  btho, xpt,xpti,unksgt, unkno, aflim, utsgc, ielem)
use constant
implicit none
!...Input
real*8, intent(in)::rc, sc
integer, intent(in)::ielem
real*8,dimension(3),            intent(in)::btho
real*8,dimension(1:2, 1:nptri), intent(in)::xpt
real*8,dimension(1:2, 1:nptri), intent(in)::xpti
real*8,dimension(1:ndegr, 1:3), intent(out)::unksgt
real*8,dimension(1:ndegr,1:nq),      intent(in)::unkno
real*8,dimension(1:nq+1),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:nq+1, 1:3),       intent(out)::utsgc
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc

integer,dimension(4, 3)::iptsg, iptsc
!...Local real array
real*8:: xvt(nvtri+1), yvt(nvtri+1)
real*8,dimension(1:ndegr)::btsg
real*8,dimension(1:ndegr,1:3)::rhsel
real*8,dimension(1:nmatr, 1:3)::matin
real*8,dimension(1:2, 1:nptri+1)::xptsg,xptisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:nptri)::shpt
real*8,dimension(1:4)::shpqs
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: postc(2, 9)
real*8:: xptc(2, 9),xptic(2, 9)
real*8:: xpqsc(2, npqua), xpqisc(2, npqua)
!...Local real
real*8:: xsc, ysc, xsci, ysci
real*8::r, s, xg,yg,dr,ds
real*8::rg,sg,rm,rp,sm,sp
real*8:: dxdr,dxds,dydr,dyds
real*8::rhomc,rhomsg
real*8::wi,djaco
real*8:: rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::f0,f1,f2,f3
real*8::masel,volel,xgaus,ygaus
!
!-xxx-real contant
data c10 / 1.0d0 /
!
dr = 0.5d0
ds = 0.5d0
!
iptsg(1, 1) = 1; iptsg(2, 1) = 4; iptsg(3, 1) = 7; iptsg(4, 1) = 6
iptsg(1, 2) = 4; iptsg(2, 2) = 2; iptsg(3, 2) = 5; iptsg(4, 2) = 7
iptsg(1, 3) = 7; iptsg(2, 3) = 5; iptsg(3, 3) = 3; iptsg(4, 3) = 6
!
postc(1, 1)=  0.25d0; postc(2, 1)=   0.d0;
postc(1, 2)=  0.75d0; postc(2, 2)=   0.d0;
postc(1, 3)=  0.75d0; postc(2, 3)= 0.25d0;
postc(1, 4)=  0.25d0; postc(2, 4)= 0.75d0;
postc(1, 5)=  0.d0;   postc(2, 5)= 0.75d0;
postc(1, 6)=  0.d0;   postc(2, 6)= 0.25d0;
postc(1, 7)=  5.d0/12.d0;   postc(2, 7)= 1.d0/6.d0;
postc(1, 8)=  5.d0/12.d0;   postc(2, 8)= 5.d0/12.d0;
postc(1, 9)=  1.d0/6.d0;    postc(2, 9)= 5.d0/12.d0;
!
iptsc(1, 1)= 1; iptsc(2, 1)= 7; iptsc(3, 1)= 9; iptsc(4, 1)= 6;
iptsc(1, 2)= 2; iptsc(2, 2)= 3; iptsc(3, 2)= 8; iptsc(4, 2)= 7;
iptsc(1, 3)= 8; iptsc(2, 3)= 4; iptsc(3, 3)= 5; iptsc(4, 3)= 9;
!
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
xvt(4) = 0.5d0; yvt(4) = 0.d0
xvt(5) = 0.5d0; yvt(5) = 0.5d0
xvt(6) = 0.d0;  yvt(6) = 0.5d0
!
xvt(7) = 1.d0/3.d0;  yvt(7) = 1.d0/3.d0
!
call ruqope(2, ngausdq, posiq, weighq)

!...Now configuration...
xptsg(1,1:6) =  xpt(1, 1:6)
xptsg(2,1:6) =  xpt(2, 1:6)

!...Initial configuration...
xptisg(1,1:6) =  xpti(1, 1:6)
xptisg(2,1:6) =  xpti(2, 1:6)

!...Get the 12 high-order curved nodes
do isc = 1, 9
r = postc(1,isc)
s = postc(2,isc)

!...  shape function
shpt(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shpt(2) = -r*(c10-2.d0*r)
shpt(3) = -s*(c10-2.d0*s)
shpt(4) = 4.d0*r*(c10-r-s)
shpt(5) = 4.d0*r*s
shpt(6) = 4.d0*s*(c10-r-s)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, nptri
xsc = xsc + shpt(ishp)*xptsg(1,ishp)
ysc = ysc + shpt(ishp)*xptsg(2,ishp)
enddo
!
xsci = 0.d0
ysci = 0.d0
!
do ishp = 1, nptri
xsci = xsci + shpt(ishp)*xptisg(1,ishp)
ysci = ysci + shpt(ishp)*xptisg(2,ishp)
enddo
!
xptc(1, isc) = xsc
xptc(2, isc) = ysc
!
xptic(1, isc) = xsci
xptic(2, isc) = ysci
enddo
!
!...Part I: Construct the nodal coordinates
!
xptsg(1, 1:6) = xpt(1, 1:6)
xptsg(2, 1:6) = xpt(2, 1:6)
!
xptsg(1, 7) = -1.d0/9.d0*xpt(1, 1)-1.d0/9.d0*xpt(1, 2)-1.d0/9.d0*xpt(1, 3) +&
4.d0/9.d0*xpt(1, 4)+4.d0/9.d0*xpt(1, 5)+4.d0/9.d0*xpt(1, 6)

xptsg(2, 7) = -1.d0/9.d0*xpt(2, 1)-1.d0/9.d0*xpt(2, 2)-1.d0/9.d0*xpt(2, 3) +&
4.d0/9.d0*xpt(2, 4)+4.d0/9.d0*xpt(2, 5)+4.d0/9.d0*xpt(2, 6)
!
xptisg(1, 1:6) = xpti(1, 1:6)
xptisg(2, 1:6) = xpti(2, 1:6)
!
xptisg(1, 7) = -1.d0/9.d0*xpti(1, 1)-1.d0/9.d0*xpti(1, 2)-1.d0/9.d0*xpti(1, 3) +&
4.d0/9.d0*xpti(1, 4)+4.d0/9.d0*xpti(1, 5)+4.d0/9.d0*xpti(1, 6)

xptisg(2, 7) = -1.d0/9.d0*xpti(2, 1)-1.d0/9.d0*xpti(2, 2)-1.d0/9.d0*xpti(2, 3) +&
4.d0/9.d0*xpti(2, 4)+4.d0/9.d0*xpti(2, 5)+4.d0/9.d0*xpti(2, 6)

!
!...Part II: Get mass matrix for four subgrids...
!

!...Mass average specific volume
rhomc = unkno(1,1)
!
do isg = 1, 3
!...Coordinates
xpqsc(1, 1:4) = xptsg(1, iptsg(1:4, isg))
xpqsc(2, 1:4) = xptsg(2, iptsg(1:4, isg))

xpqsc(1, 5:8) = xptc(1, iptsc(1:4, isg))
xpqsc(2, 5:8) = xptc(2, iptsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initialze parameters...
f0 = 0.d0
masel = 0.d0
volel = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...Get the global reference coord inside one triangle
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s

!...  shape function & its derivatives w.r.t. reference coordinates
shpqs(1) = 0.25d0*rm*sm
shpqs(2) = 0.25d0*rp*sm
shpqs(3) = 0.25d0*rp*sp
shpqs(4) = 0.25d0*rm*sp

!...Get the mapping gauss points in the standard triangle...
rg = 0.d0
sg = 0.d0

do ishp = 1, 4
rg = rg + shpqs(ishp)*xvt(iptsg(ishp, isg))
sg = sg + shpqs(ishp)*yvt(iptsg(ishp, isg))
enddo
!
!Basis function
btsg(1) = 1.d0
btsg(2) = (rg-rc)/dr
btsg(3) = (sg-sc)/ds

!...DGP2
if(npoly.eq.2)then
btsg(4) = 0.5d0*btsg(2)*btsg(2) - btho(1)
btsg(5) = 0.5d0*btsg(3)*btsg(3) - btho(2)
btsg(6) =       btsg(2)*btsg(3) - btho(3)
endif
!
rhomsg = 0.d0
do ideg = 1, mdegr
rhomsg = rhomsg + unkno(ideg,1)*btsg(ideg)
enddo
!
if(nlimi.eq.6)then
rhomsg = rhomc + aflim(1)*(rhomsg - rhomc)
endif


!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco

!...Mass and volume
masel = masel + 1.d0/rhomsg*djaco
volel = volel + djaco
!
!if(ielem.eq.1) print*,'djaco',f0,wi,dxdr,dyds,dydr,dxds,xpqsc(1:2,:)
enddo

!
if(npoly==1)then
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then
matin(16,isg) = 1.d0/f0
endif

!...The avearged specific volume at the sub-cell...
utsgc(1, isg) = volel/masel

enddo
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, 3

!...Initial coordinates
xpqisc(1, 1:4) = xptisg(1, iptsg(1:4, isg))
xpqisc(2, 1:4) = xptisg(2, iptsg(1:4, isg))

xpqisc(1, 5:8) = xptic(1, iptsc(1:4, isg))
xpqisc(2, 5:8) = xptic(2, iptsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqisc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqisc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqisc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqisc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
b(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpqisc(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqisc(2,ishp)
enddo
!
call getrhoig_triacurv2(rhoi, xpti, xgaus, ygaus)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*b(1)*djaco
!if(isg.eq.4)print*,'isgideg',rhsel(6, isg),b(6),b(2),b(3), geoq_sub(7, isg)
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

!if(isg.eq.4.and.id.eq.1) print*,'isg4',m(1, :),rhsel(:,4),unint(1)
unksgt(1, isg)= unint(1)
enddo
!
end subroutine  getrhosubgt_average
!
!...subroutine: Riemann input for hybrid curved tria using second sub-cell scheme....
!
subroutine getriem_triasubg3(iptri, geoel, gesgt, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac),intent(in)::gesgt
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),         intent(in)::afvec
real*8, dimension(1:2, 1:2, 1:npoin),        intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin),         intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin),         intent(inout)::snsigm
real*8, dimension(1:2, 1:2, 1:2, 1:3, 1:4, 1:ntria), intent(out)::munaclt
real*8, dimension(1:ndimn, 1:2,  1:3, 1:4, 1:ntria), intent(out)::munault
real*8, dimension(1:ndimn, 1:2,  1:3, 1:4, 1:ntria), intent(out)::snsigmlt

!...Local integer
integer::ie, ideg, ielem, ifa, iv, isg, ivsg, ifsg

!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvfac) :: ipf
integer,dimension(6, 4)::fntsg
integer,dimension(3, 4)::iptsg
!...local real array
real*8,dimension(1:ndegr, 1:nvtri)::bt
real*8,dimension(1:ndegr, 1:4)::btv
real*8,dimension(1:nq,1:nvtri)::unknvt
real*8,dimension(1:nq,1:4)::unsgt
real*8::aujmp(1:3, 1:4)
real*8::vnorm(1:3, 1:2, 1:6)
real*8::sigma(1:2, 1:2, 1:3)
real*8,dimension(1:2, 1:4)::murie
real*8,dimension(1:nvtri):: xvt, yvt
real*8,dimension(1:ndimn, 1:nvtri) :: xpt
real*8,dimension(1:ndimn, 1:nvtri) :: xpti
real*8,dimension(1:ndegr, 1:4)::unksgt
real*8,dimension(1:nq+1, 1:4)::utsgc
real*8,dimension(2, 3, 4)::wfgst

!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx,rhovsg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::rcsg,scsg
real*8::dxp,dyp,xc,yc
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon
!
data eps   / 1.0d-15/
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Specify some Gauss points...
!

!...Local vertex No. of gauss points in one unit
!...0 means the subcell vertex is void and inside the triangle
iptsg(1, 1) = 1; iptsg(2, 1) = 4; iptsg(3, 1) = 6;
iptsg(1, 2) = 4; iptsg(2, 2) = 5; iptsg(3, 2) = 6;
iptsg(1, 3) = 4; iptsg(2, 3) = 2; iptsg(3, 3) = 5;
iptsg(1, 4) = 6; iptsg(2, 4) = 5; iptsg(3, 4) = 3;

!
fntsg(1, 1) =  6;  fntsg(2, 1) = 1;
fntsg(3, 1) = 10;  fntsg(4, 1) =-7;
fntsg(5, 1) = -13; fntsg(6, 1) =12;

fntsg(1, 2) = 7;  fntsg(2, 2) = 14;
fntsg(3, 2) = 8;  fntsg(4, 2) = 15;
fntsg(5, 2) = 9;  fntsg(6, 2) = 13;

fntsg(1, 3) =-14; fntsg(2, 3) =10;
fntsg(3, 3) = 2; fntsg(4, 3) = 3;
fntsg(5, 3) =11; fntsg(6, 3) =-8;

fntsg(1, 4) =12; fntsg(2, 4) =-9;
fntsg(3, 4) =-15; fntsg(4, 4) =11;
fntsg(5, 4) = 4; fntsg(6, 4) =5;

!...4/6 for internal face...This is used for our work
wfgst(1, 1, 1) = 4.d0/6.d0; wfgst(2, 1, 1) = 4.d0/6.d0;
wfgst(1, 2, 1) = 4.d0/3.d0; wfgst(2, 2, 1) = 8.d0/3.d0;
wfgst(1, 3, 1) = 8.d0/3.d0; wfgst(2, 3, 1) = 4.d0/3.d0;

wfgst(1, 1, 2) = 8.d0/3.d0; wfgst(2, 1, 2) = 8.d0/3.d0;
wfgst(1, 2, 2) = 8.d0/3.d0; wfgst(2, 2, 2) = 8.d0/3.d0;
wfgst(1, 3, 2) = 8.d0/3.d0; wfgst(2, 3, 2) = 8.d0/3.d0;

wfgst(1, 1, 3) = 8.d0/3.d0; wfgst(2, 1, 3) = 4.d0/3.d0;
wfgst(1, 2, 3) = 4.d0/6.d0; wfgst(2, 2, 3) = 4.d0/6.d0;
wfgst(1, 3, 3) = 4.d0/3.d0; wfgst(2, 3, 3) = 8.d0/3.d0;

wfgst(1, 1, 4) = 4.d0/3.d0; wfgst(2, 1, 4) = 8.d0/3.d0;
wfgst(1, 2, 4) = 8.d0/3.d0; wfgst(2, 2, 4) = 4.d0/3.d0;
wfgst(1, 3, 4) = 4.d0/6.d0; wfgst(2, 3, 4) = 4.d0/6.d0;
!
!wfgst(1, 1, 1) = 4.d0/6.d0; wfgst(2, 1, 1) = 4.d0/6.d0;
!wfgst(1, 2, 1) = 4.d0/3.d0; wfgst(2, 2, 1) = 0.d0/3.d0;
!wfgst(1, 3, 1) = 0.d0/3.d0; wfgst(2, 3, 1) = 4.d0/3.d0;

!wfgst(1, 1, 2) = 0.d0/3.d0; wfgst(2, 1, 2) = 0.d0/3.d0;
!wfgst(1, 2, 2) = 0.d0/3.d0; wfgst(2, 2, 2) = 0.d0/3.d0;
!wfgst(1, 3, 2) = 0.d0/3.d0; wfgst(2, 3, 2) = 0.d0/3.d0;

!wfgst(1, 1, 3) = 0.d0/3.d0; wfgst(2, 1, 3) = 4.d0/3.d0;
!wfgst(1, 2, 3) = 4.d0/6.d0; wfgst(2, 2, 3) = 4.d0/6.d0;
!wfgst(1, 3, 3) = 4.d0/3.d0; wfgst(2, 3, 3) = 0.d0/3.d0;

!wfgst(1, 1, 4) = 4.d0/3.d0; wfgst(2, 1, 4) = 0.d0/3.d0;
!wfgst(1, 2, 4) = 0.d0/3.d0; wfgst(2, 2, 4) = 4.d0/3.d0;
!wfgst(1, 3, 4) = 4.d0/6.d0; wfgst(2, 3, 4) = 4.d0/6.d0;

!
wfgst = 0.25d0*wfgst
!
!...Part II: Loop over every quad...
!
do 350 ie = 1,ntria !...(1)ie = 1,nquad
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie

!...shape functions
dr = 0.5d0
ds = 0.5d0

!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...Vertex coordinate
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
xvt(4) = .5d0; yvt(4) = 0.d0
xvt(5) = .5d0; yvt(5) = .5d0
xvt(6) = 0.d0; yvt(6) = .5d0

!...Basis function
do iv =1 ,nvtri
bt(1, iv) = 1.d0
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bt(4, iv) = 0.5d0*bt(2, iv)*bt(2, iv) - geoel(19, ielem)
bt(5, iv) = 0.5d0*bt(3, iv)*bt(3, iv) - geoel(20, ielem)
bt(6, iv) =       bt(2, iv)*bt(3, iv) - geoel(21, ielem)
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
unknvt = 0.d0
do iv   = 1,nvtri
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
enddo
enddo

!...Get density for sub-cells...
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))

!...Get density correction
call getrhosubgt_average2(rc, sc, geoel(19:21, ielem), xpt,xpti,unksgt, unkno(:,:,ielem), aflim(:,ielem), utsgc, ielem)

!...Output for debugging
!if(ielem.eq.12)then
!print*,'Variabe',ielem,unksgt(1,:),1.d0/utsgc(1, :),1.d0/unkno(1,1,1)
!endif


!...II.1: Loop over sub-cells....
do isg = 1, 4

!...Normal vector for every quadrature point...
vnorm(1:2, 1, 1) = sign(1, fntsg(1, isg))*gesgt(1:2,abs(fntsg(1, isg)), ie); vnorm(  3, 1, 1) = gesgt(3,abs(fntsg(1, isg)), ie);
vnorm(1:2, 2, 1) = sign(1, fntsg(2, isg))*gesgt(1:2,abs(fntsg(2, isg)), ie); vnorm(  3, 2, 1) = gesgt(3,abs(fntsg(2, isg)), ie);

vnorm(1:2, 1, 2) = sign(1, fntsg(3, isg))*gesgt(1:2,abs(fntsg(3, isg)), ie); vnorm(  3, 1, 2) = gesgt(3,abs(fntsg(3, isg)), ie);
vnorm(1:2, 2, 2) = sign(1, fntsg(4, isg))*gesgt(1:2,abs(fntsg(4, isg)), ie); vnorm(  3, 2, 2) = gesgt(3,abs(fntsg(4, isg)), ie);

vnorm(1:2, 1, 3) = sign(1, fntsg(5, isg))*gesgt(1:2,abs(fntsg(5, isg)), ie); vnorm(  3, 1, 3) = gesgt(3,abs(fntsg(5, isg)), ie);
vnorm(1:2, 2, 3) = sign(1, fntsg(6, isg))*gesgt(1:2,abs(fntsg(6, isg)), ie); vnorm(  3, 2, 3) = gesgt(3,abs(fntsg(6, isg)), ie);

!...Get weighted area normal vector
do ifsg =1,2
vnorm(3, ifsg, 1) =  vnorm(  3, ifsg, 1)*wfgst(ifsg, 1, isg);
vnorm(3, ifsg, 2) =  vnorm(  3, ifsg, 2)*wfgst(ifsg, 2, isg);
vnorm(3, ifsg, 3) =  vnorm(  3, ifsg, 3)*wfgst(ifsg, 3, isg);
enddo

!...Correct the density for Reimann input
do ivsg = 1,3
!
if(ndens.eq.1)then
rhovt = 1.d0/unknvt(1, iptsg(ivsg, isg))
rhovsg = rhovt +1.d0*(unksgt(1,isg)-1.d0/utsgc(1, isg))
!rhovsg = unksgt(1,isg)

elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)

btv(1, ivsg) = 1.d0
btv(2, ivsg) = (xvt(iptsg(ivsg, isg))-rcv)/dr  !...bqv ....
btv(3, ivsg) = (yvt(iptsg(ivsg, isg))-scv)/ds
!
rhovt = 0.d0
do ideg = 1, mdegr
rhovt = rhovt + unkno(ideg,1,ielem)*btv(ideg, ivsg)
enddo
rhovsg = rhovt + (unksgt(1,isg)-rhoct)

!...Based on the sub-cell
!rcv = geoq_sub(1, isg); scv = geoq_sub(2, isg)
!bqv(1, ivsg) = 1.d0
!bqv(2, ivsg) = (xvq(ivsg)-rcv)/dr  !...bqv ....
!bqv(3, ivsg) = (yvq(ivsg)-scv)/ds

!rhovsg =0.d0
!do ideg = 1,ndegr
!rhovsg =  rhovsg + unksgq(ideg,isg)*bqv(ideg, ivsg)
!enddo

!rhovsg = unksgq(1,isg)! unkno(1,1,ielem)!
endif

uvtx = unknvt(2, iptsg(ivsg, isg))
vvtx = unknvt(3, iptsg(ivsg, isg))
evtx = unknvt(4, iptsg(ivsg, isg))
!
pvtx = max(eps, (gamlg-1.d0)*rhovsg*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unsgt(2, ivsg) = uvtx
unsgt(3 ,ivsg) = vvtx

!...Output for debugging
!if(ielem.eq.1)then
!print*,'smooth',ielem,unkno(:,1,ielem)
!print*,'Variabe',ielem,ivsg,rhovsg,unksgq(1,isg),1.d0/uqsgc(1, isg),pvtx,evtx,uvtx, vvtx
!endif

!...Limiter
!...For high-order method(>P1), if the cell includes shock, then
!...the high-order terms are chopped, and only limit the linear DG(P1) part.

if(nlimi.eq.6.and.geoel(10, ielem).gt.10.d0)then
if(ndens.eq.1)then
rhovsg = 1.d0/(rhomc + aflim(1, ielem)*(1.d0/rhovsg - rhomc))
elseif(ndens.eq.3)then
rhovsg = unksgt(1,isg) + aflim(1, ielem)*(rhovsg - unksgt(1,isg))
endif

!...Output for debugging
!print*,'Shock identification',ielem,unkno(:,:,ielem)
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bt(2, iptsg(ivsg, isg)) + duds*bt(3, iptsg(ivsg, isg))
vvtx = unkno(1,3,ielem)  + dvdr*bt(2, iptsg(ivsg, isg)) + dvds*bt(3, iptsg(ivsg, isg))
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)

!...Updtae velocity
unsgt(2, ivsg) = uvtx
unsgt(3 ,ivsg) = vvtx
endif

!...Get stress tensor at one vertex
sigma(1, 1, ivsg) = -pvtx
sigma(1, 2, ivsg) = 0.d0
sigma(2, 1, ivsg) = 0.d0
sigma(2, 2, ivsg) = -pvtx
!
!...Output for debugging
!if(ielem.eq.1)then
!print*,'smooth2',ielem,unkno(:,:,ielem)
!print*,'Variabe2',ielem,ivsg,rhovsg,pvtx,unsgq(2:3, ivsg),sigma(1, 1, ivsg)
!endif

!...Get the a_c (unit vector)
aujmp(1:2, ivsg) = vlave(1:2, ipt(iptsg(ivsg, isg))) - unsgt(2:3, ivsg)
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
do ivsg   = 1, 3

!
dux= vlave(1, ipt(iptsg(ivsg, isg)))-unsgt(2, ivsg)
duy= vlave(2, ipt(iptsg(ivsg, isg)))-unsgt(3, ivsg)
deltu = sqrt(dux**2 + duy**2)
do ifa = 1, 2
deltu = 0.d0*abs(dux*vnorm(1, ifa, ivsg) + duy*vnorm(2, ifa, ivsg))
murie(ifa, ivsg) = rhoct*sdctr + rhoct*slpdu*deltu
!murie(ifa, ivsg) = uqsgc(1, isg)*sdctr + uqsgc(1, isg)*slpdu*deltu
!...The exact shock impedance of 2 shock model
!murie(ifa, ivsg) = rhoct*slpdu*deltu/2.d0+&
!     rhoct*sqrt((slpdu*deltu/2.d0)**2 + gamlg*pctr/rhoct)
enddo
!
enddo

!...Feed the input into Riemann solver
do ivsg  = 1, 3

!...Local vertex No. of gauss points in one unit
iv = iptsg(ivsg, isg)
!
do ifa = 1, 2 !...Every corner consists of 2 faces...

!...Call Riemann solver...
call getriecoef_matrixnew(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:3, ivsg), &
unsgt(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:2, ivsg), &
!unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipt(iv)) = munacn(1:2, 1, ipt(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipt(iv)) = munacn(1:2, 2, ipt(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipt(iv)) = munacu(1:2, ipt(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipt(iv)) = snsigm(1:2, ipt(iv)) + snsigm_rie(1:2)

!...Output for debugging
!if(ipq(iv).eq.1) print*,'p36 muacn(vv) post',ie,ifa,isg,ivsg,snsigm_rie(1:2),vnorm(1:3, ifa, ivsg),sigma(1:2, 1:2, ivsg)

!...Local variable...
munaclt(1:2, 1, ifa, ivsg, isg, ie) =  munacn_rie(1:2, 1)
munaclt(1:2, 2, ifa, ivsg, isg, ie) =  munacn_rie(1:2, 2)
!
munault(1:2,    ifa, ivsg, isg, ie) =  munacu_rie(1:2)
!
snsigmlt(1:2,   ifa, ivsg, isg, ie)=  snsigm_rie(1:2)
!
enddo
enddo
!
enddo
!
350 enddo  !...(1)ie = 1,nquad

end subroutine getriem_triasubg3

!
!...Get the averaged-density within one sub-cell for triangles for second sgs...
!
subroutine  getrhosubgt_average2(rc, sc,  btho, xpt,xpti,unksgt, unkno, aflim, utsgc, ielem)
use constant
implicit none
!...Input
real*8, intent(in)::rc, sc
integer, intent(in)::ielem
real*8,dimension(3),            intent(in)::btho
real*8,dimension(1:2, 1:nptri), intent(in)::xpt
real*8,dimension(1:2, 1:nptri), intent(in)::xpti
real*8,dimension(1:ndegr, 1:4), intent(out)::unksgt
real*8,dimension(1:ndegr,1:nq),      intent(in)::unkno
real*8,dimension(1:nq+1),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:nq+1, 1:4),       intent(out)::utsgc
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc

integer,dimension(3, 4)::iptsg, iptsc
!...Local real array
real*8:: xvt(nvtri+1), yvt(nvtri+1)
real*8,dimension(1:ndegr)::btsg
real*8,dimension(1:ndegr,1:4)::rhsel
real*8,dimension(1:nmatr, 1:4)::matin
real*8,dimension(1:2, 1:nptri+1)::xptsg,xptisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:nptri)::shpt
real*8,dimension(1:3)::shps
real*8,dimension(1:nptri)::shp, dspr, dsps
real*8:: weigh(ngausd), posi(2, ngausd)
real*8:: postc(2, 9)
real*8:: xptc(2, 9),xptic(2, 9)
real*8:: xptsc(2, nptri), xptisc(2, nptri)
!...Local real
real*8:: xsc, ysc, xsci, ysci
real*8::r, s, xg,yg,dr,ds
real*8::rg,sg,rm,rp,sm,sp
real*8:: dxdr,dxds,dydr,dyds
real*8::rhomc,rhomsg
real*8::wi,djaco
real*8:: rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::f0,f1,f2,f3
real*8::masel,volel,xgaus,ygaus
!
!-xxx-real contant
data c10 / 1.0d0 /
!
dr = 0.5d0
ds = 0.5d0
!
iptsg(1, 1) = 1; iptsg(2, 1) = 4; iptsg(3, 1) = 6;
iptsg(1, 2) = 4; iptsg(2, 2) = 5; iptsg(3, 2) = 6;
iptsg(1, 3) = 4; iptsg(2, 3) = 2; iptsg(3, 3) = 5;
iptsg(1, 4) = 6; iptsg(2, 4) = 5; iptsg(3, 4) = 3;
!
postc(1, 1)=  0.25d0; postc(2, 1)=   0.d0;
postc(1, 2)=  0.75d0; postc(2, 2)=   0.d0;
postc(1, 3)=  0.75d0; postc(2, 3)= 0.25d0;
postc(1, 4)=  0.25d0; postc(2, 4)= 0.75d0;
postc(1, 5)=  0.d0;   postc(2, 5)= 0.75d0;
postc(1, 6)=  0.d0;   postc(2, 6)= 0.25d0;
postc(1, 7)=  0.25d0; postc(2, 7)= 0.25d0;
postc(1, 8)=  0.5d0;  postc(2, 8)= 0.25d0;
postc(1, 9)=  0.25d0; postc(2, 9)= 0.5d0;
!
iptsc(1, 1)= 1; iptsc(2, 1)= 7; iptsc(3, 1)= 6;
iptsc(1, 2)= 8; iptsc(2, 2)= 9; iptsc(3, 2)= 7;
iptsc(1, 3)= 2; iptsc(2, 3)= 3; iptsc(3, 3)= 8;
iptsc(1, 4)= 9; iptsc(2, 4)= 4; iptsc(3, 4)= 5;
!
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
xvt(4) = 0.5d0; yvt(4) = 0.d0
xvt(5) = 0.5d0; yvt(5) = 0.5d0
xvt(6) = 0.d0;  yvt(6) = 0.5d0
!
call rutope(2, ngausd, posi, weigh)

!...Now configuration...
xptsg(1,1:6) =  xpt(1, 1:6)
xptsg(2,1:6) =  xpt(2, 1:6)

!...Initial configuration...
xptisg(1,1:6) =  xpti(1, 1:6)
xptisg(2,1:6) =  xpti(2, 1:6)

!...Get the 12 high-order curved nodes
do isc = 1, 9
r = postc(1,isc)
s = postc(2,isc)

!...  shape function
shpt(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shpt(2) = -r*(c10-2.d0*r)
shpt(3) = -s*(c10-2.d0*s)
shpt(4) = 4.d0*r*(c10-r-s)
shpt(5) = 4.d0*r*s
shpt(6) = 4.d0*s*(c10-r-s)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, nptri
xsc = xsc + shpt(ishp)*xptsg(1,ishp)
ysc = ysc + shpt(ishp)*xptsg(2,ishp)
enddo
!
xsci = 0.d0
ysci = 0.d0
!
do ishp = 1, nptri
xsci = xsci + shpt(ishp)*xptisg(1,ishp)
ysci = ysci + shpt(ishp)*xptisg(2,ishp)
enddo
!
xptc(1, isc) = xsc
xptc(2, isc) = ysc
!
xptic(1, isc) = xsci
xptic(2, isc) = ysci
enddo
!
!...Part I: Construct the nodal coordinates
!
xptsg(1, 1:6) = xpt(1, 1:6)
xptsg(2, 1:6) = xpt(2, 1:6)
!
xptisg(1, 1:6) = xpti(1, 1:6)
xptisg(2, 1:6) = xpti(2, 1:6)

!
!...Part II: Get mass matrix for four subgrids...
!

!...Mass average specific volume
rhomc = unkno(1,1)
!
do isg = 1, 4
!...Coordinates
xptsc(1, 1:3) = xptsg(1, iptsg(1:3, isg))
xptsc(2, 1:3) = xptsg(2, iptsg(1:3, isg))

xptsc(1, 4:6) = xptc(1, iptsc(1:3, isg))
xptsc(2, 4:6) = xptc(2, iptsc(1:3, isg))
!

!...Initialze parameters...
f0 = 0.d0
masel = 0.d0
volel = 0.d0
!
do ig =1,ngausd
!
r = posi(1,ig)
s = posi(2,ig)
wi = weigh(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shps(1) = 1.d0-r-s
shps(2) = r
shps(3) = s

!...Get the mapping gauss points in the standard triangle...
rg = 0.d0
sg = 0.d0

do ishp = 1, 3
rg = rg + shps(ishp)*xvt(iptsg(ishp, isg))
sg = sg + shps(ishp)*yvt(iptsg(ishp, isg))
enddo
!
!Basis function
btsg(1) = 1.d0
btsg(2) = (rg-rc)/dr
btsg(3) = (sg-sc)/ds

!...DGP2
if(npoly.eq.2)then
btsg(4) = 0.5d0*btsg(2)*btsg(2) - btho(1)
btsg(5) = 0.5d0*btsg(3)*btsg(3) - btho(2)
btsg(6) =       btsg(2)*btsg(3) - btho(3)
endif
!
rhomsg = 0.d0
do ideg = 1, mdegr
rhomsg = rhomsg + unkno(ideg,1)*btsg(ideg)
enddo
!
if(nlimi.eq.6)then
rhomsg = rhomc + aflim(1)*(rhomsg - rhomc)
endif


!...  shape function & its derivatives w.r.t. reference coordinates
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
do ishp = 1, nptri
dxdr = dxdr + dspr(ishp)*xptsc(1,ishp)
dxds = dxds + dsps(ishp)*xptsc(1,ishp)

dydr = dydr + dspr(ishp)*xptsc(2,ishp)
dyds = dyds + dsps(ishp)*xptsc(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco

!...Mass and volume
masel = masel + 1.d0/rhomsg*djaco
volel = volel + djaco
!
!if(ielem.eq.1) print*,'djaco',f0,wi,dxdr,dyds,dydr,dxds,xpqsc(1:2,:)
enddo

!
if(npoly==1)then
matin(4,isg) = 1.d0/f0
elseif(npoly==2)then
matin(16,isg) = 1.d0/f0
endif

!...The avearged specific volume at the sub-cell...
utsgc(1, isg) = volel/masel

enddo
!
!...Part III: Get rhs for the four subgrids...
!
rhsel = 0.d0

do isg = 1, 4

!...Initial coordinates
xptisc(1, 1:3) = xptisg(1, iptsg(1:3, isg))
xptisc(2, 1:3) = xptisg(2, iptsg(1:3, isg))

xptisc(1, 4:6) = xptic(1, iptsc(1:3, isg))
xptisc(2, 4:6) = xptic(2, iptsc(1:3, isg))
!
do ig =1,ngausd
!
r = posi(1,ig)
s = posi(2,ig)
wi = weigh(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
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
do ishp = 1, nptri
dxdr = dxdr + dspr(ishp)*xptisc(1,ishp)
dxds = dxds + dsps(ishp)*xptisc(1,ishp)

dydr = dydr + dspr(ishp)*xptisc(2,ishp)
dyds = dyds + dsps(ishp)*xptisc(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!
b(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nptri
xgaus = xgaus + shp(ishp)*xptisc(1,ishp)
ygaus = ygaus + shp(ishp)*xptisc(2,ishp)
enddo
!
call getrhoig_triacurv2(rhoi, xpti, xgaus, ygaus)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*b(1)*djaco
!if(isg.eq.4)print*,'isgideg',rhsel(6, isg),b(6),b(2),b(3), geoq_sub(7, isg)
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

!if(isg.eq.4.and.id.eq.1) print*,'isg4',m(1, :),rhsel(:,4),unint(1)
unksgt(1, isg)= unint(1)
enddo
!
end subroutine  getrhosubgt_average2


!
!...Riemann input for linear tria using gauss rule similar to the Eulerian framework...
!
subroutine getriem_tria_gausseuler(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelg, 1:ntria+nbfac), intent(in)::gelag
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),            intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),             intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!
real*8, dimension(1:2, 1:2, 1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
!
real*8, dimension(1:2, 1:2, 1:ngelg, 1:ntria),       intent(out)::munaclt
real*8, dimension(1:ndimn,  1:ngelg, 1:ntria), intent(out)::munault
real*8, dimension(1:ndimn,  1:ngelg, 1:ntria), intent(out)::snsigmlt
!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvfac) :: ipf
integer, dimension(nvfac, 3) ::fglvt
integer, dimension(ngausf,3):: fglgt
!...local real array
real*8,dimension(1:ndegr)::bg, bgv
real*8,dimension(1:nq,1:ngausf)::unkng
real*8::aujmp(1:3, 1:ngausf)
real*8::vnorm(1:3, 1:ngelg)
real*8::posit(1:2, 1:ngelg), positr(1:2, 1:ngelg)
real*8::sigmg(1:2, 1:2, 1:ngausf)
real*8,dimension(1:ngausf)::murie
real*8,dimension(1:nvtri):: xvt,  yvt

real*8,dimension(1:nvtri):: rcoet
real*8,dimension(1:ndimn, 1:nvtri) :: xpt,xpht
real*8,dimension(1:ndimn, 1:nvtri) :: xpti
real*8, dimension(1:nvtri):: shp, dspr, dsps
!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
real*8::posif(ngausf),weighf(ngausf)
!
real*8::eps,c00,c05,c10,c20
real*8::rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhog,rhomg,ug,vg,eg, pg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, rg, sg,rcv,scv
real*8:: dudr, duds, dvdr, dvds
real*8::acnx, acny
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Specifiy some gauss points...
!

!...Give gaussian position and weight...
call ruqope_lobatto(1, ngausf, posif, weighf)

!...Local vertex No. of gauss points in one unit ...
fglvt(1, 1) = 1;  fglvt(2, 1) = 2; !fglvt(3, 1) = 4;
fglvt(1, 2) = 2;  fglvt(2, 2) = 3; !fglvt(3, 2) = 5;
fglvt(1, 3) = 3;  fglvt(2, 3) = 1; !fglvt(3, 3) = 6;

!...Local gauss point No. of any gauss point in one face...
fglgt(1, 1) = 1;  fglgt(2, 1) = 2;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6;

!...6-nodes tria
posit(1, 1) = 0.d0; posit(2, 1) = 0.d0;
posit(1, 2) = 1.d0; posit(2, 2) = 0.d0;
posit(1, 3) = 1.d0; posit(2, 3) = 0.d0;
posit(1, 4) = 0.d0; posit(2, 4) = 1.d0;
posit(1, 5) = 0.d0; posit(2, 5) = 1.d0;
posit(1, 6) = 0.d0; posit(2, 6) = 0.d0;

!...6-nodes tria
positr(1, 1) = 0.5d0-sqrt(3.d0)/6.d0; positr(2, 1) =0.d0;
positr(1, 2) = 0.5d0+sqrt(3.d0)/6.d0; positr(2, 2) =0.d0;
positr(1, 3) = 0.5d0+sqrt(3.d0)/6.d0; positr(2, 3) =0.5d0-sqrt(3.d0)/6.d0;
positr(1, 4) = 0.5d0-sqrt(3.d0)/6.d0; positr(2, 4) =0.5d0+sqrt(3.d0)/6.d0;
positr(1, 5) = 0.d0; positr(2, 5) = 0.5d0+sqrt(3.d0)/6.d0;
positr(1, 6) = 0.d0; positr(2, 6) = 0.5d0-sqrt(3.d0)/6.d0;

!...The coordinates of a local cell
xvt(1) =  0.d0; yvt(1) = 0.d0
xvt(2) =  1.d0; yvt(2) = 0.d0
xvt(3) =  0.d0; yvt(3) = 1.d0

dr = .5d0
ds = .5d0
!
!...Part II: Loop over every quad...
!
do 250 ie = 1,ntria !...(1)ie = 1,ntria
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie

!...shape functions
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Give the normal vector of every face...
vnorm(1:3,  1:ngelg) = gelag(1:3, 1:ngelg, ie);

!...cell averaged value...
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

!...II.1: Loop over every face...

do ifa =1, 3
!...zero out unkno
unkng = 0.d0

do ig   = 1,  2
rg = posit(1, fglgt(ig, ifa))
sg = posit(2, fglgt(ig, ifa))

bg(1) = 1.d0
bg(2) = (rg-rc)/dr
bg(3) = (sg-sc)/ds

do ideg = 1,mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq,ielem)*bg(ideg)
enddo
!
!if(ie.eq.1)print*,'triaunkno',ie,ielem,unkno(1:3,2:3,ielem)
!
if(ndens.eq.1)then
rhog  = 1.d0/unkng(1, ig)
elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bgv(1) = 1.d0
bgv(2) = (rg-rcv)/dr
bgv(3) = (sg-scv)/ds

unkng(1, ig) =0.d0
do ideg = 1,mdegr
unkng(1, ig) = unkng(1, ig) + unkno(ideg,1,ielem)*bgv(ideg)
enddo
rhog = unkng(1, ig)
endif
!
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
pg = max(eps, (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2)))
!pg = (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2))

!...Limiter
if(nlimi.eq.6)then
if(ndens.eq.1)then
rhomg = rhomc + aflim(1, ielem)*(unkng(1, ig) - rhomc)
rhog = 1.d0/rhomg
elseif(ndens.eq.3)then
rhog = 1.d0/rhomc + aflim(1, ielem)*(rhog - 1.d0/rhomc)
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
ug = unkno(1,2,ielem)  + dudr*bg(2) + duds*bg(3)
vg = unkno(1,3,ielem)  + dvdr*bg(2) + dvds*bg(3)
!
pg = pctr + aflim(4, ielem)*(pg - pctr)

!...updtae unkng(2:3,:)
unkng(2, ig) = ug
unkng(3 ,ig) = vg
endif

!...Get stress tensor at nodes
sigmg(1, 1, ig) = -pg
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pg

!...Get the a_c (unit vector)
aujmp(1:2, ig) = vlave(1:2, ipt(fglvt(ig, ifa))) - unkng(2:3, ig)
!
acnx = aujmp(1, ig)
acny = aujmp(2, ig)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig) = 1.e-11!0.d0;
else
aujmp(1:2, ig) = aujmp(1:2, ig)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig) = sqrt(acnx**2 + acny**2)

enddo !....ig

!...Get sound speed at the center...
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedence coefficient...
do ig   = 1, 2
dux= vlave(1, ipt(fglvt(ig, ifa)))-unkng(2, ig)
duy= vlave(2, ipt(fglvt(ig, ifa)))-unkng(3, ig)
deltu = sqrt(dux**2 + duy**2)
deltu = abs(dux*vnorm(1,fglgt(ig, ifa)) + duy*vnorm(2,fglgt(ig, ifa)))
murie(ig) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
enddo

!...Get the summed denominator cooefficients sum(mu*n*a_c)
do ig  = 1, 2!
!
call getriecoef_matrixnew(murie(ig), vnorm(3, fglgt(ig, ifa)), vnorm(1:2, fglgt(ig, ifa)), aujmp(1:3, ig), &
unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ig), vnorm(3, fglgq(ig, ifa)), vnorm(1:2, fglgq(ig, ifa)), aujmp(1:2, ig), &
!unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipt(fglvt(ig, ifa))) = munacn(1:2, 1, ipt(fglvt(ig, ifa))) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipt(fglvt(ig, ifa))) = munacn(1:2, 2, ipt(fglvt(ig, ifa))) + munacn_rie(1:2, 2)
!
munacu(1:2, ipt(fglvt(ig, ifa))) = munacu(1:2, ipt(fglvt(ig, ifa))) + munacu_rie(1:2)
!
snsigm(1:2, ipt(fglvt(ig, ifa))) = snsigm(1:2, ipt(fglvt(ig, ifa))) + snsigm_rie(1:2)!
!
munaclt(1:2, 1, fglgt(ig, ifa), ie) =  munacn_rie(1:2, 1)
munaclt(1:2, 2, fglgt(ig, ifa), ie) =  munacn_rie(1:2, 2)
!
munault(1:2, fglgt(ig, ifa), ie)    =  munacu_rie(1:2)
!
snsigmlt(1:2,fglgt(ig, ifa), ie)    =  snsigm_rie(1:2)
!
!if(ipt(fglvt(ig, ifa)).eq.2) print*,'ep', murie(ig),ie,ifa,ig,fglvt(ig, ifa),sigmg(1:2, 1:2, ig),vnorm(1:3, fglgt(ig, ifa)),&
!munacu(1, ipt(fglvt(ig, ifa))),munacu_rie(1)
!
enddo
!
!...Get the new snsigmlqr...
!
!...zero out unkno
unkng = 0.d0
!
do ig   = 1,  2
rg = posit(1, fglgt(ig, ifa))
sg = posit(2, fglgt(ig, ifa))
!
bg(1) = 1.d0
bg(2) = (rg-rc)/dr
bg(3) = (sg-sc)/ds
!
do ideg = 1,mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq,ielem)*bg(ideg)
enddo

!
if(ndens.eq.1)then
rhog  = 1.d0/unkng(1, ig)
elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bgv(1) = 1.d0
bgv(2) = (rg-rcv)/dr
bgv(3) = (sg-scv)/ds

unkng(1, ig) =0.d0
do ideg = 1,mdegr
unkng(1, ig) = unkng(1, ig) + unkno(ideg,1,ielem)*bgv(ideg)
enddo

rhog = unkng(1, ig)
endif
!
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
pg = max(eps, (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2)))
!pg = (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2))

!...Limiter
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomg = rhomc + aflim(1, ielem)*(unkng(1, ig) - rhomc)
rhog = 1.d0/rhomg
elseif(ndens.eq.3)then
rhog = 1.d0/rhomc + aflim(1, ielem)*(rhog - 1.d0/rhomc)
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
ug = unkno(1,2,ielem)  + dudr*bg(2) + duds*bg(3)
vg = unkno(1,3,ielem)  + dvdr*bg(2) + dvds*bg(3)
!
pg = pctr + aflim(4, ielem)*(pg - pctr)

!...updtae unkng(2:3,:)
unkng(2, ig) = ug
unkng(3 ,ig) = vg
endif

!...Get stress tensor at nodes
sigmg(1, 1, ig) = -pg
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pg
enddo

!
do ig  = 1, 2
call getriecoef_matrixnew(murie(ig), vnorm(3, fglgt(ig, ifa)), vnorm(1:2, fglgt(ig, ifa)), aujmp(1:3, ig), &
unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ig), vnorm(3, fglgq(ig, ifa)), vnorm(1:2, fglgq(ig, ifa)), aujmp(1:2, ig), &
!unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
snsigmlt(1:2,fglgt(ig, ifa), ie)    =  snsigm_rie(1:2)
enddo
!
enddo !...II.1: loop over every face
!
250 enddo  !...II: loop over every triangle

end subroutine getriem_tria_gausseuler
!
!...Riemann input for linear quad using gauss rule similar to the Eulerian framework...
!
subroutine getriem_quad_gausseuler(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),           intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),         intent(in)::afvec
!
real*8, dimension(1:2, 1:2, 1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
!
real*8, dimension(1:2, 1:2, 1:ngelgq, 1:nquad),       intent(out)::munaclq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad)::snsigmlqr
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(out)::snsigmlq
!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer, dimension(nvfac, 4)::fglvq
integer, dimension(ngausf, 4):: fglgq
!...local real array
real*8,dimension(1:ndegr)::bg, bgv
real*8,dimension(1:nq,1:ngausf)::unkng
real*8::aujmp(1:3, 1:ngausf)
real*8::vnorm(1:3, 1:ngelgq)
real*8::posiq(1:2, 1:ngelgq), posiqr(1:2, 1:ngelgq)
real*8::sigmg(1:2, 1:2, 1:ngausf)
real*8,dimension(1:ngausf)::murie
real*8,dimension(1:nvqua):: xvq,  yvq
!...arraies for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
real*8::posif(ngausf),weighf(ngausf)
!
real*8::eps,c00,c05,c10,c20
real*8::rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhog,rhomg,ug,vg,eg, pg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, rg, sg,rcv,scv
real*8:: dudr, duds, dvdr, dvds
real*8::acnx, acny
!
data eps   / 1.0d-6 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
call ruqope_lobatto(1, ngausf, posif, weighf)
!
!...Local vertex No. of gauss points in one unit ...
!
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; !fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; !fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; !fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; !fglvq(3, 4) = 8;
!
!...Local gauss point No. of any gauss point in one face...
!
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
!
!...8-nodes quad
!
posiqr(1, 1) =-0.5773502691896257645091488d0; posiqr(2, 1) =-1.d0;
posiqr(1, 2) = 0.5773502691896257645091488d0; posiqr(2, 2) =-1.d0;
posiqr(1, 3) = 1.d0; posiqr(2, 3) =-0.5773502691896257645091488d0;
posiqr(1, 4) = 1.d0; posiqr(2, 4) = 0.5773502691896257645091488d0;
posiqr(1, 5) = 0.5773502691896257645091488d0; posiqr(2, 5) = 1.d0;
posiqr(1, 6) =-0.5773502691896257645091488d0; posiqr(2, 6) = 1.d0;
posiqr(1, 7) =-1.d0; posiqr(2, 7) = 0.5773502691896257645091488d0;
posiqr(1, 8) =-1.d0; posiqr(2, 8) =-0.5773502691896257645091488d0;

!...The coordinates of a local cell

xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
dr = 1.d0
ds = 1.d0
!
!...Triangle...
!
do 250 ie = 1,nquad !...(1)ie = 1,nelem
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria

!...shape functions

rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Give the normal vector of every face...

vnorm(1:3,  1:ngelgq) = gelagq(1:3, 1:ngelgq, ie);
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

!...Give the normal vector of every face...

do ifa =1, 4

!...zero out unkno
unkng = 0.d0

do ig   = 1,  2
!
rg = posiqr(1, fglgq(ig, ifa))
sg = posiqr(2, fglgq(ig, ifa))
!
bg(1) = 1.d0
bg(2) = (rg-rc)/dr
bg(3) = (sg-sc)/ds
!
do ideg = 1,mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq,ielem)*bg(ideg)
enddo
!
if(ndens.eq.1)then
rhog  = 1.d0/unkng(1, ig)
elseif(ndens.eq.3)then
!
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bgv(1) = 1.d0
bgv(2) = (rg-rcv)/dr
bgv(3) = (sg-scv)/ds

unkng(1, ig) =0.d0
do ideg = 1,mdegr
unkng(1, ig) = unkng(1, ig) + unkno(ideg,1,ielem)*bgv(ideg)
enddo

rhog = unkng(1, ig)
endif
!
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
pg = max(eps, (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2)))
!pg = (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2))

if(nlimi.eq.6)then
pg = pctr + aflim(4, ielem)*(pg - pctr)
endif
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig) = -pg
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pg
!
enddo

!...zero out unkno
unkng = 0.d0

do ig   = 1,  2
!
rg = posiq(1, fglgq(ig, ifa))
sg = posiq(2, fglgq(ig, ifa))
!
bg(1) = 1.d0
bg(2) = (rg-rc)/dr
bg(3) = (sg-sc)/ds
!
do ideg = 1,mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq,ielem)*bg(ideg)
enddo
!
!if(ip(fagsp(ig, ifa)).eq.155) print*,'ndegr',unkno(1:3,2,ie),bg(1:3),ie,ifa,rg,sg
!
if(ndens.eq.1)then
rhog  = 1.d0/unkng(1, ig)
elseif(ndens.eq.3)then
!
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bgv(1) = 1.d0
bgv(2) = (rg-rcv)/dr
bgv(3) = (sg-scv)/ds

unkng(1, ig) =0.d0
do ideg = 1,mdegr
unkng(1, ig) = unkng(1, ig) + unkno(ideg,1,ielem)*bgv(ideg)
enddo

rhog = unkng(1, ig)
endif
!
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
!pg = max(eps, (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2)))
!pg = (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2))
!
!if(ie.eq.9) print*,'pg',rhog,unkno(1:3,1,ie),bg(1:3),rg,sg
!
!...Limiter
!
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomg = rhomc + aflim(1, ielem)*(unkng(1, ig) - rhomc)
rhog = 1.d0/rhomg
elseif(ndens.eq.3)then
rhog = 1.d0/rhomc + aflim(1, ielem)*(rhog - 1.d0/rhomc)
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
ug = unkno(1,2,ielem)  + dudr*bg(2) + duds*bg(3)
vg = unkno(1,3,ielem)  + dvdr*bg(2) + dvds*bg(3)
!
!pg = pctr + aflim(4, ielem)*(pg - pctr)

!...updtae unknv(2:3,:)

unkng(2, ig) = ug
unkng(3 ,ig) = vg
!
endif
!
!...Get stress tensor at nodes
!
!sigmg(1, 1, ig) = -pg
!sigmg(1, 2, ig) = 0.d0
!sigmg(2, 1, ig) = 0.d0
!sigmg(2, 2, ig) = -pg
!
!...Get the a_c (unit vector)
!
aujmp(1:2, ig) = vlave(1:2, ipq(fglvq(ig, ifa))) - unkng(2:3, ig)
!
acnx = aujmp(1, ig)
acny = aujmp(2, ig)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig) = 1.e-11!0.d0;
else
aujmp(1:2, ig) = aujmp(1:2, ig)/sqrt(acnx**2 + acny**2)
endif
!
aujmp(3, ig) = sqrt(acnx**2 + acny**2)
!
enddo !....ig

!...Get sound speed at the center...

sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedence coefficient...

do ig   = 1, 2
dux= vlave(1, ipq(fglvq(ig, ifa)))-unkng(2, ig)
duy= vlave(2, ipq(fglvq(ig, ifa)))-unkng(3, ig)
deltu = sqrt(dux**2 + duy**2)
deltu = 1.d0*abs(dux*vnorm(1,fglgq(ig, ifa)) + duy*vnorm(2,fglgq(ig, ifa)))
murie(ig) = rhoct*sdctr! + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
enddo
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!
do ig  = 1, 2!
!
call getriecoef_matrixnew(murie(ig), vnorm(3, fglgq(ig, ifa)), vnorm(1:2, fglgq(ig, ifa)), aujmp(1:3, ig), &
unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ig), vnorm(3, fglgq(ig, ifa)), vnorm(1:2, fglgq(ig, ifa)), aujmp(1:2, ig), &
!unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipq(fglvq(ig, ifa))) = munacn(1:2, 1, ipq(fglvq(ig, ifa))) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipq(fglvq(ig, ifa))) = munacn(1:2, 2, ipq(fglvq(ig, ifa))) + munacn_rie(1:2, 2)
!
munacu(1:2, ipq(fglvq(ig, ifa))) = munacu(1:2, ipq(fglvq(ig, ifa))) + munacu_rie(1:2)
!
snsigm(1:2,ipq(fglvq(ig, ifa))) = snsigm(1:2, ipq(fglvq(ig, ifa))) + snsigm_rie(1:2)!
!
munaclq(1:2, 1, fglgq(ig, ifa), ie) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fglgq(ig, ifa), ie) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fglgq(ig, ifa), ie)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fglgq(ig, ifa), ie)    =  snsigm_rie(1:2)
!
!if(ipq(fglvq(ig, ifa)).eq.21) print*,'epq', murie(ig),ie,ifa,ig,fglvq(ig, ifa),sigmg(1:2, 1:2, ig),vnorm(1:3, fglgq(ig, ifa)),&
!munacu(1,ipq(fglvq(ig, ifa))),munacu_rie(1),snsigm(1,ipq(fglvq(ig, ifa))),snsigm_rie(1),sigmg(1, 1, ig)
!
enddo
!
!...Get the new snsigmlqr...
!
!...zero out unkno
unkng = 0.d0
!
do ig   = 1,  2
!
rg = posiqr(1, fglgq(ig, ifa))
sg = posiqr(2, fglgq(ig, ifa))
!
bg(1) = 1.d0
bg(2) = (rg-rc)/dr
bg(3) = (sg-sc)/ds
!
do ideg = 1,mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq,ielem)*bg(ideg)
enddo
!
!if(ip(fagsp(ig, ifa)).eq.155) print*,'ndegr',unkno(1:3,2,ie),bg(1:3),ie,ifa,rg,sg
!
if(ndens.eq.1)then
rhog  = 1.d0/unkng(1, ig)
elseif(ndens.eq.3)then
!
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bgv(1) = 1.d0
bgv(2) = (rg-rcv)/dr
bgv(3) = (sg-scv)/ds

unkng(1, ig) =0.d0
do ideg = 1,mdegr
unkng(1, ig) = unkng(1, ig) + unkno(ideg,1,ielem)*bgv(ideg)
enddo

rhog = unkng(1, ig)
endif
!
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
!pg = max(eps, (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2)))
pg = (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2))
!
!...Limiter
!
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomg = rhomc + aflim(1, ielem)*(unkng(1, ig) - rhomc)
rhog = 1.d0/rhomg
elseif(ndens.eq.3)then
rhog = 1.d0/rhomc + aflim(1, ielem)*(rhog - 1.d0/rhomc)
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
ug = unkno(1,2,ielem)  + dudr*bg(2) + duds*bg(3)
vg = unkno(1,3,ielem)  + dvdr*bg(2) + dvds*bg(3)
!
pg = pctr + aflim(4, ielem)*(pg - pctr)

!...updtae unknv(2:3,:)

unkng(2, ig) = ug
unkng(3 ,ig) = vg
endif
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig) = -pg
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pg
!
enddo
!
do ig  = 1, 2!
!
call getriecoef_matrixnew(murie(ig), vnorm(3, fglgq(ig, ifa)), vnorm(1:2, fglgq(ig, ifa)), aujmp(1:3, ig), &
unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ig), vnorm(3, fglgq(ig, ifa)), vnorm(1:2, fglgq(ig, ifa)), aujmp(1:2, ig), &
!unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
snsigmlqr(1:2,fglgq(ig, ifa), ie)    =  snsigm_rie(1:2)
!
enddo

!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,nelem!
!
end subroutine getriem_quad_gausseuler
!
!...Get the nodal velocity based on gauss integration in the Eulerian framework...
!
subroutine getndvelo_lag_gausseulframe(gflag,gelag,gelagq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, ufgpt, ufgaus, fstart, fstarq, aflim, afvec, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3, 1:ngelg, 1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:ngelg, 1:ntria),   intent(out)::fstart !...Riemann forces
real*8,dimension(1:ndimn,1:ngelgq, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngelgq,1:nquad),       intent(out)::ufgaus
real*8,dimension(1:2,1:ngelg, 1:ntria),       intent(out) ::ufgpt

integer:: itime,ip
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)
integer, dimension(3, 4)::fglvq
integer, dimension(3, 3)::fglvt
!
integer, dimension(ngausf, 4):: fglgq
integer, dimension(ngausf, 3):: fglgt

!...local real array
integer::ipsin(4)
real*8::munacnsin(4,2),munacusin(2, 2),snsigmsin(2,2)
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:ndimn,1:nvqua,1:nquad)::vnulq
real*8,  dimension(1:nquad)::gqdmp
real*8::munaci(2, 2)
real*8::posif(ngausf), weighf(ngausf)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2
real*8::r,shp1,shp2,shp3
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: munacn(:,:,:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: munacu(:,:), snsigm(:,:)
real*8,allocatable:: munaclt(:,:,:,:),munaclq(:,:,:,:)
real*8,allocatable:: snsigmlt(:,:,:), munault(:,:,:)
real*8,allocatable:: snsigmlq(:,:,:), munaulq(:,:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:2, 1:2, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munaclt(1:2, 1:2, 1:ngelg, 1:ntria), munault(1:ndimn, 1:ngelg,  1:ntria),&
snsigmlt(1:ndimn, 1:ngelg,  1:ntria))
allocate (munaclq(1:2, 1:2, 1:ngelgq, 1:nquad), munaulq(1:ndimn, 1:ngelgq,  1:nquad),&
snsigmlq(1:ndimn, 1:ngelgq,  1:nquad))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
!
!...Part I: Specify some gauss points...
!

!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;

!...Local vertex No. of gauss points in one unit ...
fglvt(1, 1) = 1;  fglvt(2, 1) = 2; fglvt(3, 1) = 4;
fglvt(1, 2) = 2;  fglvt(2, 2) = 3; fglvt(3, 2) = 5;
fglvt(1, 3) = 3;  fglvt(2, 3) = 1; fglvt(3, 3) = 6;
!
!...Give gaussian position and weight...
!
call ruqope_lobatto(1, ngausf, posif, weighf)
!
!...Local gauss point No. of any gauss point in one face...
!
if(ngausf.eq.2)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8;
!
fglgt(1, 1) = 1;  fglgt(2, 1) = 2;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6;
!
elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
!
fglgt(1, 1) = 1;  fglgt(2, 1) = 2; fglgt(3, 1) = 7; fglgt(4, 1) =  8;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4; fglgt(3, 2) = 9; fglgt(4, 2) = 10;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6; fglgt(3, 3) =11; fglgt(4, 3) = 12;
!
endif
!
!...Zero out vlave
!
vlave = 0.d0
indnd = 0
!
!...Mark the boundary nodes...
!
if(ncase.eq.2)then
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!if(bface(3,ifa).eq.30)then
indnd(ipf(1:nvfac)) = 1
!endif
enddo
endif
!
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
indnd(ipf(1:nvfac)) = 1
endif
enddo

!...Get averaged velocity at nodes...
do 950 ifa = 1 , nbfac
 ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
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
!...Part II: Loop to get the information from Riemann solver
!
do iloop= 1, 1

!...Give vlave
vlave= ustar
vnulq = 0.d0

!...Zero out munacn
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0

!...Tria
if(ntria.gt.0) call getriem_tria_gausseuler(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold, aflim, afvec)

!...Quad
if(nquad.gt.0) call getriem_quad_gausseuler(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)

!...Boundary condition...
call getbc_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)

!...Periodic boundary condition for 1D isentropic Sin problem...
if(ncase.eq.12)then
 do ifa = 1, nbfac
 if(bface(3, ifa).eq.31)then !...Periodic boundary...
!
  ipsin(1:2) = bface(1:2, ifa); ipsin(3:4) = bface(1:2, bface(5, ifa))
!
  munacnsin(1,1) = munacn(1,1,ipsin(2)) + munacn(1,1,ipsin(3))
  munacnsin(2,1) = munacn(1,2,ipsin(2)) + munacn(1,2,ipsin(3))
  munacnsin(3,1) = munacn(2,1,ipsin(2)) + munacn(2,1,ipsin(3))
  munacnsin(4,1) = munacn(2,2,ipsin(2)) + munacn(2,2,ipsin(3))

  munacnsin(1,2) = munacn(1,1,ipsin(1)) + munacn(1,1,ipsin(4))
  munacnsin(2,2) = munacn(1,2,ipsin(1)) + munacn(1,2,ipsin(4))
  munacnsin(3,2) = munacn(2,1,ipsin(1)) + munacn(2,1,ipsin(4))
  munacnsin(4,2) = munacn(2,2,ipsin(1)) + munacn(2,2,ipsin(4))
!
  munacusin(1:2, 1) = munacu(1:2, ipsin(2)) + munacu(1:2, ipsin(3))
  munacusin(1:2, 2) = munacu(1:2, ipsin(1)) + munacu(1:2, ipsin(4))
!
  snsigmsin(1:2, 1) = snsigm(1:2, ipsin(2)) + snsigm(1:2, ipsin(3))
  snsigmsin(1:2, 2) = snsigm(1:2, ipsin(1)) + snsigm(1:2, ipsin(4))
!
  munacn(1,1,ipsin(2)) = munacnsin(1, 1); munacn(1,1,ipsin(3)) = munacnsin(1, 1);
  munacn(1,2,ipsin(2)) = munacnsin(2, 1); munacn(1,2,ipsin(3)) = munacnsin(2, 1);
  munacn(2,1,ipsin(2)) = munacnsin(3, 1); munacn(2,1,ipsin(3)) = munacnsin(3, 1);
  munacn(2,2,ipsin(2)) = munacnsin(4, 1); munacn(2,2,ipsin(3)) = munacnsin(4, 1);

  munacn(1,1,ipsin(1)) = munacnsin(1, 2); munacn(1,1,ipsin(4)) = munacnsin(1, 2);
  munacn(1,2,ipsin(1)) = munacnsin(2, 2); munacn(1,2,ipsin(4)) = munacnsin(2, 2);
  munacn(2,1,ipsin(1)) = munacnsin(3, 2); munacn(2,1,ipsin(4)) = munacnsin(3, 2);
  munacn(2,2,ipsin(1)) = munacnsin(4, 2); munacn(2,2,ipsin(4)) = munacnsin(4, 2);
!
  munacu(1:2, ipsin(2)) =  munacusin(1:2, 1); munacu(1:2, ipsin(3)) =  munacusin(1:2, 1);
  munacu(1:2, ipsin(1)) =  munacusin(1:2, 2); munacu(1:2, ipsin(4)) =  munacusin(1:2, 2);
!
  snsigm(1:2, ipsin(2)) =  snsigmsin(1:2, 1); snsigm(1:2, ipsin(3)) =  snsigmsin(1:2, 1);
  snsigm(1:2, ipsin(1)) =  snsigmsin(1:2, 2); snsigm(1:2, ipsin(4)) =  snsigmsin(1:2, 2);
 endif
 enddo
endif !if(ncase.eq.12)then

!....Update the velocity at the vertex...
do ipoin = 1, npoin
 if(indnd(ipoin).eq.0)then
  detma = munacn(1, 1, ipoin)*munacn(2, 2, ipoin) - munacn(2, 1, ipoin)*munacn(1, 2, ipoin)
  munaci(1, 1) = munacn(2, 2, ipoin)/detma
  munaci(1, 2) =-munacn(1, 2, ipoin)/detma
  munaci(2, 1) =-munacn(2, 1, ipoin)/detma
  munaci(2, 2) = munacn(1, 1, ipoin)/detma
!
  rhsu1 = munacu(1, ipoin) - snsigm(1, ipoin) - fpres(1, ipoin)
  rhsu2 = munacu(2, ipoin) - snsigm(2, ipoin) - fpres(2, ipoin)
!
  ustar(1, ipoin) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
  ustar(2, ipoin) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!+++++++++output for debugging++++++++++++
 ! if(ipoin.eq.21)print*,'ustar',ipoin,ustar(1:2,ipoin),munacu(1:2, ipoin),snsigm(1:2,ipoin),fpres(1:2,ipoin)
 endif
enddo

!...Get the vertex velocity at boundary....
do 900 ifa = 1 , nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)

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

!...Specify boundary velocity
 if(ncase.eq.13)then !...Saltzman
  ustar(1,ipf(1:nvfac)) = 1.d0
  ustar(2,ipf(1:nvfac)) = 0.d0
 endif
endif

!...Analytic boundary velocity for TGV
if(ncase.eq.-1)then
ustar(1, ipf(1)) = sin(pi*coord(1,ipf(1)))*cos(pi*coord(2,ipf(1)))
ustar(2, ipf(1)) =-cos(pi*coord(1,ipf(1)))*sin(pi*coord(2,ipf(1)))

ustar(1, ipf(2)) = sin(pi*coord(1,ipf(2)))*cos(pi*coord(2,ipf(2)))
ustar(2, ipf(2)) =-cos(pi*coord(1,ipf(2)))*sin(pi*coord(2,ipf(2)))
endif

900 enddo

!...Specify ufgaus and ufgpt initially...
if(iloop.eq.1)then
 if(ngausf.eq.4)then

!...tria
 do 310 ie=1, ntria !...(1) ie=1, ntria
 ipt(1:nvtri) = iptri(1:nvtri,ie)

 do ifa = 1,3
 ufgpt(1:ndimn,fglgt(1, ifa),  ie) = ustar(1:2, ipt(fglvt(1, ifa)))
 ufgpt(1:ndimn,fglgt(2, ifa),  ie) = ustar(1:2, ipt(fglvt(2, ifa)))
!
 r = -0.4472135954999579392818d0
 shp1 =  -0.5d0*(1.d0-r)*r
 shp2 =   0.5d0*(1.d0+r)*r
 shp3 =         (1.d0+r)*(1.d0-r)
!
 ufgpt(1:ndimn,fglgt(3, ifa),  ie) = shp1*ustar(1:2, ipt(fglvt(1, ifa))) + shp2*ustar(1:2, ipt(fglvt(2, ifa)))+&
 shp3*ustar(1:2, ipt(fglvt(3, ifa)))
!
 r =  0.4472135954999579392818d0
 shp1 =  -0.5d0*(1.d0-r)*r
 shp2 =   0.5d0*(1.d0+r)*r
 shp3 =         (1.d0+r)*(1.d0-r)
!
 ufgpt(1:ndimn,fglgt(4, ifa),  ie) = shp1*ustar(1:2, ipt(fglvt(1, ifa))) + shp2*ustar(1:2, ipt(fglvt(2, ifa)))+&
 shp3*ustar(1:2, ipt(fglvt(3, ifa)))
 enddo

310 enddo

!...quad
 do 300 ie=1, nquad !...(1)ifa=1,nquad
  ipq(1:nvqua) = ipqua(1:nvqua,ie)
 do ifa = 1,4
  ufgaus(1:ndimn,fglgq(1, ifa),  ie) = ustar(1:2, ipq(fglvq(1, ifa)))
  ufgaus(1:ndimn,fglgq(2, ifa),  ie) = ustar(1:2, ipq(fglvq(2, ifa)))
!
  r = posif(3)
  shp1 =   0.5d0*(1.d0-r)
  shp2 =   0.5d0*(1.d0+r)

  ufgaus(1:ndimn,fglgq(3, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))

  r = posif(4)
  shp1 =   0.5d0*(1.d0-r)
  shp2 =   0.5d0*(1.d0+r)
  ufgaus(1:ndimn,fglgq(4, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))
 enddo
300 enddo
 endif
endif

!...Get the normal velocity at the gauss point....
!...As of now, works for all-quad mesh...
 if(ngausf.eq.4)then
  call getvelo_mpt_gausseulframe(ufgaus,gelagq,intfac,ipqua,coord,geoel,unkno,indnd, &
  munaclq, munaulq, snsigmlq,afvec, aflim)
 endif

enddo
!
!...Part III:: Get the Riemann velocity and forces for face integral
!
!...III.1: Get ufgaus and ufgqt from ustar ...
!...Tria
 do 110 ie=1, ntria
  ipt(1:nvtri) = iptri(1:nvtri,ie)

 do ifa = 1,3
  ufgpt(1:ndimn,fglgt(1, ifa),  ie) = ustar(1:2, ipt(fglvt(1, ifa)))
  ufgpt(1:ndimn,fglgt(2, ifa),  ie) = ustar(1:2, ipt(fglvt(2, ifa)))
!
  if(ngausf.eq.4)then
  r = posif(2)
  shp1 =   0.5d0*(1.d0-r)
  shp2 =   0.5d0*(1.d0+r)
!if(ufgpt(1,fglgt(3, ifa),  ie).gt.1.d5)then
!ufgpt(1:ndimn,fglgt(3, ifa),  ie) = shp1*ustar(1:2, ipt(fglvt(1, ifa))) + shp2*ustar(1:2, ipt(fglvt(2, ifa)))
!endif
 r =  posif(3)
 shp1 =   0.5d0*(1.d0-r)
 shp2 =   0.5d0*(1.d0+r)
!if(ufgpt(1,fglgt(4, ifa),  ie).gt.1.d5)then
!ufgpt(1:ndimn,fglgt(4, ifa),  ie) = shp1*ustar(1:2, ipt(fglvt(1, ifa))) + shp2*ustar(1:2, ipt(fglvt(2, ifa)))
!endif
  endif
 enddo
110 enddo

!...Quad
do 100 ie=1, nquad
 ipq(1:nvqua) = ipqua(1:nvqua,ie)

 do ifa = 1,4
 ufgaus(1:ndimn,fglgq(1, ifa),  ie) = ustar(1:2, ipq(fglvq(1, ifa)))
 ufgaus(1:ndimn,fglgq(2, ifa),  ie) = ustar(1:2, ipq(fglvq(2, ifa)))
!
 if(ngausf.eq.4)then
 r = posif(2)
 shp1 =   0.5d0*(1.d0-r)
 shp2 =   0.5d0*(1.d0+r)
!if(ufgaus(1,fglgq(3, ifa),  ie).gt.1.d5)then
!ufgaus(1:ndimn,fglgq(3, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))
!endif
!
 r =  posif(3)
 shp1 =   0.5d0*(1.d0-r)
 shp2 =   0.5d0*(1.d0+r)
!if(ufgaus(1,fglgq(4, ifa),  ie).gt.1.d5)then
!ufgaus(1:ndimn,fglgq(4, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))
!endif
 endif
 enddo
100 enddo

!...III.2: Get the Riemann forces associated with every node...

!...Tria
!...Zero out fstart
fstart = 0.d0

do ie = 1, ntria
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
do ifa = 1, 3
do ig =1, ngausf
fstart(1, fglgt(ig, ifa), ie) = snsigmlt(1, fglgt(ig, ifa), ie) +&
munaclt(1,1, fglgt(ig, ifa), ie)*ufgpt(1, fglgt(ig, ifa), ie)+&
munaclt(2,1, fglgt(ig, ifa), ie)*ufgpt(2, fglgt(ig, ifa), ie)-&
munault(1, fglgt(ig, ifa), ie)
!
fstart(2, fglgt(ig, ifa), ie) = snsigmlt(2, fglgt(ig, ifa), ie) +&
munaclt(2,2,fglgt(ig, ifa), ie)*ufgpt(2, fglgt(ig, ifa), ie)+&
munaclt(1,2,fglgt(ig, ifa), ie)*ufgpt(1, fglgt(ig, ifa), ie)-&
munault(2, fglgt(ig, ifa), ie)
enddo
enddo
enddo

!...Quad
!...Zero out fstarq
fstarq = 0.d0

do ie = 1, nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
 do ifa = 1, 4
 do ig =1, ngausf
  fstarq(1, fglgq(ig, ifa), ie) = snsigmlq(1, fglgq(ig, ifa), ie) +&
  munaclq(1,1, fglgq(ig, ifa), ie)*ufgaus(1, fglgq(ig, ifa), ie)+&
  munaclq(2,1, fglgq(ig, ifa), ie)*ufgaus(2, fglgq(ig, ifa), ie)-&
  munaulq(1, fglgq(ig, ifa), ie)
!
  fstarq(2, fglgq(ig, ifa), ie) = snsigmlq(2, fglgq(ig, ifa), ie) +&
  munaclq(2,2,fglgq(ig, ifa), ie)*ufgaus(2, fglgq(ig, ifa), ie)+&
  munaclq(1,2,fglgq(ig, ifa), ie)*ufgaus(1, fglgq(ig, ifa), ie)-&
  munaulq(2, fglgq(ig, ifa), ie)
 enddo
 enddo
enddo
!
deallocate (munacn, bpres, fpres)
deallocate (munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lag_gausseulframe
!
!...Identify the gauss point location in the left cell
!
subroutine getpsfgl(ipqua, ipf, fgaus)
use constant
implicit none
integer, dimension(nvqua), intent(in):: ipqua
integer, dimension(nvfac), intent(in):: ipf
integer:: ig, nsum
integer, dimension(1:ngausf), intent(out)::fgaus
!
integer::fglgq(ngausf,4)
!
!...Specify fglgq
!
if(ngausf.eq.4)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
endif
!
!...Specify the No. of face gauss points in fglgq...
!
nsum = ngausf +3
!
do ig = 3, ngausf
!
if(ipqua(1).eq.ipf(1))then

fgaus(ig) = fglgq(ig,1)

elseif(ipf(1).eq.ipqua(2))then

fgaus(ig) = fglgq(ig,2)

elseif(ipf(1).eq.ipqua(3))then

fgaus(ig) = fglgq(ig,3)

elseif(ipf(1).eq.ipqua(4))then

fgaus(ig) = fglgq(ig,4)

endif
enddo
end subroutine getpsfgl

!
!...Identify the gauss point location in the right cell
!
subroutine getpsfgr(ipqua, ipf, fgaus)
use constant
implicit none
integer, dimension(nvqua), intent(in):: ipqua
integer, dimension(nvfac), intent(in):: ipf
integer:: ig, nsum
integer, dimension(1:ngausf), intent(out)::fgaus
!
integer::fglgq(ngausf,4)
!
!...Specify fglgq
!
if(ngausf.eq.4)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
endif
!
!...Specify the No. of face gauss points in fglgq...
!
nsum = ngausf +3
!
do ig = 3, ngausf
!
if(ipqua(1).eq.ipf(2))then

fgaus(ig) = fglgq(nsum-ig,1)

elseif(ipf(2).eq.ipqua(2))then

fgaus(ig) = fglgq(nsum-ig,2)

elseif(ipf(2).eq.ipqua(3))then

fgaus(ig) = fglgq(nsum-ig,3)

elseif(ipf(2).eq.ipqua(4))then

fgaus(ig) = fglgq(nsum-ig,4)

endif
enddo
end subroutine getpsfgr

!...Calculate the velocity at the middle point for Marie method in eulerian framework...
!
subroutine getvelo_mpt_gausseulframe(ufgaus,gelagq,intfac,ipqua,coord,geoel,unkno,indnd, &
munaclq, munaulq, snsigmlq, afvec, aflim)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
integer*4,dimension(1:npoin),                intent(in)::indnd
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngelgq,nquad),      intent(inout)::ufgaus
real*8, dimension(1:2, 1:2, 1:ngelgq, 1:nquad), intent(inout)::munaclq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(inout)::munaulq
real*8, dimension(1:ndimn, 1:ngelgq,  1:nquad), intent(inout)::snsigmlq
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(ngausf)::fgausl, fgausr
!
real*8::eps
real*8,dimension(1:2,1:ngelgq,nquad)::vlave
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::sigmg(1:2, 1:2, 1:ngausf, 1:2)
real*8::vngs(3, ngausf,2)
real*8::aujmp(1:3,1:ngausf,1:2)
real*8::murie(2)
real*8::vnorm(1:3,1:nvqua)
real*8::xvq(nvqua), yvq(nvqua), b(ndegr, ngausf)
real*8::posiq(1:2, 1:ngelgq)
real*8::posif(ngausf),weighf(ngausf)
integer,dimension(ngausf,4)::fglgq
real*8::xpf(1:2, 1:nvfac)
!
!...Riemann parameters...
!
real*8::munacn_rie(2,2), munacu_rie(2), snsigm_rie(2)
real*8::munacn(1:2, 1:2, 1:ngausf), munacu(1:2, 1:ngausf), snsigm(1:2, 1:ngausf)
real*8::munaci(2, 2)
!
!...Local real number
!
real*8::detma,dudr,duds,dvdr,dvds
real*8::pvtx,rhovt, rhomc, rhomv, rhovl, rhovr,rhsu1,rhsu2
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
real*8::othog
real*8::acnx,acny,shp1,shp2,shp3,r,delu
!
eps = 1.d-6
!
!...Give gaussian position and weight...
!
call ruqope_lobatto(1, ngausf, posif, weighf)
!
!...Zero out plnpn, ulnpn
!
dr = 1.d0
ds = 1.d0
rc = 0.0
sc = rc
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
!...The gauss quadrature points in one square....
!
if(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = posif(2); posiq(2, 9) =-1.d0;
posiq(1,10) = posif(3); posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = posif(2);
posiq(1,12) = 1.d0; posiq(2,12) = posif(3);

posiq(1,13) = posif(3); posiq(2,13) = 1.d0;
posiq(1,14) = posif(2); posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = posif(3);
posiq(1,16) =-1.d0; posiq(2,16) = posif(2);
endif
!
!...Part II: Get the gauss point velocity...
!
fgausl = 0.d0
fgausr = 0.d0
!
do 450 ifa = 1, nafac !...(1)ie = 1,nelem
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)

xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))

!...For the linear PP+

ftx = xpf(1 ,2)- xpf(1, 1)
fty = xpf(2, 2)- xpf(2, 1)

fnx = -fty/sqrt(ftx**2 + fty**2)
fny =  ftx/sqrt(ftx**2 + fty**2)
!
!...For boundary faces...
!
if(ifa.le.nbfac)then

!
rc= geoel(1, iel) !...mass center...
sc= geoel(2, iel)

!...Parameters for Left cell

rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...Identify the fgaus for left cell
!
call getpsfgl(ipqua(1:nvqua, iel), ipf, fgausl)
!
unkng = 0.d0

!...Basis fctn...

do ig =3 ,ngausf
b(1, ig) = 1.d0
b(2, ig) = (posiq(1, fgausl(ig))-rc)/dr
b(3, ig) = (posiq(2, fgausl(ig))-sc)/ds
enddo
!
!if(iel.eq.1) print*,'unng',ig,unkno(1:3,2:3,iel),b(1:3, 3)

!...Unkno at gauss points...

do ig = 3, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*b(ideg, ig)
enddo
vngs(1:3, ig, 1) = gelagq(1:3, fgausl(ig), iel) !...face normal vector
enddo

!...Limiter...

do ig = 3, ngausf
!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
!pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
pvtx = (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2))
!
!if(iel.eq.9)print*,'mptv',ig,posiq(:, fgausl(ig)),b(1:3, ig),pvtx,rhovt,unkno(1:3,1,iel)
!
if(nlimi.eq.6)then

rhomc = 1.d0/rhoct

rhomv = rhomc + aflim(1, iel)*(unkng(1, ig, 1) - rhomc)
unkng(1, ig, 1) = rhomv
rhovt = 1.d0/rhomv

dudr = afvec(1, 1, iel)*unkno(2,2,iel) +  afvec(1, 2, iel)*unkno(2,3,iel)
duds = afvec(1, 1, iel)*unkno(3,2,iel) +  afvec(1, 2, iel)*unkno(3,3,iel)
dvdr = afvec(2, 1, iel)*unkno(2,2,iel) +  afvec(2, 2, iel)*unkno(2,3,iel)
dvds = afvec(2, 1, iel)*unkno(3,2,iel) +  afvec(2, 2, iel)*unkno(3,3,iel)
!
uvtx = unkno(1,2,iel)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,iel)  + dvdr*b(2, ig) + dvds*b(3, ig)
!
pvtx = pctr + aflim(4, iel)*(pvtx - pctr)

!...updtae unknv(2:3,:)

unkng(2, ig, 1) = uvtx
unkng(3 ,ig, 1) = vvtx
!
endif
!
unkng(5, ig, 1) = pvtx
unkng(6, ig, 1) = rhoct*sdctr
unkng(7, ig, 1) = sdctr
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig, 1) = -pvtx
sigmg(1, 2, ig, 1) = 0.d0
sigmg(2, 1, ig, 1) = 0.d0
sigmg(2, 2, ig, 1) = -pvtx!
!
enddo
!
!...Get the averaged velocity at the gauss quadrature points...
!
do ig = 3, ngausf
!!
vlave(1, fgausl(ig), iel) = unkng(2, ig, 1)
vlave(2, fgausl(ig), iel) = unkng(3, ig, 1)
enddo
!
!...Impedence...
!
do ig = 3, ngausf
do ic = 1, 1
!
!aujmp(1:2, ig, ic) = vlave(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgaus(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)
!
!...Impedence
!
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) !+ unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo
!
!...Summation over corners
!
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = 3, ngausf
do ic = 1, 1
!
!call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
!if(iel.eq.1892) print*,'mptl', iel,fgausl(ig), murie(ic), vngs(1:3, ig, ic),aujmp(1:3, ig, ic)
!if(ier.eq.1892) print*,'mptr', ic, ig,ier,fgausr(ig), murie(ic), vngs(1:3, ig, ic),aujmp(1:3, ig, ic)
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
!
enddo
!
if(indnd(ipf(1))*indnd(ipf(2)).eq.0)then
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig) !- fpres(1, ipoin)
rhsu2 = munacu(2, ig) - snsigm(2, ig) !- fpres(2, ipoin)
!
!if(iel.eq.1)print*,'iel10',iel,ig,ic,fgausl(ig),unkng(2:3, ig, 1)
!
ufgaus(1, fgausl(ig), iel) = unkng(2, ig, 1)!1.d10!munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ufgaus(2, fgausl(ig), iel) = unkng(3, ig, 1)!1.d10!munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!
endif
!
!if(iel.le.10)print*,'iel1',iel,ufgaus(2, fgausl(3:4), iel)
enddo
!
!...For interior faces...
!
elseif(ifa.gt.nbfac)then
!
rc= geoel(1, iel) !...mass center...
sc= geoel(2, iel)
!
!...Parameters for Left cell
!
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...Identify the fgaus for left cell
!
call getpsfgl(ipqua(1:nvqua, iel), ipf, fgausl)
!
unkng = 0.d0
!
!...Basis fctn...
!
do ig =3 ,ngausf
b(1, ig) = 1.d0
b(2, ig) = (posiq(1, fgausl(ig))-rc)/dr
b(3, ig) = (posiq(2, fgausl(ig))-sc)/ds
enddo
!
!...Unkno at gauss points...
!
do ig = 3, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*b(ideg, ig)
enddo
vngs(1:3, ig, 1) = gelagq(1:3, fgausl(ig), iel) !...face normal vector
enddo
!
!...Limiter...
!
do ig = 3, ngausf
!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
!pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
pvtx = (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, iel)*(unkng(1, ig, 1) - rhomc)
unkng(1, ig, 1) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, iel)*unkno(2,2,iel) +  afvec(1, 2, iel)*unkno(2,3,iel)
duds = afvec(1, 1, iel)*unkno(3,2,iel) +  afvec(1, 2, iel)*unkno(3,3,iel)
dvdr = afvec(2, 1, iel)*unkno(2,2,iel) +  afvec(2, 2, iel)*unkno(2,3,iel)
dvds = afvec(2, 1, iel)*unkno(3,2,iel) +  afvec(2, 2, iel)*unkno(3,3,iel)
!
uvtx = unkno(1,2,iel)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,iel)  + dvdr*b(2, ig) + dvds*b(3, ig)
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pvtx = pctr + aflim(4, iel)*(pvtx - pctr)
!!
!...updtae unknv(2:3,:)
unkng(2, ig, 1) = uvtx
unkng(3 ,ig, 1) = vvtx
!
endif
!
unkng(5, ig, 1) = pvtx
unkng(6, ig, 1) = rhoct*sdctr
unkng(7, ig, 1) = sdctr
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig, 1) = -pvtx
sigmg(1, 2, ig, 1) = 0.d0
sigmg(2, 1, ig, 1) = 0.d0
sigmg(2, 2, ig, 1) = -pvtx!
!
enddo
!
!...Right cell
!
!
rc= geoel(1, ier) !...mass center...
sc= geoel(2, ier)
!
rhoct = 1.d0/unkno(1, 1, ier)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ier)
vctr  = unkno(1, 3, ier)
ectr  = unkno(1, 4, ier)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center..
!
!...Identify the fgaus for right cell
!
call getpsfgr(ipqua(1:nvqua, ier), ipf, fgausr)
!
!if(iel.eq.2433) print*,'fgaussr',ipqua(1:nvqua, iel),ipf(1:3),fgausl
!
do ig =3 ,ngausf
b(1, ig) = 1.d0
b(2, ig) = (posiq(1, fgausr(ig))-rc)/dr
b(3, ig) = (posiq(2, fgausr(ig))-sc)/ds
enddo
!
do ig = 3, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 2) = unkng(1:nq, ig, 2) + unkno(ideg,1:nq,ier)*b(ideg, ig)
enddo
vngs(1:3, ig, 2) = gelagq(1:3, fgausr(ig), ier) !...face normal vector
enddo
!
!...Limiter...
!
do ig = 3, ngausf
!
rhovt = 1.d0/unkng(1, ig, 2)
uvtx = unkng(2, ig, 2)
vvtx = unkng(3, ig, 2)
evtx = unkng(4, ig, 2)
!
!pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
pvtx =(gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, ier)*(unkng(1, ig, 2) - rhomc)
unkng(1, ig, 2) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, ier)*unkno(2,2,ier) +  afvec(1, 2, ier)*unkno(2,3,ier)
duds = afvec(1, 1, ier)*unkno(3,2,ier) +  afvec(1, 2, ier)*unkno(3,3,ier)
dvdr = afvec(2, 1, ier)*unkno(2,2,ier) +  afvec(2, 2, ier)*unkno(2,3,ier)
dvds = afvec(2, 1, ier)*unkno(3,2,ier) +  afvec(2, 2, ier)*unkno(3,3,ier)
!
uvtx = unkno(1,2,ier)  + dudr*b(2, ig) + duds*b(3, ig)
vvtx = unkno(1,3,ier)  + dvdr*b(2, ig) + dvds*b(3, ig)
!
!if(iel.eq.2581) print*,'ie26252628',ier,pctr,aflim(4, ier),pvtx
!
pvtx = pctr + aflim(4, ier)*(pvtx - pctr)
!!
!...updtae unknv(2:3,:)
unkng(2, ig, 2) = uvtx
unkng(3, ig, 2) = vvtx
!
endif
!
unkng(5, ig, 2) = pvtx
unkng(6, ig, 2) = rhoct*sdctr
unkng(7, ig, 2) = sdctr
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig, 2) = -pvtx
sigmg(1, 2, ig, 2) = 0.d0
sigmg(2, 1, ig, 2) = 0.d0
sigmg(2, 2, ig, 2) = -pvtx!
!
enddo
!
!...Get the averaged velocity at the gauss quadrature points...
!
do ig = 3, ngausf
!!
vlave(1, fgausl(ig), iel) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
vlave(2, fgausl(ig), iel) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
!!
vlave(1, fgausr(ig), ier) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
vlave(2, fgausr(ig), ier) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
enddo
!
!...Impedence...
!
do ig = 3, ngausf
do ic = 1, 2
!
!aujmp(1:2, ig, ic) = vlave(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgaus(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)
!
!...Impedence
!
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic)! + unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo
!
!...Summation over corners
!
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = 3, ngausf
do ic = 1, 2
!
!call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)

munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
!if(iel.eq.1892) print*,'mptl', iel,fgausl(ig), murie(ic), vngs(1:3, ig, ic),aujmp(1:3, ig, ic)
!if(ier.eq.1892) print*,'mptr', ic, ig,ier,fgausr(ig), murie(ic), vngs(1:3, ig, ic),aujmp(1:3, ig, ic)
!
if(ic.eq.1)then
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
!
!if(iel.eq.1892) print*,'mptr', ic, ig,iel, fgausl(ig),snsigmlq(1:2, fgausl(ig), iel)
!
else
!
munaclq(1:2, 1, fgausr(ig), ier) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausr(ig), ier) =  munacn_rie(1:2, 2)
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaulq(1:2, fgausr(ig), ier)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausr(ig), ier)    = snsigm_rie(1:2)
!
!if(ier.eq.1892) print*,'mptr', ic, ig,ier, fgausr(ig),snsigmlq(1:2,  fgausr(ig), ier)
endif
!
enddo
!
enddo
!
!...The normal velocity at the gauss points...
!
do ig = 3, ngausf
!
rhovl = 1.d0/unkng(1, ig, 1)
uvtxl = unkng(2, ig, 1)
vvtxl = unkng(3, ig, 1)
evtxl = unkng(4, ig, 1)
pvtxl = unkng(5, ig, 1)
mufal = unkng(6, ig, 1)
!
rhovr = 1.d0/unkng(1, ig, 2)
uvtxr = unkng(2, ig, 2)
vvtxr = unkng(3, ig, 2)
evtxr = unkng(4, ig, 2)
pvtxr = unkng(5, ig, 2)
mufar = unkng(6, ig, 2)
!
fnx = vngs(1, ig, 1) !...face normal vector for iel...
fny = vngs(2, ig, 1)
!
ftx = -fny
fty =  fnx
!
!...Mar
!
!if(indnd(ipf(1)).eq.0)then
ufgaus(1, fgausl(ig), iel) = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (pvtxr- pvtxl)/(mufal+mufar)*fnx
ufgaus(2, fgausl(ig), iel) = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (pvtxr- pvtxl)/(mufal+mufar)*fny
!
ufgaus(1, fgausr(ig), ier) = ufgaus(1, fgausl(ig), iel)
ufgaus(2, fgausr(ig), ier) = ufgaus(2, fgausl(ig), iel)
!
!if(iel.eq.1)print*,'mptvelo',ig,ifa,ufgaus(1, fgausl(ig), iel)
!
!endif
!
enddo
endif
450 enddo
!
end subroutine getvelo_mpt_gausseulframe
!
!...subroutine: Calculate the F^* N ds for all face gauss points for hybrid grids in Eulerian framework...
!
subroutine getfnds_lag_gausseulframe(gflag,gelag,gelagq,intfac,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(inout)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad),      intent(inout)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord

!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv,ishp

!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf

!...local real array
real*8,dimension(1:2, 1:2)    ::comatr !...cofactor matrix
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8,dimension(1:ndimn, 1:nptri)::xpt
real*8,dimension(1:ndimn, 1:npqua)::xpq
real*8::vnorm(1:2)
real*8,dimension(1:nvtri):: xv, yv
real*8::posit(2, ngelg)
real*8::posiq(2, ngelgq)
real*8::posif(ngausf),weighf(ngausf)
real*8, dimension(1:nptri):: shp,dspr,dsps
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
real*8, dimension(2):: rgaus
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
!...Part I: Specify some gauss points
!

!...Give gaussian position and weight...
call ruqope_lobatto(1, ngausf, posif, weighf)

!
if(ngausf.eq.2)then
!...Tria
posit(1, 1) = 0.d0; posit(2, 1) = 0.d0;
posit(1, 2) = 1.d0; posit(2, 2) = 0.d0;
posit(1, 3) = 1.d0; posit(2, 3) = 0.d0;
posit(1, 4) = 0.d0; posit(2, 4) = 1.d0;
posit(1, 5) = 0.d0; posit(2, 5) = 1.d0;
posit(1, 6) = 0.d0; posit(2, 6) = 0.d0;

!...Quad
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
!
elseif(ngausf.eq.4)then
!...Tria
posit(1, 1) = 0.d0; posit(2, 1) = 0.d0;
posit(1, 2) = 1.d0; posit(2, 2) = 0.d0;
posit(1, 3) = 1.d0; posit(2, 3) = 0.d0;
posit(1, 4) = 0.d0; posit(2, 4) = 1.d0;
posit(1, 5) = 0.d0; posit(2, 5) = 1.d0;
posit(1, 6) = 0.d0; posit(2, 6) = 0.d0;

posit(1, 7) = 0.5d0-sqrt(3.d0)/6.d0; posit(2, 7) = 0.d0;
posit(1, 8) = 0.5d0+sqrt(3.d0)/6.d0; posit(2, 8) = 0.d0;

posit(1, 9) = 0.5d0+sqrt(3.d0)/6.d0; posit(2, 9) = 0.5d0-sqrt(3.d0)/6.d0;
posit(1,10) = 0.5d0-sqrt(3.d0)/6.d0; posit(2,10) = 0.5d0+sqrt(3.d0)/6.d0;

posit(1,11) = 0.d0; posit(2,11) = 0.5d0+sqrt(3.d0)/6.d0;
posit(1,12) = 0.d0; posit(2,12) = 0.5d0-sqrt(3.d0)/6.d0;

!...Quad
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = posif(2); posiq(2, 9) =-1.d0;
posiq(1,10) = posif(3); posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = posif(2);
posiq(1,12) = 1.d0; posiq(2,12) = posif(3);

posiq(1,13) = posif(3); posiq(2,13) = 1.d0;
posiq(1,14) = posif(2); posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = posif(3);
posiq(1,16) =-1.d0; posiq(2,16) = posif(2);

endif
!
!...Part II: Boundary face...
!
do ifa =1 , nbfac
iel = intfac(1, ifa)

ipf(1:nvfac) = intfac(3:2+nvfac, ifa)

xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))

anx =   xpf(2, 2) - xpf(2, 1)
any = -(xpf(1, 2) - xpf(1, 1))

gflag(1, ifa) = anx/sqrt(anx**2 + any**2)
gflag(2, ifa) = any/sqrt(anx**2 + any**2)
gflag(3, ifa) = sqrt(anx**2 + any**2)
enddo
!
!...Part III: Get the normal for all the elements
!
!...III.1: Triangle
do 100 ie=1, ntria!...(1)ie=1, ntria

ipt(1:nvtri) = iptri(1:nvtri,ie)

!...coordinates
if(ncurv==0)then
xpt(1, 1:3) = coord(1, iptri(1:3,ie))
xpt(2, 1:3) = coord(2, iptri(1:3,ie))

xpt(1:2,4) = 0.5d0*(xpt(1:2,1) + xpt(1:2,2))
xpt(1:2,5) = 0.5d0*(xpt(1:2,2) + xpt(1:2,3))
xpt(1:2,6) = 0.5d0*(xpt(1:2,1) + xpt(1:2,3))
elseif(ncurv==1)then
xpt(1, 1:nptri) = coord(1,iptri(1:nptri, ie))
xpt(2, 1:nptri) = coord(2,iptri(1:nptri, ie))
endif
!
do ig = 1, ngelg!...(2)ig = 1,ngelg
r  = posit(1,ig)
s  = posit(2,ig)

!...Shape functions and their derivatives...
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)

dspr(1) = c10-4.d0*(c10-r-s)
dspr(2) = -1.d0 + 4.d0*r
dspr(3) = 0.d0
dspr(4) = 4.d0*(1.d0-2.d0*r-s)
dspr(5) = 4.d0*s
dspr(6) =-4.d0*s

dsps(1) = c10-4.d0*(c10-r-s)
dsps(2) = 0.d0
dsps(3) = -1.d0 + 4.d0*s
dsps(4) = -4.d0*r
dsps(5) =  4.d0*r
dsps(6) =  4.d0*(1.d0-r-2.d0*s)

dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0

do ishp = 1, nptri
dxdr = dxdr + dspr(ishp)*xpt(1,ishp)
dxds = dxds + dsps(ishp)*xpt(1,ishp)

dydr = dydr + dspr(ishp)*xpt(2,ishp)
dyds = dyds + dsps(ishp)*xpt(2,ishp)
enddo

!...Cofactor matrix for left cell
comatr(1, 1) = dyds !...yc-ya
comatr(1, 2) =-dydr !...-(yb-ya)
comatr(2, 1) =-dxds !...-(xc-xa)
comatr(2, 2) = dxdr !...xb-xa

!...Identify the local No. of one internal face for left cell...
if(ig.eq.1.or.ig.eq.2.or.ig.eq.7.or.ig.eq.8)then
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 1.d0;

idfal = ig;

elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.9.or.ig.eq.10)then
vnorm(1) = sqrt(2.d0)*0.5d0; vnorm(2) = sqrt(2.d0)*0.5d0;
larea    = sqrt(2.d0)

idfal = ig;

elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.11.or.ig.eq.12)then
vnorm(1) = -1.d0;           vnorm(2) = 0.d0;
larea    = 1.d0

idfal = ig;

endif

anx = comatr(1, 1)*vnorm(1) + comatr(1, 2)*vnorm(2)
any = comatr(2, 1)*vnorm(1) + comatr(2, 2)*vnorm(2)
vnorm(1) = anx*larea
vnorm(2) = any*larea

!...Unit vector...
farea    = sqrt(vnorm(1)**2 + vnorm(2)**2) !...farea: face area...

gelag(1, idfal, ie) = anx*larea/farea
gelag(2, idfal, ie) = any*larea/farea
gelag(3, idfal, ie) = farea
!
enddo !...ig from 1 to ngelg
!
100 enddo  !...(1)ie=1,ntria

!...III.2:Quad
do 200 ie=1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...coordinates
if(ncurv==0)then
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))
elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif
!
do ig = 1, ngelgq!...(2)ig = 1,ngelgq
r  = posiq(1,ig)
s  = posiq(2,ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
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
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 2.d0;

idfal = ig;

elseif(ig.eq.3.or.ig.eq.4.or.ig.eq.11.or.ig.eq.12)then
vnorm(1) = 1.d0; vnorm(2) = 0.d0;
larea    = 2.d0

idfal = ig;

elseif(ig.eq.5.or.ig.eq.6.or.ig.eq.13.or.ig.eq.14)then
vnorm(1) = 0.d0;           vnorm(2) = 1.d0;
larea    = 2.d0
idfal = ig;

elseif(ig.eq.7.or.ig.eq.8.or.ig.eq.15.or.ig.eq.16)then
vnorm(1) =-1.d0;           vnorm(2) = 0.d0;
larea    = 2.d0

idfal = ig;

endif
!
anx = comatr(1, 1)*vnorm(1) + comatr(1, 2)*vnorm(2)
any = comatr(2, 1)*vnorm(1) + comatr(2, 2)*vnorm(2)
vnorm(1) = anx*larea
vnorm(2) = any*larea

!...Unit vector...
farea    = sqrt(vnorm(1)**2 + vnorm(2)**2) !...farea: face area...
!
gelagq(1, idfal, ie) = anx*larea/farea
gelagq(2, idfal, ie) = any*larea/farea
gelagq(3, idfal, ie) = farea
!
enddo !...ig from 1...12
!
200 enddo  !...(1)ie=1,nquad

end subroutine getfnds_lag_gausseulframe
!
!...rhsiface for Lagrangin motion in Eulerian frmaework
!
subroutine rhsifacedg_lagtriagseulframe(iptri,  unkno, ustar, ufgpt, geoel, gelag, fstart, coord, coold,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar, coord, coold
real*8,dimension(1:ndimn,1:ngelg, 1:ntria), intent(in)  ::ufgpt
real*8,dimension(1:ndimn,1:ngelg, 1:ntria),  intent(in) ::fstart !...Riemann forces


real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::rhsel
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelg, 1:ntria+nbfac),    intent(in)::gelag

!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ielem

!...local integer array
integer,dimension(1:nvtri) :: ipt
integer, dimension(nvfac, 3)::fglvt
integer,dimension(ngausf, 3)::fglgt
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::positg(1:2, 1:ngelg)
real*8::xvt(nvtri), yvt(nvtri)
real*8::bg(3, ngausf)
!
real*8::weigh(ngausf)
real*8::weighf(ngausf), posif(ngausf)
!...local real number
real*8::eps,c00,c05,c10,c20,c13
real*8::wi
real*8:: nx, ny ,sa
real*8:: rg, sg
real*8::dr,ds,rc,sc
real*8::dwav1,dwav2
real*8::anx, any
real*8::gpnx, gpny, gpsa
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
!
!...Part I: Specify some gauss points
!

!...Give gaussian position and weight...
call ruqope_lobatto(1, ngausf, posif, weighf)

!...Local vertex No. of gauss points in one unit ...
fglvt(1, 1) = 1;  fglvt(2, 1) = 2; !fglvt(3, 1) = 4;
fglvt(1, 2) = 2;  fglvt(2, 2) = 3; !fglvt(3, 2) = 5;
fglvt(1, 3) = 3;  fglvt(2, 3) = 1; !fglvt(3, 3) = 6;

!
if(ngausf.eq.2)then
weigh(  1) = weighf(1)/2.d0
weigh(  2) = weighf(2)/2.d0

!...Local gauss point No. of any gauss point in one face...
fglgt(1, 1) = 1;  fglgt(2, 1) = 2;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6;

!...6-nodes tria
positg(1, 1) = 0.d0; positg(2, 1) = 0.d0;
positg(1, 2) = 1.d0; positg(2, 2) = 0.d0;
positg(1, 3) = 1.d0; positg(2, 3) = 0.d0;
positg(1, 4) = 0.d0; positg(2, 4) = 1.d0;
positg(1, 5) = 0.d0; positg(2, 5) = 1.d0;
positg(1, 6) = 0.d0; positg(2, 6) = 0.d0;

elseif(ngausf.eq.4)then
weigh(  1) = weighf(1)/2.d0
weigh(  2) = weighf(4)/2.d0
weigh(  3) = weighf(2)/2.d0
weigh(  4) = weighf(3)/2.d0

!...Local gauss point No. of any gauss point in one face...
fglgt(1, 1) = 1;  fglgt(2, 1) = 2; fglgt(3, 1) = 7; fglgt(4, 1) =  8;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4; fglgt(3, 2) = 9; fglgt(4, 2) = 10;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6; fglgt(3, 3) =11; fglgt(4, 3) = 12;

!
positg(1, 1) = 0.d0; positg(2, 1) = 0.d0;
positg(1, 2) = 1.d0; positg(2, 2) = 0.d0;
positg(1, 3) = 1.d0; positg(2, 3) = 0.d0;
positg(1, 4) = 0.d0; positg(2, 4) = 1.d0;
positg(1, 5) = 0.d0; positg(2, 5) = 1.d0;
positg(1, 6) = 0.d0; positg(2, 6) = 0.d0;

positg(1, 7) = 0.5d0-sqrt(3.d0)/6.d0; positg(2, 7) = 0.d0;
positg(1, 8) = 0.5d0+sqrt(3.d0)/6.d0; positg(2, 8) = 0.d0;

positg(1, 9) = 0.5d0+sqrt(3.d0)/6.d0; positg(2, 9) = 0.5d0-sqrt(3.d0)/6.d0;
positg(1,10) = 0.5d0-sqrt(3.d0)/6.d0; positg(2,10) = 0.5d0+sqrt(3.d0)/6.d0;

positg(1,11) = 0.d0; positg(2,11) = 0.5d0+sqrt(3.d0)/6.d0;
positg(1,12) = 0.d0; positg(2,12) = 0.5d0-sqrt(3.d0)/6.d0;
endif

!
dr = .5d0
ds = .5d0
!
xvt(1) =  0.d0; yvt(1) = 0.d0
xvt(2) =  1.d0; yvt(2) = 0.d0
xvt(3) =  0.d0; yvt(3) = 1.d0
!
!...Part II: Loop over the elements
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
ielem = ie

!...Zero out ulnpn, plnpn, elnpn
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
ipt(1:nvtri) = iptri(1:nvtri ,ie)
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Loop over the 3 faces
do ifa = 1, 3
do ig   =  1,  ngausf
wi  = weigh(ig)
!
rg = positg(1,fglgt(ig, ifa))
sg = positg(2,fglgt(ig, ifa))
!
bg(1,  ig) = 1.d0
bg(2,  ig) = (rg-rc)/dr
bg(3,  ig) = (sg-sc)/ds

!...The vertex constituting one cell...
gpnx = gelag(1, fglgt(ig, ifa),ie)
gpny = gelag(2, fglgt(ig, ifa),ie)
gpsa = gelag(3, fglgt(ig, ifa),ie)

!...Distribute to every corner...
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ufgpt(1, fglgt(ig, ifa),ie)*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig) +&
ufgpt(2, fglgt(ig, ifa),ie)*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstart(1, fglgt(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstart(2, fglgt(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
elnpn(1:ndegr)=elnpn(1:ndegr) +&
ufgpt(1, fglgt(ig, ifa), ie)*&
fstart(1, fglgt(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig) +&
ufgpt(2, fglgt(ig, ifa), ie)*&
fstart(2, fglgt(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
enddo
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1829) print*,'rhs iface idegr 1',ie,rhsel(1:3, 2, ielem)
!
550 enddo
!
! open(8,file='rhstxt.dat')
!   do ie = 1, nquad
!      write(8,*) ie, rhsel(1, 1:4, ie)
!   enddo
! close(8)
end subroutine rhsifacedg_lagtriagseulframe
!
!...rhsiface for Lagrangin motion in Eulerian frmaework
!
subroutine rhsifacedg_lagquadgausseulframe(ipqua,  unkno, ustar, ufgaus, geoel, gelagq, fstarq, coord, coold,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar, coord, coold!...nodal velocity
real*8,dimension(1:ndimn,1:ngelgq, 1:nquad), intent(in)  ::ufgaus
real*8,dimension(1:ndimn,1:ngelgq,1:nquad),  intent(in) ::fstarq !...Riemann forces


real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::rhsel
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ielem
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer, dimension(nvfac, 4)::fglvq
integer,dimension(ngausf, 4)::fglgq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::posiqg(1:2, 1:ngelgq)
real*8::xvq(nvqua), yvq(nvqua)
real*8::bg(3, ngausf)
!
real*8::weigh(ngausf)
real*8::weighf(ngausf), posif(ngausf)
!...local real number
real*8::eps,c00,c05,c10,c20,c13
real*8::wi
real*8:: nx, ny ,sa
real*8:: rg, sg
real*8::dr,ds,rc,sc
real*8::dwav1,dwav2
real*8::anx, any
real*8::gpnx, gpny, gpsa
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
!
!...Give gaussian position and weight...
!
call ruqope_lobatto(1, ngausf, posif, weighf)
!
!...Local vertex No. of gauss points in one unit ...
!
fglvq(1, 1) = 1;  fglvq(2, 1) = 2;! fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3;! fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4;! fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1;! fglvq(3, 4) = 8;
!
if(ngausf.eq.2)then
!
weigh(  1) = weighf(1)/2.d0
weigh(  2) = weighf(2)/2.d0
!
!...Local gauss point No. of any gauss point in one face...
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8;
!
posiqg(1, 1) =-1.d0; posiqg(2, 1) =-1.d0;
posiqg(1, 2) = 1.d0; posiqg(2, 2) =-1.d0;
posiqg(1, 3) = 1.d0; posiqg(2, 3) =-1.d0;
posiqg(1, 4) = 1.d0; posiqg(2, 4) = 1.d0;
posiqg(1, 5) = 1.d0; posiqg(2, 5) = 1.d0;
posiqg(1, 6) =-1.d0; posiqg(2, 6) = 1.d0;
posiqg(1, 7) =-1.d0; posiqg(2, 7) = 1.d0;
posiqg(1, 8) =-1.d0; posiqg(2, 8) =-1.d0;
!
elseif(ngausf.eq.4)then
!
weigh(  1) = weighf(1)/2.d0
weigh(  2) = weighf(4)/2.d0
weigh(  3) = weighf(2)/2.d0
weigh(  4) = weighf(3)/2.d0
!
!...Local gauss point No. of any gauss point in one face...
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
!
posiqg(1, 1) =-1.d0; posiqg(2, 1) =-1.d0;
posiqg(1, 2) = 1.d0; posiqg(2, 2) =-1.d0;
posiqg(1, 3) = 1.d0; posiqg(2, 3) =-1.d0;
posiqg(1, 4) = 1.d0; posiqg(2, 4) = 1.d0;
posiqg(1, 5) = 1.d0; posiqg(2, 5) = 1.d0;
posiqg(1, 6) =-1.d0; posiqg(2, 6) = 1.d0;
posiqg(1, 7) =-1.d0; posiqg(2, 7) = 1.d0;
posiqg(1, 8) =-1.d0; posiqg(2, 8) =-1.d0;

posiqg(1, 9) = posif(2); posiqg(2, 9) =-1.d0;
posiqg(1,10) = posif(3); posiqg(2,10) =-1.d0;

posiqg(1,11) = 1.d0; posiqg(2,11) = posif(2);
posiqg(1,12) = 1.d0; posiqg(2,12) = posif(3);

posiqg(1,13) = posif(3); posiqg(2,13) = 1.d0;
posiqg(1,14) = posif(2); posiqg(2,14) = 1.d0;

posiqg(1,15) =-1.d0; posiqg(2,15) = posif(3);
posiqg(1,16) =-1.d0; posiqg(2,16) = posif(2);
endif
!
!...Zero out plnpn, ulnpn
!
dr = 1.d0
ds = 1.d0
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
do 550 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
ipq(1:nvqua) = ipqua(1:nvqua ,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do ifa = 1, 4
!
do ig   =  1,  ngausf
!
wi  = weigh(ig)
!
rg = posiqg(1,fglgq(ig, ifa))
sg = posiqg(2,fglgq(ig, ifa))
!
bg(1,  ig) = 1.d0
bg(2,  ig) = (rg-rc)/dr
bg(3,  ig) = (sg-sc)/ds
!
!...The vertex constituting one cell...
!
gpnx = gelagq(1, fglgq(ig, ifa),ie)
gpny = gelagq(2, fglgq(ig, ifa),ie)
gpsa = gelagq(3, fglgq(ig, ifa),ie)
!
!...Distribute to every corner...
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ufgaus(1, fglgq(ig, ifa),ie)*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig) +&
ufgaus(2, fglgq(ig, ifa),ie)*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)
!
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
elnpn(1:ndegr)=elnpn(1:ndegr) +&
ufgaus(1, fglgq(ig, ifa), ie)*&
fstarq(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig) +&
ufgaus(2, fglgq(ig, ifa), ie)*&
fstarq(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
!if(ie.eq.9) print*,'rhs iface idegr',ie,ifa,ig,fglgq(ig, ifa),&
!fstarq(2, fglgq(ig, ifa), ie),bg(1:ndegr, ig),weigh(ig)
!if(ie==1831) print*,'rhs ifaceidegr2',ifa,ig,fglgq(ig, ifa),ulnpn(1:ndegr),ustar(1, ipq(fglvq(1, 2))),ipq(fglvq(1, 2))
!
enddo
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1829) print*,'rhs iface idegr 1',ie,rhsel(1:3, 2, ielem)
!
550 enddo
!
! open(8,file='rhstxt.dat')
!   do ie = 1, nquad
!      write(8,*) ie, rhsel(1, 1:4, ie)
!   enddo
! close(8)
!
end subroutine rhsifacedg_lagquadgausseulframe
!
!...Get the velocity at the vertex...
!
subroutine getvelo_quadvertex(ipqua, indnd, geoel,  ustar, coord,ufgaus, unkno)
use constant
implicit none
integer,dimension(1:nvqua,1:nquad),         intent(in):: ipqua
integer,dimension(1:npoin),         intent(in):: indnd
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord !...nodal velocity

real*8,  dimension(1:nquad),                 intent(in)::ufgaus(1:ndimn,1:16,1:nquad)
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
!
real*8,dimension(1:ndimn, 1:4) :: xpq
real*8,dimension(1:ndimn, 1:16) :: xpqg
real*8::posiq(1:2, 1:ngelgq)
integer, dimension(ngausf, 4):: fglgq

real*8,dimension(1:3, 1:8)::bq
real*8,dimension(1:2, 1:2)::amnul
real*8,dimension(1:2, 1:2)::dunul
real*8,dimension(1:ndimn, 1:9)::upqua
real*8::ugrad(2, 2)
integer,dimension(1:nvqua) :: ipq
!
integer::iv, ie, ielem, ij, ik,ig,ipoin,ishp
real*8::amnua, amnub, amnuc, amnud, amnuj
real*8::dxm,dym,xcv,ycv
real*8::deluq(1:ndimn, 1:8)
real*8::x5(2,2)
real*8::vevtx(2,npoin)
real*8::muvtx(npoin)
!
real*8::rhomc,rhoct,uctr,vctr,pctr,murie,ectr
real*8::rm,rp,sm,sp,r,s,rc,sc
real*8::shpq(4)
real*8::eps,c10
!
eps = 1.d-6
c10 =1.d0
!
!...Local gauss point No. of any gauss point in one face...
!
if(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.5773502691896257645091488d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.5773502691896257645091488d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.5773502691896257645091488d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.5773502691896257645091488d0;

posiq(1,13) = 0.5773502691896257645091488d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.5773502691896257645091488d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.5773502691896257645091488d0; !?
posiq(1,16) =-1.d0; posiq(2,16) =-0.5773502691896257645091488d0;

endif
!
!...Intialize array
!
vevtx = 0.d0
muvtx = 0.d0
!
do 450 ie = 1,nquad !...(1)ie = 1,nelem
!
amnul = 0.d0
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpqg = 0.d0

do ig = 1, 16
r  = posiq(1,ig)
s  = posiq(2,ig)
!
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s

!...  shape function & its derivatives w.r.t. reference coordinates

shpq(1) = 0.25d0*rm*sm
shpq(2) = 0.25d0*rp*sm
shpq(3) = 0.25d0*rp*sp
shpq(4) = 0.25d0*rm*sp

do ishp = 1, 4
xpqg(1, ig) = xpqg(1, ig) + shpq(ishp)*xpq(1,ishp)
xpqg(2, ig) = xpqg(2, ig) + shpq(ishp)*xpq(2,ishp)
enddo
enddo
!
dxm = maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dym = maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...II: Get the cell center's coordinate

rc = geoel(1, ielem);
sc =  geoel(2, ielem)
r = rc; s= sc
xcv = 0.d0; ycv = 0.d0
!
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*rm*sm
shpq(2) = 0.25d0*rp*sm
shpq(3) = 0.25d0*rp*sp
shpq(4) = 0.25d0*rm*sp
!
do ishp = 1, 4
xcv = xcv + shpq(ishp)*xpq(1,ishp)
ycv = ycv + shpq(ishp)*xpq(2,ishp)
enddo
!
do iv =1 ,8
bq(1, iv) = 1.d0
bq(2, iv) = (xpqg(1, iv+8)-xcv)/dxm
bq(3, iv) = (xpqg(2, iv+8)-ycv)/dym
enddo
!
do iv=1, 8
!
amnul(1, 1 ) = amnul(1, 1 ) + bq(2, iv)*bq(2, iv)
amnul(1, 2 ) = amnul(1, 2 ) + bq(2, iv)*bq(3, iv)

amnul(2, 1 ) = amnul(2, 1 ) + bq(3, iv)*bq(2, iv)
amnul(2, 2 ) = amnul(2, 2 ) + bq(3, iv)*bq(3, iv)
enddo
!
amnua = amnul(1, 1 )
amnub = amnul(1, 2 )
amnuc = amnul(2, 1 )
amnud = amnul(2, 2 )
!
!if(ie.eq.1831)print*,'amnua',amnul(1:2, 1:2 )
!
amnuj = amnua*amnud-amnuc*amnub
!
x5(1, 1) = amnud/amnuj
x5(1, 2) =-amnub/amnuj
x5(2, 1) = x5(1, 2)
x5(2, 2) = amnua/amnuj
!
!if(ie.eq.1831)print*,'amnuainverse',x5(1:2, 1:2 )
!
!...III: RHS of Least square
!
do iv=1, 8
upqua(1:2, iv) = ufgaus(1:2, iv+8, ie)
enddo
upqua(1:2, 9) = unkno(1,2:3,ielem)
!
do iv = 1, 8
deluq(1:2, iv) = upqua(1:2, iv) - upqua(1:2, 9)
enddo
!
dunul = 0.d0

do iv=1, 8
dunul(1, 1)  = dunul(1, 1) + bq(2, iv)*deluq(1, iv)
dunul(2, 1)  = dunul(2, 1) + bq(3, iv)*deluq(1, iv)

dunul(1, 2)  = dunul(1, 2) + bq(2, iv)*deluq(2, iv)
dunul(2, 2)  = dunul(2, 2) + bq(3, iv)*deluq(2, iv)
enddo
!
ugrad(1, 1) = x5(1, 1)*dunul(1 ,1) + x5(2, 1)*dunul(2 ,1)
ugrad(2, 1) = x5(1, 2)*dunul(1 ,1) + x5(2, 2)*dunul(2 ,1)

ugrad(1, 2) = x5(1, 1)*dunul(1 ,2) + x5(2, 1)*dunul(2 ,2)
ugrad(2, 2) = x5(1, 2)*dunul(1 ,2) + x5(2, 2)*dunul(2 ,2)
!
do iv =1 ,4
bq(1, iv) = 1.d0
bq(2, iv) = (xpq(1, iv)-xcv)/dxm
bq(3, iv) = (xpq(2, iv)-ycv)/dym
enddo
!
do iv = 1, 4
upqua(1, iv) = upqua(1, 9) + ugrad(1, 1)*bq(2, iv) + ugrad(2, 1)*bq(3, iv)
upqua(2, iv) = upqua(2, 9) + ugrad(1, 2)*bq(2, iv) + ugrad(2, 2)*bq(3, iv)
enddo
!
!...IV: Get the weighted velocity at a vertex
!
rhomc = unkno(1, 1, ielem)
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
murie = rhoct*sqrt(gamlg*pctr/rhoct)
!
vevtx(1, ipq(1:4)) = vevtx(1, ipq(1:4)) + upqua(1, 1:4)*murie
vevtx(2, ipq(1:4)) = vevtx(2, ipq(1:4)) + upqua(2, 1:4)*murie

muvtx(ipq(1:4)) = muvtx(ipq(1:4)) + murie
!
450 enddo
!
!...Get the vertex velocity...
!
do ipoin = 1, npoin
!
if(indnd(ipoin).eq.0)then
ustar(1:2, ipoin) = vevtx(1:2, ipoin)/muvtx(ipoin)
endif
!
enddo

!
end subroutine getvelo_quadvertex
!
!...Calculate the velocity at the Gauss point for novertex face...
!
subroutine getvelo_quadsubg_gs(ufgpq,dufgq,geoel,intfac,ipqua,coord,unkno,indnd, &
munaclq, munaulq, snsigmlq, afvec, aflim)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
integer*4,dimension(1:npoin),                intent(in)::indnd
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngsfq,nquad),      intent(in)::ufgpq
real*8,dimension(1:2,1:ngsfq,nquad),      intent(inout)::dufgq
real*8, dimension(1:2, 1:2, 1:ngsfq, 1:nquad), intent(inout)::munaclq
real*8, dimension(1:ndimn, 1:ngsfq,  1:nquad), intent(inout)::munaulq
real*8, dimension(1:ndimn, 1:ngsfq,  1:nquad), intent(inout)::snsigmlq
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(ngausf)::fgausl, fgausr
!
real*8::eps
real*8,dimension(1:2,1:ngelgq,nquad)::vlave
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::sigmg(1:2, 1:2, 1:ngausf, 1:2)
real*8::vngs(3, ngausf,2)
real*8::vngf(3, ngausf)
real*8::aujmp(1:3,1:ngausf,1:2)
real*8::murie(2)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, ngausf)
real*8::posiq(1:2, 1:ngsfq)
integer,dimension(nvfac, 4)::fglvq
integer,dimension(ngausf,4)::fglgq
real*8::xpf(1:2, 1:nvfac)

!...Riemann parameters...
real*8::munacn_rie(2,2), munacu_rie(2), snsigm_rie(2)
real*8::munacn(1:2, 1:2, 1:ngausf), munacu(1:2, 1:ngausf), snsigm(1:2, 1:ngausf)
real*8::munaci(2, 2)

!...Local real number
real*8::detma,dudr,duds,dvdr,dvds
real*8::pvtx,rhovt, rhomc, rhomv, rhovl, rhovr,rhsu1,rhsu2
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
real*8::othog
real*8::acnx,acny,shp1,shp2,shp3,r,delu

!...Initial velocity at gauss point
!..Initial Gauss point velocity from the interpolation of the 3-node velocity...
!
eps = 1.d-6
dufgq = 0.d0
!
dr = 1.d0
ds = 1.d0
rc = 0.0
sc = rc
!
!...Part I: Specify some gauss points
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
xvq(9) =  0.d0; yvq(9) =  0.d0

!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

!...8-nodes quad
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.4472135954999579392818d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.4472135954999579392818d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.4472135954999579392818d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.4472135954999579392818d0;

posiq(1,13) = 0.4472135954999579392818d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.4472135954999579392818d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.4472135954999579392818d0;
posiq(1,16) =-1.d0; posiq(2,16) =-0.4472135954999579392818d0;

elseif(ngausf.eq.5)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) =-0.6546536707079771437983d0; posiq(2,13) =-1.d0;
posiq(1,14) = 0.6546536707079771437983d0; posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) =-0.6546536707079771437983d0;
posiq(1,16) = 1.d0; posiq(2,16) = 0.6546536707079771437983d0;

posiq(1,17) = 0.6546536707079771437983d0; posiq(2,17) = 1.d0;
posiq(1,18) =-0.6546536707079771437983d0; posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = 0.6546536707079771437983d0;
posiq(1,20) =-1.d0; posiq(2,20) =-0.6546536707079771437983d0;
endif
!
!...Part II: Get the corrected velocity at the non-vertex gauss point...
!
fgausl = 0.d0
fgausr = 0.d0
!
do 450 ifa = 1, nafac !...(1)ifa = 1, nafac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))

!...For the linear PP+
ftx = xpf(1 ,2)- xpf(1, 1)
fty = xpf(2, 2)- xpf(2, 1)
!
fnx = -fty/sqrt(ftx**2 + fty**2)
fny =  ftx/sqrt(ftx**2 + fty**2)

!...For linear PM
ftx = xpf(1 ,3)- xpf(1, 1)
fty = xpf(2, 3)- xpf(2, 1)
!
lenmc = sqrt(ftx**2 + fty**2)
!
ftx = ftx/lenmc
fty = fty/lenmc
!
othog = abs(fnx*ftx + fny*fty)

!...Get the normal vector at the non-vertex gauss point
call getvectf(xpf, vngf)

!...Boundary face
if(ifa.le.nbfac)then

!...Parameters for Left cell
rhoct = 1.d0/unkno(1, 1, iel)
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...Identify the fgaus for left cell
call getfg_glb(ipqua(1:nvqua, iel), ipf, fgausl)
!
!if(iel.eq.593)print*,'iel593-boundary',fgausl,ipqua(1:nvqua, iel),ipf

!
unkng = 0.d0

!...Basis fctn...
do ig =nvfac+1 ,ngausf
bq(1, ig) = 1.d0
bq(2, ig) = (posiq(1, fgausl(ig))-rc)/dr
bq(3, ig) = (posiq(2, fgausl(ig))-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, iel)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, iel)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, iel)
endif
enddo

!...Unkno at gauss points...
do ig = nvfac+1, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*bq(ideg, ig)
enddo

!...The normal vector at the non-vertex gauss point
vngs(1:2, ig, 1) = vngf(1:2, ig)
vngs(  3, ig, 1) = vngf(  3, ig)
enddo

!...Limiter...
do ig = nvfac+1, ngausf
!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6.and.geoel(10, iel).gt.10.d0)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, iel)*(unkng(1, ig, 1) - rhomc)
unkng(1, ig, 1) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, iel)*unkno(2,2,iel) +  afvec(1, 2, iel)*unkno(2,3,iel)
duds = afvec(1, 1, iel)*unkno(3,2,iel) +  afvec(1, 2, iel)*unkno(3,3,iel)
dvdr = afvec(2, 1, iel)*unkno(2,2,iel) +  afvec(2, 2, iel)*unkno(2,3,iel)
dvds = afvec(2, 1, iel)*unkno(3,2,iel) +  afvec(2, 2, iel)*unkno(3,3,iel)
!
uvtx = unkno(1,2,iel)  + dudr*bq(2, ig) + duds*bq(3, ig)
vvtx = unkno(1,3,iel)  + dvdr*bq(2, ig) + dvds*bq(3, ig)

pvtx = pctr + aflim(4, iel)*(pvtx - pctr)

!...update the velocity at the non-vertex g-s
unkng(2, ig, 1) = uvtx
unkng(3 ,ig, 1) = vvtx
endif

unkng(5, ig, 1) = pvtx
unkng(6, ig, 1) = rhoct*sdctr
unkng(7, ig, 1) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig, 1) = -pvtx
sigmg(1, 2, ig, 1) = 0.d0
sigmg(2, 1, ig, 1) = 0.d0
sigmg(2, 2, ig, 1) = -pvtx!
enddo

!...Get the averaged velocity at the gauss quadrature points...
do ig = nvfac+1, ngausf
!vlave(1, fgausl(ig), iel) = unkng(2, ig, 1)
!vlave(2, fgausl(ig), iel) = unkng(3, ig, 1)
enddo

!...Impedence...
do ig = nvfac+1, ngausf
do ic = 1, 1
!
!aujmp(1:2, ig, ic) = vlave(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgpq(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)
!
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) !+ unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0

!...Excluding the non-vertex gs point
do ig = nvfac+1, ngausf
do ic = 1, 1
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
!
enddo

!...Get
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig)
rhsu2 = munacu(2, ig) - snsigm(2, ig)
!
!dufgq(1, fgausl(ig), iel) =-ufgpq(1, fgausl(ig), iel) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
!dufgq(2, fgausl(ig), iel) =-ufgpq(2, fgausl(ig), iel) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)
!
!if(iel.eq.593)print*,'dufgq',ifa,iel,ig,fgausl(ig),sigmg(1:2, 1:2, 4:5, 1),snsigm_rie(1)!/vngs(3, ig, 1)/vngs(1, ig, 1)

enddo

!...Interior face
elseif(ifa.gt.nbfac)then

!...Parameters for Left cell
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...

!...Identify the fgaus for left cell
call getfg_glb(ipqua(1:nvqua, iel), ipf, fgausl)
!
unkng = 0.d0

!...Basis fctn...
do ig =nvfac+1 ,ngausf
bq(1, ig) = 1.d0
bq(2, ig) = (posiq(1, fgausl(ig))-rc)/dr
bq(3, ig) = (posiq(2, fgausl(ig))-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, iel)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, iel)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, iel)
endif
enddo

!...Unkno at gauss points...
do ig = nvfac+1, ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*bq(ideg, ig)
enddo

!...The normal vector at the non-vertex gauss point
vngs(1:2, ig, 1) = vngf(1:2, ig)
vngs(  3, ig, 1) = vngf(  3, ig)
enddo

!...Limiter...
do ig = nvfac+1, ngausf
!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6.and.geoel(10, iel).gt.10.d0)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, iel)*(unkng(1, ig, 1) - rhomc)
unkng(1, ig, 1) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, iel)*unkno(2,2,iel) +  afvec(1, 2, iel)*unkno(2,3,iel)
duds = afvec(1, 1, iel)*unkno(3,2,iel) +  afvec(1, 2, iel)*unkno(3,3,iel)
dvdr = afvec(2, 1, iel)*unkno(2,2,iel) +  afvec(2, 2, iel)*unkno(2,3,iel)
dvds = afvec(2, 1, iel)*unkno(3,2,iel) +  afvec(2, 2, iel)*unkno(3,3,iel)
!
uvtx = unkno(1,2,iel)  + dudr*bq(2, ig) + duds*bq(3, ig)
vvtx = unkno(1,3,iel)  + dvdr*bq(2, ig) + dvds*bq(3, ig)
!
pvtx = pctr + aflim(4, iel)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig, 1) = uvtx
unkng(3 ,ig, 1) = vvtx
endif
!
unkng(5, ig, 1) = pvtx
unkng(6, ig, 1) = rhoct*sdctr
unkng(7, ig, 1) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig, 1) = -pvtx
sigmg(1, 2, ig, 1) = 0.d0
sigmg(2, 1, ig, 1) = 0.d0
sigmg(2, 2, ig, 1) = -pvtx!
enddo

!...Right cell
rhoct = 1.d0/unkno(1, 1, ier)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ier)
vctr  = unkno(1, 3, ier)
ectr  = unkno(1, 4, ier)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center..

!...Identify the fgaus for right cell
call getfg_glb(ipqua(1:nvqua, ier), ipf, fgausr)
!
do ig =nvfac+1 ,ngausf
bq(1, ig) = 1.d0
bq(2, ig) = (posiq(1, fgausr(ig))-rc)/dr
bq(3, ig) = (posiq(2, fgausr(ig))-sc)/ds
!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, ier)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, ier)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, ier)
endif
enddo
!
do ig = nvfac+1 , ngausf
do ideg =1, mdegr
unkng(1:nq, ig, 2) = unkng(1:nq, ig, 2) + unkno(ideg,1:nq,ier)*bq(ideg, ig)
enddo

!...The normal vector at the non-vertex gauss point
vngs(1:2, ig, 2) =-vngf(1:2, ig)
vngs(  3, ig, 2) = vngf(  3, ig)
enddo

!...Limiter...
do ig = nvfac+1, ngausf
!
rhovt = 1.d0/unkng(1, ig, 2)
uvtx = unkng(2, ig, 2)
vvtx = unkng(3, ig, 2)
evtx = unkng(4, ig, 2)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6.and.geoel(10, ier).gt.10.d0)then
!
rhomc = 1.d0/rhoct
!
rhomv = rhomc + aflim(1, ier)*(unkng(1, ig, 2) - rhomc)
unkng(1, ig, 2) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1, ier)*unkno(2,2,ier) +  afvec(1, 2, ier)*unkno(2,3,ier)
duds = afvec(1, 1, ier)*unkno(3,2,ier) +  afvec(1, 2, ier)*unkno(3,3,ier)
dvdr = afvec(2, 1, ier)*unkno(2,2,ier) +  afvec(2, 2, ier)*unkno(2,3,ier)
dvds = afvec(2, 1, ier)*unkno(3,2,ier) +  afvec(2, 2, ier)*unkno(3,3,ier)
!
uvtx = unkno(1,2,ier)  + dudr*bq(2, ig) + duds*bq(3, ig)
vvtx = unkno(1,3,ier)  + dvdr*bq(2, ig) + dvds*bq(3, ig)
!
pvtx = pctr + aflim(4, ier)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig, 2) = uvtx
unkng(3, ig, 2) = vvtx
!
endif
!
unkng(5, ig, 2) = pvtx
unkng(6, ig, 2) = rhoct*sdctr
unkng(7, ig, 2) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig, 2) = -pvtx
sigmg(1, 2, ig, 2) = 0.d0
sigmg(2, 1, ig, 2) = 0.d0
sigmg(2, 2, ig, 2) = -pvtx!
enddo

!...Get the averaged velocity at the gauss quadrature points...
do ig = nvfac+1, ngausf
!vlave(1, fgausl(ig), iel) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
!vlave(2, fgausl(ig), iel) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
!!
!vlave(1, fgausr(ig), ier) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
!vlave(2, fgausr(ig), ier) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
enddo
!
!...Impedence...
!
do ig = nvfac+1, ngausf
do ic = 1, 2
!
!aujmp(1:2, ig, ic) = vlave(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgpq(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)

!...Impedence
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) !+ unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = nvfac+1, ngausf
do ic = 1, 2
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)

!
if(ic.eq.1)then
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)

!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
else
!
munaclq(1:2, 1, fgausr(ig), ier) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausr(ig), ier) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausr(ig), ier)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausr(ig), ier)    = snsigm_rie(1:2)
!
endif
!
enddo

!if(iel.eq.593)print*,'dufgq',ifa,ier,ig,fgausl(ig),sigmg(1:2, 1:2, 4:5, 1),snsigm_rie(1)!/vngs(3, ig, 1)/vngs(1, ig, 1),&
!vngs(1:2, ig, 1)
!if(ier.eq.593)print*,'dufgqr',ifa,iel,ig,fgausr(ig),sigmg(1:2, 1:2, 4:5, 2),snsigm_rie(1)!/vngs(3, ig, 1)/vngs(1, ig, 2),&
!vngs(1:2, ig, 2)

!...Get the corrected velocity at the non-vertex gauss point
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig)
rhsu2 = munacu(2, ig) - snsigm(2, ig)
!
dufgq(1, fgausl(ig), iel) =-ufgpq(1, fgausl(ig), iel) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
dufgq(2, fgausl(ig), iel) =-ufgpq(2, fgausl(ig), iel) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)
!
dufgq(1, fgausr(ig), ier) =-ufgpq(1, fgausr(ig), ier) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
dufgq(2, fgausr(ig), ier) =-ufgpq(2, fgausr(ig), ier) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)
!
!if(fgausr(ig).eq.13.and.ier.eq.100)then

!print*,'bad',dufgq(1:2, fgausr(ig), ier)
!endif
!
!
enddo

endif
450 enddo
!
!print*,'dufg2',dufgq(1, 17:18, 99)
!
end subroutine getvelo_quadsubg_gs
!
!...Calculate the normal vector at the gauss point
!
subroutine getvectf(xpf, vngf)
use constant
implicit none
real*8, dimension(3, ngausf),   intent(out)::vngf
real*8, dimension(1:2, 1:nvfac), intent(in)::xpf
!...Local integer
integer::ig, ishp
integer::iv
!...Local real
real*8::dshpr(3),posig(2)
real*8::dxdr, dydr,djaco,anx,any,r
!
posig(1) = -0.6546536707079771437983d0
posig(2) =  0.6546536707079771437983d0
!
do ig = nvfac+1, ngausf
!
r  = posig(ig-nvfac)

dshpr(1) = -0.5d0 + r
dshpr(2) =  0.5d0 + r
dshpr(3) = -2.d0*r

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dshpr(ishp)*xpf(1, ishp)
dydr = dydr + dshpr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx = dydr
any =-dxdr

!...Normal vector at the gauss point
vngf(1, ig) = anx/sqrt(anx**2 + any**2)
vngf(2, ig) = any/sqrt(anx**2 + any**2)
vngf(3, ig) = djaco*2.d0

enddo

end subroutine getvectf
!
!...Find the global no. of one face gauss point
!
subroutine getfg_glb(ipqua, ipf, fgaus)
use constant
implicit none
integer, dimension(nvqua), intent(in):: ipqua
integer, dimension(nvfac), intent(in):: ipf
integer:: ig, nsum
integer, dimension(1:ngausf), intent(out)::fgaus
!
integer::fglgq(ngausf,4)

!...Give fglgq
if(ngausf.eq.5)then

fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;
endif

!...Find the No. of face gauss points in fglgq...
nsum = ngausf +4
!
do ig = nvfac+1, ngausf
!
if(ipqua(5).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(1))then
fgaus(ig) = fglgq(ig,1)
else
fgaus(ig) = fglgq(nsum-ig,1)
endif
!
elseif(ipqua(6).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(2))then
fgaus(ig) = fglgq(ig,2)
else
fgaus(ig) = fglgq(nsum-ig,2)
endif
!
elseif(ipqua(7).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(3))then
fgaus(ig) = fglgq(ig,3)
else
fgaus(ig) = fglgq(nsum-ig,3)
endif
!
elseif(ipqua(8).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(4))then
fgaus(ig) = fglgq(ig,4)
else
fgaus(ig) = fglgq(nsum-ig,4)
endif
!
endif

enddo
end subroutine getfg_glb
!
!...subroutine: Riemann input for hybrid curved quads using sub-cell scheme for 5-points integration....
!
subroutine getriem_quadsubg_gauss(ipqua, geoel, gesgq, vlave, unkno, unkgd, strnq_devtp,munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, drhosgq, coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),  intent(in)::unkgd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq, nquad)::strnq_devtp
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
real*8, dimension(1:2, 1:2, 1:2, 1:4, 1:4, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:4, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:2,  1:4, 1:4, 1:nquad), intent(out)::snsigmlq
real*8,dimension(1:4, 1:nsize), intent(out)::drhosgq

!...Local integer
integer::ie, ideg, ielem, ifa, iv, isg, ivsg, ifsg

!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(8, 4)::fnqsg
integer,dimension(4, 4)::ipqsg,ipqsg_strn
!...local real array
real*8,dimension(1:ndegr, 1:nvqua)::bq
real*8,dimension(1:ndegr, 1:4)::bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nq,1:4)::unsgq
real*8::aujmp(1:3, 1:4)
real*8::vnorm(1:3, 1:2, 1:4)
real*8::sigma(1:2, 1:2, 1:4)
real*8,dimension(1:2, 1:4)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndegr, 1:4)::unksgq
!real*8,dimension(1:ndegr, 1:4)::bqsg
real*8,dimension(1:nq+1, 1:4)::uqsgc
real*8,dimension(1:4, 1:4)::prsgq
real*8,dimension(1:4, 1:4)::bqvp
real*8,dimension(1:4)::prqz
real*8,dimension(1:7, 1:4)::geoq_sub
real*8,dimension(2, 4, 4)::wfgsq
real*8::bqsg(ndegr)


!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2), sigma_devt(1:2,1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr, eintc
real*8::rhomv,uvtx,vvtx,evtx, pvtx,rhovsg,eintv,sdv
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv,rci,sci
real*8::rcsg,scsg
real*8::dxp,dyp,xc,yc,xcrho,ycrho
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon
!
data eps   / 1.0d-6 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Specify some Gauss points...
!
!...Local vertex No. of gauss points in one unit
ipqsg_strn(1, 1) = 1; ipqsg_strn(2, 1) = 0; ipqsg_strn(3, 1) = 0; ipqsg_strn(4, 1) = 0
ipqsg_strn(1, 2) = 1; ipqsg_strn(2, 2) = 1; ipqsg_strn(3, 2) = 0; ipqsg_strn(4, 2) = 0
ipqsg_strn(1, 3) = 0; ipqsg_strn(2, 3) = 1; ipqsg_strn(3, 3) = 1; ipqsg_strn(4, 3) = 0
ipqsg_strn(1, 4) = 1; ipqsg_strn(2, 4) = 0; ipqsg_strn(3, 4) = 1; ipqsg_strn(4, 4) = 1

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
if(nfint.eq.5.or.nfint.eq.10)then
wfgsq(1, 1, 1) = 1.d0;      wfgsq(2, 1, 1) = 1.d0;
wfgsq(1, 2, 1) = .5d0;      wfgsq(2, 2, 1) = 1.d0;
wfgsq(1, 3, 1) = 0.d0;      wfgsq(2, 3, 1) = 0.d0;
wfgsq(1, 4, 1) = 1.d0;      wfgsq(2, 4, 1) = .5d0;

wfgsq(1, 1, 2) = 1.d0;      wfgsq(2, 1, 2) = .5d0;
wfgsq(1, 2, 2) = 1.d0;      wfgsq(2, 2, 2) = 1.d0;
wfgsq(1, 3, 2) = .5d0;      wfgsq(2, 3, 2) = 1.d0;
wfgsq(1, 4, 2) = 0.d0;      wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 1, 3) = 0.d0;      wfgsq(2, 1, 3) = 0.d0;
wfgsq(1, 2, 3) = 1.d0;      wfgsq(2, 2, 3) = .5d0;
wfgsq(1, 3, 3) = 1.d0;      wfgsq(2, 3, 3) = 1.d0;
wfgsq(1, 4, 3) = .5d0;      wfgsq(2, 4, 3) = 1.d0

wfgsq(1, 1, 4) = .5d0;      wfgsq(2, 1, 4) = 1.d0;
wfgsq(1, 2, 4) = 0.d0;      wfgsq(2, 2, 4) = 0.d0;
wfgsq(1, 3, 4) = 1.d0;      wfgsq(2, 3, 4) = .5d0;
wfgsq(1, 4, 4) = 1.d0;      wfgsq(2, 4, 4) = 1.d0;

elseif(nfint.eq.6)then
!xxxxxxxxx
wfgsq(1, 1, 1) = 1.d0;      wfgsq(2, 1, 1) = 1.d0;
wfgsq(1, 2, 1) = .5d0;      wfgsq(2, 2, 1) = 4.d0;
wfgsq(1, 3, 1) = 0.d0;      wfgsq(2, 3, 1) = 0.d0;
wfgsq(1, 4, 1) = 4.d0;      wfgsq(2, 4, 1) = .5d0;

wfgsq(1, 1, 2) = 4.d0;      wfgsq(2, 1, 2) = .5d0;
wfgsq(1, 2, 2) = 1.d0;      wfgsq(2, 2, 2) = 1.d0;
wfgsq(1, 3, 2) = .5d0;      wfgsq(2, 3, 2) = 4.d0;
wfgsq(1, 4, 2) = 0.d0;      wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 1, 3) = 0.d0;      wfgsq(2, 1, 3) = 0.d0;
wfgsq(1, 2, 3) = 4.d0;      wfgsq(2, 2, 3) = .5d0;
wfgsq(1, 3, 3) = 1.d0;      wfgsq(2, 3, 3) = 1.d0;
wfgsq(1, 4, 3) = .5d0;      wfgsq(2, 4, 3) = 4.d0

wfgsq(1, 1, 4) = .5d0;      wfgsq(2, 1, 4) = 4.d0;
wfgsq(1, 2, 4) = 0.d0;      wfgsq(2, 2, 4) = 0.d0;
wfgsq(1, 3, 4) = 4.d0;      wfgsq(2, 3, 4) = .5d0;
wfgsq(1, 4, 4) = 1.d0;      wfgsq(2, 4, 4) = 1.d0;
endif
!...Initilize the stress
sigma_devt  = 0.d0
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

!...The initial cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)

!...Vertex coordinate
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

!...Basis function
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
endif
enddo

!...Get density for sub-cells...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

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

eintc  = ectr - 0.5d0*(uctr**2 + vctr**2)
!
!...Call EOS
!
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)
call getrhog_initial(rhoi,  xcrho, ycrho, xcrho, ycrho)
!call getrhoig_quad(rhoi, r, s, xpqi)
!rhoi = 1.845d0
call GetEOS(nmatel,ncase,gamlg,rhoct, eintc, rhoi, pctr, sdctr, ielem)

!pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

!... u, v, e at every vertex....
unknvq = 0.d0
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
enddo

!if(ielem.eq.12)then
!print*,'smooth2',ielem,unkno(:,:,ielem)
!endif

!
!...Get density correction
!call getdens_quadsubg(rc, sc, geoel(19:21, ielem), unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc)
if(nfint.eq.5)then
!
 call getrhosubg_averg(xpq,xpqi,unksgq, geoq_sub,ielem)
 call getrhosubg_daverg(rc,sc, geoel(19:21, ielem), xpq, unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc, ielem)
elseif(nfint.eq.10)then
!
  call getrhosubcell_sms(xcrho, ycrho, xpq,xpqi,unksgq, ielem)
  call getrhosubg_daverg(rc,sc, geoel(19:21, ielem), xpq, unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc, ielem)

!
!  call getrhosubg_averg(xpq,xpqi,unksgq, geoq_sub,ielem)
!  call getrhosubg_daverg(rc,sc, geoel(19:21, ielem), xpq, unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc, ielem)
!
!!  call getrhosmsgd_averg(rci,sci,geoel(22:24, ielem),xpq,xpqi,unksgq, unkgd, geoq_sub, ielem)
!!  call getrhosms_daverg(rc,sc, rci,sci, geoel(19:21, ielem), geoel(22:24, ielem), &
! xpqi, unksgq, unkno(:,:,ielem), unkgd(:,:,ielem), aflim(:,ielem), uqsgc, ielem)
endif

!
!if(ielem.eq.1) print*,'ielem',unksgq(1,1),uqsgc(1,1),unkno(1:6,1,1)

!...II.1: Loop over sub-cells....
do isg = 1, 4

!...Record the corrected density
drhosgq(isg, ielem) = (unksgq(1,isg)-1.d0/uqsgc(1, isg))

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
!
!if(ielem.eq.1)print*,'ustar',isg,ie,vnorm(1:3, :,:)

!...Correct the density for Reimann input
do ivsg = 1,4
if(ndens.eq.1)then
rhovt = 1.d0/unknvq(1, ipqsg(ivsg, isg))
!rhovsg = unksgq(1,isg)
rhovsg = rhovt + cdrho*(unksgq(1,isg)-1.d0/uqsgc(1, isg))
!
!rcsg = geoq_sub(1, isg)
!scsg = geoq_sub(2, isg)
!Basis function
!bqsg(2) = (xvq(ivsg)-rcsg)/dr
!bqsg(3) = (yvq(ivsg)-scsg)/ds

!...DGP2
!if(npoly.eq.2)then
!bqsg(4) = 0.5d0* bqsg(2)* bqsg(2) - geoq_sub(5, isg)
!bqsg(5) = 0.5d0* bqsg(3)* bqsg(3) - geoq_sub(6, isg)
!bqsg(6) =        bqsg(2)* bqsg(3) - geoq_sub(7, isg)
!endif
!
! rhovsg = 0.d0
!do ideg = 1, mdegr
! rhovsg = rhovsg + unksgq(ideg,isg)*bqsg(ideg)
!enddo

elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bqv(1, ivsg) = 1.d0
bqv(2, ivsg) = (xvq(ipqsg(ivsg, isg))-rcv)/dr  !...bqv ....
bqv(3, ivsg) = (yvq(ipqsg(ivsg, isg))-scv)/ds
!
rhovt = 0.d0
do ideg = 1, mdegr
rhovt = rhovt + unkno(ideg,1,ielem)*bqv(ideg, ivsg)
enddo
rhovsg = rhovt + (unksgq(1,isg)-rhoct)

!...Based on the sub-cell
!rcv = geoq_sub(1, isg); scv = geoq_sub(2, isg)
!bqv(1, ivsg) = 1.d0
!bqv(2, ivsg) = (xvq(ivsg)-rcv)/dr  !...bqv ....
!bqv(3, ivsg) = (yvq(ivsg)-scv)/ds

!rhovsg =0.d0
!do ideg = 1,ndegr
!rhovsg =  rhovsg + unksgq(ideg,isg)*bqv(ideg, ivsg)
!enddo

!rhovsg = unksgq(1,isg)! unkno(1,1,ielem)!
endif

uvtx = unknvq(2, ipqsg(ivsg, isg))
vvtx = unknvq(3, ipqsg(ivsg, isg))
evtx = unknvq(4, ipqsg(ivsg, isg))

!...Derived variables at the vertex
eintv = evtx - 0.5d0*(uvtx**2 + vvtx**2)

!...Call EOS
call getrhog_initial(rhoi,  xpqi(1, ipqsg(ivsg, isg)), xpqi(2, ipqsg(ivsg, isg)), xcrho, ycrho)
call GetEOS(nmatel,ncase,gamlg,rhovsg, eintv, rhoi, pvtx, sdv, ielem)

!
!pvtx = max(eps, (gamlg-1.d0)*rhovsg*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx

!...Output for debugging
!if(ielem.eq.1.or.ielem.eq.2.or.ielem.eq.6.or.ielem.eq.7)then
!print*,'Variabes',ielem,unkno(1:3,1:4,ielem)
!print*,'Variabe for GP',ielem,isg,ivsg,rhovsg,pvtx,unsgq(2:3, ivsg),evtx,geoel(10,ielem)
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

!...Fitted pressure...
!pvtx = prsgq(ivsg, isg)

!...Updtae velocity
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx
endif

!
if(nmatel.eq.2)then
call GetStress_deviat_infsmal(xvq(ipqsg(ivsg, isg)),yvq(ipqsg(ivsg, isg)),rci,sci,geoel(22:24, ielem),&
sigma_devt, &
strnq_devtp(:,:, ipqsg(ivsg, isg),ielem), unkgd(:,:,ielem),ielem,ipqsg_strn(ivsg, isg))
!  call GetStress_deviatfem_infsmal(xvq(ipqsg(ivsg, isg)),yvq(ipqsg(ivsg, isg)),xpq,xpqi,sigma_devt,&
!                                   strnq_devtp(:,:, ipqsg(ivsg, isg),ielem), ielem)
endif

!...Get stress tensor at one vertex
sigma(1, 1, ivsg) = -pvtx + sigma_devt(1,1)
sigma(1, 2, ivsg) = 0.d0  + sigma_devt(1,2)
sigma(2, 1, ivsg) = 0.d0  + sigma_devt(2,1)
sigma(2, 2, ivsg) = -pvtx + sigma_devt(2,2)

!
!...Output for debugging
!if(ielem.eq.115)then
!print*,'Variabe2',ielem,isg,ivsg,rhovsg,pvtx,unsgq(2:3, ivsg),evtx,afvec(:, :, ielem)
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
!sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedence coefficient...
do ivsg   = 1, 4
dux= vlave(1, ipq(ipqsg(ivsg, isg)))-unsgq(2, ivsg)
duy= vlave(2, ipq(ipqsg(ivsg, isg)))-unsgq(3, ivsg)
deltu = sqrt(dux**2 + duy**2)
do ifa = 1, 2
deltu = cimpd*abs(dux*vnorm(1, ifa, ivsg) + duy*vnorm(2, ifa, ivsg))
murie(ifa, ivsg) = rhoct*sdctr + rhoct*slpdu*deltu
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
!if(ipq(iv).eq.2)then
!print*,'p36 muacn(vv) post',ie,isg,ivsg,ipq(iv),munacn(1, 1, ipq(iv)),munacn_rie(1, 1),murie(ifa, ivsg),vnorm(3, ifa, ivsg),&
!drhosgq(isg, ielem),&
!snsigm(1:2, ipq(iv)),sigma(1, 1, ivsg),unsgq(2:3, ivsg),munacu(1:2, ipq(iv)),munacu_rie(1:2)
!endif

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

end subroutine getriem_quadsubg_gauss
!
!...Get the nodal velocity based on 5-nodes Lobatto Gauss integration...
!
subroutine getndvelo_lagsubg_gaussh(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno, unkgd,strnq_devtp, ustar, ufgpq, fstart, fstarq,fstfgdu,fsqfgdu,fstfg,fsqfg,  aflim, afvec, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac), intent(in)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq, nquad)::strnq_devtp
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:4,1:3, 1:ntria),  intent(out)::fstart !...Riemann forces
real*8,dimension(1:ndimn,1:2,1:4,1:4, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:ndimn,1:ngsft, 1:ntria),   intent(out)::fstfgdu,fstfg !...Riemann forces
real*8,dimension(1:ndimn,1:ngsfq, 1:nquad),   intent(out)::fsqfgdu,fsqfg !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngsfq,1:nquad),       intent(out)::ufgpq

integer:: itime,ip
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop
integer::isg,ivsg
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(4, 4)::ipqsg
integer::indnd(npoin)
integer, dimension(3, 4)::fglvq
integer, dimension(3, 3)::fglvt
!
integer, dimension(ngausf, 4):: fglgq
integer, dimension(ngausf, 3):: fglgt

!...local real array
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:ndimn,1:nvqua,1:nquad)::vnulq
real*8,  dimension(1:nquad)::gqdmp
real*8::munaci(2, 2)
!...local real number
real*8::xp,yp,radip
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2
real*8::r,shp1,shp2,shp3
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: munacn(:,:,:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: munacu(:,:), snsigm(:,:)
real*8,allocatable:: munaclt(:,:,:,:),munaclq(:,:,:,:,:,:)
real*8,allocatable:: snsigmlt(:,:,:), munault(:,:,:)
real*8,allocatable:: snsigmlq(:,:,:,:,:), munaulq(:,:,:,:,:)
real*8,allocatable:: dufgq(:,:,:)
real*8,allocatable:: munaclq_fg(:,:,:,:),munaulq_fg(:,:,:),snsigmlq_fg(:,:,:)
real*8,allocatable:: drhosgq(:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:2, 1:2, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munaclt(1:2, 1:2, 1:ngelg, 1:ntria), munault(1:ndimn, 1:ngelg,  1:ntria),&
snsigmlt(1:ndimn, 1:ngelg,  1:ntria))
allocate (munaclq(1:2, 1:2, 1:2, 1:4, 1:4, 1:nquad), munaulq(1:ndimn, 1:2, 1:4, 1:4,  1:nquad),&
snsigmlq(1:ndimn, 1:2,  1:4, 1:4,  1:nquad))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
allocate (dufgq(1:2,1:ngsfq,1:nquad))
allocate (munaclq_fg(1:2, 1:2, 1:ngsfq, 1:nquad),munaulq_fg(1:ndimn, 1:ngsfq,  1:nquad),snsigmlq_fg(1:ndimn, 1:ngsfq,  1:nquad))
allocate(drhosgq(1:4, 1:nsize))
!
!...Part I: Specify some gauss points...
!

!...Local vertex No. of gauss points in one unit quad...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;

!...Local vertex No. of gauss points in one unit tria...
fglvt(1, 1) = 1;  fglvt(2, 1) = 2; fglvt(3, 1) = 4;
fglvt(1, 2) = 2;  fglvt(2, 2) = 3; fglvt(3, 2) = 5;
fglvt(1, 3) = 3;  fglvt(2, 3) = 1; fglvt(3, 3) = 6;

!...Local gauss point No. of any gauss point in one face...
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;
!
fglgt(1, 1) = 1;  fglgt(2, 1) = 2; fglgt(3, 1) = 7;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4; fglgt(3, 2) = 8;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6; fglgt(3, 3) = 9;
!
elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
elseif(ngausf.eq.5)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;
!
endif

!
!...Local vertex No. of gauss points in one sub-cell
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
!...Zero out vlave
!
vlave = 0.d0
indnd = 0
!
!...Mark the boundary nodes...
!
if(nmatel.eq.1)then
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

endif

!...Get averaged velocity at nodes...
!call getvlavenew(iptri, ipqua, geoel, vlave, unkno, aflim, afvec)
!
do 950 ifa = 1 , nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
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
endif!
950 enddo
!
!...Part II: Loop to get information from Riemann solver
!
do iloop= 1, 1!5

!...Give value
vlave= ustar
vnulq = 0.d0

!...Zero out munacn
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0

!...Triangle in the future...

!...Quad
if(nfint.eq.5.or.nfint.eq.10)then

if(nquad.gt.0) call getriem_quadsubg_gauss(ipqua, geoel, gesgq, vlave, unkno, unkgd, strnq_devtp, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, drhosgq, coord, coold, aflim, afvec)

elseif(nfint.eq.7)then

 if(nquad.gt.0) call getriem_quadsubg_gauss_smth(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
 munaclq, munaulq, snsigmlq, drhosgq, coord, coold, aflim, afvec)
endif
!
call getbc_lagc(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)

!...Periodic boundary condition for 1D isentropic Sin problem...
if(ncase.eq.12) call getbc_prdic(bface, munacn, munacu, snsigm)

!....Update the velocity at the end points...
do ipoin = 1,npoin

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

!if(ipoin.eq.2)print*,'ipnode',itime,ipoin,ustar(1:2,ipoin),detma,munacn(:,:,ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin)
endif
enddo

!...Get the vertex velocity at the boundary....
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
endif
endif
!
900 enddo

enddo
!
!
!print*,'velocity',ustar(1:2,2)
!
!...Update ufgpq using array ustar ...
!
do 100 ie=1, nquad !...(1)ifa=1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
do ifa = 1,4
!
!ufgaus(1:ndimn,fglgq(1, ifa),  ie) = ustar(1:2, ipq(fglvq(1, ifa)))
!ufgaus(1:ndimn,fglgq(2, ifa),  ie) = ustar(1:2, ipq(fglvq(2, ifa)))
!
r = -0.6546536707079771437983d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
!
ufgpq(1:ndimn,fglgq(4, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
shp3*ustar(1:2, ipq(fglvq(3, ifa)))

r =  0.6546536707079771437983d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
!
ufgpq(1:ndimn,fglgq(5, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
shp3*ustar(1:2, ipq(fglvq(3, ifa)))
enddo

! print*,'ie1830',ie,ifa,ufgpq(1:2,fglgq(5, 2),ie)

100 enddo

!...Get the Riemnan input for the quadrature points.
dufgq = 0.d0

call getvelo_quadsubg_gs3(ufgpq,dufgq,geoel,bface,intfac,ipqua,coord,coold,unkno,unkgd,strnq_devtp, indnd, &
munaclq_fg, munaulq_fg, snsigmlq_fg, drhosgq, afvec, aflim, itime)
!call getvelo_quadsubg_gs(ufgpq,dufgq,geoel,intfac,ipqua,coord,unkno,indnd, &
!munaclq_fg, munaulq_fg, snsigmlq_fg,  afvec, aflim)

!print*,'ie1830du',dufgq(1:2,13:20,100)
!stop

!...Zero out correction force
fsqfgdu = 0.d0

!
!...4.2: Update the Riemann forces at every node...
!...Tria
!
!...Quad
!
do ie = 1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!....Zero out the velocity at the 9th node inside one cell to avoid NaN...
ustar(1:2, ipq(9)) = 0.d0

!
do isg = 1, 4
do ivsg = 1, 4
!
iv =  ipqsg(ivsg, isg)
!
do ifa =1, 2
fstarq(1, ifa, ivsg, isg, ie) = snsigmlq(1, ifa, ivsg, isg, ie) + &
munaclq(1, 1, ifa, ivsg, isg,ie)*ustar(1, ipq(iv))+&
munaclq(2, 1, ifa, ivsg, isg,ie)*ustar(2, ipq(iv)) - munaulq(1, ifa, ivsg, isg, ie)
fstarq(2, ifa, ivsg, isg,ie) = snsigmlq(2, ifa, ivsg, isg, ie) + &
munaclq(1, 2, ifa, ivsg, isg, ie)*ustar(1, ipq(iv))+&
munaclq(2, 2, ifa, ivsg, isg, ie)*ustar(2, ipq(iv))- munaulq(2, ifa, ivsg, isg, ie)

!...Output for debugging
!if(ie.eq.1.)print*,'force vtv',ie,isg,ivsg,ifa,ipq(iv),fstarq(1, ifa, ivsg, isg, ie),&
!snsigmlq(1, ifa, ivsg, isg, ie),munaclq(1:2, 1, ifa, ivsg, isg,ie),ustar(1:2, ipq(iv)), munaulq(1, ifa, ivsg, isg, ie)
!
!if(ie.eq.1.)print*,'force vtx',ie,isg,ivsg,ifa,ipq(iv),fstarq(2, ifa, ivsg, isg, ie),&
!snsigmlq(2, ifa, ivsg, isg, ie),munaclq(1:2, 2, ifa, ivsg, isg,ie),ustar(1:2, ipq(iv)), munaulq(2, ifa, ivsg, isg, ie)
!
enddo
enddo
enddo

!...Non-vertex gauss point
do ifa = 1, 4
do ig = nvfac+1, ngausf
!
fsqfg(1, fglgq(ig, ifa), ie) = snsigmlq_fg(1, fglgq(ig, ifa), ie) +&
munaclq_fg(1,1, fglgq(ig, ifa), ie)*ufgpq(1, fglgq(ig, ifa), ie)+&
munaclq_fg(2,1, fglgq(ig, ifa), ie)*ufgpq(2, fglgq(ig, ifa), ie)-&
munaulq_fg(1, fglgq(ig, ifa), ie)
!
fsqfg(2, fglgq(ig, ifa), ie) = snsigmlq_fg(2, fglgq(ig, ifa), ie) +&
munaclq_fg(2,2,fglgq(ig, ifa), ie)*ufgpq(2, fglgq(ig, ifa), ie)+&
munaclq_fg(1,2,fglgq(ig, ifa), ie)*ufgpq(1, fglgq(ig, ifa), ie)-&
munaulq_fg(2, fglgq(ig, ifa), ie)

!...Output for debugging
!if(ie.eq.1.)print*,'force gp',ig,ifa,fglgq(ig, ifa),fsqfg(1, fglgq(ig, ifa), ie),&
!snsigmlq_fg(1, fglgq(ig, ifa), ie),&
!munaclq_fg(1,1, fglgq(ig, ifa), ie), ufgpq(1, fglgq(ig, ifa), ie),munaulq_fg(1, fglgq(ig, ifa), ie)

!if(ie.eq.1.)print*,'force gp',ig,ifa,fglgq(ig, ifa),fsqfg(2, fglgq(ig, ifa), ie),&
!snsigmlq_fg(2, fglgq(ig, ifa), ie),&
!munaclq_fg(2,2, fglgq(ig, ifa), ie), ufgpq(2, fglgq(ig, ifa), ie),munaulq_fg(2, fglgq(ig, ifa), ie)
!
enddo

!...Non-vertex gauss point
do ig = nvfac+1, ngausf
!
fsqfgdu(1, fglgq(ig, ifa), ie) = munaclq_fg(1,1, fglgq(ig, ifa), ie)*dufgq(1, fglgq(ig, ifa), ie)+&
munaclq_fg(2,1, fglgq(ig, ifa), ie)*dufgq(2, fglgq(ig, ifa), ie)
!
fsqfgdu(2, fglgq(ig, ifa), ie) = munaclq_fg(2,2,fglgq(ig, ifa), ie)*dufgq(2, fglgq(ig, ifa), ie)+&
munaclq_fg(1,2,fglgq(ig, ifa), ie)*dufgq(1, fglgq(ig, ifa), ie)

!...Output for debugging
!if(ie.eq.1.)print*,'force correct',ig,ifa,fglgq(ig, ifa),fsqfgdu(1:2, fglgq(ig, ifa), ie),&
!munaclq_fg(1,1, fglgq(ig, ifa), ie),dufgq(1:2, fglgq(ig, ifa), ie)
!
enddo
enddo

enddo
!
deallocate (munacn, bpres, fpres)
deallocate (munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
deallocate (dufgq)
deallocate (munaclq_fg,munaulq_fg,snsigmlq_fg)
deallocate (drhosgq)
end subroutine getndvelo_lagsubg_gaussh
!
!...Face integral for sub-cell scheme with both vertex and non-vertex gauss point...
!
subroutine rhsifacedg_lagsubgq_gs(ipqua, unkno, ustar, ufgpq,fstarq, fsqfg, fsqfgdu, gesgq, geoel, coord,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:2,1:ngsfq,nquad),         intent(in)::ufgpq
real*8,dimension(1:ndimn,1:2,1:4, 1:4, 1:nquad),intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndimn,1:ngsfq, 1:nquad),   intent(in)::fsqfgdu,fsqfg !...Riemann forces

real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord

!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem,isg,ivsg,ifsg
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac, 1:nvqua) :: ipfq
integer,dimension(8, 4)::fnqsg
integer,dimension(4, 4)::ipqsg
integer,dimension(ngausf,4)::fglgq
!
real*8::posiq(1:2, 1:ngsfq)
real*8,dimension(2,4,4)::wfgsqv,wfgsq
real*8,dimension(ngsfq)::wfgsqg
real*8,dimension(1:3,1:2,1:4)::vnorm
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::xvq(nvqua), yvq(nvqua),bq(1:ndegr,1:nvqua)
real*8::bg(ndegr,  ngausf)
real*8::xpf(2,nvfac)
real*8::vngf(3,ngausf)
real*8::posi(ngausf), weigh(ngausf)

!...local real number
real*8::eps,c00,c05,c10,c20,c13,c16
real*8::dr,ds,rc,sc, rg, sg, wi
real*8::gpnx,gpny,gpsa,sarea
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
data c16   / 0.1666666666666666d0 /
!
!...Part I: Specify some gauss pints
!
call ruqope_lobatto(1, ngausf, posi, weigh)
!
ipfq(1, 1) = 1; ipfq(2, 1) = 2; ipfq(3, 1) = 5
ipfq(1, 2) = 2; ipfq(2, 2) = 3; ipfq(3, 2) = 6
ipfq(1, 3) = 3; ipfq(2, 3) = 4; ipfq(3, 3) = 7
ipfq(1, 4) = 4; ipfq(2, 4) = 1; ipfq(3, 4) = 8
!
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
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) = posi(2); posiq(2,13) =-1.d0;
posiq(1,14) = posi(4); posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) = posi(2);
posiq(1,16) = 1.d0; posiq(2,16) = posi(4);

posiq(1,17) = posi(4); posiq(2,17) = 1.d0;
posiq(1,18) = posi(2); posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = posi(4);
posiq(1,20) =-1.d0; posiq(2,20) = posi(2);
!
wfgsq(1, 1, 1) = weigh(1);      wfgsq(2, 1, 1) = weigh(1);
wfgsq(1, 2, 1) = weigh(3)/1.d0; wfgsq(2, 2, 1) = weigh(3);
wfgsq(1, 3, 1) = 0.d0;          wfgsq(2, 3, 1) = 0.d0;
wfgsq(1, 4, 1) = weigh(3);      wfgsq(2, 4, 1) = weigh(3)/1.d0;

wfgsq(1, 1, 2) = weigh(3);      wfgsq(2, 1, 2) = weigh(3)/1.d0;
wfgsq(1, 2, 2) = weigh(1);      wfgsq(2, 2, 2) = weigh(1);
wfgsq(1, 3, 2) = weigh(3)/1.d0; wfgsq(2, 3, 2) = weigh(3);
wfgsq(1, 4, 2) = 0.d0;          wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 1, 3) = 0.d0;          wfgsq(2, 1, 3) = 0.d0;
wfgsq(1, 2, 3) = weigh(3);      wfgsq(2, 2, 3) = weigh(3)/1.d0;
wfgsq(1, 3, 3) = weigh(1);      wfgsq(2, 3, 3) = weigh(1);
wfgsq(1, 4, 3) = weigh(3)/1.d0; wfgsq(2, 4, 3) = weigh(3);

wfgsq(1, 1, 4) = weigh(3)/1.d0; wfgsq(2, 1, 4) = weigh(3);
wfgsq(1, 2, 4) = 0.d0;          wfgsq(2, 2, 4) = 0.d0;
wfgsq(1, 3, 4) = weigh(3);      wfgsq(2, 3, 4) = weigh(3)/1.d0;
wfgsq(1, 4, 4) = weigh(1);      wfgsq(2, 4, 4) = weigh(1);

!
wfgsq = 0.5d0*wfgsq
!
!weigh(1) = 2.d0/6.d0; weigh(3) = 8.d0/6.d0
!...weight for vertex gauss point
wfgsqv(1, 1, 1) = weigh(1);      wfgsqv(2, 1, 1) = weigh(1);
wfgsqv(1, 2, 1) = weigh(3)/2.d0; wfgsqv(2, 2, 1) = weigh(3);
wfgsqv(1, 3, 1) = 0.d0;          wfgsqv(2, 3, 1) = 0.d0;
wfgsqv(1, 4, 1) = weigh(3);      wfgsqv(2, 4, 1) = weigh(3)/2.d0;

wfgsqv(1, 1, 2) = weigh(3);      wfgsqv(2, 1, 2) = weigh(3)/2.d0;
wfgsqv(1, 2, 2) = weigh(1);      wfgsqv(2, 2, 2) = weigh(1);
wfgsqv(1, 3, 2) = weigh(3)/2.d0; wfgsqv(2, 3, 2) = weigh(3);
wfgsqv(1, 4, 2) = 0.d0;          wfgsqv(2, 4, 2) = 0.d0;

wfgsqv(1, 1, 3) = 0.d0;          wfgsqv(2, 1, 3) = 0.d0;
wfgsqv(1, 2, 3) = weigh(3);      wfgsqv(2, 2, 3) = weigh(3)/2.d0;
wfgsqv(1, 3, 3) = weigh(1);      wfgsqv(2, 3, 3) = weigh(1);
wfgsqv(1, 4, 3) = weigh(3)/2.d0; wfgsqv(2, 4, 3) = weigh(3);

wfgsqv(1, 1, 4) = weigh(3)/2.d0; wfgsqv(2, 1, 4) = weigh(3);
wfgsqv(1, 2, 4) = 0.d0;          wfgsqv(2, 2, 4) = 0.d0;
wfgsqv(1, 3, 4) = weigh(3);      wfgsqv(2, 3, 4) = weigh(3)/2.d0;
wfgsqv(1, 4, 4) = weigh(1);      wfgsqv(2, 4, 4) = weigh(1);
!
wfgsqv = 0.5d0*wfgsqv
!
!...Part II: Loop over the element
!
do 650 ie = 1,nquad
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
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
endif
enddo

!...Initialize ulnpn, plnpn, elnpn
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0

!...Contribution from vertex gauss point
do isg = 1, 4

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
vnorm(3, ifsg, 1) =  vnorm(  3, ifsg, 1)*wfgsqv(ifsg, 1, isg);
vnorm(3, ifsg, 2) =  vnorm(  3, ifsg, 2)*wfgsqv(ifsg, 2, isg);
vnorm(3, ifsg, 3) =  vnorm(  3, ifsg, 3)*wfgsqv(ifsg, 3, isg);
vnorm(3, ifsg, 4) =  vnorm(  3, ifsg, 4)*wfgsqv(ifsg, 4, isg);
enddo

!...loop over every subgrid vertex...
do ivsg = 1, 4
!
ulnpn(1:ndegr) = ulnpn(1:ndegr)  +&
ustar(1, ipq(ipqsg(ivsg, isg)))*vnorm(1, 1, ivsg)*vnorm(  3, 1, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg)) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*vnorm(2, 1, ivsg)*vnorm(  3, 1, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg)) +&
ustar(1, ipq(ipqsg(ivsg, isg)))*vnorm(1, 2, ivsg)*vnorm(  3, 2, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg)) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*vnorm(2, 2, ivsg)*vnorm(  3, 2, ivsg)*bq(1:ndegr,  ipqsg(ivsg, isg))
!
!if(ie.eq.593)print*,'rhsface-gp-v',ifa,ig,ie,ipq(ipqsg(ivsg, isg)),fstarq(1, 1:2, ivsg, isg, ie)

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
!if(ie.eq.1)print*,'rhsface-vx-bf',isg,ivsg,ipq(ipqsg(ivsg, isg)),ie,bq(1:3, ipqsg(ivsg, isg))
!if(ie.eq.1)print*,'rhsface-vertex',isg,ivsg,ipq(ipqsg(ivsg, isg)),ie,plnpn(1, 1:3), fstarq(1:2, 1:2, ivsg, isg, ie),&
!wfgsq(1:2, ivsg, isg)

enddo
!
enddo
!
!if(ie.eq.54)print*,'rhsface-vx',ifa,ig,ie,plnpn(2,1:3)


!...Contribution from non-vertex gauss point
sarea = 0.d0
!
do ifa  = 1, 4
!
xpf(1, 1:nvfac) = coord(1, ipq(ipfq(1:nvfac, ifa)))
xpf(2, 1:nvfac) = coord(2, ipq(ipfq(1:nvfac, ifa)))

!...Get the face area
call getvectf(xpf, vngf)
!
do ig   = nvfac+1, ngausf
!
wi  = weigh(2)/2.d0
!
rg = posiq(1,fglgq(ig, ifa))
sg = posiq(2,fglgq(ig, ifa))
!
bg(1,  ig) = 1.d0
bg(2,  ig) = (rg-rc)/dr
bg(3,  ig) = (sg-sc)/ds
!DGP2
if(npoly.eq.2)then
bg(4, ig) = 0.5d0*bg(2, ig)*bg(2, ig) - geoel(19, ielem)
bg(5, ig) = 0.5d0*bg(3, ig)*bg(3, ig) - geoel(20, ielem)
bg(6, ig) =       bg(2, ig)*bg(3, ig) - geoel(21, ielem)
endif

!...The vertex constituting one cell...
gpnx = vngf(1, ig)
gpny = vngf(2, ig)
gpsa = vngf(3, ig)
!
sarea = sarea + gpnx*gpsa

!if(ie.eq.593)print*,'rhsface-gp-nv',ifa,ig,ie,fsqfg(1, fglgq(ig, ifa), ie),gpnx,gpny,gpsa,sarea,&
!fsqfg(1, fglgq(ig, ifa), ie)/gpsa/gpnx
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ufgpq(1, fglgq(ig, ifa),ie)*gpnx*gpsa*bg(1:ndegr, ig)*wi +&
ufgpq(2, fglgq(ig, ifa),ie)*gpny*gpsa*bg(1:ndegr, ig)*wi
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fsqfg(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*wi
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fsqfg(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*wi
!
elnpn(1:ndegr)=elnpn(1:ndegr) +&
ufgpq(1, fglgq(ig, ifa), ie)*&
fsqfg(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*wi +&
ufgpq(2, fglgq(ig, ifa), ie)*&
fsqfg(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*wi

!...Output for debugging
!if(ie.eq.1)print*,'rhsface-gp BF',ifa,ig,ie,fglgq(ig, ifa),bg(1:3, ig)
!if(ie.eq.1)print*,'rhsface-gp',ifa,ig,ie,fglgq(ig, ifa),plnpn(1, 1:3),fsqfg(1:2, fglgq(ig, ifa), ie), wi

enddo
!
!if(rkstg.le.2)then
!...Correction force
do ig   = nvfac+1, ngausf
!
wi  = weigh(2)/2.d0

rg = posiq(1,fglgq(ig, ifa))
sg = posiq(2,fglgq(ig, ifa))
!
bg(1,  ig) = 1.d0
bg(2,  ig) = (rg-rc)/dr
bg(3,  ig) = (sg-sc)/ds
!DGP2
if(npoly.eq.2)then
bg(4, ig) = 0.5d0*bg(2, ig)*bg(2, ig) - geoel(19, ielem)
bg(5, ig) = 0.5d0*bg(3, ig)*bg(3, ig) - geoel(20, ielem)
bg(6, ig) =       bg(2, ig)*bg(3, ig) - geoel(21, ielem)
endif

!...The vertex constituting one cell...
gpnx = vngf(1, ig)
gpny = vngf(2, ig)
gpsa = vngf(3, ig)
!
!if(ie.eq.99)print*,'rhsface-gp-nv',ifa,ig,ie,vngf(1:3,ig),ufgpq(1:2, fglgq(ig, ifa),ie)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fsqfgdu(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*wi
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fsqfgdu(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*wi
!
elnpn(1:ndegr)=elnpn(1:ndegr) +&
ufgpq(1, fglgq(ig, ifa), ie)*&
fsqfgdu(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*wi +&
ufgpq(2, fglgq(ig, ifa), ie)*&
fsqfgdu(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*wi

!...Output for debugging
!if(ie.eq.1)print*,'rhsface-gp3',ifa,ig,ie,plnpn(1, 1:3),fsqfgdu(1:2, fglgq(ig, ifa), ie), wi
!
enddo

!endif

enddo
!
!...Distribute to every corner...
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)

!...Output foe debugging
!if(ie.eq.594.or.ie.eq.593.or.ie.eq.575.or.ie.eq.576)  print*,'rhs iface',ielem, ie,plnpn(1,1:3)
650 enddo
!
end subroutine rhsifacedg_lagsubgq_gs
!
!...Face integral for sub-cell scheme with both vertex and non-vertex gauss point...
!
subroutine rhsifacemem_lagsms_g(ipqua, ustar, ufgpq, gesgq0, geoel, coord,rhsgd)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:2,1:ngsfq,nquad),         intent(in)::ufgpq
real*8,dimension(1:ndegr,1:4,1:ncell),        intent(out)::rhsgd
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq0
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord

!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem,isg,ivsg,ifsg
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac, 1:nvqua) :: ipfq
integer,dimension(8, 4)::fnqsg
integer,dimension(4, 4)::ipqsg
integer,dimension(nvfac,4)::fglvq
integer,dimension(ngausf,4)::fglgq
!
real*8::posiq(1:2, 1:ngsfq)
real*8,dimension(ngausf)::wfgsqv
real*8,dimension(1:3,1:2,1:4)::vnorm
real*8, dimension(1:ndegr, 1:4) :: ulnpn
real*8::xvq(nvqua), yvq(nvqua)
real*8::bg(ndegr,  ngausf)
real*8::xpf(2,nvfac)
real*8::vngf(3,ngausf)
real*8::posi(ngausf), weigh(ngausf)

!...local real number
real*8::eps,c00,c05,c10,c20,c13,c16
real*8::dr,ds,rci,sci, rg, sg, wi
real*8::gpnx,gpny,gpsa,sarea
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
data c16   / 0.1666666666666666d0 /
!
!...Part I: Specify some gauss pints
!
call ruqope_lobatto(1, ngausf, posi, weigh)
!...Weight for every Gauss point
wfgsqv(1) = 0.5d0*weigh(1)
wfgsqv(2) = 0.5d0*weigh(5)
wfgsqv(3) = 0.5d0*weigh(3)
wfgsqv(4) = 0.5d0*weigh(2)
wfgsqv(5) = 0.5d0*weigh(4)
!
ipfq(1, 1) = 1; ipfq(2, 1) = 2; ipfq(3, 1) = 5
ipfq(1, 2) = 2; ipfq(2, 2) = 3; ipfq(3, 2) = 6
ipfq(1, 3) = 3; ipfq(2, 3) = 4; ipfq(3, 3) = 7
ipfq(1, 4) = 4; ipfq(2, 4) = 1; ipfq(3, 4) = 8

!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) = posi(2); posiq(2,13) =-1.d0;
posiq(1,14) = posi(4); posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) = posi(2);
posiq(1,16) = 1.d0; posiq(2,16) = posi(4);

posiq(1,17) = posi(4); posiq(2,17) = 1.d0;
posiq(1,18) = posi(2); posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = posi(4);
posiq(1,20) =-1.d0; posiq(2,20) = posi(2);
!
!...Part II: Loop over the element
!
do 650 ie = 1,nquad
ielem = ie + ntria

!...The vertex constituting one cell...
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
dr = 1.d0
ds = 1.d0
!
rci= geoel(7, ielem) !...mass center...
sci= geoel(8, ielem)
!
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

!...Initialize ulnpn, plnpn, elnpn
ulnpn = 0.d0

!...Loop over every face
do ifa = 1, 4
do ig   = 1, nvfac
!
rg = xvq(fglvq(ig, ifa))
sg = yvq(fglvq(ig, ifa))
!
bg(1,  ig) = 1.d0
bg(2,  ig) = (rg-rci)/dr
bg(3,  ig) = (sg-sci)/ds

!DGP2
if(npoly.eq.2)then
bg(4,  ig) = 0.5d0*bg(2,  ig)*bg(2,  ig) - geoel(22, ielem)
bg(5,  ig) = 0.5d0*bg(3,  ig)*bg(3,  ig) - geoel(23, ielem)
bg(6,  ig) =       bg(2,  ig)*bg(3,  ig) - geoel(24, ielem)
endif
!
!...The vertex constituting one cell...
gpnx = gesgq0(1, fglgq(ig, ifa), ie)
gpny = gesgq0(2, fglgq(ig, ifa), ie)
gpsa = gesgq0(3, fglgq(ig, ifa), ie)
!
ulnpn(1:ndegr, 1)  = ulnpn(1:ndegr, 1)+&
ustar(1, ipq(fglvq(ig, ifa)))*gpnx*gpsa*bg(1:ndegr, ig)*wfgsqv(ig)!*lpnpq(1, 1:ndegr, 1, iv)
!
ulnpn(1:ndegr, 2)  = ulnpn(1:ndegr, 2)+&
ustar(1, ipq(fglvq(ig, ifa)))*gpny*gpsa*bg(1:ndegr, ig)*wfgsqv(ig)
!
ulnpn(1:ndegr, 3)  = ulnpn(1:ndegr, 3)+&
ustar(2, ipq(fglvq(ig, ifa)))*gpnx*gpsa*bg(1:ndegr, ig)*wfgsqv(ig)!*lpnpq(1, 1:ndegr, 1, iv)
!
ulnpn(1:ndegr, 4)  = ulnpn(1:ndegr, 4)+&
ustar(2, ipq(fglvq(ig, ifa)))*gpny*gpsa*bg(1:ndegr, ig)*wfgsqv(ig)
!
!if(ie==1) print*,ifa,ig,ielem,ipq(fglvq(ig, ifa)), gpnx,gpsa,bg(4, ig),weigh(ig),ulnpn(4, 1),ustar(1, ipq(fglvq(ig, ifa)))
!
enddo
enddo


!...Contribution from non-vertex gauss point
sarea = 0.d0
!
do ifa  = 1, 4
!
xpf(1, 1:nvfac) = coord(1, ipq(ipfq(1:nvfac, ifa)))
xpf(2, 1:nvfac) = coord(2, ipq(ipfq(1:nvfac, ifa)))

!...Get the face area
call getvectf(xpf, vngf)
!
do ig   = nvfac+1, ngausf
!
wi  = weigh(2)/2.d0
!
rg = posiq(1,fglgq(ig, ifa))
sg = posiq(2,fglgq(ig, ifa))
!
bg(1,  ig) = 1.d0
bg(2,  ig) = (rg-rci)/dr
bg(3,  ig) = (sg-sci)/ds
!DGP2
if(npoly.eq.2)then
bg(4, ig) = 0.5d0*bg(2, ig)*bg(2, ig) - geoel(22, ielem)
bg(5, ig) = 0.5d0*bg(3, ig)*bg(3, ig) - geoel(23, ielem)
bg(6, ig) =       bg(2, ig)*bg(3, ig) - geoel(24, ielem)
endif

!...The vertex constituting one cell...
gpnx = vngf(1, ig)
gpny = vngf(2, ig)
gpsa = vngf(3, ig)
!
sarea = sarea + gpnx*gpsa

!if(ie.eq.593)print*,'rhsface-gp-nv',ifa,ig,ie,fsqfg(1, fglgq(ig, ifa), ie),gpnx,gpny,gpsa,sarea,&
!fsqfg(1, fglgq(ig, ifa), ie)/gpsa/gpnx
!
ulnpn(1:ndegr,1)  = ulnpn(1:ndegr,1)+&
ufgpq(1, fglgq(ig, ifa),ie)*gpnx*gpsa*bg(1:ndegr, ig)*wi
!
ulnpn(1:ndegr,2)  = ulnpn(1:ndegr,2)+&
ufgpq(1, fglgq(ig, ifa),ie)*gpny*gpsa*bg(1:ndegr, ig)*wi
!
ulnpn(1:ndegr,3)  = ulnpn(1:ndegr,3)+&
ufgpq(2, fglgq(ig, ifa),ie)*gpnx*gpsa*bg(1:ndegr, ig)*wi
!
ulnpn(1:ndegr,4)  = ulnpn(1:ndegr,4)+&
ufgpq(2, fglgq(ig, ifa),ie)*gpny*gpsa*bg(1:ndegr, ig)*wi

!...Output for debugging
!if(ie.eq.1)print*,'rhsface-gp BF',ifa,ig,ie,fglgq(ig, ifa),bg(1:3, ig)
!if(ie.eq.1)print*,'rhsface-gp',ifa,ig,ie,fglgq(ig, ifa),plnpn(1, 1:3),fsqfg(1:2, fglgq(ig, ifa), ie), wi

enddo

!endif

enddo
!
rhsgd(1:ndegr, 1, ielem) =  ulnpn(1:ndegr, 1)
rhsgd(1:ndegr, 2, ielem) =  ulnpn(1:ndegr, 2)
rhsgd(1:ndegr, 3, ielem) =  ulnpn(1:ndegr, 3)
rhsgd(1:ndegr, 4, ielem) =  ulnpn(1:ndegr, 4)

!...Output foe debugging
!if(ie.eq.594.or.ie.eq.593.or.ie.eq.575.or.ie.eq.576)  print*,'rhs iface',ielem, ie,plnpn(1,1:3)
650 enddo
!
end subroutine rhsifacemem_lagsms_g
!
!...Get the left and right state at one face...
!
subroutine getvar_face_general(xpqi,xpq,ipqua, ipf, fgaus, unkng, sigmg, unkno, drhosgq, geoel, afvec, &
                              aflim ,ielem)
use constant
implicit none
!...Input integer arrays
integer,dimension(1:nvqua),                intent(in)::ipqua
integer,dimension(1:nvfac),                intent(in):: ipf
integer,dimension(ngausf),                intent(out)::fgaus
!...Input real arrays
real*8,dimension(1:ndimn, 1:nvqua),        intent(in) :: xpq, xpqi
real*8, dimension(1:nq+3, 1:ngausf),       intent(out)::unkng
real*8, dimension(1:2, 1:2, 1:ngausf),     intent(out)::sigmg
real*8,dimension(1:ndegr,1:nq),            intent(in)::unkno
real*8,dimension(1:4),                     intent(in)::drhosgq
real*8,dimension(1:ngeel),                 intent(in)::geoel
real*8,dimension(1:nq+1),                  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2),                intent(in)::afvec
!
integer,                                   intent(in)::ielem

!...Local integer array
integer,dimension(ngausf,4)    ::fglgq
integer,dimension(20)    ::sglgq

!...Local real array
real*8, dimension(1:2, 1:ngsfq)::posiq
real*8, dimension(1:mdegr,1:ngausf)::bq

!...Local integer
integer::ig,ideg

!...Local real number
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr, eintc
real*8::rhovt,uvtx,vvtx,evtx,pvtx,eintv,sdv
real*8::rhoi
real*8::rhomv,rhomc
real*8::dudr,duds,dvdr,dvds
real*8::dr, ds,rc,sc
real*8::xcrho,ycrho,xgaus,ygaus
real*8::eps
!
eps = 1.d-15
!
rc = geoel(1)
sc = geoel(2)

dr = 1.d0
ds = 1.d0

!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;
!
sglgq(13) = 1; sglgq(14) = 2;
sglgq(15) = 2; sglgq(16) = 3;
sglgq(17) = 3; sglgq(18) = 4;
sglgq(19) = 4; sglgq(20) = 1;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) =-0.6546536707079771437983d0; posiq(2,13) =-1.d0;
posiq(1,14) = 0.6546536707079771437983d0; posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) =-0.6546536707079771437983d0;
posiq(1,16) = 1.d0; posiq(2,16) = 0.6546536707079771437983d0;

posiq(1,17) = 0.6546536707079771437983d0; posiq(2,17) = 1.d0;
posiq(1,18) =-0.6546536707079771437983d0; posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = 0.6546536707079771437983d0;
posiq(1,20) =-1.d0; posiq(2,20) =-0.6546536707079771437983d0;

!...Mass average value
rhoct = 1.d0/unkno(1, 1)
uctr  = unkno(1, 2)
vctr  = unkno(1, 3)
ectr  = unkno(1, 4)
!
eintc  = ectr - 0.5d0*(uctr**2 + vctr**2)

!...Get the cenetr pressure and sound speed
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)
call getrhog_initial(rhoi,  xcrho, ycrho, xcrho, ycrho)
!
call GetEOS(nmatel,ncase,gamlg,rhoct, eintc, rhoi, pctr, sdctr, ielem)

!pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )

!...Identify the fgaus for right cell
call getfg_glb(ipqua, ipf, fgaus)
!
do ig =nvfac+1 ,ngausf
bq(1, ig) = 1.d0
bq(2, ig) = (posiq(1, fgaus(ig))-rc)/dr
bq(3, ig) = (posiq(2, fgaus(ig))-sc)/ds
!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21)
endif
enddo

!...Zero out unkng
unkng = 0.d0

do ig = nvfac+1 , ngausf
do ideg =1, mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq)*bq(ideg, ig)
enddo
enddo

!...Limiter...
do ig = nvfac+1, ngausf
!
rhovt = 1.d0/unkng(1, ig) + cdrho*drhosgq(sglgq(fgaus(ig)))
uvtx = unkng(2, ig)
vvtx = unkng(3, ig)
evtx = unkng(4, ig)

!...Derived variables at the vertex
eintv = evtx - 0.5d0*(uvtx**2 + vvtx**2)

!...Call EOS
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, posiq(1, fgaus(ig)), posiq(2, fgaus(ig)), xgaus, ygaus)
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)

call GetEOS(nmatel,ncase,gamlg,rhovt, eintv, rhoi, pvtx, sdv, ielem)
!
!pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6.and.geoel(10).gt.10.d0)then
!
rhomc = 1.d0/rhoct
!
!rhomv = rhomc + aflim(1)*(unkng(1, ig) - rhomc)
rhomv = rhomc + aflim(1)*(1.d0/rhovt - rhomc)
unkng(1, ig) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1)*unkno(2,2) +  afvec(1, 2)*unkno(2,3)
duds = afvec(1, 1)*unkno(3,2) +  afvec(1, 2)*unkno(3,3)
dvdr = afvec(2, 1)*unkno(2,2) +  afvec(2, 2)*unkno(2,3)
dvds = afvec(2, 1)*unkno(3,2) +  afvec(2, 2)*unkno(3,3)
!
uvtx = unkno(1,2)  + dudr*bq(2, ig) + duds*bq(3, ig)
vvtx = unkno(1,3)  + dvdr*bq(2, ig) + dvds*bq(3, ig)
!
pvtx = pctr + aflim(4)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig) = uvtx
unkng(3, ig) = vvtx
!
endif
!
unkng(5, ig) = pvtx
unkng(6, ig) = rhoct*sdctr
unkng(7, ig) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig) = -pvtx
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pvtx!
enddo

end subroutine getvar_face_general
!
!...Get the left and right stress state at one face...
!
subroutine getvarsolid_face_general(xpqi, xpq, ipqua, ipf, fgaus, sigma_devt, unkgd, strnq_devtp,&
geoel, afvec, &
aflim ,ielem)
use constant
implicit none
!...Input integer arrays
integer,dimension(1:nvqua),                intent(in)::ipqua
integer,dimension(1:nvfac),                intent(in):: ipf
integer,dimension(ngausf),                intent(out)::fgaus
!...Input real arrays
real*8,dimension(1:ndimn, 1:nvqua),        intent(in) :: xpq, xpqi
real*8, dimension(1:2, 1:2, 1:ngausf),     intent(out)::sigma_devt
real*8,dimension(1:ndegr,1:4),             intent(in)::unkgd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq)::strnq_devtp
real*8,dimension(1:ngeel),                 intent(in)::geoel
real*8,dimension(1:nq+1),                  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2),                intent(in)::afvec
!
integer,                                   intent(in)::ielem

!...Local integer array
integer,dimension(20)    ::fgmap_q
integer,dimension(20)    ::sglgq

!...Local real array
real*8, dimension(1:2, 1:ngsfq)::posiq

!...Local integer
integer::ig

!...Local real number
real*8::rc,sc,rci,sci,rg,sg
real*8::eps
!
eps = 1.d-15
!
rc = geoel(1)
sc = geoel(2)
!...The initial cell center
rci= geoel(7)
sci= geoel(8)

!...Quadrature points mapping
fgmap_q = 0
!
fgmap_q(13) =10; fgmap_q(14) = 11
fgmap_q(15) =12; fgmap_q(16) = 13
fgmap_q(17) =14; fgmap_q(18) = 15
fgmap_q(19) =16; fgmap_q(20) = 17

!
sglgq(13) = 1; sglgq(14) = 2;
sglgq(15) = 2; sglgq(16) = 3;
sglgq(17) = 3; sglgq(18) = 4;
sglgq(19) = 4; sglgq(20) = 1;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) =-0.6546536707079771437983d0; posiq(2,13) =-1.d0;
posiq(1,14) = 0.6546536707079771437983d0; posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) =-0.6546536707079771437983d0;
posiq(1,16) = 1.d0; posiq(2,16) = 0.6546536707079771437983d0;

posiq(1,17) = 0.6546536707079771437983d0; posiq(2,17) = 1.d0;
posiq(1,18) =-0.6546536707079771437983d0; posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = 0.6546536707079771437983d0;
posiq(1,20) =-1.d0; posiq(2,20) =-0.6546536707079771437983d0;

!...Identify the fgaus for right cell
call getfg_glb(ipqua, ipf, fgaus)

!...Loop over the Quadrature points
do ig = nvfac+1, ngausf
!
rg = posiq(1, fgaus(ig))
sg = posiq(2, fgaus(ig))

!...Call stress
call GetStress_deviat_infsmal(rg,sg,rci,sci,geoel(22:24),sigma_devt(:,:,ig), strnq_devtp(:,:, fgmap_q(fgaus(ig))), unkgd, ielem, 1)
!
enddo
!
!if(ielem.eq.1) print*,'sigma',sigma_devt(:,:,4)
end subroutine getvarsolid_face_general

!
!...Get the left and right state at one face...
!
subroutine getvar_face(ipqua, ipf, fgaus, unkng, sigmg, unkno, drhosgq, geoel, afvec, aflim)
use constant
implicit none
!...Input integer arrays
integer,dimension(1:nvqua),                intent(in)::ipqua
integer,dimension(1:nvfac),                intent(in):: ipf
integer,dimension(ngausf),                intent(out)::fgaus
!...Input real arrays
real*8, dimension(1:nq+3, 1:ngausf),       intent(out)::unkng
real*8, dimension(1:2, 1:2, 1:ngausf),     intent(out)::sigmg
real*8,dimension(1:ndegr,1:nq),            intent(in)::unkno
real*8,dimension(1:4),                     intent(in)::drhosgq
real*8,dimension(1:ngeel),                 intent(in)::geoel
real*8,dimension(1:nq+1),                  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2),                intent(in)::afvec

!...Local integer array
integer,dimension(ngausf,4)    ::fglgq
integer,dimension(20)    ::sglgq

!...Local real array
real*8, dimension(1:2, 1:ngsfq)::posiq
real*8, dimension(1:mdegr,1:ngausf)::bq

!...Local integer
integer::ig,ideg

!...Local real number
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::rhovt,uvtx,vvtx,evtx,pvtx
real*8::rhomv,rhomc
real*8::dudr,duds,dvdr,dvds
real*8::dr, ds,rc,sc
real*8::eps
!
eps = 1.d-15
!
rc = geoel(1)
sc = geoel(2)

dr = 1.d0
ds = 1.d0

!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;
!
sglgq(13) = 1; sglgq(14) = 2;
sglgq(15) = 2; sglgq(16) = 3;
sglgq(17) = 3; sglgq(18) = 4;
sglgq(19) = 4; sglgq(20) = 1;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) =-0.6546536707079771437983d0; posiq(2,13) =-1.d0;
posiq(1,14) = 0.6546536707079771437983d0; posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) =-0.6546536707079771437983d0;
posiq(1,16) = 1.d0; posiq(2,16) = 0.6546536707079771437983d0;

posiq(1,17) = 0.6546536707079771437983d0; posiq(2,17) = 1.d0;
posiq(1,18) =-0.6546536707079771437983d0; posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = 0.6546536707079771437983d0;
posiq(1,20) =-1.d0; posiq(2,20) =-0.6546536707079771437983d0;

!...Mass average value
rhoct = 1.d0/unkno(1, 1)
uctr  = unkno(1, 2)
vctr  = unkno(1, 3)
ectr  = unkno(1, 4)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )

!...Identify the fgaus for right cell
call getfg_glb(ipqua, ipf, fgaus)
!
do ig =nvfac+1 ,ngausf
bq(1, ig) = 1.d0
bq(2, ig) = (posiq(1, fgaus(ig))-rc)/dr
bq(3, ig) = (posiq(2, fgaus(ig))-sc)/ds
!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21)
endif
enddo

!...Zero out unkng
unkng = 0.d0

do ig = nvfac+1 , ngausf
do ideg =1, mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq)*bq(ideg, ig)
enddo
enddo

!...Limiter...
do ig = nvfac+1, ngausf
!
rhovt = 1.d0/unkng(1, ig) + cdrho*drhosgq(sglgq(fgaus(ig)))
uvtx = unkng(2, ig)
vvtx = unkng(3, ig)
evtx = unkng(4, ig)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6.and.geoel(10).gt.10.d0)then
!
rhomc = 1.d0/rhoct
!
!rhomv = rhomc + aflim(1)*(unkng(1, ig) - rhomc)
rhomv = rhomc + aflim(1)*(1.d0/rhovt - rhomc)
unkng(1, ig) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1)*unkno(2,2) +  afvec(1, 2)*unkno(2,3)
duds = afvec(1, 1)*unkno(3,2) +  afvec(1, 2)*unkno(3,3)
dvdr = afvec(2, 1)*unkno(2,2) +  afvec(2, 2)*unkno(2,3)
dvds = afvec(2, 1)*unkno(3,2) +  afvec(2, 2)*unkno(3,3)
!
uvtx = unkno(1,2)  + dudr*bq(2, ig) + duds*bq(3, ig)
vvtx = unkno(1,3)  + dvdr*bq(2, ig) + dvds*bq(3, ig)
!
pvtx = pctr + aflim(4)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig) = uvtx
unkng(3, ig) = vvtx
!
endif
!
unkng(5, ig) = pvtx
unkng(6, ig) = rhoct*sdctr
unkng(7, ig) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig) = -pvtx
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pvtx!
enddo

end subroutine getvar_face
!
!...Get the left and right state at one face...
!
subroutine getvar_facesmth(ipqua, ipf, fgaus, unkng, sigmg, unkno, drhosgq, geoel, afvec, aflim)
use constant
implicit none
!...Input integer arrays
integer,dimension(1:nvqua),                intent(in)::ipqua
integer,dimension(1:nvfac),                intent(in):: ipf
integer,dimension(ngausf),                intent(out)::fgaus
!...Input real arrays
real*8, dimension(1:nq+3, 1:ngausf),       intent(out)::unkng
real*8, dimension(1:2, 1:2, 1:ngausf),     intent(out)::sigmg
real*8,dimension(1:ndegr,1:nq),            intent(in)::unkno
real*8,dimension(1:4),                     intent(in)::drhosgq
real*8,dimension(1:ngeel),                 intent(in)::geoel
real*8,dimension(1:nq+1),                  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2),                intent(in)::afvec

!...Local integer array
integer,dimension(ngausf,4)    ::fglgq
integer,dimension(20)    ::sglgq

!...Local real array
real*8, dimension(1:2, 1:ngsfq)::posiq
real*8, dimension(1:mdegr,1:ngausf)::bq

!...Local integer
integer::ig,ideg

!...Local real number
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::rhovt,uvtx,vvtx,evtx,pvtx
real*8::duvtx,dvvtx
real*8::rhomv,rhomc
real*8::dudr,duds,dvdr,dvds
real*8::dr, ds,rc,sc
real*8::eps
!
eps = 1.d-6
!
rc = geoel(1)
sc = geoel(2)

dr = 1.d0
ds = 1.d0

!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;
!
sglgq(13) = 1; sglgq(14) = 2;
sglgq(15) = 2; sglgq(16) = 3;
sglgq(17) = 3; sglgq(18) = 4;
sglgq(19) = 4; sglgq(20) = 1;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) =-0.6546536707079771437983d0; posiq(2,13) =-1.d0;
posiq(1,14) = 0.6546536707079771437983d0; posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) =-0.6546536707079771437983d0;
posiq(1,16) = 1.d0; posiq(2,16) = 0.6546536707079771437983d0;

posiq(1,17) = 0.6546536707079771437983d0; posiq(2,17) = 1.d0;
posiq(1,18) =-0.6546536707079771437983d0; posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = 0.6546536707079771437983d0;
posiq(1,20) =-1.d0; posiq(2,20) =-0.6546536707079771437983d0;

!...Mass average value
rhoct = 1.d0/unkno(1, 1)
uctr  = unkno(1, 2)
vctr  = unkno(1, 3)
ectr  = unkno(1, 4)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )

!...Identify the fgaus for right cell
call getfg_glb(ipqua, ipf, fgaus)
!
do ig =nvfac+1 ,ngausf
bq(1, ig) = 1.d0
bq(2, ig) = (posiq(1, fgaus(ig))-rc)/dr
bq(3, ig) = (posiq(2, fgaus(ig))-sc)/ds
!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21)
endif
enddo

!...Zero out unkng
unkng = 0.d0

do ig = nvfac+1 , ngausf
do ideg =1, mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq)*bq(ideg, ig)
enddo
enddo

!...Limiter...
do ig = nvfac+1, ngausf
!
rhovt = 1.d0/unkng(1, ig) + .8d0*drhosgq(sglgq(fgaus(ig)))
uvtx = unkng(2, ig)
vvtx = unkng(3, ig)
evtx = unkng(4, ig)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6)then
!
rhomc = 1.d0/rhoct

!...Shock cell
if(geoel(10).gt.15.d0)then
!
rhomv = rhomc + aflim(1)*(1.d0/rhovt - rhomc)
unkng(1, ig) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1)*unkno(2,2) +  afvec(1, 2)*unkno(2,3)
duds = afvec(1, 1)*unkno(3,2) +  afvec(1, 2)*unkno(3,3)
dvdr = afvec(2, 1)*unkno(2,2) +  afvec(2, 2)*unkno(2,3)
dvds = afvec(2, 1)*unkno(3,2) +  afvec(2, 2)*unkno(3,3)
!
uvtx = unkno(1,2)  + dudr*bq(2, ig) + duds*bq(3, ig)
vvtx = unkno(1,3)  + dvdr*bq(2, ig) + dvds*bq(3, ig)
!
pvtx = pctr + aflim(4)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig) = uvtx
unkng(3, ig) = vvtx

!...Transition cell
elseif(geoel(10).lt.15.d0.and.geoel(10).gt.5.d0)then
!
!rhomv = rhomc + aflim(1)*(unkng(1, ig) - rhomc)
rhomv = rhomc + aflim(1)*(1.d0/rhovt - rhomc)
unkng(1, ig) = rhomv
rhovt = 1.d0/rhomv
!
duvtx = uvtx - unkno(1,2)
dvvtx = vvtx - unkno(1,3)
!
uvtx = unkno(1,2)  + afvec(1, 1)*duvtx + afvec(1, 2)*dvvtx
vvtx = unkno(1,3)  + afvec(2, 1)*duvtx + afvec(2, 2)*dvvtx
!
pvtx = pctr + aflim(4)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig) = uvtx
unkng(3, ig) = vvtx
!
endif

endif
!
unkng(5, ig) = pvtx
unkng(6, ig) = rhoct*sdctr
unkng(7, ig) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig) = -pvtx
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pvtx!
enddo

end subroutine getvar_facesmth

!
!...Calculate the velocity at the Gauss point for novertex face(version 2)...
!
subroutine getvelo_quadsubg_gs2(ufgpq,dufgq,geoel,bface,intfac,ipqua,coord,unkno,indnd, &
munaclq, munaulq, snsigmlq, drhosgq, afvec, aflim, itime)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
integer*4,dimension(1:npoin),                intent(in)::indnd
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngsfq,nquad),      intent(in)::ufgpq
real*8,dimension(1:2,1:ngsfq,nquad),      intent(inout)::dufgq
real*8, dimension(1:2, 1:2, 1:ngsfq, 1:nquad), intent(inout)::munaclq
real*8, dimension(1:ndimn, 1:ngsfq,  1:nquad), intent(inout)::munaulq
real*8, dimension(1:ndimn, 1:ngsfq,  1:nquad), intent(inout)::snsigmlq
real*8,dimension(1:4, 1:nsize), intent(in)::drhosgq
integer*4,          intent(in)::itime
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(ngausf)::fgausl, fgausr
!
real*8::eps
real*8,dimension(1:2,1:ngelgq,nquad)::vlave
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::sigmg(1:2, 1:2, 1:ngausf, 1:2)
real*8::vngs(3, ngausf,2)
real*8::vngf(3, ngausf)
real*8::aujmp(1:3,1:ngausf,1:2)
real*8::murie(2)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, ngausf)
real*8::posiq(1:2, 1:ngsfq)
integer,dimension(nvfac, 4)::fglvq
integer,dimension(ngausf,4)::fglgq
real*8::xpf(1:2, 1:nvfac)
real*8::fpres(1:2, 1:ngausf)

!...Riemann parameters...
real*8::munacn_rie(2,2), munacu_rie(2), snsigm_rie(2)
real*8::munacn(1:2, 1:2, 1:ngausf), munacu(1:2, 1:ngausf), snsigm(1:2, 1:ngausf)
real*8::munaci(2, 2)

!...Local real number
real*8::detma,dudr,duds,dvdr,dvds
real*8::pvtx,rhovt, rhomc, rhomv, rhovl, rhovr,rhsu1,rhsu2
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
real*8::othog
real*8::r,delu
real*8::acnx,acny

!...Initial velocity at gauss point
!..Initial Gauss point velocity from the interpolation of the 3-node velocity...
!
eps = 1.d-6
dufgq = 0.d0
!
dr = 1.d0
ds = 1.d0

!....xxxxxx need modified
rc = 0.0
sc = rc
!
!...Part I: Specify some gauss points
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
xvq(9) =  0.d0; yvq(9) =  0.d0

!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

!...8-nodes quad
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.4472135954999579392818d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.4472135954999579392818d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.4472135954999579392818d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.4472135954999579392818d0;

posiq(1,13) = 0.4472135954999579392818d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.4472135954999579392818d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.4472135954999579392818d0;
posiq(1,16) =-1.d0; posiq(2,16) =-0.4472135954999579392818d0;

elseif(ngausf.eq.5)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) =-0.6546536707079771437983d0; posiq(2,13) =-1.d0;
posiq(1,14) = 0.6546536707079771437983d0; posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) =-0.6546536707079771437983d0;
posiq(1,16) = 1.d0; posiq(2,16) = 0.6546536707079771437983d0;

posiq(1,17) = 0.6546536707079771437983d0; posiq(2,17) = 1.d0;
posiq(1,18) =-0.6546536707079771437983d0; posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = 0.6546536707079771437983d0;
posiq(1,20) =-1.d0; posiq(2,20) =-0.6546536707079771437983d0;
endif
!
!...Part II: Get the corrected velocity at the non-vertex gauss point...
!
fgausl = 0.d0
fgausr = 0.d0
!
do 450 ifa = 1, nafac !...(1)ifa = 1, nafac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))

!...For the linear PP+
ftx = xpf(1 ,2)- xpf(1, 1)
fty = xpf(2, 2)- xpf(2, 1)
!
fnx = -fty/sqrt(ftx**2 + fty**2)
fny =  ftx/sqrt(ftx**2 + fty**2)

!...For linear PM
ftx = xpf(1 ,3)- xpf(1, 1)
fty = xpf(2, 3)- xpf(2, 1)
!
lenmc = sqrt(ftx**2 + fty**2)
!
ftx = ftx/lenmc
fty = fty/lenmc
!
othog = abs(fnx*ftx + fny*fty)

!...Get the normal vector at the non-vertex gauss point
call getvectf(xpf, vngf)

!...Boundary face
if(ifa.le.nbfac)then

!...The normal vector at the non-vertex gauss point
do ig = nvfac+1, ngausf
vngs(1:2, ig, 1) = vngf(1:2, ig)
vngs(  3, ig, 1) = vngf(  3, ig)
enddo

!...Get the variables at the left state
if(nfint.eq.5.or.nfint.eq.10)then

call getvar_face(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))

elseif(nfint.eq.7)then

call getvar_facesmth(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))
endif

!...Impedence...
do ig = nvfac+1, ngausf
do ic = 1, 1
!
!aujmp(1:2, ig, ic) = vlave(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgpq(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)
!aujmp(3, ig, ic) = 1.d0*abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
!
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = cimpd*abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) + unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Boundary condition imposing
call getbcgauss_lagc(bface(3, ifa), vngf, fpres, itime)

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0

!...Excluding the non-vertex gs point
do ig = nvfac+1, ngausf
do ic = 1, 1
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
!
!...Output for debugging
!if(iel.eq.1)print*,'iel',ig,ier,fgausl(ig),ie,munacn_rie(1:2, 1),aujmp(1:3, ig, 1),vngs(1:2, ig, 1)
!
enddo

!...Output for debugging
!if(iel.eq.1)print*,'iel',ig,ier,fgausl(ig),ie,unkng(6,ig,1)*vngs(3, ig, 1), ufgpq(1:2, fgausl(ig), iel),unkng(2:3, ig, 1),&
!unkng(5,ig,1),vngs(3, ig, 1),unkng(6, ig, 1),aujmp(1:3, ig, 1)

!...Get
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig) - fpres(1, ig)
rhsu2 = munacu(2, ig) - snsigm(2, ig) - fpres(2, ig)
!
dufgq(1, fgausl(ig), iel) =-ufgpq(1, fgausl(ig), iel) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
dufgq(2, fgausl(ig), iel) =-ufgpq(2, fgausl(ig), iel) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)

!...Impose the boundary condition to adjust the ustar at the boundary
if(bface(3 ,ifa).eq.22)then
if(bface(4 ,ifa).eq.221)then
dufgq(2, fgausl(ig), iel)= 0.d0
elseif(bface(4 ,ifa).eq.222)then
dufgq(1, fgausl(ig), iel)= 0.d0
endif

elseif(bface(3 ,ifa).eq.25)then
dufgq(1:2, fgausl(ig), iel)= 0.d0
elseif(bface(3 ,ifa).eq.27)then
dufgq(1:2, fgausl(ig), iel)= 0.d0
elseif(bface(3 ,ifa).eq.21)then
!dufgq(1:2, fgausl(ig), iel)= 0.d0
elseif(bface(3 ,ifa).eq.23)then
!dufgq(1:2, fgausl(ig), iel)= 0.d0
endif
!
enddo

!...Interior face
elseif(ifa.gt.nbfac)then

!...The normal vector at the non-vertex gauss point
do ig = nvfac+1, ngausf
vngs(1:2, ig, 1) = vngf(1:2, ig)
vngs(  3, ig, 1) = vngf(  3, ig)
enddo

!...Get the variables at the left state
if(nfint.eq.5.or.nfint.eq.10)then

call getvar_face(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))
elseif(nfint.eq.7)then

call getvar_facesmth(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))
endif
!...The normal vector at the non-vertex gauss point
do ig = nvfac+1 , ngausf
vngs(1:2, ig, 2) =-vngf(1:2, ig)
vngs(  3, ig, 2) = vngf(  3, ig)
enddo
!
!if(iel.eq.1)print*,'ielgs',iel,ifa,unkno(1:3,:,iel),unkng(1,4,1),unkng(4,4,1)-0.5d0*(unkng(2,4,1)**2 + unkng(3,4,1)**2),&
!unkng(1,5,1),unkng(4,5,1)-0.5d0*(unkng(2,5,1)**2 + unkng(3,5,1)**2)


!...Get the variables at the right state
if(nfint.eq.5.or.nfint.eq.10)then
call getvar_face(ipqua(1:nvqua, ier), ipf, fgausr, unkng(:,:,2), sigmg(:,:,:,2), &
unkno(:,:,ier), drhosgq(:, ier),geoel(:,ier), afvec(:,:,ier), aflim(:,ier))
elseif(nfint.eq.7)then
call getvar_facesmth(ipqua(1:nvqua, ier), ipf, fgausr, unkng(:,:,2), sigmg(:,:,:,2), &
unkno(:,:,ier), drhosgq(:, ier),geoel(:,ier), afvec(:,:,ier), aflim(:,ier))
endif

!...Get the averaged velocity at the gauss quadrature points...
do ig = nvfac+1, ngausf
!vlave(1, fgausl(ig), iel) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
!vlave(2, fgausl(ig), iel) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
!
!vlave(1, fgausr(ig), ier) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
!vlave(2, fgausr(ig), ier) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
enddo

!...Impedence...
do ig = nvfac+1, ngausf
do ic = 1, 2

!aujmp(1:2, ig, ic) = vlave(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgpq(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)

!...Impedence
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = cimpd*abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) + unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = nvfac+1, ngausf
do ic = 1, 2
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
if(ic.eq.1)then
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
else
!
munaclq(1:2, 1, fgausr(ig), ier) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausr(ig), ier) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausr(ig), ier)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausr(ig), ier)    = snsigm_rie(1:2)
endif
!
enddo

!...Get the corrected velocity at the non-vertex gauss point
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig)
rhsu2 = munacu(2, ig) - snsigm(2, ig)
!
dufgq(1, fgausl(ig), iel) =-ufgpq(1, fgausl(ig), iel) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
dufgq(2, fgausl(ig), iel) =-ufgpq(2, fgausl(ig), iel) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)
!
dufgq(1, fgausr(ig), ier) =-ufgpq(1, fgausr(ig), ier) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
dufgq(2, fgausr(ig), ier) =-ufgpq(2, fgausr(ig), ier) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)

!...Output for debugging

!if(iel.eq.1)print*,'iel',ig,ier,fgausl(ig),ie,unkng(6,ig,1)*vngs(3, ig, 1), ufgpq(1:2, fgausl(ig), iel),unkng(2:3, ig, 1),&
!unkng(5,ig,1),vngs(3, ig, 1),munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2,munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2,&
!snsigm(1:2, ig),munacu(1:2, ig)
!if(iel.eq.1)print*,'ier',ig,ier,fgausr(ig),ie,unkng(6,ig,2)*vngs(3, ig, 2), ufgpq(1:2, fgausl(ig), iel),unkng(2:3, ig, 2),&
!unkng(5,ig,2),vngs(3, ig, 2),snsigm(1:2, ig),munacu(1:2, ig)
!
enddo

endif
450 enddo
!
!print*,'dufg2',dufgq(1, 17:18, 99)
!
end subroutine getvelo_quadsubg_gs2
!
!...Calculate the velocity at the Gauss point for novertex face(version 2)...
!
subroutine getvelo_quadsubg_gs3(ufgpq,dufgq,geoel,bface,intfac,ipqua,coord,coold,unkno,unkgd,strnq_devtp,indnd, &
munaclq, munaulq, snsigmlq, drhosgq, afvec, aflim, itime)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq, nquad)::strnq_devtp
integer*4,dimension(1:npoin),                intent(in)::indnd
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngsfq,nquad),      intent(in)::ufgpq
real*8,dimension(1:2,1:ngsfq,nquad),      intent(inout)::dufgq
real*8, dimension(1:2, 1:2, 1:ngsfq, 1:nquad), intent(inout)::munaclq
real*8, dimension(1:ndimn, 1:ngsfq,  1:nquad), intent(inout)::munaulq
real*8, dimension(1:ndimn, 1:ngsfq,  1:nquad), intent(inout)::snsigmlq
real*8,dimension(1:4, 1:nsize), intent(in)::drhosgq
integer*4,          intent(in)::itime
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf,ipfp
integer,dimension(ngausf)::fgausl, fgausr
!
real*8::eps
real*8,dimension(1:2,1:ngelgq,nquad)::vlave
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::sigmg(1:2, 1:2, 1:ngausf, 1:2)
real*8::sigma_devt(1:2, 1:2, 1:ngausf,1:2)
real*8::vngs(3, ngausf,2)
real*8::vngf(3, ngausf)
real*8::aujmp(1:3,1:ngausf,1:2)
real*8::murie(2)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, ngausf)
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xpqi
real*8::posiq(1:2, 1:ngsfq)
integer,dimension(nvfac, 4)::fglvq
integer,dimension(ngausf,4)::fglgq
real*8::xpf(1:2, 1:nvfac)
real*8::fpres(1:2, 1:ngausf)

!...Riemann parameters...
real*8::munacn_rie(2,2), munacu_rie(2), snsigm_rie(2)
real*8::munacn(1:2, 1:2, 1:ngausf), munacu(1:2, 1:ngausf), snsigm(1:2, 1:ngausf)
real*8::munaci(2, 2)

!...Local real number
real*8::detma,dudr,duds,dvdr,dvds
real*8::pvtx,rhovt, rhomc, rhomv, rhovl, rhovr,rhsu1,rhsu2
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
real*8::othog
real*8::r,delu
real*8::acnx,acny
real*8::xcrho, ycrho

!...Initial velocity at gauss point
!..Initial Gauss point velocity from the interpolation of the 3-node velocity...
!
eps = 1.d-6
dufgq = 0.d0
!
dr = 1.d0
ds = 1.d0

!....xxxxxx need modified
rc = 0.0
sc = rc
!
!...Part I: Specify some gauss points
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
xvq(9) =  0.d0; yvq(9) =  0.d0

!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

!...8-nodes quad
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.4472135954999579392818d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.4472135954999579392818d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.4472135954999579392818d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.4472135954999579392818d0;

posiq(1,13) = 0.4472135954999579392818d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.4472135954999579392818d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.4472135954999579392818d0;
posiq(1,16) =-1.d0; posiq(2,16) =-0.4472135954999579392818d0;

elseif(ngausf.eq.5)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) =-0.6546536707079771437983d0; posiq(2,13) =-1.d0;
posiq(1,14) = 0.6546536707079771437983d0; posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) =-0.6546536707079771437983d0;
posiq(1,16) = 1.d0; posiq(2,16) = 0.6546536707079771437983d0;

posiq(1,17) = 0.6546536707079771437983d0; posiq(2,17) = 1.d0;
posiq(1,18) =-0.6546536707079771437983d0; posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = 0.6546536707079771437983d0;
posiq(1,20) =-1.d0; posiq(2,20) =-0.6546536707079771437983d0;
endif
!
!...Part II: Get the corrected velocity at the non-vertex gauss point...
!
fgausl = 0.d0
fgausr = 0.d0
!
do 450 ifa = 1, nafac !...(1)ifa = 1, nafac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))

!...For the linear PP+
ftx = xpf(1 ,2)- xpf(1, 1)
fty = xpf(2, 2)- xpf(2, 1)
!
fnx = -fty/sqrt(ftx**2 + fty**2)
fny =  ftx/sqrt(ftx**2 + fty**2)

!...For linear PM
ftx = xpf(1 ,3)- xpf(1, 1)
fty = xpf(2, 3)- xpf(2, 1)
!
lenmc = sqrt(ftx**2 + fty**2)
!
ftx = ftx/lenmc
fty = fty/lenmc
!
othog = abs(fnx*ftx + fny*fty)

!...Get the normal vector at the non-vertex gauss point
call getvectf(xpf, vngf)

!...Boundary face
if(ifa.le.nbfac)then
!
ipq(1:nvqua) = ipqua(1:nvqua,iel)

!
 rc = geoel(1, iel)
 sc = geoel(2, iel)
!...Get density for sub-cells...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!...Get the intial cell center
!call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

if(abs(bface(3, ifa)).ne.31)then

!...The normal vector at the non-vertex gauss point
do ig = nvfac+1, ngausf
vngs(1:2, ig, 1) = vngf(1:2, ig)
vngs(  3, ig, 1) = vngf(  3, ig)
enddo

!...Get the variables at the left state
if(nfint.eq.5.or.nfint.eq.10)then

!...Get the hydrostatic pressure
call getvar_face_general(xpqi, xpq, ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel), iel)

if(nmatel.eq.2)then
call getvarsolid_face_general(xpqi, xpq, ipqua(1:nvqua, iel), ipf, fgausl, sigma_devt(:,:,:,1),&
                               unkgd(:,:,iel), strnq_devtp(:,:,:,iel),&
                               geoel(:,iel), afvec(:,:,iel), aflim(:,iel), iel)
!
sigmg(:,:,:,1) = sigmg(:,:,:,1) + sigma_devt(:,:,:,1)
endif
!call getvar_face_general(xpqi,xpq,ipqua, ipf, fgaus, unkng, sigmg, unkno, drhosgq, geoel, afvec, aflim ,ielem)

elseif(nfint.eq.7)then

call getvar_facesmth(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))
endif

!...Impedence...
do ig = nvfac+1, ngausf
do ic = 1, 1
!
!aujmp(1:2, ig, ic) = vlave(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgpq(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)
!aujmp(3, ig, ic) = 1.d0*abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
!
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = cimpd*abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) + unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Boundary condition imposing
call getbcgauss_lagc(bface(3, ifa), vngf, fpres, itime)

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0

!...Excluding the non-vertex gs point
do ig = nvfac+1, ngausf
do ic = 1, 1
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
!
!...Output for debugging
!if(iel.eq.1)print*,'iel',ig,ier,fgausl(ig),ie,munacn_rie(1:2, 1),aujmp(1:3, ig, 1),vngs(1:2, ig, 1)
!
enddo

!...Output for debugging
!if(iel.eq.1)print*,'iel',ig,ier,fgausl(ig),ie,unkng(6,ig,1)*vngs(3, ig, 1), ufgpq(1:2, fgausl(ig), iel),unkng(2:3, ig, 1),&
!unkng(5,ig,1),vngs(3, ig, 1),unkng(6, ig, 1),aujmp(1:3, ig, 1)

!...Get
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig) - fpres(1, ig)
rhsu2 = munacu(2, ig) - snsigm(2, ig) - fpres(2, ig)
!
dufgq(1, fgausl(ig), iel) =-ufgpq(1, fgausl(ig), iel) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
dufgq(2, fgausl(ig), iel) =-ufgpq(2, fgausl(ig), iel) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)

!...Impose the boundary condition to adjust the ustar at the boundary
if(bface(3 ,ifa).eq.22)then
if(bface(4 ,ifa).eq.221)then
dufgq(2, fgausl(ig), iel)= 0.d0
elseif(bface(4 ,ifa).eq.222)then
dufgq(1, fgausl(ig), iel)= 0.d0
endif

elseif(bface(3 ,ifa).eq.25)then
dufgq(1:2, fgausl(ig), iel)= 0.d0
elseif(bface(3 ,ifa).eq.27)then
dufgq(1:2, fgausl(ig), iel)= 0.d0
elseif(bface(3 ,ifa).eq.21)then
!dufgq(1:2, fgausl(ig), iel)= 0.d0
elseif(bface(3 ,ifa).eq.23)then
!dufgq(1:2, fgausl(ig), iel)= 0.d0
endif
!
enddo

!...Period Boundary face
elseif(abs(bface(3, ifa)).eq.31)then
!
!...The normal vector at the non-vertex gauss point
do ig = nvfac+1, ngausf
vngs(1:2, ig, 1) = vngf(1:2, ig)
vngs(  3, ig, 1) = vngf(  3, ig)
enddo

!...Get the variables at the left state
if(nfint.eq.5.or.nfint.eq.10)then

call getvar_face(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))
elseif(nfint.eq.7)then

call getvar_facesmth(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))
endif
!...The normal vector at the non-vertex gauss point
do ig = nvfac+1 , ngausf
vngs(1:2, ig, 2) =-vngf(1:2, ig)
vngs(  3, ig, 2) = vngf(  3, ig)
enddo
!
!if(iel.eq.1)print*,'ielgs',iel,ifa,unkno(1:3,:,iel),unkng(1,4,1),unkng(4,4,1)-0.5d0*(unkng(2,4,1)**2 + unkng(3,4,1)**2),&
!unkng(1,5,1),unkng(4,5,1)-0.5d0*(unkng(2,5,1)**2 + unkng(3,5,1)**2)

ipfp(1) = bface(2, bface(4, ifa));
ipfp(2) = bface(1, bface(4, ifa));
ipfp(3) = bface(5, bface(4, ifa));
!...Get the variables at the right state

if(nfint.eq.5.or.nfint.eq.10)then
call getvar_face(ipqua(1:nvqua, intfac(1,bface(4, ifa))), ipfp, fgausr, unkng(:,:,2), sigmg(:,:,:,2), &
unkno(:,:,intfac(1,bface(4, ifa))), drhosgq(:, intfac(1,bface(4, ifa))),geoel(:,intfac(1,bface(4, ifa))), &
afvec(:,:,intfac(1,bface(4, ifa))), aflim(:,intfac(1,bface(4, ifa))))
elseif(nfint.eq.7)then
call getvar_facesmth(ipqua(1:nvqua, ier), ipf, fgausr, unkng(:,:,2), sigmg(:,:,:,2), &
unkno(:,:,ier), drhosgq(:, ier),geoel(:,ier), afvec(:,:,ier), aflim(:,ier))
endif

!
!print*,'period',ifa,ipfp(1:3),intfac(1,bface(4, ifa)),fgausr

!...Impedence...
do ig = nvfac+1, ngausf
do ic = 1, 2

!aujmp(1:2, ig, ic) = vlave(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgpq(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)

!...Impedence
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = cimpd*abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) + unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = nvfac+1, ngausf
do ic = 1, 2
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!...Local
if(ic.eq.1)then
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
endif
!
enddo

!...Get the corrected velocity at the non-vertex gauss point
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig)
rhsu2 = munacu(2, ig) - snsigm(2, ig)
!
dufgq(1, fgausl(ig), iel) =-ufgpq(1, fgausl(ig), iel) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
dufgq(2, fgausl(ig), iel) =-ufgpq(2, fgausl(ig), iel) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)
!
enddo

endif

!...Interior face
elseif(ifa.gt.nbfac)then
!...Left cell
!
ipq(1:nvqua) = ipqua(1:nvqua,iel)
!...Get density for sub-cells...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!...The normal vector at the non-vertex gauss point
do ig = nvfac+1, ngausf
vngs(1:2, ig, 1) = vngf(1:2, ig)
vngs(  3, ig, 1) = vngf(  3, ig)
enddo

!...Get the variables at the left state
if(nfint.eq.5.or.nfint.eq.10)then

!...Get the hydrostatic pressure
call getvar_face_general(xpqi, xpq, ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel), iel)

!
!if(iel.eq.1)then
!print*,'iel1',sigmg(:,:,4:5,1),sigma_devt(:,:,4:5,1)
!endif

if(nmatel.eq.2)then
call getvarsolid_face_general(xpqi, xpq, ipqua(1:nvqua, iel), ipf, fgausl, sigma_devt(:,:,:,1),&
unkgd(:,:,iel), strnq_devtp(:,:,:,iel),&
geoel(:,iel), afvec(:,:,iel), aflim(:,iel), iel)
!
sigmg(:,:,:,1) = sigmg(:,:,:,1) + sigma_devt(:,:,:,1)
endif

elseif(nfint.eq.7)then

call getvar_facesmth(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))
endif

!...Right cell
ipq(1:nvqua) = ipqua(1:nvqua,ier)

!
rc = geoel(1, ier)
sc = geoel(2, ier)
!...Get density for sub-cells...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!...The normal vector at the non-vertex gauss point
do ig = nvfac+1 , ngausf
vngs(1:2, ig, 2) =-vngf(1:2, ig)
vngs(  3, ig, 2) = vngf(  3, ig)
enddo
!
!if(iel.eq.1)print*,'ielgs',iel,ifa,unkno(1:3,:,iel),unkng(1,4,1),unkng(4,4,1)-0.5d0*(unkng(2,4,1)**2 + unkng(3,4,1)**2),&
!unkng(1,5,1),unkng(4,5,1)-0.5d0*(unkng(2,5,1)**2 + unkng(3,5,1)**2)


!...Get the variables at the right state
if(nfint.eq.5.or.nfint.eq.10)then
!...Get the hydrostatic pressure
call getvar_face_general(xpqi, xpq, ipqua(1:nvqua, ier), ipf, fgausr, unkng(:,:,2), sigmg(:,:,:,2), &
unkno(:,:,ier), drhosgq(:, ier), geoel(:,ier), afvec(:,:,ier), aflim(:,ier), ier)

 if(nmatel.eq.2)then
  call getvarsolid_face_general(xpqi, xpq, ipqua(1:nvqua, ier), ipf, fgausr, sigma_devt(:,:,:,2),&
                                unkgd(:,:,ier), strnq_devtp(:,:,:,ier),&
                                geoel(:,ier), afvec(:,:,ier), aflim(:,ier), ier)
!
 sigmg(:,:,:,2) = sigmg(:,:,:,2) + sigma_devt(:,:,:,2)
 endif
!

elseif(nfint.eq.7)then
call getvar_facesmth(ipqua(1:nvqua, ier), ipf, fgausr, unkng(:,:,2), sigmg(:,:,:,2), &
unkno(:,:,ier), drhosgq(:, ier),geoel(:,ier), afvec(:,:,ier), aflim(:,ier))
endif

!...Get the averaged velocity at the gauss quadrature points...
do ig = nvfac+1, ngausf
!vlave(1, fgausl(ig), iel) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
!vlave(2, fgausl(ig), iel) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
!
!vlave(1, fgausr(ig), ier) = 0.5d0*(unkng(2, ig, 2)+unkng(2, ig, 1))
!vlave(2, fgausr(ig), ier) = 0.5d0*(unkng(3, ig, 2)+unkng(3, ig, 1))
enddo

!...Impedence...
do ig = nvfac+1, ngausf
do ic = 1, 2

!aujmp(1:2, ig, ic) = vlave(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgpq(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)

!...Impedence
!delu = sqrt(aujmp(1, ig, ic)**2 + aujmp(2, ig, ic)**2)
delu = cimpd*abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) + unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Summation over corners
munacn = 0.d0
munacu = 0.d0
snsigm = 0.d0
!
do ig = nvfac+1, ngausf
do ic = 1, 2
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ig) = munacn(1:2, 1, ig) + munacn_rie(1:2, 1)
munacn(1:2, 2, ig) = munacn(1:2, 2, ig) + munacn_rie(1:2, 2)
!
munacu(1:2,   ig) = munacu(1:2,   ig) + munacu_rie(1:2)
!
snsigm(1:2,   ig) = snsigm(1:2,   ig) + snsigm_rie(1:2)
!
if(ic.eq.1)then
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
else
!
munaclq(1:2, 1, fgausr(ig), ier) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausr(ig), ier) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausr(ig), ier)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausr(ig), ier)    = snsigm_rie(1:2)
endif
!
enddo

!...Get the corrected velocity at the non-vertex gauss point
detma = munacn(1, 1, ig)*munacn(2, 2, ig) - munacn(2, 1, ig)*munacn(1, 2, ig)
munaci(1, 1) = munacn(2, 2, ig)/detma
munaci(1, 2) =-munacn(1, 2, ig)/detma
munaci(2, 1) =-munacn(2, 1, ig)/detma
munaci(2, 2) = munacn(1, 1, ig)/detma
!
rhsu1 = munacu(1, ig) - snsigm(1, ig)
rhsu2 = munacu(2, ig) - snsigm(2, ig)
!
dufgq(1, fgausl(ig), iel) =-ufgpq(1, fgausl(ig), iel) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
dufgq(2, fgausl(ig), iel) =-ufgpq(2, fgausl(ig), iel) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)
!
dufgq(1, fgausr(ig), ier) =-ufgpq(1, fgausr(ig), ier) + (munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2)
dufgq(2, fgausr(ig), ier) =-ufgpq(2, fgausr(ig), ier) + (munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2)

!...Output for debugging

!if(iel.eq.1)print*,'iel',ig,ier,fgausl(ig),ie,unkng(6,ig,1)*vngs(3, ig, 1), ufgpq(1:2, fgausl(ig), iel),unkng(2:3, ig, 1),&
!unkng(5,ig,1),vngs(3, ig, 1),munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2,munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2,&
!snsigm(1:2, ig),munacu(1:2, ig)
!if(iel.eq.1)print*,'ier',ig,ier,fgausr(ig),ie,unkng(6,ig,2)*vngs(3, ig, 2), ufgpq(1:2, fgausl(ig), iel),unkng(2:3, ig, 2),&
!unkng(5,ig,2),vngs(3, ig, 2),snsigm(1:2, ig),munacu(1:2, ig)
!
enddo

endif
450 enddo
!
!print*,'dufg2',dufgq(1, 17:18, 99)
!
end subroutine getvelo_quadsubg_gs3
!
!...subroutine: Calculate the F^* N ds for interior faces for hybrid grids with Vilar's analytical integration...
!
subroutine getfnds_lagsubgc_maire(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
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
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(1:nvfac,1:4)::ipfsc
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8::dshpr(nvfac),shp(nvfac),shpq(nvqua)
!...local real array
real*8,dimension(1:ndimn, 1:nvtri)::xp
real*8,dimension(1:ndimn, 1:nvqua)::xpq
real*8::posqc(1:2,1:4),xpqv(1:2,1:4)
real*8::xvq(nvqua),yvq(nvqua)
real*8::posi(ngausf), weigh(ngausf)
!...local real number
real*8::anx, any,xsc,ysc,dxdr,dydr
real*8::dr, ds, rc, sc,r,s,djaco
real*8::lanx,lany,wi
real*8::c16,c10
!
data c16   /0.1666666666666666d0 /
data c10   /1.d0 /
!
!...Part I: Specify some gauss pints
!
call ruqope_lobatto(1, ngausf, posi, weigh)
!
posqc(1, 1)=  0.0d0; posqc(2, 1)= -.5d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)=  0.d0;
posqc(1, 3)=  0.0d0; posqc(2, 3)=  .5d0;
posqc(1, 4)= -0.5d0; posqc(2, 4)=  0.d0;
!
ipfsc(1, 1)=5; ipfsc(2, 1)=9; ipfsc(3, 1)=1;
ipfsc(1, 2)=6; ipfsc(2, 2)=9; ipfsc(3, 2)=2;
ipfsc(1, 3)=7; ipfsc(2, 3)=9; ipfsc(3, 3)=3;
ipfsc(1, 4)=8; ipfsc(2, 4)=9; ipfsc(3, 4)=4;

!...Vertex coordinate
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
!...Boundary face...
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
!print*,ifa,iel,gflag(1:3, ifa)
!
enddo
!
!...Quads
!
do 55 ie=1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
do isv = 1, 4
r = posqc(1,isv)
s = posqc(2,isv)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpq(1,ishp)
ysc = ysc + shpq(ishp)*xpq(2,ishp)
enddo
!
xpqv(1, isv) = xsc
xpqv(2, isv) = ysc
!
enddo

!...Get the normal vector and Jacobian for curved nodes...
do ifa = 1,4

!...Edge 9 (p95)
xpf(1:2, 1) =  xpq(1:2, ipfsc(1, ifa))
xpf(1:2, 2) =  xpq(1:2, ipfsc(2, ifa))
xpf(1:2, 3) = xpqv(1:2, ipfsc(3, ifa))
!
lanx = 0.d0
lany = 0.d0
!
do ig =1, ngausf
!
r  = posi(ig)
wi = weigh(ig)
!
shp(1) =  -0.5d0*(1.d0-r)*r
shp(2) =   0.5d0*(1.d0+r)*r
shp(3) =         (1.d0+r)*(1.d0-r)
!
dshpr(1) = -0.5d0 + r
dshpr(2) =  0.5d0 + r
dshpr(3) = -2.d0*r

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dshpr(ishp)*xpf(1, ishp)
dydr = dydr + dshpr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
lanx = lanx + shp(1)*anx*wi
lany = lany + shp(1)*any*wi
enddo
!
gesgq(1, 8+ifa , ie) = lanx/sqrt(lanx**2 + lany**2)
gesgq(2, 8+ifa , ie) = lany/sqrt(lanx**2 + lany**2)
gesgq(3, 8+ifa , ie) = sqrt(lanx**2 + lany**2)

!if(ie.eq.1) print*,'iein',ie, ifa, sqrt(lanx**2 + lany**2),gesgq(1:2, 8+ifa , ie)

enddo !do ifa = 1,4
!
55 enddo  !...(1)ifa=1,nelem
!
end subroutine getfnds_lagsubgc_maire
!
!...subroutine: Calculate the F^* N dsfor all faces for hybrid grids with Maire's method...
!
subroutine getfnds_lag_simpsonh_maire(gflag,gesgq,intfac,inpoel,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
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
integer, dimension(nvfac, 4)::fglvq
integer, dimension(ngausf, 4):: fglgq

!...local real array
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8,dimension(1:ndimn, 1:nvtri)::xpt
real*8,dimension(1:ndimn, 1:nvqua)::xpq
real*8::posi(ngausf), weigh(ngausf)
real*8::vnorm(1:2)
real*8,dimension(1:nvtri):: xv, yv
real*8::posiq(2, 12)
real*8::shp(nvfac),dspr(nvfac)
real*8::lanx(3),lany(3)
real*8::dxdr,dxds,dydr,dyds
!...local real number
real*8::dwav1,dwav2,larea,farea
real*8::anx, any, djaco,wi
real*8::r,s
real*8::dr, ds, rc, sc
real*8::c16, c10
!
data c16   /0.1666666666666666d0 /
data c10   /1.0d0 /
!
!...Specify the posi for 9 simpson points...
!
call ruqope_lobatto(1, ngausf, posi, weigh)

!...Local vertex No. of gauss points in one unit quad...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;

!...Local gauss point No. of any gauss point in one face...
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;
!
! print*,'vnotmfn',gelag(1:3, 1:6, 170)
!
!
do 200 ie=1, nquad !...(1)ifa=1,nquad
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...coordinates
!
if(ncurv==0)then
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif
!
do ifa =1 ,4
!...Edge 9 (p95)
xpf(1, 1:3) =  xpq(1, fglvq(1:3, ifa))
xpf(2, 1:3) =  xpq(2, fglvq(1:3, ifa))
!
lanx = 0.d0
lany = 0.d0
!
do ig =1, ngausf
!
r  = posi(ig)
wi = weigh(ig)
!
shp(1) =  -0.5d0*(1.d0-r)*r
shp(2) =   0.5d0*(1.d0+r)*r
shp(3) =         (1.d0+r)*(1.d0-r)
!
dspr(1) = -0.5d0 + r
dspr(2) =  0.5d0 + r
dspr(3) = -2.d0*r

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dspr(ishp)*xpf(1, ishp)
dydr = dydr + dspr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx = dydr
any =-dxdr
!
lanx(1) = lanx(1) + shp(1)*anx*wi
lany(1) = lany(1) + shp(1)*any*wi

lanx(2) = lanx(2) + shp(2)*anx*wi
lany(2) = lany(2) + shp(2)*any*wi

lanx(3) = lanx(3) + shp(3)*anx*wi
lany(3) = lany(3) + shp(3)*any*wi
!
enddo
!
gelagq(1, fglgq(1, ifa), ie) = lanx(1)/sqrt(lanx(1)**2 + lany(1)**2)
gelagq(2, fglgq(1, ifa), ie) = lany(1)/sqrt(lanx(1)**2 + lany(1)**2)
gelagq(3, fglgq(1, ifa), ie) = sqrt(lanx(1)**2 + lany(1)**2)

!if(ie.eq.1) print*,'ielem',ie,sqrt(lanx(3)**2 + lany(3)**2)
!
gelagq(1, fglgq(2, ifa), ie) = lanx(2)/sqrt(lanx(2)**2 + lany(2)**2)
gelagq(2, fglgq(2, ifa), ie) = lany(2)/sqrt(lanx(2)**2 + lany(2)**2)
gelagq(3, fglgq(2, ifa), ie) = sqrt(lanx(2)**2 + lany(2)**2)
!
gelagq(1, fglgq(3, ifa), ie) = lanx(3)/sqrt(lanx(3)**2 + lany(3)**2)
gelagq(2, fglgq(3, ifa), ie) = lany(3)/sqrt(lanx(3)**2 + lany(3)**2)
gelagq(3, fglgq(3, ifa), ie) = sqrt(lanx(3)**2 + lany(3)**2)

!
enddo
!
gesgq(1:2, 1, ie) = gelagq(1:2, 1, ie);gesgq(3, 1, ie) = gelagq(3, 1, ie);
gesgq(1:2, 2, ie) = gelagq(1:2, 2, ie);gesgq(3, 2, ie) = gelagq(3, 2, ie);
gesgq(1:2, 3, ie) = gelagq(1:2, 3, ie);gesgq(3, 3, ie) = gelagq(3, 3, ie);
gesgq(1:2, 4, ie) = gelagq(1:2, 4, ie);gesgq(3, 4, ie) = gelagq(3, 4, ie);
gesgq(1:2, 5, ie) = gelagq(1:2, 5, ie);gesgq(3, 5, ie) = gelagq(3, 5, ie);
gesgq(1:2, 6, ie) = gelagq(1:2, 6, ie);gesgq(3, 6, ie) = gelagq(3, 6, ie);
gesgq(1:2, 7, ie) = gelagq(1:2, 7, ie);gesgq(3, 7, ie) = gelagq(3, 7, ie);
gesgq(1:2, 8, ie) = gelagq(1:2, 8, ie);gesgq(3, 8, ie) = gelagq(3, 8, ie);
!
gesgq(1:2, 13, ie) = gelagq(1:2, 9, ie);gesgq(3, 13, ie) = gelagq(3, 9, ie);
gesgq(1:2, 14, ie) = gelagq(1:2,10, ie);gesgq(3, 14, ie) = gelagq(3,10, ie);
gesgq(1:2, 15, ie) = gelagq(1:2,11, ie);gesgq(3, 15, ie) = gelagq(3,11, ie);
gesgq(1:2, 16, ie) = gelagq(1:2,12, ie);gesgq(3, 16, ie) = gelagq(3,12, ie);
!
200 enddo  !...(1)ifa=1,nquad
!
!do ig = 1,16
!print*,'Inside getfnds_lag',ig,gesgq(1:3,ig,1435)
!enddo
!
end subroutine getfnds_lag_simpsonh_maire
!
!...subroutine: Calculate the F^* N ds for interior faces for different moments with Vilar's analytical integration...
!
subroutine getfnds_lagsubgc_mairef(rc,sc,bqho,xpq,gesgql)
use constant
implicit none
!...Input arrays
real*8, intent(in)::rc, sc
real*8,dimension(3),            intent(in)::bqho
real*8,dimension(1:ndimn, 1:nvqua),          intent(in)::xpq
real*8,dimension(1:3,1:ngesgq,1:ndegr),      intent(inout)::gesgql
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,isv,ishp
integer::iv
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(1:nvfac,1:4)::ipfsc
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8::dshpr(nvfac),shp(nvfac),shpq(nvqua)
!...local real array
real*8,dimension(1:ndimn, 1:nvtri)::xp
real*8::posqc(1:2,1:4),xpqv(1:2,1:4)
real*8::xvq(nvqua),yvq(nvqua)
real*8::posi(ngausf), weigh(ngausf)
real*8::bq(ndegr,3),bqf(ndegr)
!...local real number
real*8::anx, any,xsc,ysc,dxdr,dydr
real*8::dr, ds, r,s,djaco
real*8::lanx,lany,wi
real*8::c16,c10
!
data c16   /0.1666666666666666d0 /
data c10   /1.d0 /
!
!...Part I: Specify some gauss pints
!
call ruqope_lobatto(1, ngausf, posi, weigh)
!
dr =1.d0
ds = 1.d0
!
posqc(1, 1)=  0.0d0; posqc(2, 1)= -.5d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)=  0.d0;
posqc(1, 3)=  0.0d0; posqc(2, 3)=  .5d0;
posqc(1, 4)= -0.5d0; posqc(2, 4)=  0.d0;
!
ipfsc(1, 1)=5; ipfsc(2, 1)=9; ipfsc(3, 1)=1;
ipfsc(1, 2)=6; ipfsc(2, 2)=9; ipfsc(3, 2)=2;
ipfsc(1, 3)=7; ipfsc(2, 3)=9; ipfsc(3, 3)=3;
ipfsc(1, 4)=8; ipfsc(2, 4)=9; ipfsc(3, 4)=4;

!...Vertex coordinate
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
!...Quads
!
do isv = 1, 4
r = posqc(1,isv)
s = posqc(2,isv)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpq(1,ishp)
ysc = ysc + shpq(ishp)*xpq(2,ishp)
enddo
!
xpqv(1, isv) = xsc
xpqv(2, isv) = ysc
!
enddo

!...Get the normal vector and Jacobian for curved nodes...
do ifa = 1,4

!...Edge 9 (p95)
xpf(1:2, 1) =  xpq(1:2, ipfsc(1, ifa))
xpf(1:2, 2) =  xpq(1:2, ipfsc(2, ifa))
xpf(1:2, 3) = xpqv(1:2, ipfsc(3, ifa))

!
bq(1, 3) = 1.d0
bq(2, 3) = (posqc(1, ifa) - rc)/dr
bq(3, 3) = (posqc(2, ifa) - sc)/ds
!
if(npoly.eq.2)then
bq(4, 3) = 0.5d0*bq(2, 3)*bq(2, 3) - bqho(1)
bq(5, 3) = 0.5d0*bq(3, 3)*bq(3, 3) - bqho(2)
bq(6, 3) =       bq(2, 3)*bq(3, 3) - bqho(3)
endif
!
bq(1, 2) = 1.d0
bq(2, 2) = 0.d0
bq(3, 2) = 0.d0
!
if(npoly.eq.2)then
bq(4, 2) = 0.5d0*bq(2, 2)*bq(2, 2) - bqho(1)
bq(5, 2) = 0.5d0*bq(3, 2)*bq(3, 2) - bqho(2)
bq(6, 2) =       bq(2, 2)*bq(3, 2) - bqho(3)
endif
!
bq(1, 1) = 1.d0
bq(2, 1) = (2.d0*posqc(1, ifa) - rc)/dr
bq(3, 1) = (2.d0*posqc(2, ifa) - sc)/ds
!
if(npoly.eq.2)then
bq(4, 1) = 0.5d0*bq(2, 1)*bq(2, 1) - bqho(1)
bq(5, 1) = 0.5d0*bq(3, 1)*bq(3, 1) - bqho(2)
bq(6, 1) =       bq(2, 1)*bq(3, 1) - bqho(3)
endif
!
do ideg = 1, ndegr
!
lanx = 0.d0
lany = 0.d0
!
do ig =1, ngausf
!
r  = posi(ig)
wi = weigh(ig)
!
shp(1) =  -0.5d0*(1.d0-r)*r
shp(2) =   0.5d0*(1.d0+r)*r
shp(3) =         (1.d0+r)*(1.d0-r)
!
bqf(ideg) = shp(1)*bq(ideg,1) + shp(2)*bq(ideg,2) + shp(3)*bq(ideg,3)
!
dshpr(1) = -0.5d0 + r
dshpr(2) =  0.5d0 + r
dshpr(3) = -2.d0*r

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dshpr(ishp)*xpf(1, ishp)
dydr = dydr + dshpr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx =-dydr
any = dxdr
!
lanx = lanx + shp(1)*anx*wi*bqf(ideg)
lany = lany + shp(1)*any*wi*bqf(ideg)
enddo
!
gesgql(1, 8+ifa , ideg) = lanx!/sqrt(lanx**2 + lany**2)
gesgql(2, 8+ifa , ideg) = lany!/sqrt(lanx**2 + lany**2)
gesgql(3, 8+ifa , ideg) = sqrt(lanx**2 + lany**2)

enddo

enddo !do ifa = 1,4
!
end subroutine getfnds_lagsubgc_mairef
!
!...subroutine: Calculate the F^* N dsfor all faces for different moments with Maire's method...
!
subroutine getfnds_lagsubgcef_mairef(rc,sc,bqho,xpq,gesgql)
use constant
implicit none
!...Input arrays
real*8, intent(in)::rc, sc
real*8,dimension(3),            intent(in)::bqho
real*8,dimension(1:ndimn, 1:nvqua),          intent(in)::xpq
real*8,dimension(1:3,1:ngesgq,1:ndegr),      intent(inout) ::gesgql
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv,ishp
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer, dimension(nvfac, 4)::fglvq
integer, dimension(ngausf, 4):: fglgq

!...local real array
real*8,dimension(1:3,1:ngelgq,1:ndegr)::gelagq
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8::posi(ngausf), weigh(ngausf)
real*8::posiq(2, 12)
real*8::shp(nvfac),dspr(nvfac)
real*8::lanx(3),lany(3)
real*8::dxdr,dxds,dydr,dyds
real*8::bq(ndegr,3),bqf(ndegr)
!...local real number
real*8::dwav1,dwav2,larea,farea
real*8::anx, any, djaco,wi
real*8::r,s
real*8::dr, ds
real*8::c16, c10
!
data c16   /0.1666666666666666d0 /
data c10   /1.0d0 /
!
!...Specify the posi for 9 simpson points...
!
call ruqope_lobatto(1, ngausf, posi, weigh)
!
dr = 1.d0
ds = 1.d0

!...Local vertex No. of gauss points in one unit quad...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;

!...Local gauss point No. of any gauss point in one face...
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

!
!...8-nodes quad
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;
!
do ifa =1 ,4
!...Edge 9 (p95)
xpf(1, 1:3) =  xpq(1, fglvq(1:3, ifa))
xpf(2, 1:3) =  xpq(2, fglvq(1:3, ifa))
!
bq(1, 3) = 1.d0
bq(2, 3) = (posiq(1, fglgq(3 ,ifa)) - rc)/dr
bq(3, 3) = (posiq(2, fglgq(3 ,ifa)) - sc)/ds
!
if(npoly.eq.2)then
bq(4, 3) = 0.5d0*bq(2, 3)*bq(2, 3) - bqho(1)
bq(5, 3) = 0.5d0*bq(3, 3)*bq(3, 3) - bqho(2)
bq(6, 3) =       bq(2, 3)*bq(3, 3) - bqho(3)
endif
!
bq(1, 2) = 1.d0
bq(2, 2) = (posiq(1, fglgq(2 ,ifa)) - rc)/dr
bq(3, 2) = (posiq(2, fglgq(2 ,ifa)) - sc)/ds
!
if(npoly.eq.2)then
bq(4, 2) = 0.5d0*bq(2, 2)*bq(2, 2) - bqho(1)
bq(5, 2) = 0.5d0*bq(3, 2)*bq(3, 2) - bqho(2)
bq(6, 2) =       bq(2, 2)*bq(3, 2) - bqho(3)
endif
!
bq(1, 1) = 1.d0
bq(2, 1) = (posiq(1, fglgq(1 ,ifa)) - rc)/dr
bq(3, 1) = (posiq(2, fglgq(1 ,ifa)) - sc)/ds
!
if(npoly.eq.2)then
bq(4, 1) = 0.5d0*bq(2, 1)*bq(2, 1) - bqho(1)
bq(5, 1) = 0.5d0*bq(3, 1)*bq(3, 1) - bqho(2)
bq(6, 1) =       bq(2, 1)*bq(3, 1) - bqho(3)
endif
!

!
do ideg = 1, ndegr
!
lanx = 0.d0
lany = 0.d0
!
do ig =1, ngausf
!
r  = posi(ig)
wi = weigh(ig)
!
shp(1) =  -0.5d0*(1.d0-r)*r
shp(2) =   0.5d0*(1.d0+r)*r
shp(3) =         (1.d0+r)*(1.d0-r)
!
dspr(1) = -0.5d0 + r
dspr(2) =  0.5d0 + r
dspr(3) = -2.d0*r
!
bqf(ideg) = shp(1)*bq(ideg,1) + shp(2)*bq(ideg,2) + shp(3)*bq(ideg,3)

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dspr(ishp)*xpf(1, ishp)
dydr = dydr + dspr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx = dydr
any =-dxdr
!
lanx(1) = lanx(1) + shp(1)*anx*wi*bqf(ideg)
lany(1) = lany(1) + shp(1)*any*wi*bqf(ideg)

lanx(2) = lanx(2) + shp(2)*anx*wi*bqf(ideg)
lany(2) = lany(2) + shp(2)*any*wi*bqf(ideg)

lanx(3) = lanx(3) + shp(3)*anx*wi*bqf(ideg)
lany(3) = lany(3) + shp(3)*any*wi*bqf(ideg)
!
!if(fglgq(2, ifa).eq.8.and.ideg.eq.3)then
!  print*,'ifa',ifa,fglgq(2, ifa),ideg,bqf(ideg),lanx(2),lany(2),anx,any,wi,shp(2),any*wi*shp(2)*bqf(ideg),bqho(3)
!  print*,'ifa',ifa,fglgq(2, ifa),ideg,bqf(ideg),lanx(2),lany(2),anx,any,wi,shp(2),any*wi*shp(2)*bqf(ideg),bqho(3)
!endif
enddo
!
gelagq(1, fglgq(1, ifa), ideg) = lanx(1)!/sqrt(lanx(1)**2 + lany(1)**2)
gelagq(2, fglgq(1, ifa), ideg) = lany(1)!/sqrt(lanx(1)**2 + lany(1)**2)
gelagq(3, fglgq(1, ifa), ideg) = sqrt(lanx(1)**2 + lany(1)**2)
!
gelagq(1, fglgq(2, ifa), ideg) = lanx(2)!/sqrt(lanx(2)**2 + lany(2)**2)
gelagq(2, fglgq(2, ifa), ideg) = lany(2)!/sqrt(lanx(2)**2 + lany(2)**2)
gelagq(3, fglgq(2, ifa), ideg) = sqrt(lanx(2)**2 + lany(2)**2)
!
!if(fglgq(2, ifa).eq.2) print*,'ideg,',ideg,
!
gelagq(1, fglgq(3, ifa), ideg) = lanx(3)!/sqrt(lanx(3)**2 + lany(3)**2)
gelagq(2, fglgq(3, ifa), ideg) = lany(3)!/sqrt(lanx(3)**2 + lany(3)**2)
gelagq(3, fglgq(3, ifa), ideg) = sqrt(lanx(3)**2 + lany(3)**2)

enddo
!
enddo
!
!
do ideg = 1, ndegr
gesgql(1:2, 1, ideg) = gelagq(1:2, 1, ideg);gesgql(3, 1, ideg) = gelagq(3, 1, ideg);
gesgql(1:2, 2, ideg) = gelagq(1:2, 2, ideg);gesgql(3, 2, ideg) = gelagq(3, 2, ideg);
gesgql(1:2, 3, ideg) = gelagq(1:2, 3, ideg);gesgql(3, 3, ideg) = gelagq(3, 3, ideg);
gesgql(1:2, 4, ideg) = gelagq(1:2, 4, ideg);gesgql(3, 4, ideg) = gelagq(3, 4, ideg);
gesgql(1:2, 5, ideg) = gelagq(1:2, 5, ideg);gesgql(3, 5, ideg) = gelagq(3, 5, ideg);
gesgql(1:2, 6, ideg) = gelagq(1:2, 6, ideg);gesgql(3, 6, ideg) = gelagq(3, 6, ideg);
gesgql(1:2, 7, ideg) = gelagq(1:2, 7, ideg);gesgql(3, 7, ideg) = gelagq(3, 7, ideg);
gesgql(1:2, 8, ideg) = gelagq(1:2, 8, ideg);gesgql(3, 8, ideg) = gelagq(3, 8, ideg);
!
gesgql(1:2, 13, ideg) = gelagq(1:2, 9, ideg);gesgql(3, 13, ideg) = gelagq(3, 9, ideg);
gesgql(1:2, 14, ideg) = gelagq(1:2,10, ideg);gesgql(3, 14, ideg) = gelagq(3,10, ideg);
gesgql(1:2, 15, ideg) = gelagq(1:2,11, ideg);gesgql(3, 15, ideg) = gelagq(3,11, ideg);
gesgql(1:2, 16, ideg) = gelagq(1:2,12, ideg);gesgql(3, 16, ideg) = gelagq(3,12, ideg)
enddo
!
!
end subroutine getfnds_lagsubgcef_mairef
!
!...subroutine: Riemann input for hybrid curved quads using sub-cell scheme for exact integration....
!
subroutine getriem_quadsubg_gauss_maire(ipqua, geoel, vlave, unkno, &
munaclq, munaulq, snsigmlq, drhosgq, gesgql_m, coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),             intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8, dimension(1:2, 1:2, 1:ndegr, 1:2, 1:4, 1:4, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn,  1:ndegr, 1:2,  1:4, 1:4, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn,  1:ndegr, 1:2,  1:4, 1:4, 1:nquad), intent(out)::snsigmlq
real*8,dimension(1:4, 1:nsize), intent(in)::drhosgq
real*8,dimension(1:3,1:ngesgq,1:ndegr, 1:nquad), intent(out)::gesgql_m

!...Local integer
integer::ie, ideg, ielem, ifa, iv, isg, ivsg, ifsg,idegf

!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(8, 4)::fnqsg
integer,dimension(4, 4)::ipqsg
!...local real array
real*8,dimension(1:3,1:ngesgq,1:ndegr)::gesgql
real*8,dimension(1:ndegr, 1:nvqua)::bq
real*8,dimension(1:ndegr, 1:4)::bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nq,1:4)::unsgq
real*8::aujmp(1:3, 1:4)
real*8::vnorm(1:3, 1:2, 1:4),vnorm1(1:3, 1:2, 1:4)
real*8::sigma(1:2, 1:2, 1:4)
real*8,dimension(1:2, 1:4)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndegr, 1:4)::unksgq
!real*8,dimension(1:ndegr, 1:4)::bqsg
real*8,dimension(1:nq+1, 1:4)::uqsgc
real*8,dimension(1:4, 1:4)::prsgq
real*8,dimension(1:4, 1:4)::bqvp
real*8,dimension(1:4)::prqz
real*8,dimension(1:7, 1:4)::geoq_sub
real*8,dimension(2, 4, 4)::wfgsq
real*8::bqsg(ndegr)

!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
real*8::indfc(2,4)
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx,rhovsg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::rcsg,scsg
real*8::dxp,dyp,xc,yc
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon
!
data eps   / 1.0d-6 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Specify some Gauss points...
!

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
wfgsq(1, 1, 1) = 1.d0;      wfgsq(2, 1, 1) = 1.d0;
wfgsq(1, 2, 1) = .5d0;      wfgsq(2, 2, 1) = 4.d0;
wfgsq(1, 3, 1) = 0.d0;      wfgsq(2, 3, 1) = 0.d0;
wfgsq(1, 4, 1) = 4.d0;      wfgsq(2, 4, 1) = .5d0;

wfgsq(1, 1, 2) = 4.d0;      wfgsq(2, 1, 2) = .5d0;
wfgsq(1, 2, 2) = 1.d0;      wfgsq(2, 2, 2) = 1.d0;
wfgsq(1, 3, 2) = .5d0;      wfgsq(2, 3, 2) = 4.d0;
wfgsq(1, 4, 2) = 0.d0;      wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 1, 3) = 0.d0;      wfgsq(2, 1, 3) = 0.d0;
wfgsq(1, 2, 3) = 4.d0;      wfgsq(2, 2, 3) = .5d0;
wfgsq(1, 3, 3) = 1.d0;      wfgsq(2, 3, 3) = 1.d0;
wfgsq(1, 4, 3) = .5d0;      wfgsq(2, 4, 3) = 4.d0

wfgsq(1, 1, 4) = .5d0;      wfgsq(2, 1, 4) = 4.d0;
wfgsq(1, 2, 4) = 0.d0;      wfgsq(2, 2, 4) = 0.d0;
wfgsq(1, 3, 4) = 4.d0;      wfgsq(2, 3, 4) = .5d0;
wfgsq(1, 4, 4) = 1.d0;      wfgsq(2, 4, 4) = 1.d0;
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

!...Vertex coordinate
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

!...Basis function
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
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

!if(ielem.eq.1)then
!print*,'smooth2',ielem,unknvq(1:nq, 8)
!endif

!...Get density for sub-cells...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
! call getrhosubg5(xpq,xpqi,unksgq, geoq_sub)
!
!if(ielem.eq.1) print*,'unksgq',unksgq(1,4)

!...Get density correction
!call getdens_quadsubg(rc, sc, geoel(19:21, ielem), unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc)

!...Zero out local gesgql
gesgql = 0.d0

!...Get the normal vector for inetrior face
call getfnds_lagsubgc_mairef(rc,sc,geoel(19:21, ielem),xpq,gesgql)

!...Get the normal vector for external face
call getfnds_lagsubgcef_mairef(rc,sc,geoel(19:21, ielem),xpq,gesgql)
!
!if(ie.eq.1)print*,'gesgql',gesgql(:,13,1)
!
gesgql_m(:,:,:,ielem) = gesgql

!...II.1: Loop over sub-cells....
do isg = 1, 4

!...Correct the density for Reimann input
do ivsg = 1,4
if(ndens.eq.1)then
rhovt = 1.d0/unknvq(1, ipqsg(ivsg, isg))
!rhovsg = unksgq(1,isg)
rhovsg = rhovt +1.d0*drhosgq(isg, ielem)!(unksgq(1,isg)-1.d0/uqsgc(1, isg))
!
rcsg = geoq_sub(1, isg)
scsg = geoq_sub(2, isg)
!Basis function
bqsg(2) = (xvq(ivsg)-rcsg)/dr
bqsg(3) = (yvq(ivsg)-scsg)/ds

!...DGP2
if(npoly.eq.2)then
bqsg(4) = 0.5d0* bqsg(2)* bqsg(2) - geoq_sub(5, isg)
bqsg(5) = 0.5d0* bqsg(3)* bqsg(3) - geoq_sub(6, isg)
bqsg(6) =        bqsg(2)* bqsg(3) - geoq_sub(7, isg)
endif
!
! rhovsg = 0.d0
do ideg = 1, mdegr
! rhovsg = rhovsg + unksgq(ideg,isg)*bqsg(ideg)
enddo

elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bqv(1, ivsg) = 1.d0
bqv(2, ivsg) = (xvq(ipqsg(ivsg, isg))-rcv)/dr  !...bqv ....
bqv(3, ivsg) = (yvq(ipqsg(ivsg, isg))-scv)/ds
!
rhovt = 0.d0
do ideg = 1, mdegr
rhovt = rhovt + unkno(ideg,1,ielem)*bqv(ideg, ivsg)
enddo
rhovsg = rhovt + (unksgq(1,isg)-rhoct)

!...Based on the sub-cell
!rcv = geoq_sub(1, isg); scv = geoq_sub(2, isg)
!bqv(1, ivsg) = 1.d0
!bqv(2, ivsg) = (xvq(ivsg)-rcv)/dr  !...bqv ....
!bqv(3, ivsg) = (yvq(ivsg)-scv)/ds

!rhovsg =0.d0
!do ideg = 1,ndegr
!rhovsg =  rhovsg + unksgq(ideg,isg)*bqv(ideg, ivsg)
!enddo

!rhovsg = unksgq(1,isg)! unkno(1,1,ielem)!
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
!if(ielem.eq.1)then
!print*,'Variabe',ielem,isg,ivsg,rhovsg,pvtx,unsgq(2:3, ivsg),unksgq(1,isg),1.d0/uqsgc(1, isg)
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
!if(ie.eq.1)print*,'ivsg',ie,ivsg,sigma(1, 1, ivsg)
!
!...Output for debugging
!if(ielem.eq.89.or.ielem.eq.90.or.ielem.eq.99.or.ielem.eq.100)then
!print*,'smooth2',ielem,unkno(:,:,ielem)
!print*,'Variabe2',ielem,ivsg,rhovsg,pvtx,unsgq(2:3, ivsg),sigma(1, 1, ivsg)
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
!
!...Different moments of the RHS
do idegf = 1, 1

!...Normal vector for every quadrature point...
vnorm1(1:2, 1, 1) = sign(1, fnqsg(1, isg))*gesgql(1:2,abs(fnqsg(1, isg)), idegf);
vnorm1(1:2, 2, 1) = sign(1, fnqsg(2, isg))*gesgql(1:2,abs(fnqsg(2, isg)), idegf);

vnorm1(1:2, 1, 2) = sign(1, fnqsg(3, isg))*gesgql(1:2,abs(fnqsg(3, isg)), idegf);
vnorm1(1:2, 2, 2) = sign(1, fnqsg(4, isg))*gesgql(1:2,abs(fnqsg(4, isg)), idegf);

vnorm1(1:2, 1, 3) = sign(1, fnqsg(5, isg))*gesgql(1:2,abs(fnqsg(5, isg)), idegf);
vnorm1(1:2, 2, 3) = sign(1, fnqsg(6, isg))*gesgql(1:2,abs(fnqsg(6, isg)), idegf);

vnorm1(1:2, 1, 4) = sign(1, fnqsg(7, isg))*gesgql(1:2,abs(fnqsg(7, isg)), idegf);
vnorm1(1:2, 2, 4) = sign(1, fnqsg(8, isg))*gesgql(1:2,abs(fnqsg(8, isg)), idegf);

enddo

!...Different moments of the RHS
do idegf = 1, ndegr
!
!
indfc = 1.d0
!...Normal vector for every quadrature point...
vnorm(1:2, 1, 1) = sign(1, fnqsg(1, isg))*gesgql(1:2,abs(fnqsg(1, isg)), idegf);
vnorm(  3, 1, 1) = gesgql(3,abs(fnqsg(1, isg)), idegf);

vnorm(1:2, 2, 1) = sign(1, fnqsg(2, isg))*gesgql(1:2,abs(fnqsg(2, isg)), idegf);
vnorm(  3, 2, 1) = gesgql(3,abs(fnqsg(2, isg)), idegf);

vnorm(1:2, 1, 2) = sign(1, fnqsg(3, isg))*gesgql(1:2,abs(fnqsg(3, isg)), idegf);
vnorm(  3, 1, 2) = gesgql(3,abs(fnqsg(3, isg)), idegf);

vnorm(1:2, 2, 2) = sign(1, fnqsg(4, isg))*gesgql(1:2,abs(fnqsg(4, isg)), idegf);
vnorm(  3, 2, 2) = gesgql(3,abs(fnqsg(4, isg)), idegf);

vnorm(1:2, 1, 3) = sign(1, fnqsg(5, isg))*gesgql(1:2,abs(fnqsg(5, isg)), idegf);
vnorm(  3, 1, 3) = gesgql(3,abs(fnqsg(5, isg)), idegf);

vnorm(1:2, 2, 3) = sign(1, fnqsg(6, isg))*gesgql(1:2,abs(fnqsg(6, isg)), idegf);
vnorm(  3, 2, 3) = gesgql(3,abs(fnqsg(6, isg)), idegf);

vnorm(1:2, 1, 4) = sign(1, fnqsg(7, isg))*gesgql(1:2,abs(fnqsg(7, isg)), idegf);
vnorm(  3, 1, 4) = gesgql(3,abs(fnqsg(7, isg)), idegf);

vnorm(1:2, 2, 4) = sign(1, fnqsg(8, isg))*gesgql(1:2,abs(fnqsg(8, isg)), idegf);
vnorm(  3, 2, 4) = gesgql(3,abs(fnqsg(8, isg)), idegf);
!
!...Get weighted area normal vector
do ifsg =1,2
if((vnorm(1, ifsg, 1)*vnorm1(1, ifsg, 1) + vnorm(2, ifsg, 1)*vnorm1(2, ifsg, 1)).lt.0.d0)then
vnorm(1:2, ifsg, 1) = -vnorm(1:2, ifsg, 1)
indfc(ifsg, 1) = -1.d0
endif

if((vnorm(1, ifsg, 2)*vnorm1(1, ifsg, 2) + vnorm(2, ifsg, 2)*vnorm1(2, ifsg, 2)).lt.0.d0)then
vnorm(1:2, ifsg, 2) = -vnorm(1:2, ifsg, 2)
indfc(ifsg, 2) = -1.d0
endif

if((vnorm(1, ifsg, 3)*vnorm1(1, ifsg, 3) + vnorm(2, ifsg, 3)*vnorm1(2, ifsg, 3)).lt.0.d0)then
vnorm(1:2, ifsg, 3) = -vnorm(1:2, ifsg, 3)
indfc(ifsg, 3) = -1.d0
endif

if((vnorm(1, ifsg, 4)*vnorm1(1, ifsg, 4) + vnorm(2, ifsg, 4)*vnorm1(2, ifsg, 4)).lt.0.d0)then
vnorm(1:2, ifsg, 4) = -vnorm(1:2, ifsg, 4)
indfc(ifsg, 4) = -1.d0
endif
enddo
!
!if(ie.eq.1.and.isg.eq.1)then
!print*,'iel',idegf,ie,isg,vnorm(1:3, 1:2, 1),indfc(1, 1)
!endif

!if(ie.eq.10.and.isg.eq.2)then
!print*,'iel',idegf,ie,isg,vnorm(1:3, 1:2, 2)
!endif

!...Get weighted area normal vector
do ifsg =1,2
vnorm(1:3, ifsg, 1) =  vnorm(  1:3, ifsg, 1)*wfgsq(ifsg, 1, isg);
vnorm(1:3, ifsg, 2) =  vnorm(  1:3, ifsg, 2)*wfgsq(ifsg, 2, isg);
vnorm(1:3, ifsg, 3) =  vnorm(  1:3, ifsg, 3)*wfgsq(ifsg, 3, isg);
vnorm(1:3, ifsg, 4) =  vnorm(  1:3, ifsg, 4)*wfgsq(ifsg, 4, isg);
enddo
!
!vnorm(3,:,:) =1.d0

!...Get impedence coefficient...
do ivsg   = 1, 4
dux= vlave(1, ipq(ipqsg(ivsg, isg)))-unsgq(2, ivsg)
duy= vlave(2, ipq(ipqsg(ivsg, isg)))-unsgq(3, ivsg)
deltu = sqrt(dux**2 + duy**2)
do ifa = 1, 2
!deltu = 1.d0*abs(dux*vnorm(1, ifa, ivsg) + duy*vnorm(2, ifa, ivsg))
murie(ifa, ivsg) = rhoct*sdctr + rhoct*slpdu*deltu
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
!
do ifa = 1, 2 !...Every corner consists of 2 faces...

!...Call Riemann solver...
call getriecoef_matrixnew2(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:3, ivsg), &
unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ifa, ivsg), vnorm(3, ifa, ivsg), vnorm(1:2, ifa, ivsg), aujmp(1:2, ivsg), &
!unsgq(2:3, ivsg), sigma(1:2, 1:2, ivsg),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!

!...Output for debugging
!if(ie.eq.1) then
!print*,'p36 muacn(vv) post',ie,isg,ivsg,ipq(iv),idegf,snsigm_rie(1:2),vnorm(1:3, ifa, ivsg),munacu_rie(1:2)
!endif

!if(ie.eq.10) then
!print*,'p36 muacn(vv) post',ie,isg,ivsg,ipq(iv),idegf,snsigm_rie(1:2),vnorm(1:3, ifa, ivsg),munacu_rie(1:2)
!endif

!...Local variable...
munaclq(1:2, 1, idegf,ifa, ivsg, isg, ie) =  munacn_rie(1:2, 1)*indfc(ifa,ivsg)
munaclq(1:2, 2, idegf, ifa, ivsg, isg, ie) =  munacn_rie(1:2, 2)*indfc(ifa,ivsg)
!
munaulq(1:2,  idegf,   ifa, ivsg, isg, ie) =  munacu_rie(1:2)*indfc(ifa,ivsg)
!
snsigmlq(1:2,  idegf,  ifa, ivsg, isg, ie)=  snsigm_rie(1:2)*indfc(ifa,ivsg)
!
!
enddo
enddo
!
enddo !...idegf
!
enddo
!
350 enddo  !...(1)ie = 1,nquad

end subroutine getriem_quadsubg_gauss_maire
!
!...Get the nodal velocity based on exact face integration...
!
subroutine getndvelo_lagsubg_gaussh_maire(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,unkgd,ustar, rhsel, aflim, afvec, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac), intent(in)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:ndegr,1:nq,1:ncell),      intent(out)::rhsel

integer:: itime,ip
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop,idegf
integer::isg,ivsg
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(4, 4)::ipqsg
integer::indnd(npoin)
integer, dimension(3, 4)::fglvq
integer, dimension(3, 3)::fglvt
!
integer, dimension(ngausf, 4):: fglgq
integer, dimension(ngausf, 3):: fglgt
integer,dimension(8, 4)::fnqsg

!...local real array
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:ndimn,1:nvqua,1:nquad)::vnulq
real*8,  dimension(1:nquad)::gqdmp
real*8::vnorm(1:3, 1:2, 1:6,1:ndegr)
real*8::munaci(2, 2)
real*8, dimension(1:ndegr)::ulnpn,elnpn
real*8, dimension(1:2, 1:ndegr)::plnpn
real*8,dimension(2, 4, 4)::wfgsq
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2
real*8::r,shp1,shp2,shp3
real*8::pc1
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: munacn(:,:,:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: munacu(:,:), snsigm(:,:)
real*8,allocatable:: munaclt(:,:,:,:),munaclq(:,:,:,:,:,:)
real*8,allocatable:: snsigmlt(:,:,:), munault(:,:,:)
real*8,allocatable:: snsigmlq(:,:,:,:,:), munaulq(:,:,:,:,:)
real*8,allocatable:: dufgq(:,:,:)
real*8,allocatable:: munaclq_m(:,:,:,:,:,:,:),munaulq_m(:,:,:,:,:,:),snsigmlq_m(:,:,:,:,:,:),fstarq_m(:,:,:,:,:,:)
real*8,allocatable:: gesgql_m(:,:,:,:)
real*8,allocatable:: drhosgq(:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:2, 1:2, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munaclt(1:2, 1:2, 1:ngelg, 1:ntria), munault(1:ndimn, 1:ngelg,  1:ntria),&
snsigmlt(1:ndimn, 1:ngelg,  1:ntria))
allocate (munaclq(1:2, 1:2, 1:2, 1:4, 1:4, 1:nquad), munaulq(1:ndimn, 1:2, 1:4, 1:4,  1:nquad),&
snsigmlq(1:ndimn, 1:2,  1:4, 1:4,  1:nquad))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
allocate (munaulq_m(1:2,1:ndegr, 1:2, 1:4, 1:4, 1:nquad),snsigmlq_m(1:2,1:ndegr, 1:2, 1:4, 1:4, 1:nquad))
allocate (fstarq_m(1:2,1:ndegr, 1:2, 1:4, 1:4, 1:nquad))
allocate (munaclq_m(1:2, 1:2,1:ndegr, 1:2, 1:4, 1:4, 1:nquad))
allocate(gesgql_m(1:3,1:ngesgq,1:ndegr, 1:nquad))
allocate(drhosgq(1:4, 1:nsize))
!
!...Part I: Specify some gauss points...
!

!...Local vertex No. of gauss points in one unit quad...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;

!...Local vertex No. of gauss points in one unit tria...
fglvt(1, 1) = 1;  fglvt(2, 1) = 2; fglvt(3, 1) = 4;
fglvt(1, 2) = 2;  fglvt(2, 2) = 3; fglvt(3, 2) = 5;
fglvt(1, 3) = 3;  fglvt(2, 3) = 1; fglvt(3, 3) = 6;

!...Local gauss point No. of any gauss point in one face...
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;
!
fglgt(1, 1) = 1;  fglgt(2, 1) = 2; fglgt(3, 1) = 7;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4; fglgt(3, 2) = 8;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6; fglgt(3, 3) = 9;
!
elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
elseif(ngausf.eq.5)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;
!
endif

!
!...Local vertex No. of gauss points in one sub-cell
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
wfgsq(1, 1, 1) = 1.d0;      wfgsq(2, 1, 1) = 1.d0;
wfgsq(1, 2, 1) = .5d0;      wfgsq(2, 2, 1) = 4.d0;
wfgsq(1, 3, 1) = 0.d0;      wfgsq(2, 3, 1) = 0.d0;
wfgsq(1, 4, 1) = 4.d0;      wfgsq(2, 4, 1) = .5d0;

wfgsq(1, 1, 2) = 4.d0;      wfgsq(2, 1, 2) = .5d0;
wfgsq(1, 2, 2) = 1.d0;      wfgsq(2, 2, 2) = 1.d0;
wfgsq(1, 3, 2) = .5d0;      wfgsq(2, 3, 2) = 4.d0;
wfgsq(1, 4, 2) = 0.d0;      wfgsq(2, 4, 2) = 0.d0;

wfgsq(1, 1, 3) = 0.d0;      wfgsq(2, 1, 3) = 0.d0;
wfgsq(1, 2, 3) = 4.d0;      wfgsq(2, 2, 3) = .5d0;
wfgsq(1, 3, 3) = 1.d0;      wfgsq(2, 3, 3) = 1.d0;
wfgsq(1, 4, 3) = .5d0;      wfgsq(2, 4, 3) = 4.d0

wfgsq(1, 1, 4) = .5d0;      wfgsq(2, 1, 4) = 4.d0;
wfgsq(1, 2, 4) = 0.d0;      wfgsq(2, 2, 4) = 0.d0;
wfgsq(1, 3, 4) = 4.d0;      wfgsq(2, 3, 4) = .5d0;
wfgsq(1, 4, 4) = 1.d0;      wfgsq(2, 4, 4) = 1.d0;
!
!...Zero out vlave
!
vlave = 0.d0
indnd = 0
!
!...Mark the boundary nodes...
!
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

!...Get averaged velocity at nodes...
!call getvlavenew(iptri, ipqua, geoel, vlave, unkno, aflim, afvec)
!
do 950 ifa = 1 , nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
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
endif!
950 enddo
!
!...Part II: Loop to get information from Riemann solver
!
do iloop= 1, 1!5

!...Give value
vlave= ustar
vnulq = 0.d0

!...Zero out munacn
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0

!...Triangle in the future...
!print*,'gesgq',gesgq(1:3,:,1)
!...Quad
!if(nquad.gt.0) call getriem_quadsubg_gauss(ipqua, geoel, gesgq, vlave, unkno, unkgd, munacn, munacu, snsigm,&
!munaclq, munaulq, snsigmlq, drhosgq, coord, coold, aflim, afvec)

!....Update the velocity at the end points...
do ipoin = 1,npoin

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

!if(ipoin.eq.2.or.ipoin.eq.22)print*,'ipnode',ustar(1:2,ipoin),detma,munacn(:,:,ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin)
endif
enddo

!
!print*,'gesgq2',gesgq(1:3,:,1)

!print*,'ustar',ustar(1,ipqua(1:8,1))

!...Get the vertex velocity at the boundary....
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
endif
!
900 enddo

enddo
!
if(nquad.gt.0) call getriem_quadsubg_gauss_maire(ipqua, geoel, vlave, unkno, &
munaclq_m, munaulq_m, snsigmlq_m, drhosgq, gesgql_m, coord, coold, aflim, afvec)
!
!...4.2: Update the Riemann forces at every node...
!...Tria
!
!...Quad
!
do ie = 1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
ustar(1:2, ipq(9))= 0.d0
!
do isg = 1, 4
do ivsg = 1, 4
!
iv =  ipqsg(ivsg, isg)
!
do ifa =1, 2
do ideg = 1, ndegr
fstarq_m(1, ideg, ifa, ivsg, isg, ie) = snsigmlq_m(1,  ideg, ifa, ivsg, isg, ie) + &
munaclq_m(1, 1,  ideg, ifa, ivsg, isg,ie)*ustar(1, ipq(iv))+&
munaclq_m(2, 1,  ideg, ifa, ivsg, isg,ie)*ustar(2, ipq(iv)) - munaulq_m(1,  ideg, ifa,ivsg, isg, ie)
fstarq_m(2,  ideg, ifa, ivsg, isg,ie) = snsigmlq_m(2,  ideg, ifa, ivsg, isg, ie) + &
munaclq_m(1, 2,  ideg, ifa, ivsg, isg, ie)*ustar(1, ipq(iv))+&
munaclq_m(2, 2,  ideg, ifa, ivsg, isg, ie)*ustar(2, ipq(iv))- munaulq_m(2,  ideg, ifa, ivsg, isg, ie)
!
enddo

enddo
enddo
enddo
!
fstarq_m(:, :, :, 3, 1, ie) = 0.d0
fstarq_m(:, :, :, 4, 2, ie) = 0.d0
fstarq_m(:, :, :, 1, 3, ie) = 0.d0
fstarq_m(:, :, :, 2, 4, ie) = 0.d0
!
!if(ie.eq.1.)then
!print*,'velo',ustar(1:2,ipq(1))
!print*,'idegf',ie, fstarq_m(1, 3, 1:2, 1, 1, ie)*0.1d0/3.3333333333333375d-002,&
!snsigmlq_m(1,  3, 1:2, 1, 1, ie)*0.1d0/3.3333333333333375d-002,munaulq_m(1,  3, 1:2, 1, 1, ie)*0.1d0/3.3333333333333333d-2
!endif

!if(ie.eq.5)then
!print*,'velo',ustar(1:2,ipq(2))
!print*,'idegf',ie, fstarq_m(2, 2, 1:2, 2, 2, ie)*0.1d0/3.3333333333333375d-002,&
!munaclq_m(1:2,2,  2, 1, 2, 2, ie),ustar(1:2,11),snsigmlq_m(2,  2, 1, 2, 2, ie),&
!munaulq_m(2,  2, 1, 2, 2, ie)
!endif

enddo
!
!...Update the RHS
!
rhsel = 0.d0
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...The vertex constituting one cell...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
pc1 = 0.d0
!
do isg = 1, 4
!
do idegf = 1, ndegr
!...Normal vector for every quadrature point...
vnorm(1:2, 1, 1, idegf) = sign(1, fnqsg(1, isg))*gesgql_m(1:2,abs(fnqsg(1, isg)), idegf, ie);
vnorm(1:2, 2, 1, idegf) = sign(1, fnqsg(2, isg))*gesgql_m(1:2,abs(fnqsg(2, isg)), idegf, ie);

vnorm(1:2, 1, 2, idegf) = sign(1, fnqsg(3, isg))*gesgql_m(1:2,abs(fnqsg(3, isg)), idegf, ie);
vnorm(1:2, 2, 2, idegf) = sign(1, fnqsg(4, isg))*gesgql_m(1:2,abs(fnqsg(4, isg)), idegf, ie);

vnorm(1:2, 1, 3, idegf) = sign(1, fnqsg(5, isg))*gesgql_m(1:2,abs(fnqsg(5, isg)), idegf, ie);
vnorm(1:2, 2, 3, idegf) = sign(1, fnqsg(6, isg))*gesgql_m(1:2,abs(fnqsg(6, isg)), idegf, ie);

vnorm(1:2, 1, 4, idegf) = sign(1, fnqsg(7, isg))*gesgql_m(1:2,abs(fnqsg(7, isg)), idegf, ie);
vnorm(1:2, 2, 4, idegf) = sign(1, fnqsg(8, isg))*gesgql_m(1:2,abs(fnqsg(8, isg)), idegf, ie);
enddo
!
!...loop over every subgrid vertex...
!
do ivsg = 1, 4
!
ulnpn(1:ndegr) = ulnpn(1:ndegr)  +&
ustar(1, ipq(ipqsg(ivsg, isg)))*vnorm(1, 1, ivsg, 1:ndegr)*wfgsq(1,ivsg,isg) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*vnorm(2, 1, ivsg, 1:ndegr)*wfgsq(1,ivsg,isg) +&
ustar(1, ipq(ipqsg(ivsg, isg)))*vnorm(1, 2, ivsg, 1:ndegr)*wfgsq(2,ivsg,isg) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*vnorm(2, 2, ivsg, 1:ndegr)*wfgsq(2,ivsg,isg)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq_m(1, 1:ndegr,1, ivsg, isg, ie) +&
fstarq_m(1, 1:ndegr,2, ivsg, isg, ie)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq_m(2, 1:ndegr, 1, ivsg, isg, ie) +&
fstarq_m(2, 1:ndegr, 2, ivsg, isg, ie)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(ipqsg(ivsg, isg)))*fstarq_m(1, 1:ndegr,1, ivsg, isg,  ie) +&
ustar(2, ipq(ipqsg(ivsg, isg)))*fstarq_m(2, 1:ndegr,1, ivsg, isg,  ie) +&
ustar(1, ipq(ipqsg(ivsg, isg)))*fstarq_m(1, 1:ndegr,2, ivsg, isg, ie)  +&
ustar(2, ipq(ipqsg(ivsg, isg)))*fstarq_m(2, 1:ndegr,2, ivsg, isg,  ie)
!
!if(ie.eq.1)then
!if(isg.eq.1.and.ivsg.eq.1)pc1 = pc1 + fstarq_m(1, 3,1, ivsg, isg, ie)+ fstarq_m(1, 3,2, ivsg, isg, ie)
!print*,'vnorm1',ie,isg,ivsg,ipq(ipqsg(ivsg, isg)),pc1,fstarq_m(1, 3,1:2, ivsg, isg, ie)
!print*,'vnorm2',ie,isg,ivsg,ipq(ipqsg(ivsg, isg)),vnorm(1:2, 1:2, ivsg, 3),wfgsq(1,ivsg,isg)
!endif

!if(ie.eq.5)then
!if(isg.eq.2.and.ivsg.eq.2)pc1 = pc1 + fstarq_m(2, 2,1, ivsg, isg, ie)+ fstarq_m(2, 2,2, ivsg, isg, ie)
!print*,'vnorm1',ie,isg,ivsg,ipq(ipqsg(ivsg, isg)),pc1,fstarq_m(2, 2,1:2, ivsg, isg, ie)
!print*,'vnorm2',ie,isg,ivsg,ipq(ipqsg(ivsg, isg)),vnorm(1:2, 1:2, ivsg, 3),wfgsq(1,ivsg,isg)
!endif

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
!if(ie.eq.1)then
! print*,'rhs iface',ielem, ie,plnpn(1, 3)
!endif

!if(ie.eq.5)then
!print*,'rhs iface',ielem, ie,plnpn(2, 2)
!endif

650 enddo
!
deallocate (munacn, bpres, fpres)
deallocate (munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
deallocate (munaclq_m, snsigmlq_m, munaulq_m, gesgql_m)
deallocate (drhosgq)
end subroutine getndvelo_lagsubg_gaussh_maire
!
!...Calculate the Riemann solutions at the non-vertex quadrature point...
!
subroutine getriem_gnvtx(ufgpq, geoel,bface,intfac,ipqua,coord,unkno, &
munaclq, munaulq, snsigmlq, drhosgq, afvec, aflim, itime)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:2,1:ngsfq,nquad),          intent(inout)::ufgpq
real*8, dimension(1:2, 1:2, 1:ngsfq, 1:nquad), intent(inout)::munaclq
real*8, dimension(1:ndimn, 1:ngsfq,  1:nquad), intent(inout)::munaulq
real*8, dimension(1:ndimn, 1:ngsfq,  1:nquad), intent(inout)::snsigmlq
real*8,dimension(1:4, 1:nsize), intent(in)::drhosgq
integer*4,          intent(in)::itime
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(ngausf)::fgausl, fgausr
!
real*8::eps
real*8,dimension(1:2,1:ngelgq,nquad)::vlave
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::sigmg(1:2, 1:2, 1:ngausf, 1:2)
real*8::vngs(3, ngausf,2)
real*8::vngf(3, ngausf)
real*8::aujmp(1:3,1:ngausf,1:2)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, ngausf)
real*8::posiq(1:2, 1:ngsfq)
integer,dimension(nvfac, 4)::fglvq
integer,dimension(ngausf,4)::fglgq
real*8::xpf(1:2, 1:nvfac)
real*8::fpres(1:2, 1:ngausf)

!...Riemann parameters...
real*8::munacn_rie(2,2), munacu_rie(2), snsigm_rie(2)

!...Local real number
real*8::pgl,pgr,ugl,ugr,vgl,vgr,mul,mur
real*8::delu
real*8::acnx,acny,ngx,ngy
!
!...Parameters specification
!
ngvtxf = nvfac
eps = 1.d-6
!
!...Part I: Specify some gauss points
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
xvq(9) =  0.d0; yvq(9) =  0.d0

!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
!...Part II: Get the variables at the left and right side of the face...
!
fgausl = 0.d0
fgausr = 0.d0
!
do 450 ifa = 1, nafac !...(1)ifa = 1, nafac

!...Face vertices
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)

!...Left and right cell
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!...
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))

!...Get the normal vector at the non-vertex gauss point
call getvectf_gnvtx(xpf, vngf)

!...Boundary face
if(ifa.le.nbfac)then

!...The normal vector at the non-vertex gauss point
do ig = ngvtxf+1, ngausf
vngs(1:2, ig, 1) = vngf(1:2, ig)
vngs(  3, ig, 1) = vngf(  3, ig)
enddo

!...Get the variables at the left state
call getvar_face_gnvtx(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))

!...Impedence...
do ig = ngvtxf+1, ngausf
do ic = 1, 1
!
aujmp(1:2, ig, ic) = ufgpq(1:2, fgausl(ig), iel) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)
!aujmp(3, ig, ic) = 1.d0*abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
!
unkng(6, ig, ic) = unkng(6, ig, ic) + unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Boundary condition imposing
!call getbcgauss_lagc(bface(3, ifa), vngf, fpres, itime)


!...Excluding the non-vertex gs point
do ig = ngvtxf+1, ngausf
do ic = 1, 1
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
!
enddo
enddo
!
!...Riemann velocity at the non-vertex quadrature point
do ig = ngvtxf+1, ngausf
!
pgl = unkng(5, ig, 1)

ugl = unkng(2, ig, 1)

vgl = unkng(3, ig, 1)

mul = unkng(6, ig, 1)

ngx = vngs(1, ig, 1)
ngy = vngs(2, ig, 1)

!...Impose the boundary condition to adjust the ustar at the boundary
if(bface(3 ,ifa).eq.22)then
if(bface(4 ,ifa).eq.221)then
pgr = pgl
ugr = ugl
vgr =-vgl
mur = mul
elseif(bface(4 ,ifa).eq.222)then
pgr = pgl
ugr =-ugl
vgr = vgl
mur = mul
endif

elseif(bface(3 ,ifa).eq.25)then
pgr = pgl
ugr = 0.d0
vgr = 0.d0
mur = mul
elseif(bface(3 ,ifa).eq.27)then
print*,'bface(3 ,ifa).eq.27 will be implemented in future for non-vertex quadrature points!'
stop
elseif(bface(3 ,ifa).eq.21)then
print*,'bface(3 ,ifa).eq.21 will be implemented in future for non-vertex quadrature points!'
stop
elseif(bface(3 ,ifa).eq.23)then
print*,'bface(3 ,ifa).eq.23 will be implemented in future for non-vertex quadrature points!'
stop
endif

!...Riemann velocity at the non-vertex quadrature points
ufgpq(1, fgausl(ig), iel) = (mul*ugl + mur*ugr)/(mul+mur) -(pgr-pgl)/(mul+mur)*ngx
ufgpq(2, fgausl(ig), iel) = (mul*vgl + mur*vgr)/(mul+mur) -(pgr-pgl)/(mul+mur)*ngy


enddo

!...Interior face
elseif(ifa.gt.nbfac)then

!...The normal vector at the non-vertex gauss point
do ig = ngvtxf+1, ngausf
vngs(1:2, ig, 1) = vngf(1:2, ig)
vngs(  3, ig, 1) = vngf(  3, ig)
enddo

!...Get the variables at the left state
call getvar_face_gnvtx(ipqua(1:nvqua, iel), ipf, fgausl, unkng(:,:,1), sigmg(:,:,:,1), &
unkno(:,:,iel), drhosgq(:, iel), geoel(:,iel), afvec(:,:,iel), aflim(:,iel))

!...The normal vector at the non-vertex gauss point
do ig = ngvtxf+1 , ngausf
vngs(1:2, ig, 2) =-vngf(1:2, ig)
vngs(  3, ig, 2) = vngf(  3, ig)
enddo
!
!if(iel.eq.1)print*,'ielgs',iel,ifa,unkno(1:3,:,iel),unkng(1,4,1),unkng(4,4,1)-0.5d0*(unkng(2,4,1)**2 + unkng(3,4,1)**2),&
!unkng(1,5,1),unkng(4,5,1)-0.5d0*(unkng(2,5,1)**2 + unkng(3,5,1)**2)


!...Get the variables at the right state
call getvar_face_gnvtx(ipqua(1:nvqua, ier), ipf, fgausr, unkng(:,:,2), sigmg(:,:,:,2), &
unkno(:,:,ier), drhosgq(:, ier),geoel(:,ier), afvec(:,:,ier), aflim(:,ier))


!...Impedence...
do ig = ngvtxf+1, ngausf
do ic = 1, 2

!aujmp(1:2, ig, ic) = vlave(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
aujmp(1:2, ig, ic) = ufgpq(1:2, fgausr(ig), ier) - unkng(2:3, ig, ic)
!
acnx = aujmp(1, ig, ic)
acny = aujmp(2, ig, ic)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig, ic) = 1.e-11!0.d0;
else
aujmp(1:2, ig, ic) = aujmp(1:2, ig, ic)/sqrt(acnx**2 + acny**2)
endif
aujmp(3, ig, ic) = sqrt(acnx**2 + acny**2)/unkng(7, ig, ic)

!...Impedence
delu = cimpd*abs(aujmp(1, ig, ic)*vngs(1, ig, ic) + aujmp(2, ig, ic)*vngs(2, ig, ic))
unkng(6, ig, ic) = unkng(6, ig, ic) + unkng(6, ig, ic)/unkng(7, ig, ic)*slpdu*delu
enddo
enddo

!...Riemann velocity at the non-vertex quadrature point
do ig = ngvtxf+1, ngausf
!
pgl = unkng(5, ig, 1)
pgr = unkng(5, ig, 2)

ugl = unkng(2, ig, 1)
ugr = unkng(2, ig, 2)

vgl = unkng(3, ig, 1)
vgr = unkng(3, ig, 2)

mul = unkng(6, ig, 1)
mur = unkng(6, ig, 2)

ngx = vngs(1, ig, 1)
ngy = vngs(2, ig, 1)

!...Riemann velocity at the non-vertex quadrature points
ufgpq(1, fgausl(ig), iel) = (mul*ugl + mur*ugr)/(mul+mur) -(pgr-pgl)/(mul+mur)*ngx
ufgpq(2, fgausl(ig), iel) = (mul*vgl + mur*vgr)/(mul+mur) -(pgr-pgl)/(mul+mur)*ngy

ufgpq(1, fgausr(ig), ier) = ufgpq(1, fgausl(ig), iel)
ufgpq(2, fgausr(ig), ier) = ufgpq(2, fgausl(ig), iel)
enddo

!...Terms for Riemann solver
do ig = ngvtxf+1, ngausf
do ic = 1, 2
!
call getriecoef_matrixnew(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:3, ig, ic), &
unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(unkng(6, ig, ic), vngs(3, ig, ic), vngs(1:2, ig, ic), aujmp(1:2, ig, ic), &
!unkng(2:3, ig, ic), sigmg(1:2, 1:2, ig, ic),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
if(ic.eq.1)then
!
munaclq(1:2, 1, fgausl(ig), iel) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausl(ig), iel) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausl(ig), iel)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausl(ig), iel)    = snsigm_rie(1:2)
else
!
munaclq(1:2, 1, fgausr(ig), ier) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fgausr(ig), ier) =  munacn_rie(1:2, 2)
!
munaulq(1:2, fgausr(ig), ier)    =  munacu_rie(1:2)
!
snsigmlq(1:2,fgausr(ig), ier)    = snsigm_rie(1:2)
endif
!
enddo
!
enddo

endif
450 enddo
!
!print*,'dufg2',dufgq(1, 17:18, 99)
!
end subroutine getriem_gnvtx
!
!...Calculate the normal vector at the non-vertex gauss point
!
subroutine getvectf_gnvtx(xpf, vngf)
use constant
implicit none
real*8, dimension(3, ngausf),   intent(out)::vngf
real*8, dimension(1:2, 1:nvfac), intent(in)::xpf
!...Local integer
integer::ig, ishp
integer::iv
!...Local real
real*8::dshpr(3),posig(ngausf)
real*8::dxdr, dydr,djaco,anx,any,r
!
if(ngausf.eq.4)then
posig(3) = -0.4472135954999579392818d0
posig(4) =  0.4472135954999579392818d0
elseif(ngausf.eq.5)then
posig(4) = -0.6546536707079771437983d0
posig(5) =  0.6546536707079771437983d0
endif
!
do ig = ngvtxf+1, ngausf
!
r  = posig(ig)

dshpr(1) = -0.5d0 + r
dshpr(2) =  0.5d0 + r
dshpr(3) = -2.d0*r

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, 3
dxdr = dxdr + dshpr(ishp)*xpf(1, ishp)
dydr = dydr + dshpr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
anx = dydr
any =-dxdr

!...Normal vector at the gauss point
vngf(1, ig) = anx/sqrt(anx**2 + any**2)
vngf(2, ig) = any/sqrt(anx**2 + any**2)
vngf(3, ig) = djaco*2.d0

enddo

end subroutine getvectf_gnvtx
!
!...Get the left and right state at one face...
!
subroutine getvar_face_gnvtx(ipqua, ipf, fgaus, unkng, sigmg, unkno, drhosgq, geoel, afvec, aflim)
use constant
implicit none
!...Input integer arrays
integer,dimension(1:nvqua),                intent(in)::ipqua
integer,dimension(1:nvfac),                intent(in):: ipf
integer,dimension(ngausf),                intent(out)::fgaus
!...Input real arrays
real*8, dimension(1:nq+3, 1:ngausf),       intent(out)::unkng
real*8, dimension(1:2, 1:2, 1:ngausf),     intent(out)::sigmg
real*8,dimension(1:ndegr,1:nq),            intent(in)::unkno
real*8,dimension(1:4),                     intent(in)::drhosgq
real*8,dimension(1:ngeel),                 intent(in)::geoel
real*8,dimension(1:nq+1),                  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2),                intent(in)::afvec

!...Local integer array
integer,dimension(ngausf,4)    ::fglgq
integer,dimension(ngausf*4)    ::sglgq

!...Local real array
real*8, dimension(1:2, 1:ngsfq)::posiq
real*8, dimension(1:mdegr,1:ngausf)::bq

!...Local integer
integer::ig,ideg

!...Local real number
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::rhovt,uvtx,vvtx,evtx,pvtx
real*8::rhomv,rhomc
real*8::dudr,duds,dvdr,dvds
real*8::dr, ds,rc,sc
real*8::eps
!
eps = 1.d-15
!
rc = geoel(1)
sc = geoel(2)

dr = 1.d0
ds = 1.d0
!
if(ngausf.eq.3)then
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

!...8-nodes quad
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

elseif(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
!
sglgq( 9) = 1; sglgq(10) = 2;
sglgq(11) = 2; sglgq(12) = 3;
sglgq(13) = 3; sglgq(14) = 4;
sglgq(15) = 4; sglgq(16) = 1;
!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = -0.4472135954999579392818d0; posiq(2, 9) =-1.d0;
posiq(1,10) =  0.4472135954999579392818d0; posiq(2,10) =-1.d0;

posiq(1,11) = 1.d0; posiq(2,11) = -0.4472135954999579392818d0;
posiq(1,12) = 1.d0; posiq(2,12) =  0.4472135954999579392818d0;

posiq(1,13) = 0.4472135954999579392818d0; posiq(2,13) = 1.d0;
posiq(1,14) =-0.4472135954999579392818d0; posiq(2,14) = 1.d0;

posiq(1,15) =-1.d0; posiq(2,15) = 0.4472135954999579392818d0;
posiq(1,16) =-1.d0; posiq(2,16) =-0.4472135954999579392818d0;
!
elseif(ngausf.eq.5)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;
!
sglgq(13) = 1; sglgq(14) = 2;
sglgq(15) = 2; sglgq(16) = 3;
sglgq(17) = 3; sglgq(18) = 4;
sglgq(19) = 4; sglgq(20) = 1;

!
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;

posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;

posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;

posiq(1,13) =-0.6546536707079771437983d0; posiq(2,13) =-1.d0;
posiq(1,14) = 0.6546536707079771437983d0; posiq(2,14) =-1.d0;

posiq(1,15) = 1.d0; posiq(2,15) =-0.6546536707079771437983d0;
posiq(1,16) = 1.d0; posiq(2,16) = 0.6546536707079771437983d0;

posiq(1,17) = 0.6546536707079771437983d0; posiq(2,17) = 1.d0;
posiq(1,18) =-0.6546536707079771437983d0; posiq(2,18) = 1.d0;

posiq(1,19) =-1.d0; posiq(2,19) = 0.6546536707079771437983d0;
posiq(1,20) =-1.d0; posiq(2,20) =-0.6546536707079771437983d0;

endif

!...Mass average value
rhoct = 1.d0/unkno(1, 1)
uctr  = unkno(1, 2)
vctr  = unkno(1, 3)
ectr  = unkno(1, 4)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )

!...Identify the fgaus for right cell
call getfg_glb2(ipqua, ipf, fgaus)
!
do ig =ngvtxf+1 ,ngausf
bq(1, ig) = 1.d0
bq(2, ig) = (posiq(1, fgaus(ig))-rc)/dr
bq(3, ig) = (posiq(2, fgaus(ig))-sc)/ds
!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21)
endif
enddo

!...Zero out unkng
unkng = 0.d0

do ig = ngvtxf+1 , ngausf
do ideg =1, mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq)*bq(ideg, ig)
enddo
enddo

!...Limiter...
do ig = ngvtxf+1, ngausf
!
rhovt = 1.d0/unkng(1, ig) + cdrho*drhosgq(sglgq(fgaus(ig)))
uvtx = unkng(2, ig)
vvtx = unkng(3, ig)
evtx = unkng(4, ig)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.6.and.geoel(10).gt.10.d0)then
!
rhomc = 1.d0/rhoct
!
!rhomv = rhomc + aflim(1)*(unkng(1, ig) - rhomc)
rhomv = rhomc + aflim(1)*(1.d0/rhovt - rhomc)
unkng(1, ig) = rhomv
rhovt = 1.d0/rhomv
!
dudr = afvec(1, 1)*unkno(2,2) +  afvec(1, 2)*unkno(2,3)
duds = afvec(1, 1)*unkno(3,2) +  afvec(1, 2)*unkno(3,3)
dvdr = afvec(2, 1)*unkno(2,2) +  afvec(2, 2)*unkno(2,3)
dvds = afvec(2, 1)*unkno(3,2) +  afvec(2, 2)*unkno(3,3)
!
uvtx = unkno(1,2)  + dudr*bq(2, ig) + duds*bq(3, ig)
vvtx = unkno(1,3)  + dvdr*bq(2, ig) + dvds*bq(3, ig)
!
pvtx = pctr + aflim(4)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig) = uvtx
unkng(3, ig) = vvtx
!
endif
!
unkng(5, ig) = pvtx
unkng(6, ig) = rhoct*sdctr
unkng(7, ig) = sdctr

!...Get stress tensor at nodes
sigmg(1, 1, ig) = -pvtx
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pvtx!
enddo

end subroutine getvar_face_gnvtx
!
!...Find the global no. of one face gauss point
!
subroutine getfg_glb2(ipqua, ipf, fgaus)
use constant
implicit none
integer, dimension(nvqua), intent(in):: ipqua
integer, dimension(nvfac), intent(in):: ipf
integer:: ig, nsum
integer, dimension(1:ngausf), intent(out)::fgaus
!
integer::fglgq(ngausf,4)

!...Give fglgq
if(ngausf.eq.4)then
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 10;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =11; fglgq(4, 2) = 12;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =13; fglgq(4, 3) = 14;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =15; fglgq(4, 4) = 16;
elseif(ngausf.eq.5)then

fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9; fglgq(4, 1) = 13;  fglgq(5, 1) = 14;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10; fglgq(4, 2) = 15;  fglgq(5, 2) = 16;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11; fglgq(4, 3) = 17;  fglgq(5, 3) = 18;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12; fglgq(4, 4) = 19;  fglgq(5, 4) = 20;
endif

!...Find the No. of face gauss points in fglgq...
nsum = 2*ngausf - 1
!
do ig = ngvtxf+1, ngausf
!
if(ipqua(5).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(1))then
fgaus(ig) = fglgq(ig,1)
else
fgaus(ig) = fglgq(nsum-ig,1)
endif
!
elseif(ipqua(6).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(2))then
fgaus(ig) = fglgq(ig,2)
else
fgaus(ig) = fglgq(nsum-ig,2)
endif
!
elseif(ipqua(7).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(3))then
fgaus(ig) = fglgq(ig,3)
else
fgaus(ig) = fglgq(nsum-ig,3)
endif
!
elseif(ipqua(8).eq.ipf(3))then
!
if(ipf(1).eq.ipqua(4))then
fgaus(ig) = fglgq(ig,4)
else
fgaus(ig) = fglgq(nsum-ig,4)
endif
!
endif

enddo
end subroutine getfg_glb2
!
!...Get the averaged-density within one sub-cell...
!
subroutine  getrhosubcell_sms(xcrho, ycrho, xpq,xpqi,unksgq, ielem)
use constant
implicit none
!...Input
integer, intent(in)::ielem
real*8,  intent(in)::xcrho,ycrho
real*8,dimension(1:2, 1:npqua), intent(in)::xpq
real*8,dimension(1:2, 1:npqua), intent(in)::xpqi
real*8,dimension(1:ndegr, 1:nq), intent(out)::unksgq
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, id, isg, ideg, iunk,isc

integer,dimension(4, 4)::ipqsg, ipqsc
!...Local real array
real*8,dimension(1:ndegr,1:4)::rhsel
real*8,dimension(1:nmatr, 1:4)::matin
real*8,dimension(1:2, 1:npqua)::xpqsg,xpqisg
real*8::b(ndegr)
real*8::unint(1)
real*8::m(ndegr, ndegr)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: posqc(2, 12)
real*8:: xpqc(2, 12),xpqic(2, 12)
real*8:: xpqsc(2, npqua), xpqisc(2, npqua)
!...Local real
real*8:: xsc, ysc, xsci, ysci
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco, volel,det
real*8:: dxdri,dxdsi,dydri,dydsi
real*8:: djacoi,voleli, rhoi
real*8:: xcel, ycel,xceli, yceli
real*8::c10
real*8::b2,b3,b4,b5,b6
real*8::bq22,bq33,bq23
real*8::f0,f1,f2,f3
real*8:: f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::masel,xgaus,ygaus
real*8::rcsg,scsg
real*8::bq(ndegr)
!
real*8,allocatable::mminv(:,:),binv(:),mmatr(:,:)
!-xxx-real contant
data c10 / 1.0d0 /
!
if(npoly==2) allocate(mminv(5,5), mmatr(5,5), binv(5))
!
ipqsg(1, 1) = 1; ipqsg(2, 1) = 5; ipqsg(3, 1) = 9; ipqsg(4, 1) = 8
ipqsg(1, 2) = 5; ipqsg(2, 2) = 2; ipqsg(3, 2) = 6; ipqsg(4, 2) = 9
ipqsg(1, 3) = 9; ipqsg(2, 3) = 6; ipqsg(3, 3) = 3; ipqsg(4, 3) = 7
ipqsg(1, 4) = 8; ipqsg(2, 4) = 9; ipqsg(3, 4) = 7; ipqsg(4, 4) = 4
!
posqc(1, 1)= -0.5d0; posqc(2, 1)= -1.d0;
posqc(1, 2)=  0.5d0; posqc(2, 2)= -1.d0;
posqc(1, 3)=  1.0d0; posqc(2, 3)= -.5d0;
posqc(1, 4)=  1.0d0; posqc(2, 4)=  .5d0;
posqc(1, 5)=  0.5d0; posqc(2, 5)=  1.d0;
posqc(1, 6)= -0.5d0; posqc(2, 6)=  1.d0;
posqc(1, 7)= -1.0d0; posqc(2, 7)=  .5d0;
posqc(1, 8)= -1.0d0; posqc(2, 8)= -.5d0;
posqc(1, 9)=  0.0d0; posqc(2, 9)= -.5d0;
posqc(1,10)=  0.5d0; posqc(2,10)=  0.d0;
posqc(1,11)=  0.0d0; posqc(2,11)=  .5d0;
posqc(1,12)= -0.5d0; posqc(2,12)=  0.d0;
!
ipqsc(1, 1)= 1; ipqsc(2, 1)= 9; ipqsc(3, 1)= 12; ipqsc(4, 1)= 8;
ipqsc(1, 2)= 2; ipqsc(2, 2)= 3; ipqsc(3, 2)= 10; ipqsc(4, 2)= 9;
ipqsc(1, 3)=10; ipqsc(2, 3)= 4; ipqsc(3, 3)=  5; ipqsc(4, 3)=11;
ipqsc(1, 4)=12; ipqsc(2, 4)=11; ipqsc(3, 4)=  6; ipqsc(4, 4)= 7;
!
call ruqope(2, ngausdq, posiq, weighq)

!...dr and ds...
dr = 1.d0
ds = 1.d0

!...Now configuration...
xpqsg(1,1:9) =  xpq(1, 1:9)
xpqsg(2,1:9) =  xpq(2, 1:9)

!...Initial configuration...
xpqisg(1,1:9) =  xpqi(1, 1:9)
xpqisg(2,1:9) =  xpqi(2, 1:9)

!...Get the 12 high-order curved nodes
do isc = 1, 12
r = posqc(1,isc)
s = posqc(2,isc)

!...  shape function
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
xsc = 0.d0
ysc = 0.d0
!
do ishp = 1, npqua
xsc = xsc + shpq(ishp)*xpqsg(1,ishp)
ysc = ysc + shpq(ishp)*xpqsg(2,ishp)
enddo
!
xsci = 0.d0
ysci = 0.d0
!
do ishp = 1, npqua
xsci = xsci + shpq(ishp)*xpqisg(1,ishp)
ysci = ysci + shpq(ishp)*xpqisg(2,ishp)
enddo
!
xpqc(1, isc) = xsc
xpqc(2, isc) = ysc
!
xpqic(1, isc) = xsci
xpqic(2, isc) = ysci
enddo
!
!print*,'test getrhosubg4'
!
!...Part I: Get cell averaged r and s for one subgrid...
!

!
!...Part II: Get high-order terms for basis function...
!

!
!...Part II: Get mass matrix for four subgrids...
!
do isg = 1, 4
!...Coordinates
xpqsc(1, 1:4) = xpqsg(1, ipqsg(1:4, isg))
xpqsc(2, 1:4) = xpqsg(2, ipqsg(1:4, isg))

xpqsc(1, 5:8) = xpqc(1, ipqsc(1:4, isg))
xpqsc(2, 5:8) = xpqc(2, ipqsc(1:4, isg))
!
xpqsc(1:2,9) = -0.25d0*(xpqsc(1:2,1) + xpqsc(1:2,2) + xpqsc(1:2,3) + xpqsc(1:2,4)) +&
0.5d0*(xpqsc(1:2,5) + xpqsc(1:2,6) + xpqsc(1:2,7) + xpqsc(1:2,8))

!...Initialze parameters...
if(npoly.eq.2)  mmatr = 0.d0
f0 = 0.d0
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqsc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqsc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqsc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqsc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
f0 = f0 + djaco
!
!if(ielem.eq.1) print*,'djaco',f0,wi,dxdr,dyds,dydr,dxds,xpqsc(1:2,:)
enddo

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

do isg = 1, 4

!...Initial coordinates
xpqisc(1, 1:4) = xpqisg(1, ipqsg(1:4, isg))
xpqisc(2, 1:4) = xpqisg(2, ipqsg(1:4, isg))

xpqisc(1, 5:8) = xpqic(1, ipqsc(1:4, isg))
xpqisc(2, 5:8) = xpqic(2, ipqsc(1:4, isg))
!
xpqisc(1:2,9) = -0.25d0*(xpqisc(1:2,1) + xpqisc(1:2,2) + xpqisc(1:2,3) + xpqisc(1:2,4)) +&
0.5d0*(xpqisc(1:2,5) + xpqisc(1:2,6) + xpqisc(1:2,7) + xpqisc(1:2,8))
!
do ig =1,ngausdq
!
r = posiq(1,ig)
s = posiq(2,ig)
wi = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)
!
!
dsprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dsprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dsprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dsprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dsprq(5) =   -r*s*(s-c10)
dsprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dsprq(7) =   -r*s*(s+c10)
dsprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dsprq(9) = -2.0d0*r*(1.d0-s*s)
!
dspsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dspsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dspsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dspsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dspsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dspsq(6) =   -r*s*(r+c10)
dspsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dspsq(8) = -r*s*(r-c10)
dspsq(9) = -2.0d0*s*(1.d0-r*r)
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, npqua
dxdr = dxdr + dsprq(ishp)*xpqisc(1,ishp)
dxds = dxds + dspsq(ishp)*xpqisc(1,ishp)

dydr = dydr + dsprq(ishp)*xpqisc(2,ishp)
dyds = dyds + dspsq(ishp)*xpqisc(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
b(1) = 1.d0

!...Get the initial density
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpqisc(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqisc(2,ishp)
enddo

!...Get the initial density
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)
!
rhsel(1, isg)=rhsel(1, isg) + rhoi*b(1)*djaco
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
end subroutine  getrhosubcell_sms
