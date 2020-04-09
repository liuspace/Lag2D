!
!...subroutine: least square for FVM through loop over faces...
!
subroutine getgraduls_lag(bface,intfac, iptri, ipqua, unkno ,cocent, coord)
use constant
implicit none
integer*4,dimension(1:nbfai,1:nbfac),      intent(in)::bface
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(inout)::unkno
real*8,dimension(1:3,1:nsize), intent(out)::cocent
real*8,dimension(1:ndimn,1:npoin), intent(in) ::coord
integer*4,dimension(1:nifai,1:nafac), intent(in) ::intfac
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8:: rhom, rhoc, uctr,vctr, pctr, ectr

integer*4::iel,ier,ifa,ie
real*8::c00
real*8::dx,dy,weigh
real*8::dxwei,dywei
real*8::du1we,du2we,du3we,du4we
real*8::rxx,ryy,rxy
real*8::rco1,rco2,rco3,rco4
real*8::ru1x,ru2x,ru3x,ru4x
real*8::ru1y,ru2y,ru3y,ru4y
real*8::rjac1
real*8::x1,y1,x2,y2,x0,y0
real*8::la,lb,lc,abp2,abm2,xg,yg
integer::ip1,ip2
real*8::unctr(1:nq+1, nsize)
real*8::fnx, fny, lnx, lx, ly, xcl, ycl, xcr, ycr, vdon,delty, deltx
real*8:: xf(2, 2)
real*8,dimension(1:2,1:nq,1:ncell)::gradu
!
real*8, allocatable::dmets(:,:)
!
!...  This sub computes the gradient using the least-square approach
!
allocate(dmets(1:3,1:nsize))
!
!...Recalculate geometry parameter..
!
call getgeoel_lag_fvm(iptri, ipqua, cocent, coord)
!
c00 = 0.0
dmets(:,:) = 0.0
!
gradu = 0.0
!
!...Give the cell-averaged density, veloccity, pressure....
!
do ie = 1, ncell
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoc = 1.d0/rhom
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))

!
unctr(1, ie) = rhom;
unctr(2, ie) = uctr;
unctr(3, ie) = vctr;
unctr(4, ie) = pctr;
unctr(5, ie) = ectr;
!
enddo
!
do ifa = 1, nbfac
!
if(bface(3, ifa).eq.21) cycle
!
!... end-elements of the faces
!
iel   = intfac(1,ifa)
ier   = intfac(2,ifa)
!
xf(1:2, 1) = coord(1:2, intfac(3, ifa))
xf(1:2, 2) = coord(1:2, intfac(4, ifa))
!
delty= xf(2, 2) - xf(2, 1)
deltx=-(xf(1, 2) - xf(1, 1))
!
fnx = delty/sqrt(delty**2 + deltx**2)
fny = deltx/sqrt(delty**2 + deltx**2)
!
xcl = cocent(1, iel)
ycl = cocent(2, iel)
!
lnx = xcl*fnx + ycl*fny
!
xcr =xcl - 2.d0*lnx*fnx
ycr =ycl - 2.d0*lnx*fny
!
cocent(1, ier) = xcr
cocent(2, ier) = ycr
!
unctr(1, ier) = unctr(1, iel)
unctr(4, ier) = unctr(4, iel)
unctr(5, ier) = unctr(5, iel)
!
uctr = unctr(2, iel)
vctr = unctr(3, iel)
!
vdon = uctr*fnx + vctr*fny
!
unctr(2, ier) = uctr - 2.d0*vdon*fnx
unctr(3, ier) = vctr - 2.d0*vdon*fny
!
!if(ifa.ge.1.and.ifa.le.30) print*,'ifa',ifa,iel,ier,xcl,ycl,xcr,ycr,unctr(2:3, ier),unctr(2:3, iel)
!
enddo
!
!print*,'center',cocent(1:2, 9:10)
!
!...Start LS
!
do 1400 ifa =11, nafac
!
if(ifa.le.nbfac)then
if(bface(3, ifa).eq.21) cycle
endif
!
!... end-elements of the faces
!
iel   = intfac(1,ifa)
ier   = intfac(2,ifa)
!
dx  = cocent(1,ier) - cocent(1,iel)
dy  = cocent(2,ier) - cocent(2,iel)
!
!...  weighting for this edge
!
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
du1we = weigh*(unctr(1,ier) - unctr(1,iel))
du2we = weigh*(unctr(2,ier) - unctr(2,iel))
du3we = weigh*(unctr(3,ier) - unctr(3,iel))
du4we = weigh*(unctr(4,ier) - unctr(4,iel))
!
!if(iel.eq.10.or.ier.eq.10) print*,'iel',iel,ier,dxwei,dywei
!
gradu(1,1,iel) = gradu(1,1,iel) + dxwei*du1we
gradu(1,2,iel) = gradu(1,2,iel) + dxwei*du2we
gradu(1,3,iel) = gradu(1,3,iel) + dxwei*du3we
gradu(1,4,iel) = gradu(1,4,iel) + dxwei*du4we
!
gradu(2,1,iel) = gradu(2,1,iel) + dywei*du1we
gradu(2,2,iel) = gradu(2,2,iel) + dywei*du2we
gradu(2,3,iel) = gradu(2,3,iel) + dywei*du3we
gradu(2,4,iel) = gradu(2,4,iel) + dywei*du4we
!
!  excluding the ghost element gradient
if(ifa.gt.nbfac) then
!
gradu(1,1,ier) = gradu(1,1,ier) + dxwei*du1we
gradu(1,2,ier) = gradu(1,2,ier) + dxwei*du2we
gradu(1,3,ier) = gradu(1,3,ier) + dxwei*du3we
gradu(1,4,ier) = gradu(1,4,ier) + dxwei*du4we
!
gradu(2,1,ier) = gradu(2,1,ier) + dywei*du1we
gradu(2,2,ier) = gradu(2,2,ier) + dywei*du2we
gradu(2,3,ier) = gradu(2,3,ier) + dywei*du3we
gradu(2,4,ier) = gradu(2,4,ier) + dywei*du4we
!
endif
!
1400 continue
!
!...
!
do 2000 ie  = 1, ncell
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
ru1x        = gradu(1,1,ie)
ru1y        = gradu(2,1,ie)

ru2x        = gradu(1,2,ie)
ru2y        = gradu(2,2,ie)
!
ru3x        = gradu(1,3,ie)
ru3y        = gradu(2,3,ie)
!
ru4x        = gradu(1,4,ie)
ru4y        = gradu(2,4,ie)
!
gradu(1,1,ie) = (ru1x*rco1 + ru1y*rco2 )*rjac1
gradu(2,1,ie) = (ru1x*rco3 + ru1y*rco4 )*rjac1
!
gradu(1,2,ie) = (ru2x*rco1 + ru2y*rco2 )*rjac1
gradu(2,2,ie) = (ru2x*rco3 + ru2y*rco4 )*rjac1
!
gradu(1,3,ie) = (ru3x*rco1 + ru3y*rco2 )*rjac1
gradu(2,3,ie) = (ru3x*rco3 + ru3y*rco4 )*rjac1
!
gradu(1,4,ie) = (ru4x*rco1 + ru4y*rco2 )*rjac1
gradu(2,4,ie) = (ru4x*rco3 + ru4y*rco4 )*rjac1
!
!if(ie.eq.10) print*,'grad',dmets(1:3, ie)
!
2000 continue
!
!...Give unkno gradu
!
do ie = 1, ncell
!
unkno(2:3, 1, ie) = gradu(1:2,1,ie)
unkno(2:3, 2, ie) = gradu(1:2,2,ie)
unkno(2:3, 3, ie) = gradu(1:2,3,ie)
unkno(2:3, 4, ie) = gradu(1:2,4,ie)
!
enddo
!
deallocate(dmets)
return
end subroutine getgraduls_lag
!
!...Find mass center in geoel for curved cell...
!
subroutine getgeoel_lag_fvm(iptri, ipqua, geoel, coord)
use constant
implicit none
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:3, 1:nsize), intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
!
!...local array...
!
real*8,dimension(1:2, 1:3)::xpin
real*8,dimension(1:2, 1:nptri)::xp
real*8,dimension(1:2, 1:npqua)::xpq
real*8,dimension(1:nptri)::shp, dspr, dsps
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weight(ngausd_geo), posit(2, ngausd_geo)
real*8:: weighq(ngausd_geoq), posiq(2, ngausd_geoq)
!...local real number
real*8:: dxdr,dxds,dydr,dyds
real*8:: eps,c00,c10,c05,c20
real*8:: r, s, djaco, volel, masel
real*8:: wi, xcent, ycent, xgaus, ygaus
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
if(ncurv==0)then
xp(1, 1:3) = coord(1, iptri(1:3,ie))
xp(2, 1:3) = coord(2, iptri(1:3,ie))
!
xp(1:2,4) = 0.5d0*(xp(1:2,1)+xp(1:2,2))
xp(1:2,5) = 0.5d0*(xp(1:2,2)+xp(1:2,3))
xp(1:2,6) = 0.5d0*(xp(1:2,1)+xp(1:2,3))
elseif(ncurv==1)then
xp(1, 1:nptri) = coord(1,iptri(1:nptri, ie))
xp(2, 1:nptri) = coord(2,iptri(1:nptri, ie))
endif!
!...
xcent = 0.d0
ycent = 0.d0
volel = 0.d0
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
do ishp = 1, nptri
dxdr = dxdr + dspr(ishp)*xp(1,ishp)
dxds = dxds + dsps(ishp)*xp(1,ishp)

dydr = dydr + dspr(ishp)*xp(2,ishp)
dyds = dyds + dsps(ishp)*xp(2,ishp)
enddo
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nptri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
volel = volel + djaco
!
!...Density distribution for different cases...
!
xcent = xcent + xgaus*djaco !...mass center at the reference element...
ycent = ycent + ygaus*djaco
!
enddo
!
geoel(1, ielem) = xcent/volel
geoel(2, ielem) = ycent/volel
geoel(2, ielem) = volel
!
!...Other components are same as those of Eulerian framework...
!
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
!print*,'ipqua',ipqua(1:8,1)
!
!...global numbering
!
ielem = ie + ntria
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
endif
!
!...
xcent = 0.d0
ycent = 0.d0
volel = 0.d0
!
do igaus =1,ngausd_geoq
!
r  = posiq(1,igaus)
s  = posiq(2,igaus)
wi  = weighq(igaus)
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
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
volel = volel + djaco
!
xcent = xcent + xgaus*djaco !...mass center at the reference element...
ycent = ycent + ygaus*djaco
!
enddo
!
geoel(1, ielem) = xcent/volel
geoel(2, ielem) = ycent/volel
geoel(3, ielem) = volel
!
enddo !...(1)ie = 1,nelem
!print*,'ipqua',ipqua(1:8,1)
end subroutine getgeoel_lag_fvm
!
!...subroutine: Calculate the nodal velocity U_p^* (mass center) for hybrid mesh with general Riemann solver...
!
subroutine getndvelo_lagfvm_curv(gflag,gelag,cocent,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord,unkno,ustar, fstar, fstarq, aflim, afvec, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3, 1:ngelg, 1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:nsize),               intent(in)::cocent
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri, 1:ntria),  intent(out)::fstar !...Riemann forces
real*8,dimension(1:ndimn,1:4,1:nvqua, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec

integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem,iloop
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)

!...local real array
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8::munaci(2, 2)
real*8,dimension(1:ndimn,1:nvqua,1:nquad):: vnulq
real*8,  dimension(1:nquad)::gqdmp
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: munacn(:,:,:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: usold(:,:), munacu(:,:), snsigm(:,:)
real*8,allocatable:: munaclt(:,:,:,:,:),munaclq(:,:,:,:,:)
real*8,allocatable:: snsigmlt(:,:,:,:), munault(:,:,:,:)
real*8,allocatable:: snsigmlq(:,:,:,:), munaulq(:,:,:,:)

!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:2, 1:2, 1:npoin))
allocate (usold(1:ndimn, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munaclt(1:2, 1:2, 1:4, 1:nvtri, 1:ntria), munault(1:ndimn, 1:4, 1:nvtri,  1:ntria),&
snsigmlt(1:ndimn, 1:4,  1:nvtri,  1:ntria))
allocate (munaclq(1:2, 1:2, 1:4, 1:nvqua, 1:nquad), munaulq(1:ndimn, 1:4, 1:nvqua,  1:nquad),&
snsigmlq(1:ndimn, 1:4,  1:nvqua,  1:nquad))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
!
!...Zero out vlave
!
vlave = 0.d0
usold = ustar
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
!
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
indnd(ipf(1:nvfac)) = 1
endif
enddo
!
!...Get averaged velocity at nodes...
!

!
do 950 ifa = 1 , nbfac
!
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

endif!
950 enddo
!
call getnullve_dampc(ipqua, geoel, cocent, coord, ustar, unkno, gqdmp)
!
vnulq = 0.d0
!
do iloop= 1,4!10
!
vlave = ustar
!
!call getnullmd_quadlp(ipqua, geoel, vnulq, ustar, coord, gqdmp)
!
!vnulq = 0.d0
!
!...Zero out munacn
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
if(nquad.gt.0) call getriem_quadfvm_matrix(ipqua, geoel, cocent, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, aflim, afvec, vnulq)
!
!...4.1: Update the Riemann forces at every node...
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
!
detma = munacn(1, 1, ipoin)*munacn(2, 2, ipoin) - munacn(2, 1, ipoin)*munacn(1, 2, ipoin)
munaci(1, 1) = munacn(2, 2, ipoin)/detma
munaci(1, 2) =-munacn(1, 2, ipoin)/detma
munaci(2, 1) =-munacn(2, 1, ipoin)/detma
munaci(2, 2) = munacn(1, 1, ipoin)/detma
!
rhsu1 = munacu(1, ipoin) - snsigm(1, ipoin) !+ fpres(1, ipoin)
rhsu2 = munacu(2, ipoin) - snsigm(2, ipoin) !+ fpres(2, ipoin)
!
ustar(1, ipoin) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ustar(2, ipoin) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!
! if(ipoin.eq.70) print*,ustar(1:2,ipoin),munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
endif
enddo
!
!call getvelo_mpt_curvfvm(ustar,geoel,cocent,gelag,gelagq,intfac,iptri,ipqua,coord,unkno,indnd,&
! aflim, afvec, vlave, vnulq)

!call  getvelo_mpt_curvfvm2(ustar,geoel,cocent,gelag,gelagq,intfac,iptri,ipqua,coord,unkno,indnd,&
!aflim, afvec, vlave, vnulq,munacn)
!
! print*,'ustar--',cos(pi/24.d0*10.d0),sin(pi/24.d0*10.d0)
!
do 900 ifa = 1 , nbfac
!
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

!...Specify the boundary velocity
if(ncase.eq.13)then !...Saltzman
ustar(2,ipf(1:nvfac)) = 0.d0
ustar(1,ipf(1:nvfac)) = 1.d0
endif
endif
!
if(ncase.eq.1)then
!
!...impose exact solution along the Boundary
!
!ustar(1, ipf(1:nvfac)) = sin(pi*coord(1,ipf(1:nvfac)))*cos(pi*coord(2,ipf(1:nvfac)))
!ustar(2, ipf(1:nvfac)) =-cos(pi*coord(1,ipf(1:nvfac)))*sin(pi*coord(2,ipf(1:nvfac)))
endif
!
!if(ipf(1).eq.100) print*,'ipf',ipf(1),ustar(1:2, ipf(1))
!
900 enddo
!
!
!
do 910 ifa = 1 ,-nafac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
ustar(1, ipf(3)) = 0.5d0*(ustar(1, ipf(1))+ustar(1, ipf(2)))
ustar(2, ipf(3)) = 0.5d0*(ustar(2, ipf(1))+ustar(2, ipf(2)))
910 enddo
!
call getnullmd_quadlp2(ipqua, geoel, vnulq, ustar, coord, gqdmp)
!
enddo !loop
!
!...Quads
!
do ie = 1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
!...shape functions
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv = 1, 4
!
do ifa =1, 4
fstarq(1, ifa, iv, ie) = snsigmlq(1, ifa, iv, ie) + &
munaclq(1, 1, ifa, iv, ie)*ustar(1, ipq(iv))+&
munaclq(2, 1, ifa, iv, ie)*ustar(2, ipq(iv)) - munaulq(1, ifa, iv, ie)
fstarq(2, ifa, iv, ie) = snsigmlq(2, ifa, iv, ie) + &
munaclq(1, 2,ifa, iv, ie)*ustar(1, ipq(iv))+&
munaclq(2, 2, ifa, iv, ie)*ustar(2, ipq(iv)) - munaulq(2, ifa, iv, ie)
!
enddo
!
!if(ie==21) print*,'fire',iv,ie,snsigmlq(1,1:2,iv,ie),munaclq(1:2,iv,ie),ipq(:),munaulq(1,1:2,iv,ie)
!if(ie.eq.21) print*,'fstarq',iv,ie,fstarq(1,1:2,iv,ie)
!
enddo
!
do iv = 5, 8
!
do ifa =1, 2
fstarq(1, ifa, iv, ie) = snsigmlq(1, ifa, iv, ie) + &
munaclq(1, 1, ifa, iv, ie)*ustar(1, ipq(iv))+&
munaclq(2, 1, ifa, iv, ie)*ustar(2, ipq(iv)) - munaulq(1, ifa, iv, ie)
fstarq(2, ifa, iv, ie) = snsigmlq(2, ifa, iv, ie) + &
munaclq(1, 2,ifa, iv, ie)*ustar(1, ipq(iv))+&
munaclq(2, 2, ifa, iv, ie)*ustar(2, ipq(iv)) - munaulq(2, ifa, iv, ie)
!
enddo
!
!do ifa =1, 2
!fstarq(1:2, ifa, 5, ie) = 0.5d0*(fstarq(1:2, 4, 1, ie) + fstarq(1:2, 3, 2, ie))
!fstarq(1:2, ifa, 6, ie) = 0.5d0*(fstarq(1:2, 4, 2, ie) + fstarq(1:2, 3, 3, ie))
!fstarq(1:2, ifa, 7, ie) = 0.5d0*(fstarq(1:2, 4, 3, ie) + fstarq(1:2, 3, 4, ie))
!fstarq(1:2, ifa, 8, ie) = 0.5d0*(fstarq(1:2, 4, 4, ie) + fstarq(1:2, 3, 1, ie))
!
!enddo
!if(ie==21) print*,'fire',iv,ie,snsigmlq(1,1:2,iv,ie),munaclq(1:2,iv,ie),ipq(:),munaulq(1,1:2,iv,ie)
!if(ie.eq.21) print*,'fstarq',iv,ie,fstarq(1,1:2,iv,ie)
!
enddo
enddo

!
deallocate (munacn, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lagfvm_curv
!
!...subroutine: Calculate the Riemann input for hybrid quad grids general Riemann solver....
!
subroutine getriem_quadfvm_matrix(ipqua, geoel, cocent, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, aflim, afvec, vnulq)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:nsize),              intent(in)::cocent
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave, coord
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:ndimn,1:nvqua,1:nquad),    intent(in):: vnulq

!
real*8, dimension(1:2, 1:2, 1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
!
real*8, dimension(1:2, 1:2, 1:4,  1:nvqua, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:4,  1:nvqua, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:4,  1:nvqua, 1:nquad), intent(out)::snsigmlq
!...Local integer
integer::ie, ideg, ielem, ifa, iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:3, 1:nvqua)::bq
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8::aujmp(1:3, 1:nvqua)
real*8::vnorm(1:3, 1:4, 1:nvqua)
real*8::sigma(1:2, 1:2, 1:nvqua)
real*8,dimension(1:nvqua)::murie
real*8,dimension(1:nvqua):: xvq, yvq
!...arraies for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8:: rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx, rhovt
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
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
rc= cocent(1, ielem) !...mass center...
sc= cocent(2, ielem)
!
xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)
bq(3, iv) = (yvq(iv)-sc)
enddo
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) =  gelagq(1:3,  4, ie);   vnorm(1:3, 2, 1) = gelagq(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 3, 1) =  gelagq(1:3, 12, ie);   vnorm(1:3, 4, 1) = gelagq(1:3, 5, ie) !...For point ip(1)

vnorm(1:3, 1, 2) = gelagq(1:3, 1, ie); vnorm(1:3, 2, 2) = gelagq(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelagq(1:3, 6, ie); vnorm(1:3, 4, 2) = gelagq(1:3, 7, ie) !...For point ip(2)

vnorm(1:3, 1, 3) = gelagq(1:3, 2, ie); vnorm(1:3, 2, 3) = gelagq(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelagq(1:3, 8, ie); vnorm(1:3, 4, 3) = gelagq(1:3, 9, ie) !...For point ip(3)

vnorm(1:3, 1, 4) = gelagq(1:3, 3, ie); vnorm(1:3, 2, 4) = gelagq(1:3, 4, ie) !...For point ip(4)
vnorm(1:3, 3, 4) = gelagq(1:3, 10, ie); vnorm(1:3, 4, 4) = gelagq(1:3, 11, ie) !...For point ip(4)

vnorm(1:3, 1, 5) = gelagq(1:3, 5, ie); vnorm(1:3, 2, 5) = gelagq(1:3, 6, ie) !...For point ip(5)
vnorm(1:3, 1, 6) = gelagq(1:3, 7, ie); vnorm(1:3, 2, 6) = gelagq(1:3, 8, ie) !...For point ip(6)
vnorm(1:3, 1, 7) = gelagq(1:3, 9, ie); vnorm(1:3, 2, 7) = gelagq(1:3, 10, ie) !...For point ip(7)
vnorm(1:3, 1, 8) = gelagq(1:3, 11, ie); vnorm(1:3, 2, 8) = gelagq(1:3, 12, ie) !...For point ip(8)
!
!...cell averaged value...
!
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unknvq = 0.d0
!
do iv   = 1,nvqua
!
!...Interpolation at the vertex...
!
do ideg = 1,mdegr
unknvq(1:nq-1, iv) = unknvq(1:nq-1, iv) + unkno(ideg,1:nq-1,ielem)*bq(ideg, iv)
enddo
!
!...Pressure...
!
unknvq(4, iv) = pctr  + unkno(2,4,ielem)*bq(2, iv) + unkno(3,4,ielem)*bq(3, iv)
!
!...Remove null mode...
!
unknvq(2:3, iv)=unknvq(2:3, iv)-vnulq(1:2,iv,ie)
!
rhovt  = 1.d0/unknvq(1, iv)
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
pvtx = unknvq(4, iv)
!
!...nlimi.eq.6:
!
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvq(1, iv) - rhomc)
rhovt = 1.d0/rhomv
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
vvtx = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!
!...updtae unknv(2:3,:)
!
unknvq(2, iv) = uvtx-vnulq(1,iv,ie)
unknvq(3 ,iv) = vvtx-vnulq(2,iv,ie)
!
!if(ielem.eq.1119.or.ielem.eq.1120)print*,'aflim',aflim(2:3,ielem),iv,uvtx,vvtx,uvtx**2+vvtx**2
!
endif
!
!...Get stress tensor at nodes
!
sigma(1, 1, iv) = -pvtx
sigma(1, 2, iv) = 0.d0
sigma(2, 1, iv) = 0.d0
sigma(2, 2, iv) = -pvtx!
!
!if(ip(iv).eq.8) print*,'velocity 8',ie, rho, uvtx,vvtx,vlave(1:2, ip(iv))
!
!...Get the a_c (unit vector)
aujmp(1:2, iv) = vlave(1:2, ipq(iv)) - unknvq(2:3, iv)
!if(ie==22) print*,'adjumpxxx22', vlave(1:2, ip(iv)) , unknv(2:3, iv), ip(iv)
acnx = aujmp(1, iv)
acny = aujmp(2, iv)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, iv) = 1.e-11!0.d0;
!print*,'point are reset', ip(iv)
else
aujmp(1:2, iv) = aujmp(1:2, iv)/sqrt(acnx**2 + acny**2)
endif
!
aujmp(3, iv) = sqrt(acnx**2 + acny**2)
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the variables at the center...
!
rhoct = 1.d0/unkno(1, 1, ielem)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ielem)
vctr  = unkno(1, 3, ielem)
ectr  = unkno(1, 4, ielem)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
aujmp(3,:)=aujmp(3,:)/sdctr
!
!...Get impedence coefficient...
!
do iv   = 1, nvqua
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ipq(iv).eq.738) print*,'murie22',murie(iv), sdctr,rhoct,deltu,vlave(1:2, ipq(iv)),unknvq(2:3, iv),dux,duy,ie
! if(ipq(iv).eq.840) then
!         print*,'murie22ang',unknvq(2,iv)**2+unknvq(3,iv)**2,ie,iv,vlave(1,ipq(iv))*unknvq(2,iv)+vlave(2,ipq(iv))*unknvq(3,iv)
! endif

enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
!!
do iv  = 1, 4
do ifa = 1, 4 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!...Call Riemann solver...
!
call getriecoef_matrixnew(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
!unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
if(ifa.le.2)then
!

!
elseif(ifa.ge.3)then
!
munacn(1:2, 1, ipq(iv)) = munacn(1:2, 1, ipq(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipq(iv)) = munacn(1:2, 2, ipq(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipq(iv)) = munacu(1:2, ipq(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipq(iv)) = snsigm(1:2, ipq(iv)) + snsigm_rie(1:2)
!
!if(ipq(iv).eq.104) print*,'p36 muacn(vv) post',ie,murie(iv),munacn_rie(1:2,1:2),&
!                         munacu_rie(1:2),unknvq(2:3, iv),vnorm(1:3, ifa, iv)
!
endif
!
!   if(ipq(iv).eq.738) print*,'p36 muacn(vv) postu',ie,munacu(1:2,ipq(iv)),munacu(1:2,738)/sqrt(munacu(1,738)**2+munacu(2,738)**2)
!
!...Local variable...
!
munaclq(1:2, 1, ifa, iv, ie) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, ifa, iv, ie) =  munacn_rie(1:2, 2)
!
munaulq(1:2, ifa, iv, ie) =  munacu_rie(1:2)
!
snsigmlq(1:2, ifa, iv, ie)=  snsigm_rie(1:2)
!
enddo

enddo
!!
do iv  = 5, 8
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!...Call Riemann solver...
!!
call getriecoef_matrixnew(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
!unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipq(iv)) = munacn(1:2, 1, ipq(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipq(iv)) = munacn(1:2, 2, ipq(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipq(iv)) = munacu(1:2, ipq(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipq(iv)) = snsigm(1:2, ipq(iv)) + snsigm_rie(1:2)
!
!if(ipq(iv).eq.2) print*,'p36 muacn(vv) post',ie,murie(iv),munacn_rie(1:2,1:2),&
!                         munacu_rie(1:2),unknvq(2:3, iv)
!   if(ipq(iv).eq.738) print*,'p36 muacn(vv) postu',ie,munacu(1:2,ipq(iv)),munacu(1:2,738)/sqrt(munacu(1,738)**2+munacu(2,738)**2)
!
!...Local variable...
!
munaclq(1:2, 1, ifa, iv, ie) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, ifa, iv, ie) =  munacn_rie(1:2, 2)
!
munaulq(1:2, ifa, iv, ie) =  munacu_rie(1:2)
!
snsigmlq(1:2, ifa, iv, ie)=  snsigm_rie(1:2)
!
enddo

enddo

350 enddo  !...(1)ie = 1,nelem!

end subroutine getriem_quadfvm_matrix
!
!...Subroutine for barth limiter based on vertex for primitive variables on curv quads....
!
subroutine barthlimit_lagfvm_curvquad(geoel, cocent, coord, coold, ustar, unkno, ipqua, bface,intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:3,1:nsize),                 intent(in)::cocent
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord, coold
real*8,dimension(1:nq+1, 1:nsize),             intent(inout)::aflim
integer,  dimension(1:nvqua,1:nquad),          intent(in):: ipqua
integer*4,dimension(1:nbfai,1:nbfac),           intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
real*8, dimension(1:nq+2, 1:npoin),           intent(in):: unmax, unmin
real*8,dimension(1:2, 1:2, 1:nsize),          intent(inout)::afvec
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
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
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:(2+nvfac), ifa)) = 1

enddo
!
eps = 1.e-6
!
!...Store the maximum and minimum values surrounding one cell...
!
do ie = 1, nquad
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)
do iq=1, nq+2
!
unmax_new(iq, ielem) = maxval(unmax(iq, ipq(1:nvqua)))
unmin_new(iq, ielem) = minval(unmin(iq, ipq(1:nvqua)))
enddo
enddo
!
!...Part 2: Impose limiter
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= cocent(1, ielem) !...mass center...
sc= cocent(2, ielem)
!
xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)
bq(3, iv) = (yvq(iv)-sc)
enddo
!
!...Cell averaged...
!
rhom = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq-1, iv) = unknvq(1:nq-1, iv) + unkno(ideg,1:nq-1,ielem)*bq(ideg, iv)
enddo
!
!...Pressure...
!
unknvq(4, iv) = pctr  + unkno(2,4,ielem)*bq(2, iv) + unkno(3,4,ielem)*bq(3, iv)
!
rhov = 1.d0/unknvq(1, iv)
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
unknvq(1, iv) = 1.d0/rhov
!
enddo
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq)  = pctr
!
do iv = 1, nvqua
do iq = 1, nq
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
!call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)
!
!call barthfctdiv(unkno(:,:,ielem),unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)

!
alfa(iq, iv) = afbar
enddo
!
!...Special treatment of boundary...
!
if(indbd(ipq(iv)).eq.1)then

if(ncase.eq.1)then
if(coord(1, ipq(iv)).lt.1.d-6.or.abs(coord(1, ipq(iv))-1.d0).lt.1.d-6) then
alfa(:, iv) = 1.d0
endif
!
if(coord(2, ipq(iv)).lt.1.d-6.or.abs(coord(2, ipq(iv))-1.d0).lt.1.d-6) then
alfa(:, iv) = 1.d0
endif
endif
endif
!
enddo
!
!...Choose the minimum one
!
do iq = 1,nq
aflim(iq, ielem) = minval(alfa(iq, 1:8))
!aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
enddo
!
enddo
!
!...Part 2.1: Impose symmetry preserving limiter for vector...
!
!call barthlimit_symprefvm_curvquad(geoel, cocent, coord, ustar, unkno, ipqua, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)

call barthlimit_symprefvm_qcquater(geoel, cocent,coord, ustar, unkno, ipqua, bface,intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
!
!...Degenarate to non-symmetric limiter...
!
do ie = 1, -nquad
!
ielem = ie + ntria
!
afvec(1, 1, ielem) = aflim(2, ielem)
afvec(1, 2, ielem) = 0.d0

afvec(2, 1, ielem) = 0.d0
afvec(2, 2, ielem) = aflim(3, ielem)
!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lagfvm_curvquad
!
!...Symmetry preserving techniques curved cell...
!
subroutine barthlimit_symprefvm_qc(geoel, cocent,coord, ustar, unkno, ipqua, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:3,1:nsize),                 intent(in)::cocent
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),           intent(in):: ipqua
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp
integer:: ielem, jelem, istor
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:nq+1) :: ummax, ummin
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
real*8,dimension(1:nsize)   ::afdu, afma
real*8,dimension(1:nsize) :: undu
real*8,dimension(1:nq+1, 1:nsize)   ::aflim
real*8::unmpx(1:nq), unmpn(1:nq)
!

real*8::eps,c00,c05,c10,c20
real*8:: uloca, vloca, ucjel, vcjel
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct, rhom
real*8:: afbar
real*8:: dudr, duds, dvdr, dvds
real*8:: etcx, etcy, xicx, xicy, fixi, fiet
real*8:: lamda1, lamda2
real*8:: mapd, mapt, matra, matrb, matrc, matrd, umap, vmap
real*8:: jf11, jf12, jf21, jf22
real*8:: fi11, fi12, fi21, fi22
real*8:: dxdr, dxds, dydr, dyds
real*8:: rm, rp, sm, sp, r, s
real*8:: drdx, drdy, dsdx, dsdy, jacom
real*8:: dudx, dudy, dvdx, dvdy, lmat1, lmat2
real*8:: dudxy, dumax, dumin
real*8:: rhomc,sdctr,macel,delu,volel
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:(2+nvfac), ifa)) = 1
enddo
!
eps = 1.e-6
mapmt = 0.d0
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
dudr = unkno(2, 2, ielem)
duds = unkno(3, 2, ielem)
dvdr = unkno(2, 3, ielem)
dvds = unkno(3, 3, ielem)
!
dudx = dudr
dudy = duds
dvdx = dvdr
dvdy = dvds
!
!...1st invariant
!
undu(ielem) = dudx +dvdy
!
!...LANL
!
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy
!
!matra = unkno(1, 2, ielem)
!matrb = unkno(1, 3, ielem)
!if(sqrt(matra**2+matrb**2).gt.1.d-8)then
!mapmt(1, 1, ielem) = matra/sqrt((matra)**2 + matrb**2)
!mapmt(1, 2, ielem) = matrb/sqrt((matra)**2 + matrb**2)
!mapmt(2, 1, ielem) =-matrb/sqrt((matra)**2 + matrb**2)
!mapmt(2, 2, ielem) = matra/sqrt((matra)**2 + matrb**2)

!else
!mapmt(1, 1, ielem) = 1.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 1.d0

!endif
!
!...Maire
!
!matra = dxdr**2 + dydr**2
!matrb = dxdr*dxds + dydr*dyds
!matrc = matrb
!matrd = dxds**2 + dyds**2
!
!...eigenvalues...
!
mapt  = matra + matrd
mapd  = matra*matrd - matrb*matrc
lamda1 = 0.5d0*mapt + sqrt(0.25d0*mapt**2-mapd)
lamda2 = 0.5d0*mapt - sqrt(0.25d0*mapt**2-mapd)
!!
if(abs(matrc).gt.1.d-6)then
!
if(matrc.gt.0.d0)then
!
! print*,'plus matrc',matrc
!
lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
mapmt(1, 2, ielem) =          matrc/lmat1
mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
mapmt(2, 2, ielem) =          matrc/lmat2
!
elseif(matrc.lt.0.d0)then
!
! print*,'minus matrc',matrc

lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
mapmt(2, 2, ielem) =          matrc/lmat1
mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
mapmt(1, 2, ielem) =          matrc/lmat2
endif
!
else
!
mapmt(1, 1, ielem) = 0.d0
mapmt(1, 2, ielem) = 0.d0
mapmt(2, 1, ielem) = 0.d0
mapmt(2, 2, ielem) = 0.d0
endif
!
enddo
!
!...Part 4: Impose limiter for mapped u and v
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
rc = cocent(1, ielem)
sc = cocent(2, ielem)
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xpq(1, iv)-rc)
bq(3, iv) = (xpq(2, iv)-sc)
enddo
!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
!
!...New mapped velocity components u and v
!
umap = uvtx*mapmt(1, 1, ielem) + vvtx*mapmt(1, 2, ielem)
vmap = uvtx*mapmt(2, 1, ielem) + vvtx*mapmt(2, 2, ielem)
!
unknvq(2, iv) = umap
unknvq(3, iv) = vmap
!
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
!
!...New mapped velocity components u and v
!
umap = uctr*mapmt(1, 1, ielem) + vctr*mapmt(1, 2, ielem)
vmap = uctr*mapmt(2, 1, ielem) + vctr*mapmt(2, 2, ielem)
!
unctr(2) = umap
unctr(3) = vmap
!
unmpx  = 1.d-10
unmpn  = 1.d10
!
ummax(2:3) = unctr(2:3)
ummin(2:3) = unctr(2:3)
!
do iv = 1, nvqua
!
!ummax(2:3) = unctr(2:3)
!ummin(2:3) = unctr(2:3)
!
do  istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
!
jelem=esuv1(istor)
!
!if(ielem.eq.3600) print*,'cell3600',iv,ipq(iv),istor,jelem
!
ucjel = unkno(1, 2, jelem)
vcjel = unkno(1, 3, jelem)
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
!
enddo
enddo
!
!do iq = 2, 3
!
!unmpx(iq) = max(ummax(iq),unmpx(iq))
!unmpn(iq) = min(ummin(iq),unmpn(iq))
!
!enddo
!
!enddo
!
do iv = 1, nvqua
!
do iq = 2, 3 !...only for vector components...
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
!call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!
call barthfctdiv(unkno(:,:,ielem),ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!
!call barthfct(unmpx(iq), unmpn(iq), unctr(iq), dunk(iq), afbar)
!
alfa(iq, iv) = afbar
!
enddo
!
if(indbd(ipq(iv)).eq.1)then

if(ncase.eq.1)then
if(coord(1, ipq(iv)).lt.1.d-6.or.abs(coord(1, ipq(iv))-1.d0).lt.1.d-6) then
!
alfa(2:3, iv) = 1.d0
!
endif
!
if(coord(2, ipq(iv)).lt.1.d-6.or.abs(coord(2, ipq(iv))-1.d0).lt.1.d-6) then
!
alfa(2:3, iv) = 1.d0
!
endif
!
endif
endif
!
enddo
!
do iq = 2, 3
aflim(iq, ielem) = minval(alfa(iq, 1:8)) !...Only consider the 4 vertices
enddo
!
enddo
!
!...Part 5: Impose limiter using the 1st invariant for mapped barth limiting coef...
!
do ie = 1, nquad
!
ielem = ie + ntria
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
dudxy = undu(ielem)
dumin = unmin_new(nq+2, ielem)
dumax = unmax_new(nq+2, ielem)
!
!if(ie.eq.6369) print*,'aflim',dumax/dudxy
!
if(dumin*dumax.gt.0.d0)then
!
if(dumax.lt.0.d0)then
afdu(ielem) = max(0.d0, min(1.d0,(dumax)/(dudxy+eps)))
elseif(dumin.gt.0.d0)then
afdu(ielem) = max(0.d0, min(1.d0,(dumin)/(dudxy+eps)))
endif
!
elseif(dumin*dumax.lt.0.d0)then
!
afdu(ielem)=0.d0
!
endif
!
geoel(10, ielem) = afdu(ielem)
!
enddo
!
!...Further technique...
!
do ie = 1, nquad
!
ielem = ie + ntria
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc !...Cell center density...
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
delu = undu(ielem)
!
volel = geoel(3, ielem)
macel = sqrt(volel)*abs(delu)/sdctr
afma(ielem) = min(1.d0,1.d0*macel)
!
afdu(ielem) =(afma(ielem))*afdu(ielem) +  (1.d0-afma(ielem))
!if(ielem.eq.1773) print*,'bad',afma(ielem),macel
enddo
!
!aflim = 1.d0
!
!...Part 4: Impose the 1st invariant limiting...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
!aflim(:,ielem) = aflim(:,ielem)*afdu(ielem)
!
enddo!
!
!...Part 5: Transfer back the limiter to the gobal frame...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!!
xicx = mapmt(1, 1, ielem); xicy = mapmt(1, 2, ielem)
etcx = mapmt(2, 1, ielem); etcy = mapmt(2, 2, ielem)
!
fixi = aflim(2, ielem) ; fiet = aflim(3, ielem)
!
fi11 = xicx**2*fixi + etcx**2*fiet
fi12 = xicx*xicy*fixi + etcx*etcy*fiet
fi21 = fi12
fi22 = xicy**2*fixi + etcy**2*fiet
!!
afvec(1, 1, ielem) = fi11
afvec(1, 2, ielem) = fi12
afvec(2, 1, ielem) = afvec(1, 2, ielem)
afvec(2, 2, ielem) = fi22
!!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_symprefvm_qc
!
!...Symmetry preserving techniques curved cell for onr 1/4 part...
!
subroutine barthlimit_symprefvm_qcquater(geoel, cocent,coord, ustar, unkno, ipqua, bface, intfac, afvec,esuv1, esuv2, &
unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:3,1:nsize),                 intent(in)::cocent
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),           intent(in):: ipqua
integer*4,dimension(1:nbfai,1:nbfac),           intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp,ivf
integer:: ielem, jelem, istor
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:nq+1) :: ummax, ummin
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
real*8,dimension(1:nsize)   ::afdu, afma
real*8,dimension(1:nsize) :: undu
real*8,dimension(1:nq+1, 1:nsize)   ::aflim
real*8::unmpx(1:nq), unmpn(1:nq)
!

real*8::eps,c00,c05,c10,c20
real*8:: uloca, vloca, ucjel, vcjel
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct, rhom
real*8:: afbar
real*8:: dudr, duds, dvdr, dvds
real*8:: etcx, etcy, xicx, xicy, fixi, fiet
real*8:: lamda1, lamda2
real*8:: mapd, mapt, matra, matrb, matrc, matrd, umap, vmap
real*8:: jf11, jf12, jf21, jf22
real*8:: fi11, fi12, fi21, fi22
real*8:: dxdr, dxds, dydr, dyds
real*8:: rm, rp, sm, sp, r, s
real*8:: drdx, drdy, dsdx, dsdy, jacom
real*8:: dudx, dudy, dvdx, dvdy, lmat1, lmat2
real*8:: dudxy, dumax, dumin
real*8:: rhomc,sdctr,macel,delu,volel
real*8:: signlim, vemag, afmag
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node

!...Specify the position of the symmetry boundary condition
if(nrz.eq.1)then
!afmag = 0.25d0*pi*(1.d0-0.5d0/(ncell/100.d0))
afmag = 39.d0*pi/160.d0
endif
!
!...Coloring the nodes at x axis
!
if(ncase.ne.1)then

do ifa =1 ,nbfac
!
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.221)then
indbd(intfac(3:(2+nvfac), ifa)) = 3
endif
endif
enddo
!
!...Coloring the nodes at y axis
!
do ifa =1 ,nbfac
!
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.222)then
do ivf = 1, nvfac
if(indbd(intfac(2+ivf, ifa)).eq.3)then
indbd(intfac(2+ivf, ifa)) = 10
else
indbd(intfac(2+ivf, ifa)) = 2
endif
enddo
endif
endif
enddo

else

do ifa =1 ,nbfac
!
indbd(intfac(3:(2+nvfac), ifa)) = 1
enddo
endif
!
eps = 1.e-6
mapmt = 0.d0
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
dudr = unkno(2, 2, ielem)
duds = unkno(3, 2, ielem)
dvdr = unkno(2, 3, ielem)
dvds = unkno(3, 3, ielem)
!
dudx = dudr
dudy = duds
dvdx = dvdr
dvdy = dvds
!
!...1st invariant
!
undu(ielem) = dudx +dvdy
!
!...LANL
!
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy
!
!matra = unkno(1, 2, ielem)
!matrb = unkno(1, 3, ielem)
if(sqrt(matra**2+matrb**2).gt.1.d-8)then
!mapmt(1, 1, ielem) = matra/sqrt((matra)**2 + matrb**2)
!mapmt(1, 2, ielem) = matrb/sqrt((matra)**2 + matrb**2)
!mapmt(2, 1, ielem) =-matrb/sqrt((matra)**2 + matrb**2)
!mapmt(2, 2, ielem) = matra/sqrt((matra)**2 + matrb**2)

else
!mapmt(1, 1, ielem) = 1.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 1.d0

endif
!
!...Maire
!
!matra = dxdr**2 + dydr**2
!matrb = dxdr*dxds + dydr*dyds
!matrc = matrb
!matrd = dxds**2 + dyds**2
!
!...eigenvalues...
!
mapt  = matra + matrd
mapd  = matra*matrd - matrb*matrc
lamda1 = 0.5d0*mapt + sqrt(0.25d0*mapt**2-mapd)
lamda2 = 0.5d0*mapt - sqrt(0.25d0*mapt**2-mapd)
!!
if(abs(matrc).gt.1.d-8)then
!
if(matrc.gt.0.d0)then
!
! print*,'plus matrc',matrc
!
lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
mapmt(1, 2, ielem) =          matrc/lmat1
mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
mapmt(2, 2, ielem) =          matrc/lmat2
!
elseif(matrc.lt.0.d0)then
!
! print*,'minus matrc',matrc

lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
mapmt(2, 2, ielem) =          matrc/lmat1
mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
mapmt(1, 2, ielem) =          matrc/lmat2
endif
!
else
!
mapmt(1, 1, ielem) = 0.d0
mapmt(1, 2, ielem) = 0.d0
mapmt(2, 1, ielem) = 0.d0
mapmt(2, 2, ielem) = 0.d0
endif
!
enddo
!
!...Part 4: Impose limiter for mapped u and v
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
rc = cocent(1, ielem)
sc = cocent(2, ielem)
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xpq(1, iv)-rc)
bq(3, iv) = (xpq(2, iv)-sc)
enddo
!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
!
!...New mapped velocity components u and v
!
umap = uvtx*mapmt(1, 1, ielem) + vvtx*mapmt(1, 2, ielem)
vmap = uvtx*mapmt(2, 1, ielem) + vvtx*mapmt(2, 2, ielem)
!
unknvq(2, iv) = umap
unknvq(3, iv) = vmap
!
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
!
!...New mapped velocity components u and v
!
umap = uctr*mapmt(1, 1, ielem) + vctr*mapmt(1, 2, ielem)
vmap = uctr*mapmt(2, 1, ielem) + vctr*mapmt(2, 2, ielem)
!
unctr(2) = umap
unctr(3) = vmap
!
ummax(2:3) = unctr(2:3)
ummin(2:3) = unctr(2:3)
!
do iv = 1, nvqua
!
!ummax(2:3) = unctr(2:3)
!ummin(2:3) = unctr(2:3)
!
do  istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
!
jelem=esuv1(istor)
!
!if(ielem.eq.3600) print*,'cell3600',iv,ipq(iv),istor,jelem
!
ucjel = unkno(1, 2, jelem)
vcjel = unkno(1, 3, jelem)
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
!
!...Symmetry BC...
!
if(ncase.ne.1)then
if(indbd(ipq(iv)).eq.3)then
!
!---1st
!
if(nrz.eq.0)then
ucjel = unkno(1, 2, jelem)
vcjel = -unkno(1, 3, jelem)
elseif(nrz.eq.1)then
signlim = sign(1.d0,unkno(1, 2, jelem))
vemag = signlim*sqrt(unkno(1, 2, jelem)**2 +unkno(1, 3, jelem)**2 )
ucjel = vemag*cos(afmag)
vcjel = vemag*sin(afmag)
endif
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)

elseif(indbd(ipq(iv)).eq.2)then
!
ucjel =-unkno(1, 2, jelem)
vcjel = unkno(1, 3, jelem)
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)

elseif(indbd(ipq(iv)).eq.10)then
!...x symmetry
ucjel = unkno(1, 2, jelem)
vcjel =-unkno(1, 3, jelem)
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
!...y symmetry
ucjel =-unkno(1, 2, jelem)
vcjel = unkno(1, 3, jelem)
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
!...0 symmetry
ucjel =-unkno(1, 2, jelem)
vcjel =-unkno(1, 3, jelem)
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)

endif
endif
!
enddo
!
enddo
!
do iv = 1, nvqua
do iq = 2, 3 !...only for vector components...
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!call barthfctdiv(unkno(:,:,ielem),ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!
alfa(iq, iv) = afbar
!
enddo
!
if(indbd(ipq(iv)).eq.1)then

if(ncase.eq.1)then
if(coord(1, ipq(iv)).lt.1.d-6.or.abs(coord(1, ipq(iv))-1.d0).lt.1.d-6) then
!
alfa(2:3, iv) = 1.d0
!
endif
!
if(coord(2, ipq(iv)).lt.1.d-6.or.abs(coord(2, ipq(iv))-1.d0).lt.1.d-6) then
!
alfa(2:3, iv) = 1.d0
!
endif
!
endif
endif
!
enddo
!
do iq = 2, 3
aflim(iq, ielem) = minval(alfa(iq, 1:8)) !...Only consider the 4 vertices
!aflim(iq, ielem) = minval(alfa(iq, 1:nvqua)) !...Only consider the 4 vertices
enddo
!
enddo
!
!...Part 5: Impose limiter using the 1st invariant for mapped barth limiting coef...
!
do ie = 1, nquad
!
ielem = ie + ntria
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
dudxy = undu(ielem)
dumin = unmin_new(nq+2, ielem)
dumax = unmax_new(nq+2, ielem)
!
!if(ie.eq.6369) print*,'aflim',dumax/dudxy
!
if(dumin*dumax.gt.0.d0)then
!
if(dumax.lt.0.d0)then
afdu(ielem) = max(0.d0, min(1.d0,(dumax)/(dudxy+eps)))
elseif(dumin.gt.0.d0)then
afdu(ielem) = max(0.d0, min(1.d0,(dumin)/(dudxy+eps)))
endif
!
elseif(dumin*dumax.lt.0.d0)then
!
afdu(ielem)=0.d0
!
endif
!
geoel(10, ielem) = afdu(ielem)
!
enddo
!
!...Further technique...
!
do ie = 1, nquad
!
ielem = ie + ntria
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc !...Cell center density...
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
delu = undu(ielem)
!
volel = geoel(3, ielem)
macel = sqrt(volel)*abs(delu)/sdctr
afma(ielem) = min(1.d0,1.d0*macel)
!
afdu(ielem) =(afma(ielem))*afdu(ielem) +  (1.d0-afma(ielem))
!if(ielem.eq.1773) print*,'bad',afma(ielem),macel
enddo
!
!aflim = 1.d0
!
!...Part 4: Impose the 1st invariant limiting...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
!aflim(:,ielem) = aflim(:,ielem)*afdu(ielem)
!
enddo!
!
!...Part 5: Transfer back the limiter to the gobal frame...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!!
xicx = mapmt(1, 1, ielem); xicy = mapmt(1, 2, ielem)
etcx = mapmt(2, 1, ielem); etcy = mapmt(2, 2, ielem)
!
fixi = aflim(2, ielem) ; fiet = aflim(3, ielem)
!
fi11 = xicx**2*fixi + etcx**2*fiet
fi12 = xicx*xicy*fixi + etcx*etcy*fiet
fi21 = fi12
fi22 = xicy**2*fixi + etcy**2*fiet
!!
afvec(1, 1, ielem) = fi11
afvec(1, 2, ielem) = fi12
afvec(2, 1, ielem) = afvec(1, 2, ielem)
afvec(2, 2, ielem) = fi22
!!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_symprefvm_qcquater
!
!...Symmetry preserving techniques curved cell FVM...
!
subroutine barthlimit_symprefvm_curvquad(geoel, cocent,coord, ustar, unkno, ipqua, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:3,1:nsize),                 intent(in)::cocent
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),           intent(in):: ipqua
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp
integer:: ielem, jelem, istor
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:nq+1) :: ummax, ummin
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
real*8,dimension(1:nsize)   ::afdu, afma
real*8,dimension(1:nsize) :: undu
real*8,dimension(1:nq+1, 1:nsize)   ::aflim
real*8::unmpx(1:nq), unmpn(1:nq)
!

real*8::eps,c00,c05,c10,c20
real*8:: uloca, vloca, ucjel, vcjel
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct, rhom
real*8:: afbar
real*8:: dudr, duds, dvdr, dvds
real*8:: etcx, etcy, xicx, xicy, fixi, fiet
real*8:: lamda1, lamda2
real*8:: mapd, mapt, matra, matrb, matrc, matrd, umap, vmap
real*8:: jf11, jf12, jf21, jf22
real*8:: fi11, fi12, fi21, fi22
real*8:: dxdr, dxds, dydr, dyds
real*8:: rm, rp, sm, sp, r, s
real*8:: drdx, drdy, dsdx, dsdy, jacom
real*8:: dudx, dudy, dvdx, dvdy, lmat1, lmat2
real*8:: dudxy, dumax, dumin
real*8:: rhomc,sdctr,macel,delu,volel
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:(2+nvfac), ifa)) = 1
enddo
!
eps = 1.e-6
mapmt = 0.d0
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
matra = unkno(1, 2, ielem)
matrb = unkno(1, 3, ielem)
if(sqrt(matra**2+matrb**2).gt.1.d-8)then
mapmt(1, 1, ielem) = matra/sqrt((matra)**2 + matrb**2)
mapmt(1, 2, ielem) = matrb/sqrt((matra)**2 + matrb**2)
mapmt(2, 1, ielem) =-matrb/sqrt((matra)**2 + matrb**2)
mapmt(2, 2, ielem) = matra/sqrt((matra)**2 + matrb**2)

else
mapmt(1, 1, ielem) = 1.d0
mapmt(1, 2, ielem) = 0.d0
mapmt(2, 1, ielem) = 0.d0
mapmt(2, 2, ielem) = 1.d0

endif
!
enddo
!
!...Part 4: Impose limiter for mapped u and v
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= cocent(1, ielem) !...mass center...
sc= cocent(2, ielem)
!
xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)
bq(3, iv) = (yvq(iv)-sc)
enddo
!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
!
!...New mapped velocity components u and v
!
umap = uvtx*mapmt(1, 1, ielem) + vvtx*mapmt(1, 2, ielem)
vmap = uvtx*mapmt(2, 1, ielem) + vvtx*mapmt(2, 2, ielem)
!
unknvq(2, iv) = umap
unknvq(3, iv) = vmap
!
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
!
!...New mapped velocity components u and v
!
umap = uctr*mapmt(1, 1, ielem) + vctr*mapmt(1, 2, ielem)
vmap = uctr*mapmt(2, 1, ielem) + vctr*mapmt(2, 2, ielem)
!
unctr(2) = umap
unctr(3) = vmap
!
unmpx  = 1.d-10
unmpn  = 1.d10
!
do iv = 1, nvqua
!
ummax(2:3) = unctr(2:3)
ummin(2:3) = unctr(2:3)
!
do  istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
!
jelem=esuv1(istor)
!
!if(ielem.eq.3600) print*,'cell3600',iv,ipq(iv),istor,jelem
!
ucjel = unkno(1, 2, jelem)
vcjel = unkno(1, 3, jelem)
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
!
enddo
!
!do iq = 2, 3
!
!unmpx(iq) = max(ummax(iq),unmpx(iq))
!unmpn(iq) = min(ummin(iq),unmpn(iq))
!
!enddo
!
!enddo
!
!do iv = 1, nvqua
!
do iq = 2, 3 !...only for vector components...
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!
!call barthfct(unmpx(iq), unmpn(iq), unctr(iq), dunk(iq), afbar)
!
alfa(iq, iv) = afbar
!
enddo
!
if(indbd(ipq(iv)).eq.1)then

if(ncase.eq.1)then
if(coord(1, ipq(iv)).lt.1.d-6.or.abs(coord(1, ipq(iv))-1.d0).lt.1.d-6) then
!
alfa(2:3, iv) = 1.d0
!
endif
!
if(coord(2, ipq(iv)).lt.1.d-6.or.abs(coord(2, ipq(iv))-1.d0).lt.1.d-6) then
!
alfa(2:3, iv) = 1.d0
!
endif
!
endif
endif
!
enddo
!
do iq = 2, 3
aflim(iq, ielem) = minval(alfa(iq, 1:4)) !...Only consider the 4 vertices
enddo
!
enddo
!
!...Part 5: Transfer back the limiter to the gobal frame...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xicx = mapmt(1, 1, ielem); xicy = mapmt(1, 2, ielem)
etcx = mapmt(2, 1, ielem); etcy = mapmt(2, 2, ielem)
!
fixi = aflim(2, ielem) ; fiet = aflim(3, ielem)
!
fi11 = xicx**2*fixi + etcx**2*fiet
fi12 = xicx*xicy*fixi + etcx*etcy*fiet
fi21 = fi12
fi22 = xicy**2*fixi + etcy**2*fiet
!!
afvec(1, 1, ielem) = fi11
afvec(1, 2, ielem) = fi12
afvec(2, 1, ielem) = afvec(1, 2, ielem)
afvec(2, 2, ielem) = fi22
!!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_symprefvm_curvquad
!
!...Face integral (mass center) for hybrid curv quad...
!
subroutine rhsifacedg_lagfvm_curvquad(ipqua, unkno, ustar,fstarq, gelagq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvqua, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(inout)::rhsel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:2, 1:nvqua) :: ipfq
integer,dimension(1:4, 1:nvqua) :: ndshp
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
!real*8,dimension(1:ndimn, 1:ndegr, 1:4, 1:nvqua)::lpnpq
real*8, dimension(1:2, 1:ndegr, 1:4):: lpnp
real*8, dimension(1:ndegr, 1:4)::clpnp
real*8::xvq(nvqua), yvq(nvqua),bq(1:ndegr,1:nvqua)
real*8::vnorm(1:3, 1:4, 1:nvqua)

!...local real number
real*8::eps,c00,c05,c10,c20,c13,c23,c16
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
data c23   / 0.6666666666666666d0 /
!data c13   / 0.5d0 /
!data c23   / 0.5d0 /
!
!...Zero out plnpn, ulnpn
!
ndshp(1, 1) = 4; ndshp(2, 1) = 8; ndshp(3, 1) = 2; ndshp(4, 1) = 5
ndshp(1, 2) = 1; ndshp(2, 2) = 5; ndshp(3, 2) = 3; ndshp(4, 2) = 6;
ndshp(1, 3) = 2; ndshp(2, 3) = 6; ndshp(3, 3) = 4; ndshp(4, 3) = 7;
ndshp(1, 4) = 3; ndshp(2, 4) = 7; ndshp(3, 4) = 1; ndshp(4, 4) = 8;

ndshp(1, 5) = 1; ndshp(2, 5) = 2
ndshp(1, 6) = 2; ndshp(2, 6) = 3
ndshp(1, 7) = 3; ndshp(2, 7) = 4
ndshp(1, 8) = 4; ndshp(2, 8) = 1
!
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) =  gelagq(1:3,  4, ie);   vnorm(1:3, 2, 1) = gelagq(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 3, 1) =  gelagq(1:3, 12, ie);   vnorm(1:3, 4, 1) = gelagq(1:3, 5, ie) !...For point ip(1)

vnorm(1:3, 1, 2) =  gelagq(1:3, 1, ie);    vnorm(1:3, 2, 2) = gelagq(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) =  gelagq(1:3, 6, ie);    vnorm(1:3, 4, 2) = gelagq(1:3, 7, ie) !...For point ip(2)

vnorm(1:3, 1, 3) =  gelagq(1:3, 2, ie);    vnorm(1:3, 2, 3) = gelagq(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) =  gelagq(1:3, 8, ie);    vnorm(1:3, 4, 3) = gelagq(1:3, 9, ie) !...For point ip(3)

vnorm(1:3, 1, 4) =  gelagq(1:3, 3, ie);    vnorm(1:3, 2, 4) = gelagq(1:3, 4, ie) !...For point ip(4)
vnorm(1:3, 3, 4) =  gelagq(1:3, 10, ie);   vnorm(1:3, 4, 4) = gelagq(1:3, 11, ie) !...For point ip(4)

vnorm(1:3, 1, 5) =  gelagq(1:3, 5, ie);    vnorm(1:3, 2, 5) = gelagq(1:3, 6, ie) !...For point ip(5)
vnorm(1:3, 1, 6) =  gelagq(1:3, 7, ie);    vnorm(1:3, 2, 6) = gelagq(1:3, 8, ie) !...For point ip(6)
vnorm(1:3, 1, 7) =  gelagq(1:3, 9, ie);    vnorm(1:3, 2, 7) = gelagq(1:3, 10, ie) !...For point ip(7)
vnorm(1:3, 1, 8) =  gelagq(1:3, 11, ie);   vnorm(1:3, 2, 8) = gelagq(1:3, 12, ie) !...For point ip(8)
!
!...The vertex constituting one cell...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, 4
!
do ideg = 1, 1
lpnp(1:2, ideg, 3) = vnorm(3, 3, iv)*vnorm(1:2, 3, iv)
lpnp(1:2, ideg, 1) = 0.d0
!
lpnp(1:2, ideg, 4) = vnorm(3, 4, iv)*vnorm(1:2, 4, iv)
lpnp(1:2, ideg, 2) = 0.d0
!
clpnp(ideg, 3) = 1.d0
clpnp(ideg, 1) = 0.d0
clpnp(ideg, 4) = 1.d0
clpnp(ideg, 2) = 0.d0
enddo
!
do ifa =1, 4
!
ulnpn(1)  = ulnpn(1) + c13*ustar(1, ipq(iv))*lpnp(1, 1, ifa) +&
c13*ustar(2, ipq(iv))*lpnp(2, 1, ifa)
plnpn(1, 1)= plnpn(1, 1) + c13*fstarq(1, ifa, iv, ie)*clpnp(1, ifa)
plnpn(2, 1)= plnpn(2, 1) + c13*fstarq(2, ifa, iv, ie)*clpnp(1, ifa)
!
elnpn(1)   = elnpn(1)+&
c13*ustar(1, ipq(iv))*fstarq(1, ifa, iv, ie)*clpnp(1, ifa) +&
c13*ustar(2, ipq(iv))*fstarq(2, ifa, iv, ie)*clpnp(1, ifa)
!
!if(ie==1831) print*,iv,ipq(iv),ielem, ulnpn(1:ndegr),ustar(1:2,ipq(iv)),&
!             lpnp(1:2,1,iv)
!
enddo
!
!if(ie==1) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
!
do iv = 5, 8
!
do ideg = 1, 1
lpnp(1:2, ideg, 1) = vnorm(3, 1, iv)*vnorm(1:2, 1, iv)
lpnp(1:2, ideg, 2) = vnorm(3, 2, iv)*vnorm(1:2, 2, iv)
!
clpnp(ideg, 1) = 1.d0
clpnp(ideg, 2) = 1.d0
enddo
!
do ifa =1, 2
!
ulnpn(1)  = ulnpn(1) + c23*ustar(1, ipq(iv))*lpnp(1, 1, ifa) +&
c23*ustar(2, ipq(iv))*lpnp(2, 1, ifa)
plnpn(1, 1)= plnpn(1, 1) + c23*fstarq(1, ifa, iv, ie)*clpnp(1, ifa)
plnpn(2, 1)= plnpn(2, 1) + c23*fstarq(2, ifa, iv, ie)*clpnp(1, ifa)
!
elnpn(1)   = elnpn(1)+&
c23*ustar(1, ipq(iv))*fstarq(1, ifa, iv, ie)*clpnp(1, ifa) +&
c23*ustar(2, ipq(iv))*fstarq(2, ifa, iv, ie)*clpnp(1, ifa)
!
!
!if(ie==1831) print*,iv,ifa,ipq(iv),ielem, ustar(1:2,ipq(iv)),&
!bq(1, 8),  bq(1, ndshp(1, 8)),lpnp(1:2, 1:3, 1)
!
enddo
!
!if(ie==1) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
!
rhsel(1, 1, ielem) =  ulnpn(1)
rhsel(1, 2, ielem) =  plnpn(1, 1)
rhsel(1, 3, ielem) =  plnpn(2, 1)
rhsel(1, 4, ielem) =  elnpn(1)
!
!if(ie==1831)  print*,'rhs iface',ielem, ie,rhsel(1:3, 1:4, ielem) !,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lagfvm_curvquad
!
!...Source term integration for curv quad grids...
!
subroutine rhsdomnsrcdg_lagfvm_curvquad(intfac, ipqua, coord, geoel,rhsel)
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
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
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
djaco = wi*(dxdr*dyds - dydr*dxds)
!
!...Gauss points...
!
xg =0.d0; yg= 0.d0
!
do ishp = 1, npqua
xg = xg + shpq(ishp)*xpq(1,ishp)
yg = yg + shpq(ishp)*xpq(2,ishp)
enddo
!
!
!...Basis function for solutions...
!
b(1) = 1.d0
!
!...source term
!
src = 0.25d0*pi*(cos(3.d0*pi*xg)*cos(pi*yg) - cos(3.d0*pi*yg)*cos(pi*xg))/(gamlg-1.d0)
!src = 3.d0/8.d0*pi*(cos(pi*xg)*cos(3.d0*pi*yg) - cos(3.d0*pi*xg)*cos(pi*yg))!/(gamlg-1.d0)
!
!    src = 0.5d0*pi/(gamlg-1.d0)*(sin(2.d0*pi*yg)*cos(pi*xg)*sin(pi*yg) - sin(2.d0*pi*xg)*sin(pi*xg)*cos(pi*yg))
!
!finally, scatter the contribution to the RHS
!
! if(ie==2) print*,'src rhs', 3.5d0/9.d0,rhsel(1, 4, 2),coorp(1, 1:3), coorp(2, 1:3)
!
rhsel(1,4,ielem)=rhsel(1,4,ielem) + src*b(1)*djaco
!
!if(ie==2) print*,'src rhs', rhsel(1, 4, 2),src,b(1),djaco,src*b(1)*djaco,xg,yg
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
end subroutine rhsdomnsrcdg_lagfvm_curvquad
!
!...Error for TGV with FVM
!
subroutine geterror_lagtgvfvm(intfac, iptri, ipqua, coord, coold, geoel, unkno)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem
!
!...local integer array
!
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvtri) :: xp
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndegr):: b,bv
real*8:: unknv(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8,dimension(1:3,1:nsize)::cocent
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s,  dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc,rm, sm,rp,sp
real*8::xg, yg, xc, yc
real*8::rho,uadv,vadv,eadv,rhom
real*8::rhoc, uctr, vctr, ectr, pctr
real*8::pres, pexa, eiadv, eiexa, rhoexa
real*8::uexa,vexa,veexa,vemag
real*8::djaco, wi, errl2
real*8::radie, radii,radie2,radii2, radig2,paras,rhoi,rhoe,htkid
real*8::rhon, rhoini,rhoct,pvtx
real*8::rcv, scv
real*8::errlp(nq)
!
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8::weigh(ngausd), posi(2,ngausd)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
call rutope(2, ngausd, posi, weigh)
call ruqope(2, ngausdq, posiq, weighq)
!
errlp = 0.d0
!
open(8,file='pres-tgv.dat')
!
!...For quads...
!
call getgeoel_lag_fvm(iptri, ipqua, cocent, coord)
!
!...Loop over elements
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Geometry parameters for reference cell...
!
rc = cocent(1, ielem)
sc = cocent(2, ielem)
!
!...Gauss loop
!
do ig = 1,ngausdq !...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
wi  = weighq(ig)
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
do ishp = 1, nvqua
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
xg =0.d0
yg = 0.d0
!
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpq(1, ishp)
yg = yg + shpq(ishp)*xpq(2, ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
b(1) = 1.d0
b(2) = (xg-rc) 
b(3) = (yg-sc) 
!
!...Jacobian transformation matrix
!
rhoct  = 1.d0/unkno(1,1,ielem) 
uctr = unkno(1,2,ielem) 
vctr = unkno(1,3,ielem) 
ectr = unkno(1,4,ielem) 
!
pres = (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
pvtx = pres + unkno(2,4,ielem)*b(2)  + unkno(3,4,ielem)*b(3)
!
pexa = 0.25d0*(cos(2.d0*pi*rc) + cos(2.d0*pi*sc)) + 10.d0
uexa   = sin(pi*rc)*cos(pi*sc)
vexa   =-cos(pi*rc)*sin(pi*sc)
rhoexa = 1.d0
!
vemag = sqrt(uctr**2 + vctr**2)
veexa = sqrt(uexa**2 + vexa**2)
!
errlp(1) = errlp(1) + abs(rhoct-rhoexa)*djaco
errlp(2) = errlp(2) + abs(vemag-veexa)*djaco
errlp(3) = errlp(3) + abs(pres - pexa)*djaco

!
enddo
!
write(8,'(i8, 17e32.16)')ielem,xpq(1,1:8),xpq(2,1:8),pres
!
650 enddo
!
close(8)
!
if(ncase.eq.1)  print*,'L2 error for Taylor-Green-Vortex',log10(errlp)
!
end subroutine geterror_lagtgvfvm
!
!...Calculate the velocity at the middle point for curved mesh with FVM...
!
subroutine getvelo_mpt_curvfvm(ustar,geoel,cocent,gelag,gelagq,intfac,iptri,ipqua,coord,unkno,indnd,&
 aflim, afvec, vlave, vnulq)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:nsize),               intent(in)::cocent
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad),      intent(in)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:2, 1:npoin),              intent(inout)::ustar
integer*4,dimension(1:npoin),                intent(in)::indnd
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),         intent(in)::afvec
real*8,dimension(1:ndimn,1:npoin),           intent(in)::vlave
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(in):: vnulq
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ielem
integer::iv
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!
integer,dimension(npoin)::icemp
!
real*8::eps
real*8::lnvp(2,nvqua),gpn(2, nvqua)
real*8::unknvt(1:nq, 1:nvtri)
real*8::unknvq(1:nq, 1:nvqua)
real*8::vnorm(3, nvqua)
real*8::xv(nvtri), yv(nvtri), bt(ndegr, nvtri) 
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, nvqua) 
real*8::xpf(1:2, 1:nvfac)
real*8::unmpt(4,2,npoin)
real*8::ft(2, 2), fn(2, 2)
real*8::r,s,rhoi,rhon
!real*8,allocatable:: ucurv(:, :)
real*8::dudr,duds,dvdr,dvds,pvtx
real*8::uvtxr,vvtxr,evtxr, pvtxr
real*8::uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu,divu,murie
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx,dux,duy
real*8::fnx,fny, rho,rhomc,rhovt,rhomv
real*8::dr, ds,rc,sc,othog,rcv,scv
real*8::ftx, fty
!
!...For quadratic mesh, only nafac high-order nodes need be recalculated at most...
!
eps = 1.d-6
!
icemp = 0
unmpt = 0.d0
!
!...Part I: Get the pres and murie at Traingle cell
!
do 250 ie = 1,ntria !...(1)ie = 1,nelem
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
!...shape functions
!
rc= cocent(1, ielem) !...mass center...
sc= cocent(2, ielem)
!
xv(1:nvtri) = coord(1, ipt(1:nvtri))
yv(1:nvtri) = coord(2, ipt(1:nvtri))
!
do iv =1 ,nvtri
!...Basis function
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc) 
bt(3, iv) = (yv(iv)-sc) 
enddo
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...
!...zero out unknv
unknvt = 0.d0
!
do iv   = 4, nvtri
!
do ideg = 1,mdegr
unknvt(1:nq-1, iv) = unknvt(1:nq-1, iv) + unkno(ideg,1:nq-1,ielem)*bt(ideg, iv)
enddo
!
unknvt(4, iv) = pctr  + unkno(2,4,ielem)*bt(2, iv) + unkno(3,4,ielem)*bt(3, iv)
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvt(1, iv)
endif
!
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
pvtx = unknvt(4, iv) 
!
!...Limiter
!
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvt(1, iv) - rhomc)
rhovt = 1.d0/rhomv
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bt(2, iv) + duds*bt(3, iv)
vvtx = unkno(1,3,ielem)  + dvdr*bt(2, iv) + dvds*bt(3, iv)
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!!
!...updtae unknv(2:3,:)
unknvt(2, iv) = uvtx
unknvt(3 ,iv) = vvtx
endif
!
!...Get impedence coefficient...
!
dux= vlave(1, ipt(iv))-unknvt(2, iv)
duy= vlave(2, ipt(iv))-unknvt(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
icemp(ipt(iv)) = icemp(ipt(iv)) + 1
unmpt(1, icemp(ipt(iv)), ipt(iv)) = pvtx
unmpt(2, icemp(ipt(iv)), ipt(iv)) = murie
unmpt(3, icemp(ipt(iv)), ipt(iv)) = unknvt(2, iv)
unmpt(4, icemp(ipt(iv)), ipt(iv)) = unknvt(3, iv)
!
enddo
250 enddo  !...(1)ie = 1,ntria!

!
!...Part II: Get the pres and murie at right cell
!
do 350 ie = 1,nquad !...(1)ie = 1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
!...shape functions
!
rc= cocent(1, ielem) !...mass center...
sc= cocent(2, ielem)
!
xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Get LN for one vertex for judging expanding or compressing
!
lnvp = 0.d0; gpn = 0.d0;
!
lnvp(1:2, 5) = gelagq(1:2,  9,ie)*gelagq(3, 9, ie)
lnvp(1:2, 6) = gelagq(1:2, 10,ie)*gelagq(3, 10,ie)
lnvp(1:2, 7) = gelagq(1:2, 11,ie)*gelagq(3, 11,ie)
lnvp(1:2, 8) = gelagq(1:2, 12,ie)*gelagq(3, 12,ie)
!
gpn(1:2, 5) = gelagq(1:2,  9,ie)
gpn(1:2, 6) = gelagq(1:2, 10,ie)
gpn(1:2, 7) = gelagq(1:2, 11,ie)
gpn(1:2, 8) = gelagq(1:2, 12,ie)
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc) 
bq(3, iv) = (yvq(iv)-sc) 
enddo
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc !...Cell center density...
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...
!...zero out unknv
unknvq = 0.d0
!
do iv   = 5, 8
!
do ideg = 1,mdegr
unknvq(1:nq-1, iv) = unknvq(1:nq-1, iv) + unkno(ideg,1:nq-1,ielem)*bq(ideg, iv)
enddo
!
unknvq(4, iv) = pctr  + unkno(2,4,ielem)*bq(2, iv) + unkno(3,4,ielem)*bq(3, iv)
!
!...Remove null mode...
!
unknvq(2:3, iv)=unknvq(2:3, iv)-vnulq(1:2,iv,ie)
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvq(1, iv)
endif

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
pvtx = unknvq(4, iv)
!
!...Limiter
!
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvq(1, iv) - rhomc)
rhovt = 1.d0/rhomv
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
vvtx = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!!
!...updtae unknv(2:3,:)
unknvq(2, iv) = uvtx-vnulq(1,iv,ie)
unknvq(3 ,iv) = vvtx-vnulq(2,iv,ie)
!
!if(ielem.eq.1119.or.ielem.eq.1120)print*,'aflim',aflim(2:3,ielem),iv,uvtx,vvtx,uvtx**2+vvtx**2
!
endif
!
!...Get impedence coefficient...
!
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
!
!dux= vlave(1, ipq(iv))-unkno(1, 2, ielem)
!duy= vlave(2, ipq(iv))-unkno(1, 3, ielem)
divu = dux*lnvp(1, iv) + duy*lnvp(2, iv)
!
if(divu.le.0.d0)then
deltu = abs(dux*gpn(1, iv) + duy*gpn(2, iv))
else
deltu = 0.d0
endif
deltu = sqrt(dux**2 + duy**2)
murie = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
icemp(ipq(iv)) = icemp(ipq(iv)) + 1
unmpt(1, icemp(ipq(iv)), ipq(iv)) = pvtx
unmpt(2, icemp(ipq(iv)), ipq(iv)) = murie
unmpt(3, icemp(ipq(iv)), ipq(iv)) = unknvq(2, iv)
unmpt(4, icemp(ipq(iv)), ipq(iv)) = unknvq(3, iv)
!
enddo
!
350 enddo  !...(1)ie = 1,nelem!
!
!...Part III: Get the velocity at middle point
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
!
ft(1, 1) = xpf(1 ,2)- xpf(1, 1)
ft(2, 1) = xpf(2, 2)- xpf(2, 1)
!
lenmc = sqrt(ft(1, 1)**2 + ft(2, 1)**2)
!
ft(1, 1) = ft(1, 1)/lenmc
ft(2, 1) = ft(2, 1)/lenmc
!
fn(1, 1) =  ft(2, 1)
fn(2, 1) = -ft(1, 1)
!
!...For linear PM
ft(1, 2) = xpf(1 ,3)- xpf(1, 1)
ft(2, 2) = xpf(2, 3)- xpf(2, 1)
!
lenmc = sqrt(ft(1, 2)**2 + ft(2, 2)**2)
!
ft(1, 2) = ft(1, 2)/lenmc
ft(2, 2) = ft(2, 2)/lenmc
!
othog = abs(ft(1, 2)*fn(1,1)+ ft(2, 2)*fn(2,1))
!
!if(ipf(3).eq.14) print*,'pt 14', lenmc, ifa
!
if(othog.lt.1.d-4)then
!
!...(ifa.le.nbfac) For boundary cells...
!
if(ifa.le.nbfac)then
!
presl = unmpt(1, 1, ipf(3))
mufal = unmpt(2, 1, ipf(3))
uvtxl = unmpt(3, 1, ipf(3))
vvtxl = unmpt(4, 1, ipf(3))
!
fnx = fn(1, 1) !...face normal vector
fny = fn(2, 1)
!
!...Ma
!
if(indnd(ipf(3)).eq.0)then
ustar(1, ipf(3)) = uvtxl + presl/mufal*fnx
ustar(2, ipf(3)) = vvtxl + presl/mufal*fny
endif
!
!...Morgan
!ustar(1, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*ftx/(mufal) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!ustar(2, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*fty/(mufal)
!
!
!...(ifa.gt.nbfac) For interior cells...
!
elseif(ifa.gt.nbfac)then
!
presl = unmpt(1, 1, ipf(3))
mufal = unmpt(2, 1, ipf(3))
uvtxl = unmpt(3, 1, ipf(3))
vvtxl = unmpt(4, 1, ipf(3))
!
presr = unmpt(1, 2, ipf(3))
mufar = unmpt(2, 2, ipf(3))
uvtxr = unmpt(3, 2, ipf(3))
vvtxr = unmpt(4, 2, ipf(3))
!
fnx = fn(1, 1) !...face normal vector
fny = fn(2, 1)
!
ftx = -fny
fty =  fnx
!...
!...Mar
!
if(indnd(ipf(3)).eq.0)then
ustar(1, ipf(3)) = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fnx
ustar(2, ipf(3)) = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fny
!
!ustar(1, ipf(3)) = (mufal*(uvtxl*fnx+vvtxl*fny) + mufar*(uvtxr*fnx+vvtxr*fny))*fnx/(mufal+mufar) &
!                    - (presr- presl)/(mufal+mufar)*fnx + &
!                  0.5d0*((uvtxl*ftx+vvtxl*fty) + (uvtxr*ftx+vvtxr*fty))*ftx
!ustar(2, ipf(3)) = (mufal*(uvtxl*fnx+vvtxl*fny) + mufar*(uvtxr*fnx+vvtxr*fny))*fny/(mufal+mufar) &
!                    - (presr- presl)/(mufal+mufar)*fny + &
!                  0.5d0*((uvtxl*ftx+vvtxl*fty) + (uvtxr*ftx+vvtxr*fty))*fty
endif
!
!ftx = -fny
!fty =  fnx
!
!...Morgan

!ustar(1, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*ftx/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!    ustar(2, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*fty/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(2, 1, iv)
!
!    if(ipf(3).eq.443) print*,'midlle velocity',ustar(1:2, ipf(3)),ipf(3),mufal+mufar,fnx,fny
endif
!
endif
!
450 enddo

!
end subroutine getvelo_mpt_curvfvm
!
!...Calculate the velocity at the middle point for curved mesh with FVM...
!
subroutine getvelo_mpt_curvfvm2(ustar,geoel,cocent,gelag,gelagq,intfac,iptri,ipqua,coord,unkno,indnd,&
aflim, afvec, vlave, vnulq,munacn)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:nsize),               intent(in)::cocent
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad),      intent(in)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:2, 1:npoin),              intent(inout)::ustar
integer*4,dimension(1:npoin),                intent(in)::indnd
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),         intent(in)::afvec
real*8,dimension(1:ndimn,1:npoin),           intent(in)::vlave
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(in):: vnulq
real*8,intent(in)::munacn(1:2, 1:2, 1:npoin)

!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ielem,ishp
integer::iv
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!
integer,dimension(npoin)::icemp
!
real*8::eps
real*8::dshpr(nvfac),dxdr,dydr,dxds,dyds,djaco,dwav1,dwav2
real*8::lnvp(2,nvqua),gpn(2, nvqua)
real*8::unknvt(1:nq, 1:nvtri)
real*8::unknvq(1:nq, 1:nvqua)
real*8::vnorm(3, nvqua)
real*8::xv(nvtri), yv(nvtri), bt(ndegr, nvtri)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, nvqua)
real*8::xpf(1:2, 1:nvfac)
real*8::unmpt(4,2,npoin)
real*8::ft(2, 2), fn(2, 2)
real*8::r,s,rhoi,rhon
!real*8,allocatable:: ucurv(:, :)
real*8::dudr,duds,dvdr,dvds,pvtx
real*8::uvtxr,vvtxr,evtxr, pvtxr
real*8::uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu,divu,murie
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx,dux,duy
real*8::fnx,fny, rho,rhomc,rhovt,rhomv
real*8::dr, ds,rc,sc,othog,rcv,scv
real*8::ftx, fty
!
!...For quadratic mesh, only nafac high-order nodes need be recalculated at most...
!
eps = 1.d-6
!
icemp = 0
unmpt = 0.d0
!
!...Part I: Get the pres and murie at Traingle cell
!
do 250 ie = 1,ntria !...(1)ie = 1,nelem
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
!...shape functions
!
rc= cocent(1, ielem) !...mass center...
sc= cocent(2, ielem)
!
xv(1:nvtri) = coord(1, ipt(1:nvtri))
yv(1:nvtri) = coord(2, ipt(1:nvtri))
!
do iv =1 ,nvtri
!...Basis function
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)
bt(3, iv) = (yv(iv)-sc)
enddo
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...
!...zero out unknv
unknvt = 0.d0
!
do iv   = 4, nvtri
!
do ideg = 1,mdegr
unknvt(1:nq-1, iv) = unknvt(1:nq-1, iv) + unkno(ideg,1:nq-1,ielem)*bt(ideg, iv)
enddo
!
unknvt(4, iv) = pctr  + unkno(2,4,ielem)*bt(2, iv) + unkno(3,4,ielem)*bt(3, iv)
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvt(1, iv)
endif
!
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
pvtx = unknvt(4, iv)
!
!...Limiter
!
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvt(1, iv) - rhomc)
rhovt = 1.d0/rhomv
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bt(2, iv) + duds*bt(3, iv)
vvtx = unkno(1,3,ielem)  + dvdr*bt(2, iv) + dvds*bt(3, iv)
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!!
!...updtae unknv(2:3,:)
unknvt(2, iv) = uvtx
unknvt(3 ,iv) = vvtx
endif
!
!...Get impedence coefficient...
!
dux= vlave(1, ipt(iv))-unknvt(2, iv)
duy= vlave(2, ipt(iv))-unknvt(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
icemp(ipt(iv)) = icemp(ipt(iv)) + 1
unmpt(1, icemp(ipt(iv)), ipt(iv)) = pvtx
unmpt(2, icemp(ipt(iv)), ipt(iv)) = murie
unmpt(3, icemp(ipt(iv)), ipt(iv)) = unknvt(2, iv)
unmpt(4, icemp(ipt(iv)), ipt(iv)) = unknvt(3, iv)
!
enddo
250 enddo  !...(1)ie = 1,ntria!

!
!...Part II: Get the pres and murie at right cell
!
do 350 ie = 1,nquad !...(1)ie = 1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
!...shape functions
!
rc= cocent(1, ielem) !...mass center...
sc= cocent(2, ielem)
!
xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Get LN for one vertex for judging expanding or compressing
!
lnvp = 0.d0; gpn = 0.d0;
!
lnvp(1:2, 5) = gelagq(1:2,  9,ie)*gelagq(3, 9, ie)
lnvp(1:2, 6) = gelagq(1:2, 10,ie)*gelagq(3, 10,ie)
lnvp(1:2, 7) = gelagq(1:2, 11,ie)*gelagq(3, 11,ie)
lnvp(1:2, 8) = gelagq(1:2, 12,ie)*gelagq(3, 12,ie)
!
gpn(1:2, 5) = gelagq(1:2,  9,ie)
gpn(1:2, 6) = gelagq(1:2, 10,ie)
gpn(1:2, 7) = gelagq(1:2, 11,ie)
gpn(1:2, 8) = gelagq(1:2, 12,ie)
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)
bq(3, iv) = (yvq(iv)-sc)
enddo
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc !...Cell center density...
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...
!...zero out unknv
unknvq = 0.d0
!
do iv   = 5, 8
!
do ideg = 1,mdegr
unknvq(1:nq-1, iv) = unknvq(1:nq-1, iv) + unkno(ideg,1:nq-1,ielem)*bq(ideg, iv)
enddo
!
unknvq(4, iv) = pctr  + unkno(2,4,ielem)*bq(2, iv) + unkno(3,4,ielem)*bq(3, iv)
!
!...Remove null mode...
!
unknvq(2:3, iv)=unknvq(2:3, iv)-vnulq(1:2,iv,ie)
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvq(1, iv)
endif

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
pvtx = unknvq(4, iv)
!
!...Limiter
!
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvq(1, iv) - rhomc)
rhovt = 1.d0/rhomv
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
uvtx = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
vvtx = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
!
!if(ie.ge.2625.and.ie.le.2628) print*,'ie26252628',ielem,pctr,aflim(4, ielem),pvtx
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!!
!...updtae unknv(2:3,:)
unknvq(2, iv) = uvtx-vnulq(1,iv,ie)
unknvq(3 ,iv) = vvtx-vnulq(2,iv,ie)
!
!if(ielem.eq.1119.or.ielem.eq.1120)print*,'aflim',aflim(2:3,ielem),iv,uvtx,vvtx,uvtx**2+vvtx**2
!
endif
!
!...Get impedence coefficient...
!
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
!
!dux= vlave(1, ipq(iv))-unkno(1, 2, ielem)
!duy= vlave(2, ipq(iv))-unkno(1, 3, ielem)
divu = dux*lnvp(1, iv) + duy*lnvp(2, iv)
!
if(divu.le.0.d0)then
deltu = abs(dux*gpn(1, iv) + duy*gpn(2, iv))
else
deltu = 0.d0
endif
deltu = sqrt(dux**2 + duy**2)
murie = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
icemp(ipq(iv)) = icemp(ipq(iv)) + 1
unmpt(1, icemp(ipq(iv)), ipq(iv)) = pvtx
unmpt(2, icemp(ipq(iv)), ipq(iv)) = murie
unmpt(3, icemp(ipq(iv)), ipq(iv)) = unknvq(2, iv)
unmpt(4, icemp(ipq(iv)), ipq(iv)) = unknvq(3, iv)
!
enddo
!
350 enddo  !...(1)ie = 1,nelem!
!
!...Part III: Get the velocity at middle point
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
!
ft(1, 1) = xpf(1 ,2)- xpf(1, 1)
ft(2, 1) = xpf(2, 2)- xpf(2, 1)
!
lenmc = sqrt(ft(1, 1)**2 + ft(2, 1)**2)
!
ft(1, 1) = ft(1, 1)/lenmc
ft(2, 1) = ft(2, 1)/lenmc
!
fn(1, 1) =  ft(2, 1)
fn(2, 1) = -ft(1, 1)
!
!...For linear PM
ft(1, 2) = xpf(1 ,3)- xpf(1, 1)
ft(2, 2) = xpf(2, 3)- xpf(2, 1)
!
lenmc = sqrt(ft(1, 2)**2 + ft(2, 2)**2)
!
ft(1, 2) = ft(1, 2)/lenmc
ft(2, 2) = ft(2, 2)/lenmc
!
othog = abs(ft(1, 2)*fn(1,1)+ ft(2, 2)*fn(2,1))
!
!if(ipf(3).eq.14) print*,'pt 14', lenmc, ifa
!
!if(othog.lt.1.d-4)then
if(munacn(1,1,ipf(3)).lt.1.d-6)then
!
!...(ifa.le.nbfac) For boundary cells...
!
if(ifa.le.nbfac)then
!
presl = unmpt(1, 1, ipf(3))
mufal = unmpt(2, 1, ipf(3))
uvtxl = unmpt(3, 1, ipf(3))
vvtxl = unmpt(4, 1, ipf(3))
!
fnx = fn(1, 1) !...face normal vector
fny = fn(2, 1)
!
!...Ma
!
if(indnd(ipf(3)).eq.0)then
ustar(1, ipf(3)) = uvtxl + presl/mufal*fnx
ustar(2, ipf(3)) = vvtxl + presl/mufal*fny
endif
!
!...Morgan
!ustar(1, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*ftx/(mufal) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!ustar(2, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*fty/(mufal)
!
!
!...(ifa.gt.nbfac) For interior cells...
!
elseif(ifa.gt.nbfac)then
!
!
dshpr(1) = -0.5d0
dshpr(2) =  0.5d0
dshpr(3) =  0.d0
!
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
!
dwav1 = dydr/djaco
dwav2 =-dxdr/djaco
!
presl = unmpt(1, 1, ipf(3))
mufal = unmpt(2, 1, ipf(3))
uvtxl = unmpt(3, 1, ipf(3))
vvtxl = unmpt(4, 1, ipf(3))
!
presr = unmpt(1, 2, ipf(3))
mufar = unmpt(2, 2, ipf(3))
uvtxr = unmpt(3, 2, ipf(3))
vvtxr = unmpt(4, 2, ipf(3))
!
fnx = dwav1 !...face normal vector
fny = dwav2
!
ftx = -fny
fty =  fnx
!...
!...Mar
!
if(indnd(ipf(3)).eq.0)then
ustar(1, ipf(3)) = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fnx
ustar(2, ipf(3)) = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fny
!
!ustar(1, ipf(3)) = (mufal*(uvtxl*fnx+vvtxl*fny) + mufar*(uvtxr*fnx+vvtxr*fny))*fnx/(mufal+mufar) &
!                    - (presr- presl)/(mufal+mufar)*fnx + &
!                  0.5d0*((uvtxl*ftx+vvtxl*fty) + (uvtxr*ftx+vvtxr*fty))*ftx
!ustar(2, ipf(3)) = (mufal*(uvtxl*fnx+vvtxl*fny) + mufar*(uvtxr*fnx+vvtxr*fny))*fny/(mufal+mufar) &
!                    - (presr- presl)/(mufal+mufar)*fny + &
!                  0.5d0*((uvtxl*ftx+vvtxl*fty) + (uvtxr*ftx+vvtxr*fty))*fty
endif
!
!ftx = -fny
!fty =  fnx
!
!...Morgan

!ustar(1, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*ftx/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!    ustar(2, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*fty/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(2, 1, iv)
!
!    if(ipf(3).eq.443) print*,'midlle velocity',ustar(1:2, ipf(3)),ipf(3),mufal+mufar,fnx,fny
endif
!
endif
!
450 enddo

!
end subroutine getvelo_mpt_curvfvm2

!
!...Get the dampening coeffcients FVM...
!
subroutine getnullve_dampc(ipqua, geoel, cocent, coord, ustar, unkno, gqdmp)
use constant
implicit none
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:3, 1:nsize),              intent(in)::cocent
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,  dimension(1:nquad),        intent(out)::gqdmp
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:3, 1:nvqua)::bq
integer,dimension(1:nvqua) :: ipq
real*8::gdampq(1:nvqua)
real*8::unknvq(1:nq, 1:nvqua)
real*8::gdux1,gdux2,gduy1,gduy2,dr,ds,rcq,scq
real*8::smthid, dmplg, du1, du2
!
integer::iv, ie, ielem,ideg
!
!...shape functions
!
do ie = 1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
rcq= cocent(1, ielem) !...mass center...
scq= cocent(2, ielem)
!
xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rcq) 
bq(3, iv) = (yvq(iv)-scq) 
enddo

!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1, nvqua
!
do ideg = 1,mdegr
unknvq(2:3, iv) = unknvq(2:3, iv) + unkno(ideg,2:3,ielem)*bq(ideg, iv)
enddo
!
gdux1 = ustar(1,ipq(iv))-unknvq(2, iv)
gduy1 = ustar(2,ipq(iv))-unknvq(3, iv)
!
gdux2 = ustar(1,ipq(iv))-unkno(1, 2, ielem)
gduy2 = ustar(2,ipq(iv))-unkno(1, 3, ielem)
!
du1 = sqrt(gdux1**2+gduy1**2)
du2 = sqrt(gdux2**2+gduy2**2)
!
if(du2.le.1.d-6)then
smthid = 0.d0
else
smthid = du1/(du2+1.d-6)
endif

gdampq(iv) = smthid!min(1.d0, .5d0, 1.d0*sqrt(gdux1**2+gduy1**2)/sqrt(gdux2**2+gduy2**2))
!
!if(ie==1) print*,'unkno',unkno(1:3,2,ie),unkno(1:3,3,ie)
!if(ie==1) print*,'basis',bq(1:3,iv),iv
!if(ie==1) print*,'basis2',ustar(1:2,ipq(iv)),unknvq(2:3, iv)
!if(ie==1) print*,'basis3',unkno(1, 2:3, ielem)

!if(ie==1) print*,'gdmp',iv,&
!sqrt(gdux1**2+gduy1**2)/sqrt(gdux2**2+gduy2**2),sqrt(gdux1**2+gduy1**2),sqrt(gdux2**2+gduy2**2)
!
enddo
!
dmplg = maxval(gdampq(1:8))
gqdmp(ie) = min(0.5d0,dmplg)
!gqdmp(ie) = minval(gdampq(1:8))! 0.8d0!minval(gdampq(1:8))
!
!if(ie==2192) print*,'gdmp2',gqdmp(ie)
!
enddo
!
!gqdmp =0.d0
!
!
!...get null mode...
!
end subroutine getnullve_dampc
!
!...Get the dampening coeffcients DGM...
!
subroutine getnullve_dampcdgp(intfac, iptri, ipqua, geoel,  coord, ustar, unkno, gqdmp)
use constant
implicit none
integer*4,dimension(1:nifai,1:nafac), intent(in) ::intfac
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,  dimension(1:nquad),        intent(out)::gqdmp
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:3, 1:nvqua)::bq
integer,dimension(1:nvqua) :: ipq
real*8::gdampq(1:nvqua)
real*8::unknvq(1:nq, 1:nvqua)
real*8,dimension(1:3, 1:nsize)::cocent
real*8,dimension(1:ndegr,1:nq,1:nsize)::unknu
real*8::gdux1,gdux2,gduy1,gduy2,dr,ds,rcq,scq
!
integer::iv, ie, ielem,ideg
!
!
!
unknu = 0.d0
unknu = unkno
!
!call getgraduls_lag(intfac, iptri, ipqua, unknu ,cocent, coord)
!
!...shape functions
!
do ie = 1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
rcq= cocent(1, ielem) !...mass center...
scq= cocent(2, ielem)
!
xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rcq) 
bq(3, iv) = (yvq(iv)-scq) 
enddo

!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1, nvqua
!
do ideg = 1,mdegr
unknvq(2:3, iv) = unknvq(2:3, iv) + unknu(ideg,2:3,ielem)*bq(ideg, iv)
enddo
!
gdux1 = ustar(1,ipq(iv))-unknvq(2, iv)
gduy1 = ustar(2,ipq(iv))-unknvq(3, iv)
!
gdux2 = ustar(1,ipq(iv))-unknu(1, 2, ielem)
gduy2 = ustar(2,ipq(iv))-unknu(1, 3, ielem)

gdampq(iv) = min(1.d0, .7d0, 15.d0*sqrt(gdux1**2+gduy1**2)/sqrt(gdux2**2+gduy2**2))
!
!if(ie==2192) print*,'gdmp',sqrt(gdux1**2+gduy1**2)/sqrt(gdux2**2+gduy2**2),iv
!
enddo
!
gqdmp(ie) = minval(gdampq(1:8))
!
!if(ie==2192) print*,'gdmp2',gqdmp(ie)
!
enddo
!
!gqdmp =0.d0
!
!
!...get null mode...
!
end subroutine getnullve_dampcdgp
!
!...Get the nodal velocity based on Simpson integration in FVM...
!
subroutine getndvelo_lagfvm_simpson(gflag,gelag,cocent, gelagq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, fstart, fstarq, aflim, afvec,vnulq, itime)
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
real*8,dimension(1:3,1:nsize),               intent(in)::cocent
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:9, 1:ntria),   intent(out)::fstart !...Riemann forces
real*8,dimension(1:ndimn,1:12, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(out)::vnulq

integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)
integer, dimension(3, 4)::fglvq, fglgq
integer, dimension(3, 3)::fglvt, fglgt

!...local real array
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,  dimension(1:nquad)::gqdmp
real*8::munaci(2, 2)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2
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
allocate (munaclt(1:2, 1:2, 1:9, 1:ntria), munault(1:ndimn, 1:9,  1:ntria),&
snsigmlt(1:ndimn, 1:9,  1:ntria))
allocate (munaclq(1:2, 1:2, 1:12, 1:nquad), munaulq(1:ndimn, 1:12,  1:nquad),&
snsigmlq(1:ndimn, 1:12,  1:nquad))
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

!...Local gauss point No. of any gauss point in one face...
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

!...Local vertex No. of gauss points in one unit ...
fglvt(1, 1) = 1;  fglvt(2, 1) = 2; fglvt(3, 1) = 4;
fglvt(1, 2) = 2;  fglvt(2, 2) = 3; fglvt(3, 2) = 5;
fglvt(1, 3) = 3;  fglvt(2, 3) = 1; fglvt(3, 3) = 6;

!...Local gauss point No. of any gauss point in one face...
fglgt(1, 1) = 1;  fglgt(2, 1) = 2; fglgt(3, 1) = 7;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4; fglgt(3, 2) = 8;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6; fglgt(3, 3) = 9;
!
!...Zero out vlave
!
vlave = 0.d0
indnd = 0
!
!...Mark the boundary nodes...
!
if(ncase.eq.2)then !...Smooth shockless Noh
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
indnd(ipf(1:nvfac)) = 1
enddo
endif

!...BC with 25 is prescribed.
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
indnd(ipf(1:nvfac)) = 1
endif
enddo

!...Get averaged velocity at nodes...
!call getvlavenew(iptri, ipqua, geoel, vlave, unkno, aflim, afvec)
 call getvlavefvm(iptri, ipqua, geoel, cocent,coord, vlave, unkno, aflim, afvec)
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
!...Part II: Loop to get the information from Riemann solver
!
!call getnullve_dampc(ipqua, geoel, cocent, coord, ustar, unkno, gqdmp)
 vnulq = 0.d0

do iloop= 1, 6!!
!
!call getnullmd_quadlp2(ipqua, geoel, vnulq, ustar, coord, gqdmp)
!...Give vlave
vlave= ustar

!call getnullmd_quadlp(ipqua, geoel, vnulq, ustar, coord, gqdmp)
!vnulq = 0.d0

!...Zero out munacn
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0

!...Quad
if(nquad.gt.0) call getriem_quadfvm_simpson(ipqua, geoel, cocent,  gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec, vnulq)

!...Boundary condition
!...The boundary condition for FVM is not implmented now.

!...Update the velocity at the vertex
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
endif
enddo

!...Get the velociy at the middle point (This is not used with |a \dot n|=1.d0)
!call getvelo_mpt_curv(ustar,geoel,gelag,gelagq,intfac,iptri,ipqua,coord,coold,unkno,indnd, aflim, afvec, vlave, vnulq)

!...Get the vertex velocity at the boundary...
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

!...Specify the boundary velocity
 if(ncase.eq.13)then !...Saltzman
 ustar(2,ipf(1:nvfac)) = 0.d0
 ustar(1,ipf(1:nvfac)) = 1.d0
 endif
endif

!...Prescribe the linear velocity at middle point
!ustar(1, ipf(3)) = 0.5d0*(ustar(1, ipf(1))+ustar(1, ipf(2)))
!ustar(2, ipf(3)) = 0.5d0*(ustar(2, ipf(1))+ustar(2, ipf(2)))

!...Impose exact solution at the boundary for TGV
 if(ncase.eq.-1)then
 ustar(1, ipf(1:nvfac)) = sin(pi*coord(1,ipf(1:nvfac)))*cos(pi*coord(2,ipf(1:nvfac)))
 ustar(2, ipf(1:nvfac)) =-cos(pi*coord(1,ipf(1:nvfac)))*sin(pi*coord(2,ipf(1:nvfac)))
 endif
900 enddo

!...Get the vnulq based on velocity filter
if(iloop.eq.1) call getnullve_dampc(ipqua, geoel, cocent, coord, ustar, unkno, gqdmp)
call getnullmd_quadlp(ipqua, geoel, vnulq, ustar, coord, gqdmp)

enddo !iloop
!
!print*,'ustar1',ustar(1,1:10)
!...Get the linear velocity at middle point...
!do 920 ifa = 1 , nafac
!ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!ustar(1, ipf(3)) = 0.5d0*(ustar(1, ipf(1))+ustar(1, ipf(2)))
!ustar(2, ipf(3)) = 0.5d0*(ustar(2, ipf(1))+ustar(2, ipf(2)))
!920 enddo
!
!print*,'ustar2',ustar(1:2,1:10)
!
!...Part III: Get the Riemann forces for face integral...
!

!...III.1: Get the Riemann forces associated with every node...

!...Quad
do ie = 1, nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
do ifa = 1, 4
do ig =1, nvfac
fstarq(1, fglgq(ig, ifa), ie) = snsigmlq(1, fglgq(ig, ifa), ie) +&
munaclq(1,1, fglgq(ig, ifa), ie)*ustar(1, ipq(fglvq(ig, ifa)))+&
munaclq(2,1, fglgq(ig, ifa), ie)*ustar(2, ipq(fglvq(ig, ifa)))-&
munaulq(1, fglgq(ig, ifa), ie)
!
fstarq(2, fglgq(ig, ifa), ie) = snsigmlq(2, fglgq(ig, ifa), ie) +&
munaclq(2,2,fglgq(ig, ifa), ie)*ustar(2, ipq(fglvq(ig, ifa)))+&
munaclq(1,2,fglgq(ig, ifa), ie)*ustar(1, ipq(fglvq(ig, ifa)))-&
munaulq(2, fglgq(ig, ifa), ie)
enddo
enddo
enddo
!
deallocate (munacn, bpres, fpres)
deallocate (munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lagfvm_simpson
!
!...Riemann input for curved quad using simpson rule in FVM...
!
subroutine getriem_quadfvm_simpson(ipqua, geoel,cocent,gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec, vnulq)
use constant
implicit none

!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:nsize),               intent(in)::cocent
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),            intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),             intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(in):: vnulq

!
real*8, dimension(1:2, 1:2, 1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
!
real*8, dimension(1:2, 1:2, 1:12, 1:nquad),       intent(out)::munaclq
real*8, dimension(1:ndimn, 1:12,  1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:12,  1:nquad), intent(out)::snsigmlq

!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig

!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(nvfac, 4)::fglvq, fglgq

!...local real array
real*8,dimension(1:ndegr)::bg
real*8,dimension(1:nq,1:nvfac)::unkng
real*8::aujmp(1:3, 1:nvfac)
real*8::vnorm(1:3, 1:12)
real*8::sigmg(1:2, 1:2, 1:nvfac)
real*8,dimension(1:nvfac)::murie
real*8,dimension(1:nvqua):: xvq,  yvq
real*8::lnvp(2,nvqua)

!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::sdimp
real*8::rhog,rhomg,ug,vg,eg, pg
real*8::dux,duy,deltu,gpnx,gpny
real*8::dr, ds, rc, sc, rg, sg,rcv,scv
real*8:: dudr, duds, dvdr, dvds
real*8::acnx, acny,divu,divgu
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I:Specify some gauss points...
!

!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;

!...Local gauss point No. of any gauss point in one face...
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

dr = 1.d0
ds = 1.d0
!
!...Part II: Loop over every triangle...
!
do 250 ie = 1,nquad !...(1)ie = 1,nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria

!...shape functions
rc= cocent(1, ielem) !...mass center...
sc= cocent(2, ielem)

xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))

!...Give the normal vector of every face...
vnorm(1:3,  1:12) = gelagq(1:3, 1:12, ie);

!...Get the LN for one vertex for identifying one cell compressed or expanded...
lnvp = 0.d0

lnvp(1:2, 1) = vnorm(1:2, 1)*vnorm(3, 1) + vnorm(1:2, 8)*vnorm(3, 8)
lnvp(1:2, 2) = vnorm(1:2, 2)*vnorm(3, 2) + vnorm(1:2, 3)*vnorm(3, 3)
lnvp(1:2, 3) = vnorm(1:2, 4)*vnorm(3, 4) + vnorm(1:2, 5)*vnorm(3, 5)
lnvp(1:2, 4) = vnorm(1:2, 6)*vnorm(3, 6) + vnorm(1:2, 7)*vnorm(3, 7)
!
lnvp(1:2, 5) = vnorm(1:2, 9)*vnorm(3, 9)
lnvp(1:2, 6) = vnorm(1:2, 10)*vnorm(3, 10)
lnvp(1:2, 7) = vnorm(1:2, 11)*vnorm(3, 11)
lnvp(1:2, 8) = vnorm(1:2, 12)*vnorm(3, 12)
!
lnvp(:,1:4) = 1.d0/6.d0*lnvp(:,1:4)
lnvp(:,5:8) = 4.d0/6.d0*lnvp(:,5:8)
!...The metrics of the divergence
divgu = 0.d0
do iv = 1, 8
divgu = divgu + vlave(1, ipq(iv))*lnvp(1, iv) + vlave(2, ipq(iv))*lnvp(2, iv)
enddo

!...cell averaged value...
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
endif

uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)

rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

!...II.1: Loop over every face...
do ifa =1, 4

!...zero out unkno
unkng = 0.d0

do ig   = 1, nvfac
rg = xvq(fglvq(ig, ifa))
sg = yvq(fglvq(ig, ifa))
!
bg(1) = 1.d0
bg(2) = (rg-rc)
bg(3) = (sg-sc)
!
do ideg = 1,mdegr
unkng(1:nq-1, ig) = unkng(1:nq-1, ig) + unkno(ideg,1:nq-1,ielem)*bg(ideg)
enddo

!...Pressure
unkng(4, ig) = pctr + unkno(2,4,ielem)*bg(2) + unkno(3,4,ielem)*bg(3)

!
if(ndens.eq.1)then
rhog  = 1.d0/unkng(1, ig)
endif
!
ug = unkng(2, ig)
vg = unkng(3, ig)
pg = unkng(4, ig)

!...Limiter
if(nlimi.eq.6)then
 if(ndens.eq.1)then
 rhomg = rhomc + aflim(1, ielem)*(unkng(1, ig) - rhomc)
 rhog = 1.d0/rhomg
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

!...Remove the null mode of velocity
unkng(2, ig) = ug -  vnulq(1,fglvq(ig, ifa),ie)
unkng(3 ,ig) = vg -  vnulq(2,fglvq(ig, ifa),ie)
endif

!...Get stress tensor at vertices
sigmg(1, 1, ig) = -pg
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pg

!...Get the a_c (unit vector)
aujmp(1:2, ig) = vlave(1:2, ipq(fglvq(ig, ifa))) - unkng(2:3, ig)
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

!...Get sound speed at the cell cenetr...
sdctr = sqrt(gamlg*pctr/rhoct)
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedance coefficient...
do ig   = 1, nvfac
 dux= vlave(1, ipq(fglvq(ig, ifa)))-unkng(2, ig)
 duy= vlave(2, ipq(fglvq(ig, ifa)))-unkng(3, ig)
 deltu = sqrt(dux**2 + duy**2)
!
 gpnx = vnorm(1, fglgq(ig, ifa))
 gpny = vnorm(2, fglgq(ig, ifa))
 dux = vlave(1, ipq(fglvq(ig, ifa)))-unkno(1, 2, ielem)
 duy = vlave(2, ipq(fglvq(ig, ifa)))-unkno(1, 3, ielem)
!
 divu = dux*lnvp(1, fglvq(ig, ifa)) + duy*lnvp(2, fglvq(ig, ifa))
if(divgu.le.0.d0)then
! deltu = abs(dux*gpnx + duy*gpny)
else
! deltu = 0.d0
endif
 !deltu = abs(dux*gpnx + duy*gpny)
 sdimp = sdctr!max(1d0, sdctr)
 murie(ig) = rhoct*sdimp !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! murie(ig) = rhoct*(slpdu*deltu*0.5d0 + sqrt((slpdu*deltu*0.5d0)**2+sdimp**2))
!
enddo
!
!print*,'sdctr',ielem,sdimp,sdctr,murie(1),rhoct,pctr

!...Get the summed denominator and numerator coefficients sum(mu*n*a_c)
do ig  = 1, nvfac

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

munacu(1:2, ipq(fglvq(ig, ifa))) = munacu(1:2, ipq(fglvq(ig, ifa))) + munacu_rie(1:2)

snsigm(1:2,ipq(fglvq(ig, ifa))) = snsigm(1:2, ipq(fglvq(ig, ifa))) + snsigm_rie(1:2)!

!
munaclq(1:2, 1, fglgq(ig, ifa), ie) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, fglgq(ig, ifa), ie) =  munacn_rie(1:2, 2)

munaulq(1:2, fglgq(ig, ifa), ie)    =  munacu_rie(1:2)

snsigmlq(1:2,fglgq(ig, ifa), ie)   = snsigm_rie(1:2)
!
enddo
!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,nquad

end subroutine getriem_quadfvm_simpson
!
!...Face integral using gauss quadrature distribution on cuvred quads for FVM...
!
subroutine rhsifacedg_lagfvmq_simpson(ipqua,  unkno, ustar, fstarq, gelagq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:12,1:nquad),       intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::rhsel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize), intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ielem
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer, dimension(3, 4)::fglvq, fglgq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::vnorm(3, 12)
real*8::xvq(nvqua), yvq(nvqua)
!
real*8::weigh(nvfac), posi(1,nvfac)
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
if(nvfac.eq.3)then
!
!...Local vertex No. of gauss points in one unit ...
!
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
posi(1, 1) = -1.d0; posi(1 ,2 )= 1.d0; posi(1 ,3 )= 0.d0
!
weigh(1) = 1.d0/6.d0; weigh(2) = 1.d0/6.d0; weigh(3) = 4.d0/6.d0
!
!...Local gauss point No. of any gauss point in one face...
!
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;
!
endif

!...Linear cell
!weigh(1) = 1.d0/2.d0; weigh(2) = 1.d0/2.d0;
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
vnorm(1:3,  1:12) = gelagq(1:3, 1:12, ie)
!
do ifa = 1, 4
!
do ig   = 1, nvfac
!
wi  = weigh(ig)
!
!if(ie.eq.1) print*,'bg',ig,ifa,infag(ig, ifa), rg, sg,bg(1:3,3-ig, infag(ig, ifa))
!
!
!...The vertex constituting one cell...
!d0
gpnx = vnorm(1, fglgq(ig, ifa))
gpny = vnorm(2, fglgq(ig, ifa))
gpsa = vnorm(3, fglgq(ig, ifa))
!
!...Distribute to every corner...
!
ulnpn(1)  = ulnpn(1)+&
ustar(1, ipq(fglvq(ig, ifa)))*gpnx*gpsa*weigh(ig) +&
ustar(2, ipq(fglvq(ig, ifa)))*gpny*gpsa*weigh(ig)
!
!
plnpn(1, 1)= plnpn(1, 1)   +&
fstarq(1, fglgq(ig, ifa), ie)*weigh(ig)
!
plnpn(2, 1)= plnpn(2, 1)  +&
fstarq(2, fglgq(ig, ifa), ie)*weigh(ig)
!
elnpn(1)   = elnpn(1)+&
ustar(1, ipq(fglvq(ig, ifa)))*fstarq(1, fglgq(ig, ifa), ie)*weigh(ig) +&
ustar(2, ipq(fglvq(ig, ifa)))*fstarq(2, fglgq(ig, ifa), ie)*weigh(ig)
!
!if(ie==1862) print*,'rhs iface idegr',ustar(1, ipq(fglvq(ig, ifa))),gpnx,gpny,gpsa,bg(1:ndegr,ig),ipq(fglvq(ig, ifa)),ig,ifa
!
enddo
!
enddo
!
rhsel(1, 1, ielem) =  ulnpn(1)
rhsel(1, 2, ielem) =  plnpn(1, 1)
rhsel(1, 3, ielem) =  plnpn(2, 1)
rhsel(1, 4, ielem) =  elnpn(1)
!
! if(ie==23) print*,'rhs iface',rhsel(1:3, 1, ie), ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)
550 enddo
!
end subroutine rhsifacedg_lagfvmq_simpson
!
!...subroutine: Calculate the averaged velocity for hybrid grids FVM...
!
subroutine getvlavefvm(iptri, ipqua, geoel, cocent,coord, vlave, unkno, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:nsize),              intent(in)::cocent
real*8,dimension(1:ndimn,1:npoin), intent(in) ::coord
real*8,dimension(1:ndimn,1:npoin),           intent(out)::vlave
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!...Local integer
integer::ie ,ielem, iv, ipoin
integer:: ideg
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:3, 1:nvtri)::bt
real*8,dimension(1:3, 1:nvqua)::bq
real*8,dimension(1:nq,1:nvtri)::unknvt
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nvtri):: xv,  yv
real*8,dimension(1:nvqua):: xvq, yvq
real*8::rcq, scq
real*8:: dudr, duds, dvdr, dvds
!...local real number
!
real*8,allocatable:: cnsup(:)
!
allocate (cnsup(1:npoin))
!
!...Zero out vlave, cnsup...cnsup:cell No surrounding one point...
!
vlave = 0.d0
cnsup = 0.d0
!
!...Quad...
!
do 300 ie = 1,nquad !...(1)ie = 1,nelem
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
!...shape functions
!
rcq= cocent(1, ielem) !...mass center...
scq= cocent(2, ielem)
!
xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))

do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rcq)
bq(3, iv) = (yvq(iv)-scq)
enddo
!
!...
!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(2:3, iv) = unknvq(2:3, iv) + unkno(ideg,2:3,ielem)*bq(ideg, iv)
!
enddo
!
if(nlimi.eq.6)then
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
unknvq(2, iv) = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
unknvq(3, iv) = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
!
endif
! if(ip(iv)==36) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
enddo
!
!...Accumulate nodal velocity...
vlave(1, ipq(1:nvqua)) = vlave(1, ipq(1:nvqua)) + unknvq(2, 1:nvqua)
vlave(2, ipq(1:nvqua)) = vlave(2, ipq(1:nvqua)) + unknvq(3, 1:nvqua)
!
!...Accumualte to get the no of cell surrounding one vertex...
!
cnsup(ipq(1:nvqua)) = cnsup(ipq(1:nvqua)) + 1.d0
!
300 enddo  !...(1)ie = 1,nelem

!
!...Get the averaged  nodal velocity...
!
do ipoin = 1, npoin
vlave(1:ndimn, ipoin) = vlave(1:ndimn, ipoin)/cnsup(ipoin)
enddo
!
!
deallocate(cnsup)
!
end subroutine getvlavefvm
!
!...Symmetry preserving techniques(vector limiting) based on face neighbours...
!
subroutine barthlimitface_vectfvm(geoel, cocent,coord, ustar, unkno, ipqua, bface, intfac, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:3,1:nsize),                 intent(in)::cocent
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),           intent(in):: ipqua
integer*4,dimension(1:nbfai,1:nbfac),           intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac

!...Local
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp,ivf
integer:: iel, ier
integer:: ielem, jelem, istor
real*8,  dimension(1:nq+1, 1:2):: unkcl, unkcr
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:nq+1) :: unctr
real*8, dimension(1:nq+1, 1:ncell) :: unmax_new, unmin_new
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
real*8,dimension(1:nq+1, 1:nsize)   ::aflim
!

real*8::eps,c00,c05,c10,c20
real*8:: uloca, vloca, ucjel, vcjel
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctrl, vctrl,uctrr,vctrr
real*8:: uctr, vctr
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct, rhom
real*8:: afbar
real*8:: dudr, duds, dvdr, dvds
real*8:: etcx, etcy, xicx, xicy, fixi, fiet
real*8:: lamda1, lamda2
real*8:: mapd, mapt, matra, matrb, matrc, matrd, umap, vmap
real*8:: jf11, jf12, jf21, jf22
real*8:: fi11, fi12, fi21, fi22
real*8:: dxdr, dxds, dydr, dyds
real*8:: rm, rp, sm, sp, r, s
real*8:: drdx, drdy, dsdx, dsdy, jacom
real*8:: dudx, dudy, dvdx, dvdy, lmat1, lmat2
real*8:: dudxy, dumax, dumin
real*8:: rhomc,sdctr,macel,delu,volel
real*8:: signlim, vemag, afmag
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Some preparing work
!

!...Coloring the boundary node
indbd = 0
eps = 1.e-6
mapmt = 0.d0
!
!...Part II: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
dudr = unkno(2, 2, ielem)
duds = unkno(3, 2, ielem)
dvdr = unkno(2, 3, ielem)
dvds = unkno(3, 3, ielem)
!
dudx = dudr
dudy = duds
dvdx = dvdr
dvdy = dvds

!...LANL
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy
!
!matra = unkno(1, 2, ielem)
!matrb = unkno(1, 3, ielem)
if(sqrt(matra**2+matrb**2).gt.1.d-8)then
!mapmt(1, 1, ielem) = matra/sqrt((matra)**2 + matrb**2)
!mapmt(1, 2, ielem) = matrb/sqrt((matra)**2 + matrb**2)
!mapmt(2, 1, ielem) =-matrb/sqrt((matra)**2 + matrb**2)
!mapmt(2, 2, ielem) = matra/sqrt((matra)**2 + matrb**2)

else
!mapmt(1, 1, ielem) = 1.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 1.d0

endif
!
!...Maire
!
!matra = dxdr**2 + dydr**2
!matrb = dxdr*dxds + dydr*dyds
!matrc = matrb
!matrd = dxds**2 + dyds**2

!...eigenvalues...
mapt  = matra + matrd
mapd  = matra*matrd - matrb*matrc
lamda1 = 0.5d0*mapt + sqrt(0.25d0*mapt**2-mapd)
lamda2 = 0.5d0*mapt - sqrt(0.25d0*mapt**2-mapd)
!!
if(abs(matrc).gt.1.d-6)then
!
if(matrc.gt.0.d0)then
!
lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
mapmt(1, 2, ielem) =          matrc/lmat1
mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
mapmt(2, 2, ielem) =          matrc/lmat2

elseif(matrc.lt.0.d0)then

lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
mapmt(2, 2, ielem) =          matrc/lmat1
mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
mapmt(1, 2, ielem) =          matrc/lmat2
endif
!
else
!
mapmt(1, 1, ielem) = 0.d0
mapmt(1, 2, ielem) = 0.d0
mapmt(2, 1, ielem) = 0.d0
mapmt(2, 2, ielem) = 0.d0
endif
!
enddo
!
!...Part III: Get the maximum and minimum around one cell
!
unmax_new = 1d-10
unmin_new = 1d10
!...III.1
do ifa = 1, nbfac
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
uctrl = unkno(1, 2,  intfac(1, ifa))
vctrl = unkno(1, 3,  intfac(1, ifa))
!
unkcl(2, 1) = uctrl*mapmt(1, 1, iel) + vctrl*mapmt(1, 2, iel)
unkcl(3, 1) = uctrl*mapmt(2, 1, iel) + vctrl*mapmt(2, 2, iel)
!
if(bface(3, ifa).eq.22)then
if(bface(4,ifa).eq.221)then
uctrr = uctrl
vctrr =-vctrl
elseif(bface(4,ifa).eq.222)then
uctrr =-uctrl
vctrr = vctrl
endif
endif
!
unkcl(2, 2) = uctrr*mapmt(1, 1, iel) + vctrr*mapmt(1, 2, iel)
unkcl(3, 2) = uctrr*mapmt(2, 1, iel) + vctrr*mapmt(2, 2, iel)
!
do iq=2, 3
unmax_new(iq, iel) = max(unkcl(iq, 1),unkcl(iq, 2), unmax_new(iq, iel))
unmin_new(iq, iel) = min(unkcl(iq, 1),unkcl(iq, 2), unmin_new(iq, iel))
enddo
enddo

!...III.2
do ifa = nbfac+1, nafac
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
uctrl = unkno(1, 2,  intfac(1, ifa))
vctrl = unkno(1, 3,  intfac(1, ifa))
!
unkcl(2, 1) = uctrl*mapmt(1, 1, iel) + vctrl*mapmt(1, 2, iel)
unkcl(3, 1) = uctrl*mapmt(2, 1, iel) + vctrl*mapmt(2, 2, iel)
!
uctrr = unkno(1, 2,  intfac(2, ifa))
vctrr = unkno(1, 3,  intfac(2, ifa))
!
unkcl(2, 2) = uctrr*mapmt(1, 1, iel) + vctrr*mapmt(1, 2, iel)
unkcl(3, 2) = uctrr*mapmt(2, 1, iel) + vctrr*mapmt(2, 2, iel)
!
do iq=2, 3
unmax_new(iq, iel) = max(unkcl(iq, 1),unkcl(iq, 2), unmax_new(iq, iel))
unmin_new(iq, iel) = min(unkcl(iq, 1),unkcl(iq, 2), unmin_new(iq, iel))
enddo
!
unkcl(2, 1) = uctrl*mapmt(1, 1, ier) + vctrl*mapmt(1, 2, ier)
unkcl(3, 1) = uctrl*mapmt(2, 1, ier) + vctrl*mapmt(2, 2, ier)
!
unkcl(2, 2) = uctrr*mapmt(1, 1, ier) + vctrr*mapmt(1, 2, ier)
unkcl(3, 2) = uctrr*mapmt(2, 1, ier) + vctrr*mapmt(2, 2, ier)
!
do iq=2, 3
unmax_new(iq, ier) = max(unkcl(iq, 1),unkcl(iq, 2), unmax_new(iq, ier))
unmin_new(iq, ier) = min(unkcl(iq, 1),unkcl(iq, 2), unmin_new(iq, ier))
enddo
enddo
!
!...Part III: Impose limiter for mapped u and v
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
rc = cocent(1, ielem)
sc = cocent(2, ielem)
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xpq(1, iv)-rc)
bq(3, iv) = (xpq(2, iv)-sc)
enddo

!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
!
!...New mapped velocity components u and v
!
umap = uvtx*mapmt(1, 1, ielem) + vvtx*mapmt(1, 2, ielem)
vmap = uvtx*mapmt(2, 1, ielem) + vvtx*mapmt(2, 2, ielem)
!
unknvq(2, iv) = umap
unknvq(3, iv) = vmap
!
enddo
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)

!...New mapped velocity components u and v
umap = uctr*mapmt(1, 1, ielem) + vctr*mapmt(1, 2, ielem)
vmap = uctr*mapmt(2, 1, ielem) + vctr*mapmt(2, 2, ielem)
!
unctr(2) = umap
unctr(3) = vmap


!...Limiting
do iv = 1, nvqua
do iq = 2, 3
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo
enddo
!
do iq = 2, 3
aflim(iq, ielem) = minval(alfa(iq, 1:8)) !...Only consider the 4 vertices
enddo
!
enddo
!
!...Part 5: Transfer back the limiter to the gobal frame...
!
do ie = 1,nquad
ielem = ie + ntria

ipq(1:nvqua) = ipqua(1:nvqua,ie)

rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
xicx = mapmt(1, 1, ielem); xicy = mapmt(1, 2, ielem)
etcx = mapmt(2, 1, ielem); etcy = mapmt(2, 2, ielem)
!
fixi = aflim(2, ielem) ; fiet = aflim(3, ielem)
!
fi11 = xicx**2*fixi + etcx**2*fiet
fi12 = xicx*xicy*fixi + etcx*etcy*fiet
fi21 = fi12
fi22 = xicy**2*fixi + etcy**2*fiet
!!
afvec(1, 1, ielem) = fi11
afvec(1, 2, ielem) = fi12
afvec(2, 1, ielem) = afvec(1, 2, ielem)
afvec(2, 2, ielem) = fi22
!
enddo
end subroutine barthlimitface_vectfvm
!
!...Subroutine for barth limiter based on face neighbours...
!
subroutine barthlimitface_fvm(geoel, cocent, coord, coold, ustar, unkno, ipqua, bface,intfac, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:3,1:nsize),                 intent(in)::cocent
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord, coold
real*8,dimension(1:nq+1, 1:nsize),             intent(inout)::aflim
integer,  dimension(1:nvqua,1:nquad),          intent(in):: ipqua
integer*4,dimension(1:nbfai,1:nbfac),           intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
real*8,dimension(1:2, 1:2, 1:nsize),          intent(inout)::afvec
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,iefa
integer:: ielem,iel,ier
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(1:ndegr, 1:nvqua),bqv(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xpqi
real*8,  dimension(1:nq+1, 1:2):: unkcl, unkcr
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy,rcv,scv
real*8:: rhov, rhoct, rhom
real*8:: afbar
real*8:: dudr, duds, dvdr, dvds
real*8:: r, s, rhoi, rhon
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:(2+nvfac), ifa)) = 1

enddo
!
eps = 1.e-6
!
!...Store the maximum and minimum values surrounding one cell based on face neighbouring...
!
unmax_new = 1d-10
unmin_new = 1d10
!
do ifa = nbfac+1, nafac
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
do iefa = 1, 2
if(ndens.eq.1)then
rhom = unkno(1, 1, intfac(iefa, ifa))
elseif(ndens.eq.2)then
rhom = 1.d0/unkno(1, 1, intfac(iefa, ifa))
elseif(ndens.eq.3)then
rhom = 1.d0/unkno(1, 1, intfac(iefa, ifa))
endif

uctr = unkno(1, 2,  intfac(iefa, ifa))
vctr = unkno(1, 3,  intfac(iefa, ifa))
ectr = unkno(1, 4,  intfac(iefa, ifa))
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!

unkcl(1, iefa)   = unkno(1, 1,  intfac(iefa, ifa))
unkcl(2:3, iefa) = unkno(1, 2:3,  intfac(iefa, ifa))
unkcl(nq, iefa) = pctr
unkcl(nq+1,iefa) = ectr
!
enddo !...iefa
!
do iq=1, nq+1
unmax_new(iq, iel) = max(unkcl(iq, 1),unkcl(iq, 2), unmax_new(iq, iel))
unmin_new(iq, iel) = min(unkcl(iq, 1),unkcl(iq, 2), unmin_new(iq, iel))
!
unmax_new(iq, ier) = max(unkcl(iq, 1),unkcl(iq, 2), unmax_new(iq, ier))
unmin_new(iq, ier) = min(unkcl(iq, 1),unkcl(iq, 2), unmin_new(iq, ier))
enddo
enddo
!
!...Part 2: Impose limiter
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= cocent(1, ielem) !...mass center...
sc= cocent(2, ielem)
!
xvq(1:nvqua) = coord(1, ipq(1:nvqua))
yvq(1:nvqua) = coord(2, ipq(1:nvqua))
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)
bq(3, iv) = (yvq(iv)-sc)
enddo
!
!...Cell averaged...
!
rhom = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq-1, iv) = unknvq(1:nq-1, iv) + unkno(ideg,1:nq-1,ielem)*bq(ideg, iv)
enddo

!...Pressure...
unknvq(4, iv) = pctr  + unkno(2,4,ielem)*bq(2, iv) + unkno(3,4,ielem)*bq(3, iv)
!
rhov = 1.d0/unknvq(1, iv)
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
unknvq(1, iv) = 1.d0/rhov
enddo
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq)  = pctr
!
do iv = 1, nvqua
do iq = 1, nq

dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo
enddo

!...Choose the minimum one
do iq = 1,nq
aflim(iq, ielem) = minval(alfa(iq, 1:8))
enddo
!
enddo
!
!...Part 2.1: Impose symmetry preserving limiter for vector...
!
call barthlimitface_vectfvm(geoel, cocent,coord, ustar, unkno, ipqua, bface,intfac, afvec)
!
end subroutine barthlimitface_fvm
