!
!...Find mass center in geoel for curved cell...
!
subroutine getgeoel_lag_mc(inpoel, iptri, ipqua, geoel, coord)
use constant
implicit none
integer*4,dimension(1:nvtri, 1:ntria),intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:ngeel, 1:nsize), intent(inout)::geoel
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
do ishp = 1, nptri
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
do ishp = 1, nptri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo
!
if(ncase.eq.1)then     !...TGV
rhog = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rhog = 1.d0
elseif(ncase.eq.3)then !...
rhog = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
!
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rhog = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rhog = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
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
!
!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then
!
rhog  = 1.d0
!
else
!
rhog  = 0.125d0
!
endif
!
elseif(ncase.eq.7)then !...sedov
!
rhog = 1.d0
!
elseif(ncase.eq.8)then !...Gresho
!
rhog = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
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
!
if(xc.lt.1.d0)then
rhog  = 1.d0
else
if(yc.gt.1.5d0)then
rhog  = 0.1d0
else
rhog = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rhog =1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rhog  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rhog = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rhog = 1.d0
!
elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rhog = 1.d0
!
endif
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
!volel = volel + djaco
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
if(ncase.eq.1)then     !...TGV
rhog = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rhog = 1.d0
elseif(ncase.eq.3)then !...
rhog = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
!
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rhog = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rhog = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
r=0.d0;s=0.d0
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo

!if(ielem.eq.1601)print*,'mass',rhog,sqrt(xgaus**2+ygaus**2)
!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then

rhog  = 1.d0
else
rhog  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
rhog = 1.d0
elseif(ncase.eq.8)then !...Gresho
!
rhog = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
r=0.d0;s=0.d0
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo
!
if(xc.lt.1.d0)then
rhog  = 1.d0
else
if(yc.gt.1.5d0)then
rhog  = 0.1d0
else
rhog = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rhog = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rhog  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rhog = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rhog = 1.d0

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rhog = 1.d0
!
endif
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
!if(ielem.eq.1) print*,'bad',djaco,xpq(1, 1:npqua)
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
!print*,'geoel',geoel(1:8,1)
!
end subroutine getgeoel_lag_mc
!
!...Find terms in geoel for high-order...
!
subroutine getgeoel_lagmc_ho(inpoel, iptri, ipqua, geoel, coord)
use constant
implicit none
integer*4,dimension(1:nvtri, 1:ntria),intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:ngeel, 1:nsize), intent(inout)::geoel
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
real*8::bt(ndegr), bq(ndegr)
!...local real number
real*8:: dxdr,dxds,dydr,dyds
real*8:: xc, yc, xcel, ycel
real*8:: eps,c00,c10,c05,c20
real*8:: r, s, djaco, volel, masel
real*8:: dr, ds, rc,sc
real*8:: xcerz, ycerz, maselrz, volelrz
real*8:: wi, xcent, ycent,xgaus, ygaus, rhog
real*8::radie, radii,radie2,radii2,radic2, radig2,sentr,rhoin,rhoex,rho0ba
real*8::preex,prein
real*8::rcoef !...Coefficient R of RZ or XY coordinate system...
real*8::areel !...2D area of the cell
real*8::xceli, yceli !...The initial cell center(volume)
real*8::bt22,bt33,bt23,bq22,bq33,bq23
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

!...Mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!
dr = 0.5d0
ds = 0.5d0

!...Cell mass
masel = geoel(4, ielem)

!...Zero out bt22,bt33,bt23
bt22 = 0.d0
bt33 = 0.d0
bt23 = 0.d0
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
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!volel = volel + djaco
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nptri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo
!
if(ncase.eq.1)then     !...TGV
rhog = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rhog = 1.d0
elseif(ncase.eq.3)then !...
rhog = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
!
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rhog = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rhog = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
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
!
!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then
!
rhog  = 1.d0
!
else
!
rhog  = 0.125d0
!
endif
!
elseif(ncase.eq.7)then !...sedov
!
rhog = 1.d0
!
elseif(ncase.eq.8)then !...Gresho
!
rhog = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
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
!
if(xc.lt.1.d0)then
rhog  = 1.d0
else
if(yc.gt.1.5d0)then
rhog  = 0.1d0
else
rhog = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rhog =1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rhog  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rhog = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rhog = 1.d0
!
elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rhog = 1.d0
!
endif
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
bt22 = bt22 + 0.5d0*rhog*bt(2)*bt(2)*djaco*rcoef
bt33 = bt33 + 0.5d0*rhog*bt(3)*bt(3)*djaco*rcoef
bt23 = bt23 + rhog*bt(2)*bt(3)*djaco*rcoef

!
if(nrz.eq.2)then
print*,'DG(P2) for AW-RZ will be implemented in the future!'
stop
endif
!
enddo
!
geoel(19, ielem) = bt22/masel
geoel(20, ielem) = bt33/masel
geoel(21, ielem) = bt23/masel
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

!...Mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!
dr = 1.d0
ds = 1.d0

!...Cell mass
masel = geoel(4, ielem)

!...Zero out bt22,bt33,bt23
bq22 = 0.d0
bq33 = 0.d0
bq23 = 0.d0
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
djaco = wi*(dxdr*dyds - dydr*dxds)
!volel = volel + djaco
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
if(ncase.eq.1)then     !...TGV
rhog = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rhog = 1.d0
elseif(ncase.eq.3)then !...
rhog = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
!
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rhog = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rhog = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
r=0.d0;s=0.d0
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo

!if(ielem.eq.1601)print*,'mass',rhog,sqrt(xgaus**2+ygaus**2)
!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then

rhog  = 1.d0
else
rhog  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
rhog = 1.d0
elseif(ncase.eq.8)then !...Gresho
!
rhog = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
r=0.d0;s=0.d0
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo
!
if(xc.lt.1.d0)then
rhog  = 1.d0
else
if(yc.gt.1.5d0)then
rhog  = 0.1d0
else
rhog = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rhog = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rhog  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rhog = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rhog = 1.d0

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rhog = 1.d0
!
endif
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
bq22 = bq22 + 0.5d0*rhog*bq(2)*bq(2)*djaco*rcoef
bq33 = bq33 + 0.5d0*rhog*bq(3)*bq(3)*djaco*rcoef
bq23 = bq23 + rhog*bq(2)*bq(3)*djaco*rcoef
!
!if(ie.eq.1)print*,'b5',ie,igaus,0.5d0*djaco*rhog*bq(3)*bq(3)/masel,0.5d0*bq(3)*bq(3)
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
geoel(19, ielem) = bq22/masel
geoel(20, ielem) = bq33/masel
geoel(21, ielem) = bq23/masel
!
if(nrz.eq.2)then
print*,'DG(P2) for AW-RZ will be implemented in the future!'
stop
endif
!
enddo !...(1)ie = 1,nquad
!
end subroutine getgeoel_lagmc_ho
!
!...Get the mass matrix for lagrangian based on mass center...
!
subroutine  getamatr_lag_mcnew(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
use constant
implicit none
!...Input
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize)::geoel
real*8,dimension(1:ndimn,1:npoin),intent(in)::coord
real*8,dimension(1:nmatr,1:ncell),intent(out)::amatr
integer,  dimension(1:nvtri,1:ntria), intent(in):: inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, ielem
!...Local real array
real*8::xp(1:2, 1:nptri)
real*8,dimension(1:2, 1:npqua)::xpq
real*8,dimension(1:nptri)::shp, dspr, dsps
real*8:: weight(ngausd), posit(2, ngausd)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
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
endif
!
!...Mass center...
!
rc= geoel(1, ielem)
sc= geoel(2, ielem)
!
dr = 0.5d0
ds = 0.5d0
!
masel = geoel(4, ielem)

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
do ishp = 1, nptri
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
if(ncase.eq.1)then     !...TGV
rho0 = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rho0 = 1.d0
elseif(ncase.eq.3)then !...
rho0 = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0)
rho0 = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rho0 = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
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
do ishp = 1, nptri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo

!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then

rho0  = 1.d0
else
rho0  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
rho0 = 1.d0
elseif(ncase.eq.8)then !...Gresho
rho0 = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
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
do ishp = 1, nptri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo
!
if(xc.lt.1.d0)then
rho0  = 1.d0
else
if(yc.gt.1.5d0)then
rho0  = 0.1d0
else
rho0 = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rho0 = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rho0  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rho0 = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rho0 = 1.d0

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rho0 = 1.d0
!
endif
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

if(npoly==1)then
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
!...Mass center...
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 1.d0
ds = 1.d0!
!
masel = geoel(4, ielem)

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
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
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
if(ncase.eq.1)then     !...TGV
rho0 = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rho0 = 1.d0
elseif(ncase.eq.3)then !...
rho0 = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rho0 = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rho0 = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
r=rc;s=sc
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo
!
!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then
rho0  = 1.d0
else
rho0  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
rho0 = 1.d0
elseif(ncase.eq.8)then !...Gresho
rho0 = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
r=rc;s=sc
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo
!
if(xc.lt.1.d0)then
rho0  = 1.d0
else
if(yc.gt.1.5d0)then
rho0  = 0.1d0
else
rho0 = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rho0 = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rho0  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rho0 = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rho0 = 1.d0

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rho0 = 1.d0
!
endif
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

if(npoly==1)then
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

end subroutine  getamatr_lag_mcnew
!
!...subroutine: Calculate the F^* N ds of all faces (mass center) only for triangle cell...
!
subroutine getfnds_lag_mc(gflag,gelag,geoel,intfac,inpoel,coord,lpnp)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(inout)::gelag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nelem+nbfac),     intent(in)::geoel
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri, 1:nelem),intent(out)::lpnp
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv
!...local integer array
integer,dimension(1:nvtri) :: ipl, ipr
!...local real array
real*8,dimension(1:2, 1:2)    ::comatr !...cofactor matrix...
real*8,dimension(1:ndimn, 1:nvtri)::coorp
real*8::b(3, nvtri)
real*8::vnorm(1:2)
real*8,dimension(1:nvtri):: xv, yv
!...local real number
real*8::dwav1,dwav2,larea,farea
real*8::anx, any
real*8::dr, ds, rc, sc
real*8::c16
!
data c16   /0.1666666666666666d0 /
!
do 100 ifa=1, nafac !...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
iel=intfac(1,ifa)
ier=intfac(2,ifa)
!
ipl(1:nvtri) = inpoel(1:nvtri,iel)
if(ifa.gt.nbfac) ipr(1:nvtri) = inpoel(1:nvtri,ier)
!
!...coordinates
!
coorp(1, 1:nvtri) = coord(1, ipl(1:nvtri))
coorp(2, 1:nvtri) = coord(2, ipl(1:nvtri))
!
!if(iel==9) print*,'coord ie==9',coorp(1:2,1),coorp(1:2,2),coorp(1:2,3),ipl(1)
!
!...Cofactor matrix for left cell
!
comatr(1, 1) = coorp(2,3) - coorp(2,1) !...yc-ya
comatr(1, 2) = coorp(2,1) - coorp(2,2) !...-(yb-ya)
comatr(2, 1) = coorp(1,1) - coorp(1,3) !...-(xc-xa)
comatr(2, 2) = coorp(1,2) - coorp(1,1) !...xb-xa
!
!if(iel==170) print*,'coord ie==170',comatr(:,:)
!
!...Identify the local No. of one internal face for left cell...
!
if(intfac(3,ifa)==ipl(1).and.intfac(4,ifa)==ipl(2))then
!
vnorm(1) = 0.d0;             vnorm(2) = -1.d0;
larea    = 1.d0;
!
idfal = 1;
!
elseif(intfac(3,ifa)==ipl(2).and.intfac(4,ifa)==ipl(3))then
!
vnorm(1) = sqrt(2.d0)*0.5d0; vnorm(2) = sqrt(2.d0)*0.5d0;
larea    = sqrt(2.d0)
!
idfal = 2;
elseif(intfac(3,ifa)==ipl(3).and.intfac(4,ifa)==ipl(1))then
!
vnorm(1) = -1.d0;           vnorm(2) = 0.d0;
larea    = 1.d0
!
idfal = 3;
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
gelag(1, idfal, iel) = anx*larea/farea
gelag(2, idfal, iel) = any*larea/farea
gelag(3, idfal, iel) = farea
!
!...Identify the local No. of one internal face for right cell...
if(ifa.gt.nbfac)then
if(intfac(3,ifa)==ipr(2).and.intfac(4,ifa)==ipr(1))then
!
idfar = 1;
!
elseif(intfac(3,ifa)==ipr(3).and.intfac(4,ifa)==ipr(2))then
!
idfar = 2;
elseif(intfac(3,ifa)==ipr(1).and.intfac(4,ifa)==ipr(3))then
!
idfar = 3;
endif
!
!...The outward unit vector of right cell is opposite to that of left cell...
gelag(1, idfar, ier) = - gelag(1, idfal, iel)
gelag(2, idfar, ier) = - gelag(2, idfal, iel)
gelag(3, idfar, ier) = farea
!
elseif(ifa.le.nbfac)then
!
gflag(1, ifa) = gelag(1, idfal, iel)
gflag(2, ifa) = gelag(2, idfal, iel)
gflag(3, ifa) = gelag(3, idfal, iel)
endif
!
100 enddo  !...(1)ifa=1,nafac
!
!  print*,'vnotmfn',gelag(1, 3, 9)
!
! print*,'Inside getfnds_lag'
!
!...Second part: Get {l_np n_np}
!
do 1000 ie = 1, nelem !...(1)ifa=1,nafac
!
! print*,'second part2',ie
!
!...shape functions
!
dr = .5d0
ds = .5d0
!rc = 1.d0/3.d0
!sc = rc
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0

do iv =1 ,nvtri
!
!print*,'iv',ivt
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!  print*,'iv',ivt
!
!...Get lpnp for different degrees of freedom...
!
!...For point(1) the first dof.
do ig = 1,ndegr
!...point 1
lpnp(1:ndimn, ig, 1, 1, ie) = c16*(2.d0*b(ig, 1) + b(ig, 3))*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
lpnp(1:ndimn, ig, 2, 1, ie) = c16*(2.d0*b(ig, 1) + b(ig, 2))*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
!
!...point 2
lpnp(1:ndimn, ig, 1, 2, ie) = c16*(2.d0*b(ig, 2) + b(ig, 1))*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
lpnp(1:ndimn, ig, 2, 2, ie) = c16*(2.d0*b(ig, 2) + b(ig, 3))*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
!
!...point 3
lpnp(1:ndimn, ig, 1, 3, ie) = c16*(2.d0*b(ig, 3) + b(ig, 2))*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
lpnp(1:ndimn, ig, 2, 3, ie) = c16*(2.d0*b(ig, 3) + b(ig, 1))*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
!
enddo
!
!if(ie==1) print*,'lnpn',lpnp(1:2,1,1,1,ie), lpnp(1:2,1,2,1,ie)!
!if(ie==1) print*,'lnpn',lpnp(1:2,1,1,2,ie), lpnp(1:2,1,2,2,ie)!
!if(ie==1) print*,'lnpn',lpnp(1:2,1,1,3,ie), lpnp(1:2,1,2,3,ie)!

1000 enddo

end subroutine getfnds_lag_mc
!
!...subroutine: Calculate the F^* N ds for hybrid linear grids...
!
subroutine getfnds_lag_mc_hybrid(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(inout)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad),      intent(inout)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,iel,ier,ie
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8,dimension(1:ndimn, 1:nvtri)::xp
real*8,dimension(1:ndimn, 1:nvqua)::xpq
!...local real number
real*8::anx, any
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
ipf(1:nvfac) = intfac(3:2+nvfac, ifa)
!
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))
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
!...Triangles...
!
do 50 ie=1, ntria !...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
!...coordinates
!
xp(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xp(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
!...Edge 1
!
anx =   xp(2, 2) - xp(2, 1)
any = -(xp(1, 2) - xp(1, 1))
!
gelag(1, 1 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 1 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 1 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 2
!
anx =   xp(2, 3) - xp(2, 2)
any = -(xp(1, 3) - xp(1, 2))
!
gelag(1, 2 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 2 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 2 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 3
!
anx =   xp(2, 1) - xp(2, 3)
any = -(xp(1, 1) - xp(1, 3))
!
gelag(1, 3 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 3 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 3 , ie) = sqrt(anx**2 + any**2)
!
50 enddo  !...(1)ifa=1,nelem
!  print*,'vnotmfn',gelag(1, 3, 9)
!
!...Quads...
!
do 150 ie=1, nquad !...(1)ifa=1,nafac
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
!...Edge 1
!
anx =   xpq(2, 2) - xpq(2, 1)
any = -(xpq(1, 2) - xpq(1, 1))
!
gelagq(1, 1 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 1 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 1 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 2
!
anx =   xpq(2, 3) - xpq(2, 2)
any = -(xpq(1, 3) - xpq(1, 2))
!
gelagq(1, 2 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 2 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 2 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 3
!
anx =   xpq(2, 4) - xpq(2, 3)
any = -(xpq(1, 4) - xpq(1, 3))
!
gelagq(1, 3 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 3 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 3 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 4
!
anx =   xpq(2, 1) - xpq(2, 4)
any = -(xpq(1, 1) - xpq(1, 4))
!
gelagq(1, 4 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 4 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 4 , ie) = sqrt(anx**2 + any**2)
!
150 enddo  !...(1)ifa=1,nelem

!
end subroutine getfnds_lag_mc_hybrid
!
!...subroutine: Calculate the nodal velocity U_p^* (mass center)...
!
subroutine getndvelo_lag_mc(gflag,gelag,geoel,bface,intfac,inpoel,coord,unkno,ustar, fstar, aflim, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3, 1:ngelg, 1:nelem+nbfac), intent(in)::gelag
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nelem+nbfac),     intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri, 1:nelem),  intent(out)::fstar !...Riemann forces
real*8,dimension(1:nq+1,1:nelem+nbfac),  intent(in)::aflim !...Limiter coef
integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2)     :: ipf
integer::indnd(npoin)

!...local real array
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:nvtri)::murie
real*8::vnorm(1:3, 1:2, 1:3)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8::aujmp(1:3, 1:nvtri)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:nvtri):: xv, yv
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::uvtx,vvtx,evtx, pvtx
real*8::rhom, rhomv
real*8::dux,duy,deltu
real*8::dwav1,dwav2
real*8::anx, any
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: cnsup(:), munacn(:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: usold(:,:), munacu(:,:), snsigm(:,:)
real*8,allocatable:: munacl(:,:,:)
real*8,allocatable:: snsigml(:,:,:,:), munaul(:,:,:,:)

!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:npoin), cnsup(1:npoin))
allocate (usold(1:ndimn, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munacl(1:2, 1:nvtri, 1:nelem), munaul(1:ndimn, 1:2, 1:nvtri,  1:nelem),&
snsigml(1:ndimn, 1:2,  1:nvtri,  1:nelem))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
!
!...Zero out vlave
!
cnsup = 0.d0
vlave = 0.d0
usold = ustar
indnd = 0
!
!
!...Mark the boundary nodes...
!
if(ncase.eq.2)then
do ifa = 1, nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
indnd(ipf(1:2)) = 1

enddo
endif
!
!
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
!
!      print*,'ifa', ifa,ipf
indnd(ipf(1:nvfac)) = 1
endif
enddo

!
do 200 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...shape functions
!
dr = .5d0
ds = .5d0
!rc = 1.d0/3.d0
!sc = rc
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0

do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...
!...zero out unknv
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
!
enddo
!
!...Impose limiter
!
if(nlimi.eq.1)then
unknv(2, iv) = unkno(1, 2, ie) + aflim(2, ie)*(unknv(2, iv) - unkno(1, 2, ie))
unknv(3, iv) = unkno(1, 3, ie) + aflim(3, ie)*(unknv(3, iv) - unkno(1, 3, ie))
endif
! if(ip(iv)==36) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
enddo
!
!...Accumulate nodal velocity...
vlave(1, ip(1:nvtri)) = vlave(1, ip(1:nvtri)) + unknv(2, 1:nvtri)
vlave(2, ip(1:nvtri)) = vlave(2, ip(1:nvtri)) + unknv(3, 1:nvtri)
!
!...Accumualte to get the no of cell surrounding one vertex...
!
cnsup(ip(1:nvtri)) = cnsup(ip(1:nvtri)) + 1.d0
!
200 enddo  !...(1)ie = 1,nelem
!
!...Get the averaged reconstructed nodal velocity...
!
do ipoin = 1, npoin
vlave(1:ndimn, ipoin) = vlave(1:ndimn, ipoin)/cnsup(ipoin)
enddo

!
!print*,'average vlave',vlave(1:2,15)
!
!
!...Zero out munacn
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
do 250 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)

!
!...shape functions
!
dr = .5d0
ds = .5d0
!rc = 1.d0/3.d0
!sc = rc
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
!...Basis function
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds

enddo
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
!
!...ndA=0.5d0*vnorm
!
vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
!...cell averaged value...
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...
!...zero out unknv
unknv = 0.d0
!
do iv   = 1,nvtri
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
!if(ip(iv)==5) print*,'adjumpunknv',unkno(2:3,2,ie)
!
rho  = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(nlimi.eq.1)then
rhomv = rhom + aflim(1, ie)*(unknv(1, iv) - rhom)
rho = 1.d0/rhomv
!
uvtx = uctr + aflim(2, ie)*(unknv(2, iv) - uctr)
vvtx = vctr + aflim(3, ie)*(unknv(3, iv) - vctr)
!
!pvtx = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
pvtx = pctr + aflim(4, ie)*(pvtx - pctr)
!
!...updtae unknv(2:3,:)
unknv(2, iv) = uvtx
unknv(3 ,iv) = vvtx
!
endif

!
!pvtx = 0.25d0*(cos(2.d0*pi*coord(1, ip(iv))) + cos(2.d0*pi*coord(2, ip(iv)))) + 1.d0
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
!aujmp(1:2, iv) = usold(1:2, ip(iv)) - unknv(2:3, iv)
aujmp(1:2, iv) = vlave(1:2, ip(iv)) - unknv(2:3, iv)
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
!if(ip(iv)==5) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the variables at the center...
!
rhoct = 1.d0/unkno(1, 1, ie)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ie)
vctr  = unkno(1, 3, ie)
ectr  = unkno(1, 4, ie)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...Get impedence coefficient...
!
do iv   = 1, nvtri
dux= vlave(1, ip(iv))-unknv(2, iv)
duy= vlave(2, ip(iv))-unknv(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ip(iv).eq.5) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ip(iv)),unknv(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
if(aujmp(3, iv)/sdctr.lt.1.d-9)then
!
!  print*,'itime',itime
munacn(ip(iv)) = munacn(ip(iv)) + murie(iv)*vnorm(3, ifa, iv)
!
munacu(1, ip(iv)) =  munacu(1, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknv(2, iv)
munacu(2, ip(iv)) =  munacu(2, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknv(3, iv)
else
munacn(ip(iv)) = munacn(ip(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
munacu(1, ip(iv)) =  munacu(1, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(2, iv)
munacu(2, ip(iv)) =  munacu(2, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(3, iv)
!
endif
!
!   if(ip(iv).eq.5) print*,'p5 muacn(vv) post',ie,murie(iv),munacn(ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p19 muacn(28) prep---',murie(iv), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                 vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
! if(ip(iv).eq.15) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Get the summed stress sum(n*A*sigma)
!
!
!  snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 1, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!
!  snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(1, 2, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!

snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) !
snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
if(aujmp(3, iv)/sdctr.lt.1.d-9)then

munacl(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)
munacl(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)

else

munacl(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)* &
abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))

munacl(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)* &
abs(vnorm(1, 2, iv)*aujmp(1, iv) + vnorm(2, 2, iv)*aujmp(2, iv))

endif
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaul(1, 1, iv, ie)    =  munacl(1, iv, ie)*unknv(2, iv)
munaul(2, 1, iv, ie)    =  munacl(1, iv, ie)*unknv(3, iv)

munaul(1:2, 2, iv, ie)    =  munacl(2, iv, ie)*unknv(2:3, iv)
!
!    snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!    snsigml(2, 1, iv, ie)= sigma(1, 2, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
!    snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!    snsigml(2, 2, iv, ie)= sigma(1, 2, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!
snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv)

snsigml(2, 1, iv, ie)= sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv)

snsigml(2, 2, iv, ie)= sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)!
enddo
!
250 enddo  !...(1)ie = 1,nelem
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!...Impose BC with pressure prescribed boundary...
!
call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)!
!...4.1: Update the Riemann forces at every node...
!
!fpres = 0.d0
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) - fpres(1, ipoin))/munacn(ipoin)
ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) - fpres(2, ipoin))/munacn(ipoin)
endif
enddo
!
!....Bd velocity
!
 !print*,'ustar--',ustar(1:2, 5),munacu(1:2,5) ,snsigm(1:2, 5), munacn(5)
!
!
if(ncase.eq.1)then
!
do 900 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
!    ustar(1, ipf(1)) = sin(pi*coord(1,ipf(1)))*cos(pi*coord(2,ipf(1)))
!    ustar(2, ipf(1)) =-cos(pi*coord(1,ipf(1)))*sin(pi*coord(2,ipf(1)))
!
!    ustar(1, ipf(2)) = sin(pi*coord(1,ipf(2)))*cos(pi*coord(2,ipf(2)))
!    ustar(2, ipf(2)) =-cos(pi*coord(1,ipf(2)))*sin(pi*coord(2,ipf(2)))
!
if(coord(1, ipf(1)).lt.1.d-6.or.abs(coord(1, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf',ipf(1)
ustar(1, ipf(1)) = 0.d0
endif
if(coord(2, ipf(1)).lt.1.d-6.or.abs(coord(2, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ipf(1)
ustar(2, ipf(1)) = 0.d0
endif
900 enddo
!
!  ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
!
!...Imposing the zero normal velocity for BC...
!
call getbcvn_lag(bface, intfac, gflag, ustar)
!
!call getbcve_exact(bface, intfac, gflag, ustar, coord, itime)
!
!   print*,'ustar',ustar(1:2, 18),ustar(1:2, 20),ustar(1:2, 21)
!  if(ustar(2, 18).gt.0.d0) then
!   print*,'Wrong movng direction'!,munacu(2, 28) - snsigm(2, 28) , fpres(2, 28),munacn(28)
!   stop
!  endif
!
!...4.2: Update the Riemann forces at every node...
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...shape functions
!
dr = .5d0
ds = .5d0
!rc = 1.d0/3.d0
!sc = rc
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
unknv = 0.d0
!
do iv = 1, nvtri
!
!...Basis function
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*ustar(1, ip(iv))- munaul(1, 1, iv, ie)
fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*ustar(2, ip(iv))- munaul(2, 1, iv, ie)
!
fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*ustar(1, ip(iv))- munaul(1, 2, iv, ie)
fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*ustar(2, ip(iv))- munaul(2, 2, iv, ie)
!
enddo
!
enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
end subroutine getndvelo_lag_mc
!
!...Face integral (mass center)...
!
subroutine rhsifacedg_lag_mc(inpoel,  unkno, ustar, fstar, lpnp, gelag, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri, 1:nelem),intent(in)::lpnp
real*8,dimension(1:ndimn,1:2,1:nvtri,1:ntria),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2, 1:nvtri) :: ipf
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri)::lpnpt
real*8::xv(3), yv(3),b(1:3,1:nvtri)
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
ipf(1, 1) = 3; ipf(2, 1) = 2
ipf(1, 2) = 1; ipf(2, 2) = 3
ipf(1, 3) = 2; ipf(2, 3) = 1

do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
!
!print*,'iv', ie
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpt(1:ndimn, ig, 1, 1) = c16*(2.d0*b(ig, 1) + b(ig, 3))*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
lpnpt(1:ndimn, ig, 2, 1) = c16*(2.d0*b(ig, 1) + b(ig, 2))*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
!
!...point 2
lpnpt(1:ndimn, ig, 1, 2) = c16*(2.d0*b(ig, 2) + b(ig, 1))*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
lpnpt(1:ndimn, ig, 2, 2) = c16*(2.d0*b(ig, 2) + b(ig, 3))*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
!
!...point 3
lpnpt(1:ndimn, ig, 1, 3) = c16*(2.d0*b(ig, 3) + b(ig, 2))*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
lpnpt(1:ndimn, ig, 2, 3) = c16*(2.d0*b(ig, 3) + b(ig, 1))*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
!
enddo

!...The vertex constituting one cell...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
!

!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvtri
!
!ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
!ustar(1, ip(iv))*lpnp(1, 1:ndegr, 1, iv, ie) +&
!ustar(2, ip(iv))*lpnp(2, 1:ndegr, 1, iv, ie) +&
!ustar(1, ip(iv))*lpnp(1, 1:ndegr, 2, iv, ie) +&
!ustar(2, ip(iv))*lpnp(2, 1:ndegr, 2, iv, ie)
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ip(iv))*lpnpt(1, 1:ndegr, 1, iv) +&
ustar(2, ip(iv))*lpnpt(2, 1:ndegr, 1, iv) +&
ustar(1, ip(iv))*lpnpt(1, 1:ndegr, 2, iv) +&
ustar(2, ip(iv))*lpnpt(2, 1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstar(1, 1, iv, ie)*c13*(2.d0*b(1:ndegr, iv) + b(1:ndegr, ipf(1,iv))) +&
fstar(1, 2, iv, ie)*c13*(2.d0*b(1:ndegr, iv) + b(1:ndegr, ipf(2,iv)))

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstar(2, 1, iv, ie)*c13*(2.d0*b(1:ndegr, iv) + b(1:ndegr, ipf(1,iv))) +&
fstar(2, 2, iv, ie)*c13*(2.d0*b(1:ndegr, iv) + b(1:ndegr, ipf(2,iv)))
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ip(iv))*fstar(1, 1, iv, ie)*c13*(2.d0*b(1:ndegr, iv) + b(1:ndegr, ipf(1,iv))) +&
ustar(2, ip(iv))*fstar(2, 1, iv, ie)*c13*(2.d0*b(1:ndegr, iv) + b(1:ndegr, ipf(1,iv))) +&
ustar(1, ip(iv))*fstar(1, 2, iv, ie)*c13*(2.d0*b(1:ndegr, iv) + b(1:ndegr, ipf(2,iv))) +&
ustar(2, ip(iv))*fstar(2, 2, iv, ie)*c13*(2.d0*b(1:ndegr, iv) + b(1:ndegr, ipf(2,iv)))
!
enddo
!
rhsel(1:ndegr, 1, ie) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ie) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ie) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ie) =  elnpn(1:ndegr)
!
! if(ie==18) print*,'rhs iface',rhsel(1, 1, ie), lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

550 enddo
!
!   open(8,file='lpnp.dat')
!    do ie = 1, nelem
!
!      do iv = 1, nvtri
!
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 1, 1, iv, ie), lpnp(1:2, 1, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 2, 1, iv, ie), lpnp(1:2, 2, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 3, 1, iv, ie), lpnp(1:2, 3, 2, iv, ie)
!      enddo
!    enddo
!   close(8)
!
!
!
! ip1 =8; ip2=37
!  print*,'dx', coord(1,ip1) ,coord(1, ip2), coord(2, ip1) -coord(2, ip2),&
!               0.5d0*(coord(2,ip1) -coord(2, ip2)),-0.5d0*(coord(1,ip1) -coord(1, ip2))
!/sqrt((coord(1,37) -coord(1, 2))**2+(coord(2,37) -coord(2, 2))**2)
end subroutine rhsifacedg_lag_mc
!
!...Domain integral (mass center)...
!
subroutine rhsdomndg_lag_mc(intfac, inpoel, coord, geoel, unkno, rhsel,aflim )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nelem+nbfac),     intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8, dimension(1:nq+1, 1:nelem+nbfac),      intent(in)::aflim
real*8,dimension(1:ndegr,1:nq,1:nelem),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv
!
!...local integer array
!
integer,dimension(1:nvtri) :: ip
!...local real array
real*8,dimension(1:ndimn, 1:nvtri) :: coorp
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc
real*8::xg, yg
real*8::rho,uadv,vadv,eadv,rhom
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr, rhomv
!
real*8,allocatable::weigh(:), posi(:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
allocate (weigh(ngausd), posi(2,ngausd))
call rutope(2, ngausd, posi, weigh)
!
!...Loop over elements
!
do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
!...Points consitituting one element...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
coorp(1, 1:nvtri) = coord(1, ip(1:nvtri))
coorp(2, 1:nvtri) = coord(2, ip(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
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
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
wi     = weigh(ig)
!
dspr(1) = -c10
dspr(2) =  c10
dspr(3) =  0.d0
!
dsps(1) = -c10
dsps(2) =  0.d0
dsps(3) =  c10
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*coorp(1,ishp)
dxds = dxds + dsps(ishp)*coorp(1,ishp)

dydr = dydr + dspr(ishp)*coorp(2,ishp)
dyds = dyds + dsps(ishp)*coorp(2,ishp)
enddo
!
djaco = 0.5d0*wi
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
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ie)*b(ideg)
enddo
!
!...Primitive variables...
!
rhom = unknod(1)
rho  = 1.d0/rhom
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
if(nlimi.eq.1)then
!
rhomc = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
rhomv = rhomc + aflim(1 ,ie)*(rhomv - rhomc)
rho = 1.d0/rhomv
!
uadv = uctr + aflim(2, ie)*(uadv - uctr)
vadv = vctr + aflim(3 ,ie)*(vadv - vctr)
!
pres = pctr + aflim(4, ie)*(pres- pctr)

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
rhsel(ideg,1:nq,ie)=rhsel(ideg,1:nq,ie) - fluxd(ideg,1:nq)*djaco
enddo
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomndg_lag_mc
!
!...Domain integral for source term (mass center)...
!
subroutine rhsdomnsrcdg_lag_mc(intfac, inpoel, coord, geoel,rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nelem+nbfac),     intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nelem),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv
!
!...local integer array
!
integer,dimension(1:nvtri) :: ip
!...local real array
real*8,dimension(1:ndimn, 1:nvtri) :: coorp
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc
real*8::xg, yg
real*8::djaco, wi,src
!
real*8,allocatable::weigh(:), posi(:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
allocate (weigh(ngausd), posi(2,ngausd))
call rutope(2, ngausd, posi, weigh)
!
!...Loop over elements
!
do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
!...Points consitituting one element...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
coorp(1, 1:nvtri) = coord(1, ip(1:nvtri))
coorp(2, 1:nvtri) = coord(2, ip(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
!rxx=
!rxy=
!ryy=
!
!...Gauss loop
!
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
wi     = weigh(ig)
!
dspr(1) = -c10
dspr(2) =  c10
dspr(3) =  0.d0
!
dsps(1) = -c10
dsps(2) =  0.d0
dsps(3) =  c10
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*coorp(1,ishp)
dxds = dxds + dsps(ishp)*coorp(1,ishp)

dydr = dydr + dspr(ishp)*coorp(2,ishp)
dyds = dyds + dsps(ishp)*coorp(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!
!...Gauss points...
!
xg = shp(1)*coorp(1, 1) + shp(2)*coorp(1, 2) + shp(3)*coorp(1, 3)
yg = shp(1)*coorp(2, 1) + shp(2)*coorp(2, 2) + shp(3)*coorp(2, 3)
!
!
!...Basis function for solutions...
!
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
!
!...source term
!
src = 0.25d0*pi*(cos(3.d0*pi*xg)*cos(pi*yg) - cos(3.d0*pi*yg)*cos(pi*xg))/(gamlg-1.d0)
!
!    src = 0.5d0*pi/(gamlg-1.d0)*(sin(2.d0*pi*yg)*cos(pi*xg)*sin(pi*yg) - sin(2.d0*pi*xg)*sin(pi*xg)*cos(pi*yg))
!
!finally, scatter the contribution to the RHS
!
! if(ie==2) print*,'src rhs', 3.5d0/9.d0,rhsel(1, 4, 2),coorp(1, 1:3), coorp(2, 1:3)
!
do ideg = 1,ndegr
rhsel(ideg,4,ie)=rhsel(ideg,4,ie) + src*b(ideg)*djaco
enddo
!
!if(ie==2) print*,'src rhs', rhsel(1, 4, 2),src,b(1),djaco,src*b(1)*djaco,xg,yg
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomnsrcdg_lag_mc
!
!...Subroutine: Set the new initial flow field (mass center)...
!
subroutine inifield_mcnew(unkno,uchar,geoel,coord,inpoel,iptri,ipqua)
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
real*8,dimension(1:2, 1:nptri)::xp
real*8:: shp(nptri), dspr(nptri), dsps(nptri)
real*8,dimension(1:2, 1:npqua)::xpq
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
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
endif
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
do ishp = 1, nptri
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
!...Shape functions for mass center...
!
shp(1) = -(c10-rc-sc)*(c10-2.d0*(c10-rc-sc))
shp(2) = -rc*(c10-2.d0*rc)
shp(3) = -sc*(c10-2.d0*sc)
shp(4) = 4.d0*rc*(c10-rc-sc)
shp(5) = 4.d0*rc*sc
shp(6) = 4.d0*sc*(c10-rc-sc)
!
!...xc, yc on physical domain corresponging to the mass center
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nptri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo
!!
dspr(1) = c10-4.d0*(c10-rc-sc)
dspr(2) = -1.d0 + 4.d0*rc
dspr(3) = 0.d0
dspr(4) = 4.d0*(1.d0-2.d0*rc-sc)
dspr(5) = 4.d0*sc
dspr(6) =-4.d0*sc

!
dsps(1) = c10-4.d0*(c10-rc-sc)
dsps(2) = 0.d0
dsps(3) = -1.d0 + 4.d0*sc
dsps(4) = -4.d0*rc
dsps(5) =  4.d0*rc
dsps(6) =  4.d0*(1.d0-rc-2.d0*sc)
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
unkno(1, 1, ielem) = rhomvel/masel;
unkno(1, 2, ielem) = rhouvel/masel;
unkno(1, 3, ielem) = rhovvel/masel;
unkno(1, 4, ielem) = rhoevel/masel
!
!if(ielem.eq.1) print*,'variable',ielem,rhomvel,masel
!
if(npoly.ge.1)then
!
call getunkct(unkct, xc, yc, xcel, ycel)
!
drhomdx = unkct(2, 1)!
drhomdy = unkct(3, 1)
du0dx   = unkct(2, 2)
du0dy   = unkct(3, 2)
dv0dx   = unkct(2, 3)
dv0dy   = unkct(3, 3)
de0dx   = unkct(2, 4)
de0dy   = unkct(3, 4)
!
unkno(2, 1, ielem) = drhomdx*dxdr + drhomdy*dydr
unkno(3, 1, ielem) = drhomdx*dxds + drhomdy*dyds
!
unkno(2, 2, ielem) = du0dx*dxdr + du0dy*dydr
unkno(3, 2, ielem) = du0dx*dxds + du0dy*dyds
!
unkno(2, 3, ielem) = dv0dx*dxdr + dv0dy*dydr
unkno(3, 3, ielem) = dv0dx*dxds + dv0dy*dyds
!
unkno(2, 4, ielem) = de0dx*dxdr + de0dy*dydr
unkno(3, 4, ielem) = de0dx*dxds + de0dy*dyds
!
!  if(ie==1)print*,'initial',du0dx, du0dy,uctr,vctr,xc,yc
!   if(ie==21.or.ie==5.or.ie==4)print*,'initial',ie,unkno(1, 3, ie),-cos(pi*xc)*sin(pi*yc),vvel,volel!+dv0dx*(0.8d0-xc) + dv0dy*(0.d0-yc)
!
!...Update the dierivative with drr and dsr...
unkno(2, 1:nq, ielem) = unkno(2, 1:nq, ielem)*dr
unkno(3, 1:nq, ielem) = unkno(3, 1:nq, ielem)*ds
!
!...DGP2
if(npoly.eq.2)then
drhomdx2 = unkct(4, 1)!
drhomdy2 = unkct(5, 1)
drhomdxy = unkct(6, 1)
!
du0dx2   = unkct(4, 2)
du0dy2   = unkct(5, 2)
du0dxy   = unkct(6, 2)
!
dv0dx2   = unkct(4, 3)
dv0dy2   = unkct(5, 3)
dv0dxy   = unkct(6, 3)
!
de0dx2   = unkct(4, 4)
de0dy2   = unkct(5, 4)
de0dxy   = unkct(6, 4)
!
unkno(4, 1, ielem) = drhomdx2*dxdr**2   + 2.d0*drhomdxy*dxdr*dydr +drhomdy2*dydr**2
unkno(5, 1, ielem) = drhomdx2*dxds**2   + 2.d0*drhomdxy*dxds*dyds +drhomdy2*dyds**2
unkno(6, 1, ielem) = drhomdx2*dxdr*dxds + drhomdxy*(dxdr*dyds + dxds*dydr) +drhomdy2*dydr*dyds
!
unkno(4, 2, ielem) = du0dx2*dxdr**2 + 2.d0*du0dxy*dxdr*dydr +du0dy2*dydr**2
unkno(5, 2, ielem) = du0dx2*dxds**2 + 2.d0*du0dxy*dxds*dyds +du0dy2*dyds**2
unkno(6, 2, ielem) = du0dx2*dxdr*dxds   + du0dxy*(dxdr*dyds + dxds*dydr) +du0dy2*dydr*dyds
!
unkno(4, 3, ielem) = dv0dx2*dxdr**2 + 2.d0*dv0dxy*dxdr*dydr +dv0dy2*dydr**2
unkno(5, 3, ielem) = dv0dx2*dxds**2 + 2.d0*dv0dxy*dxds*dyds +dv0dy2*dyds**2
unkno(6, 3, ielem) = dv0dx2*dxdr*dxds + dv0dxy*(dxdr*dyds + dxds*dydr) +dv0dy2*dydr*dyds
!
unkno(4, 4, ielem) = de0dx2*dxdr**2 + 2.d0*de0dxy*dxdr*dydr +de0dy2*dydr**2
unkno(5, 4, ielem) = de0dx2*dxds**2 + 2.d0*de0dxy*dxds*dyds +de0dy2*dyds**2
unkno(6, 4, ielem) = de0dx2*dxdr*dxds + de0dxy*(dxdr*dyds + dxds*dydr) +de0dy2*dydr*dyds

!...Update the dierivative with drr and dsr...
unkno(4, 1:nq, ielem) = unkno(4, 1:nq, ielem)*dr**2
unkno(5, 1:nq, ielem) = unkno(5, 1:nq, ielem)*ds**2
unkno(6, 1:nq, ielem) = unkno(6, 1:nq, ielem)*dr*ds
!
endif

endif !...npoly.ge.1
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
!...Reference coordinates
!
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
!
!...Avoid the centroid outof one curved cell...
!
xcel = xpq(1,9)
ycel = xpq(2,9)
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
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
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
!...rc, sc
!
r= rc
s= sc
!
!...Shape functions for mass center...
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
!...xc, yc on physical domain corresponging to the mass center
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo
!!
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
unkno(1, 1, ielem) = rhomvel/masel;
unkno(1, 2, ielem) = rhouvel/masel;
unkno(1, 3, ielem) = rhovvel/masel;
unkno(1, 4, ielem) = rhoevel/masel

if(nmatel.eq.1)then
 if(ncase.eq.3.and.ielem.gt.1)then
 ! unkno(1, 2, ielem) =-xc/sqrt(xc**2 + yc**2);
 ! unkno(1, 3, ielem) =-yc/sqrt(xc**2 + yc**2);
 endif
endif
!
!...Record the local basis
!
!geoel(5:6, ielem) = unkno(1, 2:3, ielem)
!if(ielem.eq.100) print*,'variable',ielem,unkno(1, 1:4, ielem),masel,0.5d0*(unkno(1, 2, ielem)**2+unkno(1, 3, ielem)**2),&
!(gamlg-1.d0)*(unkno(1, 4, ielem)-0.5d0*(unkno(1, 2, ielem)**2+unkno(1, 3, ielem)**2))
!
if(npoly.ge.1)then

call getunkct(unkct, xc, yc, xcel, ycel)
!
drhomdx = unkct(2, 1)!
drhomdy = unkct(3, 1)
du0dx   = unkct(2, 2)
du0dy   = unkct(3, 2)
dv0dx   = unkct(2, 3)
dv0dy   = unkct(3, 3)
de0dx   = unkct(2, 4)
de0dy   = unkct(3, 4)
!
unkno(2, 1, ielem) = drhomdx*dxdr + drhomdy*dydr
unkno(3, 1, ielem) = drhomdx*dxds + drhomdy*dyds
!
unkno(2, 2, ielem) = du0dx*dxdr + du0dy*dydr
unkno(3, 2, ielem) = du0dx*dxds + du0dy*dyds
!
unkno(2, 3, ielem) = dv0dx*dxdr + dv0dy*dydr
unkno(3, 3, ielem) = dv0dx*dxds + dv0dy*dyds
!
unkno(2, 4, ielem) = de0dx*dxdr + de0dy*dydr
unkno(3, 4, ielem) = de0dx*dxds + de0dy*dyds
!
!if(ie==1)print*,'initial',du0dx, du0dy,dxdr,dydr,dxds,dyds,unkno(1:3, 2, ielem),rc,sc,geoel(1:2,ielem)
!   if(ie==21.or.ie==5.or.ie==4)print*,'initial',ie,unkno(1, 3, ie),-cos(pi*xc)*sin(pi*yc),vvel,volel!+dv0dx*(0.8d0-xc) + dv0dy*(0.d0-yc)
!
!...Update the dierivative with drr and dsr...
!
unkno(2, 1:nq, ielem) = unkno(2, 1:nq, ielem)*dr
unkno(3, 1:nq, ielem) = unkno(3, 1:nq, ielem)*ds

!
!...DGP2
if(npoly.eq.2)then
drhomdx2 = unkct(4, 1)!
drhomdy2 = unkct(5, 1)
drhomdxy = unkct(6, 1)
!
du0dx2   = unkct(4, 2)
du0dy2   = unkct(5, 2)
du0dxy   = unkct(6, 2)
!
dv0dx2   = unkct(4, 3)
dv0dy2   = unkct(5, 3)
dv0dxy   = unkct(6, 3)
!
de0dx2   = unkct(4, 4)
de0dy2   = unkct(5, 4)
de0dxy   = unkct(6, 4)
!
unkno(4, 1, ielem) = drhomdx2*dxdr**2   + 2.d0*drhomdxy*dxdr*dydr +drhomdy2*dydr**2
unkno(5, 1, ielem) = drhomdx2*dxds**2   + 2.d0*drhomdxy*dxds*dyds +drhomdy2*dyds**2
unkno(6, 1, ielem) = drhomdx2*dxdr*dxds + drhomdxy*(dxdr*dyds + dxds*dydr) +drhomdy2*dydr*dyds
!
unkno(4, 2, ielem) = du0dx2*dxdr**2 + 2.d0*du0dxy*dxdr*dydr +du0dy2*dydr**2
unkno(5, 2, ielem) = du0dx2*dxds**2 + 2.d0*du0dxy*dxds*dyds +du0dy2*dyds**2
unkno(6, 2, ielem) = du0dx2*dxdr*dxds   + du0dxy*(dxdr*dyds + dxds*dydr) +du0dy2*dydr*dyds
!
unkno(4, 3, ielem) = dv0dx2*dxdr**2 + 2.d0*dv0dxy*dxdr*dydr +dv0dy2*dydr**2
unkno(5, 3, ielem) = dv0dx2*dxds**2 + 2.d0*dv0dxy*dxds*dyds +dv0dy2*dyds**2
unkno(6, 3, ielem) = dv0dx2*dxdr*dxds + dv0dxy*(dxdr*dyds + dxds*dydr) +dv0dy2*dydr*dyds
!
unkno(4, 4, ielem) = de0dx2*dxdr**2 + 2.d0*de0dxy*dxdr*dydr +de0dy2*dydr**2
unkno(5, 4, ielem) = de0dx2*dxds**2 + 2.d0*de0dxy*dxds*dyds +de0dy2*dyds**2
unkno(6, 4, ielem) = de0dx2*dxdr*dxds + de0dxy*(dxdr*dyds + dxds*dydr) +de0dy2*dydr*dyds

!...Update the dierivative with drr and dsr...
unkno(4, 1:nq, ielem) = unkno(4, 1:nq, ielem)*dr**2
unkno(5, 1:nq, ielem) = unkno(5, 1:nq, ielem)*ds**2
unkno(6, 1:nq, ielem) = unkno(6, 1:nq, ielem)*dr*ds
!
endif

endif !npoly.ge.1

!...Degenerated P2 initial condition
!unkno(2:6, 1:nq, ielem) = 0.d0
if(ielem.eq.1)print*,'Initial',unkno(:,1, 1770)!,unkct(:, 3)

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
end subroutine inifield_mcnew
!
!...Get the initial condition...
!
subroutine  getunknot0_lag(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
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
real*8::xp(1:2, 1:nptri)
real*8,dimension(1:2, 1:npqua)::xpq
real*8,dimension(1:nptri)::shp, dspr, dsps
real*8:: weight(ngausd), posit(2, ngausd)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
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
endif

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
do ishp = 1, nptri
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
do ishp = 1, nptri
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
do ishp = 1, nptri
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
!
if(ncase.eq.1)then     !...TGV
rho0 = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rho0 = 1.d0
elseif(ncase.eq.3)then !...
rho0 = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0)
rho0 = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rho0 = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
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
do ishp = 1, nptri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo

!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then

rho0  = 1.d0
else
rho0  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
rho0 = 1.d0
elseif(ncase.eq.8)then !...Gresho
rho0 = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
if(xc.lt.1.d0)then
rho0  = 1.d0
else
if(yc.gt.1.5d0)then
rho0  = 0.1d0
else
rho0 = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rho0 = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rho0  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rho0 = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rho0 = 1.d0

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rho0 = 1.d0
!
endif

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
!...Mass center...
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 1.d0
ds = 1.d0!

!..Physical coordinate
r=rc;s=sc
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
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
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
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

!
if(ncase.eq.1)then     !...TGV
rho0 = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rho0 = 1.d0
elseif(ncase.eq.3)then !...
rho0 = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rho0 = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rho0 = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then
rho0  = 1.d0
else
rho0  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
rho0 = 1.d0
elseif(ncase.eq.8)then !...Gresho
rho0 = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
if(xc.lt.1.d0)then
rho0  = 1.d0
else
if(yc.gt.1.5d0)then
rho0  = 0.1d0
else
rho0 = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rho0 = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rho0  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rho0 = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rho0 = 1.d0

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rho0 = 1.d0
!
endif

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

end subroutine  getunknot0_lag
!
!...Get the initial condition...
!
subroutine  getunkini_lag(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
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
integer :: ie, ig, ishp, ielem, iq, ideg
!...Local real array
real*8::xp(1:2, 1:nptri)
real*8,dimension(1:2, 1:npqua)::xpq
real*8,dimension(1:nptri)::shp, dspr, dsps
real*8:: weight(ngausd), posit(2, ngausd)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:nq)::unkin
real*8,dimension(1:ndegr, 1:nq)::rhsini
real*8,dimension(1:ndegr)::bt, bq

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
endif
!
!...Mass center...
!
rc= geoel(1, ielem)
sc= geoel(2, ielem)
!
dr = 0.5d0
ds = 0.5d0
!
volel = geoel(3, ielem)
masel = geoel(4, ielem)

!...Physical coordinate
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
do ishp = 1, nptri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo

!...Zero out the initial rhs
rhsini = 0.d0

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
do ishp = 1, nptri
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
if(ncase.eq.1)then     !...TGV
rho0 = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rho0 = 1.d0
elseif(ncase.eq.3)then !...
rho0 = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0)
rho0 = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rho0 = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
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
do ishp = 1, nptri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo

!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then

rho0  = 1.d0
else
rho0  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
rho0 = 1.d0
elseif(ncase.eq.8)then !...Gresho
rho0 = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
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
do ishp = 1, nptri
xc = xc + shp(ishp)*xp(1,ishp)
yc = yc + shp(ishp)*xp(2,ishp)
enddo
!
if(xc.lt.1.d0)then
rho0  = 1.d0
else
if(yc.gt.1.5d0)then
rho0  = 0.1d0
else
rho0 = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rho0 = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rho0  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rho0 = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rho0 = 1.d0

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rho0 = 1.d0
!
endif

!...Get the density weighted variables
call getunkin(unkin, xgaus, ygaus, ielem, volel,xc,yc)

!...Coefficient R of RZ or XY system...
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
if(nrz.eq.0.or.nrz.eq.1)then
f0 = f0 + rho0*djaco*rcoef
f1 = f1 + rho0*(xg-rc)/dr*(xg-rc)/dr*djaco*rcoef
f2 = f2 + rho0*(xg-rc)/dr*(yg-sc)/ds*djaco*rcoef
f3 = f3 + rho0*(yg-sc)/ds*(yg-sc)/ds*djaco*rcoef
!
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
!...Basis function
bt(1) = 1.d0
bt(2) = (xg-rc)/dr
bt(3) = (yg-sc)/ds
bt(4) = 0.5d0*bt(2)*bt(2) - geoel(19, ielem)
bt(5) = 0.5d0*bt(3)*bt(3) - geoel(20, ielem)
bt(6) =       bt(2)*bt(3) - geoel(21, ielem)
!
do iq = 1, nq
do ideg = 1, ndegr
rhsini(ideg, iq) = rhsini(ideg, iq) + unkin(iq)*bt(ideg)*djaco
enddo
enddo
!
!if(ie.eq.1) print*,'ielem', ie, f22,f23,f24,f25,f26

if(npoly==1)then
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
unkno(2:6, :, ielem) = 0.d0
!...Get the initial unkno
do iq = 1, nq
do ideg = 1, ndegr-1
unkno(2, iq, ielem) = unkno(2, iq, ielem) + x5(1, ideg)*rhsini(ideg+1, iq)
unkno(3, iq, ielem) = unkno(3, iq, ielem) + x5(2, ideg)*rhsini(ideg+1,iq)
unkno(4, iq, ielem) = unkno(4, iq, ielem) + x5(3, ideg)*rhsini(ideg+1,iq)
unkno(5, iq, ielem) = unkno(5, iq, ielem) + x5(4, ideg)*rhsini(ideg+1,iq)
unkno(6, iq, ielem) = unkno(6, iq, ielem) + x5(5, ideg)*rhsini(ideg+1,iq)
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
!...Mass center...
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 1.d0
ds = 1.d0!

!..Physical coordinate
r=rc;s=sc
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo

!
volel = geoel(3, ielem)
masel = geoel(4, ielem)

!...Zero out the initial rhs
rhsini = 0.d0

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
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
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
if(ncase.eq.1)then     !...TGV
rho0 = 1.d0
elseif(ncase.eq.2)then !...Shockless Noh problem...
rho0 = 1.d0
elseif(ncase.eq.3)then !...
rho0 = 1.d0
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radig2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rho0 = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rho0 = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...
!
r=rc;s=sc
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo
!
!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then
rho0  = 1.d0
else
rho0  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
rho0 = 1.d0
elseif(ncase.eq.8)then !...Gresho
rho0 = 1.d0
!
elseif(ncase.eq.9)then !...Triple point...
!
r=rc;s=sc
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
xc = 0.d0
yc = 0.d0
!
do ishp = 1, npqua
xc = xc + shpq(ishp)*xpq(1,ishp)
yc = yc + shpq(ishp)*xpq(2,ishp)
enddo
!
if(xc.lt.1.d0)then
rho0  = 1.d0
else
if(yc.gt.1.5d0)then
rho0  = 0.1d0
else
rho0 = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion...!
!
rho0 = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rho0  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rho0 = 1.d0 + 0.9999995d0*sin(pi*xgaus)

elseif(ncase.eq.13)then !...Saltzman
rho0 = 1.d0

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rho0 = 1.d0
!
endif

!...Get the density weighted variables
call getunkin(unkin, xgaus, ygaus, ielem, volel,xc,yc)

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
!...Basis function
bq(1) = 1.d0
bq(2) = (xg-rc)/dr
bq(3) = (yg-sc)/ds
bq(4) = 0.5d0*bq(2)*bq(2) - geoel(19, ielem)
bq(5) = 0.5d0*bq(3)*bq(3) - geoel(20, ielem)
bq(6) =       bq(2)*bq(3) - geoel(21, ielem)
!
do iq = 1, nq
do ideg = 1, ndegr
rhsini(ideg, iq) = rhsini(ideg, iq) + unkin(iq)*bq(ideg)*djaco
!
!if(ie.eq.1.and.iq.eq.1.and.ideg.eq.5) print*,'ielem', ie,iq, ideg,ig, xg,yg,djaco,0.5d0*bq(3)*bq(3),rhsini(ideg, iq)&
!,geoel(20, ielem)
enddo
enddo

enddo


if(npoly==1)then
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
if(ielem.eq.1)print*,'amatr',f22,f24,f44
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
!
!...Treatment for AW RZ
if(nrz.eq.2) amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)
!
!
if(ielem.eq.1)print*,rhsini(:,1)
!
unkno(2:6, :, ielem) = 0.d0
!...Get the initial unkno
do iq = 1, nq
do ideg = 1, ndegr-1
unkno(2, iq, ielem) = unkno(2, iq, ielem) + x5(1, ideg)*rhsini(ideg+1,iq)
unkno(3, iq, ielem) = unkno(3, iq, ielem) + x5(2, ideg)*rhsini(ideg+1,iq)
unkno(4, iq, ielem) = unkno(4, iq, ielem) + x5(3, ideg)*rhsini(ideg+1,iq)
unkno(5, iq, ielem) = unkno(5, iq, ielem) + x5(4, ideg)*rhsini(ideg+1,iq)
unkno(6, iq, ielem) = unkno(6, iq, ielem) + x5(5, ideg)*rhsini(ideg+1,iq)
enddo
enddo
!
enddo !...(2)ie = 1,nelem

end subroutine  getunkini_lag

!
!...Get initial unkno distribution at gauss points...
!
subroutine getunkin(unkin, xg, yg, ielem, volel,xc,yc)
use constant
implicit none
real*8, dimension(1:nq), intent(out)::unkin
real*8,  intent(in):: xg,yg, volel,xc,yc
integer, intent(in) :: ielem
real*8:: pini, uini, vini, rhoini, ieini
integer::nind
real*8::radi, radiv, radic,ratir,theta,vemag
real*8::hgre
real*8::dhgre(2),dhgr2(3)
real*8::radie, radii,radie2,radii2,radic2, radig2,sentr,rhoin,rhoex,rho0ba
real*8::prein,preex,pini1,pini2
real*8::acap,omega,omegc,s1cap,c1cap,vini1,vini2
!
!...Different materials
!
select case (nmatel)

!...Ideal gas
case (1)

if(ncase.eq.1)then     !...TGV
rhoini = 1.d0
uini   = sin(pi*xg)*cos(pi*yg)
vini   =-cos(pi*xg)*sin(pi*yg)
pini   = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 1.d0

elseif(ncase.eq.2)then !...Shockless Noh problem...
rhoini = 1.d0
uini = -xg
vini = -yg
pini = (gamlg-1.d0)

elseif(ncase.eq.3)then !...Noh implosion
rhoini  = 1.d0
uini =  -xg/sqrt(xg**2 + yg**2)
vini =  -yg/sqrt(xg**2 + yg**2)
pini  = 1.d-6
!
!pini = (gamlg-1.d0)*rhoini*pini

elseif(ncase.eq.4)then !... Kidder shell...
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2

!rhoin = 6.31d-4
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!print*,'rhoin',rhoin,sentr
!
radie2 = radie**2
radii2 = radii**2
radig2 = xg**2 + yg**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...

rhoini = rho0ba**(1.d0/(gamlg-1.d0))
uini  = 0.d0
vini  = 0.d0
pini  = sentr*(rhoini)**gamlg

elseif(ncase.eq.5)then !...Kidder ball...
rhoini  = 2.d0*exp(-(xg**2+yg**2))
uini = 0.d0
vini = 0.d0
ieini = 0.5d0
pini  = ieini*(gamlg-1.d0)*rhoini

elseif(ncase.eq.6)then !...Sod...
!if(xc.lt.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then
rhoini  = 1.d0
uini = 0.d0
vini = 0.d0
ieini = 2.5d0
pini  = ieini*(gamlg-1.d0)*rhoini
else
rhoini  = 0.125d0
uini = 0.d0
vini = 0.d0
ieini = 2.0d0
pini  = ieini*(gamlg-1.d0)*rhoini
endif

elseif(ncase.eq.7)then !...Sedov
!...Different initial conditions

!...1st: The whole domain(Box)
!...1.1:  60x60 (ielem.eq.1770.or.ielem.eq.1771.or.ielem.eq.1830.or.ielem.eq.1831)
!...1.2:120x120 (ielem.eq.7140.or.ielem.eq.7141.or.ielem.eq.7260.or.ielem.eq.7261)

!...2nd: One quadrant domain(Box)
!...2.1: ielem .eq. 1
!...XY: pini = (gamlg-1.d0)*1.d0*0.244816d0/(1.d0*volel) 
!...RZ: pini = (gamlg-1.d0)*1.d0*0.425536d0/(2.d0*pi)/(volel)

!...3rd: One domain(Polar)
!...2.1: ielem.ge.1.and.ielem.le.20
!...Quad   TV (pi/4-pi/2): pini = (gamlg-1.d0)*1.d0*0.425536d0/(2.d0*pi)/(1.8637895114842526E-005)
!...Hybrid TV (pi/4-pi/2): pini = (gamlg-1.d0)*1.d0*0.425536d0/(2.d0*pi)/(5.2998940957972727E-007)
!...Hybrid AW (0  - pi/2): pini = (gamlg-1.d0)*1.d0*0.425536d0/(2.d0*pi)/(5.7511219211513981E-007)

 if(ielem.eq.1770.or.ielem.eq.1771.or.ielem.eq.1830.or.ielem.eq.1831)then
  rhoini  = 1.d0
  uini = 0.d0
  vini = 0.d0
  ieini=658.199744d0
  pini = (gamlg-1.d0)*1.d0*0.244816d0/(1.d0*volel)
!
else
  rhoini  = 1.d0
  uini = 0.d0
  vini = 0.d0
  pini = 1.d-6
 endif
elseif(ncase.eq.16)then     !...Gresho vortex
!
nind = 6
radiv= 0.4d0
radic= sqrt(xc**2+yc**2)
radi = sqrt(xg**2+yg**2)
!
rhoini = 1.d0
!
call facos(xg, yg, theta)
!
ratir = radi/radiv
!
!...Initila velocity
!
if(radic.lt.radiv)then
vemag = 2**(2*nind)*ratir**nind*(1.d0-ratir)**nind
uini = -vemag*sin(theta)
vini =  vemag*cos(theta)
else
uini = 0.d0
vini = 0.d0
endif
!
!..Initial pressure...
!
call greshoh(xg, yg ,radic, hgre, dhgre, dhgr2)
!
pini = 5.d0 + 2**(4*nind)*hgre
!
elseif(ncase.eq.8)then     !...Gresho vortex from LLNL
!
radiv= 0.4d0
radic= sqrt(xc**2+yc**2)
radi = sqrt(xg**2+yg**2)
!
rhoini = 1.d0

!...Initila velocity
if(radic.lt.0.2d0)then
uini =  5.d0*yg
vini = -5.d0*xg
elseif(radic.lt.radiv.and.radic.gt.0.2d0)then
uini = 2.d0*yg/radi - 5.d0*yg
vini =-2.d0*xg/radi + 5.d0*xg
else
uini = 0.d0
vini = 0.d0
endif

!..Initial pressure...

if(radic.lt.0.2d0)then
pini =  5.d0 + 25.d0/2.d0*radi**2
elseif(radic.lt.radiv.and.radic.gt.0.2d0)then
pini1 = 9.d0 - 4.d0*log(0.2d0) + 25.d0/2.d0*radi**2
pini2 = 20.d0*radi - 4.d0*log(radi)
pini = (pini1-pini2)
else
pini = 3.d0 + 4.d0*log(2.d0)
endif
!
elseif(ncase.eq.9)then !...Triple point...
!
if(xc.lt.1.d0)then
rhoini  = 1.d0
uini = 0.d0
vini = 0.d0
ieini= 2.5d0
pini  = ieini*(gamlg-1.d0)*rhoini
else
if(yc.gt.1.5d0)then
rhoini  = 0.1d0
uini = 0.d0
vini = 0.d0
ieini= 2.5d0
pini  = ieini*(gamlg-1.d0)*rhoini
else
rhoini = 1.d0
uini = 0.d0
vini = 0.d0
ieini= .25d0
pini  = ieini*(gamlg-1.d0)*rhoini
endif
endif
!
elseif(ncase.eq.10)then     !...Expansion
!
rhoini = 1.d0
uini   = 0.d0
vini   = 0.d0
pini   = -(xg**2 + yg**2) + 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rhoini= 1.d0
uini =  0.d0
vini =  0.d0
ieini = 1.d-5
pini  = ieini*(gamlg-1.d0)*rhoini!
!
elseif(ncase.eq.12)then !...1D isentropic expansion...
rhoini = 1.d0 + 0.9999995d0*sin(pi*xg)
uini   = 0.d0
vini   = 0.d0
pini   = rhoini**gamlg

elseif(ncase.eq.13)then !...Saltzman
rhoini  = 1.d0
uini = 0.d0
vini = 0.d0
ieini= 1.d-6
pini  = ieini*(gamlg-1.d0)*rhoini

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rhoini  = 1.d0
uini = -xg/4.d0
vini = -yg
ieini= (3.d0*xg/8.d0)**2
pini  = ieini*(gamlg-1.d0)*rhoini
!
elseif(ncase.eq.15)then !...Shu-Osher...
!
if(xc.lt.-4.d0)then
rhoini  = 3.85714d0
uini = 2.629369d0
vini = 0.d0
pini  = 10.333333d0
else
rhoini  = 1.d0+0.2d0*sin(5.d0*xg)
uini    = 0.d0
vini    = 0.d0
pini    = 1.d0
endif
!
endif

!...Solid
case (2)
!
 if(ncase.eq.1)then !...Bending beam
 rhoini  = 2.79d0
 uini = 0.d0
 vini = 10.d0
 pini = 1d-14
 elseif(ncase.eq.2)then !...Bending beam beryllium
!
acap  = 4.3368504248678715d-003/1.d0
omega = 0.23597445253266536d0
omegc= 0.788340124143784d0
s1cap = 57.645520483362446d0
c1cap = 56.636851535303499d0
!
 rhoini  = 1.845d0
 uini = 0.d0
 vini1 = c1cap*(sinh(omegc*(xg+3.d0))+sin(omegc*(xg+3.d0)))
 vini2 = s1cap*(cosh(omegc*(xg+3.d0))+cos(omegc*(xg+3.d0)))
 vini = acap*omega*(vini1-vini2)
 pini = 1d-14

 elseif(ncase.eq.3)then !...Tayor anvil
 rhoini  = 2.785d0
 uini =-0.015d0
 vini = 0.d0
 pini = 1d-14

elseif(ncase.eq.104)then !...Be Shell
!
 radie = 10.d0
 radii =  8.d0
 radig2 = xg**2 + yg**2
 vemag = -0.04171d0*radii/sqrt(radig2)
!
 rhoini  = 1.845d0
 uini = vemag*xg/sqrt(radig2)
 vini = vemag*yg/sqrt(radig2)
 pini = 1d-14
 endif
!
end select
!
unkin(1) = 1.d0 !
unkin(2) = rhoini*uini
unkin(3) = rhoini*vini

if(nmatel.eq.1)then
unkin(4) = (pini/rhoini/(gamlg-1.d0) + 0.5d0*(uini**2 + vini**2))*rhoini
elseif(nmatel.eq.2)then
unkin(4) = (pini/rhoini/(gamlg) + 0.5d0*(uini**2 + vini**2))*rhoini
endif
!
!if(ielem==13)print*,'ie',unkin(1:4),(pini/rhoini/(gamlg-1.d0) + 0.5d0*(uini**2 + vini**2))*rhoini
!
end subroutine getunkin
!
!...Function of acos
!
subroutine facos(xg, yg ,theta)
use constant
implicit none
real*8, intent(in):: xg, yg
real*8, intent(out):: theta
!
theta = 0.d0
!
if(xg.gt.0.d0)then
theta = atan(yg/xg)
elseif(xg.lt.0.d0.and.yg.ge.0.d0)then
theta = atan(yg/xg)+ pi
elseif(xg.lt.0.d0.and.yg.lt.0.d0)then
theta = atan(yg/xg)- pi
elseif(xg.eq.0.d0.and.yg.gt.0.d0)then
theta = 0.5d0*pi
elseif(xg.eq.0.d0.and.yg.lt.0.d0)then
theta =-0.5d0*pi
endif

end subroutine facos
!
subroutine greshoh(xg, yg ,radic,hgre, dhgre, dhgr2)
implicit none
real*8, intent(in):: xg, yg, radic
real*8, intent(out):: hgre
real*8, intent(out):: dhgre(1:2), dhgr2(3)
real*8::gcoef(13)
real*8::drdx, drdy,drdx2,drdy2,drdxy, radi,radiv, ratir
integer::ih
!
gcoef( 1) =   1.d0/12.d0
gcoef( 2) = -12.d0/13.d0
gcoef( 3) =  33.d0/7.d0
gcoef( 4) = -44.d0/3.d0
gcoef( 5) = 495.d0/16.d0
gcoef( 6) =-792.d0/17.d0
gcoef( 7) = 154.d0/3.d0
gcoef( 8) =-792.d0/19.d0
gcoef( 9) =  99.d0/4.d0
gcoef(10) =-220.d0/21.d0
gcoef(11) = 3.d0
gcoef(12) =-12.d0/23.d0
gcoef(13) =  1.d0/24.d0
!
!
radiv= 0.4d0
radi = sqrt(xg**2+yg**2)
ratir = radi/radiv
!
drdx = xg/radi
drdy = yg/radi
!
drdx2 = 1.d0/radi - xg/radi**2*drdx
drdy2 = 1.d0/radi - yg/radi**2*drdy
drdxy= -xg/radi**2*drdy
!
hgre = 0.d0
dhgre = 0.d0
dhgr2 = 0.d0

!
if(radic.lt.radiv)then
do ih =1, 13
hgre     = hgre + gcoef(ih)*ratir**(ih+11)
dhgre(1) = dhgre(1) + drdx/radiv*gcoef(ih)*(ih+11)*ratir**(ih+10)
dhgre(2) = dhgre(2) + drdy/radiv*gcoef(ih)*(ih+11)*ratir**(ih+10)

dhgr2(1) = dhgr2(1) + (drdx/radiv)**2*gcoef(ih)*(ih+11)*(ih+10)*ratir**(ih+9) +&
gcoef(ih)*(ih+11)*ratir**(ih+10)*drdx2/radiv

dhgr2(2) = dhgr2(2) + (drdy/radiv)**2*gcoef(ih)*(ih+11)*(ih+10)*ratir**(ih+9) +&
gcoef(ih)*(ih+11)*ratir**(ih+10)*drdy2/radiv
dhgr2(3) = dhgr2(3) + drdy*drdx/radiv**2*gcoef(ih)*(ih+11)*(ih+10)*ratir**(ih+9) +&
gcoef(ih)*(ih+11)*ratir**(ih+10)*drdxy/radiv

enddo
else
do ih =1, 13
hgre     = hgre + gcoef(ih)
enddo
dhgre = 0.d0
dhgr2 = 0.d0
endif

end subroutine greshoh
!
!...Get the gradient at cell center...
!
subroutine getunkct(unkct, xc, yc, xct, yct)
use constant
implicit none
real*8, dimension(1:ndegr, 1:nq), intent(out)::unkct
real*8:: xc,yc,xct,yct
!
real*8::rhomc, pctr, uctr, vctr
real*8:: dp0dx, dp0dy, du0dx,du0dy,dv0dx,dv0dy
real*8:: du0dr, du0ds, dv0dr, dv0ds,de0dx,de0dy
real*8:: drho0dx, drho0dy, drhomdx, drhomdy
real*8:: drho0dx2, drho0dy2, drho0dxy
real*8:: drhomdx2, drhomdy2, drhomdxy
real*8:: du0dx2, du0dy2, du0dxy
real*8:: dv0dx2, dv0dy2, dv0dxy
real*8:: de0dx2, de0dy2, de0dxy
real*8:: dp0dx2, dp0dy2, dp0dxy
real*8:: denom
integer::nind, nind1
real*8::radi, radiv, radic,ratir,theta,vemag,drdx,drdy,dverdr,dverdr2
real*8::dtdx,dtdy
real*8::drdx2,drdy2,drdxy,dtdx2,dtdy2,dtdxy
real*8::hgre
real*8::dhgre(2),dhgr2(3)
real*8::radie, radii,radie2,radii2,radic2,sentr,rhoin,rhoex
real*8::prein,preex
real*8::rho0, rho0ba
real*8::dradxc,dradyc,rhoc,ectr
real*8::gam1,gam2,gam3
real*8::drho0bdx,drho0bdy,drho0bdx2,drho0bdy2,drho0bdxy
!
if(ncase.eq.1)then     !...TGV
!
drho0dx = 0.d0
drho0dy = 0.d0
!
rhomc = 1.d0
uctr = sin(pi*xc)*cos(pi*yc)
vctr =-cos(pi*xc)*sin(pi*yc)
pctr = 0.25d0*(cos(2.d0*pi*xc) + cos(2.d0*pi*yc)) +  1.d0
!
drhomdx = -rhomc**2*drho0dx
drhomdy = -rhomc**2*drho0dy
!
du0dx   = pi*cos(pi*xc)*cos(pi*yc)
du0dy   =-pi*sin(pi*xc)*sin(pi*yc)
!
dv0dx   = pi*sin(pi*xc)*sin(pi*yc)
dv0dy   =-pi*cos(pi*xc)*cos(pi*yc)
!
dp0dx = 0.25d0*(-2.d0*pi)*(sin(2.d0*pi*xc))
dp0dy = 0.25d0*(-2.d0*pi)*(sin(2.d0*pi*yc))
!
de0dx = (dp0dx*rhomc + pctr*drhomdx)/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
de0dy = (dp0dy*rhomc + pctr*drhomdy)/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)

!DGP2
if(npoly.eq.2)then
drho0dx2  = 0.d0
drho0dy2  = 0.d0
drho0dxy  = 0.d0
!
drhomdx2 = 0.d0
drhomdy2 = 0.d0
drhomdxy = 0.d0

du0dx2 = -pi**2*sin(pi*xc)*cos(pi*yc)
du0dy2 = -pi**2*sin(pi*xc)*cos(pi*yc)
du0dxy = -pi**2*cos(pi*xc)*sin(pi*yc)

dv0dx2 = pi**2*cos(pi*xc)*sin(pi*yc)
dv0dy2 = pi**2*cos(pi*xc)*sin(pi*yc)
dv0dxy = pi**2*sin(pi*xc)*cos(pi*yc)

dp0dx2 = -pi**2*(cos(2.d0*pi*xc))
dp0dy2 = -pi**2*(cos(2.d0*pi*yc))
dp0dxy = 0.d0

de0dx2 = (drhomdx2*pctr + 2.d0*drhomdx*dp0dx + rhomc*dp0dx2)/(gamlg-1.d0) + &
(du0dx**2 + uctr*du0dx2 + dv0dx**2 + vctr*dv0dx2)
de0dy2 = (drhomdy2*pctr + 2.d0*drhomdy*dp0dy + rhomc*dp0dy2)/(gamlg-1.d0) + &
(du0dy**2 + uctr*du0dy2 + dv0dy**2 + vctr*dv0dy2)
de0dxy = (drhomdxy*pctr + drhomdx*dp0dy+ drhomdy*dp0dx + rhomc*dp0dxy)/(gamlg-1.d0) + &
(du0dy*du0dx + uctr*du0dxy + dv0dy*dv0dx + vctr*dv0dxy)
endif
!
elseif(ncase.eq.2)then !...Shockless Noh problem...
!
drho0dx = 0.d0
drho0dy = 0.d0
!
rhomc = 1.d0
uctr = -xc
vctr = -yc
pctr = (gamlg-1.d0)
!
drhomdx = -rhomc**2*drho0dx
drhomdy = -rhomc**2*drho0dy
!
du0dx =-1.d0
du0dy = 0.d0
!
dv0dx =  0.d0
dv0dy = -1.d0
!
dp0dx = 0.d0
dp0dy = 0.d0
!
de0dx = (dp0dx*rhomc + pctr*drhomdx)/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
de0dy = (dp0dy*rhomc + pctr*drhomdy)/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)
!
!DGP2
if(npoly.eq.2)then
drho0dx2  = 0.d0
drho0dy2  = 0.d0
drho0dxy  = 0.d0
!
drhomdx2 = 2.d0*rhomc**3*drho0dx**2 - rhomc**2*drho0dx2
drhomdy2 = 2.d0*rhomc**3*drho0dy**2 - rhomc**2*drho0dy2
drhomdxy = 2.d0*rhomc**3*drho0dx*drho0dy - rhomc**2*drho0dxy

du0dx2 = 0.d0
du0dy2 = 0.d0
du0dxy = 0.d0

dv0dx2 = 0.d0
dv0dy2 = 0.d0
dv0dxy = 0.d0

dp0dx2 = 0.d0
dp0dy2 = 0.d0
dp0dxy = 0.d0

de0dx2 = (drhomdx2*pctr + 2.d0*drhomdx*dp0dx + rhomc*dp0dx2)/(gamlg-1.d0) + &
(du0dx**2 + uctr*du0dx2 + dv0dx**2 + vctr*dv0dx2)
de0dy2 = (drhomdy2*pctr + 2.d0*drhomdy*dp0dy + rhomc*dp0dy2)/(gamlg-1.d0) + &
(du0dy**2 + uctr*du0dy2 + dv0dy**2 + vctr*dv0dy2)
de0dxy = (drhomdxy*pctr + drhomdx*dp0dy+ drhomdy*dp0dx + rhomc*dp0dxy)/(gamlg-1.d0) + &
(du0dy*du0dx + uctr*du0dxy + dv0dy*dv0dx + vctr*dv0dxy)
endif
!
!print*,'deodx2',de0dx2,drhomdx2*pctr, 2.d0*drhomdx*dp0dx, rhomc*dp0dx2
!
elseif(ncase.eq.3)then !...
!
!
drho0dx = 0.d0
drho0dy = 0.d0
!
denom = sqrt(xc**2 + yc**2)
!
rhomc = 1.d0
uctr = -xc/sqrt(xc**2 + yc**2)
vctr = -yc/sqrt(xc**2 + yc**2)
pctr = 1.d-6
!pctr = (gamlg-1.d0)*pctr
!
drhomdx = -rhomc**2*drho0dx
drhomdy = -rhomc**2*drho0dy
!
du0dx = -1.d0/denom + xc**2/denom**3
du0dy = xc*yc/denom**3
!
dv0dx =  xc*yc/denom**3
dv0dy =  -1.d0/denom + yc**2/denom**3
!
dp0dx = 0.d0
dp0dy = 0.d0
!
de0dx = (dp0dx*rhomc + pctr*drhomdx)/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
de0dy = (dp0dy*rhomc + pctr*drhomdy)/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)
!
!DGP2
if(npoly.eq.2)then
drho0dx2  = 0.d0
drho0dy2  = 0.d0
drho0dxy  = 0.d0
!
drhomdx2 = 0.d0
drhomdy2 = 0.d0
drhomdxy = 0.d0

du0dx2 = 3.d0*xc/denom**3 - 3.d0*xc**3/denom**5
du0dy2 = xc/denom**3 - 3.d0*xc*yc**2/denom**5
du0dxy = yc/denom**3 - 3.d0*xc**2*yc/denom**5

dv0dx2 = yc/denom**3 - 3.d0*xc**2*yc/denom**5
dv0dy2 = 3.d0*yc/denom**3 - 3.d0*yc**3/denom**5
dv0dxy = xc/denom**3 - 3.d0*xc*yc**2/denom**5

dp0dx2 = 0.d0
dp0dy2 = 0.d0
dp0dxy = 0.d0

de0dx2 = (drhomdx2*pctr + 2.d0*drhomdx*dp0dx + rhomc*dp0dx2)/(gamlg-1.d0) + &
(du0dx**2 + uctr*du0dx2 + dv0dx**2 + vctr*dv0dx2)
de0dy2 = (drhomdy2*pctr + 2.d0*drhomdy*dp0dy + rhomc*dp0dy2)/(gamlg-1.d0) + &
(du0dy**2 + uctr*du0dy2 + dv0dy**2 + vctr*dv0dy2)
de0dxy = (drhomdxy*pctr + drhomdx*dp0dy+ drhomdy*dp0dx + rhomc*dp0dxy)/(gamlg-1.d0) + &
(du0dy*du0dx + uctr*du0dxy + dv0dy*dv0dx + vctr*dv0dxy)
endif

!
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
!rhoin = 6.31d-4
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
!sentr = 2.15d4
!
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radic2 = xc**2 + yc**2
!
rho0ba = (radie2-radic2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radic2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rho0 = rho0ba**(1.d0/(gamlg-1.d0))
!
drho0dx = 1.d0/(gamlg-1.d0)*rho0ba**(1.d0/(gamlg-1.d0)-1.d0)*&
(rhoex**(gamlg-1.d0)-rhoin**(gamlg-1.d0))/(radie**2-radii**2)*2.d0*xc
drho0dy = 1.d0/(gamlg-1.d0)*rho0ba**(1.d0/(gamlg-1.d0)-1.d0)*&
(rhoex**(gamlg-1.d0)-rhoin**(gamlg-1.d0))/(radie**2-radii**2)*2.d0*yc
!
drhomdx = -drho0dx/rho0/rho0
drhomdy = -drho0dy/rho0/rho0
!
!
du0dx = 0.d0
du0dy = 0.d0
!
dv0dx = 0.d0
dv0dy = 0.d0
!
dp0dx = sentr*gamlg*rho0**(gamlg-1.d0)*drho0dx
dp0dy = sentr*gamlg*rho0**(gamlg-1.d0)*drho0dy
!
uctr  = 0.d0
vctr  = 0.d0
pctr  = sentr*(rho0)**gamlg
!
de0dx = (dp0dx/rho0 + pctr*drhomdx)/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
de0dy = (dp0dy/rho0 + pctr*drhomdy)/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)
!
!DGP2
if(npoly.eq.2)then
!
rhoc = rho0
rhomc = 1.d0/rhoc
!
drho0bdx = (rhoex**(gamlg-1.d0)-rhoin**(gamlg-1.d0))/(radie**2-radii**2)*2.d0*xc
drho0bdy = (rhoex**(gamlg-1.d0)-rhoin**(gamlg-1.d0))/(radie**2-radii**2)*2.d0*yc

drho0bdx2 = (rhoex**(gamlg-1.d0)-rhoin**(gamlg-1.d0))/(radie**2-radii**2)*2.d0
drho0bdy2 = (rhoex**(gamlg-1.d0)-rhoin**(gamlg-1.d0))/(radie**2-radii**2)*2.d0
drho0bdxy = 0.d0

gam1 = 1.d0/(gamlg-1.d0)
gam2 = (2.d0-gamlg)/(gamlg-1.d0)
gam3 = (3.d0-2.d0*gamlg)/(gamlg-1.d0)
drho0dx2  = gam1*(gam2*rho0ba**gam3*(drho0bdx)**2 + rho0ba**gam2*drho0bdx2)
drho0dy2  = gam1*(gam2*rho0ba**gam3*(drho0bdy)**2 + rho0ba**gam2*drho0bdy2)
drho0dxy  = gam1*(gam2*rho0ba**gam3*drho0bdy*drho0bdx)
!
drhomdx2 = 2.d0*rhomc**3*drho0dx**2 - rhomc**2*drho0dx2
drhomdy2 = 2.d0*rhomc**3*drho0dy**2 - rhomc**2*drho0dy2
drhomdxy = 2.d0*rhomc**3*drho0dx*drho0dy - rhomc**2*drho0dxy

du0dx2 = 0.d0
du0dy2 = 0.d0
du0dxy = 0.d0

dv0dx2 = 0.d0
dv0dy2 = 0.d0
dv0dxy = 0.d0

dp0dx2 = sentr*gamlg*(gamlg-1.d0)*rhoc**(gamlg-2.d0)*drho0dx**2 + sentr*gamlg*rhoc**(gamlg-1.d0)*drho0dx2
dp0dy2 = sentr*gamlg*(gamlg-1.d0)*rhoc**(gamlg-2.d0)*drho0dy**2 + sentr*gamlg*rhoc**(gamlg-1.d0)*drho0dy2
dp0dxy = sentr*gamlg*(gamlg-1.d0)*rhoc**(gamlg-2.d0)*drho0dy*drho0dx + sentr*gamlg*rhoc**(gamlg-1.d0)*drho0dxy

de0dx2 = (drhomdx2*pctr + 2.d0*drhomdx*dp0dx + rhomc*dp0dx2)/(gamlg-1.d0) !+ &
!(du0dx**2 + uctr*du0dx2 + dv0dx**2 + vctr*dv0dx2)
de0dy2 = (drhomdy2*pctr + 2.d0*drhomdy*dp0dy + rhomc*dp0dy2)/(gamlg-1.d0) !+ &
!(du0dy**2 + uctr*du0dy2 + dv0dy**2 + vctr*dv0dy2)
de0dxy = (drhomdxy*pctr + drhomdx*dp0dy+ drhomdy*dp0dx + rhomc*dp0dxy)/(gamlg-1.d0) !+ &
!(du0dy*du0dx + uctr*du0dxy + dv0dy*dv0dx + vctr*dv0dxy)
endif
!
elseif(ncase.eq.5)then !...Kidder ball...
!
radic = sqrt(xc**2+yc**2)
dradxc = xc/radic
dradyc = yc/radic

rhoc = 2.d0*exp(-radic**2)
!
drho0dx = rhoc*(-2.d0*radic*dradxc)
drho0dy = rhoc*(-2.d0*radic*dradyc)
!
rhomc = 1.d0/rhoc
uctr = 0.d0
vctr = 0.d0
ectr = 0.5d0
pctr = (gamlg-1.d0)*rhoc*ectr
!
drhomdx = -rhomc**2*drho0dx
drhomdy = -rhomc**2*drho0dy
!
du0dx = 0.d0
du0dy = 0.d0
!
dv0dx =  0.d0
dv0dy =  0.d0
!
de0dx = 0.d0
de0dy = 0.d0

!
elseif(ncase.eq.6)then !...Sod...
!
drhomdx = 0.d0
drhomdy = 0.d0
!
du0dx =  0.0d0
du0dy =  0.d0
!
dv0dx = 0.d0
dv0dy = 0.d0
!
de0dx = 0.d0
de0dy = 0.d0

!DGP2
if(npoly.eq.2)then
drho0dx2  = 0.d0
drho0dy2  = 0.d0
drho0dxy  = 0.d0
!
drhomdx2 = 0.d0
drhomdy2 = 0.d0
drhomdxy = 0.d0

du0dx2 = 0.d0
du0dy2 = 0.d0
du0dxy = 0.d0

dv0dx2 = 0.d0
dv0dy2 = 0.d0
dv0dxy = 0.d0

dp0dx2 = 0.d0
dp0dy2 = 0.d0
dp0dxy = 0.d0

de0dx2 = (drhomdx2*pctr + 2.d0*drhomdx*dp0dx + rhomc*dp0dx2)/(gamlg-1.d0) + &
(du0dx**2 + uctr*du0dx2 + dv0dx**2 + vctr*dv0dx2)
de0dy2 = (drhomdy2*pctr + 2.d0*drhomdy*dp0dy + rhomc*dp0dy2)/(gamlg-1.d0) + &
(du0dy**2 + uctr*du0dy2 + dv0dy**2 + vctr*dv0dy2)
de0dxy = (drhomdxy*pctr + drhomdx*dp0dy+ drhomdy*dp0dx + rhomc*dp0dxy)/(gamlg-1.d0) + &
(du0dy*du0dx + uctr*du0dxy + dv0dy*dv0dx + vctr*dv0dxy)
endif

!
elseif(ncase.eq.7)then !...sedov
!
drhomdx = 0.d0
drhomdy = 0.d0
!
du0dx =  0.0d0
du0dy =  0.d0
!
dv0dx = 0.d0
dv0dy = 0.d0
!
de0dx = 0.d0
de0dy = 0.d0
!DGP2
if(npoly.eq.2)then
drho0dx2  = 0.d0
drho0dy2  = 0.d0
drho0dxy  = 0.d0
!
drhomdx2 = 0.d0
drhomdy2 = 0.d0
drhomdxy = 0.d0

du0dx2 = 0.d0
du0dy2 = 0.d0
du0dxy = 0.d0

dv0dx2 = 0.d0
dv0dy2 = 0.d0
dv0dxy = 0.d0

dp0dx2 = 0.d0
dp0dy2 = 0.d0
dp0dxy = 0.d0

de0dx2 = (drhomdx2*pctr + 2.d0*drhomdx*dp0dx + rhomc*dp0dx2)/(gamlg-1.d0) + &
(du0dx**2 + uctr*du0dx2 + dv0dx**2 + vctr*dv0dx2)
de0dy2 = (drhomdy2*pctr + 2.d0*drhomdy*dp0dy + rhomc*dp0dy2)/(gamlg-1.d0) + &
(du0dy**2 + uctr*du0dy2 + dv0dy**2 + vctr*dv0dy2)
de0dxy = (drhomdxy*pctr + drhomdx*dp0dy+ drhomdy*dp0dx + rhomc*dp0dxy)/(gamlg-1.d0) + &
(du0dy*du0dx + uctr*du0dxy + dv0dy*dv0dx + vctr*dv0dxy)
endif

!
elseif(ncase.eq.8)then !...Gresho
!
nind = 6
nind1 = nind-1
radic = sqrt(xct**2 + yct**2)
radiv = 0.4d0
radi  = sqrt(xc**2 + yc**2)
!
drho0dx = 0.d0
drho0dy = 0.d0
!
rhomc = 1.d0
!
drhomdx = -rhomc**2*drho0dx
drhomdy = -rhomc**2*drho0dy
!
if(radic.lt.radiv)then
!
call facos(xc, yc, theta)
ratir = radi/radiv
vemag = 2**(2*nind)*ratir**nind*(1.d0-ratir)**nind
dverdr = 2**(2*nind)*nind/radiv*(ratir-ratir**2)**nind1*(1.d0-2.d0*ratir)
dverdr2= 2**(2*nind)*nind/radiv**2*(ratir-ratir**2)**nind1*&
(nind1*(1.d0-2.d0*ratir)**2/(ratir-ratir**2)-2.d0)
!
drdx = xc/radi
drdy = yc/radi
drdx2 = 1.d0/radi - xc/radi**2*drdx
drdy2 = 1.d0/radi - yc/radi**2*drdy
drdxy= -xc/radi**2*drdy
!
dtdx =-yc/radi**2
dtdy = xc/radi**2
dtdx2 = 2.d0*yc/radi**3*drdx
dtdy2 =-2.d0*xc/radi**3*drdy
dtdxy= 1.d0/radi**2 - 2.d0*xc/radi**3*drdx

!
uctr = vemag*(-sin(theta))
vctr = vemag*( cos(theta))
!
du0dx   =-(dverdr*sin(theta)*drdx + vemag*cos(theta)*dtdx)!-dverdr*drdx*yc
du0dy   =-(dverdr*sin(theta)*drdy + vemag*cos(theta)*dtdy)
!
dv0dx   = dverdr*cos(theta)*drdx - vemag*sin(theta)*dtdx
dv0dy   = dverdr*cos(theta)*drdy - vemag*sin(theta)*dtdy

!DGP2
if(npoly.eq.2)then
drho0dx2  = 0.d0
drho0dy2  = 0.d0
drho0dxy  = 0.d0
!
drhomdx2 = 2.d0*rhomc**3*drho0dx**2 - rhomc**2*drho0dx2
drhomdy2 = 2.d0*rhomc**3*drho0dy**2 - rhomc**2*drho0dy2
drhomdxy = 2.d0*rhomc**3*drho0dx*drho0dy - rhomc**2*drho0dxy

du0dx2 =-(dverdr2*drdx*drdx*sin(theta) + dverdr*(cos(theta)*dtdx*drdx + sin(theta)*drdx2) + &
dverdr*cos(theta)*dtdx*drdx  +  vemag*(-sin(theta)*dtdx*dtdx+ cos(theta)*dtdx2))
du0dy2 =-(dverdr2*drdy*drdy*sin(theta) + dverdr*(cos(theta)*dtdy*drdy + sin(theta)*drdy2) + &
dverdr*cos(theta)*dtdy*drdy  + vemag*(-sin(theta)*dtdy*dtdy + cos(theta)*dtdy2))
du0dxy = -(dverdr2*drdx*drdy*sin(theta)+dverdr*(cos(theta)*dtdy*drdx  + sin(theta)*drdxy) + &
dverdr*cos(theta)*dtdx*drdy  + vemag*(-sin(theta)*dtdx*dtdy + cos(theta)*dtdxy))

dv0dx2 =  dverdr2*drdx*drdx*cos(theta) + dverdr*(-sin(theta)*dtdx*drdx + cos(theta)*drdx2) - &
dverdr*sin(theta)*dtdx*drdx -  vemag*( cos(theta)*dtdx*dtdx + sin(theta)*dtdx2)
dv0dy2 =  dverdr2*drdy*drdy*cos(theta) + dverdr*(-sin(theta)*dtdy*drdy + cos(theta)*drdy2) - &
dverdr*sin(theta)*dtdy*drdy -  vemag*( cos(theta)*dtdy*dtdy + sin(theta)*dtdy2)
dv0dxy =  dverdr2*drdx*drdy*cos(theta) + dverdr*(-sin(theta)*dtdy*drdx + cos(theta)*drdxy) - &
dverdr*sin(theta)*dtdx*drdy -  vemag*( cos(theta)*dtdx*dtdy + sin(theta)*dtdxy)
endif

else
!
uctr = 0.d0
vctr = 0.d0
!
du0dx   = 0.d0
du0dy   = 0.d0
!
dv0dx   = 0.d0
dv0dy   = 0.d0

!DGP2
if(npoly.eq.2)then
drho0dx2  = 0.d0
drho0dy2  = 0.d0
drho0dxy  = 0.d0
!
drhomdx2 = 2.d0*rhomc**3*drho0dx**2 - rhomc**2*drho0dx2
drhomdy2 = 2.d0*rhomc**3*drho0dy**2 - rhomc**2*drho0dy2
drhomdxy = 2.d0*rhomc**3*drho0dx*drho0dy - rhomc**2*drho0dxy

du0dx2 = 0.d0
du0dy2 = 0.d0
du0dxy = 0.d0

dv0dx2 = 0.d0
dv0dy2 = 0.d0
dv0dxy = 0.d0
endif

endif
!
call greshoh(xc, yc ,radic, hgre, dhgre, dhgr2)
!
pctr = 5.d0 + 2**(4*nind)*hgre
!
dp0dx = 2**(4.d0*nind)*dhgre(1)
dp0dy = 2**(4.d0*nind)*dhgre(2)
!
de0dx = (dp0dx*rhomc + pctr*drhomdx)/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
de0dy = (dp0dy*rhomc + pctr*drhomdy)/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)
!DGP2
if(npoly.eq.2)then
!
dp0dx2 = 2**(4.d0*nind)*dhgr2(1)
dp0dy2 = 2**(4.d0*nind)*dhgr2(2)
dp0dxy = 2**(4.d0*nind)*dhgr2(3)
!
de0dx2 = (drhomdx2*pctr + 2.d0*drhomdx*dp0dx + rhomc*dp0dx2)/(gamlg-1.d0) + &
(du0dx**2 + uctr*du0dx2 + dv0dx**2 + vctr*dv0dx2)
de0dy2 = (drhomdy2*pctr + 2.d0*drhomdy*dp0dy + rhomc*dp0dy2)/(gamlg-1.d0) + &
(du0dy**2 + uctr*du0dy2 + dv0dy**2 + vctr*dv0dy2)
de0dxy = (drhomdxy*pctr + drhomdx*dp0dy+ drhomdy*dp0dx + rhomc*dp0dxy)/(gamlg-1.d0) + &
(du0dy*du0dx + uctr*du0dxy + dv0dy*dv0dx + vctr*dv0dxy)
endif

!
elseif(ncase.eq.9)then !...Triple point..
!
drhomdx = 0.d0
drhomdy = 0.d0
!
du0dx =  0.0d0
du0dy =  0.d0
!
dv0dx = 0.d0
dv0dy = 0.d0
!
de0dx = 0.d0
de0dy = 0.d0
!
!DGP2
if(npoly.eq.2)then
drho0dx2  = 0.d0
drho0dy2  = 0.d0
drho0dxy  = 0.d0
!
drhomdx2 = 0.d0
drhomdy2 = 0.d0
drhomdxy = 0.d0

du0dx2 = 0.d0
du0dy2 = 0.d0
du0dxy = 0.d0

dv0dx2 = 0.d0
dv0dy2 = 0.d0
dv0dxy = 0.d0

dp0dx2 = 0.d0
dp0dy2 = 0.d0
dp0dxy = 0.d0

de0dx2 = 0.d0
de0dy2 = 0.d0
de0dxy = 0.d0
endif

!
elseif(ncase.eq.10)then !...Expansion...!
!
!
drho0dx = 0.d0
drho0dy = 0.d0
!
rhomc = 1.d0
uctr = 0.d0
vctr = 0.d0
pctr = 1.d0-(xc**2+yc**2)
!
drhomdx = -rhomc**2*drho0dx
drhomdy = -rhomc**2*drho0dy
!
du0dx = 0.d0
du0dy = 0.d0
!
dv0dx =  0.d0
dv0dy =  0.d0
!
dp0dx = -2.d0*xc
dp0dy = -2.d0*yc
!
de0dx = (dp0dx*rhomc + pctr*drhomdx)/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
de0dy = (dp0dy*rhomc + pctr*drhomdy)/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)
!
elseif(ncase.eq.11)then !...
!
drho0dx = 0.d0
drho0dy = 0.d0
!
denom = sqrt((xc-0.0d0)**2 + (yc-0.0d0)**2)
!
rhomc = 1.d0
uctr = 0.d0
vctr = 0.d0
pctr = (gamlg-1.d0)*1.d-5
!
drhomdx = -rhomc**2*drho0dx
drhomdy = -rhomc**2*drho0dy
!
du0dx = 0.d0
du0dy = 0.d0
!
dv0dx = 0.d0
dv0dy = 0.d0
!
dp0dx = 0.d0
dp0dy = 0.d0
!
de0dx = (dp0dx*rhomc + pctr*drhomdx)/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
de0dy = (dp0dy*rhomc + pctr*drhomdy)/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
drho0dx = 0.9999995d0*pi*cos(pi*xc)
drho0dy = 0.d0

rhoc = 1.d0 + 0.9999995d0*sin(pi*xc)
rhomc = 1.d0/rhoc
uctr = 0.d0
vctr = 0.d0
pctr = rhoc**gamlg

drhomdx = -rhomc**2*drho0dx
drhomdy = -rhomc**2*drho0dy

du0dx = 0.d0
du0dy = 0.d0

dv0dx = 0.d0
dv0dy = 0.d0

dp0dx = gamlg*rhoc**(gamlg-1.d0)*drho0dx
dp0dy = gamlg*rhoc**(gamlg-1.d0)*drho0dy

de0dx = (dp0dx*rhomc + pctr*drhomdx)/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
de0dy = (dp0dy*rhomc + pctr*drhomdy)/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)

!DGP2
if(npoly.eq.2)then
drho0dx2  =-0.9999995d0*pi**2*sin(pi*xc)
drho0dy2  = 0.d0
drho0dxy = 0.d0
!
drhomdx2 = 2.d0*rhomc**3*drho0dx**2 - rhomc**2*drho0dx2
drhomdy2 = 0.d0!2.d0*rhomc**3*drho0dy**2 - rhomc**2*drho0dy2
drhomdxy = 0.d0!2.d0*rhomc**3*drho0dx*drho0dy - rhomc**2*drho0dxy

du0dx2 = 0.d0
du0dy2 = 0.d0
du0dxy = 0.d0

dv0dx2 = 0.d0
dv0dy2 = 0.d0
dv0dxy = 0.d0

dp0dx2 = gamlg*(gamlg-1.d0)*rhoc**(gamlg-2.d0)*drho0dx**2 + gamlg*rhoc**(gamlg-1.d0)*drho0dx2
dp0dy2 = 0.d0!gamlg*(gamlg-1.d0)*rhoc**(gamlg-2.d0)*drho0dy**2 + gamlg*rhoc**(gamlg-1.d0)*drho0dy2
dp0dxy = 0.d0!gamlg*(gamlg-1.d0)*rhoc**(gamlg-2.d0)*drho0dy*drho0dx + gamlg*rhoc**(gamlg-1.d0)*drho0dxy

de0dx2 = (drhomdx2*pctr + 2.d0*drhomdx*dp0dx + rhomc*dp0dx2)/(gamlg-1.d0) !+ &
!(du0dx**2 + uctr*du0dx2 + dv0dx**2 + vctr*dv0dx2)
de0dy2 = 0.d0!(drhomdy2*pctr + 2.d0*drhomdy*dp0dy + rhomc*dp0dy2)/(gamlg-1.d0) + &
!(du0dy**2 + uctr*du0dy2 + dv0dy**2 + vctr*dv0dy2)
de0dxy = 0.d0!(drhomdxy*pctr + drhomdx*dp0dy+ drhomdy*dp0dx + rhomc*dp0dxy)/(gamlg-1.d0) + &
!(du0dy*du0dx + uctr*du0dxy + dv0dy*dv0dx + vctr*dv0dxy)
endif

elseif(ncase.eq.13)then !...Saltzman
drhomdx = 0.d0
drhomdy = 0.d0

du0dx =  0.0d0
du0dy =  0.d0

dv0dx = 0.d0
dv0dy = 0.d0

de0dx = 0.d0
de0dy = 0.d0
!
elseif(ncase.eq.14)then !......Coggeshall expansion problem
!
rhomc = 1.d0
uctr = -0.25d0*xc
vctr = -yc
pctr = (gamlg-1.d0)*(3.d0*xc/8.d0)**2
!
drhomdx = 0.d0
drhomdy = 0.d0
!
du0dx =  -0.25d0
du0dy =   0.d0
!
dv0dx =  0.d0
dv0dy = -1.d0
!
dp0dx = (gamlg-1.d0)*9.d0/32.d0*xc
dp0dy = 0.d0
!
de0dx = (dp0dx*rhomc + pctr*drhomdx)/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
de0dy = (dp0dy*rhomc + pctr*drhomdy)/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)
!
endif
!
unkct(2, 1) = drhomdx !
unkct(3, 1) = drhomdy
unkct(2, 2) = du0dx
unkct(3, 2) = du0dy
unkct(2, 3) = dv0dx
unkct(3, 3) = dv0dy
unkct(2, 4) = de0dx
unkct(3, 4) = de0dy
!
if(npoly.eq.2)then
unkct(4, 1) = drhomdx2 !
unkct(5, 1) = drhomdy2
unkct(6, 1) = drhomdxy

unkct(4, 2) = du0dx2
unkct(5, 2) = du0dy2
unkct(6, 2) = du0dxy

unkct(4, 3) = dv0dx2
unkct(5, 3) = dv0dy2
unkct(6, 3) = dv0dxy

unkct(4, 4) = de0dx2
unkct(5, 4) = de0dy2
unkct(6, 4) = de0dxy
endif
!
end subroutine getunkct
!
!...subroutine: new version of RHS in Lagrangian framework (mass center)...
!
subroutine getrhsdg_lag_mcnew(ustar,uchar,bface,unkno,unkaw,unkgd,strnq_devtp,rhsel,rhsgd,intfac,inpoel,iptri,ipqua,&
geofa,geoel,gesgt0,gesgq0,&
coord,coold,coora,&
esuv1, esuv2,ltelem,estri,esqua,amatr,itime)
use constant
implicit none

!...input arrays
integer*4,dimension(1:nbfai,nbfac)::bface
integer*4,dimension(1:nvtri,1:ntria), intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
integer::ltelem(3,ncell)
integer, dimension(1:nftri,1:ntria), intent(in)::estri
integer, dimension(1:nfqua,1:nquad), intent(in)::esqua
integer*4,dimension(1:nifai,1:nafac), intent(in)::intfac
real*8::uchar(1:nq)
real*8,dimension(1:ngeel,1:nsize)::geoel
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac),intent(in)::gesgt0
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq0
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkaw
real*8,dimension(1:ndegr,1:4,1:nsize),intent(inout)::unkgd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq, nquad)::strnq_devtp
real*8,dimension(1:ndegr,1:nq,1:ncell),intent(out)::rhsel
real*8,dimension(1:ndegr,1:4,1:ncell),intent(out)::rhsgd
real*8,dimension(1:ngefa,1:nafac)::geofa
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coold!...initial coordinates...
real*8,dimension(1:ndimn,1:npoin),  intent(inout)::coora
real*8,dimension(1:ndimn,1:npoin)::coord
real*8,dimension(nmatr,ncell),intent(inout)::amatr
real*8,dimension(1:ndimn,1:npoin)   ::ustar !...nodal velocity
real*8,dimension(1:ndegr,1:nq,1:ncell)::rhsor,rhsaw
real*8,dimension(1:ndegr,1:nq,1:nsize)::unori
real*8::m(3,3),unint(nq)
integer::itime

!...local integer
integer::ie,ideg,id,ipoin,iq,iunk,iloop,ncbad,ielem
real*8:: dudr, duds, dvdr, dvds
real*8:: dudr2,duds2,dudrs,dvdr2,dvds2,dvdrs

!...local real array
real*8, allocatable:: amatrv(:,:)
real*8, allocatable:: gflag(:,:)
real*8, allocatable:: lpnp(:,:,:,:,:)
real*8, allocatable:: gelag(:,:,:),fnew(:,:,:),fnewq(:,:,:)
real*8, allocatable:: gelagq(:,:,:)
real*8, allocatable:: gesgt(:,:,:), gesgq(:,:,:)
real*8, allocatable:: fstrz(:,:,:,:), fsqrz(:,:,:,:)
real*8, allocatable:: fstrzaw(:,:,:,:), fsqrzaw(:,:,:,:)
real*8, allocatable:: fssgt(:,:,:,:,:), fssgq(:,:,:,:,:)
real*8, allocatable:: fstar(:,:,:,:), fstarq(:,:,:,:), fcurv(:,:,:,:),fcurvq(:,:,:,:)
real*8, allocatable:: gstar(:,:,:,:,:)
real*8, allocatable:: aflim(:,:)
real*8, allocatable:: unlim(: ,: ,:)
real*8, allocatable:: unmax(: ,:), unmin(:, :)
real*8, allocatable:: afvec(: ,:, :)
real*8, allocatable:: ufgaus(:, :, :), ufgpt(:, :, :)
real*8, allocatable:: cocent(: ,:)
real*8, allocatable:: vnulq(: ,:, :)
real*8, allocatable:: presqs(: ,:, :),prests(: ,:, :) !...Pressure star at any node...
real*8, allocatable:: fstfgdu(: ,:, :),fsqfgdu(: ,:, :)
real*8, allocatable:: fstfg(: ,:, :),  fsqfg(: ,:, :)
real*8, allocatable:: ufgpq(: ,:, :)
!
allocate(amatrv(nmatr,ncell))
allocate(gflag(1:ngflg,1:nbfac))
allocate(gelag(1:3,1:ngelg,1:ntria+nbfac))
allocate(gelagq(1:3,1:ngelgq,1:nquad))
allocate(lpnp(1:ndimn, 1:ndegr, 1:2, 1:nvtri, 1:ntria))
allocate(fstar(1:ndimn,1:2,1:nvtri, 1:ntria))
allocate(fstarq(1:ndimn,1:2,1:nvqua, 1:nquad))
allocate(fcurv(1:ndimn,1:4,1:nvtri, 1:nelem))
allocate(fcurvq(1:ndimn,1:4,1:nvqua, 1:nquad))
allocate(fnew(1:ndimn,1:ngelg, 1:ntria),fnewq(1:ndimn,1:ngelgq, 1:nquad))
allocate(fstrz(1:ndimn,1:4,1:nvtri, 1:ntria), fsqrz(1:ndimn,1:4,1:nvqua, 1:nquad) )
allocate(fstrzaw(1:ndimn,1:2,1:nvtri, 1:ntria), fsqrzaw(1:ndimn,1:2,1:nvqua, 1:nquad) )
allocate(gstar(1:ndimn,1:ngausf/2,1:2,1:nvtri, 1:nelem))
allocate(aflim(1:nq+1, 1:nsize))
allocate(unlim(1:ndegr, 1:nq, 1:nelem+nbfac))
allocate(unmax(1:nq+2, 1:npoin), unmin(1:nq+2, 1:npoin))
allocate(afvec(2, 2, nsize))
allocate(ufgaus(1:ndimn,1:ngelgq,1:nquad), ufgpt(1:ndimn,1:ngelg,1:ntria))
allocate(gesgt(1:3,1:ngesgt,1:ntria+nbfac))
allocate(gesgq(1:3,1:ngesgq,1:nquad))
allocate(fssgt(1:ndimn,1:2,1:3,1:nsubt, 1:ntria))
allocate(fssgq(1:ndimn,1:2,1:4,1:nsubq, 1:nquad))
allocate(cocent(3, nsize))
allocate(vnulq(1:ndimn,1:nvqua,1:nquad))
allocate(presqs(2, 1:nvqua, 1:nquad), prests(2, 1:nvtri, 1:ntria))
allocate(fstfgdu(1:ndimn,1:ngsft, 1:ntria), fsqfgdu(1:ndimn,1:ngsfq, 1:nquad))
allocate(fstfg(1:ndimn,1:ngsft, 1:ntria), fsqfg(1:ndimn,1:ngsfq, 1:nquad))
allocate(ufgpq(1:2,1:ngsfq,1:nquad))
!
!...Initialize the rhsel
!
rhsel = 0.d0
rhsgd = 0.d0
!
!unkgd = 0.d0
!
!unkgd(1, 1, : ) = 1.d0
!unkgd(1, 2, : ) = 0.d0
!unkgd(1, 3, : ) = 0.d0
!unkgd(1, 4, : ) = 1.d0
!
!unkno(4:ndegr,:,:) = 0.d0
!unkgd(2:ndegr,:, : ) = 0.d0
!unkgd(1:ndegr,2:3, : ) = 0.d0
!
!print*,'unkno',unkgd(1:3,1:4,213),unkno(1:3,2,213),itime
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!...Part I:Linear cell
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if(ncurv==0)then
!...I.1: Riemann solver
select Case (nriem)
Case (1)  !...Burton & Morgan

!...nfint.eq.1: Anayltical integral of the face integral
if(nfint.eq.1)then

!...I.1.1: Get F^* N ds for all faces...
call getfnds_lag_mc_hybrid(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)

!...I.1.2: Get the nodal velocity U_p^* ...
if(nlimi.ne.0)then
if(nlimi.eq.6)then

!...Reduce to P0
!unkno(2:3,:,:) = 0.d0

if(npoly.ge.1)then
!...Limiting part
!call barthlimit_lag_vtxunk(unkno, iptri, ipqua, unmax, unmin)
call barthlimit_lag_vtxunknew(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

!...Tria limiting
if(ntria.gt.0)call barthlimit_lagsym_tria(geoel, coord, coold, ustar, unkno, iptri, bface, intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)

!...Quad limiting
if(nquad.gt.0)call barthlimit_lagsym_quad(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)
!

!...Riemann invariant
!call barthlimit_lag_riemaninv(geoel, coord,unkno, ipqua,&
!bface, intfac, esuv1, esuv2)
endif

!++++++++++Print for debugging++++++++++
!print*,'ustar1',ie, unkno(:,1,32),unkno(:,1,42)

!...WENO
!call wenop1_rieminvrnt_tria_shu(iptri, estri, unkno, geoel, coord, coold, esuv1, esuv2)
!call wenop1_rieminvrnt_quad_shu6(ipqua, esqua, unkno, geoel, coord, coold,esuv1, esuv2)
!call weno_char_quad(ipqua, esqua, unkno, geoel, coord, coold, esuv1, esuv2)

!...Adjust limiter
!aflim(1,:) = 1.d0

!print*,'physical',unkno(:, :, 3510)

do ie = 1,-ncell
!...'if' is used to exclude the cells near the boundary
!if(geoel(10, ie).lt.100)then
aflim(:,ie) = 1.d0
afvec(:,:,ie) = 0.d0
afvec(1,1,ie) = 1.d0
afvec(2,2,ie) = 1.d0
!endif
enddo

!...I.1.3: Get the nodal velocity and pressure...
!...(nrz.eq.0): Cartesian coordinates
if(nrz.eq.0)then
call getndvelo_lag_mc_matrixsym(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold,unkno,ustar, fstar, fstarq, aflim, afvec, itime)

!...(nrz.eq.1).or. (nrz.eq.2): R-Z coordinates...
elseif(nrz.eq.1)then
call getndvelo_lagmc_rz(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold, unkno,ustar, fstrz, fsqrz, aflim, afvec, itime)

elseif(nrz.eq.2)then
!
if(itime.gt.1.or.rkstg.gt.1)then
  call getgeoel_lagmc_rzawrev(inpoel, iptri, ipqua, geoel, coord, unkno, aflim)
  call getamatr_lagmc_rzawrev(unkno,amatr,geoel,coord,inpoel, iptri, ipqua, aflim)
endif
!
call getndvelo_lagmc_rzaw(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold, unkno,ustar, fstrzaw, fsqrzaw, aflim, afvec, itime)

!...Get the mass center at t+epsaw*dt
!if(itime.eq.1.and.rkstg.eq.1)then
!coora = coord
!unkaw = unkno
!endif
!call getgeoel_lagmc_rzawdt(inpoel, iptri, ipqua, geoel, coord, unkno, ustar, aflim)
endif

!++++++++++Print for debugging++++++++++
!print*,'ustar2',ie,ie,rhsel(1,:,1891),rhsel(1,:,1831)

!...I.1.4: Face integral for RHS
if(nrz.eq.0)then
if(ntria.gt.0) call rhsifacedg_lag_mc_hybridtria(iptri, unkno, ustar, fstar, gelag, geoel,&
rhsel)
if(nquad.gt.0) call rhsifacedg_lag_mc_hybridquad2(ipqua, unkno, ustar,fstarq, gelagq, geoel,&
rhsel)

elseif(nrz.eq.1)then
if(ntria.gt.0) call rhsifacedg_lagmc_triarz(iptri, unkno, ustar, fstrz, gelag, geoel,coord,&
rhsel)
!if(nquad.gt.0) call rhsifacedg_lagmc_quadgauss(ipqua, unkno, ustar,fstarq, gelagq, geoel,coord,&
!rhsel)
!if(nquad.gt.0) call rhsifacedg_lagquadrz_simpson(ipqua,  unkno, ustar, fstarq, gelagq, geoel, coord,&
!rhsel)
if(nquad.gt.0) call  rhsifacedg_lagmc_quadrz(ipqua, unkno, ustar,fsqrz, gelagq, geoel,coord,&
rhsel)

elseif(nrz.eq.2)then
if(ntria.gt.0) call rhsifacedg_lagmc_triarzaw(iptri, unkno, ustar, fstrzaw, gelag, geoel,&
rhsel)
if(nquad.gt.0) call rhsifacedg_lagmc_quadrzaw(ipqua, unkno, ustar,fsqrzaw, gelagq, geoel,&
rhsel)
endif
!+++++++++Print for debugging++++++++++
!print*,'ustar2',ie,rhsel(1,:,1891),rhsel(1,:,1831)

!...I.1.5: Domain integral for RHS
if(npoly.ge.1)then

if(nrz.eq.0.or.nrz.eq.1)then

if(ntria.gt.0) call rhsdomndg_lag_mc_triasym(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
if(nquad.gt.0) call rhsdomndg_lag_mc_quadsym(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )

elseif(nrz.eq.2)then

if(ntria.gt.0) call rhsdomndg_lagmc_triarzaw(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
if(nquad.gt.0) call rhsdomndg_lagmc_quadrzaw(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )

!...Get the RHS contribution from the total derivative of the basis function
!if(ntria.gt.0) call rhsdomndg_lagmc_triarzawdt(intfac, iptri, coord,coold, geoel, unkno, ustar, rhsel,aflim, afvec )
!if(nquad.gt.0) call rhsdomndg_lagmc_quadrzawdt(intfac, ipqua, coord, coold, geoel, unkno, ustar, rhsel,aflim,afvec )
endif
endif
!++++++++++Print for debugging++++++++++
!print*,'Rhs after the domian integral',rhsel(:,2,991),itime

!...I.1.6: Source term for RZ...
if(nrz.eq.1)then

if(ntria.gt.0) call rhsdomnsrcdg_lagmc_triarz(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
if(nquad.gt.0) call rhsdomnsrcdg_lagmc_quadrz(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )

!...Special treatment of soure term for R-Z in true volume from Juan and Cheng
!if(ntria.gt.0) call getsigmas_triarz(iptri, geoel, gelag, unkno,  coord, coold, aflim, afvec, ustar, prests)
!if(nquad.gt.0) call getsigmas_quadrz(ipqua, geoel, gelagq, unkno, coord, coold, aflim, afvec, ustar, presqs)

!if(ntria.gt.0) call rhsdomnsrc_lagsympre_triarz(intfac, iptri, coord,coold, geoel,prests,rhsel )
!if(nquad.gt.0) call rhsdomnsrc_lagsympre_quadrz(intfac, ipqua, coord, coold, geoel,presqs, rhsel)
elseif(nrz.eq.2)then

!...Get the RHS...
do ie = 1, ncell
rhsel(:, 1, ie) = geoel(11, ie)*rhsel(:, 1, ie)
rhsel(:, 2, ie) = geoel(11, ie)*rhsel(:, 2, ie)
rhsel(:, 3, ie) = geoel(11, ie)*rhsel(:, 3, ie)
rhsel(:, 4, ie) = geoel(11, ie)*rhsel(:, 4, ie)
enddo

!print*,'call rhsifacedg_lag src3pre',rhsel(1,4,332)/sin(23.d0/64.d0*pi),rhsel(1,4,333)/sin(25.d0/64.d0*pi)
if(ntria.gt.0) call rhsdomnsrcdg_lagmc_triarzaw(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
if(nquad.gt.0) call rhsdomnsrcdg_lagmc_quadrzaw(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )
endif
!++++++++++Print for debugging++++++++++
!print*,'Rhs after the domian source integral',rhsel(1,1,1)

endif !if(nlimi.eq.6)then
endif !if(nlimi.ne.0)then

!...nfint.eq.2: Element-wise polynomial for the deformation gradient
elseif(nfint.eq.2)then
!...I.1.1: Get F^* N ds for all faces...
call getfnds_lag_hybridgd(geoel,gflag,gesgq,gesgt,intfac,iptri,ipqua,coord,unkgd) !...Notice

!...I.1.2: Get the nodal velocity U_p^* ...
!...Degenerate to P0
!unkno(2:3,:,:) = 0.d0

!...Limiting part
!call barthlimit_lag_vtxunk(unkno, iptri, ipqua, unmax, unmin)
call barthlimit_lag_vtxunknew(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

!...Tria limiting
if(ntria.gt.0)call barthlimit_lagsym_tria(geoel, coord, coold, ustar, unkno, iptri, bface, intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)

!...Quad limiting
if(nquad.gt.0)call barthlimit_lagsym_quad(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)!

!...Riemann invariant
!call barthlimit_lag_riemaninv(geoel, coord,unkno, ipqua,&
!bface, intfac, esuv1, esuv2)

!...Adjust limiter
!aflim(1,:) = 1.d0

do ie = 1,ncell
!...'if' is used to exclude the cells near the boundary
!if(geoel(10, ie).lt.100)then
aflim(:,ie) = 1.d0
afvec(:,:,ie) = 0.d0
afvec(1,1,ie) = 1.d0
afvec(2,2,ie) = 1.d0
!endif
enddo

!...I.1.3: Get the nodal velocity and pressure...
!...(nrz.eq.0): Cartesian coordinates
if(nrz.eq.0)then
call getndvelo_lag_gd(gflag,gesgt,gesgq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold, unkno, unkgd, strnq_devtp, ustar, fstar, fstarq, aflim, afvec, itime)

endif

!++++++++++Print for debugging++++++++++
!print*,'ustar2',unkno(1:3,2,213)
!print*,'ustar3',unkgd(1:3,1:4,280),gesgq(1:3,4,280)

!...I.1.4: Face integral for RHS
if(nrz.eq.0)then
if(nquad.gt.0) call rhsifacedg_lag_hybridquadl_gd(ipqua, unkno, ustar,fstarq, gesgq, geoel,&
rhsel)
endif
!+++++++++Print for debugging++++++++++
!print*,'call ustar hybrid',rhsel(2,2,213)

!...I.1.5: Domain integral for RHS
if(nrz.eq.0)then
if(nquad.gt.0) call rhsdomndg_lag_quadl_gd(intfac, ipqua, coord, coold, geoel, unkno, unkgd, strnq_devtp,rhsel,aflim,afvec )
endif
!++++++++++Print for debugging++++++++++
!print*,'Rhs after the domian integral',rhsel(2,2,213)!,rhsel(1,1,280)

!...I.1.6: RHS for the deformation gradient...

!...I.1.6.1: RHS from face integral for the deformation gradient...
if(nquad.gt.0) call rhsifacegd_lag_hybridquadl(ipqua, unkno, ustar,fstarq, gesgq0, geoel,&
rhsgd)

!...I.1.6.2: RHS from domain integral for the deformation gradient...
if(nquad.gt.0) call rhsdomngd_quadl(intfac, ipqua, coord, coold, geoel, unkno, unkgd,rhsgd,aflim,afvec )

!...nfint.eq.3: Gauss quadrature for the face integral in the Eulerian framework
elseif(nfint.eq.3)then

!...Degenerate to P0
!unkno(2:3,:,:) = 0.d0

!...I.2.1: Get F^* N ds for all faces...
call getfnds_lag_gausseulframe(gflag,gelag,gelagq,intfac,iptri,ipqua,coord)

!...I.2.2:limiting part
!call barthlimit_lag_vtxunknew(unkno, iptri, ipqua, unmax, unmin, coord, geoel)
!...Tria limiting
!if(ntria.gt.0)call barthlimit_lagsym_tria(geoel, coord, coold, ustar, unkno, iptri, bface, intfac, aflim, afvec,&
!unmax, unmin, esuv1, esuv2)
!...Quad limiting
!if(nquad.gt.0)call barthlimit_lagsym_quad(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec,&
!unmax, unmin, esuv1, esuv2)!

!...Adjust limiter
aflim(:,:) = 1.d0
!aflim(1,:) = 0.d0
afvec(:,:,:) = 0.d0
afvec(1,1,:) = 1.d0
afvec(2,2,:) = 1.d0

!...I.2.3: Get the nodal velocity U_p^* with Riemann solver...
call getndvelo_lag_gausseulframe(gflag,gelag,gelagq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, ufgpt, ufgaus, fnew, fnewq, aflim, afvec, itime)

!++++++++++Print for debugging++++++++++
!print*,'ufgaus',ustar(1:2,1010)!,munacu(1:2, ipoin),snsigm(1:2,ipoin)

!...I.2.4:Face integral for RHS
if(ntria.gt.0) call rhsifacedg_lagtriagseulframe(iptri,  unkno, ustar, ufgpt, geoel, gelag, fnew, coord, coold,&
rhsel)
if(nquad.gt.0) call rhsifacedg_lagquadgausseulframe(ipqua,  unkno, ustar, ufgaus, geoel, gelagq, fnewq, coord, coold,&
rhsel)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',ustar(1:2,1010)!rhsel(1:3,3,1:5),rhsel(3,4,1:5)!,rhsel(1,1:4,379)

!...I.2.5:Domain integral for RHS
if(ntria.gt.0)  call rhsdomndg_lag_mc_triasym(intfac, iptri, coord, coold, geoel, unkno, rhsel,aflim,afvec )
if(nquad.gt.0)  call rhsdomndg_lag_mc_quadsym(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsdomndg',rhsel(1:3,2:3,1),unkno(1:3,2:3,1)!,rhsel(3,4,1:5)!,rhsel(1,1:4,379)

!...Linear for RSF
elseif(nfint.eq.4)then

!...I.1.1: Get F^* N ds for all faces...
call getfnds_lag_mc_hybrid(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)

!...I.1.2: Get the nodal velocity U_p^* ...

!...Degenerate to P0
!unkno(2:3,:,:) = 0.d0

if(npoly.ge.1)then

!...Limiting part
call barthlimit_lag_vtxunknew(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

!...Tria limiting
if(ntria.gt.0)call barthlimit_lagsym_tria(geoel, coord, coold, ustar, unkno, iptri, bface, intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)

!...Quad limiting
if(nquad.gt.0)call barthlimit_lagsym_quad(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)
!

!...Riemann invariant
!call barthlimit_lag_riemaninv(geoel, coord,unkno, ipqua,&
!bface, intfac, esuv1, esuv2)
endif

!...Adjust limiter
!aflim(1,:) = 1.d0

do ie = 1,-ncell
aflim(:,ie) = 1.d0
afvec(:,:,ie) = 0.d0
afvec(1,1,ie) = 1.d0
afvec(2,2,ie) = 1.d0
enddo

!...I.1.3: Get the nodal velocity and pressure...
!...(nrz.eq.0): Cartesian coordinates
if(nrz.eq.0)then
call getRiemvtx_lag_rfm(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold,unkno,ustar, fstar, fstarq, aflim, afvec, itime)

!call getndvelo_lag_mc_matrixsym(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
!coord, coold,unkno,ustar, fstar, fstarq, aflim, afvec, itime)

!...(nrz.eq.1).or. (nrz.eq.2): R-Z coordinates...
elseif(nrz.eq.1)then
print*,'nrz.eq.1 will be implemented for RSF in future!'
stop
elseif(nrz.eq.2)then
print*,'nrz.eq.2 will be implemented for RSF in future!'
stop
endif

!++++++++++Print for debugging++++++++++
!print*,'ustar2',ustar(1:2,8)

!...I.1.4: Face integral for RHS
if(nrz.eq.0)then
if(ntria.gt.0) call rhsifacedg_lag_mc_hybridtria(iptri, unkno, ustar, fstar, gelag, geoel,&
rhsel)
if(nquad.gt.0) call rhsifacedg_lag_mc_hybridquad2(ipqua, unkno, ustar,fstarq, gelagq, geoel,&
rhsel)

elseif(nrz.eq.1)then
print*,'nrz.eq.1 will be implemented for RSF in future!'
stop
elseif(nrz.eq.2)then
print*,'nrz.eq.2 will be implemented for RSF in future!'
stop
endif
!+++++++++Print for debugging++++++++++
!if(itime.gt.1000)print*,'call ustar hybrid',itime,rhsel(1,1,19)

!...I.1.5: Domain integral for RHS
if(npoly.ge.1)then

if(nrz.eq.0.or.nrz.eq.1)then

if(ntria.gt.0) call rhsdomndg_lag_mc_triasym(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
if(nquad.gt.0) call rhsdomndg_lag_mc_quadsym(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec,itime)

elseif(nrz.eq.2)then

print*,'nrz.eq.2 will be implemented for RSF in future!'
stop

endif
endif
!++++++++++Print for debugging++++++++++
!if(itime.gt.1000)print*,'Rhs after the domian integral',itime,rhsel(1,1,19)

!...I.1.6: Source term for RZ will be implemented for RSF in future!...

!...New RSF
elseif(nfint.eq.5)then

!...I.1.1: Get F^* N ds for all faces...
call getfnds_lag_mc_hybrid(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)

!...reduce to P0
!unkno(2:ndegr, :, :) = 0.d0

!...I.1.2: Get the nodal velocity U_p^* ...
if(npoly.ge.1)then
!...Limiting part
call barthlimit_lag_vtxunknew(unkno, iptri, ipqua, unmax, unmin, coord, geoel)
!...Tria limiting
if(ntria.gt.0)call barthlimit_lagsym_tria(geoel, coord, coold, ustar, unkno, iptri, bface, intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)
!...Quad limiting
if(nquad.gt.0)call barthlimit_lagsym_quad(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)
endif

!...Adjust limiter
do ie = 1,-ncell
aflim(:,ie) = 1.d0
afvec(:,:,ie) = 0.d0
afvec(1,1,ie) = 1.d0
afvec(2,2,ie) = 1.d0
enddo

!...I.1.3: Get the nodal velocity and pressure...
!...(nrz.eq.0): Cartesian coordinates
if(nrz.eq.0)then
call getRiemvtx_lag_rfmfem(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold,unkno,ustar, fstar, fstarq, aflim, afvec, rhsel, itime)

!...(nrz.eq.1).or. (nrz.eq.2): R-Z coordinates...
elseif(nrz.eq.1)then
print*,'nrz.eq.1 will be implemented for RSF in future!'
stop
elseif(nrz.eq.2)then
print*,'nrz.eq.2 will be implemented for RSF in future!'
stop
endif

!++++++++++Print for debugging++++++++++
!print*,'ustar2',ie,rhsel(1,:,1891),rhsel(1,:,1831)

!...I.1.4: Face integral for RHS has been included in subroutine getRiemvtx_lag_rfm

!...I.1.5: Domain integral for RHS
if(npoly.ge.1)then

if(nrz.eq.0.or.nrz.eq.1)then

if(ntria.gt.0) call rhsdomndg_lag_mc_triasym(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
if(nquad.gt.0) call rhsdomndg_lag_mc_quadsym(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec,itime)

elseif(nrz.eq.2)then

print*,'nrz.eq.2 will be implemented for RSF in future!'
stop

endif
endif
!++++++++++Print for debugging++++++++++
!if(itime.gt.1000)print*,'Rhs after the domian integral',itime,rhsel(1,1,19)

!...I.1.6: Source term for RZ will be implemented for RSF in future!...

endif !...nfint

!...Store the unlimited unknown as original unknown (unori)...
unori = unkno

!...I.x.7: Updating the variables' derivative...
if(nlimi.eq.1)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
unkno(2:3, 2, ie) = unkno(2:3, 2, ie)*aflim(2,ie)
unkno(2:3, 3, ie) = unkno(2:3, 3, ie)*aflim(3,ie)
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
elseif(nlimi.eq.6.and.npoly.ge.1)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)

dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)

unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
endif

!...I.x.8(RZ AW only): The part from the total derivative of the basis function for the area-weighted
if(nrz.eq.2)then

call rhsdomndg_lagmc_rzawdt(ustar,unkno,rhsel,rhsaw,intfac,inpoel,iptri,ipqua,geoel,coord,coold,amatr, aflim, afvec)
!...Store the RHS without total derivative of the basis function...
rhsor = rhsel
!
do iloop =1,1
if(ndens.eq.1)then
!
 if(ncase.ne.4)then
  do ie =1, ncell
  do iq =1, nq  !...nq variables
  rhsel(1:ndegr,iq,ie)=dtfix*rhsel(1:ndegr,iq,ie)!*7.265d-3
  enddo
  enddo
 elseif(ncase.eq.4)then
  do ie =1, ncell
  do iq =1, nq  !...nq variables
  rhsel(1:ndegr,iq,ie)=dtfix*rhsel(1:ndegr,iq,ie)*tfcus
  enddo
  enddo
 endif

!...for high-order dg...
if(npoly.ge.1)then !...(2)npoly.ge.1
 do ie=1,ncell    !...(3)ie=1,nelem
  if(npoly==1)then
   m(1,1) = amatr(4, ie)
   m(1,2) = 0.d0
   m(1,3) = 0.d0

   m(2,1) = m(1,2)
   m(2,2) = amatr(1,ie)
   m(2,3) = amatr(2,ie)

   m(3,1) = m(1,3)
   m(3,2) = amatr(2,ie)
   m(3,3) = amatr(3,ie)
endif

!...Solve the mdegr independant varaible...
do id =1,ndegr !...(4)id =1,ndegr
!...step 1
 unint = 0.d0
 do iunk = 1,ndegr
  unint(1:nq) = unint(1:nq) + m(id, iunk)*rhsel(iunk,1:nq,ie)
 enddo
!...step 2
 unkaw(id,1:nq,ie)= unori(id,1:nq,ie) +  unint(1:nq)*(epsaw)
enddo !...(4)id =1,ndegr
enddo !...(3)ie=1,nelem
endif !...(2)npoly.ge.1
!
endif !...ndens.eq.1
!++++++++++Print for debugging++++++++++
!print*,'unold',unold(1,1,1)!unold(2:3,:,:) = 0.d0

!...Updating the physical coordinates based on dx/dt=u
if(ncase.ne.4)then
do ipoin = 1, npoin
   coora(1:2, ipoin) =  coord(1:2, ipoin) + dtfix*(epsaw)*ustar(1:2, ipoin)
enddo
elseif(ncase.eq.4)then
do ipoin = 1, npoin
   coora(1:2, ipoin) = coord(1:2, ipoin) + dtfix*tfcus*(epsaw)*ustar(1:2, ipoin)
enddo
endif
!
if(ndens.eq.3)then
call getgeo_denevoldt(iptri, ipqua, geoel, coora)
call getamatr_denevoldt(amatrv,geoel,coora, iptri, ipqua)
call rhs_denevoldt(iptri, ipqua, coora, coold, geoel, unkaw, amatrv)
endif

!...Get the mass center at t+dt
call getgeoel_lagmc_rzawdt(inpoel, iptri, ipqua, geoel, coora, unkaw, ustar, aflim)

!...Zero out the rhs for the area-weighted RZ method
!rhsaw = 0.d0
!if(ntria.gt.0) call rhsdomndg_lagmc_triarzawdt(intfac, iptri, coord, coold, geoel, unori, ustar, rhsaw,aflim,afvec )
!if(nquad.gt.0) call rhsdomndg_lagmc_quadrzawdt(intfac, ipqua, coord, coold, geoel, unori, ustar, rhsaw,aflim,afvec )

!...Get the new rhsel with the contribution from the total derivate of the basis function...
do ie = 1, ncell
rhsel(:, 1, ie) = rhsor(:, 1, ie)  + geoel(11, ie)*rhsaw(:, 1, ie)
rhsel(:, 2, ie) = rhsor(:, 2, ie)  + geoel(11, ie)*rhsaw(:, 2, ie)
rhsel(:, 3, ie) = rhsor(:, 3, ie)  + geoel(11, ie)*rhsaw(:, 3, ie)
rhsel(:, 4, ie) = rhsor(:, 4, ie)  + geoel(11, ie)*rhsaw(:, 4, ie)
enddo
!
enddo!iloop
!++++++++++Print for debugging++++++++++
!print*,'Rhs after the domian source integral 2',rhsel(1,4,6)
endif !...nrz.eq.2

!...I.x.9: Domain integral for source term (TGV)
if(ncase .eq. 1) then
if(ntria.gt.0) call rhsdomnsrcdg_lag_mc_hybridtria(intfac, iptri, coord, geoel,rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_mc_hybridquad(intfac, ipqua, coord, geoel,rhsel)
endif
!
Case (2) !...Marie
!
print*,'Implement Marie Riemann solver in futrue...'
stop
!
end select
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!...Part II: Curved cell
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
elseif(ncurv==1)then
!...II.1: Riemann solver
select Case (nriem)
Case (1)  !...Burton & Morgan

!...nfint.eq.1: Analytical integration for face integral
if(nfint.eq.1)then
  print*,'This solver chocie (ncurv=1, nriem=1, nfint.eq.1) is implemented by (ncurv=1, nriem=2, nfint.eq.1) !!!'

!...nfint.eq.2: Simpson rule integration for face integral
elseif(nfint.eq.2)then

!...II.1.1: Get F^* N ds for all faces...

!call getfnds_lag_simpson(gflag,gelag,intfac,inpoel,coord)
 call getfnds_lag_simpsonh(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)

!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...Quad limiting
!!call barthlimit_lag_curvquad(geoel, coord, coold, ustar, unkno, ipqua, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
 call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...Riemann invariant limiting
!call barthlimit_lag_riemaninv(geoel, coord,unkno, ipqua,&
!bface, intfac, esuv1, esuv2)

!...Adjust limiter
!aflim(:,:) = 1.d0
!afvec(:,:,:) = 0.d0
!afvec(1,1,:) = 1.d0
!afvec(2,2,:) = 1.d0

!...II.1.3: Get the nodal velocity and pressure with Riemann solver...
!call getndvelo_lag_simpson(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fnew)
 call getndvelo_lag_simpsonh(gflag,gelag,gelagq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, fnew, fnewq, aflim, afvec, vnulq, itime)

!++++++++++Print for debugging++++++++++
!print*,'unkno after limiter-velo_lag',itime,ustar(1,ipqua(:,2))

!...II.1.4: Face integral for RHS
!call rhsifacedg_lag_simpson(inpoel,  unkno, ustar, fnew, gelag,rhsel)
 call rhsifacedg_lagtria_simpson(inpoel,  unkno, ustar, fnew, gelag, geoel,&
 rhsel)
 call rhsifacedg_lagquad_simpson(ipqua,  unkno, ustar, fnewq, gelagq, geoel,&
 rhsel)
!call rhsifacedg_lagquadsimp(ipqua,  unkno, ustar,  geoel, coord, coold, aflim, afvec,&
!rhsel)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',rhsel(2,2,1),unkno(2,2,1),unkno(:,4,1)

!...II.1.5: Domain integral for RHS

if(ntria.gt.0) call rhsdomndg_lagmc_curvtria(intfac, inpoel, coord, geoel, unkno, rhsel, aflim,afvec)
if(nquad.gt.0) call rhsdomndg_lagmc_curvquad(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec ,vnulq )

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',rhsel(2,2,1)

!...II.1.6: Domain integral for source term
if(ncase .eq. 1)then
if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!...II.1.7: Update the 1st derivative
if(nlimi.eq.6)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)

dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
!
enddo
endif

!...nfint.eq.3: Numercial Gauss integration (4) for face integral
elseif(nfint==3)then

!...II.3.1: Get F^* N ds for all faces...
call getfnds_lag_gaussh(gflag,gelag,gelagq,intfac,iptri,ipqua,coord)

!...Lmiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

!...Quad
call barthlimit_lag_curvquad(geoel, coord, coold, ustar, unkno, ipqua, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...Adjust limiter
aflim(:,:) = 1.d0
!aflim(1,:) = 0.d0
afvec(:,:,:) = 0.d0
afvec(1,1,:) = 1.d0
afvec(2,2,:) = 1.d0

!...II.3.2: Get the nodal velocity and pressure with Riemann solver...
call getndvelo_lag_gaussh(gflag,gelag,gelagq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, ufgaus, fnew, fnewq, aflim, afvec, itime)

!++++++++++Print for debugging++++++++++
!print*,'unkno after limiter-velo_lag',itime, unkno(1:3,1:4,1)

!...II.3.3: Face integral for RHS...
call rhsifacedg_lagquadgauss(ipqua,  unkno, ustar, ufgaus, geoel, gelagq, fnewq, coord, coold,&
rhsel)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',rhsel(1,1:4,22),rhsel(1,1:4,379)

!...II.3.4: Domain integral for RHS
if(nquad.gt.0)  call rhsdomndg_lagmc_curvquad(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )

!...II.3.5: Domain integral for source term
if(ncase .eq. 1)then
 if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
 if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!...II.3.6: Update the 1st derivative
if(nlimi.eq.6)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)

dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)

unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
endif

!...nfint.eq.4: Numercial integration using subgrid method
elseif(nfint.eq.4)then

!...Degenerate to P0
if(ncase.eq.13)then
unkno(2:ndegr,:,  1) = 0.d0
unkno(2:ndegr,:,100:101) = 0.d0
unkno(2:ndegr,:,200:201) = 0.d0
unkno(2:ndegr,:,300:301) = 0.d0
unkno(2:ndegr,:,400:401) = 0.d0
unkno(2:ndegr,:,500:501) = 0.d0
unkno(2:ndegr,:,600:601) = 0.d0
unkno(2:ndegr,:,700:701) = 0.d0
unkno(2:ndegr,:,800:801) = 0.d0
unkno(2:ndegr,:,900:901) = 0.d0
unkno(2:ndegr,:,1000) = 0.d0
elseif(ncase.eq.-13)then
unkno(2:ndegr,:,  1) = 0.d0
unkno(2:ndegr,:,101) = 0.d0
unkno(2:ndegr,:,201) = 0.d0
unkno(2:ndegr,:,301) = 0.d0
unkno(2:ndegr,:,401) = 0.d0
unkno(2:ndegr,:,501) = 0.d0
unkno(2:ndegr,:,601) = 0.d0
unkno(2:ndegr,:,701) = 0.d0
unkno(2:ndegr,:,801) = 0.d0
unkno(2:ndegr,:,901) = 0.d0
endif

!...Dgenerate to P1
!unkno(5:ndegr,:,  :) = 0.d0
!unkno(:,3,:) = 0.d0
!print*,'density',unkno(2:6,1,1)

!...II.4.1: Get F^* N ds for all faces...
! call getfnds_lagsubgc(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
!...xxxxcall getfnds_lag_simpsubg(gflag,gesgq, gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)
! call getfnds_lag_simpsonh3(gflag,gesgq,intfac,inpoel,iptri,ipqua,coord)

 call getfnds_lagsmsif_hybrid(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
 call getfnds_lagsms_simpsonhybrid(gflag,gesgq,gesgt,intfac,inpoel,iptri,ipqua,coord)

!...Limiting part...
if(npoly.ge.1) call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

if(npoly.eq.1)then

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!...Quad limiting
!call barthlimit_lag_curvquad(geoel, coord, coold, ustar, unkno, ipqua, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...All the cells will be limited
geoel(10, :) = 100.d0
!call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)

!...Deacitivate the limiting for smooth region of DGp1
!do ie = 1,ncell

!...Smooth cell
!if(geoel(10, ie).lt.10.d0)then
!aflim(:, ie) = 1.d0

!afvec(1, 1, ie) = 1.d0
!afvec(1, 2, ie) = 0.d0
!afvec(2, 1, ie) = 0.d0
!afvec(2, 2, ie) = 1.d0
!endif

!enddo

elseif(npoly.gt.1)then
!...Troubled-cell marker
call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)

!...Remove the troubled-cells marker
!geoel(10, :) = 00.d0

!...Limiting the troubled-cells
!call barthlimit_lag_curvquadcb_HO(geoel, coord, coold, ustar, unkno, ipqua,bface, intfac, aflim, afvec, unmax, unmin, esuv1,! esuv2)
!

!...Riemann invariant limiting
!call barthlimit_lag_riemaninv(geoel, coord,unkno, ipqua,&
!bface, intfac, esuv1, esuv2)

endif

! call  weno_char_quad_orthg2(ipqua, esqua, unkno, geoel, coord, coold, esuv1, esuv2)
!...use weno p2
!call weno_char_quad(ipqua, esqua, unkno, geoel, coord, coold, esuv1, esuv2)
 call weno_char_quad_orthgl(ipqua, esqua, unkno, geoel, coord, coold, esuv1, esuv2)

!...Remove the troubled-cells marker for WENO
geoel(10, :) = 00.d0

!...Adjust limiter
aflim(:,:) = 1.d0
afvec(:,:,:) = 0.d0
afvec(1,1,:) = 1.d0
afvec(2,2,:) = 1.d0
!

!...II.4.2: Get the nodal velocity and pressure ...
call getndvelo_lagsubg(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, fssgt, fssgq, aflim, afvec, itime)

!++++++++++Print for debugging++++++++++
!print*,'unkno after limiter-velo_lag',itime,unkno(:,2,71)

!...II.4.3: Face integral for RHS...
if(ntria.gt.0) call rhsifacedg_lagsubgt3(iptri, unkno, ustar,fssgt, gesgt, geoel,&
rhsel)

if(nquad.gt.0) call rhsifacedg_lagsubgq2(ipqua, unkno, ustar,fssgq, gesgq, geoel,&
rhsel)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',rhsel(:,2,71)
!
if(npoly.ge.1)then
!...II.4.4: Domain integral for RHS
if(ntria.gt.0) call rhsdomndg_lagmc_curvtria(intfac, inpoel, coord, geoel, unkno, rhsel, aflim,afvec)

!if(nquad.gt.0)  call rhsdomndg_lagmc_qsubg(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )
if(nquad.gt.0)  call rhsdomndg_lagmc_curvquad(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec, vnulq)

endif

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',dtfix,rhsel(:,2,71)

!...II.4.5: Domain integral for source term
if(ncase .eq. 1)then
 if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
 if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!...II.4.6: Update the 1st derivative
if(nlimi.eq.6.and.npoly.ge.1)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!
dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
endif

!...nfint.eq.5: Numercial integration using subgrid method with both vertex and non-vertex gauss...
elseif(nfint.eq.5)then

!...Degenerate to P0
!unkno(4:ndegr,:,:) = 0.d0

if(ncase.eq.3)then
unkno(2:ndegr,:,1) = 0.d0
elseif(ncase.eq.13)then
unkno(2:ndegr,:,  1) = 0.d0
unkno(2:ndegr,:,100:101) = 0.d0
unkno(2:ndegr,:,200:201) = 0.d0
unkno(2:ndegr,:,300:301) = 0.d0
unkno(2:ndegr,:,400:401) = 0.d0
unkno(2:ndegr,:,500:501) = 0.d0
unkno(2:ndegr,:,600:601) = 0.d0
unkno(2:ndegr,:,700:701) = 0.d0
unkno(2:ndegr,:,800:801) = 0.d0
unkno(2:ndegr,:,900:901) = 0.d0
unkno(2:ndegr,:,1000) = 0.d0
endif

!...II.4.1: Get F^* N ds for all faces...
call getfnds_lagsubgc(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
!call getfnds_lag_simpsubg(gflag,gesgq, gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)
call getfnds_lag_simpsonh3(gflag,gesgq,intfac,inpoel,iptri,ipqua,coord)

!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

if(npoly.eq.1)then

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!...Quad limiting
!call barthlimit_lag_curvquad(geoel, coord, coold, ustar, unkno, ipqua, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...All the cells will be limited
geoel(10, :) = 100.d0
!call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)

!...Deacitivate the limiting for smooth region of DGp1
!do ie = 1,ncell
!!...Smooth cell
!if(geoel(10, ie).lt.10.d0)then
!aflim(:, ie) = 1.d0

!afvec(1, 1, ie) = 1.d0
!afvec(1, 2, ie) = 0.d0
!afvec(2, 1, ie) = 0.d0
!afvec(2, 2, ie) = 1.d0
!endif
!enddo

elseif(npoly.gt.1)then
!...Troubled-cell marker
call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)
!call getSIHO_facejmp3(geoel,intfac,ipqua,coord,unkno)
!
!call getSIHOamlp_kim(geoel, coord,unkno, ipqua, unmax, unmin)

!...Remove the troubled-cells marker
geoel(10, :) = 0.d0

!...Limiting the troubled-cells
call barthlimit_lag_curvquadcb_HO(geoel, coord, coold, ustar, unkno, ipqua,bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!

!...Riemann invariant limiting
!call barthlimit_lag_riemaninv(geoel, coord,unkno, ipqua,&
!bface, intfac, esuv1, esuv2)

!...All the cells will be unlimited
!geoel(10, :) = 0.d0
endif
!
!print*,'limiter',aflim(:,1)

!...Adjust limiter
aflim(:,:) = 1.d0
afvec(:,:,:) = 0.d0
afvec(1,1,:) = 1.d0
afvec(2,2,:) = 1.d0
!
!++++++++++Print for debugging++++++++++
!print*,'unkno limit',itime,unkno(:,:,173)

!...II.4.2: Get the nodal velocity and pressure ...
call getndvelo_lagsubg_gaussh(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno, unkgd,strnq_devtp, ustar, ufgpq, fssgt, fssgq,fstfgdu,fsqfgdu,fstfg,fsqfg,  aflim, afvec, itime)

!++++++++++Print for debugging++++++++++
!do ie =1, ncell
!if(geoel(10,ie).lt.20)print*,'unkno after limiter-velo_lag',itime,ie,geoel(10,ie)!, ustar(1:2,2419), ustar(1:2,2218)
!enddo
!...II.4.3: Face integral for RHS...
call rhsifacedg_lagsubgq_gs(ipqua, unkno, ustar, ufgpq,fssgq, fsqfg, fsqfgdu, gesgq, geoel, coord,&
rhsel)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',itime,unkno(1:3,2,1),rhsel(1:3,2,1)

!...II.4.4: Domain integral for RHS
!if(nquad.gt.0)  call rhsdomndg_lagmc_qsubg(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )
if(nquad.gt.0)  call rhsdomndg_lagmc_curvquad(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec, vnulq)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsdomndg_lag',itime,unkno(1:6,1,1:2),rhsel(1,1,1:2)

!...II.4.5: Domain integral for source term
if(ncase .eq. 1)then
if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!...II.4.6: Update the 1st derivative
if(nlimi.eq.6)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!
dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
endif
!
!...nfint.eq.6: Analytical integration using subgrid method with both vertex and non-vertex gauss...
elseif(nfint.eq.6)then

!...Degenerate to P0
!unkno(2:3,:,:) = 0.d0

!...II.4.1: Get F^* N ds for all faces...
call getfnds_lagsubgc_maire(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
call getfnds_lag_simpsonh_maire(gflag,gesgq,intfac,inpoel,iptri,ipqua,coord)
!
!print*,'hgeole',gesgq(:,1,1)

!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

if(npoly.eq.1)then

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!...Quad limiting
!call barthlimit_lag_curvquad(geoel, coord, coold, ustar, unkno, ipqua, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...All the cells will be limited
!geoel(10, :) = 100.d0
call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)

!...Deacitivate the limiting for smooth region of DGp1
do ie = 1,ncell

!...Smooth cell
if(geoel(10, ie).lt.10.d0)then
aflim(:, ie) = 1.d0

afvec(1, 1, ie) = 1.d0
afvec(1, 2, ie) = 0.d0
afvec(2, 1, ie) = 0.d0
afvec(2, 2, ie) = 1.d0
endif

enddo

elseif(npoly.gt.1)then
!
call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)
call barthlimit_lag_curvquadcb_HO(geoel, coord, coold, ustar, unkno, ipqua,bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!

!...Riemann invariant limiting
!call barthlimit_lag_riemaninv(geoel, coord,unkno, ipqua,&
!bface, intfac, esuv1, esuv2)

!...All the cells will be limited
!geoel(10, :) = 0.d0
endif

!...Adjust limiter
!aflim(:,:) = 1.d0
!afvec(:,:,:) = 0.d0
!afvec(1,1,:) = 1.d0
!afvec(2,2,:) = 1.d0
!

!...II.4.2: Get the nodal velocity and pressure ...

call getndvelo_lagsubg_gaussh_maire(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno, unkgd, ustar, rhsel, aflim, afvec, itime)

!++++++++++Print for debugging++++++++++
print*,'unkno after limiter-velo_lag',itime,unkno(:,2,71)!, ustar(1:2,2419), ustar(1:2,2218)

!...II.4.3: Face integral for RHS...
!...This method has been included in the subroutine getndvelo_lagsubg_gaussh_maire

!++++++++++Print for debugging++++++++++

!print*,'Output RHS after sub. rhsifacedg_lag',unkno(1:3,3,1),unkno(1:3,2,5)
!print*,'Output RHS after sub. rhsifacedg_lag2',rhsel(1:3,3,1),rhsel(1:3,2,5)
!rhsel = 0.d0

!...II.4.4: Domain integral for RHS
if(nquad.gt.0)  call rhsdomndg_lagmc_curvquad(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec, vnulq)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',rhsel(1:3,3,1),rhsel(1:3,2,10)

!...II.4.5: Domain integral for source term
if(ncase .eq. 1)then
if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!...II.4.6: Update the 1st derivative
if(nlimi.eq.6)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!
dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
endif

!...nfint.eq.7: Numercial integration using subcell method with both vertex and non-vertex gauss with smooth limiting...
elseif(nfint.eq.7)then

!...Degenerate to P0
!unkno(2:3,:,:) = 0.d0

!...II.4.1: Get F^* N ds for all faces...
call getfnds_lagsubgc(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
!call getfnds_lag_simpsubg(gflag,gesgq, gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)
call getfnds_lag_simpsonh3(gflag,gesgq,intfac,inpoel,iptri,ipqua,coord)

!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

if(npoly.eq.1)then

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!...Quad limiting
!call barthlimit_lag_curvquad(geoel, coord, coold, ustar, unkno, ipqua, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...All the cells will be limited
!geoel(10, :) = 100.d0
call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)

!...Deacitivate the limiting for smooth region of DGp1
do ie = 1,ncell

!...Smooth cell
if(geoel(10, ie).lt.10.d0)then
aflim(:, ie) = 1.d0

afvec(1, 1, ie) = 1.d0
afvec(1, 2, ie) = 0.d0
afvec(2, 1, ie) = 0.d0
afvec(2, 2, ie) = 1.d0
endif

enddo

elseif(npoly.gt.1)then
!
call getSIHOsmth_Persson(intfac, ipqua, coord, coold, geoel, unkno)
call barthlimit_lag_curvquadcb_HOsmth(geoel, coord, coold, ustar, unkno, ipqua,bface, &
                                      intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...All the cells will be unlimited
!geoel(10, :) = 0.d0
endif

!...Adjust limiter
!aflim(:,:) = 1.d0
!afvec(:,:,:) = 0.d0
!afvec(1,1,:) = 1.d0
!afvec(2,2,:) = 1.d0
!
!++++++++++Print for debugging++++++++++
!print*,'limiter',itime,geoel(10,116),aflim(:,116),afvec(:,:,116)

!...II.4.2: Get the nodal velocity and pressure ...
call getndvelo_lagsubg_gaussh(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno, unkgd,strnq_devtp, ustar, ufgpq, fssgt, fssgq,fstfgdu,fsqfgdu,fstfg,fsqfg,  aflim, afvec, itime)

!++++++++++Print for debugging++++++++++
!ncbad = 0
!  open(8,file='rhstxt.dat')
!do ie =1, ncell
!if(geoel(10,ie).gt.10)then
! ncbad = ncbad+1
!   write(8,*) itime,ie, geoel(10,ie)
!endif
!enddo
!   close(8)
!...II.4.3: Face integral for RHS...
call rhsifacedg_lagsubgq_gs(ipqua, unkno, ustar, ufgpq,fssgq, fsqfg, fsqfgdu, gesgq, geoel, coord,&
rhsel)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',itime,ncbad

!...II.4.4: Domain integral for RHS
!if(nquad.gt.0)  call rhsdomndg_lagmc_qsubg(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )
if(nquad.gt.0)  call rhsdomndg_lagmc_curvquad_smth(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec, vnulq)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsdomndg_lag',itime,unkno(1,1,116),rhsel(1,1,116)

!...II.4.5: Domain integral for source term
if(ncase .eq. 1)then
if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!...II.4.6: Update the 1st derivative
if(nlimi.eq.6)then
do ie = 1, ncell
unkno(2:ndegr, 1, ie) = unkno(2:ndegr, 1, ie)*aflim(1,ie)
!
dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
dudr2 = afvec(1, 1, ie)*unkno(4,2,ie) +  afvec(1, 2, ie)*unkno(4,3,ie)
duds2 = afvec(1, 1, ie)*unkno(5,2,ie) +  afvec(1, 2, ie)*unkno(5,3,ie)
dudrs = afvec(1, 1, ie)*unkno(6,2,ie) +  afvec(1, 2, ie)*unkno(6,3,ie)

dvdr2 = afvec(2, 1, ie)*unkno(4,2,ie) +  afvec(2, 2, ie)*unkno(4,3,ie)
dvds2 = afvec(2, 1, ie)*unkno(5,2,ie) +  afvec(2, 2, ie)*unkno(5,3,ie)
dvdrs = afvec(2, 1, ie)*unkno(6,2,ie) +  afvec(2, 2, ie)*unkno(6,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds; unkno(4, 2, ie) = dudr2; unkno(5, 2, ie) = duds2; unkno(6, 2, ie) = dudrs;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds; unkno(4, 3, ie) = dvdr2; unkno(5, 3, ie) = dvds2; unkno(6, 3, ie) = dvdrs;
unkno(2:ndegr, 4, ie) = unkno(2:ndegr, 4, ie)*aflim(5,ie)
enddo
endif
!
!...Restore the smooth indicator for the shocked cell
do ie = 1,ncell
!...Shock cells
if(geoel(10, ie).gt.15.d0)then

geoel(10,ie) = 0.d0
!...Transition cells
elseif(geoel(10, ie).gt.5.d0.and.geoel(10, ie).lt.15.d0)then

geoel(10,ie) = geoel(10,ie)-10.d0
endif

enddo
!
!...nfint.eq.8: MEMP2 with Simpson's rule...
elseif(nfint.eq.8)then

!...Degenrate the high-order method

!unkgd(4:ndegr,:,:) = 0.d0
!unkno(2:3,:,:) = 0.d0
!
!...II.8.1: Get F^* N ds for all faces...
call getfnds_laghybrid_gd(geoel,gflag,gesgq,gesgt,intfac,iptri,ipqua,coord,unkgd) !...Notice

!...Option: Get the impedance weighted averaged area
!call getgesgq(intfac,ipqua,gesgq)

!...II.8.2: Get the nodal velocity U_p^* ...
!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

if(npoly.eq.1)then

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...Quad limiting
call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...All the cells will be limited
geoel(10, :) = 100.d0

elseif(npoly.gt.1)then
!...Troubled-cell marker
call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)

!...Remove the troubled-cells marker
geoel(10, :) = 0.d0

!...Limiting the troubled-cells
call barthlimit_lag_curvquadcb_HO(geoel, coord, coold, ustar, unkno, ipqua,bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

endif

do ie = 1,ncell
!...'if' is used to exclude the cells near the boundary
!if(geoel(10, ie).lt.100)then
aflim(:,ie) = 1.d0
afvec(:,:,ie) = 0.d0
afvec(1,1,ie) = 1.d0
afvec(2,2,ie) = 1.d0
!endif
enddo

!...II.8.3: Get the nodal velocity and pressure...
!...(nrz.eq.0): Cartesian coordinates
if(nrz.eq.0)then
call getndvelo_lag_simpsonh_gd(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, fnew, fnewq, aflim, afvec,vnulq, itime)
endif

!++++++++++Print for debugging++++++++++
!print*,'ustar2',unkno(2,1,8)
!print*,'ustar3',ustar(:,20)

!...II.8.4: Face integral for RHS
if(nrz.eq.0)then
if(nquad.gt.0) call rhsifacedg_lagquad_simpson_gd(ipqua,  unkno, ustar, fnewq, gesgq, geoel,&
rhsel)
endif
!+++++++++Print for debugging++++++++++
!print*,'call ustar hybrid',rhsel(2,1,8)

!...II.8.5: Domain integral for RHS
if(nrz.eq.0)then
if(nquad.gt.0) call rhsdomndg_lagquadc_gd(intfac, ipqua, coord, coold, geoel, unkno, unkgd, rhsel,aflim,afvec, vnulq )
endif
!++++++++++Print for debugging++++++++++
!print*,'Rhs after the domian integral',itime,unkgd(:,1,8)

!
!...II.8.6: RHS for the deformation gradient...
!

!...II.8.6.1: RHS from face integral for the deformation gradient...
if(nquad.gt.0) call rhsifacegd_lagquadc(ipqua, ustar, gesgq0, geoel,&
rhsgd)
!print*,'face rhs for the gradient deformation',rhsgd(:,1,8)

!...II.8.6.2: RHS from domain integral for the deformation gradient...
if(nquad.gt.0) call rhsdomngd_lagquadc(intfac, ipqua, coord, coold, geoel, unkno, unkgd,rhsgd,aflim,afvec )
!print*,'Rhs after the domian integral',itime,rhsel(1:6,2,1),rhsgd(5,1,1)
!

!...II.8.7: Domain integral for source term
if(ncase .eq. 1)then
if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!...II.8.8: Update the 1st derivative
if(nlimi.eq.6)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!
dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
endif

!...nfint.eq.9: MEMP2 with Simpson's rule using the SMS...
elseif(nfint.eq.9)then

!...Degenrate the high-order method
! unkgd(:,:,:) = 0.d0
! unkgd(1,1,:) = 1.d0
! unkgd(1,4,:) = 1.d0
! unkno(4:ndegr,:,:) = 0.d0
!
!...II.9.1: Get F^* N ds for all faces...

call getfnds_lagsmsif_hybrid(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
call getfnds_lagsms_simpsonhybrid(gflag,gesgq,gesgt,intfac,inpoel,iptri,ipqua,coord)

!...Option: Get the impedance weighted averaged area
!call getgesgq(intfac,ipqua,gesgq)

!...II.9.2: Get the nodal velocity U_p^* ...
!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

if(npoly.eq.1)then

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...Quad limiting
call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...All the cells will be limited
geoel(10, :) = 100.d0

elseif(npoly.gt.1)then
!...Troubled-cell marker
call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)

!...Remove the troubled-cells marker
geoel(10, :) = 00.d0

!...Limiting the troubled-cells
call barthlimit_lag_curvquadcb_HO(geoel, coord, coold, ustar, unkno, ipqua,bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

endif

do ie = 1,ncell
!...'if' is used to exclude the cells near the boundary
!if(geoel(10, ie).lt.100)then
aflim(:,ie) = 1.d0
afvec(:,:,ie) = 0.d0
afvec(1,1,ie) = 1.d0
afvec(2,2,ie) = 1.d0
!endif
enddo

!...II.9.3: Get the nodal velocity and pressure...
!...(nrz.eq.0): Cartesian coordinates
if(nrz.eq.0)then
call getndvelo_lagsms_gd(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,unkgd,strnq_devtp, ustar, fssgt, fssgq, aflim, afvec, itime)
endif

!++++++++++Print for debugging++++++++++
!print*,'ustar3',ustar(:,1),ustar(:,2421)
!print*,'unkno',unkno(1,1:4,1),unkno(1,1:4,2281)

!...II.9.4: Face integral for RHS
if(nrz.eq.0)then
if(nquad.gt.0) call rhsifacedg_lagsms_simpson_gd(ipqua, unkno, ustar,fssgq, gesgq, geoel,&
rhsel)
endif
!+++++++++Print for debugging++++++++++
!print*,'call ustar hybrid',rhsel(2,1,8)

!...II.9.5: Domain integral for RHS
if(nrz.eq.0)then
!if(nquad.gt.0) call rhsdomndg_lagquadc_gd(intfac, ipqua, coord, coold, geoel, unkno, unkgd, strnq_devtp, rhsel,aflim,afvec, vnulq )
if(nquad.gt.0) call rhsdomndg_lagquadc_gdmix(intfac, ipqua, coord, coold, geoel, unkno, unkgd,&
                                             strnq_devtp, rhsel,aflim,afvec, vnulq )

endif
!++++++++++Print for debugging++++++++++
!print*,'Rhs after the domian integral',itime,unkgd(:,1,8)

!
!...II.9.6: RHS for the deformation gradient...
!

!...II.9.6.1: RHS from face integral for the deformation gradient...
if(nquad.gt.0) call rhsifacegd_lagquadc(ipqua, ustar, gesgq0, geoel,&
rhsgd)
!print*,'face rhs for the gradient deformation',rhsgd(:,1,8)

!...II.8.6.2: RHS from domain integral for the deformation gradient...
if(nquad.gt.0) call rhsdomngd_lagquadc(intfac, ipqua, coord, coold, geoel, unkno, unkgd,rhsgd,aflim,afvec )
!print*,'Rhs after the domian integral',itime,rhsel(1:6,2,1),rhsgd(5,1,1)
!

!...II.9.7: Domain integral for source term
if(ncase .eq. 1)then
if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!...II.9.8: Update the 1st derivative
if(nlimi.eq.6)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!
dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
endif

!...nfint.eq.10: MEMP2 with 5-point rule using the SMS...
elseif(nfint.eq.10)then

!...Degenrate the high-order method
! unkgd(4:ndegr,:,:) = 0.d0
! unkno(4:ndegr,:,:) = 0.d0
!
if(nmatel.eq.1.and.ncase.eq.3)then
  unkno(2:ndegr,:,1) = 0.d0
endif
!...II.10.1: Get F^* N ds for all faces...

call getfnds_lagsmsif_hybrid(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
call getfnds_lagsms_simpsonhybrid(gflag,gesgq,gesgt,intfac,inpoel,iptri,ipqua,coord)

!...Option: Get the impedance weighted averaged area
!call getgesgq(intfac,ipqua,gesgq)

!...II.10.2: Get the nodal velocity U_p^* ...
!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

if(npoly.eq.1)then

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...Quad limiting
call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...All the cells will be limited
geoel(10, :) = 00.d0

elseif(npoly.gt.1)then
!...Troubled-cell marker
!call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)
call getSIHOgd_Persson(intfac, ipqua, coord, coold, geoel, unkno,unkgd)
!call getSIHOamlp_kim(geoel, coord,unkno, ipqua, unmax, unmin)

!...Remove the troubled-cells marker
geoel(10, :) = 00.d0

!...Limiting the troubled-cells
call barthlimit_lag_curvquadcb_HO(geoel, coord, coold, ustar, unkno, ipqua,bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!
do ie = 1,-nquad
ielem = ie + ntria
!...Bad cell
! if(geoel(10, ielem).gt.10.d0)then
!   unkgd(2:ndegr, :, ielem) = 0.d0
! endif
enddo
!
endif


do ie = 1,ncell
!...'if' is used to exclude the cells near the boundary
!if(geoel(10, ie).lt.100)then
aflim(:,ie) = 1.d0
afvec(:,:,ie) = 0.d0
afvec(1,1,ie) = 1.d0
afvec(2,2,ie) = 1.d0
!endif
enddo

!...II.10.3: Get the nodal velocity and pressure...
!...(nrz.eq.0): Cartesian coordinates
if(nrz.eq.0)then
call getndvelo_lagsubg_gaussh(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno, unkgd, strnq_devtp, ustar, ufgpq, fssgt, fssgq,fstfgdu,fsqfgdu,fstfg,fsqfg,  aflim, afvec, itime)
endif


!++++++++++Print for debugging++++++++++
!print*,'ustar2',unkno(1:6,1,8)
!print*,'ustar3',ustar(:,20:50)

!...II.10.4: Face integral for RHS
if(nrz.eq.0)then
if(nquad.gt.0) call rhsifacedg_lagsubgq_gs(ipqua, unkno, ustar, ufgpq,fssgq, fsqfg, fsqfgdu, gesgq, geoel, coord,&
rhsel)
endif

!+++++++++Print for debugging++++++++++
!print*,'call ustar hybrid',rhsel(2,1,8)

!...II.10.5: Domain integral for RHS
if(nrz.eq.0)then
if(nquad.gt.0) call rhsdomndg_lagquadc_gd(intfac, ipqua, coord, coold, geoel, unkno, unkgd, strnq_devtp, rhsel,aflim,afvec, vnulq )
endif
!++++++++++Print for debugging++++++++++
!print*,'Rhs after the domian integral',itime,unkgd(:,1,8)

!
!...II.10.6: RHS for the deformation gradient...
!

!...II.10.6.1: RHS from face integral for the deformation gradient...
!!if(nquad.gt.0) call rhsifacegd_lagquadc(ipqua, ustar, gesgq0, geoel,&
!!rhsgd)
if(nquad.gt.0) call rhsifacemem_lagsms_g(ipqua, ustar, ufgpq, gesgq0, geoel, coold,rhsgd)
!print*,'face rhs for the gradient deformation',rhsgd(:,1,8)

!...II.8.6.2: RHS from domain integral for the deformation gradient...
if(nquad.gt.0) call rhsdomngd_lagquadc(intfac, ipqua, coord, coold, geoel, unkno, unkgd,rhsgd,aflim,afvec )
!print*,'Rhs after the domian integral',itime,rhsel(1:6,2,1),rhsgd(5,1,1)
!

!...II.9.7: Domain integral for source term
if(ncase .eq. 1)then
if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!...II.9.8: Update the 1st derivative
if(nlimi.eq.6)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!
dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
endif

endif !...case nfint
Case (2)!...Burton & Morgan

!...nfint.eq.1: Anylytical integration for face integral
if(nfint.eq.1)then

!...1. Get F^* N ds for all faces...
!call getfnds_lagnew_curv(gflag,gelag,intfac,inpoel,coord)
 call getfnds_lagnew_curv_hybrid(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)

!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!...Quad limiting
!!call barthlimit_lag_curvquad(geoel, coord, coold, ustar, unkno, ipqua, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...Adjust limiter
aflim(:,:) = 1.d0
!aflim(1,:) = 0.d0
afvec(:,:,:) = 0.d0
afvec(1,1,:) = 1.d0
afvec(2,2,:) = 1.d0

!...2. Get the nodal velocity and pressure with Riemann solver...
!call getndvelo_lag_curv_vilar(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fcurv)
 call getndvelo_lag_mc_curv(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold, unkno,ustar, fcurv, fcurvq, aflim, afvec, itime)

!...Pseudo curved cells---
!call getndvelo_lag_mc_curvpseudo(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
!coord, coold, unkno,ustar, fcurv, fcurvq, aflim, afvec, itime)

!++++++++++Print for debugging++++++++++
!print*,'unkno after limiter-velo_lag',aflim(1:4,1)

!...3. Face integral for RHS...
if(ntria.gt.0)  call rhsifacedg_lagmc_curvtria(inpoel,  unkno, ustar, fcurv, gelag, geoel, rhsel)
if(nquad.gt.0)  call rhsifacedg_lagmc_curvquad(ipqua,   unkno, ustar,fcurvq, gelagq, geoel,rhsel)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',unkno(1:3,3,1),unkno(1:3,2,5)
!print*,'Output RHS after sub. rhsifacedg_lag2',rhsel(1:3,3,1),rhsel(1:3,2,5)

!...4. Domain integral for RHS

if(ntria.gt.0)  call rhsdomndg_lagmc_curvtria(intfac, inpoel, coord, geoel, unkno, rhsel, aflim,afvec )
if(nquad.gt.0)  call rhsdomndg_lagmc_curvquad(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec,vnulq )

!++++++++++Print for debugging++++++++++
!print*,'unkno after limiter-domn',itime, unkno(1:3,1:4,10)

!...5. Domain integral for source term
if(ncase .eq. 1)then
if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_curvquad(intfac, ipqua,  coord, geoel, rhsel)
endif

!if(ncase .eq. 1) call rhsdomnsrcdg_lag_curv3(intfac, inpoel, coord, rhsel)

!...6.Update the 1st derivative
if(nlimi.eq.6)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!
dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)

enddo
endif
!
elseif(nfint.eq.2)then
!
print*,'This solver chocie (ncurv=1, nriem=2, nfint.eq.2) is not implemented yet !!!'
!
endif

Case (3) !...curved FVM...

!...nfint.eq.1: Analytical integration for face integral
if(nfint.eq.1)then

!...1. Get F^* N ds for all faces...
call getfnds_lagnew_curv_hybrid(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)

!...Lease-square reconstruction...
call getgraduls_lag(bface,intfac, iptri, ipqua, unkno ,cocent, coord)

!unkno(2:3,:,:) = 0.d0

!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

!...Quad limiting
call barthlimit_lagfvm_curvquad(geoel, cocent, coord, coold, ustar, unkno, ipqua, bface,intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)

!...Adjust limiter
!aflim(:,:) = 1.d0
!aflim(1,:) = 0.d0
!afvec(:,:,:) = 0.d0
!afvec(1,1,:) = 1.d0
!afvec(2,2,:) = 1.d0

!...2. Get the nodal velocity and pressure ...
call getndvelo_lagfvm_curv(gflag,gelag,cocent,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord,unkno,ustar,  fcurv, fcurvq, aflim, afvec, itime)

!...3. Face integral for RHS...
if(nquad.gt.0)  call rhsifacedg_lagfvm_curvquad(ipqua, unkno, ustar,fcurvq, gelagq, geoel,&
rhsel)

!...4. Domain integral for source term
if(ncase .eq. 1)then
 if(nquad.gt.0) call rhsdomnsrcdg_lagfvm_curvquad(intfac, ipqua, coord, geoel,rhsel)
endif

!nfint.eq.2: Simpson rule integration for face integral
elseif(nfint.eq.2)then

!...1. Get F^* N ds for all faces...
call getfnds_lag_simpsonh(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)

!...Lease-square reconstruction...
call getgraduls_lag(bface,intfac, iptri, ipqua, unkno ,cocent, coord)

!...Limiting part...
call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

!...Quad limiting
call barthlimit_lagfvm_curvquad(geoel, cocent, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec,&
unmax, unmin, esuv1, esuv2)
!call barthlimitface_fvm(geoel, cocent, coord, coold, ustar, unkno, ipqua, bface,intfac, aflim, afvec)
!aflim(:,:) = 1.d0
!aflim(1,:) = 0.d0
!afvec(:,:,:) = 0.d0
!afvec(1,1,:) = 1.d0
!afvec(2,2,:) = 1.d0

!...2. Get the nodal velocity and pressure ...
call getndvelo_lagfvm_simpson(gflag,gelag,cocent, gelagq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, fnew, fnewq, aflim, afvec,vnulq, itime)

!...3. Face integral for RHS...
call rhsifacedg_lagfvmq_simpson(ipqua,  unkno, ustar, fnewq, gelagq, geoel,&
rhsel)

!...4. Domain integral for source term
if(ncase .eq. 1)then
  if(nquad.gt.0) call rhsdomnsrcdg_lagfvm_curvquad(intfac, ipqua, coord, geoel,rhsel)
endif

endif
!
end select
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!...Part III: Cubic cell
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
elseif(ncurv==2)then
!
if(ncase.eq.3)then
unkno(2:ndegr,:,1) = 0.d0
!unkno(4:6,:,:) = 0.d0
endif
!
if(nfint.eq.1)then
!...II.4.1: Get F^* N ds for all faces...
call getfnds_lagsmsif_hybrid_cubic(gflag,gesgt,gesgq,intfac,iptri,ipqua,coord)
call getfnds_lagsmsef_cubic(gflag,gesgq,gesgt,intfac,inpoel,iptri,ipqua,coord)

!...Limiting part...
if(npoly.ge.1) call barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)

if(npoly.eq.1)then

!...Tria limiting
!call barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!...Quad limiting
call barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)

!...All the cells will be limited
geoel(10, :) = 00.d0

elseif(npoly.gt.1)then

!print*,'before shock',gesgq(1:3,8,2100)

!...Troubled-cell marker
if(rkstg==1) call getSIHO_general(intfac, ipqua, coord, coold, geoel, unkno)

!call getSIHO_vorticity(intfac, ipqua, coord, coold, geoel, unkno)

!...Remove the troubled-cells marker
!geoel(10, :) = 0.d0
!
!print*,'limweird',unkno(:,1,2100)

!...Limiting the troubled-cells
 call barthlimit_lag_curvquadcb_HO(geoel, coord, coold, ustar, unkno, ipqua,bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
!...Riemann invariant limiting
!call barthlimit_lag_riemaninv(geoel, coord,unkno, ipqua,&
!bface, intfac, esuv1, esuv2)

!call barthlimit_lag_rieminvrnt(geoel, coord,unkno, ipqua,&
!bface, intfac, esuv1, esuv2)
endif

!...Adjust limiter
!aflim(:,:) = 1.d0
!afvec(:,:,:) = 0.d0
!afvec(1,1,:) = 1.d0
!afvec(2,2,:) = 1.d0
!

!...II.4.2: Get the Riemann velocity and Pressure ...
call getRiemann_lag_sms(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, fssgt, fssgq, aflim, afvec, itime)

!++++++++++Print for debugging++++++++++
!print*,'unkno after limiter-velo_lag',itime,unkno(:,1,3)

!print*,'unkno after limiter-velo_lag2',itime,fssgq(1, 1:2, :, :, ie)
!...II.4.3: Face integral for RHS...
!if(ntria.gt.0) call rhsifacedg_lagsubgt3(iptri, unkno, ustar,fssgt, gesgt, geoel,&
!rhsel)

if(nquad.gt.0) call rhsifacedg_lag_sms2(ipqua, unkno, ustar,fssgq, gesgq, geoel,&
rhsel)

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsifacedg_lag',rhsel(1,:,1)
!
if(npoly.ge.1)then
!...II.4.4: Domain integral for RHS
!if(ntria.gt.0) call rhsdomndg_lagmc_curvtria(intfac, inpoel, coord, geoel, unkno, rhsel, aflim,afvec)

!if(nquad.gt.0)  call rhsdomndg_lagmc_qsubg(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )
if(nquad.gt.0)  call rhsdomndg_lag_quadc(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec)

endif

!++++++++++Print for debugging++++++++++
!print*,'Output RHS after sub. rhsdomndg_lag',rhsel(1,:,1)

!...II.4.5: Domain integral for source term
if(ncase .eq. 1)then
!if(ntria.gt.0) call rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
if(nquad.gt.0) call rhsdomnsrcdg_lag_quadc(intfac, ipqua, coord, geoel,rhsel)
endif

!...II.4.6: Update the 1st derivative
if(nlimi.eq.6.and.npoly.ge.1)then
do ie = 1, ncell
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!
dudr = afvec(1, 1, ie)*unkno(2,2,ie) +  afvec(1, 2, ie)*unkno(2,3,ie)
duds = afvec(1, 1, ie)*unkno(3,2,ie) +  afvec(1, 2, ie)*unkno(3,3,ie)
dvdr = afvec(2, 1, ie)*unkno(2,2,ie) +  afvec(2, 2, ie)*unkno(2,3,ie)
dvds = afvec(2, 1, ie)*unkno(3,2,ie) +  afvec(2, 2, ie)*unkno(3,3,ie)
!
unkno(2, 2, ie) = dudr; unkno(3, 2, ie) = duds;
unkno(2, 3, ie) = dvdr; unkno(3, 3, ie) = dvds;
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
enddo
endif

 endif

endif
!
deallocate(amatrv)
deallocate(gflag, gelag, lpnp, fstar, fcurv, gstar,aflim,ufgaus,ufgpt)
deallocate(gelagq, fstarq, fssgt, fssgq, cocent,vnulq)
deallocate(fstrz, fsqrz, fstrzaw, fsqrzaw)
deallocate(fstfgdu,fstfg,fsqfgdu,fsqfg,ufgpq)
!
return
end subroutine getrhsdg_lag_mcnew
!
!...Subroutine for barth limiter based on vertex for primitive variables....
!
subroutine barthlimit_lag_vtx_prim(coord, ustar, unkno, inpoel, intfac, aflim)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),  intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:nq+1, 1:nelem+nbfac),             intent(out)::aflim
integer*4,dimension(1:nvtri,1:nelem),           intent(in)::inpoel
integer, intent(in)::intfac(nifai,nafac)
!
!...Local
!
integer::ip(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvtri)::alfa
real*8:: xv(3), yv(3)
real*8:: b(3, 1:nvtri)
real*8:: unmax(1:nq+1, 1:npoin), unmin(1:nq+1, 1:npoin), dunk(1:nq+1)
real*8:: unmax_new(1:nq+1, 1:nelem), unmin_new(1:nq+1, 1:nelem)
real*8,dimension(1:nq+1,  1:nvtri) ::unknv
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct, rhom, afbar
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
!
!...shape functions
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0

do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
! if(ie==1) print*,unkno(1, 1:4, ie), unknv(1:4, 1)
!
!...Part 1: Get the maximum and minimum at the vertex...
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!...
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = unkno(1, 1, ie)
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) = pctr
unctr(nq+1) = ectr
!
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvtri
do iq = 1, nq+1
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
if(unctr(iq).gt.unmax(iq, ip(iv))) then
unmax(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
if(unctr(iq).lt.unmin(iq, ip(iv))) then
unmin(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
enddo
!
!    if(ip(iv).eq.1) print*,'unmaxxxxxxp1',ie,unmax(1:nq, 1) ,unmin(1:nq, 1)
enddo
!
enddo
!
!...Get the minimum and maximum...
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
do iq=1, nq+1
   !
unmax_new(iq, ie) = maxval(unmax(iq, ip(1:nvtri)))
unmin_new(iq, ie) = minval(unmin(iq, ip(1:nvtri)))

enddo
!
enddo
!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Part 2: Impose limiter
!
do ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
!
enddo
!
rhov = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unknv(1, iv) = 1.d0/rhov
unknv(4 ,iv) = pvtx
!
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) = pctr
!
do iv = 1, nvtri
!
do iq = 1, nq
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!                                   (unmin(iq,ip(iv)) - unctr(iq)),&
!                                   (unmin(iq,ip(iv)) - unctr(iq))/dunk(iq)
if(dunk(iq).gt.0.d0)then

fiy = .9d0*(unmax(iq, ip(iv)) - unctr(iq))/dunk(iq)
!fiy = .5d0*(unmax_new(iq, ie) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
 alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
elseif(dunk(iq).lt.0.d0)then

fiy = .9d0*(unmin(iq, ip(iv)) - unctr(iq))/dunk(iq)
!fiy = .5d0*(unmin_new(iq, ie) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
!
else
!
alfa(iq, iv) = 1.d0
!
endif
!
call barthfct(unmax(iq, ip(iv)), unmin(iq, ip(iv)), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
if(indbd(ip(iv)).eq.1)then
! alfa(1:4, iv) = 1.d0
!alfa(1, iv) = 1.d0
!alfa(4, iv) = 1.d0
if(ncase .eq. 2)then
alfa(2, iv) = min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0)
alfa(3, iv) = min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0)
!
alfa(2, iv) = max(min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0),0.d0)
alfa(3, iv) = max(min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0),0.d0)
!
elseif(ncase.eq.1)then
if(coord(1, ip(iv)).lt.1.d-6.or.abs(coord(1, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max((- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).lt.1.d-6.or.abs(coord(2, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
alfa(3, iv) = min(max(( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
elseif(ncase.eq.3)then
!
if(coord(1, ip(iv)).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max(0.5d0*(- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
alfa(3, iv) = min(max(0.5d0*( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).gt.1.d-6) then
!    alfa(4, iv) = min(max(0.5d0*(1.d-6- unctr(4))/(dunk(4)),0.d0),1.d0)
endif
!
!alfa(4, iv) = 1.d0
!
elseif(ncase.eq.7)then
!
!alfa(1, iv) =1.d0

if(coord(1, ip(iv)).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max(0.3d0*(- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
alfa(3, iv) = min(max(0.3d0*( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
endif
!
endif
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 1,nq
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
enddo
!
!aflim(1,:)= 1.d0
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
!

!
!...Part 3: Correct total energy
!
do ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) = pctr
unctr(nq+1) = ectr
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
!
enddo
!
rhov = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!unknv(1, iv) = 1.d0/rhov + aflim(1, ie)*(rhom - 1.d0/rhov)
!unknv(2, iv) = uvtx + aflim(2, ie)*(uctr - uvtx)
!unknv(3, iv) = vvtx + aflim(3, ie)*(vctr - vvtx)
!unknv(4 ,iv) = pvtx + aflim(4, ie)*(pctr - pvtx)
!
unknv(1, iv) = rhom + aflim(1, ie)*(1.d0/rhov - rhom)
unknv(2, iv) = uctr + aflim(2, ie)*(uvtx - uctr)
unknv(3, iv) = vctr + aflim(3, ie)*(vvtx - vctr)
unknv(4 ,iv) = pctr + aflim(4, ie)*(pvtx - pctr)
!
unknv(5, iv) = unknv(4 ,iv)/(gamlg-1.d0)*unknv(1, iv) + 0.5d0*(unknv(2, iv)**2 + unknv(3, iv)**2)
!
enddo
!
!
do iv = 1, nvtri
!
do iq = nq+1, nq+1
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!                                   (unmin(iq,ip(iv)) - unctr(iq)),&
!                                   (unmin(iq,ip(iv)) - unctr(iq))/dunk(iq)
if(dunk(iq).gt.0.d0)then

fiy = .9d0*(unmax(iq, ip(iv)) - unctr(iq))/dunk(iq)
!fiy = .5d0*(unmax_new(iq, ie) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
elseif(dunk(iq).lt.0.d0)then

fiy = .9d0*(unmin(iq, ip(iv)) - unctr(iq))/dunk(iq)
!fiy = .5d0*(unmin_new(iq, ie) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
!
else
!
alfa(iq, iv) = 1.d0
!
endif
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = nq+1,nq+1
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
enddo
!
!aflim(1,:)= 1.d0
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
!aflim(4, :) =1.d0
!aflim(5, :) =1.d0
!
do ie = 1, nelem
!
!unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!unkno(2:3, 2, ie) = unkno(2:3, 2, ie)*aflim(2,ie)
!unkno(2:3, 3, ie) = unkno(2:3, 3, ie)*aflim(3,ie)
!unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
!
enddo

!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lag_vtx_prim
!
!...
!
subroutine barthfctdiv(unkno,umax, umin, uctr, dunk, alfa)
implicit none
real*8, dimension(1:3, 1:4),intent(in)::unkno
real*8, intent(in)::umax, umin, uctr,dunk
real*8::alfa
real*8::fiy
real*8::coef, dive
!
dive = unkno(2, 2) + unkno(3, 3)
!
if(dive.gt.0.d0)then
coef = 0.d0
else
coef = .75d0
endif
!
!
!print*,'bad'
if(dunk.gt.1.d-8)then
!
fiy = coef*(umax - uctr)/dunk
alfa = max(min(1.d0, fiy), 0.d0)
! alfa = (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)
!
elseif(dunk.lt.-1.d-8)then
fiy = coef*(umin - uctr)/dunk
alfa = max(min(1.d0, fiy), 0.d0)
! alfa = (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)
else
alfa = 1.d0
endif
!
end subroutine barthfctdiv
!
!...Barth function...
!
subroutine barthfct(umax, umin, uctr, dunk, alfa)
use constant
implicit none
real*8, intent(in)::umax, umin, uctr,dunk
real*8::alfa
real*8::fiy
real*8::coef, cthld

if(tacru.le.0.05d0)then
coef = 0.7d0
!elseif(tacru.le.0.55d0)then
!coef = 0.45d0
else
coef = .7d0
endif
!
!coef = 1.d0
cthld = 1.75d0!1.5d0!1.5d0!2.d0!1.75d0

if(dunk.gt.1.d-8)then
!
fiy = coef*(umax - uctr)/dunk
!alfa = max(min(1.d0, fiy), 0.d0)
alfa = min((fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0),1.d0)

if(fiy.gt.cthld)then
!  alfa = 1.d0
else
!  alfa = max(min(-4.d0/27.d0*fiy**3+fiy,1.d0), 0.d0)
!  alfa = max(min(-16.d0/343.d0*fiy**3-8.d0/49.d0*fiy**2+fiy,1.d0), 0.d0)
!   alfa = max(min(-0.25d0*fiy**2+fiy,1.d0), 0.d0)
endif
elseif(dunk.lt.-1.d-8)then
fiy = coef*(umin - uctr)/dunk
!alfa = max(min(1.d0, fiy), 0.d0)
alfa = min((fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0),1.d0)

if(fiy.gt.cthld)then
!  alfa = 1.d0
else
!   alfa = max(min(-4.d0/27.d0*fiy**3+fiy,1.d0), 0.d0)
!  alfa = max(min(-16.d0/343.d0*fiy**3-8.d0/49.d0*fiy**2+fiy,1.d0), 0.d0)
!   alfa = max(min(-0.25d0*fiy**2+fiy,1.d0), 0.d0)
endif
else
alfa = 1.d0
endif
!
end subroutine barthfct
!
!...subroutine: Calculate the nodal velocity U_p^* (mass center)...
!
subroutine getndvelo_lag_mc_prim(gflag,gelag,geoel,bface,intfac,inpoel,coord,unkno,ustar, fstar, aflim, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3, 1:ngelg, 1:nelem+nbfac), intent(in)::gelag
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nelem+nbfac),     intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri, 1:nelem),  intent(out)::fstar !...Riemann forces
real*8,dimension(1:nq+1,1:nelem+nbfac),  intent(in)::aflim !...Limiter coef
integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2)     :: ipf
integer::indnd(npoin)

!...local real array
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:nvtri)::murie
real*8::vnorm(1:3, 1:2, 1:3)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8::aujmp(1:3, 1:nvtri)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:nvtri):: xv, yv
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::uvtx,vvtx,evtx, pvtx
real*8::rhom, rhomv
real*8::dux,duy,deltu
real*8::dwav1,dwav2
real*8::anx, any
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: cnsup(:), munacn(:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: usold(:,:), munacu(:,:), snsigm(:,:)
real*8,allocatable:: munacl(:,:,:)
real*8,allocatable:: snsigml(:,:,:,:), munaul(:,:,:,:)

!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:npoin), cnsup(1:npoin))
allocate (usold(1:ndimn, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munacl(1:2, 1:nvtri, 1:nelem), munaul(1:ndimn, 1:2, 1:nvtri,  1:nelem),&
snsigml(1:ndimn, 1:2,  1:nvtri,  1:nelem))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
!
!...Zero out vlave
!
cnsup = 0.d0
vlave = 0.d0
usold = ustar
indnd = 0
!
!...Mark the boundary nodes...
!
if(ncase.eq.2)then
do ifa = 1, nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
indnd(ipf(1:2)) = 1

enddo
endif
!
!
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
!
!      print*,'ifa', ifa,ipf
indnd(ipf(1:nvfac)) = 1
endif
enddo

!
do 200 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...shape functions
!
dr = .5d0
ds = .5d0
!rc = 1.d0/3.d0
!sc = rc
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0

do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...
!...zero out unknv
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
!
enddo
! if(ip(iv)==36) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
enddo
!
!...Accumulate nodal velocity...
vlave(1, ip(1:nvtri)) = vlave(1, ip(1:nvtri)) + unknv(2, 1:nvtri)
vlave(2, ip(1:nvtri)) = vlave(2, ip(1:nvtri)) + unknv(3, 1:nvtri)
!
!...Accumualte to get the no of cell surrounding one vertex...
!
cnsup(ip(1:nvtri)) = cnsup(ip(1:nvtri)) + 1.d0
!
200 enddo  !...(1)ie = 1,nelem
!
!...Get the averaged reconstructed nodal velocity...
!
do ipoin = 1, npoin
vlave(1:ndimn, ipoin) = vlave(1:ndimn, ipoin)/cnsup(ipoin)
enddo

!
!print*,'average vlave',vlave(1:2,15)
!
!
!...Zero out munacn
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
do 250 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)

!
!...shape functions
!
dr = .5d0
ds = .5d0
!rc = 1.d0/3.d0
!sc = rc
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
!...Basis function
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds

enddo
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
!
!...ndA=0.5d0*vnorm
!
vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
!...cell averaged value...
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
pctr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
!
!...
!...zero out unknv
unknv = 0.d0
!
do iv   = 1,nvtri
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
rho  = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
pvtx = unknv(4, iv)
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
!aujmp(1:2, iv) = usold(1:2, ip(iv)) - unknv(2:3, iv)
aujmp(1:2, iv) = vlave(1:2, ip(iv)) - unknv(2:3, iv)
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
rhoct = 1.d0/unkno(1, 1, ie)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ie)
vctr  = unkno(1, 3, ie)
pctr  = unkno(1, 4, ie)
!pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...Get impedence coefficient...
!
do iv   = 1, nvtri
dux= vlave(1, ip(iv))-unknv(2, iv)
duy= vlave(2, ip(iv))-unknv(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ip(iv).eq.5) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ip(iv)),unknv(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
if(aujmp(3, iv)/sdctr.lt.1.d-9)then
!
!  print*,'itime',itime
munacn(ip(iv)) = munacn(ip(iv)) + murie(iv)*vnorm(3, ifa, iv)
!
munacu(1, ip(iv)) =  munacu(1, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknv(2, iv)
munacu(2, ip(iv)) =  munacu(2, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknv(3, iv)
else
munacn(ip(iv)) = munacn(ip(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
munacu(1, ip(iv)) =  munacu(1, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(2, iv)
munacu(2, ip(iv)) =  munacu(2, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(3, iv)
!
endif
!
!   if(ip(iv).eq.36) print*,'p36 muacn(vv) post',ie,murie(iv),munacn(ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p19 muacn(28) prep---',murie(iv), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                 vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
! if(ip(iv).eq.15) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Get the summed stress sum(n*A*sigma)
!
!
!  snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 1, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!
!  snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(1, 2, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!

snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) !
snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
if(aujmp(3, iv)/sdctr.lt.1.d-9)then

munacl(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)
munacl(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)

else

munacl(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)* &
abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))

munacl(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)* &
abs(vnorm(1, 2, iv)*aujmp(1, iv) + vnorm(2, 2, iv)*aujmp(2, iv))

endif
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaul(1, 1, iv, ie)    =  munacl(1, iv, ie)*unknv(2, iv)
munaul(2, 1, iv, ie)    =  munacl(1, iv, ie)*unknv(3, iv)

munaul(1:2, 2, iv, ie)    =  munacl(2, iv, ie)*unknv(2:3, iv)
!
!    snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!    snsigml(2, 1, iv, ie)= sigma(1, 2, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
!    snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!    snsigml(2, 2, iv, ie)= sigma(1, 2, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!
snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv)

snsigml(2, 1, iv, ie)= sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv)

snsigml(2, 2, iv, ie)= sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)!
enddo
!
250 enddo  !...(1)ie = 1,nelem
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!...Impose BC with pressure prescribed boundary...
!
call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)!
!...4.1: Update the Riemann forces at every node...
!
!fpres = 0.d0
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) - fpres(1, ipoin))/munacn(ipoin)
ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) - fpres(2, ipoin))/munacn(ipoin)
endif
enddo
!
!....Bd velocity
!
! print*,'ustar--',ustar(1:2, 36),munacu(1:2,36) ,snsigm(1:2, 36), munacn(36)
!
!
if(ncase.eq.1)then
!
do 900 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
!    ustar(1, ipf(1)) = sin(pi*coord(1,ipf(1)))*cos(pi*coord(2,ipf(1)))
!    ustar(2, ipf(1)) =-cos(pi*coord(1,ipf(1)))*sin(pi*coord(2,ipf(1)))
!
!    ustar(1, ipf(2)) = sin(pi*coord(1,ipf(2)))*cos(pi*coord(2,ipf(2)))
!    ustar(2, ipf(2)) =-cos(pi*coord(1,ipf(2)))*sin(pi*coord(2,ipf(2)))
!
if(coord(1, ipf(1)).lt.1.d-6.or.abs(coord(1, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf',ipf(1)
ustar(1, ipf(1)) = 0.d0
endif
if(coord(2, ipf(1)).lt.1.d-6.or.abs(coord(2, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ipf(1)
ustar(2, ipf(1)) = 0.d0
endif
900 enddo
!
!  ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
!
!...Imposing the zero normal velocity for BC...
!
call getbcvn_lag(bface, intfac, gflag, ustar)
!
!call getbcve_exact(bface, intfac, gflag, ustar, coord, itime)
!
!   print*,'ustar',ustar(1:2, 18),ustar(1:2, 20),ustar(1:2, 21)
!  if(ustar(2, 18).gt.0.d0) then
!   print*,'Wrong movng direction'!,munacu(2, 28) - snsigm(2, 28) , fpres(2, 28),munacn(28)
!   stop
!  endif
!
!...4.2: Update the Riemann forces at every node...
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...shape functions
!
dr = .5d0
ds = .5d0
!rc = 1.d0/3.d0
!sc = rc
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
unknv = 0.d0
!
do iv = 1, nvtri
!
!...Basis function
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*ustar(1, ip(iv))- munaul(1, 1, iv, ie)
fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*ustar(2, ip(iv))- munaul(2, 1, iv, ie)
!
fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*ustar(1, ip(iv))- munaul(1, 2, iv, ie)
fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*ustar(2, ip(iv))- munaul(2, 2, iv, ie)
!
enddo
!
enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
end subroutine getndvelo_lag_mc_prim
!
subroutine rhsdomndg_lag_mc_weno(intfac, inpoel, coord, geoel, unkno, rhsel,aflim )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nelem+nbfac),     intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8, dimension(1:nq, 1:nelem+nbfac),      intent(in)::aflim
real*8,dimension(1:ndegr,1:nq,1:nelem),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv
!
!...local integer array
!
integer,dimension(1:nvtri) :: ip
!...local real array
real*8,dimension(1:ndimn, 1:nvtri) :: coorp
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc
real*8::xg, yg
real*8::rho,uadv,vadv,eadv,rhom
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr, rhomv
!
real*8,allocatable::weigh(:), posi(:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
allocate (weigh(ngausd), posi(2,ngausd))
call rutope(2, ngausd, posi, weigh)
!
!...Loop over elements
!
do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
!...Points consitituting one element...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
coorp(1, 1:nvtri) = coord(1, ip(1:nvtri))
coorp(2, 1:nvtri) = coord(2, ip(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ie) !...mass center...
sc= geoel(2, ie)
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
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
wi     = weigh(ig)
!
dspr(1) = -c10
dspr(2) =  c10
dspr(3) =  0.d0
!
dsps(1) = -c10
dsps(2) =  0.d0
dsps(3) =  c10
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*coorp(1,ishp)
dxds = dxds + dsps(ishp)*coorp(1,ishp)

dydr = dydr + dspr(ishp)*coorp(2,ishp)
dyds = dyds + dsps(ishp)*coorp(2,ishp)
enddo
!
djaco = 0.5d0*wi
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
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ie)*b(ideg)
enddo
!
!...Primitive variables...
!
rhom = unknod(1)
rho  = 1.d0/rhom
uadv = unknod(2)
vadv = unknod(3)
pres = unknod(4)
!pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
if(nlimi.eq.1)then
!
rhomc = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
rhomv = rhomc + aflim(1 ,ie)*(rhomv - rhomc)
rho = 1.d0/rhomv
!
uadv = uctr + aflim(2, ie)*(uadv - uctr)
vadv = vctr + aflim(3 ,ie)*(vadv - vctr)
!
pres = pctr + aflim(4, ie)*(pres- pctr)

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
rhsel(ideg,1:nq,ie)=rhsel(ideg,1:nq,ie) - fluxd(ideg,1:nq)*djaco
enddo
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomndg_lag_mc_weno
!
!...Subroutine for barth limiter based on vertex using compatiable gradient limiting....
!
subroutine barthlimit_lag_vtx_cgl(coord, ustar, unkno, inpoel, intfac, aflim)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),  intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:nq, 1:nelem+nbfac),         intent(out)::aflim
integer*4,dimension(1:nvtri,1:nelem),           intent(in)::inpoel
integer, intent(in)::intfac(nifai,nafac)
!
!...Local
!
integer::ip(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa
real*8:: unctr(1:nq)
real*8,  dimension(1:nq, 1:nvtri)::alfa
real*8:: xv(3), yv(3)
real*8:: b(3, 1:nvtri)
real*8:: unmax(1:nq, 1:npoin), unmin(1:nq, 1:npoin), dunk(1:nq)
real*8,dimension(1:nq,  1:nvtri) ::unknv
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: eps
real*8:: rhomc, uctr, vctr, ectr, pctr
real*8:: drhomdr, drhomds, dudr, duds, dvdr, dvds, dedr, deds, dpdr, dpds
real*8:: der1, der2, dpr1, dpr2
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct,rhom
real*8:: drhodr, drhods
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
!
!...shape functions
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0

do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
! if(ie==1) print*,unkno(1, 1:4, ie), unknv(1:4, 1)
!
!...Part 1: Get the maximum and minimum at the vertex...
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!...
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = unkno(1, 1, ie)
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) =  ectr
!
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvtri
do iq = 1, nq
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
if(unctr(iq).gt.unmax(iq, ip(iv))) then
unmax(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
if(unctr(iq).lt.unmin(iq, ip(iv))) then
unmin(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
enddo
!
!    if(ip(iv).eq.1) print*,'unmaxxxxxxp1',ie,unmax(1:nq, 1) ,unmin(1:nq, 1)
enddo
!
enddo

!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Part 2: Impose limiter
!
do ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
!
enddo
!
rhov = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unknv(1, iv) = 1.d0/rhov
unknv(4 ,iv) = evtx
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) = ectr
!
do iv = 1, nvtri
!
do iq = 1, nq
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!                                   (unmin(iq,ip(iv)) - unctr(iq)),&
!                                   (unmin(iq,ip(iv)) - unctr(iq))/dunk(iq)
if(dunk(iq).gt.0.d0)then

fiy = .5d0*(unmax(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
elseif(dunk(iq).lt.0.d0)then

fiy = .5d0*(unmin(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
!
else
!
alfa(iq, iv) = 1.d0
!
endif
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
if(indbd(ip(iv)).eq.1)then
! alfa(1:4, iv) = 1.d0
!alfa(1, iv) = 1.d0
!alfa(4, iv) = 1.d0
if(ncase .eq. 2)then
alfa(2, iv) = min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0)
alfa(3, iv) = min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0)
!
alfa(2, iv) = max(min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0),0.d0)
alfa(3, iv) = max(min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0),0.d0)
!
elseif(ncase.eq.1)then
if(coord(1, ip(iv)).lt.1.d-6.or.abs(coord(1, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max((- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).lt.1.d-6.or.abs(coord(2, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
alfa(3, iv) = min(max(( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
elseif(ncase.eq.3)then
!
if(coord(1, ip(iv)).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max(0.5d0*(- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
alfa(3, iv) = min(max(0.5d0*( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).gt.1.d-6) then
!    alfa(4, iv) = min(max(0.5d0*(1.d-6- unctr(4))/(dunk(4)),0.d0),1.d0)
endif
!
endif
!
endif
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 1,nq
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
enddo
!
!aflim(1,:)= 1.d0
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
!
!
do ie = 1, nelem
!
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
unkno(2:3, 2, ie) = unkno(2:3, 2, ie)*aflim(2,ie)
unkno(2:3, 3, ie) = unkno(2:3, 3, ie)*aflim(3,ie)
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(4,ie)
!
enddo

!
!...Part 1: Get the maximum and minimum at the vertex...
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!...limite pressure....
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = unkno(1, 1, ie)
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) =  pctr
!
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvtri
do iq = 1, nq
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
if(unctr(iq).gt.unmax(iq, ip(iv))) then
unmax(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
if(unctr(iq).lt.unmin(iq, ip(iv))) then
unmin(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
enddo
!
!    if(ip(iv).eq.1) print*,'unmaxxxxxxp1',ie,unmax(1:nq, 1) ,unmin(1:nq, 1)
enddo
!
enddo

!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Part 2: Impose limiter
!
do ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
!
enddo
!
rhov = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unknv(1, iv) = 1.d0/rhov
unknv(4 ,iv) = pvtx
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) = pctr
!
do iv = 1, nvtri
!
do iq = 4, 4
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!                                   (unmin(iq,ip(iv)) - unctr(iq)),&
!                                   (unmin(iq,ip(iv)) - unctr(iq))/dunk(iq)
if(dunk(iq).gt.0.d0)then

fiy = .5d0*(unmax(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
elseif(dunk(iq).lt.0.d0)then

fiy = .5d0*(unmin(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
!
else
!
alfa(iq, iv) = 1.d0
!
endif
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!

!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 4,4
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
enddo
!
!...agian,limit pressure...
!
rhomc   = unkno(1, 1, ie)
drhomdr = unkno(2, 1, ie)
drhomds = unkno(3, 1, ie)
!
uctr = unkno(1, 2, ie)
dudr = unkno(2, 2, ie)
duds = unkno(3, 2, ie)
!
vctr = unkno(1, 3, ie)
dvdr = unkno(2, 3, ie)
dvds = unkno(3, 3, ie)
!
ectr = unkno(1, 4, ie)
dedr = unkno(2, 4, ie)
deds = unkno(3, 4, ie)
!
rhoct  = 1.d0/rhomc
pctr =  (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2))
drhodr = -drhomdr/rhomc**2
drhods = -drhomds/rhomc**2
!
dpr1 = ectr - 0.5d0*(uctr**2 + vctr**2)
dpr2 = dedr - (uctr*dudr + vctr*dvdr)
dpdr = (gamlg-1.d0)*(drhodr*dpr1 + rhoct*dpr2)
!
dpr1 = ectr - 0.5d0*(uctr**2 + vctr**2)
dpr2 = deds - (uctr*duds + vctr*dvds)
dpds = (gamlg-1.d0)*(drhods*dpr1 + rhoct*dpr2)
!
dpdr = dpdr*aflim(4, ie)
dpds = dpds*aflim(4, ie)
!
der1 = dpdr*rhomc + pctr*drhomdr
der2 = (uctr*dudr + vctr*dvdr)
!
dedr = der1/(gamlg-1.d0) + der2
!
der1 = dpds*rhomc + pctr*drhomds
der2 = (uctr*duds + vctr*dvds)
!
deds = der1/(gamlg-1.d0) + der2
!
unkno( 2, 4, ie)= dedr
unkno( 3, 4, ie)= deds
!
!aflim(1,:)= 1.d0
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lag_vtx_cgl
!
!...Subroutine for barth limiter based on vertex....
!
subroutine barthlimit_lag_vtx(coord, ustar, unkno, inpoel, intfac)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
integer*4,dimension(1:nvtri,1:nelem),           intent(in)::inpoel
integer, intent(in)::intfac(nifai,nafac)
!
!...Local
!
integer::ip(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa
real*8:: aflim(1:nq, 1:nelem), alfal(1:nelem)
real*8:: unctr(1:nq)
real*8,  dimension(1:nq, 1:nvtri)::alfa
real*8:: xv(3), yv(3)
real*8:: b(3, 1:nvtri)
real*8,  dimension(1:nq, 1:npoin)::unvtx
real*8:: unmax(1:nq, 1:npoin), unmin(1:nq, 1:npoin), dunk(1:nq)
real*8,dimension(1:nq,  1:nvtri) ::unknv
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
!
!...shape functions
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0

do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
! if(ie==1) print*,unkno(1, 1:4, ie), unknv(1:4, 1)
!
!...Part 1: Get the maximum and minimum at the vertex...
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!do ipoin = 1, npoin
!do ie =1 ,nelem
!
!ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!do iv =1 , nvtri
!if(ip(iv).eq.ipoin)then
!do iq = 1, nq
!if(unkno(1,iq,ie).gt.unmax(iq, ipoin)) unmax(iq, ipoin) = unkno(1,iq,ie)
!if(unkno(1,iq,ie).lt.unmin(iq, ipoin)) unmin(iq, ipoin) = unkno(1,iq,ie)
!enddo
!endif
!enddo
!
!enddo
!enddo

do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
unctr(1:nq) = unkno(1, 1:nq, ie)
!
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvtri
do iq = 1, nq
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
if(unctr(iq).gt.unmax(iq, ip(iv))) then
unmax(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
if(unctr(iq).lt.unmin(iq, ip(iv))) then
unmin(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
enddo
!
!    if(ip(iv).eq.1) print*,'unmaxxxxxxp1',ie,unmax(1:nq, 1) ,unmin(1:nq, 1)
enddo
!
enddo

!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Part 2: Impose limiter
!
do ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
unctr(1:nq) = unkno(1, 1:nq, ie)
!
do iv = 1, nvtri
!
do iq = 1, nq
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!                                   (unmin(iq,ip(iv)) - unctr(iq)),&
!                                   (unmin(iq,ip(iv)) - unctr(iq))/dunk(iq)
if(dunk(iq).gt.0.d0)then

fiy = .5d0*(unmax(iq, ip(iv)) - unctr(iq))/dunk(iq)
alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
!alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
elseif(dunk(iq).lt.0.d0)then

fiy = .5d0*(unmin(iq, ip(iv)) - unctr(iq))/dunk(iq)
alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
!alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
!
else
!
alfa(iq, iv) = 1.d0
!
endif
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
if(indbd(ip(iv)).eq.1)then
! alfa(1:4, iv) = 1.d0
!alfa(1, iv) = 1.d0
!alfa(4, iv) = 1.d0
if(ncase .eq. 2)then
alfa(2, iv) = min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0)
alfa(3, iv) = min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0)
!
alfa(2, iv) = max(min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0),0.d0)
alfa(3, iv) = max(min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0),0.d0)
!
elseif(ncase.eq.1)then
if(coord(1, ip(iv)).lt.1.d-6.or.abs(coord(1, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max((- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).lt.1.d-6.or.abs(coord(2, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
alfa(3, iv) = min(max(( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
endif
!
endif
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 1,nq
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
enddo
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
! aflim = 0.d0
!
!print*,'aflim',unkno(1:3, 3, 18),unkno(1:3, 3, 20),aflim(3,20)
!
do ie = 1, nelem
!
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
unkno(2:3, 2, ie) = unkno(2:3, 2, ie)*aflim(2,ie)
unkno(2:3, 3, ie) = unkno(2:3, 3, ie)*aflim(3,ie)
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(4,ie)
!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lag_vtx
!
!...Subroutine for barth limiter based on vertex using compatiable gradient limiting....
!
subroutine barthlimit_lag_vtx_cgl2(coord, ustar, unkno, inpoel, intfac, aflim)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),  intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:nq, 1:nelem+nbfac),         intent(out)::aflim
integer*4,dimension(1:nvtri,1:nelem),           intent(in)::inpoel
integer, intent(in)::intfac(nifai,nafac)
!
!...Local
!
integer::ip(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa
real*8:: unctr(1:nq)
real*8,  dimension(1:nq, 1:nvtri)::alfa
real*8:: xv(3), yv(3)
real*8:: b(3, 1:nvtri)
real*8:: unmax(1:nq, 1:npoin), unmin(1:nq, 1:npoin), dunk(1:nq)
real*8,dimension(1:nq,  1:nvtri) ::unknv
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: eps
real*8:: rhomc, uctr, vctr, ectr, pctr
real*8:: drhomdr, drhomds, dudr, duds, dvdr, dvds, dedr, deds, dpdr, dpds
real*8:: der1, der2, dpr1, dpr2
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct,rhom
real*8:: drhodr, drhods
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
!
!...shape functions
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0

do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
! if(ie==1) print*,unkno(1, 1:4, ie), unknv(1:4, 1)
!
!...Part 1: Get the maximum and minimum at the vertex...
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!...
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = unkno(1, 1, ie)
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) =  ectr
!
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvtri
do iq = 1, nq
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
if(unctr(iq).gt.unmax(iq, ip(iv))) then
unmax(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
if(unctr(iq).lt.unmin(iq, ip(iv))) then
unmin(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
enddo
!
!    if(ip(iv).eq.1) print*,'unmaxxxxxxp1',ie,unmax(1:nq, 1) ,unmin(1:nq, 1)
enddo
!
enddo

!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Part 2: Impose limiter
!
do ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
!
enddo
!
rhov = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unknv(1, iv) = 1.d0/rhov
unknv(4 ,iv) = evtx
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) = ectr
!
do iv = 1, nvtri
!
do iq = 1, nq
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!                                   (unmin(iq,ip(iv)) - unctr(iq)),&
!                                   (unmin(iq,ip(iv)) - unctr(iq))/dunk(iq)
if(dunk(iq).gt.0.d0)then

fiy = .5d0*(unmax(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
elseif(dunk(iq).lt.0.d0)then

fiy = .5d0*(unmin(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
!
else
!
alfa(iq, iv) = 1.d0
!
endif
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
if(indbd(ip(iv)).eq.1)then
! alfa(1:4, iv) = 1.d0
!alfa(1, iv) = 1.d0
!alfa(4, iv) = 1.d0
if(ncase .eq. 2)then
alfa(2, iv) = min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0)
alfa(3, iv) = min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0)
!
alfa(2, iv) = max(min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0),0.d0)
alfa(3, iv) = max(min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0),0.d0)
!
elseif(ncase.eq.1)then
if(coord(1, ip(iv)).lt.1.d-6.or.abs(coord(1, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max((- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).lt.1.d-6.or.abs(coord(2, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
alfa(3, iv) = min(max(( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
elseif(ncase.eq.3)then
!
if(coord(1, ip(iv)).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max(0.5d0*(- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
alfa(3, iv) = min(max(0.5d0*( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
if(coord(2, ip(iv)).gt.1.d-6) then
!    alfa(4, iv) = min(max(0.5d0*(1.d-6- unctr(4))/(dunk(4)),0.d0),1.d0)
endif
!
endif
!
endif
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 1,nq
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
enddo
!
!aflim(1,:)= 1.d0
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
!
!
do ie = 1, nelem
!
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
unkno(2:3, 2, ie) = unkno(2:3, 2, ie)*aflim(2,ie)
unkno(2:3, 3, ie) = unkno(2:3, 3, ie)*aflim(3,ie)
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(4,ie)
!
enddo

!
!...Part 1: Get the maximum and minimum at the vertex...
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!...limite pressure....
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = unkno(1, 1, ie)
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) =  pctr
!
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvtri
do iq = 1, nq
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
if(unctr(iq).gt.unmax(iq, ip(iv))) then
unmax(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
if(unctr(iq).lt.unmin(iq, ip(iv))) then
unmin(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
enddo
!
!    if(ip(iv).eq.1) print*,'unmaxxxxxxp1',ie,unmax(1:nq, 1) ,unmin(1:nq, 1)
enddo
!
enddo

!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Part 2: Impose limiter
!
do ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
!
enddo
!
rhov = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unknv(1, iv) = 1.d0/rhov
unknv(4 ,iv) = pvtx
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) = pctr
!
do iv = 1, nvtri
!
do iq = 4, 4
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!                                   (unmin(iq,ip(iv)) - unctr(iq)),&
!                                   (unmin(iq,ip(iv)) - unctr(iq))/dunk(iq)
if(dunk(iq).gt.0.d0)then

fiy = .5d0*(unmax(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
elseif(dunk(iq).lt.0.d0)then

fiy = .5d0*(unmin(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
!
else
!
alfa(iq, iv) = 1.d0
!
endif
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!

!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 4,4
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
enddo
!
!...again,limit density...
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ie)
uctr = unkno(1, 2, ie)
vctr = unkno(1, 3, ie)
ectr = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ie)
unctr(nq) = pctr
!
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
!
enddo
!
rhov = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
pvtx = pctr + aflim(4, ie)*(pvtx - pctr)
!
unknv(1, iv) = 1.d0/rhov
unknv(4 ,iv) = pvtx
enddo

!
do iv = 1, nvtri
!
do iq = 4, 4
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!                                   (unmin(iq,ip(iv)) - unctr(iq)),&
!                                   (unmin(iq,ip(iv)) - unctr(iq))/dunk(iq)
if(dunk(iq).gt.0.d0)then

fiy = .5d0*(unmax(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
elseif(dunk(iq).lt.0.d0)then

fiy = .5d0*(unmin(iq, ip(iv)) - unctr(iq))/dunk(iq)
!alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
alfa(iq, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
!
else
!
alfa(iq, iv) = 1.d0
!
endif
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!

!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 4,4
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
enddo
!
!aflim(1,:)= 1.d0
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lag_vtx_cgl2
!
!...subroutine: Calculate the averaged velocity for hybrid grids...
!
subroutine getvlavenew(iptri, ipqua, geoel, vlave, unkno, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
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
real*8::rc, sc, dr, ds
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
!...triangle...
!
do 200 ie = 1,ntria !...(1)ie = 1,nelem
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
!...shape functions
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0

do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...
!...zero out unknv
unknvt = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
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
unknvt(2, iv) = unkno(1,2,ielem)  + dudr*bt(2, iv) + duds*bt(3, iv)
unknvt(3, iv) = unkno(1,3,ielem)  + dvdr*bt(2, iv) + dvds*bt(3, iv)
!unknvt(2, iv) = unkno(1, 2, ielem) + aflim(2, ielem)*(unknvt(2, iv) - unkno(1, 2, ielem))
!unknvt(3, iv) = unkno(1, 3, ielem) + aflim(3, ielem)*(unknvt(3, iv) - unkno(1, 3, ielem))
endif
! if(ip(iv)==36) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
enddo
!
!...Accumulate nodal velocity...
vlave(1, ipt(1:nvtri)) = vlave(1, ipt(1:nvtri)) + unknvt(2, 1:nvtri)
vlave(2, ipt(1:nvtri)) = vlave(2, ipt(1:nvtri)) + unknvt(3, 1:nvtri)
!
!...Accumualte to get the no of cell surrounding one vertex...
!
cnsup(ipt(1:nvtri)) = cnsup(ipt(1:nvtri)) + 1.d0
!
200 enddo  !...(1)ie = 1,nelem
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

do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...
!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
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
!unknvq(2, iv) = unkno(1, 2, ielem) + aflim(2, ielem)*(unknvq(2, iv) - unkno(1, 2, ielem))
!unknvq(3, iv) = unkno(1, 3, ielem) + aflim(3, ielem)*(unknvq(3, iv) - unkno(1, 3, ielem))
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
end subroutine getvlavenew
!
!...subroutine: Calculate the averaged velocity for hybrid grids...
!
subroutine getvlave(iptri, ipqua, geoel, vlave, unkno, aflim)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(out)::vlave
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
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
real*8::rc, sc, dr, ds
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
!...triangle...
!
do 200 ie = 1,ntria !...(1)ie = 1,nelem
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
!...shape functions
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0

do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...
!...zero out unknv
unknvt = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
!
enddo
!
if(nlimi.eq.1)then
!
!unknvt(2, iv) = unkno(1, 2, ielem) + aflim(2, ielem)*(unknvt(2, iv) - unkno(1, 2, ielem))
!unknvt(3, iv) = unkno(1, 3, ielem) + aflim(3, ielem)*(unknvt(3, iv) - unkno(1, 3, ielem))
endif
! if(ip(iv)==36) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
enddo
!
!...Accumulate nodal velocity...
vlave(1, ipt(1:nvtri)) = vlave(1, ipt(1:nvtri)) + unknvt(2, 1:nvtri)
vlave(2, ipt(1:nvtri)) = vlave(2, ipt(1:nvtri)) + unknvt(3, 1:nvtri)
!
!...Accumualte to get the no of cell surrounding one vertex...
!
cnsup(ipt(1:nvtri)) = cnsup(ipt(1:nvtri)) + 1.d0
!
200 enddo  !...(1)ie = 1,nelem
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

do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...
!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
if(nlimi.eq.1)then
!
!unknvq(2, iv) = unkno(1, 2, ielem) + aflim(2, ielem)*(unknvq(2, iv) - unkno(1, 2, ielem))
!unknvq(3, iv) = unkno(1, 3, ielem) + aflim(3, ielem)*(unknvq(3, iv) - unkno(1, 3, ielem))
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
end subroutine getvlave

!
!...subroutine: Calculate the nodal velocity U_p^* (mass center) for hybrid mesh...
!
subroutine getndvelo_lag_mc_hybrid(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
                                   coord,coold, unkno,ustar, fstar, fstarq, aflim, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3, 1:ngelg, 1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coold !...initial coordinates...
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri, 1:ntria),  intent(out)::fstar !...Riemann forces
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)

!...local real array
real*8::fcx, fcy,mct
real*8,dimension(1:ndimn,1:npoin)::vlave
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: munacn(:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: usold(:,:), munacu(:,:), snsigm(:,:)
real*8,allocatable:: munaclt(:,:,:),munaclq(:,:,:)
real*8,allocatable:: snsigmlt(:,:,:,:), munault(:,:,:,:)
real*8,allocatable:: snsigmlq(:,:,:,:), munaulq(:,:,:,:)

!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:npoin))
allocate (usold(1:ndimn, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munaclt(1:2, 1:nvtri, 1:ntria), munault(1:ndimn, 1:2, 1:nvtri,  1:ntria),&
snsigmlt(1:ndimn, 1:2,  1:nvtri,  1:ntria))
allocate (munaclq(1:2, 1:nvqua, 1:nquad), munaulq(1:ndimn, 1:2, 1:nvqua,  1:nquad),&
snsigmlq(1:ndimn, 1:2,  1:nvqua,  1:nquad))
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
!
ipf(1:2) = intfac(3:4, ifa)
!
indnd(ipf(1:2)) = 1

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
 call getvlave(iptri, ipqua, geoel, vlave, unkno, aflim)
!
!...Specify the averaged velocity with BC
!
do 950 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
if(bface(3,ifa).eq.22)then
ipf(1:2) = intfac(3:4, ifa)
if(bface(4,ifa).eq.221)then
vlave(2,ipf(1:2)) = 0.d0
elseif(bface(4,ifa).eq.222)then
vlave(1,ipf(1:2)) = 0.d0
endif
elseif(bface(3,ifa).eq.25)then
vlave(2,ipf(1:2)) = 0.d0
vlave(1,ipf(1:2)) = 0.d0
endif!
950 enddo
!
!print*,'average vlave',vlave(1:2,15)
!
!
!...Zero out munacn
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
!call getriem(iptri, ipqua, geoel, gelag, gelagq, vlave, unkno, munacn, munacu, snsigm,&
!             munaclt, munault, snsigmlt, munaclq, munaulq, snsigmlq)
if(ntria.gt.0) call getriem_tria(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold,aflim)
if(nquad.gt.0) call getriem_quad(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold,aflim)
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)!
!
!...Update the nodal velocity...
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) - fpres(1, ipoin))/munacn(ipoin)
ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) - fpres(2, ipoin))/munacn(ipoin)
endif
enddo
!
!....Bd velocity
!
! print*,'ustar--',ustar(1:2, 7),munacu(1:2,7) ,snsigm(1:2, 7), munacn(7),sqrt(unkno(1,2,5)**2 + unkno(1,3,5)**2)
!print*,'ustar--',cos(pi/24.d0*5.d0),sin(pi/24.d0*5.d0),ustar(1:2,7)/sqrt(ustar(1,7)**2 + ustar(2, 7)**2)
!print*,'dustar',(ustar(1:2,7)-unkno(1,2:3,6))/sqrt((ustar(1,7)-unkno(1,2,6))**2+(ustar(2,7)-unkno(1,3,6))**2)
!print*,'dustar',sqrt((ustar(1,6)-unkno(1,2,5))**2+(ustar(2,6)-unkno(1,3,5))**2),&
!                sqrt((ustar(1,7)-unkno(1,2,5))**2+(ustar(2,7)-unkno(1,3,5))**2),&!
!sqrt((ustar(1,7)-unkno(1,2,6))**2+(ustar(2,7)-unkno(1,3,6))**2),&
!sqrt((ustar(1,8)-unkno(1,2,6))**2+(ustar(2,8)-unkno(1,3,6))**2)

!
do 900 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
 if(bface(3,ifa).eq.22)then
   ipf(1:2) = intfac(3:4, ifa)
   if(bface(4,ifa).eq.221)then
     ustar(2,ipf(1:2)) = 0.d0
   elseif(bface(4,ifa).eq.222)then
     ustar(1,ipf(1:2)) = 0.d0
   endif
elseif(bface(3,ifa).eq.25)then
 ustar(2,ipf(1:2)) = 0.d0
 ustar(1,ipf(1:2)) = 0.d0

 endif
!
!...impose exact solution along the Boundary
!
!    ustar(1, ipf(1)) = sin(pi*coord(1,ipf(1)))*cos(pi*coord(2,ipf(1)))
!    ustar(2, ipf(1)) =-cos(pi*coord(1,ipf(1)))*sin(pi*coord(2,ipf(1)))
!
!    ustar(1, ipf(2)) = sin(pi*coord(1,ipf(2)))*cos(pi*coord(2,ipf(2)))
!    ustar(2, ipf(2)) =-cos(pi*coord(1,ipf(2)))*sin(pi*coord(2,ipf(2)))
!
!do iv =1, nvfac
! if(coord(1, ipf(iv)).lt.1.d-6)then!.or.abs(coord(1, ipf(iv))-1.d0).lt.1.d-6) then
!   ustar(1, ipf(iv)) = 0.d0
! endif
! if(coord(2, ipf(iv)).lt.1.d-6)then!.or.abs(coord(2, ipf(iv))-1.d0).lt.1.d-6) then
!   ustar(2, ipf(iv)) = 0.d0
! endif
!enddo
!
900 enddo
!
!  ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
!
!...Imposing the zero normal velocity for BC...
!
! call getbcvn_lag(bface, intfac, gflag, ustar)
!
!call getbcve_exact(bface, intfac, gflag, ustar, coord, itime)
!
!   print*,'ustar',ustar(1:2, 18),ustar(1:2, 20),ustar(1:2, 21)
!  if(ustar(2, 18).gt.0.d0) then
!   print*,'Wrong movng direction'!,munacu(2, 28) - snsigm(2, 28) , fpres(2, 28),munacn(28)
!   stop
!  endif
!
!...4.2: Update the Riemann forces at every node...
!
do ie = 1, ntria
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
!...shape functions
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
fcx=0.d0
fcy = 0.d0
mct = 0.d0
!
do iv = 1, nvtri
!
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!

fstar(1, 1, iv, ie) = snsigmlt(1, 1, iv, ie) + munaclt(1, iv, ie)*ustar(1, ipt(iv))- munault(1, 1, iv, ie)
fstar(2, 1, iv, ie) = snsigmlt(2, 1, iv, ie) + munaclt(1, iv, ie)*ustar(2, ipt(iv))- munault(2, 1, iv, ie)
!
fstar(1, 2, iv, ie) = snsigmlt(1, 2, iv, ie) + munaclt(2, iv, ie)*ustar(1, ipt(iv))- munault(1, 2, iv, ie)
fstar(2, 2, iv, ie) = snsigmlt(2, 2, iv, ie) + munaclt(2, iv, ie)*ustar(2, ipt(iv))- munault(2, 2, iv, ie)

!fstar(1, 1, iv, ie) = munaclt(1, iv, ie)*ustar(1, ipt(iv))- munault(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = munaclt(1, iv, ie)*ustar(2, ipt(iv))- munault(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = munaclt(2, iv, ie)*ustar(1, ipt(iv))- munault(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = munaclt(2, iv, ie)*ustar(2, ipt(iv))- munault(2, 2, iv, ie)
!
! if(ie.eq.5.or.ie.eq.6) then
!  fcx = fstar(1, 1, iv, ie) + fstar(1, 2, iv, ie)
!  fcy = fstar(2, 1, iv, ie) + fstar(2, 2, iv, ie)
!  mct = munaclt(1, iv, ie) + munaclt(2, iv, ie)
! endif
!
!if(ie.eq.6) fcx = fcx+fstar(1, 1, iv, ie)*ustar(1,ipt(iv)) + fstar(2, 1, iv, ie)*ustar(2,ipt(iv))+&
!                  fstar(1, 2, iv, ie)*ustar(1,ipt(iv)) + fstar(2, 2, iv, ie)*ustar(2,ipt(iv))!+ fstar(1, 2, iv, ie)
!if(ie.eq.5) fcy = fcy+fstar(2, 1, iv, ie) + fstar(2, 2, iv, ie)
!
!if(ie.eq.5.or.ie.eq.6) print*,ipt(iv),iv, ie,fstar(1:2, 1, iv, ie),fstar(1:2, 2, iv, ie),&
!                              fcx/sqrt(fcx**2+fcy**2),fcy/sqrt(fcx**2+fcy**2),sqrt(fcx**2+fcy**2), mct,&
!                              unkno(1,2:3,ie)/sqrt(unkno(1,2,ie)**2 + unkno(1,3,ie)**2)!,&
!                              cos(4.5d0/24.d0*pi), sin(4.5d0/24.d0*pi)!,&
!                              fcx*unkno(1,2,ie)+fcy*unkno(1,3,ie)+fcx**2+fcy**2
!
enddo
enddo
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
do iv = 1, nvqua
!
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
fstarq(1, 1, iv, ie) = snsigmlq(1, 1, iv, ie) + munaclq(1, iv, ie)*ustar(1, ipq(iv))- munaulq(1, 1, iv, ie)
fstarq(2, 1, iv, ie) = snsigmlq(2, 1, iv, ie) + munaclq(1, iv, ie)*ustar(2, ipq(iv))- munaulq(2, 1, iv, ie)
!
fstarq(1, 2, iv, ie) = snsigmlq(1, 2, iv, ie) + munaclq(2, iv, ie)*ustar(1, ipq(iv))- munaulq(1, 2, iv, ie)
fstarq(2, 2, iv, ie) = snsigmlq(2, 2, iv, ie) + munaclq(2, iv, ie)*ustar(2, ipq(iv))- munaulq(2, 2, iv, ie)
!
!if(ie==21) print*,'fire',iv,ie,snsigmlq(1,1:2,iv,ie),munaclq(1:2,iv,ie),ipq(:),munaulq(1,1:2,iv,ie)
!if(ie.eq.21) print*,'fstarq',iv,ie,fstarq(1,1:2,iv,ie)
!
enddo
!
enddo
!
deallocate (munacn, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lag_mc_hybrid
!
!...subroutine: Calculate the Riemann input for hybrid triangle grids...
!
subroutine getriem_tria(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold, aflim)
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
!
real*8, dimension(1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
!
real*8, dimension(1:ndimn, 1:nvtri, 1:ntria),       intent(out)::munaclt
real*8, dimension(1:ndimn, 1:2, 1:nvtri,  1:ntria), intent(out)::munault
real*8, dimension(1:ndimn, 1:2, 1:nvtri,  1:ntria), intent(out)::snsigmlt
!...Local integer
integer::ie, ideg, ielem, ifa, iv
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:3, 1:nvtri)::bt
real*8,dimension(1:nq,1:nvtri)::unknvt
real*8::aujmp(1:3, 1:nvtri)
real*8::vnorm(1:3, 1:2, 1:nvtri)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8,dimension(1:nvtri)::murie
real*8,dimension(1:nvtri):: xv,  yv
real*8,dimension(1:ndimn, 1:nvtri) :: xpt
real*8,dimension(1:ndimn, 1:nvtri) :: xpti
real*8, dimension(1:nvtri):: shp, dspr, dsps
!
real*8::eps,c00,c05,c10,c20
real*8::rho, rhom, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s
real*8::acnx, acny
real*8::rhoi, rhon
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Triangle...
!
do 250 ie = 1,ntria !...(1)ie = 1,nelem
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
!...shape functions
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
!...Basis function
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds

enddo
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
!
!...ndA=0.5d0*vnorm
!
vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhom = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
!...now configuration
!
r = rc; s= sc
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
call  getrhoig_tria(rhoi, r, s, xpti)
!
call getdensity_triallnl(r, s, xpt, xpti, rhoi, rhon)
!
rhom = 1.d0/rhon
!
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...
!...zero out unknv
unknvt = 0.d0
!
do iv   = 1,nvtri
!
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
enddo
!
if(ndens.eq.1)then
rho  = 1.d0/unknvt(1, iv)
elseif(ndens.eq.2)then
!...now configuration
!
r = xv(iv); s= yv(iv)
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
call  getrhoig_tria(rhoi, r, s, xpti)
!
call getdensity_triallnl(r, s, xpt, xpti, rhoi, rhon)
!
rho = rhon
!
endif
!
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
evtx = unknvt(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!...Limiter
!
if(nlimi.eq.1)then
rhomv = rhom + aflim(1, ielem)*(unknvt(1, iv) - rhom)
rho = 1.d0/rhomv
!
uvtx = uctr + aflim(2, ielem)*(unknvt(2, iv) - uctr)
vvtx = vctr + aflim(3, ielem)*(unknvt(3, iv) - vctr)
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!
!...updtae unknv(2:3,:)
unknvt(2, iv) = uvtx
unknvt(3 ,iv) = vvtx
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
!
aujmp(1:2, iv) = vlave(1:2, ipt(iv)) - unknvt(2:3, iv)
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
!if(ipt(iv)==7) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, iv),sqrt(acnx**2 + acny**2), unknvt(2:3, iv)&
!,cos(5.d0/24.d0*pi),sin(5.d0/24.d0*pi),vlave(1:2, ipt(iv))/sqrt(vlave(1, ipt(iv))**2 + vlave(2, ipt(iv))**2)
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
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
!if(ielem.eq.5.or.ielem.eq.6) print*,'ielem56',unkno(1, 2:3, ielem),ielem
!
!...Get impedence coefficient...
!
do iv   = 1, nvtri
dux= vlave(1, ipt(iv))-unknvt(2, iv)
duy= vlave(2, ipt(iv))-unknvt(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!if(ipt(iv).eq.1.and.(ie.eq.5.or.ie.eq.6)) print*,'murie22', sdctr,ectr,pctr,sqrt(uctr**2 + vctr**2),rhoct,deltu,&
!              vlave(1:2, ipt(iv)),unknvt(2:3, iv),&
!      dux**2,duy**2,ielem

enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
if(aujmp(3, iv)/sdctr.lt.1.d-6)then
!
!  print*,'itime',itime
munacn(ipt(iv)) = munacn(ipt(iv)) + murie(iv)*vnorm(3, ifa, iv)
!
munacu(1, ipt(iv)) =  munacu(1, ipt(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknvt(2, iv)
!
munacu(2, ipt(iv)) =  munacu(2, ipt(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknvt(3, iv)
else
munacn(ipt(iv)) = munacn(ipt(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
munacu(1, ipt(iv)) =  munacu(1, ipt(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknvt(2, iv)
munacu(2, ipt(iv)) =  munacu(2, ipt(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknvt(3, iv)
!
endif
!
!   if(ipt(iv).eq.7) print*,'p36 muacn(vv) post',ie,murie(iv),munacn(ipt(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:3, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p19 muacn(28) prep---',murie(iv), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                 vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
! if(ip(iv).eq.15) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Get the summed stress sum(n*A*sigma)
!
!
!  snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 1, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!
!  snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(1, 2, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!

snsigm(1, ipt(iv)) = snsigm(1, ipt(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) !
snsigm(2, ipt(iv)) = snsigm(2, ipt(iv)) + sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
if(aujmp(3, iv)/sdctr.lt.1.d-6)then

munaclt(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)
munaclt(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)

else

munaclt(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)* &
abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))

munaclt(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)* &
abs(vnorm(1, 2, iv)*aujmp(1, iv) + vnorm(2, 2, iv)*aujmp(2, iv))

endif
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munault(1, 1, iv, ie)    =  munaclt(1, iv, ie)*unknvt(2, iv)
munault(2, 1, iv, ie)    =  munaclt(1, iv, ie)*unknvt(3, iv)

munault(1:2, 2, iv, ie)    =  munaclt(2, iv, ie)*unknvt(2:3, iv)
!
! if(ipt(iv).eq.1.and.(ie.eq.5.or.ie.eq.6))then
!print*,'p11 muacn(1) post',iv,ie,munaclt(1:2, iv, ie),murie(iv),vnorm(3,1,iv),aujmp(3,iv),&
!                       abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))
! endif

!
!    snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!    snsigml(2, 1, iv, ie)= sigma(1, 2, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
!    snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!    snsigml(2, 2, iv, ie)= sigma(1, 2, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!
snsigmlt(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv)

snsigmlt(2, 1, iv, ie)= sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
snsigmlt(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv)

snsigmlt(2, 2, iv, ie)= sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)!
enddo
!
250 enddo  !...(1)ie = 1,nelem!

end subroutine getriem_tria
!
!...subroutine: Calculate the Riemann input for hybrid quad grids...
!
subroutine getriem_quad(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::coord, coold
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
!
real*8, dimension(1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
!
real*8, dimension(1:ndimn,  1:nvqua, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:2,  1:nvqua, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:2,  1:nvqua, 1:nquad), intent(out)::snsigmlq
!...Local integer
integer::ie, ideg, ielem, ifa, iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:3, 1:nvqua)::bq
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8::aujmp(1:3, 1:nvqua)
real*8::vnorm(1:3, 1:2, 1:nvqua)
real*8::sigma(1:2, 1:2, 1:nvqua)
real*8,dimension(1:nvqua)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
!
real*8::eps,c00,c05,c10,c20
real*8::rho, rhom, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s
real*8::acnx, acny
real*8::rhoi, rhon
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
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelagq(1:3, 4, ie); vnorm(1:3, 2, 1) = gelagq(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 1, 2) = gelagq(1:3, 1, ie); vnorm(1:3, 2, 2) = gelagq(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gelagq(1:3, 2, ie); vnorm(1:3, 2, 3) = gelagq(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 1, 4) = gelagq(1:3, 3, ie); vnorm(1:3, 2, 4) = gelagq(1:3, 4, ie) !...For point ip(3)
!
!...ndA=0.5d0*vnorm
!
vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhom = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
!...now configuration
!
r = rc; s= sc
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rhom = 1.d0/rhon
!
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...
!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
!
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
if(ndens.eq.1)then
rho  = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
!...now configuration
!
r = xvq(iv); s= yvq(iv)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rho = rhon

endif
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!...Limiter
!
if(nlimi.eq.1)then
rhomv = rhom + aflim(1, ielem)*(unknvq(1, iv) - rhom)
rho = 1.d0/rhomv
!
uvtx = uctr + aflim(2, ielem)*(unknvq(2, iv) - uctr)
vvtx = vctr + aflim(3, ielem)*(unknvq(3, iv) - vctr)
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!
!...updtae unknv(2:3,:)
unknvq(2, iv) = uvtx
unknvq(3 ,iv) = vvtx
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
!rhoct = 1.d0/unkno(1, 1, ielem)         !...ct denots center of one cell; cn denotes corner of one cell.
!uctr  = unkno(1, 2, ielem)
!vctr  = unkno(1, 3, ielem)
!ectr  = unkno(1, 4, ielem)
!pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...Get impedence coefficient...
!
do iv   = 1, nvqua
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ipq(iv).eq.31) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ipq(iv)),unknvq(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, nvqua
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
if(aujmp(3, iv)/sdctr.lt.1.d-6)then
!
!  print*,'itime',itime
munacn(ipq(iv)) = munacn(ipq(iv)) + murie(iv)*vnorm(3, ifa, iv)
!
munacu(1, ipq(iv)) =  munacu(1, ipq(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknvq(2, iv)
munacu(2, ipq(iv)) =  munacu(2, ipq(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknvq(3, iv)
else
munacn(ipq(iv)) = munacn(ipq(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
munacu(1, ipq(iv)) =  munacu(1, ipq(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknvq(2, iv)
munacu(2, ipq(iv)) =  munacu(2, ipq(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknvq(3, iv)
!
endif
!
!   if(ip(iv).eq.36) print*,'p36 muacn(vv) post',ie,murie(iv),munacn(ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p19 muacn(28) prep---',murie(iv), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                 vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
! if(ip(iv).eq.15) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Get the summed stress sum(n*A*sigma)
!
!
!  snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 1, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!
!  snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(1, 2, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!

snsigm(1, ipq(iv)) = snsigm(1, ipq(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) !
snsigm(2, ipq(iv)) = snsigm(2, ipq(iv)) + sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
if(aujmp(3, iv)/sdctr.lt.1.d-6)then

munaclq(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)
munaclq(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)

else

munaclq(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)* &
abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))

munaclq(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)* &
abs(vnorm(1, 2, iv)*aujmp(1, iv) + vnorm(2, 2, iv)*aujmp(2, iv))

endif
!
! if(ipq(iv).eq.31) print*,'p11 muacn(1) post',iv,ie,munaclq(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaulq(1, 1, iv, ie)    =  munaclq(1, iv, ie)*unknvq(2, iv)
munaulq(2, 1, iv, ie)    =  munaclq(1, iv, ie)*unknvq(3, iv)

munaulq(1:2, 2, iv, ie)  =  munaclq(2, iv, ie)*unknvq(2:3, iv)
!
!    snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!    snsigml(2, 1, iv, ie)= sigma(1, 2, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
!    snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!    snsigml(2, 2, iv, ie)= sigma(1, 2, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!
snsigmlq(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv)

snsigmlq(2, 1, iv, ie)= sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
snsigmlq(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv)

snsigmlq(2, 2, iv, ie)= sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)!
enddo
!
350 enddo  !...(1)ie = 1,nelem!

end subroutine getriem_quad

!
!...subroutine: Calculate the averaged velocity for hybrid grids...
!
subroutine getriem(iptri, ipqua, geoel, gelag, gelagq, vlave, unkno, munacn, munacu, snsigm,&
                   munaclt, munault, snsigmlt, munaclq, munaulq, snsigmlq)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelg, 1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
!
real*8, dimension(1:npoin),          intent(out)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(out)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(out)::snsigm
!
real*8, dimension(1:ndimn, 1:nvtri, 1:ntria),       intent(out)::munaclt
real*8, dimension(1:ndimn, 1:2, 1:nvtri,  1:ntria), intent(out)::munault
real*8, dimension(1:ndimn, 1:2,  1:nvtri,  1:ntria), intent(out)::snsigmlt
!
real*8, dimension(1:ndimn,  1:nvqua, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:2,  1:nvqua, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:2,  1:nvqua, 1:nquad), intent(out)::snsigmlq
!...Local integer
integer::ie, ideg, ielem, ifa, iv 
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:3, 1:nvtri)::bt
real*8,dimension(1:3, 1:nvqua)::bq
real*8,dimension(1:nq,1:nvtri)::unknvt
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8::aujmp(1:3, 1:nvqua)
real*8::vnorm(1:3, 1:2, 1:nvqua)
real*8::sigma(1:2, 1:2, 1:nvqua)
real*8,dimension(1:nvqua)::murie
real*8,dimension(1:nvtri):: xv,  yv
real*8,dimension(1:nvqua):: xvq, yvq
!
real*8::eps,c00,c05,c10,c20
real*8::rho, rhom, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc
real*8::acnx, acny
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Zero out munacn
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
!...Triangle...
!
do 250 ie = 1,ntria !...(1)ie = 1,nelem
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
!...shape functions
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
!...Basis function
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds

enddo
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
!
!...ndA=0.5d0*vnorm
!
vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
!...cell averaged value...
!
rhom = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...
!...zero out unknv
unknvt = 0.d0
!
do iv   = 1,nvtri
!
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
enddo
!
rho  = 1.d0/unknvt(1, iv)
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
evtx = unknvt(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
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
!aujmp(1:2, iv) = usold(1:2, ip(iv)) - unknv(2:3, iv)
aujmp(1:2, iv) = vlave(1:2, ipt(iv)) - unknvt(2:3, iv)
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
!...Get impedence coefficient...
!
do iv   = 1, nvtri
dux= vlave(1, ipt(iv))-unknvt(2, iv)
duy= vlave(2, ipt(iv))-unknvt(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ip(iv).eq.5) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ip(iv)),unknv(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
if(aujmp(3, iv)/sdctr.lt.1.d-9)then
!
!  print*,'itime',itime
munacn(ipt(iv)) = munacn(ipt(iv)) + murie(iv)*vnorm(3, ifa, iv)
!
munacu(1, ipt(iv)) =  munacu(1, ipt(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknvt(2, iv)
munacu(2, ipt(iv)) =  munacu(2, ipt(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknvt(3, iv)
else
munacn(ipt(iv)) = munacn(ipt(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
munacu(1, ipt(iv)) =  munacu(1, ipt(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknvt(2, iv)
munacu(2, ipt(iv)) =  munacu(2, ipt(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknvt(3, iv)
!
endif
!
!   if(ip(iv).eq.36) print*,'p36 muacn(vv) post',ie,murie(iv),munacn(ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p19 muacn(28) prep---',murie(iv), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                 vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
! if(ip(iv).eq.15) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Get the summed stress sum(n*A*sigma)
!
!
!  snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 1, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!
!  snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(1, 2, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!

snsigm(1, ipt(iv)) = snsigm(1, ipt(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) !
snsigm(2, ipt(iv)) = snsigm(2, ipt(iv)) + sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
if(aujmp(3, iv)/sdctr.lt.1.d-9)then

munaclt(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)
munaclt(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)

else

munaclt(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)* &
abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))

munaclt(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)* &
abs(vnorm(1, 2, iv)*aujmp(1, iv) + vnorm(2, 2, iv)*aujmp(2, iv))

endif
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munault(1, 1, iv, ie)    =  munaclt(1, iv, ie)*unknvt(2, iv)
munault(2, 1, iv, ie)    =  munaclt(1, iv, ie)*unknvt(3, iv)

munault(1:2, 2, iv, ie)    =  munaclt(2, iv, ie)*unknvt(2:3, iv)
!
!    snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!    snsigml(2, 1, iv, ie)= sigma(1, 2, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
!    snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!    snsigml(2, 2, iv, ie)= sigma(1, 2, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!
snsigmlt(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv)

snsigmlt(2, 1, iv, ie)= sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
snsigmlt(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv)

snsigmlt(2, 2, iv, ie)= sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)!
enddo
!
250 enddo  !...(1)ie = 1,nelem!
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
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelagq(1:3, 4, ie); vnorm(1:3, 2, 1) = gelagq(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 1, 2) = gelagq(1:3, 1, ie); vnorm(1:3, 2, 2) = gelagq(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gelagq(1:3, 2, ie); vnorm(1:3, 2, 3) = gelagq(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 1, 4) = gelagq(1:3, 3, ie); vnorm(1:3, 2, 4) = gelagq(1:3, 4, ie) !...For point ip(3)
!
!...ndA=0.5d0*vnorm
!
vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
!...cell averaged value...
!
rhom = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...
!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
!
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
rho  = 1.d0/unknvq(1, iv)
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
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
!aujmp(1:2, iv) = usold(1:2, ip(iv)) - unknv(2:3, iv)
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
!...Get impedence coefficient...
!
do iv   = 1, nvqua
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ipq(iv).eq.31) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ipq(iv)),unknvq(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, nvqua
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
if(aujmp(3, iv)/sdctr.lt.1.d-13)then
!
!  print*,'itime',itime
munacn(ipq(iv)) = munacn(ipq(iv)) + murie(iv)*vnorm(3, ifa, iv)
!
munacu(1, ipq(iv)) =  munacu(1, ipq(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknvq(2, iv)
munacu(2, ipq(iv)) =  munacu(2, ipq(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*unknvq(3, iv)
else
munacn(ipq(iv)) = munacn(ipq(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
munacu(1, ipq(iv)) =  munacu(1, ipq(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknvq(2, iv)
munacu(2, ipq(iv)) =  munacu(2, ipq(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknvq(3, iv)
!
endif
!
!   if(ip(iv).eq.36) print*,'p36 muacn(vv) post',ie,murie(iv),munacn(ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p19 muacn(28) prep---',murie(iv), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                 vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
! if(ip(iv).eq.15) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Get the summed stress sum(n*A*sigma)
!
!
!  snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 1, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!
!  snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(1, 2, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
!                                          sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!

snsigm(1, ipq(iv)) = snsigm(1, ipq(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) !
snsigm(2, ipq(iv)) = snsigm(2, ipq(iv)) + sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
if(aujmp(3, iv)/sdctr.lt.1.d-13)then

munaclq(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)
munaclq(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)

else

munaclq(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)* &
abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))

munaclq(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)* &
abs(vnorm(1, 2, iv)*aujmp(1, iv) + vnorm(2, 2, iv)*aujmp(2, iv))

endif
!
! if(ipq(iv).eq.31) print*,'p11 muacn(1) post',iv,ie,munaclq(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaulq(1, 1, iv, ie)    =  munaclq(1, iv, ie)*unknvq(2, iv)
munaulq(2, 1, iv, ie)    =  munaclq(1, iv, ie)*unknvq(3, iv)

munaulq(1:2, 2, iv, ie)  =  munaclq(2, iv, ie)*unknvq(2:3, iv)
!
!    snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!    snsigml(2, 1, iv, ie)= sigma(1, 2, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
!    snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 1, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!    snsigml(2, 2, iv, ie)= sigma(1, 2, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
!                           sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!
snsigmlq(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv)

snsigmlq(2, 1, iv, ie)= sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
snsigmlq(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv)

snsigmlq(2, 2, iv, ie)= sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)!
enddo
!
350 enddo  !...(1)ie = 1,nelem!

end subroutine getriem
!
!...Face integral (mass center) for hybrid...
!
subroutine rhsifacedg_lag_mc_hybrid(inpoel, iptri, ipqua, unkno, ustar, fstar, fstarq, gelag, gelagq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri,1:ntria),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:2, 1:nvtri) :: ipf
integer,dimension(1:2, 1:nvqua) :: ipfq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri)::lpnpt
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvqua)::lpnpq
real*8::xvt(3), yvt(3),bt(1:3,1:nvtri)
real*8::xvq(nvqua), yvq(nvqua),bq(1:3,1:nvqua)

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
ipf(1, 1) = 3; ipf(2, 1) = 2
ipf(1, 2) = 1; ipf(2, 2) = 3
ipf(1, 3) = 2; ipf(2, 3) = 1
!
ipfq(1, 1) = 4; ipfq(2, 1) = 2
ipfq(1, 2) = 1; ipfq(2, 2) = 3
ipfq(1, 3) = 2; ipfq(2, 3) = 4
ipfq(1, 4) = 3; ipfq(2, 4) = 1

do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
ielem = ie
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
!
do iv =1 ,nvtri
!
!print*,'iv', ie
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds
enddo
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpt(1:ndimn, ig, 1, 1) = c16*(2.d0*bt(ig, 1) + bt(ig, 3))*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
lpnpt(1:ndimn, ig, 2, 1) = c16*(2.d0*bt(ig, 1) + bt(ig, 2))*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
!
!...point 2
lpnpt(1:ndimn, ig, 1, 2) = c16*(2.d0*bt(ig, 2) + bt(ig, 1))*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
lpnpt(1:ndimn, ig, 2, 2) = c16*(2.d0*bt(ig, 2) + bt(ig, 3))*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
!
!...point 3
lpnpt(1:ndimn, ig, 1, 3) = c16*(2.d0*bt(ig, 3) + bt(ig, 2))*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
lpnpt(1:ndimn, ig, 2, 3) = c16*(2.d0*bt(ig, 3) + bt(ig, 1))*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
!
enddo

!...The vertex constituting one cell...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
!

!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvtri
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipt(iv))*lpnpt(1, 1:ndegr, 1, iv) +&
ustar(2, ipt(iv))*lpnpt(2, 1:ndegr, 1, iv) +&
ustar(1, ipt(iv))*lpnpt(1, 1:ndegr, 2, iv) +&
ustar(2, ipt(iv))*lpnpt(2, 1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstar(1, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
fstar(1, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv)))

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstar(2, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
fstar(2, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv)))
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipt(iv))*fstar(1, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
ustar(2, ipt(iv))*fstar(2, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
ustar(1, ipt(iv))*fstar(1, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv))) +&
ustar(2, ipt(iv))*fstar(2, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv)))
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
! if(ie==18) print*,'rhs iface',rhsel(1, 1, ie), lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

550 enddo
!
!   open(8,file='lpnp.dat')
!    do ie = 1, nelem
!
!      do iv = 1, nvtri
!
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 1, 1, iv, ie), lpnp(1:2, 1, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 2, 1, iv, ie), lpnp(1:2, 2, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 3, 1, iv, ie), lpnp(1:2, 3, 2, iv, ie)
!      enddo
!    enddo
!   close(8)
!
!
!
! ip1 =8; ip2=37
!  print*,'dx', coord(1,ip1) ,coord(1, ip2), coord(2, ip1) -coord(2, ip2),&
!               0.5d0*(coord(2,ip1) -coord(2, ip2)),-0.5d0*(coord(1,ip1) -coord(1, ip2))
!/sqrt((coord(1,37) -coord(1, 2))**2+(coord(2,37) -coord(2, 2))**2)
!
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
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
!
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpq(1:ndimn, ig, 1, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 4))*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)
lpnpq(1:ndimn, ig, 2, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 2))*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)
!
!...point 2
lpnpq(1:ndimn, ig, 1, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 1))*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)
lpnpq(1:ndimn, ig, 2, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 3))*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)
!
!...point 3
lpnpq(1:ndimn, ig, 1, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 2))*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)
lpnpq(1:ndimn, ig, 2, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 4))*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)
!...point 4
lpnpq(1:ndimn, ig, 1, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 3))*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)
lpnpq(1:ndimn, ig, 2, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 1))*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)
!
enddo

!...The vertex constituting one cell...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!

!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvqua
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 1, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 1, iv) +&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 2, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
fstarq(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv)))

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
fstarq(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv)))
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(iv))*fstarq(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
ustar(2, ipq(iv))*fstarq(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
ustar(1, ipq(iv))*fstarq(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv))) +&
ustar(2, ipq(iv))*fstarq(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv)))
!
! if(ie==21) print*,iv,ie,ielem,fstarq(1, 1:2, iv,ie)
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!  print*,'rhs iface',ielem, ie,plnpn(1,1),fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lag_mc_hybrid
!
!....domain integral for hybrid linear cells
!
subroutine rhsdomndg_lag_mc_hybrid(intfac, inpoel, iptri,ipqua, coord, geoel, unkno, rhsel,aflim )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8, dimension(1:nq+1, 1:nsize),      intent(in)::aflim
real*8,dimension(1:ndegr,1:nq,1:ncell),  intent(inout)::rhsel
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
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8::weigh(ngausd), posi(2,ngausd)
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::rho,uadv,vadv,eadv,rhom
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr, rhomv
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
call rutope(2, ngausd, posi, weigh)
call ruqope(2, ngausdq, posiq, weighq)

!
!...Loop over elements
!
do 550 ie = 1,ntria!...(1)ie = 1,nelem
!
 ielem = ie
!
!...Points consitituting one element...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
!
xp(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xp(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
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
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
wi     = weigh(ig)
!
dspr(1) = -c10
dspr(2) =  c10
dspr(3) =  0.d0
!
dsps(1) = -c10
dsps(2) =  0.d0
dsps(3) =  c10
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
djaco = 0.5d0*wi
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
rhom = unknod(1)
rho  = 1.d0/rhom
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
if(nlimi.eq.1)then
!
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
rhomv = rhomc + aflim(1 ,ie)*(rhomv - rhomc)
rho = 1.d0/rhomv
!
uadv = uctr + aflim(2, ie)*(uadv - uctr)
vadv = vctr + aflim(3 ,ie)*(vadv - vctr)
!
pres = pctr + aflim(4, ie)*(pres- pctr)

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
enddo !...(2)ig = 1,ngausd
!
550 enddo
!
!...Quads
!
!
!...Loop over elements
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:4,ie)!
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
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
!
dsprq(1) = -0.25d0   *sm
dsprq(2) =  0.25d0   *sm
dsprq(3) =  0.25d0   *sp
dsprq(4) = -0.25d0   *sp
!
dspsq(1) =-0.25d0*rm
dspsq(2) =-0.25d0*rp
dspsq(3) = 0.25d0*rp
dspsq(4) = 0.25d0*rm
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 4
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
rhom = unknod(1)
rho  = 1.d0/rhom
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
if(nlimi.eq.1)then
!
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
rhomv = rhomc + aflim(1 ,ie)*(rhomv - rhomc)
rho = 1.d0/rhomv
!
uadv = uctr + aflim(2, ie)*(uadv - uctr)
vadv = vctr + aflim(3 ,ie)*(vadv - vctr)
!
pres = pctr + aflim(4, ie)*(pres- pctr)

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
enddo !...(2)ig = 1,ngausd
!
650 enddo

end subroutine rhsdomndg_lag_mc_hybrid
!
!....domain integral for hybrid linear triangle cells
!
subroutine rhsdomndg_lag_mc_hybridtria(intfac, iptri, coord, geoel, unkno, rhsel,aflim )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
real*8, dimension(1:nq+1, 1:nsize),      intent(in)::aflim
real*8,dimension(1:ndegr,1:nq,1:ncell),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem
!
!...local integer array
!
integer,dimension(1:nvtri) :: ipt
!...local real array
real*8,dimension(1:ndimn, 1:nvtri) :: xp
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8::weigh(ngausd), posi(2,ngausd)
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::rho,uadv,vadv,eadv,rhom
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr, rhomv
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
call rutope(2, ngausd, posi, weigh)

!
!...Loop over elements
!
do 550 ie = 1,ntria!...(1)ie = 1,nelem
!
 ielem = ie
!
!...Points consitituting one element...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
!
xp(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xp(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
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
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
wi     = weigh(ig)
!
dspr(1) = -c10
dspr(2) =  c10
dspr(3) =  0.d0
!
dsps(1) = -c10
dsps(2) =  0.d0
dsps(3) =  c10
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
djaco = 0.5d0*wi
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
rhom = unknod(1)
rho  = 1.d0/rhom
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!...impose imiter...
!
if(nlimi.eq.1)then
!
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
rhomv = rhom
rhomv = rhomc + aflim(1 ,ielem)*(rhomv - rhomc)
rho = 1.d0/rhomv
!
uadv = uctr + aflim(2, ielem)*(uadv - uctr)
vadv = vctr + aflim(3 ,ielem)*(vadv - vctr)
!
pres = pctr + aflim(4, ielem)*(pres- pctr)

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
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomndg_lag_mc_hybridtria
!
!....domain integral for hybrid linear quad cells
!
subroutine rhsdomndg_lag_mc_hybridquad(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8, dimension(1:nq+1, 1:nsize),      intent(in)::aflim
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
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xpqi
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::rho,uadv,vadv,eadv,rhom
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr, rhomv
real*8::rhoi, rhon
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
!
!...Loop over elements
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:4,ie)!
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
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
!
dsprq(1) = -0.25d0   *sm
dsprq(2) =  0.25d0   *sm
dsprq(3) =  0.25d0   *sp
dsprq(4) = -0.25d0   *sp
!
dspsq(1) =-0.25d0*rm
dspsq(2) =-0.25d0*rp
dspsq(3) = 0.25d0*rp
dspsq(4) = 0.25d0*rm
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 4
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
if(ndens.eq.1)then
!
rhom = unknod(1)
rho  = 1.d0/rhom
!
elseif(ndens.eq.2)then
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rhom = 1.d0/rhon
rho = rhon
!
endif
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
if(nlimi.eq.1)then
!
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
rhomv = rhom
rhomv = rhomc + aflim(1 ,ielem)*(rhomv - rhomc)
rho = 1.d0/rhomv
!
uadv = uctr + aflim(2, ielem)*(uadv - uctr)
vadv = vctr + aflim(3 ,ielem)*(vadv - vctr)
!
pres = pctr + aflim(4, ielem)*(pres- pctr)

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
enddo !...(2)ig = 1,ngausd
!
650 enddo

end subroutine rhsdomndg_lag_mc_hybridquad
!
!...Face integral (mass center) for hybrid...
!
subroutine rhsifacedg_lag_mc_hybridtria(iptri, unkno, ustar, fstar, gelag, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri,1:ntria),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:2, 1:nvtri) :: ipf
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri)::lpnpt
real*8::xvt(3), yvt(3),bt(1:3,1:nvtri)
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
ipf(1, 1) = 3; ipf(2, 1) = 2
ipf(1, 2) = 1; ipf(2, 2) = 3
ipf(1, 3) = 2; ipf(2, 3) = 1
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
ielem = ie
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
!
do iv =1 ,nvtri
!
!print*,'iv', ie
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds
enddo
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpt(1:ndimn, ig, 1, 1) = c16*(2.d0*bt(ig, 1) + bt(ig, 3))*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
lpnpt(1:ndimn, ig, 2, 1) = c16*(2.d0*bt(ig, 1) + bt(ig, 2))*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
!
!...point 2
lpnpt(1:ndimn, ig, 1, 2) = c16*(2.d0*bt(ig, 2) + bt(ig, 1))*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
lpnpt(1:ndimn, ig, 2, 2) = c16*(2.d0*bt(ig, 2) + bt(ig, 3))*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
!
!...point 3
lpnpt(1:ndimn, ig, 1, 3) = c16*(2.d0*bt(ig, 3) + bt(ig, 2))*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
lpnpt(1:ndimn, ig, 2, 3) = c16*(2.d0*bt(ig, 3) + bt(ig, 1))*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
!
enddo

!...The vertex constituting one cell...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvtri
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipt(iv))*lpnpt(1, 1:ndegr, 1, iv) +&
ustar(2, ipt(iv))*lpnpt(2, 1:ndegr, 1, iv) +&
ustar(1, ipt(iv))*lpnpt(1, 1:ndegr, 2, iv) +&
ustar(2, ipt(iv))*lpnpt(2, 1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstar(1, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
fstar(1, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv)))

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstar(2, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
fstar(2, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv)))
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipt(iv))*fstar(1, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
ustar(2, ipt(iv))*fstar(2, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
ustar(1, ipt(iv))*fstar(1, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv))) +&
ustar(2, ipt(iv))*fstar(2, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv)))
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
! if(ie==18) print*,'rhs iface',rhsel(1, 1, ie), lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

550 enddo
!
!   open(8,file='lpnp.dat')
!    do ie = 1, nelem
!
!      do iv = 1, nvtri
!
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 1, 1, iv, ie), lpnp(1:2, 1, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 2, 1, iv, ie), lpnp(1:2, 2, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 3, 1, iv, ie), lpnp(1:2, 3, 2, iv, ie)
!      enddo
!    enddo
!   close(8)
!
!
!
! ip1 =8; ip2=37
!  print*,'dx', coord(1,ip1) ,coord(1, ip2), coord(2, ip1) -coord(2, ip2),&
!               0.5d0*(coord(2,ip1) -coord(2, ip2)),-0.5d0*(coord(1,ip1) -coord(1, ip2))
!/sqrt((coord(1,37) -coord(1, 2))**2+(coord(2,37) -coord(2, 2))**2)
!
!...Quads...
!

!
end subroutine rhsifacedg_lag_mc_hybridtria
!
!...Face integral (mass center) for hybrid linear triangle using gauss integration...
!
subroutine rhsifacedg_lagmc_triagauss(iptri, unkno, ustar, fstar, gelag, geoel,coord,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri,1:ntria),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),     intent(in)::coord
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:2, 1:nvtri) :: ipf
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri)::lpnpt
real*8::xvt(3), yvt(3),bt(1:3,1:nvtri)
real*8,dimension(1:ndimn, 1:nvtri) :: xpht
real*8,dimension(1:nvtri):: rcoet
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
ipf(1, 1) = 3; ipf(2, 1) = 2
ipf(1, 2) = 1; ipf(2, 2) = 3
ipf(1, 3) = 2; ipf(2, 3) = 1
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
ielem = ie
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
!
do iv =1 ,nvtri
!
!print*,'iv', ie
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds
enddo
!
!...physical coordinate....
!
xpht(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpht(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
!...Coefficient R of RZ or XY system...
!
rcoet(1:nvtri) = 1.d0 - alfrz + alfrz*xpht(2, 1:nvtri)
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpt(1:ndimn, ig, 1, 1) = (2.d0*rcoet(1) + rcoet(3))/3.d0*0.5d0*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
lpnpt(1:ndimn, ig, 2, 1) = (2.d0*rcoet(1) + rcoet(2))/3.d0*0.5d0*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
!
!...point 2
lpnpt(1:ndimn, ig, 1, 2) = (2.d0*rcoet(2) + rcoet(1))/3.d0*0.5d0*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
lpnpt(1:ndimn, ig, 2, 2) = (2.d0*rcoet(2) + rcoet(3))/3.d0*0.5d0*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
!
!...point 3
lpnpt(1:ndimn, ig, 1, 3) = (2.d0*rcoet(3) + rcoet(2))/3.d0*0.5d0*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
lpnpt(1:ndimn, ig, 2, 3) = (2.d0*rcoet(3) + rcoet(1))/3.d0*0.5d0*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
!
enddo
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvtri
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipt(iv))*lpnpt(1, 1:ndegr, 1, iv) +&
ustar(2, ipt(iv))*lpnpt(2, 1:ndegr, 1, iv) +&
ustar(1, ipt(iv))*lpnpt(1, 1:ndegr, 2, iv) +&
ustar(2, ipt(iv))*lpnpt(2, 1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstar(1, 1, iv, ie)*bt(1:ndegr, iv)    +&
fstar(1, 2, iv, ie)*bt(1:ndegr, iv)

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstar(2, 1, iv, ie)*bt(1:ndegr, iv)   +&
fstar(2, 2, iv, ie)*bt(1:ndegr, iv)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipt(iv))*fstar(1, 1, iv, ie)*bt(1:ndegr, iv)  +&
ustar(2, ipt(iv))*fstar(2, 1, iv, ie)*bt(1:ndegr, iv)  +&
ustar(1, ipt(iv))*fstar(1, 2, iv, ie)*bt(1:ndegr, iv)  +&
ustar(2, ipt(iv))*fstar(2, 2, iv, ie)*bt(1:ndegr, iv)
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
! if(ie==18) print*,'rhs iface',rhsel(1, 1, ie), lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

550 enddo
!
!   open(8,file='lpnp.dat')
!    do ie = 1, nelem
!
!      do iv = 1, nvtri
!
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 1, 1, iv, ie), lpnp(1:2, 1, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 2, 1, iv, ie), lpnp(1:2, 2, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 3, 1, iv, ie), lpnp(1:2, 3, 2, iv, ie)
!      enddo
!    enddo
!   close(8)

!
end subroutine rhsifacedg_lagmc_triagauss
!
!...Face integral (mass center) for hybrid quad...
!
subroutine rhsifacedg_lag_mc_hybridquad(ipqua, unkno, ustar,fstarq, gelagq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:2, 1:nvqua) :: ipfq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvqua)::lpnpq
real*8::xvq(nvqua), yvq(nvqua),bq(1:3,1:nvqua)

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
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
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
!
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpq(1:ndimn, ig, 1, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 4))*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)
lpnpq(1:ndimn, ig, 2, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 2))*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)
!
!...point 2
lpnpq(1:ndimn, ig, 1, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 1))*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)
lpnpq(1:ndimn, ig, 2, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 3))*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)
!
!...point 3
lpnpq(1:ndimn, ig, 1, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 2))*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)
lpnpq(1:ndimn, ig, 2, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 4))*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)
!...point 4
lpnpq(1:ndimn, ig, 1, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 3))*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)
lpnpq(1:ndimn, ig, 2, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 1))*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)
!
enddo
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
do iv = 1, nvqua
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 1, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 1, iv) +&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 2, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
fstarq(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv)))

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
fstarq(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv)))
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(iv))*fstarq(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
ustar(2, ipq(iv))*fstarq(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
ustar(1, ipq(iv))*fstarq(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv))) +&
ustar(2, ipq(iv))*fstarq(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv)))
!
!if(ie==1) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1)  print*,'rhs iface',ielem, ie,ulnpn(1)!,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lag_mc_hybridquad
!
!...Face integral (mass center) for hybrid quad...
!
subroutine rhsifacedg_lag_mc_hybridquad2(ipqua, unkno, ustar,fstarq, gelagq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:2, 1:nvqua) :: ipfq
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
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
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
!
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
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
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpq(1:ndimn, ig, 1, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 1))*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)
lpnpq(1:ndimn, ig, 2, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 1))*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)
!
!...point 2
lpnpq(1:ndimn, ig, 1, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 2))*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)
lpnpq(1:ndimn, ig, 2, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 2))*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)
!
!...point 3
lpnpq(1:ndimn, ig, 1, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 3))*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)
lpnpq(1:ndimn, ig, 2, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 3))*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)
!...point 4
lpnpq(1:ndimn, ig, 1, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 4))*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)
lpnpq(1:ndimn, ig, 2, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 4))*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)
!
enddo
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
do iv = 1, nvqua
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 1, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 1, iv) +&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 2, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
fstarq(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv))

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
fstarq(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv))
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(iv))*fstarq(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
ustar(2, ipq(iv))*fstarq(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
ustar(1, ipq(iv))*fstarq(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
ustar(2, ipq(iv))*fstarq(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv))
!
!if(ie==1) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1)  print*,'rhs iface',ielem, ie,ulnpn(1)!,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lag_mc_hybridquad2
!
!...Face integral (mass center) for hybrid linear quad using gauss integration...
!
subroutine rhsifacedg_lagmc_quadgauss(ipqua, unkno, ustar,fstarq, gelagq, geoel,coord,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:2, 1:nvqua) :: ipfq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvqua)::lpnpq
real*8::xvq(nvqua), yvq(nvqua),bq(1:3,1:nvqua)
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:ndimn, 1:nvqua) :: xphq
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
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
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
!
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...physical coordinate....
!
xphq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xphq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Coefficient R of RZ or XY system...
!
rcoeq(1:nvqua) = 1.d0 - alfrz + alfrz*xphq(2, 1:nvqua)
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpq(1:ndimn, ig, 1, 1) = c05*(2.d0*rcoeq(1) + rcoeq(4))/3.d0*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 1)
lpnpq(1:ndimn, ig, 2, 1) = c05*(2.d0*rcoeq(1) + rcoeq(2))/3.d0*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 1)
!
!...point 2
lpnpq(1:ndimn, ig, 1, 2) = c05*(2.d0*rcoeq(2) + rcoeq(1))/3.d0*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 2)
lpnpq(1:ndimn, ig, 2, 2) = c05*(2.d0*rcoeq(2) + rcoeq(3))/3.d0*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 2)
!
!...point 3
lpnpq(1:ndimn, ig, 1, 3) = c05*(2.d0*rcoeq(3) + rcoeq(2))/3.d0*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 3)
lpnpq(1:ndimn, ig, 2, 3) = c05*(2.d0*rcoeq(3) + rcoeq(4))/3.d0*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 3)
!...point 4
lpnpq(1:ndimn, ig, 1, 4) = c05*(2.d0*rcoeq(4) + rcoeq(3))/3.d0*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 4)
lpnpq(1:ndimn, ig, 2, 4) = c05*(2.d0*rcoeq(4) + rcoeq(1))/3.d0*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 4)
!
enddo
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvqua
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 1, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 1, iv) +&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 2, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, 1, iv, ie)*bq(1:ndegr, iv)*rcoeq(iv)   +&
fstarq(1, 2, iv, ie)*bq(1:ndegr, iv)*rcoeq(iv)

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, 1, iv, ie)*bq(1:ndegr, iv)*rcoeq(iv)  +&
fstarq(2, 2, iv, ie)*bq(1:ndegr, iv)*rcoeq(iv)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(iv))*fstarq(1, 1, iv, ie)*bq(1:ndegr, iv)*rcoeq(iv)  +&
ustar(2, ipq(iv))*fstarq(2, 1, iv, ie)*bq(1:ndegr, iv)*rcoeq(iv)  +&
ustar(1, ipq(iv))*fstarq(1, 2, iv, ie)*bq(1:ndegr, iv)*rcoeq(iv)  +&
ustar(2, ipq(iv))*fstarq(2, 2, iv, ie)*bq(1:ndegr, iv)*rcoeq(iv)
!
!if(ie==1) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1)  print*,'rhs iface',ielem, ie,ulnpn(1)!,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lagmc_quadgauss
!
!...Face integral using gauss quadrature distribution on RZ quads...
!
subroutine rhsifacedg_lagquadrz_simpson(ipqua,  unkno, ustar, fstarq, gelagq, geoel, coord,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::rhsel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize), intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ielem
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer, dimension(3, 4)::fglvq, fglgq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::vnorm(3, 12)
real*8::xvq(9), yvq(9)
real*8::bg(3, 3)
real*8::rcoeq(3)
real*8::usq(2, 3), fsq(2 ,3)
real*8::xphq(2, 4)
!
real*8::weigh(3), posi(1,3)
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
!...physical coordinate....
!
xphq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xphq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...LN for every face...
!
vnorm(1:3,  1) = gelagq(1:3, 1, ie); vnorm(1:3,  2) = gelagq(1:3, 1, ie);vnorm(1:3,  9) = gelagq(1:3, 1, ie)
vnorm(1:3,  3) = gelagq(1:3, 2, ie); vnorm(1:3,  4) = gelagq(1:3, 2, ie);vnorm(1:3, 10) = gelagq(1:3, 2, ie)
vnorm(1:3,  5) = gelagq(1:3, 3, ie); vnorm(1:3,  6) = gelagq(1:3, 3, ie);vnorm(1:3, 11) = gelagq(1:3, 3, ie)
vnorm(1:3,  7) = gelagq(1:3, 4, ie); vnorm(1:3,  8) = gelagq(1:3, 4, ie);vnorm(1:3, 12) = gelagq(1:3, 4, ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do ifa = 1, 4
!
rcoeq(1) = (2.d0*xphq(2, fglvq(1, ifa)) + xphq(2, fglvq(2, ifa)))/3.d0
rcoeq(2) = (2.d0*xphq(2, fglvq(2, ifa)) + xphq(2, fglvq(1, ifa)))/3.d0
rcoeq(3) = (     xphq(2, fglvq(2, ifa)) + xphq(2, fglvq(1, ifa)))/2.d0
!
usq(1:2, 1) = ustar(1:2, ipq(fglvq(1, ifa)))
usq(1:2, 2) = ustar(1:2, ipq(fglvq(2, ifa)))
usq(1:2, 3) = 0.5d0*(usq(1:2, 1) + usq(1:2, 2))
!
fsq(1:2, 1) = fstarq(1:2, 2, fglvq(1, ifa), ie)*2.d0/rcoeq(1)
fsq(1:2, 2) = fstarq(1:2, 1, fglvq(2, ifa), ie)*2.d0/rcoeq(2)
fsq(1:2, 3) = 0.5d0*(fsq(1:2, 1) + fsq(1:2, 2))
!
!
do ig   = 1, 3
!
wi  = weigh(ig)
!
rg = xvq(fglvq(ig, ifa))
sg = yvq(fglvq(ig, ifa))
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
gpnx = vnorm(1, fglgq(ig, ifa))
gpny = vnorm(2, fglgq(ig, ifa))
gpsa = vnorm(3, fglgq(ig, ifa))
!
!...Distribute to every corner...
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
usq(1, ig)*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig)*rcoeq(ig) +&
usq(2, ig)*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)*rcoeq(ig)
!
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fsq(1, ig)*bg(1:ndegr, ig)*weigh(ig)*rcoeq(ig)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fsq(2, ig)*bg(1:ndegr, ig)*weigh(ig)*rcoeq(ig)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
usq(1, ig)*fsq(1, ig)*bg(1:ndegr, ig)*weigh(ig)*rcoeq(ig) +&
usq(2, ig)*fsq(2, ig)*bg(1:ndegr, ig)*weigh(ig)*rcoeq(ig)
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
! if(ie==23) print*,'rhs iface',rhsel(1:3, 1, ie), ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)
550 enddo
!
end subroutine rhsifacedg_lagquadrz_simpson
!
!...Source term integration for hybrid triangle grids...
!
subroutine rhsdomnsrcdg_lag_mc_hybrid(intfac, inpoel,  iptri,ipqua, coord, geoel,rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nsize),          intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:ncell),  intent(inout)::rhsel
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
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8::weigh(ngausd), posi(2,ngausd)
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
call rutope(2, ngausd, posi, weigh)
call ruqope(2, ngausdq, posiq, weighq)
!
!...Loop over elements
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
 ielem = ie
!
!...Points consitituting one element...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
!
xp(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xp(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!rxx=
!rxy=
!ryy=
!
!...Gauss loop
!
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
wi     = weigh(ig)
!
dspr(1) = -c10
dspr(2) =  c10
dspr(3) =  0.d0
!
dsps(1) = -c10
dsps(2) =  0.d0
dsps(3) =  c10
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
!...Gauss points...
!
xg =0.d0; yg= 0.d0
!
xg = shp(1)*xp(1, 1) + shp(2)*xp(1, 2) + shp(3)*xp(1, 3)
yg = shp(1)*xp(2, 1) + shp(2)*xp(2, 2) + shp(3)*xp(2, 3)
!
!
!...Basis function for solutions...
!
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
!
!...source term
!
src = 0.25d0*pi*(cos(3.d0*pi*xg)*cos(pi*yg) - cos(3.d0*pi*yg)*cos(pi*xg))/(gamlg-1.d0)
!
!    src = 0.5d0*pi/(gamlg-1.d0)*(sin(2.d0*pi*yg)*cos(pi*xg)*sin(pi*yg) - sin(2.d0*pi*xg)*sin(pi*xg)*cos(pi*yg))
!
!finally, scatter the contribution to the RHS
!
! if(ie==2) print*,'src rhs', 3.5d0/9.d0,rhsel(1, 4, 2),coorp(1, 1:3), coorp(2, 1:3)
!
do ideg = 1,ndegr
rhsel(ideg,4,ielem)=rhsel(ideg,4,ielem) + src*b(ideg)*djaco
enddo
!
!if(ie==2) print*,'src rhs', rhsel(1, 4, 2),src,b(1),djaco,src*b(1)*djaco,xg,yg
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
!
!...Loop over quads
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:4,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
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
!
dsprq(1) = -0.25d0   *sm
dsprq(2) =  0.25d0   *sm
dsprq(3) =  0.25d0   *sp
dsprq(4) = -0.25d0   *sp
!
dspsq(1) =-0.25d0*rm
dspsq(2) =-0.25d0*rp
dspsq(3) = 0.25d0*rp
dspsq(4) = 0.25d0*rm
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 4
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
xg = shpq(1)*xpq(1, 1) + shpq(2)*xpq(1, 2) + shpq(3)*xpq(1, 3) + shpq(4)*xpq(1, 4)  
yg = shpq(1)*xpq(2, 1) + shpq(2)*xpq(2, 2) + shpq(3)*xpq(2, 3) + shpq(4)*xpq(2, 4) 
!
!
!...Basis function for solutions...
!
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
!
!...source term
!
src = 0.25d0*pi*(cos(3.d0*pi*xg)*cos(pi*yg) - cos(3.d0*pi*yg)*cos(pi*xg))/(gamlg-1.d0)
!
!    src = 0.5d0*pi/(gamlg-1.d0)*(sin(2.d0*pi*yg)*cos(pi*xg)*sin(pi*yg) - sin(2.d0*pi*xg)*sin(pi*xg)*cos(pi*yg))
!
!finally, scatter the contribution to the RHS
!
! if(ie==2) print*,'src rhs', 3.5d0/9.d0,rhsel(1, 4, 2),coorp(1, 1:3), coorp(2, 1:3)
!
do ideg = 1,ndegr
rhsel(ideg,4,ielem)=rhsel(ideg,4,ielem) + src*b(ideg)*djaco
enddo
!
!if(ie==2) print*,'src rhs', rhsel(1, 4, 2),src,b(1),djaco,src*b(1)*djaco,xg,yg
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
end subroutine rhsdomnsrcdg_lag_mc_hybrid
!
!...Source term integration for hybrid triangle grids...
!
subroutine rhsdomnsrcdg_lag_mc_hybridtria(intfac, iptri,coord, geoel,rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nsize),          intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
real*8,dimension(1:ndegr,1:nq,1:ncell),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem
!
!...local integer array
!
integer,dimension(1:nvtri) :: ipt
!...local real array
real*8,dimension(1:ndimn, 1:nvtri) :: xp
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8::weigh(ngausd), posi(2,ngausd)
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
call rutope(2, ngausd, posi, weigh)
!
!...Loop over elements
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
ielem = ie
!
!...Points consitituting one element...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
!
xp(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xp(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!rxx=
!rxy=
!ryy=
!
!...Gauss loop
!
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
wi     = weigh(ig)
!
dspr(1) = -c10
dspr(2) =  c10
dspr(3) =  0.d0
!
dsps(1) = -c10
dsps(2) =  0.d0
dsps(3) =  c10
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
!...Gauss points...
!
xg =0.d0; yg= 0.d0
!
xg = shp(1)*xp(1, 1) + shp(2)*xp(1, 2) + shp(3)*xp(1, 3)
yg = shp(1)*xp(2, 1) + shp(2)*xp(2, 2) + shp(3)*xp(2, 3)
!
!
!...Basis function for solutions...
!
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
!
!...source term
!
src = 0.25d0*pi*(cos(3.d0*pi*xg)*cos(pi*yg) - cos(3.d0*pi*yg)*cos(pi*xg))/(gamlg-1.d0)
!
!    src = 0.5d0*pi/(gamlg-1.d0)*(sin(2.d0*pi*yg)*cos(pi*xg)*sin(pi*yg) - sin(2.d0*pi*xg)*sin(pi*xg)*cos(pi*yg))
!
!finally, scatter the contribution to the RHS
!
! if(ie==2) print*,'src rhs', 3.5d0/9.d0,rhsel(1, 4, 2),coorp(1, 1:3), coorp(2, 1:3)
!
do ideg = 1,ndegr
rhsel(ideg,4,ielem)=rhsel(ideg,4,ielem) + src*b(ideg)*djaco
enddo
!
!if(ie==2) print*,'src rhs', rhsel(1, 4, 2),src,b(1),djaco,src*b(1)*djaco,xg,yg
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomnsrcdg_lag_mc_hybridtria
!
!...Source term integration for hybrid quad grids...
!
subroutine rhsdomnsrcdg_lag_mc_hybridquad(intfac, ipqua, coord, geoel,rhsel)
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
ipq(1:nvqua) = ipqua(1:4,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
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
!
dsprq(1) = -0.25d0   *sm
dsprq(2) =  0.25d0   *sm
dsprq(3) =  0.25d0   *sp
dsprq(4) = -0.25d0   *sp
!
dspsq(1) =-0.25d0*rm
dspsq(2) =-0.25d0*rp
dspsq(3) = 0.25d0*rp
dspsq(4) = 0.25d0*rm
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 4
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
xg = shpq(1)*xpq(1, 1) + shpq(2)*xpq(1, 2) + shpq(3)*xpq(1, 3) + shpq(4)*xpq(1, 4)
yg = shpq(1)*xpq(2, 1) + shpq(2)*xpq(2, 2) + shpq(3)*xpq(2, 3) + shpq(4)*xpq(2, 4)
!
!
!...Basis function for solutions...
!
b(1) = 1.d0

if(npoly.ge.1)then
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
endif
!
!...source term
!
src = 0.25d0*pi*(cos(3.d0*pi*xg)*cos(pi*yg) - cos(3.d0*pi*yg)*cos(pi*xg))/(gamlg-1.d0)
!
!    src = 0.5d0*pi/(gamlg-1.d0)*(sin(2.d0*pi*yg)*cos(pi*xg)*sin(pi*yg) - sin(2.d0*pi*xg)*sin(pi*xg)*cos(pi*yg))
!
!finally, scatter the contribution to the RHS
!
! if(ie==2) print*,'src rhs', 3.5d0/9.d0,rhsel(1, 4, 2),coorp(1, 1:3), coorp(2, 1:3)
!
do ideg = 1,ndegr
rhsel(ideg,4,ielem)=rhsel(ideg,4,ielem) + src*b(ideg)*djaco
enddo
!
!if(ie==2) print*,'src rhs', rhsel(1, 4, 2),src,b(1),djaco,src*b(1)*djaco,xg,yg
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
end subroutine rhsdomnsrcdg_lag_mc_hybridquad
!
!...source domain integral for hybrid linear triangle cells
!
subroutine rhsdomnsrcdg_lagmc_triarz(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
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
integer,dimension(1:nvtri) :: ipt
!...local real array
real*8,dimension(1:ndimn, 1:nvtri) :: xp, xpi
real*8,dimension(1:ndegr):: b, dbdr, dbds,bv
real*8:: unknod(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8::weigh(ngausd), posi(2,ngausd)
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
real*8::xgaus, ygaus
real*8::rhoad,uadv,vadv,eadv,rhoma
real*8::pres
real*8::djaco, wi, rcoef
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
call rutope(2, ngausd, posi, weigh)
!
!...bv denotes the shape function based on centroid
!
!...Loop over elements
!
do 550 ie = 1,ntria!...(1)ie = 1,nelem
!
ielem = ie
!
!...Points consitituting one element...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
!
xp(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xp(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
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
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
wi     = weigh(ig)
!
dspr(1) = -c10
dspr(2) =  c10
dspr(3) =  0.d0
!
dsps(1) = -c10
dsps(2) =  0.d0
dsps(3) =  c10
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
if(ndens.eq.1)then
!
rhoma = unknod(1)
rhoad  = 1.d0/rhoma
elseif(ndens.eq.2)then
!
xpi(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpi(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
call getrhoig_tria(rhoi, r, s, xpi)
call getdensity_triallnl(r, s, xp, xpi, rhoi, rhon)
!
rhoma = 1.d0/rhon
rhoad = rhon
!
elseif(ndens.eq.3)then
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
endif

!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!...impose imiter...
!
if(nlimi.eq.1)then
!
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
rhoma = rhomc + aflim(1 ,ielem)*(rhoma - rhomc)
rhoad = 1.d0/rhoma
!
uadv = uctr + aflim(2, ielem)*(uadv - uctr)
vadv = vctr + aflim(3 ,ielem)*(vadv - vctr)
!
pres = pctr + aflim(4, ielem)*(pres- pctr)
!
elseif(nlimi.eq.6)then
!
! if(ie.ge.990.and.ie.le.1000) print*,'domn unk',ie,unkno(1, 1, ie)
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
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

endif
!
fluxd(1,2) = 0.d0
fluxd(2,2) = 0.d0
fluxd(3,2) = 0.d0
!
fluxd(1,3) = alfrz*(-pres)*b(1)
fluxd(2,3) = alfrz*(-pres)*b(2)
fluxd(3,3) = alfrz*(-pres)*b(3)
!
fluxd(1,4) = (gdshp(1, 1)*uadv + gdshp(2, 1)*vadv)*(-pres)
fluxd(2,4) = (gdshp(1, 2)*uadv + gdshp(2, 2)*vadv)*(-pres)
fluxd(3,4) = (gdshp(1, 3)*uadv + gdshp(2, 3)*vadv)*(-pres)
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg,2:3,ielem)=rhsel(ideg,2:3,ielem) - fluxd(ideg,2:3)*djaco
enddo
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomnsrcdg_lagmc_triarz
!
!....source domain integral for hybrid linear quad cells...
!
subroutine rhsdomnsrcdg_lagmc_quadrz(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )
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
real*8::xgaus, ygaus
real*8::rhoad,uadv,vadv,eadv,rhoma
real*8::pres
real*8::djaco, wi, rcoef
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
!  open(12,file='solution2.dat')
!    do ie=1,ncell
!       if(ie.ge.990.and.ie.le.1000)print*,'afterdomn',ie,unkno(1:ndegr,1,ie)
!    enddo
!  close(12)
!
!...Loop over elements
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:4,ie)!
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
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
!
dsprq(1) = -0.25d0   *sm
dsprq(2) =  0.25d0   *sm
dsprq(3) =  0.25d0   *sp
dsprq(4) = -0.25d0   *sp
!
dspsq(1) =-0.25d0*rm
dspsq(2) =-0.25d0*rp
dspsq(3) = 0.25d0*rp
dspsq(4) = 0.25d0*rm
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 4
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
if(ndens.eq.1)then
!
rhoma = unknod(1)
rhoad  = 1.d0/rhoma
elseif(ndens.eq.2)then
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rhoma = 1.d0/rhon
rhoad = rhon
!
elseif(ndens.eq.3)then
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
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
if(nlimi.eq.1)then
!
!...Added in future
!
elseif(nlimi.eq.6)then
!
! if(ie.ge.990.and.ie.le.1000) print*,'domn unk',ie,unkno(1, 1, ie)
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
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
pres = pctr + aflim(4, ielem)*(pres- pctr)
!
endif
!
fluxd(1,2) = 0.d0
fluxd(2,2) = 0.d0
fluxd(3,2) = 0.d0
!
fluxd(1,3) = alfrz*(-pres)*b(1)
fluxd(2,3) = alfrz*(-pres)*b(2)
fluxd(3,3) = alfrz*(-pres)*b(3)
!
!finally, scatter the contribution to the RHS for mementum equation
!
do ideg = 1,ndegr
rhsel(ideg,2:3,ielem)=rhsel(ideg,2:3,ielem) - fluxd(ideg,2:3)*djaco
enddo
!
enddo !...(2)ig = 1,ngausd
!
650 enddo

end subroutine rhsdomnsrcdg_lagmc_quadrz

!
!...Subroutine to calcualte the maximum and minimum nodal unknows
!for barth limiter based on vertex for primitive variables....
!
subroutine barthlimit_lag_vtxunk(unkno, iptri, ipqua, unmax, unmin)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
integer*4,dimension(1:nvtri,1:ntria),          intent(in)::iptri
integer,  dimension(1:nvqua,1:nquad),          intent(in):: ipqua
!
!...Local
!
integer::ipt(nvtri), ipq(nvqua)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa
integer:: ielem
real*8:: unctr(1:nq+1)
real*8:: unmax(1:nq+2, 1:npoin), unmin(1:nq+2, 1:npoin)
!
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rhov, rhoct, rhom
real*8:: afbar
!
eps = 1.e-6
!
! if(ie==1) print*,unkno(1, 1:4, ie), unknv(1:4, 1)
!
!...Part 1: Get the maximum and minimum at the vertex...
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!...loop over triangles' vertices..
!
do ie = 1, ntria
!
ielem = ie
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
if(ndens.eq.1)then
rhom = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhom = 1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhom = 1.d0/unkno(1, 1, ielem)
endif
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
if(ndens.eq.1)then
unctr(1)   = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
unctr(1)   = unkno(1, 1, ielem)
elseif(ndens.eq.3)then
unctr(1)   = unkno(1, 1, ielem)
endif
!
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq) = pctr
unctr(nq+1) = ectr
!
!
do iv = 1, nvtri
do iq = 1, nq+1
!
if(unctr(iq).gt.unmax(iq, ipt(iv))) then
unmax(iq, ipt(iv)) = unctr(iq)
endif
!
if(unctr(iq).lt.unmin(iq, ipt(iv))) then
unmin(iq, ipt(iv)) = unctr(iq)
endif
!
enddo
enddo
enddo
!
!...loop over quads' vertices..
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
if(ndens.eq.1)then
rhom = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhom = 1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhom = 1.d0/unkno(1, 1, ielem)
endif

uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
if(ndens.eq.1)then
unctr(1)   = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
unctr(1)   = unkno(1, 1, ielem)
elseif(ndens.eq.3)then
unctr(1)   = unkno(1, 1, ielem)
endif

unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq) = pctr
unctr(nq+1) = ectr
!
do iv = 1, nvqua
do iq = 1, nq+1
!
if(unctr(iq).gt.unmax(iq, ipq(iv))) then
unmax(iq, ipq(iv)) = unctr(iq)
endif
!
if(unctr(iq).lt.unmin(iq, ipq(iv))) then
unmin(iq, ipq(iv)) = unctr(iq)
endif
!
enddo
enddo
!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lag_vtxunk

!
!...Subroutine for barth limiter based on vertex for primitive variables on triangles....
!
subroutine barthlimit_lag_vtx_prim_hybridtria(geoel, coord, ustar, unkno, iptri, intfac, aflim, unmax, unmin)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:nq+1, 1:nsize),             intent(out)::aflim
integer*4,dimension(1:nvtri,1:ntria),          intent(in)::iptri
integer, dimension(nifai,nafac),               intent(in)::intfac
real*8, dimension(1:nq+2, 1:npoin),           intent(inout):: unmax, unmin
!
!...Local
!
integer::ipt(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ielem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvtri)::alfa
real*8:: xv(nvtri), yv(nvtri)
real*8:: bt(1:ndegr, 1:nvtri)
real*8:: unmax_new(1:nq+1, 1:nelem), unmin_new(1:nq+1, 1:nelem)
real*8,dimension(1:nq+1,  1:nvtri) ::unknvt
real*8:: dunk(1:nq+1)
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct, rhom
real*8:: afbar
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
!
!...shape functions
!
dr = .5d0
ds = .5d0
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
! if(ie==1) print*,unkno(1, 1:4, ie), unknv(1:4, 1)!
!...Get the minimum and maximum at one cell...
!
do ie = 1, ntria
!
ielem = ie
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
do iq=1, nq+1
!
unmax_new(iq, ielem) = maxval(unmax(iq, ipt(1:nvtri)))
unmin_new(iq, ielem) = minval(unmin(iq, ipt(1:nvtri)))
enddo
enddo
!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Part 2: Impose limiter
!
do ie = 1, ntria
!
ielem = ie
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...zero out unknv
!
unknvt = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
!
enddo
!
rhov = 1.d0/unknvt(1, iv)
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
evtx = unknvt(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unknvt(1, iv) = 1.d0/rhov
unknvt(4 ,iv) = pvtx
!
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
! if(ipt(1).eq.1)then

!  unmax(1:nq, 1) = unkno(1, 1:4, ielem)
!  unmin(1:nq, 1) = unkno(1, 1:4, ielem)
! endif
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq)  = pctr
!
do iv = 1, nvtri
do iq = 1, nq
!
dunk(iq) = unknvt(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!
call barthfct(unmax(iq, ipt(iv)), unmin(iq, ipt(iv)), unctr(iq), dunk(iq), afbar)
!
alfa(iq, iv) = afbar
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
if(indbd(ipt(iv)).eq.1)then

if(ncase .eq. 2)then
!
alfa(2, iv) = max(min((ustar(1, ipt(iv)) - unctr(2))/(dunk(2)),1.d0),0.d0)
alfa(3, iv) = max(min((ustar(2, ipt(iv)) - unctr(3))/(dunk(3)),1.d0),0.d0)
!
elseif(ncase.eq.1)then
if(coord(1, ipt(iv)).lt.1.d-6.or.abs(coord(1, ipt(iv))-1.d0).lt.1.d-6) then
alfa(2, iv) = min(max((- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ipt(iv)).lt.1.d-6.or.abs(coord(2, ipt(iv))-1.d0).lt.1.d-6) then
alfa(3, iv) = min(max(( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
elseif(ncase.eq.3)then
!
if(coord(1, ipt(iv)).lt.1.d-6) then
alfa(2, iv) = min(max(0.5d0*(- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ipt(iv)).lt.1.d-6) then
alfa(3, iv) = min(max(0.5d0*( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
if(coord(2, ipt(iv)).gt.1.d-6) then
!    alfa(4, iv) = min(max(0.5d0*(1.d-6- unctr(4))/(dunk(4)),0.d0),1.d0)
endif
!
elseif(ncase.eq.7)then
!
if(coord(1, ipt(iv)).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max(0.3d0*(- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ipt(iv)).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
alfa(3, iv) = min(max(0.3d0*( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
endif
!
endif
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 1,nq
aflim(iq, ielem) = minval(alfa(iq, 1:nvtri))
enddo
!
enddo
!
!...Part 3: Correct total energy
!
do ie = 1,ntria
!
ielem = ie
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq) = pctr
unctr(nq+1) = ectr
!
!...zero out unknv
!
unknvt = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
!
enddo
!
rhov = 1.d0/unknvt(1, iv)
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
evtx = unknvt(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unknvt(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
unknvt(2, iv) = uctr + aflim(2, ielem)*(uvtx - uctr)
unknvt(3, iv) = vctr + aflim(3, ielem)*(vvtx - vctr)
unknvt(4 ,iv) = pctr + aflim(4, ielem)*(pvtx - pctr)
!
unknvt(5, iv) = unknvt(4 ,iv)/(gamlg-1.d0)*unknvt(1, iv) + 0.5d0*(unknvt(2, iv)**2 + unknvt(3, iv)**2)
!
enddo
!
!
do iv = 1, nvtri
do iq = nq+1, nq+1 !...energy
!
dunk(iq) = unknvt(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!
call barthfct(unmax(iq, ipt(iv)), unmin(iq, ipt(iv)), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
enddo
!!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = nq+1,nq+1
aflim(iq, ielem) = minval(alfa(iq, 1:nvtri))
enddo
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!!
!aflim = 0.d0
!
do ie = 1, ntria
!
!aflim(1:nq+1, ie) = aflim(1:nq+1, )
!
!unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
!unkno(2:3, 2, ie) = unkno(2:3, 2, ie)*aflim(2,ie)
!unkno(2:3, 3, ie) = unkno(2:3, 3, ie)*aflim(3,ie)
!unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(5,ie)
!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lag_vtx_prim_hybridtria
!
!...Subroutine for barth limiter based on vertex for primitive variables on quads....
!
subroutine barthlimit_lag_vtx_prim_hybridquad(geoel, coord, ustar, unkno, ipqua, intfac, aflim, unmax, unmin)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:nq+1, 1:nsize),             intent(out)::aflim
integer,  dimension(1:nvqua,1:nquad),          intent(in):: ipqua
integer, dimension(nifai,nafac),               intent(in)::intfac
real*8, dimension(1:nq+2, 1:npoin),           intent(in):: unmax, unmin
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
real*8:: bq(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
!real*8:: unmax_new(1:nq+1, 1:nelem), unmin_new(1:nq+1, 1:nelem)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct, rhom
real*8:: afbar
!
!...Coloring the boundary node
!
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
!
!...shape functions
!
dr = 1.d0
ds = 1.d0
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Part 2: Impose limiter
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
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
rhov = 1.d0/unknvq(1, iv)
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unknvq(1, iv) = 1.d0/rhov
unknvq(4 ,iv) = pvtx
!
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
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
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!
!if(iq.eq.1)then
!call barthfctrho(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
!else
call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
!endif
!
alfa(iq, iv) = afbar
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
if(indbd(ipq(iv)).eq.1)then

if(ncase .eq. 2)then
!
alfa(2, iv) = max(min((ustar(1, ipq(iv)) - unctr(2))/(dunk(2)),1.d0),0.d0)
alfa(3, iv) = max(min((ustar(2, ipq(iv)) - unctr(3))/(dunk(3)),1.d0),0.d0)
!
elseif(ncase.eq.1)then
if(coord(1, ipq(iv)).lt.1.d-6.or.abs(coord(1, ipq(iv))-1.d0).lt.1.d-6) then
!
fiy = 0.5d0*(- unctr(2))/(dunk(2))
 alfa(2, iv) = min(max((- unctr(2))/(dunk(2)),0.d0),1.d0)
!alfa(2, iv) =  (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)

endif
!
if(coord(2, ipq(iv)).lt.1.d-6.or.abs(coord(2, ipq(iv))-1.d0).lt.1.d-6) then
!
 fiy = 0.5d0*( - unctr(3))/(dunk(3))
 alfa(3, iv) = min(max(( - unctr(3))/(dunk(3)),0.d0),1.d0)
!alfa(3, iv) = (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)

endif
!
elseif(ncase.eq.-3)then
!
if(coord(1, ipq(iv)).lt.1.d-6) then
alfa(2, iv) = min(max(0.4d0*(- unctr(2))/(dunk(2)),0.d0),1.d0)
endif
!
if(coord(2, ipq(iv)).lt.1.d-6) then
alfa(3, iv) = min(max(0.4d0*( - unctr(3))/(dunk(3)),0.d0),1.d0)
endif
!
if(coord(2, ipq(iv)).gt.1.d-6) then
!    alfa(4, iv) = min(max(0.5d0*(1.d-6- unctr(4))/(dunk(4)),0.d0),1.d0)
endif
!
elseif(ncase.eq.-7)then
!
if(coord(1, ipq(iv)).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
!alfa(1, iv) = 0.8d0
fiy = .6d0*(- unctr(2))/(dunk(2))
alfa(2, iv) = min(max(1.0d0*(- unctr(2))/(dunk(2)),0.d0),1.d0)
!alfa(2, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
endif
!
if(coord(2, ipq(iv)).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
!alfa(1, iv) = 0.8d0
fiy = .6d0*( - unctr(3))/(dunk(3))
alfa(3, iv) = min(max(1.0d0*( - unctr(3))/(dunk(3)),0.d0),1.d0)
!alfa(3, iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
endif
!
endif
!
endif
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 1,nq
aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
enddo
!
!aflim(2, ielem) = minval(aflim(2:3, ielem))
!aflim(3, ielem) = aflim(2, ielem)
!
enddo
!
!...Part 3: Correct total energy
!
do ie = 1,nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
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
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
rhom = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
unctr(1)   = 1.d0/rhoct
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq) = pctr
unctr(nq+1) = ectr
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
rhov = 1.d0/unknvq(1, iv)
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unknvq(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
unknvq(2, iv) = uctr + aflim(2, ielem)*(uvtx - uctr)
unknvq(3, iv) = vctr + aflim(3, ielem)*(vvtx - vctr)
unknvq(4 ,iv) = pctr + aflim(4, ielem)*(pvtx - pctr)
!
unknvq(5, iv) = unknvq(4 ,iv)/(gamlg-1.d0)*unknvq(1, iv) + 0.5d0*(unknvq(2, iv)**2 + unknvq(3, iv)**2)
!
enddo
!
!
do iv = 1, nvqua
do iq = nq+1, nq+1 !...energy
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!
call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
enddo
!!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = nq+1,nq+1
aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
enddo
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
!aflim = 1.d0
!!
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lag_vtx_prim_hybridquad
!
!
!....domain integral for hybrid linear cells
!
subroutine getdensity_triallnl(r, s, xpt, xpti, rhoi, rhon)
use constant
implicit none
!...Input arrays
!
real*8,dimension(1:ndimn, 1:nvtri), intent(in) :: xpt, xpti
real*8,intent(in):: r, s, rhoi
real*8,intent(out)::rhon
!...Local integer
!
integer::ishp
!
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8::dxdr, dxds, dydr, dyds, jacoi, jacon
real*8::eps, c00, c10, c05,c20
!...local real!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
!
dspr(1) = -c10
dspr(2) =  c10
dspr(3) =  0.d0
!
dsps(1) = -c10
dsps(2) =  0.d0
dsps(3) =  c10
!
!...New Jacobian
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*xpt(1,ishp)
dxds = dxds + dsps(ishp)*xpt(1,ishp)

dydr = dydr + dspr(ishp)*xpt(2,ishp)
dyds = dyds + dsps(ishp)*xpt(2,ishp)
enddo
!
jacon = dxdr*dyds - dxds*dydr
!
!...Initial Jacobian
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dspr(ishp)*xpti(1,ishp)
dxds = dxds + dsps(ishp)*xpti(1,ishp)

dydr = dydr + dspr(ishp)*xpti(2,ishp)
dyds = dyds + dsps(ishp)*xpti(2,ishp)
enddo
!
jacoi = dxdr*dyds - dxds*dydr
!
!...Get the new density at the (r,s) in reference cell....
!
rhon = rhoi*jacoi/jacon
!
end subroutine getdensity_triallnl
!
!....Density update for quad using llnl...
!
subroutine getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
use constant
implicit none
!...Input arrays
!
real*8,dimension(1:ndimn, 1:nvqua), intent(in) :: xpq, xpqi
real*8,intent(in):: r, s, rhoi
real*8,intent(out)::rhon
!...Local integer
!
integer::ishp
!
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::dxdr, dxds, dydr, dyds, jacoi, jacon
real*8::rp, rm, sp, sm
real*8::eps, c00, c10, c05,c20
!...local real!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Shape functions and their derivatives...
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
!
dsprq(1) = -0.25d0   *sm
dsprq(2) =  0.25d0   *sm
dsprq(3) =  0.25d0   *sp
dsprq(4) = -0.25d0   *sp
!
dspsq(1) =-0.25d0*rm
dspsq(2) =-0.25d0*rp
dspsq(3) = 0.25d0*rp
dspsq(4) = 0.25d0*rm
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 4
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
jacon = dxdr*dyds - dxds*dydr
!
!...Initial Jacobian
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 4
dxdr = dxdr + dsprq(ishp)*xpqi(1,ishp)
dxds = dxds + dspsq(ishp)*xpqi(1,ishp)

dydr = dydr + dsprq(ishp)*xpqi(2,ishp)
dyds = dyds + dspsq(ishp)*xpqi(2,ishp)
enddo
!
jacoi = dxdr*dyds - dxds*dydr
!
!...Get the new density at the (r,s) in reference cell....
!
rhon = rhoi*jacoi/jacon
!if(jacoi/jacon.lt..5d0)print*,'jaco',jacoi/jacon
!
!
end subroutine getdensity_quadllnl
!
!...Get initial unkno distribution at gauss points...
!
subroutine getrhoig_tria(rhoi, r, s, xpti)
use constant
implicit none
real*8, intent(out)::rhoi
real*8, intent(in) ::r, s
real*8,dimension(1:ndimn, 1:nvtri), intent(in) :: xpti
!
real*8::shp(nvtri)
real*8:: xg,yg,xc,yc
real*8:: rhoini
integer::ishp
real*8::radie, radii,radie2,radii2,radic2, radig2,sentr,rhoin,rhoex,rho0ba
real*8::preex,prein
!
!...Shape functions and their derivatives...
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
!
xg = 0.d0
yg = 0.d0
!
do ishp = 1, nvtri
xg = xg + shp(ishp)*xpti(1,ishp)
yg = yg + shp(ishp)*xpti(2,ishp)
enddo
!
if(ncase.eq.1)then     !...TGV
!
rhoini = 1.d0
!
elseif(ncase.eq.2)then !...Shockless Noh problem...
!
rhoini = 1.d0
!
elseif(ncase.eq.3)then !...
!
rhoini  = 1.d0
!
! print*,'pini',pini
!
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
!rhoin = 6.31d-4
!sentr = 2.15d4
!
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radig2 = xg**2 + yg**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rhoini = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rhoini= 2.d0*exp(-(xg**2+yg**2))
!
elseif(ncase.eq.6)then !...Sod...!
!
shp(1) = 1.d0/3.d0
shp(2) = 1.d0/3.d0
shp(3) = 1.d0/3.d0
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nvtri
xc = xc + shp(ishp)*xpti(1,ishp)
yc = yc + shp(ishp)*xpti(2,ishp)
enddo
!
!if(xc.le.50.d0)then
if(sqrt(xc**2+(yc-0.d0)**2).le.0.5d0)then
rhoini  = 1.d0
else
rhoini  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
!
rhoini  = 1.d0
!
elseif(ncase.eq.8)then !...Gresho...!
!
rhoini  = 1.d0
!
elseif(ncase.eq.9)then !...Triple-point...!
!
shp(1) = 1.d0/3.d0
shp(2) = 1.d0/3.d0
shp(3) = 1.d0/3.d0
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, nvtri
xc = xc + shp(ishp)*xpti(1,ishp)
yc = yc + shp(ishp)*xpti(2,ishp)
enddo
!
if(xc.lt.1.d0)then
rhoini  = 1.d0
else
if(yc.gt.1.5d0)then
rhoini  = 0.1d0
else
rhoini  = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion
!
rhoini  = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rhoini  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
rhoini = 1.d0 + 0.9999995d0*sin(pi*xg)

elseif(ncase.eq.13)then !...Saltzman
rhoini  = 1.d0

elseif(ncase.eq.14)then !...Coggeshall
rhoini  = 1.d0

else

print*,'Please specify the initital density ditribution in subroutine getrhoig_tria for Subgrid method for linear Trias!'
stop

endif

!
rhoi = rhoini

end subroutine getrhoig_tria
!
!...Get initial unkno distribution at gauss points...
!
subroutine getrhoig_quad(rhoi, r, s, xpqi)
use constant
implicit none
real*8, intent(out)::rhoi
real*8, intent(in) ::r, s
real*8,dimension(1:ndimn, 1:4), intent(in) :: xpqi
!
real*8::shpq(4)
real*8:: xc,yc,xg,yg
real*8::rp, rm, sp, sm
real*8::eps, c00, c10, c05,c20
real*8:: rhoini
real*8::radie, radii,radie2,radii2,radig2,sentr,rhoin,rhoex,rho0ba
real*8::preex,prein
!
integer::ishp
!...local real!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Shape functions and their derivatives...
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
xg = 0.d0
yg = 0.d0
!
do ishp = 1, 4
xg = xg + shpq(ishp)*xpqi(1,ishp)
yg = yg + shpq(ishp)*xpqi(2,ishp)
enddo
!
if(ncase.eq.1)then     !...TGV
!
rhoini = 1.d0
!
elseif(ncase.eq.2)then !...Shockless Noh problem...
!
rhoini = 1.d0
!
elseif(ncase.eq.3)then !...
!
rhoini  = 1.d0
!
! print*,'pini',pini
!
elseif(ncase.eq.4)then !... Kidder shell...
!
radie = 1.0d0
radii = 0.9d0
prein = 0.1d0
preex = 10.d0
rhoex = 1.d-2
!rhoin = 6.31d-4
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radig2 = xg**2 + yg**2
!
rho0ba = (radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rhoini = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rhoini= 2.d0*exp(-(xg**2+yg**2))
!
elseif(ncase.eq.6)then !...Sod...!
!
!...  Mass center
!
shpq(1) = 0.25d0
shpq(2) = 0.25d0
shpq(3) = 0.25d0
shpq(4) = 0.25d0
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, 4
xc = xc + shpq(ishp)*xpqi(1,ishp)
yc = yc + shpq(ishp)*xpqi(2,ishp)
enddo
!
! print*,'bad sg+yg',sqrt(xg**2+yg**2)
!
!if(xc.le.50.d0)then
if(sqrt(xc**2+(yc-0.0d0)**2).le.0.5d0)then
rhoini  = 1.d0
else
rhoini  = 0.125d0
endif
!
elseif(ncase.eq.7)then !...sedov
!
rhoini  = 1.d0
!
elseif(ncase.eq.8)then !...Gresho...!
!
rhoini  = 1.d0
!
elseif(ncase.eq.9)then !...Triple-point...!
!
!...  Mass center
!
shpq(1) = 0.25d0
shpq(2) = 0.25d0
shpq(3) = 0.25d0
shpq(4) = 0.25d0
!
xc = 0.d0
yc = 0.d0
!
do ishp = 1, 4
xc = xc + shpq(ishp)*xpqi(1,ishp)
yc = yc + shpq(ishp)*xpqi(2,ishp)
enddo
!
if(xc.lt.1.d0)then
rhoini  = 1.d0
else
if(yc.gt.1.5d0)then
rhoini  = 0.1d0
else
rhoini  = 1.d0
endif
endif
!
elseif(ncase.eq.10)then !...Expansion
!
rhoini  = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rhoini  = 1.d0
!
elseif(ncase.eq.12)then !...1D isentropic sin wave...
!
rhoini = 1.d0 + 0.9999995d0*sin(pi*xg)
!
elseif(ncase.eq.13)then !...Saltzman
rhoini  = 1.d0

elseif(ncase.eq.14)then !...Coggeshall
rhoini  = 1.d0
!
else

print*,'Please specify the initital density ditribution in subroutine getrhoig_tria for Subgrid method for linear Quads!'
stop
endif

rhoi = rhoini
end subroutine getrhoig_quad
!
!...subroutine: Calculate the nodal velocity U_p^* (mass center) for hybrid RZ mesh with general Riemann solver...
!
subroutine getndvelo_lagmc_rz(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold, unkno,ustar, fstrz, fsqrz, aflim, afvec, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3, 1:ngelg, 1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri, 1:ntria),  intent(out)::fstrz !...Riemann forces
real*8,dimension(1:ndimn,1:4,1:nvqua, 1:nquad),  intent(out)::fsqrz !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec

integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop, ivtest
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)

!...local real array
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8::munaci(2, 2)
real*8::routb(2),routbx(2),routby(2)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2
real*8::dtime
real*8::fverad,veradi
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
!
ipf(1:2) = intfac(3:4, ifa)
!
indnd(ipf(1:2)) = 1

enddo
endif
!
!
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
!
!     print*,'ifa', ifa,ipf
indnd(ipf(1:nvfac)) = 1
endif
enddo
!
!...Get averaged velocity at nodes...
!
call getvlavenew(iptri, ipqua, geoel, vlave, unkno, aflim, afvec)
!
!
do 950 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
if(bface(3,ifa).eq.22)then
ipf(1:2) = intfac(3:4, ifa)
if(bface(4,ifa).eq.221)then
vlave(2,ipf(1:2)) = 0.d0
elseif(bface(4,ifa).eq.222)then
vlave(1,ipf(1:2)) = 0.d0
endif
elseif(bface(3,ifa).eq.25)then
vlave(2,ipf(1:2)) = 0.d0
vlave(1,ipf(1:2)) = 0.d0

endif!
950 enddo
!
!...Use the nodal velocity from last time step...
!
do iloop= 1, 1
!
vlave= ustar
!
!vlave(:,1) = 0.d0
!
!
!print*,'average vlave limitera',aflim(1:5,1118), aflim(1:5, 1119),aflim(1:5, 1120),aflim(1:5, 1121)
!print*,'average vlave limiter',aflim(1:5,1070), aflim(1:5, 1071),aflim(1:5, 1072),aflim(1:5, 1073)
!print*,'limi6'
!
!
!...Zero out munacn
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
if(ntria.gt.0) call getriem_triarz(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold,aflim, afvec)
!
if(nquad.gt.0) call getriem_quadrz(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
!call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)!
!call getboundary_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm)
call getbc_lagmaire2(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)
!
!...Update the nodal velocity...
!
!do ipoin = 1, npoin
!if(indnd(ipoin).eq.0)then
!ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) - fpres(1, ipoin))/munacn(ipoin)
!ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) - fpres(2, ipoin))/munacn(ipoin)
!endif
!enddo
!
!...4.1: Update the Riemann forces at every node...
!
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
rhsu1 = munacu(1, ipoin) - snsigm(1, ipoin) - fpres(1, ipoin)
rhsu2 = munacu(2, ipoin) - snsigm(2, ipoin) - fpres(2, ipoin)
!
ustar(1, ipoin) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ustar(2, ipoin) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!
!if(ipoin.eq.159) print*,ustar(1:2,ipoin),detma,munacn(:,:,ipoin),munacu(1:2, ipoin),snsigm(1:2, ipoin)
endif
enddo
!
!....Bd velocity
!
!print*,'ustar--',sqrt(unkno(1,2,1:8)**2 + unkno(1,3,1:8)**2  ),unkno(1,1,1:8),unkno(1,4,1:8)
! print*,'ustar--',snsigm(1:2, 738)/sqrt(snsigm(1, 738)**2+snsigm(2, 738)**2),&
!                  munacu(1:2, 738)/sqrt(munacu(1, 738)**2+munacu(2, 738)**2)
!do ivtest = 2,5
!print*,'ustar--',ivtest, ustar(1:2,ivtest),sqrt(ustar(1,ivtest)**2 + ustar(2,ivtest)**2 ),&
!(-0.30972423147772277d0*5.9925050412607754d-002+2.5686684725619090d-004*2.5810971608771415d-004)/1.8470801575209576d-002,&

!(-0.13706546050947871d0*3.2517374968819825d-002-0.11562376925763779d0*5.0343652505702170d-002)/1.8470801575209576d-002

!print*,'ustar--',ustar(:,1)
!enddo
!
ustar(:,1) = 0.d0
!
do 900 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
if(bface(3,ifa).eq.22)then
ipf(1:2) = intfac(3:4, ifa)
if(bface(4,ifa).eq.221)then
!ustar(2,ipf(1:2)) = 0.d0
elseif(bface(4,ifa).eq.222)then
ustar(1,ipf(1:2)) = 0.d0
endif
elseif(bface(3,ifa).eq.25)then
ustar(2,ipf(1:2)) = 0.d0
ustar(1,ipf(1:2)) = 0.d0

!...Specify boundary velocity
if(ncase.eq.14)then !...Coggeshall expansion problem
dtime = (itime-1.d0)*dtfix
ustar(1, ipf(1)) = -coord(1, ipf(1))/4.d0/(1.d0-dtime);   ustar(2, ipf(1)) = -coord(2, ipf(1))/(1.d0-dtime)
ustar(1, ipf(2)) = -coord(1, ipf(2))/4.d0/(1.d0-dtime);   ustar(2, ipf(2)) = -coord(2, ipf(2))/(1.d0-dtime)

elseif(ncase.eq.11)then !...Self-similar implosion
routb(1) = sqrt(coord(2,ipf(1))**2+coord(1,ipf(1))**2)
routb(2) = sqrt(coord(2,ipf(2))**2+coord(1,ipf(2))**2)
!
routbx(1) = coord(1,ipf(1))/routb(1);routby(1)=coord(2,ipf(1))/routb(1)
routbx(2) = coord(1,ipf(2))/routb(2);routby(2)=coord(2,ipf(2))/routb(2)
!
fverad = 1.d0-0.185d0*(itime-1.d0)*dtfix-0.28d0*((itime-1.d0)*dtfix)**3
veradi = -0.6883545d0*fverad/(1.d0-fverad*(itime-1.d0)*dtfix)**(1.d0-0.6883545d0)
!
ustar(1, ipf(1)) = veradi*routbx(1);   ustar(2, ipf(1)) = veradi*routby(1);
ustar(1, ipf(2)) = veradi*routbx(2);   ustar(2, ipf(2)) = veradi*routby(2);

elseif(ncase.eq.13)then !...Saltzman
ustar(1,ipf(1:2)) = 1.d0
ustar(2,ipf(1:2)) = 0.d0
endif

endif
!
900 enddo
!endif
enddo !iloop
!
!...Imposing the zero normal velocity for BC...
!
!
!...4.2: Update the Riemann forces at every node...
!
do ie = 1, ntria
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
!...shape functions
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv = 1, nvtri
!
do ifa =1, 4
fstrz(1, ifa, iv, ie) = snsigmlt(1, ifa, iv, ie) + &
munaclt(1, 1, ifa, iv, ie)*ustar(1, ipt(iv))+&
munaclt(2, 1, ifa, iv, ie)*ustar(2, ipt(iv)) - munault(1, ifa, iv, ie)
fstrz(2, ifa, iv, ie) = snsigmlt(2, ifa, iv, ie) + &
munaclt(1, 2,ifa, iv, ie)*ustar(1, ipt(iv))+&
munaclt(2, 2, ifa, iv, ie)*ustar(2, ipt(iv)) - munault(2, ifa, iv, ie)
enddo
enddo
enddo
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
do iv = 1, nvqua
!
do ifa =1, 4
fsqrz(1, ifa, iv, ie) = snsigmlq(1, ifa, iv, ie) + &
munaclq(1, 1, ifa, iv, ie)*ustar(1, ipq(iv))+&
munaclq(2, 1, ifa, iv, ie)*ustar(2, ipq(iv)) - munaulq(1, ifa, iv, ie)
fsqrz(2, ifa, iv, ie) = snsigmlq(2, ifa, iv, ie) + &
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
enddo
!
deallocate (munacn, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lagmc_rz
!
!...subroutine: Calculate the Riemann input for hybrid tria grids(RZ) general Riemann solver....
!
subroutine getriem_triarz(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
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
real*8, dimension(1:2, 1:2, 1:4,  1:nvtri, 1:ntria),       intent(out)::munaclt
real*8, dimension(1:ndimn, 1:4, 1:nvtri,  1:ntria), intent(out)::munault
real*8, dimension(1:ndimn, 1:4, 1:nvtri,  1:ntria), intent(out)::snsigmlt
!...Local integer
integer::ie, ideg, ielem, ifa, iv
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:3, 1:nvtri)::bt, btv
real*8,dimension(1:nq,1:nvtri)::unknvt
real*8::aujmp(1:3, 1:nvtri)
real*8::vnorm(1:3, 1:4, 1:nvtri)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8,dimension(1:4, 1:nvtri)::murie
real*8,dimension(1:nvtri):: xv,  yv
real*8,dimension(1:nvtri):: rcoet
real*8,dimension(1:ndimn, 1:nvtri) :: xpt,xpht
real*8,dimension(1:ndimn, 1:nvtri) :: xpti
real*8, dimension(1:nvtri):: shp, dspr, dsps
!...arraies for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8:: dudr, duds, dvdr, dvds
real*8::acnx, acny
real*8::rhoi, rhon, unmag
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Triangle...
!
do 250 ie = 1,ntria !...(1)ie = 1,nelem
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
!...shape functions
!
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
!...Basis function
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...physical coordinate....
!
xpht(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpht(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
!...Coefficient R of RZ or XY system...
!
rcoet(1:nvtri) = 1.d0 - alfrz + alfrz*xpht(2, 1:nvtri)
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 3, 1) = gelag(1:3, 3, ie); vnorm(1:3, 4, 1) = gelag(1:3, 1, ie) !...For point ip(1)

vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 1, ie); vnorm(1:3, 4, 2) = gelag(1:3, 2, ie) !...For point ip(2)

vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 2, ie); vnorm(1:3, 4, 3) = gelag(1:3, 3, ie) !...For point ip(3)
!
!...ndA=0.5d0*vnorm
!
!vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
vnorm(3, 1, 1) = (3.d0*rcoet(1) + rcoet(3))/3.d0*0.5d0*vnorm(3, 1, 1)
vnorm(3, 2, 1) = (3.d0*rcoet(1) + rcoet(2))/3.d0*0.5d0*vnorm(3, 2, 1)
vnorm(3, 3, 1) = (1.d0*rcoet(1) + rcoet(3))/3.d0*0.5d0*vnorm(3, 3, 1)
vnorm(3, 4, 1) = (1.d0*rcoet(1) + rcoet(2))/3.d0*0.5d0*vnorm(3, 4, 1)
!
vnorm(3, 1, 2) = (3.d0*rcoet(2) + rcoet(1))/3.d0*0.5d0*vnorm(3, 1, 2)
vnorm(3, 2, 2) = (3.d0*rcoet(2) + rcoet(3))/3.d0*0.5d0*vnorm(3, 2, 2)
vnorm(3, 3, 2) = (1.d0*rcoet(2) + rcoet(1))/3.d0*0.5d0*vnorm(3, 3, 2)
vnorm(3, 4, 2) = (1.d0*rcoet(2) + rcoet(3))/3.d0*0.5d0*vnorm(3, 4, 2)
!
vnorm(3, 1, 3) = (3.d0*rcoet(3) + rcoet(2))/3.d0*0.5d0*vnorm(3, 1, 3)
vnorm(3, 2, 3) = (3.d0*rcoet(3) + rcoet(1))/3.d0*0.5d0*vnorm(3, 2, 3)
vnorm(3, 3, 3) = (1.d0*rcoet(3) + rcoet(2))/3.d0*0.5d0*vnorm(3, 3, 3)
vnorm(3, 4, 3) = (1.d0*rcoet(3) + rcoet(1))/3.d0*0.5d0*vnorm(3, 4, 3)
!
vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
!
!...now configuration
!
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
!
rhomc = 1.d0/unkno(1, 1, ielem)
!
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...
!...zero out unknv
unknvt = 0.d0
!
do iv   = 1,nvtri
!
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
enddo
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvt(1, iv)
elseif(ndens.eq.2)then
!...now configuration
!
r = xv(iv); s= yv(iv)
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
call  getrhoig_tria(rhoi, r, s, xpti)
!
call getdensity_triallnl(r, s, xpt, xpti, rhoi, rhon)
!
rhovt = rhon
!
elseif(ndens.eq.3)then
!
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
btv(1, iv) = 1.d0
btv(2, iv) = (xv(iv)-rcv)/dr
btv(3, iv) = (yv(iv)-scv)/ds
!
unknvt(1, iv) =0.d0
!
do ideg = 1,mdegr
unknvt(1, iv) = unknvt(1, iv) + unkno(ideg,1,ielem)*btv(ideg, iv)
enddo
!
rhovt = unknvt(1, iv)
!
endif
!
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
evtx = unknvt(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!...Limiter
!
if(nlimi.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvt(1, iv) - rhomc)
rhovt = 1.d0/rhomv
!
uvtx = uctr + aflim(2, ielem)*(unknvt(2, iv) - uctr)
vvtx = vctr + aflim(3, ielem)*(unknvt(3, iv) - vctr)
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!
!...updtae unknv(2:3,:)
unknvt(2, iv) = uvtx
unknvt(3 ,iv) = vvtx
!
elseif(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvt(1, iv) - rhomc)
rhovt = 1.d0/rhomv
elseif(ndens.eq.2)then
rhovt = 1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc)
!
elseif(ndens.eq.3)then
!
rhovt = 1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc)
!
endif
!
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
!
endif
!
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
!
aujmp(1:2, iv) = vlave(1:2, ipt(iv)) - unknvt(2:3, iv)
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
!if(ipt(iv)==7) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, iv),sqrt(acnx**2 + acny**2), unknvt(2:3, iv)&
!,cos(5.d0/24.d0*pi),sin(5.d0/24.d0*pi),vlave(1:2, ipt(iv))/sqrt(vlave(1, ipt(iv))**2 + vlave(2, ipt(iv))**2)
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
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
!if(ielem.eq.1.or.ielem.eq.2) print*,'ielem56',ielem,rhoct*sdctr,rhoct*sdctr*2.d0*0.5d0*&
!gelag(3, 2, ie)*(9.4801840231665494E-003+1.0347275028465509E-002)*cos(pi/32.d0)*cos(pi/32.d0),&
!rhoct*sdctr*2.d0*0.5d0*&
!gelag(3, 2, ie)*(9.4801840231665494E-003+1.0347275028465509E-002)*cos(pi/32.d0),&
!rhoct*sdctr*2.d0*0.5d0*&
!gelag(3, 2, ie)*(-9.4801840231665494E-003+1.0347275028465509E-002)*sin(pi/32.d0)
!
!...Get impedence coefficient...
!
do iv   = 1, nvtri
dux= vlave(1, ipt(iv))-unknvt(2, iv)
duy= vlave(2, ipt(iv))-unknvt(3, iv)
deltu = 10.d0*sqrt(dux**2 + duy**2)
do ifa = 1, 4
deltu = 10.d0*abs(dux*vnorm(1, ifa, iv) + duy*vnorm(2, ifa, iv))
murie(ifa, iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
enddo
enddo
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!
do iv  = 1, nvtri
do ifa = 1, 4 !...Every corner consists of 2 faces...
!
!...Call Riemann solver...
!
call getriecoef_matrixnew(murie(ifa, iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
unknvt(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ifa, iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
!unknvt(2:3, iv), sigma(1:2, 1:2, iv),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
munacn(1:2, 1, ipt(iv)) = munacn(1:2, 1, ipt(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipt(iv)) = munacn(1:2, 2, ipt(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipt(iv)) = munacu(1:2, ipt(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipt(iv)) = snsigm(1:2, ipt(iv)) + snsigm_rie(1:2)!
!
!
! if(ip(iv).eq.15) print*,'p19 muacn(28) prep---',murie(iv), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                 vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
! if(ip(iv).eq.15) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Get the summed stress sum(n*A*sigma)
!
!
! if(ipq(iv).eq.738) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ipq(iv)),ie, ifa,iv
! if(ipq(iv).eq.738) print*,'ustar--',snsigm(1:2, 738)/sqrt(snsigm(1, 738)**2+snsigm(2, 738)**2)
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Local variable...
!
munaclt(1:2, 1, ifa, iv, ie) =  munacn_rie(1:2, 1)
munaclt(1:2, 2, ifa, iv, ie) =  munacn_rie(1:2, 2)
!
munault(1:2, ifa, iv, ie) =  munacu_rie(1:2)
!
snsigmlt(1:2, ifa, iv, ie)=  snsigm_rie(1:2)
!
enddo
enddo
!
250 enddo  !...(1)ie = 1,nelem!

end subroutine getriem_triarz
!
!...subroutine: Calculate the Riemann input for hybrid quad(RZ) grids general Riemann solver....
!
subroutine getriem_quadrz(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
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
real*8,dimension(1:3, 1:nvqua)::bq,bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8::aujmp(1:3, 1:nvqua)
real*8::vnorm(1:3, 1:4, 1:nvqua)
real*8::sigma(1:2, 1:2, 1:nvqua)
real*8,dimension(1:4, 1:nvqua)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xphq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
!...arraies for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon, unmag
!
data eps   / 1.0d-09/
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
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...physical coordinate....
!
xphq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xphq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Coefficient R of RZ or XY system...
!
rcoeq(1:nvqua) = 1.d0 - alfrz + alfrz*abs(xphq(2, 1:nvqua))
!rcoeq(1:nvqua) = alfrz*xphq(2, 1:nvqua)
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelagq(1:3, 4, ie); vnorm(1:3, 2, 1) = gelagq(1:3, 1, ie);
vnorm(1:3, 3, 1) = gelagq(1:3, 4, ie); vnorm(1:3, 4, 1) = gelagq(1:3, 1, ie); !...For point ip(1)

vnorm(1:3, 1, 2) = gelagq(1:3, 1, ie); vnorm(1:3, 2, 2) = gelagq(1:3, 2, ie)
vnorm(1:3, 3, 2) = gelagq(1:3, 1, ie); vnorm(1:3, 4, 2) = gelagq(1:3, 2, ie) !...For point ip(2)

vnorm(1:3, 1, 3) = gelagq(1:3, 2, ie); vnorm(1:3, 2, 3) = gelagq(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelagq(1:3, 2, ie); vnorm(1:3, 4, 3) = gelagq(1:3, 3, ie) !...For point ip(3)

vnorm(1:3, 1, 4) = gelagq(1:3, 3, ie); vnorm(1:3, 2, 4) = gelagq(1:3, 4, ie) !...For point ip(3)
vnorm(1:3, 3, 4) = gelagq(1:3, 3, ie); vnorm(1:3, 4, 4) = gelagq(1:3, 4, ie) !...For point ip(3)

!if(ie.eq.1) print*,'rcoef',rcoeq(1:4)
!
!...ndA=0.5d0*vnorm
!
vnorm(3, 1, 1) = (3.d0*rcoeq(1) + rcoeq(4))/3.d0*0.5d0*vnorm(3, 1, 1)
vnorm(3, 2, 1) = (3.d0*rcoeq(1) + rcoeq(2))/3.d0*0.5d0*vnorm(3, 2, 1)
vnorm(3, 3, 1) = (1.d0*rcoeq(1) + rcoeq(4))/3.d0*0.5d0*vnorm(3, 3, 1)
vnorm(3, 4, 1) = (1.d0*rcoeq(1) + rcoeq(2))/3.d0*0.5d0*vnorm(3, 4, 1)
!
vnorm(3, 1, 2) = (3.d0*rcoeq(2) + rcoeq(1))/3.d0*0.5d0*vnorm(3, 1, 2)
vnorm(3, 2, 2) = (3.d0*rcoeq(2) + rcoeq(3))/3.d0*0.5d0*vnorm(3, 2, 2)
vnorm(3, 3, 2) = (1.d0*rcoeq(2) + rcoeq(1))/3.d0*0.5d0*vnorm(3, 3, 2)
vnorm(3, 4, 2) = (1.d0*rcoeq(2) + rcoeq(3))/3.d0*0.5d0*vnorm(3, 4, 2)
!
vnorm(3, 1, 3) = (3.d0*rcoeq(3) + rcoeq(2))/3.d0*0.5d0*vnorm(3, 1, 3)
vnorm(3, 2, 3) = (3.d0*rcoeq(3) + rcoeq(4))/3.d0*0.5d0*vnorm(3, 2, 3)
vnorm(3, 3, 3) = (1.d0*rcoeq(3) + rcoeq(2))/3.d0*0.5d0*vnorm(3, 3, 3)
vnorm(3, 4, 3) = (1.d0*rcoeq(3) + rcoeq(4))/3.d0*0.5d0*vnorm(3, 4, 3)
!
vnorm(3, 1, 4) = (3.d0*rcoeq(4) + rcoeq(3))/3.d0*0.5d0*vnorm(3, 1, 4)
vnorm(3, 2, 4) = (3.d0*rcoeq(4) + rcoeq(1))/3.d0*0.5d0*vnorm(3, 2, 4)
vnorm(3, 3, 4) = (1.d0*rcoeq(4) + rcoeq(3))/3.d0*0.5d0*vnorm(3, 3, 4)
vnorm(3, 4, 4) = (1.d0*rcoeq(4) + rcoeq(1))/3.d0*0.5d0*vnorm(3, 4, 4)
!
vnorm(3,:,:) = 0.5d0*vnorm(3, : ,:)
!
!...cell averaged value...
!
if(ndens.eq.1)then
!...Specific volume...
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
!
!...now configuration
!
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
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
!print*,'sound speed',sqrt( max( eps,gamlg*pctr/rhoct) )  ,ie, rhoct, uctr,vctr,pctr
!
!...
!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
!
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
!
!...now configuration
!
r = xvq(iv); s= yvq(iv)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rhovt = rhon
!
elseif(ndens.eq.3)then
!
!
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
bqv(1, iv) = 1.d0
bqv(2, iv) = (xvq(iv)-rcv)/dr
bqv(3, iv) = (yvq(iv)-scv)/ds
!
unknvq(1, iv) =0.d0
!
do ideg = 1,mdegr
unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
enddo
!
rhovt  = unknvq(1, iv)
endif

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!if(ipq(iv).eq.180)print*,'ipq',ielem,iv,ipq(iv),rhovt,evtx,0.5d0*(uvtx**2 + vvtx**2),pctr,pvtx
!
!...Limiter
!
if(nlimi.eq.1)then
!
rhomv = rhomc + aflim(1, ielem)*(unknvq(1, iv) - rhomc)
rhovt = 1.d0/rhomv
!
uvtx = uctr + aflim(2, ielem)*(unknvq(2, iv) - uctr)
vvtx = vctr + aflim(3, ielem)*(unknvq(3, iv) - vctr)
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!
!if(ie.ge.2625.and.ie.le.2628)
!print*,'ie26252628',pctr,aflim(4, ielem),pvtx
!
!...updtae unknv(2:3,:)
unknvq(2, iv) = uvtx
unknvq(3 ,iv) = vvtx
!
elseif(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvq(1, iv) - rhomc)
rhovt = 1.d0/rhomv
elseif(ndens.eq.2)then
!rhomv = 1.d0/(1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc) )
!rhovt = 1.d0/rhomv
rhovt = 1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc)
!
elseif(ndens.eq.3)then
!
rhovt = 1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc)
!
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
unknvq(2, iv) = uvtx
unknvq(3 ,iv) = vvtx
!
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
!if(ipq(iv).eq.2) print*,'velocity 8',ie, afvec(1:2,1:2,ie)
!
!...Get the a_c (unit vector)
!
aujmp(1:2, iv) = vlave(1:2, ipq(iv)) - unknvq(2:3, iv)
!if(ie==22) print*,'adjumpxxx22', vlave(1:2, ip(iv)) , unknv(2:3, iv), ip(iv)
acnx = aujmp(1, iv)
acny = aujmp(2, iv)
!
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
!rhoct = 1.d0/unkno(1, 1, ielem)         !...ct denots center of one cell; cn denotes corner of one cell.
!uctr  = unkno(1, 2, ielem)
!vctr  = unkno(1, 3, ielem)
!ectr  = unkno(1, 4, ielem)
!pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
! print*,'sound speed2',sdctr,ie ,rhoct, uctr,vctr,pctr
!
aujmp(3,:)=aujmp(3,:)/sdctr
!
!if(ielem.eq.5.or.ielem.eq.6) print*,'ielem56',ielem,rhoct*sdctr,rhoct*sdctr*2.d0*0.5d0*&
!gelagq(3, 4, ie)*(9.4801840231665494E-003+1.0347275028465509E-002)*cos(pi/32.d0)*cos(pi/32.d0)
!
!...Get impedence coefficient...
!
do iv   = 1, nvqua
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
deltu = sqrt(dux**2 + duy**2)
do ifa = 1, 4
deltu = 1.d0*abs(dux*vnorm(1, ifa, iv) + duy*vnorm(2, ifa, iv))
murie(ifa, iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
enddo
enddo
!
!if(ie==94)print*,'vnotm',murie,rhoct,sdctr,uctr,vctr,ectr!,vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, nvqua
do ifa = 1, 4 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!...Call Riemann solver...
!
call getriecoef_matrixnew(murie(ifa, iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ifa, iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
!unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
munacn(1:2, 1, ipq(iv)) = munacn(1:2, 1, ipq(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipq(iv)) = munacn(1:2, 2, ipq(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipq(iv)) = munacu(1:2, ipq(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipq(iv)) = snsigm(1:2, ipq(iv)) + snsigm_rie(1:2)!
!
!if(ipq(iv).eq.180) print*,'p36 muacn(vv) post',ipq(iv),ielem,ifa,iv,vnorm(1:3, ifa, iv),munacn(1, 1, ipq(iv)),sigma(1,1,iv)
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
!
350 enddo  !...(1)ie = 1,nelem!

end subroutine getriem_quadrz
!
!...Face integral (mass center) for hybrid linear triangle(RZ) using analytical integration...
!
subroutine rhsifacedg_lagmc_triarz(iptri, unkno, ustar, fstrz, gelag, geoel,coord,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri,1:ntria),  intent(in)::fstrz !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),     intent(in)::coord
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:2, 1:nvtri) :: ipf
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:4, 1:nvtri)::lpnpt
real*8::xvt(3), yvt(3),bt(1:3,1:nvtri)
real*8,dimension(1:ndimn, 1:nvtri) :: xpht
real*8,dimension(1:nvtri):: rcoet
real*8,dimension(1:4, 1:nvtri):: rcapt
real*8,dimension(1:3, 1:4, 1:nvtri):: bcapt
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
ipf(1, 1) = 3; ipf(2, 1) = 2
ipf(1, 2) = 1; ipf(2, 2) = 3
ipf(1, 3) = 2; ipf(2, 3) = 1
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
ielem = ie
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
!
do iv =1 ,nvtri
!
!print*,'iv', ie
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds
enddo
!
!...physical coordinate....
!
xpht(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpht(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
!...Coefficient R of RZ or XY system...
!
rcoet(1:nvtri) = 1.d0 - alfrz + alfrz*xpht(2, 1:nvtri)
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpt(1:ndimn, ig, 1, 1) = 0.5d0*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)*bt(ig, 1)
lpnpt(1:ndimn, ig, 2, 1) = 0.5d0*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)*bt(ig, 1)
lpnpt(1:ndimn, ig, 3, 1) = 0.5d0*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)*bt(ig, 3)
lpnpt(1:ndimn, ig, 4, 1) = 0.5d0*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)*bt(ig, 2)
!
!...point 2
lpnpt(1:ndimn, ig, 1, 2) = 0.5d0*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)*bt(ig, 2)
lpnpt(1:ndimn, ig, 2, 2) = 0.5d0*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)*bt(ig, 2)
lpnpt(1:ndimn, ig, 3, 2) = 0.5d0*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)*bt(ig, 1)
lpnpt(1:ndimn, ig, 4, 2) = 0.5d0*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)*bt(ig, 3)
!
!...point 3
lpnpt(1:ndimn, ig, 1, 3) = 0.5d0*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)*bt(ig, 3)
lpnpt(1:ndimn, ig, 2, 3) = 0.5d0*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)*bt(ig, 3)
lpnpt(1:ndimn, ig, 3, 3) = 0.5d0*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)*bt(ig, 2)
lpnpt(1:ndimn, ig, 4, 3) = 0.5d0*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)*bt(ig, 1)
!
enddo
!
!...Uppercase R and B
!
rcapt(1, 1) = (3.d0*rcoet(1) + rcoet(3))/6.d0
rcapt(2, 1) = (3.d0*rcoet(1) + rcoet(2))/6.d0
rcapt(3, 1) = (1.d0*rcoet(1) + rcoet(3))/6.d0
rcapt(4, 1) = (1.d0*rcoet(1) + rcoet(2))/6.d0

rcapt(1, 2) = (3.d0*rcoet(2) + rcoet(1))/6.d0
rcapt(2, 2) = (3.d0*rcoet(2) + rcoet(3))/6.d0
rcapt(3, 2) = (1.d0*rcoet(2) + rcoet(1))/6.d0
rcapt(4, 2) = (1.d0*rcoet(2) + rcoet(3))/6.d0

rcapt(1, 3) = (3.d0*rcoet(3) + rcoet(2))/6.d0
rcapt(2, 3) = (3.d0*rcoet(3) + rcoet(1))/6.d0
rcapt(3, 3) = (1.d0*rcoet(3) + rcoet(2))/6.d0
rcapt(4, 3) = (1.d0*rcoet(3) + rcoet(1))/6.d0
!
bcapt(1:3, 1, 1) = bt(1:3, 1)
bcapt(1:3, 2, 1) = bt(1:3, 1)
bcapt(1:3, 3, 1) = bt(1:3, 3)
bcapt(1:3, 4, 1) = bt(1:3, 2)

bcapt(1:3, 1, 2) = bt(1:3, 2)
bcapt(1:3, 2, 2) = bt(1:3, 2)
bcapt(1:3, 3, 2) = bt(1:3, 1)
bcapt(1:3, 4, 2) = bt(1:3, 3)

bcapt(1:3, 1, 3) = bt(1:3, 3)
bcapt(1:3, 2, 3) = bt(1:3, 3)
bcapt(1:3, 3, 3) = bt(1:3, 2)
bcapt(1:3, 4, 3) = bt(1:3, 1)
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvtri
do ifa =1 ,4
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipt(iv))*lpnpt(1, 1:ndegr, ifa, iv)*rcapt(ifa, iv) +&
ustar(2, ipt(iv))*lpnpt(2, 1:ndegr, ifa, iv)*rcapt(ifa, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstrz(1, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstrz(2, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipt(iv))*fstrz(1, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)  +&
ustar(2, ipt(iv))*fstrz(2, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)
!
enddo
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
! if(ie==18) print*,'rhs iface',rhsel(1, 1, ie), lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

550 enddo
!
!   open(8,file='lpnp.dat')
!    do ie = 1, nelem
!      do iv = 1, nvtri
!
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 1, 1, iv, ie), lpnp(1:2, 1, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 2, 1, iv, ie), lpnp(1:2, 2, 2, iv, ie)
!        write(8,*) ie, iv, inpoel(iv, ie), lpnp(1:2, 3, 1, iv, ie), lpnp(1:2, 3, 2, iv, ie)
!      enddo
!    enddo
!   close(8)
!
end subroutine rhsifacedg_lagmc_triarz

!
!...Face integral (mass center) for hybrid linear quad using analytical integration...
!
subroutine rhsifacedg_lagmc_quadrz(ipqua, unkno, ustar,fsqrz, gelagq, geoel,coord,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvqua, 1:nquad),  intent(in)::fsqrz !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:2, 1:nvqua) :: ipfq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:4, 1:nvqua)::lpnpq
real*8::xvq(nvqua), yvq(nvqua),bq(1:3,1:nvqua)
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:4, 1:nvqua):: rcapt
real*8,dimension(1:3, 1:4, 1:nvqua):: bcapt
real*8,dimension(1:ndimn, 1:nvqua) :: xphq
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
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
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
!
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...physical coordinate....
!
xphq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xphq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Coefficient R of RZ or XY system...
!
rcoeq(1:nvqua) = 1.d0 - alfrz + alfrz*abs(xphq(2, 1:nvqua))
!
!if(ie.eq.1)print*,'rcoef',rcoeq,xphq(2, :)
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpq(1:ndimn, ig, 1, 1) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 1)
lpnpq(1:ndimn, ig, 2, 1) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 1)
lpnpq(1:ndimn, ig, 3, 1) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 4)
lpnpq(1:ndimn, ig, 4, 1) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 2)
!
!...point 2
lpnpq(1:ndimn, ig, 1, 2) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 2)
lpnpq(1:ndimn, ig, 2, 2) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 2)
lpnpq(1:ndimn, ig, 3, 2) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 1)
lpnpq(1:ndimn, ig, 4, 2) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 3)
!
!...point 3
lpnpq(1:ndimn, ig, 1, 3) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 3)
lpnpq(1:ndimn, ig, 2, 3) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 3)
lpnpq(1:ndimn, ig, 3, 3) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 2)
lpnpq(1:ndimn, ig, 4, 3) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 4)
!...point 4
lpnpq(1:ndimn, ig, 1, 4) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 4)
lpnpq(1:ndimn, ig, 2, 4) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 4)
lpnpq(1:ndimn, ig, 3, 4) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 3)
lpnpq(1:ndimn, ig, 4, 4) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 1)
!
enddo
!
!...Capital r
!
rcapt(1, 1) = (3.d0*rcoeq(1) + rcoeq(4))/6.d0
rcapt(2, 1) = (3.d0*rcoeq(1) + rcoeq(2))/6.d0
rcapt(3, 1) = (1.d0*rcoeq(1) + rcoeq(4))/6.d0
rcapt(4, 1) = (1.d0*rcoeq(1) + rcoeq(2))/6.d0

rcapt(1, 2) = (3.d0*rcoeq(2) + rcoeq(1))/6.d0
rcapt(2, 2) = (3.d0*rcoeq(2) + rcoeq(3))/6.d0
rcapt(3, 2) = (1.d0*rcoeq(2) + rcoeq(1))/6.d0
rcapt(4, 2) = (1.d0*rcoeq(2) + rcoeq(3))/6.d0

rcapt(1, 3) = (3.d0*rcoeq(3) + rcoeq(2))/6.d0
rcapt(2, 3) = (3.d0*rcoeq(3) + rcoeq(4))/6.d0
rcapt(3, 3) = (1.d0*rcoeq(3) + rcoeq(2))/6.d0
rcapt(4, 3) = (1.d0*rcoeq(3) + rcoeq(4))/6.d0

rcapt(1, 4) = (3.d0*rcoeq(4) + rcoeq(3))/6.d0
rcapt(2, 4) = (3.d0*rcoeq(4) + rcoeq(1))/6.d0
rcapt(3, 4) = (1.d0*rcoeq(4) + rcoeq(3))/6.d0
rcapt(4, 4) = (1.d0*rcoeq(4) + rcoeq(1))/6.d0
!
bcapt(1:3, 1, 1) = bq(1:3, 1)
bcapt(1:3, 2, 1) = bq(1:3, 1)
bcapt(1:3, 3, 1) = bq(1:3, 4)
bcapt(1:3, 4, 1) = bq(1:3, 2)

bcapt(1:3, 1, 2) = bq(1:3, 2)
bcapt(1:3, 2, 2) = bq(1:3, 2)
bcapt(1:3, 3, 2) = bq(1:3, 1)
bcapt(1:3, 4, 2) = bq(1:3, 3)

bcapt(1:3, 1, 3) = bq(1:3, 3)
bcapt(1:3, 2, 3) = bq(1:3, 3)
bcapt(1:3, 3, 3) = bq(1:3, 2)
bcapt(1:3, 4, 3) = bq(1:3, 4)
!
bcapt(1:3, 1, 4) = bq(1:3, 4)
bcapt(1:3, 2, 4) = bq(1:3, 4)
bcapt(1:3, 3, 4) = bq(1:3, 3)
bcapt(1:3, 4, 4) = bq(1:3, 1)
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvqua
do ifa =1 ,4
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, ifa, iv)*rcapt(ifa, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, ifa, iv)*rcapt(ifa, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fsqrz(1, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fsqrz(2, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(iv))*fsqrz(1, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)  +&
ustar(2, ipq(iv))*fsqrz(2, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)
!
!if(ie==1) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1)  print*,'rhs iface',ielem, ie,ulnpn(1)!,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lagmc_quadrz
!
!...Face integral (mass center) for hybrid linear quad using analytical integration...
!
subroutine rhsifacedg_lagmc_quadrz2(ipqua, unkno, ustar,fsqrz, gelagq, geoel,coord,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvqua, 1:nquad),  intent(in)::fsqrz !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:2, 1:nvqua) :: ipfq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:4, 1:nvqua)::lpnpq
real*8::xvq(nvqua), yvq(nvqua),bq(1:3,1:nvqua)
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:4, 1:nvqua):: rcapt
real*8,dimension(1:3, 1:4, 1:nvqua):: bcapt
real*8,dimension(1:ndimn, 1:nvqua) :: xphq
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
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
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
!
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...physical coordinate....
!
xphq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xphq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Coefficient R of RZ or XY system...
!
rcoeq(1:nvqua) = 1.d0 - alfrz + alfrz*abs(xphq(2, 1:nvqua))
!
!if(ie.eq.1)print*,'rcoef',rcoeq,xphq(2, :)
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpq(1:ndimn, ig, 1, 1) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 1)/2.d0
lpnpq(1:ndimn, ig, 2, 1) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 1)/2.d0
lpnpq(1:ndimn, ig, 3, 1) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 1)/2.d0
lpnpq(1:ndimn, ig, 4, 1) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 1)/2.d0
!
!...point 2
lpnpq(1:ndimn, ig, 1, 2) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 2)/2.d0
lpnpq(1:ndimn, ig, 2, 2) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 2)/2.d0
lpnpq(1:ndimn, ig, 3, 2) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 2)/2.d0
lpnpq(1:ndimn, ig, 4, 2) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 2)/2.d0
!
!...point 3
lpnpq(1:ndimn, ig, 1, 3) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 3)/2.d0
lpnpq(1:ndimn, ig, 2, 3) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 3)/2.d0
lpnpq(1:ndimn, ig, 3, 3) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 3)/2.d0
lpnpq(1:ndimn, ig, 4, 3) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 3)/2.d0
!...point 4
lpnpq(1:ndimn, ig, 1, 4) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 4)/2.d0
lpnpq(1:ndimn, ig, 2, 4) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 4)/2.d0
lpnpq(1:ndimn, ig, 3, 4) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 4)/2.d0
lpnpq(1:ndimn, ig, 4, 4) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 4)/2.d0

!
enddo
!
!...Capital r
!
rcapt(1, 1) = (3.d0*rcoeq(1) + rcoeq(4))/3.d0
rcapt(2, 1) = (3.d0*rcoeq(1) + rcoeq(2))/3.d0
rcapt(3, 1) = (1.d0*rcoeq(1) + rcoeq(4))/3.d0
rcapt(4, 1) = (1.d0*rcoeq(1) + rcoeq(2))/3.d0

rcapt(1, 2) = (3.d0*rcoeq(2) + rcoeq(1))/3.d0
rcapt(2, 2) = (3.d0*rcoeq(2) + rcoeq(3))/3.d0
rcapt(3, 2) = (1.d0*rcoeq(2) + rcoeq(1))/3.d0
rcapt(4, 2) = (1.d0*rcoeq(2) + rcoeq(3))/3.d0

rcapt(1, 3) = (3.d0*rcoeq(3) + rcoeq(2))/3.d0
rcapt(2, 3) = (3.d0*rcoeq(3) + rcoeq(4))/3.d0
rcapt(3, 3) = (1.d0*rcoeq(3) + rcoeq(2))/3.d0
rcapt(4, 3) = (1.d0*rcoeq(3) + rcoeq(4))/3.d0

rcapt(1, 4) = (3.d0*rcoeq(4) + rcoeq(3))/3.d0
rcapt(2, 4) = (3.d0*rcoeq(4) + rcoeq(1))/3.d0
rcapt(3, 4) = (1.d0*rcoeq(4) + rcoeq(3))/3.d0
rcapt(4, 4) = (1.d0*rcoeq(4) + rcoeq(1))/3.d0
!
bcapt(1:3, 1, 1) = bq(1:3, 1)/2.d0
bcapt(1:3, 2, 1) = bq(1:3, 1)/2.d0
bcapt(1:3, 3, 1) = bq(1:3, 1)/2.d0
bcapt(1:3, 4, 1) = bq(1:3, 1)/2.d0

bcapt(1:3, 1, 2) = bq(1:3, 2)/2.d0
bcapt(1:3, 2, 2) = bq(1:3, 2)/2.d0
bcapt(1:3, 3, 2) = bq(1:3, 2)/2.d0
bcapt(1:3, 4, 2) = bq(1:3, 2)/2.d0

bcapt(1:3, 1, 3) = bq(1:3, 3)/2.d0
bcapt(1:3, 2, 3) = bq(1:3, 3)/2.d0
bcapt(1:3, 3, 3) = bq(1:3, 3)/2.d0
bcapt(1:3, 4, 3) = bq(1:3, 3)/2.d0
!
bcapt(1:3, 1, 4) = bq(1:3, 4)/2.d0
bcapt(1:3, 2, 4) = bq(1:3, 4)/2.d0
bcapt(1:3, 3, 4) = bq(1:3, 4)/2.d0
bcapt(1:3, 4, 4) = bq(1:3, 4)/2.d0
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvqua
do ifa =1 ,4
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, ifa, iv)*rcapt(ifa, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, ifa, iv)*rcapt(ifa, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fsqrz(1, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)*rcapt(ifa, iv)

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fsqrz(2, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)*rcapt(ifa, iv)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(iv))*fsqrz(1, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)*rcapt(ifa, iv)  +&
ustar(2, ipq(iv))*fsqrz(2, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)*rcapt(ifa, iv)!
!
!if(ie==1) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1)  print*,'rhs iface',ielem, ie,ulnpn(1)!,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lagmc_quadrz2
!
!...Face integral (mass center) for hybrid linear quad using analytical integration...
!
subroutine rhsifacedg_lagmc_quadrz3(ipqua, unkno, ustar,fsqrz, gelagq, geoel,coord,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvqua, 1:nquad),  intent(in)::fsqrz !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:2, 1:nvqua) :: ipfq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8,dimension(1:ndimn, 1:ndegr, 1:4, 1:nvqua)::lpnpq
real*8::xvq(nvqua), yvq(nvqua),bq(1:3,1:nvqua)
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:4, 1:nvqua):: rcapt
real*8,dimension(1:3, 1:4, 1:nvqua):: bcapt
real*8,dimension(1:ndimn, 1:nvqua) :: xphq
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
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
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
!
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
!...physical coordinate....
!
xphq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xphq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Coefficient R of RZ or XY system...
!
rcoeq(1:nvqua) = 1.d0 - alfrz + alfrz*abs(xphq(2, 1:nvqua))
!
!if(ie.eq.1)print*,'rcoef',rcoeq,xphq(2, :)
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpq(1:ndimn, ig, 1, 1) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 1)
lpnpq(1:ndimn, ig, 2, 1) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 1)
lpnpq(1:ndimn, ig, 3, 1) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 1)
lpnpq(1:ndimn, ig, 4, 1) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 1)
!
!...point 2
lpnpq(1:ndimn, ig, 1, 2) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 2)
lpnpq(1:ndimn, ig, 2, 2) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 2)
lpnpq(1:ndimn, ig, 3, 2) = c05*gelagq(1:ndimn, 1, ie)*gelagq(3, 1, ie)*bq(ig, 2)
lpnpq(1:ndimn, ig, 4, 2) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 2)
!
!...point 3
lpnpq(1:ndimn, ig, 1, 3) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 3)
lpnpq(1:ndimn, ig, 2, 3) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 3)
lpnpq(1:ndimn, ig, 3, 3) = c05*gelagq(1:ndimn, 2, ie)*gelagq(3, 2, ie)*bq(ig, 3)
lpnpq(1:ndimn, ig, 4, 3) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 3)
!...point 4
lpnpq(1:ndimn, ig, 1, 4) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 4)
lpnpq(1:ndimn, ig, 2, 4) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 4)
lpnpq(1:ndimn, ig, 3, 4) = c05*gelagq(1:ndimn, 3, ie)*gelagq(3, 3, ie)*bq(ig, 4)
lpnpq(1:ndimn, ig, 4, 4) = c05*gelagq(1:ndimn, 4, ie)*gelagq(3, 4, ie)*bq(ig, 4)
!
enddo
!
!...Capital r
!
rcapt(1, 1) = (3.d0*rcoeq(1) + rcoeq(4))/6.d0
rcapt(2, 1) = (3.d0*rcoeq(1) + rcoeq(2))/6.d0
rcapt(3, 1) = (1.d0*rcoeq(1) + rcoeq(4))/6.d0
rcapt(4, 1) = (1.d0*rcoeq(1) + rcoeq(2))/6.d0

rcapt(1, 2) = (3.d0*rcoeq(2) + rcoeq(1))/6.d0
rcapt(2, 2) = (3.d0*rcoeq(2) + rcoeq(3))/6.d0
rcapt(3, 2) = (1.d0*rcoeq(2) + rcoeq(1))/6.d0
rcapt(4, 2) = (1.d0*rcoeq(2) + rcoeq(3))/6.d0

rcapt(1, 3) = (3.d0*rcoeq(3) + rcoeq(2))/6.d0
rcapt(2, 3) = (3.d0*rcoeq(3) + rcoeq(4))/6.d0
rcapt(3, 3) = (1.d0*rcoeq(3) + rcoeq(2))/6.d0
rcapt(4, 3) = (1.d0*rcoeq(3) + rcoeq(4))/6.d0

rcapt(1, 4) = (3.d0*rcoeq(4) + rcoeq(3))/6.d0
rcapt(2, 4) = (3.d0*rcoeq(4) + rcoeq(1))/6.d0
rcapt(3, 4) = (1.d0*rcoeq(4) + rcoeq(3))/6.d0
rcapt(4, 4) = (1.d0*rcoeq(4) + rcoeq(1))/6.d0
!
bcapt(1:3, 1, 1) = bq(1:3, 1)
bcapt(1:3, 2, 1) = bq(1:3, 1)
bcapt(1:3, 3, 1) = bq(1:3, 1)
bcapt(1:3, 4, 1) = bq(1:3, 1)

bcapt(1:3, 1, 2) = bq(1:3, 2)
bcapt(1:3, 2, 2) = bq(1:3, 2)
bcapt(1:3, 3, 2) = bq(1:3, 2)
bcapt(1:3, 4, 2) = bq(1:3, 2)

bcapt(1:3, 1, 3) = bq(1:3, 3)
bcapt(1:3, 2, 3) = bq(1:3, 3)
bcapt(1:3, 3, 3) = bq(1:3, 3)
bcapt(1:3, 4, 3) = bq(1:3, 3)
!
bcapt(1:3, 1, 4) = bq(1:3, 4)
bcapt(1:3, 2, 4) = bq(1:3, 4)
bcapt(1:3, 3, 4) = bq(1:3, 4)
bcapt(1:3, 4, 4) = bq(1:3, 4)
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvqua
do ifa =1 ,4
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, ifa, iv)*rcapt(ifa, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, ifa, iv)*rcapt(ifa, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fsqrz(1, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)!*rcapt(ifa, iv)

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fsqrz(2, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)!*rcapt(ifa, iv)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(iv))*fsqrz(1, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)+&!*rcapt(ifa, iv)  +&
ustar(2, ipq(iv))*fsqrz(2, ifa, iv, ie)*bcapt(1:ndegr, ifa, iv)!*rcapt(ifa, iv)
!if(ie==1) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1)  print*,'rhs iface',ielem, ie,ulnpn(1)!,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lagmc_quadrz3
!
!...Get the mass matrix and parametrs for deformation gradient advancing...
!
subroutine  getamatr_lag_dg(amadg,geoel,coord,iptri, ipqua)
use constant
implicit none
!...Input
real*8,dimension(1:ngeel,1:nsize),intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),intent(in)::coord
real*8,dimension(1:nmatr,1:ncell),intent(out)::amadg
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!...Local integer
integer :: ie, ig, ishp, ielem
!...Local real array
real*8::xp(1:2, 1:nptri)
real*8,dimension(1:2, 1:npqua)::xpq
real*8,dimension(1:nptri)::shp, dspr, dsps
real*8:: weight(ngausd), posit(2, ngausd)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
!...Local real
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8::xc, yc
real*8:: dxdr,dxds,dydr,dyds
real*8::wi,djaco, volel,det
real*8::c10
real*8::f0,f1,f2,f3,f4
real*8::f5,f6,f7,f8
real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::b2,b3,b4,b5,b6
real*8::xgaus,ygaus
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
endif
!
!...Mass center...
!
rc= geoel(7, ielem)
sc= geoel(8, ielem)
!
dr = 0.5d0
ds = 0.5d0
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
do ishp = 1, nptri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo
!
xg = r
yg = s

!...Basis function
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
f0 = f0 + djaco
f1 = f1 + (xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + (xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + (yg-sc)/ds*(yg-sc)/ds*djaco

if(npoly==2)then
!
b4 = 0.5d0*b2*b2 - geoel(22, ielem)
b5 = 0.5d0*b3*b3 - geoel(23, ielem)
b6 =       b2*b3 - geoel(24, ielem)
!
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
enddo
!
!if(ie.eq.1) print*,'ielem', ie, f22,f23,f24,f25,f26
if(npoly==0)then
amadg(1,ielem) = 1.d0/f0

elseif(npoly==1)then
det = f1*f3-f2**2

amadg(1,ielem) = f3/det
amadg(2,ielem) = -f2/det
amadg(3,ielem) = f1/det
amadg(4,ielem) = 1.d0/f0
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
amadg(1,ielem) = x5(1,1)
amadg(2,ielem) = x5(1,2)
amadg(3,ielem) = x5(1,3)
amadg(4,ielem) = x5(1,4)
amadg(5,ielem) = x5(1,5)

amadg(6,ielem) = x5(2,2)
amadg(7,ielem) = x5(2,3)
amadg(8,ielem) = x5(2,4)
amadg(9,ielem) = x5(2,5)

amadg(10,ielem) = x5(3,3)
amadg(11,ielem) = x5(3,4)
amadg(12,ielem) = x5(3,5)

amadg(13,ielem) = x5(4,4)
amadg(14,ielem) = x5(4,5)

amadg(15,ielem) = x5(5,5)
!
amadg(16,ielem) = 1.d0/f0
!
endif

!...Treatment for AW RZ
!
enddo !...(2)ie = 1,nelem

!
!...Second loop for quads
!
do ie = 1,nquad !...(2)ie = 1,nelem
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
!...Mass center...
!
rc= geoel(7, ielem)
sc= geoel(8, ielem)
!
dr = 1.d0
ds = 1.d0!
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
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds

!...Coefficient R of RZ or XY system...
f0 = f0 + djaco
f1 = f1 + (xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + (xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + (yg-sc)/ds*(yg-sc)/ds*djaco

if(npoly==2)then
!
b4 = 0.5d0*b2*b2 - geoel(22, ielem)
b5 = 0.5d0*b3*b3 - geoel(23, ielem)
b6 =       b2*b3 - geoel(24, ielem)
!
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
!if(ie.eq.1) print*,'ielem', ie, f25,rho0*b2*b5*djaco,rho0,b2,b5,djaco,xg,yg
enddo
!

if(npoly==0)then
amadg(1,ielem) = 1.d0/f0

elseif(npoly==1)then
det = f1*f3-f2**2

amadg(1,ielem) = f3/det
amadg(2,ielem) = -f2/det
amadg(3,ielem) = f1/det
amadg(4,ielem) = 1.d0/f0
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
amadg(1,ielem) = x5(1,1)
amadg(2,ielem) = x5(1,2)
amadg(3,ielem) = x5(1,3)
amadg(4,ielem) = x5(1,4)
amadg(5,ielem) = x5(1,5)

amadg(6,ielem) = x5(2,2)
amadg(7,ielem) = x5(2,3)
amadg(8,ielem) = x5(2,4)
amadg(9,ielem) = x5(2,5)

amadg(10,ielem) = x5(3,3)
amadg(11,ielem) = x5(3,4)
amadg(12,ielem) = x5(3,5)

amadg(13,ielem) = x5(4,4)
amadg(14,ielem) = x5(4,5)

amadg(15,ielem) = x5(5,5)
!
amadg(16,ielem) = 1.d0/f0
!
endif

!...Treatment for AW RZ

!
enddo !...(2)ie = 1,nelem

end subroutine  getamatr_lag_dg
!
!...subroutine: Calculate the F^* N dsfor all faces for hybrid grids using the gadient of deformation...
!
subroutine getfnds_lag_hybridgd(geoel,gflag,gesgq,gesgt,intfac,iptri,ipqua,coord,unkgd2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac),intent(out)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(out)::gesgq
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd2
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,ie,ig,ideg,jdeg,ielem
integer::iv,ishp
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(1:nvfac, 1:3) :: ivft
integer,dimension(1:nvfac, 1:4) :: ivfq
integer,dimension(2, 3):: igt
integer,dimension(2, 4):: igq
real*8,dimension(1:ndegr,1:4,1:nsize)::unkgd
!...local real array
real*8,dimension(1:2, 1:2)        ::magd, cmagd !...cofactor matrix...
real*8,dimension(1:ndimn, 1:3)::xpf
real*8,dimension(1:ndimn, 1:nvtri)::xpt
real*8,dimension(1:ndimn, 1:nvqua)::xpq
real*8,dimension(1:ndegr)::bt, bq
real*8::vnorm(1:2)
real*8::posit(2, 9),posiq(2, 12),posif(3, 1)
real*8::dshpr(3)
real*8::dxdr,dxds,dydr,dyds,djaco
!...local real number
real*8::dwav1,dwav2,larea,farea
real*8::anx, any
real*8::r,s, rc, sc, rq, sq, rt, st
real*8::dr, ds
real*8::c16, c10
!
data c16   /0.1666666666666666d0 /
data c10   /1.0d0 /
!
!...Part I: Specify the posi for 9 simpson points...
!
!...Triangle
posit(1, 1) =  0.d0; posit(2, 1) = 0.d0;
posit(1, 2) =  1.d0; posit(2, 2) = 0.d0;
posit(1, 3) =  1.d0; posit(2, 3) = 0.d0;
posit(1, 4) =  0.d0; posit(2, 4) = 1.d0;
posit(1, 5) =  0.d0; posit(2, 5) = 1.d0;
posit(1, 6) =  0.d0; posit(2, 6) = 0.d0;
!
ivft(1, 1) = 1; ivft(2, 1) = 2;
ivft(1, 2) = 2; ivft(2, 2) = 3;
ivft(1, 3) = 3; ivft(2, 3) = 1;
!
igt(1, 1) = 1; igt(2, 1) = 2;
igt(1, 2) = 3; igt(2, 2) = 4;
igt(1, 3) = 5; igt(2, 3) = 6;

!...8-nodes quad
posiq(1, 1) =-1.d0; posiq(2, 1) =-1.d0;
posiq(1, 2) = 1.d0; posiq(2, 2) =-1.d0;
posiq(1, 3) = 1.d0; posiq(2, 3) =-1.d0;
posiq(1, 4) = 1.d0; posiq(2, 4) = 1.d0;
posiq(1, 5) = 1.d0; posiq(2, 5) = 1.d0;
posiq(1, 6) =-1.d0; posiq(2, 6) = 1.d0;
posiq(1, 7) =-1.d0; posiq(2, 7) = 1.d0;
posiq(1, 8) =-1.d0; posiq(2, 8) =-1.d0;
!
ivfq(1, 1) = 1; ivfq(2, 1) = 2;
ivfq(1, 2) = 2; ivfq(2, 2) = 3;
ivfq(1, 3) = 3; ivfq(2, 3) = 4;
ivfq(1, 4) = 4; ivfq(2, 4) = 1;
!
igq(1, 1) = 1; igq(2, 1) = 2;
igq(1, 2) = 3; igq(2, 2) = 4;
igq(1, 3) = 5; igq(2, 3) = 6;
igq(1, 4) = 7; igq(2, 4) = 8;

!...Face position
posif(1, 1) = -1.d0
posif(2, 1) =  1.d0
posif(3, 1) =  0.d0

!...Initialize the deformation gradient
unkgd = 0.d0
!
unkgd(1, 1, : ) = 1.d0
unkgd(1, 2, : ) = 0.d0
unkgd(1, 3, : ) = 0.d0
unkgd(1, 4, : ) = 1.d0
!
!...Part II: Calcualte the F^* NdS for every face of triangles...
!
dr = 0.5d0
ds = 0.5d0
!
do 100 ie=1, ntria!...(1)ie=1, ntria
!
ielem = ie
!
ipt(1:nvtri) = iptri(1:nvtri,ie)

!...The cell center
rc= geoel(7, ielem)
sc= geoel(8, ielem)
!
do ifa =1, 3
!
xpf(1, 1:2) = coord(1, ipt(ivft(1:2, ifa)))
xpf(2, 1:2) = coord(2, ipt(ivft(1:2, ifa)))
!
if(ncurv.eq.0)then
xpf(1:2, 3) = 0.5d0*(xpf(1:2, 1)+xpf(1:2, 2))
else
print*,'Implement in the future!'
endif
!
do ig = 1, 2!...(2)ig = 1,ngausd
!
r  = posif(ig,1)

!...Shape functions and their derivatives...

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
anx = dydr/djaco
any =-dxdr/djaco
!...The referece coordinate in the trignale
rt = posit(1, igt(ig, ifa))
st = posit(2, igt(ig, ifa))
!
bt(1) = 1.d0

if(npoly.ge.1)then
bt(2) = (rt - rc)/dr
bt(3) = (st - sc)/ds
endif
!
magd = 0.d0
do ideg = 1, ndegr
magd(1,1) = magd(1, 1) + unkgd(ideg, 1, ielem)*bt(ideg)
magd(1,2) = magd(1, 2) + unkgd(ideg, 2, ielem)*bt(ideg)
magd(2,1) = magd(2, 1) + unkgd(ideg, 3, ielem)*bt(ideg)
magd(2,2) = magd(2, 2) + unkgd(ideg, 4, ielem)*bt(ideg)
enddo
!
cmagd(1, 1) = magd(2,2)
cmagd(1, 2) =-magd(2,1)
cmagd(2, 1) =-magd(1,2)
cmagd(2, 2) = magd(1,1)

!...Identify the local No. of one internal face for left cell...
vnorm(1) = cmagd(1, 1)*anx*djaco + cmagd(1, 2)*any*djaco
vnorm(2) = cmagd(2, 1)*anx*djaco + cmagd(2, 2)*any*djaco
!
!...Unit vector...
farea    = sqrt(vnorm(1)**2 + vnorm(2)**2) !...farea: face area...
!
gesgt(1, igt(ig, ifa), ie) = vnorm(1)/farea
gesgt(2, igt(ig, ifa), ie) = vnorm(2)/farea
gesgt(3, igt(ig, ifa), ie) = farea*2.d0
enddo !...ig from 1...9

enddo !...ifa from 1...3
!
100 enddo  !...(1)ie=1. ntria
!
!...Part III: Calcualte the F^* NdS for every face of quads...
!
dr = 1.d0
ds = 1.d0
!
do 200 ie=1, nquad !...(1)ie=1,nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...coordinates
!...The cell center
rc= geoel(7, ielem)
sc= geoel(8, ielem)
!
do ifa =1, 4
!
xpf(1, 1:2) = coord(1, ipq(ivfq(1:2, ifa)))
xpf(2, 1:2) = coord(2, ipq(ivfq(1:2, ifa)))
!
if(ncurv.eq.0)then
xpf(1:2, 3) = 0.5d0*(xpf(1:2, 1)+xpf(1:2, 2))
else

endif
!
do ig = 1, 2!...(2)ig = 1,ngausd
!
r  = posif(ig, 1)

!...Shape functions and their derivatives...

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
anx = dydr/djaco
any =-dxdr/djaco
!...The referece coordinate in the trignale
rq = posiq(1, igq(ig, ifa))
sq = posiq(2, igq(ig, ifa))

bq(1) = 1.d0

if(npoly.ge.1)then
bq(2) = (rq - rc)/dr
bq(3) = (sq - sc)/ds
endif
!
magd = 0.d0
do ideg = 1, ndegr
magd(1,1) = magd(1, 1) + unkgd(ideg, 1, ielem)*bq(ideg)
magd(1,2) = magd(1, 2) + unkgd(ideg, 2, ielem)*bq(ideg)
magd(2,1) = magd(2, 1) + unkgd(ideg, 3, ielem)*bq(ideg)
magd(2,2) = magd(2, 2) + unkgd(ideg, 4, ielem)*bq(ideg)
enddo
!
cmagd(1, 1) = magd(2,2)
cmagd(1, 2) =-magd(2,1)
cmagd(2, 1) =-magd(1,2)
cmagd(2, 2) = magd(1,1)
!
!xxxx-output for debuggingxxx
!if(ie.eq.640.and.igq(ig, ifa).eq.4)print*,'vorm',magd(1,1),unkgd(1:3, 1, ielem),bq(1:3)
!if(ie.eq.640.and.igq(ig, ifa).eq.4)print*,'vorm',magd(1,2),unkgd(1:3, 2, ielem),bq(1:3)
!if(ie.eq.640.and.igq(ig, ifa).eq.4)print*,'vorm',magd(2,1),unkgd(1:3, 3, ielem),bq(1:3)
!if(ie.eq.640.and.igq(ig, ifa).eq.4)print*,'vorm',magd(2,2),unkgd(1:3, 4, ielem),magd
! if(ie.eq.280.and.igq(ig, ifa).eq.8)print*,'vorm',cmadg,anx,djaco

!...Identify the local No. of one internal face for left cell...
vnorm(1) = cmagd(1, 1)*anx*djaco + cmagd(1, 2)*any*djaco
vnorm(2) = cmagd(2, 1)*anx*djaco + cmagd(2, 2)*any*djaco
!
!...Unit vector...
farea    = sqrt(vnorm(1)**2 + vnorm(2)**2) !...farea: face area...
!
gesgq(1, igq(ig, ifa), ie) = vnorm(1)/farea
gesgq(2, igq(ig, ifa), ie) = vnorm(2)/farea
gesgq(3, igq(ig, ifa), ie) = farea*2.d0
enddo !...ig from 1...2

enddo !...ifa from 1...4
!
200 enddo  !...(1)ifa=1,nquad
!
end subroutine getfnds_lag_hybridgd
!
!...subroutine: Riemann input for hybrid linear quad grids using the deformation gradient....
!
subroutine getriem_quadl_gd(ipqua, geoel, gesgq, vlave, unkno, unkgd, strnq_devtp, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq, nquad)::strnq_devtp
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!
real*8, dimension(1:2, 1:2, 1:npoin),          intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin), intent(inout)::snsigm
!
real*8, dimension(1:2, 1:2, 1:2,  1:nvqua, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:2,  1:nvqua, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:2,  1:nvqua, 1:nquad), intent(out)::snsigmlq

!...Local integer
integer::ie, ideg, ielem, ifa, iv

!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf

!...local real array
real*8,dimension(1:ndegr, 1:nvqua)::bq,bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8::aujmp(1:3, 1:nvqua)
real*8::vnorm(1:3, 1:2, 1:nvqua)
real*8::sigma(1:2, 1:2, 1:nvqua)
real*8,dimension(1:nvqua)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xphq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi

!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2), sigma_devt(1:2,1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr, eintc
real*8::rhomv,uvtx,vvtx,evtx, pvtx, eintv, sdv
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv, rci,sci
real*8::xcrho, ycrho
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon
!
data eps   / 1.0d-14 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Loop every quad...
!
!...Initilize the stress
sigma_devt = 0.d0
!
do 350 ie = 1,nquad !...(1)ie = 1,nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria

!...shape functions
dr = 1.0d0
ds = 1.0d0
!...mass center
rc= geoel(1, ielem)
sc= geoel(2, ielem)
!...The initial cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)

!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
!
!DGP2
if(npoly.eq.2)then
bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
endif
enddo
!...Natural coordinate....
xphq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xphq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!...Coefficient R of RZ or XY system...
rcoeq(1:nvqua) = 1.d0 - alfrz + alfrz*xphq(2, 1:nvqua)

!...Give the normal vector of every face...
vnorm(1:3, 1, 1) = gesgq(1:3, 8, ie); vnorm(1:3, 2, 1) = gesgq(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 1, 2) = gesgq(1:3, 2, ie); vnorm(1:3, 2, 2) = gesgq(1:3, 3, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gesgq(1:3, 4, ie); vnorm(1:3, 2, 3) = gesgq(1:3, 5, ie) !...For point ip(3)
vnorm(1:3, 1, 4) = gesgq(1:3, 6, ie); vnorm(1:3, 2, 4) = gesgq(1:3, 7, ie) !...For point ip(3)

!...ndA=0.5d0*vnorm

vnorm(3, 1, 1) = (2.d0*rcoeq(1) + rcoeq(4))/3.d0*0.5d0*vnorm(3, 1, 1)
vnorm(3, 2, 1) = (2.d0*rcoeq(1) + rcoeq(2))/3.d0*0.5d0*vnorm(3, 2, 1)
!
vnorm(3, 1, 2) = (2.d0*rcoeq(2) + rcoeq(1))/3.d0*0.5d0*vnorm(3, 1, 2)
vnorm(3, 2, 2) = (2.d0*rcoeq(2) + rcoeq(3))/3.d0*0.5d0*vnorm(3, 2, 2)
!
vnorm(3, 1, 3) = (2.d0*rcoeq(3) + rcoeq(2))/3.d0*0.5d0*vnorm(3, 1, 3)
vnorm(3, 2, 3) = (2.d0*rcoeq(3) + rcoeq(4))/3.d0*0.5d0*vnorm(3, 2, 3)
!
vnorm(3, 1, 4) = (2.d0*rcoeq(4) + rcoeq(3))/3.d0*0.5d0*vnorm(3, 1, 4)
vnorm(3, 2, 4) = (2.d0*rcoeq(4) + rcoeq(1))/3.d0*0.5d0*vnorm(3, 2, 4)

!...Cell averaged value...
if(ndens.eq.1)then

!...Specific volume...
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then

!...Nodal
rhomc = 1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.3)then

!...Modal
rhomc = 1.d0/unkno(1, 1, ielem)
endif
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)

!...Derived values at the center
rhoct  = 1.d0/rhomc
eintc  = ectr - 0.5d0*(uctr**2 + vctr**2) !

!...Call EOS
!
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)
call getrhog_initial(rhoi,  xcrho, ycrho, xcrho, ycrho)
!call getrhoig_quad(rhoi, r, s, xpqi)
!rhoi = 1.845d0
call GetEOS(nmatel,ncase,gamlg,rhoct, eintc, rhoi, pctr, sdctr, ielem)
!
!print*,'rhoi',rhoi,ielem,xcrho,ycrho
!if(ielem.eq.1)print*,'ielem',ielem,rhoi,sdctr,pctr,eintc,rhoct,ectr,uctr,vctr
!
!pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

!...Sound speed at the center...
!sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )

!...Zero out unknv
unknvq = 0.d0

do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo

!...Density at the vertex
if(ndens.eq.1)then
rhovt  = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
!...Nodal
r = xvq(iv); s= yvq(iv)

xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

call getrhoig_quad(rhoi, rc, sc, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)

rhovt = rhon
elseif(ndens.eq.3)then
!...Modal
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bqv(1, iv) = 1.d0
bqv(2, iv) = (xvq(iv)-rcv)/dr
bqv(3, iv) = (yvq(iv)-scv)/ds

unknvq(1, iv) =0.d0

do ideg = 1,mdegr
unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
enddo
rhovt  = unknvq(1, iv)
endif

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)

!...Derived variables at the vertex
eintv = evtx - 0.5d0*(uvtx**2 + vvtx**2)

!...Call EOS
call getrhog_initial(rhoi,  xpqi(1, iv), xpqi(2, iv), xcrho, ycrho)
call GetEOS(nmatel,ncase,gamlg,rhovt, eintv, rhoi, pvtx, sdv, ielem)

!...Limiter
if(nlimi.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvq(1, iv) - rhomc)
rhovt = 1.d0/rhomv
!
uvtx = uctr + aflim(2, ielem)*(unknvq(2, iv) - uctr)
vvtx = vctr + aflim(3, ielem)*(unknvq(3, iv) - vctr)
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unknvq(2, iv) = uvtx
unknvq(3 ,iv) = vvtx
!
elseif(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvq(1, iv) - rhomc)
rhovt = 1.d0/rhomv
elseif(ndens.eq.2)then
!rhomv = 1.d0/(1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc) )
!rhovt = 1.d0/rhomv
rhovt = 1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc)
!
elseif(ndens.eq.3)then
rhovt = 1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc)
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
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)

!...updtae unknv(2:3,:)
unknvq(2, iv) = uvtx
unknvq(3 ,iv) = vvtx
!
endif

!...Deviatoric stress
if(nmatel.eq.2) call GetStress_deviat_infsmal(xvq(iv),yvq(iv),rci,sci,geoel(22:24, ielem),&
                                              sigma_devt, strnq_devtp(:,:,iv,ielem), unkgd(:,:,ielem),ielem)
!call GetStress_deviatfem_infsmal(xvq(iv),yvq(iv),xphq,xpqi,sigma_devt, strnq_devtp(:,:,iv,ielem),ielem)

!
!if(ielem.eq.1)then
!print*,'strn',iv,strnq_devtp(:,:,iv,ielem),sigma_devt
!endif

!...Get stress tensor at the vertex
sigma(1, 1, iv) = -pvtx + sigma_devt(1,1)
sigma(1, 2, iv) = 0.d0  + sigma_devt(1,2)
sigma(2, 1, iv) = 0.d0  + sigma_devt(2,1)
sigma(2, 2, iv) = -pvtx + sigma_devt(2,2)

!...Get the a_c (unit vector)
aujmp(1:2, iv) = vlave(1:2, ipq(iv)) - unknvq(2:3, iv)
acnx = aujmp(1, iv)
acny = aujmp(2, iv)
!
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, iv) = 1.e-11
else
aujmp(1:2, iv) = aujmp(1:2, iv)/sqrt(acnx**2 + acny**2)
endif
!
aujmp(3, iv) = sqrt(acnx**2 + acny**2)
enddo

!...Scale aujmp with sound speed...
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedence coefficient...
do iv   = 1, nvqua
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr + rhoct*slpdu*deltu
!
!if(ie.eq.1.or.ie.eq.30)then
! print*,'Shock impedence',ie,iv,unkno(:,3,ie),unknvq(3, iv)
!endif
!murie(iv) = rhoct*0.533d0 + rhoct*1.34d0*deltu
!murie(iv) = rhoct*(slpdu*deltu*0.5d0 + sqrt((slpdu*deltu*0.5d0)**2+sdctr**2))
enddo

!...Get the summed denominator and numerator: sum(mu*n*a_c)
do iv  = 1, nvqua
do ifa = 1, 2 !...Every corner consists of 2 faces...

!...Call Riemann solver...
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
snsigm(1:2, ipq(iv)) = snsigm(1:2, ipq(iv)) + snsigm_rie(1:2)!

!...Local variable...
munaclq(1:2, 1, ifa, iv, ie) =  munacn_rie(1:2, 1)
munaclq(1:2, 2, ifa, iv, ie) =  munacn_rie(1:2, 2)
!
munaulq(1:2, ifa, iv, ie) =  munacu_rie(1:2)
!
snsigmlq(1:2, ifa, iv, ie)=  snsigm_rie(1:2)
!
!if(ipq(iv).eq.1.or.ipq(iv).eq.31) print*,'epq', murie(iv),ie,ifa,iv,sigma(1:2, 1:2, iv),vnorm(1:3, ifa,iv),&
!munacu(1,ipq(iv)),munacu_rie(1),snsigm(1, ipq(iv)),snsigm_rie(1),sigma(1, 1, iv)
enddo
enddo
!
350 enddo  !...(1)ie = 1,nqaud!

end subroutine getriem_quadl_gd
!
!...subroutine: Get the nodal velocity U_p^* and pressure for hybrid meshes using the deformation gradient...
!
subroutine getndvelo_lag_gd(gflag,gesgt,gesgq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold, unkno, unkgd, strnq_devtp,ustar, fstar, fstarq, aflim, afvec, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac),intent(in)::gesgt
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq, nquad)::strnq_devtp
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri, 1:ntria),  intent(out)::fstar !...Riemann forces
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec

integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin),ipsin(4)

!...local real array
real*8::munacnsin(4,3),munacusin(2, 3),snsigmsin(2,3)
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8::munaci(2, 2)
real*8::radbf(nvfac),radbfx(nvfac), radbfy(nvfac),verad(nvfac)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2,dtime
real*8::acap,omega,omegac,c1cap,s1cap,uini,vini,vini1,vini2,xg,yg,rhoini
real*8,allocatable:: fpres(:,:)
real*8,allocatable:: munacn(:,:,:)
real*8,allocatable:: munacu(:,:), snsigm(:,:)
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
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munaclt(1:2, 1:2, 1:2, 1:nvtri, 1:ntria), munault(1:ndimn, 1:2, 1:nvtri,  1:ntria),&
snsigmlt(1:ndimn, 1:2,  1:nvtri,  1:ntria))
allocate (munaclq(1:2, 1:2, 1:2, 1:nvqua, 1:nquad), munaulq(1:ndimn, 1:2, 1:nvqua,  1:nquad),&
snsigmlq(1:ndimn, 1:2,  1:nvqua,  1:nquad))
allocate (fpres(1:2, 1:npoin))
!
!...Part I: Specify some preliminary work
!

!...Zero out vlave (averaged velocity)
vlave = 0.d0
indnd = 0

!...Mark the boundary nodes...
!...Shockless Noh...
if(nmatel.eq.1)then
 if(ncase.eq.2)then
  do ifa = 1, nbfac
   ipf(1:2) = intfac(3:4, ifa)
   indnd(ipf(1:2)) = 1
  enddo
 endif

!...Precribe time-varying boundary velocity
 do ifa = 1, nbfac
  ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
  if(bface(3, ifa).eq.25)then
   indnd(ipf(1:nvfac)) = 1
  endif
 enddo

endif

!...Get averaged velocity at nodes...
!...This is not used for present code...
!call getvlavenew(iptri, ipqua, geoel, vlave, unkno, aflim, afvec)

do 950 ifa = 1 , nbfac
ipf(1:2) = intfac(3:4, ifa)
if(bface(3,ifa).eq.22)then
ipf(1:2) = intfac(3:4, ifa)
if(bface(4,ifa).eq.221)then
vlave(2,ipf(1:2)) = 0.d0
elseif(bface(4,ifa).eq.222)then
vlave(1,ipf(1:2)) = 0.d0
endif
elseif(bface(3,ifa).eq.25)then
vlave(2,ipf(1:2)) = 0.d0
vlave(1,ipf(1:2)) = 0.d0
if(ncase.eq.13)then
vlave(2,ipf(1:2)) = 0.d0
vlave(1,ipf(1:2)) = 1.d0
endif
endif!
950 enddo
!
!...Part II: Loop to get the information from Riemann solver
!
do iloop= 1, 1

!...Give vlave
vlave= ustar

!...Zero out munacn
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0

!...Tria
!...In future
!...Quad
if(nquad.gt.0) call getriem_quadl_gd(ipqua, geoel, gesgq, vlave, unkno, unkgd,strnq_devtp, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)

!...Boundary condition
!!call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)!
!!call getboundary_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm)
!call getbc_lagmaire2(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm,itime)
!
!...Periodic boundary condition for 1D isentropic Sin problem...
if(nmatel.eq.1)then
if(ncase.eq.12)then
do ifa = 1, nbfac
if(bface(3, ifa).eq.31)then !...Periodic boundary...
!
ipsin(1:2) = bface(1:2, ifa); ipsin(3:4) = bface(1:2, bface(4, ifa))
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

endif


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

!if(ipoin.eq.1.or.ipoin.eq.31)print*,'Ustar',ipoin,munacn(1, 1, ipoin),munacu(1:2, ipoin),snsigm(1:2, ipoin),munaci(1, 1)

endif
enddo

!...Get the vertex velocity at the boundary
do 900 ifa = 1 , nbfac
ipf(1:2) = intfac(3:4, ifa)

if(bface(3,ifa).eq.22)then
ipf(1:2) = intfac(3:4, ifa)
if(bface(4,ifa).eq.221)then
ustar(2,ipf(1:2)) = 0.d0
elseif(bface(4,ifa).eq.222)then
ustar(1,ipf(1:2)) = 0.d0
endif
elseif(bface(3,ifa).eq.25)then
ustar(2,ipf(1:2)) = 0.d0
ustar(1,ipf(1:2)) = 0.d0

if(nmatel.eq.1)then
!...Specify boundary nodal velocity
if(ncase.eq.13)then !...Saltzman
ustar(2,ipf(1:2)) = 0.d0
ustar(1,ipf(1:2)) = 1.d0
endif
!
if(ncase.eq.5)then !...Kidder ball
radbf(1) = sqrt(coord(2,ipf(1))**2+coord(1,ipf(1))**2)
radbf(2) = sqrt(coord(2,ipf(2))**2+coord(1,ipf(2))**2)
!
radbfx(1) = coord(1,ipf(1));radbfy(1)=coord(2,ipf(1))
radbfx(2) = coord(1,ipf(2));radbfy(2)=coord(2,ipf(2))
!
dtime = (itime-1.d0)*dtfix
verad(1) =dtime/(1.d0+dtime**2)
verad(2) =dtime/(1.d0+dtime**2)
!
ustar(1, ipf(1)) = verad(1)*radbfx(1);   ustar(2, ipf(1)) = verad(1)*radbfy(1);
ustar(1, ipf(2)) = verad(2)*radbfx(2);   ustar(2, ipf(2)) = verad(2)*radbfy(2);
endif

endif!if(nmatel.eq.1)then

endif

!...Analytic boundary velocity for TGV
if(ncase.eq.-1)then
ustar(1, ipf(1)) = sin(pi*coord(1,ipf(1)))*cos(pi*coord(2,ipf(1)))
ustar(2, ipf(1)) =-cos(pi*coord(1,ipf(1)))*sin(pi*coord(2,ipf(1)))

ustar(1, ipf(2)) = sin(pi*coord(1,ipf(2)))*cos(pi*coord(2,ipf(2)))
ustar(2, ipf(2)) =-cos(pi*coord(1,ipf(2)))*sin(pi*coord(2,ipf(2)))
endif
900 enddo

enddo !iloop

!
!...Free boundary
if(nmatel.eq.2)then
 ! ustar(:, 19) = 0.d0
endif
!
!...Zero normal velocity for BC...
! call getbcvn_lag(bface, intfac, gflag, ustar)
! call getbcve_exact(bface, intfac, gflag, ustar, coord, itime)
!
!print*,'ustar',ustar(:,9)
!
!...Part III: Get the Riemann forces for face integral...
!

!...Tria

!...Quad
do ie = 1, nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria

!...Riemann forces
do iv = 1, nvqua
do ifa =1, 2
fstarq(1, ifa, iv, ie) = snsigmlq(1, ifa, iv, ie) + &
munaclq(1, 1, ifa, iv, ie)*ustar(1, ipq(iv))+&
munaclq(2, 1, ifa, iv, ie)*ustar(2, ipq(iv)) - munaulq(1, ifa, iv, ie)
fstarq(2, ifa, iv, ie) = snsigmlq(2, ifa, iv, ie) + &
munaclq(1, 2,ifa, iv, ie)*ustar(1, ipq(iv))+&
munaclq(2, 2, ifa, iv, ie)*ustar(2, ipq(iv)) - munaulq(2, ifa, iv, ie)
enddo
enddo

enddo
!
deallocate (munacn, fpres)
deallocate (munacu, snsigm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lag_gd
!
!...Face integral (mass center) for hybrid quad using the deformation gradient...
!
subroutine rhsifacedg_lag_hybridquadl_gd(ipqua, unkno, ustar,fstarq, gesgq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
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
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
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
!
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
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
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpq(1:ndimn, ig, 1, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 1))*gesgq(1:ndimn, 8, ie)*gesgq(3, 8, ie)
lpnpq(1:ndimn, ig, 2, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 1))*gesgq(1:ndimn, 1, ie)*gesgq(3, 1, ie)
!
!...point 2
lpnpq(1:ndimn, ig, 1, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 2))*gesgq(1:ndimn, 2, ie)*gesgq(3, 2, ie)
lpnpq(1:ndimn, ig, 2, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 2))*gesgq(1:ndimn, 3, ie)*gesgq(3, 3, ie)
!
!...point 3
lpnpq(1:ndimn, ig, 1, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 3))*gesgq(1:ndimn, 4, ie)*gesgq(3, 4, ie)
lpnpq(1:ndimn, ig, 2, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 3))*gesgq(1:ndimn, 5, ie)*gesgq(3, 5, ie)
!...point 4
lpnpq(1:ndimn, ig, 1, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 4))*gesgq(1:ndimn, 6, ie)*gesgq(3, 6, ie)
lpnpq(1:ndimn, ig, 2, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 4))*gesgq(1:ndimn, 7, ie)*gesgq(3, 7, ie)
!
enddo
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
do iv = 1, nvqua
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 1, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 1, iv) +&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 2, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
fstarq(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv))

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
fstarq(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv))
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(iv))*fstarq(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
ustar(2, ipq(iv))*fstarq(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
ustar(1, ipq(iv))*fstarq(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv)) +&
ustar(2, ipq(iv))*fstarq(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, iv))
!
!if(ie==213) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1)  print*,'rhs iface',ielem, ie,ulnpn(1)!,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lag_hybridquadl_gd
!
!....domain integral for hybrid linear quad cells using the deformation gradient
!
subroutine rhsdomndg_lag_quadl_gd(intfac, ipqua, coord, coold, geoel, unkno, unkgd,strnq_devtp,rhsel,aflim,afvec )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq, nquad)::strnq_devtp
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
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndegr):: b,bi, dbdr, dbds, bv
real*8:: unknod(1:nq)
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
real*8, dimension(1: ndimn, 1:ndimn)::sigmg,sigma_devt
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rci,sci
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::xcrho,ycrho
real*8::xgaus, ygaus
real*8::rhoad,uadv,vadv,eadv,rhoma
real*8::einta,sdadv
real*8::pres
real*8::djaco, wi, rcoef
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr
real*8:: rhoi, rhon
real*8:: a11, a12, a21, a22
!
!
data eps   / 1.0d-14 /
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
!  open(12,file='solution2.dat')
!    do ie=1,ncell
!       if(ie.ge.990.and.ie.le.1000)print*,'afterdomn',ie,unkno(1:ndegr,1,ie)
!    enddo
!  close(12)

!...Initilize the stress
sigma_devt = 0.d0
!
!...Loop over elements
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:4,ie)!
!
xpq(1, 1:4) = coold(1, ipq(1:nvqua))
xpq(2, 1:4) = coold(2, ipq(1:nvqua))
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!...The initial cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)

!...Physical cell center
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpq, rc, sc, xcrho, ycrho)
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
!
dsprq(1) = -0.25d0   *sm
dsprq(2) =  0.25d0   *sm
dsprq(3) =  0.25d0   *sp
dsprq(4) = -0.25d0   *sp
!
dspsq(1) =-0.25d0*rm
dspsq(2) =-0.25d0*rp
dspsq(3) = 0.25d0*rp
dspsq(4) = 0.25d0*rm
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 4
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo
!
djaco = wi!*(dxdr*dyds-dxds*dydr)

!...The derivative of basis function
if(npoly.eq.2)then
!
dbdr(4)= (r-rc)/dr**2
dbdr(5)= 0.d0
dbdr(6)= (s-sc)/dr/ds

dbds(4)= 0.d0
dbds(5)= (s-sc)/ds**2
dbds(6)= (r-rc)/dr/ds
endif

!...physical coordinate of gauss points....

!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0

!...Jacobian transformation matrix
!
jacbf = 0.d0
!
bi(1) = 1.d0
bi(2) = (r - rci)/dr
bi(3) = (s - sci)/ds
!
jacbf = 0.d0
do ideg = 1, ndegr
jacbf(1,1) = jacbf(1, 1) + unkgd(ideg, 1, ielem)*bi(ideg)
jacbf(1,2) = jacbf(1, 2) + unkgd(ideg, 2, ielem)*bi(ideg)
jacbf(2,1) = jacbf(2, 1) + unkgd(ideg, 3, ielem)*bi(ideg)
jacbf(2,2) = jacbf(2, 2) + unkgd(ideg, 4, ielem)*bi(ideg)
enddo
!
a11 = jacbf(1, 1)*dxdr + jacbf(1, 2)*dydr
a12 = jacbf(1, 1)*dxds + jacbf(1, 2)*dyds

a21 = jacbf(2, 1)*dxdr + jacbf(2, 2)*dydr
a22 = jacbf(2, 1)*dxds + jacbf(2, 2)*dyds

!
!...Cofactor matrix of Jacobian transformation matrix
!
jacbg(1, 1) = a22; jacbg(1, 2) =-a21
jacbg(2, 1) =-a12; jacbg(2, 2) = a11
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
!DGP2
if(npoly.eq.2)then
b(4) = 0.5d0*b(2)*b(2) - geoel(19, ielem)
b(5) = 0.5d0*b(3)*b(3) - geoel(20, ielem)
b(6) =       b(2)*b(3) - geoel(21, ielem)
endif

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
if(ndens.eq.1)then
!
rhoma = unknod(1)
rhoad  = 1.d0/rhoma
elseif(ndens.eq.2)then
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rhoma = 1.d0/rhon
rhoad = rhon
!
elseif(ndens.eq.3)then
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
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)

!...Derived variables
einta = eadv - 0.5d0*(uadv**2+vadv**2)


!...Call EOS
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpq, xg, yg, xgaus, ygaus)
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)

!rhoi=1.845d0
call GetEOS(nmatel,ncase,gamlg,rhoad, einta, rhoi, pres, sdadv, ielem)

!pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!if(ielem.eq.640.and.ig.eq.1) print*,'domn2',pres,rhoad,unkno(1:3,1,ielem),b(1:3),xg,yg,uadv,vadv,ig
!
if(nlimi.eq.1)then
!
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
rhoma = rhomc + aflim(1 ,ielem)*(rhoma - rhomc)
rhoad = 1.d0/rhoma
!
uadv = uctr + aflim(2, ielem)*(uadv - uctr)
vadv = vctr + aflim(3 ,ielem)*(vadv - vctr)
!
pres = pctr + aflim(4, ielem)*(pres- pctr)
!
!
elseif(nlimi.eq.6)then
!
! if(ie.ge.990.and.ie.le.1000) print*,'domn unk',ie,unkno(1, 1, ie)
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
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
!if(ielem.eq.9) print*,'domn',pres,aflim(4,ielem),xg,yg
!
endif
!
!if(ielem.eq.1)print*,'domn',ig,ielem
!...Get the deviatoric stress
if(nmatel.eq.2)call GetStress_deviat_infsmal(r,s,rci,sci,geoel(22:24, ielem),&
                    sigma_devt, strnq_devtp(:,:,ngstrnf+ig,ielem), unkgd(:,:,ielem),ielem)

!call GetStress_deviatfem_infsmal(r,s,coord(1:2,ipq(1:nvqua)),xpq,sigma_devt, strnq_devtp(:,:,nvqua+ig,ielem),ielem)
!...Stress tensor
sigmg(1, 1) = -pres + sigma_devt(1, 1);
sigmg(1, 2) =  0.d0 + sigma_devt(1, 2);
sigmg(2, 1) =  0.d0 + sigma_devt(2, 1);
sigmg(2, 2) = -pres + sigma_devt(2, 2);
!
fluxd(1,1) = gdshp(1, 1)*uadv + gdshp(2, 1)*vadv
fluxd(2,1) = gdshp(1, 2)*uadv + gdshp(2, 2)*vadv
fluxd(3,1) = gdshp(1, 3)*uadv + gdshp(2, 3)*vadv
!
fluxd(1,2) = gdshp(1, 1)*sigmg(1, 1) + gdshp(2, 1)*sigmg(1, 2)
fluxd(2,2) = gdshp(1, 2)*sigmg(1, 1) + gdshp(2, 2)*sigmg(1, 2)
fluxd(3,2) = gdshp(1, 3)*sigmg(1, 1) + gdshp(2, 3)*sigmg(1, 2)
!
fluxd(1,3) = gdshp(1, 1)*sigmg(2, 1) + gdshp(2, 1)*sigmg(2, 2)
fluxd(2,3) = gdshp(1, 2)*sigmg(2, 1) + gdshp(2, 2)*sigmg(2, 2)
fluxd(3,3) = gdshp(1, 3)*sigmg(2, 1) + gdshp(2, 3)*sigmg(2, 2)
!
fluxd(1,4) = gdshp(1, 1)*(sigmg(1,1)*uadv+sigmg(1,2)*vadv) + gdshp(2, 1)*(sigmg(2,1)*uadv+sigmg(2,2)*vadv)
fluxd(2,4) = gdshp(1, 2)*(sigmg(1,1)*uadv+sigmg(1,2)*vadv) + gdshp(2, 2)*(sigmg(2,1)*uadv+sigmg(2,2)*vadv)
fluxd(3,4) = gdshp(1, 3)*(sigmg(1,1)*uadv+sigmg(1,2)*vadv) + gdshp(2, 3)*(sigmg(2,1)*uadv+sigmg(2,2)*vadv)

!...High-order
do ideg = 1, -ndegr
fluxd(ideg,1) = gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv
fluxd(ideg,2) = gdshp(1, ideg)*(-pres)
fluxd(ideg,3) = gdshp(2, ideg)*(-pres)
fluxd(ideg,4) = (gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv)*(-pres)
enddo

!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) - fluxd(ideg,1:nq)*djaco*rcoef
enddo
!
!if(ielem.eq.640.and.ig.eq.1) print*,'domn',jacbf,bi(1:3),uadv,vadv,ig,afvec(:, :, ielem),uctr,vctr
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
!
!print*,'ieleme.q.2',rhsel(1:3,3,2)
!
end subroutine rhsdomndg_lag_quadl_gd
!
!...Face integral (mass center) for rhs of the deformation gradient...
!
subroutine rhsifacegd_lag_hybridquadl(ipqua, unkno, ustar,fstarq, gesgq0, geoel,&
rhsgd)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:4,1:ncell),        intent(out)::rhsgd
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq0
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
real*8, dimension(1:ndegr, 1:4) :: ulnpn
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
!...Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(7, ielem) !...mass center...
sc= geoel(8, ielem)
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
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
!
!...Get lpnp for every vertex...
!
do ig = 1,ndegr
!...point 1
lpnpq(1:ndimn, ig, 1, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 1))*gesgq0(1:ndimn, 8, ie)*gesgq0(3, 8, ie)
lpnpq(1:ndimn, ig, 2, 1) = c16*(2.d0*bq(ig, 1) + bq(ig, 1))*gesgq0(1:ndimn, 1, ie)*gesgq0(3, 1, ie)
!
!...point 2
lpnpq(1:ndimn, ig, 1, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 2))*gesgq0(1:ndimn, 2, ie)*gesgq0(3, 2, ie)
lpnpq(1:ndimn, ig, 2, 2) = c16*(2.d0*bq(ig, 2) + bq(ig, 2))*gesgq0(1:ndimn, 3, ie)*gesgq0(3, 3, ie)
!
!...point 3
lpnpq(1:ndimn, ig, 1, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 3))*gesgq0(1:ndimn, 4, ie)*gesgq0(3, 4, ie)
lpnpq(1:ndimn, ig, 2, 3) = c16*(2.d0*bq(ig, 3) + bq(ig, 3))*gesgq0(1:ndimn, 5, ie)*gesgq0(3, 5, ie)
!...point 4
lpnpq(1:ndimn, ig, 1, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 4))*gesgq0(1:ndimn, 6, ie)*gesgq0(3, 6, ie)
lpnpq(1:ndimn, ig, 2, 4) = c16*(2.d0*bq(ig, 4) + bq(ig, 4))*gesgq0(1:ndimn, 7, ie)*gesgq0(3, 7, ie)
!
enddo
!
!...The vertex constituting one cell...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, nvqua
!
ulnpn(1:ndegr, 1)  = ulnpn(1:ndegr, 1)+&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 1, iv) +&
ustar(1, ipq(iv))*lpnpq(1, 1:ndegr, 2, iv)
!
ulnpn(1:ndegr, 2)  = ulnpn(1:ndegr, 2)+&
ustar(1, ipq(iv))*lpnpq(2, 1:ndegr, 1, iv) +&
ustar(1, ipq(iv))*lpnpq(2, 1:ndegr, 2, iv)
!
ulnpn(1:ndegr, 3)  = ulnpn(1:ndegr, 3)+&
ustar(2, ipq(iv))*lpnpq(1, 1:ndegr, 1, iv) +&
ustar(2, ipq(iv))*lpnpq(1, 1:ndegr, 2, iv)
!
ulnpn(1:ndegr, 4)  = ulnpn(1:ndegr, 4)+&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 1, iv) +&
ustar(2, ipq(iv))*lpnpq(2, 1:ndegr, 2, iv)
!
!if(ie==1) print*,iv,ipq(iv),ielem, ustar(1:2,ipq(iv)),lpnpq(1:2,1,1,iv)
!
enddo
!
rhsgd(1:ndegr, 1, ielem) =  ulnpn(1:ndegr, 1)
rhsgd(1:ndegr, 2, ielem) =  ulnpn(1:ndegr, 2)
rhsgd(1:ndegr, 3, ielem) =  ulnpn(1:ndegr, 3)
rhsgd(1:ndegr, 4, ielem) =  ulnpn(1:ndegr, 4)
!
!if(ie==1)  print*,'rhs iface',ielem, ie,ulnpn(1)!,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacegd_lag_hybridquadl
!
!....domain integral for RHS of the deformation gradient
!
subroutine rhsdomngd_quadl(intfac, ipqua, coord, coold, geoel, unkno, unkgd,rhsgd,aflim,afvec )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8, dimension(1:nq+1, 1:nsize),      intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),     intent(in)::afvec
real*8,dimension(1:ndegr,1:4,1:ncell),intent(inout)::rhsgd
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
real*8::rci, sci
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::xgaus, ygaus
real*8::rhoad,uadv,vadv,eadv,rhoma
real*8::pres
real*8::djaco, wi, rcoef
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr
real*8:: rhoi, rhon
!
!
data eps   / 1.0d-14 /
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
!  open(12,file='solution2.dat')
!    do ie=1,ncell
!       if(ie.ge.990.and.ie.le.1000)print*,'afterdomn',ie,unkno(1:ndegr,1,ie)
!    enddo
!  close(12)
!
!...Loop over elements
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:4,ie)!
!
xpq(1, 1:4) = coold(1, ipq(1:nvqua))
xpq(2, 1:4) = coold(2, ipq(1:nvqua))
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!...The initial cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)
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
!
dsprq(1) = -0.25d0   *sm
dsprq(2) =  0.25d0   *sm
dsprq(3) =  0.25d0   *sp
dsprq(4) = -0.25d0   *sp
!
dspsq(1) =-0.25d0*rm
dspsq(2) =-0.25d0*rp
dspsq(3) = 0.25d0*rp
dspsq(4) = 0.25d0*rm
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, 4
dxdr = dxdr + dsprq(ishp)*xpq(1,ishp)
dxds = dxds + dspsq(ishp)*xpq(1,ishp)

dydr = dydr + dsprq(ishp)*xpq(2,ishp)
dyds = dyds + dspsq(ishp)*xpq(2,ishp)
enddo

!
djaco = wi

!...The derivative of basis function
if(npoly.eq.2)then
!
dbdr(4)= (r-rci)/dr**2
dbdr(5)= 0.d0
dbdr(6)= (s-sci)/dr/ds

dbds(4)= 0.d0
dbds(5)= (s-sci)/ds**2
dbds(6)= (r-rci)/dr/ds
endif

!...physical coordinate of gauss points....

!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0

!...Jacobian transformation matrix
!
jacbf = 0.d0
!
!...Jacobian transformation matrix
!
jacbf(1, 1) = dxdr; jacbf(1, 2) = dxds
jacbf(2, 1) = dydr; jacbf(2, 2) = dyds
!
!...Cofactor matrix of Jacobian transformation matrix
!
jacbg(1, 1) = jacbf(2, 2); jacbg(1, 2) =-jacbf(2, 1)
jacbg(2, 1) =-jacbf(1, 2); jacbg(2, 2) = jacbf(1, 1)
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
!DGP2
if(npoly.eq.2)then
b(4) = 0.5d0*b(2)*b(2) - geoel(19, ielem)
b(5) = 0.5d0*b(3)*b(3) - geoel(20, ielem)
b(6) =       b(2)*b(3) - geoel(21, ielem)
endif

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
if(ndens.eq.1)then
!
rhoma = unknod(1)
rhoad  = 1.d0/rhoma
elseif(ndens.eq.2)then
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rhoma = 1.d0/rhon
rhoad = rhon
!
elseif(ndens.eq.3)then
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
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!if(ielem.eq.9) print*,'domn2',pres,rhoad,unkno(1:3,1,ielem),b(1:3),xg,yg

if(nlimi.eq.6)then
!
! if(ie.ge.990.and.ie.le.1000) print*,'domn unk',ie,unkno(1, 1, ie)
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
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
!if(ielem.eq.9) print*,'domn',pres,aflim(4,ielem),xg,yg
!
endif
!
fluxd(1,1) = gdshp(1, 1)*uadv
fluxd(2,1) = gdshp(1, 2)*uadv
fluxd(3,1) = gdshp(1, 3)*uadv
!
fluxd(1,2) = gdshp(2, 1)*uadv
fluxd(2,2) = gdshp(2, 2)*uadv
fluxd(3,2) = gdshp(2, 3)*uadv
!
fluxd(1,3) = gdshp(1, 1)*vadv
fluxd(2,3) = gdshp(1, 2)*vadv
fluxd(3,3) = gdshp(1, 3)*vadv
!
fluxd(1,4) = gdshp(2, 1)*vadv
fluxd(2,4) = gdshp(2, 2)*vadv
fluxd(3,4) = gdshp(2, 3)*vadv

!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsgd(ideg,1:4,ielem)=rhsgd(ideg,1:4,ielem) - fluxd(ideg,1:4)*djaco*rcoef
enddo
!
!if(ielem.eq.1) print*,'domn',fluxd(2,2),djaco,-pres,ig
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
!
!print*,'ieleme.q.2',rhsel(1:3,3,2)
!
end subroutine rhsdomngd_quadl
