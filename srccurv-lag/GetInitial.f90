!
!...Initial density distribtuion for either quad or triangle..
!
subroutine getrhog_initial(rhoi,  xgaus, ygaus, xc, yc)
use constant
implicit none
real*8, intent(out)::rhoi
real*8, intent(in)::xgaus, ygaus, xc, yc
!
real*8::eps
real*8:: rhoini
real*8::radie, radii,radie2,radii2,radic2,sentr,rhoin,rhoex
real*8::prein,preex
real*8::rho0ba
!
integer::ishp
!...local real!
data eps   / 1.0d-06 /

!
!...Different materials
!
select case (nmatel)

!...Ideal gas
case (1)

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
!
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
radic2 = xgaus**2 + ygaus**2
!
rho0ba = (radie2-radic2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radic2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0) !...density at cell center...
!ie
rhoini = rho0ba**(1.d0/(gamlg-1.d0))
!
elseif(ncase.eq.5)then !...Kidder ball...
!
rhoini = 2.d0*exp(-(xgaus**2 + ygaus**2))
!
elseif(ncase.eq.6)then !...Sod...!
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
elseif(ncase.eq.9)then !...Triple point...!
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
elseif(ncase.eq.10)then !...Expansion...!
!
rhoini = 1.d0
!
elseif(ncase.eq.11)then !...Implosion Lazarus...!
!
rhoini  = 1.d0
!
elseif(ncase.eq.12)then !...Isentropic smooth flow...!
!
rhoini = 1.d0 + 0.9999995d0*sin(pi*xgaus)
!
elseif(ncase.eq.13)then !...Saltzman...!
rhoini  = 1.d0

elseif(ncase.eq.14)then !...Coggeshall expansion problem
!
rhoini = 1.d0

elseif(ncase.eq.15)then !...Coggeshall expansion problem
!
 if(xc.le.-4.d0)then
   rhoini = 3.85714d0
 else
   rhoini = 1.d0+0.2d0*sin(5.d0*xgaus)
 endif

else

print*,'Please specify the initital density ditribution in subroutine getrhoig_quadcurv for Subgrid method for curved Quads!'
stop
endif

!...Solid
case (2)
!
if(ncase.eq.1)then !...Bending beam
rhoini  = 2.79d0
elseif(ncase.eq.2)then !...Bending beam beryllium
rhoini  = 1.845d0
elseif(ncase.eq.3)then !...Tayor anvil
rhoini  = 2.785d0
elseif(ncase.eq.104)then !...Be Shell
rhoini  = 1.845d0
endif


end select

rhoi = rhoini
!
end subroutine getrhog_initial
!
!...Get the initial condition...
!
subroutine  getunkno_initial_lag(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
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

end subroutine  getunkno_initial_lag
!
!...Get the initial mass matrix for lagrangian...
!
subroutine  getamatr_initial_lag(unkno,amatr,geoel,coord,iptri, ipqua)
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
do ishp = 1, nptri
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

end subroutine  getamatr_initial_lag
!
!...Calculate geoel initially...
!
subroutine getgeoel_initial_lag(iptri, ipqua, geoel, coord)
use constant
implicit none
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

!...Physical coordinates
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
end subroutine getgeoel_initial_lag
!
!...Find terms in geoel for high-order...
!
!
!...Find terms in geoel for high-order for the gradient deformation...
!
subroutine getgeoelho_initial_lag(iptri, ipqua, geoel, coord)
use constant
implicit none
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
do ishp = 1, nptri
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
end subroutine getgeoelho_initial_lag
!
!...Get the initial cell center to determine the initial density
!
subroutine GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)
implicit none
integer, intent(in)::ncurv,ndimn,nvqua
real*8,dimension(1:ndimn, 1:nvqua), intent(in) :: xpqi
real*8,  intent(in) :: rc, sc
real*8, intent(out) :: xcrho, ycrho
!
real*8,dimension(nvqua)::shpq,dsprq,dspsq
!
integer:: ishp
real*8::r,s, rp, rm ,sp, sm
real*8::c10
!
c10 = 1.d0

!..Physical coordinate
r=rc;s=sc

!...Linear cell
if(ncurv.eq.0)then

rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
!
shpq(1) = 0.25d0*rm*sm
shpq(2) = 0.25d0*rp*sm
shpq(3) = 0.25d0*rp*sp
shpq(4) = 0.25d0*rm*sp

!...Quadratic cell
elseif(ncurv.eq.1)then

shpq(1) = 0.25d0*(c10-r)*r*s*(c10-s)
shpq(2) =-0.25d0*(c10+r)*(c10-s)*r*s
shpq(3) = 0.25d0*(c10+r)*(c10+s)*r*s
shpq(4) =-0.25d0*(c10-r)*(c10+s)*r*s
shpq(5) =-0.5d0*(c10-r**2)*(c10-s)*s
shpq(6) = 0.5d0*(c10+r)*(c10-s**2)*r
shpq(7) = 0.5d0*(c10-r**2)*(c10+s)*s
shpq(8) = -0.5d0*(c10-r)*(c10-s**2)*r
shpq(9) = (c10-r**2)*(c10-s**2)

elseif(ncurv.eq.2)then

call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

endif

!...Get the physical cell center
xcrho = 0.d0
ycrho = 0.d0
!
do ishp = 1, nvqua
xcrho = xcrho + shpq(ishp)*xpqi(1,ishp)
ycrho = ycrho + shpq(ishp)*xpqi(2,ishp)
enddo
!
end subroutine GetCellctr_quad_initial
!
!...Get the initial cell center to determine the initial density
!
subroutine GetCellctr_tria_initial (ncurv,ndimn,nvtri,xpti, rc, sc, xcrho, ycrho)
implicit none
integer, intent(in)::ncurv,ndimn,nvtri
real*8,dimension(1:ndimn, 1:nvtri), intent(in) :: xpti
real*8,  intent(in) :: rc, sc
real*8, intent(out) :: xcrho, ycrho
!
real*8,dimension(nvtri)::shpt
!
integer:: ishp
real*8::r,s
real*8::c10
!
c10 = 1.d0

!..Physical coordinate
r=rc;s=sc

!...Linear cell
if(ncurv.eq.0)then

shpt(1) = 1.d0-r-s
shpt(2) = r
shpt(3) = s

!...Quadratic cell
elseif(ncurv.eq.1)then

shpt(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shpt(2) = -r*(c10-2.d0*r)
shpt(3) = -s*(c10-2.d0*s)
shpt(4) = 4.d0*r*(c10-r-s)
shpt(5) = 4.d0*r*s
shpt(6) = 4.d0*s*(c10-r-s)
endif

!...Get the physical cell center
xcrho = 0.d0
ycrho = 0.d0
!
do ishp = 1, nvtri
xcrho = xcrho + shpt(ishp)*xpti(1,ishp)
ycrho = ycrho + shpt(ishp)*xpti(2,ishp)
enddo
!
end subroutine GetCellctr_tria_initial
