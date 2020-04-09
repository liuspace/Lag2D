subroutine wenop1_rieminvrnt_quad_shu5(ipqua, esqua, unkno, geoel, coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) :: coord, cooro
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
!
!...  local arrays
!
integer, parameter:: nsten=6,nfprj=6
integer:: ipq(nvqua)
integer:: mapfe(1:2,1:nfqua)
!
real*8,dimension(1:2, 1:nvqua)::xpq, xpqi, xpq_aj, xpqi_aj
real*8,dimension(1:ndegr,1:nq,1:nsize)::unknp,unknp2
real*8,dimension(1:ndegr,1:nq,1:nsten)::unkpf,unkf_cha
real*8,dimension(1:ndegr,1:nq,1:nsten)::unkrm,unkpe


real*8,dimension(1:5,1:nsize) ::geoph
real*8,dimension(1:4,1:nsize) ::rhsmt
real*8,dimension(1:4):: rhsmt_aj, lhsmt_aj,mtinv
real*8,dimension(1:nq,1:2)::rhs
real*8,dimension(1:4, 1:4)::qmat, qinvm
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:ndegr)::bq,bqp, bqp_aj
real*8::xpf(1:2, 1:2)
real*8, dimension(1:nq)::unknl,unknr
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: wi
real*8   weigh(1:nq, nsten)
real*8:: weigt(1:nq)
real*8   os(1:nq, nsten)
real*8:: weigl(nsten)
real*8:: weige(nq, nfprj)
real*8::lhsls(3, 2), rhsls(1:3, 1:nq),b1(nq),b2(nq)
real*8::a11, a12, a21, a22, matra,matrb
real*8::ai11, ai12, ai21, ai22
real*8::dnxl,dnyl,rhoxi,rhoet,etxi,etet,uvel,vvel
real*8::ulxi,ulet,vlxi,vlet
real*8::mapmt(2,2,ncell)
integer::jelaj(3)
real*8::lamda1,lamda2,lmat1,lmat2,mapt,mapd,matrc,matrd
real*8::unknc(4)
real*8::rhomc,uc,vc,etc,dxc,dyc
real*8::os1(nsten)
!
real*8:: c00, c05, c10, c20, epsil
real*8:: dr, ds, rc, sc, r, s
real*8:: xmc, ymc, xmc_aj, ymc_aj
real*8:: masel
real*8:: xcrho, ycrho
real*8:: xgaus, ygaus, xgausi, ygausi
real*8:: dx, dy, dx_aj, dy_aj
real*8:: detma
real*8:: djaco,dxdr,dxds,dydr,dyds
real*8:: djacoi,dxdri,dxdsi,dydri,dydsi
real*8:: dtx, dty, dlgt, dnx, dny
real*8:: dudx,dudy,dvdx,dvdy
real*8:: rhom1,uctr1,vctr1,ectr1,rhoc1,pctr1,sdctr1
real*8:: rhom2,uctr2,vctr2,ectr2,rhoc2,pctr2,sdctr2
real*8:: rhof, uf, vf, pf, sdf, vedn, rhosd
real*8:: rhogi

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
data rpowe / 2.d0    /   ! works for lilia case
!
!...Part 0: Basis parameters setup
!
!...Specify the weno
nweno =2

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
unknp = unkno

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
!
enddo !do ie = 1, nquad

!...I.2: L2 projection matrix

rhsmt= 0.d0
geoph(3:5, :) = 0.d0
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
dx =1.d0! maxval(xpq(1, 1:nvqua)) - minval(xpq(1, 1:nvqua))
dy =1.d0! maxval(xpq(2, 1:nvqua)) - minval(xpq(2, 1:nvqua))

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
rhsmt(1, ielem) = rhsmt(1, ielem) + rhogi*djacoi*bq(2)*bqp(2)
rhsmt(2, ielem) = rhsmt(2, ielem) + rhogi*djacoi*bq(3)*bqp(2)
rhsmt(3, ielem) = rhsmt(3, ielem) + rhogi*djacoi*bq(2)*bqp(3)
rhsmt(4, ielem) = rhsmt(4, ielem) + rhogi*djacoi*bq(3)*bqp(3)
!
geoph(3, ielem) = geoph(3, ielem) + rhogi*djacoi*bqp(2)*bqp(2)
geoph(4, ielem) = geoph(4, ielem) + rhogi*djacoi*bqp(2)*bqp(3)
geoph(5, ielem) = geoph(5, ielem) + rhogi*djacoi*bqp(3)*bqp(3)
!
enddo

!...Physical unkno
detma = geoph(3, ielem)*geoph(5, ielem)-geoph(4, ielem)*geoph(4, ielem)
mtinv(1) = geoph(5, ielem)
mtinv(2) =-geoph(4, ielem)
mtinv(3) = mtinv(2)
mtinv(4) = geoph(3, ielem)
!
mtinv = mtinv/detma
!
rhs(1:nq, 1) = unkno(2, 1:nq, ielem)*rhsmt(1, ielem) + unkno(3, 1:nq, ielem)*rhsmt(2, ielem)
rhs(1:nq, 2) = unkno(2, 1:nq, ielem)*rhsmt(3, ielem) + unkno(3, 1:nq, ielem)*rhsmt(4, ielem)

!...Unkowns on physical domain
unknp(1, 1:nq, ielem) = unkno(1, 1:nq, ielem)

unknp(2, 1:nq, ielem) = mtinv(1)*rhs(1:nq, 1) + mtinv(2)*rhs(1:nq, 2)
unknp(3, 1:nq, ielem) = mtinv(3)*rhs(1:nq, 1) + mtinv(4)*rhs(1:nq, 2)


!...LANL
dudx = unknp(2, 2, ielem)
dudy = unknp(3, 2, ielem)
dvdx = unknp(2, 3, ielem)
dvdy = unknp(3, 3, ielem)
!
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy
!
matra = unkno(1, 2, ielem)
matrb = unkno(1, 3, ielem)
if(sqrt(matra**2+matrb**2).gt.1.d-12)then
mapmt(1, 1, ielem) = matra/sqrt((matra)**2 + matrb**2)
mapmt(1, 2, ielem) = matrb/sqrt((matra)**2 + matrb**2)
mapmt(2, 1, ielem) =-matrb/sqrt((matra)**2 + matrb**2)
mapmt(2, 2, ielem) = matra/sqrt((matra)**2 + matrb**2)

else
mapmt(1, 1, ielem) = 0.d0
mapmt(1, 2, ielem) = 0.d0
mapmt(2, 1, ielem) = 0.d0
mapmt(2, 2, ielem) = 0.d0
endif

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
!mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
!mapmt(1, 2, ielem) =          matrc/lmat1
!mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
!mapmt(2, 2, ielem) =          matrc/lmat2
!
elseif(matrc.lt.0.d0)then
!
lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
!mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
!mapmt(2, 2, ielem) =          matrc/lmat1
!mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
!mapmt(1, 2, ielem) =          matrc/lmat2
endif

else
!mapmt(1, 1, ielem) = 0.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 0.d0
endif

mapmt(1, 1, ielem) = 1.d0
mapmt(1, 2, ielem) = 0.d0
mapmt(2, 1, ielem) = 0.d0
mapmt(2, 2, ielem) = 1.d0

enddo !do ie = 1, nquad
!
!
!...Part II: Recontruction
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...physical mass center
xmc =geoph(1, ielem)
ymc =geoph(2, ielem)

!...scale
dx =1.d0!maxval(xpq(1, 1:nvqua)) - minval(xpq(1, 1:nvqua))
dy =1.d0!maxval(xpq(2, 1:nvqua)) - minval(xpq(2, 1:nvqua))

!...  b. curvatures for the face-neighboring cells
unkpe = 0.d0

select Case (nweno)
Case (1)

!...Hermite weno
isten = 1
unkpe(:, 1:nq, isten) = unknp(:, 1:nq, ielem)
!
do ies = 1, 4
jelem = esqua(ies,ie)

!print*,'ielen',ielem,jelem,ipqua(1:nvqua,jelem),ncell

if(jelem .le. ncell) then
isten = isten + 1
!
!xpq_aj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
!xpq_aj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!
dx_aj = 1.d0!maxval(xpq_aj(1, 1:nvqua)) - minval(xpq_aj(1, 1:nvqua))
dy_aj = 1.d0!maxval(xpq_aj(2, 1:nvqua)) - minval(xpq_aj(2, 1:nvqua))
!
unkpe(1, 1:nq, isten) = unknp(1, 1:nq, ielem)
unkpe(2, 1:nq, isten) = unknp(2, 1:nq, jelem)*dx/dx_aj
unkpe(3, 1:nq, isten) = unknp(3, 1:nq, jelem)*dy/dy_aj
endif !if(jelem .le. ncell) then
enddo !do ies = 1, 4

case(2)

!...Stencil 2
isten = 1
unkpe(:, 1:nq, isten) = unknp(:, 1:nq, ielem)

call wenop1_unkstcl_quad2(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoph,&
coord, cooro, esuv1, esuv2)

case(3)

!...Stencil 3
isten = 1
unkpe(:, 1:nq, isten) = unknp(:, 1:nq, ielem)

call wenop1_unkstcl_quad(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoph, coord, cooro)

end select

!...Local system
a11 = mapmt(2, 2, ielem)
a12 =-mapmt(1, 2, ielem)
a21 =-mapmt(2, 1, ielem)
a22 = mapmt(1, 1, ielem)
!
ai11 = mapmt(1, 1, ielem)
ai12 = mapmt(1, 2, ielem)
ai21 = mapmt(2, 1, ielem)
ai22 = mapmt(2, 2, ielem)
!...Local system
do is = 1, isten
!
rhoxi = unkpe(2, 1, is)*a11 + unkpe(3, 1, is)*a21
rhoet = unkpe(2, 1, is)*a12 + unkpe(3, 1, is)*a22
!
ulxi = (ai11*unkpe(2, 2, is) + ai12*unkpe(2, 3, is))*a11 + &
(ai11*unkpe(3, 2, is) + ai12*unkpe(3, 3, is))*a21
ulet = (ai11*unkpe(2, 2, is) + ai12*unkpe(2, 3, is))*a12 + &
(ai11*unkpe(3, 2, is) + ai12*unkpe(3, 3, is))*a22
!
vlxi = (ai21*unkpe(2, 2, is) + ai22*unkpe(2, 3, is))*a11 + &
(ai21*unkpe(3, 2, is) + ai22*unkpe(3, 3, is))*a21
vlet = (ai21*unkpe(2, 2, is) + ai22*unkpe(2, 3, is))*a12 + &
(ai21*unkpe(3, 2, is) + ai22*unkpe(3, 3, is))*a22
!
etxi = unkpe(2, 4, is)*a11 + unkpe(3, 4, is)*a21
etet = unkpe(2, 4, is)*a12 + unkpe(3, 4, is)*a22
!
unkpe(2, 1, is) =rhoxi
unkpe(3, 1, is) =rhoet

unkpe(2, 2, is) =ulxi
unkpe(3, 2, is) =ulet

unkpe(2, 3, is) =vlxi
unkpe(3, 3, is) =vlet

unkpe(2, 4, is) =etxi
unkpe(3, 4, is) =etet
enddo
!
!if(ielem.eq.1542.or.ielem.eq.1543)then
!print*,'ielem',ielem,unkpe(2:3, 2:3, 1:5)
!endif
!...  b. curvatures for the face-neighboring cells
ifprj = 0
unkpf = 0.d0
unkrm = 0.d0
!
do ies = 1, 4

jelem = esqua(ies,ie)
if(jelem .le. ncell) then

!...Find the projection vector
ifprj = ifprj +1

!...normal vector
xpf(1, 1:2) = xpq(1, mapfe(1:2, ies))
xpf(2, 1:2) = xpq(2, mapfe(1:2, ies))

dtx = xpf(1, 2) - xpf(1, 1)
dty = xpf(2, 2) - xpf(2, 1)

dlgt = sqrt(dtx**2 + dty**2)

dnx = dty/dlgt
dny =-dtx/dlgt
!
dnxl = dnx*mapmt(1,1,ielem) + dny*mapmt(1,2,ielem)
dnyl = dnx*mapmt(2,1,ielem) + dny*mapmt(2,2,ielem)
!
dnx = dnxl
dny = dnyl

!
!if(ielem.eq.1)then
! print*,'bad',ielem,dnx,dny
!endif
!...left and right unknowns
unknl(1:nq) = unknp(1, 1:nq, ielem)
unknr(1:nq) = unknp(1, 1:nq, jelem)
!
uvel = unknl(2)*mapmt(1,1,ielem) + unknl(3)*mapmt(1,2,ielem)
vvel = unknl(2)*mapmt(2,1,ielem) + unknl(3)*mapmt(2,2,ielem)

unknl(2) = uvel
unknl(3) = vvel

!
uvel = unknr(2)*mapmt(1,1,ielem) + unknr(3)*mapmt(1,2,ielem)
vvel = unknr(2)*mapmt(2,1,ielem) + unknr(3)*mapmt(2,2,ielem)

unknr(2) = uvel
unknr(3) = vvel

!...Matrix
qmat = 0.d0
qinvm = 0.d0
!
qmat(1,1) = 1.d0
qmat(2,2) = 1.d0
qmat(3,3) = 1.d0
qmat(4,4) = 1.d0
!
qinvm(1,1) = 1.d0
qinvm(2,2) = 1.d0
qinvm(3,3) = 1.d0
qinvm(4,4) = 1.d0

!...Get Right eigenvectors
!
!print*,'ielem',ielem,qinvm(:,2)

!call getmatrix_prj2(qmat, qinvm, dnx, dny, unknl,unknr,ielem)

!...Riemann invariant
unkf_cha = 0.d0
do is= 1, isten
do iq= 1, nq
unkf_cha(:, 1, is) =  unkf_cha(:, 1, is) + qinvm(1, iq)*unkpe(:,iq, is)
unkf_cha(:, 2, is) =  unkf_cha(:, 2, is) + qinvm(2, iq)*unkpe(:,iq, is)
unkf_cha(:, 3, is) =  unkf_cha(:, 3, is) + qinvm(3, iq)*unkpe(:,iq, is)
unkf_cha(:, 4, is) =  unkf_cha(:, 4, is) + qinvm(4, iq)*unkpe(:,iq, is)
enddo
enddo

!...Smooth indicator
os = 0.d0

!...Average
rhomc = unknp(1, 1, ielem)
etc = unknp(1, 4, ielem)
uc = unknl(2)*mapmt(1,1,ielem) + unknl(3)*mapmt(1,2,ielem)
vc = unknl(2)*mapmt(2,1,ielem) + unknl(3)*mapmt(2,2,ielem)
!
unknc(1) = rhomc
unknc(2) = uc
unknc(3) = vc
unknc(4) = etc

!...scale
dxc =maxval(xpq(1, 1:nvqua)) - minval(xpq(1, 1:nvqua))
dyc =maxval(xpq(2, 1:nvqua)) - minval(xpq(2, 1:nvqua))

!
dnxl = dxc*mapmt(1,1,ielem) + dyc*mapmt(1,2,ielem)
dnyl = dxc*mapmt(2,1,ielem) + dyc*mapmt(2,2,ielem)
!
dxc = 1.d0!dnxl
dyc = 1.d0!dnyl

!
do is= 1, isten
!os(:, is) = ((unkf_cha(2, :, is)*dxc)**2 + (unkf_cha(3, :, is)*dyc)**2)
os(:, is) = ((unkf_cha(2, :, is)*dxc)**2 +&
            (unkf_cha(3, :, is)*dyc)**2)!/(abs(unknc(:))+epsil)
enddo

!...Mometum
!os1(:) = 0.25d0*(os(1, :)+os(2, :)+os(3, :)+os(4, :))
os1(:) = max(os(1, :),os(2, :),os(3, :),os(4, :))


do is= 1, isten
os(1:nq, is) = os1(is)
enddo

!...Linear weight
weigl(1)=0.5d0;
weigl(2:isten)= 0.5d0/(isten-1.d0)

!weigl = 1.d0

weigt = 0.d0
do is= 1, isten
do iq =1,nq
weigh(iq, is) = weigl(is)/(epsil+os(iq,is))**rpowe
enddo
enddo
!
do iq =1,nq
do is= 1, isten
weigt(iq) = weigt(iq) + weigh(iq,is)
enddo
enddo
!
do is= 1, isten
weigh(1:nq,is) = weigh(1:nq,is)/weigt(1:nq)
enddo

!...Reconstructed face unknown
do is= 1, isten
do iq = 1, nq
unkrm(:,iq,ifprj) = unkrm(:,iq, ifprj) + weigh(iq,is)*unkf_cha(:, iq, is)
enddo
enddo

!
do iq = 1, nq
unkpf(:, 1, ifprj) =  unkpf(:, 1, ifprj) + qmat(1, iq)*unkrm(:,iq, ifprj)
unkpf(:, 2, ifprj) =  unkpf(:, 2, ifprj) + qmat(2, iq)*unkrm(:,iq, ifprj)
unkpf(:, 3, ifprj) =  unkpf(:, 3, ifprj) + qmat(3, iq)*unkrm(:,iq, ifprj)
unkpf(:, 4, ifprj) =  unkpf(:, 4, ifprj) + qmat(4, iq)*unkrm(:,iq, ifprj)
enddo
!
endif !if(jelem .le. ncell) then
enddo

!... compute the oscillation indicators
weige=1.d0/ifprj

!... finally, we can reconstruct the physical polynomails for one cell
unknp2(2:3, 1:nq ,ielem) = c00
do ifa = 1, ifprj
unknp2(2, 1:nq, ielem) =unknp2(2, 1:nq, ielem) + weige(1:nq, ifa)*unkpf(2,1:nq,ifa)
unknp2(3, 1:nq, ielem) =unknp2(3, 1:nq, ielem) + weige(1:nq, ifa)*unkpf(3,1:nq,ifa)
enddo
!
a11 = mapmt(1, 1, ielem)
a12 = mapmt(1, 2, ielem)
a21 = mapmt(2, 1, ielem)
a22 = mapmt(2, 2, ielem)
!
ai11 = a22
ai12 =-a12
ai21 =-a21
ai22 = a11
!
rhoxi = unknp2(2, 1, ielem)*a11 + unknp2(3, 1, ielem)*a21
rhoet = unknp2(2, 1, ielem)*a12 + unknp2(3, 1, ielem)*a22
!
ulxi = (ai11*unknp2(2, 2, ielem) + ai12*unknp2(2, 3, ielem))*a11 + &
(ai11*unknp2(3, 2, ielem) + ai12*unknp2(3, 3, ielem))*a21
ulet = (ai11*unknp2(2, 2, ielem) + ai12*unknp2(2, 3, ielem))*a12 + &
(ai11*unknp2(3, 2, ielem) + ai12*unknp2(3, 3, ielem))*a22
!
vlxi = (ai21*unknp2(2, 2, ielem) + ai22*unknp2(2, 3, ielem))*a11 + &
(ai21*unknp2(3, 2, ielem) + ai22*unknp2(3, 3, ielem))*a21
vlet = (ai21*unknp2(2, 2, ielem) + ai22*unknp2(2, 3, ielem))*a12 + &
(ai21*unknp2(3, 2, ielem) + ai22*unknp2(3, 3, ielem))*a22
!
etxi = unknp2(2, 4, ielem)*a11 + unknp2(3, 4, ielem)*a21
etet = unknp2(2, 4, ielem)*a12 + unknp2(3, 4, ielem)*a22
!
unknp2(2, 1, ielem) =rhoxi
unknp2(3, 1, ielem) =rhoet

unknp2(2, 2, ielem) = ulxi
unknp2(3, 2, ielem) =  ulet

unknp2(2, 3, ielem) =  vlxi
unknp2(3, 3, ielem) = vlet
!
!unknp2(2:3, 2:3, ielem) = 0.d0

unknp2(2, 4, ielem)  =etxi
unknp2(3, 4, ielem)  =etet


!...  end of the loop over the quad cells
enddo
!

!print*,'unkno', unkno(2:3, 4, 5051),unknp(2:3, 4, 5051)

!
!...Get the final reference polynomial
!
do ie = 1, nquad
!
ielem = ie + ntria

detma = rhsmt(1, ielem)*rhsmt(4, ielem) - rhsmt(2, ielem)*rhsmt(3, ielem)
mtinv(1) = rhsmt(4, ielem)
mtinv(2) =-rhsmt(2, ielem)
mtinv(3) =-rhsmt(3, ielem)
mtinv(4) = rhsmt(1, ielem)
!
mtinv = mtinv/detma
!
rhs(1:nq, 1) = unknp2(2, 1:nq, ielem)*geoph(3, ielem) + unknp2(3, 1:nq, ielem)*geoph(4, ielem)
rhs(1:nq, 2) = unknp2(2, 1:nq, ielem)*geoph(4, ielem) + unknp2(3, 1:nq, ielem)*geoph(5, ielem)
!
unkno(2, 1:nq, ielem) = mtinv(1)*rhs(1:nq, 1) + mtinv(2)*rhs(1:nq, 2)
unkno(3, 1:nq, ielem) = mtinv(3)*rhs(1:nq, 1) + mtinv(4)*rhs(1:nq, 2)
!
enddo
!
!unkno(2:3,2:3,:) = 0.d0
!
!
return
end subroutine wenop1_rieminvrnt_quad_shu5
