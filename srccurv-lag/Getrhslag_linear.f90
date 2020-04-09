!
!...subroutine: Riemann free solver to calculate the u* and p* for linear quad....
!
subroutine getriem_quad_rsf(ipqua, geoel, gelagq, vlave, unkno, ustar,rhust,&
sgmst,lhsrm,lhsrmu,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!
real*8, dimension(1:2, 1:npoin),          intent(in)::ustar
real*8, dimension(1:2, 1:npoin),          intent(out)::rhust
real*8, dimension(1:2, 1:2, 1:npoin),     intent(out)::sgmst
real*8, dimension(1:npoin),               intent(out)::lhsrm,lhsrmu

!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig,ishp
!...local integer array
integer,dimension(1:nvqua) :: ipq

!...local real array
real*8,dimension(1:ndegr):: b,bv
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xphq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:nq)::unknod
real*8,dimension(1:ndimn, 1:ndimn)::sigmg,sgmas
real*8,dimension(2,ngausdq)::posiq
real*8,dimension(ngausdq)  ::weighq

!
real*8::eps,c00,c05,c10,c20
real*8::wi
real*8::rcoef,xg,yg,xgaus,ygaus
real*8::dxdr,dxds,dydr,dyds,djaco
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhoma,rhoad,uadv,vadv,eadv,pres,sdg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::dudr, duds, dvdr, dvds
real*8::rhoi, rhon
real*8::darea
!
data eps   / 1.0d-6 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /

!
!...Part 0: Preliminary setup...
!

!...Initialize sgmas
sgmas = 0.d0
lhsrm = 0.d0
lhsrmu = 0.d0
sgmst = 0.d0
rhust = 0.d0

darea = 0.d0

!...Quadrature points
call ruqope(2, ngausdq, posiq, weighq)

!...Geometry parameters for reference cell...
dr = 1.d0
ds = 1.d0
!
!...Part I: Loop every quad...
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria

!...Points consitituting one element...
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
!...Gauss loop
!
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
djaco = wi*(dxdr*dyds - dydr*dxds)

!...physical coordinate of gauss points....
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo

!...Coefficient R of RZ or XY system...
rcoef = 1.d0 - alfrz + alfrz*ygaus

!...Gauss points...
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

!...Solution at the Gauss points...
unknod = 0.d0
do ideg =1,mdegr
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo

!...Primitive variables...
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
if(ielem.eq.-1)then
print*,'ielem',ig,eadv,b,unkno(:,4,ielem)
endif
!
if(nlimi.eq.6)then

! if(ie.ge.990.and.ie.le.1000) print*,'domn unk',ie,unkno(1, 1, ie)

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
sdg = sqrt(gamlg*pres/rhoad)

!...Stress tensor
sigmg(1, 1) = -pres; sigmg(1, 2) =  0.d0;
sigmg(2, 1) = 0.d0;  sigmg(2, 2) = -pres;
!
!if(ielem.eq.991)print*,'volume integral'
call getAS(xg,yg, unkno(:,:,ielem), geoel(:,ielem) ,xpq, xpqi, sgmas, ielem)
!
!sgmas(1,1) = rhoct*sqrt(gamlg*pctr/rhoct)*(uadv-)
!
sigmg(1, 1) = -pres !+ sgmas(1, 1)
sigmg(1, 2) = 0.d0  !+ sgmas(1, 2)
sigmg(2, 1) = 0.d0  !+ sgmas(2, 1)
sigmg(2, 2) = -pres !+ sgmas(2, 2)!

!...lhs for stress and velovity
lhsrm(ipq(1:nvqua)) =  lhsrm(ipq(1:nvqua)) + shpq(1:nvqua)*djaco
lhsrmu(ipq(1:nvqua)) =  lhsrmu(ipq(1:nvqua)) + rhoad*sdg*shpq(1:nvqua)*djaco

!...rhs for nodal stress
sgmst(1,1,ipq(1:nvqua)) = sgmst(1,1,ipq(1:nvqua)) + sigmg(1,1)*shpq(1:nvqua)*djaco
sgmst(1,2,ipq(1:nvqua)) = sgmst(1,2,ipq(1:nvqua)) + sigmg(1,2)*shpq(1:nvqua)*djaco
sgmst(2,1,ipq(1:nvqua)) = sgmst(2,1,ipq(1:nvqua)) + sigmg(2,1)*shpq(1:nvqua)*djaco
sgmst(2,2,ipq(1:nvqua)) = sgmst(2,2,ipq(1:nvqua)) + sigmg(2,2)*shpq(1:nvqua)*djaco

!...rhs for nodal velocity
rhust(1,ipq(1:nvqua)) = rhust(1,ipq(1:nvqua)) + rhoad*sdg*uadv*shpq(1:nvqua)*djaco
rhust(2,ipq(1:nvqua)) = rhust(2,ipq(1:nvqua)) + rhoad*sdg*vadv*shpq(1:nvqua)*djaco
!
!if(ielem.eq.1)darea=darea+djaco
!
enddo !...(2)ig = 1,ngausd
!
!
650 enddo
!

end subroutine getriem_quad_rsf
!
!...subroutine: Get the nodal velocity U_p^* and pressure for hybrid meshes with RSF...
!
subroutine getRiemvtx_lag_rfm(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold, unkno,ustar, fstar, fstarq, aflim, afvec, itime)
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
integer,dimension(npoin)   ::idxnd

!...local real array
real*8,dimension(1:3,1:2,1:nvqua)::anmq
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8::radbf(nvfac),radbfx(nvfac), radbfy(nvfac),verad(nvfac)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::anx,any
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2,dtime
!
real*8,allocatable::sgmst(:,:,:),lhsrm(:),lhsrmu(:),rhust(:,:)

!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /

!
!...Part I: Specify some preliminary work
!
allocate(sgmst(2,2,npoin),lhsrm(npoin),lhsrmu(npoin),rhust(2,npoin))

!...Zero out vlave (averaged velocity)
vlave = 0.d0
idxnd = 0

!...Mark the boundary nodes...
!...Shockless Noh...
if(ncase.eq.2)then
 do ifa = 1, nbfac
  ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
  idxnd(ipf(1:nvfac)) = 1
 enddo
endif
!
!...Part II: Loop to get the information from Riemann solver
!
do iloop= 1, 1

!...Give vlave
vlave= ustar

!...Tria

!...Quad
if(nquad.gt.0) call getriem_quad_rsf(ipqua, geoel, gelagq, vlave, unkno, ustar, rhust,&
sgmst,lhsrm,lhsrmu,coord, coold, aflim, afvec)

do ipoin = 1, npoin
  sgmst(1:2,1:2, ipoin) = sgmst(1:2,1:2, ipoin)/lhsrm(ipoin)
enddo
call getriem_quad_velo_rsf2(ipqua, geoel, gelagq, vlave, unkno, ustar,rhust,&
sgmst,lhsrmu,coord, coold, aflim, afvec)

!...Boundary condition
!print*,'ipoin8',sgmst(1,1, 2),sgmst(1,1, 8),lhsrm(2),lhsrm(8)

!...Periodic boundary condition for 1D isentropic Sin problem...

!...Update the velocity and stress at the vertex...
do ipoin = 1, npoin
 if(idxnd(ipoin).eq.0)then
   ustar(1:2, ipoin)     = rhust(1:2, ipoin)/lhsrmu(ipoin)
 endif
!sgmst(1:2,1:2, ipoin) = sgmst(1:2,1:2, ipoin)/lhsrm(ipoin)
ustar(1,ipoin) = sin(pi*coord(1,ipoin))*cos(pi*coord(2,ipoin))
ustar(2,ipoin) =-cos(pi*coord(1,ipoin))*sin(pi*coord(2,ipoin))
!
sgmst(1,1, ipoin) =-0.25d0*(cos(2.d0*pi*coord(1,ipoin))+cos(2.d0*pi*coord(2,ipoin)))-1.d0
sgmst(2,2, ipoin) =-0.25d0*(cos(2.d0*pi*coord(1,ipoin))+cos(2.d0*pi*coord(2,ipoin)))-1.d0

enddo

!print*,'ipoin8',sgmst(1,1, 1:2),sgmst(1,1, 7:8)

!
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
!...Part III: Get the Riemann forces for face integral...
!
!print*,'normal vector',gelagq(:,:,1)
!...Tria
!...'Traingle will be implemented in future!'

!...Quad
do ie = 1, nquad

ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!...Give the normal vector of every face...
anmq(1:3, 1, 1) = gelagq(1:3, 4, ie); anmq(1:3, 2, 1) = gelagq(1:3, 1, ie)
anmq(1:3, 1, 2) = gelagq(1:3, 1, ie); anmq(1:3, 2, 2) = gelagq(1:3, 2, ie)
anmq(1:3, 1, 3) = gelagq(1:3, 2, ie); anmq(1:3, 2, 3) = gelagq(1:3, 3, ie)
anmq(1:3, 1, 4) = gelagq(1:3, 3, ie); anmq(1:3, 2, 4) = gelagq(1:3, 4, ie)
!

!...Riemann forces
do iv = 1, nvqua
do ifa =1, 2

anx = anmq(1, ifa, iv)*anmq(3, ifa, iv)*0.5d0
any = anmq(2, ifa, iv)*anmq(3, ifa, iv)*0.5d0

fstarq(1, ifa, iv, ie) = sgmst(1, 1, ipq(iv))*anx + sgmst(1, 2, ipq(iv))*any
fstarq(2, ifa, iv, ie) = sgmst(2, 1, ipq(iv))*anx + sgmst(2, 2, ipq(iv))*any
enddo
enddo
!

enddo

!...Free the allocatable arrays
deallocate(sgmst, lhsrm)

end subroutine getRiemvtx_lag_rfm
!
!...subroutine: Riemann free solver to calculate the u* and p* for linear quad....
!
subroutine getriem_quad_rsf2(ipqua, geoel, gelagq, vlave, unkno, ustar, rhust,&
sgmst,lhsrm,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!
real*8, dimension(1:2, 1:npoin),          intent(in)::ustar
real*8, dimension(1:2, 1:npoin),          intent(out)::rhust
real*8, dimension(1:2, 1:2, 1:npoin),     intent(out)::sgmst
real*8, dimension(1:npoin),               intent(out)::lhsrm

!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig,ishp
!...local integer array
integer,dimension(1:nvqua) :: ipq

!...local real array
real*8,dimension(1:ndegr):: b,bv
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xphq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:nq)::unknod
real*8,dimension(1:ndimn, 1:ndimn)::sigmg,sgmas
real*8,dimension(2,ngausdq)::posiq
real*8,dimension(ngausdq)  ::weighq

!
real*8::eps,c00,c05,c10,c20
real*8::wi
real*8::rcoef,xg,yg,xgaus,ygaus
real*8::dxdr,dxds,dydr,dyds,djaco
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhoma,rhoad,uadv,vadv,eadv,pres
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::dudr, duds, dvdr, dvds
real*8::rhoi, rhon
real*8::darea
!
data eps   / 1.0d-6 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /

!
!...Part 0: Preliminary setup...
!

!...Initialize sgmas
sgmas = 0.d0
lhsrm = 0.d0
sgmst = 0.d0
rhust = 0.d0

darea = 0.d0

!...Quadrature points
call ruqope(2, ngausdq, posiq, weighq)

!...Geometry parameters for reference cell...
dr = 1.d0
ds = 1.d0
!
!...Part I: Loop every quad...
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria

!...Points consitituting one element...
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
!...Gauss loop
!
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
djaco = wi*(dxdr*dyds - dydr*dxds)

!...Gauss points...
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

!...Solution at the Gauss points...
unknod = 0.d0
do ideg =1,mdegr
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo

!...Primitive variables...
if(ndens.eq.1)then
!
rhoma = unknod(1)
rhoad  = 1.d0/rhoma
!
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))

!...Stress tensor
sigmg(1, 1) = -pres; sigmg(1, 2) =  0.d0;
sigmg(2, 1) = 0.d0;  sigmg(2, 2) = -pres;
!
!if(ielem.eq.991)print*,'volume integral'
!
!
sigmg(1, 1) = -pres !+ sgmas(1, 1)
sigmg(1, 2) = 0.d0  !+ sgmas(1, 2)
sigmg(2, 1) = 0.d0  !+ sgmas(2, 1)
sigmg(2, 2) = -pres !+ sgmas(2, 2)!

!...lhs for stress and velovity
lhsrm(ipq(1:nvqua)) =  lhsrm(ipq(1:nvqua)) + shpq(1:nvqua)*djaco

!...rhs for nodal stress
sgmst(1,1,ipq(1:nvqua)) = sgmst(1,1,ipq(1:nvqua)) + sigmg(1,1)*shpq(1:nvqua)*djaco
sgmst(1,2,ipq(1:nvqua)) = sgmst(1,2,ipq(1:nvqua)) + sigmg(1,2)*shpq(1:nvqua)*djaco
sgmst(2,1,ipq(1:nvqua)) = sgmst(2,1,ipq(1:nvqua)) + sigmg(2,1)*shpq(1:nvqua)*djaco
sgmst(2,2,ipq(1:nvqua)) = sgmst(2,2,ipq(1:nvqua)) + sigmg(2,2)*shpq(1:nvqua)*djaco

!...rhs for nodal velocity
rhust(1,ipq(1:nvqua)) = rhust(1,ipq(1:nvqua)) + uadv*shpq(1:nvqua)*djaco
rhust(2,ipq(1:nvqua)) = rhust(2,ipq(1:nvqua)) + vadv*shpq(1:nvqua)*djaco
!!
!if(ielem.eq.1)darea=darea+djaco
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
!
!print*,'area',darea

end subroutine getriem_quad_rsf2
!
!...subroutine: Riemann free solver to calculate the u*  for linear quad....
!
subroutine getriem_quad_velo_rsf(ipqua, geoel, gelagq, vlave, unkno, ustar,rhust,&
sgmst,lhsrmu,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!
real*8, dimension(1:2, 1:npoin),          intent(in)::ustar
real*8, dimension(1:2, 1:npoin),          intent(out)::rhust
real*8, dimension(1:2, 1:2, 1:npoin),     intent(in)::sgmst
real*8, dimension(1:npoin),               intent(out)::lhsrmu

!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig,ishp
!...local integer array
integer,dimension(1:nvqua) :: ipq

!...local real array
real*8,dimension(1:ndegr):: b,bv
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xphq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:nq)::unknod
real*8,dimension(1:ndimn, 1:ndimn)::sigmg
real*8,dimension(1:nvqua,1:ndimn)         ::dissg
real*8,dimension(2,ngausdq)::posiq
real*8,dimension(ngausdq)  ::weighq
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
real*8, dimension(1:ndimn, 1:nvqua):: gdshp
!
real*8::eps,c00,c05,c10,c20
real*8::wi
real*8::rcoef,xg,yg,xgaus,ygaus
real*8::dxdr,dxds,dydr,dyds,djaco
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhoma,rhoad,uadv,vadv,eadv,pres,sdg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::dudr, duds, dvdr, dvds
real*8::rhoi, rhon
real*8::darea
!
data eps   / 1.0d-6 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /

!
!...Part 0: Preliminary setup...
!

!...Initialize sgmas
lhsrmu = 0.d0
rhust = 0.d0

darea = 0.d0

!...Quadrature points
call ruqope(2, ngausdq, posiq, weighq)

!...Geometry parameters for reference cell...
dr = 1.d0
ds = 1.d0
!
!...Part I: Loop every quad...
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria

!...Points consitituting one element...
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)


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
djaco = wi*(dxdr*dyds-dxds*dydr)

!...physical coordinate of gauss points....
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo

!...Riemann stress at Gauss point
sigmg = 0.d0
do ishp = 1, nvqua
sigmg(1, 1) = sigmg(1, 1) + shpq(ishp)*sgmst(1, 1, ipq(ishp))
sigmg(2, 2) = sigmg(2, 2) + shpq(ishp)*sgmst(2, 2, ipq(ishp))
enddo

!...Coefficient R of RZ or XY system...
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

!...Calculate G dot dbdx or dbdy
gdshp(1, 1:nvqua) = jacbg(1, 1)*dsprq(1:nvqua) + jacbg(1, 2)*dspsq(1:nvqua)
gdshp(2, 1:nvqua) = jacbg(2, 1)*dsprq(1:nvqua) + jacbg(2, 2)*dspsq(1:nvqua)


!...Gauss points...
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

!...Solution at the Gauss points...
unknod = 0.d0
do ideg =1,mdegr
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo

!...Primitive variables...
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
if(ielem.eq.-1)then
print*,'ielem',ig,eadv,b,unkno(:,4,ielem)
endif
!
if(nlimi.eq.6)then

! if(ie.ge.990.and.ie.le.1000) print*,'domn unk',ie,unkno(1, 1, ie)

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
sdg = sqrt(gamlg*pres/rhoad)

!...Stress tensor
!...lhs for stress and velovity
lhsrmu(ipq(1:nvqua)) =  lhsrmu(ipq(1:nvqua)) + rhoad*sdg*shpq(1:nvqua)*djaco

!...dissipation
dissg(1:nvqua,1) = gdshp(1,1:nvqua)*(sigmg(1, 1) + pres)*sqrt(geoel(3, ielem))
dissg(1:nvqua,2) = gdshp(2,1:nvqua)*(sigmg(2, 2) + pres)*sqrt(geoel(3, ielem))


!...rhs for nodal velocity
rhust(1,ipq(1:nvqua)) = rhust(1,ipq(1:nvqua)) + rhoad*sdg*uadv*shpq(1:nvqua)*djaco+dissg(1:nvqua,1)*wi
rhust(2,ipq(1:nvqua)) = rhust(2,ipq(1:nvqua)) + rhoad*sdg*vadv*shpq(1:nvqua)*djaco+dissg(1:nvqua,2)*wi
!
!if(ielem.eq.1)darea=darea+djaco
!
enddo !...(2)ig = 1,ngausd
!
!
650 enddo
!

end subroutine getriem_quad_velo_rsf

!
!...subroutine: Riemann free solver to calculate the u*  for linear quad....
!
subroutine getriem_quad_velo_rsf2(ipqua, geoel, gelagq, vlave, unkno, ustar,rhust,&
sgmst,lhsrmu,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::vlave
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!
real*8, dimension(1:2, 1:npoin),          intent(in)::ustar
real*8, dimension(1:2, 1:npoin),          intent(out)::rhust
real*8, dimension(1:2, 1:2, 1:npoin),     intent(in)::sgmst
real*8, dimension(1:npoin),               intent(out)::lhsrmu

!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig,ishp
!...local integer array
integer,dimension(1:nvqua) :: ipq

!...local real array
real*8,dimension(1:ndegr):: b,bv
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xphq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:nq)::unknod
real*8,dimension(1:ndimn, 1:ndimn)::sigmg
real*8,dimension(1:nvqua,1:ndimn)         ::dissg
real*8,dimension(2,ngausdq)::posiq
real*8,dimension(ngausdq)  ::weighq
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
real*8, dimension(1:ndimn, 1:nvqua):: gdshp
!
real*8::eps,c00,c05,c10,c20
real*8::wi
real*8::rcoef,xg,yg,xgaus,ygaus
real*8::dxdr,dxds,dydr,dyds,djaco
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhoma,rhoad,uadv,vadv,eadv,pres,sdg
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::dudr, duds, dvdr, dvds
real*8::rhoi, rhon
real*8::darea
!
data eps   / 1.0d-6 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /

!
!...Part 0: Preliminary setup...
!

!...Initialize sgmas
lhsrmu = 0.d0
rhust = 0.d0

darea = 0.d0

!...Quadrature points
call ruqope(2, ngausdq, posiq, weighq)

!...Geometry parameters for reference cell...
dr = 1.d0
ds = 1.d0
!
!...Part I: Loop every quad...
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem
!
ielem = ie + ntria

!...Points consitituting one element...
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)


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
djaco = wi*(dxdr*dyds-dxds*dydr)

!...physical coordinate of gauss points....
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo

!...Riemann stress at Gauss point
sigmg = 0.d0
do ishp = 1, nvqua
sigmg(1, 1) = sigmg(1, 1) + shpq(ishp)*sgmst(1, 1, ipq(ishp))
sigmg(2, 2) = sigmg(2, 2) + shpq(ishp)*sgmst(2, 2, ipq(ishp))
enddo

!...Coefficient R of RZ or XY system...
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

!...Calculate G dot dbdx or dbdy
gdshp(1, 1:nvqua) = jacbg(1, 1)*dsprq(1:nvqua) + jacbg(1, 2)*dspsq(1:nvqua)
gdshp(2, 1:nvqua) = jacbg(2, 1)*dsprq(1:nvqua) + jacbg(2, 2)*dspsq(1:nvqua)


!...Gauss points...
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

!...Solution at the Gauss points...
unknod = 0.d0
do ideg =1,mdegr
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo

!...Primitive variables...
if(ndens.eq.1)then
!
rhoma = unknod(1)
rhoad  = 1.d0/rhoma
!
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))

!
sdg = sqrt(gamlg*pres/rhoad)

!...Stress tensor
!...lhs for stress and velovity
lhsrmu(ipq(1:nvqua)) =  lhsrmu(ipq(1:nvqua)) + rhoad*sdg*shpq(1:nvqua)*djaco

!...dissipation
dissg(1:nvqua,1) = gdshp(1,1:nvqua)*(sigmg(1, 1) + pres)*sqrt(geoel(3, ielem))
dissg(1:nvqua,2) = gdshp(2,1:nvqua)*(sigmg(2, 2) + pres)*sqrt(geoel(3, ielem))


!...rhs for nodal velocity
rhust(1,ipq(1:nvqua)) = rhust(1,ipq(1:nvqua)) + rhoad*sdg*uadv*shpq(1:nvqua)*djaco!+dissg(1:nvqua,1)*wi
rhust(2,ipq(1:nvqua)) = rhust(2,ipq(1:nvqua)) + rhoad*sdg*vadv*shpq(1:nvqua)*djaco!+dissg(1:nvqua,2)*wi
!
!if(ielem.eq.1)darea=darea+djaco
!
enddo !...(2)ig = 1,ngausd
!
!
650 enddo
!

end subroutine getriem_quad_velo_rsf2
