!-xxxxxxxx---The details of array geoel---xxxxxxxx
!...1: rc 2: sc 3:Volume 4:Mass 5:rcv 6:scv
!...7: rcv(initital) 8: scv(initial)
!...9:  10:
!...11: overbar(r) for RZ 12: area for RZ
!...13: xi_c for AW RZ(dt) 14: eta_c for AW RZ(dt)
!
!...subroutine: Calculate P star for source term....
!
subroutine getsigmas_quadrz(ipqua, geoel, gelagq,  unkno, coord, coold, aflim, afvec, ustar, presqs)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in):: ustar
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!
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
real*8,dimension(1:nvqua)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xphq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(2, nvqua, nquad) :: presqs

!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8::acnx, acny
real*8:: dudr, duds, dvdr, dvds
real*8::rhoi, rhon
real*8::dnx,dny
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
!...physical coordinate....
!
xphq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xphq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
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
!
!if(ie.eq.1)print*,'ieeq1',vnorm(3,2,1),rcoeq(1),rcoeq(2),ipq(1:4),xphq(2, 1:nvqua)
!
!...cell averaged value...
!
if(ndens.eq.1)then
!...Specific volume...
rhomc = unkno(1, 1, ielem)
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
elseif(ndens.eq.3)then
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
!...Limiter
!
if(nlimi.eq.6)then
!
if(ndens.eq.1)then
rhomv = rhomc + aflim(1, ielem)*(unknvq(1, iv) - rhomc)
rhovt = 1.d0/rhomv
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
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the variables at the center...
!
sdctr = sqrt( max( 1d-6,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!
!...Get impedence coefficient...
!
do iv   = 1, nvqua
murie(iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
! if(ipq(iv).eq.738) print*,'murie22',murie(iv), sdctr,rhoct,deltu,vlave(1:2, ipq(iv)),unknvq(2:3, iv),dux,duy,ie
! if(ipq(iv).eq.840) then
!         print*,'murie22ang',unknvq(2,iv)**2+unknvq(3,iv)**2,ie,iv,vlave(1,ipq(iv))*unknvq(2,iv)+vlave(2,ipq(iv))*unknvq(3,iv)
! endif
!
enddo
!
!if(ie==94)print*,'vnotm',murie,rhoct,sdctr,uctr,vctr,ectr!,vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, nvqua
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!...Call Riemann solver...
!
dux = ustar(1, ipq(iv)) - unknvq(2, iv)
duy = ustar(2, ipq(iv)) - unknvq(3, iv)
!
dnx = vnorm(1, ifa, iv)
dny = vnorm(2, ifa, iv)

presqs(ifa, iv, ie) = -(sigma(1, 1, iv) + murie(iv)*(dux*dnx + duy*dny))
!
enddo

enddo
!
350 enddo  !...(1)ie = 1,nelem!

end subroutine getsigmas_quadrz
!
!...subroutine: Calculate P star for source term(triangle)....
!
subroutine getsigmas_triarz(iptri, geoel, gelag,  unkno,  coord, coold, aflim, afvec, ustar, prests)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3, 1:ngelg, 1:ntria+nbfac), intent(in)::gelag
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndimn,1:npoin),            intent(in):: ustar
real*8,dimension(1:ndimn,1:npoin),            intent(in)::coord, coold
real*8,dimension(1:nq+1,1:nsize),             intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!
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
real*8,dimension(1:nvtri)::murie
real*8,dimension(1:nvtri):: xv,  yv
real*8,dimension(1:ndimn, 1:nvtri) :: xpt,xpht
real*8,dimension(1:ndimn, 1:nvtri) :: xpti
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8,dimension(2, nvtri, ntria):: prests
!...arraies for Riemann solver
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhomv,uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv
real*8:: dudr, duds, dvdr, dvds
real*8::acnx, acny
real*8::rhoi, rhon,dnx,dny
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
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 3, 1) = gelag(1:3, 3, ie); vnorm(1:3, 4, 1) = gelag(1:3, 1, ie) !...For point ip(1)

vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 1, ie); vnorm(1:3, 4, 2) = gelag(1:3, 2, ie) !...For point ip(2)

vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 2, ie); vnorm(1:3, 4, 3) = gelag(1:3, 3, ie) !...For point ip(3)
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
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
if(nlimi.eq.6)then
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
!...Get stress tensor at nodes
!
sigma(1, 1, iv) = -pvtx
sigma(1, 2, iv) = 0.d0
sigma(2, 1, iv) = 0.d0
sigma(2, 2, iv) = -pvtx!
!
!if(ip(iv).eq.8) print*,'velocity 8',ie, rho, uvtx,vvtx,vlave(1:2, ip(iv))
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the variables at the center...
!
sdctr = sqrt( max( 1d-6,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!if(ielem.eq.5.or.ielem.eq.6) print*,'ielem56',unkno(1, 2:3, ielem),ielem
!
!...Get impedence coefficient...
!
do iv   = 1, nvtri
murie(iv) = rhoct*sdctr
enddo
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!
do iv  = 1, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
!...Call Riemann solver...
!
dux = ustar(1, ipt(iv)) - unknvt(2, iv)
duy = ustar(2, ipt(iv)) - unknvt(3, iv)
!
dnx = vnorm(1, ifa, iv)
dny = vnorm(2, ifa, iv)

prests(ifa, iv, ie) = -(sigma(1, 1, iv) + murie(iv)*(dux*dnx + duy*dny))
!
enddo
enddo
!
250 enddo  !...(1)ie = 1,nelem!

end subroutine getsigmas_triarz
!
!....source domain integral for hybrid linear quad cells(sympre)...
!
subroutine rhsdomnsrc_lagsympre_quadrz(intfac, ipqua, coord, coold, geoel,presqs, rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),     intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8,dimension(2, nvqua, nquad), intent(in) :: presqs
real*8,dimension(1:ndegr,1:nq,1:ncell),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem,ifa
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
real*8::rprz(2, 4), zprz(2, 4), xirz(2 ,4)
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
real*8:: rhoi, rhon, area
real*8:: pstar1, pstar2, pstar
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
area = 0.d0
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
rprz(1, 1) = 1.d0/3.d0*(2.d0*coord(1, ipq(1)) + coord(1, ipq(4)))
rprz(2, 1) = 1.d0/3.d0*(2.d0*coord(1, ipq(1)) + coord(1, ipq(2)))

rprz(1, 2) = 1.d0/3.d0*(2.d0*coord(1, ipq(2)) + coord(1, ipq(1)) )
rprz(2, 2) = 1.d0/3.d0*(2.d0*coord(1, ipq(2)) + coord(1, ipq(3)) )

rprz(1, 3) = 1.d0/3.d0*(2.d0*coord(1, ipq(3)) + coord(1, ipq(2)) )
rprz(2, 3) = 1.d0/3.d0*(2.d0*coord(1, ipq(3)) + coord(1, ipq(4)) )

rprz(1, 4) = 1.d0/3.d0*(2.d0*coord(1, ipq(4)) + coord(1, ipq(3)) )
rprz(2, 4) = 1.d0/3.d0*(2.d0*coord(1, ipq(4)) + coord(1, ipq(1)) )
!
zprz(1, 1) = 1.d0/3.d0*(2.d0*coord(2, ipq(1)) + coord(2, ipq(4)))
zprz(2, 1) = 1.d0/3.d0*(2.d0*coord(2, ipq(1)) + coord(2, ipq(2)))

zprz(1, 2) = 1.d0/3.d0*(2.d0*coord(2, ipq(2)) + coord(2, ipq(1)) )
zprz(2, 2) = 1.d0/3.d0*(2.d0*coord(2, ipq(2)) + coord(2, ipq(3)) )

zprz(1, 3) = 1.d0/3.d0*(2.d0*coord(2, ipq(3)) + coord(2, ipq(2)) )
zprz(2, 3) = 1.d0/3.d0*(2.d0*coord(2, ipq(3)) + coord(2, ipq(4)) )

zprz(1, 4) = 1.d0/3.d0*(2.d0*coord(2, ipq(4)) + coord(2, ipq(3)) )
zprz(2, 4) = 1.d0/3.d0*(2.d0*coord(2, ipq(4)) + coord(2, ipq(1)) )
!
do iv = 1, 4
do ifa = 1,2
xirz(ifa, iv) = sqrt(rprz(ifa, iv)**2 + zprz(ifa, iv)**2)
enddo
enddo
!
pstar1 = xirz(2, 1)*presqs(2, 1, ie) + xirz(1, 2)*presqs(1, 2, ie) +&
xirz(2, 3)*presqs(2, 3, ie) + xirz(1, 4)*presqs(1, 4, ie)
!
pstar2 = xirz(2, 1) + xirz(1, 2) +&
xirz(2, 3) + xirz(1, 4)
pstar = pstar1/pstar2
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
!
!finally, scatter the contribution to the RHS for mementum equation
!
rhsel(1,3,ielem)=rhsel(1,3,ielem) + pstar*djaco
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
end subroutine rhsdomnsrc_lagsympre_quadrz
!
!....source domain integral for hybrid linear quad cells(sympre)...
!
subroutine rhsdomnsrc_lagsympre_triarz(intfac, iptri, coord,coold, geoel,prests,rhsel )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
real*8,dimension(2, nvtri, ntria), intent(in) :: prests
real*8,dimension(1:ndegr,1:nq,1:ncell),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem,ifa
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
real*8::rprz(2, 3), zprz(2, 3), xirz(2 ,3)
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
real*8:: pstar1, pstar2, pstar
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

!...The derivatives of basis function...
!...Here dbdx means dbd(xsi), dbdy means dbd(eta)
!
rprz(1, 1) = 1.d0/3.d0*(2.d0*coord(1, ipt(1)) + coord(1, ipt(3)))
rprz(2, 1) = 1.d0/3.d0*(2.d0*coord(1, ipt(1)) + coord(1, ipt(2)))

rprz(1, 2) = 1.d0/3.d0*(2.d0*coord(1, ipt(2)) + coord(1, ipt(1)) )
rprz(2, 2) = 1.d0/3.d0*(2.d0*coord(1, ipt(2)) + coord(1, ipt(3)) )

rprz(1, 3) = 1.d0/3.d0*(2.d0*coord(1, ipt(3)) + coord(1, ipt(2)) )
rprz(2, 3) = 1.d0/3.d0*(2.d0*coord(1, ipt(3)) + coord(1, ipt(1)) )

!
zprz(1, 1) = 1.d0/3.d0*(2.d0*coord(2, ipt(1)) + coord(2, ipt(3)))
zprz(2, 1) = 1.d0/3.d0*(2.d0*coord(2, ipt(1)) + coord(2, ipt(2)))

zprz(1, 2) = 1.d0/3.d0*(2.d0*coord(2, ipt(2)) + coord(2, ipt(1)) )
zprz(2, 2) = 1.d0/3.d0*(2.d0*coord(2, ipt(2)) + coord(2, ipt(3)) )

zprz(1, 3) = 1.d0/3.d0*(2.d0*coord(2, ipt(3)) + coord(2, ipt(2)) )
zprz(2, 3) = 1.d0/3.d0*(2.d0*coord(2, ipt(3)) + coord(2, ipt(1)) )
!
do iv = 1, 3
do ifa = 1,2
xirz(ifa, iv) = sqrt(rprz(ifa, iv)**2 + zprz(ifa, iv)**2)
enddo
enddo
!
pstar1 = xirz(2, 1)*prests(2, 1, ie) +xirz(1, 1)*prests(1, 1, ie) + &
xirz(1, 2)*prests(1, 2, ie) +&
xirz(2, 3)*prests(2, 3, ie)
!
pstar2 = xirz(2, 1) + xirz(1, 2) +&
xirz(2, 3) + xirz(1, 1)
pstar = pstar1/pstar2
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
!finally, scatter the contribution to the RHS
!
rhsel(1,3,ielem)=rhsel(1,3,ielem) + pstar*djaco
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomnsrc_lagsympre_triarz
!
!...subroutine: Calculate the nodal velocity U_p^* (mass center) for hybrid RZ mesh (Area-weighted) with general Riemann solver...
!
subroutine getndvelo_lagmc_rzaw(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold, unkno,ustar, fstrzaw, fsqrzaw, aflim, afvec, itime)
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
real*8,dimension(1:ndimn,1:2,1:nvtri, 1:ntria),  intent(out)::fstrzaw !...Riemann forces
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(out)::fsqrzaw !...Riemann forces
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
allocate (munaclt(1:2, 1:2, 1:2, 1:nvtri, 1:ntria), munault(1:ndimn, 1:2, 1:nvtri,  1:ntria),&
snsigmlt(1:ndimn, 1:2,  1:nvtri,  1:ntria))
allocate (munaclq(1:2, 1:2, 1:2, 1:nvqua, 1:nquad), munaulq(1:ndimn, 1:2, 1:nvqua,  1:nquad),&
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

do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
indnd(ipf(1:nvfac)) = 1
endif
enddo
!
!...Get averaged velocity at nodes...
!
call getvlavenew(iptri, ipqua, geoel, vlave, unkno, aflim, afvec)
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
!...Zero out munacn
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
if(ntria.gt.0) call getriem_triarzaw(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold,aflim, afvec)
!
if(nquad.gt.0) call getriem_quadrzaw(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)
!
!...Third part: Impose the boundary condition
!
call getbc_lagmaire_rzaw(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)
!
!...Update the nodal velocity...
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
!if(ipoin.eq.2)print*,'ipoin',ustar(1:2, ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin),fpres(1:2, ipoin)
!
endif
enddo
!
!ustar(:,1) = 0.d0
!
!print*,'badin'
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

!...Specify boundary velocity
if(ncase.eq.14)then !...Coggeshall expansion problem
  dtime = (itime-1.d0)*dtfix
!  if(dtime.le.0.5d0)then
  ustar(1, ipf(1)) = -coord(1, ipf(1))/4.d0/(1.d0-dtime);   ustar(2, ipf(1)) = -coord(2, ipf(1))/(1.d0-dtime)
  ustar(1, ipf(2)) = -coord(1, ipf(2))/4.d0/(1.d0-dtime);   ustar(2, ipf(2)) = -coord(2, ipf(2))/(1.d0-dtime)
!  endif

!print*,'itime',dtime,ustar(1:2, ipf(1))

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
!
900 enddo

enddo !iloop
!
!print*,'ustar',ustar(:,421),coord(1:2,421)/(1.d0-dtime)
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
do ifa =1, 2
fstrzaw(1, ifa, iv, ie) = snsigmlt(1, ifa, iv, ie) + &
munaclt(1, 1, ifa, iv, ie)*ustar(1, ipt(iv))+&
munaclt(2, 1, ifa, iv, ie)*ustar(2, ipt(iv)) - munault(1, ifa, iv, ie)
fstrzaw(2, ifa, iv, ie) = snsigmlt(2, ifa, iv, ie) + &
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
do ifa =1, 2
fsqrzaw(1, ifa, iv, ie) = snsigmlq(1, ifa, iv, ie) + &
munaclq(1, 1, ifa, iv, ie)*ustar(1, ipq(iv))+&
munaclq(2, 1, ifa, iv, ie)*ustar(2, ipq(iv)) - munaulq(1, ifa, iv, ie)
fsqrzaw(2, ifa, iv, ie) = snsigmlq(2, ifa, iv, ie) + &
munaclq(1, 2,ifa, iv, ie)*ustar(1, ipq(iv))+&
munaclq(2, 2, ifa, iv, ie)*ustar(2, ipq(iv)) - munaulq(2, ifa, iv, ie)
!
enddo
!
!if(ie==21) print*,'fire',iv,ie,snsigmlq(1,1:2,iv,ie),munaclq(1:2,iv,ie),ipq(:),munaulq(1,1:2,iv,ie)
!
enddo
!
enddo
!
deallocate (munacn, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lagmc_rzaw
!
!...subroutine: Calculate the Riemann input for hybrid tria grids general Riemann solver(area-weighted RZ)....
!
subroutine getriem_triarzaw(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
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
real*8, dimension(1:2, 1:2, 1:2,  1:nvtri, 1:ntria),       intent(out)::munaclt
real*8, dimension(1:ndimn, 1:2, 1:nvtri,  1:ntria), intent(out)::munault
real*8, dimension(1:ndimn, 1:2, 1:nvtri,  1:ntria), intent(out)::snsigmlt
!...Local integer
integer::ie, ideg, ielem, ifa, iv
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:3, 1:nvtri)::bt, btv
real*8,dimension(1:nq,1:nvtri)::unknvt
real*8::aujmp(1:3, 1:nvtri)
real*8::vnorm(1:3, 1:2, 1:nvtri)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8,dimension(1:2,1:nvtri)::murie
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
real*8::rhoi, rhon
!
data eps   / 1.0d-08 /
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
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
!
!...ndA=0.5d0*vnorm
!
!vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
vnorm(3, 1, 1) = 0.5d0*vnorm(3, 1, 1)
vnorm(3, 2, 1) = 0.5d0*vnorm(3, 2, 1)
!
vnorm(3, 1, 2) = 0.5d0*vnorm(3, 1, 2)
vnorm(3, 2, 2) = 0.5d0*vnorm(3, 2, 2)
!
vnorm(3, 1, 3) = 0.5d0*vnorm(3, 1, 3)
vnorm(3, 2, 3) = 0.5d0*vnorm(3, 2, 3)
!
!...cell averaged value...
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
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
!...zero out unknv
!
unknvt = 0.d0
!
do iv   = 1,nvtri
 do ideg = 1,mdegr
   unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
 enddo
!
 if(ndens.eq.1)then
  rhovt  = 1.d0/unknvt(1, iv)
 elseif(ndens.eq.2)then
  r = xv(iv); s= yv(iv)

  xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
  xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
  xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
  xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
  call  getrhoig_tria(rhoi, r, s, xpti)
  call getdensity_triallnl(r, s, xpt, xpti, rhoi, rhon)
!
  rhovt = rhon
 elseif(ndens.eq.3)then
  rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
  btv(1, iv) = 1.d0
  btv(2, iv) = (xv(iv)-rcv)/dr
  btv(3, iv) = (yv(iv)-scv)/ds
!
 unknvt(1, iv) =0.d0

 do ideg = 1,mdegr
  unknvt(1, iv) = unknvt(1, iv) + unkno(ideg,1,ielem)*btv(ideg, iv)
 enddo

 rhovt = unknvt(1, iv)
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
if(nlimi.eq.6)then
 if(ndens.eq.1)then
   rhomv = rhomc + aflim(1, ielem)*(unknvt(1, iv) - rhomc)
   rhovt = 1.d0/rhomv
 elseif(ndens.eq.2)then
   rhovt = 1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc)
 elseif(ndens.eq.3)then
   rhovt = 1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc)
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
 pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)
!...updtae unknv(2:3,:)
 unknvt(2, iv) = uvtx
 unknvt(3 ,iv) = vvtx
endif
!
!...Get stress tensor at nodes
!
sigma(1, 1, iv) = -pvtx
sigma(1, 2, iv) = 0.d0
sigma(2, 1, iv) = 0.d0
sigma(2, 2, iv) = -pvtx!
!
!...Get the a_c (unit vector)
!
aujmp(1:2, iv) = vlave(1:2, ipt(iv)) - unknvt(2:3, iv)
acnx = aujmp(1, iv)
acny = aujmp(2, iv)
!
 if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
  aujmp(1:2, iv) = 1.e-11
 else
  aujmp(1:2, iv) = aujmp(1:2, iv)/sqrt(acnx**2 + acny**2)
 endif

  aujmp(3, iv) = sqrt(acnx**2 + acny**2)
enddo
!
!...Get the variables at the center...
!
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
aujmp(3,:)=aujmp(3,:)/sdctr
!
!...Get impedence coefficient...
!
do iv   = 1, nvtri
dux= vlave(1, ipt(iv))-unknvt(2, iv)
duy= vlave(2, ipt(iv))-unknvt(3, iv)
deltu = sqrt(dux**2 + duy**2)
do ifa = 1,2
deltu = 0.d0*abs(dux*vnorm(1,ifa, iv) + duy*vnorm(2,ifa,iv))
murie(ifa, iv) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u

!if(ipt(iv).eq.2) print*,'iptv',ie,rhoct*sdctr,rhoct,sdctr,unkno(1, 1, ielem)
enddo
!
enddo
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!
do iv  = 1, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
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
!if(ipt(iv).eq.2) print*,'The output for Riemann solver',sigma(1,1,iv),munacu_rie(1:2), murie(ifa,iv),unknvt(2:3,iv),&
!                    vnorm(1:3,ifa,iv),snsigm(1:2, ipt(iv)),ie, ifa,iv

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

end subroutine getriem_triarzaw
!
!...subroutine: Calculate the Riemann input for hybrid quad grids general Riemann solver (area-weighted RZ)....
!
subroutine getriem_quadrzaw(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
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
real*8, dimension(1:2, 1:2, 1:2,  1:nvqua, 1:nquad),      intent(out)::munaclq
real*8, dimension(1:ndimn, 1:2,  1:nvqua, 1:nquad), intent(out)::munaulq
real*8, dimension(1:ndimn, 1:2,  1:nvqua, 1:nquad), intent(out)::snsigmlq
!...Local integer
integer::ie, ideg, ielem, ifa, iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!...local real array
real*8,dimension(1:3, 1:nvqua)::bq,bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8::aujmp(1:3, 1:nvqua)
real*8::vnorm(1:3, 1:2, 1:nvqua)
real*8::sigma(1:2, 1:2, 1:nvqua)
real*8,dimension(1:2, 1:nvqua)::murie
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
real*8::ax,bx
real*8::rhoi, rhon
!
data eps   / 1.0d-08 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
ax= 0.d0
bx = 0.d0
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
rcoeq(1:nvqua) = 1.d0 - alfrz + alfrz*xphq(2, 1:nvqua)
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
vnorm(3, 1, 1) = 0.5d0*vnorm(3, 1, 1)
vnorm(3, 2, 1) = 0.5d0*vnorm(3, 2, 1)
!
vnorm(3, 1, 2) = 0.5d0*vnorm(3, 1, 2)
vnorm(3, 2, 2) = 0.5d0*vnorm(3, 2, 2)
!
vnorm(3, 1, 3) = 0.5d0*vnorm(3, 1, 3)
vnorm(3, 2, 3) = 0.5d0*vnorm(3, 2, 3)
!
vnorm(3, 1, 4) = 0.5d0*vnorm(3, 1, 4)
vnorm(3, 2, 4) = 0.5d0*vnorm(3, 2, 4)
!
!...cell averaged value...
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

!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
 do ideg = 1,mdegr
  unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
 enddo
!
if(ndens.eq.1)then
 rhovt  = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
 r = xvq(iv); s= yvq(iv)

 xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
 xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
 xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
 xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
 call getrhoig_quad(rhoi, r, s, xpqi)!
 call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
 rhovt = rhon
elseif(ndens.eq.3)then
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
 rhovt  = unknvq(1, iv)
endif
!
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!...Limiter
!
if(nlimi.eq.6)then
 if(ndens.eq.1)then
  rhomv = rhomc + aflim(1, ielem)*(unknvq(1, iv) - rhomc)
  rhovt = 1.d0/rhomv
 elseif(ndens.eq.2)then
 rhovt = 1.d0/rhomc + aflim(1, ielem)*(rhovt - 1.d0/rhomc)
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
endif
!
!...Get stress tensor at nodes
!
sigma(1, 1, iv) = -pvtx
sigma(1, 2, iv) = 0.d0
sigma(2, 1, iv) = 0.d0
sigma(2, 2, iv) = -pvtx!
!
!...Get the a_c (unit vector)
!
aujmp(1:2, iv) = vlave(1:2, ipq(iv)) - unknvq(2:3, iv)
acnx = aujmp(1, iv)
acny = aujmp(2, iv)
!
 if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
  aujmp(1:2, iv) = 1.e-11
 else
  aujmp(1:2, iv) = aujmp(1:2, iv)/sqrt(acnx**2 + acny**2)
 endif

aujmp(3, iv) = sqrt(acnx**2 + acny**2)
enddo
!
!...Get the variables at the center...
!
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
do ifa = 1,2
deltu = 0.d0*abs(dux*vnorm(1,ifa, iv) + duy*vnorm(2,ifa,iv))
murie(ifa, iv) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
enddo
enddo
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!
do iv  = 1, nvqua
do ifa = 1, 2 !...Every corner consists of 2 faces...
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
 !if(ipq(iv).eq.2) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(1,1,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ipq(iv)),ie, ifa,iv
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

end subroutine getriem_quadrzaw
!
!...Face integral (mass center) for hybrid quad in R-Z (area-weighted)...
!
subroutine rhsifacedg_lagmc_quadrzaw(ipqua, unkno, ustar,fsqrzaw, gelagq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvqua, 1:nquad),  intent(in)::fsqrzaw !...Riemann forces
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
fsqrzaw(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
fsqrzaw(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv)))

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fsqrzaw(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
fsqrzaw(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv)))
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(iv))*fsqrzaw(1, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
ustar(2, ipq(iv))*fsqrzaw(2, 1, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(1,iv))) +&
ustar(1, ipq(iv))*fsqrzaw(1, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv))) +&
ustar(2, ipq(iv))*fsqrzaw(2, 2, iv, ie)*c13*(2.d0*bq(1:ndegr, iv) + bq(1:ndegr, ipfq(2,iv)))
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
!if(ie==1)  print*,'rhs iface',ielem, ie,ulnpn(1)!,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie)
!

650 enddo
!
end subroutine rhsifacedg_lagmc_quadrzaw
!
!...Face integral (mass center) for hybrid triangle in RZ (area-weighted)...
!
subroutine rhsifacedg_lagmc_triarzaw(iptri, unkno, ustar, fstrzaw, gelag, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri,1:ntria),  intent(in)::fstrzaw !...Riemann forces
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
fstrzaw(1, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
fstrzaw(1, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv)))

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstrzaw(2, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
fstrzaw(2, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv)))
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipt(iv))*fstrzaw(1, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
ustar(2, ipt(iv))*fstrzaw(2, 1, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(1,iv))) +&
ustar(1, ipt(iv))*fstrzaw(1, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv))) +&
ustar(2, ipt(iv))*fstrzaw(2, 2, iv, ie)*c13*(2.d0*bt(1:ndegr, iv) + bt(1:ndegr, ipf(2,iv)))
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
end subroutine rhsifacedg_lagmc_triarzaw
!
!....domain integral for hybrid linear triangle cells for R-Z (Area-weighted)
!
subroutine rhsdomndg_lagmc_triarzaw(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
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
djaco = 0.5d0*wi
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
enddo
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
  rhoma = unknod(1)
  rhoad  = 1.d0/rhoma
elseif(ndens.eq.2)then

  xpi(1, 1:nvtri) = coold(1, ipt(1:nvtri))
  xpi(2, 1:nvtri) = coold(2, ipt(1:nvtri))

  call getrhoig_tria(rhoi, r, s, xpi)
  call getdensity_triallnl(r, s, xp, xpi, rhoi, rhon)

  rhoma = 1.d0/rhon
  rhoad = rhon

elseif(ndens.eq.3)then
 rcv = geoel(5, ielem); scv = geoel(6, ielem)

 bv(1) = 1.d0
 bv(2) = (xg-rcv)/dr
 bv(3) = (yg-scv)/ds

 unknod(1) =0.d0
 do ideg = 1,mdegr
  unknod(1) = unknod(1) + unkno(ideg,1,ielem)*bv(ideg)
 enddo
 rhoma = 1.d0/unknod(1)
 rhoad  = unknod(1)
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!...Impose imiter...
!
if(nlimi.eq.1)then

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

elseif(nlimi.eq.6)then
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
end subroutine rhsdomndg_lagmc_triarzaw
!
!....domain integral for hybrid linear quad cells for R-Z (Area-weighted)
!
subroutine rhsdomndg_lagmc_quadrzaw(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )
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
djaco = wi
!
!...physical coordinate of gauss points....
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, 4
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
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
 rhoma = unknod(1)
 rhoad  = 1.d0/rhoma

elseif(ndens.eq.2)then
 xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
 xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
 call getrhoig_quad(rhoi, r, s, xpqi)!
 call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
 rhoma = 1.d0/rhon
 rhoad = rhon

elseif(ndens.eq.3)then
!
 rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
 bv(1) = 1.d0
 bv(2) = (xg-rcv)/dr
 bv(3) = (yg-scv)/ds

 unknod(1) =0.d0

 do ideg = 1,mdegr
 unknod(1) = unknod(1) + unkno(ideg,1,ielem)*bv(ideg)
 enddo

 rhoma = 1.d0/unknod(1)
 rhoad  = unknod(1)

endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!if(ielem.eq.9) print*,'domn2',pres,rhoad,unkno(1:3,1,ielem),b(1:3),xg,yg
!
if(nlimi.eq.1)then
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
!finally, collect the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) - fluxd(ideg,1:nq)*djaco
enddo
!
!if(ielem.eq.9) print*,'domn',-fluxd(3,3)*2.d0,gdshp(2, 3),-pres,ig
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
!
end subroutine rhsdomndg_lagmc_quadrzaw
!
!....Source term domain integral  for hybrid linear triangle cells for R-Z (Area-weighted)
!
subroutine rhsdomnsrcdg_lagmc_triarzaw(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
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
djaco = 0.5d0*wi*(dxdr*dyds - dxds*dydr)
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgaus = xgaus + shp(ishp)*xp(1,ishp)
ygaus = ygaus + shp(ishp)*xp(2,ishp)
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
 rhoma = unknod(1)
 rhoad  = 1.d0/rhoma
elseif(ndens.eq.2)then
!
 xpi(1, 1:nvtri) = coold(1, ipt(1:nvtri))
 xpi(2, 1:nvtri) = coold(2, ipt(1:nvtri))

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
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!...impose imiter...
!
if(nlimi.eq.6)then
!
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
uadv = unkno(1,2,ielem)  !+ dudr*b(2) + duds*b(3)
vadv = unkno(1,3,ielem)  !+ dvdr*b(2) + dvds*b(3)
!
pres = pctr + aflim(4, ielem)*(pres- pctr)
!
endif
!
fluxd(1,1) = vadv*b(1)
fluxd(2,1) = vadv*b(2)
fluxd(3,1) = vadv*b(3)
!
fluxd(1,2) = 0.d0
fluxd(2,2) = 0.d0
fluxd(3,2) = 0.d0
!
fluxd(1,3) = 0.d0
fluxd(2,3) = 0.d0
fluxd(3,3) = 0.d0
!
fluxd(1,4) = vadv*(-pres)*b(1)
fluxd(2,4) = vadv*(-pres)*b(2)
fluxd(3,4) = vadv*(-pres)*b(3)
!
!finally, collect the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) + fluxd(ideg,1:nq)*djaco
enddo
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomnsrcdg_lagmc_triarzaw
!
!....source term domain integral for hybrid linear quad cells for R-Z (Area-weighted)
!
subroutine rhsdomnsrcdg_lagmc_quadrzaw(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec )
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
real*8:: veloy
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
!...Zero veloy
!
veloy = 0.d0
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
djaco = wi*(dxdr*dyds - dxds*dydr)
!
!...physical coordinate of gauss points....
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, 4
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
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
endif
!
uadv = unknod(2)
vadv = unknod(3)
eadv = unknod(4)
pres = (gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
if(nlimi.eq.6)then
!
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
uadv = unkno(1,2,ielem)  !+ dudr*b(2) + duds*b(3)
vadv = unkno(1,3,ielem)  !+ dvdr*b(2) + dvds*b(3)
!
pres = pctr + aflim(4, ielem)*(pres- pctr)
!
endif
!
fluxd(1,1) = vadv*b(1)
fluxd(2,1) = vadv*b(2)
fluxd(3,1) = vadv*b(3)
!
fluxd(1,2) = 0.d0
fluxd(2,2) = 0.d0
fluxd(3,2) = 0.d0
!
fluxd(1,3) = 0.d0
fluxd(2,3) = 0.d0
fluxd(3,3) = 0.d0
!
fluxd(1,4) = vadv*(-pres)*b(1)
fluxd(2,4) = vadv*(-pres)*b(2)
fluxd(3,4) = vadv*(-pres)*b(3)
!
!finally, collect the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) + fluxd(ideg,1:nq)*djaco
enddo
!
!veloy = veloy + vadv*djaco
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
!
end subroutine rhsdomnsrcdg_lagmc_quadrzaw
!
!...Get the mass matrix for lagrangian based on mass center...
!
subroutine  getamatr_lagmc_rzaw(unkno,amatr,geoel,coord,inpoel, iptri, ipqua, aflim)
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
real*8,dimension(1:nq+1, 1:nsize),intent(in)::aflim
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
real*8::rhom, rho0, rho0m
real*8::wi,djaco, volel,det
real*8::c10
real*8::f0,f1,f2,f3,f4
real*8::f5,f6,f7,f8
real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::b2,b3,b4,b5,b6
real*8::radie, radii,radie2,radii2,radic2, radig2,sentr,rhoin,rhoex,rho0ba
real*8::masel,xgaus,ygaus
real*8::rcoef
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
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 0.5d0
ds = 0.5d0!
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
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
rho0m = unkno(1, 1, ielem) + unkno(2, 1, ielem)*b2 + unkno(3, 1, ielem)*b3
!
if(nlimi.eq.6)then
rho0m = unkno(1, 1, ielem) + aflim(1, ielem)*(rho0m - unkno(1, 1, ielem))
endif
!
rho0 = 1.d0/rho0m
!
!...Coefficient R of RZ or XY system...
!
!
f0 = f0 + rho0*djaco
f1 = f1 + rho0*(xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + rho0*(xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + rho0*(yg-sc)/ds*(yg-sc)/ds*djaco

if(npoly==2)then
f22 = f22 + rho0*b2*b2*djaco
f23 = f23 + rho0*b2*b3*djaco
f24 = f24 + rho0*b2*b4*djaco
f25 = f25 + rho0*b2*b5*djaco
f26 = f26 + rho0*b2*b6*djaco

f33 = f33 + rho0*b3*b3*djaco
f34 = f34 + rho0*b3*b4*djaco
f35 = f35 + rho0*b3*b5*djaco
f36 = f36 + rho0*b3*b6*djaco

f44 = f44 + rho0*b4*b4*djaco
f45 = f45 + rho0*b4*b5*djaco
f46 = f46 + rho0*b4*b6*djaco

f55 = f55 + rho0*b5*b5*djaco
f56 = f56 + rho0*b5*b6*djaco

f66 = f66 + rho0*b6*b6*djaco
endif
!   if(ie==400) print*,'gauss',xg, yg, f0
enddo
!
!print*,'ielem', ie, f0, volel

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
!...
!   if(ie==2) print*,'ie',mmatr(:,:)
x5 = 0.d0
b55 = 0.d0
call getinvmat(5, mmatr, x5, b55)
!
!   if(ie==2) print*,'ie',x5(:,:)
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
amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)
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
!
rho0m = unkno(1, 1, ielem) + unkno(2, 1, ielem)*b2 + unkno(3, 1, ielem)*b3
!
if(nlimi.eq.6)then
rho0m = unkno(1, 1, ielem) + aflim(1, ielem)*(rho0m - unkno(1, 1, ielem))
endif
!
rho0 = 1.d0/rho0m
!
!...Coefficient R of RZ or XY system...
!
f0 = f0 + rho0*djaco
f1 = f1 + rho0*(xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + rho0*(xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + rho0*(yg-sc)/ds*(yg-sc)/ds*djaco

if(npoly==2)then
f22 = f22 + rho0*b2*b2*djaco
f23 = f23 + rho0*b2*b3*djaco
f24 = f24 + rho0*b2*b4*djaco
f25 = f25 + rho0*b2*b5*djaco
f26 = f26 + rho0*b2*b6*djaco

f33 = f33 + rho0*b3*b3*djaco
f34 = f34 + rho0*b3*b4*djaco
f35 = f35 + rho0*b3*b5*djaco
f36 = f36 + rho0*b3*b6*djaco

f44 = f44 + rho0*b4*b4*djaco
f45 = f45 + rho0*b4*b5*djaco
f46 = f46 + rho0*b4*b6*djaco

f55 = f55 + rho0*b5*b5*djaco
f56 = f56 + rho0*b5*b6*djaco

f66 = f66 + rho0*b6*b6*djaco
endif
!   if(ie==400) print*,'gauss',xg, yg, f0
enddo
!
!print*,'ielem', ie, f0, volel

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
!...
!   if(ie==2) print*,'ie',mmatr(:,:)
x5 = 0.d0
b55 = 0.d0
call getinvmat(5, mmatr, x5, b55)
!
!   if(ie==2) print*,'ie',x5(:,:)
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
amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)
!
!if(ielem.eq.332.or.ielem.eq.333)print*,'amatr+rzaw',ielem,unkno(1:3,1,ielem),amatr(4, ielem)*geoel(11, ielem),&
!geoel(11, ielem)/sin(23.d0/64.d0*pi),&
!geoel(11, ielem)/sin(25.d0/64.d0*pi)
!
enddo !...(2)ie = 1,nelem

end subroutine  getamatr_lagmc_rzaw
!
!...Get the mass matrix for lagrangian based on the updated mass center (revised version)...
!
subroutine  getamatr_lagmc_rzawrev(unkno,amatr,geoel,coord,inpoel, iptri, ipqua, aflim)
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
real*8,dimension(1:nq+1, 1:nsize),intent(in)::aflim
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp, ielem, ideg
!...Local real array
real*8::xp(1:2, 1:nptri)
real*8,dimension(1:2, 1:npqua)::xpq
real*8,dimension(1:nptri)::shp, dspr, dsps
real*8:: weight(ngausd), posit(2, ngausd)
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:mdegr)::bt, btv
real*8,dimension(1:mdegr)::bq, bqv
!...Local real
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8::rcv, scv
real*8::xc, yc
real*8:: dxdr,dxds,dydr,dyds
real*8::rhog,rhomg,rhomc
real*8::wi,djaco, volel,det
real*8::c10
real*8::f0,f1,f2,f3,f4
real*8::f5,f6,f7,f8
real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::radie, radii,radie2,radii2,radic2, radig2,sentr,rhoin,rhoex,rho0ba
real*8::masel,xgaus,ygaus
real*8::rcoef
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
!...is stored in the last term of amatr ...
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
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 0.5d0
ds = 0.5d0!
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
!
bt(1) = 1.d0
bt(2) = (r - rc)/dr
bt(3) = (s - sc)/ds
!
!...Density calculations with different density evolutions
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
rhomg = 0.d0
do ideg = 1,mdegr
rhomg = rhomg + unkno(ideg,1,ielem)*bt(ideg)
enddo

rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
rcv = geoel(5, ielem); scv = geoel(6, ielem)

btv(1) = 1.d0
btv(2) = (r-rcv)/dr
btv(3) = (s-scv)/ds

rhog = 0.d0
do ideg = 1,mdegr
rhog = rhog + unkno(ideg,1,ielem)*btv(ideg)
enddo

endif
!
if(nlimi.eq.6)then
if(ndens.eq.1)then
rhomg = rhomc + aflim(1 ,ielem)*(rhomg - rhomc)
rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
!
elseif(ndens.eq.3)then
rhog = 1.d0/rhomc + aflim(1 ,ielem)*(rhog - 1.d0/rhomc)
endif
endif
!
!...Coefficient R of RZ or XY system...
!
f0 = f0 + rhog*djaco
f1 = f1 + rhog*(xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + rhog*(xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + rhog*(yg-sc)/ds*(yg-sc)/ds*djaco
!
!if(ielem.eq.1)print*,'f0',f0,f1,f2,f3,ig,btv(3)

if(npoly==2)then
f22 = f22 + rhog*bt(2)*bt(2)*djaco
f23 = f23 + rhog*bt(2)*bt(3)*djaco
f24 = f24 + rhog*bt(2)*bt(4)*djaco
f25 = f25 + rhog*bt(2)*bt(5)*djaco
f26 = f26 + rhog*bt(2)*bt(6)*djaco

f33 = f33 + rhog*bt(3)*bt(3)*djaco
f34 = f34 + rhog*bt(3)*bt(4)*djaco
f35 = f35 + rhog*bt(3)*bt(5)*djaco
f36 = f36 + rhog*bt(3)*bt(6)*djaco

f44 = f44 + rhog*bt(4)*bt(4)*djaco
f45 = f45 + rhog*bt(4)*bt(5)*djaco
f46 = f46 + rhog*bt(4)*bt(6)*djaco

f55 = f55 + rhog*bt(5)*bt(5)*djaco
f56 = f56 + rhog*bt(5)*bt(6)*djaco

f66 = f66 + rhog*bt(6)*bt(6)*djaco
endif
!
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

!...Get the inverse matrix

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
!1:nmatr
amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)
!amatr(1:nmatr-1, ielem) = amatr(1:nmatr-1, ielem)/geoel(11, ielem)
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
bq(1) = 1.d0
bq(2) = (r - rc)/dr
bq(3) = (s - sc)/ds
!
!...Density calculations with different density evolutions
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
rhomg = 0.d0
do ideg = 1,mdegr
rhomg = rhomg + unkno(ideg,1,ielem)*bq(ideg)
enddo

rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bqv(1) = 1.d0
bqv(2) = (r-rcv)/dr
bqv(3) = (s-scv)/ds

rhog = 0.d0
do ideg = 1,mdegr
rhog = rhog + unkno(ideg,1,ielem)*bqv(ideg)
enddo

endif
!
if(nlimi.eq.6)then
if(ndens.eq.1)then
rhomg = rhomc + aflim(1 ,ielem)*(rhomg - rhomc)
rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
!
elseif(ndens.eq.3)then
rhog = 1.d0/rhomc + aflim(1 ,ielem)*(rhog - 1.d0/rhomc)
endif
endif
!
!...Coefficient R of RZ or XY system...
!
f0 = f0 + rhog*djaco
f1 = f1 + rhog*(xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + rhog*(xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + rhog*(yg-sc)/ds*(yg-sc)/ds*djaco

if(npoly==2)then
f22 = f22 + rhog*bq(2)*bq(2)*djaco
f23 = f23 + rhog*bq(2)*bq(3)*djaco
f24 = f24 + rhog*bq(2)*bq(4)*djaco
f25 = f25 + rhog*bq(2)*bq(5)*djaco
f26 = f26 + rhog*bq(2)*bq(6)*djaco

f33 = f33 + rhog*bq(3)*bq(3)*djaco
f34 = f34 + rhog*bq(3)*bq(4)*djaco
f35 = f35 + rhog*bq(3)*bq(5)*djaco
f36 = f36 + rhog*bq(3)*bq(6)*djaco

f44 = f44 + rhog*bq(4)*bq(4)*djaco
f45 = f45 + rhog*bq(4)*bq(5)*djaco
f46 = f46 + rhog*bq(4)*bq(6)*djaco

f55 = f55 + rhog*bq(5)*bq(5)*djaco
f56 = f56 + rhog*bq(5)*bq(6)*djaco

f66 = f66 + rhog*bq(6)*bq(6)*djaco
endif
!
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

!...Get the inverse matrix

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
amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)
!amatr(1:nmatr-1, ielem) = amatr(1:nmatr-1, ielem)/geoel(11, ielem)
!
!if(ielem.eq.332.or.ielem.eq.333)print*,'amatr+rzaw',ielem,unkno(1:3,1,ielem),amatr(4, ielem)*geoel(11, ielem)
!
enddo !...(2)ie = 1,nelem

end subroutine  getamatr_lagmc_rzawrev
!
!...Find mass center in geoel for curved cell...
!
subroutine getgeoel_lagmc_rzaw(inpoel, iptri, ipqua, geoel, coord, unkno, aflim)
use constant
implicit none
integer*4,dimension(1:nvtri, 1:ntria),intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:ngeel, 1:nsize), intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:nq+1, 1:nsize),intent(in)::aflim
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
real*8,dimension(1:3)::bt
real*8,dimension(1:3)::bq

!...local real number
real*8:: dxdr,dxds,dydr,dyds
real*8:: xc, yc, xcel, ycel
real*8:: eps,c00,c10,c05,c20
real*8:: r, s, djaco, volel, masel
real*8:: wi, xcent, ycent,xgaus, ygaus
real*8:: rhomc, rhogc, rhog
real*8::rcoef !...Coefficient R of RZ or XY coordinate system...
real*8::areel !...2D area of the cell
!
real*8:: dr, ds, rc ,sc
!
integer::ielem, igaus, ishp,ideg, iloop
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
!
do iloop = 1,5
!...
xcent = 0.d0
ycent = 0.d0
xcel  = 0.d0
ycel  = 0.d0
volel = 0.d0
masel = 0.d0
areel = 0.d0
!
!...Mass averaged center for Area-weighted RZ
!
rc = geoel(1, ielem)
sc = geoel(2, ielem)
!
dr = 0.5d0
ds = 0.5d0
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
!...Restore r, s
!
r  = posit(1,igaus)
s  = posit(2,igaus)
!
bt(1) = 1.d0
bt(2) = (r - rc)/dr
bt(3) = (s - sc)/ds
!
rhogc = 0.d0
!
do ideg = 1,mdegr
rhogc = rhogc + unkno(ideg,1,ielem)*bt(ideg)
enddo
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
endif
!
if(nlimi.eq.6)then
rhogc = rhomc + aflim(1, ielem)*(rhogc - rhomc)
endif
!
if(ndens.eq.1)then
rhog = 1.d0/rhogc
endif
!
!...Mass...
!
masel = masel + 1.d0/rhomc*djaco
!
xcent = xcent + rhog*r*djaco!...mass center at the reference element...
ycent = ycent + rhog*s*djaco
!
xcel = xcel + 1.d0/rhomc*r*djaco
ycel = ycel + 1.d0/rhomc*s*djaco
!
enddo
!
geoel(1, ielem) = xcent/masel
geoel(2, ielem) = ycent/masel
!
geoel(15, ielem) = xcel/masel
geoel(16, ielem) = ycel/masel
!
enddo !iloop(Iterate to get the 2nd order mass cenetr)
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
!...Mass averaged center for Area-weighted RZ
!
do iloop = 1,5
!
xcent = 0.d0
ycent = 0.d0
xcel = 0.d0
ycel = 0.d0
volel = 0.d0
masel = 0.d0
areel = 0.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 1.d0
ds = 1.d0!
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
!...Restore r, s
r  = posiq(1,igaus)
s  = posiq(2,igaus)
!
bq(1) = 1.d0
bq(2) = (r - rc)/dr
bq(3) = (s - sc)/ds
!
rhogc = 0.d0
!
do ideg = 1,mdegr
rhogc = rhogc + unkno(ideg,1,ielem)*bq(ideg)
enddo
!
!...cell averaged value...
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
endif
!
if(iloop.ne.1)then
if(nlimi.eq.6)then
rhogc = rhomc + aflim(1, ielem)*(rhogc - rhomc)
endif
else
if(nlimi.eq.6)then
rhogc = rhomc !+ aflim(1, ielem)*(rhogc - rhomc)
endif
endif
!
if(ndens.eq.1)then
rhog = 1.d0/rhogc
endif
!
!...Mass...
!
masel = masel + 1.d0/rhomc*djaco
!
!if(ielem.eq.1) print*,'bad',djaco,xpq(1, 1:npqua)
!
xcent = xcent + rhog*r*djaco !...mass center at the reference element...
ycent = ycent + rhog*s*djaco
!
xcel = xcel + 1.d0/rhomc*r*djaco
ycel = ycel + 1.d0/rhomc*s*djaco
!
!if(ie.eq.1) print*,'ieeq1',djaco,r,s
!
enddo
!
geoel(1, ielem) = xcent/masel
geoel(2, ielem) = ycent/masel
!
geoel(15, ielem) = xcel/masel
geoel(16, ielem) = ycel/masel
!
enddo !iloop
!
enddo !...(1)ie = 1,nelem
!
!print*,'geoel',geoel(1:8,1)
!
end subroutine getgeoel_lagmc_rzaw
!
!...Find mass center in geoel for curved cell...
!
subroutine getgeoel_lagmc_rzawrev(inpoel, iptri, ipqua, geoel, coord, unkno, aflim)
use constant
implicit none
integer*4,dimension(1:nvtri, 1:ntria),intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:ngeel, 1:nsize), intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:nq+1, 1:nsize),intent(in)::aflim
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
real*8,dimension(1:3)::bt,btv
real*8,dimension(1:3)::bq,bqv

!...local real number
real*8:: dxdr,dxds,dydr,dyds
real*8:: xc, yc, xcel, ycel
real*8:: eps,c00,c10,c05,c20
real*8:: r, s, djaco, volel, masel
real*8:: wi, xcent, ycent,xgaus, ygaus
real*8:: rhomc, rhogc, rhog
real*8:: rhomg,rcv,scv
real*8::rcoef !...Coefficient R of RZ or XY coordinate system...
real*8::areel !...2D area of the cell
!
real*8:: dr, ds, rc ,sc
!
integer::ielem, igaus, ishp,ideg, iloop
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
!
do iloop = 1,5
!...
xcent = 0.d0
ycent = 0.d0
xcel  = 0.d0
ycel  = 0.d0
volel = 0.d0
masel = 0.d0
areel = 0.d0
!
!...Mass averaged center for Area-weighted RZ
!
rc = geoel(1, ielem)
sc = geoel(2, ielem)
!
dr = 0.5d0
ds = 0.5d0
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
!...Restore r, s
!
r  = posit(1,igaus)
s  = posit(2,igaus)
!
bt(1) = 1.d0
bt(2) = (r - rc)/dr
bt(3) = (s - sc)/ds
!
!...Density calculations with different density evolutions
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
rhomg = 0.d0
do ideg = 1,mdegr
rhomg = rhomg + unkno(ideg,1,ielem)*bt(ideg)
enddo

if(iloop.eq.1)then
rhomg  = rhomc
endif

rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
rcv = geoel(5, ielem); scv = geoel(6, ielem)

btv(1) = 1.d0
btv(2) = (r - rcv)/dr
btv(3) = (s - scv)/ds

rhog = 0.d0
do ideg = 1,mdegr
rhog = rhog + unkno(ideg,1,ielem)*btv(ideg)
enddo

if(iloop.eq.1)then
rhog  = 1.d0/rhomc
endif
endif
!
if(nlimi.eq.6)then
if(ndens.eq.1)then
rhomg = rhomc + aflim(1 ,ielem)*(rhomg - rhomc)
rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
!
elseif(ndens.eq.3)then
rhog = 1.d0/rhomc + aflim(1 ,ielem)*(rhog - 1.d0/rhomc)
endif
endif
!
!...Mass...
!
masel = masel + 1.d0/rhomc*djaco
!
xcent = xcent + rhog*r*djaco!...mass center at the reference element...
ycent = ycent + rhog*s*djaco
!
xcel = xcel + 1.d0/rhomc*r*djaco
ycel = ycel + 1.d0/rhomc*s*djaco
!
enddo
!
geoel(1, ielem) = xcent/masel
geoel(2, ielem) = ycent/masel
!
geoel(15, ielem) = xcel/masel
geoel(16, ielem) = ycel/masel
!
enddo !iloop(Iterate to get the 2nd order mass cenetr)
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
!...Mass averaged center for Area-weighted RZ
!
do iloop = 1,5
!
xcent = 0.d0
ycent = 0.d0
xcel = 0.d0
ycel = 0.d0
volel = 0.d0
masel = 0.d0
areel = 0.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 1.d0
ds = 1.d0!
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
!...Restore r, s
r  = posiq(1,igaus)
s  = posiq(2,igaus)
!
bq(1) = 1.d0
bq(2) = (r - rc)/dr
bq(3) = (s - sc)/ds
!
!...Density calculations with different density evolutions
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
rhomg = 0.d0
do ideg = 1,mdegr
rhomg = rhomg + unkno(ideg,1,ielem)*bq(ideg)
enddo

if(iloop.eq.1)then
rhomg  = rhomc
endif

rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
!
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bqv(1) = 1.d0
bqv(2) = (r-rcv)/dr
bqv(3) = (s-scv)/ds

rhog = 0.d0
do ideg = 1,mdegr
rhog = rhog + unkno(ideg,1,ielem)*bqv(ideg)
enddo

if(iloop.eq.1)then
rhog  = 1.d0/rhomc
endif
endif
!
if(nlimi.eq.6)then
if(ndens.eq.1)then
rhomg = rhomc + aflim(1 ,ielem)*(rhomg - rhomc)
rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
!
elseif(ndens.eq.3)then
rhog = 1.d0/rhomc + aflim(1 ,ielem)*(rhog - 1.d0/rhomc)
endif
endif
!
!...Mass...
!
masel = masel + 1.d0/rhomc*djaco
!
!if(ielem.eq.1) print*,'bad',djaco,xpq(1, 1:npqua)
!
xcent = xcent + rhog*r*djaco !...mass center at the reference element...
ycent = ycent + rhog*s*djaco
!
xcel = xcel + 1.d0/rhomc*r*djaco
ycel = ycel + 1.d0/rhomc*s*djaco
!
!if(ie.eq.1) print*,'ieeq1',djaco,r,s
!
enddo
!
geoel(1, ielem) = xcent/masel
geoel(2, ielem) = ycent/masel
!
geoel(15, ielem) = xcel/masel
geoel(16, ielem) = ycel/masel
!
enddo !iloop
!
enddo !...(1)ie = 1,nelem
!
!print*,'geoel',geoel(1:8,1)
!
end subroutine getgeoel_lagmc_rzawrev
!
!....domain integral for hybrid linear triangle cells for R-Z (Area-weighted)for the total derivative of the basis function
!
subroutine rhsdomndg_lagmc_triarzawdt(intfac, iptri, coord,coold, geoel, unkno, ustar, rhsel,aflim, afvec )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),     intent(in)::ustar !...nodal velocity
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
real*8,dimension(1:ndimn, 1:nvtri) :: xp, xpi, xpn
real*8,dimension(1:ndegr):: b, dbdr, dbds,bv
real*8,dimension(1:ndegr):: bdt, dbdt
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
real*8:: drhomdr, drhomds, dedr, deds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rcdt,scdt,rc2,sc2
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
data eps   / 1.0d-6 /
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
!...New points in dt
!
xpn(1, 1:nvtri) = xp(1, 1:nvtri) + eps*dtfix*ustar(1, ipt(1:nvtri))
xpn(2, 1:nvtri) = xp(2, 1:nvtri) + eps*dtfix*ustar(2, ipt(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
!
rc= geoel(15, ielem) !...mass center with the first-order density distribution...
sc= geoel(16, ielem)
!
rc2= geoel(1, ielem) !...mass center with the second-order density distribution...
sc2= geoel(2, ielem)
!
rcdt= geoel(13, ielem) !...mass center after dt...
scdt= geoel(14, ielem)
!
!...Get the total time derivative of the basis funciton: dphi/dt
!
if(ncase.ne.4)then
dbdt(1) = 0.d0
dbdt(2) = (rc-rcdt)/dr/(epsaw*dtfix)
dbdt(3) = (sc-scdt)/ds/(epsaw*dtfix)
else
dbdt(1) = 0.d0
dbdt(2) = (rc-rcdt)/dr/(epsaw*dtfix*6.726842539779d-3)
dbdt(3) = (sc-scdt)/ds/(epsaw*dtfix*6.726842539779d-3)
endif

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
djaco = 0.5d0*wi*(dxdr*dyds - dxds*dydr)
!
!...Gauss points...
!
xg = r
yg = s
!
!...Basis function for solutions...
!
b(1) = 1.d0
b(2) = (xg-rc2)/dr
b(3) = (yg-sc2)/ds
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
drhomdr  = unkno(2, 1, ielem)*aflim(1,ielem)
drhomds  = unkno(3, 1, ielem)*aflim(1,ielem)
!
!dudr = unkno(2,2,ielem)
!duds = unkno(3,2,ielem)
!dvdr = unkno(2,3,ielem)
!dvds = unkno(3,3,ielem)
!
dedr = unkno(2, 4, ielem)*aflim(5,ielem)
deds = unkno(3, 4, ielem)*aflim(5,ielem)
!
endif
!
fluxd(1,1) = rhoad*b(1)*drhomdr*dbdt(2) + rhoad*b(1)*drhomds*dbdt(3)
fluxd(2,1) = rhoad*b(2)*drhomdr*dbdt(2) + rhoad*b(2)*drhomds*dbdt(3)
fluxd(3,1) = rhoad*b(3)*drhomdr*dbdt(2) + rhoad*b(3)*drhomds*dbdt(3)
!
fluxd(1,2) = rhoad*b(1)*dudr*dbdt(2) + rhoad*b(1)*duds*dbdt(3)
fluxd(2,2) = rhoad*b(2)*dudr*dbdt(2) + rhoad*b(2)*duds*dbdt(3)
fluxd(3,2) = rhoad*b(3)*dudr*dbdt(2) + rhoad*b(3)*duds*dbdt(3)
!
fluxd(1,3) = rhoad*b(1)*dvdr*dbdt(2) + rhoad*b(1)*dvds*dbdt(3)
fluxd(2,3) = rhoad*b(2)*dvdr*dbdt(2) + rhoad*b(2)*dvds*dbdt(3)
fluxd(3,3) = rhoad*b(3)*dvdr*dbdt(2) + rhoad*b(3)*dvds*dbdt(3)
!
fluxd(1,4) = rhoad*b(1)*dedr*dbdt(2) + rhoad*b(1)*deds*dbdt(3)
fluxd(2,4) = rhoad*b(2)*dedr*dbdt(2) + rhoad*b(2)*deds*dbdt(3)
fluxd(3,4) = rhoad*b(3)*dedr*dbdt(2) + rhoad*b(3)*deds*dbdt(3)
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
end subroutine rhsdomndg_lagmc_triarzawdt
!
!....domain integral for hybrid linear quad cells for R-Z (Area-weighted) for the total derivative of the basis function
!
subroutine rhsdomndg_lagmc_quadrzawdt(intfac, ipqua, coord, coold, geoel, unkno, ustar, rhsel,aflim,afvec )
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),     intent(in)::ustar !...nodal velocity
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
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xpqn
real*8,dimension(1:ndegr):: b, dbdr, dbds, bv
real*8,dimension(1:ndegr):: bdt,dbdt
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
real*8:: drhomdr, drhomds, dedr, deds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rcdt,scdt,rc2,sc2
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
data eps   / 1.d-6/
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
call ruqope(2, ngausdq, posiq, weighq)
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
!...New points in dt
!
xpqn(1, 1:4) = xpq(1, 1:4) + eps*dtfix*ustar(1, ipq(1:nvqua))
xpqn(2, 1:4) = xpq(2, 1:4) + eps*dtfix*ustar(2, ipq(1:nvqua))
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(15, ielem) !...mass center with the first-order density distribution...
sc= geoel(16, ielem)
!
rc2= geoel(1, ielem) !...mass center with the second-order density distribution...
sc2= geoel(2, ielem)
!
rcdt= geoel(13, ielem) !...mass center after dt...
scdt= geoel(14, ielem)
!
!...Get the total time derivative of the basis funciton: dphi/dt
!
if(ncase.ne.4)then
dbdt(1) = 0.d0
dbdt(2) = (rc-rcdt)/dr/(epsaw*dtfix)
dbdt(3) = (sc-scdt)/ds/(epsaw*dtfix)
else
dbdt(1) = 0.d0
dbdt(2) = (rc-rcdt)/dr/(epsaw*dtfix*6.726842539779d-3)
dbdt(3) = (sc-scdt)/ds/(epsaw*dtfix*6.726842539779d-3)
endif
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
djaco = wi*(dxdr*dyds - dxds*dydr)
!
!...Gauss points...
!
xg = r
yg = s
!
!...Basis function for solutions...
!
b(1) = 1.d0
b(2) = (xg-rc2)/dr
b(3) = (yg-sc2)/ds
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
!
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
drhomdr  = unkno(2, 1, ielem)*aflim(1,ielem)
drhomds  = unkno(3, 1, ielem)*aflim(1,ielem)
!
!dudr = unkno(2,2,ielem)
!duds = unkno(3,2,ielem)
!dvdr = unkno(2,3,ielem)
!dvds = unkno(3,3,ielem)
!
dedr = unkno(2, 4, ielem)*aflim(5,ielem)
deds = unkno(3, 4, ielem)*aflim(5,ielem)

!
!if(ielem.eq.9) print*,'domn',pres,aflim(4,ielem),xg,yg
!
endif
!
fluxd(1,1) = rhoad*b(1)*drhomdr*dbdt(2) + rhoad*b(1)*drhomds*dbdt(3)
fluxd(2,1) = rhoad*b(2)*drhomdr*dbdt(2) + rhoad*b(2)*drhomds*dbdt(3)
fluxd(3,1) = rhoad*b(3)*drhomdr*dbdt(2) + rhoad*b(3)*drhomds*dbdt(3)
!
fluxd(1,2) = rhoad*b(1)*dudr*dbdt(2) + rhoad*b(1)*duds*dbdt(3)
fluxd(2,2) = rhoad*b(2)*dudr*dbdt(2) + rhoad*b(2)*duds*dbdt(3)
fluxd(3,2) = rhoad*b(3)*dudr*dbdt(2) + rhoad*b(3)*duds*dbdt(3)
!
fluxd(1,3) = rhoad*b(1)*dvdr*dbdt(2) + rhoad*b(1)*dvds*dbdt(3)
fluxd(2,3) = rhoad*b(2)*dvdr*dbdt(2) + rhoad*b(2)*dvds*dbdt(3)
fluxd(3,3) = rhoad*b(3)*dvdr*dbdt(2) + rhoad*b(3)*dvds*dbdt(3)
!
fluxd(1,4) = rhoad*b(1)*dedr*dbdt(2) + rhoad*b(1)*deds*dbdt(3)
fluxd(2,4) = rhoad*b(2)*dedr*dbdt(2) + rhoad*b(2)*deds*dbdt(3)
fluxd(3,4) = rhoad*b(3)*dedr*dbdt(2) + rhoad*b(3)*deds*dbdt(3)
!
!fluxd(:,2:4) = 0.d0
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) - fluxd(ideg,1:nq)*djaco
enddo
!
!if(ielem.eq.4000) print*,'domn',ielem,rhsel(1,1,ielem),fluxd(1,1)*djaco,rcdt,rc,(eps*dtfix)
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
!
!print*,'ieleme.q.2',rhsel(1:3,3,2)
!
end subroutine rhsdomndg_lagmc_quadrzawdt
!
!...Find mass center in geoel for curved cell...
!
subroutine getgeoel_lagmc_rzawdt(inpoel, iptri, ipqua, geoel, coord, unkno, ustar, aflim)
use constant
implicit none
integer*4,dimension(1:nvtri, 1:ntria),intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:ngeel, 1:nsize), intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),     intent(in)::ustar !...nodal velocity
real*8,dimension(1:nq+1, 1:nsize),intent(in)::aflim
!
!...local array...
!
real*8,dimension(1:2, 1:3)::xpin
real*8,dimension(1:2, 1:nptri)::xp, xpn
real*8,dimension(1:2, 1:npqua)::xpq, xpqn
real*8,dimension(1:nptri)::shp, dspr, dsps
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weight(ngausd_geo), posit(2, ngausd_geo)
real*8:: weighq(ngausd_geoq), posiq(2, ngausd_geoq)
real*8,dimension(1:3)::bt
real*8,dimension(1:3)::bq

!...local real number
real*8:: dxdr,dxds,dydr,dyds
real*8:: xc, yc, xcel, ycel
real*8:: eps,c00,c10,c05,c20
real*8:: r, s, djaco, volel, masel
real*8:: wi, xcent, ycent,xgaus, ygaus
real*8:: rhomc, rhomg, rhog
real*8::rcoef !...Coefficient R of RZ or XY coordinate system...
real*8::areel !...2D area of the cell
!
real*8:: dr, ds, rc ,sc
!
integer::ielem, igaus, ishp,ideg,iloop
integer::ie
!
data eps / 1.0d-02 /
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
!
!...New points in dt
!
xpn(1, 1:3) = xp(1, 1:3) + eps*dtfix*ustar(1, iptri(1:3, ie))
xpn(2, 1:3) = xp(2, 1:3) + eps*dtfix*ustar(2, iptri(1:3, ie))
!
xpn(1:2,4) = 0.5d0*(xpn(1:2,1)+xpn(1:2,2))
xpn(1:2,5) = 0.5d0*(xpn(1:2,2)+xpn(1:2,3))
xpn(1:2,6) = 0.5d0*(xpn(1:2,1)+xpn(1:2,3))
!...
xcent = 0.d0
ycent = 0.d0
xcel  = 0.d0
ycel  = 0.d0
volel = 0.d0
masel = 0.d0
areel = 0.d0
!
!...Mass averaged center for Area-weighted RZ
!
rc = geoel(13, ielem)
sc = geoel(14, ielem)
!
dr = 0.5d0
ds = 0.5d0
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
!...Restore r, s
!
r  = posit(1,igaus)
s  = posit(2,igaus)
!
bt(1) = 1.d0
bt(2) = (r - rc)/dr
bt(3) = (s - sc)/ds
!
!...Density calculations with different density evolutions
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
rhomg  = rhomc
rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
rhog = 1.d0/rhomc
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
rhog = 1.d0/rhomc
endif
!
!...Mass...
!
masel = masel + 1.d0/rhomc*djaco
!
xcent = xcent + rhog*r*djaco!...mass center at the reference element...
ycent = ycent + rhog*s*djaco
!
enddo
!
geoel(13, ielem) = xcent/masel
geoel(14, ielem) = ycent/masel
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
!...New points in dt
!
xpqn(1, 1:4) = xpq(1, 1:4) + eps*dtfix*ustar(1, ipqua(1:4, ie))
xpqn(2, 1:4) = xpq(2, 1:4) + eps*dtfix*ustar(2, ipqua(1:4, ie))
!
xpqn(1:2,5) = 0.5d0*(xpqn(1:2,1)+xpqn(1:2,2))
xpqn(1:2,6) = 0.5d0*(xpqn(1:2,2)+xpqn(1:2,3))
xpqn(1:2,7) = 0.5d0*(xpqn(1:2,3)+xpqn(1:2,4))
xpqn(1:2,8) = 0.5d0*(xpqn(1:2,4)+xpqn(1:2,1))
xpqn(1:2,9) = 0.5d0*(xpqn(1:2,5)+xpqn(1:2,7))
!
!...Mass averaged center for Area-weighted RZ
!
do iloop = 1,1
!
xcent = 0.d0
ycent = 0.d0
xcel = 0.d0
ycel = 0.d0
volel = 0.d0
masel = 0.d0
areel = 0.d0
!
!if(iloop.eq.1)then
!rc= geoel(1, ielem) !...mass center...
!sc= geoel(2, ielem)
!else
rc= geoel(13, ielem) !...mass center...
sc= geoel(14, ielem)
!endif
!
dr = 1.d0
ds = 1.d0!
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
!...Restore r, s
r  = posiq(1,igaus)
s  = posiq(2,igaus)
!
bq(1) = 1.d0
bq(2) = (r - rc)/dr
bq(3) = (s - sc)/ds
!
!...Density calculations with different density evolutions
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
rhomg  = rhomc
rhog = 1.d0/rhomg
elseif(ndens.eq.2)then
rhomc = 1.d0/unkno(1, 1, ielem)
rhog = 1.d0/rhomc
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
rhog = 1.d0/rhomc
endif
!
!...Mass...
!
masel = masel + 1.d0/rhomc*djaco
!
xcent = xcent + rhog*r*djaco !...mass center at the reference element...
ycent = ycent + rhog*s*djaco
!
enddo
!
geoel(13, ielem) = xcent/masel
geoel(14, ielem) = ycent/masel
!
enddo!iloop
!
enddo !...(1)ie = 1,nelem
!
!print*,'geoel',geoel(1:8,1)
!
end subroutine getgeoel_lagmc_rzawdt
!
!...Find geometry parameters in geoel (17 ,18) for the 3rd density evolution at t+dt...
!
subroutine getgeo_denevoldt(iptri, ipqua, geoel, coord)
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
do ishp = 1, nptri
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
do ishp = 1, nptri
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
if(nrz.eq.2)then
geoel(17, ielem) = rcvrz/areel
geoel(18, ielem) = scvrz/areel
endif
!
!    print*,'old geoel',ielem, geoel(1:5, ielem)
!    print*,'new geoel',ielem, xcent/volel, ycent/volel, volel, maxval(xp(1, 1:3))-minval(xp(1, 1:3)), &
!                         maxval(xp(2, 1:3))- minval(xp(2, 1:3))
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
!...Zero out some variables...
!
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
if(nrz.eq.2)then
geoel(17, ielem) = rcvrz/areel
geoel(18, ielem) = scvrz/areel
endif
!
!    print*,'new geoel',ielem, xcent/volel, ycent/volel, volel, maxval(xp(1, 1:3))-minval(xp(1, 1:3)), &
!                         maxval(xp(2, 1:3))- minval(xp(2, 1:3))
!
enddo !...(1)ie = 1,nelem

end subroutine getgeo_denevoldt
!
!...Get the mass matrix based on volume center(density evolution) at t+dt...
!
subroutine  getamatr_denevoldt(amatr,geoel,coord, iptri, ipqua)
use constant
implicit none
!...Input
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
real*8:: dxdr,dxds,dydr,dyds
real*8::rhom, rho0
real*8::wi,djaco,det
real*8::c10
real*8::f0,f1,f2,f3,f4
real*8::f5,f6,f7,f8
real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::b2,b3,b4,b5,b6
real*8::radie, radii,radie2,radii2,radic2, radig2,paras,rhoi,rhoe
real*8::masel,xgaus,ygaus
real*8::rcoef
real*8::volel,areel
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
!...volum center...
!
rc= geoel(17, ielem) !...mass center...
sc= geoel(18, ielem)
!
!
dr = 0.5d0
ds = 0.5d0!
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
!
volel = 0.d0
areel = 0.d0
!
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
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
volel = volel + djaco*rcoef !The cell true volume
areel = areel + djaco !The cell area
!
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
if(nrz.eq.0.or.nrz.eq.1)then
f0 = f0 + djaco*rcoef
f1 = f1 + (xg-rc)/dr*(xg-rc)/dr*djaco*rcoef
f2 = f2 + (xg-rc)/dr*(yg-sc)/ds*djaco*rcoef
f3 = f3 + (yg-sc)/ds*(yg-sc)/ds*djaco*rcoef
!
elseif(nrz.eq.2)then
f0 = f0 + djaco
f1 = f1 + (xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + (xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + (yg-sc)/ds*(yg-sc)/ds*djaco
endif
!

if(npoly==2)then
!
if(nrz.eq.0.or.nrz.eq.1)then
f22 = f22 + b2*b2*djaco*rcoef
f23 = f23 + b2*b3*djaco*rcoef
f24 = f24 + b2*b4*djaco*rcoef
f25 = f25 + b2*b5*djaco*rcoef
f26 = f26 + b2*b6*djaco*rcoef

f33 = f33 + b3*b3*djaco*rcoef
f34 = f34 + b3*b4*djaco*rcoef
f35 = f35 + b3*b5*djaco*rcoef
f36 = f36 + b3*b6*djaco*rcoef

f44 = f44 + b4*b4*djaco*rcoef
f45 = f45 + b4*b5*djaco*rcoef
f46 = f46 + b4*b6*djaco*rcoef

f55 = f55 + b5*b5*djaco*rcoef
f56 = f56 + b5*b6*djaco*rcoef

f66 = f66 + b6*b6*djaco*rcoef
!
elseif(nrz.eq.2)then
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

endif
!   if(ie==400) print*,'gauss',xg, yg, f0
enddo
!
!print*,'ielem', ie, f0, volel

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
!...
!   if(ie==2) print*,'ie',mmatr(:,:)
x5 = 0.d0
b55 = 0.d0
call getinvmat(5, mmatr, x5, b55)
!
!   if(ie==2) print*,'ie',x5(:,:)
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
if(nrz.eq.2) amatr(:, ielem) = amatr(:, ielem)/(volel/areel)!geoel(11, ielem)
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
rc= geoel(17, ielem) !...mass center...
sc= geoel(18, ielem)
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
!
volel = 0.d0
areel = 0.d0
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
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
!
volel = volel + djaco*rcoef !The cell true volume
areel = areel + djaco !The cell area
!
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
if(nrz.eq.0.or.nrz.eq.1)then
f0 = f0 + djaco*rcoef
f1 = f1 + (xg-rc)/dr*(xg-rc)/dr*djaco*rcoef
f2 = f2 + (xg-rc)/dr*(yg-sc)/ds*djaco*rcoef
f3 = f3 + (yg-sc)/ds*(yg-sc)/ds*djaco*rcoef
elseif(nrz.eq.2)then
f0 = f0 + djaco
f1 = f1 + (xg-rc)/dr*(xg-rc)/dr*djaco
f2 = f2 + (xg-rc)/dr*(yg-sc)/ds*djaco
f3 = f3 + (yg-sc)/ds*(yg-sc)/ds*djaco
endif

if(npoly==2)then
if(nrz.eq.0.or.nrz.eq.1)then
f22 = f22 + b2*b2*djaco*rcoef
f23 = f23 + b2*b3*djaco*rcoef
f24 = f24 + b2*b4*djaco*rcoef
f25 = f25 + b2*b5*djaco*rcoef
f26 = f26 + b2*b6*djaco*rcoef

f33 = f33 + b3*b3*djaco*rcoef
f34 = f34 + b3*b4*djaco*rcoef
f35 = f35 + b3*b5*djaco*rcoef
f36 = f36 + b3*b6*djaco*rcoef

f44 = f44 + b4*b4*djaco*rcoef
f45 = f45 + b4*b5*djaco*rcoef
f46 = f46 + b4*b6*djaco*rcoef

f55 = f55 + b5*b5*djaco*rcoef
f56 = f56 + b5*b6*djaco*rcoef

f66 = f66 + rho0*b6*b6*djaco*rcoef
elseif(nrz.eq.2)then
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

endif
!   if(ie==400) print*,'gauss',xg, yg, f0
enddo
!
!print*,'ielem', ie, f0, volel

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
!...
!   if(ie==2) print*,'ie',mmatr(:,:)
x5 = 0.d0
b55 = 0.d0
call getinvmat(5, mmatr, x5, b55)
!
!   if(ie==2) print*,'ie',x5(:,:)
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
if(nrz.eq.2) amatr(:, ielem) = amatr(:, ielem)/(volel/areel)!geoel(11, ielem)
!
enddo !...(2)ie = 1,nelem

end subroutine  getamatr_denevoldt
!
!....rhs for density evoltuion at t+dt
!
subroutine rhs_denevoldt(iptri, ipqua, coord, coold, geoel, unkno, amatr)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
real*8,dimension(1:nmatr,1:ncell),     intent(in)::amatr
real*8,dimension(1:ndimn,1:npoin),     intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!
real*8,dimension(1:ndegr,1:ncell)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem,iunk,id
!
!...local integer array
!
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8::m(ndegr, ndegr)
real*8,dimension(1:ndimn, 1:nvtri) :: xpti
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndegr):: b, dbdr, dbds
!
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8::weigh(ngausd), posi(2,ngausd)

real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8::unint(1)

!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::rhoad,uadv,vadv,eadv,rhoma
real*8::pres
real*8::djaco, wi
real*8:: rhoi
real*8::xgaus, ygaus
real*8::rcoef
real*8::areel,volel
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
!...Quads
!
!  open(12,file='solution2.dat')
!    do ie=1,ncell
!       if(ie.ge.990.and.ie.le.1000)print*,'afterdomn',ie,unkno(1:ndegr,1,ie)
!    enddo
!  close(12)
!
rhsel = 0.d0
!
!...Loop over elements-triangles
!
do 550 ie = 1,ntria!...(1)ie = 1,nelem
!
ielem = ie
!
!...Points consitituting one element...
!
ipt(1:nvtri) = iptri(1:nvtri, ie)
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
!
rc= geoel(7, ielem) !...Initial cell center(volume)...
sc= geoel(8, ielem)

!...Zero out sth

volel = 0.d0
areel = 0.d0
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
dxdr = dxdr + dspr(ishp)*xpti(1,ishp)
dxds = dxds + dsps(ishp)*xpti(1,ishp)

dydr = dydr + dspr(ishp)*xpti(2,ishp)
dyds = dyds + dsps(ishp)*xpti(2,ishp)
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
xgaus = xgaus + shp(ishp)*xpti(1,ishp)
ygaus = ygaus + shp(ishp)*xpti(2,ishp)
enddo
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
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
volel = volel + djaco*rcoef
areel = areel + djaco
!
call  getrhoig_tria(rhoi, r, s, xpti)
!
!finally, scatter the contribution to the RHS
!
if(nrz.eq.0.or.nrz.eq.1)then
do ideg = 1,ndegr
rhsel(ideg, ielem)=rhsel(ideg, ielem) + rhoi*b(ideg)*djaco*rcoef
enddo
elseif(nrz.eq.2)then
do ideg = 1,ndegr
rhsel(ideg, ielem)=rhsel(ideg, ielem) + rhoi*b(ideg)*djaco
enddo
endif
!
enddo !...(2)ig = 1,ngausd
!
if(nrz.eq.2) rhsel(:, ielem) = rhsel(:, ielem)*volel/areel
!
550 enddo
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
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(7, ielem) !...mass center...
sc= geoel(8, ielem)

!...Zero out sth

volel = 0.d0
areel = 0.d0
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
dxdr = dxdr + dsprq(ishp)*xpqi(1,ishp)
dxds = dxds + dspsq(ishp)*xpqi(1,ishp)

dydr = dydr + dsprq(ishp)*xpqi(2,ishp)
dyds = dyds + dspsq(ishp)*xpqi(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
!...Density distribution for different cases...
!
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, 4
xgaus = xgaus + shpq(ishp)*xpqi(1,ishp)
ygaus = ygaus + shpq(ishp)*xpqi(2,ishp)
enddo
!
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
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
volel = volel + djaco*rcoef
areel = areel + djaco
!
!...Solution at the Gauss points...
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
!
!finally, scatter the contribution to the RHS
!
if(nrz.eq.0.or.nrz.eq.1)then
do ideg = 1,ndegr
rhsel(ideg, ielem)=rhsel(ideg, ielem) + rhoi*b(ideg)*djaco*rcoef
enddo
elseif(nrz.eq.2)then
do ideg = 1,ndegr
rhsel(ideg, ielem)=rhsel(ideg, ielem) + rhoi*b(ideg)*djaco
enddo
endif
!
enddo !...(2)ig = 1,ngausd
!
if(nrz.eq.2) rhsel(:, ielem) = rhsel(:, ielem)*volel/areel
!
650 enddo

!
!...Update density...
!

if(npoly.ge.1)then !...(2)npoly.ge.1
!
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

elseif(npoly==2)then

m(1,1) = amatr(16, ie)
m(1,2) = 0.d0
m(1,3) = 0.d0
m(1,4) = 0.d0
m(1,5) = 0.d0
m(1,6) = 0.d0

m(2,1) = m(1,2)
m(2,2) = amatr(1,ie)
m(2,3) = amatr(2,ie)
m(2,4) = amatr(3,ie)
m(2,5) = amatr(4,ie)
m(2,6) = amatr(5,ie)

m(3,1) = m(1,3)
m(3,2) = m(2,3)
m(3,3) = amatr(6,ie)
m(3,4) = amatr(7,ie)
m(3,5) = amatr(8,ie)
m(3,6) = amatr(9,ie)

m(4,1) = m(1,4)
m(4,2) = m(2,4)
m(4,3) = m(3,4)
m(4,4) = amatr(10,ie)
m(4,5) = amatr(11,ie)
m(4,6) = amatr(12,ie)

m(5,1) = m(1,5)
m(5,2) = m(2,5)
m(5,3) = m(3,5)
m(5,4) = m(4,5)
m(5,5) = amatr(13,ie)
m(5,6) = amatr(14,ie)

m(6,1) = m(1,6)
m(6,2) = m(2,6)
m(6,3) = m(3,6)
m(6,4) = m(4,6)
m(6,5) = m(5,6)
m(6,6) = amatr(15,ie)
endif
!
!...solve the mdegr independant varaible...
!
do id =1,ndegr !...(4)id =1,ndegr
!...step 1
unint = 0.d0
do iunk = 1,ndegr
unint(1) = unint(1) + m(id, iunk)*rhsel(iunk,ie)
enddo
!...step 2
unkno(id, 1, ie)= unint(1)
enddo !...(4)id =1,ndegr
enddo !...(3)ie=1,nelem
!
endif !...(2)npoly.ge.1

end subroutine rhs_denevoldt
!
!...subroutine: Get the RHS from the total derivative of the basis function for AW RZ
!
subroutine rhsdomndg_lagmc_rzawdt(ustar,unkno,rhsel,rhsaw,intfac,inpoel,iptri,ipqua,geoel,coord,coold,amatr, aflim, afvec)
use constant
implicit none

!...input arrays
integer*4,dimension(1:nvtri,1:ntria), intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer*4,dimension(1:nifai,1:nafac), intent(in)::intfac
real*8,dimension(1:ngeel,1:nsize), intent(inout)::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndegr,1:nq,1:ncell),intent(in)::rhsel
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coold!...initial coordinates...
real*8,dimension(1:ndimn,1:npoin),intent(in)::coord
real*8,dimension(nmatr,ncell),    intent(in)::amatr
real*8,dimension(1:ndimn,1:npoin),   intent(in)::ustar !...nodal velocity
real*8, dimension(1:nq+1, 1:nsize),  intent(in)::aflim
real*8, dimension(2, 2, nsize),      intent(in)::afvec
real*8,dimension(1:ndegr,1:nq,1:ncell),intent(inout)::rhsaw

!...Local inetger
integer:: iloop, ie, iq, id ,ipoin, iunk

!...Local real arrays
real*8::m(3,3),unint(nq)
!
real*8, allocatable:: coora(:,:)
real*8, allocatable:: amatrv(:,:)
real*8, allocatable:: rhsor(:,:, :)
real*8, allocatable:: unkaw(:,:, :),unori(:,:, :)
!
allocate(coora(1:ndimn,1:npoin))
allocate(amatrv(nmatr,ncell))
allocate(unkaw(1:ndegr,1:nq,1:nsize), unori(1:ndegr,1:nq,1:nsize))!
allocate(rhsor(1:ndegr,1:nq,1:ncell))

!...Store the RHS without total derivative of the basis function...
rhsor = rhsel
unori = unkno

!...Iterative loop
do iloop =1, 1

if(ndens.eq.1)then
 if(ncase.ne.4)then
 do ie =1, ncell
 do iq =1, nq  !...nq variables
 rhsor(1:ndegr,iq,ie)=dtfix*rhsor(1:ndegr,iq,ie)!*7.265d-3
 enddo
 enddo
elseif(ncase.eq.4)then
 do ie =1, ncell
 do iq =1, nq  !...nq variables
 rhsor(1:ndegr,iq,ie)=dtfix*rhsor(1:ndegr,iq,ie)*tfcus
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
unint(1:nq) = unint(1:nq) + m(id, iunk)*rhsor(iunk,1:nq,ie)
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
rhsaw = 0.d0
if(ntria.gt.0) call rhsdomndg_lagmc_triarzawdt(intfac, iptri, coord, coold, geoel, unori, ustar, rhsaw,aflim,afvec )
if(nquad.gt.0) call rhsdomndg_lagmc_quadrzawdt(intfac, ipqua, coord, coold, geoel, unori, ustar, rhsaw,aflim,afvec )

enddo !iloop

deallocate(amatrv, coora, unkaw, unori, rhsor)
!
end subroutine rhsdomndg_lagmc_rzawdt

