!
!...Riemann input for curved quad using simpson rule using the Gradient deformation...
!
subroutine getriem_quad_simpson_gd(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec, vnulq)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq
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
real*8,dimension(1:ndegr)::bg, bgv
real*8,dimension(1:nq,1:nvfac)::unkng
real*8::aujmp(1:3, 1:nvfac)
real*8::vnorm(1:3, 1:12)
real*8::sigmg(1:2, 1:2, 1:nvfac)
real*8,dimension(1:nvfac)::murie
real*8,dimension(1:nvqua):: xvq,  yvq
real*8::lnvp(2,nvqua)
!...arraies for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
!
real*8::eps,c00,c05,c10,c20
real*8::rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::rhog,rhomg,ug,vg,eg, pg
real*8::dux,duy,deltu,gpnx,gpny
real*8::dr, ds, rc, sc, rg, sg,rcv,scv
real*8:: dudr, duds, dvdr, dvds
real*8::acnx, acny,divu
real*8::sdimp
!
data eps   / 1.0d-14/
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Part I: Specify some Gauss points
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

!...9-node quadrilateral
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
!...Part II: Loop over every quad...
!
do 250 ie = 1,nquad !...(1)ie = 1,nquad
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria

!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...Give the normal vector of every face...
vnorm(1:3,  1:12) = gesgq(1:3, 1:12, ie);

!...Get the LN for one vertex used for judging cell compressing or expanding...
lnvp = 0.d0

lnvp(1:2, 1) = vnorm(1:2, 1)*vnorm(3, 1) + vnorm(1:2, 8)*vnorm(3, 8)
lnvp(1:2, 2) = vnorm(1:2, 2)*vnorm(3, 2) + vnorm(1:2, 3)*vnorm(3, 3)
lnvp(1:2, 3) = vnorm(1:2, 4)*vnorm(3, 4) + vnorm(1:2, 5)*vnorm(3, 5)
lnvp(1:2, 4) = vnorm(1:2, 6)*vnorm(3, 6) + vnorm(1:2, 7)*vnorm(3, 7)
lnvp(1:2, 5) = vnorm(1:2, 9)*vnorm(3, 9)
lnvp(1:2, 6) = vnorm(1:2, 10)*vnorm(3, 10)
lnvp(1:2, 7) = vnorm(1:2, 11)*vnorm(3, 11)
lnvp(1:2, 8) = vnorm(1:2, 12)*vnorm(3, 12)

!...cell averaged value...
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

!...Part II.1: Loop over every face...

do ifa =1, 4

!...zero out unkno
unkng = 0.d0

!...Loop over every gauss point
do ig   = 1, nvfac
rg = xvq(fglvq(ig, ifa))
sg = yvq(fglvq(ig, ifa))
!
bg(1) = 1.d0
bg(2) = (rg-rc)/dr
bg(3) = (sg-sc)/ds

!...DGP2
if(npoly.eq.2)then
bg(4) = 0.5d0*bg(2)*bg(2) - geoel(19, ielem)
bg(5) = 0.5d0*bg(3)*bg(3) - geoel(20, ielem)
bg(6) =       bg(2)*bg(3) - geoel(21, ielem)
endif
!
do ideg = 1,mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq,ielem)*bg(ideg)
enddo

!...Different density
if(ndens.eq.1)then
rhog  = 1.d0/unkng(1, ig)
elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)

bgv(1) = 1.d0
bgv(2) = (rg-rcv)/dr
bgv(3) = (sg-scv)/ds

!...DGP2
!if(npoly.eq.2)then
!bgv(4) = 0.5d0*bgv(2)*bgv(2) - geoel(19, ielem)
!bgv(5) = 0.5d0*bgv(3)*bgv(3) - geoel(20, ielem)
!bgv(6) =       bgv(2)*bgv(3) - geoel(21, ielem)
!endif
unkng(1, ig) =0.d0
do ideg = 1,mdegr
unkng(1, ig) = unkng(1, ig) + unkno(ideg,1,ielem)*bgv(ideg)
enddo
rhog = unkng(1, ig)
endif

!...Velocity, total enrgy, pressure at some Gauss point
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
pg = max(eps, (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2)))

!...Limiter
if(nlimi.eq.6.and.geoel(10,ielem).gt.10.d0)then
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
pg = pctr + aflim(4, ielem)*(pg - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig) = ug! -  vnulq(1,fglvq(ig, ifa),ie)
unkng(3 ,ig) = vg! -  vnulq(2,fglvq(ig, ifa),ie)
endif

!...Get stress tensor at nodes
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

enddo !...end of Loop over every gauss point

!...Get sound speed at the center...
sdctr = sqrt( gamlg*pctr/rhoct)
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedence coefficient (Loop over every gauss point)...
do ig   = 1, nvfac

!...Original impedence
!dux= vlave(1, ipq(fglvq(ig, ifa)))-unkng(2, ig)
!duy= vlave(2, ipq(fglvq(ig, ifa)))-unkng(3, ig)
!deltu = 1.d0*sqrt(dux**2 + duy**2)
!
gpnx = vnorm(1, fglgq(ig, ifa))
gpny = vnorm(2, fglgq(ig, ifa))
dux = vlave(1, ipq(fglvq(ig, ifa)))-unkno(1, 2, ielem)
duy = vlave(2, ipq(fglvq(ig, ifa)))-unkno(1, 3, ielem)
deltu = abs(dux*gpnx + duy*gpny)
!

!...Treatmen for compression or expansion
!divu = dux*lnvp(1, fglvq(ig, ifa)) + duy*lnvp(2, fglvq(ig, ifa))
!if(divu.gt.0.d0)then
!deltu = 0.d0
!endif

deltu = abs(dux*gpnx + duy*gpny)
sdimp = sdctr!max(1d-8,sdctr)
murie(ig) = rhoct*sdimp + cimpd*rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!murie(ig) = rhoct*(slpdu*deltu*0.5d0 + sqrt((slpdu*deltu*0.5d0)**2+sdimp**2))
enddo

!...Feed the above into Riemann solver(Loop over every gauss point)
do ig  = 1, nvfac

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
snsigmlq(1:2,fglgq(ig, ifa), ie)   = snsigm_rie(1:2)

!...Output for debugging
! if(ipq(fglvq(ig, ifa)).eq.20) print*,'ep54 ',ie,ifa,ig,ipq(fglvq(ig, ifa)),&
!munacn_rie(1, 1),murie(ig),vnorm(1:3, fglgq(ig, ifa)),&
!sigmg(1, 1, ig),unkng(2:3, ig)
!
enddo
!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,ntria!

end subroutine getriem_quad_simpson_gd
!
!...Get the nodal velocity based on simpson integral rule using the Gradient deformation...
!
subroutine getndvelo_lag_simpsonh_gd(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno,ustar, fstart, fstarq, aflim, afvec,vnulq, itime)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac),intent(in)::gesgt
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
!...Part I: Specify some Gauss points
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

!...Zero out vlave
vlave = 0.d0
indnd = 0

!...Mark the boundary nodes...
if(ncase.eq.2)then !...2D shockless Noh
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
!call getnullve_gqdmpquadc(ipqua, geoel, ustar, unkno, gqdmp)
!
vnulq = 0.d0
!
!...Part II: Loop to get the info from Riemann solver
!
do iloop= 1,1!7!!!

!...Give vlave
vlave= ustar

!call getnullmd_quadlp(ipqua, geoel, vnulq, ustar, coord, gqdmp)

!...Zero out munacn
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0

!...Tria
!if(ntria.gt.0) call getriem_tria_simpson(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
!munaclt, munault, snsigmlt, coord, coold,aflim, afvec)

!...Quad
if(nquad.gt.0) call getriem_quad_simpson_gd(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec, vnulq)

!...Boundary condition...
!call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)
!call getboundary_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm)

!...Periodic boundary condition for 1D isentropic Sin problem...
!if(ncase.eq.12) call getbc_prdic(bface, munacn, munacu, snsigm)

!...Update the velocity at the vertex...
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
!
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
!
!if(ipoin.eq.20) print*,'ep',ustar(1:2,ipoin),detma,munacn(1,1,ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin),&
!munacu(1, ipoin)/munacn(1,1,ipoin),gesgq(1:3,5,1),gesgq(1:3,11,1)
endif
enddo

!...Velocity at the surface vertex
!if(ntria.gt.0) call getvelo_mpt_marie(ustar,gelag,intfac,inpoel,coord,unkno,indnd)
!if(nquad.gt.0) call getvelo_mpt_mariequad(ustar,gelagq,intfac,ipqua,coord,unkno,indnd)

!call getvelo_mpt_curv(ustar,geoel,gelag,gelagq,intfac,iptri,ipqua,coord,coold,unkno,indnd, aflim, afvec, vlave, vnulq)
!call getvelo_mpt(ustar,gelag,intfac,iptri,coord,munacn,vlave,unkno)
! call getvelo_mpt_curv2(ustar,geoel,gelag,gelagq,intfac,iptri,ipqua,coord,unkno,indnd, aflim, afvec, vlave, vnulq,munacn)
!
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
!
if(ncase.eq.13)then
ustar(2,ipf(1:nvfac)) = 0.d0
ustar(1,ipf(1:nvfac)) = 1.d0
endif
endif

!...Linearize the moving velocity
!ustar(1, ipf(3)) = 0.5d0*(ustar(1, ipf(1))+ustar(1, ipf(2)))
!ustar(2, ipf(3)) = 0.5d0*(ustar(2, ipf(1))+ustar(2, ipf(2)))
!
if(ncase.eq.-1)then
!...impose exact solution along the Boundary
ustar(1, ipf(1:nvfac)) = sin(pi*coord(1,ipf(1:nvfac)))*cos(pi*coord(2,ipf(1:nvfac)))
ustar(2, ipf(1:nvfac)) =-cos(pi*coord(1,ipf(1:nvfac)))*sin(pi*coord(2,ipf(1:nvfac)))
endif
!
!if(ipf(1).eq.100) print*,'ipf',ipf(1),ustar(1:2, ipf(1))
!
900 enddo
!endif
!
!call getvelofit_mpt3(intfac, ipqua, geoel, ustar, vnulq, coord)

!...Velocity filter

!...1) Get the dampening coefficient
!if(iloop.eq.1) call getcoef_nuvequadc(ipqua, geoel, ustar, unkno, gqdmp)

!...2) Get the null-mode velocity
!call getnullmd_quadlp(ipqua, geoel, vnulq, ustar, coord, gqdmp)
!
enddo !iloop

!...Imposing the zero normal velocity for BC...
! call getbcvn_lag(bface, intfac, gflag, ustar)
! call getbcve_exact(bface, intfac, gflag, ustar, coord, itime)
!
!...4.2: Update the Riemann forces at every node...
!...Tria
do ie = 1, ntria
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!!
do ifa = 1, 3
!
do ig =1, nvfac
!
!...Basis function!
!
fstart(1, fglgt(ig, ifa), ie) = snsigmlt(1, fglgt(ig, ifa), ie) +&
munaclt(1,1, fglgt(ig, ifa), ie)*ustar(1, ipt(fglvt(ig, ifa)))+&
munaclt(2,1, fglgt(ig, ifa), ie)*ustar(2, ipt(fglvt(ig, ifa)))-&
munault(1, fglgt(ig, ifa), ie)
!
fstart(2, fglgt(ig, ifa), ie) = snsigmlt(2, fglgt(ig, ifa), ie) +&
munaclt(2,2,fglgt(ig, ifa), ie)*ustar(2, ipt(fglvt(ig, ifa)))+&
munaclt(1,2,fglgt(ig, ifa), ie)*ustar(1, ipt(fglvt(ig, ifa)))-&
munault(2, fglgt(ig, ifa), ie)
!
!if(ie.eq.23.and.fglgt(ig, ifa).eq.1) print*,fstart(1:2,fglgt(ig, ifa),23),snsigmlt(1:2, fglgt(ig, ifa), ie),&
!munaclt(1:2,1:2, fglgt(ig, ifa), ie),munault(1:2, fglgt(ig, ifa), ie)
!
enddo
enddo
enddo
!
!...Quad
!
do ie = 1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!!
do ifa = 1, 4
do ig =1, nvfac
!
!...Basis function!
!
fstarq(1, fglgq(ig, ifa), ie) = snsigmlq(1, fglgq(ig, ifa), ie) +&
munaclq(1,1, fglgq(ig, ifa), ie)*ustar(1, ipq(fglvq(ig, ifa)))+&
munaclq(2,1, fglgq(ig, ifa), ie)*ustar(2, ipq(fglvq(ig, ifa)))-&
munaulq(1, fglgq(ig, ifa), ie)
!
fstarq(2, fglgq(ig, ifa), ie) = snsigmlq(2, fglgq(ig, ifa), ie) +&
munaclq(2,2,fglgq(ig, ifa), ie)*ustar(2, ipq(fglvq(ig, ifa)))+&
munaclq(1,2,fglgq(ig, ifa), ie)*ustar(1, ipq(fglvq(ig, ifa)))-&
munaulq(2, fglgq(ig, ifa), ie)
!
!if(ie.eq.1.and.ifa.eq.1)print*,ie,ig,fglgq(ig, ifa),snsigmlq(1, fglgq(ig, ifa), ie),munaclq(1,1, fglgq(ig, ifa), ie),&
!ustar(1, ipq(fglvq(ig, ifa))),&
!munaclq(2,1, fglgq(ig, ifa), ie),ustar(2, ipq(fglvq(ig, ifa))),&
!munaulq(1, fglgq(ig, ifa), ie),munaclq(1,1, fglgq(ig, ifa), ie)*ustar(1, ipq(fglvq(ig, ifa)))
!
enddo
enddo
enddo
!
deallocate (munacn, bpres, fpres)
deallocate (munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lag_simpsonh_gd
!
!...subroutine: Calculate the nodal velocity U_p^* (mass center) for quadratic meshes (MEM) using SMS..
!
subroutine getndvelo_lagsms_gd(gflag,gesgt,gesgq,geoel,bface,intfac,iptri,ipqua,&
coord, coold, unkno, unkgd,strnq_devtp,ustar, fstart, fstarq, aflim, afvec, itime)
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
if(ntria.gt.0) print*,'Triangular cells will be impelmented in future!'

!...Quad
if(nquad.gt.0) call getriem_quadsms_simpson_gd(ipqua, geoel, gesgq, vlave, unkno, unkgd, strnq_devtp, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)

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
rhsu1 = munacu(1, ipoin) - snsigm(1, ipoin) - fpres(1, ipoin)
rhsu2 = munacu(2, ipoin) - snsigm(2, ipoin) - fpres(2, ipoin)
!
ustar(1, ipoin) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ustar(2, ipoin) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2

!if(ipoin.eq.1.or.ipoin.eq.2421) print*,'ep',ustar(1:2,ipoin),detma,munacn(1,1,ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin)
endif
enddo

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
ustar(:,ipq(nvqua)) = 0.d0
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
end subroutine getndvelo_lagsms_gd
!
!...subroutine: Calculate the F^* N dsfor all faces for hybrid curved grids using the gadient of deformation...
!
subroutine getfnds_laghybrid_gd(geoel,gflag,gesgq,gesgt,intfac,iptri,ipqua,coord,unkgd2)
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
integer,dimension(3, 3):: igt
integer,dimension(3, 4):: igq
!...local real array
real*8,dimension(1:2, 1:2)        ::magd, cmagd !...cofactor matrix...
real*8,dimension(1:ndimn, 1:3)::xpf
real*8,dimension(1:ndimn, 1:nvtri)::xpt
real*8,dimension(1:ndimn, 1:nvqua)::xpq
real*8,dimension(1:ndegr)::bti, bqi
real*8::vnorm(1:2)
real*8::posit(2, ngesgt),posiq(2, ngesgq),posif(nvfac, 1)
real*8::dshpr(3)
real*8::dxdr,dxds,dydr,dyds,djaco
real*8,dimension(1:ndegr,1:4,1:nsize)::unkgd
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
igt(1, 1) = 1; igt(2, 1) = 2; igt(3, 1) = 7;
igt(1, 2) = 3; igt(2, 2) = 4; igt(3, 2) = 8;
igt(1, 3) = 5; igt(2, 3) = 6; igt(3, 3) = 9;
!
if(ncurv.eq.1)then
posit(1, 7) =  .5d0; posit(2, 7) = 0.d0;
posit(1, 8) =  .5d0; posit(2, 8) = .5d0;
posit(1, 9) =  0.d0; posit(2, 9) = .5d0;
!
ivft(3, 1) = 4; ivft(3, 2) = 5; ivft(3, 3) = 6;

endif

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
igq(1, 1) = 1; igq(2, 1) = 2; igq(3, 1) = 9;
igq(1, 2) = 3; igq(2, 2) = 4; igq(3, 2) =10;
igq(1, 3) = 5; igq(2, 3) = 6; igq(3, 3) =11;
igq(1, 4) = 7; igq(2, 4) = 8; igq(3, 4) =12;
!
!
if(ncurv.eq.1)then
posiq(1, 9) = 0.d0; posiq(2, 9) =-1.d0;
posiq(1,10) = 1.d0; posiq(2,10) = 0.d0;
posiq(1,11) = 0.d0; posiq(2,11) = 1.d0;
posiq(1,12) =-1.d0; posiq(2,12) = 0.d0;
!
ivfq(3, 1) = 5; ivfq(3, 2) = 6; ivfq(3, 3) = 7; ivfq(3, 4) = 8;
endif

!...Face position
posif(1, 1) = -1.d0
posif(2, 1) =  1.d0

if(ncurv.eq.1)then
posif(3, 1) =  0.d0
endif

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
xpf(1, 1:nvfac) = coord(1, ipt(ivft(1:nvfac, ifa)))
xpf(2, 1:nvfac) = coord(2, ipt(ivft(1:nvfac, ifa)))
!
if(ncurv.eq.0)then
xpf(1:2, 3) = 0.5d0*(xpf(1:2, 1)+xpf(1:2, 2))
elseif(ncurv.eq.1)then

else
print*,'Implement in the future!'
endif
!
do ig = 1, nvfac!...(2)ig = 1,ngausd
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
bti(1) = 1.d0

if(npoly.ge.1)then
 bti(2) = (rt - rc)/dr
 bti(3) = (st - sc)/ds

!DGP2
 if(npoly.eq.2)then
 bti(4) = 0.5d0*bti(2)*bti(2) - geoel(22, ielem)
 bti(5) = 0.5d0*bti(3)*bti(3) - geoel(23, ielem)
 bti(6) =       bti(2)*bti(3) - geoel(24, ielem)
 endif
endif
!
magd = 0.d0
do ideg = 1, ndegr
magd(1,1) = magd(1, 1) + unkgd(ideg, 1, ielem)*bti(ideg)
magd(1,2) = magd(1, 2) + unkgd(ideg, 2, ielem)*bti(ideg)
magd(2,1) = magd(2, 1) + unkgd(ideg, 3, ielem)*bti(ideg)
magd(2,2) = magd(2, 2) + unkgd(ideg, 4, ielem)*bti(ideg)
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
xpf(1, 1:nvfac) = coord(1, ipq(ivfq(1:nvfac, ifa)))
xpf(2, 1:nvfac) = coord(2, ipq(ivfq(1:nvfac, ifa)))
!
if(ncurv.eq.0)then
xpf(1:2, 3) = 0.5d0*(xpf(1:2, 1)+xpf(1:2, 2))
else

endif
!
do ig = 1, nvfac!...(2)ig = 1,ngausd
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

bqi(1) = 1.d0

 if(npoly.ge.1)then
  bqi(2) = (rq - rc)/dr
  bqi(3) = (sq - sc)/ds
!DGP2
  if(npoly.eq.2)then
  bqi(4) = 0.5d0*bqi(2)*bqi(2) - geoel(22, ielem)
  bqi(5) = 0.5d0*bqi(3)*bqi(3) - geoel(23, ielem)
  bqi(6) =       bqi(2)*bqi(3) - geoel(24, ielem)
 endif
endif
!
magd = 0.d0
do ideg = 1, ndegr
magd(1,1) = magd(1, 1) + unkgd(ideg, 1, ielem)*bqi(ideg)
magd(1,2) = magd(1, 2) + unkgd(ideg, 2, ielem)*bqi(ideg)
magd(2,1) = magd(2, 1) + unkgd(ideg, 3, ielem)*bqi(ideg)
magd(2,2) = magd(2, 2) + unkgd(ideg, 4, ielem)*bqi(ideg)

!if(ie.eq.1.and.ifa.eq.3.and.abs(unkgd(5, 3, ielem)).gt.1e-14)then
!print*,'ideg',magd(2,1),unkgd(ideg, 3, ielem)*bq(ideg),unkgd(ideg, 3, ielem),bq(ideg)
!endif
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
!
!if(ie.eq.1.and.ifa.eq.3)then
!print*,'geosq',ifa,ig,igq(ig, ifa),vnorm(1:2),cmagd(:, :),&
!unkgd(1:6, 3, ielem)
!endif

enddo !...ig from 1...2

enddo !...ifa from 1...4
!
200 enddo  !...(1)ifa=1,nquad
!
end subroutine getfnds_laghybrid_gd
!
!...Face integral using gauss quadrature distribution on cuvred quads using the gradient deformation...
!
subroutine rhsifacedg_lagquad_simpson_gd(ipqua,  unkno, ustar, fstarq, gesgq, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:12,1:nquad),       intent(in)::fstarq !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::rhsel
real*8,dimension(1:3, 1:ngesgq, 1:nquad),    intent(in)::gesgq
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
real*8::bg(ndegr, nvfac)
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

!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;

posi(1, 1) = -1.d0; posi(1 ,2 )= 1.d0; posi(1 ,3 )= 0.d0
weigh(1) = 1.d0/6.d0; weigh(2) = 1.d0/6.d0; weigh(3) = 4.d0/6.d0

!...Local gauss point No. of any gauss point in one face...
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;
!
endif

!...Zero out plnpn, ulnpn
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
!...Loop over every quad
!
do 550 ie = 1,nquad !...(1)ie = 1,nquad
!
ielem = ie + ntria

!...Zero out ulnpn, plnpn, elnpn
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0

ipq(1:nvqua) = ipqua(1:nvqua ,ie)
vnorm(1:3,  1:12) = gesgq(1:3, 1:12, ie)

!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...Loop over every face
do ifa = 1, 4
do ig   = 1, nvfac
wi  = weigh(ig)
!
rg = xvq(fglvq(ig, ifa))
sg = yvq(fglvq(ig, ifa))
!
bg(1,  ig) = 1.d0
bg(2,  ig) = (rg-rc)/dr
bg(3,  ig) = (sg-sc)/ds

!DGP2
if(npoly.eq.2)then
bg(4,  ig) = 0.5d0*bg(2,  ig)*bg(2,  ig) - geoel(19, ielem)
bg(5,  ig) = 0.5d0*bg(3,  ig)*bg(3,  ig) - geoel(20, ielem)
bg(6,  ig) =       bg(2,  ig)*bg(3,  ig) - geoel(21, ielem)
endif

!...The vertex constituting one cell...
gpnx = vnorm(1, fglgq(ig, ifa))
gpny = vnorm(2, fglgq(ig, ifa))
gpsa = vnorm(3, fglgq(ig, ifa))

!...Distribute to every corner...
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipq(fglvq(ig, ifa)))*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig) +&
ustar(2, ipq(fglvq(ig, ifa)))*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipq(fglvq(ig, ifa)))*fstarq(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig) +&
ustar(2, ipq(fglvq(ig, ifa)))*fstarq(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)

!...Output for debugging
!if(ie==1) print*,'rhs iface idegr',ipq(fglvq(ig, ifa)),ig,ifa,fglgq(ig, ifa),fstarq(1, fglgq(ig, ifa), ie),weigh(ig)

enddo
!
enddo
!
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
550 enddo
!
end subroutine rhsifacedg_lagquad_simpson_gd
!
!...Face integral (mass center) for hybrid quad (MEM) wit SMS using Simpson's rule...
!
subroutine rhsifacedg_lagsms_simpson_gd(ipqua, unkno, ustar,fstarq, gesgq, geoel,&
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
!if(ie==1)  print*,'rhs iface',ielem, ie, plnpn(1, 1:ndegr),fstarq(1, 1, ivsg, isg, ie)*bq(1:ndegr, ipqsg(ivsg, isg))*wfgsq(1, ivsg, isg)
650 enddo
!
end subroutine rhsifacedg_lagsms_simpson_gd
!
!....domain integral for hybrid curv quad cells using the gradient deformation
!
subroutine rhsdomndg_lagquadc_gd(intfac, ipqua, coord, coold, geoel, unkno, unkgd,strnq_devtp, rhsel,aflim,afvec, vnulq )
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
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(in)::vnulq
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
real*8,dimension(1:ndimn, 1:npqua) :: xpqi
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:ndegr):: b, bi,dbdr, dbds, bv
real*8:: unknod(1:nq)
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
real*8, dimension(1: ndimn):: vgnul
real*8, dimension(1: ndimn, 1:ndimn)::sigmg
!
real*8::sigma_devt(2, 2)!,strnq_devtp(2, 2)
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rci,sci
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::rhoad,uadv,vadv,eadv,rhoma,einta, sdadv
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr
real*8:: rhoi, rhon
real*8:: a11, a12, a21, a22
real*8:: xcrho ,ycrho, xgaus,ygaus
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
call ruqope(2, ngausdq, posiq, weighq)
!
!...Loop over quad
!
!...Initilize the stress
sigma_devt  = 0.d0
!strnq_devtp = 0.d0

do 650 ie = 1,nquad!...(1)ie = 1,nelem

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

!...Points consitituting one element...
if(ncurv==0)then
xpqi(1, 1:4) = coold(1, ipqua(1:4,ie))
xpqi(2, 1:4) = coold(2, ipqua(1:4,ie))
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))
!


elseif(ncurv==1)then
xpqi(1, 1:npqua) = coold(1,ipqua(1:npqua, ie))
xpqi(2, 1:npqua) = coold(2,ipqua(1:npqua, ie))

xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif

!...Geometry parameters for reference cell...
dr = 1.d0
ds = 1.d0

!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...The initial cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)

!...Physical cell center
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi(1:2,1:nvqua), rc, sc, xcrho, ycrho)

!...The derivatives of basis function...
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
!
jacbf = 0.d0
!
bi(1) = 1.d0
bi(2) = (r - rci)/dr
bi(3) = (s - sci)/ds
!DGP2
if(npoly.eq.2)then
bi(4) = 0.5d0*bi(2)*bi(2) - geoel(22, ielem)
bi(5) = 0.5d0*bi(3)*bi(3) - geoel(23, ielem)
bi(6) =       bi(2)*bi(3) - geoel(24, ielem)
endif
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

!...Cofactor matrix of Jacobian transformation matrix
jacbg(1, 1) = a22; jacbg(1, 2) =-a21
jacbg(2, 1) =-a12; jacbg(2, 2) = a11

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

!...Derived variables
einta = eadv - 0.5d0*(uadv**2+vadv**2)

!...Call EOS
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi(1:2,1:nvqua), xg, yg, xgaus, ygaus)
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)

!rhoi=1.845d0
call GetEOS(nmatel,ncase,gamlg,rhoad, einta, rhoi, pres, sdadv, ielem)

!pres = max(eps,(gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2)))

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
!...Null mode at gauss points....
vgnul = 0.d0
do ishp = 1, npqua
vgnul(1:2)= vgnul(1:2) + shpq(ishp)*vnulq(1:2, ishp, ie)
enddo
uadv = uadv !- vgnul(1)
vadv = vadv !- vgnul(2)
pres = pctr + aflim(4, ielem)*(pres- pctr)
!
endif
!
if(nmatel.eq.2)then
!  call GetStress_deviat_finite
!  call GetStress_deviat_infsmal(r,s,rci,sci,geoel(22:24, ielem),sigma_devt,&
!                                             strnq_devtp(:,:,ngstrnf+ig,ielem), unkgd(:,:,ielem),ielem, 1)
  call GetStress_deviatfem_infsmal(r,s,coord(1:2,ipqua(1:nvqua, ie)),coold(1:2,ipqua(1:nvqua, ie)),&
                                 sigma_devt, strnq_devtp(:,:,ngstrnf+ig,ielem), ielem, 1)
endif
!
!...Stress tensor
sigmg(1, 1) = -pres + sigma_devt(1, 1);
sigmg(1, 2) =  0.d0 + sigma_devt(1, 2);
sigmg(2, 1) =  0.d0 + sigma_devt(2, 1);
sigmg(2, 2) = -pres + sigma_devt(2, 2);

!...High-order
!do ideg = 1,-ndegr
!fluxd(ideg,1) = gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv
!fluxd(ideg,2) = gdshp(1, ideg)*(-pres)
!fluxd(ideg,3) = gdshp(2, ideg)*(-pres)
!fluxd(ideg,4) = (gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv)*(-pres)
!enddo

do ideg = 1, ndegr
fluxd(ideg,1) = gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv
fluxd(ideg,2) = gdshp(1, ideg)*sigmg(1, 1) + gdshp(2, ideg)*sigmg(1, 2)
fluxd(ideg,3) = gdshp(1, ideg)*sigmg(2, 1) + gdshp(2, ideg)*sigmg(2, 2)
fluxd(ideg,4) = gdshp(1, ideg)*(sigmg(1,1)*uadv+sigmg(1,2)*vadv) + gdshp(2, ideg)*(sigmg(2,1)*uadv+sigmg(2,2)*vadv)
enddo

!finally, scatter the contribution to the RHS
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) - fluxd(ideg,1:nq)*djaco
enddo

!...Output for debugging
!if(ie==8) print*,'rhs iface idegr',ie,ig,gdshp(1:2,2),uadv,vadv,djaco
!
enddo !...(2)ig = 1,ngausd
!
650 enddo

end subroutine rhsdomndg_lagquadc_gd
!
!....domain integral for hybrid curv quad cells using the gradient deformation
!...Jacobian matrix from FEM, Stress from
!
subroutine rhsdomndg_lagquadc_gdmix(intfac, ipqua, coord, coold, geoel, unkno, unkgd,strnq_devtp, rhsel,aflim,afvec, vnulq )
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
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(in)::vnulq
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
real*8,dimension(1:ndimn, 1:npqua) :: xpqi
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:ndegr):: b, bi,dbdr, dbds, bv
real*8:: unknod(1:nq)
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
real*8, dimension(1: ndimn):: vgnul
real*8, dimension(1: ndimn, 1:ndimn)::sigmg
!
real*8::sigma_devt(2, 2)!,strnq_devtp(2, 2)
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rci,sci
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::rhoad,uadv,vadv,eadv,rhoma,einta, sdadv
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr
real*8:: rhoi, rhon
real*8:: a11, a12, a21, a22
real*8:: xcrho ,ycrho, xgaus,ygaus
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
call ruqope(2, ngausdq, posiq, weighq)
!
!...Loop over quad
!
!...Initilize the stress
sigma_devt  = 0.d0
!strnq_devtp = 0.d0

do 650 ie = 1,nquad!...(1)ie = 1,nelem

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

!...Points consitituting one element...
if(ncurv==0)then
xpqi(1, 1:4) = coold(1, ipqua(1:4,ie))
xpqi(2, 1:4) = coold(2, ipqua(1:4,ie))
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))
!


elseif(ncurv==1)then
xpqi(1, 1:npqua) = coold(1,ipqua(1:npqua, ie))
xpqi(2, 1:npqua) = coold(2,ipqua(1:npqua, ie))

xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif

!...Geometry parameters for reference cell...
dr = 1.d0
ds = 1.d0

!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...The initial cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)

!...Physical cell center
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi(1:2,1:nvqua), rc, sc, xcrho, ycrho)

!...The derivatives of basis function...
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
!
jacbf = 0.d0
!
bi(1) = 1.d0
bi(2) = (r - rci)/dr
bi(3) = (s - sci)/ds
!DGP2
if(npoly.eq.2)then
bi(4) = 0.5d0*bi(2)*bi(2) - geoel(22, ielem)
bi(5) = 0.5d0*bi(3)*bi(3) - geoel(23, ielem)
bi(6) =       bi(2)*bi(3) - geoel(24, ielem)
endif

!
a11 = dxdr
a12 = dxds

a21 = dydr
a22 = dyds

!...Cofactor matrix of Jacobian transformation matrix
jacbg(1, 1) = a22; jacbg(1, 2) =-a21
jacbg(2, 1) =-a12; jacbg(2, 2) = a11

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

!...Derived variables
einta = eadv - 0.5d0*(uadv**2+vadv**2)

!...Call EOS
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi(1:2,1:nvqua), xg, yg, xgaus, ygaus)
call getrhog_initial(rhoi,  xgaus, ygaus, xcrho, ycrho)

!rhoi=1.845d0
call GetEOS(nmatel,ncase,gamlg,rhoad, einta, rhoi, pres, sdadv, ielem)

!pres = max(eps,(gamlg-1.d0)*rhoad*(eadv - 0.5d0*(uadv**2 + vadv**2)))

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
!...Null mode at gauss points....
vgnul = 0.d0
do ishp = 1, npqua
vgnul(1:2)= vgnul(1:2) + shpq(ishp)*vnulq(1:2, ishp, ie)
enddo
uadv = uadv !- vgnul(1)
vadv = vadv !- vgnul(2)
pres = pctr + aflim(4, ielem)*(pres- pctr)
!
endif
!
if(nmatel.eq.2)then
!  call GetStress_deviat_finite  GetStress_deviat_infsmal
call GetStress_deviat_infsmal(r,s,rci,sci,geoel(22:24, ielem),sigma_devt,&
strnq_devtp(:,:,ngstrnf+ig,ielem), unkgd(:,:,ielem),ielem, 1)
!call GetStress_deviatfem_infsmal(r,s,coord(1:2,ipqua(1:nvqua, ie)),coold(1:2,ipqua(1:nvqua, ie)),&
!sigma_devt, strnq_devtp(:,:,ngstrnf+ig,ielem), ielem, 1)
endif
!
!...Stress tensor
sigmg(1, 1) = -pres + sigma_devt(1, 1);
sigmg(1, 2) =  0.d0 + sigma_devt(1, 2);
sigmg(2, 1) =  0.d0 + sigma_devt(2, 1);
sigmg(2, 2) = -pres + sigma_devt(2, 2);

!...High-order
!do ideg = 1,-ndegr
!fluxd(ideg,1) = gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv
!fluxd(ideg,2) = gdshp(1, ideg)*(-pres)
!fluxd(ideg,3) = gdshp(2, ideg)*(-pres)
!fluxd(ideg,4) = (gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv)*(-pres)
!enddo

do ideg = 1, ndegr
fluxd(ideg,1) = gdshp(1, ideg)*uadv + gdshp(2, ideg)*vadv
fluxd(ideg,2) = gdshp(1, ideg)*sigmg(1, 1) + gdshp(2, ideg)*sigmg(1, 2)
fluxd(ideg,3) = gdshp(1, ideg)*sigmg(2, 1) + gdshp(2, ideg)*sigmg(2, 2)
fluxd(ideg,4) = gdshp(1, ideg)*(sigmg(1,1)*uadv+sigmg(1,2)*vadv) + gdshp(2, ideg)*(sigmg(2,1)*uadv+sigmg(2,2)*vadv)
enddo

!finally, scatter the contribution to the RHS
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) - fluxd(ideg,1:nq)*djaco
enddo

!...Output for debugging
!if(ie==8) print*,'rhs iface idegr',ie,ig,gdshp(1:2,2),uadv,vadv,djaco
!
enddo !...(2)ig = 1,ngausd
!
650 enddo

end subroutine rhsdomndg_lagquadc_gdmix
!
!...Face integral (mass center) for rhs of the deformation gradient on curvilinear meshes...
!
subroutine rhsifacegd_lagquadc(ipqua, ustar, gesgq0, geoel,&
rhsgd)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndegr,1:4,1:ncell),        intent(out)::rhsgd
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq0
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvqua) :: ipq
real*8, dimension(1:ndegr, 1:4) :: ulnpn
real*8::xvq(nvqua), yvq(nvqua),bg(1:ndegr,1:nvfac)
real*8::weigh(nvfac)
integer, dimension(3, 4)::fglvq,fglgq

!...local real number
real*8::eps,c00,c05,c10,c20,c13,c16
real*8::dr,ds,rc,sc
real*8::dwav1,dwav2
real*8::anx, any
real*8::gpnx, gpny,gpsa,rg,sg
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
data c16   / 0.1666666666666666d0 /

!
!...Part I: Preliminary data preparation
!
!...Local vertex No. of gauss points in one unit ...
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;

!
weigh(1) = 1.d0/6.d0; weigh(2) = 1.d0/6.d0; weigh(3) = 4.d0/6.d0

!...Local gauss point No. of any gauss point in one face...
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;

!...The coordinates for vertices
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
if(ncurv.eq.1)then
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
!
xvq(9) =  0.d0; yvq(9) =  0.d0
endif
!
!...Part II: Loop over every quads...
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
rc= geoel(7, ielem) !...mass center...
sc= geoel(8, ielem)

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
bg(2,  ig) = (rg-rc)/dr
bg(3,  ig) = (sg-sc)/ds

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
ustar(1, ipq(fglvq(ig, ifa)))*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig)!*lpnpq(1, 1:ndegr, 1, iv)
!
ulnpn(1:ndegr, 2)  = ulnpn(1:ndegr, 2)+&
ustar(1, ipq(fglvq(ig, ifa)))*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)
!
ulnpn(1:ndegr, 3)  = ulnpn(1:ndegr, 3)+&
ustar(2, ipq(fglvq(ig, ifa)))*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig)!*lpnpq(1, 1:ndegr, 1, iv)
!
ulnpn(1:ndegr, 4)  = ulnpn(1:ndegr, 4)+&
ustar(2, ipq(fglvq(ig, ifa)))*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)
!
!if(ie==1) print*,ifa,ig,ielem,ipq(fglvq(ig, ifa)), gpnx,gpsa,bg(4, ig),weigh(ig),ulnpn(4, 1),ustar(1, ipq(fglvq(ig, ifa)))
!
enddo
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
end subroutine rhsifacegd_lagquadc
!
!...Get thr RHS for the deformation gradient in curved cells
!
subroutine rhsdomngd_lagquadc(intfac, ipqua, coord, coold, geoel, unkno, unkgd,rhsgd,aflim,afvec )
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
real*8,dimension(1:ndimn, 1:npqua) :: xpq
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
data eps   / 1.0d-6 /
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
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
!...Points consitituting one element...
if(ncurv==0)then
xpq(1, 1:4) = coold(1, ipqua(1:4,ie))
xpq(2, 1:4) = coold(2, ipqua(1:4,ie))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coold(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coold(2,ipqua(1:npqua, ie))
endif

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

!...Coefficient R of RZ or XY system...
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
!if(ie.eq.1)print*,'jacbg',jacbg(1,1:2),dbdr(4),dbds(4),r,s
!
!...Calculate G dot dbdx or dbdy
!
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
!if(ielem.eq.9) print*,'domn2',pres,rhoad,unkno(1:3,1,ielem),b(1:3),xg,yg

if(nlimi.eq.6.and.geoel(10,ielem).gt.10.d0)then
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
!...High-order
do ideg = 1, ndegr
fluxd(ideg,1) = gdshp(1, ideg)*uadv
fluxd(ideg,2) = gdshp(2, ideg)*uadv
fluxd(ideg,3) = gdshp(1, ideg)*vadv
fluxd(ideg,4) = gdshp(2, ideg)*vadv
enddo

!if(ielem.eq.213) print*,'domn1',gdshp(2,1:3),uadv,ig,rhsgd(1:3,2,ielem)

!finally, scatter the contribution to the RHS
do ideg = 1,ndegr
rhsgd(ideg,1:4,ielem)=rhsgd(ideg,1:4,ielem) - fluxd(ideg,1:4)*djaco*rcoef
enddo
!
!if(ielem.eq.1) print*,'domn',gdshp(1,4),uadv,djaco,ig,rhsgd(4,1,ielem),fluxd(4, 1)
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
!
!print*,'ieleme.q.2',rhsel(1:3,3,2)
!
end subroutine rhsdomngd_lagquadc
!
!...Find terms in geoel for high-order for the gradient deformation...
!
subroutine getgeoel_lag_hogd(iptri, ipqua, geoel, coord)
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
rhog = 1.d0 + 0.995d0*sin(pi*xgaus)

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
rhog = 1.d0 + 0.995d0*sin(pi*xgaus)

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
end subroutine getgeoel_lag_hogd
!
!...subroutine: Riemann solver for curved quad (MEM) with SMS using Simpson rule....
!
subroutine getriem_quadsms_simpson_gd(ipqua, geoel, gesgq, vlave, unkno, unkgd, strnq_devtp, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq,coord, coold, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
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

!...Local integer
integer::ie, ideg, ielem, ifa, iv, isg, ivsg, ifsg

!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(8, 4)::fnqsg
integer,dimension(4, 4)::ipqsg, ipqsg_strn
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
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2), sigma_devt(1:2,1:2)
!real*8,dimension(1:ndimn,1:ndimn)::strnq_devtp
!
real*8::eps,c00,c05,c10,c20
real*8::rhovt, rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr, eintc
real*8::rhomv,uvtx,vvtx,evtx, pvtx,rhovsg,eintv,sdv
real*8::dux,duy,deltu
real*8::dr, ds, rc, sc, r, s,rcv,scv, rci,sci
real*8::xcrho, ycrho
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
!
wfgsq = 0.25d0*wfgsq
!
!...Part II: Loop over every quad...
!

!...Initilize the stress
 sigma_devt  = 0.d0
! strnq_devtp = 0.d0
!print*,'sigma',sigma_devt
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
!
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

!
!call getrhosubg_averg(xpq,xpqi,unksgq, geoq_sub,ielem)

call getrhosubcell_sms(xcrho, ycrho, xpq,xpqi,unksgq, ielem)

!...Get density correction
!call getdens_quadsubg(rc, sc, geoel(19:21, ielem), unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc)
call getrhosubg_daverg(rc,sc, geoel(19:21, ielem), xpq, unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc, ielem)

!...Output for debugging
!if(ielem.eq.1)then
!print*,'Variabe',ielem,ivsg,rhovsg,unksgq(1,1),1.d0/uqsgc(1, 1)
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

!...SMS for ndens=1
if(ndens.eq.1)then

rhovt = 1.d0/unknvq(1, ipqsg(ivsg, isg))
rhovsg = rhovt +cdrho*(unksgq(1,isg)-1.d0/uqsgc(1, isg))

elseif(ndens.eq.3)then

Print*,'SMS for ndens=3 will be implemented in future'
stop

endif

uvtx = unknvq(2, ipqsg(ivsg, isg))
vvtx = unknvq(3, ipqsg(ivsg, isg))
evtx = unknvq(4, ipqsg(ivsg, isg))

!...Derived variables at the vertex
eintv = evtx - 0.5d0*(uvtx**2 + vvtx**2)

!...Call EOS
call getrhog_initial(rhoi,  xpqi(1, ipqsg(ivsg, isg)), xpqi(2, ipqsg(ivsg, isg)), xcrho, ycrho)
call GetEOS(nmatel,ncase,gamlg,rhovsg, eintv, rhoi, pvtx, sdv, ielem)

!pvtx = max(eps, (gamlg-1.d0)*rhovsg*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx

!...Output for debugging
!if(ielem.eq.8)then
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

!...Updtae velocity
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx
endif

!
if(nmatel.eq.2)then
!  call GetStress_deviat_finite GetStress_deviat_infsmal
  call GetStress_deviat_infsmal(xvq(ipqsg(ivsg, isg)),yvq(ipqsg(ivsg, isg)),rci,sci,geoel(22:24, ielem),&
                                             sigma_devt, &
                                             strnq_devtp(:,:, ipqsg(ivsg, isg),ielem), unkgd(:,:,ielem),ielem,ipqsg_strn(ivsg, isg))
!  call GetStress_deviatfem_infsmal(xvq(ipqsg(ivsg, isg)),yvq(ipqsg(ivsg, isg)),xpq,xpqi,sigma_devt,&
!                                   strnq_devtp(:,:, ipqsg(ivsg, isg),ielem), ielem,ipqsg_strn(ivsg, isg))
endif

!...Get stress tensor at one vertex
!sigma(1, 1, ivsg) = -pvtx
!sigma(1, 2, ivsg) = 0.d0
!sigma(2, 1, ivsg) = 0.d0
!sigma(2, 2, ivsg) = -pvtx

sigma(1, 1, ivsg) = -pvtx + sigma_devt(1,1)
sigma(1, 2, ivsg) = 0.d0  + sigma_devt(1,2)
sigma(2, 1, ivsg) = 0.d0  + sigma_devt(2,1)
sigma(2, 2, ivsg) = -pvtx + sigma_devt(2,2)
!
!...Output for debugging
!if(ielem.eq.1.or.ielem.eq.2281)then
!print*,'smooth2',ielem,isg,ivsg,ipqsg(ivsg, isg),sigma_devt
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
!if(ipq(iv).eq.20) print*,'p36 muacn(vv) post',ie,ifa,isg,ivsg,munacn_rie(1, 1),murie(ifa, ivsg),vnorm(1:3, ifa, ivsg),&
!sigma(1, 1, ivsg),unsgq(2:3, ivsg)

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

end subroutine getriem_quadsms_simpson_gd
!
!...Calculate the velocity at the Gauss point for novertex face...
!
subroutine getgesgq(intfac,ipqua,gesgq)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(inout)::gesgq
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(ngausf)::fgausl, fgausr
integer,dimension(nvfac, 4)::fglvq
integer,dimension(ngausf,4)::fglgq
real*8::vnorm(3,3)
!
!...Initial velocity at gauss point
!..Initial Gauss point velocity from the interpolation of the 3-node velocity...
!
!
!...Part I: Specify some gauss points


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

endif
!
!...Part II: Get the corrected velocity at the non-vertex gauss point...
!
fgausl = 0.d0
fgausr = 0.d0
!
do 450 ifa = nbfac+1, nafac !...(1)ifa = 1, nafac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)

!...Boundary face
!...Identify the fgaus for left cell
call getfg_glb_gd(ipqua(1:nvqua, iel), ipf, fgausl)

!...Identify the fgaus for right cell
call getfg_glb_gd(ipqua(1:nvqua, ier), ipf, fgausr)
!
vnorm = 0.d0

do ig =1, nvfac
vnorm(1:2 ,ig) = 0.5d0*(gesgq(1:2,fgausl(ig), iel) - gesgq(1:2,fgausr(ig), ier))
vnorm(3 ,ig)= 0.5d0*(gesgq(3,fgausl(ig), iel) + gesgq(3,fgausr(ig), ier))
enddo
!
do ig =1, nvfac
gesgq(1:2,fgausl(ig), iel) = vnorm(1:2 ,ig)
gesgq(3,  fgausl(ig), iel) = vnorm(3 ,  ig)
!
gesgq(1:2,fgausr(ig), ier) =-vnorm(1:2 ,ig)
gesgq(3,  fgausr(ig), ier) = vnorm(3 ,  ig)
enddo
!

450 enddo
!
!print*,'dufg2',dufgq(1, 17:18, 99)
!
end subroutine getgesgq
!
!...Find the global no. of one face gauss point
!
subroutine getfg_glb_gd(ipqua, ipf, fgaus)
use constant
implicit none
integer, dimension(nvqua), intent(in):: ipqua
integer, dimension(nvfac), intent(in):: ipf
integer:: ig, nsum
integer, dimension(1:ngausf), intent(out)::fgaus
!
integer::fglgq(3,4)

!...Give fglgq
!...Local gauss point No. of any gauss point in one face...
fglgq(1, 1) = 1;  fglgq(2, 1) = 2; fglgq(3, 1) = 9;
fglgq(1, 2) = 3;  fglgq(2, 2) = 4; fglgq(3, 2) =10;
fglgq(1, 3) = 5;  fglgq(2, 3) = 6; fglgq(3, 3) =11;
fglgq(1, 4) = 7;  fglgq(2, 4) = 8; fglgq(3, 4) =12;


!...Find the No. of face gauss points in fglgq...
nsum = 3
!
do ig = 1, 3
!
if(ipqua(5).eq.ipf(3))then
!
if(ig.eq.3)then
fgaus(ig) = fglgq(ig,1)
else
if(ipf(1).eq.ipqua(1))then
fgaus(ig) = fglgq(ig,1)
else
fgaus(ig) = fglgq(nsum-ig,1)
endif
endif
!
elseif(ipqua(6).eq.ipf(3))then

if(ig.eq.3)then
fgaus(ig) = fglgq(ig,2)
else
if(ipf(1).eq.ipqua(2))then
fgaus(ig) = fglgq(ig,2)
else
fgaus(ig) = fglgq(nsum-ig,2)
endif
endif
!
elseif(ipqua(7).eq.ipf(3))then
!
if(ig.eq.3)then
fgaus(ig) = fglgq(ig,3)
else
if(ipf(1).eq.ipqua(3))then
fgaus(ig) = fglgq(ig,3)
else
fgaus(ig) = fglgq(nsum-ig,3)
endif
endif
!
elseif(ipqua(8).eq.ipf(3))then
!
if(ig.eq.3)then
fgaus(ig) = fglgq(ig,4)
else
if(ipf(1).eq.ipqua(4))then
fgaus(ig) = fglgq(ig,4)
else
fgaus(ig) = fglgq(nsum-ig,4)
endif
endif
!
endif

enddo
end subroutine getfg_glb_gd
