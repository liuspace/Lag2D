!
!...Subroutine for barth limiter based on vertex for primitive variables on quads with symmetry preserving....
!
subroutine barthlimit_lagsym_quad(geoel, coord, coold, ustar, unkno, ipqua, bface,intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord, coold
real*8,dimension(1:nq+1, 1:nsize),             intent(inout)::aflim
integer,  dimension(1:nvqua,1:nquad),          intent(in):: ipqua
integer*4,dimension(1:nbfai,1:nbfac),          intent(in)::bface
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
!...Part I: Some preparing work
!

!...Coloring the boundary node

indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6

!...shape functions

dr = 1.d0
ds = 1.d0
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
!++++++++++Print for debugging++++++++++
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
!

!...Get the max and min surrounding one cell

do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

 do iq=1, nq+2
  unmax_new(iq, ielem) = maxval(unmax(iq, ipq(1:nvqua)))
  unmin_new(iq, ielem) = minval(unmin(iq, ipq(1:nvqua)))
 enddo

enddo
!
!...Part 2: Impose limiter
!
do ie = 1, nquad

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

!...zero out unknv

unknvq = 0.d0
do iv   = 1,nvqua
 do ideg = 1,mdegr
  unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
 enddo

 if(ndens.eq.1)then
  rhov = 1.d0/unknvq(1, iv)
 elseif(ndens.eq.2)then
  r = xvq(iv); s = yvq(iv)

  xpq(1, 1:4) = coord(1, ipq(1:nvqua))
  xpq(2, 1:4) = coord(2, ipq(1:nvqua))
  xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
  xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

  call getrhoig_quad(rhoi, r, s, xpqi)!
  call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)

  rhov = rhon
 elseif(ndens.eq.3)then
  rcv = geoel(5, ielem); scv = geoel(6, ielem)

  bqv(1, iv) = 1.d0
  bqv(2, iv) = (xvq(iv)-rcv)/dr
  bqv(3, iv) = (yvq(iv)-scv)/ds

  unknvq(1, iv) =0.d0
  do ideg = 1,mdegr
   unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
  enddo
  rhov = unknvq(1, iv)
 endif !if(ndens.eq.1)then

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
 if(ndens.eq.1)then
  unknvq(1, iv) = 1.d0/rhov
 elseif(ndens.eq.2)then
  unknvq(1, iv) = rhov
 elseif(ndens.eq.3)then
  unknvq(1, iv) = rhov
 endif

 unknvq(4 ,iv) = pvtx !do iv   = 1,nvqua
enddo
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
  unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
  unctr(1)   = rhoct
elseif(ndens.eq.3)then
  unctr(1)   = rhoct
endif
  unctr(2:3) = unkno(1, 2:3, ielem)
  unctr(nq)  = pctr
!
 do iv = 1, nvqua
  do iq = 1, nq
   dunk(iq) = unknvq(iq, iv) - unctr(iq)
   call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
!  call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)

   alfa(iq, iv) = afbar
  enddo !do iq = 1, nq

!...Special treatment of Boundary condition

  if(indbd(ipq(iv)).eq.1)then

  !...Special treatment of Boundary condition for any case...
  !...This makes the symmetry BC different from the whole domain

  !alfa(1, iv) = 1.d0
  !alfa(4, iv) = 1.d0
!
  if(ncase.eq.1)then
   if(coord(1, ipq(iv)).lt.1.d-6.or.abs(coord(1, ipq(iv))-1.d0).lt.1.d-6) then
    !alfa(1, iv) = 1.d0
    !alfa(4, iv) = 1.d0
    !call  barthfct(unmax(2, ipq(iv)), unmin(2, ipq(iv)), unctr(2), dunk(2), afbar)
    !fiy = (- unctr(2))/(dunk(2))
    alfa(2, iv) = min(max((- unctr(2))/(dunk(2)),0.d0),1.d0)
    !alfa(2, iv) =  (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)
   endif
!
   if(coord(2, ipq(iv)).lt.1.d-6.or.abs(coord(2, ipq(iv))-1.d0).lt.1.d-6) then
    !alfa(1, iv) = 1.d0
    !alfa(4, iv) = 1.d0
    !call  barthfct(unmax(3, ipq(iv)), unmin(3, ipq(iv)), unctr(3), dunk(3), afbar)
    ! fiy = ( - unctr(3))/(dunk(3))
    alfa(3, iv) = min(max(( - unctr(3))/(dunk(3)),0.d0),1.d0)
    !alfa(3, iv) = (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)
   endif
  endif !if(ncase.eq.1)then
 endif !if(indbd(ipq(iv)).eq.1)then
!
enddo ! do iv = 1, nvqua

!...Get teh minimum value for the cell

 do iq = 1,nq
  aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
 enddo
!
enddo  !do ie = 1, nquad
!
!...Part 2.1: Impose symmetry preserving limiter for velocity...
!
!call barthlimit_sympre_phy3(geoel, coord, ustar, unkno, ipqua, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
 call  barthlimit_sympre_quadlb(geoel, coord, ustar, unkno, ipqua, bface,intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
!
!...Degenerate to non-symmetric...
!
 do ie = 1, -nquad
  ielem = ie + ntria

  afvec(1, 1, ielem) = aflim(2, ielem)
  afvec(1, 2, ielem) = 0.d0
  afvec(2, 1, ielem) = 0.d0
  afvec(2, 2, ielem) = aflim(3, ielem)
enddo
!
!...The limitation of the first invariant...
!
!call barthlimit_sympre_inva(geoel, coord, ustar, unkno, ipqua, intfac, afvec)
!call barthlimit_sympre_inva2(geoel, coord, ustar, unkno, ipqua, intfac, afvec, aflim)
!
!...Part 3: Correct total energy
!
do ie = 1,nquad

  ielem = ie + ntria
  ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
  rc= geoel(1, ielem) !...mass center...
  sc= geoel(2, ielem)

!...Shape function

  do iv =1 ,nvqua
   bq(1, iv) = 1.d0
   bq(2, iv) = (xvq(iv)-rc)/dr
   bq(3, iv) = (yvq(iv)-sc)/ds
  enddo

!...Cell average of inverse density, velocity and total energy

  if(ndens.eq.1)then
   rhom = unkno(1, 1, ielem)
  elseif(ndens.eq.2)then
   rhom = 1.d0/unkno(1, 1, ielem)
  elseif(ndens.eq.3)then
   rhom = 1.d0/unkno(1, 1, ielem)
  endif
!
 uctr = unkno(1, 2, ielem)
 vctr = unkno(1, 3, ielem)
 ectr = unkno(1, 4, ielem)
!
 rhoct  = 1.d0/rhom
 pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
 if(ndens.eq.1)then
  unctr(1)   = 1.d0/rhoct
 elseif(ndens.eq.2)then
  unctr(1)   = rhoct
 elseif(ndens.eq.3)then
  unctr(1)   = rhoct
 endif

 unctr(2:3) = unkno(1, 2:3, ielem)
 unctr(nq) = pctr
 unctr(nq+1) = ectr

!...zero out unknv

 unknvq = 0.d0
 do iv   = 1,nvqua
  do ideg = 1,mdegr
   unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
  enddo
!
 if(ndens.eq.1)then
  rhov = 1.d0/unknvq(1, iv)
 elseif(ndens.eq.2)then
  r = xvq(iv); s = yvq(iv)

  xpq(1, 1:4) = coord(1, ipq(1:nvqua))
  xpq(2, 1:4) = coord(2, ipq(1:nvqua))
  xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
  xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

  call getrhoig_quad(rhoi, r, s, xpqi)!
  call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)

  rhov = rhon
 elseif(ndens.eq.3)then
  rcv = geoel(5, ielem); scv = geoel(6, ielem)

  bqv(1, iv) = 1.d0
  bqv(2, iv) = (xvq(iv)-rcv)/dr
  bqv(3, iv) = (yvq(iv)-scv)/ds

  unknvq(1, iv) =0.d0
!
  do ideg = 1,mdegr
   unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
  enddo
  rhov = unknvq(1, iv)
 endif !if(ndens.eq.1)then

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
 if(ndens.eq.1)then
  unknvq(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
 elseif(ndens.eq.2)then
  unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
 elseif(ndens.eq.3)then
  unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
 endif
!
 dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
 duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
 dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
 dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)

!  unknvq(2, iv) = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
!  unknvq(3, iv) = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
 unknvq(4 ,iv) = pctr + aflim(4, ielem)*(pvtx - pctr)
 unknvq(5, iv) = unknvq(4 ,iv)/(gamlg-1.d0)*unknvq(1, iv) + 0.5d0*(unknvq(2, iv)**2 + unknvq(3, iv)**2)
enddo  !do iv   = 1,nvqua

!...Get the limiting coefficient for the total energy...

 do iv = 1, nvqua
  do iq = nq+1, nq+1 !...energy
   dunk(iq) = unknvq(iq, iv) - unctr(iq)

   call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
!  call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)
   alfa(iq, iv) = afbar
  enddo

!...Special treatment for the boundary nodes...

  if(indbd(ipq(iv)).eq.1)then
! alfa(5, iv) = 1.d0
  endif
 enddo ! do iv = 1, nvqua

!...Get teh minimum value for the whole cell

 do iq = nq+1,nq+1
  aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
 enddo

enddo !do ie = 1,nquad

!
!...Call the smooth indicator
!
!call getSI_divg(geoel, coord, unkno, ipqua, bface,intfac,esuv1, esuv2, unmax_new, unmin_new)
!
do ie = 1, -nquad
!
 if(log10(geoel(10 ,ie)).lt.-2.8d0)then
!
 aflim(1, ie) = geoel(10 ,ie)*aflim(1, ie) + (1.d0-geoel(10 ,ie))
 aflim(3, ie) = geoel(10 ,ie)*aflim(3, ie) + (1.d0-geoel(10 ,ie))
 aflim(5, ie) = geoel(10 ,ie)*aflim(5, ie) + (1.d0-geoel(10 ,ie))
!
 afvec(1, 1, ie) = geoel(10 ,ie)*afvec(1, 1, ie) + 1.d0-geoel(10 ,ie)
 afvec(1, 2, ie) = geoel(10 ,ie)*afvec(1, 2, ie)
 afvec(2, 1, ie) = geoel(10 ,ie)*afvec(2, 1, ie) + 1.d0-geoel(10 ,ie)
 afvec(2, 2, ie) = geoel(10 ,ie)*afvec(2, 2, ie)
endif

enddo
!
end subroutine barthlimit_lagsym_quad
!
!...subroutine: Get the nodal velocity U_p^* and pressure for hybrid meshes with general Riemann solver...
!
subroutine getndvelo_lag_mc_matrixsym(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
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
!usold = ustar
indnd = 0

!...Mark the boundary nodes...
!...Shockless Noh...
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
if(ntria.gt.0) call getriem_tria_matrixsym(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold,aflim, afvec)

!...Quad
if(nquad.gt.0) call getriem_quad_matrixsym(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)

!...Boundary condition
!call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)!
!call getboundary_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm)
 call getbc_lagmaire2(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm,itime)
!
!...Periodic boundary condition for 1D isentropic Sin problem...
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
!
!if(ipoin.eq.3630)print*,'Ustar',ipoin,munaci(1, 1),rhsu1,rhsu2
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

!...Zero normal velocity for BC...
! call getbcvn_lag(bface, intfac, gflag, ustar)
! call getbcve_exact(bface, intfac, gflag, ustar, coord, itime)
!
!...Part III: Get the Riemann forces for face integral...
!

!...Tria
do ie = 1, ntria
 ipt(1:nvtri) = iptri(1:nvtri,ie)
 ielem = ie

!...Riemann forces
 do iv = 1, nvtri
 do ifa =1, 2
 fstar(1, ifa, iv, ie) = snsigmlt(1, ifa, iv, ie) + &
                        munaclt(1, 1, ifa, iv, ie)*ustar(1, ipt(iv))+&
                        munaclt(2, 1, ifa, iv, ie)*ustar(2, ipt(iv)) - munault(1, ifa, iv, ie)
 fstar(2, ifa, iv, ie) = snsigmlt(2, ifa, iv, ie) + &
                        munaclt(1, 2,ifa, iv, ie)*ustar(1, ipt(iv))+&
                        munaclt(2, 2, ifa, iv, ie)*ustar(2, ipt(iv)) - munault(2, ifa, iv, ie)
 enddo
 enddo
enddo

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
!
 enddo
 enddo

enddo
!
deallocate (munacn, fpres)
deallocate (munacu, snsigm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lag_mc_matrixsym
!
!...subroutine: Riemann input for hybrid linear quad grids with general Riemann solver....
!
subroutine getriem_quad_matrixsym(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
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
real*8,dimension(1:ndegr, 1:nvqua)::bq,bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8::aujmp(1:3, 1:nvqua)
real*8::vnorm(1:3, 1:2, 1:nvqua)
real*8::sigma(1:2, 1:2, 1:nvqua)
real*8,dimension(1:2, 1:nvqua)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:nvqua):: rcoeq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xphq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi

!...Arrays for Riemann solver
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2),sgmas(2, 2)
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
data eps   / 1.0d-6 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
ax= 0.d0
bx = 0.d0
!
!...Part I: Loop every quad...
!
!...Initialize sigmas
sgmas = 0.d0
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

!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
 do iv =1 ,nvqua
 bq(1, iv) = 1.d0
!
 if(npoly.ge.1)then
  bq(2, iv) = (xvq(iv)-rc)/dr
  bq(3, iv) = (yvq(iv)-sc)/ds
!
!DGP2
  if(npoly.eq.2)then
   bq(4, iv) = 0.5d0*bq(2, iv)*bq(2, iv) - geoel(19, ielem)
   bq(5, iv) = 0.5d0*bq(3, iv)*bq(3, iv) - geoel(20, ielem)
   bq(6, iv) =       bq(2, iv)*bq(3, iv) - geoel(21, ielem)
 endif
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
vnorm(1:3, 1, 1) = gelagq(1:3, 4, ie); vnorm(1:3, 2, 1) = gelagq(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 1, 2) = gelagq(1:3, 1, ie); vnorm(1:3, 2, 2) = gelagq(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gelagq(1:3, 2, ie); vnorm(1:3, 2, 3) = gelagq(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 1, 4) = gelagq(1:3, 3, ie); vnorm(1:3, 2, 4) = gelagq(1:3, 4, ie) !...For point ip(3)

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
!
rhoct  = 1.d0/rhomc !...Cell center density...
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))

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

call getrhoig_quad(rhoi, r, s, xpqi)!
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
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))

!...Limiter
if(nlimi.eq.1.and.npoly.ge.1)then
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
elseif(nlimi.eq.6.and.npoly.ge.1)then
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

!...Get the artifical stress from LLNL
!if(ielem.eq.1)print*,'sgmas',xvq(iv),yvq(iv)
!call getAS_llnl(xvq(iv) ,yvq(iv), unkno(:,:,ielem), geoel(:,ielem) ,xphq, xpqi, sgmas, ielem)

!if(ielem.eq.2.or.ielem.eq.31)then
! print*,'sgmas',ielem,iv,sgmas
!endif
!
!if(abs(sgmas(1,1)).ge.1e-5)print*,'ielem',ielem,rkstg
!...Get stress tensor at the vertex
sigma(1, 1, iv) = -pvtx + sgmas(1, 1)
sigma(1, 2, iv) = 0.d0  + sgmas(1, 2)
sigma(2, 1, iv) = 0.d0  + sgmas(2, 1)
sigma(2, 2, iv) = -pvtx + sgmas(2, 2)!
!
!if(ielem.eq.921)print*,'sgmas2',iv,sigma(:,:,iv)

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

!...Sound speed at the center...
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedence coefficient...
do iv   = 1, nvqua
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
do ifa = 1, 2
! deltu = sqrt(dux**2 + duy**2)
 deltu = abs(dux*vnorm(1, ifa, iv) + duy*vnorm(2, ifa, iv))
 murie(ifa, iv) = rhoct*sdctr + cimpd*rhoct*slpdu*deltu
!murie(iv) = rhoct*(slpdu*deltu*0.5d0 + sqrt((slpdu*deltu*0.5d0)**2+sdctr**2))
enddo
enddo

!...Get the summed denominator and numerator: sum(mu*n*a_c)
do iv  = 1, nvqua
do ifa = 1, 2 !...Every corner consists of 2 faces...

!...Call Riemann solver...
 call getriecoef_matrixnew(murie(ifa, iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
                         unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
                         munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ifa, iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
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
!if(ipq(iv).eq.3630) print*,'epq',ipq(iv), murie(ifa, iv),ie,ifa,iv,sigma(1:2, 1:2, iv),vnorm(1:3, ifa,iv),&
!munacu(1:2,ipq(iv)),unknvq(2:3, iv),sigma(1,1,iv)
enddo
enddo
!
350 enddo  !...(1)ie = 1,nqaud!

end subroutine getriem_quad_matrixsym
!
!...Coefficients for Riemann-like problem (Burton & Morgan)matrix similar to Marie...
!
subroutine getriecoef_matrixnew(murie, area, vnorm, aujmp, velo, sigma, munacn, munacu, snsigm)
use constant
implicit none
real*8, intent(in):: murie
real*8, intent(in):: area
real*8, intent(in):: vnorm(ndimn), aujmp(3), velo(ndimn)
real*8, intent(in):: sigma(ndimn, ndimn)
!
real*8, intent(out):: munacn(1:ndimn, 1:ndimn)
real*8, intent(out):: munacu(ndimn), snsigm(ndimn)
!
real*8:: vnau
!
if(aujmp(3).lt.1.d12)then
vnau   =  1.d0
munacn(1, 1) = murie*area*abs(vnau); munacn(2, 1) = 0.d0; 
munacn(1, 2) = 0.d0;                 munacn(2, 2) = murie*area*abs(vnau);
!
!
munacu(1) = munacn(1, 1)*velo(1) + munacn(2, 1)*velo(2) 
munacu(2) = munacn(2, 1)*velo(1) + munacn(2, 2)*velo(2)
!
snsigm(1) = sigma(1, 1)*area*vnorm(1) + sigma(1, 2)*area*vnorm(2)!
snsigm(2) = sigma(2, 1)*area*vnorm(1) + sigma(2, 2)*area*vnorm(2) 
else
        !
vnau   =  vnorm(1)*aujmp(1) + vnorm(2)*aujmp(2)
munacn(1, 1) = murie*area*abs(vnau); munacn(2, 1) = 0.d0; 
munacn(1, 2) = 0.d0;                 munacn(2, 2) = murie*area*abs(vnau);
!
!
munacu(1) = munacn(1, 1)*velo(1) + munacn(2, 1)*velo(2) 
munacu(2) = munacn(2, 1)*velo(1) + munacn(2, 2)*velo(2)
!
snsigm(1) = sigma(1, 1)*area*vnorm(1) + sigma(1, 2)*area*vnorm(2)!
snsigm(2) = sigma(2, 1)*area*vnorm(1) + sigma(2, 2)*area*vnorm(2) 

endif
!
end subroutine getriecoef_matrixnew
!
!...Coefficients for Riemann-like problem (Burton & Morgan)matrix similar to Marie...
!
subroutine getriecoef_matrixnew2(murie, area, vnorm, aujmp, velo, sigma, munacn, munacu, snsigm)
use constant
implicit none
real*8, intent(in):: murie
real*8, intent(in):: area
real*8, intent(in):: vnorm(ndimn), aujmp(3), velo(ndimn)
real*8, intent(in):: sigma(ndimn, ndimn)
!
real*8, intent(out):: munacn(1:ndimn, 1:ndimn)
real*8, intent(out):: munacu(ndimn), snsigm(ndimn)
!
real*8:: vnau
!
vnau   =  1.d0
munacn(1, 1) = murie*area*abs(vnau); munacn(2, 1) = 0.d0;
munacn(1, 2) = 0.d0;                 munacn(2, 2) = murie*area*abs(vnau);
!
!
munacu(1) = munacn(1, 1)*velo(1) + munacn(2, 1)*velo(2)
munacu(2) = munacn(2, 1)*velo(1) + munacn(2, 2)*velo(2)
!
snsigm(1) = sigma(1, 1)*vnorm(1) + sigma(1, 2)*vnorm(2)!
snsigm(2) = sigma(2, 1)*vnorm(1) + sigma(2, 2)*vnorm(2)
!
end subroutine getriecoef_matrixnew2
!
!....domain integral for hybrid linear quad cells
!
subroutine rhsdomndg_lag_mc_quadsym(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec, itime )
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
integer,intent(in)::itime
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
real*8::sigmg(2,2),sgmas(2,2)
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
ipq(1:nvqua) = ipqua(1:4,ie)!
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:4) = coold(1, ipq(1:nvqua))
xpqi(2, 1:4) = coold(2, ipq(1:nvqua))
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

!...Stress tensor
sigmg(1, 1) = -pres; sigmg(1, 2) =  0.d0;
sigmg(2, 1) = 0.d0;  sigmg(2, 2) = -pres;
!
!if(ielem.eq.991)print*,'volume integral'
call getAS_llnl(xg,yg, unkno(:,:,ielem), geoel(:,ielem) ,xpq, xpqi, sgmas, ielem)
!
sigmg(1, 1) = -pres + sgmas(1, 1)
sigmg(1, 2) = 0.d0  + sgmas(1, 2)
sigmg(2, 1) = 0.d0  + sgmas(2, 1)
sigmg(2, 2) = -pres + sgmas(2, 2)!
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
!if(itime.gt.1000.and.ielem.eq.19) print*,'domn',ig,fluxd(1,1),gdshp(1, 1)*uadv+gdshp(2, 1)*vadv
!
enddo !...(2)ig = 1,ngausd
!
650 enddo
!
!if(itime.gt.1000) print*,'domn',rhsel(1,1,19)
!
end subroutine rhsdomndg_lag_mc_quadsym
!
!...Symmetry preserving techniques...
!
subroutine barthlimit_sympre_phy(geoel, coord, ustar, unkno, ipqua, intfac, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),          intent(in):: ipqua
integer, dimension(nifai,nafac),               intent(in)::intfac
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp
integer:: ielem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:nq+1, 1:npoin) :: unmax, unmin
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
real*8,dimension(1:nq+1, 1:nsize)   ::aflim
real*8, dimension(1:nq+1, 1:ncell):: unmax_new, unmin_new
!

real*8::eps,c00,c05,c10,c20
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
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
mapmt = 0.d0
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
!...Part 1: Mapping matrix
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
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
!...Get dudx, dudy, dvdx, dvdy
!
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
!
!...LANL
!
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy
!
matra = unkno(1, 2, ielem)
matrb = unkno(1, 3, ielem)
!
if(sqrt(matra**2+matrb**2).gt.1.d-8)then
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
!mapt  = matra + matrd
!mapd  = matra*matrd - matrb*matrc
!lamda1 = 0.5d0*mapt + sqrt(0.25d0*mapt**2-mapd)
!lamda2 = 0.5d0*mapt - sqrt(0.25d0*mapt**2-mapd)
!
!if(abs(matrc).gt.1.d-10)then
!
!lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
!lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
!mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
!mapmt(1, 2, ielem) =          matrc/lmat1
!mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
!mapmt(2, 2, ielem) =          matrc/lmat2
!!
!else
!
!mapmt(1, 1, ielem) = 1.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 1.d0
!endif
!
!
enddo
!
!...Part 2: Find the mapped maxi and mini
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!...loop over quads' vertices..
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
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
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvqua
do iq = 2, 3 !...only velocity...
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
if(unctr(iq).gt.unmax(iq, ipq(iv))) then
unmax(iq, ipq(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
if(unctr(iq).lt.unmin(iq, ipq(iv))) then
unmin(iq, ipq(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
enddo
enddo
!
enddo
!
!
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua, ie)
!
do iq=2, 3
!
unmax_new(iq, ielem) = maxval(unmax(iq, ipq(1:nvqua)))
unmin_new(iq, ielem) = minval(unmin(iq, ipq(1:nvqua)))
enddo
enddo
!
!...Part 3: Impose limiter
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
!
do iv = 1, nvqua
do iq = 2, 3 !...only for vector components...
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
!call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)
!
alfa(iq, iv) = afbar
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 2, 3
aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
enddo
!
enddo
!
!...Part 3: Transform back the limiter...
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
!
afvec(1, 1, ielem) = xicx**2*fixi + etcx**2*fiet
afvec(1, 2, ielem) = xicx*xicy*fixi + etcx*etcy*fiet
afvec(2, 1, ielem) = afvec(1, 2, ielem)
afvec(2, 2, ielem) = xicy**2*fixi + etcy**2*fiet
!!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_sympre_phy
!
!...Find new volume in geoel for curved cell...
!
subroutine getcellvol_lag(iptri, ipqua, geoel, coord)
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
!...Other components are same as those of Eulerian framework...
!
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
!    print*,'old geoel',ielem, geoel(1:5, ielem)
!    print*,'new geoel',ielem, xcent/volel, ycent/volel, volel, maxval(xp(1, 1:3))-minval(xp(1, 1:3)), &
!                         maxval(xp(2, 1:3))- minval(xp(2, 1:3))
!
enddo !...(1)ie = 1,nelem

end subroutine getcellvol_lag
!
!...Find new volume in geoel for linear cell...
!
subroutine getcellvollinear_lag(iptri, ipqua, geoel, coord)
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
real*8:: rm, rp, sm, sp


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
if(ncurv==1)then
xp(1, 1:3) = coord(1, iptri(1:3,ie))
xp(2, 1:3) = coord(2, iptri(1:3,ie))
!
xp(1:2,4) = 0.5d0*(xp(1:2,1)+xp(1:2,2))
xp(1:2,5) = 0.5d0*(xp(1:2,2)+xp(1:2,3))
xp(1:2,6) = 0.5d0*(xp(1:2,1)+xp(1:2,3))
elseif(ncurv==0)then
xp(1, 1:nptri) = coord(1,iptri(1:nptri, ie))
xp(2, 1:nptri) = coord(2,iptri(1:nptri, ie))
endif!
!...
rcv = 0.d0
scv = 0.d0
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
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
volel = volel + djaco
rcv = rcv + r*djaco
scv = scv + s*djaco
!
enddo
!
geoel(3, ielem) = volel
geoel(5, ielem) = rcv/volel
geoel(6, ielem) = scv/volel
!geoel(5, ielem) = geoel(1, ielem)
!geoel(6, ielem) = geoel(2, ielem)
!
!...Other components are same as those of Eulerian framework...
!
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
!...global numbering
!
ielem = ie + ntria
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!
!...
rcv = 0.d0
scv = 0.d0
volel = 0.d0
!
do igaus =1,ngausd_geoq
!
r  = posiq(1,igaus)
s  = posiq(2,igaus)
wi  = weighq(igaus)
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
volel = volel + djaco
rcv = rcv + r*djaco
scv = scv + s*djaco
!
enddo
!
geoel(3, ielem) = volel
geoel(5, ielem) = rcv/volel
geoel(6, ielem) = scv/volel
!geoel(5, ielem) = geoel(1, ielem)
!geoel(6, ielem) = geoel(2, ielem)
!
!    print*,'old geoel',ielem, geoel(1:5, ielem)
!    print*,'new geoel',ielem, xcent/volel, ycent/volel, volel, maxval(xp(1, 1:3))-minval(xp(1, 1:3)), &
!                         maxval(xp(2, 1:3))- minval(xp(2, 1:3))
!
enddo !...(1)ie = 1,nelem

end subroutine getcellvollinear_lag
!
!....rhs updating density for hybrid linear quad cells
!
subroutine rhsdomndg_lag_density(ipqua, coord, coold, geoel, unkno, amatr)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
real*8,dimension(1:nmatr,1:ncell),intent(in)::amatr
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
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
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8::m(ndegr, ndegr)
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds
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
real*8:: rhoi, rhon
real*8:: rcvq,scvq
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
rhsel = 0.d0
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
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(5, ielem) !...mass center...
sc= geoel(6, ielem)
!
!rcvq = 0.d0
!scvq = 0.d0
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
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rhoma = 1.d0/rhon
rhoad = rhon
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg, ielem)=rhsel(ideg, ielem) + rhoi*b(ideg)*djaco
enddo
!
!rcvq = rcvq + r*djaco
!scvq = scvq + s*djaco
!

!
enddo !...(2)ig = 1,ngausd
!
!if(ielem.eq.1)print*,'rcvq, scvq',rcvq/geoel(3,ielem),scvq/geoel(4,ielem)
!
650 enddo
!
!...Update density...
!
if(npoly.ge.1)then !...(2)npoly.ge.1
!
do ie=1,ncell    !...(3)ie=1,nelem
!
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
!
!...step 1
!
!if(ie==2) print*,'rhsel',rhsel(2:3,1,ie),rhsel(2:3,2,ie),rhsel(2:3,3,ie),rhsel(1,4,ie),ie
!
do id =1,ndegr !...(4)id =1,ndegr
!
!...step 1
!
unint = 0.d0
do iunk = 1,ndegr
unint(1) = unint(1) + m(id, iunk)*rhsel(iunk,ie)
enddo
!
!...step 2
!
unkno(id, 1, ie)= unint(1)
!
enddo !...(4)id =1,ndegr
enddo !...(3)ie=1,nelem
!
endif !...(2)npoly.ge.1

end subroutine rhsdomndg_lag_density
!
!....rhs updating density for hybrid linear quad cells
!
subroutine rhsdensity_lag(iptri, ipqua, coord, coold, geoel, unkno, amatr)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
real*8,dimension(1:nmatr,1:ncell),intent(in)::amatr
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
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

end subroutine rhsdensity_lag 
!
!...Get the mass matrix for lagrangian based on volume center(density evolution)...
!
subroutine  getamatr_lagdensity(amatr,geoel,coord, iptri, ipqua)
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
real*8::wi,djaco, volel,det
real*8::c10
real*8::f0,f1,f2,f3,f4
real*8::f5,f6,f7,f8
real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::b2,b3,b4,b5,b6
real*8::radie, radii,radie2,radii2,radic2, radig2,paras,rhoi,rhoe
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
!...volum center...
!
rc= geoel(5, ielem) !...mass center...
sc= geoel(6, ielem)
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
rc= geoel(5, ielem) !...mass center...
sc= geoel(6, ielem)
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
!...Coefficient R of RZ or XY system...
!
rcoef = 1.d0 - alfrz + alfrz*ygaus
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
if(nrz.eq.2) amatr(:, ielem) = amatr(:, ielem)/geoel(11, ielem)
!
enddo !...(2)ie = 1,nelem

end subroutine  getamatr_lagdensity
!
!...Subroutine for barth limiter based on vertex for primitive variables on quads with symmetry preserving....
!
subroutine barthlimit_lagsym_quad2(geoel, coord, coold, ustar, unkno, ipqua, intfac, aflim, afvec, unmax, unmin)
use constant
implicit none
integer, parameter::nvbar=8
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord, coold
real*8,dimension(1:nq+1, 1:nsize),             intent(out)::aflim
integer,  dimension(1:nvqua,1:nquad),          intent(in):: ipqua
integer, dimension(nifai,nafac),               intent(in)::intfac
real*8, dimension(1:nq+1, 1:npoin),           intent(in):: unmax, unmin
real*8, dimension(1:nq+1, 1:ncell):: unmax_new, unmin_new
real*8,dimension(1:2, 1:2, 1:nsize),          intent(out)::afvec
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa
integer:: ielem 
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvbar)::alfa
real*8:: xvq(nvbar), yvq(nvbar)
real*8:: bq(1:ndegr, 1:nvbar),bqv(1:ndegr, 1:nvbar)
real*8:: dunk(1:nq+1)
!real*8:: unmax_new(1:nq+1, 1:nelem), unmin_new(1:nq+1, 1:nelem)
real*8,dimension(1:nq+1,  1:nvbar) ::unknvq
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
xvq(5) = -0.5773502691896257645091488d0; yvq(5) = -0.5773502691896257645091488d0
xvq(6) = -0.5773502691896257645091488d0; yvq(6) =  0.5773502691896257645091488d0
xvq(7) =  0.5773502691896257645091488d0; yvq(7) = -0.5773502691896257645091488d0
xvq(8) =  0.5773502691896257645091488d0; yvq(8) =  0.5773502691896257645091488d0
!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua, ie)
!
do iq=1, nq+1
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
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvbar
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
do iv   = 1,nvbar
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
if(ndens.eq.1)then
rhov = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
!
r = xvq(iv); s = yvq(iv)
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rhov = rhon
!
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
rhov = unknvq(1, iv)
!
endif
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvq(1, iv) = 1.d0/rhov
elseif(ndens.eq.2)then
unknvq(1, iv) = rhov
elseif(ndens.eq.3)then
unknvq(1, iv) = rhov
endif
unknvq(4 ,iv) = pvtx
!
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
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
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq)  = pctr
!
do iv = 1, nvbar
do iq = 1, nq
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!
!if(iq.eq.1)then
!call barthfctrho(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
!else
!call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)
!endif
!
alfa(iq, iv) = afbar
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = 1,nq
aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
enddo
!
!
!aflim(2, ielem) = minval(aflim(2:3, ielem))
!aflim(3, ielem) = aflim(2, ielem)
!
enddo
!
!...Part 2.1: Impose symmetry preserving limiter...
!
call barthlimit_sympre_phy(geoel, coord, ustar, unkno, ipqua, intfac, afvec)
!
!...Degenarate to non-symmetric...
!
do ie = 1, -nquad
!
ielem = ie + ntria
!
afvec(1, 1, ielem) = aflim(2, ielem)
afvec(1, 2, ielem) = 0.d0

afvec(2, 1, ielem) = 0.d0
afvec(2, 2, ielem) = aflim(3, ielem)

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
do iv =1 ,nvbar
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
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
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq) = pctr
unctr(nq+1) = ectr
!

!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1,nvbar
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
if(ndens.eq.1)then
rhov = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
!
r = xvq(iv); s = yvq(iv)
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quad(rhoi, r, s, xpqi)!
call getdensity_quadllnl(r, s, xpq, xpqi, rhoi, rhon)
!
rhov = rhon
!
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
rhov = unknvq(1, iv)
!
endif

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvq(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
elseif(ndens.eq.2)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
elseif(ndens.eq.3)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
unknvq(2, iv) = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
unknvq(3, iv) = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
unknvq(4 ,iv) = pctr + aflim(4, ielem)*(pvtx - pctr)
!
unknvq(5, iv) = unknvq(4 ,iv)/(gamlg-1.d0)*unknvq(1, iv) + 0.5d0*(unknvq(2, iv)**2 + unknvq(3, iv)**2)
!
enddo
!
!
do iv = 1, nvbar
do iq = nq+1, nq+1 !...energy
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!
!call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)

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
enddo
!
aflim = 1.d0
afvec = 0.d0
afvec(1,1,:) = 1.d0
afvec(2,2,:) = 1.d0
!
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lagsym_quad2
!
!...Barth limiter based on Characteristic variables...
!
subroutine barthlimit_lag_shu2(unkno, ipqua, intfac)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
integer, intent(in) :: intfac(nifai,nafac)
!
!...Local
!
integer::ip(nvqua)
integer::indbd(npoin)
integer:: ie, iv, iest, iq, ideg,ifa
real*8:: aflim(1:nq, 1:ncell)
real*8:: unctr(1:nq)
real*8,  dimension(1:nq, 1:nvqua)::alfa
real*8:: xv(4), yv(4)
real*8:: b(3, 1:nvqua)
real*8,  dimension(1:nq, 1:npoin)::unvtx
real*8:: unmax(1:nq, 1:npoin), unmin(1:nq, 1:npoin), dunk(1:nq)
real*8,dimension(1:nq,  1:nvqua) ::unknv
real*8,dimension(1:ndegr,1:nq,1:nsize)::uncha  !...Characteristic variable...
!
real*8::  uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy, ndxc, ndyc, tdxc, tdyc
real*8:: rhom, drhomdr, drhomds, dudr,duds, dvdr, dvds,dedr, deds
real*8:: rhoct, sdctr
!
eps = 1.e-6
!
indbd = 0
!
do ifa = 1, nbfac
indbd(intfac(3:4, ifa)) =1
enddo
!
!
!...shape functions
!
dr = 1.d0
ds = 1.d0
rc = 0.d0
sc = 0.d0
!
xv(1) = -1.d0; yv(1) = -1.d0
xv(2) =  1.d0; yv(2) = -1.d0
xv(3) =  1.d0; yv(3) =  1.d0
xv(4) = -1.d0; yv(4) =  1.d0

do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...Get the successive momentum of characteristic variables...
!
do ie =1 ,nquad
!
rhom    = unkno(1, 1, ie)
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
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
if(sqrt(uctr**2 + vctr**2).gt.1.d13)then
!
ndxc = uctr/sqrt(uctr**2 + vctr**2)
ndyc = vctr/sqrt(uctr**2 + vctr**2)
!
tdxc = -ndyc
tdyc =  ndxc
!
else
!
ndxc =uctr
ndyc =vctr
tdxc =-ndyc
tdyc = ndxc

endif
!
uncha(1, 1, ie) = ectr - (uctr*uctr + vctr*vctr) + pctr*rhom
uncha(2, 1, ie) = dedr - (uctr*dudr + vctr*dvdr) + pctr*drhomdr
uncha(3, 1, ie) = deds - (uctr*duds + vctr*dvds) + pctr*drhomds
!
uncha(1, 2, ie) = rhom    - (uctr*ndxc + vctr*ndyc)/rhoct/sdctr
uncha(2, 2, ie) = drhomdr - (dudr*ndxc + dvdr*ndyc)/rhoct/sdctr
uncha(3, 2, ie) = drhomds - (duds*ndxc + dvds*ndyc)/rhoct/sdctr
!
uncha(1, 3, ie) = rhom    + (uctr*ndxc + vctr*ndyc)/rhoct/sdctr
uncha(2, 3, ie) = drhomdr + (dudr*ndxc + dvdr*ndyc)/rhoct/sdctr
uncha(3, 3, ie) = drhomds + (duds*ndxc + dvds*ndyc)/rhoct/sdctr
!
uncha(1, 4, ie) = (uctr*tdxc + vctr*tdyc)
uncha(2, 4, ie) = (dudr*tdxc + dvdr*tdyc)
uncha(3, 4, ie) = (duds*tdxc + dvds*tdyc)

enddo
!
! if(ie==1) print*,unkno(1, 1:4, ie), unknv(1:4, 1)
!
!...Get the maximum at nodes...
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
do ie = 1, nquad
!
ip(1:nvqua) = ipqua(1:nvqua,ie)
!
unctr(1:nq) = uncha(1, 1:nq, ie)
!
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvqua
do iq = 1, nq
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
if(unctr(iq).gt.unmax(iq, ip(iv))) then
unmax(iq, ip(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmaxpost',iv, iq, unctr(iq),unmax(iq, ip(iv))
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
!
enddo
!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Impose limiter
!
do ie = 1,nquad
!
ip(1:nvqua) = ipqua(1:nvqua,ie)
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + uncha(ideg,1:nq,ie)*b(ideg, iv)
enddo
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
unctr(1:nq) = uncha(1, 1:nq, ie)
!
do iv = 1, nvqua
!
do iq = 1, nq
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1.and.iv==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq,ip(iv)) - unctr(iq)),(unmax(iq,ip(iv)) - unctr(iq))/dunk(iq),&
!                                   (unmin(iq,ip(iv)) - unctr(iq)),&
!                                   (unmin(iq,ip(iv)) - unctr(iq))/dunk(iq)
if(dunk(iq).gt.1.d-13)then

fiy = .8d0*(unmax(iq, ip(iv)) - unctr(iq))/dunk(iq)
alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)

elseif(dunk(iq).lt.-1.d-13)then

fiy = .8d0*(unmin(iq, ip(iv)) - unctr(iq))/dunk(iq)
alfa(iq, iv) = max(min(1.d0, fiy), 0.d0)
!
else
!
alfa(iq, iv) = 1.d0
endif
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==2) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
!if(indbd(ip(iv)).eq.1) alfa(:, iv) = 1.d0
!
enddo
!
!...Get the minimum alfa...
!
do iq = 1,nq
aflim(iq, ie) = minval(alfa(iq, 1:nvqua))
enddo
!
enddo
!
! aflim = .0d0
!
do ie = 1, nquad
!
uncha(2:3, 1, ie) = uncha(2:3, 1, ie)*aflim(1,ie)
uncha(2:3, 2, ie) = uncha(2:3, 2, ie)*aflim(2,ie)
uncha(2:3, 3, ie) = uncha(2:3, 3, ie)*aflim(3,ie)
uncha(2:3, 4, ie) = uncha(2:3, 4, ie)*aflim(4,ie)
!
enddo
!
!...Recover the original variables...
!
do ie =1 ,nquad
!
rhom = unkno(1, 1, ie)
!
uctr  = unkno(1, 2, ie)
!
vctr  = unkno(1, 3, ie)

ectr  = unkno(1, 4, ie)

!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
if(sqrt(uctr**2 + vctr**2).gt.1.d13)then
!
ndxc = uctr/sqrt(uctr**2 + vctr**2)
ndyc = vctr/sqrt(uctr**2 + vctr**2)
!
tdxc = -ndyc
tdyc =  ndxc
!
else
!
!
ndxc =uctr
ndyc =vctr
tdxc =-ndyc
tdyc = ndxc
endif
!!
do ideg = 2, ndegr
!
unkno(ideg, 1, ie) = 0.5d0*(uncha(ideg, 2,ie) + uncha(ideg, 3,ie))

unkno(ideg, 2, ie) = 0.5d0*rhoct*sdctr*(uncha(ideg, 3,ie) - uncha(ideg, 2,ie))*ndxc +&
uncha(ideg, 4,ie)*tdxc

unkno(ideg, 3, ie) = 0.5d0*rhoct*sdctr*(uncha(ideg, 3,ie) - uncha(ideg, 2,ie))*ndyc +&
uncha(ideg, 4,ie)*tdyc

unkno(ideg, 4, ie) = uncha(ideg, 1,ie) + &
0.5d0*rhoct*sdctr*(uncha(ideg, 3,ie) - uncha(ideg, 2,ie))*(uctr*ndxc+vctr*ndyc) + &
uncha(ideg, 4,ie)*(uctr*tdxc+vctr*tdyc) -&
0.5d0*pctr*(uncha(ideg, 2,ie) + uncha(ideg, 3,ie))

!if(ie.eq.1) print*,'ie==1',unkno(ideg, 2, ie),ideg, uncha(ideg,2:4,ie),ndxc,tdxc
enddo
!
enddo
!
! print*,'aflim', aflim(1:4,1)
end subroutine barthlimit_lag_shu2
!
!...Transfer the reference derivatibe to the physical domain
!
subroutine getgradphy_curv(gradr,xpq, r, s)
!
use constant
implicit none
!
real*8,  dimension(1:2, 1:nvqua), intent(in)::xpq
real*8,  dimension(1:2, 1:nq),   intent(inout)::gradr

real*8,  intent(in):: r, s
!
integer:: iq, ishp
!
real*8,  dimension(1:2, 1:nq)::grade
real*8,  dimension(1:nvqua)::shpq,dsprq, dspsq
real*8:: drdx, drdy, dsdx, dsdy, jacom
real*8:: dxdr, dxds, dydr, dyds
real*8:: rm, rp, sm, sp
real*8:: c00, c05, c10, c20, epsil
!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /
!
if(ncurv.eq.0)then
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
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
elseif(ncurv.eq.1)then
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
elseif(ncurv.eq.2)then
 call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
endif
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

!...Get dudx, dudy, dvdx, dvdy
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom

!
do iq = 1, nq
grade(1, iq) = gradr(1, iq)*drdx + gradr(2, iq)*dsdx
grade(2, iq) = gradr(1, iq)*drdy + gradr(2, iq)*dsdy
enddo
!
gradr = grade
!
end subroutine getgradphy_curv
!
!...Symmetry preserving techniques...
!
subroutine barthlimit_sympre_inva(geoel, coord, ustar, unkno, ipqua, intfac, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(inout)::afvec
integer, dimension(1:nvqua,1:nquad),          intent(in):: ipqua
integer, dimension(nifai,nafac),               intent(in)::intfac
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp
integer:: ielem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:nq+1, 1:npoin) :: unmax, unmin
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
real*8,dimension(1:nsize)   ::aflim, afma
real*8,dimension(1:nsize) :: undu
real*8, dimension(1:nq+1, 1:ncell):: unmax_new, unmin_new

!

real*8::eps,c00,c05,c10,c20
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct
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
real*8:: dudx, dudy, dvdx, dvdy
real*8:: dudxy, dumax, dumin
!
real*8:: rhomc,sdctr,delu,volel,macel
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
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
mapmt = 0.d0
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
!...Part 1: Mapping matrix
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
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
!...Get dudx, dudy, dvdx, dvdy
!
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
!
undu(ielem) = dudx + dvdy!
!
enddo
!
!...Part 2: Find the mapped maxi and mini
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!...loop over quads' vertices..
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
unctr(4) = undu(ielem)
!
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvqua
do iq = 4, 4 !...only velocity...
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!

!
if(unctr(iq).gt.unmax(iq, ipq(iv))) then
unmax(iq, ipq(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
if(unctr(iq).lt.unmin(iq, ipq(iv))) then
unmin(iq, ipq(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
enddo
enddo
!
enddo
!
!
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua, ie)
!
do iq=1, nq+1
!
unmax_new(iq, ielem) = maxval(unmax(iq, ipq(1:nvqua)))
unmin_new(iq, ielem) = minval(unmin(iq, ipq(1:nvqua)))
enddo
enddo

!
!...Part 3: Impose limiter
!
do ie = 1, nquad
!
ielem = ie + ntria
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
dudxy = undu(ielem)
dumin = unmin_new(4, ielem)
dumax = unmax_new(4, ielem)
!
!if(ie.eq.6369) print*,'aflim',dumax/dudxy
!
!aflim(ielem) = max(0.d0, min(1.d0,(dumax)/(dudxy)))
!
 if(dumin*dumax.gt.0.d0)then
!
  if(dumax.lt.0.d0)then
    aflim(ielem) = max(0.d0, min(1.d0,(dumax)/(dudxy+eps)))
  elseif(dumin.gt.0.d0)then
    aflim(ielem) = max(0.d0, min(1.d0,(dumin)/(dudxy+eps)))
  endif
!
 elseif(dumin*dumax.le.0.d0)then
!
  aflim(ielem)=0.d0
!
 endif
!
 geoel(10, ielem) = aflim(ielem)
!
! if(ie.eq.2) print*,'aflim',aflim(2),dumax,dumin,eps
!
enddo
!
!print*,'delu lim', aflim(1:60)
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
afma(ielem) = min(1.d0, .0001d0*macel)
!
!if(ielem.eq.2) print*,'bad',afma(ielem),(afma(ielem)),aflim(ielem)
!
aflim(ielem) =(afma(ielem))*aflim(ielem) +  (1.d0-afma(ielem))
!
enddo
!
!print*,'bad',aflim(1:ncell)
!
!aflim = 1.d0
!
!...Part 3: Transform back the limiter...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
afvec(:,:,ielem) = afvec(:,:,ielem)*aflim(ielem)
!!
!!
enddo!
!print*,'aflimccc',afvec(:,:,2),aflim(2),afma(2)
!
end subroutine barthlimit_sympre_inva
!
!...Symmetry preserving techniques...
!
subroutine barthlimit_sympre_inva2(geoel, coord, ustar, unkno, ipqua, intfac, afvec, aflim)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(inout)::afvec
real*8,dimension(1:nq+1, 1:nsize),             intent(inout)::aflim
integer, dimension(1:nvqua,1:nquad),          intent(in):: ipqua
integer, dimension(nifai,nafac),               intent(in)::intfac
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp
integer:: ielem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:nq+1, 1:npoin) :: unmax, unmin
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
real*8,dimension(1:nsize)   ::afdu, afma, afshk
real*8,dimension(1:nsize) :: undu
real*8, dimension(1:nq+1, 1:ncell):: unmax_new, unmin_new

!

real*8::eps,c00,c05,c10,c20
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr
real*8:: rc, sc, dr, ds, fiy
real*8:: rhov, rhoct
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
real*8:: dudx, dudy, dvdx, dvdy
real*8:: dudxy, dumax, dumin, sdmax
!
real*8:: rhomc,sdctr,delu,volel,macel
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
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
mapmt = 0.d0
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
!...Part 1: Mapping matrix
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
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
!...Get dudx, dudy, dvdx, dvdy
!
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
!
undu(ielem) = dudx + dvdy!
!
enddo
!
!...Part 2: Find the mapped maxi and mini
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!...loop over quads' vertices..
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
unctr(4) = undu(ielem)
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.1)then
rhomc = 1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.1)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif

uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhomc !...Cell center density...
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
unctr(5) = sqrt((uctr**2 + vctr**2))/sdctr
!
! if(ie==5244) print*,'unctr1',unctr(5),sdctr,sqrt((uctr**2 + vctr**2))
!
do iv = 1, nvqua
do iq = 4, 5 !...only velocity...
!
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!

!
if(unctr(iq).gt.unmax(iq, ipq(iv))) then
unmax(iq, ipq(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
if(unctr(iq).lt.unmin(iq, ipq(iv))) then
unmin(iq, ipq(iv)) = unctr(iq)
endif
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
enddo
enddo
!
enddo
!
!
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua, ie)
!
do iq=1, nq+1
!
unmax_new(iq, ielem) = maxval(unmax(iq, ipq(1:nvqua)))
unmin_new(iq, ielem) = minval(unmin(iq, ipq(1:nvqua)))
enddo
enddo

!
!...Part 3: Impose limiter
!
do ie = 1, nquad
!
ielem = ie + ntria
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
dudxy = undu(ielem)
dumin = unmin_new(4, ielem)
dumax = unmax_new(4, ielem)
!
!if(ie.eq.6369) print*,'aflim',dumax/dudxy
!
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
!print*,'delu lim', aflim(1:60)
!
!...Further technique...
!
do ie = 1, nquad
!
ielem = ie + ntria
!
if(ndens.eq.1)then
rhomc = unkno(1, 1, ielem)
elseif(ndens.eq.1)then
rhomc = 1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.1)then
rhomc = 1.d0/unkno(1, 1, ielem)
endif
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
afma(ielem) = min(1.d0,.2d0*macel)
!
!afma(ielem)=0.d0
!
!afdu(ielem) =(1.d0-afma(ielem))*afdu(ielem) +  afma(ielem)
 afdu(ielem) =  afma(ielem)*afdu(ielem) +  1.d0-afma(ielem)
!!
!if(ielem.eq.1) print*,'bad',afma(ielem),macel
enddo
!
!...Part 3: Transform back the limiter...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
afvec(:,:,ielem) = afvec(:,:,ielem)*afdu(ielem)
!afvec(1,1,ielem) = afvec(1,1,ielem)*(afdu(ielem)) + (1.d0-afdu(ielem))
!afvec(2,2,ielem) = afvec(2,2,ielem)*(afdu(ielem)) + (1.d0-afdu(ielem))
!afvec(1,2,ielem) = afvec(1,2,ielem)*(afdu(ielem))
!afvec(2,1,ielem) = afvec(2,1,ielem)*(afdu(ielem))
enddo!
!
!...Further technique...
!
do ie = 1, nquad
!
ielem = ie + ntria
!
!print*,'mach',sdmax,ielem
!
sdmax = unmax_new(5, ielem)
!
!if(sdmax.lt.1.d0)print*,'bad',ielem,sdmax
!
if(sdmax.le.0.8d0)then
!
afshk(ielem) = 1.d0

!elseif(sdmax.gt.0.8d0.sdmax.lt.0.85d0)then
!
!afshk(ielem) = 1.d0

elseif(sdmax.gt.0.80d0)then
afshk(ielem) = 0.d0
endif
!!
!if(ielem.eq.1) print*,'bad',afma(ielem),macel
enddo
!
!...Part 3: Transform back the limiter...
!
do ie = 1,-nquad
!
ielem = ie + ntria
!
aflim(1,ielem) = aflim(1,ielem)*(1.d0-afshk(ielem)) + (afshk(ielem))
aflim(4,ielem) = aflim(4,ielem)*(1.d0-afshk(ielem)) + (afshk(ielem))
aflim(5,ielem) = aflim(5,ielem)*(1.d0-afshk(ielem)) + (afshk(ielem))
!
enddo!

!
!...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
!
!afdu(ielem) =(afma(ielem))*afdu(ielem) +  (1.d0-afma(ielem))
!
!aflim(1,ielem) = aflim(1,ielem)*(afma(ielem)) + (1.d0-afma(ielem))
!aflim(4,ielem) = aflim(4,ielem)*(afma(ielem)) + (1.d0-afma(ielem))
!aflim(5,ielem) = aflim(5,ielem)*(afma(ielem)) + (1.d0-afma(ielem))
!
!afvec(1,1,ielem) = afvec(1,1,ielem)*(afma(ielem)) + (1.d0-afma(ielem))
!afvec(2,2,ielem) = afvec(2,2,ielem)*(afma(ielem)) + (1.d0-afma(ielem))
!afvec(1,2,ielem) = afvec(1,2,ielem)*(afma(ielem))
!afvec(2,1,ielem) = afvec(2,1,ielem)*(afma(ielem))
!
enddo!

!
!aflim = 1.d0
!print*,'aflimccc',aflim(1)
!
end subroutine barthlimit_sympre_inva2
!
!...Symmetry preserving techniques...
!
subroutine barthlimit_sympre_phy2(geoel, coord, ustar, unkno, ipqua, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),          intent(in):: ipqua
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
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
mapmt = 0.d0
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
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
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
!...Get dudx, dudy, dvdx, dvdy
!
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
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
!mapmt(1, 1, ielem) = 0.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 0.d0

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
!
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
do iq = 2, 3 !...only for vector components...
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!
alfa(iq, iv) = afbar
!
enddo
!
enddo
!
do iq = 2, 3
aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
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
rhom = unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhom = 1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhom = 1.d0/unkno(1, 1, ielem)
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
afma(ielem) = min(1.d0, .03d0*macel)
!
afdu(ielem) =(afma(ielem))*afdu(ielem) +  (1.d0-afma(ielem))
!if(ielem.eq.1773) print*,'bad',afma(ielem),macel
enddo
!
!aflim = 1.d0
!
!...Part 4: Transform back the limiter...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
!aflim(:,ielem) = aflim(:,ielem)*afdu(ielem)
!
enddo!
!
!...Part 6: Transfer back the limiter to the gobal frame...
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
end subroutine barthlimit_sympre_phy2
!
!...Symmetry preserving techniques for genral test cases...
!
subroutine barthlimit_sympre_phy3(geoel, coord, ustar, unkno, ipqua, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),          intent(in):: ipqua
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
indbd(intfac(3:4, ifa)) = 1
enddo
!
eps = 1.e-6
mapmt = 0.d0
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
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
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
!...Get dudx, dudy, dvdx, dvdy
!
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
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
!mapmt(1, 1, ielem) = 0.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 0.d0

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
!
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
do iq = 2, 3 !...only for vector components...
!
dunk(iq) = unknvq(iq, iv) - unctr(iq)
!
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
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
aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
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
afma(ielem) = min(1.d0, .03d0*macel)
!
afdu(ielem) =(afma(ielem))*afdu(ielem) +  (1.d0-afma(ielem))
!if(ielem.eq.1773) print*,'bad',afma(ielem),macel
enddo
!
!aflim = 1.d0
!
!...Part 4: Transform back the limiter...
!
do ie = 1,nquad
!
ielem = ie + ntria
!
!aflim(:,ielem) = aflim(:,ielem)*afdu(ielem)
!
enddo!
!
!...Part 6: Transfer back the limiter to the gobal frame...
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
end subroutine barthlimit_sympre_phy3
!
!...Subroutine to calcualte the maximum and minimum nodal unknows
!for barth limiter based on vertex for primitive variables and 1st invariant on curved cell....
!
subroutine barthlimit_lagcurv(unkno, iptri, ipqua, unmax, unmin, coord, geoel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
integer*4,dimension(1:nvtri,1:ntria),          intent(in)::iptri
integer,  dimension(1:nvqua,1:nquad),          intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::coord
!
!...Local
!
integer::ipt(nvtri), ipq(nvqua)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp
integer:: ielem
real*8:: unctr(1:nq+2)
real*8:: unmax(1:nq+2, 1:npoin), unmin(1:nq+2, 1:npoin)
!
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rhov, rhoct, rhom
real*8:: afbar
real*8:: c10
real*8:: dudr, duds, dvdr, dvds
real*8:: drdx, drdy, dsdx, dsdy, jacom
real*8:: dudx, dudy, dvdx, dvdy
real*8:: dxdr, dxds, dydr, dyds
real*8:: dr, ds, r, s, rc, sc, rm, rp ,sm ,sp
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8,  dimension(1:2, 1:nvtri)::xpt
real*8,  dimension(1:nvtri)::dspr, dsps
real*8,  dimension(1:mdegr, 1:nvqua)::bq
real*8,  dimension(1:nvqua)::xvq, yvq
real*8,  dimension(1,1:nq,1:nvqua)::unknvq
!
data c10   / 1.0d0 /
!
eps = 1.e-6
!
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0

if(npoly.eq.2)then
xvq(5) =  0.d0; yvq(5) = -1.d0
xvq(6) =  1.d0; yvq(6) =  0.d0
xvq(7) =  0.d0; yvq(7) =  1.d0
xvq(8) = -1.d0; yvq(8) =  0.d0
xvq(9) =  0.d0; yvq(9) =  0.d0
endif
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
!...Get the 1st invariant...
!
dr = 0.5d0
ds = 0.5d0
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
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
dxdr = dxdr + dspr(ishp)*xpt(1,ishp)
dxds = dxds + dsps(ishp)*xpt(1,ishp)

dydr = dydr + dspr(ishp)*xpt(2,ishp)
dyds = dyds + dsps(ishp)*xpt(2,ishp)
enddo
!
!...Get dudx, dudy, dvdx, dvdy
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy

!..1st invariant of strain tensor...
unctr(nq+2) = dudx + dvdy!

!..Characteristic length
unctr(nq+2) = sqrt(geoel(3, ielem))
!
do iv = 1, nvtri
do iq = 1, nq+2
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
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dr = 1.d0
ds = 1.d0
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

!...Basis function
do iv = 1, nvqua
!...Basis function for solutions...
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
!...Get the 1st invariant
!
dr = 1.d0
ds = 1.d0
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
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
!...Get dudx, dudy, dvdx, dvdy
!
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy

!..1st invariant of strain tensor...
unctr(nq+2) = dudx + dvdy!

!..Characteristic length
unctr(nq+2) = sqrt(geoel(3, ielem))

!...zero out unknv
unknvq = 0.d0

do iv = 1, nvqua
 !...Store the interpolation value at the Gauss point.
 if(npoly.eq.2)then
 !...Get the P1 projection at the vertex
 do ideg = 1, mdegr
  unknvq(1, 1:nq, iv) = unknvq(1, 1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
 enddo
 unctr(2) = unknvq(1, 1, iv)
 endif
!
do iq = 1, nq+2
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
end subroutine barthlimit_lagcurv
!
!...Subroutine for barth limiter based on vertex for primitive variables on curv quads....
!
subroutine barthlimit_lag_curvquad(geoel, coord, coold, ustar, unkno, ipqua, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord, coold
real*8,dimension(1:nq+1, 1:nsize),             intent(inout)::aflim
integer,  dimension(1:nvqua,1:nquad),          intent(in):: ipqua
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
!...shape functions
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
!...Store the maximum and minimum values surrounding one cell...
!
do ie = 1, nquad
!
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua, ie)
!
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
if(ndens.eq.1)then
rhov = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
!
r = xvq(iv); s = yvq(iv)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quadcurv(rhoi, xpqi)!
call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)
!
rhov = rhon

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
rhov = unknvq(1, iv)
!
endif
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvq(1, iv) = 1.d0/rhov
elseif(ndens.eq.2)then
unknvq(1, iv) = rhov
elseif(ndens.eq.3)then
unknvq(1, iv) = rhov
endif
unknvq(4 ,iv) = pvtx
!
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
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
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif
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
!if(iq.ne.4)then
call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
!else
!call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)
!endif
!
alfa(iq, iv) = afbar
!
!  if(ie==56) print*,'dunk alfa',iq,iv,alfa(iq,iv),dunk(iq)
!
enddo
!
!...Special treatment of boundary...
!
if(indbd(ipq(iv)).eq.1)then

if(ncase.eq.1)then
if(coord(1, ipq(iv)).lt.1.d-6.or.abs(coord(1, ipq(iv))-1.d0).lt.1.d-6) then
!
!fiy = (- unctr(2))/(dunk(2))
!alfa(2, iv) = min(max((- unctr(2))/(dunk(2)),0.d0),1.d0)
!alfa(2, iv) =  (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)
!
!...Setting as 1..
!
alfa(:, iv) = 1.d0
endif
!
if(coord(2, ipq(iv)).lt.1.d-6.or.abs(coord(2, ipq(iv))-1.d0).lt.1.d-6) then
!
! fiy = ( - unctr(3))/(dunk(3))
!alfa(3, iv) = min(max(( - unctr(3))/(dunk(3)),0.d0),1.d0)
!alfa(3, iv) = (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)
!
!...Setting as 1..
!
alfa(:, iv) = 1.d0
endif
!
endif
endif
!
enddo
!
!...Choose the minimum one
!
do iq = 1,nq
if(iq.ne.4)then
aflim(iq, ielem) = minval(alfa(iq, 1:4))
else
aflim(iq, ielem) = minval(alfa(iq, 1:8))
endif
enddo
!
enddo
!
!print*,'aflim56',aflim(2:3,56)
!
!...Part 2.1: Impose symmetry preserving limiter...
!
call barthlimit_sympre_curv(geoel, coord, ustar, unkno, ipqua, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
!
!...Degenarate to non-symmetric...
!
do ie = 1, -nquad
!
ielem = ie + ntria
!
afvec(1, 1, ielem) = aflim(2, ielem)
afvec(1, 2, ielem) = 0.d0

afvec(2, 1, ielem) = 0.d0
afvec(2, 2, ielem) = aflim(3, ielem)

enddo
!
!call barthlimit_sympre_inva(geoel, coord, ustar, unkno, ipqua, intfac, afvec)
!
!call barthlimit_sympre_inva2(geoel, coord, ustar, unkno, ipqua, intfac, afvec, aflim)
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
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif
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
if(ndens.eq.1)then
rhov = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
!
r = xvq(iv); s = yvq(iv)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quadcurv(rhoi, xpqi)!
call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)
!
rhov = rhon

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
rhov = unknvq(1, iv)
!
endif

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvq(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
elseif(ndens.eq.2)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
elseif(ndens.eq.3)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
unknvq(2, iv) = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
unknvq(3, iv) = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
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
 call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
!call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)

alfa(iq, iv) = afbar
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
enddo
!
!...Define boundary total energy...
!
if(indbd(ipq(iv)).eq.1)then
alfa(5, iv) = 1.d0
endif
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = nq+1,nq+1
aflim(iq, ielem) = minval(alfa(iq, 1:4))
enddo
!
enddo
!
!print*,'aflimccc',unkno(1:3, 3, 18),unkno(1:3, 3, 20)
!
end subroutine barthlimit_lag_curvquad
!
!...Subroutine for barth limiter based on vertex for primitive variables on curvd quads with symmetry preserving....
!
subroutine barthlimit_lag_curvquadcb(geoel, coord, coold, ustar, unkno, ipqua, &
                                     bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
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
!...Part I: Some preparing work
!

!...Coloring the boundary node
indbd = 0  !...indbd represents index of boundary node
do ifa =1 ,nbfac
indbd(intfac(3:(2+nvfac), ifa)) = 1
enddo
!
eps = 1.e-6

!...shape functions
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

!...Store the maximum and minimum values surrounding one cell...
do ie = 1, nquad
 ielem = ie + ntria
 ipq(1:nvqua) = ipqua(1:nvqua, ie)
 do iq=1, nq+2
  unmax_new(iq, ielem) = maxval(unmax(iq, ipq(1:nvqua)))
  unmin_new(iq, ielem) = minval(unmin(iq, ipq(1:nvqua)))
 enddo
enddo
!
!...Part II: Impose limiter
!
do ie = 1, nquad

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

!...zero out unknv
unknvq = 0.d0
do iv   = 1,nvqua
  do ideg = 1,mdegr
   unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
  enddo

 if(ndens.eq.1)then
  rhov = 1.d0/unknvq(1, iv)
 elseif(ndens.eq.2)then
  r = xvq(iv); s = yvq(iv)
!
  xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
  xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
  xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
  xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

  call getrhoig_quadcurv(rhoi, xpqi)!
  call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)

  rhov = rhon
elseif(ndens.eq.3)then
 rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
 bqv(1, iv) = 1.d0
 bqv(2, iv) = (xvq(iv)-rcv)/dr
 bqv(3, iv) = (yvq(iv)-scv)/ds

 unknvq(1, iv) =0.d0
 do ideg = 1,mdegr
 unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
 enddo
 rhov = unknvq(1, iv)
endif

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
 unknvq(1, iv) = 1.d0/rhov
elseif(ndens.eq.2)then
 unknvq(1, iv) = rhov
elseif(ndens.eq.3)then
 unknvq(1, iv) = rhov
endif

unknvq(4 ,iv) = pvtx
enddo !do iv   = 1,nvqua
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
 unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
 unctr(1)   = rhoct
elseif(ndens.eq.3)then
 unctr(1)   = rhoct
endif
 unctr(2:3) = unkno(1, 2:3, ielem)
 unctr(nq)  = pctr
!
do iv = 1, nvqua
 do iq = 1, nq
  dunk(iq) = unknvq(iq, iv) - unctr(iq)
  call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)

  alfa(iq, iv) = afbar
 enddo

!...Special treatment of boundary...
if(indbd(ipq(iv)).eq.1)then

!...Set special boundary
!  alfa(:, iv) = 1.d0
!
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
enddo ! do iv = 1, nvqua

!...Get the minimum value for one cell
do iq = 1,nq
   aflim(iq, ielem) = minval(alfa(iq, 1:4))
enddo
!
enddo
!
!...Part 2.1: Impose symmetry preserving limiter for velocity...
!
call barthlimit_symprecb(geoel, coord, ustar, unkno, ipqua, bface, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
!
!...Degenarate to non-symmetric...
!
do ie = 1, -nquad
!
ielem = ie + ntria
!
afvec(1, 1, ielem) = aflim(2, ielem)
afvec(1, 2, ielem) = 0.d0

afvec(2, 1, ielem) = 0.d0
afvec(2, 2, ielem) = aflim(3, ielem)
enddo
!
!...The limitation of the first invariant...
!
!call barthlimit_sympre_inva(geoel, coord, ustar, unkno, ipqua, intfac, afvec)
!call barthlimit_sympre_inva2(geoel, coord, ustar, unkno, ipqua, intfac, afvec, aflim)
!
!...Part 3: Correct total energy
!
do ie = 1,nquad

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Shape function
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo

!...Cell average of inverse density, velocity and total energy.
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
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif

unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq) = pctr
unctr(nq+1) = ectr

!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
if(ndens.eq.1)then
rhov = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
r = xvq(iv); s = yvq(iv)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quadcurv(rhoi, xpqi)!
call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)
!
rhov = rhon
elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
bqv(1, iv) = 1.d0
bqv(2, iv) = (xvq(iv)-rcv)/dr
bqv(3, iv) = (yvq(iv)-scv)/ds
!
unknvq(1, iv) =0.d0
do ideg = 1,mdegr
unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
enddo

rhov = unknvq(1, iv)
endif !if(ndens.eq.1)then

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvq(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
elseif(ndens.eq.2)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
elseif(ndens.eq.3)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
unknvq(2, iv) = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
unknvq(3, iv) = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
unknvq(4 ,iv) = pctr + aflim(4, ielem)*(pvtx - pctr)
unknvq(5, iv) = unknvq(4 ,iv)/(gamlg-1.d0)*unknvq(1, iv) + 0.5d0*(unknvq(2, iv)**2 + unknvq(3, iv)**2)
!
enddo

!...Get the limiting coefficient for the total energy
do iv = 1, nvqua
 do iq = nq+1, nq+1 !...energy
  dunk(iq) = unknvq(iq, iv) - unctr(iq)
  call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
! call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)

  alfa(iq, iv) = afbar
 enddo

!...Special treatment for the boundary nodes...
! if(indbd(ipq(iv)).eq.1)then
!  alfa(5, iv) = 1.d0
! endif
enddo !do iv = 1, nvqua

!...Get the minimum value for one cell
do iq = nq+1,nq+1
 aflim(iq, ielem) = minval(alfa(iq, 1:4))
enddo

enddo !...Do ie = 1, nquad
!
end subroutine barthlimit_lag_curvquadcb
!
!...Symmetry preserving techniques curved cell including symmetry boundary condition...
!
subroutine barthlimit_symprecb(geoel, coord, ustar, unkno, ipqua, bface, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),           intent(in):: ipqua
integer*4,dimension(1:nbfai,nbfac),           intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
!
!...Local
!
integer:: ipq(nvqua)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp, ivf
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
elseif(bface(3,ifa).eq.25.and.ncase.eq.13)then
indbd(intfac(3:(2+nvfac), ifa)) = 4
endif
enddo

!...Coloring the nodes at y axis
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
!...For basis function scaling...
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
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
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
!...Get dudx, dudy, dvdx, dvdy
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
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
matra = unkno(1, 2, ielem)
matrb = unkno(1, 3, ielem)
if(sqrt(matra**2+matrb**2).gt.1.d-8)then
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
!...Symmetry BC...
!
if(ncase.ne.1)then
if(indbd(ipq(iv)).eq.3)then
!
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

elseif(indbd(ipq(iv)).eq.4)then
!
ucjel = 1.d0
vcjel = 0.d0
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)

elseif(indbd(ipq(iv)).eq.2)then

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

!...Get the max and min around one cell
!do iq = 2, 3
!unmpx(iq) = max(ummax(iq),unmpx(iq))
!unmpn(iq) = min(ummin(iq),unmpn(iq))
!enddo
!
!enddo
!
!...Vector limiting
!do iv = 1, nvqua
do iq = 2, 3
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!call barthfct(unmpx(iq), unmpn(iq), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo

!...Special treatment about the boundary
!if(indbd(ipq(iv)).ne.0)then
!alfa(2:3, iv) = 1.d0
!endif
if(indbd(ipq(iv)).eq.1)then
if(ncase.eq.1)then
 alfa(2:3, iv) = 1.d0
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
afma(ielem) = min(1.d0, .03d0*macel)
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
end subroutine barthlimit_symprecb
!
!...Symmetry preserving techniques curved cell...
!
subroutine barthlimit_sympre_curv(geoel, coord, ustar, unkno, ipqua, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
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
!...For basis function scaling...
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
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
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
!...Get dudx, dudy, dvdx, dvdy
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
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
afma(ielem) = min(1.d0, .3d0*macel)
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
end subroutine barthlimit_sympre_curv
!
!...Symmetry preserving techniques for linear grid(1/4)...
!
subroutine barthlimit_sympre_quadlb(geoel, coord, ustar, unkno, ipqua, bface,intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),          intent(in):: ipqua
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
real*8:: vemag,afmag,signlim
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
if(nrz.eq.1.or.nrz.eq.2)then
afmag = 0.25d0*pi*(1.d0-0.5d0/(ncell/100.d0))
!afmag = 0.125d0*3.d0*pi*(1.d0-0.5d0/3.d0/(ncell/98.d0))
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
!
!
eps = 1.e-6
mapmt = 0.d0
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
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
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
!...Get dudx, dudy, dvdx, dvdy
!
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy

!...1st invariant
undu(ielem) = dudx +dvdy

!...Local velocity
matra = unkno(1, 2, ielem)
matrb = unkno(1, 3, ielem)

!...
!matra = geoel(5, ielem)
!matrb = geoel(6, ielem)
!
if(sqrt(matra**2+matrb**2).gt.1.d-10)then
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

!...LANL Principal strain
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy

!...Maire deformation gradient
!matra = dxdr**2 + dydr**2
!matrb = dxdr*dxds + dydr*dyds
!matrc = matrb
!matrd = dxds**2 + dyds**2

!...eigenvalues of the matrix...
mapt  = matra + matrd
mapd  = matra*matrd - matrb*matrc
lamda1 = 0.5d0*mapt + sqrt(0.25d0*mapt**2-mapd)
lamda2 = 0.5d0*mapt - sqrt(0.25d0*mapt**2-mapd)
!
if(abs(matrc).gt.1.d-8)then
 if(matrc.gt.0.d0)then

 lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
 lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
 mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
 mapmt(1, 2, ielem) =          matrc/lmat1
 mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
 mapmt(2, 2, ielem) =          matrc/lmat2

 elseif(matrc.lt.0.d0)then
!
 lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
 lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)

 mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
 mapmt(2, 2, ielem) =          matrc/lmat1
 mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
 mapmt(1, 2, ielem) =          matrc/lmat2
 endif

else
 mapmt(1, 1, ielem) = 0.d0
 mapmt(1, 2, ielem) = 0.d0
 mapmt(2, 1, ielem) = 0.d0
 mapmt(2, 2, ielem) = 0.d0
endif !if(abs(matrc).gt.1.d-8)then

enddo !...ielem
!
!...Part 4: Impose limiter for mapped u and v
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

!...zero out unknv

 unknvq = 0.d0
 do iv   = 1,nvqua
  do ideg = 1,mdegr
   unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
  enddo
!
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)

!...New mapped velocity components u and v

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
!
!ummax(2:3) = unctr(2:3)
!ummin(2:3) = unctr(2:3)
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

!...Symmetry BC...

if(ncase.ne.1)then
if(indbd(ipq(iv)).eq.3)then

!---1st
if(nrz.eq.0)then
ucjel = unkno(1, 2, jelem)
vcjel = -unkno(1, 3, jelem)
elseif(nrz.eq.1)then
signlim = sign(1.d0,unkno(1, 2, jelem))
vemag = signlim*sqrt(unkno(1, 2, jelem)**2 +unkno(1, 3, jelem)**2 )
ucjel = vemag*cos(afmag)
vcjel = vemag*sin(afmag)

ucjel = unkno(1, 2, jelem)
vcjel = -unkno(1, 3, jelem)
!...Area-weighted only requires theta = 0

elseif(nrz.eq.2)then
ucjel = unkno(1, 2, jelem)
vcjel = -unkno(1, 3, jelem)
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
!enddo
!do iv =1 ,nvqua
!
do iq = 2, 3 !...only for vector components...
 dunk(iq) = unknvq(iq, iv) - unctr(iq)
 call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
 alfa(iq, iv) = afbar
enddo

!...Some special treatment for the coeficient on the boundary...

!if(indbd(ipq(iv)).ne.0)then
!alfa(2:3, iv) = 1.d0
!endif
!
 if(indbd(ipq(iv)).eq.1)then
 if(ncase.eq.1)then
   if(coord(1, ipq(iv)).lt.1.d-6.or.abs(coord(1, ipq(iv))-1.d0).lt.1.d-6) then
     alfa(2:3, iv) = 1.d0
   endif
   if(coord(2, ipq(iv)).lt.1.d-6.or.abs(coord(2, ipq(iv))-1.d0).lt.1.d-6) then
    alfa(2:3, iv) = 1.d0
   endif
 endif
endif
!
enddo !...iv
!
 do iq = 2, 3
  aflim(iq, ielem) = minval(alfa(iq, 1:nvqua))
 enddo
!
enddo  !...ielem
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
afma(ielem) = min(1.d0, .03d0*macel)

!if(ielem.eq.1500) print*,'bad',afma(ielem),afdu(ielem)

afdu(ielem) =(afma(ielem))*afdu(ielem) +  (1.d0-afma(ielem))
geoel(10, ielem) = afdu(ielem)
!
enddo

!
!aflim = 1.d0
!
!...Part 4: Transform back the limiter...
!
!do ie = 1,nquad
! ielem = ie + ntria
! aflim(:,ielem) = aflim(:,ielem)*afdu(ielem)
!enddo!
!
!...Part 6: Transfer back the limiter to the gobal frame...
!
do ie = 1,nquad
ielem = ie + ntria
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
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
enddo
!
end subroutine barthlimit_sympre_quadlb
!
!...Barth limiter based on Characteristic variables...
!
subroutine barthlimit_lag_rieminvrnt(geoel, coord,unkno, ipqua,&
bface, intfac, esuv1, esuv2)
use constant
implicit none

!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
integer*4,dimension(1:nvqua,1:nquad),     intent(in)::ipqua
integer, intent(in) :: intfac(nifai,nafac)
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
integer*4,dimension(1:nbfai,nbfac)::bface

!...Local arrays
integer::ipq(nvqua)
integer::idxbd(npoin)
integer:: ie, iv, iest, iq, ideg,ifa,ielem,ishp,jelem,istor,ivf,idirc
real*8:: aflim(1:nq, 1:ncell)
real*8:: aflimp2(1:nq, 1:ncell)
real*8:: afdirp2(1:2, 1:nq)
real*8:: unctr(1:3, 1:nq)
real*8,  dimension(1:nq, 1:nvqua)::alfa
real*8:: xpq(1:2,1:nvqua)
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(ndegr, 1:nvqua)
real*8,  dimension(1:nq, 1:npoin)::unvtx
real*8:: ummax(1:nq), ummin(1:nq), dunk(1:nq)
real*8,dimension(1:nq,  1:nvqua) ::unknv
real*8,dimension(1:ndegr,1:nq,1:nsize)::uncha  !...Characteristic variable...
real*8,dimension(1:3, 1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8,dimension(1:2, 1:nq) ::gradr
!
real*8:: c10
real*8::  uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy, dnxc, dnyc, dtxc, dtyc
real*8:: rm, sm, rp,sp,r,s
real*8:: dxdr,dxds,dydr,dyds
real*8:: drdx, drdy, dsdx, dsdy
real*8:: rhomc,rhomv
real*8:: dnudx,dnudy,dedx,dedy
real*8:: dudx,dudy,dvdx,dvdy
real*8:: rhomcje,jacom,ecjel,ucjel,vcjel
real*8:: drhomdr, drhomds, drhomdr2, drhomds2, drhomdrs
real*8:: dudr, duds, dudr2, duds2, dudrs
real*8:: dvdr, dvds, dvdr2, dvds2, dvdrs
real*8:: dedr, deds, dedr2, deds2, dedrs
real*8:: dnudje, dudjel, dvdjel, dedjel
real*8:: uloca,vloca,rhomloc,eloca
real*8:: rhoct, sdctr
real*8:: matra,matrb,matrc,matrd,afbar
real*8:: umap, vmap,umpj,vmpj
real*8:: lamda1,lamda2,mapt,mapd,lmat1,lmat2
!
eps = 1.e-6
c10 =1.d0
idxbd = 0

!...Coloring the nodes at x axis

do ifa =1 ,nbfac
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.221)then
idxbd(intfac(3:(2+nvfac), ifa)) = 3
endif
endif
enddo
!
do ifa =1 ,nbfac
!
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.222)then
do ivf = 1, nvfac
if(idxbd(intfac(2+ivf, ifa)).eq.3)then
idxbd(intfac(2+ivf, ifa)) = 10
else
idxbd(intfac(2+ivf, ifa)) = 2
endif
enddo
endif
endif
enddo

!...Vertex in one quad
dr = 1.d0
ds = 1.d0
rc = 0.d0
sc = 0.d0
!
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
xvq(9) =  0.d0; yvq(9) =  0.d0
endif
!
!...Part I: Mapping matrix and store the physical derivative
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...mass center...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem)
sc= geoel(2, ielem)
!
matra = unkno(1, 2, ielem)
matrb = unkno(1, 3, ielem)
!
if(sqrt(matra**2+matrb**2).gt.1.d-8)then
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
!
!...Maire
!
!matra = dxdr**2 + dydr**2
!matrb = dxdr*dxds + dydr*dyds
!matrc = matrb
!matrd = dxds**2 + dyds**2
!
!...LANL
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy
!
!...eigenvalues...
!
mapt  = matra + matrd
mapd  = matra*matrd - matrb*matrc
lamda1 = 0.5d0*mapt + sqrt(0.25d0*mapt**2-mapd)
lamda2 = 0.5d0*mapt - sqrt(0.25d0*mapt**2-mapd)
!
if(abs(matrc).gt.1.d-8)then
if(matrc.gt.0.d0)then

lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)

! mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
! mapmt(1, 2, ielem) =          matrc/lmat1
! mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
! mapmt(2, 2, ielem) =          matrc/lmat2

elseif(matrc.lt.0.d0)then
!
lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)

! mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
! mapmt(2, 2, ielem) =          matrc/lmat1
! mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
! mapmt(1, 2, ielem) =          matrc/lmat2
endif

else
! mapmt(1, 1, ielem) = 0.d0
! mapmt(1, 2, ielem) = 0.d0
! mapmt(2, 1, ielem) = 0.d0
! mapmt(2, 2, ielem) = 0.d0
endif !if(abs(matrc).gt.1.d-8)then
!
enddo
!
!...Part II:Get the characteristic variables...
!
do ie = 1, nquad

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Vertex coordinate
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct = 1.d0/rhomc
pctr = (gamlg-1.d0)*rhoct*(ectr-0.5d0*(uctr**2+vctr**2))
sdctr = sqrt(gamlg*pctr/rhoct)
!
dnxc = mapmt(1, 1, ielem)
dnyc = mapmt(1, 2, ielem)
!
dtxc = mapmt(2, 1, ielem)
dtyc = mapmt(2, 2, ielem)

!...zero out unknv
unknvq = 0.d0

!...Get the mapped varaible at the vertex
do iv   = 1,nvqua
do ideg = 1,3
unknvq(1, 1:nq, iv) = unknvq(1, 1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
rhomv = unknvq(1, 1, iv)
uvtx = unknvq(1, 2, iv)
vvtx = unknvq(1, 3, iv)
evtx = unknvq(1, 4, iv)

!...New mapped velocity components u and v
umap = uctr
vmap = vctr
!
unknvq(1, 1, iv) = evtx - (umap*uvtx + vmap*vvtx) + pctr*rhomv
unknvq(1, 2, iv) = rhomv  - (uvtx*dnxc + vvtx*dnyc)/rhoct/sdctr
unknvq(1, 3, iv) = rhomv  + (uvtx*dnxc + vvtx*dnyc)/rhoct/sdctr
unknvq(1, 4, iv) = (uvtx*dtxc + vvtx*dtyc)
enddo


!...New mapped velocity components u and v
umap = uctr
vmap = vctr

!...Riemann invariant
!...Average
unctr(1, 1) = ectr - (uctr*umap + vctr*vmap) + pctr*rhomc
unctr(1, 2) = rhomc- (uctr*dnxc + vctr*dnyc)/rhoct/sdctr
unctr(1, 3) = rhomc+ (uctr*dnxc + vctr*dnyc)/rhoct/sdctr
unctr(1, 4) = (uctr*dtxc + vctr*dtyc)

!...Average
do iv = 1, nvqua

!...Pre-store some
ummax(1:nq) = unctr(1, 1:nq)
ummin(1:nq) = unctr(1, 1:nq)
!
do  istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
!
rhomcje= unkno(1, 1, jelem)
ucjel = unkno(1, 2, jelem)
vcjel = unkno(1, 3, jelem)
ecjel = unkno(1, 4, jelem)
!
umpj = ucjel
vmpj = vcjel
!
rhomloc= ecjel  - (uctr*umpj + vctr*vmpj) + pctr*rhomcje
uloca = rhomcje - (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
vloca = rhomcje + (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
eloca = (umpj*dtxc + vmpj*dtyc)
!
ummax(1) = max(ummax(1), rhomloc)
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)
ummax(4) = max(ummax(4), eloca)

ummin(1) = min(ummin(1), rhomloc)
ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
ummin(4) = min(ummin(4), eloca)

!...Boundary node
if(idxbd(ipq(iv)).eq.3)then
!
rhomcje= unkno(1, 1, jelem)
ucjel = unkno(1, 2, jelem)
vcjel =-unkno(1, 3, jelem)
ecjel = unkno(1, 4, jelem)
!
umpj = ucjel
vmpj = vcjel
!
rhomloc= ecjel  - (uctr*umpj + vctr*vmpj) + pctr*rhomcje
uloca = rhomcje - (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
vloca = rhomcje + (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
eloca = (umpj*dtxc + vmpj*dtyc)
!
ummax(1) = max(ummax(1), rhomloc)
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)
ummax(4) = max(ummax(4), eloca)

ummin(1) = min(ummin(1), rhomloc)
ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
ummin(4) = min(ummin(4), eloca)

elseif(idxbd(ipq(iv)).eq.2)then
!
rhomcje= unkno(1, 1, jelem)
ucjel =-unkno(1, 2, jelem)
vcjel = unkno(1, 3, jelem)
ecjel = unkno(1, 4, jelem)
!
umpj = ucjel
vmpj = vcjel
!
rhomloc= ecjel  - (uctr*umpj + vctr*vmpj) + pctr*rhomcje
uloca = rhomcje - (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
vloca = rhomcje + (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
eloca = (umpj*dtxc + vmpj*dtyc)
!
ummax(1) = max(ummax(1), rhomloc)
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)
ummax(4) = max(ummax(4), eloca)

ummin(1) = min(ummin(1), rhomloc)
ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
ummin(4) = min(ummin(4), eloca)

elseif(idxbd(ipq(iv)).eq.10)then
!
rhomcje= unkno(1, 1, jelem)
ucjel = unkno(1, 2, jelem)
vcjel =-unkno(1, 3, jelem)
ecjel = unkno(1, 4, jelem)
!
umpj = ucjel
vmpj = vcjel
!
rhomloc= ecjel  - (uctr*umpj + vctr*vmpj) + pctr*rhomcje
uloca = rhomcje - (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
vloca = rhomcje + (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
eloca = (umpj*dtxc + vmpj*dtyc)
!
ummax(1) = max(ummax(1), rhomloc)
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)
ummax(4) = max(ummax(4), eloca)

ummin(1) = min(ummin(1), rhomloc)
ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
ummin(4) = min(ummin(4), eloca)

!
rhomcje= unkno(1, 1, jelem)
ucjel =-unkno(1, 2, jelem)
vcjel = unkno(1, 3, jelem)
ecjel = unkno(1, 4, jelem)
!
umpj = ucjel
vmpj = vcjel
!
rhomloc= ecjel  - (uctr*umpj + vctr*vmpj) + pctr*rhomcje
uloca = rhomcje - (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
vloca = rhomcje + (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
eloca = (umpj*dtxc + vmpj*dtyc)
!
ummax(1) = max(ummax(1), rhomloc)
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)
ummax(4) = max(ummax(4), eloca)

ummin(1) = min(ummin(1), rhomloc)
ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
ummin(4) = min(ummin(4), eloca)

!
rhomcje= unkno(1, 1, jelem)
ucjel =-unkno(1, 2, jelem)
vcjel =-unkno(1, 3, jelem)
ecjel = unkno(1, 4, jelem)
!
umpj = ucjel
vmpj = vcjel
!
rhomloc= ecjel  - (uctr*umpj + vctr*vmpj) + pctr*rhomcje
uloca = rhomcje - (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
vloca = rhomcje + (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
eloca = (umpj*dtxc + vmpj*dtyc)
!
ummax(1) = max(ummax(1), rhomloc)
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)
ummax(4) = max(ummax(4), eloca)

ummin(1) = min(ummin(1), rhomloc)
ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
ummin(4) = min(ummin(4), eloca)

endif
enddo
!
do iq = 1, nq !...only for vector components...
dunk(iq) = unknvq(1, iq, iv) - unctr(1, iq)
call barthfct(ummax(iq), ummin(iq), unctr(1, iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo
!
enddo !...iv

!...Get the coefficient...
do iq = 1, nq
aflim(iq, ielem) = minval(alfa(iq, 1:4))
enddo
!
enddo  !...ielem
!
!...Get the successive momentum of characteristic variables...
!
do ie =1 ,nquad
!
ielem = ie + ntria
!
rhomc    = unkno(1, 1, ie)
drhomdr  = unkno(2, 1, ie)
drhomds  = unkno(3, 1, ie)
!
uctr  = unkno(1, 2, ie)
dudr  = unkno(2, 2, ie)
duds  = unkno(3, 2, ie)
!
vctr  = unkno(1, 3, ie)
dvdr  = unkno(2, 3, ie)
dvds  = unkno(3, 3, ie)
!
ectr  = unkno(1, 4, ie)
dedr  = unkno(2, 4, ie)
deds  = unkno(3, 4, ie)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
dnxc = mapmt(1, 1, ielem)
dnyc = mapmt(1, 2, ielem)
!
dtxc = mapmt(2, 1, ielem)
dtyc = mapmt(2, 2, ielem)
!
uncha(1, 1, ie) = ectr  - (uctr*uctr  + vctr*vctr ) + pctr*rhomc
uncha(2, 1, ie) = dedr  - (uctr*dudr  + vctr*dvdr ) + pctr*drhomdr
uncha(3, 1, ie) = deds  - (uctr*duds  + vctr*dvds ) + pctr*drhomds
!
uncha(1, 2, ie) = rhomc    - (uctr*dnxc  + vctr*dnyc )/rhoct/sdctr
uncha(2, 2, ie) = drhomdr  - (dudr*dnxc  + dvdr*dnyc )/rhoct/sdctr
uncha(3, 2, ie) = drhomds  - (duds*dnxc  + dvds*dnyc )/rhoct/sdctr
!
uncha(1, 3, ie) = rhomc    + (uctr*dnxc  + vctr*dnyc)/rhoct/sdctr
uncha(2, 3, ie) = drhomdr  + (dudr*dnxc  + dvdr*dnyc)/rhoct/sdctr
uncha(3, 3, ie) = drhomds  + (duds*dnxc  + dvds*dnyc)/rhoct/sdctr
!
uncha(1, 4, ie) = (uctr*dtxc  + vctr*dtyc)
uncha(2, 4, ie) = (dudr*dtxc  + dvdr*dtyc)
uncha(3, 4, ie) = (duds*dtxc  + dvds*dtyc)

enddo
!
!...Part III: Project back the global system
!
do ie = 1, nquad
!
uncha(2:3, 1, ie) = uncha(2:3, 1, ie)*aflim(1,ie)
uncha(2:3, 2, ie) = uncha(2:3, 2, ie)*aflim(2,ie)
uncha(2:3, 3, ie) = uncha(2:3, 3, ie)*aflim(3,ie)
uncha(2:3, 4, ie) = uncha(2:3, 4, ie)*aflim(4,ie)
!
enddo
!
!print*,'sod before',unkno(6,1,399),uncha(6, 2,399) + uncha(6, 3,399)
!
!...Recover the original variables...
!
do ie =1 ,nquad
!
ielem = ie + ntria
!
if(geoel(10, ielem).lt.10.d0) cycle
!
rhomc = unkno(1, 1, ie)
uctr  = unkno(1, 2, ie)
vctr  = unkno(1, 3, ie)
ectr  = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
dnxc = mapmt(1, 1, ielem)
dnyc = mapmt(1, 2, ielem)
!
dtxc = mapmt(2, 1, ielem)
dtyc = mapmt(2, 2, ielem)
!
do ideg = 2, ndegr
unkno(ideg, 1, ie) = 0.5d0*(uncha(ideg, 2,ie) + uncha(ideg, 3,ie))
unkno(ideg, 2, ie) = 0.5d0*rhoct*sdctr*(uncha(ideg, 3,ie) - uncha(ideg, 2,ie))*dnxc +&
uncha(ideg, 4,ie)*dtxc
unkno(ideg, 3, ie) = 0.5d0*rhoct*sdctr*(uncha(ideg, 3,ie) - uncha(ideg, 2,ie))*dnyc +&
uncha(ideg, 4,ie)*dtyc

unkno(ideg, 4, ie) = uncha(ideg, 1,ie) + &
0.5d0*rhoct*sdctr*(uncha(ideg, 3,ie) - uncha(ideg,2,ie))*(uctr*dnxc+vctr*dnyc) + &
uncha(ideg, 4,ie)*(uctr*dtxc+vctr*dtyc) -&
0.5d0*pctr*(uncha(ideg, 2,ie) + uncha(ideg, 3,ie))
enddo
!
if(npoly.gt.1) unkno(4:ndegr, :, ie) = 0.d0
!
enddo
!
!print*,'sod after',unkno(6,1,399)
!
end subroutine barthlimit_lag_rieminvrnt
!
!...Barth limiter based on Characteristic variables...
!
subroutine barthlimit_lag_riemaninv(geoel, coord,unkno, ipqua,&
bface, intfac, esuv1, esuv2)
use constant
implicit none

!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
integer*4,dimension(1:nvqua,1:nquad),     intent(in)::ipqua
integer, intent(in) :: intfac(nifai,nafac)
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
integer*4,dimension(1:nbfai,nbfac)::bface

!...Local arrays
integer::ipq(nvqua)
integer::indbd(npoin)
integer:: ie, iv, iest, iq, ideg,ifa,ielem,ishp,jelem,istor,ivf,idirc
real*8:: aflim(1:nq, 1:ncell)
real*8:: aflimp2(1:nq, 1:ncell)
real*8:: afdirp2(1:2, 1:nq)
real*8:: unctr(1:3, 1:nq)
real*8,  dimension(1:nq, 1:nvqua)::alfa
real*8:: xpq(1:2,1:nvqua)
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(ndegr, 1:nvqua)
real*8,  dimension(1:nq, 1:npoin)::unvtx
real*8:: ummax(1:nq), ummin(1:nq), dunk(1:nq)
real*8,dimension(1:nq,  1:nvqua) ::unknv
real*8,dimension(1:ndegr,1:nq,1:nsize)::uncha  !...Characteristic variable...
real*8,dimension(1:3, 1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8, dimension(2,1:nq,1:ncell) :: unkphy
real*8, dimension(1:3,1:2,1:nq) :: unkp2
real*8,dimension(1:2, 1:nq) ::gradr
!
real*8:: c10
real*8::  uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy, dnxc, dnyc, dtxc, dtyc
real*8:: rm, sm, rp,sp,r,s
real*8:: dxdr,dxds,dydr,dyds
real*8:: drdx, drdy, dsdx, dsdy
real*8:: rhomc,rhomv
real*8:: dnudx,dnudy,dedx,dedy
real*8:: dudx,dudy,dvdx,dvdy
real*8:: rhomcje,jacom,ecjel,ucjel,vcjel
real*8:: drhomdr, drhomds, drhomdr2, drhomds2, drhomdrs
real*8:: dudr, duds, dudr2, duds2, dudrs
real*8:: dvdr, dvds, dvdr2, dvds2, dvdrs
real*8:: dedr, deds, dedr2, deds2, dedrs
real*8:: dnudje, dudjel, dvdjel, dedjel
real*8:: uloca,vloca,rhomloc,eloca
real*8:: rhoct, sdctr
real*8:: matra,matrb,matrc,matrd,afbar
real*8:: umap, vmap,umpj,vmpj
real*8:: lamda1,lamda2,mapt,mapd,lmat1,lmat2
!
eps = 1.e-6
c10 =1.d0
indbd = 0
!
unkphy = 0.d0

!...Coloring the nodes at x axis
if(ncase.ne.1)then
do ifa =1 ,nbfac

if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.221)then
indbd(intfac(3:(2+nvfac), ifa)) = 3
endif
endif
enddo
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
indbd(intfac(3:(2+nvfac), ifa)) = 1
enddo
endif

!...Vertex in one quad
dr = 1.d0
ds = 1.d0
rc = 0.d0
sc = 0.d0
!
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
xvq(9) =  0.d0; yvq(9) =  0.d0
endif
!
!...Part I: Mapping matrix and store the physical derivative
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...mass center...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem)
sc= geoel(2, ielem)
!
unkphy(1, 1:nq, ielem) = unkno(2, 1:nq, ielem)/dr
unkphy(2, 1:nq, ielem) = unkno(3, 1:nq, ielem)/ds

!...Get the phycial derivative
call getgradphy_curv(unkphy(1:2, 1:nq, ielem),xpq, rc, sc)

matra = unkno(1, 2, ielem)
matrb = unkno(1, 3, ielem)
!
if(sqrt(matra**2+matrb**2).gt.1.d-8)then
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
!
!...Maire
!
!matra = dxdr**2 + dydr**2
!matrb = dxdr*dxds + dydr*dyds
!matrc = matrb
!matrd = dxds**2 + dyds**2
!
!...LANL
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy
!
!...eigenvalues...
!
mapt  = matra + matrd
mapd  = matra*matrd - matrb*matrc
lamda1 = 0.5d0*mapt + sqrt(0.25d0*mapt**2-mapd)
lamda2 = 0.5d0*mapt - sqrt(0.25d0*mapt**2-mapd)
!
if(abs(matrc).gt.1.d-8)then
if(matrc.gt.0.d0)then

lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)

! mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
! mapmt(1, 2, ielem) =          matrc/lmat1
! mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
! mapmt(2, 2, ielem) =          matrc/lmat2

elseif(matrc.lt.0.d0)then
!
lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)

! mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
! mapmt(2, 2, ielem) =          matrc/lmat1
! mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
! mapmt(1, 2, ielem) =          matrc/lmat2
endif

else
! mapmt(1, 1, ielem) = 0.d0
! mapmt(1, 2, ielem) = 0.d0
! mapmt(2, 1, ielem) = 0.d0
! mapmt(2, 2, ielem) = 0.d0
endif !if(abs(matrc).gt.1.d-8)then
!
enddo
!
!...Part II:Get the characteristic variables...
!
do ie = 1, nquad

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Vertex coordinate
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
!
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct = 1.d0/rhomc
pctr = (gamlg-1.d0)*rhoct*(ectr-0.5d0*(uctr**2+vctr**2))
sdctr = sqrt(gamlg*pctr/rhoct)
!
dnxc = mapmt(1, 1, ielem)
dnyc = mapmt(1, 2, ielem)
!
dtxc = mapmt(2, 1, ielem)
dtyc = mapmt(2, 2, ielem)

!...zero out unknv
unknvq = 0.d0

!...Get the mapped varaible at the vertex
do iv   = 1,nvqua
do ideg = 1,3
unknvq(1, 1:nq, iv) = unknvq(1, 1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
rhomv = unknvq(1, 1, iv)
uvtx = unknvq(1, 2, iv)
vvtx = unknvq(1, 3, iv)
evtx = unknvq(1, 4, iv)

!...New mapped velocity components u and v
umap = uctr
vmap = vctr
!
unknvq(1, 1, iv) = evtx - (umap*uvtx + vmap*vvtx) + pctr*rhomv
unknvq(1, 2, iv) = rhomv  - (uvtx*dnxc + vvtx*dnyc)/rhoct/sdctr
unknvq(1, 3, iv) = rhomv  + (uvtx*dnxc + vvtx*dnyc)/rhoct/sdctr
unknvq(1, 4, iv) = (uvtx*dtxc + vvtx*dtyc)
enddo


!...Get the variable derivative at the vertex
if(npoly.eq.2)then

!...direction r
unkp2(1, 1, 1:nq) =  unkno(2,1:nq,ielem)/dr
unkp2(2, 1, 1:nq) =  unkno(4,1:nq,ielem)/dr
unkp2(3, 1, 1:nq) =  unkno(6,1:nq,ielem)/dr

!...direction s
unkp2(1, 2, 1:nq) =  unkno(3,1:nq,ielem)/ds
unkp2(2, 2, 1:nq) =  unkno(6,1:nq,ielem)/ds
unkp2(3, 2, 1:nq) =  unkno(5,1:nq,ielem)/ds

!...
do iv   = 1,nvqua
 do ideg = 1,3
  unknvq(2, 1:nq, iv) = unknvq(2, 1:nq, iv) + unkp2(ideg, 1, 1:nq)*bq(ideg, iv)
  unknvq(3, 1:nq, iv) = unknvq(3, 1:nq, iv) + unkp2(ideg, 2, 1:nq)*bq(ideg, iv)
 enddo

!
gradr(1, 1:nq) = unknvq(2, 1:nq, iv)
gradr(2, 1:nq) = unknvq(3, 1:nq, iv)
!...Get the phycial derivative
call getgradphy_curv(gradr,xpq, xvq(iv), yvq(iv))

!
unknvq(2, 1:nq, iv) =  gradr(1, 1:nq)
unknvq(3, 1:nq, iv) =  gradr(2, 1:nq)
!
dnudx = unknvq(2, 1, iv)
dnudy = unknvq(3, 1, iv)

dudx = unknvq(2, 2, iv)
dudy = unknvq(3, 2, iv)

dvdx = unknvq(2, 3, iv)
dvdy = unknvq(3, 3, iv)

dedx = unknvq(2, 4, iv)
dedy = unknvq(3, 4, iv)

!...The physical derivative of the Riemann invariant
unknvq(2, 1, iv) = dedx - (dudx*uctr + dvdx*vctr) + pctr*dnudx
unknvq(2, 2, iv) = dnudx  - (dudx*dnxc + dvdx*dnyc)/rhoct/sdctr
unknvq(2, 3, iv) = dnudx  + (dudx*dnxc + dvdx*dnyc)/rhoct/sdctr
unknvq(2, 4, iv) = (dudx*dtxc + dvdx*dtyc)

unknvq(3, 1, iv) = dedy - (dudy*uctr + dvdy*vctr) + pctr*dnudy
unknvq(3, 2, iv) = dnudy  - (dudy*dnxc + dvdy*dnyc)/rhoct/sdctr
unknvq(3, 3, iv) = dnudy  + (dudy*dnxc + dvdy*dnyc)/rhoct/sdctr
unknvq(3, 4, iv) = (dudy*dtxc + dvdy*dtyc)
enddo
endif

!...New mapped velocity components u and v
umap = uctr
vmap = vctr

!...Riemann invariant
!...Average
unctr(1, 1) = ectr - (uctr*umap + vctr*vmap) + pctr*rhomc
unctr(1, 2) = rhomc- (uctr*dnxc + vctr*dnyc)/rhoct/sdctr
unctr(1, 3) = rhomc+ (uctr*dnxc + vctr*dnyc)/rhoct/sdctr
unctr(1, 4) = (uctr*dtxc + vctr*dtyc)

!...Derivative
dnudx = unkphy(1, 1, ielem)
dnudy = unkphy(2, 1, ielem)

dudx = unkphy(1, 2, ielem)
dudy = unkphy(2, 2, ielem)

dvdx = unkphy(1, 3, ielem)
dvdy = unkphy(2, 3, ielem)

dedx = unkphy(1, 4, ielem)
dedy = unkphy(2, 4, ielem)
!
unctr(2, 1) = dedx - (dudx*umap + dvdx*vmap) + pctr*dnudx
unctr(2, 2) = dnudx  - (dudx*dnxc + dvdx*dnyc)/rhoct/sdctr
unctr(2, 3) = dnudx  + (dudx*dnxc + dvdx*dnyc)/rhoct/sdctr
unctr(2, 4) = (dudx*dtxc + dvdx*dtyc)

unctr(3, 1) = dedy - (dudy*umap + dvdy*vmap) + pctr*dnudy
unctr(3, 2) = dnudy  - (dudy*dnxc + dvdy*dnyc)/rhoct/sdctr
unctr(3, 3) = dnudy  + (dudy*dnxc + dvdy*dnyc)/rhoct/sdctr
unctr(3, 4) = (dudy*dtxc + dvdy*dtyc)

!...Gthe limiting coefficient

if(npoly.eq.2)then
!...Derivative
do idirc = 1, 2
!...Loop over every node
do iv = 1, nvqua

!...Pre-store some
ummax(1:nq) = unctr(1+idirc, 1:nq)
ummin(1:nq) = unctr(1+idirc, 1:nq)
!
do  istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
!
dnudje = unkphy(idirc, 1, jelem)
dudjel = unkphy(idirc, 2, jelem)
dvdjel = unkphy(idirc, 3, jelem)
dedjel = unkphy(idirc, 4, jelem)
!
rhomloc= dedjel - (umap*dudjel + vmap*dvdjel) + pctr*dnudje
uloca = dnudje - (dudjel*dnxc + dvdjel*dnyc)/rhoct/sdctr
vloca = dnudje + (dudjel*dnxc + dvdjel*dnyc)/rhoct/sdctr
eloca = (dudjel*dtxc + dvdjel*dtyc)
!
ummax(1) = max(ummax(1), rhomloc)
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)
ummax(4) = max(ummax(4), eloca)

ummin(1) = min(ummin(1), rhomloc)
ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
ummin(4) = min(ummin(4), eloca)
enddo
!
do iq = 1, nq !...only for vector components...
dunk(iq) = unknvq(1+idirc,iq, iv) - unctr(1+idirc, iq)
call barthfct(ummax(iq), ummin(iq), unctr(1+idirc, iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo
!
enddo !...iv

!...Get the limiting coefficient for every direction
 do iq = 1, nq
 afdirp2(idirc, iq) = minval(alfa(iq, 1:4))
 enddo
enddo !do idirc = 1, 2

!...Get the smallest coefficient forthird derivative...
do iq = 1, nq
aflimp2(iq, ielem) = minval(afdirp2(:, iq))
enddo

endif

!...Average
do iv = 1, nvqua

!...Pre-store some
ummax(1:nq) = unctr(1, 1:nq)
ummin(1:nq) = unctr(1, 1:nq)
!
do  istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
!
rhomcje= unkno(1, 1, jelem)
ucjel = unkno(1, 2, jelem)
vcjel = unkno(1, 3, jelem)
ecjel = unkno(1, 4, jelem)
!
umpj = ucjel
vmpj = vcjel
!
rhomloc= ecjel  - (uctr*umpj + vctr*vmpj) + pctr*rhomcje
uloca = rhomcje - (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
vloca = rhomcje + (umpj*dnxc + vmpj*dnyc)/rhoct/sdctr
eloca = (umpj*dtxc + vmpj*dtyc)
!
ummax(1) = max(ummax(1), rhomloc)
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)
ummax(4) = max(ummax(4), eloca)

ummin(1) = min(ummin(1), rhomloc)
ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
ummin(4) = min(ummin(4), eloca)
enddo
!
do iq = 1, nq !...only for vector components...
dunk(iq) = unknvq(1, iq, iv) - unctr(1, iq)
call barthfct(ummax(iq), ummin(iq), unctr(1, iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo
!
enddo !...iv

!...Get the coefficient...
do iq = 1, nq
aflim(iq, ielem) = minval(alfa(iq, 1:4))
enddo

!...Get the final  one
if(npoly.eq.2)then
do iq = 1, nq
aflim(iq, ielem) = max(aflim(iq, ielem), aflimp2(iq, ielem))
enddo
endif
!
enddo  !...ielem
!
!...Get the successive momentum of characteristic variables...
!
do ie =1 ,nquad
!
ielem = ie + ntria
!
rhomc    = unkno(1, 1, ie)
drhomdr  = unkno(2, 1, ie)
drhomds  = unkno(3, 1, ie)
!
uctr  = unkno(1, 2, ie)
dudr  = unkno(2, 2, ie)
duds  = unkno(3, 2, ie)
!
vctr  = unkno(1, 3, ie)
dvdr  = unkno(2, 3, ie)
dvds  = unkno(3, 3, ie)
!
ectr  = unkno(1, 4, ie)
dedr  = unkno(2, 4, ie)
deds  = unkno(3, 4, ie)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
dnxc = mapmt(1, 1, ielem)
dnyc = mapmt(1, 2, ielem)
!
dtxc = mapmt(2, 1, ielem)
dtyc = mapmt(2, 2, ielem)
!
uncha(1, 1, ie) = ectr  - (uctr*uctr  + vctr*vctr ) + pctr*rhomc
uncha(2, 1, ie) = dedr  - (uctr*dudr  + vctr*dvdr ) + pctr*drhomdr
uncha(3, 1, ie) = deds  - (uctr*duds  + vctr*dvds ) + pctr*drhomds
!
uncha(1, 2, ie) = rhomc    - (uctr*dnxc  + vctr*dnyc )/rhoct/sdctr
uncha(2, 2, ie) = drhomdr  - (dudr*dnxc  + dvdr*dnyc )/rhoct/sdctr
uncha(3, 2, ie) = drhomds  - (duds*dnxc  + dvds*dnyc )/rhoct/sdctr
!
uncha(1, 3, ie) = rhomc    + (uctr*dnxc  + vctr*dnyc)/rhoct/sdctr
uncha(2, 3, ie) = drhomdr  + (dudr*dnxc  + dvdr*dnyc)/rhoct/sdctr
uncha(3, 3, ie) = drhomds  + (duds*dnxc  + dvds*dnyc)/rhoct/sdctr
!
uncha(1, 4, ie) = (uctr*dtxc  + vctr*dtyc)
uncha(2, 4, ie) = (dudr*dtxc  + dvdr*dtyc)
uncha(3, 4, ie) = (duds*dtxc  + dvds*dtyc)

!
if(npoly.eq.2)then
drhomdr2 = unkno(4, 1, ielem)
drhomds2 = unkno(5, 1, ielem)
drhomdrs = unkno(6, 1, ielem)
!
dudr2 = unkno(4, 2, ielem)
duds2 = unkno(5, 2, ielem)
dudrs = unkno(6, 2, ielem)
!
dvdr2 = unkno(4, 3, ielem)
dvds2 = unkno(5, 3, ielem)
dvdrs = unkno(6, 3, ielem)
!
dedr2 = unkno(4, 4, ielem)
deds2 = unkno(5, 4, ielem)
dedrs = unkno(6, 4, ielem)
!
uncha(4, 1, ie) = dedr2 - (uctr*dudr2 + vctr*dvdr2) + pctr*drhomdr2
uncha(5, 1, ie) = deds2 - (uctr*duds2 + vctr*dvds2) + pctr*drhomds2
uncha(6, 1, ie) = dedrs - (uctr*dudrs + vctr*dvdrs) + pctr*drhomdrs
!
uncha(4, 2, ie) = drhomdr2 - (dudr2*dnxc + dvdr2*dnyc)/rhoct/sdctr
uncha(5, 2, ie) = drhomds2 - (duds2*dnxc + dvds2*dnyc)/rhoct/sdctr
uncha(6, 2, ie) = drhomdrs - (dudrs*dnxc + dvdrs*dnyc)/rhoct/sdctr
!
uncha(4, 3, ie) = drhomdr2 + (dudr2*dnxc + dvdr2*dnyc)/rhoct/sdctr
uncha(5, 3, ie) = drhomds2 + (duds2*dnxc + dvds2*dnyc)/rhoct/sdctr
uncha(6, 3, ie) = drhomdrs + (dudrs*dnxc + dvdrs*dnyc)/rhoct/sdctr
!
uncha(4, 4, ie) = (dudr2*dtxc + dvdr2*dtyc)
uncha(5, 4, ie) = (duds2*dtxc + dvds2*dtyc)
uncha(6, 4, ie) = (dudrs*dtxc + dvdrs*dtyc)
!
!if(ie.eq.1)then
!print*,'sod cha',
!endif


endif
enddo
!
!...Part III: Project back the global system
!
do ie = 1, nquad
!
uncha(2:3, 1, ie) = uncha(2:3, 1, ie)*aflim(1,ie)
uncha(2:3, 2, ie) = uncha(2:3, 2, ie)*aflim(2,ie)
uncha(2:3, 3, ie) = uncha(2:3, 3, ie)*aflim(3,ie)
uncha(2:3, 4, ie) = uncha(2:3, 4, ie)*aflim(4,ie)
!
if(npoly.eq.2)then
uncha(4:6, 1, ie) = uncha(4:6, 1, ie)*aflimp2(1,ie)
uncha(4:6, 2, ie) = uncha(4:6, 2, ie)*aflimp2(2,ie)
uncha(4:6, 3, ie) = uncha(4:6, 3, ie)*aflimp2(3,ie)
uncha(4:6, 4, ie) = uncha(4:6, 4, ie)*aflimp2(4,ie)
endif
!
enddo
!
!print*,'sod before',unkno(6,1,399),uncha(6, 2,399) + uncha(6, 3,399)
!
!...Recover the original variables...
!
do ie =1 ,nquad
!
ielem = ie + ntria
!
if(geoel(10, ielem).lt.10.d0) cycle
!
rhomc = unkno(1, 1, ie)
uctr  = unkno(1, 2, ie)
vctr  = unkno(1, 3, ie)
ectr  = unkno(1, 4, ie)
!
rhoct  = 1.d0/rhomc
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
dnxc = mapmt(1, 1, ielem)
dnyc = mapmt(1, 2, ielem)
!
dtxc = mapmt(2, 1, ielem)
dtyc = mapmt(2, 2, ielem)
!
do ideg = 2, ndegr
unkno(ideg, 1, ie) = 0.5d0*(uncha(ideg, 2,ie) + uncha(ideg, 3,ie))
unkno(ideg, 2, ie) = 0.5d0*rhoct*sdctr*(uncha(ideg, 3,ie) - uncha(ideg, 2,ie))*dnxc +&
uncha(ideg, 4,ie)*dtxc
unkno(ideg, 3, ie) = 0.5d0*rhoct*sdctr*(uncha(ideg, 3,ie) - uncha(ideg, 2,ie))*dnyc +&
uncha(ideg, 4,ie)*dtyc

unkno(ideg, 4, ie) = uncha(ideg, 1,ie) + &
0.5d0*rhoct*sdctr*(uncha(ideg, 3,ie) - uncha(ideg,2,ie))*(uctr*dnxc+vctr*dnyc) + &
uncha(ideg, 4,ie)*(uctr*dtxc+vctr*dtyc) -&
0.5d0*pctr*(uncha(ideg, 2,ie) + uncha(ideg, 3,ie))
enddo
enddo
!
!print*,'sod after',unkno(6,1,399)
!
end subroutine barthlimit_lag_riemaninv
!
!...Smooth indicator (Kim)...
!
subroutine getSIHO_kim(geoel, coord,unkno, ipqua, unmax, unmin)
use constant
implicit none

!...Input arrays
real*8,dimension(1:ngeel,1:nsize),    intent(inout) ::geoel
real*8,dimension(1:ndimn,1:npoin),       intent(in) ::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
integer*4,dimension(1:nvqua,1:nquad),     intent(in)::ipqua
real*8, dimension(1:nq+2, 1:npoin),       intent(in):: unmax, unmin

!...Local arrays
integer::ipq(nvqua)
integer::indpt(4)
integer:: ie, iv, iest, iq, ideg,ifa,ielem,ivf
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(ndegr, 1:nvqua)
real*8,dimension(1:3, 1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt

!
real*8:: c10
real*8::  uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds
real*8:: rhomc,rhomv
real*8:: rhoct, sdctr
real*8:: matra,matrb,matrc,matrd,afbar
real*8:: kxtrem, lchara
real*8:: rhomvmax, rhomvmin
real*8:: rhomhvmax, rhomhvmin, rhomvp1
real*8:: pnfilt,pnslop
!
eps = 1.e-6
c10 =1.d0
!
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
xvq(9) =  0.d0; yvq(9) =  0.d0
endif
!
dr =1.d0
ds =1.d0
!
!...Part II:Get the characteristic variables...
!
do ie = 1, nquad

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Characteric length and coefficient
kxtrem = 10.d0

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

!...Mass average
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct = 1.d0/rhomc
pctr = (gamlg-1.d0)*rhoct*(ectr-0.5d0*(uctr**2+vctr**2))
sdctr = sqrt(gamlg*pctr/rhoct)

!...zero out indpt
indpt = 0

!...zero out unknv
unknvq = 0.d0

!...Get the P1 projection at the vertex
do iv   = 1, 4!nvqua
do ideg = 1, mdegr
unknvq(1, 1:nq, iv) = unknvq(1, 1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
rhomv = unknvq(1, 1, iv)
uvtx = unknvq(1, 2, iv)
vvtx = unknvq(1, 3, iv)
evtx = unknvq(1, 4, iv)

!...Maximum and minimum rhom around one vertex
rhomvmax = unmax(1, ipq(iv))
rhomvmin = unmin(1, ipq(iv))
!...Maximum and minimum rhomh around one vertex
rhomhvmax = unmax(2, ipq(iv))
rhomhvmin = unmin(2, ipq(iv))

!...P1-projected MLP condition
if(rhomhvmin.ge.rhomvmin.and.rhomhvmax.le.rhomvmax)then
indpt(iv) = 1
endif

!...Detect the local extrema
if(indpt(iv).eq.0)then

!...Local characteristic length
lchara = unmax(nq+2, ipq(iv))

!...Local extrema detection
!if((rhomvmax - rhomvmin).le.kxtrem*lchara**2)then
!indpt(iv) = 1
!endif

!...The rhomv with Pn
rhomv = 0.d0
do ideg = 1,mdegr
rhomv = rhomv + unkno(ideg,1,ielem)*bq(ideg, iv)
enddo

!...The rhomv with projected P1
rhomvp1 = 0.d0
do ideg = 1, 3
rhomvp1 = rhomvp1 + unkno(ideg,1,ielem)*bq(ideg, iv)
enddo

pnslop = rhomvp1 - rhomc
pnfilt = rhomv - rhomvp1

!if(abs(rhomv-rhomc).le.max(1.d-3*rhomc, geoel(3, ielem)))then
if(abs(rhomv-rhomc).le.max(1.d-3*rhomc, 1.d-10))then
 indpt(iv) = 1
else

!...Local maximum
if(pnslop.gt.0.d0.and.pnfilt.lt.0.d0.and.rhomv.gt.rhomvmin)then
 indpt(iv) = 1
endif
!...Local minimum
 if(pnslop.lt.0.d0.and.pnfilt.gt.0.d0.and.rhomv.lt.rhomvmax)then
  indpt(iv) = 1
 endif
endif
!
endif
!
enddo !...iv

!...Mark the troubled cell...
if(minval(indpt(1:4)).eq.1)then
geoel(10, ielem) = 0.d0
else
geoel(10, ielem) = 100.d0 !...Means troubled cell
endif
!
enddo  !...ielem
!
end subroutine getSIHO_kim
!
!...Smooth indicator (Kim) from augmented MLP...
!
subroutine getSIHOamlp_kim(geoel, coord,unkno, ipqua, unmax, unmin)
use constant
implicit none

!...Input arrays
real*8,dimension(1:ngeel,1:nsize),    intent(inout) ::geoel
real*8,dimension(1:ndimn,1:npoin),       intent(in) ::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
integer*4,dimension(1:nvqua,1:nquad),     intent(in)::ipqua
real*8, dimension(1:nq+2, 1:npoin),       intent(in):: unmax, unmin

!...Local arrays
integer::ipq(nvqua)
integer::indpt(4)
integer:: ie, iv, iest, iq, ideg,ifa,ielem,ivf
real*8:: xvq(nvqua), yvq(nvqua)
real*8:: bq(ndegr, 1:nvqua)
real*8,dimension(1:3, 1:nq+1,  1:nvqua) ::unknvq

!
real*8:: c10
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds
real*8:: rhomc,rhomv
real*8:: rhoct, sdctr
real*8:: kxtrem, lchara
real*8:: rhomvmax, rhomvmin
real*8:: rhomhvmax, rhomhvmin, rhomvp1
real*8:: pnfilt,pnslop
real*8:: drhom
!
eps = 1.e-6
c10 =1.d0
!
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
xvq(9) =  0.d0; yvq(9) =  0.d0
endif
!
dr =1.d0
ds =1.d0

!...Characteric length and coefficient
kxtrem = 10.d0
!
!...Part II:Get the characteristic variables...
!
do ie = 1, nquad

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

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

!...Mass average
rhomc = unkno(1, 1, ielem)
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct = 1.d0/rhomc
pctr = (gamlg-1.d0)*rhoct*(ectr-0.5d0*(uctr**2+vctr**2))
sdctr = sqrt(gamlg*pctr/rhoct)

!...zero out indpt
indpt = 0

!...zero out unknv
unknvq = 0.d0

!...Get the interpolations at the vertex
do iv   = 1, 4
do ideg = 1, mdegr
unknvq(1, 1:nq, iv) = unknvq(1, 1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
rhomv = unknvq(1, 1, iv)

!...Maximum and minimum averaged rhom around one vertex
rhomvmax = unmax(1, ipq(iv))
rhomvmin = unmin(1, ipq(iv))

!...Maximum and minimum interpolation rhomh around one vertex
rhomhvmax = unmax(2, ipq(iv))
rhomhvmin = unmin(2, ipq(iv))
!
drhom = max(1d-4*(rhomvmax-rhomvmin),1d-3)

!... MLP condition
if((rhomhvmin+drhom).ge.rhomvmin.and.rhomhvmax.le.(rhomvmax+drhom))then
!if((rhomv).ge.rhomvmin.and.rhomv.le.(rhomvmax))then
indpt(iv) = 1
endif

!
enddo !...iv

!...Mark the troubled cell...
if(minval(indpt(1:4)).eq.1)then
geoel(10, ielem) = 0.d0
else
geoel(10, ielem) = 100.d0 !...Means troubled cell
endif
!
enddo  !...ielem
!
end subroutine getSIHOamlp_kim
!
!...high order subroutine for barth limiter on curved quads with symmetry preserving....
!
subroutine barthlimit_lag_curvquadcb_HO(geoel, coord, coold, ustar, unkno, ipqua, &
bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
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

!...Store the maximum and minimum values surrounding one cell...
do ie = 1, nquad
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)
do iq=1, nq+2
unmax_new(iq, ielem) = maxval(unmax(iq, ipq(1:nvqua)))
unmin_new(iq, ielem) = minval(unmin(iq, ipq(1:nvqua)))
enddo
enddo
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
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo

if(ndens.eq.1)then
rhov = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
r = xvq(iv); s = yvq(iv)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

call getrhoig_quadcurv(rhoi, xpqi)!
call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)

rhov = rhon
elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
bqv(1, iv) = 1.d0
bqv(2, iv) = (xvq(iv)-rcv)/dr
bqv(3, iv) = (yvq(iv)-scv)/ds

unknvq(1, iv) =0.d0
do ideg = 1,3 !mdegr
unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
enddo
rhov = unknvq(1, iv)
endif

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvq(1, iv) = 1.d0/rhov
elseif(ndens.eq.2)then
unknvq(1, iv) = rhov
elseif(ndens.eq.3)then
unknvq(1, iv) = rhov
endif

unknvq(4 ,iv) = pvtx
enddo !do iv   = 1,nvqua
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
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq)  = pctr

do iv = 1, nvqua

!...
do iq = 1, nq
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)

alfa(iq, iv) = afbar
enddo

!...Special treatment of boundary...
if(indpt(ipq(iv)).eq.1)then
!...Set special boundary
!  alfa(:, iv) = 1.d0

if(ncase.eq.1)then
if(coord(1, ipq(iv)).lt.1.d-6.or.abs(coord(1, ipq(iv))-1.d0).lt.1.d-6) then
alfa(:, iv) = 1.d0
endif
if(coord(2, ipq(iv)).lt.1.d-6.or.abs(coord(2, ipq(iv))-1.d0).lt.1.d-6) then
alfa(:, iv) = 1.d0
endif
endif
endif
!
enddo ! do iv = 1, nvqua

!...Get the minimum value for one cell
!...Exclude pressure
do iq = 1,nq
if(geoel(10, ielem).gt.10.d0)then
aflim(iq, ielem) = minval(alfa(iq, 1:4))
else
aflim(iq, ielem) = 1.d0
endif
enddo
!
enddo
!
!...Part 2.1: Impose symmetry preserving limiter for velocity...
!
!if(ncurv.le.1)then
 call barthlimit_symprecb_HO(geoel, coord, ustar, unkno, ipqua, bface, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
!else
! call      barthlimit_vector(geoel, coord, ustar, unkno, ipqua, bface, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
!endif
!
!...Part 3: Correct total energy
!
do ie = 1,nquad

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Skip the smooth cell
if(geoel(10, ielem).lt.10.d0) cycle

!...Shape function
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo

!...Cell average of inverse density, velocity and total energy.
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
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif

unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq) = pctr
unctr(nq+1) = ectr

!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1, 3!mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
if(ndens.eq.1)then
rhov = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
r = xvq(iv); s = yvq(iv)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quadcurv(rhoi, xpqi)!
call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)
!
rhov = rhon
elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
bqv(1, iv) = 1.d0
bqv(2, iv) = (xvq(iv)-rcv)/dr
bqv(3, iv) = (yvq(iv)-scv)/ds
!
unknvq(1, iv) =0.d0
do ideg = 1, 3!mdegr
unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
enddo

rhov = unknvq(1, iv)
endif !if(ndens.eq.1)then

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvq(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
elseif(ndens.eq.2)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
elseif(ndens.eq.3)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
!unknvq(2, iv) = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
!unknvq(3, iv) = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
unknvq(4 ,iv) = pctr + aflim(4, ielem)*(pvtx - pctr)
!unknvq(4 ,iv) = pvtx
unknvq(5, iv) = unknvq(4 ,iv)/(gamlg-1.d0)*unknvq(1, iv) + 0.5d0*(unknvq(2, iv)**2 + unknvq(3, iv)**2)
!
enddo

!...Get the limiting coefficient for the total energy
do iv = 1, nvqua
do iq = nq+1, nq+1 !...energy
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
! call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)

alfa(iq, iv) = afbar
enddo

!if(ie.eq.1)print*,'Energy limiter',iv, unmax(nq+1, ipq(iv)), unmin(nq+1, ipq(iv)), unctr(nq+1), dunk(nq+1), afbar

!...Special treatment for the boundary nodes...
! if(indpt(ipq(iv)).eq.1)then
!alfa(5, iv) = 1.d0
! endif
enddo !do iv = 1, nvqua


!if(ielem.eq.89)print*,'limiter',ielem,unmax(nq+1, ipq(3)), unmin(nq+1, ipq(3))

!...Get the minimum value for one cell
do iq = nq+1,nq+1
if(geoel(10, ielem).gt.10.d0)then
aflim(iq, ielem) = minval(alfa(iq, 1:4))
else
aflim(iq, ielem) = 1.d0
endif
enddo
enddo !...Do ie = 1, nquad
!
!...Part 4: Annihilate the HO moments of DGp2
!
do ie = 1,nquad

ielem = ie + ntria

!...Bad cell
 if(geoel(10, ielem).gt.10.d0)then
  unkno(4:6, :, ielem) = 0.d0
 else
  aflim(:, ielem) = 1.d0

  afvec(1, 1, ielem) = 1.d0
  afvec(1, 2, ielem) = 0.d0
  afvec(2, 1, ielem) = 0.d0
  afvec(2, 2, ielem) = 1.d0
 endif

enddo

end subroutine barthlimit_lag_curvquadcb_HO
!
!...High-order symmetry preserving techniques curved cell including symmetry boundary condition...
!
subroutine barthlimit_symprecb_HO(geoel, coord, ustar, unkno, ipqua, bface, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),           intent(in):: ipqua
integer*4,dimension(1:nbfai,nbfac),           intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
!
!...Local
!
integer:: ipq(nvqua)
integer:: indpt(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp, ivf
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
indpt = 0  !...indpt represents index of boundary node
!
!...Coloring the nodes at x axis
!
if(ncase.ne.1)then

do ifa =1 ,nbfac
!
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.221)then
indpt(intfac(3:(2+nvfac), ifa)) = 3
endif
endif
enddo

!...Coloring the nodes at y axis
do ifa =1 ,nbfac
!
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.222)then
do ivf = 1, nvfac
if(indpt(intfac(2+ivf, ifa)).eq.3)then
indpt(intfac(2+ivf, ifa)) = 10
else
indpt(intfac(2+ivf, ifa)) = 2
endif
enddo
endif
endif
enddo

else

do ifa =1 ,nbfac
!
indpt(intfac(3:(2+nvfac), ifa)) = 1
enddo
endif
!
eps = 1.e-6
mapmt = 0.d0
!
!...For basis function scaling...
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
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
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
!...Get dudx, dudy, dvdx, dvdy
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
!
!...LANL
!
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy
!
matra = unkno(1, 2, ielem)
matrb = unkno(1, 3, ielem)
if(sqrt(matra**2+matrb**2).gt.1.d-8)then
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
do ideg = 1,3!mdegr
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
!...Symmetry BC...
!
if(ncase.ne.1)then
if(indpt(ipq(iv)).eq.3)then
!
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

elseif(indpt(ipq(iv)).eq.2)then

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

elseif(indpt(ipq(iv)).eq.10)then
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

!...Get the max and min around one cell
!do iq = 2, 3
!unmpx(iq) = max(ummax(iq),unmpx(iq))
!unmpn(iq) = min(ummin(iq),unmpn(iq))
!enddo
!
!enddo
!
!...Vector limiting
!do iv = 1, nvqua
do iq = 2, 3
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!call barthfct(unmpx(iq), unmpn(iq), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo

!...Special treatment about the boundary
!if(indpt(ipq(iv)).ne.0)then
!alfa(2:3, iv) = 1.d0
!endif
if(indpt(ipq(iv)).eq.1)then
if(ncase.eq.1)then
alfa(2:3, iv) = 1.d0
endif
endif
!
enddo
!
do iq = 2, 3
if(geoel(10, ielem).gt.10.d0)then
aflim(iq, ielem) = minval(alfa(iq, 1:4)) !...Only consider the 4 vertices
else
aflim(iq, ielem) = 1.d0
endif
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
if(geoel(10, ielem).gt.10.d0)then
afvec(1, 1, ielem) = fi11
afvec(1, 2, ielem) = fi12
afvec(2, 1, ielem) = afvec(1, 2, ielem)
afvec(2, 2, ielem) = fi22
else
afvec(1, 1, ielem) = 1.d0
afvec(1, 2, ielem) = 0.d0
afvec(2, 1, ielem) = 0.d0
afvec(2, 2, ielem) = 1.d0
endif
!!
enddo
!
end subroutine barthlimit_symprecb_HO
!
!...Vector limiting for high-order ...
!
subroutine barthlimit_vector(geoel, coord, ustar, unkno, ipqua, bface, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),           intent(in):: ipqua
integer*4,dimension(1:nbfai,nbfac),           intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
!
!...Local
!
integer:: ipq(nvqua)
integer:: idxpt(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp, ivf
integer:: ielem, jelem, istor
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvqua)::alfa
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::shpq,dsprq, dspsq
real*8,  dimension(1:nvqua):: rvq, svq
real*8:: bq(1:ndegr, 1:nvqua)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvqua) ::unknvq
real*8, dimension(1:nq+1) :: ummax, ummin
real*8, dimension(1:2, 1:2,1:ncell) :: mapmt
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
!...Part I: Preliminary seutp
!

!...idxpt represents index of boundary node
idxpt = 0

!...Coloring the nodes at x axis
if(ncase.ne.1)then

do ifa =1 ,nbfac
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.221)then
idxpt(intfac(3:(2+nvfac), ifa)) = 3
endif
endif
enddo
!...Coloring the nodes at y axis
do ifa =1 ,nbfac
!
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.222)then
do ivf = 1, nvfac
if(idxpt(intfac(2+ivf, ifa)).eq.3)then
idxpt(intfac(2+ivf, ifa)) = 10
else
idxpt(intfac(2+ivf, ifa)) = 2
endif
enddo
endif
endif
enddo

else

do ifa =1 ,nbfac
idxpt(intfac(3:(2+nvfac), ifa)) = 1
enddo
endif
!
eps = 1.e-6
mapmt = 0.d0

!...Basis function scaling...
dr = 1.d0
ds = 1.d0

!...Vertex reference coordinates
rvq(1) = -1.d0; svq(1) = -1.d0
rvq(2) =  1.d0; svq(2) = -1.d0
rvq(3) =  1.d0; svq(3) =  1.d0
rvq(4) = -1.d0; svq(4) =  1.d0


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
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc

!...Get the points for FEM
if(ncurv==2) call getcoord_fe(ncurv, nvfac, nvqua, xpq)

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
!...Get dudx, dudy, dvdx, dvdy
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
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
!mapmt(1, 1, ielem) = 0.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 0.d0
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
lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
mapmt(2, 2, ielem) =          matrc/lmat1
mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
mapmt(1, 2, ielem) =          matrc/lmat2
endif

else
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
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,4
bq(1, iv) = 1.d0
bq(2, iv) = (rvq(iv)-rc)/dr
bq(3, iv) = (svq(iv)-sc)/ds
enddo

!...zero out unknv
unknvq = 0.d0
do iv   = 1,4
do ideg = 1,3!mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)

!...New mapped velocity components u and v
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
!
unmpx  = 1.d-10
unmpn  = 1.d10
!
do iv = 1, 4
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
!...Symmetry BC...
!
if(ncase.ne.1)then
if(idxpt(ipq(iv)).eq.3)then
!
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

elseif(idxpt(ipq(iv)).eq.2)then

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

elseif(idxpt(ipq(iv)).eq.10)then
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

!...Get the max and min around one cell
!do iq = 2, 3
!unmpx(iq) = max(ummax(iq),unmpx(iq))
!unmpn(iq) = min(ummin(iq),unmpn(iq))
!enddo
!
!enddo
!
!...Vector limiting
!do iv = 1, nvqua
do iq = 2, 3
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!call barthfct(unmpx(iq), unmpn(iq), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo

!...Special treatment about the boundary
!if(indpt(ipq(iv)).ne.0)then
!alfa(2:3, iv) = 1.d0
!endif
if(idxpt(ipq(iv)).eq.1)then
if(ncase.eq.1)then
alfa(2:3, iv) = 1.d0
endif
endif
!
enddo
!
do iq = 2, 3
if(geoel(10, ielem).gt.10.d0)then
aflim(iq, ielem) = minval(alfa(iq, 1:4)) !...Only consider the 4 vertices
else
aflim(iq, ielem) = 1.d0
endif
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
if(geoel(10, ielem).gt.10.d0)then
afvec(1, 1, ielem) = fi11
afvec(1, 2, ielem) = fi12
afvec(2, 1, ielem) = afvec(1, 2, ielem)
afvec(2, 2, ielem) = fi22
else
afvec(1, 1, ielem) = 1.d0
afvec(1, 2, ielem) = 0.d0
afvec(2, 1, ielem) = 0.d0
afvec(2, 2, ielem) = 1.d0
endif
!!
enddo
!
end subroutine barthlimit_vector
!
!....Smooth indicator from LLNL
!
subroutine getSIHO_vorticity(intfac, ipqua, coord, coold, geoel, unkno)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(inout)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua

!...Local integer
integer::ie,ig,ideg,ishp,iv,ielem

!...local integer array
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndegr):: b
real*8, dimension(1:nq):: unknod,unkp1,unvar1,unvar2
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)

!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8:: drdx, drdy, dsdx, dsdy
real*8:: dudx,dudy,dvdx,dvdy
real*8::dr,ds,rc,sc, rcv, scv
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::djaco, wi, smthv, jacom
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

!
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds

!....Zero out unvar2, unvar1
!...Numerator and denumerator of the smooth indicator

unvar2 = 0.d0
unvar1 = 0.d0

!...Gauss loop
do ig = 1,1!ngausdq !...(2)ig = 1,ngausd
!
r  = rc!posiq(1,ig)
s  = sc!posiq(2,ig)
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

!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy


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
!
smthv = abs(dudx+dvdy)/(abs(dudx*dvdy-dvdx*dudy)+1d-6)
!
!unvar1(1:nq) = unvar1(1:nq) + (unkno(1,1:nq,ielem))**2*djaco
!
enddo !...(2)ig = 1,ngausd

!...Coarse smooth indicator
if(log10(smthv).le.1.d0)then
!if(log10(unvar2(1)/unvar1(1)).le.(-1.5d0*log10(2.d0**4)-4.d0))then
geoel(10, ielem) = 0.d0
print*,'vorticity', ielem,log10(smthv)
endif
!
!if(ielem.eq.2100) print*,'ielem',log10(unvar2(1)/unvar1(1)),unvar2(1),unvar1(1)
!
650 enddo

end subroutine getSIHO_vorticity
!
!....Smooth indicator from machine learning
!
subroutine getSIHO_general_ml(intfac, ipqua, coord, coold, geoel, unkno)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(inout)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua

!...Local integer
integer::ie,ig,ideg,ishp,iv,ielem, imt, jmt

!...local integer array
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8, dimension(1:nq):: unknod,unkp1,unvar1,unvar2
real*8, dimension(ngausdq):: unvarg
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)

real*8, dimension(8)::outp0
real*8, dimension(6)::outp1
real*8, dimension(4)::outp2
real*8, dimension(1)::outp3

!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::djaco, wi
!
!
data eps   / 1.0d-6 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /

!...Give gaussian position and weight...
call ruqope(2, ngausdq, posiq, weighq)
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
if(ncurv.eq.2) call getcoord_fe(ncurv, nvfac, nvqua, xpq)

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

!....Zero out unvar2, unvar1
!...Numerator and denumerator of the smooth indicator

unvar2 = 0.d0
unvar1 = 0.d0

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
djaco = wi*(dxdr*dyds - dydr*dxds)

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

!...Record the sepcific vlaue at quadrature point
unvarg(ig) = unknod(1)
!

enddo !...(2)ig = 1,ngausd

!...ANN

!...1st layer
outp0 = 0.d0
!
do jmt = 1, 8
do imt = 1, 25
outp0(jmt) = outp0(jmt) + weigh0t(jmt, imt)*unvarg(imt)
enddo
enddo
!
outp0(:) = 1.d0/(1.d0+exp(-outp0(:)))

!...2nd layer
outp1 = 0.d0
!
do jmt = 1, 6
do imt = 1, 8
outp1(jmt) = outp1(jmt) + weigh1t(jmt, imt)*outp0(imt)
enddo
enddo

!
outp1(:) = 1.d0/(1.d0+exp(-outp1(:)))

!...3rd layer
outp2 = 0.d0
!
do jmt = 1, 4
do imt = 1, 6
outp2(jmt) = outp2(jmt) + weigh2t(jmt, imt)*outp1(imt)
enddo
enddo

!
outp2(:) = 1.d0/(1.d0+exp(-outp2(:)))

!...3rd layer
outp3 = 0.d0
!
do jmt = 1, 1
do imt = 1, 4
outp3(jmt) = outp3(jmt) + weigh3t(jmt, imt)*outp2(imt)
enddo
enddo

!
outp3(:) = 1.d0/(1.d0+exp(-outp3(:)))

!...Record
geoel(10, ielem) = outp3(1)

!
650 enddo

end subroutine getSIHO_general_ml
!
!....Smooth indicator from Persson
!
subroutine getSIHO_general(intfac, ipqua, coord, coold, geoel, unkno)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(inout)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua

!...Local integer
integer::ie,ig,ideg,ishp,iv,ielem

!...local integer array
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8, dimension(1:nq):: unknod,unkp1,unvar1,unvar2
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)

!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::djaco, wi
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

!....Zero out unvar2, unvar1
!...Numerator and denumerator of the smooth indicator

unvar2 = 0.d0
unvar1 = 0.d0

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
djaco = wi*(dxdr*dyds - dydr*dxds)

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

!...Solution at P1
unkp1 = 0.d0

if(npoly.eq.1)then
do ideg =1, 1!
unkp1(1:nq) = unkp1(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo
elseif(npoly.eq.2)then
do ideg =1, 3!
unkp1(1:nq) = unkp1(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo
endif

unvar2(1:nq) = unvar2(1:nq) + (unkp1(1:nq)-unknod(1:nq))**2*djaco
unvar1(1:nq) = unvar1(1:nq) + (unknod(1:nq))**2*djaco
!
!if(ielem.eq.2100) print*,'ielemig',ig,unkno(:,1,ielem)

!
!unvar1(1:nq) = unvar1(1:nq) + (unkno(1,1:nq,ielem))**2*djaco
!
enddo !...(2)ig = 1,ngausd

!...Coarse smooth indicator
geoel(10, ielem) = 0.d0
if(log10(unvar2(1)/(unvar1(1))).le.(-3.d0*log10(4.d0)-4.5d0))then
!if(log10(unvar2(1)/unvar1(1)).le.(-1.5d0*log10(2.d0**4)-4.d0))then
geoel(10, ielem) = 0.d0
else
geoel(10, ielem) = 100.d0
endif
!
!if(ielem.eq.2100) print*,'ielem',log10(unvar2(1)/unvar1(1)),unvar2(1),unvar1(1)
!
650 enddo

end subroutine getSIHO_general
!
!....Smooth indicator from Persson
!
subroutine getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(inout)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua

!...Local integer
integer::ie,ig,ideg,ishp,iv,ielem

!...local integer array
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds, bv
real*8, dimension(1:nq):: unknod,unkp1,unvar1,unvar2
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)

!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::djaco, wi
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
!...Loop over quad
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

!...Points consitituting one element...
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

!....Zero out unvar2, unvar1
!...Numerator and denumerator of the smooth indicator

unvar2 = 0.d0
unvar1 = 0.d0

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
djaco = wi*(dxdr*dyds - dydr*dxds)

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

!...Solution at P1
unkp1 = 0.d0

if(npoly.eq.1)then
 do ideg =1, 1!
 unkp1(1:nq) = unkp1(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
 enddo
elseif(npoly.eq.2)then
 do ideg =1, 3!
 unkp1(1:nq) = unkp1(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
 enddo
endif

unvar2(1:nq) = unvar2(1:nq) + (unkp1(1:nq)-unknod(1:nq))**2*djaco
unvar1(1:nq) = unvar1(1:nq) + (unknod(1:nq))**2*djaco
!
enddo !...(2)ig = 1,ngausd

!...Coarse smooth indicator
geoel(10, ielem) = 0.d0
if(log10(unvar2(1)/(unvar1(1))).le.(-3.d0*log10(4.d0)-4.d0))then
!if(log10(unvar2(1)/unvar1(1)).le.(-1.5d0*log10(2.d0**4)-4.d0))then
geoel(10, ielem) = 0.d0
else
geoel(10, ielem) = 100.d0
endif
!
!if(ielem.eq.1489) print*,'ielem',log10(unvar2(1)/unvar1(1))
!
650 enddo

end subroutine getSIHO_Persson
!
!....Smooth indicator from Persson for the deformation gradient
!
subroutine getSIHOgd_Persson(intfac, ipqua, coord, coold, geoel, unkno,unkgd)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(inout)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua

!...Local integer
integer::ie,ig,ideg,ishp,iv,ielem

!...local integer array
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:ndegr):: b, bi, dbdr, dbds, bv
real*8, dimension(1:nq):: unknod,unkp1,unvar1,unvar2
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8, dimension(1: ndimn, 1:ndimn)::jacbf


!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv,rci,sci
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::djaco, wi
real*8:: a11, a12, a21, a22
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
!...Loop over quad
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

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

!...Geometry parameters for reference cell...
dr = 1.d0
ds = 1.d0
!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...The initial cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)

!...The derivatives of basis function...
!...Here dbdx means dbd(xsi), dbdy means dbd(eta)
dbdr(1)= 0.d0
dbdr(2)= 1.d0/dr
dbdr(3)= 0.d0

dbds(1)= 0.d0
dbds(2)= 0.d0
dbds(3)= 1.0/ds

!....Zero out unvar2, unvar1
!...Numerator and denumerator of the smooth indicator

unvar2 = 0.d0
unvar1 = 0.d0

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

!...Gauss points...
xg = r
yg = s

!
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
do ideg = 1, 1!ndegr
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
djaco = wi*(a11*a22-a12*a21)
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

!...Solution at P1
unkp1 = 0.d0

if(npoly.eq.1)then
do ideg =1, 1!
unkp1(1:nq) = unkp1(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo
elseif(npoly.eq.2)then
do ideg =1, 3!
unkp1(1:nq) = unkp1(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo
endif

unvar2(1:nq) = unvar2(1:nq) + (unkp1(1:nq)-unknod(1:nq))**2*djaco
unvar1(1:nq) = unvar1(1:nq) + (unknod(1:nq))**2*djaco
!
enddo !...(2)ig = 1,ngausd

!...Coarse smooth indicator
geoel(10, ielem) = 0.d0
if(log10(unvar2(1)/(unvar1(1))).le.(-3.d0*log10(4.d0)-4.d0))then
!if(log10(unvar2(1)/unvar1(1)).le.(-1.5d0*log10(2.d0**4)-4.d0))then
geoel(10, ielem) = 0.d0
else
geoel(10, ielem) = 100.d0
endif
!
!if(ielem.eq.1489) print*,'ielem',log10(unvar2(1)/unvar1(1))
!
650 enddo

end subroutine getSIHOgd_Persson
!
!...High order barth limiter on curvilinear meshes with smooth indicator transition....
!
subroutine barthlimit_lag_curvquadcb_HOsmth(geoel, coord, coold, ustar, unkno, ipqua, &
bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
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
real*8:: smthid, fi11,fi12,fi21,fi22
!
!...Part I: Some preparing work
!

!... Annihilate the HO moments of DGp2
do ie = 1,nquad

ielem = ie + ntria
!...Shock cell
if(geoel(10, ielem).gt.15.d0)then
unkno(4:6, :, ielem) = 0.d0
endif
enddo


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

!...Store the maximum and minimum values surrounding one cell...
do ie = 1, nquad
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)
do iq=1, nq+2
unmax_new(iq, ielem) = maxval(unmax(iq, ipq(1:nvqua)))
unmin_new(iq, ielem) = minval(unmin(iq, ipq(1:nvqua)))
enddo
enddo
!
!...Part II: Impose limiter
!
do ie = 1, nquad

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Skip the smooth cells
if(geoel(10, ielem).lt.5.d0) cycle
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
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

!...zero out unknv
unknvq = 0.d0
do iv   = 1,nvqua
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo

if(ndens.eq.1)then
rhov = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
r = xvq(iv); s = yvq(iv)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))

call getrhoig_quadcurv(rhoi, xpqi)!
call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)

rhov = rhon
!
print*,'Please chech the P2 density evolution for ndens=2 '
stop

elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
bqv(1, iv) = 1.d0
bqv(2, iv) = (xvq(iv)-rcv)/dr
bqv(3, iv) = (yvq(iv)-scv)/ds

unknvq(1, iv) =0.d0
do ideg = 1,mdegr
unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
enddo
rhov = unknvq(1, iv)

print*,'Please chech the P2 density evolution for ndens=2 '
stop

endif

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvq(1, iv) = 1.d0/rhov
elseif(ndens.eq.2)then
unknvq(1, iv) = rhov
elseif(ndens.eq.3)then
unknvq(1, iv) = rhov
endif

unknvq(4 ,iv) = pvtx
enddo !do iv   = 1,nvqua
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
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif
unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq)  = pctr
!
do iv = 1, nvqua
do iq = 1, nq
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)

alfa(iq, iv) = afbar
enddo

!...Special treatment of boundary...
if(indpt(ipq(iv)).eq.1)then

!...Set special boundary
!  alfa(:, iv) = 1.d0
!
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
enddo ! do iv = 1, nvqua

!...Get the minimum value for one cell
do iq = 1,nq
aflim(iq, ielem) = minval(alfa(iq, 1:4))
enddo
!
enddo
!
!...Part 2.1: Impose symmetry preserving limiter for velocity...
!
call barthlimit_symprecb_HOsmth(geoel, coord, ustar, unkno, ipqua, bface, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
!
!...Part 3: Correct total energy
!
do ie = 1,nquad

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Skip the smooth cells
if(geoel(10, ielem).lt.5.d0) cycle

!...Shape function
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
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

!...Cell average of inverse density, velocity and total energy.
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
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif

unctr(2:3) = unkno(1, 2:3, ielem)
unctr(nq) = pctr
unctr(nq+1) = ectr

!...zero out unknv
unknvq = 0.d0
!
do iv   = 1,nvqua
do ideg = 1, mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
!
enddo
!
if(ndens.eq.1)then
rhov = 1.d0/unknvq(1, iv)
elseif(ndens.eq.2)then
r = xvq(iv); s = yvq(iv)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quadcurv(rhoi, xpqi)!
call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)
!
rhov = rhon
elseif(ndens.eq.3)then
rcv = geoel(5, ielem); scv = geoel(6, ielem)
!
bqv(1, iv) = 1.d0
bqv(2, iv) = (xvq(iv)-rcv)/dr
bqv(3, iv) = (yvq(iv)-scv)/ds
!
unknvq(1, iv) =0.d0
do ideg = 1, mdegr
unknvq(1, iv) = unknvq(1, iv) + unkno(ideg,1,ielem)*bqv(ideg, iv)
enddo

rhov = unknvq(1, iv)
endif !if(ndens.eq.1)then

uvtx = unknvq(2, iv)
vvtx = unknvq(3, iv)
evtx = unknvq(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvq(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
elseif(ndens.eq.2)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
elseif(ndens.eq.3)then
unknvq(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
!unknvq(2, iv) = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
!unknvq(3, iv) = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
unknvq(4 ,iv) = pctr + aflim(4, ielem)*(pvtx - pctr)
unknvq(5, iv) = unknvq(4 ,iv)/(gamlg-1.d0)*unknvq(1, iv) + 0.5d0*(unknvq(2, iv)**2 + unknvq(3, iv)**2)
!
enddo

!...Get the limiting coefficient for the total energy
do iv = 1, nvqua
do iq = nq+1, nq+1 !...energy
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(unmax(iq, ipq(iv)), unmin(iq, ipq(iv)), unctr(iq), dunk(iq), afbar)
! call barthfct(unmax_new(iq, ielem), unmin_new(iq, ielem), unctr(iq), dunk(iq), afbar)

alfa(iq, iv) = afbar
enddo

!...Special treatment for the boundary nodes...
! if(indpt(ipq(iv)).eq.1)then
!alfa(5, iv) = 1.d0
! endif
enddo !do iv = 1, nvqua
!
!if(ielem.eq.89)print*,'limiter',ielem,unmax(nq+1, ipq(3)), unmin(nq+1, ipq(3))
!...Get the minimum value for one cell
do iq = nq+1,nq+1
aflim(iq, ielem) = minval(alfa(iq, 1:4))
enddo
enddo !...Do ie = 1, nquad
!
!...Part 4: Annihilate the HO moments of DGp2
!
do ie = 1,nquad

ielem = ie + ntria

!...Shock cell
if(geoel(10, ielem).gt.15.d0)then
unkno(4:6, :, ielem) = 0.d0

!...Smooth cell
elseif(geoel(10, ielem).lt.5.d0)then
aflim(:, ielem) = 1.d0

afvec(1, 1, ielem) = 1.d0
afvec(1, 2, ielem) = 0.d0
afvec(2, 1, ielem) = 0.d0
afvec(2, 2, ielem) = 1.d0

!...Transition cells
else

smthid = geoel(10,ielem)-10.d0
fi11 = smthid  + (1.d0 - smthid)*afvec(1, 1, ielem)
fi12 =           (1.d0 - smthid)*afvec(1, 2, ielem)

fi21 =           (1.d0 - smthid)*afvec(2, 1, ielem)
fi22 = smthid  + (1.d0 - smthid)*afvec(2, 2, ielem)

!
afvec(1, 1, ielem) = fi11
afvec(1, 2, ielem) = fi12
afvec(2, 1, ielem) = fi21
afvec(2, 2, ielem) = fi22

!
aflim(1, ielem) = smthid + (1.d0 - smthid)*aflim(1, ielem)
aflim(4, ielem) = smthid + (1.d0 - smthid)*aflim(4, ielem)
aflim(5, ielem) = smthid + (1.d0 - smthid)*aflim(5, ielem)
!

endif

enddo

end subroutine barthlimit_lag_curvquadcb_HOsmth
!
!...High order barth limiter (vector) on curvilinear meshes with smooth indicator transition....
!
subroutine barthlimit_symprecb_HOsmth(geoel, coord, ustar, unkno, ipqua, bface, intfac, afvec,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer, dimension(1:nvqua,1:nquad),           intent(in):: ipqua
integer*4,dimension(1:nbfai,nbfac),           intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
!
!...Local
!
integer:: ipq(nvqua)
integer:: indpt(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp, ivf
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
indpt = 0  !...indpt represents index of boundary node
!
!...Coloring the nodes at x axis
!
if(ncase.ne.1)then

do ifa =1 ,nbfac
!
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.221)then
indpt(intfac(3:(2+nvfac), ifa)) = 3
endif
endif
enddo

!...Coloring the nodes at y axis
do ifa =1 ,nbfac
!
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.222)then
do ivf = 1, nvfac
if(indpt(intfac(2+ivf, ifa)).eq.3)then
indpt(intfac(2+ivf, ifa)) = 10
else
indpt(intfac(2+ivf, ifa)) = 2
endif
enddo
endif
endif
enddo

else

do ifa =1 ,nbfac
!
indpt(intfac(3:(2+nvfac), ifa)) = 1
enddo
endif
!
eps = 1.e-6
mapmt = 0.d0
!
!...For basis function scaling...
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
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Skip the smooth cells
if(geoel(10, ielem).lt.5.d0) cycle
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
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
!...Get dudx, dudy, dvdx, dvdy
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy
!
!...LANL
!
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy
!
matra = unkno(1, 2, ielem)
matrb = unkno(1, 3, ielem)
if(sqrt(matra**2+matrb**2).gt.1.d-8)then
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
!
enddo
!
!...Part 4: Impose limiter for mapped u and v
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Skip the smooth cells
if(geoel(10, ielem).lt.5.d0) cycle
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
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
!...Symmetry BC...
!
if(ncase.ne.1)then
if(indpt(ipq(iv)).eq.3)then
!
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

elseif(indpt(ipq(iv)).eq.2)then

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

elseif(indpt(ipq(iv)).eq.10)then
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

!...Vector limiting
do iq = 2, 3
dunk(iq) = unknvq(iq, iv) - unctr(iq)
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo

!...Special treatment about the boundary
!if(indpt(ipq(iv)).ne.0)then
!alfa(2:3, iv) = 1.d0
!endif
if(indpt(ipq(iv)).eq.1)then
if(ncase.eq.1)then
alfa(2:3, iv) = 1.d0
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

!...Skip the smooth cells
if(geoel(10, ielem).lt.5.d0) cycle
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
!
afvec(1, 1, ielem) = fi11
afvec(1, 2, ielem) = fi12
afvec(2, 1, ielem) = afvec(1, 2, ielem)
afvec(2, 2, ielem) = fi22
!!
enddo
!
end subroutine barthlimit_symprecb_HOsmth
!
!...Smooth indicator transition....
!
subroutine getSIHOsmth_Persson(intfac, ipqua, coord, coold, geoel, unkno)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(inout)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua

!...Local integer
integer::ie,ig,ideg,ishp,iv,ielem

!...local integer array
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds, bv
real*8, dimension(1:nq):: unknod,unkp1,unvar1,unvar2
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2,ngausdq)

!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::djaco, wi
real*8::ksup,supi,sup0,khalf,lhalf,ksup2
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
!...Loop over quad
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

!...Points consitituting one element...
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

!....Zero out unvar2, unvar1
!...Numerator and denumerator of the smooth indicator

unvar2 = 0.d0
unvar1 = 0.d0

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
djaco = wi*(dxdr*dyds - dydr*dxds)

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

!...Solution at P1
unkp1 = 0.d0

if(npoly.eq.1)then
do ideg =1, 1!
unkp1(1:nq) = unkp1(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo
elseif(npoly.eq.2)then
do ideg =1, 3!
unkp1(1:nq) = unkp1(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo
endif

unvar2(1:nq) = unvar2(1:nq) + (unkp1(1:nq)-unknod(1:nq))**2*djaco
unvar1(1:nq) = unvar1(1:nq) + (unknod(1:nq))**2*djaco
!
enddo !...(2)ig = 1,ngausd

!...Coarse smooth indicator
geoel(10, ielem) = 1.d0
!
supi = log10(unvar2(1)/unvar1(1))
sup0 = -3.d0*log10(4.d0)
ksup = 4.d0
ksup2 = 3.d0

!...Smooth region
if(supi.le.(sup0-ksup))then
!if(log10(unvar2(1)/unvar1(1)).le.(-1.5d0*log10(2.d0**4)-4.d0))then
geoel(10, ielem) = 1.d0

!...Shock region
elseif(supi.gt.(sup0-ksup2))then
geoel(10, ielem) = 100.d0
!
!print*,'limit',ielem,geoel(10, ielem)

!...Transition region
else
khalf = .5d0*(sup0-ksup+sup0-ksup2)
lhalf = .5d0*(ksup-ksup2)
geoel(10, ielem) = 10.d0 + 0.5d0*(1.d0-sin(0.5d0*pi*(supi-khalf)/lhalf))
endif
!
650 enddo

end subroutine getSIHOsmth_Persson
!
!...subroutine: Riemann input for hybrid curved quads using sub-cell scheme for 5-points integration with smooth indicator transition ....
!
subroutine getriem_quadsubg_gauss_smth(ipqua, geoel, gesgq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, drhosgq, coord, coold, aflim, afvec)
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
real*8,dimension(1:4, 1:nsize), intent(out)::drhosgq

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
if(nfint.eq.5.or.nfint.eq.7.or.nfint.eq.10)then
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

!if(ielem.eq.12)then
!print*,'smooth2',ielem,unkno(:,:,ielem)
!endif

!...Get density for sub-cells...
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
call getrhosubg_averg(xpq,xpqi,unksgq, geoq_sub)

!...Get density correction
!call getdens_quadsubg(rc, sc, geoel(19:21, ielem), unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc)
call getrhosubg_daverg(rc,sc, geoel(19:21, ielem), xpq, unksgq, unkno(:,:,ielem), aflim(:,ielem), uqsgc)

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
rhovsg = rhovt +.8d0*(unksgq(1,isg)-1.d0/uqsgc(1, isg))
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
!if(ielem.eq.115)then
!print*,'Variabe1',ielem,isg,ivsg,rhovsg,pvtx,unsgq(2:3, ivsg),evtx,geoel(10,ielem)
!endif

!...Limiter
!...For high-order method(>P1), if the cell includes shock, then
!...the high-order terms are chopped, and only limit the linear DG(P1) part.

if(nlimi.eq.6)then

!...Shock cell
if(geoel(10, ielem).gt.15.d0)then

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
!uvtx = unkno(1,2,ielem)  + afvec(1, 1, ielem)*(uvtx - uctr) + afvec(1, 2, ielem)*(vvtx - vctr)
!vvtx = unkno(1,3,ielem)  + afvec(2, 1, ielem)*(uvtx - uctr) + afvec(2, 2, ielem)*(vvtx - vctr)
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)

!...Updtae velocity
!unsgq(2, ivsg) = unkno(1,2,ielem)  + afvec(1, 1, ielem)*(uvtx - uctr) + afvec(1, 2, ielem)*(vvtx - vctr)
!unsgq(3 ,ivsg) = unkno(1,3,ielem)  + afvec(2, 1, ielem)*(uvtx - uctr) + afvec(2, 2, ielem)*(vvtx - vctr)

!...Non-shock cell
!else

!uvtx = unkno(1,2,ielem)  + afvec(1, 1, ielem)*(uvtx - unkno(1,2,ielem)) + afvec(1, 2, ielem)*(vvtx - unkno(1,3,ielem))
!vvtx = unkno(1,3,ielem)  + afvec(2, 1, ielem)*(uvtx - unkno(1,2,ielem)) + afvec(2, 2, ielem)*(vvtx - unkno(1,3,ielem))
!
!pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)

!...Updtae velocity
unsgq(2, ivsg) = uvtx
unsgq(3 ,ivsg) = vvtx

!endif

!...Transition cells
elseif(geoel(10,ielem).lt.15.d0.and.geoel(10,ielem).gt.5.d0)then

if(ndens.eq.1)then
rhovsg = 1.d0/(rhomc + aflim(1, ielem)*(1.d0/rhovsg - rhomc))
elseif(ndens.eq.3)then
rhovsg = unksgq(1,isg) + aflim(1, ielem)*(rhovsg - unksgq(1,isg))
endif
!
pvtx = pctr + aflim(4, ielem)*(pvtx - pctr)

!...Updtae velocity
unsgq(2, ivsg) = unkno(1,2,ielem)  + afvec(1, 1, ielem)*(uvtx - uctr) + afvec(1, 2, ielem)*(vvtx - vctr)
unsgq(3 ,ivsg) = unkno(1,3,ielem)  + afvec(2, 1, ielem)*(uvtx - uctr) + afvec(2, 2, ielem)*(vvtx - vctr)

endif

endif

!...Get stress tensor at one vertex
sigma(1, 1, ivsg) = -pvtx
sigma(1, 2, ivsg) = 0.d0
sigma(2, 1, ivsg) = 0.d0
sigma(2, 2, ivsg) = -pvtx
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
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
aujmp(3,:)=aujmp(3,:)/sdctr

!...Get impedence coefficient...
do ivsg   = 1, 4
dux= vlave(1, ipq(ipqsg(ivsg, isg)))-unsgq(2, ivsg)
duy= vlave(2, ipq(ipqsg(ivsg, isg)))-unsgq(3, ivsg)
deltu = sqrt(dux**2 + duy**2)
do ifa = 1, 2
deltu = 4.d0*abs(dux*vnorm(1, ifa, ivsg) + duy*vnorm(2, ifa, ivsg))
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
!if(ipq(iv).eq.149) print*,'p36 muacn(vv) post',ie,isg,ivsg,murie(ifa, ivsg),sigma(1, 1, ivsg),unsgq(2:3, ivsg)

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

end subroutine getriem_quadsubg_gauss_smth
!
!....domain integral for hybrid curv quad cells with smooth indicator transition
!
subroutine rhsdomndg_lagmc_curvquad_smth(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec, vnulq )
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
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: dudr, duds, dvdr, dvds
real*8::dr,ds,rc,sc, rcv, scv
real*8::rm,sm,rp,sp
real*8::xg, yg
real*8::rhoad,uadv,vadv,eadv,rhoma
real*8::duadv,dvadv
real*8::pres
real*8::djaco, wi
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
!...Loop over quad
!
do 650 ie = 1,nquad!...(1)ie = 1,nelem

ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

!...Points consitituting one element...
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

!...Limiting
!if(nlimi.eq.6.and.geoel(10,ielem).gt.10.d0)then
if(nlimi.eq.6)then

!...Cell average
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

!...Shcok cells
if(geoel(10,ielem).gt.15.d0)then
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)

uadv = unkno(1,2,ielem)  + dudr*b(2) + duds*b(3)
vadv = unkno(1,3,ielem)  + dvdr*b(2) + dvds*b(3)
!duadv = uadv - unkno(1,2,ielem)
!dvadv = vadv - unkno(1,3,ielem)
!
!uadv = unkno(1,2,ielem)  + afvec(1, 1, ielem)*duadv + afvec(1, 2, ielem)*dvadv
!vadv = unkno(1,3,ielem)  + afvec(2, 1, ielem)*duadv + afvec(2, 2, ielem)*dvadv
!
!if(geoel(10,ielem).lt.10.d0)then
!print*,'ielem',ielem,geoel(10,ielem)
!if(abs(uadv-unknod(2)).ge.1d-14.or.abs(vadv-unknod(3)).ge.1d-14)then
!print*,'non equal',ielem,uadv,unknod(2),vadv,unknod(3), afvec(2, 1, ielem)*duadv + afvec(2, 2, ielem)*dvadv,dvadv
!endif
!endif
!
pres = pctr + aflim(4, ielem)*(pres- pctr)

!...Transition cells
elseif(geoel(10,ielem).lt.15.d0.and.geoel(10,ielem).gt.5.d0)then
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)

!
!print*,'ielem',ielem,geoel(10,ielem),unkno(4:6,:,ielem),uadv,vadv,dudr*b(2) + duds*b(3),&
!afvec(1, 1, ielem)*(uadv - unkno(1,2,ielem)) + afvec(1, 2, ielem)*(vadv - unkno(1,3,ielem)),&
!dvdr*b(2) + dvds*b(3),&
!afvec(2, 1, ielem)*(uadv - unkno(1,2,ielem)) + afvec(2, 2, ielem)*(vadv - unkno(1,3,ielem))

!uadv = unkno(1,2,ielem)  + dudr*b(2) + duds*b(3)
!vadv = unkno(1,3,ielem)  + dvdr*b(2) + dvds*b(3)
duadv = uadv - unkno(1,2,ielem)
dvadv = vadv - unkno(1,3,ielem)
!
uadv = unkno(1,2,ielem)  + afvec(1, 1, ielem)*duadv + afvec(1, 2, ielem)*dvadv
vadv = unkno(1,3,ielem)  + afvec(2, 1, ielem)*duadv + afvec(2, 2, ielem)*dvadv
!
!if(geoel(10,ielem).lt.10.d0)then
!print*,'ielem',ielem,geoel(10,ielem)
!if(abs(uadv-unknod(2)).ge.1d-14.or.abs(vadv-unknod(3)).ge.1d-14)then
!print*,'non equal',ielem,uadv,unknod(2),vadv,unknod(3), afvec(2, 1, ielem)*duadv + afvec(2, 2, ielem)*dvadv,dvadv
!endif
!endif
!
!
!if(ielem.eq.173) print*,'duadv',duadv,dvadv,&
!unkno(1,2,ielem)  + dudr*b(2) + duds*b(3),&
!unkno(1,3,ielem)  + dvdr*b(2) + dvds*b(3),&
!uadv,vadv,6.

pres = pctr + aflim(4, ielem)*(pres- pctr)
endif
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
!if(ie==1772) print*,'rhs iface idegr',ie,ig,fluxd(1:3,4),pres,djaco
!
enddo !...(2)ig = 1,ngausd
!
650 enddo

end subroutine rhsdomndg_lagmc_curvquad_smth
!
!...Calculate the velocity at the Gauss point for novertex face...
!
subroutine getSIHO_facejmp(geoel,intfac,ipqua,coord,unkno)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic,ishp
integer::iv,fnoel
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(nvfac,4) :: fglvq
!
real*8::eps
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, ngausf)
real*8::posif(1:ngausf),weighf(ngausf)
real*8::xpf(1:2, 1:nvfac)
real*8::dshpr(3)

!...Local real number
real*8::detma,dudr,duds,dvdr,dvds
real*8::pvtx,rhovt, rhomc, rhomv, rhovl, rhovr,rhsu1,rhsu2
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx
real*8::fnx,fny, ftx, fty, rho
real*8::rf,s
real*8::dr, ds,rc,sc
real*8::othog
real*8::acnx,acny,shp1,shp2,shp3,r,delu
real*8::dpres,averp
real*8::dxdr,dydr
real*8::djaco,wi
!
real*8,allocatable::sfael(:)
!
allocate(sfael(ncell))
sfael = 0.d0
geoel(10,:) = 0.d0
!
call ruqope_lobatto(1, ngausf, posif, weighf)

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
!...Part II: Get the corrected velocity at the non-vertex gauss point...
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

!...Boundary face
if(ifa.le.nbfac)then

!...Identify the fgaus for left cell
call getfgfj_glb(ipqua(1:nvqua, iel), ipf, fnoel)
!
do ig = 1, ngausf

!...Get the normal vector for the GPs
rf = posif(ig)
wi = weighf(ig)
!
dshpr(1) = -0.5d0 + rf
dshpr(2) =  0.5d0 + rf
dshpr(3) = -2.d0*rf

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


!...Get the variables at the GPs

!...Get the reference coordinates for the GPs
r = (rf + 1.d0)/2.d0*(xvq(fglvq(2, fnoel))-xvq(fglvq(1, fnoel))) + xvq(fglvq(1, fnoel))
s = (rf + 1.d0)/2.d0*(yvq(fglvq(2, fnoel))-yvq(fglvq(1, fnoel))) + yvq(fglvq(1, fnoel))

!...Zero out unkng for the Gauss points
unkng = 0.d0

!...Basis fctn...
bq(1, ig) = 1.d0
bq(2, ig) = (r-rc)/dr
bq(3, ig) = (s-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, iel)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, iel)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, iel)
endif

!...Unkno at gauss points...

do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*bq(ideg, ig)
enddo

!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unkng(nq+1, ig, 1) = pvtx

!
sfael(iel) = sfael(iel) + djaco*wi
geoel(10, iel) = geoel(10 ,iel) + 0.d0
!
!if(iel.eq.900) print*,'jump',geoel(10,iel)
!
enddo

!...Interior face
elseif(ifa.gt.nbfac)then

!...Zero out unkng for the Gauss points
unkng = 0.d0

!...Identify the fgaus for left cell
call getfgfj_glb(ipqua(1:nvqua, iel), ipf, fnoel)
!
do ig = 1, ngausf

!...Get the normal vector for the GPs
rf = posif(ig)
wi = weighf(ig)
!
dshpr(1) = -0.5d0 + rf
dshpr(2) =  0.5d0 + rf
dshpr(3) = -2.d0*rf

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

!...Get the variables at the GPs

!...Get the reference coordinates for the GPs
r = (rf + 1.d0)/2.d0*(xvq(fglvq(2, fnoel))-xvq(fglvq(1, fnoel))) + xvq(fglvq(1, fnoel))
s = (rf + 1.d0)/2.d0*(yvq(fglvq(2, fnoel))-yvq(fglvq(1, fnoel))) + yvq(fglvq(1, fnoel))

!...Basis fctn...
bq(1, ig) = 1.d0
bq(2, ig) = (r-rc)/dr
bq(3, ig) = (s-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, iel)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, iel)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, iel)
endif

!...Unkno at gauss points...

do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*bq(ideg, ig)
enddo

!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!if(ier.eq.900)print*,'jump',ier,iel,unkng(nq+1, ig, 2),unkng(nq+1, ig, 1),ig,uvtx,vvtx,evtx,rhovt,pvtx
!
unkng(nq+1, ig, 1) = pvtx
!
!if(ier.eq.900)print*,'jump3',ier,iel,unkng(nq+1, 1, 1),ig
!
sfael(iel) = sfael(iel) + djaco*wi
!
enddo
!

!if(ier.eq.900)print*,'jump2',ier,iel,unkng(nq+1, :, 1)

!...Right cell
!...Identify the fgaus for left cell
call getfgfj_glb(ipqua(1:nvqua, ier), ipf, fnoel)
!
do ig = 1, ngausf

!...Get the normal vector for the GPs
rf = posif(ig)
wi = weighf(ig)
!
dshpr(1) = -0.5d0 + rf
dshpr(2) =  0.5d0 + rf
dshpr(3) = -2.d0*rf

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

!...Get the variables at the GPs

!...Get the reference coordinates for the GPs
r = (rf + 1.d0)/2.d0*(xvq(fglvq(1, fnoel))-xvq(fglvq(2, fnoel))) + xvq(fglvq(2, fnoel))
s = (rf + 1.d0)/2.d0*(yvq(fglvq(1, fnoel))-yvq(fglvq(2, fnoel))) + yvq(fglvq(2, fnoel))


!...Basis fctn...
bq(1, ig) = 1.d0
bq(2, ig) = (r-rc)/dr
bq(3, ig) = (s-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, ier)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, ier)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, ier)
endif

!...Unkno at gauss points...

do ideg =1, mdegr
unkng(1:nq, ig, 2) = unkng(1:nq, ig ,2) + unkno(ideg,1:nq,ier)*bq(ideg, ig)
enddo

!
rhovt = 1.d0/unkng(1, ig, 2)
uvtx = unkng(2, ig, 2)
vvtx = unkng(3, ig, 2)
evtx = unkng(4, ig, 2)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unkng(nq+1, ig, 2) = pvtx
!
sfael(ier) = sfael(ier) + djaco*wi
!
dpres = abs(unkng(nq+1, ig, 2)-unkng(nq+1, ig, 1))
averp = 0.5d0*abs(unkng(nq+1, ig, 2) + unkng(nq+1, ig, 1))
!
geoel(10,iel) = geoel(10, iel) + dpres/averp*djaco*wi
geoel(10,ier) = geoel(10, ier) + dpres/averp*djaco*wi
!
!if(ier.eq.900)print*,'jump',ier,iel,unkng(nq+1, ig, 2),unkng(nq+1, ig, 1),geoel(10, ier),ig

enddo


endif
450 enddo
!
!...Part II: Designate shocked cell
!
do ie = 1, ncell
!
!if(ie.eq.900) print*,'pres',ie,geoel(10,ie)

if(log10(geoel(10,ie)/sfael(ie)).ge.-1.d0)then
geoel(10, ie) = 100
else
geoel(10 ,ie) = 0
endif

enddo
!
deallocate(sfael)
!
end subroutine getSIHO_facejmp
!
!...Calculate the velocity at the Gauss point for novertex face...
!
subroutine getSIHO_facejmp2(geoel,intfac,ipqua,coord,unkno)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ic,ishp
integer::iv,fnoel,fnoer
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(nvfac,4) :: fglvq
!
real*8::eps
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, ngausf)
real*8::posif(1:ngausf),weighf(ngausf)
real*8::xpf(1:2, 1:nvfac)
real*8::dshpr(3)

!...Local real number
real*8::detma,dudr,duds,dvdr,dvds
real*8::pvtx,rhovt, rhomc, rhomv, rhovl, rhovr,rhsu1,rhsu2
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx
real*8::fnx,fny, ftx, fty, rho
real*8::rf,s
real*8::dr, ds,rc,sc,rl,sl,rr,sr
real*8::rcl,scl,rcr,scr
real*8::othog
real*8::acnx,acny,shp1,shp2,shp3,r,delu
real*8::dpres,averp
real*8::dxdr,dydr
real*8::djaco,wi
!
real*8,allocatable::sfael(:)
!
allocate(sfael(ncell))
sfael = 0.d0
geoel(10,:) = 0.d0
!
call ruqope_lobatto(1, ngausf, posif, weighf)

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
!...Part II: Get the corrected velocity at the non-vertex gauss point...
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

!...Boundary face
if(ifa.le.nbfac)then

!...Identify the fgaus for left cell
call getfgfj_glb(ipqua(1:nvqua, iel), ipf, fnoel)
!
do ig = 1, ngausf

!...Get the normal vector for the GPs
rf = posif(ig)
wi = weighf(ig)
!
dshpr(1) = -0.5d0 + rf
dshpr(2) =  0.5d0 + rf
dshpr(3) = -2.d0*rf

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


!...Get the variables at the GPs

!...Get the reference coordinates for the GPs
r = (rf + 1.d0)/2.d0*(xvq(fglvq(2, fnoel))-xvq(fglvq(1, fnoel))) + xvq(fglvq(1, fnoel))
s = (rf + 1.d0)/2.d0*(yvq(fglvq(2, fnoel))-yvq(fglvq(1, fnoel))) + yvq(fglvq(1, fnoel))

!...Zero out unkng for the Gauss points
unkng = 0.d0

!...Basis fctn...
bq(1, ig) = 1.d0
bq(2, ig) = (r-rc)/dr
bq(3, ig) = (s-sc)/ds

!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, iel)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, iel)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, iel)
endif

!...Unkno at gauss points...

do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*bq(ideg, ig)
enddo

!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unkng(nq+1, ig, 1) = pvtx

!
sfael(iel) = sfael(iel) + djaco*wi
geoel(10, iel) = geoel(10 ,iel) + 0.d0
!
!if(iel.eq.900) print*,'jump',geoel(10,iel)
!
enddo

!...Interior face
elseif(ifa.gt.nbfac)then

!...Zero out unkng for the Gauss points
unkng = 0.d0

!...Identify the fgaus for left and right cell
call getfgfj_glb(ipqua(1:nvqua, iel), ipf, fnoel)
call getfgfj_glb(ipqua(1:nvqua, ier), ipf, fnoer)

!...Left and right mass center
rcl = geoel(1, iel); scl = geoel(2, iel)
rcr = geoel(1, ier); scr = geoel(2, ier)

!
do ig = 1, ngausf

!...Get the normal vector for the GPs
rf = posif(ig)
wi = weighf(ig)
!
dshpr(1) = -0.5d0 + rf
dshpr(2) =  0.5d0 + rf
dshpr(3) = -2.d0*rf

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

!...Get the variables at the GPs

!...Get the reference coordinates for the GPs
rl = (rf + 1.d0)/2.d0*(xvq(fglvq(2, fnoel))-xvq(fglvq(1, fnoel))) + xvq(fglvq(1, fnoel))
sl = (rf + 1.d0)/2.d0*(yvq(fglvq(2, fnoel))-yvq(fglvq(1, fnoel))) + yvq(fglvq(1, fnoel))

!...Basis fctn...
bq(1, ig) = 1.d0
bq(2, ig) = (rl-rcl)/dr
bq(3, ig) = (sl-scl)/ds

!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, iel)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, iel)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, iel)
endif

!...Unkno at gauss points...

do ideg =1, mdegr
unkng(1:nq, ig, 1) = unkng(1:nq, ig ,1) + unkno(ideg,1:nq,iel)*bq(ideg, ig)
enddo

!
rhovt = 1.d0/unkng(1, ig, 1)
uvtx = unkng(2, ig, 1)
vvtx = unkng(3, ig, 1)
evtx = unkng(4, ig, 1)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!if(ier.eq.900)print*,'jump',ier,iel,unkng(nq+1, ig, 2),unkng(nq+1, ig, 1),ig,uvtx,vvtx,evtx,rhovt,pvtx
!
unkng(nq+1, ig, 1) = pvtx

!...Get the reference coordinates for the GPs
rr = (rf + 1.d0)/2.d0*(xvq(fglvq(1, fnoer))-xvq(fglvq(2, fnoer))) + xvq(fglvq(2, fnoer))
sr = (rf + 1.d0)/2.d0*(yvq(fglvq(1, fnoer))-yvq(fglvq(2, fnoer))) + yvq(fglvq(2, fnoer))

!...Basis fctn...
bq(1, ig) = 1.d0
bq(2, ig) = (rr-rcr)/dr
bq(3, ig) = (sr-scr)/ds

!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, ier)
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, ier)
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, ier)
endif

!...Unkno at gauss points...

do ideg =1, mdegr
unkng(1:nq, ig, 2) = unkng(1:nq, ig ,2) + unkno(ideg,1:nq,ier)*bq(ideg, ig)
enddo

!
rhovt = 1.d0/unkng(1, ig, 2)
uvtx = unkng(2, ig, 2)
vvtx = unkng(3, ig, 2)
evtx = unkng(4, ig, 2)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!if(iel.eq.145)print*,iel,ier,ig,rcr,scr,rr,sr,geoel(19:21,ier)
!
unkng(nq+1, ig, 2) = pvtx
!
!if(ier.eq.900)print*,'jump3',ier,iel,unkng(nq+1, 1, 1),ig
!
sfael(iel) = sfael(iel) + djaco*wi
sfael(ier) = sfael(ier) + djaco*wi
!
!
dpres = abs(unkng(nq+1, ig, 2)-unkng(nq+1, ig, 1))
averp = 0.5d0*abs(unkng(nq+1, ig, 2) + unkng(nq+1, ig, 1))
!
geoel(10,iel) = geoel(10, iel) + dpres/averp*djaco*wi
geoel(10,ier) = geoel(10, ier) + dpres/averp*djaco*wi
!
!if(iel.eq.145)print*,iel,ier,dpres,averp,djaco*wi,unkng(nq+1, ig, 2),unkng(nq+1, ig, 1)
!if(ier.eq.145)print*,iel,ier,dpres,averp,djaco*wi
!
enddo

endif
450 enddo
!
!...Part II: Designate shocked cell
!
do ie = 1, ncell
!
if(log10(geoel(10,ie)/sfael(ie)).ge.0.d0)then
geoel(10, ie) = 100
else
geoel(10 ,ie) = 0
endif

enddo
!
deallocate(sfael)
!
end subroutine getSIHO_facejmp2
!
!...Calculate the jump over the face for smooth indicator...
!
subroutine getSIHO_facejmp3(geoel,intfac,ipqua,coord,unkno)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(inout)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
!...Local integer
integer::ifa,iel,ier,ie,idfal,ig,ideg,jdeg,ic,ishp,isf
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(nvfac,4) :: fglvq
!
real*8::eps
real*8::unkng(1:nq+3, 1:ngausf, 1:2)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, ngausf)
real*8::posif(1:ngausf),weighf(ngausf)
real*8::xpf(1:2, 1:nvfac)
real*8::dshpr(3)
real*8::rca(2),sca(2)
integer::fnoe(2)

!...Local real number
real*8::pvtx,rhovt,uvtx,vvtx,evtx
real*8::rf,dr,ds,r,s
real*8::dpres,averp
real*8::dxdr,dydr,djaco,wi
!
real*8,allocatable::sfael(:)
!
allocate(sfael(ncell))
!
eps = 1.d-6
!
sfael = 0.d0
geoel(10,:) = 0.d0
!
call ruqope_lobatto(1, ngausf, posif, weighf)
!
dr = 1.d0
ds = 1.d0
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
!...Part II: Loop over the faces to get the jump...
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

!...Left and right mass center
rca(1) = geoel(1, iel); sca(1) = geoel(2, iel)
rca(2) = geoel(1, ier); sca(2) = geoel(2, ier)

!...Boundary face
if(ifa.le.nbfac)then
!
do ig = 1, ngausf

!...Gauss points and weights
rf = posif(ig)
wi = weighf(ig)
!
dshpr(1) = -0.5d0 + rf
dshpr(2) =  0.5d0 + rf
dshpr(3) = -2.d0*rf

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
sfael(iel) = sfael(iel) + djaco*wi
geoel(10, iel) = geoel(10 ,iel) + 0.d0
!
enddo

!...Interior face
elseif(ifa.gt.nbfac)then

!...Zero out unkng for the Gauss points
unkng = 0.d0

!...Identify the fgaus for left and right cell
call getfgfj_glb(ipqua(1:nvqua, iel), ipf, fnoe(1))
call getfgfj_glb(ipqua(1:nvqua, ier), ipf, fnoe(2))

!
do ig = 1, ngausf

!...Gauss points and weights
rf = posif(ig)
wi = weighf(ig)
!
dshpr(1) = -0.5d0 + rf
dshpr(2) =  0.5d0 + rf
dshpr(3) = -2.d0*rf

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

!...Get the variables at the GPs for the left(right) cell
do isf = 1, 2
!...Get the reference coordinates for the GPs
if(isf.eq.1)then
r = (rf + 1.d0)/2.d0*(xvq(fglvq(2, fnoe(isf)))-xvq(fglvq(1, fnoe(isf)))) + xvq(fglvq(1, fnoe(isf)))
s = (rf + 1.d0)/2.d0*(yvq(fglvq(2, fnoe(isf)))-yvq(fglvq(1, fnoe(isf)))) + yvq(fglvq(1, fnoe(isf)))
else
r = (rf + 1.d0)/2.d0*(xvq(fglvq(1, fnoe(isf)))-xvq(fglvq(2, fnoe(isf)))) + xvq(fglvq(2, fnoe(isf)))
s = (rf + 1.d0)/2.d0*(yvq(fglvq(1, fnoe(isf)))-yvq(fglvq(2, fnoe(isf)))) + yvq(fglvq(2, fnoe(isf)))
endif

!...Basis fctn...
bq(1, ig) = 1.d0
bq(2, ig) = (r-rca(isf))/dr
bq(3, ig) = (s-sca(isf))/ds

!DGP2
if(npoly.eq.2)then
bq(4, ig) = 0.5d0*bq(2, ig)*bq(2, ig) - geoel(19, intfac(isf, ifa))
bq(5, ig) = 0.5d0*bq(3, ig)*bq(3, ig) - geoel(20, intfac(isf, ifa))
bq(6, ig) =       bq(2, ig)*bq(3, ig) - geoel(21, intfac(isf, ifa))
endif

!...Unkno at gauss points...
do ideg =1, mdegr
unkng(1:nq, ig, isf) = unkng(1:nq, ig ,isf) + unkno(ideg,1:nq,intfac(isf, ifa))*bq(ideg, ig)
enddo
!
rhovt = 1.d0/unkng(1, ig, isf)
uvtx = unkng(2, ig, isf)
vvtx = unkng(3, ig, isf)
evtx = unkng(4, ig, isf)
!
pvtx = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
unkng(nq+1, ig, isf) = pvtx
!
!if(iel.eq.145.and.isf.eq.2)print*,iel,ier,ig,rca(isf),sca(isf),r,s,geoel(19:21,intfac(isf, ifa))

enddo
!
sfael(iel) = sfael(iel) + djaco*wi
sfael(ier) = sfael(ier) + djaco*wi
!
dpres = abs(unkng(nq+1, ig, 2)-unkng(nq+1, ig, 1))
averp = 0.5d0*abs(unkng(nq+1, ig, 2) + unkng(nq+1, ig, 1))
!
geoel(10,iel) = geoel(10, iel) + dpres/averp*djaco*wi
geoel(10,ier) = geoel(10, ier) + dpres/averp*djaco*wi
!
!if(iel.eq.145)print*,iel,ier,dpres,averp,djaco*wi,unkng(nq+1, ig, 2),unkng(nq+1, ig, 1)
!if(ier.eq.145)print*,iel,ier,dpres,averp,djaco*wi
!
enddo

endif
450 enddo
!
!...Part III: Designate shocked cell
!
do ie = 1, ncell
!
if(log10(geoel(10,ie)/sfael(ie)).ge.0.d0)then
geoel(10, ie) = 100
else
geoel(10 ,ie) = 0
endif

enddo
!
deallocate(sfael)
!
end subroutine getSIHO_facejmp3
!
!...Find the global no. of one face gauss point
!
subroutine getfgfj_glb(ipqua, ipf, fnoel)
use constant
implicit none
integer, dimension(nvqua), intent(in):: ipqua
integer, dimension(nvfac), intent(in):: ipf
integer, intent(out)::fnoel

!...Find the No. of face in fglgq...
!
if(ipqua(5).eq.ipf(3))then
!
fnoel = 1
!
elseif(ipqua(6).eq.ipf(3))then
!
fnoel = 2
!
elseif(ipqua(7).eq.ipf(3))then
!
fnoel = 3
!
elseif(ipqua(8).eq.ipf(3))then
!
fnoel = 4
!
endif

end subroutine getfgfj_glb
!
!...Get the smooth indicator from the flow divergence
!
subroutine getSI_divg(geoel, coord, unkno, ipqua, bface,intfac,esuv1, esuv2, unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),        intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in)::coord
integer, dimension(1:nvqua,1:nquad),           intent(in):: ipqua
integer*4,dimension(1:nbfai,1:nbfac),          intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)

!...Local
integer:: ipq(nvqua)
integer:: ie, ipoin,ishp,ielem
real*8,  dimension(1:2, 1:nvqua)::xpq
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8,dimension(1:nsize)   ::alfdu
real*8,dimension(1:nsize) :: unkdu
!

real*8::eps,c00,c05,c10,c20
real*8:: uctr, vctr, ectr, pctr
real*8:: rc, sc, dr, ds
real*8:: rhoct
real*8:: dudr, duds, dvdr, dvds
real*8:: dxdr, dxds, dydr, dyds, jacom
real*8:: rm, rp, sm, sp, r, s
real*8:: drdx, drdy, dsdx, dsdy
real*8:: dudx, dudy, dvdx, dvdy
real*8:: dudxy, dumax, dumin
real*8:: rhomc,sdctr,macel,delu,volel
real*8:: alfam
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
dr = 1.d0
ds = 1.d0
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, nquad
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
dudr = unkno(2, 2, ielem)/dr
duds = unkno(3, 2, ielem)/ds
dvdr = unkno(2, 3, ielem)/dr
dvds = unkno(3, 3, ielem)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rc; s=sc
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
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
!...Get dudx, dudy, dvdx, dvdy
!
!
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy

!...1st invariant
unkdu(ielem) = dudx +dvdy
enddo !...ielem

!
!...Part II: Impose limiter using the 1st invariant for mapped barth limiting coef...
!
do ie = 1, nquad
!
ielem = ie + ntria

dudxy = unkdu(ielem)
dumin = unmin_new(nq+2, ielem)
dumax = unmax_new(nq+2, ielem)

!if(ie.eq.6369) print*,'aflim',dumax/dudxy

if(dumin*dumax.gt.0.d0)then

if(dumax.lt.0.d0)then
alfdu(ielem) = max(0.d0, min(1.d0,(dumax)/(dudxy+eps)))
elseif(dumin.gt.0.d0)then
alfdu(ielem) = max(0.d0, min(1.d0,(dumin)/(dudxy+eps)))
endif

elseif(dumin*dumax.lt.0.d0)then
alfdu(ielem)=0.d0
endif
!
enddo
!
!...Part III:Further technique...
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
delu = unkdu(ielem)
!
volel = geoel(3, ielem)
macel = sqrt(volel)*abs(delu)/sdctr

alfam = min(1.d0, .03d0*macel)
!
alfdu(ielem) = alfam!*alfdu(ielem) +  (1.d0-alfam)
geoel(10, ielem) = alfdu(ielem)
!
enddo

end subroutine getSI_divg
!
!...Get the artifical stress from LLNL
!
subroutine getAS_llnl(rg ,sg, unkno, geoel ,xpq, xpqi, sgmas, ielem)
use constant
implicit none
!...Input arrays
integer, intent(in):: ielem
real*8, intent(in):: rg ,sg
real*8,dimension(1:ngeel),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq),         intent(in)::unkno
real*8,dimension(1:2, 1:2),            intent(out)::sgmas
real*8,dimension(1:2, 1:nvqua),        intent(in) ::xpq, xpqi
!...Local
integer:: ie, ipoin,ishp, ideg, idrec
real*8,  dimension(1:nvqua)::dsprq, dspsq
real*8, dimension(1:2, 1:2) :: mapmt
real*8, dimension(ndegr)    ::bg
real*8, dimension(nq)       ::unkng
real*8:: lamda(2)
real*8:: m(2,2),matsym(2,2)
!

real*8::eps,c00,c05,c10,c20
real*8:: rc, sc, dr, ds
real*8:: dudr, duds, dvdr, dvds
real*8:: dxdr, dxds, dydr, dyds, jacom
real*8:: dxdri, dxdsi, dydri, dydsi, jacomi
real*8:: a11,a12,a21,a22
real*8:: rm, rp, sm, sp, r, s
real*8:: drdx, drdy, dsdx, dsdy
real*8:: drdxi, drdyi, dsdxi, dsdyi
real*8:: dudx, dudy, dvdx, dvdy
real*8:: mapd, mapt, matra, matrb, matrc, matrd
real*8:: rhog, ug, vg, eg, pg, sdg
real*8:: unkdu
real*8:: lengs,muvis,psi0,psi1,q1,q2,rhoi

!
data eps   / 1.0d-010 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
dr = 1.d0
ds = 1.d0
!
rc = geoel(1)
sc = geoel(2)
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor..
!
dudr = unkno(2, 2)/dr
duds = unkno(3, 2)/ds
dvdr = unkno(2, 3)/dr
dvds = unkno(3, 3)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rg; s=sg
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s
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
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, 4
dxdri = dxdri + dsprq(ishp)*xpqi(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqi(1,ishp)

dydri = dydri + dsprq(ishp)*xpqi(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqi(2,ishp)
enddo
!
jacomi = dxdri*dydsi - dydri*dxdsi
drdxi = dydsi/jacomi
drdyi = -dxdsi/jacomi
dsdxi = -dydri/jacomi
dsdyi = dxdri/jacomi

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

!...Get dudx, dudy, dvdx, dvdy
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy


!...Store the divergen
unkdu = dudx + dvdy
!

!...Get the eigenvector
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy

!...eigenvalues of the matrix...
mapt  = matra + matrd
mapd  = matra*matrd - matrb*matrc
lamda(1) = 0.5d0*mapt + sqrt(0.25d0*mapt**2-mapd)
lamda(2) = 0.5d0*mapt - sqrt(0.25d0*mapt**2-mapd)
!
matsym(1, 1) = matra
matsym(1, 2) = matrb
matsym(2, 1) = matrc
matsym(2, 2) = matrd
!
call getEigensystem(matsym, mapmt, lamda,ielem)
!
if(abs(matrc).gt.1.d-8)then
!
!mapmt(1, 1) = (lamda(1)-matrd)/sqrt((lamda(1)-matrd)**2+matrc**2)
!mapmt(1, 2) =          matrc/sqrt((lamda(1)-matrd)**2+matrc**2)
!mapmt(2, 1) = (lamda(2)-matrd)/sqrt((lamda(2)-matrd)**2+matrc**2)
!mapmt(2, 2) =          matrc/sqrt((lamda(2)-matrd)**2+matrc**2)

else
!mapmt(1, 1) = 1.d0
!mapmt(1, 2) = 0.d0
!mapmt(2, 1) = 0.d0
!mapmt(2, 2) = 1.d0
endif !if(abs(matrc).gt.1.d-8)then

!
!...Part II: The general viscosity coefficient
!
bg(1) = 1.d0
bg(2) = 0.d0!(rg-rc)/dr
bg(3) = 0.d0!(sg-sc)/ds
!
unkng = 0.d0
!
do ideg = 1,mdegr
unkng(1:nq) = unkng(1:nq) + unkno(ideg,1:nq)*bg(ideg)
enddo
!
if(ndens.eq.1)then
rhog = 1.d0/unkng(1)
elseif(ndens.eq.2)then
!
rhog = unkng(1)
!

endif
ug = unkng(2)
vg = unkng(3)
eg = unkng(4)
!
pg = max(eps, (gamlg-1.d0)*rhog*(eg - 0.5d0*(ug**2 + vg**2)))
!
sdg = sqrt( max( eps,gamlg*pg/rhog) )

!
a11= dxdr*drdxi + dxds*dsdxi
a12= dxdr*drdyi + dxds*dsdyi
a21= dydr*drdxi + dyds*dsdxi
a22= dydr*drdyi + dyds*dsdyi

lengs = sqrt(geoel(4))*abs(a11*a22-a12*a21)
q2 = 1.d0/3.d0
q1 = 4.d0/2.d0
psi0 = min(abs(unkdu)/(abs(dudx*dvdy-dudy*dvdx)+eps),1.d0)


!...Zero out mpos
m = 0
do idrec = 1, 2
if(lamda(1).le.0.d0)then
psi1 = 1.d0
else
psi1 = 1.d0
endif
!
muvis = rhog*(q2*lengs**2*abs(lamda(1))+q1*psi0*psi1*lengs*sdg)!*lamda(idrec)
m(1,1) = m(1, 1) + muvis*mapmt(idrec,1)*mapmt(idrec,1)
m(1,2) = m(1, 2) + muvis*mapmt(idrec,1)*mapmt(idrec,2)
m(2,1) = m(2, 1) + muvis*mapmt(idrec,2)*mapmt(idrec,1)
m(2,2) = m(2, 2) + muvis*mapmt(idrec,2)*mapmt(idrec,2)

enddo
!
!sgmas(1,1) = m(1,1)*matra + m(1,2)*matrc
!sgmas(1,2) = m(1,1)*matrb + m(1,2)*matrd
!sgmas(2,1) = m(2,1)*matra + m(2,2)*matrc
!sgmas(2,2) = m(2,1)*matrb + m(2,2)*matrd
!
!...Set artificial vicosity equal to zero
sgmas(1,1) = muvis*matra
sgmas(1,2) = muvis*matrb
sgmas(2,1) = muvis*matrc
sgmas(2,2) = muvis*matrd
!
sgmas = 0.d0
!
if(ielem.eq.-2.or.ielem.eq.-31)then
print*,'in1+derivative',ielem,dudr,duds,dvdr,dvds
print*,'in2+Jacobian',ielem,dxdr,dxds,dydr,dyds
print*,'matrix',matra,matrb,matrc,matrd,muvis
print*,'eigen',ielem,q2,lengs**2,lamda(1:2)
endif
!
end subroutine getAS_llnl
!
!...Get the artifical stress for RSF
!
subroutine getAS(rg ,sg, unkno, geoel ,xpq, xpqi, sgmas, ielem)
use constant
implicit none
!...Input arrays
integer, intent(in):: ielem
real*8, intent(in):: rg ,sg
real*8,dimension(1:ngeel),             intent(in) ::geoel
real*8,dimension(1:ndegr,1:nq),         intent(in)::unkno
real*8,dimension(1:2, 1:2),            intent(out)::sgmas
real*8,dimension(1:2, 1:nvqua),        intent(in) ::xpq, xpqi
!...Local
integer:: ie, ipoin,ishp, ideg, idrec
real*8,  dimension(1:nvqua)::shpq,dsprq, dspsq
real*8, dimension(1:2, 1:2) :: mapmt
real*8, dimension(ndegr)    ::bg
real*8, dimension(nq)       ::unkng
real*8:: lamda(2)
real*8:: m(2,2),matsym(2,2)
!

real*8::eps,c00,c05,c10,c20
real*8:: rc, sc, dr, ds
real*8:: dudr, duds, dvdr, dvds
real*8:: dxdr, dxds, dydr, dyds, jacom
real*8:: dxdri, dxdsi, dydri, dydsi, jacomi
real*8:: a11,a12,a21,a22
real*8:: rm, rp, sm, sp, r, s
real*8:: drdx, drdy, dsdx, dsdy
real*8:: drdxi, drdyi, dsdxi, dsdyi
real*8:: dudx, dudy, dvdx, dvdy
real*8:: mapd, mapt, matra, matrb, matrc, matrd
real*8:: rhoct, uctr, vctr, ectr, pctr, sdctr
real*8:: rhog, ug, vg, eg, pg, sdg
real*8:: unkdu
real*8:: muvis

!
data eps   / 1.0d-010 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
dr = 1.d0
ds = 1.d0
!
rc = geoel(1)
sc = geoel(2)
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor..
!
dudr = unkno(2, 2)/dr
duds = unkno(3, 2)/ds
dvdr = unkno(2, 3)/dr
dvds = unkno(3, 3)/ds
!
!...Find dxdr,dxds,dydr,dyds
!
r = rg; s=sg
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
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
jacomi = dxdri*dydsi - dydri*dxdsi
drdxi = dydsi/jacomi
drdyi = -dxdsi/jacomi
dsdxi = -dydri/jacomi
dsdyi = dxdri/jacomi

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

!...Get dudx, dudy, dvdx, dvdy
jacom = dxdr*dyds - dydr*dxds
drdx = dyds/jacom
drdy = -dxds/jacom
dsdx = -dydr/jacom
dsdy = dxdr/jacom
!
dudx = dudr*drdx + duds*dsdx
dudy = dudr*drdy + duds*dsdy
dvdx = dvdr*drdx + dvds*dsdx
dvdy = dvdr*drdy + dvds*dsdy


!...Store the divergen
unkdu = dudx + dvdy
!

!...Get the eigenvector
matra = dudx
matrb = 0.5d0*(dudy + dvdx)
matrc = matrb
matrd = dvdy

!...Part II: The general viscosity coefficient

rhoct = 1.d0/unkno(1,1)
uctr = unkno(1,2)
vctr = unkno(1,3)
ectr = unkno(1,4)
!
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )
!
muvis =rhoct*sdctr*sqrt(geoel(3))/.5d0
!...Set artificial vicosity equal to zero
sgmas(1,1) = muvis*matra
sgmas(1,2) = muvis*matrb
sgmas(2,1) = muvis*matrc
sgmas(2,2) = muvis*matrd
!
if(ielem.eq.-1)then
 print*,'muvis',muvis,dudr,duds,dvdr,dvds
 print*,'mu2',dxdr,dxds,dydr,dyds
 print*,'mu3',matra,matrb,matrc,matrd
 print*,'physics',pctr,rhoct,sdctr,sqrt(geoel(3))
endif
!
!sgmas = 0.d0
!
end subroutine getAS
!
!...WENO limiter for P1 on Riemann invariant...
!
subroutine wenop1_rieminvrnt_quad_shu2(ipqua, esqua, unkno, geoel, coord, cooro, esuv1, esuv2)
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
real*8,dimension(1:ndegr,1:nq,1:nsize)::unknp
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
real*8::a11, a12, a21, a22
integer::jelaj(3)

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
dx = maxval(xpq(1, 1:nvqua)) - minval(xpq(1, 1:nvqua))
dy = maxval(xpq(2, 1:nvqua)) - minval(xpq(2, 1:nvqua))

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
if(jelem .le. ncell) then
isten = isten + 1
!
xpq_aj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpq_aj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!
dx_aj = maxval(xpq_aj(1, 1:nvqua)) - minval(xpq_aj(1, 1:nvqua))
dy_aj = maxval(xpq_aj(2, 1:nvqua)) - minval(xpq_aj(2, 1:nvqua))
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
!

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

!...left and right unknowns
unknl(1:nq) = unknp(1, 1:nq, ielem)
unknr(1:nq) = unknp(1, 1:nq, jelem)

!...Get Right eigenvectors
call getmatrix_prj(qmat, qinvm, dnx, dny, unknl, unknr, ielem)

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

do is= 1, isten
os(:, is) = unkf_cha(2, :, is)**2 + unkf_cha(3, :, is)**2
enddo

!...Linear weight
weigl(1)=0.95d0;
weigl(2:isten)= .05d0/(isten-1.d0)

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
unknp(2:3, 1:nq ,ielem) = c00
do ifa = 1, ifprj
unknp(2, 1:nq, ielem) =unknp(2, 1:nq, ielem) + weige(1:nq, ifa)*unkpf(2,1:nq,ifa)
unknp(3, 1:nq, ielem) =unknp(3, 1:nq, ielem) + weige(1:nq, ifa)*unkpf(3,1:nq,ifa)
enddo
!...  end of the loop over the quad cells
enddo
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
rhs(1:nq, 1) = unknp(2, 1:nq, ielem)*geoph(3, ielem) + unknp(3, 1:nq, ielem)*geoph(4, ielem)
rhs(1:nq, 2) = unknp(2, 1:nq, ielem)*geoph(4, ielem) + unknp(3, 1:nq, ielem)*geoph(5, ielem)
!
unkno(2, 1:nq, ielem) = mtinv(1)*rhs(1:nq, 1) + mtinv(2)*rhs(1:nq, 2)
unkno(3, 1:nq, ielem) = mtinv(3)*rhs(1:nq, 1) + mtinv(4)*rhs(1:nq, 2)
!
enddo
!
return
end subroutine wenop1_rieminvrnt_quad_shu2

!
!...Rigt vector matrix
!
subroutine getmatrix_prj(qmat, qinvm, dnx, dny, unknl,unknr, ielem)
use constant
implicit none
!...Input arrays
real*8,dimension(1:nq),   intent(in) ::unknl, unknr
real*8,dimension(1:4,1:4),intent(out)::qmat, qinvm
real*8,                    intent(in)::dnx,dny
integer, intent(in)::ielem

!...local array
real*8:: rhom1,uctr1,vctr1,ectr1,rhoc1,pctr1,sdctr1
real*8:: rhom2,uctr2,vctr2,ectr2,rhoc2,pctr2,sdctr2
real*8:: rhof, uf, vf, pf, sdf, vedn, rhosd

!...Riemann invariant
rhom1 = unknl(1)
uctr1 = unknl(2)
vctr1 = unknl(3)
ectr1 = unknl(4)
rhoc1 = 1.d0/rhom1
pctr1 = (gamlg-1.d0)*rhoc1*(ectr1-0.5d0*(uctr1**2 + vctr1**2))
sdctr1 = sqrt(gamlg*pctr1*rhom1)

rhom2 = unknr(1)
uctr2 = unknr(2)
vctr2 = unknr(3)
ectr2 = unknr(4)
rhoc2 = 1.d0/rhom2
pctr2 = (gamlg-1.d0)*rhoc2*(ectr2-0.5d0*(uctr2**2 + vctr2**2))
sdctr2 = sqrt(gamlg*pctr2*rhom2)

!
rhof = 0.5d0*(rhoc1 + rhoc2)
sdf = 0.5d0*(sdctr1 + sdctr2)
pf  = 0.5d0*(pctr1 + pctr2)
uf  = 0.5d0*(uctr1 + uctr2)
vf  = 0.5d0*(vctr1 + vctr2)
rhosd = rhof*sdf

!...Inverse matrix
qinvm = 0.d0
!
qinvm(1, 1) = -(rhosd)**2
qinvm(1, 2) = -(rhosd)*dnx
qinvm(1, 3) = -(rhosd)*dny

qinvm(2, 2) =  dny
qinvm(2, 3) = -dnx

qinvm(3, 1) =  pf
qinvm(3, 2) = -uf
qinvm(3, 3) = -vf
qinvm(3, 4) = 1.d0

qinvm(4, 1) = -(rhosd)**2
qinvm(4, 2) = (rhosd)*dnx
qinvm(4, 3) = (rhosd)*dny

!...Transfer back from Reimann invariant
qmat = 0.d0
vedn = uf*dnx + vf*dny
!
!if(ielem.eq.5)print*,'rhosd',rhof,pctr1,rhom1, sdctr2
!
qmat(1, 1) = -.5d0/(rhosd)**2
qmat(1, 4) = -.5d0/(rhosd)**2

qmat(2, 1) = -0.5d0*dnx/rhosd
qmat(2, 2) =  dny
qmat(2, 4) =  0.5d0*dnx/rhosd

qmat(3, 1) = -0.5d0*dny/rhosd
qmat(3, 2) = -dnx
qmat(3, 4) =  0.5d0*dny/rhosd

qmat(4, 1) = 0.5d0*pf/(rhosd)**2 - 0.5d0*vedn/rhosd
qmat(4, 2) = uf*dny - vf*dnx
qmat(4, 3) = 1.d0
qmat(4, 4) = 0.5d0*pf/(rhosd)**2 + 0.5d0*vedn/rhosd

end subroutine getmatrix_prj
!
!...Rigt vector matrix
!
subroutine getmatrix_prj2(qmat, qinvm, dnx, dny, unknl,unknr,ielem)
use constant
implicit none
!...Input arrays
real*8,dimension(1:nq),   intent(in) ::unknl, unknr
real*8,dimension(1:4,1:4),intent(out)::qmat, qinvm
real*8,                    intent(in)::dnx,dny
integer,   intent(in) ::ielem

!...local array
real*8:: rhom1,uctr1,vctr1,ectr1,rhoc1,pctr1,sdctr1, hctr1
real*8:: rhom2,uctr2,vctr2,ectr2,rhoc2,pctr2,sdctr2, hctr2
real*8:: rhof, uf, vf, pf, sdf, vedn, rhosd,vnm,hf, rhosq
real*8:: gam1
real*8:: eps
!
eps=1.d-6

!...Riemann invariant
rhom1 = unknl(1)
uctr1 = unknl(2)
vctr1 = unknl(3)
ectr1 = unknl(4)
rhoc1 = 1.d0/rhom1
pctr1 = max(eps,(gamlg-1.d0)*rhoc1*(ectr1-0.5d0*(uctr1**2 + vctr1**2)))
sdctr1 = sqrt(gamlg*pctr1*rhom1)

rhom2 = unknr(1)
uctr2 = unknr(2)
vctr2 = unknr(3)
ectr2 = unknr(4)
rhoc2 = 1.d0/rhom2
pctr2 = max(eps,(gamlg-1.d0)*rhoc2*(ectr2-0.5d0*(uctr2**2 + vctr2**2)))
sdctr2 = sqrt(gamlg*pctr2*rhom2)

!Simple average
rhof = 0.5d0*(rhoc1 + rhoc2)
sdf = 0.5d0*(sdctr1 + sdctr2)
pf  = 0.5d0*(pctr1 + pctr2)
uf  = 0.5d0*(uctr1 + uctr2)
vf  = 0.5d0*(vctr1 + vctr2)

!Roe average
rhof = sqrt(rhoc1*rhoc2)
rhosq = sqrt(rhoc1) + sqrt(rhoc2)
uf  = (sqrt(rhoc1)*uctr1 + sqrt(rhoc2)*uctr2)/rhosq
vf  = (sqrt(rhoc1)*vctr1 + sqrt(rhoc2)*vctr2)/rhosq
hctr1 = max(eps,gamlg*ectr1 - (gamlg-1.d0)*0.5d0*(uctr1**2 + vctr1**2))
hctr2 = max(eps,gamlg*ectr2 - (gamlg-1.d0)*0.5d0*(uctr2**2 + vctr2**2))
hf = (sqrt(rhoc1)*hctr1 + sqrt(rhoc2)*hctr2)/rhosq
sdf = sqrt(max(eps, (gamlg-1.d0)*(hf-0.5d0*(uf**2+vf**2))))
rhosd = rhof*sdf
pf = sdf**2*rhof/gamlg

!...Constant
gam1 = (gamlg-1.d0)/gamlg
vnm  = vf*dnx - uf*dny

!...Inverse matrix
qinvm = 0.d0
!
qinvm(1, 1) = -.5d0/gamlg
qinvm(1, 2) =  .5d0*dnx/rhosd - 0.5d0*uf/pf*gam1
qinvm(1, 3) =  .5d0*dny/rhosd - 0.5d0*vf/pf*gam1
qinvm(1, 4) =  .5d0/pf*gam1

qinvm(2, 1) =  1.d0/gamlg
qinvm(2, 2) = -dny*vnm/pf - uf/pf/gamlg
qinvm(2, 3) =  dnx*vnm/pf - vf/pf/gamlg
qinvm(2, 4) =  1.d0/pf/gamlg

qinvm(3, 1) =  0.d0
qinvm(3, 2) =  dny/pf
qinvm(3, 3) = -dnx/pf
qinvm(3, 4) =  0.d0

qinvm(4, 1) =  .5d0/gamlg
qinvm(4, 2) =  .5d0*dnx/rhosd + 0.5d0*uf/pf*gam1
qinvm(4, 3) =  .5d0*dny/rhosd + 0.5d0*vf/pf*gam1
qinvm(4, 4) = -.5d0/pf*gam1

!...Transfer back from Reimann invariant
qmat = 0.d0
vedn = uf*dnx + vf*dny
!
qmat(1, 1) = -1.d0
qmat(1, 2) =  gamlg - 1.d0
qmat(1, 3) =  (gamlg - 1.d0)*vnm
qmat(1, 4) =  1.d0

qmat(2, 1) =  rhosd*dnx
qmat(2, 2) =  0.d0
qmat(2, 3) =  pf*dny
qmat(2, 4) =  rhosd*dnx

qmat(3, 1) =  rhosd*dny
qmat(3, 2) =  0.d0
qmat(3, 3) = -pf*dnx
qmat(3, 4) =  rhosd*dny

qmat(4, 1) =  pf + rhosd*vedn
qmat(4, 2) =  pf
qmat(4, 3) =  0.d0
qmat(4, 4) = -pf + rhosd*vedn

!if(ielem.eq.71)print*,'inside',hf,0.5d0*(uf**2+vf**2)

end subroutine getmatrix_prj2

!
!...WENO limiter for P1 on Riemann invariant...
!
subroutine wenop1_unkstcl_quad(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoph, coord, cooro)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:5,1:nsize), intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),             intent(in) :: coord, cooro


integer, intent(in)::nsten, ielem
integer, intent(out)::isten

integer:: ipq(nvqua)
integer:: mapfe(1:2,1:nfqua)
!
real*8,dimension(1:2, 1:nvqua)::xpq, xpqi,xpq_aj, xpqi_aj

real*8,dimension(1:4):: rhsmt_aj, lhsmt_aj,mtinv
real*8,dimension(1:nq,1:2)::rhs
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:ndegr)::bq,bqp, bqp_aj

real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: wi
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
real*8:: rhom1,uctr1,vctr1,ectr1,rhoc1,pctr1,sdctr1
real*8:: rhom2,uctr2,vctr2,ectr2,rhoc2,pctr2,sdctr2
real*8:: rhof, uf, vf, pf, sdf, vedn, rhosd
real*8:: rhogi
!
integer:: ie, ies, iq, is, ntrsf,  jelem, irm, igaus, ishp
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

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...Scaling parameters
dr   = 1.d0
ds   = 1.d0
!
 ie = ielem-ntria
!
!...Part II: Recontruction
!
!...physical mass center
xmc =geoph(1, ielem)
ymc =geoph(2, ielem)

xpq(1, 1:nvqua) = coord(1, ipqua(1:nvqua,ielem))
xpq(2, 1:nvqua) = coord(2, ipqua(1:nvqua,ielem))
!
xpqi(1, 1:nvqua) = cooro(1, ipqua(1:nvqua,ielem))
xpqi(2, 1:nvqua) = cooro(2, ipqua(1:nvqua,ielem))

!...scale
dx = 1.d0!maxval(xpq(1, 1:nvqua)) - minval(xpq(1, 1:nvqua))
dy = 1.d0!maxval(xpq(2, 1:nvqua)) - minval(xpq(2, 1:nvqua))

!...  b. curvatures for the face-neighboring cells
!isten = 1
!unkpe = 0.d0
!
do ies = 1, 4

jelem = esqua(ies,ie)
if(jelem .le. ncell) then
!
isten = isten + 1
!
xpq_aj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpq_aj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!
xpqi_aj(1, 1:nvqua) = cooro(1, ipqua(1:nvqua,jelem))
xpqi_aj(2, 1:nvqua) = cooro(2, ipqua(1:nvqua,jelem))

!...physical mass center
xmc_aj =geoph(1, jelem)
ymc_aj =geoph(2, jelem)
!
dx_aj = 1.d0!maxval(xpq_aj(1, 1:nvqua)) - minval(xpq_aj(1, 1:nvqua))
dy_aj = 1.d0!maxval(xpq_aj(2, 1:nvqua)) - minval(xpq_aj(2, 1:nvqua))

!...Initialize
rhsmt_aj = 0.d0
lhsmt_aj = 0.d0

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
dxdri = dxdri + dsprq(ishp)*xpqi_aj(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqi_aj(1,ishp)

dydri = dydri + dsprq(ishp)*xpqi_aj(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqi_aj(2,ishp)
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
xgausi = xgausi + shpq(ishp)*xpqi_aj(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqi_aj(2,ishp)

xgaus  = xgaus + shpq(ishp)*xpq_aj(1,ishp)
ygaus  = ygaus + shpq(ishp)*xpq_aj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqp(1) = 1.d0
bqp(2) = (xgaus-xmc)/dx
bqp(3) = (ygaus-ymc)/dy

bqp_aj(1) = 1.d0
bqp_aj(2) = (xgaus-xmc_aj)/dx_aj
bqp_aj(3) = (ygaus-ymc_aj)/dy_aj

!...Matrix
rhsmt_aj(1) = rhsmt_aj(1) + rhogi*djacoi*bqp_aj(2)*bqp_aj(2)
rhsmt_aj(2) = rhsmt_aj(2) + rhogi*djacoi*bqp_aj(3)*bqp_aj(2)
rhsmt_aj(3) = rhsmt_aj(3) + rhogi*djacoi*bqp_aj(2)*bqp_aj(3)
rhsmt_aj(4) = rhsmt_aj(4) + rhogi*djacoi*bqp_aj(3)*bqp_aj(3)
!
lhsmt_aj(1) =  lhsmt_aj(1) + rhogi*djacoi*bqp(2)*bqp_aj(2)
lhsmt_aj(2) =  lhsmt_aj(2) + rhogi*djacoi*bqp(3)*bqp_aj(2)
lhsmt_aj(3) =  lhsmt_aj(3) + rhogi*djacoi*bqp(2)*bqp_aj(3)
lhsmt_aj(4) =  lhsmt_aj(4) + rhogi*djacoi*bqp(3)*bqp_aj(3)
enddo

!...Matrix inverse
detma = lhsmt_aj(1)*lhsmt_aj(4) - lhsmt_aj(2)*lhsmt_aj(3)
mtinv(1) = lhsmt_aj(4)
mtinv(2) =-lhsmt_aj(2)
mtinv(3) =-lhsmt_aj(3)
mtinv(4) = lhsmt_aj(1)
!
mtinv = mtinv/detma

!
unkpe(1, 1:nq, isten) = unknp(1, 1:nq, ielem)

rhs(1:nq, 1) = rhsmt_aj(1)*unknp(2, 1:nq, jelem) + rhsmt_aj(2)*unknp(3, 1:nq, jelem)
rhs(1:nq, 2) = rhsmt_aj(3)*unknp(2, 1:nq, jelem) + rhsmt_aj(4)*unknp(3, 1:nq, jelem)

unkpe(2, 1:nq, isten) = mtinv(1)*rhs(1:nq, 1) +  mtinv(2)*rhs(1:nq, 2)
unkpe(3, 1:nq, isten) = mtinv(3)*rhs(1:nq, 1) +  mtinv(4)*rhs(1:nq, 2)
endif
enddo
!
return
end subroutine wenop1_unkstcl_quad
!
!...WENO limiter for P1 using stencil 2...
!
subroutine wenop1_unkstcl_quad2(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoph,&
coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:5,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)

integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(5)::jelaj
integer,  dimension(1:nvqua):: ipq
!
real*8,dimension(1:2, 1:nvqua)::xpq
real*8, dimension(1:3, 1:2)::lhsmt
real*8, dimension(1:3,1:nq)::rhsmt
real*8,dimension(1:4):: mtinv, lhsls
real*8,dimension(1:nq,1:2)::rhsls

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dx, dy, xmc, ymc
real*8:: detma
!
integer:: ie, ies, jelem,  ishp, ieaj, istor,iv,neaj

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Scaling parameter
dx = 1.d0!maxval(xpq(1, 1:nvqua)) - minval(xpq(1, 1:nvqua))
dy = 1.d0!maxval(xpq(2, 1:nvqua)) - minval(xpq(2, 1:nvqua))

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Loop over the 4 vertices for one quad
do iv = 1, 4

jelaj = 0
neaj = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)

!
!if(ielem.eq.900)then
!print*,'iv',iv,istor,jelem,esqua(:,ielem)
!endif

if(jelem.ne.ielem)then
neaj = neaj +1
jelaj(neaj) = jelem
endif
enddo !istor

!...For neaj lsee that 2, there is no stencil.
if(neaj.eq.3)then

!...Stencil No.
isten = isten + 1

!...LHS matrix
lhsmt(1, 1) = geoph(1, jelaj(1)) - xmc; lhsmt(1, 2) = geoph(2, jelaj(1)) - ymc
lhsmt(2, 1) = geoph(1, jelaj(2)) - xmc; lhsmt(2, 2) = geoph(2, jelaj(2)) - ymc
lhsmt(3, 1) = geoph(1, jelaj(3)) - xmc; lhsmt(3, 2) = geoph(2, jelaj(3)) - ymc

!...Scale
lhsmt(:, 1) = lhsmt(:, 1)/dx
lhsmt(:, 2) = lhsmt(:, 2)/dy

!
rhsmt(1, 1:nq) = unknp(1, 1:nq, jelaj(1)) - unknp(1, 1:nq, ielem)
rhsmt(2, 1:nq) = unknp(1, 1:nq, jelaj(2)) - unknp(1, 1:nq, ielem)
rhsmt(3, 1:nq) = unknp(1, 1:nq, jelaj(3)) - unknp(1, 1:nq, ielem)
!
lhsls = 0.d0
rhsls = 0.d0
!
do ieaj=1,3
lhsls(1) = lhsls(1) + lhsmt(ieaj, 1)**2
lhsls(2) = lhsls(2) + lhsmt(ieaj, 1)*lhsmt(ieaj, 2)
lhsls(3) = lhsls(3) + lhsmt(ieaj, 1)*lhsmt(ieaj, 2)
lhsls(4) = lhsls(4) + lhsmt(ieaj, 2)**2
!
rhsls(1:nq, 1) = rhsls(1:nq, 1) + lhsmt(ieaj, 1)*rhsmt(ieaj, 1:nq)
rhsls(1:nq, 2) = rhsls(1:nq, 2) + lhsmt(ieaj, 2)*rhsmt(ieaj, 1:nq)
enddo

!...Inverse matrix
detma = lhsls(1)*lhsls(4) - lhsls(2)*lhsls(3)
mtinv(1) = lhsls(4)
mtinv(2) =-lhsls(2)
mtinv(3) =-lhsls(3)
mtinv(4) = lhsls(1)

mtinv = mtinv /detma

unkpe(1, 1:nq, isten) = unknp(1, 1:nq, ielem)
unkpe(2, 1:nq, isten) = mtinv(1)*rhsls(1:nq, 1) + mtinv(2)*rhsls(1:nq, 2)
unkpe(3, 1:nq, isten) = mtinv(3)*rhsls(1:nq, 1) + mtinv(4)*rhsls(1:nq, 2)
endif
enddo !do iv = 1, 4
!
return
end subroutine wenop1_unkstcl_quad2
!
!...WENO limiter for P1 using stencil 2...
!
subroutine wenop1_unkstcl_tria2(isten, nsten, ielem, iptri, estri, unkpe, unknp, geoph,&
coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer, dimension(1:nftri,1:ntria),  intent(in)::estri
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(out)::unkpe
real*8,dimension(1:5,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)

integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(8)::jelaj
integer,  dimension(1:nvtri):: ipt
!
real*8,dimension(1:2, 1:nvtri)::xpt
real*8, dimension(1:8, 1:2)::lhsmt
real*8, dimension(1:8,1:nq)::rhsmt
real*8,dimension(1:4):: mtinv, lhsls
real*8,dimension(1:nq,1:2)::rhsls

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dx, dy, xmc, ymc
real*8:: detma
!
integer:: ie, ies, jelem,  ishp, ieaj, istor,iv,neaj

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Local triangle No.
ie = ielem

!...Local vertices
ipt(1:nvtri) = iptri(1:nvtri,ie)

!...Coordinates
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))

!...Scaling parameter
dx = maxval(xpt(1, 1:nvtri)) - minval(xpt(1, 1:nvtri))
dy = maxval(xpt(2, 1:nvtri)) - minval(xpt(2, 1:nvtri))

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Loop over the 4 vertices for one quad
do iv = 1, 3

jelaj = 0
neaj = 0
do istor=esuv2(ipt(iv))+1,esuv2(ipt(iv)+1)
jelem=esuv1(istor)
!
!if(ielem.eq.5) print*,'ielem5',jelem,istor,iv

if(jelem.ne.ielem)then
neaj = neaj +1
jelaj(neaj) = jelem
endif
enddo !istor

!...For neaj lsee that 2, there is no stencil.
if(neaj.ge.3)then

!...Stencil No.
isten = isten + 1

!...LHS matrix
do ieaj=1,neaj
lhsmt(ieaj, 1) = geoph(1, jelaj(ieaj)) - xmc;
lhsmt(ieaj, 2) = geoph(2, jelaj(ieaj)) - ymc
enddo

!...Scale
lhsmt(:, 1) = lhsmt(:, 1)/dx
lhsmt(:, 2) = lhsmt(:, 2)/dy

!
do ieaj=1,neaj
rhsmt(ieaj, 1:nq) = unknp(1, 1:nq, jelaj(ieaj)) - unknp(1, 1:nq, ielem)
enddo
!
lhsls = 0.d0
rhsls = 0.d0
!
do ieaj=1,neaj
lhsls(1) = lhsls(1) + lhsmt(ieaj, 1)**2
lhsls(2) = lhsls(2) + lhsmt(ieaj, 1)*lhsmt(ieaj, 2)
lhsls(3) = lhsls(3) + lhsmt(ieaj, 1)*lhsmt(ieaj, 2)
lhsls(4) = lhsls(4) + lhsmt(ieaj, 2)**2
!
rhsls(1:nq, 1) = rhsls(1:nq, 1) + lhsmt(ieaj, 1)*rhsmt(ieaj, 1:nq)
rhsls(1:nq, 2) = rhsls(1:nq, 2) + lhsmt(ieaj, 2)*rhsmt(ieaj, 1:nq)
enddo

!...Inverse matrix
detma = lhsls(1)*lhsls(4) - lhsls(2)*lhsls(3)
mtinv(1) = lhsls(4)
mtinv(2) =-lhsls(2)
mtinv(3) =-lhsls(3)
mtinv(4) = lhsls(1)

mtinv = mtinv /detma

unkpe(1, 1:nq, isten) = unknp(1, 1:nq, ielem)
unkpe(2, 1:nq, isten) = mtinv(1)*rhsls(1:nq, 1) + mtinv(2)*rhsls(1:nq, 2)
unkpe(3, 1:nq, isten) = mtinv(3)*rhsls(1:nq, 1) + mtinv(4)*rhsls(1:nq, 2)
endif
enddo !do iv = 1, 4
!
return
end subroutine wenop1_unkstcl_tria2

!
!...WENO limiter for P1 on Riemann invariant...
!
subroutine wenop1_rieminvrnt_quad_shu3(ipqua, esqua, unkno, geoel, coord, cooro, esuv1, esuv2)
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
mapmt(1, 1, ielem) = 1.d0
mapmt(1, 2, ielem) = 0.d0
mapmt(2, 1, ielem) = 0.d0
mapmt(2, 2, ielem) = 1.d0
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

!mapmt(1, 1, ielem) = 1.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 1.d0

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

do is= 1, isten
os(:, is) = unkf_cha(2, :, is)**2 + unkf_cha(3, :, is)**2
enddo

!...Linear weight
weigl(1)=0.8d0;
weigl(2:isten)= 0.2d0/(isten-1.d0)

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
end subroutine wenop1_rieminvrnt_quad_shu3
!
!...WENO limiter for P1 on Riemann invariant...
!
subroutine wenop1_rieminvrnt_tria_shu(iptri, estri, unkno, geoel, coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer, dimension(1:nftri,1:ntria),  intent(in)::estri
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) :: coord, cooro
real*8,dimension(1:ngeel,1:nsize),             intent(in) ::geoel
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
!
!...  local arrays
!
integer, parameter:: nsten=8,nfprj=8
integer:: ipt(nvtri)
integer:: mapfe(1:2,1:nftri)
!
real*8,dimension(1:2, 1:nvtri)::xpt, xpti, xpt_aj, xpti_aj
real*8,dimension(1:ndegr,1:nq,1:nsize)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten)::unkpf,unkf_cha
real*8,dimension(1:ndegr,1:nq,1:nsten)::unkrm,unkpe


real*8,dimension(1:5,1:nsize) ::geoph
real*8,dimension(1:4,1:nsize) ::rhsmt
real*8,dimension(1:4):: rhsmt_aj, lhsmt_aj,mtinv
real*8,dimension(1:nq,1:2)::rhs
real*8,dimension(1:4, 1:4)::qmat, qinvm
real*8,dimension(1:nvtri)::shpt, dsprt, dspst
real*8,dimension(1:ndegr)::bt,btp, btp_aj
real*8::xpf(1:2, 1:2)
real*8, dimension(1:nq)::unknl,unknr
real*8:: weight(ngausd), posit(2, ngausd)
real*8:: wi
real*8   weigh(1:nq, nsten)
real*8:: weigt(1:nq)
real*8   os(1:nq, nsten)
real*8:: weigl(nsten)
real*8:: weige(nq, nfprj)
real*8::lhsls(3, 2), rhsls(1:3, 1:nq),b1(nq),b2(nq)
real*8::a11, a12, a21, a22
integer::jelaj(3)

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
data rpowe / 10.d0    /   ! works for lilia case
!
!...Part 0: Basis parameters setup
!
!...Specify the weno
nweno =2

!...Find weight and position for gauss points...
call rutope(2, ngausd, posit, weight)

!...Scaling parameters
dr   = .5d0
ds   = .5d0

!...Mapping array
mapfe(1, 1) = 1; mapfe(2, 1) = 2;
mapfe(1, 2) = 2; mapfe(2, 2) = 3;
mapfe(1, 3) = 3; mapfe(2, 3) = 1;
!
!...Part I: Get the L2 projection matrix
!
unknp = unkno

!...I.1:Get the physical mass center xmc, ymc
do ie = 1, ntria
!
ielem = ie
!
!if(ie.eq.1)print*,'tyia',estri(:,ie)
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))

xpti(1, 1:nvtri) = cooro(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = cooro(2, ipt(1:nvtri))
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Cell cenetr for density
!call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpti, rc, sc, xcrho, ycrho)

!...mass
masel = geoel(4, ielem)

!...physical mass center
xmc = 0.d0
ymc = 0.d0
!
do igaus =1,ngausd
!
r  = posit(1,igaus)
s  = posit(2,igaus)
wi  = weight(igaus)

!...  shape function & its derivatives w.r.t. reference coordinates
!call getshapfct_tria(ncurv,nvtri,shpt, dsprt, dspst, r, s)

shpt(1) = 1.d0-r-s
shpt(2) = r
shpt(3) = s
!
dsprt(1) = -c10
dsprt(2) =  c10
dsprt(3) =  0.d0
!
dspst(1) = -c10
dspst(2) =  0.d0
dspst(3) =  c10

!...Current domain
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvtri
dxdr = dxdr + dsprt(ishp)*xpt(1,ishp)
dxds = dxds + dspst(ishp)*xpt(1,ishp)

dydr = dydr + dsprt(ishp)*xpt(2,ishp)
dyds = dyds + dspst(ishp)*xpt(2,ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)

!...Initial domain
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, nvtri
dxdri = dxdri + dsprt(ishp)*xpti(1,ishp)
dxdsi = dxdsi + dspst(ishp)*xpti(1,ishp)

dydri = dydri + dsprt(ishp)*xpti(2,ishp)
dydsi = dydsi + dspst(ishp)*xpti(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)*0.5d0

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgausi = xgausi + shpt(ishp)*xpti(1,ishp)
ygausi = ygausi + shpt(ishp)*xpti(2,ishp)

xgaus  = xgaus + shpt(ishp)*xpt(1,ishp)
ygaus  = ygaus + shpt(ishp)*xpt(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

rhogi =1.d0

!...physical mass center
xmc = xmc + rhogi*djacoi*xgaus
ymc = ymc + rhogi*djacoi*ygaus
enddo

!...Physical mass center
geoph(1, ielem) = xmc/masel
geoph(2, ielem) = ymc/masel
enddo !do ie = 1, nquad
!...I.2: L2 projection matrix

rhsmt= 0.d0
geoph(3:5, :) = 0.d0
!
do ie = 1, ntria
!
ielem = ie
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))

xpti(1, 1:nvtri) = cooro(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = cooro(2, ipt(1:nvtri))
!...mass center...
rc= geoel(1, ielem)
sc= geoel(2, ielem)

!...physical mass center
xmc =geoph(1, ielem)
ymc =geoph(2, ielem)

!...scale
dx = maxval(xpt(1, 1:nvtri)) - minval(xpt(1, 1:nvtri))
dy = maxval(xpt(2, 1:nvtri)) - minval(xpt(2, 1:nvtri))

!...Cell cenetr for density
!call GetCellctr_quad_initial (ncurv,ndimn,nvtri,xpti, rc, sc, xcrho, ycrho)

!...mass
masel = geoel(4, ielem)

!
do igaus =1,ngausd
!
r  = posit(1,igaus)
s  = posit(2,igaus)
wi  = weight(igaus)

!...  shape function & its derivatives w.r.t. reference coordinates
!call getshapfct_quad(ncurv,nvtri,shpt, dsprt, dspst, r, s)

shpt(1) = 1.d0-r-s
shpt(2) = r
shpt(3) = s
!
dsprt(1) = -c10
dsprt(2) =  c10
dsprt(3) =  0.d0
!
dspst(1) = -c10
dspst(2) =  0.d0
dspst(3) =  c10

!...Initial domain
dxdri = 0.d0
dxdsi = 0.d0
dydri = 0.d0
dydsi = 0.d0
!
do ishp = 1, nvtri
dxdri = dxdri + dsprt(ishp)*xpti(1,ishp)
dxdsi = dxdsi + dspst(ishp)*xpti(1,ishp)

dydri = dydri + dsprt(ishp)*xpti(2,ishp)
dydsi = dydsi + dspst(ishp)*xpti(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)*0.5d0

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, nvtri
xgausi = xgausi + shpt(ishp)*xpti(1,ishp)
ygausi = ygausi + shpt(ishp)*xpti(2,ishp)

xgaus  = xgaus + shpt(ishp)*xpt(1,ishp)
ygaus  = ygaus + shpt(ishp)*xpt(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

rhogi = 1.d0

!...Basis functions
bt(1) = 1.d0
bt(2) = (r-rc)/dr
bt(3) = (s-sc)/ds

btp(1) = 1.d0
btp(2) = (xgaus-xmc)/dx
btp(3) = (ygaus-ymc)/dy

!...Matrix
rhsmt(1, ielem) = rhsmt(1, ielem) + rhogi*djacoi*bt(2)*btp(2)
rhsmt(2, ielem) = rhsmt(2, ielem) + rhogi*djacoi*bt(3)*btp(2)
rhsmt(3, ielem) = rhsmt(3, ielem) + rhogi*djacoi*bt(2)*btp(3)
rhsmt(4, ielem) = rhsmt(4, ielem) + rhogi*djacoi*bt(3)*btp(3)
!
geoph(3, ielem) = geoph(3, ielem) + rhogi*djacoi*btp(2)*btp(2)
geoph(4, ielem) = geoph(4, ielem) + rhogi*djacoi*btp(2)*btp(3)
geoph(5, ielem) = geoph(5, ielem) + rhogi*djacoi*btp(3)*btp(3)
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

enddo !do ie = 1, nquad
!
!...Part II: Recontruction
!
do ie = 1,ntria
!
ielem = ie
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))

!...physical mass center
xmc =geoph(1, ielem)
ymc =geoph(2, ielem)

!...scale
dx = maxval(xpt(1, 1:nvtri)) - minval(xpt(1, 1:nvtri))
dy = maxval(xpt(2, 1:nvtri)) - minval(xpt(2, 1:nvtri))

!...  b. curvatures for the face-neighboring cells
unkpe = 0.d0

select Case (nweno)
Case (1)

case(2)

!...Stencil 2
isten = 1
unkpe(:, 1:nq, isten) = unknp(:, 1:nq, ielem)

!call wenop1_unkstcl_quad2(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoph,&
!coord, cooro, esuv1, esuv2)

call wenop1_unkstcl_tria2(isten, nsten, ielem, iptri, estri, unkpe, unknp, geoph,&
coord, cooro, esuv1, esuv2)

case(3)

end select

!if(ielem.eq.5)print*,'tria',ifprj,isten,unkpe(2,:,1),unknp(2, :, ielem)

!...  b. curvatures for the face-neighboring cells
ifprj = 0
unkpf = 0.d0
unkrm = 0.d0
!
do ies = 1, 3

jelem = estri(ies,ie)
if(jelem .le. ncell) then

!...Find the projection vector
ifprj = ifprj +1

!
!if(ie.eq.5) print*,'triabefore',ies,jelem,ncell,ifprj

!...normal vector
xpf(1, 1:2) = xpt(1, mapfe(1:2, ies))
xpf(2, 1:2) = xpt(2, mapfe(1:2, ies))

dtx = xpf(1, 2) - xpf(1, 1)
dty = xpf(2, 2) - xpf(2, 1)

dlgt = sqrt(dtx**2 + dty**2)

dnx = dty/dlgt
dny =-dtx/dlgt

!...left and right unknowns
unknl(1:nq) = unknp(1, 1:nq, ielem)
unknr(1:nq) = unknp(1, 1:nq, jelem)

!...Get Right eigenvectors
call getmatrix_prj(qmat, qinvm, dnx, dny, unknl, unknr,ielem)


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

!if(ielem.eq.5)print*,'tria',ifprj,isten,unkpe(2,:,1),unkpe(2,:,2),unkpe(2,:,3)

!...Smooth indicator
os = 0.d0

do is= 1, isten
os(:, is) = unkf_cha(2, :, is)**2 + unkf_cha(3, :, is)**2
enddo

!...Linear weight
weigl(1)=0.00d0;
weigl(2:isten)= 1.d0/(isten-1.d0)

!weigl = 1.d0
!weigl(2:isten)= 0.d0

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
unknp(2:3, 1:nq ,ielem) = c00
do ifa = 1, ifprj
unknp(2, 1:nq, ielem) =unknp(2, 1:nq, ielem) + weige(1:nq, ifa)*unkpf(2,1:nq,ifa)
unknp(3, 1:nq, ielem) =unknp(3, 1:nq, ielem) + weige(1:nq, ifa)*unkpf(3,1:nq,ifa)
enddo
!
!if(ielem.eq.5)print*,'iel54',ielem,weige(2, 1:2),unkpf(2,2,1:2)
!...  end of the loop over the quad cells
enddo

!print*,'tria',unknp(:,2, 5)
!
!...Get the final reference polynomial
!
do ie = 1, ntria
!
ielem = ie

detma = rhsmt(1, ielem)*rhsmt(4, ielem) - rhsmt(2, ielem)*rhsmt(3, ielem)
mtinv(1) = rhsmt(4, ielem)
mtinv(2) =-rhsmt(2, ielem)
mtinv(3) =-rhsmt(3, ielem)
mtinv(4) = rhsmt(1, ielem)
!
mtinv = mtinv/detma
!
rhs(1:nq, 1) = unknp(2, 1:nq, ielem)*geoph(3, ielem) + unknp(3, 1:nq, ielem)*geoph(4, ielem)
rhs(1:nq, 2) = unknp(2, 1:nq, ielem)*geoph(4, ielem) + unknp(3, 1:nq, ielem)*geoph(5, ielem)
!
unkno(2, 1:nq, ielem) = mtinv(1)*rhs(1:nq, 1) + mtinv(2)*rhs(1:nq, 2)
unkno(3, 1:nq, ielem) = mtinv(3)*rhs(1:nq, 1) + mtinv(4)*rhs(1:nq, 2)
!
enddo
!
return
end subroutine wenop1_rieminvrnt_tria_shu
!
!...WENO limiter for P1 on Riemann invariant...
!
subroutine wenop1_rieminvrnt_quad_shu6(ipqua, esqua, unkno, geoel, coord, cooro, esuv1, esuv2)
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
integer, parameter:: nsten=10,nfprj=10
integer:: ipq(nvqua)
integer:: mapfe(1:2,1:nfqua)
!
real*8,dimension(1:2, 1:nvqua)::xpq, xpqi, xpq_aj, xpqi_aj
real*8,dimension(1:ndegr,1:nq,1:nsize)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten)::unkf_phy,unksc, unksl,unksf_cha, unkf_cha
real*8,dimension(1:ndegr,1:nq)::unk_phy


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
integer:: ie, ies, iq, is, isten, ntrsf, ielem, jelem, ifa, igaus, ishp,ifprj,neaj,ieaj,istor,iv,id
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

!
!print*,'esqua',esqua(:,1),esqua(:,900),esqua(:,600)

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
dx =1.d0
dy =1.d0

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

!if(ielem.eq.42)print*,'ielem',rhs(1,1:2),unkno(2:3, 1, ielem),rhsmt(:,ielem)

!...Local system
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
if(sqrt(matra**2+matrb**2).gt.1.d-10)then
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

!...No local system
!mapmt(1, 1, ielem) = 1.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 1.d0

enddo !do ie = 1, nquad
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
dx =1.d0
dy =1.d0

!...  b. curvatures for the face-neighboring cells
unksc = 0.d0 !...Unknown for Cartesian stencile

select Case (nweno)

Case (1)

!...Hermite weno
isten = 1
unksc(:, 1:nq, isten) = unknp(:, 1:nq, ielem)
!
do ies = 1, 4
jelem = esqua(ies,ie)

if(jelem .le. ncell) then
isten = isten + 1
unksc(1, 1:nq, isten) = unknp(1, 1:nq, ielem)
unksc(2, 1:nq, isten) = unknp(2, 1:nq, jelem)
unksc(3, 1:nq, isten) = unknp(3, 1:nq, jelem)
endif
enddo !do ies = 1, 4

!call wenop1_unkstcl_quad2(isten, nsten, ielem, ipqua, esqua, unksc, unknp, geoph,&
!coord, cooro, esuv1, esuv2)

!...Stencil 2
case(2)
isten = 1
unksc(:, 1:nq, isten) = unknp(:, 1:nq, ielem)

call wenop1_unkstcl_quad2(isten, nsten, ielem, ipqua, esqua, unksc, unknp, geoph,&
coord, cooro, esuv1, esuv2)

case(3)

!...Stencil 3
isten = 1
unksc(:, 1:nq, isten) = unknp(:, 1:nq, ielem)

!...Hermite
call wenop1_unkstcl_quad(isten, nsten, ielem, ipqua, esqua, unksc, unknp, geoph, coord, cooro)

!...Lagrangian
call wenop1_unkstcl_quad2(isten, nsten, ielem, ipqua, esqua, unksc, unknp, geoph,&
coord, cooro, esuv1, esuv2)

end select

!
!if(ielem.eq.1153.or.ielem.eq.1154.or.ielem.eq.1176)then
!print*,'Bef',ielem,unksc(2:3,1,1:isten)
!endif

!...Local system
call getunkno_local(0, isten, mapmt(:,:,ielem), unksc(:,:,1:isten))

!
!if(ielem.eq.1153.or.ielem.eq.1154.or.ielem.eq.1176)then
! print*,'bad',ielem,unksc(2:3,1,1:isten)
!endif

!...Update local stencile values
unksl = unksc

!...  b. curvatures for the face-neighboring cells
ifprj = 0
unkf_cha = 0.d0
unkf_phy = 0.d0

!...Loop over the faces
do ies = 1, 1!4

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
call getvector_local(mapmt(:,:,ielem), dnx,dny)

!...left and right unknowns
unknl(1:nq) = unknp(1, 1:nq, ielem)
unknr(1:nq) = unknp(1, 1:nq, jelem)
!
call getvector_local(mapmt(:,:,ielem), unknl(2), unknl(3))
call getvector_local(mapmt(:,:,ielem), unknr(2), unknr(3))

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
!call getmatrix_prj2(qmat, qinvm, dnx, dny, unknl,unknr,ielem)

!...Riemann invariant
unksf_cha = 0.d0
do is= 1, isten
do iq= 1, nq
do id = 1, ndegr
unksf_cha(id, 1, is) =  unksf_cha(id, 1, is) + qinvm(1, iq)*unksl(id,iq, is)
unksf_cha(id, 2, is) =  unksf_cha(id, 2, is) + qinvm(2, iq)*unksl(id,iq, is)
unksf_cha(id, 3, is) =  unksf_cha(id, 3, is) + qinvm(3, iq)*unksl(id,iq, is)
unksf_cha(id, 4, is) =  unksf_cha(id, 4, is) + qinvm(4, iq)*unksl(id,iq, is)
enddo
enddo
enddo

!...Smooth indicator
os = 0.d0

!...scale
dxc = maxval(xpq(1, 1:nvqua)) - minval(xpq(1, 1:nvqua))
dyc = maxval(xpq(2, 1:nvqua)) - minval(xpq(2, 1:nvqua))

!
call getvector_local(mapmt(:,:,ielem), dxc,dyc)

dxc = 1.d0
dyc = 1.d0

!...Smooth indicator for every stencle
do is= 1, isten
do iq= 1, nq
!os(iq, is) = ((unkf_cha(2, iq, is)*dxc)**2 + (unkf_cha(3, iq, is)*dyc)**2)
os(iq, is) = ((unksf_cha(2, iq, is)*dxc)**2 +&
(unksf_cha(3, iq, is)*dyc)**2)!/(abs(unksf_cha(1, iq, 1))+epsil)**1
enddo
enddo

!...Mometum
!os1(:) = 0.25d0*(os(1, :)+os(2, :)+os(3, :)+os(4, :))
os1(:) = max(os(1, :),os(2, :),os(3, :),os(4, :))


do is= 1, isten
os(1:nq, is) = os1(is)
enddo

!...Linear weight for central and biased stencle
weigl(1)=0.2d0;
weigl(2:isten)= 0.8d0/(isten-1.d0)


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
unkf_cha(:,iq,ifprj) = unkf_cha(:,iq, ifprj) + weigh(iq,is)*unksf_cha(:, iq, is)
enddo
enddo

!++++++++++Print for debugging++++++++++
!if(ielem.eq.1153.or.ielem.eq.1154.or.ielem.eq.1176)then
!print*,'Weight',ielem,weigh(1,1:isten)
!endif

!...Only update the dirivatives
do iq = 1, nq
do id = 2,ndegr
unkf_phy(id, 1, ifprj) =  unkf_phy(id, 1, ifprj) + qmat(1, iq)*unkf_cha(id,iq, ifprj)
unkf_phy(id, 2, ifprj) =  unkf_phy(id, 2, ifprj) + qmat(2, iq)*unkf_cha(id,iq, ifprj)
unkf_phy(id, 3, ifprj) =  unkf_phy(id, 3, ifprj) + qmat(3, iq)*unkf_cha(id,iq, ifprj)
unkf_phy(id, 4, ifprj) =  unkf_phy(id, 4, ifprj) + qmat(4, iq)*unkf_cha(id,iq, ifprj)
enddo
enddo

endif !if(jelem .le. ncell) then

!...End of loop over the faces
enddo

!... Compute the average weight for every face-based polynomial
weige=1.d0/ifprj

!... finally, we can reconstruct the physical polynomails for one cell
unk_phy = c00
do ifa = 1, ifprj
do iq = 1, nq
unk_phy(2, iq) =unk_phy(2, iq) + weige(iq, ifa)*unkf_phy(2,iq,ifa)
unk_phy(3, iq) =unk_phy(3, iq) + weige(iq, ifa)*unkf_phy(3,iq,ifa)
enddo
enddo

!++++++++++Print for debugging++++++++++
!if(ielem.eq.1153.or.ielem.eq.1154.or.ielem.eq.1176)then
!print*,'Final',ielem,unk_phy(2:3,1)
!endif

!...Local system
call getunkno_local(1, 1, mapmt(:,:,ielem), unk_phy)


!...Transfer back the reference configuration

detma = rhsmt(1, ielem)*rhsmt(4, ielem) - rhsmt(2, ielem)*rhsmt(3, ielem)
mtinv(1) = rhsmt(4, ielem)
mtinv(2) =-rhsmt(2, ielem)
mtinv(3) =-rhsmt(3, ielem)
mtinv(4) = rhsmt(1, ielem)
!
mtinv = mtinv/detma
!
rhs(1:nq, 1) = unk_phy(2, 1:nq)*geoph(3, ielem) + unk_phy(3, 1:nq)*geoph(4, ielem)
rhs(1:nq, 2) = unk_phy(2, 1:nq)*geoph(4, ielem) + unk_phy(3, 1:nq)*geoph(5, ielem)
!
unkno(2, 1:nq, ielem) = mtinv(1)*rhs(1:nq, 1) + mtinv(2)*rhs(1:nq, 2)
unkno(3, 1:nq, ielem) = mtinv(3)*rhs(1:nq, 1) + mtinv(4)*rhs(1:nq, 2)


!...  end of the loop over the quad cells
enddo
!
return
end subroutine wenop1_rieminvrnt_quad_shu6

!
!...Rigt vector matrix
!
subroutine getunkno_local(idtrf, isten, mapmt, unkps)
use constant
implicit none
!...Input arrays
integer,   intent(in) ::idtrf, isten
real*8,dimension(1:ndegr, 1:nq, 1:isten),intent(inout) ::unkps
real*8,dimension(1:2,1:2),intent(in)::mapmt

!...local array
real*8:: a11, a12, a21, a22
real*8:: ai11, ai12, ai21, ai22
integer:: is
real*8::rhoxi,rhoet,etxi,etet
real*8::ulxi,ulet,vlxi,vlet
real*8::ucel, vcel

!...Local system
if(idtrf.eq.0)then
a11 = mapmt(2, 2)
a12 =-mapmt(1, 2)
a21 =-mapmt(2, 1)
a22 = mapmt(1, 1)
!
ai11 = mapmt(1, 1)
ai12 = mapmt(1, 2)
ai21 = mapmt(2, 1)
ai22 = mapmt(2, 2)
elseif(idtrf.eq.1)then
ai11 = mapmt(2, 2)
ai12 =-mapmt(1, 2)
ai21 =-mapmt(2, 1)
ai22 = mapmt(1, 1)
!
a11 = mapmt(1, 1)
a12 = mapmt(1, 2)
a21 = mapmt(2, 1)
a22 = mapmt(2, 2)
endif

!...Local stencle
do is = 1, isten
!
rhoxi = unkps(2, 1, is)*a11 + unkps(3, 1, is)*a21
rhoet = unkps(2, 1, is)*a12 + unkps(3, 1, is)*a22
!
ulxi = (ai11*unkps(2, 2, is) + ai12*unkps(2, 3, is))*a11 + &
       (ai11*unkps(3, 2, is) + ai12*unkps(3, 3, is))*a21
ulet = (ai11*unkps(2, 2, is) + ai12*unkps(2, 3, is))*a12 + &
       (ai11*unkps(3, 2, is) + ai12*unkps(3, 3, is))*a22
!
vlxi = (ai21*unkps(2, 2, is) + ai22*unkps(2, 3, is))*a11 + &
       (ai21*unkps(3, 2, is) + ai22*unkps(3, 3, is))*a21
vlet = (ai21*unkps(2, 2, is) + ai22*unkps(2, 3, is))*a12 + &
       (ai21*unkps(3, 2, is) + ai22*unkps(3, 3, is))*a22
!
etxi = unkps(2, 4, is)*a11 + unkps(3, 4, is)*a21
etet = unkps(2, 4, is)*a12 + unkps(3, 4, is)*a22
!
unkps(1, 1, is) =unkps(1, 1, is)
unkps(2, 1, is) =rhoxi
unkps(3, 1, is) =rhoet

ucel = unkps(1, 2, is)*ai11 + unkps(1, 3, is)*ai12
vcel = unkps(1, 2, is)*ai21 + unkps(1, 3, is)*ai22

unkps(1, 2, is) =ucel
unkps(2, 2, is) =ulxi
unkps(3, 2, is) =ulet

unkps(1, 3, is) =vcel
unkps(2, 3, is) =vlxi
unkps(3, 3, is) =vlet

unkps(1, 4, is) =unkps(1, 4, is)
unkps(2, 4, is) =etxi
unkps(3, 4, is) =etet
enddo

end subroutine getunkno_local

!
!...Vector projection
!
subroutine getvector_local(mapmt, dnx,dny)
use constant
implicit none
!...Input arrays
real*8,intent(inout) ::dnx,dny
real*8,dimension(1:2,1:2),intent(in)::mapmt

!...local array
real*8:: a11, a12, a21, a22
real*8:: dnxl, dnyl
integer:: is
real*8::rhoxi,rhoet,etxi,etet
real*8::ulxi,ulet,vlxi,vlet
real*8:: gam1

!...Local system
a11 = mapmt(1, 1)
a12 = mapmt(1, 2)
a21 = mapmt(2, 1)
a22 = mapmt(2, 2)

!...Local vector
dnxl = dnx*a11 + dny*a12
dnyl=  dnx*a21 + dny*a22
!
dnx = dnxl
dny = dnyl
end subroutine getvector_local
