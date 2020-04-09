!
!...Barthkimiter for curved triangle...
!
subroutine barthlimit_lag_curvtria(geoel, coord, coold, ustar, unkno, iptri, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),           intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord, coold
real*8,dimension(1:nq+1, 1:nsize),             intent(out)::aflim
integer*4,dimension(1:nvtri,1:ntria),          intent(in)::iptri
integer, dimension(nifai,nafac),               intent(in)::intfac
real*8, dimension(1:nq+2, 1:npoin),           intent(in):: unmax, unmin
real*8,dimension(1:2, 1:2, 1:nsize),          intent(out)::afvec
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
!
!...Local
!
integer::ipt(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ielem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvtri)::alfa
real*8:: xv(nvtri), yv(nvtri)
real*8:: bt(1:ndegr, 1:nvtri),btv(1:ndegr, 1:nvtri)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
real*8,dimension(1:nq+1,  1:nvtri) ::unknvt
real*8,dimension(1:ndimn, 1:nvtri) :: xpt, xpti
real*8:: dunk(1:nq+1)
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
indbd(intfac(3:(2+nvfac),ifa)) = 1
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
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
!...Get the minimum and maximum at one cell...
!
do ie = 1, ntria
!
ielem = ie
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
do iq=1, nq+2
!
unmax_new(iq, ielem) = maxval(unmax(iq, ipt(1:nvtri)))
unmin_new(iq, ielem) = minval(unmin(iq, ipt(1:nvtri)))
enddo
enddo
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
!
if(ndens.eq.1)then
rhov = 1.d0/unknvt(1, iv)
elseif(ndens.eq.2)then
!
r = xv(iv); s = yv(iv)
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
call getrhoig_tria(rhoi, r, s, xpti)
call getdensity_triallnl(r, s, xpt, xpti, rhoi, rhon)
!
rhov = rhon
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
rhov = unknvt(1, iv)
!
endif
!
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
evtx = unknvt(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!
if(ndens.eq.1)then
unknvt(1, iv) = 1.d0/rhov
elseif(ndens.eq.2)then
unknvt(1, iv) = rhov
elseif(ndens.eq.3)then
unknvt(1, iv) = rhov
endif
unknvt(4 ,iv) = pvtx
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
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
rhoct  = 1.d0/rhom
pctr = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
!
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
alfa(:, iv) = 1.d0
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
elseif(ncase.eq.-7)then
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
aflim(iq, ielem) = minval(alfa(iq, 1:3))
enddo
!
enddo
!
!...Call the symetry-preserving subroutine...
!
call barthlimit_sympre_triac(geoel, coord, ustar, unkno, iptri, intfac, afvec, esuv1, esuv2,unmax_new, unmin_new)
!
!...Degenarate to non-symmetric...
!
do ie = 1, -ntria
!
ielem = ie  
!
afvec(1, 1, ielem) = aflim(2, ielem)
afvec(1, 2, ielem) = 0.d0

afvec(2, 1, ielem) = 0.d0
afvec(2, 2, ielem) = aflim(3, ielem)

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
!
if(ndens.eq.1)then
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif
!
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
!
if(ndens.eq.1)then
rhov = 1.d0/unknvt(1, iv)
elseif(ndens.eq.2)then
!
r = xv(iv); s = yv(iv)
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
call getrhoig_tria(rhoi, r, s, xpti)
call getdensity_triallnl(r, s, xpt, xpti, rhoi, rhon)
!
rhov = rhon
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
rhov = unknvt(1, iv)
!
endif
!
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
evtx = unknvt(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvt(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
elseif(ndens.eq.2)then
unknvt(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
elseif(ndens.eq.3)then
unknvt(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
unknvt(2, iv) = unkno(1,2,ielem)  + dudr*bt(2, iv) + duds*bt(3, iv)
unknvt(3, iv) = unkno(1,3,ielem)  + dvdr*bt(2, iv) + dvds*bt(3, iv)
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
enddo
!!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
do iq = nq+1,nq+1
aflim(iq, ielem) = minval(alfa(iq, 1:3))
enddo
!
!
enddo
!
end subroutine barthlimit_lag_curvtria
!
!...Subroutine for barth limiter based on vertex for primitive variables on trias with symmetry preserving....
!
subroutine barthlimit_lagsym_tria(geoel, coord, coold, ustar, unkno, iptri, bface, intfac, aflim, afvec, unmax, unmin, esuv1, esuv2)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),           intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord, coold
real*8,dimension(1:nq+1, 1:nsize),             intent(out)::aflim
integer*4,dimension(1:nvtri,1:ntria),          intent(in)::iptri
integer*4,dimension(1:nbfai,nbfac),            intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
real*8, dimension(1:nq+2, 1:npoin),           intent(in):: unmax, unmin
real*8,dimension(1:2, 1:2, 1:nsize),          intent(out)::afvec
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
!
!...Local
!
integer::ipt(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ielem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvtri)::alfa
real*8:: xv(nvtri), yv(nvtri)
real*8:: bt(1:ndegr, 1:nvtri),btv(1:ndegr, 1:nvtri)
real*8:: unmax_new(1:nq+2, 1:ncell), unmin_new(1:nq+2, 1:ncell)
real*8,dimension(1:nq+1,  1:nvtri) ::unknvt
real*8,dimension(1:ndimn, 1:nvtri) :: xpt, xpti
real*8:: dunk(1:nq+1)
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
dr = .5d0
ds = .5d0
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
!...Get the minimum and maximum at one cell...
!
do ie = 1, ntria
!
ielem = ie
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
do iq=1, nq+2
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
!
if(ndens.eq.1)then
rhov = 1.d0/unknvt(1, iv)
elseif(ndens.eq.2)then
!
r = xv(iv); s = yv(iv)
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
call getrhoig_tria(rhoi, r, s, xpti)
call getdensity_triallnl(r, s, xpt, xpti, rhoi, rhon)
!
rhov = rhon
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
rhov = unknvt(1, iv)
!
endif
!
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
evtx = unknvt(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
!
if(ndens.eq.1)then
unknvt(1, iv) = 1.d0/rhov
elseif(ndens.eq.2)then
unknvt(1, iv) = rhov
elseif(ndens.eq.3)then
unknvt(1, iv) = rhov
endif
unknvt(4 ,iv) = pvtx
!
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
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
!
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
elseif(ncase.eq.-7)then
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
!...Call the symetry-preserving subroutine(lb:linear boundary)...
!
call barthlimit_sympre_trialb(geoel, coord, ustar, unkno, iptri, bface,intfac, afvec, esuv1, esuv2,unmax_new, unmin_new)
!
!...Degenarate to non-symmetric...
!
do ie = 1, -ntria
!
ielem = ie  
!
afvec(1, 1, ielem) = 0.d0!aflim(2, ielem)
afvec(1, 2, ielem) = 0.d0

afvec(2, 1, ielem) = 0.d0
afvec(2, 2, ielem) = 0.d0!aflim(3, ielem)

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
!
if(ndens.eq.1)then
unctr(1)   = 1.d0/rhoct
elseif(ndens.eq.2)then
unctr(1)   = rhoct
elseif(ndens.eq.3)then
unctr(1)   = rhoct
endif
!
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
!
if(ndens.eq.1)then
rhov = 1.d0/unknvt(1, iv)
elseif(ndens.eq.2)then
!
r = xv(iv); s = yv(iv)
!
xpt(1, 1:nvtri) = coord(1, ipt(1:nvtri))
xpt(2, 1:nvtri) = coord(2, ipt(1:nvtri))
!
xpti(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xpti(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
call getrhoig_tria(rhoi, r, s, xpti)
call getdensity_triallnl(r, s, xpt, xpti, rhoi, rhon)
!
rhov = rhon
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
rhov = unknvt(1, iv)
!
endif
!
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
evtx = unknvt(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rhov*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
if(ndens.eq.1)then
unknvt(1, iv) = rhom + aflim(1, ielem)*(1.d0/rhov - rhom)
elseif(ndens.eq.2)then
unknvt(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
elseif(ndens.eq.3)then
unknvt(1, iv) = 1.d0/(1.d0/rhom + aflim(1, ielem)*(rhov - 1.d0/rhom))
endif
!
dudr = afvec(1, 1, ielem)*unkno(2,2,ielem) +  afvec(1, 2, ielem)*unkno(2,3,ielem)
duds = afvec(1, 1, ielem)*unkno(3,2,ielem) +  afvec(1, 2, ielem)*unkno(3,3,ielem)
dvdr = afvec(2, 1, ielem)*unkno(2,2,ielem) +  afvec(2, 2, ielem)*unkno(2,3,ielem)
dvds = afvec(2, 1, ielem)*unkno(3,2,ielem) +  afvec(2, 2, ielem)*unkno(3,3,ielem)
!
!unknvt(2, iv) = unkno(1,2,ielem)  + dudr*bt(2, iv) + duds*bt(3, iv)
!unknvt(3, iv) = unkno(1,3,ielem)  + dvdr*bt(2, iv) + dvds*bt(3, iv)
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

call barthfct(unmax(iq, ipt(iv)), unmin(iq, ipt(iv)), unctr(iq), dunk(iq), afbar)
alfa(iq, iv) = afbar
enddo
!!
enddo
!
do iq = nq+1,nq+1
aflim(iq, ielem) = minval(alfa(iq, 1:nvtri))
enddo
!
!
enddo
!
end subroutine barthlimit_lagsym_tria
!
!
!...subroutine: Calculate the Riemann input for hybrid tria grids general Riemann solver....
!
subroutine getriem_tria_matrixsym(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
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
real*8,dimension(1:nvtri)::murie
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
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
!
!...ndA=0.5d0*vnorm
!
!vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
!
vnorm(3, 1, 1) = (2.d0*rcoet(1) + rcoet(3))/3.d0*0.5d0*vnorm(3, 1, 1)
vnorm(3, 2, 1) = (2.d0*rcoet(1) + rcoet(2))/3.d0*0.5d0*vnorm(3, 2, 1)
!
vnorm(3, 1, 2) = (2.d0*rcoet(2) + rcoet(1))/3.d0*0.5d0*vnorm(3, 1, 2)
vnorm(3, 2, 2) = (2.d0*rcoet(2) + rcoet(3))/3.d0*0.5d0*vnorm(3, 2, 2)
!
vnorm(3, 1, 3) = (2.d0*rcoet(3) + rcoet(2))/3.d0*0.5d0*vnorm(3, 1, 3)
vnorm(3, 2, 3) = (2.d0*rcoet(3) + rcoet(1))/3.d0*0.5d0*vnorm(3, 2, 3)
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
!if(ielem.eq.5.or.ielem.eq.6) print*,'ielem56',unkno(1, 2:3, ielem),ielem
!
!...Get impedence coefficient...
!
do iv   = 1, nvtri
dux= vlave(1, ipt(iv))-unknvt(2, iv)
duy= vlave(2, ipt(iv))-unknvt(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr + cimpd*rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
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
!...Call Riemann solver...
!
call getriecoef_matrixnew(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
                          unknvt(2:3, iv), sigma(1:2, 1:2, iv),&
                         munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
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
!   if(ipq(iv).eq.738) print*,'p36 muacn(vv) post',ie,murie(iv),munacn(ipq(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:3, iv)
!   if(ipq(iv).eq.738) print*,'p36 muacn(vv) postu',ie,munacu(1:2,ipq(iv)),munacu(1:2,738)/sqrt(munacu(1,738)**2+munacu(2,738)**2)
!   if(ipq(iv).eq.738) print*,'p36 muacn(vv) postmuna',ie,abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
   
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

end subroutine getriem_tria_matrixsym
!
!....domain integral for hybrid linear triangle cells
!
subroutine rhsdomndg_lag_mc_triasym(intfac, iptri, coord,coold, geoel, unkno, rhsel,aflim, afvec )
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
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem) - fluxd(ideg,1:nq)*djaco*rcoef
enddo
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomndg_lag_mc_triasym
!
!...Symmetry preserving techniques...
!
subroutine barthlimit_sympre_tria(geoel, coord, ustar, unkno, iptri, intfac, afvec, esuv1, esuv2,unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8, dimension(1:nq+2, 1:ncell):: unmax_new, unmin_new
!
!...Local
!
integer:: ipt(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp,istor
integer:: ielem,jelem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvtri)::alfa
real*8,  dimension(1:2, 1:nvtri)::xpt
real*8,  dimension(1:nvtri)::dspr, dsps
real*8:: xvt(nvtri), yvt(nvtri)
real*8:: bt(1:ndegr, 1:nvtri)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvtri) ::unknvt
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
dr = 0.5d0
ds = 0.5d0
!
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)
!
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, ntria
!
ielem = ie 
ipt(1:nvtri) = iptri(1:nvtri,ie)
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
dxdr = dxdr + dspr(ishp)*xpt(1,ishp)
dxds = dxds + dsps(ishp)*xpt(1,ishp)

dydr = dydr + dspr(ishp)*xpt(2,ishp)
dyds = dyds + dsps(ishp)*xpt(2,ishp)
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
!..1st invariant of strain tensor...
!
undu(ielem) = dudx + dvdy!
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
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds
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
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
!
!...New mapped velocity components u and v
!
umap = uvtx*mapmt(1, 1, ielem) + vvtx*mapmt(1, 2, ielem)
vmap = uvtx*mapmt(2, 1, ielem) + vvtx*mapmt(2, 2, ielem)
!
unknvt(2, iv) = umap
unknvt(3, iv) = vmap
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
do iv = 1, nvtri
!
ummax(2:3) = unctr(2:3)
ummin(2:3) = unctr(2:3)
!
do  istor=esuv2(ipt(iv))+1,esuv2(ipt(iv)+1)
!
jelem=esuv1(istor)
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
dunk(iq) = unknvt(iq, iv) - unctr(iq)
!
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!
alfa(iq, iv) = afbar
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
enddo
!
do iq = 2, 3
aflim(iq, ielem) = minval(alfa(iq, 1:nvtri))
enddo
!
enddo
!
!...Part 5: Impose limiter using the 1st invariant for mapped barth limiting coef...
!
do ie = 1, ntria
!
ielem = ie  
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
do ie = 1, ntria
!
ielem = ie  
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
afma(ielem) = min(1.d0, .2d0*macel)
!
afdu(ielem) =(afma(ielem))*afdu(ielem) +  (1.d0-afma(ielem))
!if(ielem.eq.1) print*,'bad',afma(ielem),macel
enddo
!
!aflim = 1.d0
!
!...Part 4: Transform back the limiter...
!
do ie = 1,ntria
!
ielem = ie  
!
!aflim(:,ielem) = aflim(:,ielem)*afdu(ielem)
!
enddo!
!
!...Part 6: Transfer back the limiter to the gobal frame...
!
do ie = 1,ntria
!
ielem = ie  
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
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
end subroutine barthlimit_sympre_tria
!
!...Symmetry preserving techniques for curved triangle...
!
subroutine barthlimit_sympre_triac(geoel, coord, ustar, unkno, iptri, intfac, afvec, esuv1, esuv2,unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8, dimension(1:nq+2, 1:ncell):: unmax_new, unmin_new
!
!...Local
!
integer:: ipt(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp,istor
integer:: ielem,jelem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvtri)::alfa
real*8,  dimension(1:2, 1:nvtri)::xpt
real*8,  dimension(1:nvtri)::dspr, dsps
real*8:: xvt(nvtri), yvt(nvtri)
real*8:: bt(1:ndegr, 1:nvtri)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvtri) ::unknvt
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
indbd(intfac(3:(2+nvfac), ifa)) = 1
enddo
!
eps = 1.e-6
mapmt = 0.d0
!
!...shape functions
!
dr = 0.5d0
ds = 0.5d0
!
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
xvt(4) = 0.5d0; yvt(4) = 0.d0
xvt(5) = 0.5d0; yvt(5) = 0.5d0
xvt(6) = 0.d0;  yvt(6) = 0.5d0
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, ntria
!
ielem = ie 
ipt(1:nvtri) = iptri(1:nvtri,ie)
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
!..1st invariant of strain tensor...
!
undu(ielem) = dudx + dvdy!
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
!mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
!mapmt(1, 2, ielem) =          matrc/lmat1
!mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
!mapmt(2, 2, ielem) =          matrc/lmat2
!
elseif(matrc.lt.0.d0)then
!
! print*,'minus matrc',matrc

lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
!mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
!mapmt(2, 2, ielem) =          matrc/lmat1
!mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
!mapmt(1, 2, ielem) =          matrc/lmat2
endif
!
else
!
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
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds
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
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
!
!...New mapped velocity components u and v
!
umap = uvtx*mapmt(1, 1, ielem) + vvtx*mapmt(1, 2, ielem)
vmap = uvtx*mapmt(2, 1, ielem) + vvtx*mapmt(2, 2, ielem)
!
unknvt(2, iv) = umap
unknvt(3, iv) = vmap
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
do iv = 1, nvtri
!
ummax(2:3) = unctr(2:3)
ummin(2:3) = unctr(2:3)
!
do  istor=esuv2(ipt(iv))+1,esuv2(ipt(iv)+1)
!
jelem=esuv1(istor)
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
dunk(iq) = unknvt(iq, iv) - unctr(iq)
!
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!
alfa(iq, iv) = afbar
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
enddo
!
do iq = 2, 3
aflim(iq, ielem) = minval(alfa(iq, 1:3))
enddo
!
enddo
!
!...Part 5: Impose limiter using the 1st invariant for mapped barth limiting coef...
!
do ie = 1, ntria
!
ielem = ie  
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
!geoel(10, ielem) = afdu(ielem)
!
enddo
!
!...Further technique...
!
do ie = 1, ntria
!
ielem = ie  
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
afma(ielem) = min(1.d0, .2d0*macel)
!
afdu(ielem) =(afma(ielem))*afdu(ielem) +  (1.d0-afma(ielem))
!if(ielem.eq.1) print*,'bad',afma(ielem),macel
enddo
!
!aflim = 1.d0
!
!...Part 4: Transform back the limiter...
!
do ie = 1,ntria
!
ielem = ie  
!
!aflim(:,ielem) = aflim(:,ielem)*afdu(ielem)
!
enddo!
!
!...Part 6: Transfer back the limiter to the gobal frame...
!
do ie = 1,ntria
!
ielem = ie  
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
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
end subroutine barthlimit_sympre_triac
!
!...Subroutine to calcualte the maximum and minimum nodal unknows
!for barth limiter based on vertex for primitive variables and 1st invariant....
!
subroutine barthlimit_lag_vtxunknew(unkno, iptri, ipqua, unmax, unmin, coord, geoel)
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
!
data c10   / 1.0d0 /
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
!
!..1st invariant of strain tensor...
!
unctr(nq+2) = dudx + dvdy!
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
!...Get the 1st invariant
!
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
!..1st invariant of strain tensor...
!
unctr(nq+2) = dudx + dvdy!
!
do iv = 1, nvqua
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
end subroutine barthlimit_lag_vtxunknew
!
!...Symmetry preserving techniques for linear triagle with symmetry boudnary(not whole domain)...
!
subroutine barthlimit_sympre_trialb(geoel, coord, ustar, unkno, iptri, bface, intfac, afvec, esuv1, esuv2,unmax_new, unmin_new)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),             intent(inout) ::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
real*8,dimension(1:2, 1:2, 1:nsize),           intent(out)::afvec
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer, dimension(nifai,nafac),               intent(in)::intfac
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8, dimension(1:nq+2, 1:ncell):: unmax_new, unmin_new
!
!...Local
!
integer:: ipt(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa,ishp,istor,ivf
integer:: ielem,jelem
real*8:: unctr(1:nq+1)
real*8,  dimension(1:nq+1, 1:nvtri)::alfa
real*8,  dimension(1:2, 1:nvtri)::xpt
real*8,  dimension(1:nvtri)::dspr, dsps
real*8:: xvt(nvtri), yvt(nvtri)
real*8:: bt(1:ndegr, 1:nvtri)
real*8:: dunk(1:nq+1)
real*8,dimension(1:nq+1,  1:nvtri) ::unknvt
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
real*8:: vemag,ucvem,vcvem,afmag,signlim
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Zero out arrays
!
indbd = 0  !...indbd represents index of boundary node
eps = 1.e-6
mapmt = 0.d0

!...Specify the position of the symmetry boundary condition
if(nrz.eq.1.or.nrz.eq.2)then
  afmag = 0.25d0*pi*(1.d0-0.5d0/(ncell/100.d0))
endif
!
!...Coloring the boundary node
!
!---X axis---
if(ncase.ne.1)then !...Not TGV case...
do ifa =1 ,nbfac
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.221)then
indbd(intfac(3:(2+nvfac), ifa)) = 3 !...Nodes at axis X
endif
endif
enddo
!---Y axis---
do ifa =1 ,nbfac
if(bface(3,ifa).eq.22)then
if(bface(4,ifa).eq.222)then
do ivf = 1, nvfac
if(indbd(intfac(2+ivf, ifa)).eq.3)then
indbd(intfac(2+ivf, ifa)) = 10  !...Intersection node for 2 symmetry boundary faces
else
indbd(intfac(2+ivf, ifa)) = 2   !...Nodes at axis Y
endif
enddo
endif
endif
enddo
else !... TGV case...
do ifa =1 ,nbfac
indbd(intfac(3:(2+nvfac), ifa)) = 1
enddo
endif
!
!...shape functions
!
dr = 0.5d0
ds = 0.5d0
!
xvt(1) = 0.d0; yvt(1) = 0.d0
xvt(2) = 1.d0; yvt(2) = 0.d0
xvt(3) = 0.d0; yvt(3) = 1.d0
!
!...Part 1: Mapping matrix and the 1st invariant of strain tensor...
!
do ie = 1, ntria
!
ielem = ie
ipt(1:nvtri) = iptri(1:nvtri,ie)
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
!
!..1st invariant of strain tensor...
!
undu(ielem) = dudx + dvdy!
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
!
!...Maire deformation gadient
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
!mapmt(1, 1, ielem) = (lamda1-matrd)/lmat1
!mapmt(1, 2, ielem) =          matrc/lmat1
!mapmt(2, 1, ielem) = (lamda2-matrd)/lmat2
!mapmt(2, 2, ielem) =          matrc/lmat2
!
elseif(matrc.lt.0.d0)then
!
! print*,'minus matrc',matrc

lmat1 = sqrt((lamda1-matrd)**2 + matrc**2)
lmat2 = sqrt((lamda2-matrd)**2 + matrc**2)
!
!mapmt(2, 1, ielem) = (lamda1-matrd)/lmat1
!mapmt(2, 2, ielem) =          matrc/lmat1
!mapmt(1, 1, ielem) = (lamda2-matrd)/lmat2
!mapmt(1, 2, ielem) =          matrc/lmat2
endif
!
else
!
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
bt(2, iv) = (xvt(iv)-rc)/dr
bt(3, iv) = (yvt(iv)-sc)/ds
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
uvtx = unknvt(2, iv)
vvtx = unknvt(3, iv)
!
!...New mapped velocity components u and v
!
umap = uvtx*mapmt(1, 1, ielem) + vvtx*mapmt(1, 2, ielem)
vmap = uvtx*mapmt(2, 1, ielem) + vvtx*mapmt(2, 2, ielem)
!
unknvt(2, iv) = umap
unknvt(3, iv) = vmap
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
do iv = 1, nvtri
!
ummax(2:3) = unctr(2:3)
ummin(2:3) = unctr(2:3)
!
do  istor=esuv2(ipt(iv))+1,esuv2(ipt(iv)+1)
!
jelem=esuv1(istor)
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
!...Special treatment of Symmetry BC...
!
if(ncase.ne.1)then
if(indbd(ipt(iv)).eq.3)then
!
!---1st
!
if(nrz.eq.0)then
 ucjel = unkno(1, 2, jelem)
 vcjel =-unkno(1, 3, jelem)
elseif(nrz.eq.1)then
 signlim = sign(1.d0,unkno(1, 2, jelem))
 vemag =signlim*sqrt(unkno(1, 2, jelem)**2 +unkno(1, 3, jelem)**2 )
 ucjel = vemag*cos(afmag)
 vcjel = vemag*sin(afmag)

 ucjel = unkno(1, 2, jelem)
 vcjel =-unkno(1, 3, jelem)

!...Area-weighted only requires theta = 0
elseif(nrz.eq.2)then
ucjel = unkno(1, 2, jelem)
vcjel =-unkno(1, 3, jelem)
endif
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
!
elseif(indbd(ipt(iv)).eq.2)then
!
!...1st
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

elseif(indbd(ipt(iv)).eq.10)then

!...1st x symmetry
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

!...2nd y symmetry
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

!...3rd symmetry
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

!...Special treatment for RZ symmetry boundary
if(nrz.eq.1)then
 if(ielem.eq.1)then
!---Velocity magnitude
signlim = sign(1.d0,unkno(1, 2, jelem))
vemag = signlim*sqrt(unkno(1, 2, jelem)**2 +unkno(1, 3, jelem)**2 )

!...4th
ucvem = vemag*cos(afmag)
vcvem = vemag*sin(afmag)

!...For 1/4 domain 
ucvem = unkno(1, 2, jelem)
vcvem =-unkno(1, 3, jelem)
!
ucjel = ucvem
vcjel = vcvem
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)

!...5th
ucjel = ucvem
vcjel = -vcvem
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)

!...6th
!
ucjel = -ucvem
vcjel =  vcvem
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)

!...7th
ucjel = -ucvem
vcjel = -vcvem
!
uloca = ucjel*mapmt(1, 1, ielem) + vcjel*mapmt(1, 2, ielem)
vloca = ucjel*mapmt(2, 1, ielem) + vcjel*mapmt(2, 2, ielem)
!
ummax(2) = max(ummax(2), uloca)
ummax(3) = max(ummax(3), vloca)

ummin(2) = min(ummin(2), uloca)
ummin(3) = min(ummin(3), vloca)
!
 endif
endif
!
endif
endif
!
enddo
!
do iq = 2, 3 !...only for vector components...
!
dunk(iq) = unknvt(iq, iv) - unctr(iq)
!
call barthfct(ummax(iq), ummin(iq), unctr(iq), dunk(iq), afbar)
!
alfa(iq, iv) = afbar
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
!
enddo
!
do iq = 2, 3
aflim(iq, ielem) = minval(alfa(iq, 1:nvtri))
enddo
!
enddo
!
!...Part 5: Impose limiter using the 1st invariant for mapped barth limiting coef...
!
do ie = 1, ntria
!
ielem = ie
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
!geoel(10, ielem) = afdu(ielem)
!
enddo
!
!...Further technique...
!
do ie = 1, ntria
!
ielem = ie
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
afma(ielem) = min(1.d0, .2d0*macel)
!
afdu(ielem) =(afma(ielem))*afdu(ielem) +  (1.d0-afma(ielem))
!if(ielem.eq.1) print*,'bad',afma(ielem),macel
enddo
!
!aflim = 1.d0
!
!...Part 4: Transform back the limiter...
!
do ie = 1,ntria
!
ielem = ie
!
!aflim(:,ielem) = aflim(:,ielem)*afdu(ielem)
!
enddo!
!
!...Part 6: Transfer back the limiter to the gobal frame...
!
do ie = 1,ntria
!
ielem = ie
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
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
end subroutine barthlimit_sympre_trialb
