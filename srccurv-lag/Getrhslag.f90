!
!...subroutine: Calculate the F^* N dsfor all faces...
!
subroutine getfnds_lag(gflag,gelag,intfac,inpoel,coord,lpnp)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(inout)::gelag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
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
if(ifa.gt.nbfac)ipr(1:nvtri) = inpoel(1:nvtri,ier)
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
 rc = 1.d0/3.d0
 sc = rc
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

end subroutine getfnds_lag
!
!...subroutine: Calculate the nodal velocity U_p^*...
!
subroutine getndvelo_lag(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fstar, itime)
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
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri, 1:nelem),  intent(out)::fstar !...Riemann forces
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
rc = 1.d0/3.d0
sc = rc
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
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
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
 fpres = 0.d0
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
rc = 1.d0/3.d0
sc = rc
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
 fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
 fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
 !
 fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
 fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*ustar(1, ip(iv))- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*ustar(2, ip(iv))- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*ustar(1, ip(iv))- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*ustar(2, ip(iv))- munaul(2, 2, iv, ie)
!
  enddo
!
 enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
end subroutine getndvelo_lag
!
!...Face integral...
!
subroutine rhsifacedg_lag(inpoel,  unkno, ustar, fstar, lpnp, gelag,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri, 1:nelem),intent(in)::lpnp
real*8,dimension(1:ndimn,1:2,1:nvtri,1:nelem),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nelem),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(in)::gelag
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2, 1:nvtri) :: ipf
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::xv(3), yv(3),b(1:3,1:nvtri)
!...local real number
real*8::eps,c00,c05,c10,c20,c13
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
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
!
!print*,'iv',ivt
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
ipf(1, 1) = 3; ipf(2, 1) = 2
ipf(1, 2) = 1; ipf(2, 2) = 3
ipf(1, 3) = 2; ipf(2, 3) = 1

do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
  ip(1:nvtri) = inpoel(1:nvtri, ie)
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
                     ustar(1, ip(iv))*lpnp(1, 1:ndegr, 1, iv, ie) +&
                     ustar(2, ip(iv))*lpnp(2, 1:ndegr, 1, iv, ie) +&
                     ustar(1, ip(iv))*lpnp(1, 1:ndegr, 2, iv, ie) +&
                     ustar(2, ip(iv))*lpnp(2, 1:ndegr, 2, iv, ie) 
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
end subroutine rhsifacedg_lag
!
!...Domain integral...
!
subroutine rhsdomndg_lag(intfac, inpoel, coord, unkno, rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
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
  rc = 1.d0/3.d0
  sc = rc
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
end subroutine rhsdomndg_lag
!
!...Domain integral for source term...
!
subroutine rhsdomnsrcdg_lag(intfac, inpoel, coord, rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
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
rc = 1.d0/3.d0
sc = rc
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
end subroutine rhsdomnsrcdg_lag
!
!...subroutine: Calculate the mass matrix for DGM in Lagrangian framework, npoly=1 for DGP1, and npoly=2 for DGP2...
!
subroutine  getamatr_lag_linear(unkno,amatr,geoel,coord,inpoel)
use constant
implicit none
!...Input
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel
real*8,dimension(1:ndimn,1:npoin),intent(in)::coord
real*8,dimension(1:nmatr,1:nelem),intent(out)::amatr
integer:: inpoel(1:nvtri,1:nelem)
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig
!...Local real array
real*8::xp(1:3), yp(1:3)
!...Local real
real*8::r, s, rc,sc,xg,yg,dr,ds
real*8::rhom, rho0
real*8::wi,djac, volel
real*8::f0,f1,f2,f3,f4
real*8::f5,f6,f7,f8
real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
real*8::b2,b3,b4,b5,b6
!
real*8,allocatable::x5(:,:),b55(:),mmatr(:,:)
real*8,allocatable::weigh(:), posi(:,:)
!
if(npoly==2) allocate(x5(5,5), mmatr(5,5), b55(5))

allocate(weigh(ngausd), posi(2,ngausd))
call rutope(ndimn, ngausd, posi, weigh)
!
!...get amatr...
!...Note: The first term of mass matrix, mass in one cell,
!...is stored in the last term of amatr for convenience...
!
do ie = 1,nelem !...(2)ie = 1,nelem
!
rc = 1.d0/3.d0
sc = rc

dr = 0.5d0
ds = 0.5d0!
volel = geoel(3, ie)

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
r = posi(1,ig)
s = posi(2,ig)
!
wi = weigh(ig)*volel
!
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
!rhom = 0.d0
rhom = unkno(1,1,ie) + unkno(2,1,ie)*b2 + unkno(3,1,ie)*b3
rho0 = 1.d0/rhom
!
f0 = f0 + rho0*wi
f1 = f1 + rho0*(xg-rc)/dr*(xg-rc)/dr*wi
f2 = f2 + rho0*(xg-rc)/dr*(yg-sc)/ds*wi
f3 = f3 + rho0*(yg-sc)/ds*(yg-sc)/ds*wi

if(npoly==2)then
f22 = f22 + rho0*b2*b2*wi
f23 = f23 + rho0*b2*b3*wi
f24 = f24 + rho0*b2*b4*wi
f25 = f25 + rho0*b2*b5*wi
f26 = f26 + rho0*b2*b6*wi

f33 = f33 + rho0*b3*b3*wi
f34 = f34 + rho0*b3*b4*wi
f35 = f35 + rho0*b3*b5*wi
f36 = f36 + rho0*b3*b6*wi

f44 = f44 + rho0*b4*b4*wi
f45 = f45 + rho0*b4*b5*wi
f46 = f46 + rho0*b4*b6*wi

f55 = f55 + rho0*b5*b5*wi
f56 = f56 + rho0*b5*b6*wi

f66 = f66 + rho0*b6*b6*wi
endif
!   if(ie==1) print*,'gauss',xg, yg
enddo

if(npoly==1)then
djac = f1*f3-f2**2

amatr(1,ie) = f3/djac
amatr(2,ie) = -f2/djac
amatr(3,ie) = f1/djac
amatr(4,ie) = 1.d0/f0
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
amatr(1,ie) = x5(1,1)
amatr(2,ie) = x5(1,2)
amatr(3,ie) = x5(1,3)
amatr(4,ie) = x5(1,4)
amatr(5,ie) = x5(1,5)

amatr(6,ie) = x5(2,2)
amatr(7,ie) = x5(2,3)
amatr(8,ie) = x5(2,4)
amatr(9,ie) = x5(2,5)

amatr(10,ie) = x5(3,3)
amatr(11,ie) = x5(3,4)
amatr(12,ie) = x5(3,5)

amatr(13,ie) = x5(4,4)
amatr(14,ie) = x5(4,5)

amatr(15,ie) = x5(5,5)
!
amatr(16,ie) = 1.d0/f0
!
endif
enddo !...(2)ie = 1,nelem

deallocate(weigh, posi)
end subroutine  getamatr_lag_linear
!
!...subroutine: Calculate the mass matrix for DGM in Lagrangian framework, npoly=1 for DGP1, and npoly=2 for DGP2...
!
subroutine  getamatr_lag(unkno,amatr,geoel,coord,inpoel)
use constant
implicit none
!...Input
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel
real*8,dimension(1:ndimn,1:npoin),intent(in)::coord
real*8,dimension(1:nmatr,1:nelem),intent(out)::amatr
integer:: inpoel(1:nvtri,1:nelem)
integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
!...Local integer
integer :: ie, ig, ishp
!...Local real array
real*8::xp(1:2, 1:nptri)
real*8,dimension(1:nptri)::shp, dspr, dsps
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
!
real*8,allocatable::x5(:,:),b55(:),mmatr(:,:)
real*8,allocatable::weigh(:), posi(:,:)
!
data c10 / 1.0d0 /
!
if(npoly==2) allocate(x5(5,5), mmatr(5,5), b55(5))

allocate(weigh(ngausd), posi(2,ngausd))
call rutope(ndimn, ngausd, posi, weigh)
!
!...get amatr...
!...Note: The first term of mass matrix, mass in one cell,
!...is stored in the last term of amatr for convenience...
!
do ie = 1,nelem !...(2)ie = 1,nelem
!
if(ncurv==0)then
xp(1, 1:3) = coord(1, inpoel(1:3,ie))
xp(2, 1:3) = coord(2, inpoel(1:3,ie))
!
xp(1:2,4) = 0.5d0*(xp(1:2,1)+xp(1:2,2))
xp(1:2,5) = 0.5d0*(xp(1:2,2)+xp(1:2,3))
xp(1:2,6) = 0.5d0*(xp(1:2,1)+xp(1:2,3))
elseif(ncurv==1)then
xp(1, 1:nptri) = coord(1,inpoel(1:nptri, ie))
xp(2, 1:nptri) = coord(2,inpoel(1:nptri, ie))
endif
!
rc = 1.d0/3.d0
sc = rc
!
dr = 0.5d0
ds = 0.5d0!
volel = geoel(3, ie)

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
r = posi(1,ig)
s = posi(2,ig)
wi = weigh(ig)!*volel
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
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
!rhom = 0.d0
rhom = unkno(1,1,ie) + unkno(2,1,ie)*b2 + unkno(3,1,ie)*b3
rho0 = 1.d0/rhom
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
!   if(ie==1) print*,'gauss',xg, yg
enddo
!
!print*,'ielem', ie, f0, volel

if(npoly==1)then
det = f1*f3-f2**2

amatr(1,ie) = f3/det
amatr(2,ie) = -f2/det
amatr(3,ie) = f1/det
amatr(4,ie) = 1.d0/f0
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
amatr(1,ie) = x5(1,1)
amatr(2,ie) = x5(1,2)
amatr(3,ie) = x5(1,3)
amatr(4,ie) = x5(1,4)
amatr(5,ie) = x5(1,5)

amatr(6,ie) = x5(2,2)
amatr(7,ie) = x5(2,3)
amatr(8,ie) = x5(2,4)
amatr(9,ie) = x5(2,5)

amatr(10,ie) = x5(3,3)
amatr(11,ie) = x5(3,4)
amatr(12,ie) = x5(3,5)

amatr(13,ie) = x5(4,4)
amatr(14,ie) = x5(4,5)

amatr(15,ie) = x5(5,5)
!
amatr(16,ie) = 1.d0/f0
!
endif
enddo !...(2)ie = 1,nelem

deallocate(weigh, posi)
end subroutine  getamatr_lag
!
!...Subroutine: Solver for DG in Lagrangian frame...
!
subroutine getsolvdg_lag(intfac,inpoel,iptri,ipqua,bface,cooro,coord,coora,geofa,geoel,gesgt0,gesgq0,uchar,&
                         unkno,unold,unkaw,unkgd,unogd,strnq_devtp,&
                         ltelem,estri,esqua,itime,amatr,amagd,ustar,coold, esuv1, esuv2)
use constant
implicit none
!...Input
integer*4,dimension(1:nvtri,1:ntria), intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer*4,dimension(1:nifai,1:nafac), intent(in)::intfac
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coold !...initial coordinates...
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
real*8,dimension(1:ndimn,1:npoin)::cooro
real*8,dimension(1:ngefa,1:nafac)::geofa
real*8,dimension(1:ngeel,1:nsize)::geoel
real*8,dimension(1:3,1:ngesgt,1:ntria+nbfac),intent(in)::gesgt0
real*8,dimension(1:3,1:ngesgq,1:nquad),      intent(in)::gesgq0
real*8,dimension(1:nq)::uchar
real*8,dimension(1:ndegr,1:nq,1:nsize)::unold
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
integer*4,dimension(1:3,1:ncell)::ltelem
integer, dimension(1:nftri,1:ntria), intent(in)::estri
integer, dimension(1:nfqua,1:nquad), intent(in)::esqua
real*8,dimension(nmatr,ncell), intent(in)::amatr
real*8,dimension(1:ndimn,1:npoin)     ::ustar !...nodal velocity
integer*4,dimension(1:nbfai,1:nbfac)::bface
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
integer::itime
real*8,dimension(1:nmatr,1:ncell),     intent(in)::amagd
real*8,dimension(1:ndegr,1:4,1:nsize), intent(in)::unkgd
real*8,dimension(1:ndegr,1:4,1:nsize)::unogd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq, nquad)::strnq_devtp
!...2:Local
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvtri) :: ipt

!...2.1:local real array
real*8,dimension(1:ndimn,1:npoin)::coora
real*8,dimension(1:ndegr,1:nq,1:nsize)::unkaw
real*8,dimension(nmatr,ncell)::amatrd
real*8,dimension(nmati,ncell)::matin
real*8,dimension(1:ndimn, 1:nvtri) :: xpt
real*8,dimension(1:ndimn, 1:nvtri) :: xpti
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndimn, 1:nvqua) :: xpqi
real*8,dimension(1:ndegr,1:nq,1:ncell)::rhsel
real*8,dimension(1:ndegr,1:4,1:ncell) ::rhsgd
real*8,dimension(1:nsize)::volum
real*8,dimension(1:ncell)::dvdt
real*8::m(ndegr, ndegr)
real*8::mgd(ndegr, ndegr)
real*8::gdint(4)
real*8::unint(1:nq)
real*8::alfa
real*8:: rhoi, rhon, rc, sc
real*8:: eps
!...2.2:local integer
integer*4::istag,ie,iq,ig,id,iunk,ideg
integer  ::ielem
integer  ::ipoin
integer  ::im, jm, ijmat
!
eps =1.d-5
!
!...Update the cell center(volume) and cell volume for modal density evolution
!
 if(ncurv.le.1)then
   call getcellvol_lag(iptri, ipqua, geoel, coord)
 elseif(ncurv.gt.1)then
   call getcell_evol_lag(iptri, ipqua, geoel, coord)
 endif
!!call getcellvollinear_lag(iptri, ipqua, geoel, coord)

!if(nrz.eq.2)then
!call getgeoel_lagmc_rzaw(inpoel, iptri, ipqua, geoel, coord, unkno)
!call getamatr_lagmc_rzaw(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
!endif
!
!...Update density with modal(ndens=3) and nodal density(ndens=2) evolution
!
if(ndens.eq.2)then
!...Triangle
 do ie = 1, ntria
   ielem = ie
   unold(1, 1, ielem) = geoel(4,ielem)/geoel(3, ielem)!...density at mass center
 enddo
!...Quad
 do ie = 1, nquad
   ielem = ie + ntria
   unold(1, 1, ielem) = geoel(4,ielem)/geoel(3, ielem) !...density at mass center
 enddo

!...Zero our griaident of density
  unold(2:3, 1, :) = 0.d0

elseif(ndens.eq.3)then
!
 if(ncurv.eq.0)then
   call getamatr_lagdensity(amatrd,geoel,coord, iptri, ipqua)
  !call rhsdomndg_lag_density(ipqua, coord, coold, geoel, unold, amatrd)
   call rhsdensity_lag(iptri, ipqua, coord, coold, geoel, unold, amatrd)
 elseif(ncurv.eq.1)then
   call getamatr_lagdensitycurv(amatrd,matin,geoel,coord, iptri, ipqua)
   call rhsdensity_lagcurv(iptri, ipqua, coord, coold, geoel, unold, matin)
 endif

endif
!
!...Set up time step
!
epsaw = 1.d-1

!...Set cfl No. dynamicly
if(itime.le.2000)then !30300
!cfl=0.002d0
else
!cfl=0.02d0
endif
!

if(nrz.ge.1.and.ncase.eq.7)then !...33860
if(itime.le.2000)then !30300
dtfix = 1d-6
epsaw = 1.d0
ndump=1000
elseif(itime.gt.2000.and.itime.le.6000)then
dtfix = 2.5d-6
epsaw = 4.d-1
ndump=2000
elseif(itime.gt.6000.and.itime.le.10000)then
dtfix = 5.d-6
epsaw = 2.d-1
ndump=2000
else
dtfix = 2.d-5
epsaw = 5.d-2
dtfix = 5.d-5
epsaw = 2.d-2
ndump=160
endif
!
elseif(nrz.ge.1.and.ncase.eq.11)then !11

if(itime.le.7500)then
dtfix = 1.d-4
epsaw = 1.d-3
else
dtfix = 2.5d-6
epsaw = 2.d-2
endif
!
!elseif(nrz.ge.1.and.ncase.eq.4)then !11
elseif(ncase.eq.-4)then !11

if(itime.le.100000)then
!dtfix = 1.d-4
dtfix = 2.d-6
else
dtfix = 2.d-6
endif
!
elseif(ncase.eq.13)then !13
if(itime.le.7500)then
dtfix = 1.d-4
elseif(itime.gt.7500.and.itime.le.17500)then
dtfix = 2.d-5
else
dtfix = 1d-5
endif
!
!if(itime.le.9000)then
!dtfix = 1.d-4
!else
!dtfix = 1.d-5
!endif
!
elseif(ncase.eq.-14)then !14
if(itime.le.1000)then
dtfix = 5.d-4
epsaw = 1.d-2
else
dtfix = 1.d-4
endif
endif
!

!
!...Time loop...
!
do 2500 istag=1, nstag !...(1)istag=1,nstag
!
  alfa=rklag(nstag-1, istag)
  rkstg = istag
!
call getrhsdg_lag_mcnew(ustar,uchar,bface,unold,unkaw,unogd,strnq_devtp,rhsel,rhsgd,intfac,inpoel,iptri,ipqua,&
                     geofa,geoel,gesgt0,gesgq0,&
                     cooro,coold,coora,&
                     esuv1, esuv2,ltelem,estri,esqua,amatr,itime)
!
!
!if(itime.gt.1000)then
!print*,'gerhs112',rhsel(1:3,1,19)
!endif

!...Old subroutine
! call getrhsdg_lag(ustar,uchar,bface,unold,rhsel,intfac,inpoel,geofa,geoel,cooro,ltelem,amatr,itime)

!...Determine the time step size dtfix
!...ndt.eq.1: fixed cfl No.; ndt.eq.0: fixed time step size
if(ndt.eq.1)then

if(istag.eq.1)then
dvdt(:) = rhsel(1, 1, :)
call getdelt_lag(ustar,bface,unold,intfac,iptri,ipqua,geoel,coord,dvdt,itime)
!
if((tacru + dtfix).gt.tend)then
dtfix = tend-tacru
endif
!...Accruded time
tacru = tacru + dtfix

endif

else

if(istag.eq.1)then
!...Accruded time
tacru = tacru + dtfix
endif

endif
!
if(ncase.eq.8.and.tacru.gt..6d0) cfl = 0.005d0



!
if(ncase.eq.4.and.nmatel.eq.1)then
  do ie =1, ncell
  do iq =1, nq
   rhsel(1:ndegr,iq,ie)=dtfix*rhsel(1:ndegr,iq,ie)*tfcus
  enddo
  enddo
!...Rhs for deformation gradient
  do ie =1, ncell
  do iq =1, 4 !...nq variables
   rhsgd(1:ndegr,iq,ie)=dtfix*rhsgd(1:ndegr,iq,ie)*tfcus
  enddo
  enddo

else
  do ie =1, ncell
  do iq =1, nq  !...nq variables
   rhsel(1:ndegr,iq,ie)=dtfix*rhsel(1:ndegr,iq,ie)!*7.265d-3
  enddo
  enddo
!...Rhs for deformation gradient
  do ie =1, ncell
  do iq =1, 4 !...nq variables
   rhsgd(1:ndegr,iq,ie)=dtfix*rhsgd(1:ndegr,iq,ie)!*7.265d-3
  enddo
  enddo
endif

!

!...The mass matrix for high-order dg...
if(npoly.ge.0)then !...(2)npoly.ge.1
!
do ie=1,ncell    !...(3)ie=1,nelem
!
if(npoly.eq.0)then !...(2)npoly.eq.0
!
m(1,1) = amatr(1, ie)
!...Mass matrix for the gradient of deformation
mgd(1,1) = amagd(1, ie)

elseif(npoly==1)then
m(1,1) = amatr(4, ie)
m(1,2) = 0.d0
m(1,3) = 0.d0

m(2,1) = m(1,2)
m(2,2) = amatr(1,ie)
m(2,3) = amatr(2,ie)

m(3,1) = m(1,3)
m(3,2) = amatr(2,ie)
m(3,3) = amatr(3,ie)

!...Mass matrix for the gradient of deformation
mgd(1,1) = amagd(4, ie)
mgd(1,2) = 0.d0
mgd(1,3) = 0.d0

mgd(2,1) = mgd(1,2)
mgd(2,2) = amagd(1,ie)
mgd(2,3) = amagd(2,ie)

mgd(3,1) = mgd(1,3)
mgd(3,2) = amagd(2,ie)
mgd(3,3) = amagd(3,ie)

elseif(npoly==2)then

m = 0.d0
m(1, 1) = amatr(16 ,ie)
!
ijmat = 0
do im = 2, ndegr
do jm = im, ndegr
ijmat = ijmat +1
m(im, jm) = amatr(ijmat, ie)
m(jm, im) = amatr(ijmat, ie)
enddo
enddo

mgd = 0.d0
mgd(1, 1) = amagd(16 ,ie)
!
ijmat = 0
do im = 2, ndegr
do jm = im, ndegr
ijmat = ijmat +1
mgd(im, jm) = amagd(ijmat, ie)
mgd(jm, im) = amagd(ijmat, ie)
enddo
enddo

endif
!
!...solve the mdegr independant varaible...
!
do id =1,ndegr !...(4)id =1,ndegr
!...step 1
unint = 0.d0
do iunk = 1,ndegr
unint(1:nq) = unint(1:nq) + m(id, iunk)*rhsel(iunk,1:nq,ie)
enddo
!
!if(ie.eq.1)print*,'rhs',ie,unint(1),rhsel(1:3,1,ie),m(1, 1:3)
!
!...step 2
!if(nrz.eq.2)then
!unkaw(id,1:nq,ie)= alfa*unkno(id,1:nq,ie) + (1.d0-alfa)*(unold(id,1:nq,ie) + unint(1:nq)*(1.d0+eps))
!endif
!
!if(ie.eq.1832.and.id.eq.2)then
!print*,'ielem1832before',unold(id,nq,ie)
!endif
!
unold(id,1:nq,ie)= alfa*unkno(id,1:nq,ie) + (1.d0-alfa)*(unold(id,1:nq,ie) + unint(1:nq))

!
!...Deformation gradient
!...step 1
gdint = 0.d0
do iunk = 1,ndegr
gdint(1:4) = gdint(1:4) + mgd(id, iunk)*rhsgd(iunk,1:4,ie)
enddo

!...Output foe debugging
!if(ie.eq.1.and.id.eq.5)then
!print*,'ielem1832',id,unogd(id,3,ie),alfa*unkgd(id,3,ie),(1.d0-alfa)*unogd(id,3,ie),(1.d0-alfa)*gdint(3),&
!(1.d0-alfa),mgd(1:6,6),m(1:6,6)
!endif

!...step 2
unogd(id,1:4,ie)= alfa*unkgd(id,1:4,ie) + (1.d0-alfa)*(unogd(id,1:4,ie) + gdint(1:4))

!
enddo !...(4)id =1,ndegr
enddo !...(3)ie=1,nelem
endif !...(2)npoly.ge.1
!
!unold(5:6,:,:) = 0.d0
!unold(3,:,:) = 0.d0
!unold(:,3,:) = 0.d0

!...Updating the physical coordinates based on dx/dt=u
if(ncase.eq.4.and.nmatel.eq.1)then
!
!if(nrz.eq.2)then
!do ipoin = 1, npoin
!coora(1:2, ipoin) = alfa*coord(1:2, ipoin) + (1.d0-alfa)*(cooro(1:2, ipoin) +&
!dtfix*6.726842539779d-3*(1.d0+eps)*ustar(1:2, ipoin))
!enddo
!endif
!
do ipoin = 1, npoin
cooro(1:2, ipoin) = alfa*coord(1:2, ipoin) + (1.d0-alfa)*(cooro(1:2, ipoin) + dtfix*ustar(1:2, ipoin)*tfcus)
enddo

else
!
! if(nrz.eq.2)then
!  do ipoin = 1, npoin
!   coora(1:2, ipoin) = alfa*coord(1:2, ipoin) + (1.d0-alfa)*(cooro(1:2, ipoin) + dtfix*(1.d0+eps)*ustar(1:2, ipoin))
!  enddo
! endif
!
do ipoin = 1, npoin
cooro(1:2, ipoin) = alfa*coord(1:2, ipoin) + (1.d0-alfa)*(cooro(1:2, ipoin) + dtfix*ustar(1:2, ipoin))
enddo

endif


!...Special treatment for the 9th node on quadratic quads...
if(ncurv.eq.1)then
do ie = 1, nquad
xpq(1, 1:8) = cooro(1,ipqua(1:8, ie))
xpq(2, 1:8) = cooro(2,ipqua(1:8, ie))
!
cooro(1, ipqua(9, ie)) = -0.25d0*(xpq(1, 1) + xpq(1, 2) + xpq(1, 3) + xpq(1, 4)) +&
                           0.5d0*(xpq(1, 5) + xpq(1, 6) + xpq(1, 7) + xpq(1, 8))
cooro(2, ipqua(9, ie)) = -0.25d0*(xpq(2, 1) + xpq(2, 2) + xpq(2, 3) + xpq(2, 4)) +&
                           0.5d0*(xpq(2, 5) + xpq(2, 6) + xpq(2, 7) + xpq(2, 8))
enddo

!...Impose some bound for the maximum variation of the curved node (Juan and Shu) ...
!call getmpt_modify(intfac, cooro)

endif

!...Update the cell center(volume) and cell volume for modal density evolution
 if(ncurv.le.1)then
  call getcellvol_lag(iptri, ipqua, geoel, cooro)
 elseif(ncurv.gt.1)then
  call getcell_evol_lag(iptri, ipqua, geoel, cooro)
 endif
!
!if(nrz.eq.2)then
!call getgeoel_lagmc_rzaw(inpoel, iptri, ipqua, geoel, cooro, unold)
!call getamatr_lagmc_rzaw(unold,amatr,geoel,cooro,inpoel, iptri, ipqua)
!endif
!
!  call getcellvollinear_lag(iptri, ipqua, geoel, coord)
!
!...Update the density using LLNL
!
if(ndens.eq.2)then
!
do ie = 1, ntria
   ielem = ie
   unold(1, 1, ielem) = geoel(4,ielem)/geoel(3, ielem) !...density at mass center
enddo
!
do ie = 1, nquad
   ielem = ie + ntria
   unold(1, 1, ielem) = geoel(4,ielem)/geoel(3, ielem) !...density at mass center
enddo
!...Zero our gradident of density
! unold(2:3, 1, :) = 0.d0

elseif(ndens.eq.3)then
 if(ncurv.eq.0)then
  call getamatr_lagdensity(amatrd,geoel,cooro, iptri, ipqua)
! call rhsdomndg_lag_density(ipqua, cooro, coold, geoel, unold, amatrd)
  call rhsdensity_lag(iptri, ipqua, cooro, coold, geoel, unold, amatrd)
 elseif(ncurv.eq.1)then
  call getamatr_lagdensitycurv(amatrd,matin,geoel,cooro, iptri, ipqua)
  call rhsdensity_lagcurv(iptri, ipqua, cooro, coold, geoel, unold, matin)
 endif
endif

!...Output for debugging
!print*,'unkno',unold(1,2,:)
! unold(2:3, 1, :) = 0.d0
!
! open(12,file='solution.dat')
!   do ie=1,ncell
!      if(ie.ge.990.and.ie.le.1000)print*,'after',ie,unold(1:ndegr,1,ie)
!   enddo
! close(12)
!print*,'face', intfac(1:4,94),itime
! stop
!
2500 enddo  !...(1)istag=1,nstag
!
!...Output the solution field after calculation...
!
end subroutine getsolvdg_lag
!
!...subroutine: Calculate the F^* N dsfor all faces...
!
subroutine getfnds_lagnew(gflag,gelag,intfac,inpoel,coord,lpnp)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:ngelg,1:nvtri,1:nelem+nbfac), intent(inout)::gelag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri, 1:nelem),intent(out)::lpnp
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv
!...local integer array
integer,dimension(1:nvtri) :: ip
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
do 50 ie=1, nelem !...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!

!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...coordinates
!
coorp(1, 1:nvtri) = coord(1, ip(1:nvtri))
coorp(2, 1:nvtri) = coord(2, ip(1:nvtri))
!
!...Edge 1
!
anx =   coorp(2, 2) - coorp(2, 1)
any = -(coorp(1, 2) - coorp(1, 1))
!
gelag(1, 1 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 1 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 1 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 2
!
anx =   coorp(2, 3) - coorp(2, 2)
any = -(coorp(1, 3) - coorp(1, 2))
!
gelag(1, 2 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 2 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 2 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 3
!
anx =   coorp(2, 1) - coorp(2, 3)
any = -(coorp(1, 1) - coorp(1, 3))
!
gelag(1, 3 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 3 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 3 , ie) = sqrt(anx**2 + any**2)
!
50 enddo  !...(1)ifa=1,nelem
!  print*,'vnotmfn',gelag(1, 3, 9)
!
! print*,'Inside getfnds_lag'
!
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
!
!print*,'iv',ivt
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...Second part: Get {l_np n_np}
!
do 1100 ie = 1, nelem !...(1)ifa=1,nafac
!
! print*,'second part2',ie
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
!lpnp(1:ndimn, ig, 1, 1, ie) = 0.5d0*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
!lpnp(1:ndimn, ig, 2, 1, ie) = 0.5d0*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
!
!...point 2
!lpnp(1:ndimn, ig, 1, 2, ie) = 0.5d0*gelag(1:ndimn, 1, ie)*gelag(3, 1, ie)
!lpnp(1:ndimn, ig, 2, 2, ie) = 0.5d0*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
!
!...point 3
!lpnp(1:ndimn, ig, 1, 3, ie) = 0.5d0*gelag(1:ndimn, 2, ie)*gelag(3, 2, ie)
!lpnp(1:ndimn, ig, 2, 3, ie) = 0.5d0*gelag(1:ndimn, 3, ie)*gelag(3, 3, ie)
!
enddo
!
!if(ie==1) print*,'lnpn',lpnp(1:2,1,1,1,ie), lpnp(1:2,1,2,1,ie)!
!if(ie==1) print*,'lnpn',lpnp(1:2,1,1,2,ie), lpnp(1:2,1,2,2,ie)!
!if(ie==1) print*,'lnpn',lpnp(1:2,1,1,3,ie), lpnp(1:2,1,2,3,ie)!

1100 enddo

end subroutine getfnds_lagnew
!
!...subroutine getustar_slnoh set up the boundary nodal velocity for shockless Noh test problems...
!
subroutine getustar_slnoh(ustar, bface, intfac, coord)
use constant
implicit none
integer*4,dimension(1:nbfai,1:nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),          intent(in)::intfac
real*8,dimension(1:ndimn,1:npoin),             intent(in)::coord
!
real*8,dimension(1:ndimn,1:npoin),           intent(out)::ustar(1:ndimn,1:npoin)
!
!...Local
!
integer,dimension(1:nvfac)::ipf
integer:: ifa
!
 ustar=0.d0
!00
do ifa = 1, nbfac
!
ipf(1:nvfac) = intfac(3:2+nvfac, ifa) !...Specifying the boundary nodes...
!
ustar(1:2, ipf(1:nvfac)) = -coord(1:2, ipf(1:nvfac))
!
enddo
end subroutine getustar_slnoh
!
!...subroutine getustar_slnoh set up the boundary nodal velocity for shockless Noh test problems...
!
subroutine getustar_slnohmc(ustar, bface, intfac, coord)
use constant
implicit none
integer*4,dimension(1:nbfai,1:nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),          intent(in)::intfac
real*8,dimension(1:ndimn,1:npoin),             intent(in)::coord
!
real*8,dimension(1:ndimn,1:npoin),           intent(out)::ustar(1:ndimn,1:npoin)

!...Local
integer,dimension(1:nvfac)::ipf
integer:: ifa, ipoin
real*8:: xp, yp

!...Zero out ustar
ustar=0.d0

if(ncase.eq.2)then !...Shockless noh
do ipoin = 1, npoin
  xp = coord(1, ipoin)
  yp = coord(2, ipoin)

  ustar(1, ipoin)   = -xp
  ustar(2, ipoin)   = -yp
enddo

elseif(ncase.eq.1)then !...TGV
do ipoin = 1, npoin
  xp = coord(1, ipoin)
  yp = coord(2, ipoin)

  ustar(1, ipoin)   = sin(pi*xp)*cos(pi*yp)
  ustar(2, ipoin)   =-cos(pi*xp)*sin(pi*yp)
enddo

elseif(ncase.eq.3)then !...Noh
do ipoin = 1, npoin
   xp = coord(1, ipoin)
   yp = coord(2, ipoin)

   ustar(1, ipoin)   = -xp/sqrt(xp**2 + yp**2)
   ustar(2, ipoin)   = -yp/sqrt(xp**2 + yp**2)
enddo

elseif(ncase.eq.5)then !...Kidder ball
do ipoin = 1, npoin
  ustar(1, ipoin)   = 0.d0
  ustar(2, ipoin)   = 0.d0
enddo

elseif(ncase.eq.6)then !...Sod
do ipoin = 1, npoin
  ustar(1, ipoin)   = 0.d0
  ustar(2, ipoin)   = 0.d0
enddo

elseif(ncase.eq.7)then !...Sedov
do ipoin = 1, npoin
  xp = coord(1, ipoin)
  yp = coord(2, ipoin)

  ustar(1, ipoin)   = 0.d0
  ustar(2, ipoin)   = 0.d0
enddo

elseif(ncase.eq.13)then !...Saltzman
do ipoin = 1, npoin
   xp = coord(1, ipoin)
  if(abs(xp).lt.1.d-6)then
   ustar(1, ipoin)   = 1.d0
   ustar(2, ipoin)   = 0.d0
  endif
enddo

endif
end subroutine getustar_slnohmc
!
!
!
subroutine callbarthlimiter(unkno, inpoel, intfac)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(inout)::unkno
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
integer, intent(in)::intfac(nifai,nafac)
!
!...Local
!
integer:: iq
!
do iq =1 ,nq

call barthlimiter(unkno(1:ndegr,iq,1:nelem+nbfac), inpoel, intfac)
enddo
!
! print*,'aflim', aflim(1:4,1)
end subroutine callbarthlimiter


!
subroutine barthlimiter(unkno, inpoel, intfac)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1,1:nelem+nbfac),intent(inout)::unkno
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
integer, intent(in)::intfac(nifai,nafac)
!
!...Local
!
integer:: ip(nvtri)
integer:: indbd(npoin)
integer:: ie, iv, iest, iq, ideg, ipoin,ifa
real*8:: aflim(1:nelem), alfal(1:nelem)
real*8:: unctr
real*8,  dimension(1:nvtri)::alfa
real*8:: xv(3), yv(3)
real*8:: b(3, 1:nvtri)
real*8:: unmax(1:npoin), unmin(1:npoin), dunk(nvtri)
real*8,dimension(1:nvtri) ::unknv
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy
!
indbd = 0
!
do ifa =1 ,nbfac
!
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
!...Get the maximum at nodes...
!
unmax(:) = -1.d10
unmin(:) =  1.d10
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
unctr = unkno(1, 1, ie)
!
! if(ie==4) print*,'unctr1',ip,unctr
!
do iv = 1, nvtri
!   if(ie==2) print*,'unmax',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
unmax(ip(iv)) = max(unctr, unmax(ip(iv)))
!
!   if(ie==2) print*,'unmaxpost',iv, iq, unctr(iq),unmax(iq, ip(iv))
!
!   if(ie==2) print*,'unmin',iv, iq, unctr(iq),unmin(iq, ip(iv))
!
unmin(ip(iv)) = min(unctr, unmin(ip(iv)))
!
!   if(ie==2) print*,'unminpost',iv, iq, unctr(iq),unmin(iq, ip(iv))
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
do ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...zero out unknv
!
 dunk = 0.d0
!
do iv   = 1,nvtri
do ideg = 2,mdegr
dunk(iv) = dunk(iv) + unkno(ideg,1,ie)*b(ideg, iv)
enddo
enddo
!
! if(ie==1) print*,'unknv', unknv(1:nq, 1)
!
unctr = unkno(1, 1, ie)
!
do iv = 1, nvtri  !...iv = 1, nvtri
!
if(dunk(iv).gt.0.d0)then
!
fiy = (unmax(ip(iv)) - unctr)/dunk(iv)
alfa(iv) = max(min(1.d0, fiy), 0.d0)
!alfa(iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
elseif(dunk(iv).lt.0.d0)then

fiy = (unmin(ip(iv)) - unctr)/dunk(iv)
alfa(iv) = max(min(1.d0, fiy), 0.d0)
!alfa(iv) = min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0))
!
else
!
alfa(iv) = 1.d0
!
endif
!
!  if(ip(iv)==1) print*,'alfa',ie,iq,iv,alfa(iq, iv)
!
!  if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo  !...iv = 1, nvtri
!
!...Get the minimum alfa...
!
aflim(ie) = minval(alfa(1:nvtri))

!
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
! aflim = 0.d0
!
do ie = 1, nelem
!
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(ie)
!
enddo
!
! print*,'aflim', aflim(1:4,1)
end subroutine barthlimiter
!
!...Limiter
!
subroutine barthlimit_lag_mairev3(coord, ustar, unkno, inpoel, intfac)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),             intent(in) ::ustar, coord
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
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
indbd = 0
!
do ifa =1 ,nbfac
!
indbd(intfac(3:4, ifa)) = 1
!
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
!...Get the maximum at nodes...
!
unmax(:, :) = -1.d10
unmin(:, :) =  1.d10
!
!do ipoin = 1, npoin
!
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

!
enddo

!
! print*,'maximum and minimum1',unmax(1:nq, 1), unmin(1:nq, 1)
! print*,'maximum and minimum122',unmax(1:nq, 122), unmin(1:nq, 122)

!
!...Impose limiter
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
!alfa(1:4, iv) = 1.d0
!alfa(1, iv) = 1.d0
!alfa(4, iv) = 1.d0
!alfa(2, iv) = min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0)
!alfa(3, iv) = min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0)
if(ncase.eq.2)then
 alfa(2, iv) = max(min((ustar(1, ip(iv)) - unctr(2))/(dunk(2)),1.d0),0.d0)
 alfa(3, iv) = max(min((ustar(2, ip(iv)) - unctr(3))/(dunk(3)),1.d0),0.d0)
endif
!
if(ncase .eq. 3)then
if(coord(1, ip(iv)).lt.1.d-6)then!or.abs(coord(1, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf',ip(iv),ie,iv,(- unctr(2))/(dunk(2))
    alfa(2, iv) = min(max(0.5d0*(- unctr(2))/(dunk(2)),0.d0),1.d0)
 !   alfa(3, iv) = alfa(2, iv)
endif
if(coord(2, ip(iv)).lt.1.d-6)then!.or.abs(coord(2, ip(iv))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ip(iv)
    alfa(3, iv) = min(max(0.5d0*( - unctr(3))/(dunk(3)),0.d0),1.d0)
!    alfa(2, iv) = alfa(3, iv)
endif
endif
!
endif
!
enddo
!
!if(ie.eq.20) print*,'aflim',alfa(3,1:nvtri)
!
!...Get the minimum alfa...
!
do iq = 1,nq
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
enddo
!
!
!if(indbd(ie).eq.1) aflim(:, ie) = 0.d0
!
enddo
!
!aflim(1, :) = 0.d0
!aflim(2:3, :) = 0.d0
!aflim(4, :) = 0.d0
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
end subroutine barthlimit_lag_mairev3
!
!...subroutine: Calculate the nodal velocity U_p^* using Legendre-Guass quadrature...
!
subroutine getndvelo_lag_gauss_legendre(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fstar)
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
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri, 1:nelem),intent(out)::fstar !...Riemann forces
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)
integer::fagsp(ngausf, 3)
integer::gssid(ngausf)

!...local real array
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:nq,1:ngausf)::unkng
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8::vnorm(1:3, 1:3)
real*8,dimension(1:ngausf)::murie
real*8::sigmg(1:2, 1:2, 1:ngausf)
real*8::aujmp(1:2, 1:ngausf)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:ndegr)::bg
real*8,dimension(1:nvtri):: xv, yv
real*8::weigh(ngausf), posi(1,ngausf)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::ug, vg, pg, eg
real*8:: r1, r2, s2,s1, rg, sg, r
real*8:: shp1, shp2, wi
real*8::uvtx,vvtx,evtx, pvtx
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
!...Gauss quadrature...
!
call rutope(1, ngausf, posi, weigh)
!
!...Make comparison with Lobatto integration...
!
!posi(1, 1) = -1.d0
!posi(1, 2) =  1.d0
!
 if(ngausf ==2)then
!
    fagsp(1, 1) = 1;  fagsp(2, 1) = 2;
    fagsp(1, 2) = 2;  fagsp(2, 2) = 3;
    fagsp(1, 3) = 3;  fagsp(2, 3) = 1;
!
    gssid(1) = 2;  gssid(2) = 1;
!
elseif(ngausf.eq.4)then
!
  fagsp(1, 1) = 1;  fagsp(2, 1) = 1; fagsp(3, 1) = 2;  fagsp(4, 1) = 2;
  fagsp(1, 2) = 2;  fagsp(2, 2) = 2; fagsp(3, 2) = 3;  fagsp(4, 2) = 3;
  fagsp(1, 3) = 3;  fagsp(2, 3) = 3; fagsp(3, 3) = 1;  fagsp(4, 3) = 1;
!
  gssid(1) = 2; gssid(2) = 2; gssid(3) = 1; gssid(4) = 1;
!
 endif
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
 ipf(1:nvfac) = intfac(3:2+nvfac, ifa)
 indnd(ipf(1:nvfac)) = 1
enddo
endif
!
!...Basis function for vertices...
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!

do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo

!
do 200 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...zero out unknv..
!
unknv = 0.d0
!
do iv   = 1,nvtri
  do ideg = 1,mdegr
     unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
  enddo
! if(ip(iv)==36) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
enddo
!
!...Accumulate nodal velocity...
!
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
vnorm(1:3,  1) = gelag(1:3, 1, ie);
vnorm(1:3,  2) = gelag(1:3, 2, ie);
vnorm(1:3,  3) = gelag(1:3, 3, ie);
!
!...ndA=0.5d0*vnorm
!
vnorm(3, :) = 0.5d0*vnorm(3, :)
!
!...Give the normal vector of every face...
!
do ifa =1, 3
!
r1 = xv(fagsp(1, ifa))
s1 = yv(fagsp(1, ifa))
r2 = xv(fagsp(2, ifa))
s2 = yv(fagsp(2, ifa))
!...
!...zero out unkng
unkng = 0.d0
!
do ig   = 1, ngausf !(1) Gauss point loop...
!
r   = posi(1, ig)
wi  = weigh(ig)
!
shp1 = 0.5d0*(1.d0 - r)
shp2 = 0.5d0*(1.d0 + r)
!
rg = r1*shp1 + r2*shp2
sg = s1*shp1 + s2*shp2
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
rho  = 1.d0/unkng(1, ig)
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
pg = max(eps, (gamlg-1.d0)*rho*(eg - 0.5d0*(ug**2 + vg**2)))
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig) = -pg
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pg
!
!
!if(ip(iv).eq.8) print*,'velocity 8',ie, rho, uvtx,vvtx,vlave(1:2, ip(iv))
!
!...Get the a_c (unit vector)
!aujmp(1:2, iv) = usold(1:2, ip(iv)) - unknv(2:3, iv)
aujmp(1:2, ig) = vlave(1:2, ip(fagsp(ig, ifa))) - unkng(2:3, ig)
!
acnx = aujmp(1, ig)
acny = aujmp(2, ig)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig) = 1.e-11!0.d0;
!print*,'point are reset', ip(iv)
else
aujmp(1:2, ig) = aujmp(1:2, ig)/sqrt(acnx**2 + acny**2)
endif
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo !(1) Gauss point loop...
!
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
do ig   = 1, ngausf
dux= vlave(1, ip(fagsp(ig, ifa)))-unkng(2, ig)
duy= vlave(2, ip(fagsp(ig, ifa)))-unkng(3, ig)
deltu = sqrt(dux**2 + duy**2)
murie(ig) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!if(ip(fagsp(ig, ifa)).eq.155) print*,'murie22', sdctr,rhoct,vlave(1:2, ip(fagsp(ig, ifa))),unkng(2:3, ig),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do ig  = 1, ngausf !(3) ig loop...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
munacn(ip(fagsp(ig, ifa))) = munacn(ip(fagsp(ig, ifa))) + murie(ig)*vnorm(3, ifa)*weigh(ig)* &
abs(gelag(1, ifa, ie)*aujmp(1, ig) + gelag(2, ifa, ie)*aujmp(2, ig))
!
!if(ip(fagsp(ig, ifa)).eq.155) print*,'p36 muacn(vv) ',ie,murie(ig),munacn(ip(fagsp(ig, ifa))),&
!           vnorm(3, ifa),vnorm(1:2, ifa),aujmp(1:2, ig)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p19 muacn(28) prep---',murie(iv), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                 vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
munacu(1, ip(fagsp(ig, ifa))) =  munacu(1, ip(fagsp(ig, ifa))) +&
murie(ig)*vnorm(3, ifa)*weigh(ig)*abs(gelag(1, ifa, ie)*aujmp(1, ig) + gelag(2, ifa, ie)*aujmp(2, ig))*unkng(2, ig)
munacu(2, ip(fagsp(ig, ifa))) =  munacu(2, ip(fagsp(ig, ifa))) +&
murie(ig)*vnorm(3, ifa)*weigh(ig)*abs(gelag(1, ifa, ie)*aujmp(1, ig) + gelag(2, ifa, ie)*aujmp(2, ig))*unkng(3, ig)
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
snsigm(1, ip(fagsp(ig, ifa))) = snsigm(1, ip(fagsp(ig, ifa))) + sigmg(1, 1, ig)*vnorm(3, ifa)*weigh(ig)*gelag(1, ifa, ie) !
snsigm(2, ip(fagsp(ig, ifa))) = snsigm(2, ip(fagsp(ig, ifa))) + sigmg(2, 2, ig)*vnorm(3, ifa)*weigh(ig)*gelag(2, ifa, ie)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
enddo !(3) ig loop...
!
!
!
do ig =1, ngausf
!
munacl(gssid(ig), fagsp(ig, ifa), ie) =  murie(ig)*vnorm(3, ifa)*weigh(ig)* &
abs(gelag(1, ifa, ie)*aujmp(1, ig) + gelag(2, ifa, ie)*aujmp(2, ig))
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaul(1, gssid(ig), fagsp(ig, ifa), ie)    =  munacl(gssid(ig), fagsp(ig, ifa), ie)*unkng(2, ig)
munaul(2, gssid(ig), fagsp(ig, ifa), ie)    =  munacl(gssid(ig), fagsp(ig, ifa), ie)*unkng(3, ig)
!
snsigml(1, gssid(ig), fagsp(ig, ifa), ie)= sigmg(1, 1, ig)*vnorm(3, ifa)*weigh(ig)*gelag(1, ifa, ie)

snsigml(2, gssid(ig), fagsp(ig, ifa), ie)= sigmg(2, 2, ig)*vnorm(3, ifa)*weigh(ig)*gelag(2, ifa, ie)
!
enddo!
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
!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,nelem
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
!...Zero out the bnorm
!
!
!
!...Fourth part: Solve the nodal velocity...
!
!...deactivate the pseudo pressure....
!
fpres = 0.d0
!
!...4.1: Update the Riemann forces at every node...
!
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
endif
enddo
!
!....Bd velocity
!
! print*,'ustar--',ustar(1:2, 155),munacu(1:2,155) ,snsigm(1:2, 155), munacn(155)
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
ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
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
!!
do iv = 1, nvtri
!
!...Basis function!
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
end subroutine getndvelo_lag_gauss_legendre
!
!...Face integral using gauss legendre quadrature distribution...
!
subroutine rhsifacedg_lag_gauss_legendre(inpoel,  unkno, ustar, fstar, lpnp, gelag,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri, 1:nelem),intent(in)::lpnp
real*8,dimension(1:ndimn,1:2,1:nvtri,1:nelem),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nelem),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac),    intent(in)::gelag
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa
integer::ip1,ip2
integer:: infag(2, 3)
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2, 1:nvtri) :: ipf
integer::gssid(1:ngausf)
integer::fagsp(ngausf, 3)
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::xv(3), yv(3),b(1:3,1:nvtri)
real*8,dimension(1:3, 1:2, 1:nvtri)::bg
!
real*8::weigh(ngausf), posi(1,ngausf)
!...local real number
real*8::eps,c00,c05,c10,c20,c13
real*8:: shp1, shp2, wi
real*8:: r1, s1, r2, s2, rg, sg, r
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
!
!...Gauss quadrature...
!
call rutope(1, ngausf, posi, weigh)
!
! posi(1, 1) = -1.d0
! posi(1 ,2 )= 1.d0
!
 if(ngausf.eq.2)then
!
   fagsp(1, 1) = 1;  fagsp(2, 1) = 2;
   fagsp(1, 2) = 2;  fagsp(2, 2) = 3;
   fagsp(1, 3) = 3;  fagsp(2, 3) = 1;
!
   gssid(1) = 2; gssid(2) = 1
!
 elseif(ngausf.eq.4)then
!
  fagsp(1, 1) = 1;  fagsp(2, 1) = 1; fagsp(3, 1) = 2;  fagsp(4, 1) = 2;
  fagsp(1, 2) = 2;  fagsp(2, 2) = 2; fagsp(3, 2) = 3;  fagsp(4, 2) = 3;
  fagsp(1, 3) = 3;  fagsp(2, 3) = 3; fagsp(3, 3) = 1;  fagsp(4, 3) = 1;
  gssid(1) = 2; gssid(2) = 2; gssid(3) = 1; gssid(4) = 1;
!
 endif
!
!...Basis function for all vertices...
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
   b(1, iv) = 1.d0
   b(2, iv) = (xv(iv)-rc)/dr
   b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!
ipf(1, 1) = 3; ipf(2, 1) = 2
ipf(1, 2) = 1; ipf(2, 2) = 3
ipf(1, 3) = 2; ipf(2, 3) = 1

do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
do ifa = 1, 3
!
r1 = xv(fagsp(1, ifa))
s1 = yv(fagsp(1, ifa))
r2 = xv(fagsp(2, ifa))
s2 = yv(fagsp(2, ifa))
!
do ig   = 1, ngausf
!
r   = posi(1, ig)
wi  = weigh(ig)
!
shp1 = 0.5d0*(1.d0 - r)
shp2 = 0.5d0*(1.d0 + r)
!
rg = r1*shp1 + r2*shp2
sg = s1*shp1 + s2*shp2
!
bg(1, gssid(ig), fagsp(ig, ifa)) = 1.d0
bg(2, gssid(ig), fagsp(ig, ifa)) = (rg-rc)/dr
bg(3, gssid(ig), fagsp(ig, ifa)) = (sg-sc)/ds
!
!if(ie.eq.1) print*,'bg',ig,ifa,infag(ig, ifa), rg, sg,bg(1:3,3-ig, infag(ig, ifa))
!
enddo
enddo
!
!...The vertex constituting one cell...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
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
ustar(1, ip(iv))*lpnp(1, 1, 1, iv, ie)*bg(1:ndegr, 1, iv) +&
ustar(2, ip(iv))*lpnp(2, 1, 1, iv, ie)*bg(1:ndegr, 1, iv) +&
ustar(1, ip(iv))*lpnp(1, 1, 2, iv, ie)*bg(1:ndegr, 2, iv) +&
ustar(2, ip(iv))*lpnp(2, 1, 2, iv, ie)*bg(1:ndegr, 2, iv)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstar(1, 1, iv, ie)*bg(1:ndegr, 1, iv) +&
fstar(1, 2, iv, ie)*bg(1:ndegr, 2, iv)

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstar(2, 1, iv, ie)*bg(1:ndegr, 1, iv)  +&
fstar(2, 2, iv, ie)*bg(1:ndegr, 2, iv)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ip(iv))*fstar(1, 1, iv, ie)*bg(1:ndegr, 1, iv) +&
ustar(2, ip(iv))*fstar(2, 1, iv, ie)*bg(1:ndegr, 1, iv) +&
ustar(1, ip(iv))*fstar(1, 2, iv, ie)*bg(1:ndegr, 2, iv) +&
ustar(2, ip(iv))*fstar(2, 2, iv, ie)*bg(1:ndegr, 2, iv)
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
end subroutine rhsifacedg_lag_gauss_legendre
!
!...subroutine: Calculate the F^* N dsfor all faces...
!
subroutine getfnds_lagnew_curv(gflag,gelag,intfac,inpoel,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(inout)::gelag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv
!...local integer array
integer,dimension(1:nvtri) :: ip
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
!...Allocatable array...
real*8, allocatable::cordl(:, :)
!
data c16   /0.1666666666666666d0 /
!
allocate(cordl(1:ndimn, 1:npoin))
!
cordl = coord
!
!...Part 1: Get the control point coordinates...
!
do 40 ie=1, nelem !...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...coordinates
!
cordl(1:2, ip(4)) = 0.5d0*(4.d0*coord(1:2, ip(4)) - coord(1:2, ip(1)) - coord(1:2, ip(2))) !
cordl(1:2, ip(5)) = 0.5d0*(4.d0*coord(1:2, ip(5)) - coord(1:2, ip(2)) - coord(1:2, ip(3))) !...Control points...
cordl(1:2, ip(6)) = 0.5d0*(4.d0*coord(1:2, ip(6)) - coord(1:2, ip(3)) - coord(1:2, ip(1))) !
!
40 enddo
!
!
do 50 ie=1, nelem !...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...coordinates
!
coorp(1, 1:nvtri) = cordl(1, ip(1:nvtri))
coorp(2, 1:nvtri) = cordl(2, ip(1:nvtri))
!
!...Edge 1
!
anx =   coorp(2, 2) - coorp(2, 1)
any = -(coorp(1, 2) - coorp(1, 1))
!
gelag(1, 1 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 1 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 1 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 2
!
anx =   coorp(2, 3) - coorp(2, 2)
any = -(coorp(1, 3) - coorp(1, 2))
!
gelag(1, 2 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 2 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 2 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 3
!
anx =   coorp(2, 1) - coorp(2, 3)
any = -(coorp(1, 1) - coorp(1, 3))
!
gelag(1, 3 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 3 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 3 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 4
!
anx =   coorp(2, 4) - coorp(2, 1)
any = -(coorp(1, 4) - coorp(1, 1))
!
gelag(1, 4 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 4 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 4 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 5
!
anx =   coorp(2, 2) - coorp(2, 4)
any = -(coorp(1, 2) - coorp(1, 4))
!
gelag(1, 5 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 5 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 5 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 6
!
anx =   coorp(2, 5) - coorp(2, 2)
any = -(coorp(1, 5) - coorp(1, 2))
!
gelag(1, 6 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 6 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 6 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 7
!
anx =   coorp(2, 3) - coorp(2, 5)
any = -(coorp(1, 3) - coorp(1, 5))
!
gelag(1, 7 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 7 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 7 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 8
!
anx =   coorp(2, 6) - coorp(2, 3)
any = -(coorp(1, 6) - coorp(1, 3))
!
gelag(1, 8 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 8 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 8 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 9
!
anx =   coorp(2, 1) - coorp(2, 6)
any = -(coorp(1, 1) - coorp(1, 6))
!
gelag(1, 9 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 9 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 9 , ie) = sqrt(anx**2 + any**2)
!
50 enddo  !...(1)ifa=1,nelem
!
! print*,'vnotmfn',gelag(1:3, 1:9, 50)
!
! print*,'Inside getfnds_lag'
!
deallocate(cordl)
!
end subroutine getfnds_lagnew_curv
!
!...Face integral for curved...
!
subroutine rhsifacedg_lagmc_curvtria(inpoel,  unkno, ustar, fstar, gelag, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri,1:ntria),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:ncell),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa,ielem
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2, 1:nvtri) :: ipf
integer,dimension(1:4, 1:nvtri) :: ndshp
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8, dimension(1:2, 1:ndegr, 1:4):: lpnp
real*8, dimension(1:ndegr, 1:4)::clpnp
real*8::xv(nvtri), yv(nvtri),b(1:ndegr,1:nvtri)
real*8::vnorm(1:3, 1:4, 1:nvtri)
!...local real number
real*8::eps,c00,c05,c10,c20,c13,c23
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
data c23   / 0.6666666666666666d0 /
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
!rc = 1.d0/3.d0
!sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
!...Nodes for one vertex
!
ndshp(1, 1) = 3; ndshp(2, 1) = 6; ndshp(3, 1) = 2; ndshp(4, 1) = 4
ndshp(1, 2) = 1; ndshp(2, 2) = 4; ndshp(3, 2) = 3; ndshp(4, 2) = 5;
ndshp(1, 3) = 2; ndshp(2, 3) = 5; ndshp(3, 3) = 1; ndshp(4, 3) = 6;
ndshp(1, 4) = 1; ndshp(2, 4) = 2
ndshp(1, 5) = 2; ndshp(2, 5) = 3
ndshp(1, 6) = 3; ndshp(2, 6) = 1
!
do 550 ie = 1,ntria!...(1)ie = 1,nelem
!
ielem = ie
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
!rc = 1.d0/3.d0
!sc = rc
!
!...Basis function...
!
do iv =1 , nvtri !...low order nodes...
!
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...The vertex constituting one cell...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 3, 1) = gelag(1:3, 9, ie); vnorm(1:3, 4, 1) = gelag(1:3, 4, ie) !...For point ip(1)
!
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 5, ie); vnorm(1:3, 4, 2) = gelag(1:3, 6, ie) !...For point ip(2)
!
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 7, ie); vnorm(1:3, 4, 3) = gelag(1:3, 8, ie) !...For point ip(3)
!
vnorm(1:3, 1, 4) = gelag(1:3, 4, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(3)
!
vnorm(1:3, 1, 5) = gelag(1:3, 6, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(3)
!
vnorm(1:3, 1, 6) = gelag(1:3, 8, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(3)
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, 3
!
do ideg = 1, ndegr
lpnp(1:2, ideg, 3) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(2, iv)))*vnorm(3, 3, iv)*vnorm(1:2, 3, iv)
lpnp(1:2, ideg, 1) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(1, iv)))*vnorm(3, 1, iv)*vnorm(1:2, 1, iv)
!
lpnp(1:2, ideg, 4) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(4, iv)))*vnorm(3, 4, iv)*vnorm(1:2, 4, iv)
lpnp(1:2, ideg, 2) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(3, iv)))*vnorm(3, 2, iv)*vnorm(1:2, 2, iv)
!
clpnp(ideg, 3) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(2, iv)))
clpnp(ideg, 1) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(1, iv)))
clpnp(ideg, 4) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(4, iv)))
clpnp(ideg, 2) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(3, iv)))
enddo

do ifa =1, 4
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr) + c13*ustar(1, ip(iv))*lpnp(1, 1:ndegr, ifa) +&
c13*ustar(2, ip(iv))*lpnp(2, 1:ndegr, ifa)
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr) + c13*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa)
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr) + c13*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
c13*ustar(1, ip(iv))*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa) +&
c13*ustar(2, ip(iv))*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
enddo
!
enddo
!
!...High-order nodes...
!
do iv = 4, nvtri
!
do ideg = 1, ndegr
lpnp(1:2, ideg, 1) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(1, iv)))*vnorm(3, 1, iv)*vnorm(1:2, 1, iv)
lpnp(1:2, ideg, 2) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(2, iv)))*vnorm(3, 2, iv)*vnorm(1:2, 2, iv)
!
clpnp(ideg, 1) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(1, iv)))
clpnp(ideg, 2) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(2, iv)))

enddo

do ifa =1, 2
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr) + c23*ustar(1, ip(iv))*lpnp(1, 1:ndegr, ifa) +&
c23*ustar(2, ip(iv))*lpnp(2, 1:ndegr, ifa)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr) + c23*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa)
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr) + c23*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
c23*ustar(1, ip(iv))*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa) +&
c23*ustar(2, ip(iv))*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
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
end subroutine rhsifacedg_lagmc_curvtria
!
!...Face integral (mass center) for hybrid curv quad...
!
subroutine rhsifacedg_lagmc_curvquad(ipqua, unkno, ustar,fstarq, gelagq, geoel,&
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
real*8::pc1
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
data c13   / 0.3333333333333333d0 /
data c16   / 0.1666666666666666d0 /
data c23   / 0.6666666666666666d0 /
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
!print*,'iv',iv
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
pc1 =0.d0
!
!...Distribute to every corner...
!
do iv = 1, 4
!
do ideg = 1, ndegr
lpnp(1:2, ideg, 3) = 0.1d0*(6.d0*bq(ideg, iv) + 4.d0*bq(ideg, ndshp(2, iv)))*vnorm(3, 3, iv)*vnorm(1:2, 3, iv)
lpnp(1:2, ideg, 1) = 0.1d0*(     bq(ideg, iv) -      bq(ideg, ndshp(1, iv)))*vnorm(3, 1, iv)*vnorm(1:2, 1, iv)
!
lpnp(1:2, ideg, 4) = 0.1d0*(6.d0*bq(ideg, iv) + 4.d0*bq(ideg, ndshp(4, iv)))*vnorm(3, 4, iv)*vnorm(1:2, 4, iv)
lpnp(1:2, ideg, 2) = 0.1d0*(     bq(ideg, iv) -      bq(ideg, ndshp(3, iv)))*vnorm(3, 2, iv)*vnorm(1:2, 2, iv)
!
 clpnp(ideg, 3) = 0.1d0*(6.d0*bq(ideg, iv) + 4.d0*bq(ideg, ndshp(2, iv)))
 clpnp(ideg, 1) = 0.1d0*(     bq(ideg, iv) -      bq(ideg, ndshp(1, iv)))
 clpnp(ideg, 4) = 0.1d0*(6.d0*bq(ideg, iv) + 4.d0*bq(ideg, ndshp(4, iv)))
 clpnp(ideg, 2) = 0.1d0*(     bq(ideg, iv) -      bq(ideg, ndshp(3, iv)))
enddo
!
do ifa =1, 4
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr) + c13*ustar(1, ipq(iv))*lpnp(1, 1:ndegr, ifa) +&
 c13*ustar(2, ipq(iv))*lpnp(2, 1:ndegr, ifa)
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr) + c13*fstarq(1, ifa, iv, ie)*clpnp(1:ndegr, ifa)
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr) + c13*fstarq(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
 c13*ustar(1, ipq(iv))*fstarq(1, ifa, iv, ie)*clpnp(1:ndegr, ifa) +&
 c13*ustar(2, ipq(iv))*fstarq(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
!if(ie==1.and.iv.eq.1)then
!print*,iv,ifa, ipq(iv),ielem,c13*clpnp(2, ifa)
!endif

!if(ie.eq.5.and.iv.eq.2)then
!print*,iv,ifa, ipq(iv),ielem, lpnp(1:2, 3, 1), lpnp(1:2, 3, 3),lpnp(1:2, 3, 2) , lpnp(1:2, 3, 4)
!endif
!
enddo
!
enddo
!
do iv = 5, 8
!
do ideg = 1, ndegr
lpnp(1:2, ideg, 1) = 0.2d0*(4.d0*bq(ideg, iv) +   bq(ideg, ndshp(1, iv)))*vnorm(3, 1, iv)*vnorm(1:2, 1, iv)
lpnp(1:2, ideg, 2) = 0.2d0*(4.d0*bq(ideg, iv) +   bq(ideg, ndshp(2, iv)))*vnorm(3, 2, iv)*vnorm(1:2, 2, iv)
!
 clpnp(ideg, 1) = 0.2d0*(4.d0*bq(ideg, iv) +   bq(ideg, ndshp(1, iv)))
 clpnp(ideg, 2) = 0.2d0*(4.d0*bq(ideg, iv) +   bq(ideg, ndshp(2, iv)))
enddo
!
do ifa =1, 2
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr) + c23*ustar(1, ipq(iv))*lpnp(1, 1:ndegr, ifa) +&
 c23*ustar(2, ipq(iv))*lpnp(2, 1:ndegr, ifa)
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr) + c23*fstarq(1, ifa, iv, ie)*clpnp(1:ndegr, ifa)
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr) + c23*fstarq(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
 c23*ustar(1, ipq(iv))*fstarq(1, ifa, iv, ie)*clpnp(1:ndegr, ifa) +&
 c23*ustar(2, ipq(iv))*fstarq(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
!elnpn(1:ndegr)   = elnpn(1:ndegr)+&
! c23*0.5d0*(ustar(1, ipq(ndshp(1,iv)))*fstarq(1, 4, ndshp(1,iv), ie)+ustar(1, ipq(ndshp(2,iv)))*&
! fstarq(1, 3, ndshp(2,iv), ie))*clpnp(1:ndegr, ifa) +&
!c23*0.5d0*(ustar(2, ipq(ndshp(1,iv)))*fstarq(2, 4, ndshp(1,iv), ie)+ustar(2, ipq(ndshp(2,iv)))*&
! fstarq(2, 3, ndshp(2,iv), ie))*clpnp(1:ndegr, ifa)
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
rhsel(1:ndegr, 1, ielem) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ielem) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ielem) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ielem) =  elnpn(1:ndegr)
!
!if(ie==1831)  print*,'rhs iface',ielem, ie,rhsel(1:3, 1:4, ielem) !,fstarq(1,1:2,)!, lpnp(1:2, 1, 1, 1, ie),lpnp(1:2, 1, 2, 1, ie),lpnp(1:2, 1, 1, 2, ie),&
!                               lpnp(1:2, 1, 2, 2, ie),&
!                            lpnp(1:2, 1, 1, 3, ie),lpnp(1:2, 1, 2, 3, ie),ustar(1:2,ip(1)), &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

650 enddo
!
end subroutine rhsifacedg_lagmc_curvquad
!
!...Domain integral for source term for curved triangle...
!
subroutine rhsdomnsrcdg_lag_curvtria(intfac, inpoel, coord, geoel, rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
real*8,dimension(1:ngeel,1:nsize),          intent(in)::geoel
real*8,dimension(1:ndegr,1:nq,1:ncell),  intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem
!
!...local integer array
!
integer,dimension(1:nvtri) :: ip
!...local real array
real*8,dimension(1:ndimn, 1:nptri) :: coorp
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nptri):: shp, dspr, dsps
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
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
ielem = ie
!
!...Points consitituting one element...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
if(ncurv==0)then
coorp(1, 1:3) = coord(1, inpoel(1:3,ie))
coorp(2, 1:3) = coord(2, inpoel(1:3,ie))
!
coorp(1:2,4) = 0.5d0*(coorp(1:2,1) + coorp(1:2,2))
coorp(1:2,5) = 0.5d0*(coorp(1:2,2) + coorp(1:2,3))
coorp(1:2,6) = 0.5d0*(coorp(1:2,1) + coorp(1:2,3))
elseif(ncurv==1)then
coorp(1, 1:nptri) = coord(1,inpoel(1:nptri, ie))
coorp(2, 1:nptri) = coord(2,inpoel(1:nptri, ie))
endif
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
wi     = weigh(ig)
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
xg = 0.d0
yg = 0.d0
!
do ishp = 1, nptri
xg = xg + shp(ishp)*coorp(1,ishp)
yg = yg + shp(ishp)*coorp(2,ishp)
enddo
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

!   src = 0.5d0*pi/(gamlg-1.d0)*(sin(2.d0*pi*yg)*cos(pi*xg)*sin(pi*yg) - sin(2.d0*pi*xg)*sin(pi*xg)*cos(pi*yg))
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
end subroutine rhsdomnsrcdg_lag_curvtria
!
!...Source term integration for curv quad grids...
!
subroutine rhsdomnsrcdg_lag_curvquad(intfac, ipqua, coord, geoel,rhsel)
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

if(npoly.ge.1)then
 b(2) = (r-rc)/dr
 b(3) = (s-sc)/ds

!DGP2
 if(npoly.eq.2)then
 b(4) = 0.5d0*b(2)*b(2) - geoel(19, ielem)
 b(5) = 0.5d0*b(3)*b(3) - geoel(20, ielem)
 b(6) =       b(2)*b(3) - geoel(21, ielem)
 endif
endif

!...source term
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
end subroutine rhsdomnsrcdg_lag_curvquad

!
!...Face integral for curved...
!
subroutine rhsifacedg_lag_curv3(inpoel,  unkno, ustar, fstar, gelag,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri,1:nelem),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nelem),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(in)::gelag
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2, 1:nvtri) :: ipf
integer,dimension(1:4, 1:nvtri) :: ndshp
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8, dimension(1:2, 1:ndegr, 1:4):: lpnp
real*8, dimension(1:ndegr, 1:4)::clpnp
real*8::xv(nvtri), yv(nvtri),b(1:3,1:nvtri)
real*8::vnorm(1:3, 1:4, 1:nvtri)
!...local real number
real*8::eps,c00,c05,c10,c20,c13,c23
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
data c23   / 0.6666666666666666d0 /
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
!...Nodes for one vertex
!
ndshp(1, 1) = 3; ndshp(2, 1) = 6; ndshp(3, 1) = 2; ndshp(4, 1) = 4
ndshp(1, 2) = 1; ndshp(2, 2) = 4; ndshp(3, 2) = 3; ndshp(4, 2) = 5;
ndshp(1, 3) = 2; ndshp(2, 3) = 5; ndshp(3, 3) = 1; ndshp(4, 3) = 6;
ndshp(1, 4) = 1; ndshp(2, 4) = 2
ndshp(1, 5) = 2; ndshp(2, 5) = 3
ndshp(1, 6) = 3; ndshp(2, 6) = 1

!
do iv =1 , nvtri !...low order nodes...
!
!print*,'iv',ivt
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!

do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 3, 1) = gelag(1:3, 9, ie); vnorm(1:3, 4, 1) = gelag(1:3, 4, ie) !...For point ip(1)
!
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 5, ie); vnorm(1:3, 4, 2) = gelag(1:3, 6, ie) !...For point ip(2)
!
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 7, ie); vnorm(1:3, 4, 3) = gelag(1:3, 8, ie) !...For point ip(3)
!
vnorm(1:3, 1, 4) = gelag(1:3, 4, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(3)
!
vnorm(1:3, 1, 5) = gelag(1:3, 6, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(3)
!
vnorm(1:3, 1, 6) = gelag(1:3, 8, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(3)
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, 3
!
do ideg = 1, ndegr
lpnp(1:2, ideg, 3) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(2, iv)))*vnorm(3, 3, iv)*vnorm(1:2, 3, iv)
lpnp(1:2, ideg, 1) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(1, iv)))*vnorm(3, 1, iv)*vnorm(1:2, 1, iv)
!
lpnp(1:2, ideg, 4) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(4, iv)))*vnorm(3, 4, iv)*vnorm(1:2, 4, iv)
lpnp(1:2, ideg, 2) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(3, iv)))*vnorm(3, 2, iv)*vnorm(1:2, 2, iv)
!
clpnp(ideg, 3) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(2, iv)))
clpnp(ideg, 1) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(1, iv)))
clpnp(ideg, 4) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(4, iv)))
clpnp(ideg, 2) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(3, iv)))
enddo

do ifa =1, 4
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr) + c13*ustar(1, ip(iv))*lpnp(1, 1:ndegr, ifa) +&
                                   c13*ustar(2, ip(iv))*lpnp(2, 1:ndegr, ifa)
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr) + c13*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa)
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr) + c13*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
                                   c13*ustar(1, ip(iv))*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa) +&
                                   c13*ustar(2, ip(iv))*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
enddo
!
enddo
!
!...High-order nodes...
!
do iv = 4, nvtri
!
do ideg = 1, ndegr
lpnp(1:2, ideg, 1) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(1, iv)))*vnorm(3, 1, iv)*vnorm(1:2, 1, iv)
lpnp(1:2, ideg, 2) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(2, iv)))*vnorm(3, 2, iv)*vnorm(1:2, 2, iv)
!
clpnp(ideg, 1) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(1, iv)))
clpnp(ideg, 2) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(2, iv)))

enddo

do ifa =1, 2
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr) + c23*ustar(1, ip(iv))*lpnp(1, 1:ndegr, ifa) +&
                                   c23*ustar(2, ip(iv))*lpnp(2, 1:ndegr, ifa)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr) + c23*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa)
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr) + c23*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
                                  c23*ustar(1, ip(iv))*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa) +&
                                  c23*ustar(2, ip(iv))*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
enddo
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
end subroutine rhsifacedg_lag_curv3
!
!...Curved domain integral ...
!
subroutine rhsdomndg_lagmc_curvtria(intfac, inpoel, coord, geoel, unkno, rhsel, aflim, afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
real*8, dimension(1:nq+1, 1:nsize),      intent(in)::aflim
real*8,dimension(1:2, 1:2, 1:nsize),     intent(in)::afvec
real*8,dimension(1:ndegr,1:nq,1:ncell),   intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem
!
!...local integer array
!
integer,dimension(1:nvtri) :: ip
!...local real array
real*8,dimension(1:ndimn, 1:nptri) :: coorp
real*8,dimension(1:ndegr):: b, dbdr, dbds,bv
real*8:: unknod(1:nq)
real*8, dimension(1:nptri):: shp, dspr, dsps
real*8, dimension(1:ndimn, 1:ndegr):: gdshp
real*8, dimension(1:ndegr, 1:nq)::fluxd
real*8, dimension(1: ndimn, 1:ndimn)::jacbf, jacbg
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8:: rcv, scv
real*8::dr,ds,rc,sc
real*8:: dudr, duds, dvdr, dvds
real*8::xg, yg
real*8::uadv,vadv,eadv,rhoma,rhoad
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr
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
!...Loop over trias
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
ielem = ie

!...Points consitituting one element...
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
if(ncurv==0)then
coorp(1, 1:3) = coord(1, inpoel(1:3,ie))
coorp(2, 1:3) = coord(2, inpoel(1:3,ie))
!
coorp(1:2,4) = 0.5d0*(coorp(1:2,1) + coorp(1:2,2))
coorp(1:2,5) = 0.5d0*(coorp(1:2,2) + coorp(1:2,3))
coorp(1:2,6) = 0.5d0*(coorp(1:2,1) + coorp(1:2,3))
elseif(ncurv==1)then
coorp(1, 1:nptri) = coord(1,inpoel(1:nptri, ie))
coorp(2, 1:nptri) = coord(2,inpoel(1:nptri, ie))
endif

!...Geometry parameters for reference cell...
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
wi     = weigh(ig)
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
dxdr = dxdr + dspr(ishp)*coorp(1,ishp)
dxds = dxds + dsps(ishp)*coorp(1,ishp)

dydr = dydr + dspr(ishp)*coorp(2,ishp)
dyds = dyds + dsps(ishp)*coorp(2,ishp)
enddo
!
djaco = 0.5d0*wi

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
do ideg =1,mdegr
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo

!...Primitive variables...
if(ndens.eq.1)then
rhoma = unknod(1)
rhoad  = 1.d0/rhoma
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

!...Limiting
if(nlimi.eq.6.and.geoel(10,ielem).gt.10.d0)then
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
!
if(ndens.eq.1)then
rhoma = rhomc + aflim(1 ,ielem)*(rhoma - rhomc)
rhoad = 1.d0/rhoma
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
!
enddo !...(2)ig = 1,ngausd
!
550 enddo
end subroutine rhsdomndg_lagmc_curvtria

!
!....domain integral for hybrid curv quad cells
!
subroutine rhsdomndg_lagmc_curvquad(intfac, ipqua, coord, coold, geoel, unkno, rhsel,aflim,afvec, vnulq )
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
real*8::pres
real*8::djaco, wi
real*8::rhomc, rhoct, pctr, uctr, vctr, ectr
real*8:: rhoi, rhon
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
call ruqope_lobatto(2, ngausdq, posiq, weighq)
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
!if(ie==1) print*,'rhs iface idegr',ie,ig,gdshp(1:2,4:5),uadv,vadv,djaco
!
enddo !...(2)ig = 1,ngausd
!
650 enddo

end subroutine rhsdomndg_lagmc_curvquad
!
!...Curved domain integral ...
!
subroutine rhsdomndg_lag_curv3(intfac, inpoel, coord, unkno, rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nelem),   intent(inout)::rhsel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv
!
!...local integer array
!
integer,dimension(1:nvtri) :: ip
!...local real array
real*8,dimension(1:ndimn, 1:nptri) :: coorp
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nptri):: shp, dspr, dsps
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
if(ncurv==0)then
coorp(1, 1:3) = coord(1, inpoel(1:3,ie))
coorp(2, 1:3) = coord(2, inpoel(1:3,ie))
!
coorp(1:2,4) = 0.5d0*(coorp(1:2,1) + coorp(1:2,2))
coorp(1:2,5) = 0.5d0*(coorp(1:2,2) + coorp(1:2,3))
coorp(1:2,6) = 0.5d0*(coorp(1:2,1) + coorp(1:2,3))
elseif(ncurv==1)then
coorp(1, 1:nptri) = coord(1,inpoel(1:nptri, ie))
coorp(2, 1:nptri) = coord(2,inpoel(1:nptri, ie))
endif

!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
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
wi     = weigh(ig)
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
end subroutine rhsdomndg_lag_curv3
!
!...Domain integral for source term for curved element...
!
subroutine rhsdomnsrcdg_lag_curv3(intfac, inpoel, coord, rhsel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
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
real*8,dimension(1:ndimn, 1:nptri) :: coorp
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8:: unknod(1:nq)
real*8, dimension(1:nptri):: shp, dspr, dsps
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
if(ncurv==0)then
coorp(1, 1:3) = coord(1, inpoel(1:3,ie))
coorp(2, 1:3) = coord(2, inpoel(1:3,ie))
!
coorp(1:2,4) = 0.5d0*(coorp(1:2,1) + coorp(1:2,2))
coorp(1:2,5) = 0.5d0*(coorp(1:2,2) + coorp(1:2,3))
coorp(1:2,6) = 0.5d0*(coorp(1:2,1) + coorp(1:2,3))
elseif(ncurv==1)then
coorp(1, 1:nptri) = coord(1,inpoel(1:nptri, ie))
coorp(2, 1:nptri) = coord(2,inpoel(1:nptri, ie))
endif
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
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
wi     = weigh(ig)
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
xg = 0.d0
yg = 0.d0
!
do ishp = 1, nptri
xg = xg + shp(ishp)*coorp(1,ishp)
yg = yg + shp(ishp)*coorp(2,ishp)
enddo
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

!   src = 0.5d0*pi/(gamlg-1.d0)*(sin(2.d0*pi*yg)*cos(pi*xg)*sin(pi*yg) - sin(2.d0*pi*xg)*sin(pi*xg)*cos(pi*yg))
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
end subroutine rhsdomnsrcdg_lag_curv3
!
!...subroutine: Calculate the nodal velocity U_p^* using gauss quadrature distribution...
!
subroutine getndvelo_lag_gauss_labatto(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fstar)
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
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:ngausf/2,1:2,1:nvtri, 1:nelem),intent(out)::fstar !...Riemann forces
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2)     :: ipf
integer::indnd(npoin)
integer:: infag(2, 3)
integer::fagsp(ngausf, 3)
integer::gssid(ngausf)
integer::gpord(ngausf)

!...local real array
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:nq,1:ngausf)::unkng
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8::vnorm(1:3, 1:3)
real*8,dimension(1:ngausf)::murie
real*8::sigmg(1:2, 1:2, 1:ngausf)
real*8::aujmp(1:2, 1:ngausf)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:ndegr)::bg
real*8,dimension(1:nvtri):: xv, yv
real*8::weigh(ngausf), posi(1,ngausf)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::ug, vg, pg, eg
real*8:: r1, r2, s2,s1, rg, sg, r
real*8:: shp1, shp2, wi
real*8::uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dwav1,dwav2
real*8::anx, any
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: cnsup(:), munacn(:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: usold(:,:), munacu(:,:), snsigm(:,:)
real*8,allocatable:: munacl(:,:,:,:)
real*8,allocatable:: snsigml(:,:,:,:,:), munaul(:,:,:,:,:)

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
allocate (munacl(1:ngausf/2, 1:2, 1:nvtri, 1:nelem), munaul(1:ndimn, 1:ngausf/2, 1:2, 1:nvtri,  1:nelem),&
snsigml(1:ndimn, 1:ngausf/2, 1:2,  1:nvtri,  1:nelem))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
!
!...Gauss quadrature...
!
call ruqope_lobatto(1, ngausf, posi, weigh)
!
!posi(1, 1) = -1.d0
!posi(1, 2) =  1.d0
!
if(ngausf ==2)then
fagsp(1, 1) = 1;  fagsp(2, 1) = 2;
fagsp(1, 2) = 2;  fagsp(2, 2) = 3;
fagsp(1, 3) = 3;  fagsp(2, 3) = 1;
!
gssid(1) = 2;  gssid(2) = 1;
!
gpord(1) = 1; gpord(2) = 1!
!
elseif(ngausf.eq.4)then
!
fagsp(1, 1) = 1;  fagsp(2, 1) = 1; fagsp(3, 1) = 2;  fagsp(4, 1) = 2;
fagsp(1, 2) = 2;  fagsp(2, 2) = 2; fagsp(3, 2) = 3;  fagsp(4, 2) = 3;
fagsp(1, 3) = 3;  fagsp(2, 3) = 3; fagsp(3, 3) = 1;  fagsp(4, 3) = 1;
!
gssid(1) = 2; gssid(2) = 2; gssid(3) = 1; gssid(4) = 1;
!
gpord(1) = 1; gpord(2) = 2; gpord(3) = 2; gpord(4) = 1;
!
endif
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
do 200 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
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
!

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
vnorm(1:3,  1) = gelag(1:3, 1, ie);
vnorm(1:3,  2) = gelag(1:3, 2, ie);
vnorm(1:3,  3) = gelag(1:3, 3, ie);
!
!...ndA=0.5d0*vnorm
!
vnorm(3, :) = 0.5d0*vnorm(3, :)
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
!
do iv =1 ,nvtri
!...Basis function
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds

enddo
!
infag (1 ,1) = 1; infag (2 ,1) = 2;
infag (1 ,2) = 2; infag (2 ,2) = 3;
infag (1 ,3) = 3; infag (2 ,3) = 1;
!
!...Give the normal vector of every face...
!
do ifa =1, 3
!
r1 = xv(infag(1, ifa))
s1 = yv(infag(1, ifa))
r2 = xv(infag(2, ifa))
s2 = yv(infag(2, ifa))
!...
!...zero out unkng
unkng = 0.d0
!
do ig   = 1, ngausf
!
r   = posi(1, ig)
wi  = weigh(ig)
!
shp1 = 0.5d0*(1.d0 - r)
shp2 = 0.5d0*(1.d0 + r)
!
rg = r1*shp1 + r2*shp2
sg = s1*shp1 + s2*shp2
!
bg(1) = 1.d0
bg(2) = (rg-rc)/dr
bg(3) = (sg-sc)/ds
!
do ideg = 1,mdegr
unkng(1:nq, ig) = unkng(1:nq, ig) + unkno(ideg,1:nq,ie)*bg(ideg)
enddo
!
rho  = 1.d0/unkng(1, ig)
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
pg = max(eps, (gamlg-1.d0)*rho*(eg - 0.5d0*(ug**2 + vg**2)))
!
!pvtx = 0.25d0*(cos(2.d0*pi*coord(1, ip(iv))) + cos(2.d0*pi*coord(2, ip(iv)))) + 1.d0
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig) = -pg
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pg
!
!
!if(ip(iv).eq.8) print*,'velocity 8',ie, rho, uvtx,vvtx,vlave(1:2, ip(iv))
!
!...Get the a_c (unit vector)
!aujmp(1:2, iv) = usold(1:2, ip(iv)) - unknv(2:3, iv)
aujmp(1:2, ig) = vlave(1:2, ip(fagsp(ig, ifa))) - unkng(2:3, ig)
!
acnx = aujmp(1, ig)
acny = aujmp(2, ig)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig) = 1.e-11!0.d0;
!print*,'point are reset', ip(iv)
else
aujmp(1:2, ig) = aujmp(1:2, ig)/sqrt(acnx**2 + acny**2)
endif
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
!
!
enddo !....ig
!
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
do ig   = 1, ngausf
dux= vlave(1, ip(fagsp(ig, ifa)))-unkng(2, ig)
duy= vlave(2, ip(fagsp(ig, ifa)))-unkng(3, ig)
deltu = sqrt(dux**2 + duy**2)
murie(ig) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ip(iv).eq.5) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ip(iv)),unknv(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do ig  = 1, ngausf
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
munacn(ip(fagsp(ig, ifa))) = munacn(ip(fagsp(ig, ifa))) + murie(ig)*vnorm(3, ifa)*weigh(ig)* &
abs(gelag(1, ifa, ie)*aujmp(1, ig) + gelag(2, ifa, ie)*aujmp(2, ig))
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
munacu(1, ip(fagsp(ig, ifa))) =  munacu(1, ip(fagsp(ig, ifa))) +&
murie(ig)*vnorm(3, ifa)*weigh(ig)*abs(gelag(1, ifa, ie)*aujmp(1, ig) + gelag(2, ifa, ie)*aujmp(2, ig))*unkng(2, ig)
munacu(2, ip(fagsp(ig, ifa))) =  munacu(2, ip(fagsp(ig, ifa))) +&
murie(ig)*vnorm(3, ifa)*weigh(ig)*abs(gelag(1, ifa, ie)*aujmp(1, ig) + gelag(2, ifa, ie)*aujmp(2, ig))*unkng(3, ig)
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
snsigm(1, ip(fagsp(ig, ifa))) = snsigm(1, ip(fagsp(ig, ifa))) + sigmg(1, 1, ig)*vnorm(3, ifa)*weigh(ig)*gelag(1, ifa, ie) !
snsigm(2, ip(fagsp(ig, ifa))) = snsigm(2, ip(fagsp(ig, ifa))) + sigmg(2, 2, ig)*vnorm(3, ifa)*weigh(ig)*gelag(2, ifa, ie)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
!
do ig =1, ngausf
!
munacl(gpord(ig), gssid(ig), fagsp(ig, ifa), ie) =  murie(ig)*vnorm(3, ifa)*weigh(ig)* &
abs(gelag(1, ifa, ie)*aujmp(1, ig) + gelag(2, ifa, ie)*aujmp(2, ig))
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaul(1, gpord(ig), gssid(ig), fagsp(ig, ifa), ie)    =  munacl(gpord(ig), gssid(ig), fagsp(ig, ifa), ie)*unkng(2, ig)
munaul(2, gpord(ig), gssid(ig), fagsp(ig, ifa), ie)    =  munacl(gpord(ig), gssid(ig), fagsp(ig, ifa), ie)*unkng(3, ig)
!
snsigml(1, gpord(ig), gssid(ig), fagsp(ig, ifa), ie)= sigmg(1, 1, ig)*vnorm(3, ifa)*weigh(ig)*gelag(1, ifa, ie)

snsigml(2, gpord(ig), gssid(ig), fagsp(ig, ifa), ie)= sigmg(2, 2, ig)*vnorm(3, ifa)*weigh(ig)*gelag(2, ifa, ie)
!
enddo!
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
!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,nelem
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
!...Zero out the bnorm
!
!
!
!...Fourth part: Solve the nodal velocity...
!
!...deactivate the pseudo pressure....
!
fpres = 0.d0
!
!...4.1: Update the Riemann forces at every node...
!
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
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
ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
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
!!
do ifa = 1, 3
!
do ig =1, ngausf
!
!...Basis function!
!
  fstar(1, gpord(ig), gssid(ig), fagsp(ig, ifa), ie) = snsigml(1, gpord(ig), gssid(ig), fagsp(ig, ifa), ie) +&
                           munacl(gpord(ig), gssid(ig), fagsp(ig, ifa), ie)*ustar(1, ip(fagsp(ig, ifa)))-&
                                                       munaul(1, gpord(ig), gssid(ig), fagsp(ig, ifa), ie)
!
  fstar(2, gpord(ig), gssid(ig), fagsp(ig, ifa), ie) = snsigml(2, gpord(ig),gssid(ig), fagsp(ig, ifa), ie) +&
                           munacl(gpord(ig), gssid(ig), fagsp(ig, ifa), ie)*ustar(2, ip(fagsp(ig, ifa)))-&
                                                       munaul(2, gpord(ig),  gssid(ig), fagsp(ig, ifa), ie)
!
enddo
!
enddo
!
enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
end subroutine getndvelo_lag_gauss_labatto
!
!...iface for gauss labatto
!
!
!...Face integral using gauss quadrature distribution...
!
subroutine rhsifacedg_lag_gauss_labatto(inpoel,  unkno, ustar, fstar, lpnp, gelag,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri, 1:nelem),intent(in)::lpnp
real*8,dimension(1:ndimn,1:ngausf/2,1:2,1:nvtri,1:nelem),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nelem),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(in)::gelag
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa
integer::ip1,ip2
integer:: infag(2, 3)
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2, 1:nvtri) :: ipf
integer::gssid(1:ngausf)
integer::fagsp(ngausf, 3)
integer::gpord(ngausf)
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::xv(3), yv(3),b(1:3,1:nvtri)
real*8,dimension(1:3, 1:ngausf/2, 1:2, 1:nvtri)::bg
!
real*8::weigh(ngausf), posi(1,ngausf), weig(ngausf/2)
!...local real number
real*8::eps,c00,c05,c10,c20,c13
real*8:: shp1, shp2, wi
real*8:: r1, s1, r2, s2, rg, sg, r
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
!
!...Gauss quadrature...
!
!
call ruqope_lobatto(1, ngausf, posi, weigh)
!
! posi(1, 1) = -1.d0
! posi(1 ,2 )= 1.d0
!
if(ngausf.eq.2)then
!
fagsp(1, 1) = 1;  fagsp(2, 1) = 2;
fagsp(1, 2) = 2;  fagsp(2, 2) = 3;
fagsp(1, 3) = 3;  fagsp(2, 3) = 1;
!
gssid(1) = 2; gssid(2) = 1
!
gpord(1) = 1; gpord(2) = 1!
!
weig(1:ngausf/2) = weigh(1:ngausf/2)
!
elseif(ngausf.eq.4)then
!
fagsp(1, 1) = 1;  fagsp(2, 1) = 1; fagsp(3, 1) = 2;  fagsp(4, 1) = 2;
fagsp(1, 2) = 2;  fagsp(2, 2) = 2; fagsp(3, 2) = 3;  fagsp(4, 2) = 3;
fagsp(1, 3) = 3;  fagsp(2, 3) = 3; fagsp(3, 3) = 1;  fagsp(4, 3) = 1;
!
gssid(1) = 2; gssid(2) = 2; gssid(3) = 1; gssid(4) = 1;
!
gpord(1) = 1; gpord(2) = 2; gpord(3) = 2; gpord(4) = 1;
!
weig(1:ngausf/2) = weigh(1:ngausf/2)
!
endif
!
infag (1 ,1) = 1; infag (2 ,1) = 2;
infag (1 ,2) = 2; infag (2 ,2) = 3;
infag (1 ,3) = 3; infag (2 ,3) = 1;
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
!
!print*,'iv',ivt
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!

ipf(1, 1) = 3; ipf(2, 1) = 2
ipf(1, 2) = 1; ipf(2, 2) = 3
ipf(1, 3) = 2; ipf(2, 3) = 1

do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
do ifa = 1, 3
!
!
r1 = xv(infag(1, ifa))
s1 = yv(infag(1, ifa))
r2 = xv(infag(2, ifa))
s2 = yv(infag(2, ifa))
!
do ig   = 1, ngausf
!
r   = posi(1, ig)
wi  = weigh(ig)
!
shp1 = 0.5d0*(1.d0 - r)
shp2 = 0.5d0*(1.d0 + r)
!
!
rg = r1*shp1 + r2*shp2
sg = s1*shp1 + s2*shp2
!
bg(1, gpord(ig), gssid(ig), fagsp(ig, ifa)) = 1.d0
bg(2, gpord(ig), gssid(ig), fagsp(ig, ifa)) = (rg-rc)/dr
bg(3, gpord(ig), gssid(ig), fagsp(ig, ifa)) = (sg-sc)/ds
!
!if(ie.eq.1) print*,'bg',ig,ifa,infag(ig, ifa), rg, sg,bg(1:3,3-ig, infag(ig, ifa))
!
enddo
enddo
!
!...The vertex constituting one cell...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
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
do ig = 1, ngausf/2
!
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ip(iv))*lpnp(1, 1, 1, iv, ie)*bg(1:ndegr, ig, 1, iv)*weig(ig) +&
ustar(2, ip(iv))*lpnp(2, 1, 1, iv, ie)*bg(1:ndegr, ig, 1, iv)*weig(ig) +&
ustar(1, ip(iv))*lpnp(1, 1, 2, iv, ie)*bg(1:ndegr, ig, 2, iv)*weig(ig) +&
ustar(2, ip(iv))*lpnp(2, 1, 2, iv, ie)*bg(1:ndegr, ig, 2, iv)*weig(ig)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstar(1, ig, 1, iv, ie)*bg(1:ndegr, ig, 1, iv) +&
fstar(1, ig, 2, iv, ie)*bg(1:ndegr, ig, 2, iv)

!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstar(2, ig, 1, iv, ie)*bg(1:ndegr, ig, 1, iv)  +&
fstar(2, ig, 2, iv, ie)*bg(1:ndegr, ig, 2, iv)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ip(iv))*fstar(1, ig, 1, iv, ie)*bg(1:ndegr, ig, 1, iv) +&
ustar(2, ip(iv))*fstar(2, ig, 1, iv, ie)*bg(1:ndegr, ig, 1, iv) +&
ustar(1, ip(iv))*fstar(1, ig, 2, iv, ie)*bg(1:ndegr, ig, 2, iv) +&
ustar(2, ip(iv))*fstar(2, ig, 2, iv, ie)*bg(1:ndegr, ig, 2, iv)
!
enddo
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
end subroutine rhsifacedg_lag_gauss_labatto
!
!...subroutine: Calculate the nodal velocity U_p^* for curved cell...
!
subroutine getndvelo_lag_curv4(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fstar)
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
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri, 1:nelem),  intent(out)::fstar !...Riemann forces
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)
!
!...local real array
!
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:nvtri)::murie
real*8::vnorm(1:3, 1:4, 1:nvtri)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8::aujmp(1:3, 1:nvtri)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:nvtri):: xv, yv
real*8::usnom(2), ustng(2)
!...local real number
real*8::unorm, tfx, tfy, nfx, nfy,signm, utang
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dwav1,dwav2
real*8::anx, any
!...Allocatable arrays
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: cnsup(:), munacn(:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: usold(:,:), munacu(:,:), snsigm(:,:)
real*8,allocatable:: munacl(:,:,:)
real*8,allocatable:: snsigml(:,:,:,:), munaul(:,:,:,:)
!
!...For high-order nodes...
!
real*8,allocatable::muutg(:, :), mutag(:)

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
allocate (munacl(1:4, 1:nvtri, 1:nelem), munaul(1:ndimn, 1:4, 1:nvtri,  1:nelem),&
snsigml(1:ndimn, 1:4,  1:nvtri,  1:nelem))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
allocate (muutg(2, npoin), mutag(npoin))
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
ipf(1:3) = intfac(3:5, ifa)
indnd(ipf(1:3)) = 1
enddo
endif
!
!...Part I: Calculate the averaged velocity...
!
do 200 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...shape functions
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0;  yv(1) = 0.d0
xv(2) = 1.d0;  yv(2) = 0.d0
xv(3) = 0.d0;  yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
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
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
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
vnorm(1:3, 3, 1) = gelag(1:3, 9, ie); vnorm(1:3, 4, 1) = gelag(1:3, 4, ie) !...For point ip(1)
!
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 5, ie); vnorm(1:3, 4, 2) = gelag(1:3, 6, ie) !...For point ip(2)
!
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 7, ie); vnorm(1:3, 4, 3) = gelag(1:3, 8, ie) !...For point ip(3)
!
vnorm(1:3, 1, 4) = gelag(1:3, 4, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(3)
!
vnorm(1:3, 1, 5) = gelag(1:3, 6, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(3)
!
vnorm(1:3, 1, 6) = gelag(1:3, 8, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(3)

!
!...ndA=0.5d0*vnorm
!
!vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
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
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
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
!aujmp(3, iv) = sqrt(acnx**2 + acny**2)
!print*,'point are reset', ip(iv)
else
aujmp(1:2, iv) = aujmp(1:2, iv)/sqrt(acnx**2 + acny**2)
!aujmp(3, iv) = sqrt(acnx**2 + acny**2)
endif
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
murie(iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ip(iv).eq.5) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ip(iv)),unknv(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, 3 !...Low-order nodes...
do ifa = 1, 4 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!  if(abs(aujmp(3, iv)).lt.1.d-2) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
!
  if(ifa.le.2)then
!
  else
!
munacn(ip(iv)) = munacn(ip(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
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
munacu(1, ip(iv)) =  munacu(1, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(2, iv)
munacu(2, ip(iv)) =  munacu(2, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(3, iv)
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
endif
!
!...Local variable...
!
munacl(ifa, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
!
munaul(1:2, ifa, iv, ie)    =  munacl(ifa, iv, ie)*unknv(2:3, iv)
!
snsigml(1, ifa, iv, ie)= sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv)

snsigml(2, ifa, iv, ie)= sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
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
!     snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv)

!     snsigml(2, 1, iv, ie)= sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
!     snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv)

!     snsigml(2, 2, iv, ie)= sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)!
enddo
!
!...Get the summed denominator cooefficients sum(mu*n*a_c) for high-order nodes...
!
do iv  = 4, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!!
! if(ip(iv).eq.281) print*,'p19 muacn(281) pre++', munacn(281),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
munacn(ip(iv)) = munacn(ip(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
!   if(ip(iv).eq.281) print*,'p281 muacn(vv) post',ie,murie(iv),munacn(ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.281) print*,'p19 muacn(28) prep---',murie(iv), munacn(281), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
munacu(1, ip(iv)) =  munacu(1, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(2, iv)
munacu(2, ip(iv)) =  munacu(2, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(3, iv)
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
munacl(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)* &
abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))
!
munacl(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)* &
abs(vnorm(1, 2, iv)*aujmp(1, iv) + vnorm(2, 2, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
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
!
!...Zero out the bnorm
!
!
!...
!...Zero out boundary pressure...
!
!
!...Fourth part: Solve the nodal velocity...
!
!...deactivate the pseudo pressure....
!
fpres = 0.d0
!
!...4.1: Update the Riemann forces at every node...
!
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
!if(idcrd(ipoin).eq.0)then
ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
!elseif(idcrd(ipoin).eq.1)then
!  if(muacn(ipoin).gt.1.d-9)then
!   ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
!   ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
!  endif
!endif
endif
enddo
!
!...Deal with the high-order nodes...
!
!
muutg = 0.d0
mutag = 0.d0
!
do 350 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
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
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
do iv =1 ,nvtri
!...Basis function
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds

enddo
!
!...Impedence...
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
murie(iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!if(ip(iv).eq.53) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ip(iv)),unknv(2, iv),unkno(1,2,ie),ie
enddo
!
!...
!...zero out unknv
unknv = 0.d0
!
do iv   = 4,nvtri
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
enddo
!
!
do iv  = 4, nvtri
muutg(1:2, ip(iv)) = muutg(1:2, ip(iv)) + murie(iv)*unknv(2:3, iv)
mutag(ip(iv)) = mutag(ip(iv)) + murie(iv)
enddo
!
350 enddo
!
!
!
do 450 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 4) = gelag(1:3, 1, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(3)
!
vnorm(1:3, 1, 5) = gelag(1:3, 2, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(3)
!
vnorm(1:3, 1, 6) = gelag(1:3, 3, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(3)
!
do iv = 4, nvtri
!
!if(ip(iv).eq.212.or.ip(iv).eq.216.or.ip(iv).eq.51) print*,'bad point',ie,ip(iv),munacn(ip(iv))
!
if(munacn(ip(iv)).lt.1.d-6)then
!
!
nfx = vnorm(1, 1, iv); nfy = vnorm(2, 1, iv)
!
tfx = nfy ; tfy = -nfx
!
unorm = munacu(1, ip(iv))*nfx + munacu(2, ip(iv))*nfy
signm = snsigm(1, ip(iv))*nfx + snsigm(2, ip(iv))*nfy
!
usnom(1) = (unorm*nfx - signm*nfx)/(munacn(ip(iv)) + 1.d-9)
usnom(2) = (unorm*nfy - signm*nfy)/(munacn(ip(iv)) + 1.d-9)
!
usnom = 0.d0
!
utang = muutg(1, ip(iv))*tfx + muutg(2, ip(iv))*tfy
!
ustng(1) = utang*tfx/mutag(ip(iv))
ustng(2) = utang*tfy/mutag(ip(iv))
!
ustar(1 ,ip(iv)) = usnom(1) + ustng(1)
ustar(2, ip(iv)) = usnom(2) + ustng(2)
!
endif
!
!if(ip(iv).eq.281) print*, 'ustar', usnom, ustng
enddo
!
450 enddo

!
!....Bd velocity
!
! print*,'ustar--',ustar(1:2, 281),munacu(1:2,281) ,snsigm(1:2,281), munacn(281)
!
!
if(ncase.eq.1)then
!
do 900 ifa = 1 , nbfac
!
ipf(1:3) = intfac(3:5, ifa)
!!
!   ustar(1, ipf(3)) = sin(pi*coord(1,ipf(3)))*cos(pi*coord(2,ipf(3)))
!    ustar(2, ipf(3)) =-cos(pi*coord(1,ipf(3)))*sin(pi*coord(2,ipf(3)))
!
!
if(coord(1, ipf(1)).lt.1.d-6.or.abs(coord(1, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf',ipf(1)
ustar(1, ipf(1)) = 0.d0
endif
if(coord(2, ipf(1)).lt.1.d-6.or.abs(coord(2, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ipf(1)
ustar(2, ipf(1)) = 0.d0
endif
if(coord(1, ipf(3)).lt.1.d-6.or.abs(coord(1, ipf(3))-1.d0).lt.1.d-6) then
!    print*,'ipf',ipf(1)
ustar(1, ipf(3)) = 0.d0
endif
if(coord(2, ipf(3)).lt.1.d-6.or.abs(coord(2, ipf(3))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ipf(1)
ustar(2, ipf(3)) = 0.d0
endif
900 enddo
!
ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
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
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
unknv = 0.d0
!
!...Low order nodes...
!
do iv = 1, 3
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
fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
fstar(1, 3, iv, ie) = snsigml(1, 3, iv, ie) + munacl(3, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 3, iv, ie) = snsigml(2, 3, iv, ie) + munacl(3, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
fstar(1, 4, iv, ie) = snsigml(1, 4, iv, ie) + munacl(4, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 4, iv, ie) = snsigml(2, 4, iv, ie) + munacl(4, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*ustar(1, ip(iv))- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*ustar(2, ip(iv))- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*ustar(1, ip(iv))- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*ustar(2, ip(iv))- munaul(2, 2, iv, ie)
!
enddo
!
!...High-order nodes
!
do iv = 4, nvtri
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
fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*ustar(1, ip(iv))- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*ustar(2, ip(iv))- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*ustar(1, ip(iv))- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*ustar(2, ip(iv))- munaul(2, 2, iv, ie)
!
enddo
enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
deallocate (muutg, mutag)
end subroutine getndvelo_lag_curv4!
!
!...subroutine: Calculate the nodal velocity U_p^* for curved cell...
!
subroutine getndvelo_lag_curv7(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fstar)
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
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri, 1:nelem),  intent(out)::fstar !...Riemann forces
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin
integer::iel,ier
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)

!...local real array
real*8::munacn_rie, munacu_rie(1:2), snsigm_rie(1:2)
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:nvtri)::murie
real*8::vnorm(1:3, 1:4, 1:nvtri)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8::aujmp(1:2, 1:nvtri)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:nvtri):: xv, yv
!
!...local real number
!
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::uvtx,vvtx,evtx, pvtx
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
allocate (munacl(1:4, 1:nvtri, 1:nelem), munaul(1:ndimn, 1:4, 1:nvtri,  1:nelem),&
snsigml(1:ndimn, 1:4,  1:nvtri,  1:nelem))
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
!...Mark the boundary nodes for case 2...
!
if(ncase.eq.2)then
do ifa = 1, nbfac
   ipf(1:nvfac) = intfac(3:2+nvfac, ifa)
   indnd(ipf(1:3)) = 1
enddo
endif
!
!...shape functions
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0;  yv(1) = 0.d0
xv(2) = 1.d0;  yv(2) = 0.d0
xv(3) = 0.d0;  yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...Part I: Get the averaged nodal velocity...
!
do 200 ie = 1,nelem !...(1)ie = 1,nelem
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
! if(ip(iv)==36) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
enddo
!
!...Accumulate nodal velocity...
!
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
!...Part II:Get the related coefficients for Riemann-like
!
!
!...Zero out munacn, munacu, snsigm...
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
do 250 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 3, 1) = gelag(1:3, 9, ie); vnorm(1:3, 4, 1) = gelag(1:3, 4, ie) !...For point ip(1)
!
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 5, ie); vnorm(1:3, 4, 2) = gelag(1:3, 6, ie) !...For point ip(2)
!
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 7, ie); vnorm(1:3, 4, 3) = gelag(1:3, 8, ie) !...For point ip(3)
!
vnorm(1:3, 1, 4) = gelag(1:3, 4, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(4)
!
vnorm(1:3, 1, 5) = gelag(1:3, 6, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(5)
!
vnorm(1:3, 1, 6) = gelag(1:3, 8, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(6)
!
!...Loop 1: Calculate the tensor stress and a_c at vertices...
!
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
!...Stress tensor...
!
rho  = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
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
!
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
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
!
!
!...Loop 2: Get the impedence at the vertices...
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
murie(iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ip(iv).eq.5) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ip(iv)),unknv(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, 3 !...Low-order nodes...
do ifa = 1, 4 !...Every corner consists of 4 linear faces...
!
!...Call subroutine....
!
call getriecoef(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), unknv(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
if(ifa.le.2)then
!

!
elseif(ifa.ge.3)then
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
munacn(ip(iv)) = munacn(ip(iv)) + munacn_rie
!
munacu(1:2, ip(iv)) = munacu(1:2, ip(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ip(iv)) = snsigm(1:2, ip(iv)) + snsigm_rie(1:2)!
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
endif
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Local variable...
!
munacl(ifa, iv, ie) =  munacn_rie
!
munaul(1:2, ifa, iv, ie) =  munacu_rie(1:2)
!
snsigml(1:2, ifa, iv, ie)=  snsigm_rie(1:2)
!
enddo
!
enddo
!
!...Get the summed denominator cooefficients sum(mu*n*a_c) for high-order nodes...
!
do iv  = 4, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.281) print*,'p19 muacn(281) pre++', munacn(281),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
call getriecoef(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), unknv(2:3, iv), sigma(1:2, 1:2, iv), &
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
munacn(ip(iv)) = munacn(ip(iv)) + munacn_rie
!
munacu(1:2, ip(iv)) = munacu(1:2, ip(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ip(iv)) = snsigm(1:2, ip(iv)) + snsigm_rie(1:2)!
!
!
!
!   if(ip(iv).eq.281) print*,'p281 muacn(vv) post',ie,murie(iv),munacn(ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.281) print*,'p19 muacn(28) prep---',murie(iv), munacn(281), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))

! if(ip(iv).eq.15) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munacl(ifa, iv, ie) =  munacn_rie
!
munaul(1:2, ifa, iv, ie) =  munacu_rie(1:2)
!
snsigml(1:2, ifa, iv, ie)=  snsigm_rie(1:2)
!
enddo

enddo
!
250 enddo  !...(1)ie = 1,nelem
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
!...Zero out the bnorm
!
!
!...
!...Zero out boundary pressure...
!
!
!...Fourth part: Solve the nodal velocity...
!
!...deactivate the pseudo pressure....
!
fpres = 0.d0
!
!...4.1: Update the Riemann forces at every node...
!
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
endif
enddo
!
!
call getvelo_mpt(ustar,gelag,intfac,inpoel,coord,munacn,vlave,unkno)
!
!....Bd velocity
!
! print*,'ustar--',ustar(1:2, 281),munacu(1:2,281) ,snsigm(1:2,281), munacn(281)
!
!
if(ncase.eq.1)then
!
do 900 ifa = 1 , nbfac
!
ipf(1:nvfac) = intfac(3:2+nvfac, ifa)
!
!    ustar(1, ipf(3)) = sin(pi*coord(1,ipf(3)))*cos(pi*coord(2,ipf(3)))
!    ustar(2, ipf(3)) =-cos(pi*coord(1,ipf(3)))*sin(pi*coord(2,ipf(3)))
!
!
if(coord(1, ipf(1)).lt.1.d-6.or.abs(coord(1, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf',ipf(1)
ustar(1, ipf(1)) = 0.d0
endif
if(coord(2, ipf(1)).lt.1.d-6.or.abs(coord(2, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ipf(1)
ustar(2, ipf(1)) = 0.d0
endif
if(coord(1, ipf(3)).lt.1.d-6.or.abs(coord(1, ipf(3))-1.d0).lt.1.d-6) then
!    print*,'ipf',ipf(1)
ustar(1, ipf(3)) = 0.d0
endif
if(coord(2, ipf(3)).lt.1.d-6.or.abs(coord(2, ipf(3))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ipf(1)
ustar(2, ipf(3)) = 0.d0
endif
900 enddo
!
!...Fix the 4 vertices....
!
ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
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
unknv = 0.d0
!
!...Low order nodes...
!
do iv = 1, 3
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
do ifa = 1, 4
!
fstar(1, ifa, iv, ie) = snsigml(1, ifa, iv, ie) + munacl(ifa, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
fstar(2, ifa, iv, ie) = snsigml(2, ifa, iv, ie) + munacl(ifa, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
enddo
!
enddo
!
!...High-order nodes
!
do iv = 4, nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
enddo
enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
end subroutine getndvelo_lag_curv7
!
!
!
!
subroutine getriecoef(murie, area, vnorm, aujmp, velo, sigma, munacn, munacu, snsigm)
use constant
implicit none
real*8, intent(in):: murie
real*8, intent(in):: area
real*8, intent(in):: vnorm(ndimn), aujmp(ndimn), velo(ndimn)
real*8, intent(in):: sigma(ndimn, ndimn)
!
real*8, intent(out):: munacn
real*8, intent(out):: munacu(ndimn), snsigm(ndimn)
!
!...
!
real*8:: vnau
!
vnau   =  vnorm(1)*aujmp(1) + vnorm(2)*aujmp(2)
munacn =  murie*area*abs(vnau)
!
munacu(1:ndimn) = munacn*velo(1:ndimn)
!
snsigm(1) = sigma(1, 1)*area*vnorm(1) + sigma(1, 2)*area*vnorm(2)!
snsigm(2) = sigma(2, 1)*area*vnorm(1) + sigma(2, 2)*area*vnorm(2)
!
end subroutine getriecoef
!
!...subroutine: Calculate the nodal velocity U_p^* for curved cell...
!
subroutine getndvelo_lag_curv5(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fstar)
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
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri, 1:nelem),  intent(out)::fstar !...Riemann forces
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:3)     :: ipf
integer::indnd(npoin)
integer::idcrv(npoin)

!...local real array
real*8::muutg(1:2, 1:npoin), mutag(npoin)
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:nvtri)::murie
real*8::vnorm(1:3, 1:4, 1:nvtri)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8::aujmp(1:2, 1:nvtri)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:nvtri):: xv, yv
real*8::usnom(2), ustng(2)
!...local real number
real*8::unorm, tfx, tfy, nfx, nfy,signm, utang
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::uvtx,vvtx,evtx, pvtx
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
allocate (munacl(1:4, 1:nvtri, 1:nelem), munaul(1:ndimn, 1:4, 1:nvtri,  1:nelem),&
snsigml(1:ndimn, 1:4,  1:nvtri,  1:nelem))
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
ipf(1:3) = intfac(3:5, ifa)
indnd(ipf(1:3)) = 1
enddo
endif
!
!...Coloring the high-order nodes...
!
idcrv = 0
!
do ie =1, nelem
ip(1:nvtri) = inpoel(1:nvtri,ie)
idcrv(ip(4:nvtri)) = 1
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
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0;  yv(1) = 0.d0
xv(2) = 1.d0;  yv(2) = 0.d0
xv(3) = 0.d0;  yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
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
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
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
vnorm(1:3, 3, 1) = gelag(1:3, 9, ie); vnorm(1:3, 4, 1) = gelag(1:3, 4, ie) !...For point ip(1)
!
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 5, ie); vnorm(1:3, 4, 2) = gelag(1:3, 6, ie) !...For point ip(2)
!
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 7, ie); vnorm(1:3, 4, 3) = gelag(1:3, 8, ie) !...For point ip(3)
!
vnorm(1:3, 1, 4) = gelag(1:3, 4, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(3)
!
vnorm(1:3, 1, 5) = gelag(1:3, 6, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(3)
!
vnorm(1:3, 1, 6) = gelag(1:3, 8, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(3)

!
!...ndA=0.5d0*vnorm
!
!vnorm(3, :, :) = 0.5d0*vnorm(3, :, :)
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
evtx = unknv(4, iv)
!
pvtx = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
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
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!...Validation...
!
sigma(1, 1, 4) = 0.5d0*(sigma(1, 1, 1) + sigma(1, 1, 2))
sigma(2, 2, 4) = 0.5d0*(sigma(2, 2, 1) + sigma(2, 2, 2))
!
sigma(1, 1, 5) = 0.5d0*(sigma(1, 1, 2) + sigma(1, 1, 3))
sigma(2, 2, 5) = 0.5d0*(sigma(2, 2, 2) + sigma(2, 2, 3))
!
sigma(1, 1, 6) = 0.5d0*(sigma(1, 1, 1) + sigma(1, 1, 3))
sigma(2, 2, 6) = 0.5d0*(sigma(2, 2, 1) + sigma(2, 2, 3))

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
!if(ip(iv).eq.23) print*,'murie22', sdctr,rhoct,deltu,vlave(1:2, ip(iv)),unknv(2:3, iv),unkno(1:3,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, 3 !...Low-order nodes...
do ifa = 1, 4 !...Every corner consists of 2 faces...
!
if(ifa.le.2)then
!

!
elseif(ifa.ge.3)then
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
munacn(ip(iv)) = munacn(ip(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
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
munacu(1, ip(iv)) =  munacu(1, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(2, iv)
munacu(2, ip(iv)) =  munacu(2, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(3, iv)
!if(ip(iv).eq.23) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv,&
!vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
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
!
! if(ip(iv).eq.23) print*,'p19 muacn(28) post-snsigmacurv5',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv&
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
endif

!
!...Local variable...
!
munacl(ifa, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
!
munaul(1:2, ifa, iv, ie)    =  munacl(ifa, iv, ie)*unknv(2:3, iv)
!
snsigml(1, ifa, iv, ie)= sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv)

snsigml(2, ifa, iv, ie)= sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!

!
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
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
!     snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv)

!     snsigml(2, 1, iv, ie)= sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
!     snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv)

!     snsigml(2, 2, iv, ie)= sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)!
enddo
!
!...Get the summed denominator cooefficients sum(mu*n*a_c) for high-order nodes...
!
do iv  = 4, nvtri
do ifa = 1, 1 !...Every corner consists of 2 faces...
!
if(iv.eq.4) vnorm(1:3, ifa ,iv) = gelag(1:3, 1, ie)
if(iv.eq.5)  vnorm(1:3, ifa ,iv) = gelag(1:3, 3, ie)
if(iv.eq.6)  vnorm(1:3, ifa ,iv) = gelag(1:3, 3, ie)
!
! if(ip(iv).eq.281) print*,'p19 muacn(281) pre++', munacn(281),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
munacn(ip(iv)) = munacn(ip(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
!   if(ip(iv).eq.281) print*,'p281 muacn(vv) post',ie,murie(iv),munacn(ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.281) print*,'p19 muacn(28) prep---',murie(iv), munacn(281), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
munacu(1, ip(iv)) =  munacu(1, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(2, iv)
munacu(2, ip(iv)) =  munacu(2, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(3, iv)
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
munacl(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)* &
abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))
!
munacl(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)* &
abs(vnorm(1, 2, iv)*aujmp(1, iv) + vnorm(2, 2, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
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
!
!...Zero out the bnorm
!
!
!...
!...Zero out boundary pressure...
!
!
!...Fourth part: Solve the nodal velocity...
!
!...deactivate the pseudo pressure....
!
fpres = 0.d0
!
!...4.1: Update the Riemann forces at every node...
!
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
!if(idcrd(ipoin).eq.0)then
ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
!elseif(idcrd(ipoin).eq.1)then
!  if(muacn(ipoin).gt.1.d-9)then
!   ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
!   ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
!  endif
!endif
endif
enddo
!
!...Deal with the high-order nodes...
!
!
muutg = 0.d0
mutag = 0.d0
!
do 350 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
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
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
do iv =1 ,nvtri
!...Basis function
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds

enddo
!
!...Impedence...
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
!...
!...zero out unknv
unknv = 0.d0
!
do iv   = 4,nvtri
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
enddo
!
!
do iv  = 4, nvtri
muutg(1:2, ip(iv)) = muutg(1:2, ip(iv)) + murie(iv)*unknv(2:3, iv)
mutag(ip(iv)) = mutag(ip(iv)) + murie(iv)
enddo
!
350 enddo
!
!
!
do 450 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 4) = gelag(1:3, 4, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(3)
!
vnorm(1:3, 1, 5) = gelag(1:3, 6, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(3)
!
vnorm(1:3, 1, 6) = gelag(1:3, 8, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(3)
!
do iv = 4, nvtri
!
if(munacn(ip(iv)).lt.1.d-9)then
!
nfx = vnorm(1, 1, iv); nfy = vnorm(2, 1, iv)
!
tfx = nfy ; tfy = -nfx
!
unorm = munacu(1, ip(iv))*nfx + munacu(2, ip(iv))*nfy
signm = snsigm(1, ip(iv))*nfx + snsigm(2, ip(iv))*nfy
!
usnom(1) = (unorm*nfx - signm*nfx)/(munacn(ip(iv)) + 1.d-9)
usnom(2) = (unorm*nfy - signm*nfy)/(munacn(ip(iv)) + 1.d-9)
!
usnom = 0.d0
!
utang = muutg(1, ip(iv))*tfx + muutg(2, ip(iv))*tfy
!
ustng(1) = utang*tfx/mutag(ip(iv))
ustng(2) = utang*tfy/mutag(ip(iv))
!
ustar(1 ,ip(iv)) = usnom(1) + ustng(1)
ustar(2, ip(iv)) = usnom(2) + ustng(2)
!
endif
!
!if(ip(iv).eq.281) print*, 'ustar', usnom, ustng
enddo
!
450 enddo

!
!....Bd velocity
!
!print*,'ustar--',ustar(1:2, 23),munacu(1:2,23) ,snsigm(1:2,23), munacn(23)
!
!
if(ncase.eq.1)then
!
do 900 ifa = 1 , nbfac
!
ipf(1:3) = intfac(3:5, ifa)
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
if(coord(1, ipf(3)).lt.1.d-6.or.abs(coord(1, ipf(3))-1.d0).lt.1.d-6) then
!    print*,'ipf',ipf(1)
ustar(1, ipf(3)) = 0.d0
endif
if(coord(2, ipf(3)).lt.1.d-6.or.abs(coord(2, ipf(3))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ipf(1)
ustar(2, ipf(3)) = 0.d0
endif
900 enddo
!
ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
!
!...validation...
!
do 460 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
ustar(1:2,ip(4)) = 0.5d0*(ustar(1:2,ip(1)) + ustar(1:2,ip(2)))
ustar(1:2,ip(5)) = 0.5d0*(ustar(1:2,ip(2)) + ustar(1:2,ip(3)))
ustar(1:2,ip(6)) = 0.5d0*(ustar(1:2,ip(1)) + ustar(1:2,ip(3)))
!
!if(ip(iv).eq.281) print*, 'ustar', usnom, ustng
!
460 enddo
!
!  print*,'ustar',ustar(1:2, 23)!,ustar(1:2, 20),ustar(1:2, 21)
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
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
unknv = 0.d0
!
!...Low order nodes...
!
do iv = 1, 3
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
fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
fstar(1, 3, iv, ie) = snsigml(1, 3, iv, ie) + munacl(3, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 3, iv, ie) = snsigml(2, 3, iv, ie) + munacl(3, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
fstar(1, 4, iv, ie) = snsigml(1, 4, iv, ie) + munacl(4, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 4, iv, ie) = snsigml(2, 4, iv, ie) + munacl(4, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*ustar(1, ip(iv))- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*ustar(2, ip(iv))- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*ustar(1, ip(iv))- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*ustar(2, ip(iv))- munaul(2, 2, iv, ie)
!
enddo
!
!...High-order nodes
!
do iv = 4, nvtri
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
fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*ustar(1, ip(iv))- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*ustar(2, ip(iv))- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*ustar(1, ip(iv))- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*ustar(2, ip(iv))- munaul(2, 2, iv, ie)
!
enddo
!
!...Linear validation...
!
fstar(1:2, 1, 4, ie) = 0.5d0*(fstar(1:2, 4, 1, ie) + fstar(1:2, 3, 2, ie))
fstar(1:2, 2, 4, ie) = fstar(1:2, 1, 4, ie)
!
fstar(1:2, 1, 5, ie) = 0.5d0*(fstar(1:2, 4, 2, ie) + fstar(1:2, 3, 3, ie))
fstar(1:2, 2, 5, ie) = fstar(1:2, 1, 5, ie)
!
fstar(1:2, 1, 6, ie) = 0.5d0*(fstar(1:2, 3, 1, ie) + fstar(1:2, 4, 3, ie))
fstar(1:2, 2, 6, ie) = fstar(1:2, 1, 6, ie)

!
!fstar(1:2, 1, 4, ie) = fstar(1:2, 4, 1, ie)
!fstar(1:2, 2, 4, ie) = fstar(1:2, 3, 2, ie)
!
!fstar(1:2, 1, 5, ie) = fstar(1:2, 4, 2, ie)
!fstar(1:2, 2, 5, ie) = fstar(1:2, 3, 3, ie)
!
!fstar(1:2, 1, 6, ie) = fstar(1:2, 4, 3, ie)
!fstar(1:2, 2, 6, ie) = fstar(1:2, 3, 1, ie)

!
enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
end subroutine getndvelo_lag_curv5!
!
!...Face integral for curved...
!
subroutine rhsifacedg_lag_curv5(inpoel,  unkno, ustar, fstar, gelag,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri,1:nelem),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nelem),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac),    intent(in)::gelag
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:4, 1:nvtri) :: ndshp
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8, dimension(1:2, 1:ndegr, 1:4):: lpnp
real*8, dimension(1:ndegr, 1:4)::clpnp
real*8::xv(nvtri), yv(nvtri),b(1:3,1:nvtri)
real*8::vnorm(1:3, 1:4, 1:nvtri)
!...local real number
real*8::eps,c00,c05,c10,c20,c13,c23
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
data c23   / 0.6666666666666666d0 /
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
!...Nodes for one vertex
!
ndshp(1, 1) = 3; ndshp(2, 1) = 6; ndshp(3, 1) = 2; ndshp(4, 1) = 4
ndshp(1, 2) = 1; ndshp(2, 2) = 4; ndshp(3, 2) = 3; ndshp(4, 2) = 5;
ndshp(1, 3) = 2; ndshp(2, 3) = 5; ndshp(3, 3) = 1; ndshp(4, 3) = 6;
ndshp(1, 4) = 1; ndshp(2, 4) = 2
ndshp(1, 5) = 2; ndshp(2, 5) = 3
ndshp(1, 6) = 3; ndshp(2, 6) = 1

!
do iv =1 , nvtri !...low order nodes...
!
!print*,'iv',ivt
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!

do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 3, 1) = gelag(1:3, 9, ie); vnorm(1:3, 4, 1) = gelag(1:3, 4, ie) !...For point ip(1)
!
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 5, ie); vnorm(1:3, 4, 2) = gelag(1:3, 6, ie) !...For point ip(2)
!
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 7, ie); vnorm(1:3, 4, 3) = gelag(1:3, 8, ie) !...For point ip(3)
!
vnorm(1:3, 1, 4) = gelag(1:3, 4, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(3)
!
vnorm(1:3, 1, 5) = gelag(1:3, 6, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(3)
!
vnorm(1:3, 1, 6) = gelag(1:3, 8, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(3)
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
!...Distribute to every corner...
!
do iv = 1, 3
!
do ideg = 1, ndegr
lpnp(1:2, ideg, 3) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(2, iv)))*vnorm(3, 3, iv)*vnorm(1:2, 3, iv)
lpnp(1:2, ideg, 1) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(1, iv)))*vnorm(3, 1, iv)*vnorm(1:2, 1, iv)
!
lpnp(1:2, ideg, 4) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(4, iv)))*vnorm(3, 4, iv)*vnorm(1:2, 4, iv)
lpnp(1:2, ideg, 2) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(3, iv)))*vnorm(3, 2, iv)*vnorm(1:2, 2, iv)
!
clpnp(ideg, 3) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(2, iv)))
clpnp(ideg, 1) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(1, iv)))
clpnp(ideg, 4) = 0.1d0*(6.d0*b(ideg, iv) + 4.d0*b(ideg, ndshp(4, iv)))
clpnp(ideg, 2) = 0.1d0*(     b(ideg, iv) -      b(ideg, ndshp(3, iv)))
enddo

do ifa =1, 4
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr) + c13*ustar(1, ip(iv))*lpnp(1, 1:ndegr, ifa) +&
c13*ustar(2, ip(iv))*lpnp(2, 1:ndegr, ifa)
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr) + c13*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa)

plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr) + c13*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
c13*ustar(1, ip(iv))*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa) +&
c13*ustar(2, ip(iv))*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
!if(ie==50) print*,'rhs iface interplus', ip(iv),iv, plnpn(1, 1:3),fstar(1, ifa, iv, ie),clpnp(1, ifa),&
!                                         c13*fstar(1, ifa, iv, ie)
!
enddo
!

!
enddo
!
do iv = 4, nvtri
!
do ideg = 1, ndegr
lpnp(1:2, ideg, 1) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(1, iv)))*vnorm(3, 1, iv)*vnorm(1:2, 1, iv)
lpnp(1:2, ideg, 2) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(2, iv)))*vnorm(3, 2, iv)*vnorm(1:2, 2, iv)
!
clpnp(ideg, 1) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(1, iv)))
clpnp(ideg, 2) = 0.2d0*(4.d0*b(ideg, iv) +   b(ideg, ndshp(2, iv)))

enddo

do ifa =1, 2
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr) + c23*ustar(1, ip(iv))*lpnp(1, 1:ndegr, ifa) +&
c23*ustar(2, ip(iv))*lpnp(2, 1:ndegr, ifa)
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr) + c23*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa)
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr) + c23*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
!elnpn(1:ndegr)   = elnpn(1:ndegr)+&
!c23*ustar(1, ip(iv))*fstar(1, ifa, iv, ie)*clpnp(1:ndegr, ifa) +&
!c23*ustar(2, ip(iv))*fstar(2, ifa, iv, ie)*clpnp(1:ndegr, ifa)
!
!
!... validation for linear problems....
!
if(iv.eq.4)then
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
c23*0.5d0*(ustar(1, ip(1))*fstar(1, 4, 1, ie) + ustar(1, ip(2))*fstar(1, 3, 2, ie))*clpnp(1:ndegr, ifa) +&
c23*0.5d0*(ustar(2, ip(1))*fstar(2, 4, 1, ie) + ustar(2, ip(2))*fstar(2, 3, 2, ie))*clpnp(1:ndegr, ifa)
!
elseif(iv.eq.5)then
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
c23*0.5d0*(ustar(1, ip(2))*fstar(1, 4, 2, ie) + ustar(1, ip(3))*fstar(1, 3, 3, ie))*clpnp(1:ndegr, ifa) +&
c23*0.5d0*(ustar(2, ip(2))*fstar(2, 4, 2, ie) + ustar(2, ip(3))*fstar(2, 3, 3, ie))*clpnp(1:ndegr, ifa)
!
elseif(iv.eq.6)then

elnpn(1:ndegr)   = elnpn(1:ndegr)+&
c23*0.5d0*(ustar(1, ip(3))*fstar(1, 4, 3, ie) + ustar(1, ip(1))*fstar(1, 3, 1, ie))*clpnp(1:ndegr, ifa) +&
c23*0.5d0*(ustar(2, ip(3))*fstar(2, 4, 3, ie) + ustar(2, ip(1))*fstar(2, 3, 1, ie))*clpnp(1:ndegr, ifa)
endif
!
!...
!
!if(ie==50) print*,'rhs iface interplus', ip(iv), iv, elnpn(1:3),fstar(1, ifa, iv, ie),clpnp(1, ifa), &
!                                         c23*fstar(1, ifa, iv, ie)
!
enddo
!
enddo
!
rhsel(1:ndegr, 1, ie) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ie) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ie) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ie) =  elnpn(1:ndegr)
!
!if(ie==50) print*,'rhs iface low order',rhsel(1:3, 4, ie),ustar(1:2,ip(1)), &
!ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:6)
!
! if(ie==83) print*,'rhs iface',rhsel(1, 1, ie), lpnp(1:2, 1,  1:2),ustar(1:2,ip(4)), &
!                               ustar(1:2,ip(5)),ustar(1:2,ip(6)), ip(1:6)

550 enddo
!
!   open(8,file='lpnp.dat')
!    do ie = 1, nelem8
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
end subroutine rhsifacedg_lag_curv5
!
!...subroutine: Calculate the nodal velocity U_p^*...
!
subroutine getndvelo_lag_linear(gflag,gelag,bface,intfac,inpoel,ltelem,coord,unkno,ustar, fstar)
use constant
implicit none
!...Input arrays
real*8,dimension(1:3, 1:ngelg, 1:nelem+nbfac), intent(in)::gelag
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
integer, dimension(3,nelem),                 intent(in)::ltelem
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri, 1:nelem),  intent(out)::fstar !...Riemann forces
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvfac) :: ipf
integer,dimension(1:3)     :: estri
integer::indnd(npoin)

!...local real array
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:nvtri)::murie
real*8::vnorm(1:3, 1:2, 1:3)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8::aujmp(1:2, 1:nvtri)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:nvtri):: xv, yv
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::uvtx,vvtx,evtx, pvtx
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
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
indnd(ipf(1:nvfac)) = 1
enddo
endif
!
!...Shape functions for reference triagnle...
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
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
!...Part 1: Get the averaged velocity at nodes...
!
do 200 ie = 1,nelem !...(1)ie = 1,nelem
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
! if(ip(iv)==15) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
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
!...Part II:Calculate the nodal coefficients for RIemann-like...
!
!...Zero out munacn , munacu, snsigm...
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
do 250 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...Find elements surounding elements
!
estri(1:3) = ltelem(1:3, ie)
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
!...zero out unknv...
!
unknv = 0.d0
!
!...Get the numerical value at the vertex...
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
enddo
!
!...Get the Riemann-like coefficients...
!
do iv   = 1,nvtri
!
rho  = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
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
!
!aujmp(1:2, iv) = usold(1:2, ip(iv)) - unknv(2:3, iv)
aujmp(1:2, iv) = vlave(1:2, ip(iv)) - unknv(2:3, iv)
!if(ie==22) print*,'adjumpxxx22', vlave(1:2, ip(iv)) , unknv(2:3, iv), ip(iv)
acnx = aujmp(1, iv)
acny = aujmp(2, iv)
if(sqrt(acnx**2 + acny**2).lt.1.e-12)then
!aujmp(1:2, iv) = 0.d0;
aujmp(1:2, iv) = 1.d-12
!print*,'point are reset', ip(iv), aujmp(1:2, iv),ie,iv
else
aujmp(1:2, iv) = aujmp(1:2, iv)/sqrt(acnx**2 + acny**2)
endif
!if(ip(iv)==6878) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, iv),sqrt(acnx**2 + acny**2)
!vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
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
if(ip(iv).eq.23) print*,'murie22', sdctr,rhoct,deltu,vlave(1:2, ip(iv)),unknv(2:3, iv),unkno(1:3,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.6878) print*,'p19 muacn(28) pre++', munacn(ip(iv)),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
!  if(sqrt(aujmp(1, iv)**2 + aujmp(2, iv)**2).lt.1.e-12)then
!     aujmp(1:2, iv) = vnorm(1:2, ifa, iv)
!  endif
!
munacn(ip(iv)) = munacn(ip(iv)) + murie(iv)*vnorm(3, ifa, iv)* &
abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
!   if(ip(iv).eq.6878) print*,'p9741 muacn(vv) post',ie,murie(iv),munacn(ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p19 muacn(28) prep---',murie(iv), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                 vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))
munacu(1, ip(iv)) =  munacu(1, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(2, iv)
munacu(2, ip(iv)) =  munacu(2, ip(iv)) +&
murie(iv)*vnorm(3, ifa, iv)*abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))*unknv(3, iv)
if(ip(iv).eq.23) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv,&
vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Get the summed stress sum(n*A*sigma)
!
!
snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
sigma(2, 1, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!
snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(1, 2, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) + &
sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!

! snsigm(1, ip(iv)) = snsigm(1, ip(iv)) + sigma(1, 1, iv)*vnorm(3, ifa, iv)*vnorm(1, ifa, iv) !
! snsigm(2, ip(iv)) = snsigm(2, ip(iv)) + sigma(2, 2, iv)*vnorm(3, ifa, iv)*vnorm(2, ifa, iv)
!if(ip(iv).eq.23) print*,'p19 muacn(23) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv&
!            ,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
!  if(sqrt(aujmp(1, iv)**2 + aujmp(2, iv)**2).lt.1.e-12)then
!     aujmp(1:2, iv) = vnorm(1:2, 1, iv)
!  endif
!
munacl(1, iv, ie) =  murie(iv)*vnorm(3, 1, iv)* &
abs(vnorm(1, 1, iv)*aujmp(1, iv) + vnorm(2, 1, iv)*aujmp(2, iv))
!
!  if(sqrt(aujmp(1, iv)**2 + aujmp(2, iv)**2).lt.1.e-12)then
!     aujmp(1:2, iv) = vnorm(1:2, 2, iv)
!  endif
!
munacl(2, iv, ie) =  murie(iv)*vnorm(3, 2, iv)* &
abs(vnorm(1, 2, iv)*aujmp(1, iv) + vnorm(2, 2, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
!
munaul(1, 1, iv, ie)    =  munacl(1, iv, ie)*unknv(2, iv)
munaul(2, 1, iv, ie)    =  munacl(1, iv, ie)*unknv(3, iv)

munaul(1:2, 2, iv, ie)    =  munacl(2, iv, ie)*unknv(2:3, iv)
!
snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
sigma(2, 1, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
snsigml(2, 1, iv, ie)= sigma(1, 2, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv) + &
sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
sigma(2, 1, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
snsigml(2, 2, iv, ie)= sigma(1, 2, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv) + &
sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)
!

!
!     snsigml(1, 1, iv, ie)= sigma(1, 1, iv)*vnorm(3, 1, iv)*vnorm(1, 1, iv)

!     snsigml(2, 1, iv, ie)= sigma(2, 2, iv)*vnorm(3, 1, iv)*vnorm(2, 1, iv)
!
!     snsigml(1, 2, iv, ie)= sigma(1, 1, iv)*vnorm(3, 2, iv)*vnorm(1, 2, iv)

!     snsigml(2, 2, iv, ie)= sigma(2, 2, iv)*vnorm(3, 2, iv)*vnorm(2, 2, iv)!
enddo
!
250 enddo  !...(1)ie = 1,nelem
!

!
!...Third part: Impose the boundary condition
!
!...Impose BC with pressure prescribed boundary...
!
!call getbcfc_lag(bface, intfac, gflag, fpres)
!
 fpres= 0.d0
!
!...3.1: Update the Riemann forces at every node...
!
do ipoin = 1, npoin
!
!   print*,'good',ipoin
!
if(indnd(ipoin).eq.0)then !Excluding the boundaries with prescribed velocity..
ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
endif
!
enddo
!
print*,'ustar--',ustar(1:2, 23),munacu(1:2,23) ,snsigm(1:2, 23), munacn(23)
!
!...Special treatment for Taylor-Green vortex boundary nodes...
!
if(ncase.eq.1)then
!
 do 900 ifa = 1 , nbfac
!
   ipf(1:2) = intfac(3:4, ifa)
!
!!    ustar(1, ipf(1)) = sin(pi*coord(1,ipf(1)))*cos(pi*coord(2,ipf(1)))
!!    ustar(2, ipf(1)) =-cos(pi*coord(1,ipf(1)))*sin(pi*coord(2,ipf(1)))
!
!!    ustar(1, ipf(2)) = sin(pi*coord(1,ipf(2)))*cos(pi*coord(2,ipf(2)))
!!    ustar(2, ipf(2)) =-cos(pi*coord(1,ipf(2)))*sin(pi*coord(2,ipf(2)))
!
   if(coord(1, ipf(1)).lt.1.d-6.or.abs(coord(1, ipf(1))-1.d0).lt.1.d-6) then
! !    print*,'ipf',ipf(1)
      ustar(1, ipf(1)) = 0.d0
   endif
  if(coord(2, ipf(1)).lt.1.d-6.or.abs(coord(2, ipf(1))-1.d0).lt.1.d-6) then
! !    print*,'ipf2',ipf(1)
     ustar(2, ipf(1)) = 0.d0
   endif
900 enddo
endif
!
!...Imposing the zero normal velocity for BC...
!
!call getbcvn_lag(bface, intfac, gflag, ustar)
!
! print*,'ustar',ustar(1:2, 1:4)
!
!
!   print*,'ustar',ustar(1:2, 18),ustar(1:2, 20),ustar(1:2, 21)
!  if(ustar(2, 18).gt.0.d0) then
!   print*,'Wrong movng direction'!,munacu(2, 28) - snsigm(2, 28) , fpres(2, 28),munacn(28)
!   stop
!  endif

!
!...Fourth part: Update the Riemann forces at every node...
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...Zero out unknown vector at nodes...
!
unknv = 0.d0
!
do iv = 1, nvtri
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 1, iv, ie)
fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 1, iv, ie)
!
fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv)) !- munaul(1, 2, iv, ie)
fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(2, 2, iv, ie)
!
!fstar(1, 1, iv, ie) = snsigml(1, 1, iv, ie) + munacl(1, iv, ie)*ustar(1, ip(iv))- munaul(1, 1, iv, ie)
!fstar(2, 1, iv, ie) = snsigml(2, 1, iv, ie) + munacl(1, iv, ie)*ustar(2, ip(iv))- munaul(2, 1, iv, ie)
!
!fstar(1, 2, iv, ie) = snsigml(1, 2, iv, ie) + munacl(2, iv, ie)*ustar(1, ip(iv))- munaul(1, 2, iv, ie)
!fstar(2, 2, iv, ie) = snsigml(2, 2, iv, ie) + munacl(2, iv, ie)*ustar(2, ip(iv))- munaul(2, 2, iv, ie)
!
enddo
!
enddo
!
deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
end subroutine getndvelo_lag_linear
!
!...Face integral...
!
subroutine rhsifacedg_lag_linear(inpoel,  unkno, ustar, fstar, lpnp, gelag,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn, 1:ndegr, 1:2, 1:nvtri, 1:nelem),intent(in)::lpnp
real*8,dimension(1:ndimn,1:2,1:nvtri,1:nelem),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nelem),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(in)::gelag
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa
integer::ip1,ip2
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2, 1:nvtri) :: ipf
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::xv(3), yv(3),b(1:3,1:nvtri)
!...local real number
real*8::eps,c00,c05,c10,c20,c13
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
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
do iv =1 ,nvtri
!
!print*,'iv',ivt
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
ipf(1, 1) = 3; ipf(2, 1) = 2
ipf(1, 2) = 1; ipf(2, 2) = 3
ipf(1, 3) = 2; ipf(2, 3) = 1

do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
!...The vertex constituting one cell...
!
ip(1:nvtri) = inpoel(1:nvtri, ie)
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
ustar(1, ip(iv))*lpnp(1, 1:ndegr, 1, iv, ie) +&
ustar(2, ip(iv))*lpnp(2, 1:ndegr, 1, iv, ie) +&
ustar(1, ip(iv))*lpnp(1, 1:ndegr, 2, iv, ie) +&
ustar(2, ip(iv))*lpnp(2, 1:ndegr, 2, iv, ie)
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
if(ie==50) print*,'rhs iface inter', ip(iv), elnpn(1:3), ustar(1, ip(iv))*fstar(1, 1, iv, ie),iv
!
enddo
!
rhsel(1:ndegr, 1, ie) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ie) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ie) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ie) =  elnpn(1:ndegr)
!
if(ie==50) print*,'rhs iface',rhsel(1:3, 4, ie),ustar(1:2,ip(1)), &
ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)

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
end subroutine rhsifacedg_lag_linear
!
!
!
subroutine getriecoef_vilar(murie, area, vnorm, aujmp, velo, sigma, munacn, munacu, snsigm)
use constant
implicit none
real*8, intent(in):: murie
real*8, intent(in):: area
real*8, intent(in):: vnorm(1:ndimn), aujmp(1:ndimn), velo(1:ndimn)
real*8, intent(in):: sigma(1:ndimn, 1:ndimn)
!
real*8, intent(out):: munacn(1:ndimn, 1:ndimn)
real*8, intent(out):: munacu(1:ndimn), snsigm(1:ndimn)
!
!...
!
real*8:: vnau,nx,ny
!
nx = vnorm(1)
ny = vnorm(2)
!
munacn(1, 1) = nx*nx; munacn(2, 1) = nx*ny;
munacn(1, 2) = nx*ny; munacn(2, 2) = ny*ny;
!
munacn = munacn*area*murie
!
munacu(1) = munacn(1, 1)*velo(1) + munacn(2, 1)*velo(2)
munacu(2) = munacn(2, 1)*velo(1) + munacn(2, 2)*velo(2)
!
!print*,'munacu sub',munacu, munacn, area,murie,velo
!print*,'mu',murie*(velo(1)*nx + velo(2)*ny)*area*nx
!
snsigm(1) = sigma(1, 1)*area*vnorm(1) + sigma(1, 2)*area*vnorm(2)!
snsigm(2) = sigma(2, 1)*area*vnorm(1) + sigma(2, 2)*area*vnorm(2)
!
end subroutine getriecoef_vilar
!
!...Simpson method...
!

subroutine getndvelo_lag_simpson(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fstar)
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
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:9, 1:nelem),intent(out)::fstar !...Riemann forces
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin, iel, ier
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)
integer::infag(2,3)
integer::fagsp(nvfac, 3),cegsp(nvfac, 3)
integer::gssid(nvfac)

!...local real array
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:nq,1:nvfac)::unkng
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8::vnorm(1:3, 1:9)
real*8,dimension(1:nvfac)::murie
real*8::sigmg(1:2, 1:2, 1:nvfac)
real*8::aujmp(1:2, 1:nvfac)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:ndegr)::bg
real*8,dimension(1:nvtri):: xv, yv
real*8::weigh(ngausf), posi(1,ngausf)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::ug, vg, pg, eg
real*8:: r1, r2, s2,s1, rg, sg, r
real*8:: shp1, shp2, wi
real*8::uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dwav1,dwav2
real*8::anx, any
real*8::nfx, nfy, tfx, tfy,rhsu1,rhsu2,detma
real*8::xpf(1:2,1:3),xcf(1:2)
real*8::len1,len2,tfx1,tfy1
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2),munaci(2,2)
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: cnsup(:), munacn(:,:,:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: usold(:,:), munacu(:,:), snsigm(:,:)
real*8,allocatable:: munacl(:,:,:,:)
real*8,allocatable:: snsigml(:,:,:), munaul(:,:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:2,1:2,1:npoin), cnsup(1:npoin))
allocate (usold(1:ndimn, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munacl(1:2,1:2,1:9, 1:nelem), munaul(1:ndimn, 1:9,  1:nelem),&
snsigml(1:ndimn, 1:9,  1:nelem))
allocate (bnorm(1:3, 1:npoin))
allocate (bpres(1:npoin))
allocate (fpres(1:2, 1:npoin))
!
!...Gauss quadrature...
!
!call ruqope_lobatto(1, ngausf, posi, weigh)
!
!posi(1, 1) = -1.d0
!posi(1, 2) =  1.d0
!
!
if(nvfac.eq.2)then
!
fagsp(1, 1) = 1;  fagsp(2, 1) = 2;
fagsp(1, 2) = 2;  fagsp(2, 2) = 3;
fagsp(1, 3) = 3;  fagsp(2, 3) = 1;
!
cegsp(1, 1) = 1;  cegsp(2, 1) = 2;
cegsp(1, 2) = 3;  cegsp(2, 2) = 4;
cegsp(1, 3) = 5;  cegsp(2, 3) = 6;
!
!print*,'nvfac', nvfac
!
elseif(nvfac.eq.3)then
!
!print*,'nvfac', nvfac
!
fagsp(1, 1) = 1;  fagsp(2, 1) = 2; fagsp(3, 1) = 4;
fagsp(1, 2) = 2;  fagsp(2, 2) = 3; fagsp(3, 2) = 5;
fagsp(1, 3) = 3;  fagsp(2, 3) = 1; fagsp(3, 3) = 6;
!
gssid(1) = 2; gssid(2) = 1; gssid(3)= 1
!
!
cegsp(1, 1) = 1;  cegsp(2, 1) = 2; cegsp(3, 1) = 7;
cegsp(1, 2) = 3;  cegsp(2, 2) = 4; cegsp(3, 2) = 8;
cegsp(1, 3) = 5;  cegsp(2, 3) = 6; cegsp(3, 3) = 9;
!
endif
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
ipf(1:nvfac) = intfac(3:2+nvfac, ifa)
!
indnd(ipf(1:nvfac)) = 1

enddo
endif

!
do 200 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
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
if(ncurv.eq.1)then
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
endif
!

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
if(ncurv.eq.1)then
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
endif
!
do iv =1 ,nvtri
!...Basis function
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds

enddo
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
vnorm(1:3,  1:9) = gelag(1:3, 1:9, ie);
!

!
!...Give the normal vector of every face...
!
do ifa =1, 3
!
!...zero out unkng
unkng = 0.d0
!
do ig   = 1, nvfac
!
!r   = posi(1, ig)
!wi  = weigh(ig)
!
rg = xv(fagsp(ig, ifa))
sg = yv(fagsp(ig, ifa))
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
rho  = 1.d0/unkng(1, ig)
ug = unkng(2, ig)
vg = unkng(3, ig)
eg = unkng(4, ig)
!
pg = max(eps, (gamlg-1.d0)*rho*(eg - 0.5d0*(ug**2 + vg**2)))
!
!pvtx = 0.25d0*(cos(2.d0*pi*coord(1, ip(iv))) + cos(2.d0*pi*coord(2, ip(iv)))) + 1.d0
!
!...Get stress tensor at nodes
!
sigmg(1, 1, ig) = -pg
sigmg(1, 2, ig) = 0.d0
sigmg(2, 1, ig) = 0.d0
sigmg(2, 2, ig) = -pg
!
!
!if(ip(iv).eq.8) print*,'velocity 8',ie, rho, uvtx,vvtx,vlave(1:2, ip(iv))
!
!...Get the a_c (unit vector)
!aujmp(1:2, iv) = usold(1:2, ip(iv)) - unknv(2:3, iv)
aujmp(1:2, ig) = vlave(1:2, ip(fagsp(ig, ifa))) - unkng(2:3, ig)
!
acnx = aujmp(1, ig)
acny = aujmp(2, ig)
if(sqrt(acnx**2 + acny**2).lt.1.e-11)then
aujmp(1:2, ig) = 1.e-11!0.d0;
!print*,'point are reset', ip(iv)
else
aujmp(1:2, ig) = aujmp(1:2, ig)/sqrt(acnx**2 + acny**2)
endif
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
!
!
enddo !....ig
!
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
do ig   = 1, nvfac
dux= vlave(1, ip(fagsp(ig, ifa)))-unkng(2, ig)
duy= vlave(2, ip(fagsp(ig, ifa)))-unkng(3, ig)
deltu = sqrt(dux**2 + duy**2)
murie(ig) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!if(ip(fagsp(ig, ifa)).eq.155) print*,'murie22', sdctr,rhoct,vlave(1:2, ip(fagsp(ig, ifa))),unkng(2:3, ig),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do ig  = 1, nvfac
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
call getriecoef_matrixnew(murie(ig), vnorm(3, cegsp(ig, ifa)), vnorm(1:2, cegsp(ig, ifa)), aujmp(1:2, ig), &
unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ig), vnorm(3, cegsp(ig, ifa)), vnorm(1:2, cegsp(ig, ifa)), aujmp(1:2, ig), &
!unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ip(fagsp(ig, ifa))) = munacn(1:2, 1, ip(fagsp(ig, ifa))) + munacn_rie(1:2, 1)
munacn(1:2, 2, ip(fagsp(ig, ifa))) = munacn(1:2, 2, ip(fagsp(ig, ifa))) + munacn_rie(1:2, 2)
!
munacu(1:2, ip(fagsp(ig, ifa))) = munacu(1:2, ip(fagsp(ig, ifa))) + munacu_rie(1:2)
!
snsigm(1:2,ip(fagsp(ig, ifa))) = snsigm(1:2, ip(fagsp(ig, ifa))) + snsigm_rie(1:2)!
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
munacl(1:2, 1, cegsp(ig, ifa), ie) =  munacn_rie(1:2, 1)
munacl(1:2, 2, cegsp(ig, ifa), ie) =  munacn_rie(1:2, 2)
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munaul(1:2, cegsp(ig, ifa), ie)    =  munacu_rie(1:2)
!
snsigml(1:2, cegsp(ig, ifa), ie)   = snsigm_rie(1:2)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,nelem
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
!...Zero out the bnorm
!
!
!
!...Fourth part: Solve the nodal velocity...
!
!...deactivate the pseudo pressure....
!
fpres = 0.d0
!
!...4.1: Update the Riemann forces at every node...
!
!
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
endif
enddo
!
!...Handle the midlle point...
!
call getvelo_mpt(ustar,gelag,intfac,inpoel,coord,munacn,vlave,unkno)
!
!
!print*,'ustar--',ustar(1:2, 45),munacu(1:2,45) ,snsigm(1:2, 45), munacn(45)
!
!
if(ncase.eq.1)then
!
do 900 ifa = 1 , nbfac
!
ipf(1:nvfac) = intfac(3:2+nvfac, ifa)
!
 ustar(1, ipf(1:nvfac)) = sin(pi*coord(1,ipf(1:nvfac)))*cos(pi*coord(2,ipf(1:nvfac)))
 ustar(2, ipf(1:nvfac)) =-cos(pi*coord(1,ipf(1:nvfac)))*sin(pi*coord(2,ipf(1:nvfac)))
!

900 enddo
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
!
!  print*,'ustar',ustar(1:2, 18),ustar(1:2, 20),ustar(1:2, 21)
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
!!
do ifa = 1, 3
!
do ig =1, nvfac
!
!...Basis function!
!
fstar(1, cegsp(ig, ifa), ie) = snsigml(1, cegsp(ig, ifa), ie) +&
munacl(1,1, cegsp(ig, ifa), ie)*ustar(1, ip(fagsp(ig, ifa)))+&
munacl(2,1, cegsp(ig, ifa), ie)*ustar(2, ip(fagsp(ig, ifa)))-&
munaul(1, cegsp(ig, ifa), ie)
!
fstar(2, cegsp(ig, ifa), ie) = snsigml(2, cegsp(ig, ifa), ie) +&
munacl(2,2,cegsp(ig, ifa), ie)*ustar(2, ip(fagsp(ig, ifa)))+&
munacl(1,2,cegsp(ig, ifa), ie)*ustar(1, ip(fagsp(ig, ifa)))-&
munaul(2, cegsp(ig, ifa), ie)
!
enddo
!
enddo
!
enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)

end subroutine getndvelo_lag_simpson
!
!...Face integral using gauss quadrature distribution...
!
subroutine rhsifacedg_lag_simpson(inpoel,  unkno, ustar, fstar, gelag,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:9,1:nelem),  intent(in)::fstar !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nelem),        intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(in)::gelag
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa
integer::ip1,ip2
integer:: infag(2, 3)
!...local integer array
integer,dimension(1:nvtri) :: ip
integer::fagsp(nvfac, 3), cegsp(nvfac, 3)
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::vnorm(3, 9)
real*8::xv(nvtri), yv(nvtri),b(1:3,1:nvtri)
real*8::bg(3, nvfac)
!
real*8::weigh(nvfac), posi(1,nvfac)
!...local real number
real*8::eps,c00,c05,c10,c20,c13
real*8:: shp1, shp2, wi,umid,vmid
real*8:: nx, ny ,sa
real*8:: r1, s1, r2, s2, rg, sg, r
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
!
!
if(nvfac.eq.2)then
!
fagsp(1, 1) = 1;  fagsp(2, 1) = 2;
fagsp(1, 2) = 2;  fagsp(2, 2) = 3;
fagsp(1, 3) = 3;  fagsp(2, 3) = 1;
!
!
cegsp(1, 1) = 1;  cegsp(2, 1) = 2;
cegsp(1, 2) = 3;  cegsp(2, 2) = 4;
cegsp(1, 3) = 5;  cegsp(2, 3) = 6;
!
weigh(1) = .5d0; weigh(2) = .5d0
!
elseif(nvfac.eq.3)then
!
fagsp(1, 1) = 1;  fagsp(2, 1) = 2; fagsp(3, 1) = 4;
fagsp(1, 2) = 2;  fagsp(2, 2) = 3; fagsp(3, 2) = 5;
fagsp(1, 3) = 3;  fagsp(2, 3) = 1; fagsp(3, 3) = 6;
!!
!
posi(1, 1) = -1.d0; posi(1 ,2 )= 1.d0; posi(1 ,3 )= 0.d0
!
weigh(1) = 1.d0/6.d0; weigh(2) = 1.d0/6.d0; weigh(3) = 4.d0/6.d0
!
cegsp(1, 1) = 1;  cegsp(2, 1) = 2; cegsp(3, 1) = 7;
cegsp(1, 2) = 3;  cegsp(2, 2) = 4; cegsp(3, 2) = 8;
cegsp(1, 3) = 5;  cegsp(2, 3) = 6; cegsp(3, 3) = 9;
!
endif
!
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
if(ncurv.eq.1)then
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
endif
!
do iv =1 ,nvtri
!
!print*,'iv',ivt
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo


do 550 ie = 1,nelem !...(1)ie = 1,nelem
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
ip(1:nvtri) = inpoel(1:nvtri ,ie)
!
vnorm(1:3,  1:9) = gelag(1:3, 1:9, ie);
!
do ifa = 1, 3
!
!
do ig   = 1, nvfac
!
r   = posi(1, ig)
wi  = weigh(ig)
!
rg = xv(fagsp(ig, ifa))
sg = yv(fagsp(ig, ifa))
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
gpnx = vnorm(1, cegsp(ig, ifa))
gpny = vnorm(2, cegsp(ig, ifa))
gpsa = vnorm(3, cegsp(ig, ifa))
!
!...Distribute to every corner...
!
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ip(fagsp(ig, ifa)))*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig) +&
ustar(2, ip(fagsp(ig, ifa)))*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)
!
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstar(1, cegsp(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstar(2, cegsp(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ip(fagsp(ig, ifa)))*fstar(1, cegsp(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig) +&
ustar(2, ip(fagsp(ig, ifa)))*fstar(2, cegsp(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
!if(ie==23) print*,'rhs iface idegr',ustar(1:2, ip(fagsp(ig, ifa))),gpnx,gpny,gpsa,bg(1:ndegr,ig),ip(fagsp(ig, ifa)),ig,ifa
!
enddo
!
enddo
!
rhsel(1:ndegr, 1, ie) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ie) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ie) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ie) =  elnpn(1:ndegr)
!
! if(ie==23) print*,'rhs iface',rhsel(1:3, 1, ie), ustar(1:2,ip(1)), &
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
end subroutine rhsifacedg_lag_simpson
!
!...Simpson...
!
!
!...subroutine: Calculate the F^* N dsfor all faces...
!
subroutine getfnds_lag_simpson(gflag,gelag,intfac,inpoel,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(inout)::gelag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv,ishp
!...local integer array
integer,dimension(1:nvtri) :: ipl, ipr
!...local real array
real*8,dimension(1:2, 1:2)    ::comatr !...cofactor matrix...
real*8,dimension(1:ndimn, 1:nptri)::coorp
real*8::b(3, nvtri)
real*8::vnorm(1:2)
real*8,dimension(1:nvtri):: xv, yv
real*8::posi(2, 9)
real*8::shp(nptri),dspr(nptri),dsps(nptri)
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
posi(1, 1) = 0.d0; posi(2, 1) = 0.d0;
posi(1, 2) = 1.d0; posi(2, 2) = 0.d0;
posi(1, 3) = 1.d0; posi(2, 3) = 0.d0;
posi(1, 4) = 0.d0; posi(2, 4) = 1.d0;
posi(1, 5) = 0.d0; posi(2, 5) = 1.d0;
posi(1, 6) = 0.d0; posi(2, 6) = 0.d0;
posi(1, 7) = .5d0; posi(2, 7) = 0.d0;
posi(1, 8) = .5d0; posi(2, 8) = .5d0;
posi(1, 9) = 0.d0; posi(2, 9) = .5d0;

!
do 100 ie=1, nelem !...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
!...First step: calcualte the F^* NdS for every face of all the cells...
!
!
ipl(1:nvtri) = inpoel(1:nvtri,ie)
!
!...coordinates
!
if(ncurv==0)then
coorp(1, 1:3) = coord(1, inpoel(1:3,ie))
coorp(2, 1:3) = coord(2, inpoel(1:3,ie))
!
coorp(1:2,4) = 0.5d0*(coorp(1:2,1) + coorp(1:2,2))
coorp(1:2,5) = 0.5d0*(coorp(1:2,2) + coorp(1:2,3))
coorp(1:2,6) = 0.5d0*(coorp(1:2,1) + coorp(1:2,3))
elseif(ncurv==1)then
coorp(1, 1:nptri) = coord(1,inpoel(1:nptri, ie))
coorp(2, 1:nptri) = coord(2,inpoel(1:nptri, ie))
endif
!
do ig = 1, 9!...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
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
dxdr = dxdr + dspr(ishp)*coorp(1,ishp)
dxds = dxds + dsps(ishp)*coorp(1,ishp)

dydr = dydr + dspr(ishp)*coorp(2,ishp)
dyds = dyds + dsps(ishp)*coorp(2,ishp)
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
!  print*,'Inside getfnds_lag'
!
end subroutine getfnds_lag_simpson
!
!...Calculate the velocity at the middle point...
!
subroutine getvelo_mpt(ustar,gelag,intfac,inpoel,coord,munacn,vlave,unkno)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(in)::gelag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:2,1:2,1:npoin),                   intent(in)::munacn
real*8,dimension(1:2, 1:npoin),              intent(in)::vlave
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:2, 1:npoin),              intent(inout)::ustar
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvfac) :: ipf
!
real*8::eps
real*8::unknv(1:nq, 1:nvtri)
real*8::vnorm(3, nvtri)
real*8::xv(nvtri), yv(nvtri), b(ndegr, nvtri)
!real*8,allocatable:: ucurv(:, :)
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx,dux,duy
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
!
!...For quadratic mesh, only nafac high-order nodes need be recalculated at most...
!
!allocate ucurv(2, nafac)
!
!...
!
eps = 1.d-6
!
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
if(ncurv.eq.1)then
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
endif
!
do iv =1 ,nvtri
!
!print*,'iv',ivt
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo

!
do 450 ifa = 1, nafac !...(1)ie = 1,nelem
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
!...For gauss integration...
!
vnorm(1:3, 4) = gelag(1:3, 7, iel)
vnorm(1:3, 5) = gelag(1:3, 8, iel)
vnorm(1:3, 6) = gelag(1:3, 9, iel)
!
!...For
!
!vnorm(1:3, 4) = gelag(1:3, 1, iel)
!vnorm(1:3, 5) = gelag(1:3, 2, iel)
!vnorm(1:3, 6) = gelag(1:3, 3, iel)
!
!if(ipf(3).eq.14) print*,'pt 14', lenmc, ifa
!
!if(munacn(1,1,ipf(3)).lt.1.d-1)then
!
! print*,'criterion', lenmc, munacn(ipf(3))
!
!...(ifa.le.nbfac) For boundary cells...
!
if(ifa.le.nbfac)then
!
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
unknv = 0.d0
!
do iv = 4, nvtri
!
if(inpoel(iv,iel).eq.ipf(3))then
do ideg =1, mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,iel)*b(ideg, iv)
enddo
rho  = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
presl = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
dux = vlave(1, ipf(3)) - uvtx
duy = vlave(2, ipf(3)) - vvtx
deltu = sqrt(dux**2 + duy**2)
!
mufal = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
fnx = vnorm(1, iv) !...face normal vector
fny = vnorm(2, iv)
!
ftx = fny
fty = -fnx
!
!...Ma
ustar(1, ipf(3)) = uvtx + presl/mufal*vnorm(1, iv)
ustar(2, ipf(3)) = vvtx + presl/mufal*vnorm(2, iv)
!
!...Morgan
!ustar(1, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*ftx/(mufal) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!ustar(2, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*fty/(mufal)
!
endif
enddo
!
!...(ifa.gt.nbfac) For interior cells...
!
elseif(ifa.gt.nbfac)then
!
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
unknv = 0.d0
!
do iv = 4, nvtri
if(inpoel(iv,iel).eq.ipf(3))then
do ideg =1, mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,iel)*b(ideg, iv)
enddo
rhol  = 1.d0/unknv(1, iv)
uvtxl = unknv(2, iv)
vvtxl = unknv(3, iv)
evtxl = unknv(4, iv)
presl = max(eps, (gamlg-1.d0)*rhol*(evtxl - 0.5d0*(uvtxl**2 + vvtxl**2)))
!
dux = vlave(1, ipf(3)) - uvtxl
duy = vlave(2, ipf(3)) - vvtxl
deltu = sqrt(dux**2 + duy**2)
!
mufal = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
fnx = vnorm(1, iv) !...face normal vector
fny = vnorm(2, iv)
!
ftx = fny
fty = -fnx
endif
enddo
!
rhoct = 1.d0/unkno(1, 1, ier)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ier)
vctr  = unkno(1, 3, ier)
ectr  = unkno(1, 4, ier)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center..
!
unknv = 0.d0
!
do iv = 4, nvtri
!
if(inpoel(iv,ier).eq.ipf(3))then
do ideg =1, mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ier)*b(ideg, iv)
enddo
rhor  = 1.d0/unknv(1, iv)
uvtxr = unknv(2, iv)
vvtxr = unknv(3, iv)
evtxr = unknv(4, iv)
presr = max(eps, (gamlg-1.d0)*rhor*(evtxr - 0.5d0*(uvtxr**2 + vvtxr**2)))
!
dux = vlave(1, ipf(3)) - uvtxr
duy = vlave(2, ipf(3)) - vvtxr
deltu = sqrt(dux**2 + duy**2)
!
mufar = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
endif
enddo
!...
!...Mar
ustar(1, ipf(3)) = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fnx
ustar(2, ipf(3)) = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fny
!
!...Morgan
!    ustar(1, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*ftx/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!    ustar(2, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*fty/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(2, 1, iv)
!
!    if(ipf(3).eq.190) print*,'midlle velocity',ustar(1:2, ipf(3)),ipf(3),mufal+mufar,tfx,tfy
endif
!
!endif
!
450 enddo
!
end subroutine getvelo_mpt
!
subroutine getndvelo_lag_curv_vilar(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fstar)
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
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:4,1:nvtri, 1:nelem),  intent(out)::fstar !...Riemann forces
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin
integer::iel,ier
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:3)     :: ipf
integer::indnd(npoin)

!...local real array
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
real*8::munaci(2, 2)
real*8::muutg(1:2, 1:2, 1:npoin), mutag(npoin)
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:nvtri)::murie
real*8::vnorm(1:3, 1:4, 1:nvtri)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8::aujmp(1:2, 1:nvtri)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:nvtri):: xv, yv
real*8::usnom(2), ustng(2)
real*8::xpf(1:2,1:3),xcf(1:2)
!...local real number
real*8::unorm, tfx, tfy, nfx, nfy,signm, utang
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dwav1,dwav2
real*8::anx, any, len1, len2, tfx1, tfy1
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufac,mufal,mufar,rhsu1,rhsu2,detma
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: cnsup(:), munacn(:,:,:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: usold(:,:), munacu(:,:), snsigm(:,:)
real*8,allocatable:: munacl(:,:,:,:,:)
real*8,allocatable:: snsigml(:,:,:,:), munaul(:,:,:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:2,1:2,1:npoin), cnsup(1:npoin))
allocate (usold(1:ndimn, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munacl(1:2,1:2,1:4, 1:nvtri, 1:nelem), munaul(1:ndimn, 1:4, 1:nvtri,  1:nelem),&
snsigml(1:ndimn, 1:4,  1:nvtri,  1:nelem))
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
!...Mark the boundary nodes for case 2...
!
if(ncase.eq.2)then
do ifa = 1, nbfac
ipf(1:3) = intfac(3:5, ifa)
indnd(ipf(1:3)) = 1
enddo
endif
!
!...shape functions
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0;  yv(1) = 0.d0
xv(2) = 1.d0;  yv(2) = 0.d0
xv(3) = 0.d0;  yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...Part I: Get the averaged nodal velocity...
!
do 200 ie = 1,nelem !...(1)ie = 1,nelem
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
! if(ip(iv)==36) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
enddo
!
!...Accumulate nodal velocity...
!
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
!...Part II:Get the related coefficients for Riemann-like
!
!
!...Zero out munacn, munacu, snsigm...
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
!...Deal with the high-order nodes...
!
muutg = 0.d0
mutag = 0.d0
!
vlave = ustar
!
do 250 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!
!...Give the normal vector of every face...
!
vnorm(1:3, 1, 1) = gelag(1:3, 3, ie); vnorm(1:3, 2, 1) = gelag(1:3, 1, ie) !...For point ip(1)
vnorm(1:3, 3, 1) = gelag(1:3, 9, ie); vnorm(1:3, 4, 1) = gelag(1:3, 4, ie) !...For point ip(1)
!
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 5, ie); vnorm(1:3, 4, 2) = gelag(1:3, 6, ie) !...For point ip(2)
!
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 7, ie); vnorm(1:3, 4, 3) = gelag(1:3, 8, ie) !...For point ip(3)
!
vnorm(1:3, 1, 4) = gelag(1:3, 4, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(4)
!
vnorm(1:3, 1, 5) = gelag(1:3, 6, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(5)
!
vnorm(1:3, 1, 6) = gelag(1:3, 8, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(6)
!
!...Loop 1: Calculate the tensor stress and a_c at vertices...
!
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
!...Stress tensor...
!
rho  = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
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
!
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
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
!
!
!...Loop 2: Get the impedence at the vertices...
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
murie(iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ip(iv).eq.5) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ip(iv)),unknv(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, 3 !...Low-order nodes...
do ifa = 1, 4 !...Every corner consists of 2 faces...
!
!
!
call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), unknv(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
if(ifa.le.2)then
!

!
elseif(ifa.ge.3)then
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
munacn(1:2, 1, ip(iv)) = munacn(1:2, 1, ip(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ip(iv)) = munacn(1:2, 2, ip(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ip(iv)) = munacu(1:2, ip(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ip(iv)) = snsigm(1:2, ip(iv)) + snsigm_rie(1:2)!
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
endif
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Local variable...
!
munacl(1:2, 1, ifa, iv, ie) =  munacn_rie(1:2, 1)
munacl(1:2, 2, ifa, iv, ie) =  munacn_rie(1:2, 2)

!
munaul(1:2, ifa, iv, ie) =  munacu_rie(1:2)
!
snsigml(1:2, ifa, iv, ie)=  snsigm_rie(1:2)
!
enddo
!
enddo
!
!...Get the summed denominator cooefficients sum(mu*n*a_c) for high-order nodes...
!
do iv  = 4, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.281) print*,'p19 muacn(281) pre++', munacn(281),ie,iv,ifa
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!

call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), unknv(2:3, iv), sigma(1:2, 1:2, iv), &
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
munacn(1:2, 1, ip(iv)) = munacn(1:2, 1, ip(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ip(iv)) = munacn(1:2, 2, ip(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ip(iv)) = munacu(1:2, ip(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ip(iv)) = snsigm(1:2, ip(iv)) + snsigm_rie(1:2)!
!
!   if(ip(iv).eq.14) print*,'p281 muacn(vv)',ie,murie(iv),munacn(1:2,1:2,ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!   munacl(1, iv, ie) = munacl(1, iv, ie) + murie(iv)*vnorm(3, ifa, iv)* &
!                    abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!   munacl(1, iv, ie) =  murie(iv)*vnorm(3, ifa, iv)* &
!                        abs(vnorm(1, ifa, iv)*aujmp(1, iv) + vnorm(2, ifa, iv)*aujmp(2, iv))
!
! if(ip(iv).eq.281) print*,'p19 muacn(28) prep---',murie(iv), munacn(281), munacu(1:2,ip(iv)),vnorm(3, ifa, iv),&
!                                                vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknv(2:3,iv),&
!                                                 vlave(1:2, ip(iv))

! if(ip(iv).eq.15) print*,'p19 muacn(28) postxxxx',murie(iv), munacu(1:2,ip(iv)),ie, ifa,iv!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
munacl(1:2, 1, ifa, iv, ie) =  munacn_rie(1:2, 1)
munacl(1:2, 2, ifa, iv, ie) =  munacn_rie(1:2, 2)
!
munaul(1:2, ifa, iv, ie) =  munacu_rie(1:2)
!
snsigml(1:2, ifa, iv, ie)=  snsigm_rie(1:2)
!
enddo
!
!...Coefficients for middle points...
!
!  muutg(1:2, ip(iv)) = muutg(1:2, ip(iv)) + murie(iv)*unknv(2:3, iv)
!  mutag(ip(iv)) = mutag(ip(iv)) + murie(iv)
!

enddo
!
250 enddo  !...(1)ie = 1,nelem
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
!...Zero out the bnorm
!
!
!...
!...Zero out boundary pressure...
!
!
!...Fourth part: Solve the nodal velocity...
!
!...deactivate the pseudo pressure....
!
fpres = 0.d0
!
!...4.1: Update the Riemann forces at every node...
!
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
!if(idcrd(ipoin).eq.0)then
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
!elseif(idcrd(ipoin).eq.1)then
!  if(muacn(ipoin).gt.1.d-9)then
!   ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
!   ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
!  endif
!endif
!
! if(ipoin.eq.70) print*,ustar(1:2,ipoin),munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
endif
enddo
!
!...Calculate the midlle point velocity....
!
call getvelo_mpt_marie(ustar,gelag,intfac,inpoel,coord,unkno,indnd)
!
!....Bd velocity
!
! print*,'ustar--',ustar(1:2, 14),munacu(1:2,14) ,snsigm(1:2,14), munacn(1:2,1:2,14)
!
!
if(ncase.eq.1)then
!
do 900 ifa = 1 , nbfac
!
ipf(1:3) = intfac(3:5, ifa)
!
ustar(1, ipf(1)) = sin(pi*coord(1,ipf(1)))*cos(pi*coord(2,ipf(1)))
ustar(2, ipf(1)) =-cos(pi*coord(1,ipf(1)))*sin(pi*coord(2,ipf(1)))
!
ustar(1, ipf(2)) = sin(pi*coord(1,ipf(2)))*cos(pi*coord(2,ipf(2)))
ustar(2, ipf(2)) =-cos(pi*coord(1,ipf(2)))*sin(pi*coord(2,ipf(2)))
!
ustar(1, ipf(3)) = sin(pi*coord(1,ipf(3)))*cos(pi*coord(2,ipf(3)))
ustar(2, ipf(3)) =-cos(pi*coord(1,ipf(3)))*sin(pi*coord(2,ipf(3)))
!
if(coord(1, ipf(1)).lt.1.d-6.or.abs(coord(1, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf',ipf(1)
ustar(1, ipf(1)) = 0.d0
endif
if(coord(2, ipf(1)).lt.1.d-6.or.abs(coord(2, ipf(1))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ipf(1)
ustar(2, ipf(1)) = 0.d0
endif
if(coord(1, ipf(3)).lt.1.d-6.or.abs(coord(1, ipf(3))-1.d0).lt.1.d-6) then
!    print*,'ipf',ipf(1)
ustar(1, ipf(3)) = 0.d0
endif
if(coord(2, ipf(3)).lt.1.d-6.or.abs(coord(2, ipf(3))-1.d0).lt.1.d-6) then
!    print*,'ipf2',ipf(1)
ustar(2, ipf(3)) = 0.d0
endif
900 enddo
!
!ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
!
!   print*,'ustar',ustar(1:2, 14)!,ustar(1:2, 20),ustar(1:2, 21)
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
unknv = 0.d0
!
!...Low order nodes...
!
do iv = 1, 3
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
do ifa = 1, 4
!
fstar(1, ifa, iv, ie) = snsigml(1, ifa, iv, ie) + &
munacl(1, 1, ifa, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv))+&
munacl(2, 1, ifa, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(1, 1, iv, ie)
fstar(2, ifa, iv, ie) = snsigml(2, ifa, iv, ie) + &
munacl(1,2,ifa, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv))+&
munacl(2, 2, ifa, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(1, 1, iv, ie)
enddo
!
enddo
!
!...High-order nodes
!
do iv = 4, nvtri
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
do ifa = 1, 2
!
fstar(1, ifa, iv, ie) = snsigml(1, ifa, iv, ie) + &
munacl(1, 1, ifa, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv))+&
munacl(2, 1, ifa, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(1, 1, iv, ie)
fstar(2, ifa, iv, ie) = snsigml(2, ifa, iv, ie) + &
munacl(2,1,ifa, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv))+&
munacl(2, 2, ifa, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(1, 1, iv, ie)
enddo
!
enddo
enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
end subroutine getndvelo_lag_curv_vilar
!
!...Calculate the velocity at the middle point for Marie method...
!
subroutine getvelo_mpt_marie(ustar,gelag,intfac,inpoel,coord,unkno,indnd)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:3,1:ngelg,1:nelem+nbfac), intent(in)::gelag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:2, 1:npoin),              intent(inout)::ustar
integer*4,dimension(1:npoin),                intent(in)::indnd
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvfac) :: ipf
!
real*8::eps
real*8::unknv(1:nq, 1:nvtri)
real*8::vnorm(3, nvtri)
real*8::xv(nvtri), yv(nvtri), b(ndegr, nvtri)
real*8::xpf(1:2, 1:nvfac)
!real*8,allocatable:: ucurv(:, :)
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx,dux,duy
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
real*8::othog
!
!...For quadratic mesh, only nafac high-order nodes need be recalculated at most...
!
eps = 1.d-6
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
if(ncurv.eq.1)then
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
endif
!
do iv =1 ,nvtri
!
!print*,'iv',ivt
!
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo

!
do 450 ifa = 1, nafac !...(1)ie = 1,nelem
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
!...For gauss integration...
!
!vnorm(1:3, 4) = gelag(1:3, 7, iel)
!vnorm(1:3, 5) = gelag(1:3, 8, iel)
!vnorm(1:3, 6) = gelag(1:3, 9, iel)
!
!...For
!
vnorm(1:3, 4) = gelag(1:3, 1, iel)
vnorm(1:3, 5) = gelag(1:3, 2, iel)
vnorm(1:3, 6) = gelag(1:3, 3, iel)
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
!if(ipf(3).eq.14) print*,'pt 14', lenmc, ifa
!
if(othog.lt.1.d-8)then
!
! print*,'criterion', lenmc, munacn(ipf(3))
!
!...(ifa.le.nbfac) For boundary cells...
!
if(ifa.le.nbfac)then
!
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
unknv = 0.d0
!
do iv = 4, nvtri
!
if(inpoel(iv,iel).eq.ipf(3))then
do ideg =1, mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,iel)*b(ideg, iv)
enddo
rho  = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
presl = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
mufal = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
fnx = vnorm(1, iv) !...face normal vector
fny = vnorm(2, iv)
!
ftx = fny
fty = -fnx
!
!...Ma
!
if(indnd(ipf(3)).eq.0)then
ustar(1, ipf(3)) = uvtx + presl/mufal*vnorm(1, iv)
ustar(2, ipf(3)) = vvtx + presl/mufal*vnorm(2, iv)
endif
!
!...Morgan
!ustar(1, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*ftx/(mufal) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!ustar(2, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*fty/(mufal)
!
endif
enddo
!
!...(ifa.gt.nbfac) For interior cells...
!
elseif(ifa.gt.nbfac)then
!
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
unknv = 0.d0
!
do iv = 4, nvtri
if(inpoel(iv,iel).eq.ipf(3))then
do ideg =1, mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,iel)*b(ideg, iv)
enddo
rhol  = 1.d0/unknv(1, iv)
uvtxl = unknv(2, iv)
vvtxl = unknv(3, iv)
evtxl = unknv(4, iv)
presl = max(eps, (gamlg-1.d0)*rhol*(evtxl - 0.5d0*(uvtxl**2 + vvtxl**2)))
!
!
mufal = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
fnx = vnorm(1, iv) !...face normal vector
fny = vnorm(2, iv)
!
ftx = fny
fty = -fnx
endif
enddo
!
rhoct = 1.d0/unkno(1, 1, ier)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ier)
vctr  = unkno(1, 3, ier)
ectr  = unkno(1, 4, ier)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center..
!
unknv = 0.d0
!
do iv = 4, nvtri
!
if(inpoel(iv,ier).eq.ipf(3))then
do ideg =1, mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ier)*b(ideg, iv)
enddo
rhor  = 1.d0/unknv(1, iv)
uvtxr = unknv(2, iv)
vvtxr = unknv(3, iv)
evtxr = unknv(4, iv)
presr = max(eps, (gamlg-1.d0)*rhor*(evtxr - 0.5d0*(uvtxr**2 + vvtxr**2)))
!
!
mufar = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
endif
enddo
!...
!...Mar
!
if(indnd(ipf(3)).eq.0)then
ustar(1, ipf(3)) = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fnx
ustar(2, ipf(3)) = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fny
endif
!
!...Morgan
!    ustar(1, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*ftx/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!    ustar(2, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*fty/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(2, 1, iv)
!
!    if(ipf(3).eq.190) print*,'midlle velocity',ustar(1:2, ipf(3)),ipf(3),mufal+mufar,tfx,tfy
endif
!
endif
!
450 enddo
!
end subroutine getvelo_mpt_marie
!
!...Calculate the velocity at the middle point for Marie method...
!
subroutine getvelo_mpt_mariequad(ustar,gelagq,intfac,ipqua,coord,unkno,indnd)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:3, 1:ngelgq, 1:nquad),    intent(in)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:2, 1:npoin),              intent(inout)::ustar
integer*4,dimension(1:npoin),                intent(in)::indnd
!...Local integer
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg
integer::iv
!...local integer array
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!
real*8::eps
real*8::unknv(1:nq, 1:nvqua)
real*8::vnorm(3, nvqua)
real*8::xvq(nvqua), yvq(nvqua), b(ndegr, nvqua)
real*8::xpf(1:2, 1:nvfac)
!real*8,allocatable:: ucurv(:, :)
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar
real*8::deltu
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx,dux,duy
real*8::fnx,fny, ftx, fty, rho
real*8::dr, ds,rc,sc
real*8::othog
!
!...For quadratic mesh, only nafac high-order nodes need be recalculated at most...
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
!
xvq(9) =  0.d0; yvq(9) =  0.d0
!
do iv =1 ,nvqua
!
!print*,'iv',ivt
!
b(1, iv) = 1.d0
b(2, iv) = (xvq(iv)-rc)/dr
b(3, iv) = (yvq(iv)-sc)/ds
enddo

!
do 450 ifa = 1, nafac !...(1)ie = 1,nelem
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
!...For gauss integration...
!
!vnorm(1:3, 4) = gelag(1:3, 7, iel)
!vnorm(1:3, 5) = gelag(1:3, 8, iel)
!vnorm(1:3, 6) = gelag(1:3, 9, iel)
!
!...For
!
vnorm(1:3, 5) = gelagq(1:3, 1, iel)
vnorm(1:3, 6) = gelagq(1:3, 2, iel)
vnorm(1:3, 7) = gelagq(1:3, 3, iel)
vnorm(1:3, 8) = gelagq(1:3, 4, iel)
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
!if(ipf(3).eq.14) print*,'pt 14', lenmc, ifa
!
if(othog.lt.1.d-8)then
!
! print*,'criterion', lenmc, munacn(ipf(3))
!
!...(ifa.le.nbfac) For boundary cells...
!
if(ifa.le.nbfac)then
!
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
unknv = 0.d0
!
do iv = 5, 8
!
if(ipqua(iv,iel).eq.ipf(3))then
do ideg =1, mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,iel)*b(ideg, iv)
enddo
rho  = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
presl = max(eps, (gamlg-1.d0)*rho*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))
!
mufal = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
fnx = vnorm(1, iv) !...face normal vector
fny = vnorm(2, iv)
!
ftx = fny
fty = -fnx
!
!...Ma
!
if(indnd(ipf(3)).eq.0)then
ustar(1, ipf(3)) = uvtx + presl/mufal*vnorm(1, iv)
ustar(2, ipf(3)) = vvtx + presl/mufal*vnorm(2, iv)
endif
!
!...Morgan
!ustar(1, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*ftx/(mufal) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!ustar(2, ipf(3)) = (mufal*uvtx*ftx + mufal*vvtx*fty)*fty/(mufal)
!
endif
enddo
!
!...(ifa.gt.nbfac) For interior cells...
!
elseif(ifa.gt.nbfac)then
!
rhoct = 1.d0/unkno(1, 1, iel)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, iel)
vctr  = unkno(1, 3, iel)
ectr  = unkno(1, 4, iel)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
unknv = 0.d0
!
do iv = 5, 8
if(ipqua(iv,iel).eq.ipf(3))then
do ideg =1, mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,iel)*b(ideg, iv)
enddo
rhol  = 1.d0/unknv(1, iv)
uvtxl = unknv(2, iv)
vvtxl = unknv(3, iv)
evtxl = unknv(4, iv)
presl = max(eps, (gamlg-1.d0)*rhol*(evtxl - 0.5d0*(uvtxl**2 + vvtxl**2)))
!
!
mufal = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
fnx = vnorm(1, iv) !...face normal vector
fny = vnorm(2, iv)
!
ftx = fny
fty = -fnx
endif
enddo
!
rhoct = 1.d0/unkno(1, 1, ier)         !...ct denots center of one cell; cn denotes corner of one cell.
uctr  = unkno(1, 2, ier)
vctr  = unkno(1, 3, ier)
ectr  = unkno(1, 4, ier)
pctr  = max(eps, (gamlg-1.d0)*rhoct*(ectr - 0.5d0*(uctr**2 + vctr**2)))
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center..
!
unknv = 0.d0
!
do iv = 5, 8
!
if(ipqua(iv,ier).eq.ipf(3))then
do ideg =1, mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ier)*b(ideg, iv)
enddo
rhor  = 1.d0/unknv(1, iv)
uvtxr = unknv(2, iv)
vvtxr = unknv(3, iv)
evtxr = unknv(4, iv)
presr = max(eps, (gamlg-1.d0)*rhor*(evtxr - 0.5d0*(uvtxr**2 + vvtxr**2)))
!
!
mufar = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
endif
enddo
!...
!...Mar
!
if(indnd(ipf(3)).eq.0)then
ustar(1, ipf(3)) = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fnx
ustar(2, ipf(3)) = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fny
endif
!
!...Morgan
!    ustar(1, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*ftx/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(1, 1, iv)
!    ustar(2, ipf(3)) = ((mufal*uvtxl + mufar*uvtxr)*ftx + (mufal*vvtxl + mufar*vvtxr)*fty)*fty/(mufal+mufar) !- (presr- presl)/(mufal+mufar)*vnorm(2, 1, iv)
!
!    if(ipf(3).eq.190) print*,'midlle velocity',ustar(1:2, ipf(3)),ipf(3),mufal+mufar,tfx,tfy
endif
!
endif
!
450 enddo
!
end subroutine getvelo_mpt_mariequad

!
!...Subroutine to deal with BC with prescribed pressure and normal velocity
!
subroutine getbcfc_lag(bface, intfac, gflag,  fpres, coord, ustar,  itime)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),           intent(out)::ustar
real*8,dimension(1:ndimn,1:npoin),          intent(out)::fpres
integer:: itime
!
!...Local variable...
!
real*8:: pres0, bnxpr, bnypr, bnx, bny
real*8:: presi, prese
real*8:: dtime, htkid, taokid
real*8:: coorp(2, nvfac) !nvfac: no of vertices of one face
integer::ifa
integer::ipf(nvfac)
!
!...Zero our fpres and pres0...
!
fpres = 0.d0
!
pres0 = 0.d0
!
if(ncase.eq.3)then
!
pres0 = 1.d-6  !...Noh test problem
!
!  print*,'gflag',gflag(1:3, 12)
!
do 1300 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
if(bface(3, ifa).eq.21)then
!
bnxpr =  gflag(1, ifa)*gflag(3,ifa)*0.5d0
bnypr =  gflag(2, ifa)*gflag(3,ifa)*0.5d0
!
fpres(1, ipf(1)) = fpres(1, ipf(1)) + bnxpr*pres0
fpres(2, ipf(1)) = fpres(2, ipf(1)) + bnypr*pres0
!
fpres(1, ipf(2)) = fpres(1, ipf(2)) + bnxpr*pres0
fpres(2, ipf(2)) = fpres(2, ipf(2)) + bnypr*pres0
endif
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1300 enddo
!
elseif(ncase.eq.4)then
!
presi = 0.1d0
prese = 10.d0
!
!  if(rkstg.eq.1)then
dtime = (itime)*dtfix
!  elseif(rkstg.eq.2)then
!   dtime = (itime+1)*dtfix
!  endif
!
htkid = (sqrt(1.d0-(dtime)**2))**(-4)
!
!  print*,'gflag',itime,dtime,itime,dtfix
!
do 1350 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
!print*,'gflag',ifa,ipf(1:2),bface(1:2, ifa), bface(3, ifa)

!
if(bface(3, ifa).eq.21)then !...inner circle
!
bnxpr =  gflag(1, ifa)*gflag(3,ifa)*0.5d0
bnypr =  gflag(2, ifa)*gflag(3,ifa)*0.5d0
!
fpres(1, ipf(1)) = fpres(1, ipf(1)) + bnxpr*presi*htkid
fpres(2, ipf(1)) = fpres(2, ipf(1)) + bnypr*presi*htkid
!
fpres(1, ipf(2)) = fpres(1, ipf(2)) + bnxpr*presi*htkid
fpres(2, ipf(2)) = fpres(2, ipf(2)) + bnypr*presi*htkid
!
!print*,'gflag',ifa,ipf(1), fpres(1:2, ipf(1)),gflag(1:3, ifa), htkid, bnxpr*presi*htkid,bnypr*presi*htkid
!print*,'gflagipf2',ifa,ipf(2), fpres(1:2, ipf(2)),gflag(1:3, ifa), bnxpr*presi*htkid,bnypr*presi*htkid
!
elseif(bface(3, ifa).eq.23)then !...outer circle
!
bnxpr =  gflag(1, ifa)*gflag(3,ifa)*0.5d0
bnypr =  gflag(2, ifa)*gflag(3,ifa)*0.5d0
!
fpres(1, ipf(1)) = fpres(1, ipf(1)) + bnxpr*prese*htkid
fpres(2, ipf(1)) = fpres(2, ipf(1)) + bnypr*prese*htkid
!
fpres(1, ipf(2)) = fpres(1, ipf(2)) + bnxpr*prese*htkid
fpres(2, ipf(2)) = fpres(2, ipf(2)) + bnypr*prese*htkid
endif
!
!print*,'ifa',fpres(1:2, 5)
!
1350 enddo
!
elseif(ncase.eq.5)then !...kidder ball
!
!if(rkstg.eq.1)then
   dtime = (itime)*dtfix
!elseif(rkstg.eq.2)then
!   dtime = (itime+1)*dtfix
!endif
!
do 1450 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
!   print*,'gflag',ifa,ipf(1:2),gflag(1:2, ifa),bface(3, ifa)
!
if(bface(3, ifa).eq.25)then !...inner circle
!
coorp(1, 1:nvfac) = coord(1, ipf(1:nvfac))
coorp(2, 1:nvfac) = coord(2, ipf(1:nvfac))
!
ustar(1, ipf(1:nvfac)) = coorp(1, 1:nvfac)*(dtime-1.d0)/(1.d0+(dtime-1.d0)**2)
ustar(2, ipf(1:nvfac)) = coorp(2, 1:nvfac)*(dtime-1.d0)/(1.d0+(dtime-1.d0)**2)
!
!
endif
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1450 enddo
!
elseif(ncase.eq.7)then !...sedov
!
do 1550 ifa = 1 , nbfac
!
ipf(1:nvfac) = intfac(3:2+nvfac, ifa)
!
!   print*,'gflag',ifa,ipf(1:2),gflag(1:2, ifa),bface(3, ifa)
!
if(bface(3, ifa).eq.25)then !...inner circle
!
ustar(1, ipf(1:nvfac)) = 0.d0
ustar(2, ipf(1:nvfac)) = 0.d0
!
endif
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1550 enddo

!
endif
!
end subroutine getbcfc_lag
!
!...Subroutine to deal with BC with prescribed pressure and normal velocity
!
subroutine getbc_lagmaire(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),          intent(out)::fpres
real*8,intent(in)::munacn(1:2, 1:2, 1:npoin)
real*8,intent(in)::snsigm(1:ndimn, 1:npoin)
real*8,intent(in)::munacu(1:ndimn, 1:npoin)
!
!...Local variable...
!
real*8:: pres0, bnxpr, bnypr, bnx, bny
real*8:: presi, prese
real*8:: dtime, htkid, taokid
real*8:: coorp(2, nvfac) !nvfac: no of vertices of one face
integer::ifa, ip
integer::ipf(nvfac)
real*8::bnpre(2, npoin)
real*8::munaci(2,2)
real*8::lnx,lny,coefl,coefr,rhsx,rhsy,detma
real*8:: rcoef(2)
real*8:: xph(2, 2)
!
!...Zero our fpres and pres0...
!
fpres = 0.d0
pres0 = 0.d0
!
!if(ncase.eq.1)then
!
bnpre = 0.d0
!
!  print*,'gflag',gflag(1:3, 12)
!
do 1300 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
xph(2, 1:2) = coord(2, ipf(1:2))
!
rcoef(1) = 1.d0 - alfrz + alfrz*xph(2, 1)
rcoef(2) = 1.d0 - alfrz + alfrz*xph(2, 2)
!
if(bface(4, ifa).eq.221)then
bnxpr =  0.d0
bnypr =  -1.d0*gflag(3,ifa)*0.5d0
elseif(bface(4, ifa).eq.222)then
bnxpr =  -1.d0*gflag(3,ifa)*0.5d0
bnypr =  0.d0
endif
!
bnpre(1, ipf(1)) = bnpre(1, ipf(1)) + bnxpr!*(2.d0*rcoef(1) + rcoef(2))/3.d0
bnpre(2, ipf(1)) = bnpre(2, ipf(1)) + bnypr!*(2.d0*rcoef(1) + rcoef(2))/3.d0

bnpre(1, ipf(2)) = bnpre(1, ipf(2)) + bnxpr!*(2.d0*rcoef(2) + rcoef(1))/3.d0
bnpre(2, ipf(2)) = bnpre(2, ipf(2)) + bnypr!*(2.d0*rcoef(2) + rcoef(1))/3.d0
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1300 enddo
!
do 1350 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
!
do ip = 1, 2
detma = munacn(1, 1, ipf(ip))*munacn(2, 2, ipf(ip)) - munacn(2, 1, ipf(ip))*munacn(1, 2, ipf(ip))
munaci(1, 1) = munacn(2, 2, ipf(ip))/detma
munaci(1, 2) =-munacn(1, 2, ipf(ip))/detma
munaci(2, 1) =-munacn(2, 1, ipf(ip))/detma
munaci(2, 2) = munacn(1, 1, ipf(ip))/detma
!
!print*,'ipf',ifa, detma, ip
!
lnx = bnpre(1, ipf(ip))
lny = bnpre(2, ipf(ip))
!
coefl = (munaci(1,1)*lnx + munaci(1,2)*lny)*lnx + (munaci(2,1)*lnx + munaci(2,2)*lny)*lny
!
rhsx = munacu(1, ipf(ip))  - snsigm(1, ipf(ip))
rhsy = munacu(2, ipf(ip))  - snsigm(2, ipf(ip))
coefr = (munaci(1,1)*rhsx + munaci(1,2)*rhsy)*lnx + (munaci(2,1)*rhsx + munaci(2,2)*rhsy)*lny
!
fpres(1, ipf(ip)) = coefr/max(1.d-10,coefl)*lnx
fpres(2, ipf(ip)) = coefr/max(1.d-10,coefl)*lny
enddo
!
!print*,'ifa',ipf(1:2),ifa,fpres(1,:)
!
1350 enddo
!
!endif
!
end subroutine getbc_lagmaire
!
!...Subroutine to deal with BC with prescribed pressure and normal velocity prescribed(Maire)
!
subroutine getbc_lagmaire2(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),          intent(out)::fpres
real*8,intent(in)::munacn(1:2, 1:2, 1:npoin)
real*8,intent(in)::snsigm(1:ndimn, 1:npoin)
real*8,intent(in)::munacu(1:ndimn, 1:npoin)
integer, intent(in):: itime
!
!...Local variable...
!
real*8:: pres0, bnxpr, bnypr, bnx, bny
real*8:: presi, prese
real*8:: dtime, htkid, taokid
real*8:: coorp(2, nvfac) !nvfac: no of vertices of one face
integer::ifa, ip
integer::ipf(nvfac)
real*8:: bsnv(2, npoin), bsnpf(2, npoin)
real*8::munaci(2,2)
real*8::lnx,lny,coefl,coefr,rhsx,rhsy,detma
real*8:: bsnx, bsny, prex, prin
real*8:: rcoef(2)
real*8:: xph(2, 2)
!
!...Zero our fpres and pres0...
!
fpres = 0.d0
pres0 = 1.d-6
!
bsnv = 0.d0  !...area*n for velocity boundary
bsnpf = 0.d0 !...force for pressure boundary


!
!...Part I: The pressure force for any boundary node...
!
do 1300 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
xph(1, 1:2) = coord(1, ipf(1:2))
xph(2, 1:2) = coord(2, ipf(1:2))
!
rcoef(1) = 1.d0 - alfrz + alfrz*xph(2, 1)
rcoef(2) = 1.d0 - alfrz + alfrz*xph(2, 2)
!
bsnx = 0.d0
bsny = 0.d0
!
if(bface(3, ifa).eq.22)then !...Normal velocity boundary
!
if(bface(4, ifa).eq.221)then
!bsnx = gflag(1,ifa)*gflag(3,ifa)*0.5d0
!bsny=  gflag(2,ifa)*gflag(3,ifa)*0.5d0
bsnx = sqrt(2.d0)/2.d0*gflag(3,ifa)*0.5d0
bsny =-sqrt(2.d0)/2.d0*gflag(3,ifa)*0.5d0
!
bsnx = 0.d0*gflag(3,ifa)*0.5d0
bsny =-1.d0*gflag(3,ifa)*0.5d0
!
if(ncase.eq.1.or.ncase.eq.2.or.ncase.eq.12)then
if(xph(2,1).lt.1.d-6.and.xph(2,2).lt.1.d-6)then
bsnx = 0.d0*gflag(3,ifa)*0.5d0
bsny =-1.d0*gflag(3,ifa)*0.5d0
elseif(abs(xph(2,1)-1.d0).lt.1.d-6.and.abs(xph(2,2)-1.d0).lt.1.d-6)then
bsnx = 0.d0*gflag(3,ifa)*0.5d0
bsny = 1.d0*gflag(3,ifa)*0.5d0
endif
endif
!
!bsnx = sin(3.d0/8.d0*pi)*gflag(3,ifa)*0.5d0
!bsny =-cos(3.d0/8.d0*pi)*gflag(3,ifa)*0.5d0
elseif(bface(4, ifa).eq.222)then
!bsnx = gflag(1,ifa)*gflag(3,ifa)*0.5d0
!bsny=  gflag(2,ifa)*gflag(3,ifa)*0.5d0
bsnx = -1.0*gflag(3,ifa)*0.5d0
bsny=   0.d0*gflag(3,ifa)*0.5d0

if(ncase.eq.1.)then
if(xph(1,1).lt.1.d-6.and.xph(1,2).lt.1.d-6)then
bsnx = -1.d0*gflag(3,ifa)*0.5d0
bsny =  0.d0*gflag(3,ifa)*0.5d0
elseif(abs(xph(1,1)-1.d0).lt.1.d-6.and.abs(xph(1,2)-1.d0).lt.1.d-6)then
bsnx =  1.d0*gflag(3,ifa)*0.5d0
bsny =  0.d0*gflag(3,ifa)*0.5d0
endif
endif

endif
!
bsnv(1, ipf(1)) = bsnv(1, ipf(1)) + bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0
bsnv(2, ipf(1)) = bsnv(2, ipf(1)) + bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0

bsnv(1, ipf(2)) = bsnv(1, ipf(2)) + bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0
bsnv(2, ipf(2)) = bsnv(2, ipf(2)) + bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0
!
elseif(bface(3, ifa).eq.21)then !...Pressure described
!
bsnx = gflag(1,ifa)*gflag(3,ifa)*0.5d0
bsny=  gflag(2,ifa)*gflag(3,ifa)*0.5d0
!
if(ncase.eq.3)then
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0
!
elseif(ncase.eq.4)then !...Kidder shell
!if(rkstg.eq.1)then
dtime = (itime-1.d0)*dtfix
!elseif(rkstg.eq.2)then
!dtime = (itime)*dtfix
!endif
prin = 0.1d0
prex = 10.d0
!
pres0 = prin*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0

elseif(ncase.eq.14)then !...Coggshell

dtime = (itime-1.d0)*dtfix
pres0 = (gamlg-1.d0)*(1.d0-dtime)**(-2.25d0)*(3.d0/8.d0/(1.d0-dtime))**2

bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0*coord(1, ipf(1))**2
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0*coord(1, ipf(1))**2

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0*coord(1, ipf(2))**2
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0*coord(1, ipf(2))**2

endif
!
elseif(bface(3, ifa).eq.23)then
!
bsnx =  gflag(1,ifa)*gflag(3,ifa)*0.5d0
bsny =  gflag(2,ifa)*gflag(3,ifa)*0.5d0
!
if(ncase.eq.4)then !...Kidder shell
!
!if(rkstg.eq.1)then
dtime = (itime-1.d0)*dtfix
!elseif(rkstg.eq.2)then
!dtime = (itime)*dtfix
!endif

prin = 0.1d0
prex = 10.d0
!
pres0 = prex*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0
!
endif
!
endif
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1300 enddo
!
!...Part II: Get the pressure from the normal velocity prescription
!
do 1350 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
if(bface(3, ifa).eq.22)then
!
do ip = 1, 2
!
detma = munacn(1, 1, ipf(ip))*munacn(2, 2, ipf(ip)) - munacn(2, 1, ipf(ip))*munacn(1, 2, ipf(ip))
munaci(1, 1) = munacn(2, 2, ipf(ip))/detma
munaci(1, 2) =-munacn(1, 2, ipf(ip))/detma
munaci(2, 1) =-munacn(2, 1, ipf(ip))/detma
munaci(2, 2) = munacn(1, 1, ipf(ip))/detma
!
!print*,'ipf',ifa, detma, ip
!
lnx = bsnv(1, ipf(ip))
lny = bsnv(2, ipf(ip))
!
coefl = (munaci(1,1)*lnx + munaci(1,2)*lny)*lnx + (munaci(2,1)*lnx + munaci(2,2)*lny)*lny
!
rhsx = munacu(1, ipf(ip))  - snsigm(1, ipf(ip)) - bsnpf(1, ipf(ip))
rhsy = munacu(2, ipf(ip))  - snsigm(2, ipf(ip)) - bsnpf(2, ipf(ip))
!
coefr = (munaci(1,1)*rhsx + munaci(1,2)*rhsy)*lnx + (munaci(2,1)*rhsx + munaci(2,2)*rhsy)*lny
!
fpres(1, ipf(ip)) = coefr/max(1.d-16,coefl)*lnx
fpres(2, ipf(ip)) = coefr/max(1.d-16,coefl)*lny
enddo
!
endif
!
!print*,'ifa',ipf(1:2),ifa,fpres(1,:)
!
1350 enddo
!
!...Part III: Get the whole force...
!
fpres(:, :) = fpres(:, :) + bsnpf(:, :)
!
!print*,'ifa',fpres(1:2,500),bsnpf(1:2, 500)
!
end subroutine getbc_lagmaire2
!
!...Subroutine to deal with BC with prescribed pressure and normal velocity(Maire)
!
subroutine getbc_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),          intent(out)::fpres
real*8,intent(in)::munacn(1:2, 1:2, 1:npoin)
real*8,intent(in)::snsigm(1:ndimn, 1:npoin)
real*8,intent(in)::munacu(1:ndimn, 1:npoin)
integer, intent(in):: itime
!
!...Local variable...
!
real*8:: pres0, bnxpr, bnypr, bnx, bny
real*8:: presi, prese
real*8:: dtime, htkid, taokid
real*8:: coorp(2, nvfac) !nvfac: no of vertices of one face
integer::ifa, ip
integer::ipf(nvfac)
real*8:: bsnv(2, npoin), bsnpf(2, npoin)
real*8::munaci(2,2)
real*8::lnx,lny,coefl,coefr,rhsx,rhsy,detma
real*8:: bsnx, bsny, prex, prin
real*8:: rcoef(2)
real*8:: xph(2, 2)
!
real*8:: farea
!
!...Zero our fpres and pres0...
!
fpres = 0.d0
pres0 = 1.d-6
!
bsnv = 0.d0  !...area*n for velocity boundary
bsnpf = 0.d0 !...force for pressure boundary
!
if(nriem.eq.1)then
if(nfint.eq.1)then
farea = .5d0
elseif(nfint.eq.3)then
farea = 1.d0
endif
endif

!
!...Part I: The pressure force for any boundary node...
!
do 1300 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
xph(1, 1:2) = coord(1, ipf(1:2))
xph(2, 1:2) = coord(2, ipf(1:2))
!
rcoef(1) = 1.d0 - alfrz + alfrz*xph(2, 1)
rcoef(2) = 1.d0 - alfrz + alfrz*xph(2, 2)
!
bsnx = 0.d0
bsny = 0.d0
!
if(bface(3, ifa).eq.22)then !...Normal velocity boundary
!
if(bface(4, ifa).eq.221)then
!bsnx = gflag(1,ifa)*gflag(3,ifa)*farea
!bsny=  gflag(2,ifa)*gflag(3,ifa)*farea
bsnx = sqrt(2.d0)/2.d0*gflag(3,ifa)*farea
bsny =-sqrt(2.d0)/2.d0*gflag(3,ifa)*farea
!
bsnx = 0.d0*gflag(3,ifa)*farea
bsny =-1.d0*gflag(3,ifa)*farea
!
if(ncase.eq.1.or.ncase.eq.2.or.ncase.eq.12)then
if(xph(2,1).lt.1.d-6.and.xph(2,2).lt.1.d-6)then
bsnx = 0.d0*gflag(3,ifa)*farea
bsny =-1.d0*gflag(3,ifa)*farea
elseif(abs(xph(2,1)-1.d0).lt.1.d-6.and.abs(xph(2,2)-1.d0).lt.1.d-6)then
bsnx = 0.d0*gflag(3,ifa)*farea
bsny = 1.d0*gflag(3,ifa)*farea
endif
endif
!
!bsnx = sin(3.d0/8.d0*pi)*gflag(3,ifa)*farea
!bsny =-cos(3.d0/8.d0*pi)*gflag(3,ifa)*farea
elseif(bface(4, ifa).eq.222)then
!bsnx = gflag(1,ifa)*gflag(3,ifa)*farea
!bsny=  gflag(2,ifa)*gflag(3,ifa)*farea
bsnx = -1.0*gflag(3,ifa)*farea
bsny=   0.d0*gflag(3,ifa)*farea

if(ncase.eq.1.)then
if(xph(1,1).lt.1.d-6.and.xph(1,2).lt.1.d-6)then
bsnx = -1.d0*gflag(3,ifa)*farea
bsny =  0.d0*gflag(3,ifa)*farea
elseif(abs(xph(1,1)-1.d0).lt.1.d-6.and.abs(xph(1,2)-1.d0).lt.1.d-6)then
bsnx =  1.d0*gflag(3,ifa)*farea
bsny =  0.d0*gflag(3,ifa)*farea
endif
endif

endif
!
bsnv(1, ipf(1)) = bsnv(1, ipf(1)) + bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0
bsnv(2, ipf(1)) = bsnv(2, ipf(1)) + bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0

bsnv(1, ipf(2)) = bsnv(1, ipf(2)) + bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0
bsnv(2, ipf(2)) = bsnv(2, ipf(2)) + bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0
!
elseif(bface(3, ifa).eq.21)then !...Pressure described
!
bsnx = gflag(1,ifa)*gflag(3,ifa)*farea
bsny=  gflag(2,ifa)*gflag(3,ifa)*farea
!
if(ncase.eq.3)then
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0
!
elseif(ncase.eq.4)then !...Kidder shell
!if(rkstg.eq.1)then
dtime = (itime-1.d0)*dtfix + crklg(rkstg)*dtfix
!elseif(rkstg.eq.2)then
!dtime = (itime)*dtfix
!endif
prin = 0.1d0
prex = 10.d0
!
pres0 = prin*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0

elseif(ncase.eq.14)then !...Coggshell

dtime = (itime-1.d0)*dtfix
pres0 = (gamlg-1.d0)*(1.d0-dtime)**(-2.25d0)*(3.d0/8.d0/(1.d0-dtime))**2

bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0*coord(1, ipf(1))**2
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0*coord(1, ipf(1))**2

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0*coord(1, ipf(2))**2
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0*coord(1, ipf(2))**2

endif
!
elseif(bface(3, ifa).eq.23)then
!
bsnx =  gflag(1,ifa)*gflag(3,ifa)*farea
bsny =  gflag(2,ifa)*gflag(3,ifa)*farea
!
if(ncase.eq.4)then !...Kidder shell
!
!if(rkstg.eq.1)then
dtime = (itime-1.d0)*dtfix + crklg(rkstg)*dtfix
!elseif(rkstg.eq.2)then
!dtime = (itime)*dtfix
!endif

prin = 0.1d0
prex = 10.d0
!
pres0 = prex*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx*(2.d0*rcoef(1) + rcoef(2))/3.d0
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny*(2.d0*rcoef(1) + rcoef(2))/3.d0

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx*(2.d0*rcoef(2) + rcoef(1))/3.d0
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny*(2.d0*rcoef(2) + rcoef(1))/3.d0
!
endif
!
endif
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1300 enddo
!
!...Part II: Get the pressure from the normal velocity prescription
!
do 1350 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
if(bface(3, ifa).eq.22)then
!
do ip = 1, 2
!
detma = munacn(1, 1, ipf(ip))*munacn(2, 2, ipf(ip)) - munacn(2, 1, ipf(ip))*munacn(1, 2, ipf(ip))
munaci(1, 1) = munacn(2, 2, ipf(ip))/detma
munaci(1, 2) =-munacn(1, 2, ipf(ip))/detma
munaci(2, 1) =-munacn(2, 1, ipf(ip))/detma
munaci(2, 2) = munacn(1, 1, ipf(ip))/detma
!
!print*,'ipf',ifa, detma, ip
!
lnx = bsnv(1, ipf(ip))
lny = bsnv(2, ipf(ip))
!
coefl = (munaci(1,1)*lnx + munaci(1,2)*lny)*lnx + (munaci(2,1)*lnx + munaci(2,2)*lny)*lny
!
rhsx = munacu(1, ipf(ip))  - snsigm(1, ipf(ip)) - bsnpf(1, ipf(ip))
rhsy = munacu(2, ipf(ip))  - snsigm(2, ipf(ip)) - bsnpf(2, ipf(ip))
!
coefr = (munaci(1,1)*rhsx + munaci(1,2)*rhsy)*lnx + (munaci(2,1)*rhsx + munaci(2,2)*rhsy)*lny
!
fpres(1, ipf(ip)) = coefr/max(1.d-16,coefl)*lnx
fpres(2, ipf(ip)) = coefr/max(1.d-16,coefl)*lny
enddo
!
endif
!
!print*,'ifa',ipf(1:2),ifa,fpres(1,:)
!
1350 enddo
!
!...Part III: Get the whole force...
!
fpres(:, :) = fpres(:, :) + bsnpf(:, :)
!
!print*,'ifa',fpres(1:2,500),bsnpf(1:2, 500)
!
end subroutine getbc_lag

!
!...Subroutine to deal with BC with prescribed pressure and normal velocity
!
subroutine getboundary_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),          intent(out)::fpres
real*8,intent(in)::munacn(1:2, 1:2, 1:npoin)
real*8,intent(in)::snsigm(1:ndimn, 1:npoin)
real*8,intent(in)::munacu(1:ndimn, 1:npoin)
!
!...Local variable...
!
real*8:: pres0, bnxpr, bnypr, bnx, bny
real*8:: presi, prese
real*8:: dtime, htkid, taokid
real*8:: coorp(2, nvfac) !nvfac: no of vertices of one face
integer::ifa, ip
integer::ipf(nvfac)
real*8::bnpre(2, npoin)
real*8::munaci(2,2)
real*8::lnx,lny,coefl,coefr,rhsx,rhsy,detma
!
!...Zero our fpres and pres0...
!
fpres = 0.d0
pres0 = 0.d0
!
if(ncase.eq.1)then
!
bnpre = 0.d0
!
!  print*,'gflag',gflag(1:3, 12)
!
do 1300 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
bnxpr =  gflag(1, ifa)*gflag(3,ifa)*0.5d0
bnypr =  gflag(2, ifa)*gflag(3,ifa)*0.5d0
!
bnpre(1, ipf(1)) = bnpre(1, ipf(1)) + bnxpr
bnpre(2, ipf(1)) = bnpre(2, ipf(1)) + bnypr

bnpre(1, ipf(2)) = bnpre(1, ipf(2)) + bnxpr
bnpre(2, ipf(2)) = bnpre(2, ipf(2)) + bnypr
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1300 enddo
!
do 1350 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
!
do ip = 1, 2
detma = munacn(1, 1, ipf(ip))*munacn(2, 2, ipf(ip)) - munacn(2, 1, ipf(ip))*munacn(1, 2, ipf(ip))
munaci(1, 1) = munacn(2, 2, ipf(ip))/detma
munaci(1, 2) =-munacn(1, 2, ipf(ip))/detma
munaci(2, 1) =-munacn(2, 1, ipf(ip))/detma
munaci(2, 2) = munacn(1, 1, ipf(ip))/detma
!
lnx = bnpre(1, ipf(ip))
lny = bnpre(2, ipf(ip))
!
coefl = (munaci(1,1)*lnx + munaci(1,2)*lny)*lnx + (munaci(2,1)*lnx + munaci(2,2)*lny)*lny
!
rhsx = munacu(1, ipf(ip))  - snsigm(1, ipf(ip))
rhsy = munacu(2, ipf(ip))  - snsigm(2, ipf(ip))
coefr = (munaci(1,1)*rhsx + munaci(1,2)*rhsy)*lnx + (munaci(2,1)*rhsx + munaci(2,2)*rhsy)*lny
!
fpres(1, ipf(ip)) = coefr/coefl*lnx
fpres(2, ipf(ip)) = coefr/coefl*lny
enddo
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1350 enddo
!
endif
!
end subroutine getboundary_lag
!
!...subroutine to deal with norma velocity prescribed...
!
subroutine getbcvn_lag(bface, intfac, gflag, ustar)
use constant
implicit none
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ndimn,1:npoin),        intent(inout)::ustar
!
real*8:: bnx, bny, vdon
integer::ifa
integer::ipf(2)
!
!  print*,'gflag',gflag(1:3, 12)
!
do 1600 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
if(bface(3, ifa).eq.22)then
!
bnx = gflag(1, ifa)
bny = gflag(2, ifa)
!
vdon = bnx*ustar(1, ipf(1)) + bny*ustar(2, ipf(1))
!
ustar(1, ipf(1)) = ustar(1, ipf(1)) - vdon*bnx
ustar(2, ipf(1)) = ustar(2, ipf(1)) - vdon*bny
!
bnx = gflag(1, ifa)
bny = gflag(2, ifa)
!
vdon = bnx*ustar(1, ipf(2)) + bny*ustar(2, ipf(2))
!
ustar(1, ipf(2)) = ustar(1, ipf(2)) - vdon*bnx
ustar(2, ipf(2)) = ustar(2, ipf(2)) - vdon*bny
!
if(ifa.eq.94) print*,'getbcfc',bnx,bny,ustar(1:2, ipf(1)),ipf(1)
!
endif
1600 enddo

end subroutine getbcvn_lag
!
!...Subroutine to deal with BC with prescribed pressure and normal velocity
!
subroutine getbcve_exact(bface, intfac, gflag,  ustar, coord, itime)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ngflg,1:nbfac),           intent(in)::gflag
real*8,dimension(1:ndimn,1:npoin),          intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),          intent(inout)::ustar
integer:: itime
!
!...Local variable...
!
real*8:: pres0, bnxpr, bnypr, bnx, bny
real*8:: presi, prese
real*8:: dtime, htkid, taokid
real*8:: xp1, yp1, xp2, yp2
integer::ifa
integer::ipf(2)
!
!...bface(3, ifa)==25 denotes the Dirichlet BC...
!
pres0 = 0.d0
!
if(ncase.eq.3)then
!
pres0 = 1.d-6  !...Noh test problem
!
!  print*,'gflag',gflag(1:3, 12)
!
do 1300 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
if(bface(3, ifa).eq.21)then
!
!
endif
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1300 enddo
!
elseif(ncase.eq.4)then
!
presi = 0.1d0
prese = 10.d0
!
!  if(rkstg.eq.1)then
dtime = (itime)*dtfix
!  elseif(rkstg.eq.2)then
!   dtime = (itime+1)*dtfix
!  endif
!
htkid = -dtime/(sqrt(1.d0-(dtime)**2))/7.265d-3
!
!  print*,'gflag',itime,dtime,itime,dtfix
!
do 1350 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
!   print*,'gflag',ifa,ipf(1:2),gflag(1:2, ifa),bface(3, ifa)
!
!if(bface(3, ifa).eq.22)then !...inner circle
!
xp1 = coord(1, ipf(1)); yp1 = coord(2, ipf(1))
xp2 = coord(1, ipf(2)); yp2 = coord(2, ipf(2))
!
ustar(1, ipf(1)) = xp1*htkid
ustar(2, ipf(1)) = yp1*htkid
!
ustar(1, ipf(2)) = xp2*htkid
ustar(2, ipf(2)) = yp2*htkid
!
!endif
!
!
!elseif(bface(3, ifa).eq.23)then !...outer circle
!
!
!xp1 = coord(1, ipf(1)); yp1 = coord(2, ipf(1))
!xp2 = coord(1, ipf(2)); yp2 = coord(2, ipf(2))
!
!ustar(1, ipf(1)) = xp1*htkid
!ustar(2, ipf(1)) = yp1*htkid
!
!ustar(1, ipf(2)) = xp2*htkid
!ustar(2, ipf(2)) = yp2*htkid
!
!endif
!
!print*,'ifa',ipf(1:2),ifa,bnorm(1:2,ipf(1)),bnorm(1:2,ipf(2))
!
1350 enddo
!
endif
!
end subroutine getbcve_exact
!
!...subroutine: Calculate the nodal velocity U_p^* for linear cell with vilar...
!
subroutine getndvelo_lag_linear_vilar(gflag,gelag,bface,intfac,inpoel,coord,unkno,ustar, fstar, itime)
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
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:2,1:nvtri, 1:nelem),  intent(out)::fstar !...Riemann forces
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin, itime
integer::iel,ier
!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:2)     :: ipf
integer::indnd(npoin)

!...local real array
real*8::munacn_rie(2, 2), munacu_rie(1:2), snsigm_rie(1:2)
real*8::munaci(2, 2)
real*8::muutg(1:2, 1:2, 1:npoin), mutag(npoin)
real*8,dimension(1:nq,1:nvtri)::unknv
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:nvtri)::murie
real*8::vnorm(1:3, 1:2, 1:nvtri)
real*8::sigma(1:2, 1:2, 1:nvtri)
real*8::aujmp(1:2, 1:nvtri)
real*8,dimension(1:3, 1:nvtri)::b
real*8,dimension(1:nvtri):: xv, yv
real*8::usnom(2), ustng(2)
real*8::xpf(1:2,1:3),xcf(1:2)
!...local real number
real*8::unorm, tfx, tfy, nfx, nfy,signm, utang
real*8::eps,c00,c05,c10,c20
real*8::dr,ds,farea,larea,rc,sc,acnx,acny
real*8::bnx, bny
real*8::rho, rhoct, sdctr, pctr, uctr, vctr, ectr
real*8::uvtx,vvtx,evtx, pvtx
real*8::dux,duy,deltu
real*8::dwav1,dwav2
real*8::anx, any, len1, len2, tfx1, tfy1
real*8::uvtxr,vvtxr,evtxr, pvtxr,uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufac,mufal,mufar,rhsu1,rhsu2,detma
real*8,allocatable:: bnorm(:,:), fpres(:,:)
real*8,allocatable:: cnsup(:), munacn(:,:,:), bpres(:) !...Count no surrounding one vertex....
real*8,allocatable:: usold(:,:), munacu(:,:), snsigm(:,:)
real*8,allocatable:: munacl(:,:,:,:,:)
real*8,allocatable:: snsigml(:,:,:,:), munaul(:,:,:,:)
!
data eps   / 1.0d-06 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
allocate (munacn(1:2,1:2,1:npoin), cnsup(1:npoin))
allocate (usold(1:ndimn, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munacl(1:2,1:2,1:2, 1:nvtri, 1:nelem), munaul(1:ndimn, 1:2, 1:nvtri,  1:nelem),&
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
!...Mark the boundary nodes for case 2...
!
if(ncase.eq.2)then
do ifa = 1, nbfac
ipf(1:2) = intfac(3:4, ifa)
indnd(ipf(1:2)) = 1
enddo
endif
!
!...shape functions
!
dr = .5d0
ds = .5d0
rc = 1.d0/3.d0
sc = rc
!
xv(1) = 0.d0;  yv(1) = 0.d0
xv(2) = 1.d0;  yv(2) = 0.d0
xv(3) = 0.d0;  yv(3) = 1.d0
!
do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xv(iv)-rc)/dr
b(3, iv) = (yv(iv)-sc)/ds
enddo
!
!...Part I: Get the averaged nodal velocity...
!
do 200 ie = 1,nelem !...(1)ie = 1,nelem
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
! if(ip(iv)==36) print*,'average p21',unknv(2:3, iv),ip(iv),ie, unkno(1, 2:3, ie)
enddo
!
!...Accumulate nodal velocity...
!
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
!...Part II:Get the related coefficients for Riemann-like
!
!
!...Zero out munacn, munacu, snsigm...
!
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
!...Deal with the high-order nodes...
!
muutg = 0.d0
mutag = 0.d0
!
do 250 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
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
!...Loop 1: Calculate the tensor stress and a_c at vertices...
!
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
!...Stress tensor...
!
rho  = 1.d0/unknv(1, iv)
uvtx = unknv(2, iv)
vvtx = unknv(3, iv)
evtx = unknv(4, iv)
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
!
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
!if(ip(iv)==36) print*,'adjumpxxx9471', acnx,acny,aujmp(1:2, ip(iv)),sqrt(acnx**2 + acny**2),&
!                     vlave(1:2, ip(iv)) , unknv(2:3, iv), aujmp(1:2,iv), ip(iv),ie,iv
enddo
!
!if(ie==3) print*,'vnotmxxx',vnorm(3,1,1),gelag(1, 3, 9),aujmp(1:2,1)
!
!
!...Loop 2: Get the impedence at the vertices...
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
murie(iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
! if(ip(iv).eq.5) print*,'murie22', sdctr,rhoct,deltu,vlave(1, ip(iv)),unknv(2, iv),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, 3 !...Low-order nodes...
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
!
!
! if(ip(iv).eq.40)then
call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), unknv(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
! endif
!
!
! if(ip(iv).eq.40) print*,'p19 muacn(28) pre++', munacn(1:2,1:2,ip(iv)),ie,iv,ifa,munacu(1:2,ip(iv))
!
!  if(abs(aujmp(1, iv))+abs(aujmp(2, iv)).lt.1.d-7) then
!   aujmp(1, iv)=vnorm(1, ifa, iv); aujmp(2, iv)=vnorm(2, ifa, iv)
!  endif
!
munacn(1:2, 1, ip(iv)) = munacn(1:2, 1, ip(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ip(iv)) = munacn(1:2, 2, ip(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ip(iv)) = munacu(1:2, ip(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ip(iv)) = snsigm(1:2, ip(iv)) + snsigm_rie(1:2)!
!
!   if(ip(iv).eq.40) print*,'p36 muacn(vv) post',ie,murie(iv),munacn(1:2,1:2,ip(iv)),vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),&
!                            unknv(2:3,iv),&
!                            sigma(1:2, 1:2, iv),munacn_rie, munacu_rie(1:2), snsigm_rie(1:2)
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
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
!...Local variable...
!
munacl(1:2, 1, ifa, iv, ie) =  munacn_rie(1:2, 1)
munacl(1:2, 2, ifa, iv, ie) =  munacn_rie(1:2, 2)

!
munaul(1:2, ifa, iv, ie) =  munacu_rie(1:2)
!
snsigml(1:2, ifa, iv, ie)=  snsigm_rie(1:2)
!
enddo
!
enddo
!
250 enddo  !...(1)ie = 1,nelem
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
!...Zero out the bnorm
!
!
!...
!...Zero out boundary pressure...
!
!
!...Fourth part: Solve the nodal velocity...
!
!...deactivate the pseudo pressure....
!
fpres = 0.d0
!
 call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar,  itime)
!
!call getbcvnforce_lag(bface, intfac, gflag, munacu, snsigm, fpres)
!
!...4.1: Update the Riemann forces at every node...
!
!
do ipoin = 1, npoin
if(indnd(ipoin).eq.0)then
!if(idcrd(ipoin).eq.0)then
detma = munacn(1, 1, ipoin)*munacn(2, 2, ipoin) - munacn(2, 1, ipoin)*munacn(1, 2, ipoin)
munaci(1, 1) = munacn(2, 2, ipoin)/detma
munaci(1, 2) =-munacn(1, 2, ipoin)/detma
munaci(2, 1) =-munacn(2, 1, ipoin)/detma
munaci(2, 2) = munacn(1, 1, ipoin)/detma
!
rhsu1 = munacu(1, ipoin)  - snsigm(1, ipoin) - fpres(1, ipoin)
rhsu2 = munacu(2, ipoin)  - snsigm(2, ipoin) - fpres(2, ipoin)
!
ustar(1, ipoin) = munaci(1, 1)*rhsu1 + munaci(1, 2)*rhsu2
ustar(2, ipoin) = munaci(2, 1)*rhsu1 + munaci(2, 2)*rhsu2
!elseif(idcrd(ipoin).eq.1)then
!  if(muacn(ipoin).gt.1.d-9)then
!   ustar(1, ipoin) = (munacu(1, ipoin) - snsigm(1, ipoin) + fpres(1, ipoin))/munacn(ipoin)
!   ustar(2, ipoin) = (munacu(2, ipoin) - snsigm(2, ipoin) + fpres(2, ipoin))/munacn(ipoin)
!  endif
!endif
!
! if(ipoin.eq.40) print*,munaci, munacn(1:2,1:2,ipoin)
endif
enddo
!
!....Bd velocity
!
! print*,'ustar--',ustar(1:2, 14),munacu(1:2,14) ,snsigm(1:2,14), munacn(1:2,1:2,14)
!
!
if(ncase.eq.1)then
!
do 900 ifa = 1 , nbfac
!
ipf(1:2) = intfac(3:4, ifa)
!
ustar(1, ipf(1)) = sin(pi*coord(1,ipf(1)))*cos(pi*coord(2,ipf(1)))
ustar(2, ipf(1)) =-cos(pi*coord(1,ipf(1)))*sin(pi*coord(2,ipf(1)))
!
ustar(1, ipf(2)) = sin(pi*coord(1,ipf(2)))*cos(pi*coord(2,ipf(2)))
ustar(2, ipf(2)) =-cos(pi*coord(1,ipf(2)))*sin(pi*coord(2,ipf(2)))
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
ustar(1:2, 1:4) = 0.d0
!
!ustar(1:2, 1) = 0.d0;ustar(1:2, 81) = 0.d0;ustar(1:2, 6481) = 0.d0;ustar(1:2, 6561) = 0.d0
!
endif
!
!   print*,'ustar',ustar(1:2, 14)!,ustar(1:2, 20),ustar(1:2, 21)
!  if(ustar(2, 18).gt.0.d0) then
!   print*,'Wrong movng direction'!,munacu(2, 28) - snsigm(2, 28) , fpres(2, 28),munacn(28)
!   stop
!  endif
!
 call getbcvn_lag(bface, intfac, gflag, ustar)
!
! call getbcve_exact(bface, intfac, gflag, ustar, coord, itime)
!
!...4.2: Update the Riemann forces at every node...
!
do ie = 1, nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
unknv = 0.d0
!
!...Low order nodes...
!
do iv = 1, 3
!
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
!
do ifa = 1, 2
!
fstar(1, ifa, iv, ie) = snsigml(1, ifa, iv, ie) +&
munacl(1, 1, ifa, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv))+&
munacl(2, 1, ifa, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(1, 1, iv, ie)
fstar(2, ifa, iv, ie) = snsigml(2, ifa, iv, ie) + &
munacl(1,2,ifa, iv, ie)*(ustar(1, ip(iv))-unknv(2,iv))+&
munacl(2, 2, ifa, iv, ie)*(ustar(2, ip(iv))-unknv(3,iv)) !- munaul(1, 1, iv, ie)
enddo
!
enddo
!
enddo
!

deallocate (munacn, cnsup, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
end subroutine getndvelo_lag_linear_vilar
!
!...subroutine: Calculate the F^* N dsfor all curved faces...
!
subroutine getfnds_lagnew_curv_hybrid(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngflg,1:nbfac),         intent(inout)::gflag  !...Geometry of face in lagrangian
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(inout)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad),      intent(inout)::gelagq
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
!...Allocatable array...
real*8, allocatable::cordl(:, :)
!
data c16   /0.1666666666666666d0 /
!
allocate(cordl(1:ndimn, 1:npoin))
!
cordl = coord
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
!...Part 1: Get the control point coordinates...
!
do 40 ie=1, ntria !...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...coordinates
!
cordl(1:2, ip(4)) = 0.5d0*(4.d0*coord(1:2, ip(4)) - coord(1:2, ip(1)) - coord(1:2, ip(2))) !
cordl(1:2, ip(5)) = 0.5d0*(4.d0*coord(1:2, ip(5)) - coord(1:2, ip(2)) - coord(1:2, ip(3))) !...Control points...
cordl(1:2, ip(6)) = 0.5d0*(4.d0*coord(1:2, ip(6)) - coord(1:2, ip(3)) - coord(1:2, ip(1))) !
!
40 enddo
!
!...qauds
!
do 45 ie=1, nquad !...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...coordinates
!
cordl(1:2, ipq(5)) = 0.5d0*(4.d0*coord(1:2, ipq(5)) - coord(1:2, ipq(1)) - coord(1:2, ipq(2))) !
cordl(1:2, ipq(6)) = 0.5d0*(4.d0*coord(1:2, ipq(6)) - coord(1:2, ipq(2)) - coord(1:2, ipq(3))) !...Control points...
cordl(1:2, ipq(7)) = 0.5d0*(4.d0*coord(1:2, ipq(7)) - coord(1:2, ipq(3)) - coord(1:2, ipq(4))) !
cordl(1:2, ipq(8)) = 0.5d0*(4.d0*coord(1:2, ipq(8)) - coord(1:2, ipq(1)) - coord(1:2, ipq(4))) !
!
45 enddo
!
!...Trias
!
do 50 ie=1, ntria!...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...coordinates
!
xp(1, 1:nvtri) = cordl(1, ip(1:nvtri))
xp(2, 1:nvtri) = cordl(2, ip(1:nvtri))
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
!...Edge 4
!
anx =   xp(2, 4) - xp(2, 1)
any = -(xp(1, 4) - xp(1, 1))
!
gelag(1, 4 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 4 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 4 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 5
!
anx =   xp(2, 2) - xp(2, 4)
any = -(xp(1, 2) - xp(1, 4))
!
gelag(1, 5 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 5 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 5 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 6
!
anx =   xp(2, 5) - xp(2, 2)
any = -(xp(1, 5) - xp(1, 2))
!
gelag(1, 6 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 6 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 6 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 7
!
anx =   xp(2, 3) - xp(2, 5)
any = -(xp(1, 3) - xp(1, 5))
!
gelag(1, 7 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 7 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 7 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 8
!
anx =   xp(2, 6) - xp(2, 3)
any = -(xp(1, 6) - xp(1, 3))
!
gelag(1, 8 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 8 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 8 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 9
!
anx =   xp(2, 1) - xp(2, 6)
any = -(xp(1, 1) - xp(1, 6))
!
gelag(1, 9 , ie) = anx/sqrt(anx**2 + any**2)
gelag(2, 9 , ie) = any/sqrt(anx**2 + any**2)
gelag(3, 9 , ie) = sqrt(anx**2 + any**2)
!
50 enddo  !...(1)ifa=1,nelem
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
xpq(1, 1:nvqua) = cordl(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = cordl(2, ipq(1:nvqua))
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
!...Edge 5
!
anx =   xpq(2, 5) - xpq(2, 1)
any = -(xpq(1, 5) - xpq(1, 1))
!
gelagq(1, 5 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 5 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 5 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 6
!
anx =   xpq(2, 2) - xpq(2, 5)
any = -(xpq(1, 2) - xpq(1, 5))
!
gelagq(1, 6 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 6 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 6 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 7
!
anx =   xpq(2, 6) - xpq(2, 2)
any = -(xpq(1, 6) - xpq(1, 2))
!
gelagq(1, 7 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 7 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 7 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 8
!
anx =   xpq(2, 3) - xpq(2, 6)
any = -(xpq(1, 3) - xpq(1, 6))
!
gelagq(1, 8 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 8 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 8 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 9
!
anx =   xpq(2, 7) - xpq(2, 3)
any = -(xpq(1, 7) - xpq(1,3))
!
gelagq(1, 9 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 9 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 9 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 10
!
anx =   xpq(2, 4) - xpq(2, 7)
any = -(xpq(1, 4) - xpq(1, 7))
!
gelagq(1, 10 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 10 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 10 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 11
!
anx =   xpq(2, 8) - xpq(2, 4)
any = -(xpq(1, 8) - xpq(1, 4))
!
gelagq(1, 11 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 11 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 11 , ie) = sqrt(anx**2 + any**2)
!
!...Edge 12
!
anx =   xpq(2, 1) - xpq(2, 8)
any = -(xpq(1, 1) - xpq(1, 8))
!
gelagq(1, 12 , ie) = anx/sqrt(anx**2 + any**2)
gelagq(2, 12 , ie) = any/sqrt(anx**2 + any**2)
gelagq(3, 12 , ie) = sqrt(anx**2 + any**2)
!
55 enddo  !...(1)ifa=1,nelem
!
! print*,'vnotmfn',gelag(1:3, 1:9, 50)
!
! print*,'Inside getfnds_lag'
!
deallocate(cordl)
!
end subroutine getfnds_lagnew_curv_hybrid
!
!
!...subroutine: Calculate the Riemann input for hybrid curv tria grids general Riemann solver....
!
subroutine getriem_tria_curv(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold, aflim, afvec,vnult)
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
real*8,dimension(1:ndimn,1:nvtri,1:ntria),   intent(in):: vnult
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
real*8,dimension(1:nvtri)::murie
real*8,dimension(1:nvtri):: xv,  yv
real*8,dimension(1:ndimn, 1:nvtri) :: xpt
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
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
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
vnorm(1:3, 3, 1) = gelag(1:3, 9, ie); vnorm(1:3, 4, 1) = gelag(1:3, 4, ie) !...For point ip(1)
!
vnorm(1:3, 1, 2) = gelag(1:3, 1, ie); vnorm(1:3, 2, 2) = gelag(1:3, 2, ie) !...For point ip(2)
vnorm(1:3, 3, 2) = gelag(1:3, 5, ie); vnorm(1:3, 4, 2) = gelag(1:3, 6, ie) !...For point ip(2)
!
vnorm(1:3, 1, 3) = gelag(1:3, 2, ie); vnorm(1:3, 2, 3) = gelag(1:3, 3, ie) !...For point ip(3)
vnorm(1:3, 3, 3) = gelag(1:3, 7, ie); vnorm(1:3, 4, 3) = gelag(1:3, 8, ie) !...For point ip(3)
!
vnorm(1:3, 1, 4) = gelag(1:3, 4, ie); vnorm(1:3, 2, 4) = gelag(1:3, 5, ie) !...For point ip(4)
!
vnorm(1:3, 1, 5) = gelag(1:3, 6, ie); vnorm(1:3, 2, 5) = gelag(1:3, 7, ie) !...For point ip(5)
!
vnorm(1:3, 1, 6) = gelag(1:3, 8, ie); vnorm(1:3, 2, 6) = gelag(1:3, 9, ie) !...For point ip(6)

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
!...Remove null mode...
!
unknvt(2:3, iv)=unknvt(2:3, iv)-vnult(1:2,iv,ie)
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
unknvt(2, iv) = uvtx-vnult(1,iv,ie)
unknvt(3 ,iv) = vvtx-vnult(2,iv,ie)
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
do iv  = 1, 3
do ifa = 1, 4 !...Every corner consists of 2 faces...
!
!...Call Riemann solver...
!
!call getriecoef_matrixnew(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
!unknvt(2:3, iv), sigma(1:2, 1:2, iv),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
unknvt(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
if(ifa.le.2)then
!

!
elseif(ifa.ge.3)then
!
munacn(1:2, 1, ipt(iv)) = munacn(1:2, 1, ipt(iv)) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipt(iv)) = munacn(1:2, 2, ipt(iv)) + munacn_rie(1:2, 2)
!
munacu(1:2, ipt(iv)) = munacu(1:2, ipt(iv)) + munacu_rie(1:2)
!
snsigm(1:2, ipt(iv)) = snsigm(1:2, ipt(iv)) + snsigm_rie(1:2)!
!
endif
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
!!
do iv  = 4, nvtri
do ifa = 1, 2 !...Every corner consists of 2 faces...
!
!...Call Riemann solver...
!
!call getriecoef_matrixnew(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
!unknvt(2:3, iv), sigma(1:2, 1:2, iv),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
unknvt(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
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
!if(ipt(iv).eq.119) print*,'p19 muacn(28) prep---',murie(iv), munacn_rie(:,:),munacu(1:2,ipt(iv)),vnorm(3, ifa, iv),&
!                                                vnorm(1:2, ifa, iv),aujmp(1:2, iv),unknvt(2:3,iv)
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

end subroutine getriem_tria_curv
!
!...subroutine: Calculate the Riemann input for hybrid curved quad grids general Riemann solver....
!
subroutine getriem_quad_curv(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq,coord, coold, aflim, afvec,vnulq)
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
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(in):: vnulq
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
real*8,dimension(1:ndegr, 1:nvqua)::bq,bqv
real*8,dimension(1:nq,1:nvqua)::unknvq
real*8::aujmp(1:3, 1:nvqua)
real*8::vnorm(1:3, 1:4, 1:nvqua)
real*8::sigma(1:2, 1:2, 1:nvqua)
real*8,dimension(1:nvqua)::murie
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
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
data eps   / 1.0d-14/
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
!...Remove null mode...
!
unknvq(2:3, iv)=unknvq(2:3, iv)-vnulq(1:2,iv,ie)

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
call getrhoig_quadcurv(rhoi, xpqi)!
call getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)
!call getrhoig_quad(rhoi,r,s, xpqi(1:2,1:4))!
!call getdensity_quadllnl(r, s, xpq(1:2,1:4), xpqi(1:2,1:4), rhoi, rhon)
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
!uvtx = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
!vvtx = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
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
!...Get impedence coefficient...
!
do iv   = 1, nvqua
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie(iv) = rhoct*sdctr !+ rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!
! if(ipq(iv).eq.1922) print*,'murie22',murie(iv), sdctr,rhoct,deltu,vlave(1:2, ipq(iv)),unknvq(2:3, iv),dux,duy,ie
! if(ipq(iv).eq.840) then
!         print*,'murie22ang',unknvq(2,iv)**2+unknvq(3,iv)**2,ie,iv,vlave(1,ipq(iv))*unknvq(2,iv)+vlave(2,ipq(iv))*unknvq(3,iv)
! endif
!
enddo
!
!if(ie==94) print*,'vnotm',murie,rhoct,sdctr,uctr,vctr,ectr!,vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do iv  = 1, 4
do ifa = 1, 4 !...Every corner consists of 2 faces...
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
!...Call Riemann solver...
!
!call getriecoef_matrixnew(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
!unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
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
!   if(ipq(iv).eq.1923) print*,'p36 muacn(vv) postmuna',ie

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
!if(ipq(iv).eq.1922) print*,'Riem',aujmp(1:3, iv),munacu_rie(1:2),snsigm_rie(1:2),ie, ifa,iv,ipq(iv)
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
!call getriecoef_matrixnew(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:3, iv), &
!unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
call getriecoef_vilar(murie(iv), vnorm(3, ifa, iv), vnorm(1:2, ifa, iv), aujmp(1:2, iv), &
unknvq(2:3, iv), sigma(1:2, 1:2, iv),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
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

end subroutine getriem_quad_curv
!
!...subroutine: Calculate the nodal velocity U_p^* (mass center) for hybrid mesh with general Riemann solver...
!
subroutine getndvelo_lag_mc_curv(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
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
real*8,dimension(1:ndimn,1:4,1:nvtri, 1:ntria),  intent(out)::fstar !...Riemann forces
real*8,dimension(1:ndimn,1:4,1:nvqua, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec

integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)
integer::ipsin(6)

!...local real array
real*8::munacnsin(4,3),munacusin(2, 3),snsigmsin(2,3)
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:ndimn,1:nvqua,1:nquad)::vnulq
real*8,  dimension(1:nquad)::gqdmp
real*8,dimension(1:ndimn,1:nvtri,1:ntria)::vnult
real*8,  dimension(1:ntria)::gtdmp
real*8::munaci(2, 2)
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
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
indnd(ipf(1:nvfac)) = 1
enddo
endif
!
!
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
!
indnd(ipf(1:nvfac)) = 1
endif
enddo
!
!...Get averaged velocity at nodes...
!
!call getvlavenew(iptri, ipqua, geoel, vlave, unkno, aflim, afvec)
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
!...Get the dampening coefficients...
!
!call getnullve_gqdmpquadc(ipqua, geoel, ustar, unkno, gqdmp)
!
!...Use the nodal velocity from last time step...
!
vnulq = 0.d0; vnult= 0.d0
!
do iloop= 1, 1!5
!
vlave= ustar

!...Zero out munacn
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0
!
if(ntria.gt.0) call getriem_tria_curv(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold,aflim, afvec, vnult)
!
if(nquad.gt.0) call getriem_quad_curv(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec,vnulq)
!
!...Periodic boundary condition for 1D isentropic Sin problem...
if(ncase.eq.12)then
do ifa = 1, nbfac
if(bface(3, ifa).eq.31)then !...Periodic boundary...
!
ipsin(1:2) = bface(1:2, ifa); ipsin(3:4) = bface(1:2, bface(4, ifa))
ipsin(5)   = bface(5, ifa);   ipsin(6)   = bface(5, bface(4, ifa))
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
munacnsin(1,3) = munacn(1,1,ipsin(5)) + munacn(1,1,ipsin(6))
munacnsin(2,3) = munacn(1,2,ipsin(5)) + munacn(1,2,ipsin(6))
munacnsin(3,3) = munacn(2,1,ipsin(5)) + munacn(2,1,ipsin(6))
munacnsin(4,3) = munacn(2,2,ipsin(5)) + munacn(2,2,ipsin(6))
!
munacusin(1:2, 1) = munacu(1:2, ipsin(2)) + munacu(1:2, ipsin(3))
munacusin(1:2, 2) = munacu(1:2, ipsin(1)) + munacu(1:2, ipsin(4))
munacusin(1:2, 3) = munacu(1:2, ipsin(5)) + munacu(1:2, ipsin(6))
!
snsigmsin(1:2, 1) = snsigm(1:2, ipsin(2)) + snsigm(1:2, ipsin(3))
snsigmsin(1:2, 2) = snsigm(1:2, ipsin(1)) + snsigm(1:2, ipsin(4))
snsigmsin(1:2, 3) = snsigm(1:2, ipsin(5)) + snsigm(1:2, ipsin(6))
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
munacn(1,1,ipsin(5)) = munacnsin(1, 3); munacn(1,1,ipsin(6)) = munacnsin(1, 3);
munacn(1,2,ipsin(5)) = munacnsin(2, 3); munacn(1,2,ipsin(6)) = munacnsin(2, 3);
munacn(2,1,ipsin(5)) = munacnsin(3, 3); munacn(2,1,ipsin(6)) = munacnsin(3, 3);
munacn(2,2,ipsin(5)) = munacnsin(4, 3); munacn(2,2,ipsin(6)) = munacnsin(4, 3);
!
munacu(1:2, ipsin(2)) =  munacusin(1:2, 1); munacu(1:2, ipsin(3)) =  munacusin(1:2, 1);
munacu(1:2, ipsin(1)) =  munacusin(1:2, 2); munacu(1:2, ipsin(4)) =  munacusin(1:2, 2);
munacu(1:2, ipsin(5)) =  munacusin(1:2, 3); munacu(1:2, ipsin(6)) =  munacusin(1:2, 3);
!
snsigm(1:2, ipsin(2)) =  snsigmsin(1:2, 1); snsigm(1:2, ipsin(3)) =  snsigmsin(1:2, 1);
snsigm(1:2, ipsin(1)) =  snsigmsin(1:2, 2); snsigm(1:2, ipsin(4)) =  snsigmsin(1:2, 2);
snsigm(1:2, ipsin(5)) =  snsigmsin(1:2, 3); snsigm(1:2, ipsin(6)) =  snsigmsin(1:2, 3);
endif
enddo
endif !if(ncase.eq.12)then
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
!call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)
!call getboundary_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm)
 call getbc_lagmaire(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm)
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
!if(ipoin.eq.8) print*,'ep',ustar(1:2,ipoin),detma,munacn(:,:,ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin)
endif
enddo
!
!print*,ustar(1:2,1923),unkno(1,1:4,1831)!,detma,munacn(:,:,ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin),indnd(7594)
!
!....Bd velocity
!
!print*,'ustar--',indnd(65),ustar(1:2, 65),itime!/sqrt(ustar(1, 738)**2+ustar(2, 738)**2),munacu(1:2,36) ,snsigm(1:2, 36), munacn(36)

!if(ntria.gt.0) call getvelo_mpt_marie(ustar,gelag,intfac,inpoel,coord,unkno,indnd)
!if(nquad.gt.0) call getvelo_mpt_mariequad(ustar,gelagq,intfac,ipqua,coord,unkno,indnd)
!
!print*,'ustar',ustar(1,:)

!call getvelo_mpt_curv(ustar,geoel,gelag,gelagq,intfac,iptri,ipqua,coord,coold,unkno,indnd, aflim, afvec, vlave, vnulq)
call getvelo_mpt_curv2(ustar,geoel,gelag,gelagq,intfac,iptri,ipqua,coord,unkno,indnd, aflim, afvec, vlave, vnulq,munacn)
!print*,'ustar--',ustar(1:2,1922)
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

!...Specified boundary nodal velocity

if(ncase.eq.13)then !...Saltzman
ustar(2,ipf(1:nvfac)) = 0.d0
ustar(1,ipf(1:nvfac)) = 1.d0
endif

endif
!
!...impose exact solution along the Boundary for TGV case
!
if(ncase.eq.-1)then
ustar(1, ipf(1:nvfac)) = sin(pi*coord(1,ipf(1:nvfac)))*cos(pi*coord(2,ipf(1:nvfac)))
ustar(2, ipf(1:nvfac)) =-cos(pi*coord(1,ipf(1:nvfac)))*sin(pi*coord(2,ipf(1:nvfac)))
endif
!
!...Linearize the mid-point velocity at the boundary...
!
!ustar(1, ipf(3)) = 0.5d0*(ustar(1, ipf(1))+ustar(1, ipf(2)))
!ustar(2, ipf(3)) = 0.5d0*(ustar(2, ipf(1))+ustar(2, ipf(2)))
!
900 enddo

!

if(ntria.gt.0)then
! if(iloop.eq.1) call getcoef_nuvetriac(iptri, geoel, ustar, unkno, gtdmp)
! call getnullmd_trialp(iptri, geoel, vnult, ustar, coord, gtdmp)
endif

if(nquad.gt.0)then
! if(iloop.eq.1) call getcoef_nuvequadc(ipqua, geoel, ustar, unkno, gqdmp)
! call getnullmd_quadlp(ipqua, geoel, vnulq, ustar, coord, gqdmp)
endif
!
!print*,'vnulq',vnulq(1:2,1:8,1)
!
!do 910 ifa = 1 , -nafac
!ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
!ustar(1, ipf(3)) = 0.5d0*(ustar(1, ipf(1))+ustar(1, ipf(2)))
!ustar(2, ipf(3)) = 0.5d0*(ustar(2, ipf(1))+ustar(2, ipf(2)))
!910 enddo
!
enddo !iloop
!
!...Imposing the zero normal velocity for BC...
!
! call getbcvn_lag(bface, intfac, gflag, ustar)
!
!call getbcve_exact(bface, intfac, gflag, ustar, coord, itime)
!
!print*,'ustar',ustar(1, 37:50)!,ustar(1:2, 20),ustar(1:2, 21)
!
!ustar(:,1) = 0.d0
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
do iv = 1, 3
!
do ifa =1, 4
fstar(1, ifa, iv, ie) = snsigmlt(1, ifa, iv, ie) + &
munaclt(1, 1, ifa, iv, ie)*ustar(1, ipt(iv))+&
munaclt(2, 1, ifa, iv, ie)*ustar(2, ipt(iv)) - munault(1, ifa, iv, ie)
fstar(2, ifa, iv, ie) = snsigmlt(2, ifa, iv, ie) + &
munaclt(1, 2,ifa, iv, ie)*ustar(1, ipt(iv))+&
munaclt(2, 2, ifa, iv, ie)*ustar(2, ipt(iv)) - munault(2, ifa, iv, ie)
!
enddo
enddo
!
do iv = 4, nvtri
!
do ifa = 1, 2
fstar(1, ifa, iv, ie) = snsigmlt(1, ifa, iv, ie) + &
munaclt(1, 1, ifa, iv, ie)*ustar(1, ipt(iv))+&
munaclt(2, 1, ifa, iv, ie)*ustar(2, ipt(iv)) - munault(1, ifa, iv, ie)
fstar(2, ifa, iv, ie) = snsigmlt(2, ifa, iv, ie) + &
munaclt(1, 2,ifa, iv, ie)*ustar(1, ipt(iv))+&
munaclt(2, 2, ifa, iv, ie)*ustar(2, ipt(iv)) - munault(2, ifa, iv, ie)
!
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
!if(ie.eq.1.or.ie.eq.5)print*,'velo',ie,ustar(1:2,ipq(6))
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
!if(ie==1.and.iv.eq.1)then
!print*,'fire',iv,ie, ifa, fstarq(1, 3:4, iv, ie),snsigmlq(1,3:4,iv,ie),munaclq(1:2,1,3:4,iv,ie),munaulq(1,3:4,iv,ie)
!endif

!if(ie.eq.5.and.iv.eq.2)then
!print*,'fire',iv,ie, ifa, fstarq(1, 1:4, iv, ie)!snsigmlq(1,1:2,iv,ie),munaclq(1:2,iv,ie),ipq(:),munaulq(1,1:2,iv,ie)
!endif

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
!
enddo
!
!print*,'genndve'
!
deallocate (munacn, bpres, fpres)
deallocate (usold, munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lag_mc_curv
!
!...subroutine: Calculate the nodal velocity U_p^* (mass center) with pseduo curved mesh...
!
subroutine getndvelo_lag_mc_curvpseudo(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
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
real*8,dimension(1:ndimn,1:4,1:nvtri, 1:ntria),  intent(out)::fstar !...Riemann forces
real*8,dimension(1:ndimn,1:4,1:nvqua, 1:nquad),  intent(out)::fstarq !...Riemann forces
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec

integer:: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer::indnd(npoin)

!...local real array
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8,dimension(1:ndimn,1:nvqua,1:nquad):: vnulq
real*8::munaci(2, 2)
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
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
indnd(ipf(1:nvfac)) = 1
enddo
endif
!
!
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
if(bface(3, ifa).eq.25)then
!
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
!...Use the nodal velocity from last time step...
!
do iloop= 1, 4
!
!vlave= ustar
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
!if(ntria.gt.0) call getriem_tria_curv(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
!munaclt, munault, snsigmlt, coord, coold,aflim, afvec)
!
if(nquad.gt.0) call getriem_quad_curv(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec,vnulq)
!
!
!...Third part: Impose the boundary condition
!
!...Right now, only prescribed normal velocity is given...
!
call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)!
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
!if(ipoin.eq.104) print*,ustar(1:2,ipoin),detma,munacn(:,:,ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin)
endif
enddo
!
!....Bd velocity
!
! print*,'ustar--',indnd(53)!,ustar(1:2, 738)/sqrt(ustar(1, 738)**2+ustar(2, 738)**2),munacu(1:2,36) ,snsigm(1:2, 36), munacn(36)
! print*,'ustar--',snsigm(1:2, 738)/sqrt(snsigm(1, 738)**2+snsigm(2, 738)**2),&
!                  munacu(1:2, 738)/sqrt(munacu(1, 738)**2+munacu(2, 738)**2)

!if(ntria.gt.0) call getvelo_mpt_marie(ustar,gelag,intfac,inpoel,coord,unkno,indnd)
if(nquad.gt.0) call getvelo_mpt_mariequad(ustar,gelagq,intfac,ipqua,coord,unkno,indnd)
! print*,'ustar--',ustar(1,ipqua(2:3,881)),ustar(2,ipqua(2:3,881))
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

!...Specified boundary nodal velocity

if(ncase.eq.13)then !...Saltzman
ustar(2,ipf(1:2)) = 0.d0
ustar(1,ipf(1:2)) = 1.d0
endif

endif
!
if(ncase.eq.1)then
!
!...impose exact solution along the Boundary
!
ustar(1, ipf(1)) = sin(pi*coord(1,ipf(1)))*cos(pi*coord(2,ipf(1)))
ustar(2, ipf(1)) =-cos(pi*coord(1,ipf(1)))*sin(pi*coord(2,ipf(1)))
!
ustar(1, ipf(2)) = sin(pi*coord(1,ipf(2)))*cos(pi*coord(2,ipf(2)))
ustar(2, ipf(2)) =-cos(pi*coord(1,ipf(2)))*sin(pi*coord(2,ipf(2)))
!
ustar(1, ipf(3)) = sin(pi*coord(1,ipf(3)))*cos(pi*coord(2,ipf(3)))
ustar(2, ipf(3)) =-cos(pi*coord(1,ipf(3)))*sin(pi*coord(2,ipf(3)))
!
endif
!
!if(ipf(1).eq.100) print*,'ipf',ipf(1),ustar(1:2, ipf(1))
!
900 enddo
!endif
enddo !iloop
!
!...Imposing the zero normal velocity for BC...
!
! call getbcvn_lag(bface, intfac, gflag, ustar)
!
!call getbcve_exact(bface, intfac, gflag, ustar, coord, itime)
!
!  print*,'ustar',ustar(1:2, 18),ustar(1:2, 20),ustar(1:2, 21)
!  if(ustar(2, 18).gt.0.d0) then
!   print*,'Wrong movng direction'!,munacu(2, 28) - snsigm(2, 28) , fpres(2, 28),munacn(28)
!   stop
!  endif
!ustar(:,1) = 0.d0
!
!...4.2: Update the Riemann forces at every node...
!
do 990 ifa = 1 , -nafac
  ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
  ustar(1:2, ipf(3)) = 0.5d0*(ustar(1:2, ipf(1)) + ustar(1:2, ipf(2)))
990 enddo
!
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
!if(ie==21) print*,'fire',iv,ie,snsigmlq(1,1:2,iv,ie),munaclq(1:2,iv,ie),ipq(:),munaulq(1,1:2,iv,ie)
!if(ie.eq.21) print*,'fstarq',iv,ie,fstarq(1,1:2,iv,ie)
!
enddo

!
do ifa =1, 2
!fstarq(1:2, ifa, 5, ie) = 0.5d0*(fstarq(1:2, 4, 1, ie) + fstarq(1:2, 3, 2, ie))
!fstarq(1:2, ifa, 6, ie) = 0.5d0*(fstarq(1:2, 4, 2, ie) + fstarq(1:2, 3, 3, ie))
!fstarq(1:2, ifa, 7, ie) = 0.5d0*(fstarq(1:2, 4, 3, ie) + fstarq(1:2, 3, 4, ie))
!fstarq(1:2, ifa, 8, ie) = 0.5d0*(fstarq(1:2, 4, 4, ie) + fstarq(1:2, 3, 1, ie))
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
end subroutine getndvelo_lag_mc_curvpseudo
!
!...Calculate the velocity at the middle point for curved mesh...
!
subroutine getvelo_mpt_curv(ustar,geoel,gelag,gelagq,intfac,iptri,ipqua,coord,coold,unkno,indnd, aflim, afvec, vlave, vnulq)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:3,1:ngelgq,1:nquad),      intent(in)::gelagq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
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
real*8::xv(nvtri), yv(nvtri), bt(ndegr, nvtri),btv(ndegr, nvtri)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, nvqua), bqv(ndegr, nvqua)
real*8::xpq(1:2,1:nvqua), xpqi(1:2,1:nvqua)
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
!...Reference triangle...
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
if(ncurv.eq.1)then
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
endif
!
!...Reference quad
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
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvtri
!...Basis function
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds
enddo
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
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...
!...zero out unknv
unknvt = 0.d0
!
do iv   = 4, nvtri
!
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
enddo
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvt(1, iv)
elseif(ndens.eq.3)then
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
murie = rhoct*sdctr! + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
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
dr = 1.0d0
ds = 1.0d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
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
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
enddo
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
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...
!...zero out unknv
unknvq = 0.d0
!
do iv   = 5, 8
!
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
!...Remove null mode...
!
unknvq(2:3, iv)=unknvq(2:3, iv)-vnulq(1:2,iv,ie)
!
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvq(1, iv)
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

!call getrhoig_quad(rhoi,r,s, xpqi(1:2,1:4))!
!call getdensity_quadllnl(r, s, xpq(1:2,1:4), xpqi(1:2,1:4), rhoi, rhon)
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
if(othog.lt.1.d-6)then
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
end subroutine getvelo_mpt_curv
!
!...Calculate the velocity at the middle point for curved mesh...
!
subroutine getvelo_mpt_curv2(ustar,geoel,gelag,gelagq,intfac,iptri,ipqua,coord,unkno,indnd, aflim, afvec, vlave, vnulq,munacn)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
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
integer::ifa,iel,ier,ie,idfal,idfar,ig,ideg,jdeg,ielem
integer::iv,ishp
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
!
integer,dimension(npoin)::icemp
!
real*8::eps
real*8::unknvt(1:nq, 1:nvtri)
real*8::unknvq(1:nq, 1:nvqua)
real*8::vnorm(3, nvqua)
real*8::xv(nvtri), yv(nvtri), bt(ndegr, nvtri),btv(ndegr, nvtri)
real*8::xvq(nvqua), yvq(nvqua), bq(ndegr, nvqua), bqv(ndegr, nvqua)
real*8::xpf(1:2, 1:nvfac)
real*8::unmpt(5,2,npoin)
real*8::ft(2, 2), fn(2, 2)
!real*8,allocatable:: ucurv(:, :)
real*8::dudr,duds,dvdr,dvds,pvtx
real*8::uvtxr,vvtxr,evtxr, pvtxr
real*8::uvtxl,vvtxl,evtxl, pvtxl,rhol,rhor,presl,presr,lenmc,mufal,mufar,rhovl,rhovr
real*8::deltu,murie
real*8::rhoct,uctr,vctr,ectr,pctr,sdctr
real*8::uvtx,vvtx,evtx,dux,duy
real*8::fnx,fny, rho,rhomc,rhovt,rhomv
real*8::dr, ds,rc,sc,othog,rcv,scv
real*8::dshpr(3),dxdr,dydr,djaco,dwav1,dwav2
real*8::ftx,fty,uhll,vhll,vdon,vtan
real*8::csoul, csour,pstar,vdonl,vdonr
!
!...For quadratic mesh, only nafac high-order nodes need be recalculated at most...
!
eps = 1.d-6
!
!...Reference triangle...
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
if(ncurv.eq.1)then
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
endif
!
!...Reference quad
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
dr = .5d0
ds = .5d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do iv =1 ,nvtri
!...Basis function
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds
enddo
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
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) ) !...sound speed at the center...
!
!...
!...zero out unknv
unknvt = 0.d0
!
do iv   = 4, nvtri
!
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
enddo
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvt(1, iv)
elseif(ndens.eq.3)then
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
unmpt(5, icemp(ipt(iv)), ipt(iv)) = rhovt
!
!if(ipt(iv).eq.71)print*,'midppp',unknvt(2:3, iv),ie,unkno(2:3,2,ielem),unkno(2:3,3,ielem)
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
dr = 1.0d0
ds = 1.0d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
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
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
!...Remove null mode...
!
unknvq(2:3, iv)=unknvq(2:3, iv)-vnulq(1:2,iv,ie)
!
if(ndens.eq.1)then
rhovt  = 1.d0/unknvq(1, iv)
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
!uvtx = unkno(1,2,ielem)  + dudr*bq(2, iv) + duds*bq(3, iv)
!vvtx = unkno(1,3,ielem)  + dvdr*bq(2, iv) + dvds*bq(3, iv)
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
icemp(ipq(iv)) = icemp(ipq(iv)) + 1
dux= vlave(1, ipq(iv))-unknvq(2, iv)
duy= vlave(2, ipq(iv))-unknvq(3, iv)
deltu = sqrt(dux**2 + duy**2)
murie = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
unmpt(1, icemp(ipq(iv)), ipq(iv)) = max(eps, (gamlg-1.d0)*rhovt*(evtx - 0.5d0*(uvtx**2 + vvtx**2)))!pvtx
unmpt(2, icemp(ipq(iv)), ipq(iv)) = murie
unmpt(3, icemp(ipq(iv)), ipq(iv)) = unknvq(2, iv)
unmpt(4, icemp(ipq(iv)), ipq(iv)) = unknvq(3, iv)
unmpt(5, icemp(ipq(iv)), ipq(iv)) = rhovt
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
!if(othog.lt.1.d6)then
if(munacn(1,1,ipf(3)).lt.1.d6)then
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
ustar(1, ipf(3)) = uvtxl !+ presl/mufal*fnx
ustar(2, ipf(3)) = vvtxl !+ presl/mufal*fny
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
rhovl = unmpt(5, 1, ipf(3))
csoul = sqrt(gamlg*presl/rhovl)
!
presr = unmpt(1, 2, ipf(3))
mufar = unmpt(2, 2, ipf(3))
uvtxr = unmpt(3, 2, ipf(3))
vvtxr = unmpt(4, 2, ipf(3))
rhovr = unmpt(5, 2, ipf(3))
csour = sqrt(gamlg*presr/rhovr)
!
fnx = fn(1, 1) !...face normal vector
fny = fn(2, 1)
!
ftx = -fny
fty = fnx
!
!uhll =  (sqrt(rhovl)*uvtxl + sqrt(rhovr)*uvtxr)/(sqrt(rhovl) + sqrt(rhovr))
!vhll =  (sqrt(rhovl)*vvtxl + sqrt(rhovr)*vvtxr)/(sqrt(rhovl) + sqrt(rhovr))
!
!vdon = uhll*fnx + vhll*fny
!
!vdonl = uvtxl*dwav1 + vvtxl*dwav2
!vdonr = uvtxr*dwav1 + vvtxr*dwav2
!
!call getvelo_mpt_dukowicz(vdonl,vdonr,rhovl,rhovr,presl,presr,csoul,csour,gamlg,pstar,vdon)
!
!vtan = 0.5d0*(uvtxl*ftx+vvtxl*fty + uvtxr*ftx + vvtxr*fty)
!...
!...Mar
!
if(indnd(ipf(3)).eq.0)then
ustar(1, ipf(3)) = (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fnx
ustar(2, ipf(3)) = (mufal*vvtxl + mufar*vvtxr)/(mufal+mufar) - (presr- presl)/(mufal+mufar)*fny
!
!    if(ipf(3).eq.71) print*,'midlle velocity',ustar(1:2, ipf(3)),ipf(3),ustar(1:2,ipf(1)),ustar(1:2,ipf(2)),&
!                              (mufal*uvtxl + mufar*uvtxr)/(mufal+mufar),mufal,uvtxl , mufar,uvtxr,mufal,mufar
endif
!
endif
endif
!
450 enddo
!
end subroutine getvelo_mpt_curv2
!
!...Get the null mode of velocity in curved quad...
!
subroutine getnullmd_quadc(ipqua, geoel, vnulq, ustar, gqdmp)
use constant
implicit none
integer,dimension(1:nvqua,1:nquad),         intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(out):: vnulq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,  dimension(1:nquad),                 intent(in)::gqdmp
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:6, 1:nvqua)::bq
real*8,dimension(1:5, 1:5)::amnul
real*8,dimension(1:5, 1:2)::dunul
real*8,dimension(1:ndimn, 1:nvqua)::upqua
real*8::ugrad(5, 2)
integer,dimension(1:nvqua) :: ipq
!
integer::iv, ie, ielem, ij, ik
real*8::ama, amb, amc, amd, amjac
real*8::dr, ds, rcq, scq
real*8::deluq(1:ndimn, 1:nvqua)
real*8::x5(5,5), mmatr(5,5), b55(5)
!
!...shape functions
!
dr = 1.0d0
ds = 1.0d0
!
rcq = 0.d0
scq = 0.d0
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
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rcq)/dr
bq(3, iv) = (yvq(iv)-scq)/ds
bq(4, iv) = bq(2, iv)**2
bq(5, iv) = bq(3, iv)**2
bq(6, iv) = bq(2, iv)*bq(3, iv)
enddo
!
!...zero vnulq
!
vnulq = 0.d0
amnul = 0.d0
!
do iv=1, 8
!
amnul(1, 1) = amnul(1, 1) + bq(2, iv)*bq(2, iv)
amnul(1, 2) = amnul(1, 2) + bq(2, iv)*bq(3, iv)
amnul(1, 3) = amnul(1, 3) + bq(2, iv)*bq(4, iv)
amnul(1, 4) = amnul(1, 4) + bq(2, iv)*bq(5, iv)
amnul(1, 5) = amnul(1, 5) + bq(2, iv)*bq(6, iv)

amnul(2, 1) = amnul(2, 1) + bq(3, iv)*bq(2, iv)
amnul(2, 2) = amnul(2, 2) + bq(3, iv)*bq(3, iv)
amnul(2, 3) = amnul(2, 3) + bq(3, iv)*bq(4, iv)
amnul(2, 4) = amnul(2, 4) + bq(3, iv)*bq(5, iv)
amnul(2, 5) = amnul(2, 5) + bq(3, iv)*bq(6, iv)

amnul(3, 1) = amnul(3, 1) + bq(4, iv)*bq(2, iv)
amnul(3, 2) = amnul(3, 2) + bq(4, iv)*bq(3, iv)
amnul(3, 3) = amnul(3, 3) + bq(4, iv)*bq(4, iv)
amnul(3, 4) = amnul(3, 4) + bq(4, iv)*bq(5, iv)
amnul(3, 5) = amnul(3, 5) + bq(4, iv)*bq(6, iv)

amnul(4, 1) = amnul(4, 1) + bq(5, iv)*bq(2, iv)
amnul(4, 2) = amnul(4, 2) + bq(5, iv)*bq(3, iv)
amnul(4, 3) = amnul(4, 3) + bq(5, iv)*bq(4, iv)
amnul(4, 4) = amnul(4, 4) + bq(5, iv)*bq(5, iv)
amnul(4, 5) = amnul(4, 5) + bq(5, iv)*bq(6, iv)

amnul(5, 1) = amnul(5, 1) + bq(6, iv)*bq(2, iv)
amnul(5, 2) = amnul(5, 2) + bq(6, iv)*bq(3, iv)
amnul(5, 3) = amnul(5, 3) + bq(6, iv)*bq(4, iv)
amnul(5, 4) = amnul(5, 4) + bq(6, iv)*bq(5, iv)
amnul(5, 5) = amnul(5, 5) + bq(6, iv)*bq(6, iv)

enddo
!
!...Inverse matrix
!
call getinvmat(5,amnul, x5, b55)
!
!
!print*,'amunul',amnul
!
!
do 350 ie = 1,nquad !...(1)ie = 1,nelem
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
upqua(1, 1:8) = ustar(1, ipq(1:8))
upqua(2, 1:8) = ustar(2, ipq(1:8))
!
upqua(1, 9) = -0.25d0*(upqua(1, 1) + upqua(1, 2) + upqua(1, 3) + upqua(1, 4)) +&
0.5d0*(upqua(1, 5) + upqua(1, 6) + upqua(1, 7) + upqua(1, 8))
upqua(2, 9) = -0.25d0*(upqua(2, 1) + upqua(2, 2) + upqua(2, 3) + upqua(2, 4)) +&
0.5d0*(upqua(2, 5) + upqua(2, 6) + upqua(2, 7) + upqua(2, 8))
!
do iv = 1, 8
deluq(1:2, iv) = upqua(1:2, iv) - upqua(1:2, 9)
enddo
!
dunul = 0.d0
!
do iv=1, 8
!
dunul(1, 1)  = dunul(1, 1) + bq(2, iv)*deluq(1, iv)
dunul(2, 1)  = dunul(2, 1) + bq(3, iv)*deluq(1, iv)
dunul(3, 1)  = dunul(3, 1) + bq(4, iv)*deluq(1, iv)
dunul(4, 1)  = dunul(4, 1) + bq(5, iv)*deluq(1, iv)
dunul(5, 1)  = dunul(5, 1) + bq(6, iv)*deluq(1, iv)
!
dunul(1, 2)  = dunul(1, 2) + bq(2, iv)*deluq(2, iv)
dunul(2, 2)  = dunul(2, 2) + bq(3, iv)*deluq(2, iv)
dunul(3, 2)  = dunul(3, 2) + bq(4, iv)*deluq(2, iv)
dunul(4, 2)  = dunul(4, 2) + bq(5, iv)*deluq(2, iv)
dunul(5, 2)  = dunul(5, 2) + bq(6, iv)*deluq(2, iv)
enddo
!
ugrad = 0.d0
!
do ij = 1, 5
do ik = 1, 5
ugrad(ij, 1) = ugrad(ij, 1) + x5(ik, ij)*dunul(ik ,1)
ugrad(ij, 2) = ugrad(ij, 2) + x5(ik, ij)*dunul(ik ,2)
enddo
enddo
!
!ucell = 0.d0
!
do iv = 1, 8
upqua(1, iv) = upqua(1, 9) + ugrad(1, 1)*bq(2, iv) + ugrad(2, 1)*bq(3, iv)  !+&
!ugrad(3, 1)*bq(4, iv)   + ugrad(4, 1)*bq(5, iv)   + ugrad(5, 1)*bq(6, iv)
upqua(2, iv) = upqua(2, 9) + ugrad(1, 2)*bq(2, iv) + ugrad(2, 2)*bq(3, iv) ! +&
!ugrad(3, 2)*bq(4, iv)   + ugrad(4, 2)*bq(5, iv)   + ugrad(5, 2)*bq(6, iv)
enddo
!
vnulq(1, 1:8, ie) = ustar(1, ipq(1:8)) - upqua(1, 1:8)
vnulq(2, 1:8, ie) = ustar(2, ipq(1:8)) - upqua(2, 1:8)
!
!...Imposing the dampening coefficient
!
vnulq(1:2, 1:8, ie) = vnulq(1:2, 1:8, ie)*gqdmp(ie)
!
350 enddo
!
!...get null mode...
!
end subroutine getnullmd_quadc
!
!...Get the dampening coeffcients...
!
subroutine getnullve_gqdmpquadc(ipqua, geoel, ustar, unkno, gqdmp)
use constant
implicit none
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,  dimension(1:nquad),        intent(out)::gqdmp
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:3, 1:nvqua)::bq
integer,dimension(1:nvqua) :: ipq
real*8::gdampq(1:nvqua)
real*8::unknvq(1:nq, 1:nvqua)
real*8::gdux1,gdux2,gduy1,gduy2,dr,ds,rcq,scq
!
integer::iv, ie, ielem,ideg
!
!...shape functions
!
dr = 1.0d0
ds = 1.0d0
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
do ie = 1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
rcq= geoel(1, ielem) !...mass center...
scq= geoel(2, ielem)
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rcq)/dr
bq(3, iv) = (yvq(iv)-scq)/ds
enddo

!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1, nvqua
!
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
enddo
!
gdux1 = ustar(1,ipq(iv))-unknvq(2, iv)
gduy1 = ustar(2,ipq(iv))-unknvq(3, iv)
!
gdux2 = ustar(1,ipq(iv))-unkno(1, 2, ielem)
gduy2 = ustar(2,ipq(iv))-unkno(1, 3, ielem)

gdampq(iv) = min(1.d0, .5d0, 1.0d0*sqrt(gdux1**2+gduy1**2)/sqrt(gdux2**2+gduy2**2))
!
!if(ie==2192) print*,'gdmp',sqrt(gdux1**2+gduy1**2)/sqrt(gdux2**2+gduy2**2),iv
!
enddo
!
gqdmp(ie) = 0.5d0!minval(gdampq(1:8))
!
!if(ie==2192) print*,'gdmp2',gqdmp(ie)
!
enddo
!
!gqdmp =1.d0
!
!
!...get null mode...
!
end subroutine getnullve_gqdmpquadc
!
!...Mpt velocity fitting...
!
subroutine getvelofit_mpt3(intfac, ipqua, geoel, ustar, vnulq, coord)
use constant
implicit none
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,dimension(1:nvqua,1:nquad),         intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord !...nodal velocity
real*8,dimension(1:ndimn,1:nvqua,1:nquad),           intent(out)::vnulq
!
integer,parameter::neqls=6
real*8,dimension(1:ndimn, 1:neqls) :: xpq
real*8,dimension(1:3, 1:neqls)::bq
real*8,dimension(1:3, 1:3)::amnul
real*8,dimension(1:3, 1:2)::dunul
real*8::ugrad(3, 2)
integer::ipm(neqls), ipf(3)
integer:: ipq(1:nvqua)
!
integer::iv, ie, ielem, ij, ik,iel,ier,ifa,nsten
real*8::amnua, amnub, amnuc, amnud, amnuj
real*8::dxm,dym,xcf,ycf
real*8::deluq(1:ndimn, 1:neqls)
real*8::x6(3,3),b6(3)
real*8,dimension(1:ndimn,1:npoin)::usfit

!
!...Initilize usfit
!
usfit  = ustar
!
do ifa = nbfac+1, nafac
!
iel = intfac(1, ifa)
ier = intfac(2, ifa)
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
!...Get the stencile...
!
nsten = 2

ipm(1:2) = ipf(1:2)
!
do iv = 1, 4
!
if(ipf(1).ne.ipqua(iv, iel).and.ipf(2).ne.ipqua(iv, iel).and.ipf(3).ne.ipqua(iv, iel))then
nsten = nsten + 1
ipm(nsten) = ipqua(iv, iel)
endif
!
if(ipf(1).ne.ipqua(iv, ier).and.ipf(2).ne.ipqua(iv, ier).and.ipf(3).ne.ipqua(iv, ier))then
nsten = nsten + 1
ipm(nsten) = ipqua(iv, ier)
endif
enddo
!
!...Get basis function...
!
xpq(1, 1:neqls) = coord(1, ipm(1:neqls))
xpq(2, 1:neqls) = coord(2, ipm(1:neqls))
!
dxm = maxval(xpq(1, 1:neqls)) - minval(xpq(1, 1:neqls))
dym = maxval(xpq(2, 1:neqls)) - minval(xpq(2, 1:neqls))

!...Face center

xcf = coord(1, ipf(3))
ycf = coord(2, ipf(3))
!

do iv =1 ,neqls
bq(1, iv) = 1.d0
bq(2, iv) = (xpq(1, iv)-xcf)/dxm
bq(3, iv) = (xpq(2, iv)-ycf)/dym
enddo
!
!...Get LHS matrix...
!
amnul = 0.d0
!
do iv=1,neqls
!
amnul(1, 1) = amnul(1, 1) + bq(1, iv)*bq(1, iv)
amnul(1, 2) = amnul(1, 2) + bq(1, iv)*bq(2, iv)
amnul(1, 3) = amnul(1, 3) + bq(1, iv)*bq(3, iv)

amnul(2, 1) = amnul(2, 1) + bq(2, iv)*bq(1, iv)
amnul(2, 2) = amnul(2, 2) + bq(2, iv)*bq(2, iv)
amnul(2, 3) = amnul(2, 3) + bq(2, iv)*bq(3, iv)

amnul(3, 1) = amnul(3, 1) + bq(3, iv)*bq(1, iv)
amnul(3, 2) = amnul(3, 2) + bq(3, iv)*bq(2, iv)
amnul(3, 3) = amnul(3, 3) + bq(3, iv)*bq(3, iv)
!if(ipf(3).eq.996)print*,'usfit2',bq(1:6, iv)
enddo
!
!if(ipf(3).eq.996)print*,'usfit3',amnul
!
!...Inverse matrix
!
call getinvmat(3,amnul, x6, b6)
!
!...Get RHS....
!
do iv = 1, neqls
deluq(1:2, iv) = ustar(1:2, ipm(iv))
enddo
!
dunul = 0.d0
!
do iv=1, neqls
!
dunul(1, 1)  = dunul(1, 1) + bq(1, iv)*deluq(1, iv)
dunul(2, 1)  = dunul(2, 1) + bq(2, iv)*deluq(1, iv)
dunul(3, 1)  = dunul(3, 1) + bq(3, iv)*deluq(1, iv)
!
dunul(1, 2)  = dunul(1, 2) + bq(1, iv)*deluq(2, iv)
dunul(2, 2)  = dunul(2, 2) + bq(2, iv)*deluq(2, iv)
dunul(3, 2)  = dunul(3, 2) + bq(3, iv)*deluq(2, iv)

enddo
!
ugrad = 0.d0
!
do ij = 1, 3
do ik = 1, 3
ugrad(ij, 1) = ugrad(ij, 1) + x6(ik, ij)*dunul(ik ,1)
ugrad(ij, 2) = ugrad(ij, 2) + x6(ik, ij)*dunul(ik ,2)
enddo
enddo
!
!...Update the ustar for middle point...
!
usfit(1:2, ipf(3)) =  ugrad(1, 1:2)
!
!if(ipf(3).eq.996)print*,'usfit',x6(1:6,4)
!
enddo
!
!
do ie = 1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua, ie)
!
vnulq(1, 1:nvqua, ie) = vnulq(1, 1:nvqua, ie) + ustar(1, ipq(1:nvqua)) - usfit(1, ipq(1:nvqua))
vnulq(2, 1:nvqua, ie) = vnulq(2, 1:nvqua, ie) + ustar(2, ipq(1:nvqua)) - usfit(2, ipq(1:nvqua))

enddo

end subroutine getvelofit_mpt3!
!
!...Get the dampening coeffcients of null mode velocity on curved quad...
!
subroutine getcoef_nuvequadc(ipqua, geoel, ustar, unkno, gqdmp)
use constant
implicit none
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,  dimension(1:nquad),        intent(out)::gqdmp
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:3, 1:nvqua)::bq
integer,dimension(1:nvqua) :: ipq
real*8::gdampq(1:nvqua)
real*8::unknvq(1:nq, 1:nvqua)
real*8::gdux1,gdux2,gduy1,gduy2,dr,ds,rcq,scq
real*8::dmplg,smthid,du1,du2
!
integer::iv, ie, ielem,ideg
!
!...shape functions
!
dr = 1.0d0
ds = 1.0d0
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
do ie = 1, nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
rcq= geoel(1, ielem) !...mass center...
scq= geoel(2, ielem)
!
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rcq)/dr
bq(3, iv) = (yvq(iv)-scq)/ds
enddo
!
!...zero out unknv
!
unknvq = 0.d0
!
do iv   = 1, nvqua
!
do ideg = 1,mdegr
unknvq(1:nq, iv) = unknvq(1:nq, iv) + unkno(ideg,1:nq,ielem)*bq(ideg, iv)
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

gdampq(iv) = smthid
!
!if(ie==2192) print*,'gdmp',sqrt(gdux1**2+gduy1**2)/sqrt(gdux2**2+gduy2**2),iv
!
enddo
!
dmplg = maxval(gdampq(1:8))
gqdmp(ie) = min(0.5d0,dmplg)
!
!if(ie==2192) print*,'gdmp2',gqdmp(ie)
!
enddo
!
!gqdmp =1.d0
!
!
!...get null mode...
!
end subroutine getcoef_nuvequadc
!
!...Get the dampening coeffcients of null mode velocity on curved tria...
!
subroutine getcoef_nuvetriac(iptri, geoel, ustar, unkno, gtdmp)
use constant
implicit none
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,  dimension(1:ntria),        intent(out)::gtdmp
real*8,dimension(1:nvtri):: xv, yv
real*8,dimension(1:3, 1:nvtri)::bt
integer,dimension(1:nvtri) :: ipt
real*8::gdampt(1:nvtri)
real*8::unknvt(1:nq, 1:nvtri)
real*8::gdux1,gdux2,gduy1,gduy2,dr,ds,rct,sct
real*8::dmplg,smthid,du1,du2
!
integer::iv, ie, ielem,ideg
!
!...shape functions
!
dr = 0.5d0
ds = 0.5d0
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
do ie = 1, ntria
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
rct= geoel(1, ielem) !...mass center...
sct= geoel(2, ielem)
!
do iv =1 ,nvtri
!...Left cell + intfac(3,ifa)
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rct)/dr
bt(3, iv) = (yv(iv)-sct)/ds
enddo
!
!...zero out unknv
!
unknvt = 0.d0
!
do iv   = 1, nvtri
!
do ideg = 1,mdegr
unknvt(1:nq, iv) = unknvt(1:nq, iv) + unkno(ideg,1:nq,ielem)*bt(ideg, iv)
enddo
!
gdux1 = ustar(1,ipt(iv))-unknvt(2, iv)
gduy1 = ustar(2,ipt(iv))-unknvt(3, iv)
!
gdux2 = ustar(1,ipt(iv))-unkno(1, 2, ielem)
gduy2 = ustar(2,ipt(iv))-unkno(1, 3, ielem)
!
du1 = sqrt(gdux1**2+gduy1**2)
du2 = sqrt(gdux2**2+gduy2**2)
!
if(du2.le.1.d-6)then
smthid = 0.d0
else
smthid = du1/(du2+1.d-6)
endif

gdampt(iv) = smthid
!
!if(ie==2192) print*,'gdmp',sqrt(gdux1**2+gduy1**2)/sqrt(gdux2**2+gduy2**2),iv
!
enddo
!
dmplg = maxval(gdampt(1:6))
gtdmp(ie) = min(0.5d0,dmplg)
!
!if(ie==2192) print*,'gdmp2',gtdmp(:)
!
enddo
!
!gqdmp =1.d0
!
!
!...get null mode...
!
end subroutine getcoef_nuvetriac

!
!...Get the null mode of velocity on reference domain without averaged velocity given...
!
subroutine getnullmd_quadl(ipqua, geoel, vnulq, ustar, gqdmp)
use constant
implicit none
integer,dimension(1:nvqua,1:nquad),         intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(out):: vnulq
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::ustar !...nodal velocity
real*8,  dimension(1:nquad),                 intent(in)::gqdmp
!
integer, parameter::ncbf = 2
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:3, 1:nvqua)::bq
real*8,dimension(1:2, 1:2)::amnul
real*8,dimension(1:2, 1:2)::dunul
real*8,dimension(1:ndimn, 1:nvqua)::upqua
real*8::ugrad(2, 2)
integer,dimension(1:nvqua) :: ipq
integer,dimension(npoin)::indnd
!
integer::iv, ie, ielem, ij, ik
real*8::amnua, amnub, amnuc, amnud, amnuj
real*8::dr, ds, rcq, scq
real*8::deluq(1:ndimn, 1:nvqua)
real*8::x5(2,2), mmatr(2,2), b55(2)
!
!...shape functions
!
dr = 1.0d0
ds = 1.0d0
!
rcq = 0.d0
scq = 0.d0
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
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rcq)/dr
bq(3, iv) = (yvq(iv)-scq)/ds
enddo
!
!...zero vnulq
!
vnulq = 0.d0
amnul = 0.d0
!
do iv=1, 8
!
amnul(1, 1) = amnul(1, 1) + bq(2, iv)*bq(2, iv)
amnul(1, 2) = amnul(1, 2) + bq(2, iv)*bq(3, iv)

amnul(2, 1) = amnul(2, 1) + bq(3, iv)*bq(2, iv)
amnul(2, 2) = amnul(2, 2) + bq(3, iv)*bq(3, iv)

enddo
!
!
!...Inverse matrix
!
!call getinvmat(2,amnul, x5, b55)
!
amnua = amnul(1, 1)
amnub = amnul(1, 2)
amnuc = amnul(2, 1)
amnud = amnul(2, 2)
!
amnuj = amnua*amnud-amnuc*amnub
!
x5(1, 1) = amnud/amnuj
x5(1, 2) =-amnub/amnuj
x5(2, 1) = x5(1, 2)
x5(2, 2) = amnua/amnuj
!
!print*,'amunul',x5
!
!
do 350 ie = 1,nquad !...(1)ie = 1,nelem
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
upqua(1, 1:8) = ustar(1, ipq(1:8))
upqua(2, 1:8) = ustar(2, ipq(1:8))
!
!if(ie.eq.1950) print*,'i1950',upqua(1,1:8),upqua(2,1:8)
!
upqua(1, 9) = -0.25d0*(upqua(1, 1) + upqua(1, 2) + upqua(1, 3) + upqua(1, 4)) +&
0.50d0*(upqua(1, 5) + upqua(1, 6) + upqua(1, 7) + upqua(1, 8))
upqua(2, 9) = -0.25d0*(upqua(2, 1) + upqua(2, 2) + upqua(2, 3) + upqua(2, 4)) +&
0.5d0*(upqua(2, 5) + upqua(2, 6) + upqua(2, 7) + upqua(2, 8))
!
do iv = 1, 8
deluq(1:2, iv) = upqua(1:2, iv) - upqua(1:2, 9)
enddo
!
!if(ie.eq.1950) print*,'i19502',deluq(1, 1:8)
!
dunul = 0.d0
!
do iv=1, 8
!
dunul(1, 1)  = dunul(1, 1) + bq(2, iv)*deluq(1, iv)
dunul(2, 1)  = dunul(2, 1) + bq(3, iv)*deluq(1, iv)
!
dunul(1, 2)  = dunul(1, 2) + bq(2, iv)*deluq(2, iv)
dunul(2, 2)  = dunul(2, 2) + bq(3, iv)*deluq(2, iv)
enddo
!
!if(ie.eq.1950) print*,'i1950g',dunul(1:2,1)
!
ugrad = 0.d0
!
do ij = 1, 2
do ik = 1, 2
ugrad(ij, 1) = ugrad(ij, 1) + x5(ik, ij)*dunul(ik ,1)
ugrad(ij, 2) = ugrad(ij, 2) + x5(ik, ij)*dunul(ik ,2)
enddo
enddo
!
!if(ie.eq.1950) print*,'i1950g',ugrad(1:2,1),ugrad(1:2,2),upqua(1:2,9)
!
!ucell = 0.d0
!
do iv = 1, 8
upqua(1, iv) = upqua(1, 9) + ugrad(1, 1)*bq(2, iv) + ugrad(2, 1)*bq(3, iv)
upqua(2, iv) = upqua(2, 9) + ugrad(1, 2)*bq(2, iv) + ugrad(2, 2)*bq(3, iv)
enddo
!
vnulq(1, 1:8, ie) = ustar(1, ipq(1:8)) - upqua(1, 1:8)
vnulq(2, 1:8, ie) = ustar(2, ipq(1:8)) - upqua(2, 1:8)
!
!if(ie.eq.1950) print*,'delu', upqua(1,1:8),ustar(1,ipq(1:8))
!
!...Imposing the dampening coefficient
!
vnulq(1:2, 1:8, ie) = vnulq(1:2, 1:8, ie)*gqdmp(ie)

!
350 enddo
!
!...get null mode...
!
end subroutine getnullmd_quadl
!
!...Get the null mode of velocity on physical domain with 8 vertices velocity fitting...
!
subroutine getnullmd_quadlp(ipqua, geoel, vnulq, ustar, coord, gqdmp)
use constant
implicit none
integer,dimension(1:nvqua,1:nquad),         intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(inout):: vnulq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord !...nodal velocity

real*8,  dimension(1:nquad),                 intent(in)::gqdmp
!
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:3, 1:nvqua)::bq
real*8,dimension(1:2, 1:2)::amnul
real*8,dimension(1:2, 1:2)::dunul
real*8,dimension(1:ndimn, 1:nvqua)::upqua
real*8::ugrad(2, 2)
integer,dimension(1:nvqua) :: ipq
integer,dimension(npoin)::indnd
!
integer::iv, ie, ielem, ij, ik
real*8::amnua, amnub, amnuc, amnud, amnuj
real*8::dxm,dym,xcv,ycv
real*8::deluq(1:ndimn, 1:nvqua)
real*8::x5(2,2)
!
do 450 ie = 1,nquad !...(1)ie = 1,nelem
!
amnul = 0.d0
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
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
dxm = maxval(xpq(1, 1:npqua)) - minval(xpq(1, 1:npqua))
dym = maxval(xpq(2, 1:npqua)) - minval(xpq(2, 1:npqua))
!
xcv = xpq(1,9)
ycv = xpq(2,9)
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xpq(1, iv)-xcv)/dxm
bq(3, iv) = (xpq(2, iv)-ycv)/dym
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
!...Part 2
!
upqua(1, 1:8) = ustar(1, ipq(1:8))
upqua(2, 1:8) = ustar(2, ipq(1:8))
!
!if(ie.eq.1950) print*,'i1950',upqua(1,1:8),upqua(2,1:8)
!
upqua(1, 9) = -0.25d0*(upqua(1, 1) + upqua(1, 2) + upqua(1, 3) + upqua(1, 4)) +&
0.50d0*(upqua(1, 5) + upqua(1, 6) + upqua(1, 7) + upqua(1, 8))
upqua(2, 9) = -0.25d0*(upqua(2, 1) + upqua(2, 2) + upqua(2, 3) + upqua(2, 4)) +&
0.5d0*(upqua(2, 5) + upqua(2, 6) + upqua(2, 7) + upqua(2, 8))
!
!upqua(1, 9) = 0.125d0*(upqua(1, 1) + upqua(1, 2) + upqua(1, 3) + upqua(1, 4)) +&
!0.1250d0*(upqua(1, 5) + upqua(1, 6) + upqua(1, 7) + upqua(1, 8))
!upqua(2, 9) = 0.125d0*(upqua(2, 1) + upqua(2, 2) + upqua(2, 3) + upqua(2, 4)) +&
!0.125d0*(upqua(2, 5) + upqua(2, 6) + upqua(2, 7) + upqua(2, 8))
!
do iv = 1, 8
deluq(1:2, iv) = upqua(1:2, iv) - upqua(1:2, 9)
enddo
!
dunul = 0.d0
!
do iv=1, 8
!
dunul(1, 1)  = dunul(1, 1) + bq(2, iv)*deluq(1, iv)
dunul(2, 1)  = dunul(2, 1) + bq(3, iv)*deluq(1, iv)
!
dunul(1, 2)  = dunul(1, 2) + bq(2, iv)*deluq(2, iv)
dunul(2, 2)  = dunul(2, 2) + bq(3, iv)*deluq(2, iv)
enddo
!
!if(ie.eq.1831) print*,'i1950g',deluq(1, :)
!
ugrad(1, 1) = x5(1, 1)*dunul(1 ,1) + x5(2, 1)*dunul(2 ,1)
ugrad(2, 1) = x5(1, 2)*dunul(1 ,1) + x5(2, 2)*dunul(2 ,1)

ugrad(1, 2) = x5(1, 1)*dunul(1 ,2) + x5(2, 1)*dunul(2 ,2)
ugrad(2, 2) = x5(1, 2)*dunul(1 ,2) + x5(2, 2)*dunul(2 ,2)
!
!if(ie.eq.1831) print*,'i1950g',ugrad(1:2,1)
!
!ucell = 0.d0
!
do iv = 1, 8
upqua(1, iv) = upqua(1, 9) + ugrad(1, 1)*bq(2, iv) + ugrad(2, 1)*bq(3, iv)
upqua(2, iv) = upqua(2, 9) + ugrad(1, 2)*bq(2, iv) + ugrad(2, 2)*bq(3, iv)
enddo
!
!vnulq(1, 1:8, ie) = ustar(1, ipq(1:8)) - upqua(1, 1:8)
!vnulq(2, 1:8, ie) = ustar(2, ipq(1:8)) - upqua(2, 1:8)
!
vnulq(1, 1:8, ie) = vnulq(1, 1:8, ie) + (ustar(1, ipq(1:8)) - upqua(1, 1:8))*gqdmp(ie)
vnulq(2, 1:8, ie) = vnulq(2, 1:8, ie) + (ustar(2, ipq(1:8)) - upqua(2, 1:8))*gqdmp(ie)
!
!if(ie.eq.1950) print*,'delu', upqua(1,1:8),ustar(1,ipq(1:8))
!
!...Imposing the dampening coefficient
!
!vnulq(1:2, 1:8, ie) = vnulq(1:2, 1:8, ie)*gqdmp(ie)

!
450 enddo
!
!...get null mode...
!
end subroutine getnullmd_quadlp
!
!...Get the null mode of velocity on physical domain with 6 vertices velocity fitting...
!
subroutine getnullmd_trialp(iptri, geoel, vnult, ustar, coord, gtdmp)
use constant
implicit none
integer,dimension(1:nvtri,1:ntria),         intent(in):: iptri
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:nvtri,1:ntria),   intent(inout):: vnult
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord !...nodal velocity

real*8,  dimension(1:ntria),                 intent(in)::gtdmp
!
real*8,dimension(1:ndimn, 1:nptri) :: xpt
real*8,dimension(1:3, 1:nvtri)::bt
real*8,dimension(1:2, 1:2)::amnul
real*8,dimension(1:2, 1:2)::dunul
real*8,dimension(1:ndimn, 1:nvtri)::uptri
real*8::ugrad(2, 2)
real*8::uptrc(2)
integer,dimension(1:nvtri) :: ipt
integer,dimension(npoin)::indnd
!
integer::iv, ie, ielem, ij, ik
real*8::amnua, amnub, amnuc, amnud, amnuj
real*8::dxm,dym,xcv,ycv
real*8::delut(1:ndimn, 1:nvtri)
real*8::x5(2,2)
!
do 450 ie = 1,ntria !...(1)ie = 1,nelem
!
amnul = 0.d0
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
ielem = ie
!
xpt(1, 1:nptri) = coord(1,iptri(1:nptri, ie))
xpt(2, 1:nptri) = coord(2,iptri(1:nptri, ie))
!
dxm = maxval(xpt(1, 1:nptri)) - minval(xpt(1, 1:nptri))
dym = maxval(xpt(2, 1:nptri)) - minval(xpt(2, 1:nptri))
!
!xcv = 1.d0/3.d0*(xpt(1, 1)+xpt(1, 2)+xpt(1, 3))
!ycv = 1.d0/3.d0*(xpt(2, 1)+xpt(2, 2)+xpt(2, 3))

xcv = -1.d0/9.d0*(xpt(1, 1)+xpt(1, 2)+xpt(1, 3)) + 4.d0/9.d0*(xpt(1, 4)+xpt(1, 5)+xpt(1, 6))
ycv = -1.d0/9.d0*(xpt(2, 1)+xpt(2, 2)+xpt(2, 3)) + 4.d0/9.d0*(xpt(2, 4)+xpt(2, 5)+xpt(2, 6))
!
do iv =1 ,nvtri
bt(1, iv) = 1.d0
bt(2, iv) = (xpt(1, iv)-xcv)/dxm
bt(3, iv) = (xpt(2, iv)-ycv)/dym
enddo
!
do iv=1, 6
!
amnul(1, 1 ) = amnul(1, 1 ) + bt(2, iv)*bt(2, iv)
amnul(1, 2 ) = amnul(1, 2 ) + bt(2, iv)*bt(3, iv)

amnul(2, 1 ) = amnul(2, 1 ) + bt(3, iv)*bt(2, iv)
amnul(2, 2 ) = amnul(2, 2 ) + bt(3, iv)*bt(3, iv)
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
!...Part 2
!
uptri(1, 1:6) = ustar(1, ipt(1:6))
uptri(2, 1:6) = ustar(2, ipt(1:6))
!
!if(ie.eq.1950) print*,'i1950',upqua(1,1:8),upqua(2,1:8)
!
uptrc(1) = -1.d0/9.d0*(uptri(1, 1) + uptri(1, 2) + uptri(1, 3)) +&
4.d0/9.d0*(uptri(1, 4) + uptri(1, 5) + uptri(1, 6))
uptrc(2) = -1.d0/9.d0*(uptri(2, 1) + uptri(2, 2) + uptri(2, 3)) +&
4.d0/9.d0*(uptri(2, 4) + uptri(2, 5) + uptri(2, 6))
!
!upqua(1, 9) = 0.125d0*(upqua(1, 1) + upqua(1, 2) + upqua(1, 3) + upqua(1, 4)) +&
!0.1250d0*(upqua(1, 5) + upqua(1, 6) + upqua(1, 7) + upqua(1, 8))
!upqua(2, 9) = 0.125d0*(upqua(2, 1) + upqua(2, 2) + upqua(2, 3) + upqua(2, 4)) +&
!0.125d0*(upqua(2, 5) + upqua(2, 6) + upqua(2, 7) + upqua(2, 8))
!
do iv = 1, 6
delut(1:2, iv) = uptri(1:2, iv) - uptrc(1:2)
enddo
!
dunul = 0.d0
!
do iv=1, 6
!
dunul(1, 1)  = dunul(1, 1) + bt(2, iv)*delut(1, iv)
dunul(2, 1)  = dunul(2, 1) + bt(3, iv)*delut(1, iv)
!
dunul(1, 2)  = dunul(1, 2) + bt(2, iv)*delut(2, iv)
dunul(2, 2)  = dunul(2, 2) + bt(3, iv)*delut(2, iv)
enddo
!
!if(ie.eq.1831) print*,'i1950g',deluq(1, :)
!
ugrad(1, 1) = x5(1, 1)*dunul(1 ,1) + x5(2, 1)*dunul(2 ,1)
ugrad(2, 1) = x5(1, 2)*dunul(1 ,1) + x5(2, 2)*dunul(2 ,1)

ugrad(1, 2) = x5(1, 1)*dunul(1 ,2) + x5(2, 1)*dunul(2 ,2)
ugrad(2, 2) = x5(1, 2)*dunul(1 ,2) + x5(2, 2)*dunul(2 ,2)
!
!if(ie.eq.1831) print*,'i1950g',ugrad(1:2,1)
!
!ucell = 0.d0
!
do iv = 1, 6
uptri(1, iv) = uptrc(1) + ugrad(1, 1)*bt(2, iv) + ugrad(2, 1)*bt(3, iv)
uptri(2, iv) = uptrc(2) + ugrad(1, 2)*bt(2, iv) + ugrad(2, 2)*bt(3, iv)
enddo
!
!vnulq(1, 1:8, ie) = ustar(1, ipq(1:8)) - upqua(1, 1:8)
!vnulq(2, 1:8, ie) = ustar(2, ipq(1:8)) - upqua(2, 1:8)
!
vnult(1, 1:6, ie) = vnult(1, 1:6, ie) + (ustar(1, ipt(1:6)) - uptri(1, 1:6))*gtdmp(ie)
vnult(2, 1:6, ie) = vnult(2, 1:6, ie) + (ustar(2, ipt(1:6)) - uptri(2, 1:6))*gtdmp(ie)
!
!if(ie.eq.1950) print*,'delu', upqua(1,1:8),ustar(1,ipq(1:8))
!
!...Imposing the dampening coefficient
!
!vnulq(1:2, 1:8, ie) = vnulq(1:2, 1:8, ie)*gqdmp(ie)

!
450 enddo
!
!print*,'vnult',vnult(1, 1:6, :)
!
!...get null mode...
!
end subroutine getnullmd_trialp

!
!...Get the null mode of velocity on physical domain with 4 vertices velocity fitting...
!
subroutine getnullmd_quadlp2(ipqua, geoel, vnulq, ustar, coord, gqdmp)
use constant
implicit none
integer,dimension(1:nvqua,1:nquad),         intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(inout):: vnulq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord !...nodal velocity

real*8,  dimension(1:nquad),                 intent(in)::gqdmp
!
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:3, 1:nvqua)::bq
real*8,dimension(1:2, 1:2)::amnul
real*8,dimension(1:2, 1:2)::dunul
real*8,dimension(1:ndimn, 1:nvqua)::upqua
real*8::ugrad(2, 2)
integer,dimension(1:nvqua) :: ipq
integer,dimension(npoin)::indnd
!
integer::iv, ie, ielem, ij, ik
real*8::amnua, amnub, amnuc, amnud, amnuj
real*8::dxm,dym,xcv,ycv
real*8::deluq(1:ndimn, 1:nvqua)
real*8::x5(2,2)
!
do 450 ie = 1,nquad !...(1)ie = 1,nelem
!
amnul = 0.d0
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
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
dxm = maxval(xpq(1, 1:npqua)) - minval(xpq(1, 1:npqua))
dym = maxval(xpq(2, 1:npqua)) - minval(xpq(2, 1:npqua))
!
xcv = xpq(1,9)
ycv = xpq(2,9)
!
do iv =1 ,nvqua
bq(1, iv) = 1.d0
bq(2, iv) = (xpq(1, iv)-xcv)/dxm
bq(3, iv) = (xpq(2, iv)-ycv)/dym
enddo
!
do iv=1, 4
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
!...Part 2
!
upqua(1, 1:8) = ustar(1, ipq(1:8))
upqua(2, 1:8) = ustar(2, ipq(1:8))
!
!if(ie.eq.1950) print*,'i1950',upqua(1,1:8),upqua(2,1:8)
!
upqua(1, 9) = -0.25d0*(upqua(1, 1) + upqua(1, 2) + upqua(1, 3) + upqua(1, 4)) +&
0.50d0*(upqua(1, 5) + upqua(1, 6) + upqua(1, 7) + upqua(1, 8))
upqua(2, 9) = -0.25d0*(upqua(2, 1) + upqua(2, 2) + upqua(2, 3) + upqua(2, 4)) +&
0.5d0*(upqua(2, 5) + upqua(2, 6) + upqua(2, 7) + upqua(2, 8))
!
!upqua(1, 9) = 0.125d0*(upqua(1, 1) + upqua(1, 2) + upqua(1, 3) + upqua(1, 4)) +&
!0.1250d0*(upqua(1, 5) + upqua(1, 6) + upqua(1, 7) + upqua(1, 8))
!upqua(2, 9) = 0.125d0*(upqua(2, 1) + upqua(2, 2) + upqua(2, 3) + upqua(2, 4)) +&
!0.125d0*(upqua(2, 5) + upqua(2, 6) + upqua(2, 7) + upqua(2, 8))
!
do iv = 1, 4
deluq(1:2, iv) = upqua(1:2, iv) - upqua(1:2, 9)
enddo
!
dunul = 0.d0
!
do iv=1, 4
!
dunul(1, 1)  = dunul(1, 1) + bq(2, iv)*deluq(1, iv)
dunul(2, 1)  = dunul(2, 1) + bq(3, iv)*deluq(1, iv)
!
dunul(1, 2)  = dunul(1, 2) + bq(2, iv)*deluq(2, iv)
dunul(2, 2)  = dunul(2, 2) + bq(3, iv)*deluq(2, iv)
enddo
!
!if(ie.eq.1831) print*,'i1950g',deluq(1, :)
!
ugrad(1, 1) = x5(1, 1)*dunul(1 ,1) + x5(2, 1)*dunul(2 ,1)
ugrad(2, 1) = x5(1, 2)*dunul(1 ,1) + x5(2, 2)*dunul(2 ,1)

ugrad(1, 2) = x5(1, 1)*dunul(1 ,2) + x5(2, 1)*dunul(2 ,2)
ugrad(2, 2) = x5(1, 2)*dunul(1 ,2) + x5(2, 2)*dunul(2 ,2)
!
!if(ie.eq.1831) print*,'i1950g',ugrad(1:2,1)
!
!ucell = 0.d0
!
do iv = 1, 8
upqua(1, iv) = upqua(1, 9) + ugrad(1, 1)*bq(2, iv) + ugrad(2, 1)*bq(3, iv)
upqua(2, iv) = upqua(2, 9) + ugrad(1, 2)*bq(2, iv) + ugrad(2, 2)*bq(3, iv)
enddo
!
!vnulq(1, 1:8, ie) = ustar(1, ipq(1:8)) - upqua(1, 1:8)
!vnulq(2, 1:8, ie) = ustar(2, ipq(1:8)) - upqua(2, 1:8)
!
vnulq(1, 1:8, ie) = vnulq(1, 1:8, ie) + (ustar(1, ipq(1:8)) - upqua(1, 1:8))*gqdmp(ie)
vnulq(2, 1:8, ie) = vnulq(2, 1:8, ie) + (ustar(2, ipq(1:8)) - upqua(2, 1:8))*gqdmp(ie)
!
vnulq(1:2, 9, ie) = 0.d0
!
!if(ie.eq.1950) print*,'delu', upqua(1,1:8),ustar(1,ipq(1:8))
!
!...Imposing the dampening coefficient
!
!vnulq(1:2, 1:8, ie) = vnulq(1:2, 1:8, ie)*gqdmp(ie)

!
450 enddo
!
!...get null mode...
!
end subroutine getnullmd_quadlp2
!
!...Get the null mode of velocity on reference domain without averaged velocity determined...
!
subroutine getnullmd_quadl2(ipqua, geoel, vnulq, ustar, gqdmp)
use constant
implicit none
integer,dimension(1:nvqua,1:nquad),         intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in):: geoel
real*8,dimension(1:ndimn,1:nvqua,1:nquad),   intent(out):: vnulq
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,  dimension(1:nquad),                 intent(in)::gqdmp
!
integer, parameter::ncbf = 2
real*8,dimension(1:nvqua):: xvq, yvq
real*8,dimension(1:3, 1:nvqua)::bq
real*8,dimension(1:3, 1:3)::amnul
real*8,dimension(1:3, 1:2)::dunul
real*8,dimension(1:ndimn, 1:nvqua)::upqua
real*8::ugrad(3, 2)
integer,dimension(1:nvqua) :: ipq
!
integer::iv, ie, ielem, ij, ik
real*8::amnua, amnub, amnuc, amnud, amnuj
real*8::dr, ds, rcq, scq
real*8::deluq(1:ndimn, 1:nvqua)
real*8::x5(3,3), mmatr(3,3), b55(3)
!
!...shape functions
!
dr = 1.0d0
ds = 1.0d0
!
rcq = 0.d0
scq = 0.d0
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
do iv =1 ,nvqua
!...Left cell + intfac(3,ifa)
bq(1, iv) = 1.d0
bq(2, iv) = (xvq(iv)-rcq)/dr
bq(3, iv) = (yvq(iv)-scq)/ds
enddo
!
!...zero vnulq
!
vnulq = 0.d0
amnul = 0.d0
!
do iv=1, 8
!
amnul(1, 1) = amnul(1, 1) + bq(1, iv)*bq(1, iv)
amnul(1, 2) = amnul(1, 2) + bq(1, iv)*bq(2, iv)
amnul(1, 3) = amnul(1, 3) + bq(1, iv)*bq(3, iv)

amnul(2, 2) = amnul(2, 2) + bq(2, iv)*bq(2, iv)
amnul(2, 3) = amnul(2, 3) + bq(2, iv)*bq(3, iv)

amnul(3, 3) = amnul(2, 2) + bq(3, iv)*bq(3, iv)
!
enddo
!
amnul(2 ,1) = amnul(1, 2)
amnul(3 ,1) = amnul(1, 3)
amnul(3 ,2) = amnul(2, 3)

!
!
!...Inverse matrix
!
call getinvmat(3,amnul, x5, b55)
!
!print*,'amunul',x5
!
!
do 350 ie = 1,nquad !...(1)ie = 1,nelem
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!
upqua(1, 1:8) = ustar(1, ipq(1:8))
upqua(2, 1:8) = ustar(2, ipq(1:8))
!
!if(ie.eq.1950) print*,'i1950',upqua(1,1:8),upqua(2,1:8)
!
do iv = 1, 8
deluq(1:2, iv) = upqua(1:2, iv) !- upqua(1:2, 9)
enddo
!
!if(ie.eq.1950) print*,'i19502',deluq(1, 1:8)
!
dunul = 0.d0
!
do iv=1, 8
!
dunul(1, 1)  = dunul(1, 1) + bq(1, iv)*deluq(1, iv)
dunul(2, 1)  = dunul(2, 1) + bq(2, iv)*deluq(1, iv)
dunul(3, 1)  = dunul(3, 1) + bq(3, iv)*deluq(1, iv)
!
dunul(1, 2)  = dunul(1, 2) + bq(1, iv)*deluq(2, iv)
dunul(2, 2)  = dunul(2, 2) + bq(2, iv)*deluq(2, iv)
dunul(3, 2)  = dunul(3, 2) + bq(3, iv)*deluq(2, iv)

enddo
!
!if(ie.eq.1950) print*,'i1950g',dunul(1:2,1)
!
ugrad = 0.d0
!
do ij = 1, 3
do ik = 1, 3
ugrad(ij, 1) = ugrad(ij, 1) + x5(ik, ij)*dunul(ik ,1)
ugrad(ij, 2) = ugrad(ij, 2) + x5(ik, ij)*dunul(ik ,2)
enddo
enddo
!
!if(ie.eq.1950) print*,'i1950g',ugrad(1:2,1),ugrad(1:2,2),upqua(1:2,9)
!
!ucell = 0.d0
!
do iv = 1, 8
upqua(1, iv) = ugrad(1, 1) + ugrad(2, 1)*bq(2, iv) + ugrad(3, 1)*bq(3, iv)
upqua(2, iv) = ugrad(1, 2) + ugrad(2, 2)*bq(2, iv) + ugrad(3, 2)*bq(3, iv)
enddo
!
vnulq(1, 1:8, ie) = ustar(1, ipq(1:8)) - upqua(1, 1:8)
vnulq(2, 1:8, ie) = ustar(2, ipq(1:8)) - upqua(2, 1:8)
!
!if(ie.eq.1950) print*,'delu', upqua(1,1:8),ustar(1,ipq(1:8))
!
!...Imposing the dampening coefficient
!
vnulq(1:2, 1:8, ie) = vnulq(1:2, 1:8, ie)*gqdmp(ie)

!
350 enddo

!
!...get null mode...
!
end subroutine getnullmd_quadl2

!
!...Dukowicz solver for middle point velocity...
!
subroutine getvelo_mpt_dukowicz(WL,WR,RHOL,RHOR,PL,PR,SSL,SSr,gamlg,p12,w12)
real*8, intent(in):: WL,WR,RHOL,RHOR,PL,PR,SSL,SSr,gamlg
real*8, intent(out)::p12, w12
real*8::wmin,wmax,plmin,prmin,al,ar
real*8::a,b,c,d,bl,br,dd
!
al = (gamlg+1.d0)/2.d0
ar = (gamlg+1.d0)/2.d0
!
WMIN=WR-0.5d0*SSR/AR
WMAX=WL+0.5d0*SSL/AL
PLMIN=PL-0.25d0*RHOL*SSL**2/AL
PRMIN=PR-0.25d0*RHOR*SSR**2/AR
BL=RHOL*AL
BR=RHOR*AR
A=(BR-BL)*(PRMIN-PLMIN)
B=BR*WMIN**2-BL*WMAX**2
C=BR*WMIN-BL*WMAX
D=BR*BL*(WMIN-WMAX)**2
! ***
! ***CASE A: W12-WMIN>O,Wl2-WMAX<O
! ***
DD=sqrt(AMAX1(0.d0,D-A))
W12=(B+PRMIN-PLMIN)/(C-SIGN(DD, WMAX-WMIN))
IF(W12-WMIN.GE.0.d0.AND.W12-WMAX.LE.0.) goto 10
! ***
! ***CASE B: W12-WMIN<O,Wl2-WMAX>O
! ***
dd=SQRT(AMAX1(0.d0, D+A))
W12=(B-PRMIN+PLMIN)/(C-SIGN(DD, WMAX-WMIN))
IF(W12-WMIN.LE.0..and.W12-WMAX.GE.0.) goto 10
A=(BL+BR)*(PLMIN-PRMIN)
B=bL*WMAX+BR*WMIN
C=1./(BL+BR)
! ***
! ***CASE C: W12-WMIN>O,Wl2-WMAx>O
! ***
DD=sqrt(AMAX1(0.d0,A-D))
W12=(B+DD)*C
IF(W12-WMIN.GE.0..AND.W12-WMAX.GE.0.) goto 10
! ***
! *** CASE D: W12-WMIN<O,Wl2-WMAx<O
! ***
DD=sqrt(AMAX1(0.d0,-A-D))
W12=(B-DD)*C

10 P12=0.5d0*(PLMIN+PRMIN+BR*ABS(W12-WMIN)*(Wl2-WMIN)-BL*ABS(W12-WMAX)*(W12-WMAX))

end subroutine getvelo_mpt_dukowicz

!
!...Subrotuine adjust the curved nodes...
!
subroutine getmpt_modify(intfac, coord)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::coord
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
!
!...Local integer
!
integer::ifa, ishp
!
!...local integer array
!
integer,dimension(1:nvfac) :: ipf
real*8, dimension(1:nvfac) :: dspr1, dspr2
!
real*8, dimension(1:ndimn, 1:nvfac)::xpf, xpm
real*8::cpara, cradc, cradm
real*8::dxl,dyl,dxm,dym
real*8::fmx, fmy,fnx,fny,len12,lenm
real*8::dxdr1,dxdr2,dydr1,dydr2
!
!...Adjustable parameter...
!
cpara = 0.1d0 !.sedov
!cpara = 1.d-3 !.sedov
!cpara = 0.01
!
!...Loop over elements
!
do 550 ifa = 1, nafac !...(1)ie = 1,nelem
!
!...Points consitituting one element...
!
ipf(1:nvfac) = intfac(3:(nvfac+2), ifa)
!
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))
!
!...Get the pre-set node xpm...
!
dxl = xpf(1, 2)-xpf(1, 1)
dyl = xpf(2, 2)-xpf(2,1)
len12 = sqrt(dxl**2 + dyl**2)
!
fnx = dyl/len12
fny =-dxl/len12
!
xpm = xpf
xpm(1:2, 3) = 0.5d0*(xpf(1:2, 1) + xpf(1:2, 2))
!
dxm = xpf(1,3) - xpm(1, 3)
dym = xpf(2,3) - xpm(2, 3)
!
lenm = sqrt(dxm**2 + dym**2)
!
fmx = dxm/lenm
fmy = dym/lenm
!
if((fnx*fmx+fny*fmy).lt.0.d0)then
!
xpm(1, 3) = xpm(1, 3) - cpara*len12*fnx
xpm(2, 3) = xpm(2, 3) - cpara*len12*fny
else
xpm(1, 3) = xpm(1, 3) + cpara*len12*fnx
xpm(2, 3) = xpm(2, 3) + cpara*len12*fny
endif
!
!....Get the curvature radius...
!
dspr1(1) = -0.5d0
dspr1(2) =  0.5d0
dspr1(3) =  0.d0
!
dspr2(1) =  1.d0
dspr2(2) =  1.d0
dspr2(3) = -2.d0
!
dxdr1 = 0.d0
dydr1 = 0.d0
dxdr2 = 0.d0
dydr2 = 0.d0
!
do ishp = 1, 3
dxdr1 = dxdr1 + dspr1(ishp)*xpf(1, ishp)
dydr1 = dydr1 + dspr1(ishp)*xpf(2, ishp)
!
dxdr2 = dxdr2 + dspr2(ishp)*xpf(1, ishp)
dydr2 = dydr2 + dspr2(ishp)*xpf(2, ishp)
enddo
!
cradc = abs((dxdr1**2+dydr1**2)**(1.5d0)/(dxdr1*dydr2 - dydr1*dxdr2))
!
!
dxdr1 = 0.d0
dydr1 = 0.d0
dxdr2 = 0.d0
dydr2 = 0.d0
!
do ishp = 1, 3
dxdr1 = dxdr1 + dspr1(ishp)*xpm(1, ishp)
dydr1 = dydr1 + dspr1(ishp)*xpm(2, ishp)
!
dxdr2 = dxdr2 + dspr2(ishp)*xpm(1, ishp)
dydr2 = dydr2 + dspr2(ishp)*xpm(2, ishp)
enddo
!
cradm = abs((dxdr1**2+dydr1**2)**(1.5d0)/(dxdr1*dydr2 - dydr1*dxdr2))
!
!if(ipf(3).eq.1683)print*,'bad',cradc,cradm,(dxdr1*dydr2 - dydr1*dxdr2)
!
if(cradc.le.cradm)then
!print*,'bad',cradc,cradm
coord(1:2, ipf(3)) = xpm(1:2, 3)
!coord(1:2, ipf(3)) = 0.5d0*(coord(1:2, ipf(1)) + coord(1:2, ipf(2))  )
endif
!
!
550 enddo
end subroutine getmpt_modify
!
!...Get the mass matrix for lagrangian based on mass center...
!
subroutine  getamatr_lagdensitycurv(amatr,matin,geoel,coord, iptri, ipqua)
use constant
implicit none
!...Input
real*8,dimension(1:ngeel,1:nsize)::geoel
real*8,dimension(1:ndimn,1:npoin),intent(in)::coord
real*8,dimension(1:nmatr,1:ncell),intent(out)::amatr
real*8,dimension(1:nmati,1:ncell),intent(out)::matin
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
mmatr = 0.d0
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
!   if(ie==400) print*,'gauss',xg, yg, f0
enddo
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
matin(1,ielem) = x5(1,1)
matin(2,ielem) = x5(1,2)
matin(3,ielem) = x5(1,3)
matin(4,ielem) = x5(2,2)
matin(5,ielem) = x5(2,3)
matin(6,ielem) = x5(3,3)
!
endif
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
!
mmatr = 0.d0
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
xg = r
yg = s
!
b2 = (xg-rc)/dr
b3 = (yg-sc)/ds
!
!
mmatr(1,1) = mmatr(1,1) + djaco
mmatr(2,1) = mmatr(2,1) + b2*djaco
mmatr(3,1) = mmatr(3,1) + b3*djaco

mmatr(2,2) = mmatr(2,2) + b2*b2*djaco
mmatr(3,2) = mmatr(3,2) + b2*b3*djaco

mmatr(3,3) = mmatr(3,3) + b3*b3*djaco
!   if(ie==400) print*,'gauss',xg, yg, f0
enddo
!
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
matin(1,ielem) = x5(1,1)
matin(2,ielem) = x5(1,2)
matin(3,ielem) = x5(1,3)
matin(4,ielem) = x5(2,2)
matin(5,ielem) = x5(2,3)
matin(6,ielem) = x5(3,3)
!
endif
enddo !...(2)ie = 1,nelem

end subroutine  getamatr_lagdensitycurv
!
!
!....rhs updating density for hybrid curvd quad cells
!
subroutine rhsdensity_lagcurv(iptri, ipqua, coord, coold, geoel, unkno, matin)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
real*8,dimension(1:nmati,1:ncell),intent(in)::matin
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
real*8,dimension(1:ndimn, 1:nptri) :: xpti
real*8,dimension(1:ndimn, 1:npqua) :: xpqi
real*8,dimension(1:ndegr):: b, dbdr, dbds
!
real*8, dimension(1:nptri):: shp, dspr, dsps
real*8::weigh(ngausd), posi(2,ngausd)

real*8, dimension(1:npqua):: shpq, dsprq, dspsq
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
!
if(ncurv==0)then
xpti(1, 1:3) = coold(1, iptri(1:3,ie))
xpti(2, 1:3) = coold(2, iptri(1:3,ie))
!
xpti(1:2,4) = 0.5d0*(xpti(1:2,1)+xpti(1:2,2))
xpti(1:2,5) = 0.5d0*(xpti(1:2,2)+xpti(1:2,3))
xpti(1:2,6) = 0.5d0*(xpti(1:2,1)+xpti(1:2,3))
elseif(ncurv==1)then
xpti(1, 1:nptri) = coold(1,iptri(1:nptri, ie))
xpti(2, 1:nptri) = coold(2,iptri(1:nptri, ie))
endif
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
!
rc= geoel(7, ielem) !...mass center...
sc= geoel(8, ielem)
!
!...Gauss loop
!
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
wi     = weigh(ig)
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
dxdr = dxdr + dspr(ishp)*xpti(1,ishp)
dxds = dxds + dsps(ishp)*xpti(1,ishp)

dydr = dydr + dspr(ishp)*xpti(2,ishp)
dyds = dyds + dsps(ishp)*xpti(2,ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
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
call  getrhoig_triacurv(rhoi, xpti)
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg, ielem)=rhsel(ideg, ielem) + rhoi*b(ideg)*djaco
enddo
!
enddo !...(2)ig = 1,ngausd
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
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
!
if(ncurv==0)then
xpqi(1, 1:4) = coold(1, ipqua(1:4,ie))
xpqi(2, 1:4) = coold(2, ipqua(1:4,ie))
!
xpqi(1:2,5) = 0.5d0*(xpqi(1:2,1)+xpqi(1:2,2))
xpqi(1:2,6) = 0.5d0*(xpqi(1:2,2)+xpqi(1:2,3))
xpqi(1:2,7) = 0.5d0*(xpqi(1:2,3)+xpqi(1:2,4))
xpqi(1:2,8) = 0.5d0*(xpqi(1:2,4)+xpqi(1:2,1))
xpqi(1:2,9) = 0.5d0*(xpqi(1:2,5)+xpqi(1:2,7))

elseif(ncurv==1)then
xpqi(1, 1:npqua) = coold(1,ipqua(1:npqua, ie))
xpqi(2, 1:npqua) = coold(2,ipqua(1:npqua, ie))
endif
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
!
rc= geoel(7, ielem) !...mass center...
sc= geoel(8, ielem)
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
do ishp = 1, npqua
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
!xpqi(1, 1:npqua) = coold(1, ipq(1:nvqua))
!xpqi(2, 1:npqua) = coold(2, ipq(1:nvqua))
!
call getrhoig_quadcurv(rhoi, xpqi)!
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg, ielem)=rhsel(ideg, ielem) + rhoi*b(ideg)*djaco
enddo
!
enddo !...(2)ig = 1,ngausd
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
m(1,1) = matin(1, ie)
m(1,2) = 0.d0!matin(2, ie)
m(1,3) = 0.d0!matin(3, ie)

m(2,1) = m(1,2)
m(2,2) = matin(4, ie)
m(2,3) = matin(5, ie)

m(3,1) = m(1,3)
m(3,2) = m(2,3)
m(3,3) = matin(6, ie)

endif
!
!...solve the mdegr independant varaible...
!
!
!...step 1
!
!if(ie==1771) print*,'rhsel',rhsel(2:3,ie),ie
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

end subroutine rhsdensity_lagcurv
!
!...Initial density distribtuion for curved quad..
!
subroutine getrhoig_quadcurv(rhoi, xpqi)
use constant
implicit none
real*8, intent(out)::rhoi
real*8,dimension(1:ndimn, 1:npqua), intent(in) :: xpqi
!
real*8::shpq(npqua)
real*8:: xc,yc,xg,yg
real*8::r,s
real*8::eps, c00, c10, c05,c20
real*8:: rhoini
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
elseif(ncase.eq.5)then !...Kidder ball...

!
elseif(ncase.eq.6)then !...Sod...!
!
!...  Mass center
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
elseif(ncase.eq.9)then !...Triple point...!
!
!...  Mass center
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
elseif(ncase.eq.13)then !...Saltzman...!
rhoini  = 1.d0

else

print*,'Please specify the initital density ditribution in subroutine getrhoig_quadcurv for Subgrid method for curved Quads!'
stop
endif

rhoi = rhoini
!
end subroutine getrhoig_quadcurv
!
!...Initial density distribtuion for curved triangle...
!
subroutine getrhoig_triacurv2(rhoi, xpti, xgaus, ygaus)
use constant
implicit none
real*8, intent(out)::rhoi
real*8,dimension(1:ndimn, 1:nptri), intent(in) :: xpti
real*8, intent(in)::xgaus, ygaus
!
real*8::shp(nptri)
real*8:: xg,yg,xc,yc, r,s
real*8::eps, c00, c10, c05,c20
real*8:: rhoini
real*8::radie, radii,radie2,radii2,radic2,sentr,rhoin,rhoex
real*8::prein,preex
real*8::rho0ba
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
elseif(ncase.eq.4)then !... Kidder shell...
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
elseif(ncase.eq.6)then !...Sod...!
!
r = 1.d0/3.d0; s=1.d0/3.d0
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
xc = xc + shp(ishp)*xpti(1,ishp)
yc = yc + shp(ishp)*xpti(2,ishp)
enddo
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
elseif(ncase.eq.9)then !...Triple point...!
!
r = 1.d0/3.d0; s=1.d0/3.d0
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
xc = xc + shp(ishp)*xpti(1,ishp)
yc = yc + shp(ishp)*xpti(2,ishp)
enddo
!
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
elseif(ncase.eq.12)then !...Isentropic smooth flow...!
!
rhoini = 1.d0 + 0.9999995d0*sin(pi*xgaus)
elseif(ncase.eq.13)then !...Saltzman...!
rhoini  = 1.d0

else

print*,'Please specify the initital density ditribution in subroutine getrhoig_tricurv for Subgrid method for curved Trias!'
stop
endif
!
rhoi = rhoini
!
end subroutine getrhoig_triacurv2
!
!...Initial density distribtuion for curved quad(2)..
!
subroutine getrhoig_quadcurv2(rhoi, xpqi, xgaus, ygaus)
use constant
implicit none
real*8, intent(out)::rhoi
real*8,dimension(1:ndimn, 1:npqua), intent(in) :: xpqi
real*8, intent(in)::xgaus, ygaus
!
real*8::shpq(npqua)
real*8:: xc,yc,xg,yg
real*8::r,s
real*8::eps, c00, c10, c05,c20
real*8:: rhoini
real*8::radie, radii,radie2,radii2,radic2,sentr,rhoin,rhoex
real*8::prein,preex
real*8::rho0ba
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
elseif(ncase.eq.6)then !...Sod...!
!
!...  Mass center
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
elseif(ncase.eq.9)then !...Triple point...!
!
!...  Mass center
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
elseif(ncase.eq.12)then !...Isentropic smooth flow...!
!
rhoini = 1.d0 + 0.9999995d0*sin(pi*xgaus)
!
elseif(ncase.eq.13)then !...Saltzman...!
rhoini  = 1.d0

elseif(ncase.eq.14)then
rhoini  = 1.d0
elseif(ncase.eq.15)then !...Shu-Osher

!...  Mass center
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
xc = xc + shpq(ishp)*xpqi(1,ishp)
yc = yc + shpq(ishp)*xpqi(2,ishp)
enddo
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

rhoi = rhoini
!
end subroutine getrhoig_quadcurv2

!
!....Density update for curved quad using llnl...
!
subroutine getdensity_quadllnl_curv(r, s, xpq, xpqi, rhoi, rhon)
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
jacon = dxdr*dyds - dxds*dydr
!
!...Initial Jacobian
!
dxdr = 0.d0
dxds = 0.d0
dydr = 0.d0
dyds = 0.d0
!
do ishp = 1, nvqua
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
end subroutine getdensity_quadllnl_curv
!
!...Initial density distribtuion for curved triangle...
!
subroutine getrhoig_triacurv(rhoi, xpti)
use constant
implicit none
real*8, intent(out)::rhoi
real*8,dimension(1:ndimn, 1:nptri), intent(in) :: xpti
!
real*8::shp(nptri)
real*8:: xg,yg,xc,yc, r,s
real*8::eps, c00, c10, c05,c20
real*8:: rhoini
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
elseif(ncase.eq.5)then !...Kidder ball...

!
elseif(ncase.eq.6)then !...Sod...!
!
!
r = 0.d0; s=0.d0
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
xc = xc + shp(ishp)*xpti(1,ishp)
yc = yc + shp(ishp)*xpti(2,ishp)
enddo
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
elseif(ncase.eq.9)then !...Triple point...!
!
r = 0.d0; s=0.d0
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
xc = xc + shp(ishp)*xpti(1,ishp)
yc = yc + shp(ishp)*xpti(2,ishp)
enddo
!
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
elseif(ncase.eq.13)then !...Saltzman...!
rhoini  = 1.d0

else

print*,'Please specify the initital density ditribution in subroutine getrhoig_quadcurv for Subgrid method for curved Trias!'
stop
endif
!
rhoi = rhoini
!
end subroutine getrhoig_triacurv
!
!...subroutine: Calculate the F^* N dsfor all faces for hybrid grids...
!
subroutine getfnds_lag_simpsonh(gflag,gelag,gelagq,intfac,inpoel,iptri,ipqua,coord)
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
!if(ie==1) print*,'coord ie==170',dxdr,dxds,dydr,dyds,ig
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
!if(ie.eq.1)then
!print*,'face',idfal,ie,farea,comatr(:, :),xpq(1, 1:9),xpq(2,1:9)
!endif
!
enddo !...ig from 1...9
!
200 enddo  !...(1)ifa=1,nquad
!
!  print*,'Inside getfnds_lag'
!
end subroutine getfnds_lag_simpsonh
!
!...Riemann input for curved triangle using simpson rule...
!
subroutine getriem_tria_simpson(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
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
real*8, dimension(1:2, 1:2, 1:9, 1:ntria),       intent(out)::munaclt
real*8, dimension(1:ndimn, 1:9,  1:ntria), intent(out)::munault
real*8, dimension(1:ndimn, 1:9,  1:ntria), intent(out)::snsigmlt
!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvfac) :: ipf
integer,dimension(nvfac, 3)::fglvt, fglgt
!...local real array
real*8,dimension(1:ndegr)::bg, bgv
real*8,dimension(1:nq,1:nvfac)::unkng
real*8::aujmp(1:3, 1:nvfac)
real*8::vnorm(1:3, 1:9)
real*8::sigmg(1:2, 1:2, 1:nvfac)
real*8,dimension(1:nvfac)::murie
real*8,dimension(1:nvtri):: xv,  yv
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
fglvt(1, 1) = 1;  fglvt(2, 1) = 2; fglvt(3, 1) = 4;
fglvt(1, 2) = 2;  fglvt(2, 2) = 3; fglvt(3, 2) = 5;
fglvt(1, 3) = 3;  fglvt(2, 3) = 1; fglvt(3, 3) = 6;
!
!...Local gauss point No. of any gauss point in one face...
!
fglgt(1, 1) = 1;  fglgt(2, 1) = 2; fglgt(3, 1) = 7;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4; fglgt(3, 2) = 8;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6; fglgt(3, 3) = 9;
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
dr = .5d0
ds = .5d0
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
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
!...Give the normal vector of every face...
!
vnorm(1:3,  1:9) = gelag(1:3, 1:9, ie);
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
do ifa =1, 3
!
!...zero out unkno
!
unkng = 0.d0
!
do ig   = 1, nvfac
!
rg = xv(fglvt(ig, ifa))
sg = yv(fglvt(ig, ifa))
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
aujmp(1:2, ig) = vlave(1:2, ipt(fglvt(ig, ifa))) - unkng(2:3, ig)
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
do ig   = 1, nvfac
dux= vlave(1, ipt(fglvt(ig, ifa)))-unkng(2, ig)
duy= vlave(2, ipt(fglvt(ig, ifa)))-unkng(3, ig)
deltu = sqrt(dux**2 + duy**2)
murie(ig) = rhoct*sdctr + rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
!if(ip(fagsp(ig, ifa)).eq.155) print*,'murie22', sdctr,rhoct,vlave(1:2, ip(fagsp(ig, ifa))),unkng(2:3, ig),unkno(1,2,ie),ie
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do ig  = 1, nvfac
!
! if(ip(iv).eq.5) print*,'p19 muacn(28) pre++', munacn(5),ie,iv,ifa
!
call getriecoef_matrixnew(murie(ig), vnorm(3, fglgt(ig, ifa)), vnorm(1:2, fglgt(ig, ifa)), aujmp(1:3, ig), &
unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ig), vnorm(3, fglgt(ig, ifa)), vnorm(1:2, fglgt(ig, ifa)), aujmp(1:2, ig), &
!unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
munacn(1:2, 1, ipt(fglvt(ig, ifa))) = munacn(1:2, 1, ipt(fglvt(ig, ifa))) + munacn_rie(1:2, 1)
munacn(1:2, 2, ipt(fglvt(ig, ifa))) = munacn(1:2, 2, ipt(fglvt(ig, ifa))) + munacn_rie(1:2, 2)
!
munacu(1:2, ipt(fglvt(ig, ifa))) = munacu(1:2, ipt(fglvt(ig, ifa))) + munacu_rie(1:2)
!
snsigm(1:2,ipt(fglvt(ig, ifa))) = snsigm(1:2, ipt(fglvt(ig, ifa))) + snsigm_rie(1:2)!
! if(ip(iv).eq.15) print*,'p19 muacn(28) post-snsigmaxxxx',sigma(:,:,iv),vnorm(1:3,ifa,iv),snsigm(1:2, ip(iv)),ie, ifa,iv
!,vnorm(3, ifa, iv),vnorm(1:2, ifa, iv),aujmp(1:2, iv)
!
munaclt(1:2, 1, fglgt(ig, ifa), ie) =  munacn_rie(1:2, 1)
munaclt(1:2, 2, fglgt(ig, ifa), ie) =  munacn_rie(1:2, 2)
!
! if(ie.eq.23.and.fglgt(ig, ifa).eq.1) print*,'p11 muacn(1) post',munacn_rie(1:2, 1),aujmp(3,ig),ig,vnorm(1:3, fglgt(ig, ifa))
!
munault(1:2, fglgt(ig, ifa), ie)    =  munacu_rie(1:2)
!
snsigmlt(1:2,fglgt(ig, ifa), ie)   = snsigm_rie(1:2)
!
! if(ipt(fglvt(ig, ifa)).eq.54) print*,'p11 muacn(1) prep--munacl',murie(ig),aujmp(1:2, ig),sigmg(1:2, 1:2, ig),snsigm_rie(1:2),ie,ig
enddo
!
!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,nelem!

end subroutine getriem_tria_simpson
!
!...Riemann input for curved quad using simpson rule...
!
subroutine getriem_quad_simpson(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec, vnulq)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
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
data eps   / 1.0d-6/
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
vnorm(1:3,  1:12) = gelagq(1:3, 1:12, ie);

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
pg = pctr + aflim(4, ielem)*(pg - pctr)

!...updtae unknv(2:3,:)
unkng(2, ig) = ug -  vnulq(1,fglvq(ig, ifa),ie)
unkng(3 ,ig) = vg -  vnulq(2,fglvq(ig, ifa),ie)
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
dux= vlave(1, ipq(fglvq(ig, ifa)))-unkng(2, ig)
duy= vlave(2, ipq(fglvq(ig, ifa)))-unkng(3, ig)
deltu = 1.d0*sqrt(dux**2 + duy**2)
!
gpnx = vnorm(1, fglgq(ig, ifa))
gpny = vnorm(2, fglgq(ig, ifa))
dux = vlave(1, ipq(fglvq(ig, ifa)))-unkno(1, 2, ielem)
duy = vlave(2, ipq(fglvq(ig, ifa)))-unkno(1, 3, ielem)
!
divu = dux*lnvp(1, fglvq(ig, ifa)) + duy*lnvp(2, fglvq(ig, ifa))

if(divu.le.0.d0)then
!deltu = abs(dux*gpnx + duy*gpny)
else
!deltu = 0.d0
endif
deltu = abs(dux*gpnx + duy*gpny)
sdimp = sdctr!max(1d-8,sdctr)
murie(ig) = rhoct*sdimp  + cdrho*rhoct*slpdu*deltu !...slpdu denotes the slope of delt u
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
! if(ipq(fglvq(ig, ifa)).eq.87) print*,'ep54 ',ie,ifa,ig,ipq(fglvq(ig, ifa)),murie(ig),fglvq(ig, ifa),&
!sigmg(1, 1, ig),unkng(2:3, ig)&
!,vnorm(1:3, fglgq(ig, ifa)),fglgq(ig, ifa),munacn_rie(1,1),aujmp(1:2,ig)
!
enddo
!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,ntria!

end subroutine getriem_quad_simpson
!
!...Get the nodal velocity based on simpson integral rule...
!
subroutine getndvelo_lag_simpsonh(gflag,gelag,gelagq,geoel,bface,intfac,iptri,ipqua,&
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
if(ntria.gt.0) call getriem_tria_simpson(iptri, geoel, gelag, vlave, unkno, munacn, munacu, snsigm,&
munaclt, munault, snsigmlt, coord, coold,aflim, afvec)

!...Quad
if(nquad.gt.0) call getriem_quad_simpson(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec, vnulq)

!...Boundary condition...
!call getbcfc_lag(bface, intfac, gflag, fpres,  coord, ustar, itime)
!call getboundary_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm)

!...Periodic boundary condition for 1D isentropic Sin problem...
if(ncase.eq.12) call getbc_prdic(bface, munacn, munacu, snsigm)

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
!if(ipoin.eq.87) print*,'ep',ustar(1:2,ipoin),detma,munacn(1,1,ipoin),snsigm(1:2, ipoin),munacu(1:2, ipoin)
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
!if(ie.eq.1.and.ifa.eq.1) print*,ifa,ig,snsigmlq(1, fglgq(ig, ifa), ie),&
!munaclq(1,1, fglgq(ig, ifa), ie),ustar(1, ipq(fglvq(ig, ifa))),&
!munaclq(2,1, fglgq(ig, ifa), ie),ustar(2, ipq(fglvq(ig, ifa))),&
!munaulq(1, fglgq(ig, ifa), ie)
!
enddo
enddo
enddo
!
deallocate (munacn, bpres, fpres)
deallocate (munacu, snsigm, bnorm)
deallocate (munaclt, snsigmlt, munault)
deallocate (munaclq, snsigmlq, munaulq)
end subroutine getndvelo_lag_simpsonh
!
!...Face integral using gauss quadrature distribution on cuvred triangle...
!
subroutine rhsifacedg_lagtria_simpson(inpoel,  unkno, ustar, fstart, gelag, geoel,&
rhsel)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:ndimn,1:9,1:ntria),       intent(in)::fstart !...Riemann forces
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(out)::rhsel
real*8,dimension(1:3,1:ngelg,1:ntria+nbfac), intent(in)::gelag
real*8,dimension(1:ngeel,1:nsize), intent(in)::geoel
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ielem
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer, dimension(3, 3)::fglvt, fglgt
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::vnorm(3, 9)
real*8::xv(nvtri), yv(nvtri)
real*8::bg(3, nvfac)
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
fglvt(1, 1) = 1;  fglvt(2, 1) = 2; fglvt(3, 1) = 4;
fglvt(1, 2) = 2;  fglvt(2, 2) = 3; fglvt(3, 2) = 5;
fglvt(1, 3) = 3;  fglvt(2, 3) = 1; fglvt(3, 3) = 6;
!
posi(1, 1) = -1.d0; posi(1 ,2 )= 1.d0; posi(1 ,3 )= 0.d0
!
weigh(1) = 1.d0/6.d0; weigh(2) = 1.d0/6.d0; weigh(3) = 4.d0/6.d0
!
!...Local gauss point No. of any gauss point in one face...
!
fglgt(1, 1) = 1;  fglgt(2, 1) = 2; fglgt(3, 1) = 7;
fglgt(1, 2) = 3;  fglgt(2, 2) = 4; fglgt(3, 2) = 8;
fglgt(1, 3) = 5;  fglgt(2, 3) = 6; fglgt(3, 3) = 9;
!
endif
!
!...Zero out plnpn, ulnpn
!
dr = .5d0
ds = .5d0
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
if(ncurv.eq.1)then
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
endif
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
ielem = ie
!
!...Initialize ulnpn, plnpn, elnpn
!
ulnpn = 0.d0
plnpn = 0.d0
elnpn = 0.d0
!
ipt(1:nvtri) = inpoel(1:nvtri ,ie)
vnorm(1:3,  1:9) = gelag(1:3, 1:9, ie)
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
do ifa = 1, 3
!
do ig   = 1, nvfac
!
wi  = weigh(ig)
!
rg = xv(fglvt(ig, ifa))
sg = yv(fglvt(ig, ifa))
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
gpnx = vnorm(1, fglgt(ig, ifa))
gpny = vnorm(2, fglgt(ig, ifa))
gpsa = vnorm(3, fglgt(ig, ifa))
!
!...Distribute to every corner...
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustar(1, ipt(fglvt(ig, ifa)))*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig) +&
ustar(2, ipt(fglvt(ig, ifa)))*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)
!
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstart(1, fglgt(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstart(2, fglgt(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
elnpn(1:ndegr)   = elnpn(1:ndegr)+&
ustar(1, ipt(fglvt(ig, ifa)))*fstart(1, fglgt(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig) +&
ustar(2, ipt(fglvt(ig, ifa)))*fstart(2, fglgt(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
!if(ie==54) print*,'rhs iface idegr',ustar(1:2, ipt(fglvt(ig, ifa))),fstart(1:2, fglgt(ig, ifa), ie),ig,ifa,ipt(fglvt(ig, ifa))
!
enddo
!
enddo
!
rhsel(1:ndegr, 1, ie) =  ulnpn(1:ndegr)
rhsel(1:ndegr, 2, ie) =  plnpn(1, 1:ndegr)
rhsel(1:ndegr, 3, ie) =  plnpn(2, 1:ndegr)
rhsel(1:ndegr, 4, ie) =  elnpn(1:ndegr)
!
!if(ie==23) print*,'rhs iface',rhsel(1:3, 4, ie), fstart(1:2,7,ie)!, &
!                               ustar(1:2,ip(2)),ustar(1:2,ip(3)), ip(1:3)
550 enddo
!
!   open(8,file='rhstxt.dat')
!    do ie = 1, ntria
!   write(8,*) ie, rhsel(1, 1:4, ie)
!    enddo
!   close(8)
!
end subroutine rhsifacedg_lagtria_simpson
!
!...Face integral using gauss quadrature distribution on cuvred quads...
!
subroutine rhsifacedg_lagquad_simpson(ipqua,  unkno, ustar, fstarq, gelagq, geoel,&
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
vnorm(1:3,  1:12) = gelagq(1:3, 1:12, ie)

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
!if(ie==1) print*,'rhs iface idegr',ipq(fglvq(ig, ifa)),ig,ifa,ustar(1:2, ipq(fglvq(ig, ifa)))

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
end subroutine rhsifacedg_lagquad_simpson
!
!...subroutine: Calculate the F^* N dsfor sub grids...
!
subroutine getgeosg_lag(geoel,geoels,geosgq,intfac,iptri,ipqua,coord)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria),        intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel, 1:nsize),          intent(in)::geoel
!
integer,dimension(1:4,1:4,1:nquad)::ipqsg
real*8,dimension(1:2,1:4,1:nquad),        intent(out)::geoels
real*8,dimension(1:3,1:4,1:nquad),        intent(out)::geosgq
!...Local integer
integer::ifa,iel,ier,ie,isq,ielem
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(1:4,1:4) :: ivqsg
!...local real array
real*8,dimension(1:ndimn, 1:nvfac)::xpf
real*8,dimension(1:ndimn, 1:nvtri)::xp
real*8,dimension(1:ndimn, 1:nvqua)::xpq,xsq
!...local real number
real*8::anx, any,volq
real*8::c16
!
data c16   /0.1666666666666666d0 /
!  print*,'vnotmfn',gelag(1, 3, 9)
!
ivqsg(1, 1) = 1; ivqsg(2, 1) = 5;
ivqsg(3, 1) = 9; ivqsg(4, 1) = 8;
!
!...Cell 2
!
ivqsg(1, 2) = 5; ivqsg(2, 2) = 2;
ivqsg(3, 2) = 6; ivqsg(4, 2) = 9;
!
!...Cell 3
!
ivqsg(1, 3) = 9; ivqsg(2, 3) = 6;
ivqsg(3, 3) = 3; ivqsg(4, 3) = 7;
!
!...Cell 4
!
ivqsg(1, 4) = 8; ivqsg(2, 4) = 9;
ivqsg(3, 4) = 7; ivqsg(2, 4) = 4;
!
!...Quads...
!
do 150 ie=1, nquad !...(1)ifa=1,nafac
!
!  print*,'ifa',ifa
!
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
!...coordinates
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Cell 1
!
ipqsg(1:4, 1, ie) = ipq(ivqsg(1:4, 1))
ipqsg(1:4, 2, ie) = ipq(ivqsg(1:4, 2))
ipqsg(1:4, 3, ie) = ipq(ivqsg(1:4, 3))
ipqsg(1:4, 4, ie) = ipq(ivqsg(1:4, 4))
!
!...The mass of every subgrid...
!...1--4
geoels(1, 1:4, ie) = 0.25d0*geoel(4, ielem)
!
!...The volume of every subgrid...
!...5--8
do isq = 1,4
xsq(1, 1:4) = xpq(1, ivqsg(1:4, isq))
xsq(2, 1:4) = xpq(2, ivqsg(1:4, isq))
!
xsq(1:2,5) = 0.5d0*(xsq(1:2,1)+xsq(1:2,2))
xsq(1:2,6) = 0.5d0*(xsq(1:2,2)+xsq(1:2,3))
xsq(1:2,7) = 0.5d0*(xsq(1:2,3)+xsq(1:2,4))
xsq(1:2,8) = 0.5d0*(xsq(1:2,4)+xsq(1:2,1))
xsq(1:2,9) = 0.5d0*(xsq(1:2,5)+xsq(1:2,7))
!
call getvol_sgq(xsq, volq)
!
geoels(2, isq, ie) = volq
enddo
!
!...LN....
!9--12
!
anx =   xpq(2, 9) - xpq(2, 5)
any = -(xpq(1, 9) - xpq(1, 5))
!
geosgq(1, 1, ie) = anx/sqrt(anx**2 + any**2)
geosgq(2, 1, ie) = any/sqrt(anx**2 + any**2)
geosgq(3, 1, ie) = sqrt(anx**2 + any**2)
!
anx =   xpq(2, 9) - xpq(2, 6)
any = -(xpq(1, 9) - xpq(1, 6))
!
geosgq(1, 2, ie) = anx/sqrt(anx**2 + any**2)
geosgq(2, 2, ie) = any/sqrt(anx**2 + any**2)
geosgq(3, 2, ie) = sqrt(anx**2 + any**2)
!
anx =   xpq(2, 9) - xpq(2, 7)
any = -(xpq(1, 9) - xpq(1, 7))
!
geosgq(1, 3, ie) = anx/sqrt(anx**2 + any**2)
geosgq(2, 3, ie) = any/sqrt(anx**2 + any**2)
geosgq(3, 3, ie) = sqrt(anx**2 + any**2)
!
anx =   xpq(2, 9) - xpq(2, 8)
any = -(xpq(1, 9) - xpq(1, 8))
!
geosgq(1, 4, ie) = anx/sqrt(anx**2 + any**2)
geosgq(2, 4, ie) = any/sqrt(anx**2 + any**2)
geosgq(3, 4, ie) = sqrt(anx**2 + any**2)
!
150 enddo  !...(1)ifa=1,nelem

!
end subroutine getgeosg_lag
!
!...Find new volume of subgrid for curved cell...
!
subroutine getvol_sgq(xpq, volq)
use constant
implicit none
real*8,dimension(1:2, 1:npqua), intent(in):: xpq
real*8, intent(out):: volq
real*8 :: rcvq, scvq
!
!...local array...
!
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausd_geoq), posiq(2, ngausd_geoq)
!...local real number
real*8:: dxdr,dxds,dydr,dyds
real*8:: eps,c00,c10,c05,c20
real*8:: r, s, djaco, volel, masel,wi
real*8:: rcv, scv
!
integer:: igaus, ishp
integer::ie
!
data eps / 1.0d-06 /
data c00 / 0.0d0 /
data c10 / 1.0d0 /
data c05 / 0.5d0 /
data c20 / 2.0d0 /
!
call ruqope(2, ngausd_geoq, posiq, weighq)

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
volel = volel + djaco
rcv = rcv + r*djaco
scv = scv + s*djaco
!
enddo
!
volq = volel
rcvq = rcv/volel
scvq = scv/volel
!geoel(5, ielem) = geoel(1, ielem

end subroutine getvol_sgq
!
!...Face integral using high-order gauss quadrature distribution on cuvred quads...
!
subroutine rhsifacedg_lagquadsimp(ipqua,  unkno, ustar,  geoel, coord, coold, aflim, afvec,&
rhsel)
use constant
implicit none
!...Input arrays
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar, coord, coold!...nodal velocity
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(out)::rhsel
real*8,dimension(1:ngeel,1:nsize), intent(in)::geoel
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ielem
!...local integer array
real*8, dimension(1:ndimn, 1:16,  1:nquad)::ustarq
real*8,dimension(1:ndimn,1:16,1:nquad) ::fstarq !...Riemann forces
real*8,dimension(1:3, 1:16,  1:nquad)::vnormq
integer,dimension(1:nvqua) :: ipq
integer, dimension(nvfac, 4)::fglvq
integer,dimension(ngausf, 4)::fglgq
real*8, dimension(1:ndegr) :: ulnpn, elnpn
real*8, dimension(1:ndimn, 1:ndegr) :: plnpn
real*8::posiq(1:2, 1:16)
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
if(ngausf.eq.4)then
!
!...Local vertex No. of gauss points in one unit ...
!
fglvq(1, 1) = 1;  fglvq(2, 1) = 2; fglvq(3, 1) = 5;
fglvq(1, 2) = 2;  fglvq(2, 2) = 3; fglvq(3, 2) = 6;
fglvq(1, 3) = 3;  fglvq(2, 3) = 4; fglvq(3, 3) = 7;
fglvq(1, 4) = 4;  fglvq(2, 4) = 1; fglvq(3, 4) = 8;
!
posi(1, 1) = -1.d0; posi(1 ,2)= 1.d0;
posi(1 ,3)= -0.4472135954999579392818d0;
posi(1 ,4)=  0.4472135954999579392818d0;
!
weigh(  1) = 0.166666666666666666666d0/2.d0
weigh(  2) = 0.166666666666666666666d0/2.d0
weigh(  3) = 0.833333333333333333333d0/2.d0
weigh(  4) = 0.833333333333333333333d0/2.d0

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
!...Get high-order quadratrature point...
!
call getvelo_quad(ipqua, geoel,fstarq, coord, coold,vnormq,ustarq,ustar,unkno, aflim,afvec)
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
gpnx = vnormq(1, fglgq(ig, ifa),ie)
gpny = vnormq(2, fglgq(ig, ifa),ie)
gpsa = vnormq(3, fglgq(ig, ifa),ie)
!
!...Distribute to every corner...
!
ulnpn(1:ndegr)  = ulnpn(1:ndegr)+&
ustarq(1, fglgq(ig, ifa),ie)*gpnx*gpsa*bg(1:ndegr, ig)*weigh(ig) +&
ustarq(2, fglgq(ig, ifa),ie)*gpny*gpsa*bg(1:ndegr, ig)*weigh(ig)
!
!
plnpn(1, 1:ndegr)= plnpn(1, 1:ndegr)   +&
fstarq(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
plnpn(2, 1:ndegr)= plnpn(2, 1:ndegr)  +&
fstarq(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
elnpn(1:ndegr)=elnpn(1:ndegr) +&
ustarq(1, fglgq(ig, ifa), ie)*&
fstarq(1, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig) +&
ustarq(2, fglgq(ig, ifa), ie)*&
fstarq(2, fglgq(ig, ifa), ie)*bg(1:ndegr, ig)*weigh(ig)
!
!if(ie==1831) print*,'rhs iface idegr',ifa,ig,fglgq(ig, ifa),ulnpn(1:ndegr),ustarq(1, fglgq(ig, ifa),ie),gpnx,gpny,gpsa
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
550 enddo
!
end subroutine rhsifacedg_lagquadsimp
!
!...The nodal velocity for high-order quardrature nodes..
!
subroutine getvelo_quad(ipqua, geoel,fstarq, coord, coold,vnormq, ustarq, ustar,unkno,aflim,afvec)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),            intent(in)::coord, coold, ustar
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:nq+1,1:nsize),  intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),          intent(in)::afvec
!
real*8, dimension(1:2, 1:2, 1:16, 1:nquad)::munaclq_g
real*8, dimension(1:ndimn, 1:16,  1:nquad)::munaulq_g
real*8, dimension(1:ndimn, 1:16,  1:nquad)::snsigmlq_g
real*8, dimension(1:ndimn, 1:16,  1:nquad), intent(out)::ustarq
real*8, dimension(1:ndimn, 1:16,  1:nquad), intent(out)::fstarq
real*8, dimension(1:3, 1:16,  1:nquad), intent(out)::vnormq
!...Local integer
integer::ie, ideg, ielem, ifa, iv,ig,idfal,ishp
!...local integer array
real*8, dimension(1:ndimn, 1:16,  1:nquad)::vlave
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(nvfac, 4)::fglvq
integer,dimension(ngausf, 4)::fglgq
!...local real array
real*8::dsprq(9),dspsq(9), shpq(9)
real*8::xpq(1:ndimn, 1:9)
real*8::comatr(1:2, 1:2)
real*8,dimension(1:ndegr)::bg, bgv
real*8,dimension(1:nq,1:ngausf)::unkng
real*8::aujmp(1:3, 1:ngausf)
real*8::vnormg(1:3, 1:16), vnorm(1:2)
real*8::posiq(1:2, 1:16)
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
real*8::acnx, acny,anx,any
real*8::dxdr, dxds,dydr,dyds,farea,larea, r, s
real*8::shp1,shp2, shp3
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
!
endif
!
!...Get the velocity...
!
do 100 ie=1, nquad !...(1)ifa=1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
do ifa = 1,4
!
ustarq(1:ndimn,fglgq(1, ifa),  ie) = ustar(1:2, ipq(fglvq(1, ifa)))
ustarq(1:ndimn,fglgq(2, ifa),  ie) = ustar(1:2, ipq(fglvq(2, ifa)))
!
r = -0.4472135954999579392818d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
!
ustarq(1:ndimn,fglgq(3, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
shp3*ustar(1:2, ipq(fglvq(3, ifa)))
!
r =  0.4472135954999579392818d0
shp1 =  -0.5d0*(1.d0-r)*r
shp2 =   0.5d0*(1.d0+r)*r
shp3 =         (1.d0+r)*(1.d0-r)
ustarq(1:ndimn,fglgq(4, ifa),  ie) = shp1*ustar(1:2, ipq(fglvq(1, ifa))) + shp2*ustar(1:2, ipq(fglvq(2, ifa)))+&
shp3*ustar(1:2, ipq(fglvq(3, ifa)))
enddo
!
!if(ie.eq.1831)print*,'ustaq',ustar(1, ipq(fglvq(3, 1))),ipq(fglvq(3, 1))
!
100 enddo
!
vlave = ustarq
!
!...Get the norm at the gauss...
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
!
do ig = 1, 16!...(2)ig = 1,ngausd
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
vnormq(1, idfal, ie) = anx*larea/farea
vnormq(2, idfal, ie) = any*larea/farea
vnormq(3, idfal, ie) = farea
!
enddo !...ig from 1...9
!
!print*,'fglgqxx3',ie,fglgq
!
200 enddo  !...(1)ifa=1,nquad
!

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
!
!...Part II: Triangle...
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
vnormg(1:3,  1:16) = vnormq(1:3, 1:16, ie);
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
!print*,'fglgqxx',fglgq
!
do ifa =1, 4
!
!...zero out unkno
!
unkng = 0.d0
!
do ig   = 1, ngausf
!
!print*,'fglgq',ig,ifa,fglgq(ig, ifa)
rg = posiq(1, fglgq(ig, ifa))
sg = posiq(2, fglgq(ig, ifa))
!
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
aujmp(1:2, ig) = vlave(1:2, fglgq(ig, ifa), ie) - unkng(2:3, ig)
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
!if(ie.eq.1892) print*,ie,ig,ifa, vlave(1:2, fglgq(ig, ifa), ie), unkng(2:3, ig)
!
enddo !....ig
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
do ig   = 1, ngausf
murie(ig) = rhoct*sdctr   !...slpdu denotes the slope of delt u
enddo
!
!if(ie==3) print*,'vnotm',vnorm(1:3,1,3)!,gelag(1, 3, 9),aujmp(1:2,1)
!
!...Get the summed denominator cooefficients sum(mu*n*a_c)
!!
do ig  = 1, ngausf
!
!
call getriecoef_matrixnew(murie(ig), vnormg(3, fglgq(ig, ifa)), vnormg(1:2, fglgq(ig, ifa)), aujmp(1:3, ig), &
unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!call getriecoef_vilar(murie(ig), vnormg(3, fglgq(ig, ifa)), vnormg(1:2, fglgq(ig, ifa)), aujmp(1:2, ig), &
!unkng(2:3, ig), sigmg(1:2, 1:2, ig),&
!munacn_rie, munacu_rie(1:2), snsigm_rie(1:2))
!
!
! if(ipq(fglvq(ig, ifa)).eq.62) print*,'p19 muacn(28) post-snsigmaxxxx',murie(ig),ie
!
munaclq_g(1:2, 1, fglgq(ig, ifa), ie) =  munacn_rie(1:2, 1)
munaclq_g(1:2, 2, fglgq(ig, ifa), ie) =  munacn_rie(1:2, 2)
!

!
! if(ip(iv).eq.15) print*,'p11 muacn(1) post',munacl(1:2,iv,ie),murie(iv),aujmp(1:2, iv),vnormg(1:3, 1, iv),ie,iv
!
munaulq_g(1:2, fglgq(ig, ifa), ie)    =  munacu_rie(1:2)
!
snsigmlq_g(1:2,fglgq(ig, ifa), ie)   = snsigm_rie(1:2)
!
!if(ie.eq.1892) print*,ie,ifa,ig,fglgq(ig, ifa),snsigmlq_g(1:2, fglgq(ig, ifa), ie)!murie(ig),vnormg(1:3, fglgq(ig, ifa)),aujmp(1:3, ig)
enddo
!
! if(ip(iv).eq.15) print*,'p11 muacn(1) prep--munacl',murie(iv),aujmp(1:2, iv),vnorm(1:3, 1, iv),ie,iv
!
!
enddo ! ifa
!
250 enddo  !...(1)ie = 1,nelem!
!
!...Get the force...
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
fstarq(1, fglgq(ig, ifa), ie) = snsigmlq_g(1, fglgq(ig, ifa), ie) +&
munaclq_g(1,1, fglgq(ig, ifa), ie)*ustarq(1, fglgq(ig, ifa), ie)+&
munaclq_g(2,1, fglgq(ig, ifa), ie)*ustarq(2, fglgq(ig, ifa), ie)-&
munaulq_g(1, fglgq(ig, ifa), ie)
!
!if(ie.eq.1892) print*,'fstarq',ifa,ig,snsigmlq_g(1, fglgq(ig, ifa), ie),munaclq_g(1,1, fglgq(ig, ifa), ie),&
!munaulq_g(1, fglgq(ig, ifa), ie),fstarq(1, fglgq(ig, ifa), ie),ustarq(1:2, fglgq(ig, ifa), ie)
!
fstarq(2, fglgq(ig, ifa), ie) = snsigmlq_g(2, fglgq(ig, ifa), ie) +&
munaclq_g(2,2,fglgq(ig, ifa), ie)*ustarq(2, fglgq(ig, ifa), ie)+&
munaclq_g(1,2,fglgq(ig, ifa), ie)*ustarq(1, fglgq(ig, ifa), ie)-&
munaulq_g(2, fglgq(ig, ifa), ie)
!
enddo
enddo
enddo

end subroutine getvelo_quad
!
!...subroutine: Get time step size for Lagrangian hydrodynamics...
!
subroutine getdelt_lag(ustar,bface,unkno,intfac,iptri,ipqua,geoel,coord,dvdt,itime)
use constant
implicit none

!...input arrays
integer*4,dimension(1:nbfai,nbfac)::bface
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer*4,dimension(1:nifai,1:nafac), intent(in)::intfac
real*8,   dimension(1:ngeel,1:nsize), intent(in)::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),  intent(in)::ustar !...nodal velocity
real*8,dimension(1:ncell),          intent(in)::dvdt
!
integer, intent(in)::itime

!...local integer array
integer,dimension(1:nvqua) :: ipq

!...local integer
integer::ie, ielem

!...local real array
real*8:: xpq(1:2, nvqua)
real*8:: dledgq(4), dledgt(3)

!...local real
real*8::eps
real*8:: dtcfl, dtvol, dlen, dlenda
real*8:: rhomc, rhoct, sdctr, pctr, uctr, vctr, ectr
!
data eps   / 1.0d-06 /

!...Set big No. for dtcfl,dtvol
dtcfl = 1d10
dtvol = 1d10

!...I: Get the time step size for cfl
do 650 ie = 1,nquad
!
ielem = ie + ntria

!...Points consitituting one element...
ipq(1:nvqua) = ipqua(1:nvqua,ie)!

xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
if(ncurv.eq.0)then
  dledgq(1) = sqrt((xpq(1, 1)-xpq(1, 2))**2 + (xpq(2, 1)-xpq(2, 2))**2)
  dledgq(2) = sqrt((xpq(1, 2)-xpq(1, 3))**2 + (xpq(2, 2)-xpq(2, 3))**2)
  dledgq(3) = sqrt((xpq(1, 3)-xpq(1, 4))**2 + (xpq(2, 3)-xpq(2, 4))**2)
  dledgq(4) = sqrt((xpq(1, 4)-xpq(1, 1))**2 + (xpq(2, 4)-xpq(2, 1))**2)
elseif(ncurv.ge.1)then
 dledgq(1) = sqrt((xpq(1, 1)-xpq(1, 5))**2 + (xpq(2, 1)-xpq(2, 5))**2) +&
             sqrt((xpq(1, 5)-xpq(1, 2))**2 + (xpq(2, 5)-xpq(2, 2))**2)
 dledgq(2) = sqrt((xpq(1, 2)-xpq(1, 6))**2 + (xpq(2, 2)-xpq(2, 6))**2) +&
             sqrt((xpq(1, 6)-xpq(1, 3))**2 + (xpq(2, 6)-xpq(2, 3))**2)
 dledgq(3) = sqrt((xpq(1, 3)-xpq(1, 7))**2 + (xpq(2, 3)-xpq(2, 7))**2) +&
             sqrt((xpq(1, 7)-xpq(1, 4))**2 + (xpq(2, 7)-xpq(2, 4))**2)
 dledgq(4) = sqrt((xpq(1, 4)-xpq(1, 8))**2 + (xpq(2, 4)-xpq(2, 8))**2) +&
             sqrt((xpq(1, 8)-xpq(1, 1))**2 + (xpq(2, 8)-xpq(2, 1))**2)
endif
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

!...Sound speed at the center...
sdctr = sqrt( max( eps,gamlg*pctr/rhoct) )

!...Characteristic length
dlen = minval(dledgq)!sqrt(geoel(3, ielem))
!
dtcfl = min(dtcfl, dlen/sdctr)
!
!if(ielem.eq.1)then
! print*,'dlen',dlen,dledgq
! print*,'x',xpq(1,1:8)
! print*,'y',xpq(2,1:8)
!endif
!
650 enddo

!...II: Get the time step size limiting the variation of the cell volume
do 750 ie = 1,nquad
ielem = ie + ntria
!print*,'ie',ielem,dtvol,geoel(3, ielem),abs(dvdt(ielem))
dtvol = min(dtvol, geoel(3, ielem)/abs(dvdt(ielem)))
750 enddo

!print*,'dt',cfl*dtcfl, 0.1d0*dtvol, 1.01d0*dtfix

!...III: Get the minimum time tsep size...
if(itime.eq.1)then
dtfix = min(cfl*dtcfl, 0.1d0*dtvol)
else
dtfix = min(cfl*dtcfl, 0.1d0*dtvol, 1.01d0*dtfix)
endif
!
!print*,'dt',0.05d0*dtcfl, 0.1d0*dtvol, 1.01d0*dtfix!,geoel(3,:)

end subroutine getdelt_lag
!
!...Get the edge length for tria or quad
!
subroutine GetLength_edge(xpf, lnged)
use constant
real*8, dimension(2, nvfac), intent(in)::xpf
real*8, intent(out)::lnged

!...Local array
real*8,dimension(1:nvfac)::shp, dshpr
real*8::weigh(ngausf), posi(1,ngausf)

!...Local real or integer
integer::ishp, ig
real*8 ::r, wi,dxdr,dydr,djaco

!...Get position and weight
call rutope(1, ngausf, posi, weigh)

!...Initialize lnged
lnged = 0.d0

!
if(ncurv.eq.0)then
lnged = sqrt((xpf(1, 2)- xpf(1, 1))**2 + (xpf(2, 2)- xpf(2, 1))**2)
elseif(ncurv.eq.1)then

do ig =1,ngausf !...(4)ig=1,ngausf
!
r  = posi(1, ig)
wi  = weigh(ig)
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
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, nvfac
dxdr = dxdr + dshpr(ishp)*xpf(1, ishp)
dydr = dydr + dshpr(ishp)*xpf(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
lnged = lnged + djaco*wi
enddo

else

print*,'Undefined ncurv>1 in subroutine GetLength_edge!'
stop

endif

end subroutine GetLength_edge
