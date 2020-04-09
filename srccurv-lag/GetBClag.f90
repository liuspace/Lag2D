!
!...Periodic boundary condition for 1D isentropic Sin problem...
!
subroutine getbc_prdic(bface, munacn, munacu, snsigm)
use constant
implicit none

!...Input arrays
integer*4,dimension(1:nbfai,nbfac),      intent(in)::bface
real*8, dimension(1:2, 1:2, 1:npoin), intent(inout)::munacn
real*8, dimension(1:ndimn, 1:npoin),  intent(inout)::munacu
real*8, dimension(1:ndimn, 1:npoin),  intent(inout)::snsigm

!...local real array
real*8::muanpr(1:2,1:2,3),muaupr(2, 3),snsgpr(2,3)

!...Local integer array
integer:: iprdic(1:3, 2)

!...Local integer
integer:: ifa, ipr

do ifa = 1, nbfac
!
if(bface(3, ifa).eq.31)then !...Periodic boundary...
!
iprdic(1:2, 1) = bface(1:2, ifa);
iprdic(  3, 1) = bface(  5, ifa);
iprdic(1:2, 2) = bface(1:2, bface(4, ifa))
iprdic(  3, 2) = bface(5, bface(4, ifa))
!
do ipr = 1, 3
!
muanpr(:,:,ipr) = munacn(:,:,iprdic(ipr,1)) + munacn(:,:,iprdic(ipr,2))
muaupr(1:2,ipr) = munacu(1:2, iprdic(ipr,1)) + munacu(1:2, iprdic(ipr,2))
snsgpr(1:2,ipr) = snsigm(1:2, iprdic(ipr,1)) + snsigm(1:2, iprdic(ipr,2))
enddo
!
do ipr = 1, 3
!
munacn(:,:,iprdic(ipr,1)) = muanpr(:,:,ipr); munacn(:,:,iprdic(ipr,2)) = muanpr(:,:, ipr);
!
munacu(1:2,iprdic(ipr,1)) = muaupr(1:2, ipr);
munacu(1:2,iprdic(ipr,2)) = muaupr(1:2, ipr);
!
snsigm(1:2,iprdic(ipr,1)) = snsgpr(1:2, ipr);
snsigm(1:2,iprdic(ipr,2)) = snsgpr(1:2, ipr);
enddo
!
endif
enddo

end subroutine getbc_prdic
!
!...Subroutine to deal with BC with pressure and normal velocity prescribed (Maire) in R-Z (Area-weighted)
!
subroutine getbc_lagmaire_rzaw(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)
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
bsnx = 0.d0
bsny = 0.d0
!
if(bface(3, ifa).eq.22)then !...Normal velocity boundary
!
if(bface(4, ifa).eq.221)then
bsnx = gflag(1,ifa)*gflag(3,ifa)*0.5d0
bsny=  gflag(2,ifa)*gflag(3,ifa)*0.5d0
!bsnx = sqrt(2.d0)/2.d0*gflag(3,ifa)*0.5d0
!bsny =-sqrt(2.d0)/2.d0*gflag(3,ifa)*0.5d0
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
bsnv(1, ipf(1)) = bsnv(1, ipf(1)) + bsnx
bsnv(2, ipf(1)) = bsnv(2, ipf(1)) + bsny

bsnv(1, ipf(2)) = bsnv(1, ipf(2)) + bsnx
bsnv(2, ipf(2)) = bsnv(2, ipf(2)) + bsny
!
elseif(bface(3, ifa).eq.21)then !...Pressure described
!
bsnx = gflag(1,ifa)*gflag(3,ifa)*0.5d0
bsny=  gflag(2,ifa)*gflag(3,ifa)*0.5d0
!
if(ncase.eq.3)then
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny
!
elseif(ncase.eq.4)then
!
dtime = (itime-1.d0)*dtfix + crklg(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prin*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny
!
endif
!
elseif(bface(3, ifa).eq.23)then
!
bsnx =  gflag(1,ifa)*gflag(3,ifa)*0.5d0
bsny =  gflag(2,ifa)*gflag(3,ifa)*0.5d0
!
if(ncase.eq.4)then !...Kidder shell
!
dtime = (itime-1.d0)*dtfix + crklg(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prex*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny
!
endif
!
endif
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
1350 enddo
!
!...Part III: Get the whole force...
!
fpres(:, :) = fpres(:, :) + bsnpf(:, :)
!
end subroutine getbc_lagmaire_rzaw
!
!...Subroutine to deal with Simpson rule BC with prescribed pressure and normal velocity(Maire)
!
subroutine getbccurv_lag(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ngflg,1:nbfac),           intent(inout)::gflag
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
real*8:: bsnx, bsny, prex, prin,anx,any
real*8:: rcoef(2)
real*8:: xph(2, 2)
!
real*8:: farea(nvfac)
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
elseif(nfint.eq.4)then
farea(1) = 1.d0/6.d0
farea(2) = 1.d0/6.d0
farea(3) = 4.d0/6.d0
endif
endif

!
!...Part I: The pressure force for any boundary node...
!
do 1300 ifa = 1 , nbfac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
xph(1, 1:2) = coord(1, ipf(1:2))
xph(2, 1:2) = coord(2, ipf(1:2))
!
anx =   xph(2, 2) - xph(2, 1)
any = -(xph(1, 2) - xph(1, 1))
!
gflag(1, ifa) = anx/sqrt(anx**2 + any**2)
gflag(2, ifa) = any/sqrt(anx**2 + any**2)
gflag(3, ifa) = sqrt(anx**2 + any**2)
!
rcoef(1) = 1.d0 - alfrz + alfrz*xph(2, 1)
rcoef(2) = 1.d0 - alfrz + alfrz*xph(2, 2)
!
bsnx = 0.d0
bsny = 0.d0
!
if(bface(3, ifa).eq.21)then !...Pressure described
!
bsnx = gflag(1,ifa)*gflag(3,ifa)
bsny=  gflag(2,ifa)*gflag(3,ifa)
!
if(ncase.eq.3)then
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx*farea(1)
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny*farea(1)

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx*farea(2)
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny*farea(2)

bsnpf(1, ipf(3)) = bsnpf(1, ipf(3)) + pres0*bsnx*farea(3)
bsnpf(2, ipf(3)) = bsnpf(2, ipf(3)) + pres0*bsny*farea(3)
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
!
!...Part III: Get the whole force...
!
fpres(:, :) = fpres(:, :) + bsnpf(:, :)
!
!print*,'ifa',fpres(1:2,500),bsnpf(1:2, 500)
!
end subroutine getbccurv_lag
!
!...Subroutine to deal with Simpson rule BC with prescribed pressure and normal velocity(Maire)
!
subroutine getbc_lagc(bface, intfac, gflag,  fpres, coord, munacn, munacu, snsigm, itime)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ngflg,1:nbfac),           intent(inout)::gflag
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
real*8:: dxdr, dydr, djaco, r
real*8:: coorp(2, nvfac) !nvfac: no of vertices of one face
integer::ifa, ip, ishp, ig
integer::ipf(nvfac)
real*8:: bsnv(2, npoin), bsnpf(2, npoin)
real*8::munaci(2,2)
real*8::lnx,lny,coefl,coefr,rhsx,rhsy,detma
real*8:: prex, prin,anx,any
real*8:: rcoef(2)
real*8:: xph(2, nvfac), shp(nvfac), dshpr(nvfac),posi(1, nvfac)
real*8:: bsnx(nvfac), bsny(nvfac)
!
real*8:: farea(nvfac)
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
elseif(nfint.eq.4.or.nfint.eq.9)then
farea(1) = 1.d0/6.d0
farea(2) = 1.d0/6.d0
farea(3) = 4.d0/6.d0
!
posi(1, 1) =-1.d0
posi(1, 2) = 1.d0
posi(1, 3) = 0.d0
elseif(nfint.eq.5.or.nfint.eq.10)then
farea(1) = 1.d0
farea(2) = 1.d0
farea(3) = 1.d0
!
posi(1, 1) =-1.d0
posi(1, 2) = 1.d0
posi(1, 3) = 0.d0
!
endif
endif

!
!...Part I: The pressure force for any boundary node...
!
do 1300 ifa = 1 , nbfac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
xph(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xph(2, 1:nvfac) = coord(2, ipf(1:nvfac))
!
!
do ig = 1, nvfac
!
r = posi(1, ig)
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
dxdr = dxdr + dshpr(ishp)*xph(1, ishp)
dydr = dydr + dshpr(ishp)*xph(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
!
bsnx(ig) = dydr*2.d0
bsny(ig) =-dxdr*2.d0
enddo
!
rcoef(1) = 1.d0 - alfrz + alfrz*xph(2, 1)
rcoef(2) = 1.d0 - alfrz + alfrz*xph(2, 2)
!
!bsnx = 0.d0
!bsny = 0.d0
!
if(bface(3, ifa).eq.21)then !...Pressure described
!
!bsnx = gflag(1,ifa)*gflag(3,ifa)
!bsny=  gflag(2,ifa)*gflag(3,ifa)
!
if(ncase.eq.3)then
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx(1)*farea(1)
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny(1)*farea(1)

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx(2)*farea(2)
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny(2)*farea(2)

bsnpf(1, ipf(3)) = bsnpf(1, ipf(3)) + pres0*bsnx(3)*farea(3)
bsnpf(2, ipf(3)) = bsnpf(2, ipf(3)) + pres0*bsny(3)*farea(3)
!
elseif(ncase.eq.4)then
!
dtime = (itime-1.d0)*dtfix + crklg2(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prin*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx(1)*farea(1)
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny(1)*farea(1)

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx(2)*farea(2)
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny(2)*farea(2)

bsnpf(1, ipf(3)) = bsnpf(1, ipf(3)) + pres0*bsnx(3)*farea(3)
bsnpf(2, ipf(3)) = bsnpf(2, ipf(3)) + pres0*bsny(3)*farea(3)
!
endif
!
!
elseif(bface(3, ifa).eq.23)then
!
if(ncase.eq.4)then !...Kidder shell
!
dtime = (itime-1.d0)*dtfix + crklg2(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prex*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx(1)*farea(1)
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny(1)*farea(1)

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx(2)*farea(2)
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny(2)*farea(2)

bsnpf(1, ipf(3)) = bsnpf(1, ipf(3)) + pres0*bsnx(3)*farea(3)
bsnpf(2, ipf(3)) = bsnpf(2, ipf(3)) + pres0*bsny(3)*farea(3)
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
!
!...Part III: Get the whole force...
!
fpres(:, :) = fpres(:, :) + bsnpf(:, :)
!
!print*,'ifa',fpres(1:2,500),bsnpf(1:2, 500)
!
end subroutine getbc_lagc
!
!...Subroutine to deal with high-order quadrature BC with prescribed pressure and normal velocity(Maire)
!
subroutine getbcgauss_lagc(bcindx, vngf, fpres, itime)
use constant
implicit none

!...Input variable...
integer, intent(in):: bcindx, itime
real*8, dimension(3, ngausf),        intent(in)::vngf
real*8,dimension(1:ndimn,1:ngausf), intent(out)::fpres

!
!...Local variable...
!
real*8:: pres0
real*8:: prin, prex
real*8:: dtime
!
integer:: ig
!
!...Part I: Zero our fpres and pres0...
!
fpres = 0.d0
pres0 = 1.d-6
!
!...Part II: Impose pressure boundary...
!
do ig = nvfac+1, ngausf
!
if(bcindx.eq.21)then !...Pressure described
!
if(ncase.eq.3)then
!
fpres(1, ig) =   pres0*vngf(1, ig)*vngf(3, ig)
fpres(2, ig) =   pres0*vngf(2, ig)*vngf(3, ig)
!
elseif(ncase.eq.4)then
!
dtime = (itime-1.d0)*dtfix + crklg2(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prin*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
fpres(1, ig) =   pres0*vngf(1, ig)*vngf(3, ig)
fpres(2, ig) =   pres0*vngf(2, ig)*vngf(3, ig)
!
endif
!
elseif(bcindx.eq.23)then
!
if(ncase.eq.4)then !...Kidder shell
!
dtime = (itime-1.d0)*dtfix + crklg2(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prex*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
fpres(1, ig) =   pres0*vngf(1, ig)*vngf(3, ig)
fpres(2, ig) =   pres0*vngf(2, ig)*vngf(3, ig)

endif
!
endif
!
enddo
!
end subroutine getbcgauss_lagc
!
!...Subroutine to deal with Simpson rule BC with prescribed pressure and normal velocity suing general methods
!
subroutine getbc_lagc_general(bface, intfac,  fpres, coord, munacn, munacu, snsigm, itime)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),          intent(out)::fpres
real*8,intent(in)::munacn(1:2, 1:2, 1:npoin)
real*8,intent(in)::snsigm(1:ndimn, 1:npoin)
real*8,intent(in)::munacu(1:ndimn, 1:npoin)
integer, intent(in):: itime
!
!...Local variable...
!
real*8:: pres0
real*8:: dtime, htkid, taokid
real*8:: dxdr, dydr, djaco, r
integer::ifa, ip, ishp, ig
integer::ipf(nvfac)
real*8:: bsnv(2, npoin), bsnpf(2, npoin)
real*8::munaci(2,2)
real*8::lnx,lny,coefl,coefr,rhsx,rhsy,detma
real*8:: prex, prin
real*8:: rcoef(2)
real*8:: xph(2, nvfac), shp(nvfac), dshpr(nvfac),posi(1, nvfac)
real*8:: bsnx(nvfac), bsny(nvfac)
!
real*8:: farea(nvfac)
!
!...Zero our fpres and pres0...
!
fpres = 0.d0
pres0 = 1.d-14
!
bsnv = 0.d0  !...area*n for velocity boundary
bsnpf = 0.d0 !...force for pressure boundary
!
if(nriem.eq.1)then
if(nfint.eq.1)then
farea = .5d0
elseif(nfint.eq.3)then
farea = 1.d0
elseif(nfint.eq.4.or.nfint.eq.9)then
farea(1) = 1.d0/6.d0
farea(2) = 1.d0/6.d0
farea(3) = 4.d0/6.d0
!
posi(1, 1) =-1.d0
posi(1, 2) = 1.d0
posi(1, 3) = 0.d0
elseif(nfint.eq.5.or.nfint.eq.10)then
farea(1) = 1.d0
farea(2) = 1.d0
farea(3) = 1.d0
!
posi(1, 1) =-1.d0
posi(1, 2) = 1.d0
posi(1, 3) = 0.d0
!
endif
endif

!
!...Part I: The pressure force for any boundary node...
!
do 1300 ifa = 1 , nbfac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
xph(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xph(2, 1:nvfac) = coord(2, ipf(1:nvfac))
!
!
do ig = 1, nvfac
!
r = posi(1, ig)
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
dxdr = dxdr + dshpr(ishp)*xph(1, ishp)
dydr = dydr + dshpr(ishp)*xph(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
!
bsnx(ig) = dydr*2.d0
bsny(ig) =-dxdr*2.d0
enddo
!
rcoef(1) = 1.d0 - alfrz + alfrz*xph(2, 1)
rcoef(2) = 1.d0 - alfrz + alfrz*xph(2, 2)

!
if(bface(3, ifa).eq.21)then !...Pressure described
!
if(ncase.eq.3)then
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx(1)*farea(1)
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny(1)*farea(1)

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx(2)*farea(2)
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny(2)*farea(2)

bsnpf(1, ipf(3)) = bsnpf(1, ipf(3)) + pres0*bsnx(3)*farea(3)
bsnpf(2, ipf(3)) = bsnpf(2, ipf(3)) + pres0*bsny(3)*farea(3)
!
elseif(ncase.eq.4)then
!
dtime = (itime-1.d0)*dtfix + crklg2(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prin*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx(1)*farea(1)
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny(1)*farea(1)

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx(2)*farea(2)
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny(2)*farea(2)

bsnpf(1, ipf(3)) = bsnpf(1, ipf(3)) + pres0*bsnx(3)*farea(3)
bsnpf(2, ipf(3)) = bsnpf(2, ipf(3)) + pres0*bsny(3)*farea(3)
!
endif
!
!
elseif(bface(3, ifa).eq.23)then
!
if(ncase.eq.4)then !...Kidder shell
!
dtime = (itime-1.d0)*dtfix + crklg2(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prex*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx(1)*farea(1)
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny(1)*farea(1)

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx(2)*farea(2)
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny(2)*farea(2)

bsnpf(1, ipf(3)) = bsnpf(1, ipf(3)) + pres0*bsnx(3)*farea(3)
bsnpf(2, ipf(3)) = bsnpf(2, ipf(3)) + pres0*bsny(3)*farea(3)
endif
!
elseif(bface(3, ifa).eq.22)then
!
if(bface(4, ifa).eq.221)then
!
bsnv(1, ipf(1)) = bsnv(1, ipf(1)) +  0.d0
bsnv(2, ipf(1)) = bsnv(2, ipf(1)) -1.d0*bsny(1)*farea(1)
!
bsnv(1, ipf(2)) = bsnv(1, ipf(2)) +  0.d0
bsnv(2, ipf(2)) = bsnv(2, ipf(2)) -1.d0*bsny(2)*farea(2)
!
bsnv(1, ipf(3)) = bsnv(1, ipf(3)) +  0.d0
bsnv(2, ipf(3)) = bsnv(2, ipf(3)) -1.d0*bsny(3)*farea(3)

elseif(bface(4, ifa).eq.222)then
!
bsnv(1, ipf(1)) = bsnv(1, ipf(1)) -1.d0*bsnx(1)*farea(1)
bsnv(2, ipf(1)) = bsnv(2, ipf(1)) +0.d0
!
bsnv(1, ipf(2)) = bsnv(1, ipf(2)) -1.d0*bsnx(2)*farea(2)
bsnv(2, ipf(2)) = bsnv(2, ipf(2)) +0.d0
!
bsnv(1, ipf(3)) = bsnv(1, ipf(3)) -1.d0*bsnx(3)*farea(3)
bsnv(2, ipf(3)) = bsnv(2, ipf(3)) +0.d0
endif
!
!
elseif(bface(3, ifa).eq.29.and.nmatel.eq.2)then !...Pressure described
!
if(ncase.eq.2)then
!
bsnpf(1, ipf(1)) = bsnpf(1, ipf(1)) + pres0*bsnx(1)*farea(1)
bsnpf(2, ipf(1)) = bsnpf(2, ipf(1)) + pres0*bsny(1)*farea(1)

bsnpf(1, ipf(2)) = bsnpf(1, ipf(2)) + pres0*bsnx(2)*farea(2)
bsnpf(2, ipf(2)) = bsnpf(2, ipf(2)) + pres0*bsny(2)*farea(2)

bsnpf(1, ipf(3)) = bsnpf(1, ipf(3)) + pres0*bsnx(3)*farea(3)
bsnpf(2, ipf(3)) = bsnpf(2, ipf(3)) + pres0*bsny(3)*farea(3)
!
!
endif
!
endif
!
!print*,'ifa',ipf(1:2),ifa
!
1300 enddo
!
!...Part II: Get the pressure from the normal velocity prescription
!
do 1350 ifa = 1 , nbfac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
if(bface(3, ifa).eq.22)then
!
do ip = 1, nvfac
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
1350 enddo
!
!print*,'pres'
!
!...Part III: Get the whole force...
!
fpres(:, :) = fpres(:, :) + bsnpf(:, :)
!
!print*,'ifa',fpres(1:2,500),bsnpf(1:2, 500)
!
end subroutine getbc_lagc_general
!
!...Subroutine to deal with Simpson rule BC with prescribed pressure and normal velocity suing general methods
!
subroutine getbc_lagc_general_curv(bface, intfac,  fpres, coord, munacn, munacu, snsigm, itime)
use constant
implicit none
!
!...Input variable...
!
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ndimn,1:npoin),          intent(out)::fpres
real*8,intent(in)::munacn(1:2, 1:2, 1:npoin)
real*8,intent(in)::snsigm(1:ndimn, 1:npoin)
real*8,intent(in)::munacu(1:ndimn, 1:npoin)
integer, intent(in):: itime
!
!...Local variable...
!
real*8:: pres0
real*8:: dtime, htkid, taokid
real*8:: dxdr, dydr, djaco, r
integer::ifa, ip, ishp, ig
integer::ipf(nvfac)
real*8:: bsnv(2, npoin), bsnpf(2, npoin)
real*8::munaci(2,2)
real*8::lnx,lny,coefl,coefr,rhsx,rhsy,detma
real*8:: prex, prin
real*8:: rcoef(2)
real*8:: xph(2, nvfac), shpf(nvfac), dshprf(nvfac),posi(1, nvfac)
real*8:: bsnx(nvfac), bsny(nvfac)
!
real*8:: farea(nvfac)
!
!...Zero our fpres and pres0...
!
fpres = 0.d0
pres0 = 1.d-6
!
bsnv = 0.d0  !...area*n for velocity boundary
bsnpf = 0.d0 !...force for pressure boundary
!
if(ncurv.le.1)then

if(nriem.eq.1)then
if(nfint.eq.1)then
farea = .5d0
elseif(nfint.eq.3)then
farea = 1.d0
elseif(nfint.eq.4.or.nfint.eq.9)then
farea(1) = 1.d0/6.d0
farea(2) = 1.d0/6.d0
farea(3) = 4.d0/6.d0
!
posi(1, 1) =-1.d0
posi(1, 2) = 1.d0
posi(1, 3) = 0.d0
elseif(nfint.eq.5.or.nfint.eq.10)then
farea(1) = 1.d0
farea(2) = 1.d0
farea(3) = 1.d0
!
posi(1, 1) =-1.d0
posi(1, 2) = 1.d0
posi(1, 3) = 0.d0
endif
endif

elseif(ncurv.eq.2)then
farea(1) = 1.d0
farea(2) = 1.d0
farea(3) = 2.d0
farea(4) = 2.d0
!!
posi(1, 1) =-1.d0
posi(1, 2) = 1.d0
posi(1, 3) =-sqrt(5.d0)/5.d0
posi(1, 4) = sqrt(5.d0)/5.d0
endif
!
!...Part I: The pressure force for any boundary node...
!
do 1300 ifa = 1 , nbfac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
xph(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xph(2, 1:nvfac) = coord(2, ipf(1:nvfac))
!
!
do ig = 1, nvfac
!
r = posi(1, ig)
!
call getshapfct_edge(ncurv,nvfac,shpf, dshprf, r)
!
!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, nvfac
dxdr = dxdr + dshprf(ishp)*xph(1, ishp)
dydr = dydr + dshprf(ishp)*xph(2, ishp)
enddo
!
djaco = sqrt(dxdr**2 + dydr**2)
!
bsnx(ig) = dydr*2.d0
bsny(ig) =-dxdr*2.d0
enddo
!
rcoef(1) = 1.d0 - alfrz + alfrz*xph(2, 1)
rcoef(2) = 1.d0 - alfrz + alfrz*xph(2, 2)

!
if(bface(3, ifa).eq.21)then !...Pressure described
!
if(ncase.eq.3)then
!
do ig = 1, nvfac
bsnpf(1, ipf(ig)) = bsnpf(1, ipf(ig)) + pres0*bsnx(ig)*farea(ig)
bsnpf(2, ipf(ig)) = bsnpf(2, ipf(ig)) + pres0*bsny(ig)*farea(ig)
enddo
!
elseif(ncase.eq.4)then
!
dtime = (itime-1.d0)*dtfix + crklg2(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prin*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
do ig = 1, nvfac
bsnpf(1, ipf(ig)) = bsnpf(1, ipf(ig)) + pres0*bsnx(ig)*farea(ig)
bsnpf(2, ipf(ig)) = bsnpf(2, ipf(ig)) + pres0*bsny(ig)*farea(ig)
enddo
!
endif
!
!
elseif(bface(3, ifa).eq.23)then
!
if(ncase.eq.4)then !...Kidder shell
!
dtime = (itime-1.d0)*dtfix + crklg2(rkstg)*dtfix

prin = 0.1d0
prex = 10.d0
!
pres0 = prex*sqrt(1.d0-dtime**2)**(-2.d0*gamlg/(gamlg-1.d0))
!
do ig = 1, nvfac
bsnpf(1, ipf(ig)) = bsnpf(1, ipf(ig)) + pres0*bsnx(ig)*farea(ig)
bsnpf(2, ipf(ig)) = bsnpf(2, ipf(ig)) + pres0*bsny(ig)*farea(ig)
enddo

endif
!
elseif(bface(3, ifa).eq.22)then
!
if(bface(4, ifa).eq.221)then
!
do ig = 1, nvfac
bsnv(1, ipf(ig)) = bsnv(1, ipf(ig)) + 0.d0
bsnv(2, ipf(ig)) = bsnv(2, ipf(ig)) - 1.d0*bsny(ig)*farea(ig)
enddo

elseif(bface(4, ifa).eq.222)then
!
do ig = 1, nvfac
bsnv(1, ipf(ig)) = bsnv(1, ipf(ig)) -1.d0*bsnx(ig)*farea(ig)
bsnv(2, ipf(ig)) = bsnv(2, ipf(ig)) + 0.d0
enddo

endif
!
!
elseif(bface(3, ifa).eq.29.and.nmatel.eq.2)then !...Pressure described
!
if(ncase.eq.2)then
!
do ig = 1, nvfac
bsnpf(1, ipf(ig)) = bsnpf(1, ipf(ig)) + pres0*bsnx(ig)*farea(ig)
bsnpf(2, ipf(ig)) = bsnpf(2, ipf(ig)) + pres0*bsny(ig)*farea(ig)
enddo
!
endif
!
endif
!
!print*,'ifa',ipf(1:2),ifa
!
1300 enddo
!
!...Part II: Get the pressure from the normal velocity prescription
!
do 1350 ifa = 1 , nbfac
!
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
!
if(bface(3, ifa).eq.22)then
!
do ip = 1, nvfac
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
1350 enddo
!
!print*,'pres'
!
!...Part III: Get the whole force...
!
fpres(:, :) = fpres(:, :) + bsnpf(:, :)
!
!print*,'ifa',fpres(1:2,500),bsnpf(1:2, 500)
!
end subroutine getbc_lagc_general_curv
