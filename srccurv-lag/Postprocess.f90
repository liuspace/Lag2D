!
!...subroutine: post processing...
!
subroutine getpostp1(unkno,inpoel,coord,geoel)
use constant
implicit none
real*8,dimension(1:ndegr,1:nq,1:nsize)::unkno
integer*4,dimension(1:nvtri,1:ncell)::inpoel
real*8,dimension(1:ndimn,1:npoin)::coord
real*8,dimension(1:ngeel,1:nsize)::geoel
integer*4::i,ie,ip1,ip2,ip3,id,ideg
real*8::xcent,ycent,dist,invdis,vxave,vyave,weight,pave
real*8::ux,uy,re,ua,presi,aspd,dx,dy,xc,yc,xg,yg,ui,dpi,qi,qave
real*8::b(6)
real*8::rhoi,rhoui,rhovi,rhoei
real*8,dimension(1:npoin)::pres,vxpoin,vypoin
real*8,dimension(1:npoin)::rho,rhou,rhov
!
 id=1
!
 rho = 0.d0
 rhou = 0.d0
 rhov = 0.d0
!calculate the point value
  do i=1,npoin
    weight=0.d0
    vxave=0.d0
    qave=0.d0
    pave=0.d0
   do ie=1,ncell
    ip1=inpoel(1,ie)
    ip2=inpoel(2,ie)
    ip3=inpoel(3,ie)
!
   dx = geoel(4,ie)*0.5d0
   dy = geoel(5,ie)*0.5d0
!
   if((ip1.eq.i).or.(ip2.eq.i).or.(ip3.eq.i))then
    xcent=(coord(1,ip1)+coord(1,ip2)+coord(1,ip3))/3.d0
    ycent=(coord(2,ip1)+coord(2,ip2)+coord(2,ip3))/3.d0
    dist=sqrt((coord(1,i)-xcent)**2+(coord(2,i)-ycent)**2)
    invdis=1.d0/(dist)
!   
    xc = xcent
    yc = ycent
!
    xg = coord(1,i)
    yg = coord(2,i)
!
     b(1) = 1.d0
     b(2) = (xg-xc)/dx
     b(3) = (yg-yc)/dy
     b(4) = 0.5d0*b(2)**2-geoel(11,ie)
     b(5) = 0.5d0*b(3)**2-geoel(12,ie)
     b(6) = b(2)*b(3) - geoel(13,ie)
!
    rhoi  = 0.d0
    rhoui = 0.d0
    rhovi = 0.d0
    rhoei = 0.d0
!
   do ideg = 1,ndegr
    rhoi = rhoi  + unkno(ideg,1,ie)*b(ideg)
    rhoui= rhoui + unkno(ideg,2,ie)*b(ideg)
    rhovi= rhovi + unkno(ideg,3,ie)*b(ideg)
    rhoei= rhoei + unkno(ideg,4,ie)*b(ideg)
   enddo
!
    rho(i)=rho(i)+rhoi*invdis
    rhou(i)=rhou(i)+rhoui/rhoi*invdis
    rhov(i)=rhov(i)+rhovi/rhoi*invdis
    weight=invdis+weight
   endif
   enddo  
   rho(i) =rho(i)/weight
   rhou(i)=rhou(i)/weight
   rhov(i)=rhov(i)/weight
  enddo
 !output the value
   call getoutput(inpoel,coord,rho,rhou,rhov)

 return
! close(11)
end subroutine getpostp1
!
!...subroutine: output tecplot format...
!
subroutine getoutput(inpoel,coord,pres,vxpoin,vypoin)
use constant
integer*4,dimension(1:nvtri,1:nelem)::inpoel
real*8,dimension(1:ndimn,1:npoin)::coord
real*8,dimension(1:npoin)::pres,vxpoin,vypoin
real*8::ux,uy
integer*4::ip

  open(11,file='output.dat')
  !  print*,'hi'
  write(11,*)'TITLE="Example"'
  write(11,*)'VARIABLES="X","Y","<greek><b>f</b></greek>","u","v"'
  write(11,*)'ZONE T="P_1", DATAPACKING=BLOCK, NODES=',npoin,',ELEMENTS=',nelem, ',ZONETYPE=FETRIANGLE'
  write(11,*)coord(1,1:npoin)
  write(11,*)
  write(11,*)coord(2,1:npoin)
  write(11,*)
  write(11,*)pres(1:npoin)
  write(11,*)
  write(11,*)vxpoin(1:npoin)
  write(11,*)
  write(11,*)vypoin(1:npoin)
  write(11,*)
  do ie=1,nelem
    write(11,*)inpoel(1:3,ie)
  enddo
end subroutine getoutput
!
!...Subroutine: Check the mesh quality (triangle)...
!
subroutine getoutput_cell(inpoel,coord,geoel,unkno)
use constant
implicit none
!...Real array
integer*4,dimension(1:nvtri,1:nelem)::inpoel
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unkno
!...Integer array
real*8,dimension(1:ndimn,1:npoin), intent(in)::coord
real*8,dimension(1:ngeel,1:nsize), intent(in)::geoel
!
!...Local
integer:: ielem, ip, ie
real*8 :: eps
real*8 :: rho, rhou, rhov, rhoe, uadv, vadv, rhom
!
!...Real array
real*8, allocatable::pres(:), mach(:)
!
  data   eps / 1.0d-06 /
!
!...Step 1: Find pressure, Mach...
!
  allocate(pres(1:ncell), mach(1:ncell))
!
  do ielem = 1,ncell
!
      rho  = unkno(1, 1, ielem)
      rhou = unkno(1, 2, ielem)
      rhov = unkno(1, 3, ielem)
      rhoe = unkno(1, 4, ielem)
!
      rhom = 1.d0/rho
      uadv = rhou*rhom
      vadv = rhov*rhom
     if(nmeth==1)then
      pres(ielem) = (gamma-1.d0)*(rhoe-0.5d0*rho*(uadv*uadv+vadv*vadv)) 
     elseif(nmeth==2)then
      pres(ielem) = (gamlg-1.d0)*(rhoe-0.5d0*rho*(uadv*uadv+vadv*vadv)) 
     endif  
      mach(ielem) = sqrt(uadv*uadv+vadv*vadv)/sqrt(max(eps,gamma*pres(ielem)/rho))
  enddo 
!
!...Step 2: Output tecplot format...
!
  open(11,file='output_cell.dat')
  write(11,*)'VARIABLES = "X","Y","Density","U","V","P","Ma",'
  write(11,*)'ZONE'
  write(11,*)'     T="2D"'
  write(11,*)'     DATAPACKING=BLOCK'
  write(11,*)'     ZONETYPE=FEtriangle'
  write(11,*)'     N=',npoin
  write(11,*)'     E=',nelem
  write(11,*)'     VARLOCATION=([3-7]=CELLCENTERED)'
!
!...X, Y coordinates...
  do ip =1,npoin
     write(11,*)coord(1,ip)
  enddo
  do ip=1,npoin
    write(11,*)coord(2,ip)
  enddo
!
!...Density
  do ie =1,nelem
    write(11,*)unkno(1, 1, ie)
  enddo
!
!...U
  do ie =1,nelem
    write(11,*)unkno(1, 1, ie)/unkno(1, 1, ie)
  enddo
!
!...V	
  do ie =1,nelem
    write(11,*)unkno(1, 2, ie)/unkno(1, 1, ie)
  enddo
!
!...Pressure
  do ie =1,nelem
    write(11,*)pres(ie)!pres(ie)
  enddo
!
!...Mach
  do ie =1,nelem
    write(11,*)mach(ie)!mach(ie)
  enddo
    write(11,*)
  do ie=1,nelem
    write(11,*)inpoel(1:3,ie) 
  enddo
 close(11)
!
!...
  deallocate(pres, mach)
!
end subroutine getoutput_cell
!
!...Find face center in geoel...
!
 subroutine geterror(inpoel, iptri, ipqua, geoel, coord, unkno)
 use constant
 implicit none
 integer*4,dimension(1:nvtri, 1:ntria),intent(in)::inpoel
 integer*4,dimension(1:nvtri,1:ntria)::iptri
 integer*4,dimension(1:nvqua,1:nquad)::ipqua
!
  real*8,dimension(1:ngeel, 1:ncell), intent(in)::geoel
  real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
  real*8,dimension(1:ndegr,1:nq,1:nsize),  intent(in)::unkno
!
!...local array...
!
  real*8,dimension(1:2, 1:3)::xpin
  real*8,dimension(1:2, 1:nptri)::xp
  real*8,dimension(1:2, 1:npqua)::xpq
  real*8,dimension(1:nptri)::shp, dspr, dsps
  real*8,dimension(1:npqua)::shpq, dsprq, dspsq
  real*8:: weight(ngausd), posit(2, ngausd)
  real*8:: weighq(ngausdq), posiq(2, ngausdq)
  real*8::b(6)
!...local real number
  real*8:: dxdr,dxds,dydr,dyds 
  real*8:: eps,c00,c10,c05,c20
  real*8:: r, s, djaco, volel
  real*8:: wi, xc, yc, xg, yg, dx, dy
  real*8::rho, rhou, rhov, rhoe, pres, difun, errl2
!
  integer::ielem, igaus, ishp, ideg
  integer:: ie
!
  data eps / 1.0d-06 / 
  data c00 / 0.0d0 / 
  data c10 / 1.0d0 / 
  data c05 / 0.5d0 / 
  data c20 / 2.0d0 / 
!
  errl2 = 0.d0
!
!...Find weight and position for gauss points...
  call rutope(2, ngausd, posit, weight)
  call ruqope(2, ngausdq, posiq, weighq)
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
  endif
!
    xc = geoel(1, ielem)
    yc = geoel(2, ielem)

    dx = geoel(4, ielem)*0.5d0
    dy = geoel(5, ielem)*0.5d0

    volel = geoel(3, ielem)
!
  do igaus =1,ngausd
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
!
     xg = 0.d0
     yg = 0.d0     
!
     do ishp = 1, nptri  
      xg = xg + shp(ishp)*xp(1,ishp) 
      yg = yg + shp(ishp)*xp(2,ishp) 
     enddo  
!
     b(1) = 1.d0
     b(2) = (xg-xc)/dx
     b(3) = (yg-yc)/dy
     b(4) = 0.5d0*b(2)**2-geoel(6,ielem)
     b(5) = 0.5d0*b(3)**2-geoel(7,ielem)
     b(6) = b(2)*b(3) - geoel(8,ielem)
!
     rho = 0.d0; rhou = 0.d0; rhov = 0.d0; rhoe = 0.d0;
!     
     do ideg = 1, ndegr  
      rho = rho + unkno(ideg, 1, ielem)*b(ideg)
      rhou= rhou+ unkno(ideg, 2, ielem)*b(ideg)
      rhov= rhov+ unkno(ideg, 3, ielem)*b(ideg)
      rhoe= rhoe+ unkno(ideg, 4, ielem)*b(ideg)
     enddo  
!
     pres = (gamma-1.d0)*(rhoe-0.5d0*(rhou*rhou+rhov*rhov)/rho)
     difun = (pres*gamma)*(1.d0/rho)**gamma - 1.d0
     errl2 = errl2 + difun*difun*djaco
! 
  enddo
 enddo !...(1)ie = 1,nelem
!
!...2nd loop for quads..
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
endif!
!
xc = geoel(1, ielem)
yc = geoel(2, ielem)

dx = geoel(4, ielem)*0.5d0
dy = geoel(5, ielem)*0.5d0

volel = geoel(3, ielem)
!
do igaus =1,ngausdq
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
xg = 0.d0
yg = 0.d0
!
do ishp = 1, npqua
xg = xg + shpq(ishp)*xpq(1,ishp)
yg = yg + shpq(ishp)*xpq(2,ishp)
enddo
!
b(1) = 1.d0
b(2) = (xg-xc)/dx
b(3) = (yg-yc)/dy
b(4) = 0.5d0*b(2)**2-geoel(6,ielem)
b(5) = 0.5d0*b(3)**2-geoel(7,ielem)
b(6) = b(2)*b(3) - geoel(8,ielem)
!
rho = 0.d0; rhou = 0.d0; rhov = 0.d0; rhoe = 0.d0;
!
do ideg = 1, ndegr
rho = rho + unkno(ideg, 1, ielem)*b(ideg)
rhou= rhou+ unkno(ideg, 2, ielem)*b(ideg)
rhov= rhov+ unkno(ideg, 3, ielem)*b(ideg)
rhoe= rhoe+ unkno(ideg, 4, ielem)*b(ideg)
enddo
!
pres = (gamma-1.d0)*(rhoe-0.5d0*(rhou*rhou+rhov*rhov)/rho)
difun = (pres*gamma)*(1.d0/rho)**gamma - 1.d0
errl2 = errl2 + difun*difun*djaco
!
enddo
enddo !...(1)ie = 1,nelem
!
    print*,'L2 error',0.5d0*log10(errl2),0.5d0*log10(100.0)
!

 end subroutine geterror
!
!...Find face center in geoel...
!
 subroutine geterror_curv2(inpoel, geoel, coord, unkno)
 use constant
 implicit none
 integer*4,dimension(1:nvtri, 1:nelem),intent(in)::inpoel
!
  real*8,dimension(1:ngeel, 1:nelem), intent(in)::geoel
  real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
  real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),  intent(in)::unkno
!
!...local array...
!
  real*8,dimension(1:2, 1:3)::xpin
  real*8,dimension(1:2, 1:nptri)::xp 
  real*8,dimension(1:nptri)::shp, dspr, dsps
  real*8:: weigh(ngausd), posi(2, ngausd)
  real*8::b(6)
!...local real number
  real*8:: dxdr,dxds,dydr,dyds 
  real*8:: eps,c00,c10,c05,c20
  real*8:: r, s, t, djaco, volel
  real*8:: wi, xc, yc, xg, yg, dx, dy
  real*8::rho, rhou, rhov, rhoe, pres, difun, errl2
!
  integer::ielem, igaus, ishp, ideg
!
  data eps / 1.0d-06 / 
  data c00 / 0.0d0 / 
  data c10 / 1.0d0 / 
  data c05 / 0.5d0 / 
  data c20 / 2.0d0 / 
!
  errl2 = 0.d0
!
!...Find weight and position for gauss points...
  call rutope(2, ngausd, posi, weigh)
!
!...1st loop to find center and volulme...
!
 do ielem = 1, nelem !...(1)ie = 1,nelem
!
   if(ncurv==0)then
    xp(1, 1:3) = coord(1, inpoel(1:3,ielem))
    xp(2, 1:3) = coord(2, inpoel(1:3,ielem))
!
    xp(1:2,4) = 0.5d0*(xp(1:2,1)+xp(1:2,2))
    xp(1:2,5) = 0.5d0*(xp(1:2,2)+xp(1:2,3))
    xp(1:2,6) = 0.5d0*(xp(1:2,1)+xp(1:2,3))    
   elseif(ncurv==2)then
    xp(1, 1:nptri) = coord(1,inpoel(1:nptri, ielem))
    xp(2, 1:nptri) = coord(2,inpoel(1:nptri, ielem))
   endif
!
    xc = geoel(1, ielem)
    yc = geoel(2, ielem)

    dx = geoel(4, ielem)*0.5d0
    dy = geoel(5, ielem)*0.5d0

    volel = geoel(3, ielem)
!
  do igaus =1,ngausd
!
     r  = posi(1,igaus)
     s  = posi(2,igaus)
     t  = 1.d0-r-s
    wi  = weigh(igaus)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
     shp(1) =       t*(-1.d0+3.d0*t)*(-2.d0+3.d0*t)
     shp(2) =       r*(-1.d0+3.d0*r)*(-2.d0+3.d0*r)
     shp(3) =       s*(-1.d0+3.d0*s)*(-2.d0+3.d0*s)
     shp(4) = 9.d0*t*r*(-1.d0+3.d0*t)
     shp(5) = 9.d0*r*s*(-1.d0+3.d0*r)
     shp(6) = 9.d0*t*s*(-1.d0+3.d0*s)
     shp(7) = 9.d0*t*r*(-1.d0+3.d0*r)
     shp(8) = 9.d0*r*s*(-1.d0+3.d0*s)
     shp(9) = 9.d0*t*s*(-1.d0+3.d0*t)
     shp(10)=54.d0*r*t*s 
!
     shp(:) = 0.5d0*shp(:)
    !
    !
     dspr(1) = -2.d0+18.d0*t-27.d0*t**2
     dspr(2) =  2.d0-18.d0*r+27.d0*r**2  
     dspr(3) =  0.d0
     dspr(4) =  9.d0*t*(-1.d0+3.d0*t-6.d0*r)+9.d0*r 
     dspr(5) =  9.d0*s*(-1.d0+6.d0*r)
     dspr(6) = -9.d0*s*(-1.d0+3.d0*s)
     dspr(7) =  9.d0*r*(1.d0+6.d0*t-3.d0*r)-9.d0*t 
     dspr(8) =  9.d0*s*(-1.d0+3.d0*s)
     dspr(9) = -9.d0*s*(-1.d0+6.d0*t)
     dspr(10)= 54.d0*s*(t-r)
    ! 
     dspr = dspr*0.5d0
     !
     dsps(1) = -2.d0+18.d0*t-27.d0*t**2
     dsps(2) =  0.d0 
     dsps(3) =  2.d0-18.d0*s+27.d0*s**2
     dsps(4) = -9.d0*r*(-1.d0+6.d0*t) 
     dsps(5) =  9.d0*r*(-1.d0+3.d0*r)
     dsps(6) =  9.d0*s*(1.d0+6.d0*t-3.d0*s)-9.d0*t 
     dsps(7) = -9.d0*r*(-1.d0+3.d0*r)
     dsps(8) =  9.d0*r*(-1.d0+6.d0*s) 
     dsps(9) =  9.d0*t*(-1.d0+3.d0*t-6.d0*s)+9.d0*s  
     dsps(10)= 54.d0*r*(t-s) 
   !
     dsps = dsps*0.5d0
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
     xg = 0.d0
     yg = 0.d0     
!
     do ishp = 1, nptri  
      xg = xg + shp(ishp)*xp(1,ishp) 
      yg = yg + shp(ishp)*xp(2,ishp) 
     enddo  
!
     b(1) = 1.d0
     b(2) = (xg-xc)/dx
     b(3) = (yg-yc)/dy
     b(4) = 0.5d0*b(2)**2-geoel(6,ielem)
     b(5) = 0.5d0*b(3)**2-geoel(7,ielem)
     b(6) = b(2)*b(3) - geoel(8,ielem)
!
     rho = 0.d0; rhou = 0.d0; rhov = 0.d0; rhoe = 0.d0;
!     
     do ideg = 1, ndegr  
      rho = rho + unkno(ideg, 1, ielem)*b(ideg)
      rhou= rhou+ unkno(ideg, 2, ielem)*b(ideg)
      rhov= rhov+ unkno(ideg, 3, ielem)*b(ideg)
      rhoe= rhoe+ unkno(ideg, 4, ielem)*b(ideg)
     enddo  
!
     pres = (gamma-1.d0)*(rhoe-0.5d0*(rhou*rhou+rhov*rhov)/rho)
     difun = (pres*gamma)*(1.d0/rho)**gamma - 1.d0
     errl2 = errl2 + difun*difun*djaco
! 
  enddo
 enddo !...(1)ie = 1,nelem
!
!
    print*,'L2 error',0.5d0*log10(errl2),0.5d0*log10(100.0)
!
 end subroutine geterror_curv2
!
!...L2 Error called during time loop... 
!
 subroutine geterror_resi(inpoel, geoel, coord, unkno, errl2)
 use constant
 implicit none
 integer*4,dimension(1:nvtri, 1:nelem),intent(in)::inpoel
!
  real*8,dimension(1:ngeel, 1:nelem), intent(in)::geoel
  real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
  real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),  intent(in)::unkno
!
  real*8, intent(out):: errl2
!
!...local array...
!
  real*8,dimension(1:2, 1:3)::xpin
  real*8,dimension(1:2, 1:nptri)::xp 
  real*8,dimension(1:nptri)::shp, dspr, dsps
  real*8:: weigh(ngausd), posi(2, ngausd)
  real*8::b(6)
!...local real number
  real*8:: dxdr,dxds,dydr,dyds 
  real*8:: eps,c00,c10,c05,c20
  real*8:: r, s, djaco, volel
  real*8:: wi, xc, yc, xg, yg, dx, dy
  real*8::rho, rhou, rhov, rhoe, pres, difun
!
  integer::ielem, igaus, ishp, ideg
!
  data eps / 1.0d-06 / 
  data c00 / 0.0d0 / 
  data c10 / 1.0d0 / 
  data c05 / 0.5d0 / 
  data c20 / 2.0d0 / 
!
  errl2 = 0.d0
!
!...Find weight and position for gauss points...
  call rutope(2, ngausd, posi, weigh)
!
!...1st loop to find center and volulme...
!
 do ielem = 1, nelem !...(1)ie = 1,nelem
!
   if(ncurv==0)then
    xp(1, 1:3) = coord(1, inpoel(1:3,ielem))
    xp(2, 1:3) = coord(2, inpoel(1:3,ielem))
!
    xp(1:2,4) = 0.5d0*(xp(1:2,1)+xp(1:2,2))
    xp(1:2,5) = 0.5d0*(xp(1:2,2)+xp(1:2,3))
    xp(1:2,6) = 0.5d0*(xp(1:2,1)+xp(1:2,3))    
   elseif(ncurv==1)then
    xp(1, 1:nptri) = coord(1,inpoel(1:nptri, ielem))
    xp(2, 1:nptri) = coord(2,inpoel(1:nptri, ielem))
   endif
!
    xc = geoel(1, ielem)
    yc = geoel(2, ielem)

    dx = geoel(4, ielem)*0.5d0
    dy = geoel(5, ielem)*0.5d0

    volel = geoel(3, ielem)
!
  do igaus =1,ngausd
!
     r  = posi(1,igaus)
     s  = posi(2,igaus)
    wi  = weigh(igaus)
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
     xg = 0.d0
     yg = 0.d0     
!
     do ishp = 1, nptri  
      xg = xg + shp(ishp)*xp(1,ishp) 
      yg = yg + shp(ishp)*xp(2,ishp) 
     enddo  
!
     b(1) = 1.d0
     b(2) = (xg-xc)/dx
     b(3) = (yg-yc)/dy
     b(4) = 0.5d0*b(2)**2-geoel(6,ielem)
     b(5) = 0.5d0*b(3)**2-geoel(7,ielem)
     b(6) = b(2)*b(3) - geoel(8,ielem)
!
     rho = 0.d0; rhou = 0.d0; rhov = 0.d0; rhoe = 0.d0;
!     
     do ideg = 1, ndegr  
      rho = rho + unkno(ideg, 1, ielem)*b(ideg)
      rhou= rhou+ unkno(ideg, 2, ielem)*b(ideg)
      rhov= rhov+ unkno(ideg, 3, ielem)*b(ideg)
      rhoe= rhoe+ unkno(ideg, 4, ielem)*b(ideg)
     enddo  
!
     pres = (gamma-1.d0)*(rhoe-0.5d0*(rhou*rhou+rhov*rhov)/rho)
     difun = (pres*gamma)*(1.d0/rho)**gamma - 1.d0
     errl2 = errl2 + difun*difun*djaco
! 
  enddo
 enddo !...(1)ie = 1,nelem
!
!
!    print*,'L2 error',0.5d0*log10(errl2),0.5d0*log10(100.0)
!

 end subroutine geterror_resi
!
!...Subroutine: Check the mesh quality (triangle)...
!
subroutine getoutput_celllag(inpoel,ipqua, coord,geoel,unkno)
use constant
implicit none
!...Real array
integer*4,dimension(1:nvtri,1:ntria)::inpoel
integer*4,dimension(1:nvqua,1:nquad)::ipqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unkno
!...Integer array
real*8,dimension(1:ndimn,1:npoin), intent(in)::coord
real*8,dimension(1:ngeel,1:nsize), intent(in)::geoel
!
!...Local
integer:: ielem, ip, ie
real*8 :: eps
real*8 :: rho,uadv, vadv, eadv
!
!...Real array
real*8, allocatable::pres(:)
!
  data   eps / 1.0d-06 /
!
!...Step 1: Find pressure, Mach...
!
  allocate(pres(1:ncell))
!
  do ielem = 1,ncell
!
     if(ndens.eq.1)then
      rho  = 1.d0/unkno(1, 1, ielem)
     elseif(ndens.eq.2)then
      rho  = unkno(1, 1, ielem)
     elseif(ndens.eq.3)then
      rho  = unkno(1, 1, ielem)
     endif
      uadv = unkno(1, 2, ielem)
      vadv = unkno(1, 3, ielem)
      eadv = unkno(1, 4, ielem)
!
      pres(ielem) = (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
!

  enddo 
!
!...Step 2: Output tecplot format...
!
  open(11,file='output_cell.dat')
  write(11,*)'VARIABLES = "X","Y","Density","U","V","P",'
  write(11,*)'ZONE'
  write(11,*)'     T="2D"'
  write(11,*)'     DATAPACKING=BLOCK'
  !write(11,*)'     ZONETYPE=FEtriangle'
  write(11,*)'     ZONETYPE=FEQuadrilateral'
  write(11,*)'     N=',npoin
  write(11,*)'     E=',ncell
  write(11,*)'     VARLOCATION=([3-6]=CELLCENTERED)'
!
!...X, Y coordinates...
  do ip =1,npoin
     write(11,*)coord(1,ip)
  enddo
  do ip=1,npoin
    write(11,*)coord(2,ip)
  enddo
!
!...Density
  if(ndens.eq.1)then
  do ie =1,ncell
    write(11,*)1.d0/unkno(1, 1, ie)
  enddo
  elseif(ndens.eq.2)then
  do ie =1,ncell
    write(11,*)unkno(1, 1, ie)
  enddo
  elseif(ndens.eq.3)then
  do ie =1,ncell
    write(11,*)unkno(1, 1, ie)
  enddo
  endif
!
!...U
  do ie =1,ncell
    write(11,*)unkno(1, 2, ie) 
  enddo
!
!...V	
  do ie =1,ncell
    write(11,*)unkno(1, 3, ie) 
  enddo
!
!...Pressure
  do ie =1,ncell
    write(11,*)pres(ie)!pres(ie)
  enddo
!
    write(11,*)
!
  do ie=1,ntria
    write(11,*)inpoel(1:3,ie),inpoel(3,ie) 
  enddo
!
  do ie=1,nquad
    write(11,*)ipqua(1:4,ie) 
  enddo
 close(11)
!
!...
  deallocate(pres)
!
end subroutine getoutput_celllag
!
!...Lagrangian hybrid meshes with .vtu format
!
subroutine getoutput_cellpara_hybrid(inpoel,coord,geoel,unkno,ujacb,iptri, ipqua)
use constant
implicit none
integer::cntr,cntr2,dump,ie,ip
integer::mtype(2)
integer,intent(in)::inpoel(nvtri,ntria)
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
real*8, intent(in)::coord(ndimn,npoin)
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ncell),             intent(in)::ujacb
!
real*8:: rho, uadv, vadv, eadv

open(22,file='output_cell_para.vtu',status='unknown')
write(22,"(T1,A)")'<?xml version="1.0"?>'
write(22,"(T1,A)")'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
write(22,"(T1,A)")'<UnstructuredGrid>'
write(22,*)'<Piece NumberOfPoints="',npoin,'" NumberOfCells="',ncell,'" >'
!
!...Set mesh type
!
if(ncurv==0)then
mtype(1) = 5
mtype(2) = 9
elseif(ncurv==1)then
mtype(1) = 22
mtype(2) = 23
endif
!
!...Output cell-based data...
!
write(22,"(T1,A)")'<CellData Scalars="scalars">'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Density_cell" NumberOfComponents="1" format="ascii">'
if(ndens.eq.1)then
 do ie = 1,ncell
   write(22,'(3e20.8)') 1.d0/unkno(1,1,ie)
 end do
else
 do ie = 1,ncell
   write(22,'(3e20.8)') unkno(1,1,ie)
 end do
endif
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Velocity_cell" NumberOfComponents="3" format="ascii">'
do ie = 1,ncell
write(22,'(3e20.8)') unkno(1,2,ie),unkno(1,3,ie),0.0
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Pressure_cell" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
!
if(ndens.eq.1)then
rho  = 1.d0/unkno(1, 1, ie)
else
rho  = unkno(1, 1, ie)
endif
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Internal energy_cell" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
!
if(ndens.eq.1)then
rho  = 1.d0/unkno(1, 1, ie)
else
rho  = unkno(1, 1, ie)
endif
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
write(22,"(T1,A)")'</DataArray>'

!...Limiter marker
write(22,"(T1,A)")'<DataArray type="Float32" Name="Limiter_inva" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
!
write(22,'(3e20.8)') geoel(10, ie)
end do
write(22,"(T1,A)")'</DataArray>'

!...Minimum Jacobian
write(22,"(T1,A)")'<DataArray type="Float32" Name="Jacobian_min" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
write(22,'(3e20.8)') ujacb(ie)
end do
write(22,"(T1,A)")'</DataArray>'


write(22,"(T1,A)")'</CellData>'
!
!...Output nodal coordinates...
!
write(22,"(T1,A)")'<Points>'
write(22,"(T1,A)")'<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
do ip = 1,npoin
write(22,'(3e14.6)')coord(1,ip),coord(2,ip),0.d0
end do      !i
write(22,"(T1,A)")'</DataArray>'
write(22,"(T1,A)")'</Points>'
!
!print*,'ip',coord(1:2,3630)
!
!...Output cell connectivity...
!
write(22,"(T1,A)")'<Cells>'
write(22,"(T1,A)")'<DataArray type="Int64" Name="connectivity" format="ascii">'
!
do ie=1,ntria
write(22,'(T1,6I6)')iptri(1:nvtri, ie) - 1!inpoel(1,ie)-1,inpoel(2,ie)-1,inpoel(3,ie)-1
cntr = cntr+1
end do
!
do ie=1,nquad
write(22,'(T1,9I6)')ipqua(1:nvqua, ie) - 1!inpoel(1,ie)-1,inpoel(2,ie)-1,inpoel(3,ie)-1
!print*,'ie',ie
cntr = cntr+1
end do
write(22,"(T1,A)")'</DataArray>'
!
!...Output offsets...
!
write(22,"(T1,A)")'<DataArray type="Int64" Name="offsets" format="ascii">'
cntr2 = 0
cntr = 0

do ie=1,ntria
cntr = cntr + nvtri
write(22,'(T1,1I8)')cntr
end do
!cntr = cntr2
!
do ie=1,nquad
cntr = cntr + nvqua
write(22,'(T1,1I8)')cntr
end do
cntr = cntr2
write(22,"(T1,A)")'</DataArray>'
!
!...Output the mesh type...
!
write(22,"(T1,A)")'<DataArray type="UInt8" Name="types" format="ascii">'
!-------- tria
do ie=1,ntria
write(22,'(T1,1I8)') mtype(1)
end do
do ie=1,nquad
write(22,'(T1,1I8)') mtype(2)
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'</Cells>'

write(22,"(T1,A)")'</Piece>'
write(22,"(T1,A)")'</UnstructuredGrid>'
write(22,"(T1,A)")'</VTKFile>'

close (22)

write(*,'(1A17,1A35)') 'Data written to:'!,trim(adjustl(filename2))
end subroutine getoutput_cellpara_hybrid
!
!...Lagrangian hybrid meshes
!
subroutine getfile_animation_hybrid(inpoel,coord,geoel,unkno,iptri, ipqua)
use constant
implicit none
integer::cntr,cntr2,dump,ie,ip
integer::mtype(2)
character(len=80):: fname
integer,intent(in)::inpoel(nvtri,ntria)
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
real*8, intent(in)::coord(ndimn,npoin)
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
!
real*8:: rho, uadv, vadv, eadv
!
!
nofile = nofile +1
!
!print*,nofile
!
! write local restart file
!
call rtname(nofile,fname)
!
!print*,fname
!
!write(filename2,'(1I1,1A1,1A50)')dgp,'.',filename
!filename3 = trim(adjustl(filename2)) // '.'//'vtu'
!open(22,file='./output/'//trim(adjustl(filename3)),status='unknown')
open(22,file=fname,status='unknown')
write(22,"(T1,A)")'<?xml version="1.0"?>'
write(22,"(T1,A)")'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
write(22,"(T1,A)")'<UnstructuredGrid>'
write(22,*)'<Piece NumberOfPoints="',npoin,'" NumberOfCells="',ncell,'" >'
!
!print*,'bad',fname
!
!...Output nodal data...
!
if(ncurv==0)then
mtype(1) = 5
mtype(2) = 9
elseif(ncurv==1)then
mtype(1) = 22
mtype(2) = 23
endif
!
!...Output cell-based data...
!
write(22,"(T1,A)")'<CellData Scalars="scalars">'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Density_cell" NumberOfComponents="1" format="ascii">'
if(ndens.eq.1)then
do ie = 1,ncell
write(22,'(3e20.8)') 1.d0/unkno(1,1,ie)
end do
else
do ie = 1,ncell
write(22,'(3e20.8)') unkno(1,1,ie)
end do
endif
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Velocity_cell" NumberOfComponents="3" format="ascii">'
do ie = 1,ncell
write(22,'(3e20.8)') unkno(1,2,ie),unkno(1,3,ie),0.0
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Pressure_cell" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
!
if(ndens.eq.1)then
rho  = 1.d0/unkno(1, 1, ie)
else
rho  = unkno(1, 1, ie)
endif
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Internal energy_cell" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
!
if(ndens.eq.1)then
rho  = 1.d0/unkno(1, 1, ie)
else
rho  = unkno(1, 1, ie)
endif
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
write(22,"(T1,A)")'</DataArray>'


write(22,"(T1,A)")'<DataArray type="Float32" Name="Limiter_inva" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
!
write(22,'(3e20.8)') geoel(10, ie)
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'</CellData>'
!
!...Output nodal coordinates...
!
write(22,"(T1,A)")'<Points>'
write(22,"(T1,A)")'<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
do ip = 1,npoin
write(22,'(3e14.6)')coord(1,ip),coord(2,ip),0.d0
end do      !i
write(22,"(T1,A)")'</DataArray>'
write(22,"(T1,A)")'</Points>'
!
!...Output cell connectivity...
!
write(22,"(T1,A)")'<Cells>'
write(22,"(T1,A)")'<DataArray type="Int64" Name="connectivity" format="ascii">'
!
do ie=1,ntria
write(22,'(T1,6I6)')iptri(1:nvtri, ie) - 1!inpoel(1,ie)-1,inpoel(2,ie)-1,inpoel(3,ie)-1
cntr = cntr+1
end do
!
do ie=1,nquad
write(22,'(T1,9I6)')ipqua(1:nvqua, ie) - 1!inpoel(1,ie)-1,inpoel(2,ie)-1,inpoel(3,ie)-1
!print*,'ie',ie
cntr = cntr+1
end do
write(22,"(T1,A)")'</DataArray>'
!
!...Output offsets...
!
write(22,"(T1,A)")'<DataArray type="Int64" Name="offsets" format="ascii">'
cntr2 = 0
cntr = 0

do ie=1,ntria
cntr = cntr + nvtri
write(22,'(T1,1I8)')cntr
end do
!cntr = cntr2
!
do ie=1,nquad
cntr = cntr + nvqua
write(22,'(T1,1I8)')cntr
end do
cntr = cntr2
write(22,"(T1,A)")'</DataArray>'
!
!...Output the mesh type...
!
write(22,"(T1,A)")'<DataArray type="UInt8" Name="types" format="ascii">'
!-------- tria
do ie=1,ntria
write(22,'(T1,1I8)') mtype(1)
end do
do ie=1,nquad
write(22,'(T1,1I8)') mtype(2)
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'</Cells>'

write(22,"(T1,A)")'</Piece>'
write(22,"(T1,A)")'</UnstructuredGrid>'
write(22,"(T1,A)")'</VTKFile>'

close (22)

!write(*,'(1A17,1A35)') 'Data written to:'!,trim(adjustl(filename2))
end subroutine getfile_animation_hybrid
!
!...Lagrangian hybrid meshes for vtk legacy
!
subroutine getoutput_cellvtk(inpoel,coord,geoel,unkno,iptri, ipqua)
use constant
implicit none
integer::cntr,cntr2,dump,ie,ip
integer::mtype(2)
integer,intent(in)::inpoel(nvtri,ntria)
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
real*8, intent(in)::coord(ndimn,npoin)
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
!
real*8:: rho, uadv, vadv, eadv

open(22,file='output.vtk',status='unknown')
write(22,"(T1,A)")'# vtk DataFile Version 2.0'
write(22,"(T1,A)")'Unstructured Grid Example'
write(22,"(T1,A)")'ASCII'
write(22,*)
write(22,"(T1,A)")'DATASET UNSTRUCTURED_GRID'
write(22,*)'Points',npoin,' float '
!
!...Cell type...
!
if(ncurv==0)then
mtype(1) = 5
mtype(2) = 9
elseif(ncurv==1)then
mtype(1) = 22
mtype(2) = 23
endif
!
!...Output nodal coordinates...
!
do ip = 1,npoin
write(22,*)coord(1,ip),coord(2,ip),0.d0
end do

write(22,*)
write(22,*)'CELLS',ncell,(nvtri+1)*ntria + (nvqua+1)*nquad
!
do ie=1,ntria
write(22,'(T1,7I10)')nvtri, iptri(1:nvtri, ie) - 1
end do
!
do ie=1,nquad
write(22,'(T1,10I10)')nvqua, ipqua(1:nvqua, ie) - 1
end do
!

write(22,*)
write(22,*)'CELL_TYPES',ncell
!
do ie=1,ntria
write(22,*)mtype(1)
end do
!
do ie=1,nquad
write(22,*)mtype(2)
end do
!
!...Scalar
!
write(22,*)
write(22,*)
write(22,*)'CELL_DATA',ncell
write(22,*)'SCALARS density float'
write(22,*)'LOOKUP_TABLE default'
!-xxx-Density
if(ndens.eq.1)then
do ie=1,ncell
write(22,*)1.d0/unkno(1,1,ie)
end do
else
do ie=1,ncell
write(22,*)unkno(1,1,ie)
end do
endif

write(22,*)
write(22,*)'SCALARS internalenergy float'
write(22,*)'LOOKUP_TABLE default'
!-xxx-IE
if(ndens.eq.1)then
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
else
do ie=1,ncell
rho =  unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
endif

write(22,*)
write(22,*)'SCALARS pressure float'
write(22,*)'LOOKUP_TABLE default'
!-xxx-pres
if(ndens.eq.1)then
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
else
do ie=1,ncell
rho =  unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
endif

write(22,*)
write(22,*)'SCALARS Badcell float'
write(22,*)'LOOKUP_TABLE default'
!-xxx-shock dtector
do ie=1,ncell
write(22,'(1e20.8)') geoel(10, ie)
end do
!
!...Vectors
!
write(22,*)
write(22,*)'VECTORS vectors float'
!-xxx-Velocity
do ie=1,ncell
write(22,'(3e32.16)')unkno(1,2:3,ie),0.d0
end do

write(*,'(1A17,1A35)') 'Data written !'
end subroutine getoutput_cellvtk
!
!...Lagrangian hybrid meshes for vtk legacy for high-order cells
!
subroutine getoutput_cellvtk_ho(inpoel,coord,geoel,unkno,iptri, ipqua)
use constant
implicit none
integer::cntr,cntr2,dump,ie,ip
integer::mtype(2)
integer,intent(in)::inpoel(nvtri,ntria)
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
integer, dimension(nvqua)::ipq,idmpp
real*8, intent(in)::coord(ndimn,npoin)
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
real*8, dimension(1:nvqua)::shpq,dsprq, dspsq
real*8, dimension(1:ndimn, 1:nvqua)::xpq
!
integer, allocatable::ipqdh(:,:),iptah(:,:)
real*8,  allocatable::coorh(:,:)
real*8,  allocatable::rqh(:), sqh(:)
real*8:: rho, uadv, vadv, eadv
real*8::xh, yh
!
integer::npinh,nvqdh,nvtah,ielem,ishp
integer::iph,nivqh, ipoin
!
!...Cell type...
!
if(ncurv==0)then
mtype(1) = 5
mtype(2) = 9
elseif(ncurv==1)then
mtype(1) = 22
mtype(2) = 23
!
npinh = npoin
elseif(ncurv==2)then !...VTK_LAGRANGE_QUADRILATERAL = 70,
mtype(1) = 69
mtype(2) = 70

!...Mapping array
idmpp(1:4) =  (/1,2,3,4/)
idmpp(5:8) =  (/5,9,6,10/)
idmpp(9:12) = (/11,7,12,8/)
!
npinh = npoin + ntria*1 + nquad*4
nvtah = nvtri
nvqdh = 16
nivqh = 4
!
allocate(rqh(nivqh), sqh(nivqh))
!
rqh(1) = -1.d0/3.d0; rqh(2) =  1.d0/3.d0; rqh(3) =  1.d0/3.d0; rqh(4) = -1.d0/3.d0
sqh(1) = -1.d0/3.d0; sqh(2) = -1.d0/3.d0; sqh(3) =  1.d0/3.d0; sqh(4) =  1.d0/3.d0

endif
!
allocate(coorh(ndimn, npinh))
allocate(iptah(nvtah, ntria), ipqdh(nvqdh, nquad))

!...New array for quad
ipqdh(1:nvqua, :)  = ipqua(idmpp(1:nvqua), :)

do ie = 1, nquad
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

!...Points consitituting one element...
xpq(1, 1:nvqua) = coord(1,ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2,ipq(1:nvqua))

!...Get the points for FEM
call getcoord_fe(ncurv, nvfac, nvqua, xpq)
!
do iph = 1, nivqh

!...Point No.
ipoin = npoin + (ie-1)*4 + iph
ipqdh(nvqua+iph, ie) = ipoin

!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, rqh(iph), sqh(iph))

xh = 0.d0
yh = 0.d0
do ishp = 1, nvqua
xh = xh + shpq(ishp)*xpq(1,ishp)
yh = yh + shpq(ishp)*xpq(2,ishp)
enddo
!
coorh(1,ipoin) = xh
coorh(2,ipoin) = yh
enddo
!
coorh(1,ipq(1:nvqua)) = xpq(1, 1:nvqua)
coorh(2,ipq(1:nvqua)) = xpq(2, 1:nvqua)

enddo

open(22,file='output_ho.vtk',status='unknown')
write(22,"(T1,A)")'# vtk DataFile Version 4.2'
write(22,"(T1,A)")'vtk output'
write(22,"(T1,A)")'ASCII'
write(22,"(T1,A)")'DATASET UNSTRUCTURED_GRID'
write(22,*)'POINTS ',npinh,' double '
!
!...Output nodal coordinates...
!
do ip = 1,npinh
write(22,*)coorh(1,ip),coorh(2,ip),0.d0
end do

write(22,*)
write(22,*)'CELLS',ncell,(nvtah+1)*ntria + (nvqdh+1)*nquad
!
do ie=1,ntria
write(22,'(T1,7I10)')nvtah, iptah(1:nvtah, ie) - 1
end do
!
do ie=1,nquad
write(22,'(T1,25I10)')nvqdh, ipqdh(1:nvqdh, ie) - 1
end do
!

write(22,*)
write(22,*)'CELL_TYPES',ncell
!
do ie=1,ntria
write(22,*)mtype(1)
end do
!
do ie=1,nquad
write(22,*)mtype(2)
end do
!
!...Scalar
!
write(22,*)'CELL_DATA',ncell
write(22,*)'FIELD FieldData 4'

write(22,*)'Density 1',ncell,' double '
!...Desnity
if(ndens.eq.1)then
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') rho
end do
else
do ie=1,ncell
rho =  unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') rho
end do
endif

!...Pres
write(22,*)'Pressure 1',ncell,' double '
!-xxx-IE
if(ndens.eq.1)then
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
else
do ie=1,ncell
rho =  unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
endif

!...Troubled-cell inicator
write(22,*)'Troubled-cell 1',ncell,' double '
do ie=1,ncell
write(22,'(3e20.8)') geoel(10 ,ie)
end do

!...Velocity
write(22,*)'Velocity 2',ncell,' double '
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)')uadv, vadv
end do

write(*,'(1A17,1A35)') 'Data written for a file!'
end subroutine getoutput_cellvtk_ho
!
!...Lagrangian hybrid meshes for vtk legacy for high-order cells
!
subroutine getfile_animation_cellvtk_ho(coord,geoel,unkno,iptri, ipqua)
use constant
implicit none
integer::cntr,cntr2,dump,ie,ip
integer::mtype(2)
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
integer, dimension(nvqua)::ipq,idmpp
real*8, intent(in)::coord(ndimn,npoin)
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
real*8, dimension(1:nvqua)::shpq,dsprq, dspsq
real*8, dimension(1:ndimn, 1:nvqua)::xpq
!
integer, allocatable::ipqdh(:,:),iptah(:,:)
real*8,  allocatable::coorh(:,:)
real*8,  allocatable::rqh(:), sqh(:)
real*8:: rho, uadv, vadv, eadv
real*8::xh, yh
!
integer::npinh,nvqdh,nvtah,ielem,ishp
integer::iph,nivqh, ipoin
character(len=80):: fname
!
nofile = nofile +1
!
call rtnamevtk(nofile,fname)
!
!...Cell type...
!
if(ncurv==0)then
mtype(1) = 5
mtype(2) = 9
elseif(ncurv==1)then
mtype(1) = 22
mtype(2) = 23
!
npinh = npoin
elseif(ncurv==2)then !...VTK_LAGRANGE_QUADRILATERAL = 70,
mtype(1) = 69
mtype(2) = 70

!...Mapping array
idmpp(1:4) =  (/1,2,3,4/)
idmpp(5:8) =  (/5,9,6,10/)
idmpp(9:12) = (/11,7,12,8/)
!
npinh = npoin + ntria*1 + nquad*4
nvtah = nvtri
nvqdh = 16
nivqh = 4
!
allocate(rqh(nivqh), sqh(nivqh))
!
rqh(1) = -1.d0/3.d0; rqh(2) =  1.d0/3.d0; rqh(3) =  1.d0/3.d0; rqh(4) = -1.d0/3.d0
sqh(1) = -1.d0/3.d0; sqh(2) = -1.d0/3.d0; sqh(3) =  1.d0/3.d0; sqh(4) =  1.d0/3.d0

endif
!
allocate(coorh(ndimn, npinh))
allocate(iptah(nvtah, ntria), ipqdh(nvqdh, nquad))

!...New array for quad
ipqdh(1:nvqua, :)  = ipqua(idmpp(1:nvqua), :)

do ie = 1, nquad
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)

!...Points consitituting one element...
xpq(1, 1:nvqua) = coord(1,ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2,ipq(1:nvqua))

!...Get the points for FEM
call getcoord_fe(ncurv, nvfac, nvqua, xpq)
!
do iph = 1, nivqh

!...Point No.
ipoin = npoin + (ie-1)*4 + iph
ipqdh(nvqua+iph, ie) = ipoin

!...  shape function & its derivatives w.r.t. reference coordinates
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, rqh(iph), sqh(iph))

xh = 0.d0
yh = 0.d0
do ishp = 1, nvqua
xh = xh + shpq(ishp)*xpq(1,ishp)
yh = yh + shpq(ishp)*xpq(2,ishp)
enddo
!
coorh(1,ipoin) = xh
coorh(2,ipoin) = yh
enddo
!
coorh(1,ipq(1:nvqua)) = xpq(1, 1:nvqua)
coorh(2,ipq(1:nvqua)) = xpq(2, 1:nvqua)

enddo

!
open(22,file=fname,status='unknown')
write(22,"(T1,A)")'# vtk DataFile Version 4.2'
write(22,"(T1,A)")'vtk output'
write(22,"(T1,A)")'ASCII'
write(22,*)
write(22,"(T1,A)")'DATASET UNSTRUCTURED_GRID'
write(22,*)'POINTS',npinh,'double'
!
!...Output nodal coordinates...
!
do ip = 1,npinh
write(22,*)coorh(1,ip),coorh(2,ip),0.d0
end do

write(22,*)
write(22,*)'CELLS',ncell,(nvtah+1)*ntria + (nvqdh+1)*nquad
!
do ie=1,ntria
write(22,'(T1,7I10)')nvtah, iptah(1:nvtah, ie) - 1
end do
!
do ie=1,nquad
write(22,'(T1,25I10)')nvqdh, ipqdh(1:nvqdh, ie) - 1
end do
!

write(22,*)
write(22,*)'CELL_TYPES',ncell
!
do ie=1,ntria
write(22,*)mtype(1)
end do
!
do ie=1,nquad
write(22,*)mtype(2)
end do
!
!...Scalar
!
write(22,*)'CELL_DATA',ncell
write(22,*)'FIELD FieldData 4'

write(22,*)'Density 1',ncell,' double '
!...Desnity
if(ndens.eq.1)then
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') rho
end do
else
do ie=1,ncell
rho =  unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') rho
end do
endif

!...Pres
write(22,*)'Pressure 1',ncell,' double '
!-xxx-IE
if(ndens.eq.1)then
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
else
do ie=1,ncell
rho =  unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
endif

!...Troubled-cell inicator
write(22,*)'Troubled-cell 1',ncell,' double '
do ie=1,ncell
write(22,'(3e20.8)') geoel(10 ,ie)
end do

!...Velocity
write(22,*)'Velocity 2',ncell,' double '
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)')uadv, vadv
end do

close(22)

end subroutine getfile_animation_cellvtk_ho

!
!...Lagrangian hybrid meshes for MEM
!
subroutine getoutput_particlevtk(coord,geoel,unkno,iptri, ipqua)
use constant
implicit none
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
real*8, intent(in)::coord(ndimn,npoin)
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
!
real*8::shpq(nvqua)
real*8::weighq(ngausdq_pst), posiq(2,ngausdq_pst)
real*8::xgq(1:2,ngausdq_pst,nquad)
real*8::unkng(1:nq,ngausdq_pst,nquad)
real*8::b(ndegr)
!
integer::ie,ip, ig, ielem, ishp, ideg
!
real*8:: rho, uadv, vadv, eadv,rhom
real*8:: r, s, rm, sm, rp, sp
real*8:: rc, sc, dr, ds
!
real*8:: c10
!
c10 = 1.d0
!
call ruqope(2, ngausdq_pst, posiq, weighq)
!
!...Part I: Get the particle positions and variables
!
xgq = 0.d0
unkng = 0.d0
!
do ie = 1, nquad
!
ielem = ie + ntria

!
dr = 1.d0
ds = 1.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Gauss loop
do ig = 1,ngausdq_pst!...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
!
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*rm*sm
shpq(2) = 0.25d0*rp*sm
shpq(3) = 0.25d0*rp*sp
shpq(4) = 0.25d0*rm*sp

!...Get the points
do ishp = 1, 4
xgq(1, ig, ie) = xgq(1, ig, ie) + shpq(ishp)*coord(1,ipqua(ishp, ie))
xgq(2, ig, ie) = xgq(2, ig, ie) + shpq(ishp)*coord(2,ipqua(ishp, ie))
enddo

!...Get the variable values

!...Basis function for solutions...
b(1) = 1.d0

if(npoly.ge.1)then
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
endif
!
rhom = 0.d0
uadv = 0.d0
vadv = 0.d0
eadv = 0.d0
!
do ideg = 1, 1!ndegr
rhom = rhom + unkno(ideg, 1, ielem)*b(ideg)
uadv = uadv + unkno(ideg, 2, ielem)*b(ideg)
vadv = vadv + unkno(ideg, 3, ielem)*b(ideg)
eadv = eadv + unkno(ideg, 4, ielem)*b(ideg)
enddo
!
unkng(1, ig, ielem) = 1.d0/rhom
unkng(2, ig, ielem) = uadv
unkng(3, ig, ielem) = vadv
unkng(4, ig, ielem) = (gamlg-1.d0)*1.d0/rhom*(eadv-0.5d0*(uadv**2 + vadv**2))

enddo
!
enddo
!
!...Part II: output for vtk
!
open(22,file='output-particle.vtk',status='unknown')
write(22,"(T1,A)")'# vtk DataFile Version 2.0'
write(22,"(T1,A)")'Unstructured Grid Example'
write(22,"(T1,A)")'ASCII'
write(22,*)
write(22,"(T1,A)")'DATASET UNSTRUCTURED_GRID'
write(22,*)'Points',nquad*ngausdq_pst,' float '
!
do ie=1,nquad
do ig = 1,ngausdq_pst
write(22,*)xgq(1:2,ig,ie),0.d0
enddo
end do
!
write(22,*)
write(22,*)'CELL_TYPES',ncell*ngausdq_pst
!
do ie=1,nquad
do ig = 1,ngausdq_pst
write(22,*)1
enddo
end do
!
!...Scalar
!
write(22,*)
write(22,*)
write(22,*)'POINT_DATA',ncell*ngausdq_pst
write(22,*)'SCALARS density float'
write(22,*)'LOOKUP_TABLE default'

!-xxx-Density
if(ndens.eq.1)then
do ie=1,nquad
ielem = ie +ntria
do ig = 1,ngausdq_pst
write(22,*)unkng(1, ig, ielem)
enddo
enddo
endif
!
write(22,*)
write(22,*)'SCALARS pressure float'
write(22,*)'LOOKUP_TABLE default'

!-xxx-Density
if(ndens.eq.1)then
do ie=1,nquad
ielem = ie +ntria
do ig = 1,ngausdq_pst
write(22,*)unkng(4, ig, ielem)
enddo
enddo
endif

write(22,*)
write(22,*)'SCALARS velocity float'
write(22,*)'LOOKUP_TABLE default'

!-xxx-Density
if(ndens.eq.1)then
do ie=1,nquad
ielem = ie +ntria
do ig = 1,ngausdq_pst
write(22,*)sqrt(unkng(2, ig, ielem)**2 + unkng(3, ig, ielem)**2)
enddo
enddo
endif


!...Vectors
write(22,*)
write(22,*)'VECTORS vectors float'
!-xxx-Velocity
do ie=1,nquad
ielem = ie +ntria
do ig = 1,ngausdq_pst
write(22,'(3e32.16)')unkng(2:3,ig,ielem),0.d0
end do
enddo

!
write(*,'(1A17,1A35)') 'Particle data written for MEM!'
end subroutine getoutput_particlevtk
!
!...Output files (.vtk) for animation on Lagrangian hybrid meshes
!
subroutine getfile_animation_hybridvtk(inpoel,coord,geoel,unkno,iptri, ipqua)
use constant
implicit none
integer::cntr,cntr2,dump,ie,ip
integer::mtype(2)
character(len=80):: fname
integer,intent(in)::inpoel(nvtri,ntria)
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
real*8, intent(in)::coord(ndimn,npoin)
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
!
real*8:: rho, uadv, vadv, eadv
!
nofile = nofile +1
!
call rtnamevtk(nofile,fname)
!
open(22,file=fname,status='unknown')
write(22,"(T1,A)")'# vtk DataFile Version 2.0'
write(22,"(T1,A)")'Unstructured Grid Example'
write(22,"(T1,A)")'ASCII'
write(22,*)
write(22,"(T1,A)")'DATASET UNSTRUCTURED_GRID'
write(22,*)'Points',npoin,' float '
!
!...Cell type...
!
if(ncurv==0)then
mtype(1) = 5
mtype(2) = 9
elseif(ncurv==1)then
mtype(1) = 22
mtype(2) = 23
endif
!
!...Output nodal coordinates...
!
do ip = 1,npoin
write(22,'(3e20.8)')coord(1,ip),coord(2,ip),0.d0
end do

write(22,*)
write(22,*)'CELLS',ncell,(nvtri+1)*ntria + (nvqua+1)*nquad
!
do ie=1,ntria
write(22,'(T1,7I10)')nvtri, iptri(1:nvtri, ie) - 1
end do
!
do ie=1,nquad
write(22,'(T1,10I10)')nvqua, ipqua(1:nvqua, ie) - 1
end do
!

write(22,*)
write(22,*)'CELL_TYPES',ncell
!
do ie=1,ntria
write(22,*)mtype(1)
end do
!
do ie=1,nquad
write(22,*)mtype(2)
end do
!
!...Scalar
!
write(22,*)
write(22,*)
write(22,*)'CELL_DATA',ncell
write(22,*)'SCALARS density float'
write(22,*)'LOOKUP_TABLE default'
!-xxx-Density
if(ndens.eq.1)then
do ie=1,ncell
write(22,*)1.d0/unkno(1,1,ie)
end do
else
do ie=1,ncell
write(22,*)unkno(1,1,ie)
end do
endif

write(22,*)
write(22,*)'SCALARS internalenergy float'
write(22,*)'LOOKUP_TABLE default'
!-xxx-IE
if(ndens.eq.1)then
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
else
do ie=1,ncell
rho =  unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
endif

write(22,*)
write(22,*)'SCALARS pressure float'
write(22,*)'LOOKUP_TABLE default'
!-xxx-pres
if(ndens.eq.1)then
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
else
do ie=1,ncell
rho =  unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
endif

write(22,*)
write(22,*)'SCALARS shockmarker float'
write(22,*)'LOOKUP_TABLE default'
!-xxx-pres
do ie=1,ncell
rho = 1.d0/unkno(1,1,ie)
uadv = unkno(1, 2, ie)
vadv = unkno(1, 3, ie)
eadv = unkno(1, 4, ie)
write(22,'(3e20.8)') geoel(10,ie)
end do
!
!...Vectors
!
write(22,*)
write(22,*)'VECTORS vectors float'
!-xxx-Velocity
do ie=1,ncell
write(22,'(3e32.16)')unkno(1,2:3,ie),0.d0
end do

close (22)

!write(*,'(1A17,1A35)') 'Data written to:'!,trim(adjustl(filename2))
end subroutine getfile_animation_hybridvtk
!
!...Output animation files for MEM
!
subroutine getfile_animation_particlevtk(coord,geoel,unkno,iptri, ipqua)
use constant
implicit none
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
real*8, intent(in)::coord(ndimn,npoin)
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
!
character(len=80):: fname
!
real*8::shpq(nvqua)
real*8::weighq(ngausdq_pst), posiq(2,ngausdq_pst)
real*8::xgq(1:2,ngausdq_pst,nquad)
real*8::unkng(1:nq,ngausdq_pst,nquad)
real*8::b(ndegr)
!
integer::ie,ip, ig, ielem, ishp, ideg
!
real*8:: rho, uadv, vadv, eadv,rhom
real*8:: r, s, rm, sm, rp, sp
real*8:: rc, sc, dr, ds
!
real*8:: c10
!
c10 = 1.d0
!
call ruqope(2, ngausdq_pst, posiq, weighq)
!
!...Part I: Get the particle positions and variables
!
xgq = 0.d0
unkng = 0.d0
!
do ie = 1, nquad
!
ielem = ie + ntria

!
dr = 1.d0
ds = 1.d0
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...Gauss loop
do ig = 1,ngausdq_pst!...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
!
rp = c10 + r
rm = c10 - r
sp = c10 + s
sm = c10 - s

!...  shape function & its derivatives w.r.t. reference coordinates
shpq(1) = 0.25d0*rm*sm
shpq(2) = 0.25d0*rp*sm
shpq(3) = 0.25d0*rp*sp
shpq(4) = 0.25d0*rm*sp

!...Get the points
do ishp = 1, 4
xgq(1, ig, ie) = xgq(1, ig, ie) + shpq(ishp)*coord(1,ipqua(ishp, ie))
xgq(2, ig, ie) = xgq(2, ig, ie) + shpq(ishp)*coord(2,ipqua(ishp, ie))
enddo

!...Get the variable values

!...Basis function for solutions...
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds

!
rhom = 0.d0
uadv = 0.d0
vadv = 0.d0
eadv = 0.d0
!
do ideg = 1, 3! ndegr
rhom = rhom + unkno(ideg, 1, ielem)*b(ideg)
uadv = uadv + unkno(ideg, 2, ielem)*b(ideg)
vadv = vadv + unkno(ideg, 3, ielem)*b(ideg)
eadv = eadv + unkno(ideg, 4, ielem)*b(ideg)
enddo
!
unkng(1, ig, ielem) = 1.d0/rhom
unkng(2, ig, ielem) = uadv
unkng(3, ig, ielem) = vadv
unkng(4, ig, ielem) = (gamlg-1.d0)*1.d0/rhom*(eadv-0.5d0*(uadv**2 + vadv**2))

enddo
!
enddo
!
!...Part II: output for vtk
!
nofile = nofile +1
!
call rtnamevtk(nofile,fname)
!
open(22,file=fname,status='unknown')
write(22,"(T1,A)")'# vtk DataFile Version 2.0'
write(22,"(T1,A)")'Unstructured Grid Example'
write(22,"(T1,A)")'ASCII'
write(22,*)
write(22,"(T1,A)")'DATASET UNSTRUCTURED_GRID'
write(22,*)'Points',nquad*ngausdq_pst,' float '
!
do ie=1,nquad
do ig = 1,ngausdq_pst
write(22,*)xgq(1:2,ig,ie),0.d0
enddo
end do
!
write(22,*)
write(22,*)'CELL_TYPES',ncell*ngausdq_pst
!
do ie=1,nquad
do ig = 1,ngausdq_pst
write(22,*)1
enddo
end do
!
!...Scalar
!
write(22,*)
write(22,*)
write(22,*)'POINT_DATA',ncell*ngausdq_pst
write(22,*)'SCALARS density float'
write(22,*)'LOOKUP_TABLE default'

!-xxx-Density
if(ndens.eq.1)then
do ie=1,nquad
ielem = ie +ntria
do ig = 1,ngausdq_pst
write(22,*)unkng(1, ig, ielem)
enddo
enddo
endif
!
write(22,*)
write(22,*)'SCALARS pressure float'
write(22,*)'LOOKUP_TABLE default'

!-xxx-Density
if(ndens.eq.1)then
do ie=1,nquad
ielem = ie +ntria
do ig = 1,ngausdq_pst
write(22,*)unkng(4, ig, ielem)
enddo
enddo
endif

write(22,*)
write(22,*)'SCALARS velocity float'
write(22,*)'LOOKUP_TABLE default'

!-xxx-Density
if(ndens.eq.1)then
do ie=1,nquad
ielem = ie +ntria
do ig = 1,ngausdq_pst
write(22,*)sqrt(unkng(2, ig, ielem)**2 + unkng(3, ig, ielem)**2)
enddo
enddo
endif


!...Vectors
write(22,*)
write(22,*)'VECTORS vectors float'
!-xxx-Velocity
do ie=1,nquad
ielem = ie +ntria
do ig = 1,ngausdq_pst
write(22,'(3e32.16)')unkng(2:3,ig,ielem),0.d0
end do
enddo
!
close (22)
!
end subroutine getfile_animation_particlevtk

!
!...Output files to plot curved mesh using GNUPLOT
!
subroutine getoutput_cellgeo(inpoel,coord,geoel,unkno,iptri, ipqua)
use constant
implicit none
integer::cntr,cntr2,dump,ie,ip,iv
integer::mtype(2)
integer,intent(in)::inpoel(nvtri,ntria)
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
real*8, intent(in)::coord(ndimn,npoin)
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
!
real*8:: rho, uadv, vadv, eadv

!
open(22,file='output.geo',status='unknown')
!-xxx-cell
do ie=1,nquad
write(22,*)'#cell'
do iv = 1, nvqua
write(22,'(2e32.16)')coord(1:2, ipqua(1, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(5, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(5, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(9, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(9, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(2, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(2, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(6, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(6, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(10, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(10, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(3, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(3, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(7, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(7, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(11, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(11, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(4, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(4, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(8, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(8, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(12, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(12, ie))
write(22,'(2e32.16)')coord(1:2, ipqua(1, ie))
enddo

write(22,*)
write(22,*)
end do


end subroutine getoutput_cellgeo

!
!...Output data for paraview for eulerian framework...
!
subroutine getoutput_cellpara_euler(inpoel,coord,unkno,iptri, ipqua)
use constant
implicit none
integer::cntr,cntr2,dump,ie,ip
integer::mtype(2)
integer,intent(in)::inpoel(nvtri,ntria)
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
real*8, intent(in)::coord(ndimn,npoin)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
!
real*8:: rho, uadv, vadv, eadv
!write(filename2,'(1I1,1A1,1A50)')dgp,'.',filename
!filename3 = trim(adjustl(filename2)) // '.'//'vtu'
!open(22,file='./output/'//trim(adjustl(filename3)),status='unknown')
open(22,file='output_cell_para.vtu',status='unknown')
write(22,"(T1,A)")'<?xml version="1.0"?>'
write(22,"(T1,A)")'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
write(22,"(T1,A)")'<UnstructuredGrid>'
write(22,*)'<Piece NumberOfPoints="',npoin,'" NumberOfCells="',ncell,'" >'
!
!...Output nodal data...
!
if(ncurv==0)then
mtype(1) = 5
mtype(2) = 9
elseif(ncurv==1)then
mtype(1) = 22
mtype(2) = 23
endif
!
!...Output cell-based data...
!
write(22,"(T1,A)")'<CellData Scalars="scalars">'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Density_cell" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
write(22,'(3e20.8)') unkno(1,1,ie)
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Velocity_cell" NumberOfComponents="3" format="ascii">'
do ie = 1,ncell
write(22,'(3e20.8)') unkno(1,2,ie)/unkno(1,1,ie),unkno(1,3,ie)/unkno(1,1,ie),0.0
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Pressure_cell" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
!
rho  = unkno(1, 1, ie)
uadv = unkno(1, 2, ie)/rho
vadv = unkno(1, 3, ie)/rho
eadv = unkno(1, 4, ie)/rho
write(22,'(3e20.8)') (gamma-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'<DataArray type="Float32" Name="Mach_cell" NumberOfComponents="1" format="ascii">'
do ie = 1,ncell
!
rho  = unkno(1, 1, ie)
uadv = unkno(1, 2, ie)/rho
vadv = unkno(1, 3, ie)/rho
eadv = unkno(1, 4, ie)/rho
write(22,'(3e20.8)') sqrt(uadv**2 + vadv**2)/sqrt((gamma-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))*gamma/rho)
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'</CellData>'
!
!...Output nodal coordinates...
!
write(22,"(T1,A)")'<Points>'
write(22,"(T1,A)")'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
do ip = 1,npoin
! if(coord(1, ip).lt.0.d0) print*,'bad point',ip
write(22,'(3e14.6)')coord(1,ip),coord(2,ip),0.d0
end do      !i
write(22,"(T1,A)")'</DataArray>'
write(22,"(T1,A)")'</Points>'
!
!...Output cell connectivity...
!
write(22,"(T1,A)")'<Cells>'
write(22,"(T1,A)")'<DataArray type="Int32" Name="connectivity" format="ascii">'
!
do ie=1,ntria
write(22,'(T1,6I6)')iptri(1:nvtri, ie) - 1!inpoel(1,ie)-1,inpoel(2,ie)-1,inpoel(3,ie)-1
cntr = cntr+1
end do
!
do ie=1,nquad
write(22,'(T1,9I6)')ipqua(1:nvqua, ie) - 1!inpoel(1,ie)-1,inpoel(2,ie)-1,inpoel(3,ie)-1
cntr = cntr+1
end do
write(22,"(T1,A)")'</DataArray>'
!
!...Output offsets...
!
write(22,"(T1,A)")'<DataArray type="Int32" Name="offsets" format="ascii">'
cntr2 = 0
cntr = 0

do ie=1,ntria
cntr = cntr + nvtri
write(22,'(T1,1I8)')cntr
end do
!cntr = cntr2
!
do ie=1,nquad
cntr = cntr + nvqua
write(22,'(T1,1I8)')cntr
end do
cntr = cntr2
write(22,"(T1,A)")'</DataArray>'
!
!...Output the mesh type...
!
write(22,"(T1,A)")'<DataArray type="UInt8" Name="types" format="ascii">'
!-------- tria
do ie=1,ntria
write(22,'(T1,1I8)') mtype(1)
end do
do ie=1,nquad
write(22,'(T1,1I8)') mtype(2)
end do
write(22,"(T1,A)")'</DataArray>'

write(22,"(T1,A)")'</Cells>'

write(22,"(T1,A)")'</Piece>'
write(22,"(T1,A)")'</UnstructuredGrid>'
write(22,"(T1,A)")'</VTKFile>'

close (22)

write(*,'(1A17,1A35)') 'Data written to:'!,trim(adjustl(filename2))
end subroutine getoutput_cellpara_euler
!
!...Get the nodal values...
!
subroutine getpostnd_lag(unkno,inpoel,coord,geoel,unknp)
use constant
implicit none
!...Input array
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nelem+nbfac),     intent(in)::geoel
!Loacal
real*8,dimension(1:nvtri)::xv, yv
real*8,dimension(1:ndegr,1:nvtri)::b
real*8,dimension(1:npoin, 1:nq), intent(out)::unknp
real*8,  dimension(1:npoin)::cnsup
integer, dimension(1:nvtri):: ip

integer::iv, ideg,ipoin, ie
real*8::rho, uadv, vadv, pres,eadv
real*8::rc,sc,dr,ds
!
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
dr = 0.5d0
ds = 0.5d0
!
rc = 1.d0/3.d0
sc = rc
!
unknp = 0.d0
cnsup = 0.d0
!
!...Find No. of cells surrunding one node...dgflo2d.domn
!
do 300 ie = 1,nelem !...(1)ie = 1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!
!...Accumualte to get the no of cell surrounding one vertex...
!
cnsup(ip(1:nvtri)) = cnsup(ip(1:nvtri)) + 1.d0
300 enddo
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
!calculate the point value
do ipoin=1,npoin
!
!
do ie=1,nelem
!
ip(1:nvtri) = inpoel(1:nvtri,ie)
!

!
do iv = 1, nvtri
!
if(ip(iv).eq.ipoin)then
!
do ideg = 1,ndegr
unknp(ipoin, 1:nq) = unknp(ipoin, 1:nq) + unkno(ideg, 1:nq, ie)*b(ideg,iv)
enddo
endif
enddo
!
enddo
!
enddo !ipoin=1,npoin
!
!...Get nodal value
!
do ipoin = 1, npoin
unknp(ipoin, 1:nq) = unknp(ipoin, 1:nq)/cnsup(ipoin)
!
rho = 1.d0/unknp(ipoin, 1)
uadv = unknp(ipoin, 2)
vadv = unknp(ipoin, 3)
eadv = unknp(ipoin, 4)
!
pres = (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
!
unknp(ipoin, 1) = 1.d0/unknp(ipoin, 1)
unknp(ipoin, 2) = unknp(ipoin, 2)
unknp(ipoin, 3) = unknp(ipoin, 3)
unknp(ipoin, 4) = pres
!
enddo
!
!output the valu
return
! close(11)
end subroutine getpostnd_lag
!
!...Find volume at final stage in geoel for curved cell...
!
subroutine getgeoel_lag_post(iptri, ipqua, geoel, coord)
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
real*8:: r, s, djaco, volel
real*8:: wi
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
!
enddo
!
geoel(3, ielem) = volel
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
!
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
!
enddo
!
geoel(3, ielem) = volel
!
!    print*,'old geoel',ielem, geoel(1:5, ielem)
!    print*,'new geoel',ielem, xcent/volel, ycent/volel, volel, maxval(xp(1, 1:3))-minval(xp(1, 1:3)), &
!                         maxval(xp(2, 1:3))- minval(xp(2, 1:3))
!
enddo !...(1)ie = 1,nelem

end subroutine getgeoel_lag_post
!
!...Get the nodal values...
!
subroutine getpostnd_lag_hybrid(unkno,inpoel,iptri, ipqua, coord,geoel,unknp)
use constant
implicit none
!...Input array
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
integer*4,dimension(1:nvtri,1:ntria),  intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria),  intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad),  intent(in):: ipqua
real*8,dimension(1:ndimn,1:npoin),     intent(in)::coord
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
!Loacal
real*8,dimension(1:nvtri)::xv, yv
real*8,dimension(1:nvqua)::xvq, yvq
real*8,dimension(1:ndegr,1:nvtri)::bt
real*8,dimension(1:ndegr, 1:nvqua)::bq
real*8,dimension(1:npoin, 1:nq), intent(out)::unknp
real*8,  dimension(1:npoin)::vlsup
integer, dimension(1:nvtri):: ipt
integer, dimension(1:nvqua):: ipq

integer::iv, ideg,ipoin, ie, ielem
real*8::rho, uadv, vadv, pres,eadv,volel
real*8::rc,sc,dr,ds
!
xv(1) = 0.d0; yv(1) = 0.d0
xv(2) = 1.d0; yv(2) = 0.d0
xv(3) = 0.d0; yv(3) = 1.d0
!
 if(ncurv.eq.1)then
!
xv(4) = 0.5d0; yv(4) = 0.d0
xv(5) = 0.5d0; yv(5) = 0.5d0
xv(6) = 0.d0;  yv(6) = 0.5d0
!
 endif
!
xvq(1) = -1.d0; yvq(1) = -1.d0
xvq(2) =  1.d0; yvq(2) = -1.d0
xvq(3) =  1.d0; yvq(3) =  1.d0
xvq(4) = -1.d0; yvq(4) =  1.d0
!
unknp = 0.d0
vlsup = 0.d0
!
!calculate the point value
!
do ipoin=1,npoin
!
!...Triangle...
!
do ie=1,ntria
!
ielem = ie
!
dr = 0.5d0
ds = 0.5d0
!
rc = geoel(1, ielem)
sc = geoel(2, ielem)
volel = geoel(3, ielem)
!
ipt(1:nvtri) = iptri(1:nvtri,ie)
!
do iv = 1, nvtri
!
bt(1, iv) = 1.d0
bt(2, iv) = (xv(iv)-rc)/dr
bt(3, iv) = (yv(iv)-sc)/ds
!
if(ipt(iv).eq.ipoin)then
!
do ideg = 1,ndegr
unknp(ipoin, 1:nq) = unknp(ipoin, 1:nq) + unkno(ideg, 1:nq, ielem)*bt(ideg,iv)*volel/(nvtri*1.d0)
enddo
vlsup(ipoin) = vlsup(ipoin) + volel/(nvtri*1.d0)
endif
enddo
!
enddo
!
!...Quads...
!
do ie=1,nquad
!
ielem = ie + ntria
!
dr = 1.d0
ds = 1.d0
!
rc = geoel(1, ielem)
sc = geoel(2, ielem)
volel = geoel(3, ielem)
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
do iv = 1, nvqua
!
bq(1, iv) = 1.d0
!
if(npoly.ge.1)then
bq(2, iv) = (xvq(iv)-rc)/dr
bq(3, iv) = (yvq(iv)-sc)/ds
endif
!
if(ipq(iv).eq.ipoin)then
!
do ideg = 1,ndegr
unknp(ipoin, 1:nq) = unknp(ipoin, 1:nq) + unkno(ideg, 1:nq, ielem)*bq(ideg,iv)*volel/(nvqua*1.d0)
enddo
vlsup(ipoin) = vlsup(ipoin) + volel/(nvqua*1.d0)

endif
enddo
!
enddo
!
enddo !ipoin=1,npoin
!
!...Get nodal value
!
do ipoin = 1, npoin
unknp(ipoin, 1:nq) = unknp(ipoin, 1:nq)/vlsup(ipoin)
!
rho = 1.d0/unknp(ipoin, 1)
uadv = unknp(ipoin, 2)
vadv = unknp(ipoin, 3)
eadv = unknp(ipoin, 4)
!
pres = (gamlg-1.d0)*rho*(eadv-0.5d0*(uadv*uadv+vadv*vadv))
!
unknp(ipoin, 1) = rho
unknp(ipoin, 2) = uadv
unknp(ipoin, 3) = vadv
unknp(ipoin, 4) = pres
!
enddo
!
!output the valu
return
! close(11)
end subroutine getpostnd_lag_hybrid
!
!...L2 error of linear hybrid meshes in Lagrangian
!
subroutine geterror_lag_hybrid(intfac, inpoel,iptri, ipqua, coord, coold, geoel, unkno)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
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
real*8,dimension(1:ndimn, 1:nvtri) :: xp,xo
real*8,dimension(1:ndimn, 1:nvqua) :: xpq,xoq
real*8,dimension(1:ndegr):: b,bv
real*8:: unknv(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s,  dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc,rm, sm,rp,sp
real*8::xg, yg, xc, yc, xgo, ygo
real*8::rho,uadv,vadv,eadv,rhom,rho0ba,rho0
real*8::rhoc, uctr, vctr, ectr, pctr
real*8::pres, pexa, eiadv, eiexa, rhoexa
real*8::djaco, wi
real*8::radie, radii,radie2,radii2, radig2,sentr,rhoin,rhoex,htkid
real*8::htkid2,rhoexa2,pexa2,tmend2
real*8::rhon, rhoini
real*8::rcv, scv
real*8::rexp,veexa,tmend
real*8::uex,uexa,vexa,vemag
real*8::prein,preex
real*8::errorl2(3)
!
real*8::errlp(nq),errl2(nq)!...Error array
!
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8::weigh(ngausd), posi(2,ngausd)
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
open(8,file='unradi.dat')  
!
errorl2 = 0.d0
errl2 = 0.d0
errlp = 0.d0
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
xo(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xo(2, 1:nvtri) = coold(2, ipt(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
rc = geoel(1, ielem)
sc = geoel(2, ielem)
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
xg =0.d0
yg = 0.d0
!
do ishp = 1, nvtri
xg = xg + shp(ishp)*xp(1, ishp)
yg = yg + shp(ishp)*xp(2, ishp)
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
!
unknv = 0.d0
!
do ideg = 1,mdegr
unknv(1:nq) = unknv(1:nq) + unkno(ideg,1:nq,ielem)*b(ideg)
enddo
!
!...Jacobian transformation matrix
!
if(ndens.eq.1)then
 rho  = 1.d0/unknv(1)
elseif(ndens.eq.2)then
!
call  getrhoig_tria(rhoini, r, s, xo)
!
call getdensity_triallnl(r, s, xp, xo, rhoini, rhon)
!
 rho = rhon
!
elseif(ndens.eq.3)then
rho  = unknv(1)
endif
uadv = unknv(2)
vadv = unknv(3)
eadv = unknv(4)

!
if(ncase.eq.1)then
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
pexa = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 1.d0
!
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
elseif(ncase.eq.2)then
!
eiadv = eadv - 0.5d0*(uadv**2 + vadv**2)
eiexa = (1.d0-0.6d0)**(-2.d0*(gamlg-1.d0))
!
!    rhoexa = (1.d0-0.4d0)**(-2.d0)
!
errl2 = errl2 + (eiadv - eiexa)**2*djaco
!    errl2 = errl2 + (rhoexa - rho)**2*djaco
!
elseif(ncase.eq.4)then
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!...some parameters
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
xgo =0.d0
ygo = 0.d0
!
xgo = shp(1)*xo(1, 1) + shp(2)*xo(1, 2) + shp(3)*xo(1, 3)
ygo = shp(1)*xo(2, 1) + shp(2)*xo(2, 2) + shp(3)*xo(2, 3)
!
radig2 = xgo**2 + ygo**2
!
rho0ba =(radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0)
!
!rho = rho0ba**(1.d0/(gamlg-1.d0))
!
htkid = (sqrt(1.d0-(.5d0)**2))**(-2.d0*gamlg/(gamlg-1.d0))
pexa= sentr*(rho)**gamlg*htkid
!
!errl2 = errl2 + (pres -pexa)**2*djaco
errl2 = errl2 + (pres/sentr/rho**gamlg -1.d0)**2*djaco*yg
!
!
endif
!
!...Primitive variables...
!
enddo !...(2)ig = 1,ngausd
!
if(ncase.eq.7.or.ncase.eq.6.or.ncase.eq.3.or.ncase.eq.8)then
!
if(ndens.eq.1)then
rhoc=1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhoc=unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhoc=unkno(1, 1, ielem)
endif
uctr=unkno(1, 2, ielem)
vctr=unkno(1, 3, ielem)
ectr=unkno(1, 4, ielem)
!
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
r = rc; s =sc

shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
!
xc =0.d0
yc = 0.d0
!
do ishp = 1, nvtri
xc = xc + shp(ishp)*xp(1, ishp)
yc = yc + shp(ishp)*xp(2, ishp)
enddo
!
!if(ncase.eq.6) write(8,'(i8, 4e32.16)')ielem,sqrt(xc**2 + (yc-0.d0)**2),rhoc,pctr,pctr/rhoc/(gamlg-1.d0)
if(ncase.eq.6.or.ncase.eq.7.or.ncase.eq.3.or.ncase.eq.8) write(8,'(i8, 8e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,&
                                           pctr/rhoc/(gamlg-1.d0),xc,yc,uctr,vctr
!
endif
!
550 enddo
!
!...For quads...
!
!
!...Loop over elements
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
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
xoq(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xoq(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
rc = geoel(1, ielem)
sc = geoel(2, ielem)
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
xg =0.d0
yg = 0.d0
!
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpq(1, ishp)
yg = yg + shpq(ishp)*xpq(2, ishp)
enddo
!
djaco = wi*(dxdr*dyds - dydr*dxds)
!
b(1) = 1.d0

if(npoly.ge.1)then
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
endif

!DGP2
if(npoly.eq.2)then
b(4) = 0.5d0*b(2)*b(2) - geoel(19, ielem)
b(5) = 0.5d0*b(3)*b(3) - geoel(20, ielem)
b(6) =       b(2)*b(3) - geoel(21, ielem)
endif
!
unknv = 0.d0
!
do ideg = 1,mdegr
unknv(1:nq) = unknv(1:nq) + unkno(ideg,1:nq,ielem)*b(ideg)
enddo
!
!...Jacobian transformation matrix
!
if(ndens.eq.1)then
rho  = 1.d0/unknv(1)
elseif(ndens.eq.2)then
!
call getrhoig_quad(rhoini, r, s, xoq)!
call getdensity_quadllnl(r, s, xpq, xoq, rhoini, rhon)
!
rho = rhon
!
elseif(ndens.eq.3)then
rcv = geoel(5, ielem)
scv = geoel(6, ielem)
!
!
bv(1) = 1.d0
bv(2) = (r-rcv)/dr
bv(3) = (s-scv)/ds
!
unknv(1) = 0.d0
!
do ideg = 1,mdegr
unknv(1) = unknv(1) + unkno(ideg,1,ielem)*bv(ideg)
enddo
!
rho  = unknv(1)
!
endif
uadv = unknv(2)
vadv = unknv(3)
eadv = unknv(4)
!
if(ncase.eq.1)then
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
pexa = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 1.d0
uexa   = sin(pi*xg)*cos(pi*yg)
vexa   =-cos(pi*xg)*sin(pi*yg)
rhoexa = 1.d0
!
vemag = sqrt(uadv**2 + vadv**2)
veexa = sqrt(uexa**2 + vexa**2)
!
errlp(1) = errlp(1) + abs(rho-rhoexa)**2*djaco
errlp(2) = errlp(2) + abs(vemag-veexa)**2*djaco
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
!
elseif(ncase.eq.2)then
!
eiadv = eadv - 0.5d0*(uadv**2 + vadv**2)
rhoexa = (1.d0-0.6d0)**(-2.d0)
uex = -xg/(1.d0-0.6d0)
eiexa = (1.d0-0.6d0)**(-2.d0*(gamlg-1.d0))
!
!    rhoexa = (1.d0-0.4d0)**(-2.d0)
!
errl2 = errl2 + (eiadv - eiexa)**2*djaco
!
errorl2(1) = errorl2(1) + (rho - rhoexa)**2*djaco
errorl2(2) = errorl2(2) + (uex - uadv)**2*djaco
errorl2(3) = errorl2(3) + (eadv -eiexa-uadv**2/2.d0-vadv**2/2.d0)**2*djaco
!
elseif(ncase.eq.4)then
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!...some parameters
!
radie = 1.0d0
radii = 0.9d0
!rhoin = 6.31d-4
rhoex = 1.d-2
!sentr = 2.15d4
prein = 0.1d0
preex = 10.d0
tmend = 0.8d0
tmend2 = 0.5d0
!
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
!
xgo =0.d0
ygo = 0.d0
!
do ishp = 1, nvqua
xgo = xgo + shpq(ishp)*xoq(1, ishp)
ygo = ygo + shpq(ishp)*xoq(2, ishp)
enddo
!
radig2 = xgo**2 + ygo**2
!
rho0ba =(radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
     (radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0)
rho0 = rho0ba**(1.d0/(gamlg-1.d0))

htkid = sqrt(1.d0-(tmend)**2)
htkid2 = sqrt(1.d0-(tmend2)**2)
!
rhoexa = htkid**(-2.d0/(gamlg-1.d0))*rho0
pexa= sentr*(htkid)**(-2.d0*gamlg/(gamlg-1.d0))*(rho0)**gamlg
!
rhoexa2 = htkid2**(-2.d0/(gamlg-1.d0))*rho0
pexa2= sentr*(htkid2)**(-2.d0*gamlg/(gamlg-1.d0))*(rho0)**gamlg
!
 errlp(1) = errlp(1) + (rhoexa -rho)**2*djaco*yg
 errlp(2) = errlp(2) + (pres/sentr/rho**gamlg -1.d0)**2*djaco*yg!geoel(11, ielem)!yg
 errlp(3) = errlp(3) + (pexa -pres)**2*djaco*yg

 errl2(1) = errl2(1) + (rhoexa2 -rho)**2*djaco*yg
 errl2(2) = errl2(2) + (pres/sentr/rho**gamlg -1.d0)**2*djaco*yg!geoel(11, ielem)!yg
 errl2(3) = errl2(3) + (pexa2 -pres)**2*djaco*yg
!
!endif
!
!
elseif(ncase.eq.10)then !...Expansion
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
tmend = 0.4d0
rexp = sqrt(1.d0 + 2.d0*tmend**2)
rhoexa = 1.d0/rexp**3
veexa = 2.d0*tmend/rexp**2*sqrt(xg**2 + yg**2)
pexa = (1.d0-(xg**2 + yg**2)/rexp**2)/rexp**5

if(geoel(10, ielem).gt.100)then
errlp(1) = errlp(1) + (rhoexa-rho)**2*djaco*yg
errlp(2) = errlp(2) + (rhoexa*veexa-rho*sqrt(uadv**2 + vadv**2))**2*djaco*yg
errlp(3) = errlp(3) + (pexa-pres)**2*djaco*yg
endif
!
elseif(ncase.eq.12)then !...1D isentropic expansion...
!
call solve_isentropic(500,xg,0.1d0,rhoexa,uex,pexa)
!
errl2 = errl2 + (rho - rhoex)**2*djaco
errorl2(1) = errorl2(1) + (rho - rhoexa)**2*djaco
errorl2(2) = errorl2(2) + (uex - uadv)**2*djaco
errorl2(3) = errorl2(3) + (eadv -pexa/(gamlg-1.d0)/rhoexa-0.5d0*uex**2)**2*djaco
!
!print*,'Focusing time for Kidder',sqrt((gamlg-1.d0)/sentr/gamlg/2.d0*(radie2-radii2)/(rhoex**(gamlg-1.d0)-rhoin**(gamlg-1.d0)))
!
elseif(ncase.eq.14)then !...Coggeshall problem...
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
tmend = 0.5d0
rhoexa = 1.d0/(1.d0-tmend)**2.25d0

!if(geoel(10, ielem).gt.100)then
errlp(1) = errlp(1) + (rhoexa-rho)**2*djaco*geoel(11, ielem)!yg
!errlp(2) = errlp(2) + (rhoexa*veexa-rho*sqrt(uadv**2 + vadv**2))**2*djaco*yg
!errlp(3) = errlp(3) + (pexa-pres)**2*djaco*yg
!endif
!
endif
!
!...Primitive variables...
!
enddo !...(2)ig = 1,ngausd
!
if(ncase.eq.7.or.ncase.eq.6.or.ncase.eq.3.or.ncase.eq.4.or.ncase.eq.8.or.ncase.eq.10 &
   .or.ncase.eq.11.or.ncase.eq.12.or.ncase.eq.13.or.ncase.eq.14)then
!
if(ndens.eq.1)then
rhoc=1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhoc=unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhoc=unkno(1, 1, ielem)
endif
uctr=unkno(1, 2, ielem)
vctr=unkno(1, 3, ielem)
ectr=unkno(1, 4, ielem)
!
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
r = rc; s =sc

r = 0.d0; s =0.d0
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
xc =0.d0
yc = 0.d0
!
do ishp = 1, nvqua
xc = xc + shpq(ishp)*xpq(1, ishp)
yc = yc + shpq(ishp)*xpq(2, ishp)
enddo
!
!if(rhoc.gt.1.d0)print*,'bad cell',ielem,sqrt(xc**2+yc**2),rhoc,unkno(1, 1, ielem)
!if(ncase.eq.6) write(8,'(i8, 4e32.16)')ielem,sqrt(xc**2 + (yc-0.0d0)**2),rhoc,pctr,pctr/rhoc/(gamlg-1.d0)

if(ncase.eq.6.or.ncase.eq.7.or.ncase.eq.8.or.ncase.eq.3.or.ncase.eq.4.or.ncase.eq.10.or.ncase.eq.11.or.ncase.eq.14)then
  write(8,'(i8, 9e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,&
                                           pctr/rhoc/(gamlg-1.d0),xc,yc,sqrt(uctr**2+vctr**2),uctr,vctr
endif

if(ncase.eq.12) then
  write(8,'(i8, 7e32.16)')ielem,xc,rhoc,pctr,&
pctr/rhoc/(gamlg-1.d0),sqrt(uctr**2+vctr**2),uctr,vctr
!write(8,'(4e32.16)')xc,rhoexa,uex,pexa
endif

if(ncase.eq.13) write(8,'(i8, 8e32.16)')ielem,sqrt(xc**2),rhoc,pctr,&
pctr/rhoc/(gamlg-1.d0),xc,yc,uctr,vctr
!
endif
!
650 enddo
!
close(8)
!
if(ncase.eq.1)  print*,'L2 error for Taylor-Green-Vortex',0.5d0*log10(errlp),sqrt(errlp)
if(ncase.eq.2)  print*,'L2 error for Shockless Noh',0.5d0*log10(errl2),sqrt(errl2)!0.5d0*log10(errorl2),sqrt(errorl2)
if(ncase.eq.4)  print*,'L2 error for Kidder shell',0.5d0*log10(errlp),sqrt(errlp)
if(ncase.eq.4)  print*,'L2 error for Kidder shell0.5',0.5d0*log10(errl2),sqrt(errl2)
if(ncase.eq.10)  print*,'L2 error for Kidder shell',0.5d0*log10(errlp(1:3)),sqrt(errlp(1:3))
if(ncase.eq.12)  print*,'L2 error for 1D isentropic sin wave flow',0.5d0*log10(errorl2),sqrt(errorl2)
if(ncase.eq.14)  print*,'L2 error for Coggeshall',0.5d0*log10(errlp)
!
end subroutine geterror_lag_hybrid
!
!...L2 error of curved hybrid meshes in Lagrangian
!
subroutine geterror_lag_curvhybrid(intfac, inpoel, ipqua, coord, coold, unkno, geoel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem,nind
!
!...local integer array
!
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpq,xoq
real*8,dimension(1:ndimn, 1:nvtri) :: coorp
real*8,dimension(1:ndegr):: b,bv
real*8:: unknv(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq_pst), posiq(2,ngausdq_pst)
real*8::dhgre(2), dhgr2(3)

!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s,  dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc,xc,yc
real*8::xg, yg, xgo, ygo
real*8::rho,uadv,vadv,eadv,rhom
real*8::pres, pexa, eiadv, eiexa
real*8::djaco, wi, errl2
real*8::rhoexa,uexa,vexa,veexa,vemag,uex
real*8::rhoc, uctr, vctr, ectr, pctr
real*8::radie, radii,radie2,radii2, radig2,paras,rhoi,rhoe,htkid
real*8::preex, prein, rho0, rho0ba, rhoex, rhoin, sentr, tmend
real*8::rhon, rhoini
real*8::rcv, scv
real*8::hgre,radi,radic,radiv,ratir
!
real*8:: errlp(nq)
!
real*8,allocatable::weigh(:), posi(:,:)
!
data eps   / 1.0d-14 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
allocate (weigh(ngausd), posi(2,ngausd))
call rutope(2, ngausd, posi, weigh)
call ruqope(2, ngausdq_pst, posiq, weighq)
!
errl2 = 0.d0
errlp = 0.d0
!
open(8,file='unradicurv.dat')
!
!...Loop over elements
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
!...Points consitituting one element...
!
ielem = ie
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
coorp(1, 1:nvtri) = coord(1, ip(1:nvtri))
coorp(2, 1:nvtri) = coord(2, ip(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
!...Gauss loop
!
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
wi     = weigh(ig)
!
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
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
!DGP2
if(npoly.eq.2)then
b(4) = 0.5d0*b(2)*b(2) - geoel(19, ielem)
b(5) = 0.5d0*b(3)*b(3) - geoel(20, ielem)
b(6) =       b(2)*b(3) - geoel(21, ielem)
endif
!
unknv = 0.d0
!
do ideg = 1,mdegr
unknv(1:nq) = unknv(1:nq) + unkno(ideg,1:nq,ie)*b(ideg)
enddo
!
!...Jacobian transformation matrix
!
rho  = 1.d0/unknv(1)
uadv = unknv(2)
vadv = unknv(3)
eadv = unknv(4)
!
if(ncase.eq.1)then
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
pexa = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 1.d0
!
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
elseif(ncase.eq.2)then
!
eiadv = eadv - 0.5d0*(uadv**2 + vadv**2)
eiexa = (1.d0-0.6d0)**(-2.d0*(gamlg-1.d0))
!
!    rhoexa = (1.d0-0.6d0)**(-2.d0)
!
errl2 = errl2 + (eiadv - eiexa)**2*djaco
!
!print*,'errl2',errl2,(eiadv - eiexa)**2*djaco,(eiadv - eiexa)**2,djaco,ie
!    errl2 = errl2 + (rhoexa - rho)**2*djaco
!
endif
!
!...Primitive variables...
!
enddo !...(2)ig = 1,ngausd
!
if(ncase.eq.7.or.ncase.eq.6.or.ncase.eq.3.or.ncase.eq.8)then
!
if(ndens.eq.1)then
rhoc=1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhoc=unkno(1, 1, ielem)
endif
uctr=unkno(1, 2, ielem)
vctr=unkno(1, 3, ielem)
ectr=unkno(1, 4, ielem)
!
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
r = rc; s =sc
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
xc =0.d0
yc = 0.d0
!
do ishp = 1, nvtri
xc = xc + shp(ishp)*coorp(1, ishp)
yc = yc + shp(ishp)*coorp(2, ishp)
enddo
!
!if(rhoc.gt.1.d0)print*,'bad cell',ielem,sqrt(xc**2+yc**2),rhoc,unkno(1, 1, ielem)
if(ncase.eq.6) write(8,'(i8, 4e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,pctr/rhoc/(gamlg-1.d0)
if(ncase.eq.7.or.ncase.eq.3.or.ncase.eq.8) write(8,'(i8, 7e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,&
                                           pctr/rhoc/(gamlg-1.d0),sqrt(uctr**2+vctr**2),uctr,vctr
endif
!
550 enddo
!
!
!...For quads...
!
!
!...Loop over elements
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xoq(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xoq(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
rc = geoel(1, ielem)
sc = geoel(2, ielem)
!
!...Gauss loop
!
do ig = 1,ngausdq_pst !...(2)ig = 1,ngausd
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
if(djaco/geoel(4, ielem).le.1.d-6) print*,'Negative Jacobian for curved cells',djaco,ielem
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
!
unknv = 0.d0
!
do ideg = 1,mdegr
unknv(1:nq) = unknv(1:nq) + unkno(ideg,1:nq,ielem)*b(ideg)
enddo
!
!...Jacobian transformation matrix
!
if(ndens.eq.1)then
rho  = 1.d0/unknv(1)
elseif(ndens.eq.2)then
!
call getrhoig_quadcurv(rhoini, xoq)!
call getdensity_quadllnl_curv(r, s, xpq, xoq, rhoini, rhon)
!
rho = rhon
!
elseif(ndens.eq.3)then
rcv = geoel(5, ielem)
scv = geoel(6, ielem)
!
!
bv(1) = 1.d0
bv(2) = (r-rcv)/dr
bv(3) = (s-scv)/ds
!
unknv(1) = 0.d0
!
do ideg = 1,mdegr
unknv(1) = unknv(1) + unkno(ideg,1,ielem)*bv(ideg)
enddo
!
rho  = unknv(1)
!
endif
uadv = unknv(2)
vadv = unknv(3)
eadv = unknv(4)
!
if(ncase.eq.1)then
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
pexa = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 1.d0
uexa   = sin(pi*xg)*cos(pi*yg)
vexa   =-cos(pi*xg)*sin(pi*yg)
rhoexa = 1.d0
!
vemag = sqrt(uadv**2 + vadv**2)
veexa = sqrt(uexa**2 + vexa**2)
!
errlp(1) = errlp(1) + abs(rho-rhoexa)**2*djaco
errlp(2) = errlp(2) + abs(vemag-veexa)**2*djaco
!errlp(2) = errlp(2) + abs(uadv-uexa)**2*djaco
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
!
!if(ielem.eq.313)print*,'TGV error',errlp(3),abs(pres - pexa)**2*djaco,ielem,unkno(:,:,313)
!
elseif(ncase.eq.2)then
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
eiadv = eadv - 0.5d0*(uadv**2 + vadv**2)
eiexa = (1.d0-0.6d0)**(-2.d0*(gamlg-1.d0))
rhoexa = (1.d0-0.6d0)**(-2.d0)
!
pexa = (gamlg-1.d0)*rhoexa*eiadv
!
errlp(1) = errlp(1) + abs(rho-rhoexa)**2*djaco
errlp(2) = errlp(2) + abs(eiadv - eiexa)**2*djaco
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
!
!if(ielem.eq.1) print*,'errlp',errlp(2),unkno(1:6,4,ielem),ielem
!
elseif(ncase.eq.4)then
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!...some parameters
!
radie = 1.0d0
radii = 0.9d0
!rhoin = 6.31d-4
rhoex = 1.d-2
!sentr = 2.15d4
prein = 0.1d0
preex = 10.d0
tmend = 0.5d0
!
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
!
xgo =0.d0
ygo = 0.d0
!
do ishp = 1, nvqua
xgo = xgo + shpq(ishp)*xoq(1, ishp)
ygo = ygo + shpq(ishp)*xoq(2, ishp)
enddo
!
radig2 = xgo**2 + ygo**2
!
rho0ba =(radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0)
rho0 = rho0ba**(1.d0/(gamlg-1.d0))

htkid = sqrt(1.d0-(tmend)**2)
!
rhoexa = htkid**(-2.d0/(gamlg-1.d0))*rho0
pexa= sentr*(htkid)**(-2.d0*gamlg/(gamlg-1.d0))*(rho0)**gamlg
!
errlp(1) = errlp(1) + (rhoexa -rho)**2*djaco
errlp(2) = errlp(2) + (pres/sentr/rho**gamlg -1.d0)**2*djaco
!errlp(3) = errlp(3) + (pexa -pres)**2*djaco
errlp(3) = errlp(3) + (pexa -pres)**2*djaco
!...Gresho
elseif(ncase.eq.8)then
!
nind = 6.d0
radiv= 0.4d0
radi = sqrt(xg**2+yg**2)
radic = radi
!
ratir = radi/radiv
!
vemag = sqrt(uadv**2 + vadv**2)
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))

if(radi.lt.radiv)then
veexa = 2**(2*nind)*ratir**nind*(1.d0-ratir)**nind
else
veexa = 0.d0
endif
!
call greshoh(xg, yg ,radic, hgre, dhgre, dhgr2)
!
rhoexa =1.d0
pexa = 5.d0 + 2**(4*nind)*hgre
!
errlp(1) = errlp(1) + abs(1.d0/rho-1.d0/rhoexa)**2*djaco
errlp(2) = errlp(2) + abs(vemag-veexa)**2*djaco
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
errlp(4) = errlp(4) + abs(eadv -pexa/(gamlg-1.d0)/rhoexa-0.5d0*veexa**2)**2*djaco
!
elseif(ncase.eq.12)then !...1D isentropic expansion...
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
call solve_isentropic(1000,xg,0.d0,rhoexa,uex,pexa)
!
errlp(1) = errlp(1) + abs(1.d0/rho-1.d0/rhoexa)**2*djaco
errlp(2) = errlp(2) + abs(uex - uadv)**2*djaco
errlp(3) = errlp(3) + abs(eadv -pexa/(gamlg-1.d0)/rhoexa-0.5d0*uex**2)**2*djaco
errlp(4) = errlp(3) + abs(pres -pexa)**2*djaco
!
endif
!
enddo !...(2)ig = 1,ngausd
!
!...Output cell centered values to the file ''
!
if(ncase.eq.7.or.ncase.eq.6.or.ncase.eq.3.or.ncase.eq.4.or.ncase.eq.8.or.ncase.eq.12.or.ncase.eq.13.or.ncase.eq.15)then
!
if(ndens.eq.1)then
rhoc=1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhoc=unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhoc=unkno(1, 1, ielem)
endif
uctr=unkno(1, 2, ielem)
vctr=unkno(1, 3, ielem)
ectr=unkno(1, 4, ielem)
!
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
!r = rc; s =sc
r = 0.d0; s =0.d0
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
xc =0.d0
yc = 0.d0
!
do ishp = 1, nvqua
xc = xc + shpq(ishp)*xpq(1, ishp)
yc = yc + shpq(ishp)*xpq(2, ishp)
enddo
!
if(ncase.eq.6.or.ncase.eq.7.or.ncase.eq.3.or.ncase.eq.4.or.ncase.eq.10.or.ncase.eq.11.or.ncase.eq.14.or.ncase.eq.8)then
 write(8,'(i8, 9e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,&
                         pctr/rhoc/(gamlg-1.d0),xc,yc,sqrt(uctr**2+vctr**2),uctr,vctr
endif

!...Isentropic smooth flow and Saltzman
if(ncase.eq.12.or.ncase.eq.13.or.ncase.eq.15)then
 write(8,'(i8, 7e32.16)')ielem,xc,rhoc,pctr,&
                         pctr/rhoc/(gamlg-1.d0),sqrt(uctr**2+vctr**2),uctr,vctr
endif
!
endif
!
650 enddo
!
close(8)

!...Show the L2 error on the screen
 print*,'L2 error',0.5d0*log10(errlp),sqrt(errlp)
!
end subroutine geterror_lag_curvhybrid
!
!...L2 error of curved hybrid meshes in Lagrangian
!
subroutine geterror_lag_curvhybrid_general(intfac, inpoel, ipqua, coord, coold, unkno, geoel)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
!
!...Local integer
!
integer::ie,ig,ideg,ishp,iv,ielem,nind
!
!...local integer array
!
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpq,xoq
real*8,dimension(1:ndimn, 1:nvtri) :: coorp
real*8,dimension(1:ndegr):: b,bv
real*8:: unknv(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq_pst), posiq(2,ngausdq_pst)
real*8::dhgre(2), dhgr2(3)

!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s,  dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc,xc,yc
real*8::xg, yg, xgo, ygo
real*8::rho,uadv,vadv,eadv,rhom
real*8::pres, pexa, eiadv, eiexa
real*8::djaco, wi, errl2
real*8::rhoexa,uexa,vexa,veexa,vemag,uex
real*8::rhoc, uctr, vctr, ectr, pctr
real*8::radie, radii,radie2,radii2, radig2,paras,rhoi,rhoe,htkid
real*8::preex, prein, rho0, rho0ba, rhoex, rhoin, sentr, tmend
real*8::rhon, rhoini
real*8::rcv, scv
real*8::hgre,radi,radic,radiv,ratir
!
real*8:: errlp(nq)
!
real*8,allocatable::weigh(:), posi(:,:)
!
data eps   / 1.0d-14 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
allocate (weigh(ngausd), posi(2,ngausd))
call rutope(2, ngausd, posi, weigh)
call ruqope(2, ngausdq_pst, posiq, weighq)
!
errl2 = 0.d0
errlp = 0.d0
!
open(8,file='unradicurv.dat')
!
!...Loop over elements
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
!...Points consitituting one element...
!
ielem = ie
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
coorp(1, 1:nvtri) = coord(1, ip(1:nvtri))
coorp(2, 1:nvtri) = coord(2, ip(1:nvtri))
!
!...Geometry parameters for reference cell...
!
dr = .5d0
ds = .5d0
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)
!
!...Gauss loop
!
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
wi     = weigh(ig)
!
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
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
!DGP2
if(npoly.eq.2)then
b(4) = 0.5d0*b(2)*b(2) - geoel(19, ielem)
b(5) = 0.5d0*b(3)*b(3) - geoel(20, ielem)
b(6) =       b(2)*b(3) - geoel(21, ielem)
endif
!
unknv = 0.d0
!
do ideg = 1,mdegr
unknv(1:nq) = unknv(1:nq) + unkno(ideg,1:nq,ie)*b(ideg)
enddo
!
!...Jacobian transformation matrix
!
rho  = 1.d0/unknv(1)
uadv = unknv(2)
vadv = unknv(3)
eadv = unknv(4)
!
if(ncase.eq.1)then
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
pexa = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 1.d0
!
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
elseif(ncase.eq.2)then
!
eiadv = eadv - 0.5d0*(uadv**2 + vadv**2)
eiexa = (1.d0-0.6d0)**(-2.d0*(gamlg-1.d0))
!
!    rhoexa = (1.d0-0.6d0)**(-2.d0)
!
errl2 = errl2 + (eiadv - eiexa)**2*djaco
!
!print*,'errl2',errl2,(eiadv - eiexa)**2*djaco,(eiadv - eiexa)**2,djaco,ie
!    errl2 = errl2 + (rhoexa - rho)**2*djaco
!
endif
!
!...Primitive variables...
!
enddo !...(2)ig = 1,ngausd
!
if(ncase.eq.7.or.ncase.eq.6.or.ncase.eq.3.or.ncase.eq.8)then
!
if(ndens.eq.1)then
rhoc=1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhoc=unkno(1, 1, ielem)
endif
uctr=unkno(1, 2, ielem)
vctr=unkno(1, 3, ielem)
ectr=unkno(1, 4, ielem)
!
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
r = rc; s =sc
!
shp(1) = -(c10-r-s)*(c10-2.d0*(c10-r-s))
shp(2) = -r*(c10-2.d0*r)
shp(3) = -s*(c10-2.d0*s)
shp(4) = 4.d0*r*(c10-r-s)
shp(5) = 4.d0*r*s
shp(6) = 4.d0*s*(c10-r-s)
!
xc =0.d0
yc = 0.d0
!
do ishp = 1, nvtri
xc = xc + shp(ishp)*coorp(1, ishp)
yc = yc + shp(ishp)*coorp(2, ishp)
enddo
!
!if(rhoc.gt.1.d0)print*,'bad cell',ielem,sqrt(xc**2+yc**2),rhoc,unkno(1, 1, ielem)
if(ncase.eq.6) write(8,'(i8, 4e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,pctr/rhoc/(gamlg-1.d0)
if(ncase.eq.7.or.ncase.eq.3.or.ncase.eq.8) write(8,'(i8, 7e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,&
pctr/rhoc/(gamlg-1.d0),sqrt(uctr**2+vctr**2),uctr,vctr
endif
!
550 enddo
!
!
!...For quads...
!
!
!...Loop over elements
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
xoq(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xoq(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
!...Get the points for FEM
call getcoord_fe(ncurv, nvfac, nvqua, xpq)
!
!...Get the points for FEM
call getcoord_fe(ncurv, nvfac, nvqua, xoq)
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
rc = geoel(1, ielem)
sc = geoel(2, ielem)
!
!...Gauss loop
!
do ig = 1,ngausdq_pst !...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
wi  = weighq(ig)
!
!...  shape function & its derivatives w.r.t. reference coordinates
!
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
!if(djaco/geoel(4, ielem).le.1.d-6) print*,'Negative Jacobian for curved cells',djaco,ielem
!
!...Gauss points...
!
xg =0.d0; yg= 0.d0
!
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpq(1,ishp)
yg = yg + shpq(ishp)*xpq(2,ishp)
enddo
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
!
unknv = 0.d0
!
do ideg = 1,mdegr
unknv(1:nq) = unknv(1:nq) + unkno(ideg,1:nq,ielem)*b(ideg)
enddo
!
!...Jacobian transformation matrix
!
if(ndens.eq.1)then
rho  = 1.d0/unknv(1)
elseif(ndens.eq.2)then
!
call getrhoig_quadcurv(rhoini, xoq)!
call getdensity_quadllnl_curv(r, s, xpq, xoq, rhoini, rhon)
!
rho = rhon
!
elseif(ndens.eq.3)then
rcv = geoel(5, ielem)
scv = geoel(6, ielem)
!
!
bv(1) = 1.d0
bv(2) = (r-rcv)/dr
bv(3) = (s-scv)/ds
!
unknv(1) = 0.d0
!
do ideg = 1,mdegr
unknv(1) = unknv(1) + unkno(ideg,1,ielem)*bv(ideg)
enddo
!
rho  = unknv(1)
!
endif
uadv = unknv(2)
vadv = unknv(3)
eadv = unknv(4)
!
if(ncase.eq.1)then
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
pexa = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 1.d0
uexa   = sin(pi*xg)*cos(pi*yg)
vexa   =-cos(pi*xg)*sin(pi*yg)
rhoexa = 1.d0
!
vemag = sqrt(uadv**2 + vadv**2)
veexa = sqrt(uexa**2 + vexa**2)
!
errlp(1) = errlp(1) + abs(rho-rhoexa)**1*djaco
!errlp(2) = errlp(2) + abs(vemag-veexa)**2*djaco
errlp(2) = errlp(2) + abs(uadv-uexa)**2*djaco
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
errlp(4) = errlp(4) + abs(pexa/(gamlg-1.d0)+0.5d0*veexa**2 - eadv)**2*djaco
!
!if(ielem.eq.313)print*,'TGV error',errlp(3),abs(pres - pexa)**2*djaco,ielem,unkno(:,:,313)
!
elseif(ncase.eq.2)then
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
eiadv = eadv - 0.5d0*(uadv**2 + vadv**2)
eiexa = (1.d0-0.6d0)**(-2.d0*(gamlg-1.d0))
rhoexa = (1.d0-0.6d0)**(-2.d0)
!
pexa = (gamlg-1.d0)*rhoexa*eiadv
!
errlp(1) = errlp(1) + abs(rho-rhoexa)**2*djaco
errlp(2) = errlp(2) + abs(eiadv - eiexa)**2*djaco
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
!
!if(ielem.eq.1) print*,'errlp',errlp(2),unkno(1:6,4,ielem),ielem
!
elseif(ncase.eq.4)then
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!...some parameters
!
radie = 1.0d0
radii = 0.9d0
!rhoin = 6.31d-4
rhoex = 1.d-2
!sentr = 2.15d4
prein = 0.1d0
preex = 10.d0
tmend = 0.5d0
!
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
!
xgo =0.d0
ygo = 0.d0
!
do ishp = 1, nvqua
xgo = xgo + shpq(ishp)*xoq(1, ishp)
ygo = ygo + shpq(ishp)*xoq(2, ishp)
enddo
!
radig2 = xgo**2 + ygo**2
!
rho0ba =(radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0)
rho0 = rho0ba**(1.d0/(gamlg-1.d0))

htkid = sqrt(1.d0-(tmend)**2)
!
rhoexa = htkid**(-2.d0/(gamlg-1.d0))*rho0
pexa= sentr*(htkid)**(-2.d0*gamlg/(gamlg-1.d0))*(rho0)**gamlg
!
errlp(1) = errlp(1) + (rhoexa -rho)**2*djaco
errlp(2) = errlp(2) + (pres/sentr/rho**gamlg -1.d0)**2*djaco
!errlp(3) = errlp(3) + (pexa -pres)**2*djaco
errlp(3) = errlp(3) + (pexa -pres)**2*djaco
!...Gresho
elseif(ncase.eq.8)then
!
nind = 6.d0
radiv= 0.4d0
radi = sqrt(xg**2+yg**2)
radic = radi
!
ratir = radi/radiv
!
vemag = sqrt(uadv**2 + vadv**2)
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))

if(radi.lt.radiv)then
veexa = 2**(2*nind)*ratir**nind*(1.d0-ratir)**nind
else
veexa = 0.d0
endif
!
call greshoh(xg, yg ,radic, hgre, dhgre, dhgr2)
!
rhoexa =1.d0
pexa = 5.d0 + 2**(4*nind)*hgre
!
errlp(1) = errlp(1) + abs(1.d0/rho-1.d0/rhoexa)**2*djaco
errlp(2) = errlp(2) + abs(vemag-veexa)**2*djaco
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
errlp(4) = errlp(4) + abs(eadv -pexa/(gamlg-1.d0)/rhoexa-0.5d0*veexa**2)**2*djaco
!
elseif(ncase.eq.12)then !...1D isentropic expansion...
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
call solve_isentropic(1000,xg,0.d0,rhoexa,uex,pexa)
!
errlp(1) = errlp(1) + abs(1.d0/rho-1.d0/rhoexa)**2*djaco
errlp(2) = errlp(2) + abs(uex - uadv)**2*djaco
errlp(3) = errlp(3) + abs(eadv -pexa/(gamlg-1.d0)/rhoexa-0.5d0*uex**2)**2*djaco
errlp(4) = errlp(3) + abs(pres -pexa)**2*djaco
!
endif
!
enddo !...(2)ig = 1,ngausd
!
!...Output cell centered values to the file ''
!
if(ncase.eq.7.or.ncase.eq.6.or.ncase.eq.3.or.ncase.eq.4.or.ncase.eq.8.or.ncase.eq.12.or.ncase.eq.13.or.ncase.eq.15)then
!
if(ndens.eq.1)then
rhoc=1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhoc=unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhoc=unkno(1, 1, ielem)
endif
uctr=unkno(1, 2, ielem)
vctr=unkno(1, 3, ielem)
ectr=unkno(1, 4, ielem)
!
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
!r = rc; s =sc
r = 0.d0; s =0.d0
!
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)
!
xc =0.d0
yc = 0.d0
!
do ishp = 1, nvqua
xc = xc + shpq(ishp)*xpq(1, ishp)
yc = yc + shpq(ishp)*xpq(2, ishp)
enddo
!
if(ncase.eq.6.or.ncase.eq.7.or.ncase.eq.3.or.ncase.eq.4.or.ncase.eq.10.or.ncase.eq.11.or.ncase.eq.14.or.ncase.eq.8)then
write(8,'(i8, 10e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,&
pctr/rhoc/(gamlg-1.d0),xc,yc,sqrt(uctr**2+vctr**2),uctr,vctr,geoel(10,ielem)
endif

!...Isentropic smooth flow and Saltzman
if(ncase.eq.12.or.ncase.eq.13.or.ncase.eq.15)then
write(8,'(i8, 7e32.16)')ielem,xc,rhoc,pctr,&
pctr/rhoc/(gamlg-1.d0),sqrt(uctr**2+vctr**2),uctr,vctr
endif
!
endif
!
650 enddo
!
close(8)

!...Show the L2 error on the screen
print*,'L2 error',0.5d0*log10(errlp),sqrt(errlp)
!
end subroutine geterror_lag_curvhybrid_general
!
!...Subroutine to get the minimum Jacobian for every cell
!
subroutine getJacobian_lagcurv(iptri, ipqua, coord, ujacb)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::iptri
integer,  dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,   dimension(1:ncell),                intent(out)::ujacb
!
!...Local integer
!
integer::ie,ig,ishp,ielem
!
!...local integer array
!
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:npqua) :: xpq
real*8,dimension(1:ndimn, 1:nptri) :: xpt
real*8, dimension(1:nptri):: shp, dspr, dsps
real*8, dimension(1:npqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq_pst), posiq(2,ngausdq_pst)
real*8::jtria(ngausd), jquad(ngausdq_pst)
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
real*8::djaco, wi
!
real*8,allocatable::weigh(:), posi(:,:)
!
data eps   / 1.0d-14 /
data c00   / 0.0d0 /
data c10   / 1.0d0 /
data c05   / 0.5d0 /
data c20   / 2.0d0 /
!
!...Give gaussian position and weight...
!
allocate (weigh(ngausd), posi(2,ngausd))
call rutope(2, ngausd, posi, weigh)
call ruqope(2, ngausdq_pst, posiq, weighq)
!
!...Part I: Loop over Triangles
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
!...Points consitituting one element...
!
ielem = ie
ip(1:nvtri) = iptri(1:nvtri, ie)
!
if(ncurv==0)then
xpt(1, 1:3) = coord(1, iptri(1:3,ie))
xpt(2, 1:3) = coord(2, iptri(1:3,ie))
!
xpt(1:2,4) = 0.5d0*(xpt(1:2,1)+xpt(1:2,2))
xpt(1:2,5) = 0.5d0*(xpt(1:2,2)+xpt(1:2,3))
xpt(1:2,6) = 0.5d0*(xpt(1:2,1)+xpt(1:2,3))
elseif(ncurv==1)then
xpt(1, 1:nptri) = coord(1,iptri(1:nptri, ie))
xpt(2, 1:nptri) = coord(2,iptri(1:nptri, ie))
endif

!...Gauss loop
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
wi     = weigh(ig)
!
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
enddo
!
djaco = 0.5d0*wi*(dxdr*dyds - dydr*dxds)
!
jtria(ig) = djaco

enddo !...(2)ig = 1,ngausd

!...Store the minimum Jacobian for triangles
ujacb(ielem) = minval(jtria)
!
550 enddo
!
!...Part II: Quads...
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria

!...Points consitituting one element...
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

if(ncurv==0)then
xpq(1, 1:4) = coord(1, ipq(1:4))
xpq(2, 1:4) = coord(2, ipq(1:4))
!
xpq(1:2,5) = 0.5d0*(xpq(1:2,1)+xpq(1:2,2))
xpq(1:2,6) = 0.5d0*(xpq(1:2,2)+xpq(1:2,3))
xpq(1:2,7) = 0.5d0*(xpq(1:2,3)+xpq(1:2,4))
xpq(1:2,8) = 0.5d0*(xpq(1:2,4)+xpq(1:2,1))
xpq(1:2,9) = 0.5d0*(xpq(1:2,5)+xpq(1:2,7))

elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipq(1:npqua))
xpq(2, 1:npqua) = coord(2,ipq(1:npqua))
endif

!...Gauss loop
do ig = 1,ngausdq_pst !...(2)ig = 1,ngausd
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
jquad(ig) = djaco
!
enddo !...(2)ig = 1,ngausd

!...Store the minimum Jacobian for triangles
ujacb(ielem) = minval(jquad)
!
650 enddo
!
end subroutine getJacobian_lagcurv
!
!...L2 error of linear hybrid meshes in Lagrangian using the deformation gradeint
!
subroutine geterror_lag_hybrid_gd(intfac,iptri, ipqua, coord, coold, geoel, unkno, unkgd)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize),       intent(in)::unkgd
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
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
real*8,dimension(1:ndimn, 1:nvtri) :: xp,xo
real*8,dimension(1:ndimn, 1:nvqua) :: xpq,xoq
real*8,dimension(1:ndegr):: b,bv
real*8:: unknv(1:nq)
real*8, dimension(1:nvtri):: shp, dspr, dsps
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s,  dxdr, dxds, dydr, dyds
real*8::dr,ds,rc,sc,rm, sm,rp,sp
real*8::xg, yg, xc, yc, xgo, ygo
real*8::rho,uadv,vadv,eadv,rhom,rho0ba,rho0
real*8::rhoc, uctr, vctr, ectr, pctr
real*8::pres, pexa, eiadv, eiexa, rhoexa
real*8::djaco, wi
real*8::radie, radii,radie2,radii2, radig2,sentr,rhoin,rhoex,htkid
real*8::htkid2,rhoexa2,pexa2,tmend2
real*8::rhon, rhoini
real*8::rcv, scv
real*8::rexp,veexa,tmend
real*8::uex,uexa,vexa,vemag
real*8::prein,preex
real*8::errorl2(3)
real*8::a11,a12,a21,a22
real*8::rci,sci
real*8::jacbf(2, 2)
real*8::bi(ndegr)
!
real*8::errlp(nq),errl2(nq)!...Error array
!
real*8::weighq(ngausdq_pst), posiq(2,ngausdq_pst)
real*8::weigh(ngausd), posi(2,ngausd)
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
call ruqope(2, ngausdq_pst, posiq, weighq)
!
open(8,file='unradi.dat')
!
errorl2 = 0.d0
errl2 = 0.d0
errlp = 0.d0
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
xo(1, 1:nvtri) = coold(1, ipt(1:nvtri))
xo(2, 1:nvtri) = coold(2, ipt(1:nvtri))

!...Geometry parameters for reference cell...
dr = .5d0
ds = .5d0
rc = geoel(1, ielem)
sc = geoel(2, ielem)

!...The initial cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)

!...Gauss loop
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
xg =0.d0
yg = 0.d0
!
do ishp = 1, nvtri
xg = xg + shp(ishp)*xp(1, ishp)
yg = yg + shp(ishp)*xp(2, ishp)
enddo
!...Deformation gradient transformation matrix

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
djaco = 0.5d0*wi*(a11*a22 - a12*a21)

!
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds
!
unknv = 0.d0
!
do ideg = 1,mdegr
unknv(1:nq) = unknv(1:nq) + unkno(ideg,1:nq,ielem)*b(ideg)
enddo
!
!...Jacobian transformation matrix
!
if(ndens.eq.1)then
rho  = 1.d0/unknv(1)
elseif(ndens.eq.2)then
!
call  getrhoig_tria(rhoini, r, s, xo)
!
call getdensity_triallnl(r, s, xp, xo, rhoini, rhon)
!
rho = rhon
!
elseif(ndens.eq.3)then
rho  = unknv(1)
endif
uadv = unknv(2)
vadv = unknv(3)
eadv = unknv(4)

!
if(ncase.eq.1)then
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
pexa = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 10.d0
!
errl2 = errl2 + (pres -pexa)**2*djaco
elseif(ncase.eq.2)then
!
eiadv = eadv - 0.5d0*(uadv**2 + vadv**2)
eiexa = (1.d0-0.6d0)**(-2.d0*(gamlg-1.d0))
!
!    rhoexa = (1.d0-0.4d0)**(-2.d0)
!
errl2 = errl2 + (eiadv - eiexa)**2*djaco
!    errl2 = errl2 + (rhoexa - rho)**2*djaco
!
elseif(ncase.eq.4)then
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!...some parameters
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
xgo =0.d0
ygo = 0.d0
!
xgo = shp(1)*xo(1, 1) + shp(2)*xo(1, 2) + shp(3)*xo(1, 3)
ygo = shp(1)*xo(2, 1) + shp(2)*xo(2, 2) + shp(3)*xo(2, 3)
!
radig2 = xgo**2 + ygo**2
!
rho0ba =(radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0)
!
!rho = rho0ba**(1.d0/(gamlg-1.d0))
!
htkid = (sqrt(1.d0-(.5d0)**2))**(-2.d0*gamlg/(gamlg-1.d0))
pexa= sentr*(rho)**gamlg*htkid
!
!errl2 = errl2 + (pres -pexa)**2*djaco
errl2 = errl2 + (pres/sentr/rho**gamlg -1.d0)**2*djaco*yg
!
!
endif
!
!...Primitive variables...
!
enddo !...(2)ig = 1,ngausd

!...Output density scatter
if(ndens.eq.1)then
rhoc=1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhoc=unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhoc=unkno(1, 1, ielem)
endif
uctr=unkno(1, 2, ielem)
vctr=unkno(1, 3, ielem)
ectr=unkno(1, 4, ielem)
!
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
r = rc; s =sc

shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
!
xc =0.d0
yc = 0.d0
!
do ishp = 1, nvtri
xc = xc + shp(ishp)*xp(1, ishp)
yc = yc + shp(ishp)*xp(2, ishp)
enddo
!
if(ncase.eq.6.or.ncase.eq.7.or.ncase.eq.3.or.ncase.eq.8) write(8,'(i8, 8e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,&
pctr/rhoc/(gamlg-1.d0),xc,yc,uctr,vctr
!
550 enddo
!
!...For quads...
!
!
!...Loop over elements
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:4) = coord(1, ipq(1:nvqua))
xpq(2, 1:4) = coord(2, ipq(1:nvqua))
!
xoq(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xoq(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
!...Geometry parameters for reference cell...
!
dr = 1.d0
ds = 1.d0
rc = geoel(1, ielem)
sc = geoel(2, ielem)

!...The initial cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)
!
!...Gauss loop
!
do ig = 1,ngausdq_pst !...(2)ig = 1,ngausd
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
dxdr = dxdr + dsprq(ishp)*xoq(1,ishp)
dxds = dxds + dspsq(ishp)*xoq(1,ishp)

dydr = dydr + dsprq(ishp)*xoq(2,ishp)
dyds = dyds + dspsq(ishp)*xoq(2,ishp)
enddo
!
xg = 0.d0
yg = 0.d0
!
do ishp = 1, nvqua
xg = xg + shpq(ishp)*xpq(1, ishp)
yg = yg + shpq(ishp)*xpq(2, ishp)
enddo

!...Deformation gradient transformation matrix

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
djaco = wi*(a11*a22 - a12*a21)
!
b(1) = 1.d0
b(2) = (r-rc)/dr
b(3) = (s-sc)/ds

!DGP2
if(npoly.eq.2)then
b(4) = 0.5d0*b(2)*b(2) - geoel(19, ielem)
b(5) = 0.5d0*b(3)*b(3) - geoel(20, ielem)
b(6) =       b(2)*b(3) - geoel(21, ielem)
endif
!
unknv = 0.d0
!
do ideg = 1,mdegr
unknv(1:nq) = unknv(1:nq) + unkno(ideg,1:nq,ielem)*b(ideg)
enddo
!
!...Jacobian transformation matrix
!
if(ndens.eq.1)then
rho  = 1.d0/unknv(1)
elseif(ndens.eq.2)then
!
call getrhoig_quad(rhoini, r, s, xoq)!
call getdensity_quadllnl(r, s, xpq, xoq, rhoini, rhon)
!
rho = rhon
!
elseif(ndens.eq.3)then
rcv = geoel(5, ielem)
scv = geoel(6, ielem)
!
!
bv(1) = 1.d0
bv(2) = (r-rcv)/dr
bv(3) = (s-scv)/ds
!
unknv(1) = 0.d0
!
do ideg = 1,mdegr
unknv(1) = unknv(1) + unkno(ideg,1,ielem)*bv(ideg)
enddo
!
rho  = unknv(1)
!
endif
uadv = unknv(2)
vadv = unknv(3)
eadv = unknv(4)
!
if(ncase.eq.1)then
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
pexa = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 1.d0
uexa   = sin(pi*xg)*cos(pi*yg)
vexa   =-cos(pi*xg)*sin(pi*yg)
rhoexa = 1.d0
!
vemag = sqrt(uadv**2 + vadv**2)
veexa = sqrt(uexa**2 + vexa**2)
!
errlp(1) = errlp(1) + abs(rho-rhoexa)**2*djaco
errlp(2) = errlp(2) + abs(vemag-veexa)**2*djaco
errlp(3) = errlp(3) + abs(pres - pexa)**2*djaco
!
elseif(ncase.eq.2)then
!
eiadv = eadv - 0.5d0*(uadv**2 + vadv**2)
rhoexa = (1.d0-0.4d0)**(-2.d0)
uex = -xg/(1.d0-0.4d0)
eiexa = (1.d0-0.4d0)**(-2.d0*(gamlg-1.d0))
!
!    rhoexa = (1.d0-0.4d0)**(-2.d0)
!
errl2 = errl2 + (eiadv - eiexa)**2*djaco
!
errlp(1) = errlp(1) + (rho - rhoexa)**2*djaco
errlp(2) = errlp(2) + (uex - uadv)**2*djaco
errlp(3) = errlp(3) + (eadv -eiexa-uadv**2/2.d0-vadv**2/2.d0)**2*djaco
!
elseif(ncase.eq.4)then
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
!...some parameters
!
radie = 1.0d0
radii = 0.9d0
!rhoin = 6.31d-4
rhoex = 1.d-2
!sentr = 2.15d4
prein = 0.1d0
preex = 10.d0
tmend = 0.8d0
tmend2 = 0.5d0
!
rhoin = rhoex*(prein/preex)**(1.d0/gamlg)
sentr = preex/rhoex**gamlg!2.15d4
!
radie2 = radie**2
radii2 = radii**2
!
xgo =0.d0
ygo = 0.d0
!
do ishp = 1, nvqua
xgo = xgo + shpq(ishp)*xoq(1, ishp)
ygo = ygo + shpq(ishp)*xoq(2, ishp)
enddo
!
radig2 = xgo**2 + ygo**2
!
rho0ba =(radie2-radig2)/(radie2-radii2)*rhoin**(gamlg-1.d0) +&
(radig2-radii2)/(radie2-radii2)*rhoex**(gamlg-1.d0)
rho0 = rho0ba**(1.d0/(gamlg-1.d0))

htkid = sqrt(1.d0-(tmend)**2)
htkid2 = sqrt(1.d0-(tmend2)**2)
!
rhoexa = htkid**(-2.d0/(gamlg-1.d0))*rho0
pexa= sentr*(htkid)**(-2.d0*gamlg/(gamlg-1.d0))*(rho0)**gamlg
!
rhoexa2 = htkid2**(-2.d0/(gamlg-1.d0))*rho0
pexa2= sentr*(htkid2)**(-2.d0*gamlg/(gamlg-1.d0))*(rho0)**gamlg
!
errlp(1) = errlp(1) + (rhoexa -rho)**2*djaco*yg
errlp(2) = errlp(2) + (pres/sentr/rho**gamlg -1.d0)**2*djaco*yg!geoel(11, ielem)!yg
errlp(3) = errlp(3) + (pexa -pres)**2*djaco*yg

errl2(1) = errl2(1) + (rhoexa2 -rho)**2*djaco*yg
errl2(2) = errl2(2) + (pres/sentr/rho**gamlg -1.d0)**2*djaco*yg!geoel(11, ielem)!yg
errl2(3) = errl2(3) + (pexa2 -pres)**2*djaco*yg
!
!endif
!
!
elseif(ncase.eq.10)then !...Expansion
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
tmend = 0.4d0
rexp = sqrt(1.d0 + 2.d0*tmend**2)
rhoexa = 1.d0/rexp**3
veexa = 2.d0*tmend/rexp**2*sqrt(xg**2 + yg**2)
pexa = (1.d0-(xg**2 + yg**2)/rexp**2)/rexp**5

if(geoel(10, ielem).gt.100)then
errlp(1) = errlp(1) + (rhoexa-rho)**2*djaco*yg
errlp(2) = errlp(2) + (rhoexa*veexa-rho*sqrt(uadv**2 + vadv**2))**2*djaco*yg
errlp(3) = errlp(3) + (pexa-pres)**2*djaco*yg
endif
!
elseif(ncase.eq.12)then !...1D isentropic expansion...
!
call solve_isentropic(500,xg,0.1d0,rhoexa,uex,pexa)
!
errl2 = errl2 + (rho - rhoex)**2*djaco
errorl2(1) = errorl2(1) + (rho - rhoexa)**2*djaco
errorl2(2) = errorl2(2) + (uex - uadv)**2*djaco
errorl2(3) = errorl2(3) + (eadv -pexa/(gamlg-1.d0)/rhoexa-0.5d0*uex**2)**2*djaco
!
!print*,'Focusing time for Kidder',sqrt((gamlg-1.d0)/sentr/gamlg/2.d0*(radie2-radii2)/(rhoex**(gamlg-1.d0)-rhoin**(gamlg-1.d0)))
!
elseif(ncase.eq.14)then !...Coggeshall problem...
!
pres = (gamlg-1.d0)*rho*(eadv - 0.5d0*(uadv**2 + vadv**2))
!
tmend = 0.5d0
rhoexa = 1.d0/(1.d0-tmend)**2.25d0

!if(geoel(10, ielem).gt.100)then
errlp(1) = errlp(1) + (rhoexa-rho)**2*djaco*geoel(11, ielem)!yg
!errlp(2) = errlp(2) + (rhoexa*veexa-rho*sqrt(uadv**2 + vadv**2))**2*djaco*yg
!errlp(3) = errlp(3) + (pexa-pres)**2*djaco*yg
!endif
!
endif
!
!...Primitive variables...
!
enddo !...(2)ig = 1,ngausd

!...Output density scatter

if(ndens.eq.1)then
rhoc=1.d0/unkno(1, 1, ielem)
elseif(ndens.eq.2)then
rhoc=unkno(1, 1, ielem)
elseif(ndens.eq.3)then
rhoc=unkno(1, 1, ielem)
endif
uctr=unkno(1, 2, ielem)
vctr=unkno(1, 3, ielem)
ectr=unkno(1, 4, ielem)
!
pctr = (gamlg-1.d0)*rhoc*(ectr - 0.5d0*(uctr**2 + vctr**2))
!
r = rc; s =sc

r = 0.d0; s =0.d0
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
xc =0.d0
yc = 0.d0
!
do ishp = 1, nvqua
xc = xc + shpq(ishp)*xpq(1, ishp)
yc = yc + shpq(ishp)*xpq(2, ishp)
enddo
!
if(ncase.eq.6.or.ncase.eq.7.or.ncase.eq.8.or.ncase.eq.3.or.ncase.eq.4.or.ncase.eq.10.or.ncase.eq.11.or.ncase.eq.14)then
write(8,'(i8, 8e32.16)')ielem,sqrt(xc**2 + yc**2),rhoc,pctr,&
pctr/rhoc/(gamlg-1.d0),xc,yc,uctr,vctr
endif

if(ncase.eq.12) then
write(8,'(i8, 7e32.16)')ielem,xc,rhoc,pctr,&
pctr/rhoc/(gamlg-1.d0),sqrt(uctr**2+vctr**2),uctr,vctr
!write(8,'(4e32.16)')xc,rhoexa,uex,pexa
endif

if(ncase.eq.13) write(8,'(i8, 8e32.16)')ielem,sqrt(xc**2),rhoc,pctr,&
pctr/rhoc/(gamlg-1.d0),xc,yc,uctr,vctr
!
650 enddo
!
close(8)
!
print*,'L2 error',ncase, 0.5d0*log10(errlp),sqrt(errlp)
!
end subroutine geterror_lag_hybrid_gd

!
!...Correct the negative cells in Lagrangian
!
subroutine getmesh_lagcrect(intfac, inpoel, ipqua, coord)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(inout)::coord
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:ntria),        intent(in)::inpoel
integer,  dimension(1:nvqua,1:nquad),        intent(in):: ipqua

!...Local integer
integer::ie,ig,ideg,ishp,iv,ielem,nind

!...local integer array
integer,dimension(1:nvtri) :: ip
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nsize) :: indxe

!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndimn, 1:nvtri) :: coorp
real*8, dimension(1:nvtri):: dspr, dsps
real*8, dimension(1:nvqua):: dsprq, dspsq
real*8::weighq(ngausdq_pst), posiq(2,ngausdq_pst)


!...local real
real*8::eps,c00,c05,c10,c20
real*8::r, s, dxdr, dxds, dydr, dyds
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
call ruqope(2, ngausdq_pst, posiq, weighq)
!
indxe = 0
!
!...Loop over elements
!
do 550 ie = 1,ntria !...(1)ie = 1,nelem
!
!...Points consitituting one element...
!
ielem = ie
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
coorp(1, 1:nvtri) = coord(1, ip(1:nvtri))
coorp(2, 1:nvtri) = coord(2, ip(1:nvtri))
!
!...Gauss loop
!
do ig = 1,ngausd !...(2)ig = 1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
wi = weigh(ig)
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
if(djaco.le.-1.d-6)then
print*,'Negative Jacobian for curved trias',djaco,ielem
indxe(ielem) = -1
goto 540
endif
!
enddo !...(2)ig = 1,ngausd
!
540 continue
!
550 enddo
!
!...For quads...
!
!...Loop over elements
do 650 ie = 1,nquad !...(1)ie = 1,nelem
!
ielem = ie + ntria
!
!...Points consitituting one element...
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)!
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))
!
!...Gauss loop
!
do ig = 1,ngausdq_pst !...(2)ig = 1,ngausd
!
r  = posiq(1,ig)
s  = posiq(2,ig)
wi = weighq(ig)
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
if(djaco.le.0.d0)then
print*,'Negative Jacobian for curved cells',djaco,wi,dxdr*dyds - dydr*dxds,ielem
indxe(ielem) = -1
goto 640
endif
!
enddo !...(2)ig = 1,ngausd
!
640 continue
!
650 enddo
!
!...Correct the negative cell.
!
do ie = 1, ntria
ielem = ie
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
if(indxe(ielem).lt.0)then
coord(1:2, ip(4)) = 0.5d0*(coord(1:2, ip(1))+coord(1:2, ip(2)))
coord(1:2, ip(5)) = 0.5d0*(coord(1:2, ip(2))+coord(1:2, ip(3)))
coord(1:2, ip(6)) = 0.5d0*(coord(1:2, ip(3))+coord(1:2, ip(1)))
endif
!
enddo

!...Quads
do ie = 1, nquad
ielem = ie + ntria
ipq(1:nvqua) = ipqua(1:nvqua, ie)
!
if(indxe(ielem).lt.-2)then
coord(1:2,ipq(5)) = 0.5d0*(coord(1:2,ipq(1))+coord(1:2,ipq(2)))
coord(1:2,ipq(6)) = 0.5d0*(coord(1:2,ipq(2))+coord(1:2,ipq(3)))
coord(1:2,ipq(7)) = 0.5d0*(coord(1:2,ipq(3))+coord(1:2,ipq(4)))
coord(1:2,ipq(8)) = 0.5d0*(coord(1:2,ipq(4))+coord(1:2,ipq(1)))
!
coord(1:2, ipq(9)) = -0.25d0*(coord(1:2, ipq(1)) + coord(1:2, ipq(2)) + coord(1:2, ipq(3)) + coord(1:2, ipq(4))) +&
0.5d0*(coord(1:2, ipq(5)) + coord(1:2, ipq(6)) + coord(1:2, ipq(7)) + coord(1:2, ipq(8)))

endif
!
enddo

end subroutine getmesh_lagcrect
!
!
!
subroutine getcellbd(geoel,coord,iptri,ipqua)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
real*8,dimension(1:ngeel,1:nsize),           intent(inout)::geoel
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8, dimension(1:nvqua):: shpq
!
integer,dimension(1:nvqua) :: ipq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
!
integer::ie, ielem,ishp
real*8:: r, s, rc, sc, rp, sp, rm, sm, xc, yc
real*8:: c10
!
c10 = 1.d0
!
geoel(10,:)=0.d0
!
!...Loop over elements
!
do 650 ie = 1,nquad !...(1)ie = 1,nelem
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
rc = geoel(1, ielem)
sc = geoel(2, ielem)
!
r = rc; s =sc
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
xc =0.d0
yc = 0.d0
!
do ishp = 1, nvqua
xc = xc + shpq(ishp)*xpq(1, ishp)
yc = yc + shpq(ishp)*xpq(2, ishp)
enddo
!
!print*,'bad',sqrt(xc**2+yc**2)
!
!if(sqrt(xc**2+yc**2).gt.0.92d0.and.sqrt(xc**2+yc**2).lt.0.98d0)then
!if(atan(yc/xc).lt.19.d0*pi/40.d0.and.atan(yc/xc).gt.11.d0*pi/40.d0)then

if(sqrt(xc**2+yc**2).gt.0.3d0.and.sqrt(xc**2+yc**2).lt.0.7d0)then
if(atan(yc/xc).lt.9.d0*pi/20.d0.and.atan(yc/xc).gt.6.d0*pi/20.d0)then
geoel(10 ,ielem) = 1000
endif
endif

650 enddo

end subroutine getcellbd
!
!...The naming rule for a set of output files(vtu)...
!
subroutine rtname(id,fname)
!
! this sub creates a file name for a restart file
!
implicit none
!
integer, intent(in) :: id
character(len=80), intent(out) :: fname
!
! create the file name based on the id number
!
if(id.lt.   10.and.id.ge.   0) write(fname,1) id
if(id.lt.  100.and.id.ge.  10) write(fname,2) id
if(id.lt. 1000.and.id.ge. 100) write(fname,3) id
if(id.lt.10000.and.id.ge.1000) write(fname,4) id
!
! format statements
!
1 format('output.lag.000',i1,'.vtu')
2 format('output.lag.00' ,i2,'.vtu')
3 format('output.lag.0'  ,i3,'.vtu')
4 format('output.lag.'   ,i4,'.vtu')
!
return
!
end subroutine rtname
!
!...The naming rule for a set of output files(vtk)...
!
subroutine rtnamevtk(id,fname)
!
! this sub creates a file name for a restart file
!
implicit none
!
integer, intent(in) :: id
character(len=80), intent(out) :: fname
!
! create the file name based on the id number
!
if(id.lt.   10.and.id.ge.   0) write(fname,1) id
if(id.lt.  100.and.id.ge.  10) write(fname,2) id
if(id.lt. 1000.and.id.ge. 100) write(fname,3) id
if(id.lt.10000.and.id.ge.1000) write(fname,4) id
!
! format statements
!
1 format('output.lag.000',i1,'.vtk')
2 format('output.lag.00' ,i2,'.vtk')
3 format('output.lag.0'  ,i3,'.vtk')
4 format('output.lag.'   ,i4,'.vtk')
!
return
!
end subroutine rtnamevtk
!
!...Isentropic smooth solution (Sin wave) functions
!
subroutine solve_isentropic(niter,xi,tend,rhoex,uexa,pexa)

integer, intent(in) :: niter
real*8,  intent(in) :: xi, tend

integer :: iter
real*8  :: fx0, fpx0, x0, deltx, j1, j2,c30,c20,c10,c05
real*8  :: uexa, rhoex, pexa
real*8  :: fisen,dfisen
!
c30 = 3.d0
c20 = 2.d0
c10 = 1.d0
c05= 0.5d0
pi =4.d0*atan(1.d0)

!... 1st Riemann invariant
!... Initial guess
x0 = xi

!... Newton loop to find x_0
do iter = 1,niter

!... Function and derivative evaluation
fx0  = fisen(c10,x0,xi,tend)
fpx0 = dfisen(c10,x0,tend)

!... Newton update
deltx = fx0/fpx0
x0 = x0 - deltx

!... Convergence
if (dabs(deltx)<(1.d-12*dabs(x0))) then
exit
end if

enddo !iter

!...Solution 1
j1 = dsqrt(c30) * (c10 + 0.9999995d+00*dsin(pi*x0))
!j1 = dsqrt(c30) * (c10 + 0.100d+00*dsin(c20*pi*x0))

!... 2nd Riemann invariant
!... Initial guess
x0 = xi

!... Newton loop to find x_0
do iter = 1,niter

!... Function and derivative evaluation
fx0  = fisen(-c10,x0,xi,tend)
fpx0 = dfisen(-c10,x0,tend)

!... Newton update
deltx = fx0/fpx0
x0    = x0 - deltx

!... Convergence
if (dabs(deltx)<(1.d-12*dabs(x0))) then
exit
end if

enddo !iter

!...Solution 2
j2 = - dsqrt(c30) * (c10 + 0.9999995d+00*dsin(pi*x0))
!j2 = - dsqrt(c30) * (c10 + 0.100d+00*dsin(2.d0*pi*x0))

!...Primitive variables
uexa  = c05 * (j1 + j2)
rhoex = c05 * (j1 - j2)/(dsqrt(3.d0))
pexa  = rhoex**3

end subroutine solve_isentropic
!
!...The numerator for Newton update
!
function fisen(eig,x0,xi,t)

real*8, intent(in) :: eig, x0, xi, t
real*8 :: sqrt3, fisen

pi =4.d0*atan(1.d0)
sqrt3 = dsqrt(3.d0)
fisen = x0 + eig*(sqrt3 * 0.9999995d0 * dsin(pi*x0) * t) + eig*(sqrt3*t) - xi
!fisentr = x0 + eig*(sqrt3 * 0.100d0 * dsin(2.d0*pi*x0) * t) + eig*(sqrt3*t) - xi

end function fisen
!
!...The denominator for Newton update
!
function dfisen(eig,x0,t)

real*8, intent(in) :: eig, x0, t
real*8 :: sqrt3, dfisen
c10 =1.d0
pi =4.d0*atan(1.d0)
sqrt3 = dsqrt(3.d0)
dfisen = c10 + eig*(sqrt3 * 0.9999995d0 * pi * dcos(pi*x0) * t)
!fpisentr = c10 + eig*(sqrt3 * 0.100d0 * 2.d0*pi * dcos(2.d0*pi*x0) * t)

end function dfisen
!
!....Smooth indicator from machine learning
!
subroutine output_csv(intfac, ipqua, coord, coold, geoel, unkno)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ngeel,1:nsize),     intent(inout)::geoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua

!...Local integer
integer::ie,ig,ideg,ishp,iv,ielem, imt, jmt, iecsv

!...local integer array
integer,dimension(1:nvqua) :: ipq
!...local real array
real*8,dimension(1:ndimn, 1:nvqua) :: xpq
real*8,dimension(1:ndegr):: b, dbdr, dbds
real*8, dimension(1:nq):: unknod,unkp1,unvar1,unvar2
real*8, dimension(ngausdq_pst, ncell):: unvarg
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::weighq(ngausdq_pst), posiq(2,ngausdq_pst)

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

!...Initial
iecsv = 0
!...Give gaussian position and weight...
call ruqope(2, ngausdq_pst, posiq, weighq)
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
do ig = 1,ngausdq_pst !...(2)ig = 1,ngausd
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
!
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

unvarg(ig, ielem) = unknod(1)
!!
enddo !...(2)ig = 1,ngausd
!
geoel(10, ielem) = log10(unvar2(1)/unvar1(1))
!if(ielem.eq.10.or.ielem.eq.11)then
! print*,'ielem',log10(unvar2(1)),log10(unvar1(1))
!print*,'unkno',unkno(:,1,ielem)
!endif
650 enddo

!...Output training data csv fromat.
open(12,file='X.data')
do ie=1,ncell
if(geoel(10, ie).le.-4.5d0.and.geoel(10, ie).ge.-6.5d0)then
write(12,'(f32.16, 24(",", f32.16))')unvarg(:, ie)
endif
enddo
close(12)

!...Y dat
open(12,file='Y.data')
do ie=1,ncell
if(geoel(10, ie).le.-4.5d0.and.geoel(10, ie).ge.-6.5d0)then
write(12,'(f32.16)')geoel(10, ie)
!write(12,'(i8, f32.16)')ie, geoel(10, ie)
endif
enddo
close(12)

!...Recover the geoel
do ie=1,ncell
if(geoel(10, ie).le.(-3.d0*log10(4.d0)-3.d0))then
  geoel(10 ,ie) = 0
else
  geoel(10 ,ie) = 100
endif
enddo


end subroutine output_csv
!
