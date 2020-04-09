!
!...subroutine: RHS for BR2 ...
! 
subroutine getrhsdg(uchar,bface,unkno,rhsel,intfac,inpoel,iptri,ipqua,geofa,geoel,coord,ltelem,fsuel,amatr,itime)
use constant
implicit none
 integer*4,dimension(1:nbfai,nbfac)::bface 
 integer*4,dimension(1:nvtri,1:ntria)::inpoel
 integer*4,dimension(1:nvtri,1:ntria)::iptri
 integer*4,dimension(1:nvqua,1:nquad)::ipqua
 real*8::uchar(1:nq),b(6)
 integer, intent(in)::ltelem(3,ncell), fsuel(3,ncell)
 real*8,dimension(1:ngeel,1:nsize)::geoel
 real*8,dimension(1:ndegr,1:nq,1:nsize)::unkno
 real*8,dimension(1:ndegr,1:nq,1:nsize)::unpri
 real*8,dimension(1:ndegr,1:nq,1:ncell),intent(out)::rhsel
 real*8,dimension(1:ngefa,1:nafac)::geofa
 real*8,dimension(1:ngaus,1:nq,1:nbfac)::unint 
 real*8,dimension(1:ndimn,1:npoin)::coord 
 real*8,intent(in)::amatr(nmatr,ncell)
 integer*4,dimension(1:nifai,1:nafac)::intfac
 integer::itime,ie,ideg
 real*8::dx,dy,p,q
!
!...get the ghost state value...
!
!    call gtghstat(unint,uchar,bface,intfac,geofa,coord)
if(ncurv==0)then
!
  if(nreco == 1)then
!
!     if(nreco.eq.1) call getgraduls_euler(bface, intfac,geoel ,unkno ,coord ,geofa, uchar)

  if(nreco.eq.1) call getgraduls_euler_primitive(bface, intfac,geoel ,unkno ,coord ,geofa, uchar, unpri)
!
!      print*,'good'
!
  if(nlimi.eq.1) call barthlimit_fvm(inpoel, intfac, ltelem, fsuel, geoel, unpri, coord)
!
!...1. face integral for boundary face...
!
  call rhsinvbface_fvm_prim(geofa,geoel,bface,uchar,unpri,intfac,inpoel,rhsel,coord)
!
!...2. face integral for internal face...
!
  call rhsinviface_fvm_prim(geofa,geoel,bface,unpri,intfac,inpoel,rhsel,coord)

else

!
!...1. face integral for boundary face...
!
    call rhsinvbfacedg(geofa,geoel,bface,uchar,unkno,intfac,inpoel,rhsel,coord)
!
!...2. face integral for internal face...
!
    call rhsinvifacedg(geofa,geoel,bface,unkno,intfac,inpoel,rhsel,coord)
!
!...3. domain integral...
!
    if(ntria.gt.0) call rhsinvdomndg_tria(geoel,unkno,intfac,rhsel,coord,inpoel,iptri)
    if(nquad.gt.0) call rhsinvdomndg_quad(geoel,unkno,intfac,rhsel,coord,ipqua)
 endif
!
   elseif(ncurv==1)then
!
 !  if(nlimi.eq.1) call barthlimit(intfac, ltelem, inpoel, unkno ,geoel,coord)
!
!...1. face integral for boundary face...
!

    call rhsinvbfacedg_curv(geofa,geoel,bface,uchar,unkno,intfac,inpoel,rhsel,coord)
!    print*,'11111rhs',rhsel(1:3,2,1)
!
!...2. face integral for internal face...
!
!    print*,'22222'
    call rhsinvifacedg_curv(geofa,geoel,bface,unkno,intfac,inpoel,rhsel,coord)
!    print*,'22222rhs',rhsel(1:3,2,1)
!
!...3. domain integral...
!
!    print*,'33333'
    if(ntria.gt.0)  call rhsinvdomndg_tria_curv(geoel,unkno,intfac,rhsel,coord,iptri)
    if(nquad.gt.0)  call rhsinvdomndg_quad_curv(geoel,unkno,intfac,rhsel,coord,ipqua)

 !  print*,'33333rhs',rhsel(1:3,1,1)
   elseif(ncurv==2)then
!
!...1. face integral for boundary face...
!
!    print*,'11111rhs',unkno(1:4,1,1)
    call rhsinvbfacedg_curv2(geofa,geoel,bface,uchar,unkno,intfac,inpoel,rhsel,coord)
!    print*,'11111rhs',rhsel(1:4,1,1)
!
!...2. face integral for internal face...
!
!    print*,'22222'
    call rhsinvifacedg_curv2(geofa,geoel,bface,unkno,intfac,inpoel,rhsel,coord)
!    print*,'22222rhs',rhsel(1:4,1,1)
!
!...3. domain integral...
!
!    print*,'33333'
 if(ntria.gt.0)  call rhsinvdomndg_tria_curv2(geoel,unkno,intfac,rhsel,coord,inpoel)
!    print*,'33333rhs',rhsel(1:4,1,1)

   endif
!
   return
end subroutine getrhsdg
!
!...subroutine: get the ghost state value...
!
subroutine gtghstat(unint,uchar,bface,intfac,geofa,coord)
 use constant
implicit none
 real*8,intent(inout):: unint(1:ngaus,1:nq,nbfac)
 real*8,dimension(1:ngefa,1:nafac),intent(in)::geofa
 real*8,dimension(1:ndimn,1:npoin),intent(in)::coord 
 integer*4,dimension(1:nbfai,1:nbfac)::bface
 real*8,dimension(1:nq)::uchar
 integer*4::ifa,iel,ier,ig
 integer*4,dimension(1:nifai,1:nafac)::intfac
 real*8::rnx,rny,rnn,rho,uadv1,vadv1,vn
 real*8::rtx,rty,rtt
 real*8::a,b,c,d
 real*8::roinf,uinf,uxinf,uyinf,pinf,roeinf,ainf,uninf,utinf
 real*8::ug,roeg,pg,ugx,ugy,rhog,ugt,ugn
 real*8::p2,u2
 real*8::rhobd,uxbd,uybd,roebd,unbd,utbd,pbd,aspd1,amach1,ubd
 real*8::xg, yg
 real*8::x1,x2,y1,y2,xi ,shp1,shp2
 real*8,allocatable::weigh(:), posi(:,:)
!
 allocate (weigh(ngausf), posi(1,ngausf))

 call rutope(1, ngausf, posi, weigh)

  do 100 ifa=1,nbfac
 
    iel=intfac(1,ifa)
    ier=intfac(2,ifa)
!
    x1 = coord(1,intfac(3,ifa))
    y1 = coord(2,intfac(3,ifa))

    x2 = coord(1,intfac(4,ifa))
    y2 = coord(2,intfac(4,ifa))
!
!...for inviscid...
!
    if(bface(3,ifa).eq.2)then

    elseif(bface(3,ifa).eq.9)then
!
      rnx = geofa(1,ifa)
      rny = geofa(2,ifa)

    if(npoly==0)then
!
!...FVM, npoly=0...
!
      do ig=1,ngausf
       xi =posi(1, ig)
       shp1 = 0.5d0*(1.d0-xi)
       shp2 = 0.5d0*(1.d0+xi)
 
       xg = x1*shp1 + x2*shp2
       yg = y1*shp1 + y2*shp2

      if(nmeth==1)then
       unint(ig,1,ifa) =  (sinh(pi*xg)*sin(pi*yg)+sinh(pi*yg)*sin(pi*xg))/sinh(pi)
       unint(ig,2,ifa) =  (cosh(pi*xg)*sin(pi*yg)+sinh(pi*yg)*cos(pi*xg))*pi/sinh(pi)
       unint(ig,3,ifa) =  (sinh(pi*xg)*cos(pi*yg)+cosh(pi*yg)*sin(pi*xg))*pi/sinh(pi)
!
!...for russian
!        unint(ig,1,ifa) =  2.d0*cos(pi*xg)*sin(2.d0*pi*yg)+2.d0
!        unint(ig,2,ifa) = -2.d0*pi*sin(pi*xg)*sin(2.d0*pi*yg)
!        unint(ig,3,ifa) =  4.d0*pi*cos(pi*xg)*cos(2.d0*pi*yg)
! 
      endif
   
     enddo
     
     elseif(npoly.ge.1)then
!
!...DGM, npoly=1...
!
      do ig=1,ngausf

       xi =posi(1, ig)
       shp1 = 0.5d0*(1.d0-xi)
       shp2 = 0.5d0*(1.d0+xi)
 
       xg = x1*shp1 + x2*shp2
       yg = y1*shp1 + y2*shp2
!
!...for hiro..
        unint(ig,1,ifa) =  (sinh(pi*xg)*sin(pi*yg)+sinh(pi*yg)*sin(pi*xg))/sinh(pi)
        unint(ig,2,ifa) =  (cosh(pi*xg)*sin(pi*yg)+sinh(pi*yg)*cos(pi*xg))*pi/sinh(pi)
        unint(ig,3,ifa) =  (sinh(pi*xg)*cos(pi*yg)+cosh(pi*yg)*sin(pi*xg))*pi/sinh(pi)
!
!...for russian
!        unint(ig,1,ifa) =  2.d0*cos(pi*xg)*sin(2.d0*pi*yg)+2.d0
!        unint(ig,2,ifa) = -2.d0*pi*sin(pi*xg)*sin(2.d0*pi*yg)
!        unint(ig,3,ifa) =  4.d0*pi*cos(pi*xg)*cos(2.d0*pi*yg)
! 
     enddo
     endif

    endif
!...end of BCs
100 enddo 
 end subroutine gtghstat
!
!...subroutine: inviscid flow --face integral for internal face...
!
subroutine rhsinvbfacedg(geofa,geoel,bface,uchar,unkno,intfac,inpoel,rhsel,coord)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nsize)::geoel 
 real*8,intent(in)::uchar(1:nq)
 real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
 real*8,dimension(1:ndegr,1:nq,1:ncell),intent(out)::rhsel
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer::ifa,iel,ier,ie,ifb,k,m,ig,ideg,jdeg
 integer*4,dimension(1:nifai,1:nafac)::intfac
 integer*4,dimension(1:nbfai,nbfac)::bface 
 integer*4,dimension(1:nvtri,1:ntria)::inpoel
 real*8,dimension(1:ngefa,1:nafac)::geofa
 real*8,dimension(1:nq)::unkno1,unkno2
!
!...local array
 real*8::flux(1:nq)
 real*8::bl(6) 
!...local real number
 real*8::xcl,ycl,dxl,dyl,rnx,rny
 real*8::x1,x2,y1,y2,xi,shp1,shp2,xg,yg,wi
 real*8::rho1,rhou1,rhov1,rhoe1,rho2,rhou2,rhov2,rhoe2
 real*8::vn,uadv2,vadv2,dwav1,dwav2
 real*8::rhom1,uadv1,vadv1,pres1
 real*8::vnorm(ndimn),unkno_rie(1:nq,1:2)

 real*8,allocatable::weigh(:), posi(:,:)

 allocate (weigh(ngausf), posi(1,ngausf))
 call rutope(1, ngausf, posi, weigh)
!
!...zero out rhs...
!
   rhsel = 0.d0
!
!...Start loop of boundary faces for RHS....
!
    do 250 ifa=1,nbfac !...(1)ifa=1,nafac
!
     iel=intfac(1,ifa)
     ier=intfac(2,ifa)
   !
     x1 = coord(1,intfac(3,ifa))
     y1 = coord(2,intfac(3,ifa))

     x2 = coord(1,intfac(4,ifa))
     y2 = coord(2,intfac(4,ifa))
   !
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
!
     rnx = geofa(1,ifa)
     rny = geofa(2,ifa) 
!
!...2nd gauss loop...
!
     do ig =1,ngausf !...(4)ig=1,ngausf
!
     xi =posi(1, ig)
 
     shp1 = 0.5d0*(1.d0-xi)
     shp2 = 0.5d0*(1.d0+xi)
 
     xg = x1*shp1 + x2*shp2
     yg = y1*shp1 + y2*shp2

     bl(1) = 1.d0
     bl(2) = (xg-xcl)/dxl
     bl(3) = (yg-ycl)/dyl
     bl(4) = 0.5d0*bl(2)**2-geoel(6,iel)
     bl(5) = 0.5d0*bl(3)**2-geoel(7,iel)
     bl(6) = bl(2)*bl(3)   -geoel(8,iel)
!
!...
!...zero out unkno1, unkno2
      unkno1 = 0.d0
      unkno2 = 0.d0
!
     do ideg = 1,ndegr
       unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
     enddo
!
!...Specify unkno2 for different boundary condition...
    if(bface(3,ifa).eq.2)then
!   
!...Inviscid boundary wall
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!
!...  derived variables
      rhom1 = 1.d0/rho1
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
!
!      rnx = -xg/sqrt(xg**2+yg**2)
!      rny = -yg/sqrt(xg**2+yg**2)
!
!...  unknows for the ghost state
      rho2  = rho1
      rhoe2 = rhoe1
      vn    = uadv1*rnx+ vadv1*rny 
      uadv2 = uadv1 - 2.d0*vn*rnx
      vadv2 = vadv1 - 2.d0*vn*rny
      rhou2 = rho2*uadv2
      rhov2 = rho2*vadv2
!
      unkno2(1) = rho2
      unkno2(2) = rhou2
      unkno2(3) = rhov2
      unkno2(4) = rhoe2
    elseif(bface(3,ifa).eq.4)then
!    
!...Farfield 
      unkno2(1:4) = uchar(1:4)
!
!...Interpolation...
!
    elseif(bface(3,ifa).eq.5)then
!
       unkno2(1:4) = uchar(1:4)
!
!...Extrapolation...
!
    elseif(bface(3,ifa).eq.6)then
!
       unkno2(1:4) = unkno1(1:4)
!
!
     endif
!
!...call Riemann solver...
     dwav1 = rnx
     dwav2 = rny
     vnorm(1) = dwav1
     vnorm(2) = dwav2
!
     unkno_rie(1:nq,1) = unkno1(1:nq) 
     unkno_rie(1:nq,2) = unkno2(1:nq)   
!
     call riemann_2dhllc(unkno_rie, vnorm, flux)
!
!...adding to the RHS...
    do ideg =1 ,ndegr 
       rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
    enddo  
!    
   enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nbfac

end subroutine rhsinvbfacedg
!
!...subroutine: inviscid flow --face integral for internal face...
!
subroutine rhsinvifacedg(geofa,geoel,bface,unkno,intfac,inpoel,rhsel,coord)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nsize)::geoel 
 real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
 real*8,dimension(1:ndegr,1:nq,1:ncell),intent(out)::rhsel
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer::ifa,iel,ier,ie,ifb,k,m,ig,ideg,jdeg
 integer*4,dimension(1:nifai,1:nafac)::intfac
 integer*4,dimension(1:nbfai,nbfac)::bface 
 integer*4,dimension(1:nvtri,1:ntria)::inpoel
 real*8,dimension(1:ngefa,1:nafac)::geofa
 real*8,dimension(1:nq)::unkno1,unkno2
!...local array
 real*8::unkno_rie(1:nq,1:2)
 real*8::vnorm(1:2)
 real*8::flux(1:nq)
 real*8::bl(6),br(6)
!...local real number
 real*8::xcl,ycl,xcr,ycr,dxl,dyl,dxr,dyr,rnx,rny
 real*8::x1,x2,y1,y2,xi,shp1,shp2,xg,yg,wi
 real*8::dwav1,dwav2

 real*8,allocatable::weigh(:), posi(:,:)

 allocate (weigh(ngausf), posi(1,ngausf))
 call rutope(1, ngausf, posi, weigh)
!
!
!
    do 250 ifa=nbfac+1,nafac !...(1)ifa=1,nafac
!
     iel=intfac(1,ifa)
     ier=intfac(2,ifa)
   !
     x1 = coord(1,intfac(3,ifa))
     y1 = coord(2,intfac(3,ifa))

     x2 = coord(1,intfac(4,ifa))
     y2 = coord(2,intfac(4,ifa))
   !
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
   !
     xcr = geoel(1,ier)
     ycr = geoel(2,ier)
     dxr = geoel(4,ier)*0.5d0
     dyr = geoel(5,ier)*0.5d0
!
     rnx = geofa(1,ifa)
     rny = geofa(2,ifa) 
!
!...2nd gauss loop...
!
     do ig =1,ngausf !...(4)ig=1,ngausf
!
     xi =posi(1, ig)
 
     shp1 = 0.5d0*(1.d0-xi)
     shp2 = 0.5d0*(1.d0+xi)
 
     xg = x1*shp1 + x2*shp2
     yg = y1*shp1 + y2*shp2

     bl(1) = 1.d0
     bl(2) = (xg-xcl)/dxl
     bl(3) = (yg-ycl)/dyl
     bl(4) = 0.5d0*bl(2)**2-geoel(6,iel)
     bl(5) = 0.5d0*bl(3)**2-geoel(7,iel)
     bl(6) = bl(2)*bl(3) - geoel(8,iel)
!
     br(1) = 1.d0
     br(2) = (xg-xcr)/dxr
     br(3) = (yg-ycr)/dyr
     br(4) = 0.5d0*br(2)**2-geoel(6,ier)
     br(5) = 0.5d0*br(3)**2-geoel(7,ier)
     br(6) = br(2)*br(3) - geoel(8,ier)
!
!...
!...zero out unkno1, unkno2
      unkno1 = 0.d0
      unkno2 = 0.d0
!
     do ideg = 1,mdegr
       unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
       unkno2(1:nq) = unkno2(1:nq) + unkno(ideg,1:nq,ier)*br(ideg)
     enddo
!
!...call Riemann solver...
!
     dwav1 = rnx
     dwav2 = rny
     vnorm(1) = dwav1
     vnorm(2) = dwav2
!
     unkno_rie(1:nq,1) = unkno1(1:nq) 
     unkno_rie(1:nq,2) = unkno2(1:nq)   
!
     call riemann_2dhllc(unkno_rie, vnorm, flux)
!
!...adding to the RHS...
    do ideg =1 ,ndegr 
       rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
       rhsel(ideg,1:nq,ier)=rhsel(ideg,1:nq,ier) + flux(1:nq)*br(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
    enddo  
!    
   enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nafac

end subroutine rhsinvifacedg
!
!...subroutine: domain integral for BR2...
!
subroutine rhsinvdomndg_tria(geoel,unkno,intfac,rhsel,coord,inpoel,iptri)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nsize)::geoel 
 real*8,dimension(1:ndegr,1:nq,1:nsize)::unkno
 real*8,dimension(1:ndegr,1:nq,1:ncell),intent(inout)::rhsel
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer*4,dimension(1:nvtri,1:ntria)::inpoel
 integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
 real*8::shp1,shp2,shp3
 real*8::volel,wi,c10,p,q,u
 real*8::dbdx(6), dbdy(6), b(6)
 real*8::fluxd(1:ndegr,1:nq)
 real*8::unknod(1:nq)
 real*8::xc,yc,b2,b3,xg,yg,dx,dy
 real*8::xp(3), yp(3)
 real*8::rho,rhou,rhov,rhoe,rhom,uadv,vadv,ua,pres
 real*8::dp,dq
 integer*4::ig,ie,ideg
 real*8,allocatable::weigh(:), posi(:,:)
 integer*4,dimension(1:nifai,1:nafac)::intfac
!
 allocate (weigh(ngausd), posi(2,ngausd))
 call rutope(2, ngausd, posi, weigh)
!
!  print*,'domn',posi
    dbdx = 0.d0
    dbdy = 0.d0

    do 3000 ie=1,ntria !...(1)ie=1,nelem
      xc = geoel(1,ie)
      yc = geoel(2,ie)
!
      dx = geoel(4,ie)*0.5d0
      dy = geoel(5,ie)*0.5d0
!
     dbdx(1)= 0.d0
     dbdx(2)= 1.d0/dx
     dbdx(3)= 0.d0

     dbdy(1)= 0.d0
     dbdy(2)= 0.d0
     dbdy(3)= 1.0/dy
!
     xp(1:3) = coord(1,iptri(1:3,ie))
     yp(1:3) = coord(2,iptri(1:3,ie))
!
     volel=geoel(3,ie) !the area of the elemnt
!
    do ig=1,ngausd !...(2)ngausd
!
     shp1=posi(1,ig)
     shp2=posi(2,ig)
     shp3=1.d0-shp1-shp2
     wi=weigh(ig)*volel 

     xg = xp(1)*shp1 + xp(2)*shp2 + xp(3)*shp3 
     yg = yp(1)*shp1 + yp(2)*shp2 + yp(3)*shp3 

     b(1) = 1.d0
     b(2) = (xg-xc)/dx
     b(3) = (yg-yc)/dy
     b(4) = 0.5d0*b(2)**2-geoel(6,ie)
     b(5) = 0.5d0*b(3)**2-geoel(7,ie)
     b(6) = b(2)*b(3) - geoel(8,ie)
!
     dbdx(4)= (xg-xc)/dx/dx
     dbdx(5)= 0.d0
     dbdx(6)= (yg-yc)/dy/dx         

     dbdy(4)= 0.d0
     dbdy(5)= (yg-yc)/dy/dy
     dbdy(6)= (xg-xc)/dx/dy       
!  
     unknod = 0.d0
!
     do ideg =1,mdegr
        unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ie)*b(ideg)
     enddo
!
     rho  = unknod(1)
     rhou = unknod(2)
     rhov = unknod(3)
     rhoe = unknod(4)
!
     rhom = 1.d0/rho
     uadv = rhou*rhom
     vadv = rhov*rhom
     pres = (gamma-1.d0)*(rhoe-0.5d0*rho*(uadv*uadv+vadv*vadv))   
!
     fluxd(1,1) = rhou*dbdx(1) + rhov*dbdy(1) 
     fluxd(2,1) = rhou*dbdx(2) + rhov*dbdy(2) 
     fluxd(3,1) = rhou*dbdx(3) + rhov*dbdy(3) 
!
     if(npoly.eq.2)then
      fluxd(4,1) = rhou*dbdx(4) + rhov*dbdy(4) 
      fluxd(5,1) = rhou*dbdx(5) + rhov*dbdy(5) 
      fluxd(6,1) = rhou*dbdx(6) + rhov*dbdy(6) 
     endif

     fluxd(1,2) = (uadv* rhou + pres)*dbdx(1) + vadv* rhou*dbdy(1)
     fluxd(2,2) = (uadv* rhou + pres)*dbdx(2) + vadv* rhou*dbdy(2)
     fluxd(3,2) = (uadv* rhou + pres)*dbdx(3) + vadv* rhou*dbdy(3)
     if(npoly.eq.2)then
      fluxd(4,2) = (uadv* rhou + pres)*dbdx(4) + vadv* rhou*dbdy(4)
      fluxd(5,2) = (uadv* rhou + pres)*dbdx(5) + vadv* rhou*dbdy(5)
      fluxd(6,2) = (uadv* rhou + pres)*dbdx(6) + vadv* rhou*dbdy(6)
     endif

     fluxd(1,3) = uadv* rhov*dbdx(1) + (vadv* rhov + pres)*dbdy(1)
     fluxd(2,3) = uadv* rhov*dbdx(2) + (vadv* rhov + pres)*dbdy(2)
     fluxd(3,3) = uadv* rhov*dbdx(3) + (vadv* rhov + pres)*dbdy(3)
     if(npoly.eq.2)then
      fluxd(4,3) = uadv* rhov*dbdx(4) + (vadv* rhov + pres)*dbdy(4)
      fluxd(5,3) = uadv* rhov*dbdx(5) + (vadv* rhov + pres)*dbdy(5)
      fluxd(6,3) = uadv* rhov*dbdx(6) + (vadv* rhov + pres)*dbdy(6)
     endif

     fluxd(1,4) = uadv*(rhoe+pres)*dbdx(1) + vadv*(rhoe+pres)*dbdy(1)
     fluxd(2,4) = uadv*(rhoe+pres)*dbdx(2) + vadv*(rhoe+pres)*dbdy(2)
     fluxd(3,4) = uadv*(rhoe+pres)*dbdx(3) + vadv*(rhoe+pres)*dbdy(3)
     if(npoly.eq.2)then
      fluxd(4,4) = uadv*(rhoe+pres)*dbdx(4) + vadv*(rhoe+pres)*dbdy(4)
      fluxd(5,4) = uadv*(rhoe+pres)*dbdx(5) + vadv*(rhoe+pres)*dbdy(5)
      fluxd(6,4) = uadv*(rhoe+pres)*dbdx(6) + vadv*(rhoe+pres)*dbdy(6)
     endif
!
!
!finally, scatter the contribution to the RHS
!    
    do ideg = 1,ndegr
     rhsel(ideg,1:nq,ie)=rhsel(ideg,1:nq,ie)+fluxd(ideg,1:nq)*wi
    enddo
!
   enddo !...(2)ngausd

3000 enddo !...(1)ie=1,nelem

    return
end subroutine rhsinvdomndg_tria
!
!...subroutine: domain integral for linear quad..
!
subroutine rhsinvdomndg_quad(geoel,unkno,intfac,rhsel,coord,ipqua)
use constant
implicit none
real*8,dimension(1:ngeel,1:nsize),intent(in)::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndegr,1:nq,1:ncell),intent(inout)::rhsel
real*8,dimension(1:ndimn,1:npoin)::coord
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer*4,dimension(1:nifai,1:nafac)::intfac
!
!...local array
real*8::volel,wi
real*8:: r, s, djaco
real*8:: rm, sm, rp, sp
real*8:: dxdr,dxds,dydr,dyds
real*8:: eps,c00,c10,c05,c20
real*8::xc,yc,xg,yg,dx,dy
real*8::dbdx(6), dbdy(6), b(6)
real*8::fluxd(1:ndegr,1:nq)
real*8::unknod(1:nq)
real*8,dimension(1:2, 1:4)::xpq
real*8,dimension(1:4)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8::rho,rhou,rhov,rhoe,rhom,uadv,vadv,pres
integer*4::ig,ideg,ielem,ishp
integer::ie
!
data eps / 1.0d-06 /
data c00 / 0.0d0 /
data c10 / 1.0d0 /
data c05 / 0.5d0 /
data c20 / 2.0d0 /
!
call ruqope(2, ngausdq, posiq, weighq)
!
!  print*,'domn',posi
dbdx = 0.d0
dbdy = 0.d0

do 3000 ie = 1,nquad !...(1)ie=1,nelem
!
!...global numbering
!
ielem = ie + ntria
!
xpq(1, 1:4) = coord(1, ipqua(1:4,ie))
xpq(2, 1:4) = coord(2, ipqua(1:4,ie))
!...
xc    = geoel(1, ielem)
yc    = geoel(2, ielem)
volel = geoel(3, ielem)
dx    = geoel(4, ielem)*0.5d0
dy    = geoel(5, ielem)*0.5d0
!
dbdx(1)= 0.d0
dbdx(2)= 1.d0/dx
dbdx(3)= 0.d0

dbdy(1)= 0.d0
dbdy(2)= 0.d0
dbdy(3)= 1.0/dy
!
do ig=1,ngausdq !...(2)ngausd
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
xg = 0.d0
yg = 0.d0
!
do ishp = 1, 4
xg = xg + shpq(ishp)*xpq(1,ishp)
yg = yg + shpq(ishp)*xpq(2,ishp)
enddo
!
b(1) = 1.d0
b(2) = (xg-xc)/dx
b(3) = (yg-yc)/dy
b(4) = 0.5d0*b(2)**2-geoel(6,ielem)
b(5) = 0.5d0*b(3)**2-geoel(7,ielem)
b(6) = b(2)*b(3)   - geoel(8,ielem)
!
dbdx(4)= (xg-xc)/dx/dx
dbdx(5)= 0.d0
dbdx(6)= (yg-yc)/dy/dx

dbdy(4)= 0.d0
dbdy(5)= (yg-yc)/dy/dy
dbdy(6)= (xg-xc)/dx/dy
!
unknod = 0.d0
!
do ideg =1,mdegr
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo
!
rho  = unknod(1)
rhou = unknod(2)
rhov = unknod(3)
rhoe = unknod(4)
!
rhom = 1.d0/rho
uadv = rhou*rhom
vadv = rhov*rhom
pres = (gamma-1.d0)*(rhoe-0.5d0*rho*(uadv*uadv+vadv*vadv))
!
fluxd(1,1) = rhou*dbdx(1) + rhov*dbdy(1)
fluxd(2,1) = rhou*dbdx(2) + rhov*dbdy(2)
fluxd(3,1) = rhou*dbdx(3) + rhov*dbdy(3)
!
if(npoly.eq.2)then
fluxd(4,1) = rhou*dbdx(4) + rhov*dbdy(4)
fluxd(5,1) = rhou*dbdx(5) + rhov*dbdy(5)
fluxd(6,1) = rhou*dbdx(6) + rhov*dbdy(6)
endif

fluxd(1,2) = (uadv* rhou + pres)*dbdx(1) + vadv* rhou*dbdy(1)
fluxd(2,2) = (uadv* rhou + pres)*dbdx(2) + vadv* rhou*dbdy(2)
fluxd(3,2) = (uadv* rhou + pres)*dbdx(3) + vadv* rhou*dbdy(3)
if(npoly.eq.2)then
fluxd(4,2) = (uadv* rhou + pres)*dbdx(4) + vadv* rhou*dbdy(4)
fluxd(5,2) = (uadv* rhou + pres)*dbdx(5) + vadv* rhou*dbdy(5)
fluxd(6,2) = (uadv* rhou + pres)*dbdx(6) + vadv* rhou*dbdy(6)
endif

fluxd(1,3) = uadv* rhov*dbdx(1) + (vadv* rhov + pres)*dbdy(1)
fluxd(2,3) = uadv* rhov*dbdx(2) + (vadv* rhov + pres)*dbdy(2)
fluxd(3,3) = uadv* rhov*dbdx(3) + (vadv* rhov + pres)*dbdy(3)
if(npoly.eq.2)then
fluxd(4,3) = uadv* rhov*dbdx(4) + (vadv* rhov + pres)*dbdy(4)
fluxd(5,3) = uadv* rhov*dbdx(5) + (vadv* rhov + pres)*dbdy(5)
fluxd(6,3) = uadv* rhov*dbdx(6) + (vadv* rhov + pres)*dbdy(6)
endif

fluxd(1,4) = uadv*(rhoe+pres)*dbdx(1) + vadv*(rhoe+pres)*dbdy(1)
fluxd(2,4) = uadv*(rhoe+pres)*dbdx(2) + vadv*(rhoe+pres)*dbdy(2)
fluxd(3,4) = uadv*(rhoe+pres)*dbdx(3) + vadv*(rhoe+pres)*dbdy(3)
if(npoly.eq.2)then
fluxd(4,4) = uadv*(rhoe+pres)*dbdx(4) + vadv*(rhoe+pres)*dbdy(4)
fluxd(5,4) = uadv*(rhoe+pres)*dbdx(5) + vadv*(rhoe+pres)*dbdy(5)
fluxd(6,4) = uadv*(rhoe+pres)*dbdx(6) + vadv*(rhoe+pres)*dbdy(6)
endif
!
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem)+fluxd(ideg,1:nq)*djaco
enddo
!
enddo !...(2)ngausd

3000 enddo !...(1)ie=1,nelem

return
end subroutine rhsinvdomndg_quad

!
!...subroutine: inviscid flow --face integral for internal curved face...
!
subroutine rhsinvbfacedg_curv(geofa,geoel,bface,uchar,unkno,intfac,inpoel,rhsel,coord)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nsize),intent(in)::geoel
 real*8,intent(in)::uchar(1:nq)
 real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
 real*8,dimension(1:ndegr,1:nq,1:ncell),intent(inout)::rhsel
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer*4,dimension(1:nifai,1:nafac)::intfac
 integer*4,dimension(1:nbfai,nbfac)::bface 
 integer*4,dimension(1:nvtri,1:ntria)::inpoel
 real*8,dimension(1:ngefa,1:nafac)::geofa
!
!...local array
 real*8,dimension(1:nq)::unkno1,unkno2
 real*8::xpin(2, 2)
 real*8::xp(2, nptfa)
 real*8,dimension(1:nptfa)::shp, dshpr
 real*8::weigh(ngausf), posi(1,ngausf)
 real*8::flux(1:nq)
 real*8::bl(6) 
!...local real number
 integer::ifa,iel,ier,ie,ifb,k,m,ig,ishp,ideg,jdeg
 real*8::xcl,ycl,dxl,dyl,rnx,rny
 real*8:: r, djaco, dxdr, dydr 
 real*8::xg,yg,wi
 real*8::rho1,rhou1,rhov1,rhoe1,rho2,rhou2,rhov2,rhoe2
 real*8::vn,dwav1,dwav2
 real*8::rhom1,uadv1,vadv1,pres1,temp1,qadv1
 real*8::rhom2,uadv2,vadv2,pres2,temp2,qadv2
!
 real*8::cbdry,entrb,vbdry
 real*8::csinf,pinf,roinf,uinf,vinf,vninf
 real*8::csou1,mach1,vdon1
 real*8::rinv1,rinv2
 real*8::vnorm(ndimn),unkno_rie(1:nq,1:2)
!
 call rutope(1, ngausf, posi, weigh)
!
! print*,'gauss',posi(1,ngausf)
!
!...zero out rhs...
!
   rhsel = 0.d0
!
    do 250 ifa=1,nbfac !...(1)ifa=1,nafac
!
     iel=intfac(1,ifa)
     ier=intfac(2,ifa)
!
   if(ncurv==0)then
!
     xpin(1, 1:2) = coord(1, intfac(3:4, ifa)) !...One face constitutes of 'nptfa' points...
     xpin(2, 1:2) = coord(2, intfac(3:4, ifa)) 
!
     xp(1, 1:2) = xpin(1, 1:2)
     xp(2, 1:2) = xpin(2, 1:2)
     xp(1, 3) = 0.5d0*(xpin(1, 1) + xpin(1, 2))
     xp(2, 3) = 0.5d0*(xpin(2, 1) + xpin(2, 2))
   elseif(ncurv==1)then
!
     xp(1, 1:nptfa) = coord(1, intfac(3:(2+nptfa), ifa))  
     xp(2, 1:nptfa) = coord(2, intfac(3:(2+nptfa), ifa)) 
   endif
   !
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
!
!...2nd gauss loop...
!
     do ig =1, ngausf !...(4)ig=1,ngausf
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
    do ishp = 1, nptfa
     dxdr = dxdr + dshpr(ishp)*xp(1, ishp)
     dydr = dydr + dshpr(ishp)*xp(2, ishp)
    enddo 
!
    djaco = sqrt(dxdr**2 + dydr**2)
!
    dwav1 = dydr/djaco
    dwav2 =-dxdr/djaco
!
    rnx = dwav1
    rny = dwav2
!
     xg = 0.d0
     yg = 0.d0
!
   do ishp = 1, nptfa
     xg = xg + shp(ishp)*xp(1, ishp) 
     yg = yg + shp(ishp)*xp(2, ishp) 
    enddo
!
!    print*,'xg,yg',rnx,rny,ifa
!
     bl(1) = 1.d0
     bl(2) = (xg-xcl)/dxl
     bl(3) = (yg-ycl)/dyl
     bl(4) = 0.5d0*bl(2)**2-geoel(6,iel)
     bl(5) = 0.5d0*bl(3)**2-geoel(7,iel)
     bl(6) = bl(2)*bl(3)   -geoel(8,iel)
!
!...zero out unkno1, unkno2
      unkno1 = 0.d0
      unkno2 = 0.d0
!
     do ideg = 1,ndegr
       unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
     enddo

 !   print*,'unkno1', unkno(1:3, 2, iel),iel, dxl,dyl
!
!...Specify unkno2 for different boundary condition...
    if(bface(3,ifa).eq.2)then
!   print*,'sym ifa',ifa,intfac(3:4,ifa)
!...Inviscid boundary wall
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!
!...  derived variables
      rhom1 = 1.d0/rho1
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
!
!      rnx = -xg/sqrt(xg**2+yg**2)
!      rny = -yg/sqrt(xg**2+yg**2)
!
!...  unknows for the ghost state
      rho2  = rho1
      rhoe2 = rhoe1
      vn    = uadv1*rnx+ vadv1*rny 
      uadv2 = uadv1 - 2.d0*vn*rnx
      vadv2 = vadv1 - 2.d0*vn*rny
      rhou2 = rho2*uadv2
      rhov2 = rho2*vadv2
!
      unkno2(1) = rho2
      unkno2(2) = rhou2
      unkno2(3) = rhov2
      unkno2(4) = rhoe2
!
!...Variation for invisicd wall...
!
    elseif(bface(3,ifa).eq.8)then
!   print*,'sym ifa',ifa,intfac(3:4,ifa)
!...Inviscid boundary wall
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!
!...  derived variables
      rhom1 = 1.d0/rho1
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
      qadv1 = uadv1**2 + vadv1**2
      pres1 = max(1.e-6,(gamma-1.d0)*(rhoe1-0.5d0*rho1*qadv1))
! 
      flux(1)  = 0.d0
      flux(2)  = pres1*dwav1
      flux(3)  = pres1*dwav2
      flux(4)  = 0.d0
!
!...Riemann solver for farfield...
!
    elseif(bface(3,ifa).eq.4)then
!   print*,'far',ifa,intfac(3:4,ifa)
!...Farfield 
      unkno2(1:4) = uchar(1:4)
!
!...Riemann invariant for farfield...
!
  elseif(bface(3,ifa).eq.7)then
!   print*,'sym ifa',ifa,intfac(3:4,ifa)
!...Inviscid boundary wall
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!
!...  derived variables
      rhom1 = 1.d0/rho1
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
      qadv1 = uadv1**2 + vadv1**2
      pres1 = max(1.e-6,(gamma-1.d0)*(rhoe1-0.5d0*rho1*qadv1))
      csou1 = sqrt( max(1.e-6,gamma*pres1*rhom1) )
      vdon1 = uadv1*rnx+ vadv1*rny 
! 
      mach1 = vdon1/csou1
!
!...Farfield variables...
!
      roinf = uchar(1)
      uinf  = uchar(2)/roinf
      vinf  = uchar(3)/roinf
      pinf  = max(1.e-6,(gamma-1.d0)*(uchar(4)-0.5d0*roinf*(uinf**2 + vinf**2)))
      csinf = sqrt( max(1.e-6,gamma*pinf/roinf) )
      vninf = uinf*rnx + vinf*rny

!...Subsonic inflow...
     if(mach1.gt.-1.d0.and.mach1.lt.0.d0)then
!
!...Riemann invariant...
!
       rinv1 = vdon1 + 2.d0*csou1/(gamma-1.d0)
       rinv2 = vninf - 2.d0*csinf/(gamma-1.d0)
!
!...vbdry:velocity at the boundary...
!
       vbdry = 0.5d0*(rinv1 + rinv2)
       cbdry = 0.25d0*(gamma-1.d0)*(rinv2-rinv1)
!
!...entrb: entropy at boundary
!
       entrb = csinf**2/gamma/roinf**(gamma-1.d0)
!
       rho2  = (cbdry**2/gamma/entrb)**(1.d0/(gamma-1.d0))
       uadv2 = uinf + (vbdry-vninf)*rnx
       vadv2 = vinf + (vbdry-vninf)*rny
       pres2 = rho2*cbdry**2/gamma
       rhou2 = rho2*uadv2
       rhov2 = rho2*vadv2
       rhoe2 = pres2/(gamma-1.d0) + 0.5d0*rho2*(uadv2**2 + vadv2**2)
!
      flux(1)  = vbdry* rho2
      flux(2)  = vbdry* rhou2 + pres2*dwav1
      flux(3)  = vbdry* rhov2 + pres2*dwav2
      flux(4)  = vbdry*(rhoe2 + pres2)
!
!...Subsonic outflow...
!
     elseif(mach1.gt.0.d0.and.mach1.lt.1.d0)then
!
!...Riemann invariant...
!
       rinv1 = vdon1 + 2.d0*csou1/(gamma-1.d0)
       rinv2 = vninf - 2.d0*csinf/(gamma-1.d0)
!
!...vbdry:velocity at the boundary...
!
       vbdry = 0.5d0*(rinv1 + rinv2)
       cbdry = 0.25d0*(gamma-1.d0)*(rinv2-rinv1)
!
!...entrb: entropy at boundary
!
       entrb = csou1**2/gamma/rho1**(gamma-1.d0)
!
       rho2  = (cbdry**2/gamma/entrb)**(1.d0/(gamma-1.d0))
       uadv2 = uadv1 + (vbdry-vdon1)*rnx
       vadv2 = vadv1 + (vbdry-vdon1)*rny
       pres2 = rho2*cbdry**2/gamma
       rhou2 = rho2*uadv2
       rhov2 = rho2*vadv2
       rhoe2 = pres2/(gamma-1.d0) + 0.5d0*rho2*(uadv2**2 + vadv2**2)
!
      flux(1)  = vbdry* rho2
      flux(2)  = vbdry* rhou2 + pres2*dwav1
      flux(3)  = vbdry* rhov2 + pres2*dwav2
      flux(4)  = vbdry*(rhoe2 + pres2)
!
      endif
!
!...Isothermal inviscid wall...
    elseif(bface(3,ifa).eq.3)then
!
!   print*,'sym ifa',ifa,intfac(3:4,ifa)
!...Inviscid boundary wall
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!
!...  derived variables
      rhom1 = 1.d0/rho1
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
!
      qadv1 = uadv1**2 + vadv1**2
      pres1 = max(1.e-6,(gamma-1.d0)*(rhoe1-0.5d0*rho1*qadv1))
      temp1 = pres1*rhom1*gamma
!
      pres2 = pres1
      temp2 = 2.d0-temp1   
!
!...  unknows for the ghost state
      rho2  = pres2/temp2*gamma
      vn    = uadv1*rnx+ vadv1*rny 
      uadv2 = uadv1 - 2.d0*vn*rnx
      vadv2 = vadv1 - 2.d0*vn*rny
      qadv2 = uadv2**2 + vadv2**2      
      rhou2 = rho2*uadv2
      rhov2 = rho2*vadv2
      rhoe2 = pres2/(gamma-1.d0) + 0.5d0*rho2*qadv2
!
      unkno2(1) = rho2
      unkno2(2) = rhou2
      unkno2(3) = rhov2
      unkno2(4) = rhoe2
!
!...Inflow...
!
    elseif(bface(3,ifa).eq.5)then
!
       unkno2(1:4) = uchar(1:4)
!
!...Outflow...
!
    elseif(bface(3,ifa).eq.6)then
!
       unkno2 = unkno1
!
     endif
!
!...call Riemann solver...
!     dwav1 = rnx
!     dwav2 = rny
     vnorm(1) = dwav1
     vnorm(2) = dwav2
!
     unkno_rie(1:nq,1) = unkno1(1:nq) 
     unkno_rie(1:nq,2) = unkno2(1:nq)   
!
    if(bface(3,ifa).ne.7.and.bface(3,ifa).ne.8)then
!     if(bface(3,ifa).eq.2)then
!      call riemann_2dlf(unkno_rie, vnorm, flux)
!     else
      call riemann_2dhllc(unkno_rie, vnorm, flux)
!     endif
    endif
!
!...adding to the RHS...
    do ideg =1 ,ndegr 
       rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*djaco*wi
    enddo  
!    
   enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nbfac

end subroutine rhsinvbfacedg_curv
!
!...subroutine: inviscid flow --face integral for curved internal face...
!
subroutine rhsinvifacedg_curv(geofa,geoel,bface,unkno,intfac,inpoel,rhsel,coord)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nsize)::geoel
 real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
 real*8,dimension(1:ndegr,1:nq,1:ncell),intent(inout)::rhsel
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer::ifa,iel,ier,ie,ifb,k,m,ig,ideg,jdeg
 integer*4,dimension(1:nifai,1:nafac)::intfac
 integer*4,dimension(1:nbfai,nbfac)::bface 
 integer*4,dimension(1:nvtri,1:ntria)::inpoel
 real*8,dimension(1:ngefa,1:nafac)::geofa
!
!...local array
 real*8,dimension(1:nq)::unkno1,unkno2
 real*8::unkno_rie(1:nq,1:2)
 real*8::xpin(2, 2)
 real*8::xp(2, nptfa)
 real*8,dimension(1:nptfa)::shp, dshpr
 real*8::weigh(ngausf), posi(1,ngausf)
 real*8::vnorm(1:2)
 real*8::flux(1:nq)
 real*8::bl(6),br(6)
!
!...local real number
 real*8::xcl,ycl,xcr,ycr,dxl,dyl,dxr,dyr,dxdr,dydr
 real*8:: r, djaco
 real*8::xg,yg,wi
 real*8::dwav1,dwav2
!
 integer ishp
!
 call rutope(1, ngausf, posi, weigh)
!
!
    do 250 ifa=nbfac+1,nafac !...(1)ifa=1,nafac
!
     iel=intfac(1,ifa)
     ier=intfac(2,ifa)
!
!    print*,'ifa',ifa,iel,ier
   !
  if(ncurv==0)then
!
     xpin(1, 1:2) = coord(1, intfac(3:4, ifa)) !...One face constitutes of 'nptfa' points...
     xpin(2, 1:2) = coord(2, intfac(3:4, ifa)) 
!
     xp(1, 1:2) = xpin(1, 1:2)
     xp(2, 1:2) = xpin(2, 1:2)
     xp(1, 3) = 0.5d0*(xpin(1, 1) + xpin(1, 2))
     xp(2, 3) = 0.5d0*(xpin(2, 1) + xpin(2, 2))
   elseif(ncurv==1)then
     xp(1, 1:nptfa) = coord(1, intfac(3:(2+nptfa), ifa))  
     xp(2, 1:nptfa) = coord(2, intfac(3:(2+nptfa), ifa)) 
   endif
   !
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
   !
     xcr = geoel(1,ier)
     ycr = geoel(2,ier)
     dxr = geoel(4,ier)*0.5d0
     dyr = geoel(5,ier)*0.5d0
!
!...2nd gauss loop...
!
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
    do ishp = 1, nptfa
     dxdr = dxdr + dshpr(ishp)*xp(1, ishp)
     dydr = dydr + dshpr(ishp)*xp(2, ishp)
    enddo 
!
    djaco = sqrt(dxdr**2 + dydr**2)
!
    dwav1 = dydr/djaco
    dwav2 =-dxdr/djaco
!
     xg = 0.d0
     yg = 0.d0
!
   do ishp = 1, nptfa
     xg = xg + shp(ishp)*xp(1, ishp) 
     yg = yg + shp(ishp)*xp(2, ishp) 
    enddo

     bl(1) = 1.d0
     bl(2) = (xg-xcl)/dxl
     bl(3) = (yg-ycl)/dyl
     bl(4) = 0.5d0*bl(2)**2-geoel(6,iel)
     bl(5) = 0.5d0*bl(3)**2-geoel(7,iel)
     bl(6) = bl(2)*bl(3)   -geoel(8,iel)
!
     br(1) = 1.d0
     br(2) = (xg-xcr)/dxr
     br(3) = (yg-ycr)/dyr
     br(4) = 0.5d0*br(2)**2-geoel(6,ier)
     br(5) = 0.5d0*br(3)**2-geoel(7,ier)
     br(6) = br(2)*br(3)   -geoel(8,ier)
!
!...
!...zero out unkno1, unkno2
      unkno1 = 0.d0
      unkno2 = 0.d0
!
     do ideg = 1,mdegr
       unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
       unkno2(1:nq) = unkno2(1:nq) + unkno(ideg,1:nq,ier)*br(ideg)
     enddo
!
!...call Riemann solver...
!
     vnorm(1) = dwav1
     vnorm(2) = dwav2
!
     unkno_rie(1:nq,1) = unkno1(1:nq) 
     unkno_rie(1:nq,2) = unkno2(1:nq)   
!
     call riemann_2dhllc(unkno_rie, vnorm, flux)
!
!...adding to the RHS...
    do ideg =1 ,ndegr 
!       rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
!       rhsel(ideg,1:nq,ier)=rhsel(ideg,1:nq,ier) + flux(1:nq)*br(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
    enddo
!
    do ideg =1 ,ndegr 
       rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*wi*djaco
       rhsel(ideg,1:nq,ier)=rhsel(ideg,1:nq,ier) + flux(1:nq)*br(ideg)*wi*djaco
    enddo    
!    
   enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nafac

end subroutine rhsinvifacedg_curv
!
!...subroutine: domain integral for curved triangle..
!
subroutine rhsinvdomndg_tria_curv(geoel,unkno,intfac,rhsel,coord,iptri)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nsize),intent(in)::geoel
 real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
 real*8,dimension(1:ndegr,1:nq,1:ncell),intent(inout)::rhsel
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
 integer*4,dimension(1:nifai,1:nafac)::intfac
!
!...local array
 real*8::volel,wi
 real*8:: r, s, djaco
 real*8:: dxdr,dxds,dydr,dyds 
 real*8:: eps,c00,c10,c05,c20
 real*8::xc,yc,xg,yg,dx,dy
 real*8::dbdx(6), dbdy(6), b(6)
 real*8::fluxd(1:ndegr,1:nq)
 real*8::unknod(1:nq)
 real*8,dimension(1:2, 1:nptri)::xp 
 real*8,dimension(1:nptri)::shp, dspr, dsps
 real*8:: weigh(ngausd), posi(2, ngausd)
 real*8::rho,rhou,rhov,rhoe,rhom,uadv,vadv,pres
 integer*4::ig,ideg,ielem,ishp
 integer:: ie
!
  data eps / 1.0d-06 / 
  data c00 / 0.0d0 / 
  data c10 / 1.0d0 / 
  data c05 / 0.5d0 / 
  data c20 / 2.0d0 / 
!
 call rutope(2, ngausd, posi, weigh)
!
!  print*,'domn',posi
    dbdx = 0.d0
    dbdy = 0.d0

   do 3000 ie =1,ntria!...(1)ie=1,nelem
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
!...
    xc    = geoel(1, ielem)
    yc    = geoel(2, ielem)
    volel = geoel(3, ielem)
    dx    = geoel(4, ielem)*0.5d0
    dy    = geoel(5, ielem)*0.5d0
!
    dbdx(1)= 0.d0
    dbdx(2)= 1.d0/dx
    dbdx(3)= 0.d0

    dbdy(1)= 0.d0
    dbdy(2)= 0.d0
    dbdy(3)= 1.0/dy
!
    do ig=1,ngausd !...(2)ngausd
!
     r  = posi(1,ig)
     s  = posi(2,ig)
    wi  = weigh(ig)
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
     b(6) = b(2)*b(3)   - geoel(8,ielem)
!
     dbdx(4)= (xg-xc)/dx/dx
     dbdx(5)= 0.d0
     dbdx(6)= (yg-yc)/dy/dx         

     dbdy(4)= 0.d0
     dbdy(5)= (yg-yc)/dy/dy
     dbdy(6)= (xg-xc)/dx/dy       
!  
     unknod = 0.d0
!
     do ideg =1,mdegr
        unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
     enddo
!
     rho  = unknod(1)
     rhou = unknod(2)
     rhov = unknod(3)
     rhoe = unknod(4)
!
     rhom = 1.d0/rho
     uadv = rhou*rhom
     vadv = rhov*rhom
     pres = (gamma-1.d0)*(rhoe-0.5d0*rho*(uadv*uadv+vadv*vadv))   
!
     fluxd(1,1) = rhou*dbdx(1) + rhov*dbdy(1) 
     fluxd(2,1) = rhou*dbdx(2) + rhov*dbdy(2) 
     fluxd(3,1) = rhou*dbdx(3) + rhov*dbdy(3) 
!
     if(npoly.eq.2)then
      fluxd(4,1) = rhou*dbdx(4) + rhov*dbdy(4) 
      fluxd(5,1) = rhou*dbdx(5) + rhov*dbdy(5) 
      fluxd(6,1) = rhou*dbdx(6) + rhov*dbdy(6) 
     endif

     fluxd(1,2) = (uadv* rhou + pres)*dbdx(1) + vadv* rhou*dbdy(1)
     fluxd(2,2) = (uadv* rhou + pres)*dbdx(2) + vadv* rhou*dbdy(2)
     fluxd(3,2) = (uadv* rhou + pres)*dbdx(3) + vadv* rhou*dbdy(3)
     if(npoly.eq.2)then
      fluxd(4,2) = (uadv* rhou + pres)*dbdx(4) + vadv* rhou*dbdy(4)
      fluxd(5,2) = (uadv* rhou + pres)*dbdx(5) + vadv* rhou*dbdy(5)
      fluxd(6,2) = (uadv* rhou + pres)*dbdx(6) + vadv* rhou*dbdy(6)
     endif

     fluxd(1,3) = uadv* rhov*dbdx(1) + (vadv* rhov + pres)*dbdy(1)
     fluxd(2,3) = uadv* rhov*dbdx(2) + (vadv* rhov + pres)*dbdy(2)
     fluxd(3,3) = uadv* rhov*dbdx(3) + (vadv* rhov + pres)*dbdy(3)
     if(npoly.eq.2)then
      fluxd(4,3) = uadv* rhov*dbdx(4) + (vadv* rhov + pres)*dbdy(4)
      fluxd(5,3) = uadv* rhov*dbdx(5) + (vadv* rhov + pres)*dbdy(5)
      fluxd(6,3) = uadv* rhov*dbdx(6) + (vadv* rhov + pres)*dbdy(6)
     endif

     fluxd(1,4) = uadv*(rhoe+pres)*dbdx(1) + vadv*(rhoe+pres)*dbdy(1)
     fluxd(2,4) = uadv*(rhoe+pres)*dbdx(2) + vadv*(rhoe+pres)*dbdy(2)
     fluxd(3,4) = uadv*(rhoe+pres)*dbdx(3) + vadv*(rhoe+pres)*dbdy(3)
     if(npoly.eq.2)then
      fluxd(4,4) = uadv*(rhoe+pres)*dbdx(4) + vadv*(rhoe+pres)*dbdy(4)
      fluxd(5,4) = uadv*(rhoe+pres)*dbdx(5) + vadv*(rhoe+pres)*dbdy(5)
      fluxd(6,4) = uadv*(rhoe+pres)*dbdx(6) + vadv*(rhoe+pres)*dbdy(6)
     endif
!
!
!finally, scatter the contribution to the RHS
!    
    do ideg = 1,ndegr
     rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem)+fluxd(ideg,1:nq)*djaco
    enddo
!
   enddo !...(2)ngausd

3000 enddo !...(1)ie=1,nelem

    return
end subroutine rhsinvdomndg_tria_curv
!
!...subroutine: domain integral for curved quad..
!
subroutine rhsinvdomndg_quad_curv(geoel,unkno,intfac,rhsel,coord,ipqua)
use constant
implicit none
real*8,dimension(1:ngeel,1:nsize),intent(in)::geoel
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndegr,1:nq,1:ncell),intent(inout)::rhsel
real*8,dimension(1:ndimn,1:npoin)::coord
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer*4,dimension(1:nifai,1:nafac)::intfac
!
!...local array
real*8::volel,wi
real*8:: r, s, djaco
real*8:: dxdr,dxds,dydr,dyds
real*8:: eps,c00,c10,c05,c20
real*8::xc,yc,xg,yg,dx,dy
real*8::dbdx(6), dbdy(6), b(6)
real*8::fluxd(1:ndegr,1:nq)
real*8::unknod(1:nq)
real*8,dimension(1:2, 1:npqua)::xpq
real*8,dimension(1:npqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8::rho,rhou,rhov,rhoe,rhom,uadv,vadv,pres
integer*4::ig,ideg,ielem,ishp
integer::ie
!
data eps / 1.0d-06 /
data c00 / 0.0d0 /
data c10 / 1.0d0 /
data c05 / 0.5d0 /
data c20 / 2.0d0 /
!
call ruqope(2, ngausdq, posiq, weighq)
!
!  print*,'domn',posi
dbdx = 0.d0
dbdy = 0.d0

do 3000 ie = 1,nquad !...(1)ie=1,nelem
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
!
elseif(ncurv==1)then
xpq(1, 1:npqua) = coord(1,ipqua(1:npqua, ie))
xpq(2, 1:npqua) = coord(2,ipqua(1:npqua, ie))
endif!
!...
xc    = geoel(1, ielem)
yc    = geoel(2, ielem)
volel = geoel(3, ielem)
dx    = geoel(4, ielem)*0.5d0
dy    = geoel(5, ielem)*0.5d0
!
dbdx(1)= 0.d0
dbdx(2)= 1.d0/dx
dbdx(3)= 0.d0

dbdy(1)= 0.d0
dbdy(2)= 0.d0
dbdy(3)= 1.0/dy
!
do ig=1,ngausdq !...(2)ngausd
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
b(6) = b(2)*b(3)   - geoel(8,ielem)
!
dbdx(4)= (xg-xc)/dx/dx
dbdx(5)= 0.d0
dbdx(6)= (yg-yc)/dy/dx

dbdy(4)= 0.d0
dbdy(5)= (yg-yc)/dy/dy
dbdy(6)= (xg-xc)/dx/dy
!
unknod = 0.d0
!
do ideg =1,mdegr
unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
enddo
!
rho  = unknod(1)
rhou = unknod(2)
rhov = unknod(3)
rhoe = unknod(4)
!
rhom = 1.d0/rho
uadv = rhou*rhom
vadv = rhov*rhom
pres = (gamma-1.d0)*(rhoe-0.5d0*rho*(uadv*uadv+vadv*vadv))
!
fluxd(1,1) = rhou*dbdx(1) + rhov*dbdy(1)
fluxd(2,1) = rhou*dbdx(2) + rhov*dbdy(2)
fluxd(3,1) = rhou*dbdx(3) + rhov*dbdy(3)
!
if(npoly.eq.2)then
fluxd(4,1) = rhou*dbdx(4) + rhov*dbdy(4)
fluxd(5,1) = rhou*dbdx(5) + rhov*dbdy(5)
fluxd(6,1) = rhou*dbdx(6) + rhov*dbdy(6)
endif

fluxd(1,2) = (uadv* rhou + pres)*dbdx(1) + vadv* rhou*dbdy(1)
fluxd(2,2) = (uadv* rhou + pres)*dbdx(2) + vadv* rhou*dbdy(2)
fluxd(3,2) = (uadv* rhou + pres)*dbdx(3) + vadv* rhou*dbdy(3)
if(npoly.eq.2)then
fluxd(4,2) = (uadv* rhou + pres)*dbdx(4) + vadv* rhou*dbdy(4)
fluxd(5,2) = (uadv* rhou + pres)*dbdx(5) + vadv* rhou*dbdy(5)
fluxd(6,2) = (uadv* rhou + pres)*dbdx(6) + vadv* rhou*dbdy(6)
endif

fluxd(1,3) = uadv* rhov*dbdx(1) + (vadv* rhov + pres)*dbdy(1)
fluxd(2,3) = uadv* rhov*dbdx(2) + (vadv* rhov + pres)*dbdy(2)
fluxd(3,3) = uadv* rhov*dbdx(3) + (vadv* rhov + pres)*dbdy(3)
if(npoly.eq.2)then
fluxd(4,3) = uadv* rhov*dbdx(4) + (vadv* rhov + pres)*dbdy(4)
fluxd(5,3) = uadv* rhov*dbdx(5) + (vadv* rhov + pres)*dbdy(5)
fluxd(6,3) = uadv* rhov*dbdx(6) + (vadv* rhov + pres)*dbdy(6)
endif

fluxd(1,4) = uadv*(rhoe+pres)*dbdx(1) + vadv*(rhoe+pres)*dbdy(1)
fluxd(2,4) = uadv*(rhoe+pres)*dbdx(2) + vadv*(rhoe+pres)*dbdy(2)
fluxd(3,4) = uadv*(rhoe+pres)*dbdx(3) + vadv*(rhoe+pres)*dbdy(3)
if(npoly.eq.2)then
fluxd(4,4) = uadv*(rhoe+pres)*dbdx(4) + vadv*(rhoe+pres)*dbdy(4)
fluxd(5,4) = uadv*(rhoe+pres)*dbdx(5) + vadv*(rhoe+pres)*dbdy(5)
fluxd(6,4) = uadv*(rhoe+pres)*dbdx(6) + vadv*(rhoe+pres)*dbdy(6)
endif
!
!
!finally, scatter the contribution to the RHS
!
do ideg = 1,ndegr
rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem)+fluxd(ideg,1:nq)*djaco
enddo
!
enddo !...(2)ngausd

3000 enddo !...(1)ie=1,nelem

return
end subroutine rhsinvdomndg_quad_curv
!
!...subroutine: inviscid flow --face integral for internal cubic curved face...
!
subroutine rhsinvbfacedg_curv2(geofa,geoel,bface,uchar,unkno,intfac,inpoel,rhsel,coord)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nelem+nbfac),intent(in)::geoel 
 real*8,intent(in)::uchar(1:nq)
 real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
 real*8,dimension(1:ndegr,1:nq,1:nelem),intent(out)::rhsel
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer*4,dimension(1:nifai,1:nafac)::intfac
 integer*4,dimension(1:nbfai,nbfac)::bface 
 integer*4,dimension(1:nvtri,1:nelem)::inpoel
 real*8,dimension(1:ngefa,1:nafac)::geofa
!
!...local array
 real*8,dimension(1:nq)::unkno1,unkno2
 real*8::xpin(2, 2)
 real*8::xp(2, nptfa)
 real*8,dimension(1:nptfa)::shp, dshpr
 real*8::weigh(ngausf), posi(1,ngausf)
 real*8::flux(1:nq)
 real*8::bl(6) 
!...local real number
 integer::ifa,iel,ier,ie,ifb,k,m,ig,ishp,ideg,jdeg
 real*8::xcl,ycl,dxl,dyl,rnx,rny
 real*8:: r, djaco, dxdr, dydr 
 real*8::xg,yg,wi
 real*8::rho1,rhou1,rhov1,rhoe1,rho2,rhou2,rhov2,rhoe2
 real*8::vn,uadv2,vadv2,dwav1,dwav2
 real*8::rhom1,uadv1,vadv1,pres1
 real*8::vnorm(ndimn),unkno_rie(1:nq,1:2)
!
 call rutope(1, ngausf, posi, weigh)
!
! print*,'gauss',posi(1,ngausf)
!
!...zero out rhs...
!
   rhsel = 0.d0
!
    do 250 ifa=1,nbfac !...(1)ifa=1,nafac
!
     iel=intfac(1,ifa)
     ier=intfac(2,ifa)
!
   if(ncurv==0)then
!
     xpin(1, 1:2) = coord(1, intfac(3:4, ifa)) !...One face constitutes of 'nptfa' points...
     xpin(2, 1:2) = coord(2, intfac(3:4, ifa)) 
!
     xp(1, 1:2) = xpin(1, 1:2)
     xp(2, 1:2) = xpin(2, 1:2)
     xp(1, 3) = 0.5d0*(xpin(1, 1) + xpin(1, 2))
     xp(2, 3) = 0.5d0*(xpin(2, 1) + xpin(2, 2))
   elseif(ncurv==2)then
!
     xp(1, 1:nptfa) = coord(1, intfac(3:(2+nptfa), ifa))  
     xp(2, 1:nptfa) = coord(2, intfac(3:(2+nptfa), ifa)) 
   endif
   !
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
!
!...2nd gauss loop...
!
     do ig =1, ngausf !...(4)ig=1,ngausf
!
     r  = posi(1, ig) 
    wi  = weigh(ig)
!
    shp(1) =  9.d0/16.d0*(1.d0-r)*(r+1.d0/3.d0)*(r-1.d0/3.d0)
    shp(2) = -9.d0/16.d0*(1.d0+r)*(1.d0/3.d0-r)*(r+1.d0/3.d0)   
    shp(3) =  27.d0/16.d0*(r-1.d0)*(r+1.d0)*(r-1.d0/3.d0)     
    shp(4) = -27.d0/16.d0*(r-1.d0)*(r+1.d0)*(r+1.d0/3.d0) 
!
    dshpr(1) = 9.d0/16.d0*(1.d0/9.d0-3.d0*r**2+2.d0*r)
    dshpr(2) =-9.d0/16.d0*(1.d0/9.d0-3.d0*r**2-2.d0*r) 
    dshpr(3) = 27.d0/16.d0*(-1.d0+3.d0*r**2-2.d0/3.d0*r)
    dshpr(4) =-27.d0/16.d0*(-1.d0+3.d0*r**2+2.d0/3.d0*r)
!
!...Jacobian determinant...
    dxdr = 0.d0
    dydr = 0.d0
!
    do ishp = 1, nptfa
     dxdr = dxdr + dshpr(ishp)*xp(1, ishp)
     dydr = dydr + dshpr(ishp)*xp(2, ishp)
    enddo 
!
    djaco = sqrt(dxdr**2 + dydr**2)
!
    dwav1 = dydr/djaco
    dwav2 =-dxdr/djaco
!
    rnx = dwav1
    rny = dwav2
!
     xg = 0.d0
     yg = 0.d0
!
   do ishp = 1, nptfa
     xg = xg + shp(ishp)*xp(1, ishp) 
     yg = yg + shp(ishp)*xp(2, ishp) 
    enddo
!
     bl(1) = 1.d0
     bl(2) = (xg-xcl)/dxl
     bl(3) = (yg-ycl)/dyl
     bl(4) = 0.5d0*bl(2)**2-geoel(6,iel)
     bl(5) = 0.5d0*bl(3)**2-geoel(7,iel)
     bl(6) = bl(2)*bl(3)   -geoel(8,iel)
!
!...zero out unkno1, unkno2
      unkno1 = 0.d0
      unkno2 = 0.d0
!
     do ideg = 1,ndegr
       unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
     enddo

 !   print*,'unkno1', unkno(1:3, 2, iel),iel, dxl,dyl
!
!...Specify unkno2 for different boundary condition...
    if(bface(3,ifa).eq.2)then
!   print*,'sym ifa',ifa,intfac(3:4,ifa)
!...Inviscid boundary wall
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!
!...  derived variables
      rhom1 = 1.d0/rho1
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
!
!      rnx = -xg/sqrt(xg**2+yg**2)
!      rny = -yg/sqrt(xg**2+yg**2)
!
!...  unknows for the ghost state
      rho2  = rho1
      rhoe2 = rhoe1
      vn    = uadv1*rnx+ vadv1*rny 
      uadv2 = uadv1 - 2.d0*vn*rnx
      vadv2 = vadv1 - 2.d0*vn*rny
      rhou2 = rho2*uadv2
      rhov2 = rho2*vadv2
!
      unkno2(1) = rho2
      unkno2(2) = rhou2
      unkno2(3) = rhov2
      unkno2(4) = rhoe2
    elseif(bface(3,ifa).eq.4)then
!   print*,'far',ifa,intfac(3:4,ifa)
!...Farfield 
      unkno2(1:4) = uchar(1:4)
!
!...Inflow...
!
    elseif(bface(3,ifa).eq.5)then
!
       unkno2(1:4) = uchar(1:4)
!
!...Outflow...
!
    elseif(bface(3,ifa).eq.6)then
!
       unkno2 = unkno1
!
     endif
!
!...call Riemann solver...
!     dwav1 = rnx
!     dwav2 = rny
     vnorm(1) = dwav1
     vnorm(2) = dwav2
!
     unkno_rie(1:nq,1) = unkno1(1:nq) 
     unkno_rie(1:nq,2) = unkno2(1:nq)   
!
!    print*,'left xg',nbfac,iel,ifa,unkno_rie(1:nq,1)
!    print*,'right xg',nbfac,ier,ifa,unkno_rie(1:nq,2)
!
     call riemann_2dhllc(unkno_rie, vnorm, flux)
!
!...adding to the RHS...
    do ideg =1 ,ndegr 
       rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*djaco*wi
    enddo  
!    
   enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nbfac

end subroutine rhsinvbfacedg_curv2
!
!...subroutine: inviscid flow --face integral for cubic curved internal face...
!
subroutine rhsinvifacedg_curv2(geofa,geoel,bface,unkno,intfac,inpoel,rhsel,coord)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel 
 real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
 real*8,dimension(1:ndegr,1:nq,1:nelem),intent(out)::rhsel
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer::ifa,iel,ier,ie,ifb,k,m,ig,ideg,jdeg
 integer*4,dimension(1:nifai,1:nafac)::intfac
 integer*4,dimension(1:nbfai,nbfac)::bface 
 integer*4,dimension(1:nvtri,1:nelem)::inpoel
 real*8,dimension(1:ngefa,1:nafac)::geofa
!
!...local array
 real*8,dimension(1:nq)::unkno1,unkno2
 real*8::unkno_rie(1:nq,1:2)
 real*8::xpin(2, 2)
 real*8::xp(2, nptfa)
 real*8,dimension(1:nptfa)::shp, dshpr
 real*8::weigh(ngausf), posi(1,ngausf)
 real*8::vnorm(1:2)
 real*8::flux(1:nq)
 real*8::bl(6),br(6)
!
!...local real number
 real*8::xcl,ycl,xcr,ycr,dxl,dyl,dxr,dyr,dxdr,dydr
 real*8:: r, djaco
 real*8::xg,yg,wi
 real*8::dwav1,dwav2
!
 integer ishp
!
 call rutope(1, ngausf, posi, weigh)
!
!
    do 250 ifa=nbfac+1,nafac !...(1)ifa=1,nafac
!
     iel=intfac(1,ifa)
     ier=intfac(2,ifa)
   !
  if(ncurv==0)then
!
     xpin(1, 1:2) = coord(1, intfac(3:4, ifa)) !...One face constitutes of 'nptfa' points...
     xpin(2, 1:2) = coord(2, intfac(3:4, ifa)) 
!
     xp(1, 1:2) = xpin(1, 1:2)
     xp(2, 1:2) = xpin(2, 1:2)
     xp(1, 3) = 0.5d0*(xpin(1, 1) + xpin(1, 2))
     xp(2, 3) = 0.5d0*(xpin(2, 1) + xpin(2, 2))
   elseif(ncurv==2)then
     xp(1, 1:nptfa) = coord(1, intfac(3:(2+nptfa), ifa))  
     xp(2, 1:nptfa) = coord(2, intfac(3:(2+nptfa), ifa)) 
   endif
   !
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
   !
     xcr = geoel(1,ier)
     ycr = geoel(2,ier)
     dxr = geoel(4,ier)*0.5d0
     dyr = geoel(5,ier)*0.5d0
!
!...2nd gauss loop...
!
     do ig =1,ngausf !...(4)ig=1,ngausf
!
     r  = posi(1, ig) 
    wi  = weigh(ig)
!
    shp(1) =  9.d0/16.d0*(1.d0-r)*(r+1.d0/3.d0)*(r-1.d0/3.d0)
    shp(2) = -9.d0/16.d0*(1.d0+r)*(1.d0/3.d0-r)*(r+1.d0/3.d0)   
    shp(3) =  27.d0/16.d0*(r-1.d0)*(r+1.d0)*(r-1.d0/3.d0)     
    shp(4) = -27.d0/16.d0*(r-1.d0)*(r+1.d0)*(r+1.d0/3.d0) 
!
    dshpr(1) = 9.d0/16.d0*(1.d0/9.d0-3.d0*r**2+2.d0*r)
    dshpr(2) =-9.d0/16.d0*(1.d0/9.d0-3.d0*r**2-2.d0*r) 
    dshpr(3) = 27.d0/16.d0*(-1.d0+3.d0*r**2-2.d0/3.d0*r)
    dshpr(4) =-27.d0/16.d0*(-1.d0+3.d0*r**2+2.d0/3.d0*r)
!
!...Jacobian determinant...
    dxdr = 0.d0
    dydr = 0.d0
!
    do ishp = 1, nptfa
     dxdr = dxdr + dshpr(ishp)*xp(1, ishp)
     dydr = dydr + dshpr(ishp)*xp(2, ishp)
    enddo 
!
    djaco = sqrt(dxdr**2 + dydr**2)
!
    dwav1 = dydr/djaco
    dwav2 =-dxdr/djaco
!
     xg = 0.d0
     yg = 0.d0
!
   do ishp = 1, nptfa
     xg = xg + shp(ishp)*xp(1, ishp) 
     yg = yg + shp(ishp)*xp(2, ishp) 
    enddo

     bl(1) = 1.d0
     bl(2) = (xg-xcl)/dxl
     bl(3) = (yg-ycl)/dyl
     bl(4) = 0.5d0*bl(2)**2-geoel(6,iel)
     bl(5) = 0.5d0*bl(3)**2-geoel(7,iel)
     bl(6) = bl(2)*bl(3)   -geoel(8,iel)
!
     br(1) = 1.d0
     br(2) = (xg-xcr)/dxr
     br(3) = (yg-ycr)/dyr
     br(4) = 0.5d0*br(2)**2-geoel(6,ier)
     br(5) = 0.5d0*br(3)**2-geoel(7,ier)
     br(6) = br(2)*br(3)   -geoel(8,ier)
!
!...
!...zero out unkno1, unkno2
      unkno1 = 0.d0
      unkno2 = 0.d0
!
     do ideg = 1,mdegr
       unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
       unkno2(1:nq) = unkno2(1:nq) + unkno(ideg,1:nq,ier)*br(ideg)
     enddo
!
!...call Riemann solver...
!
     vnorm(1) = dwav1
     vnorm(2) = dwav2
!
     unkno_rie(1:nq,1) = unkno1(1:nq) 
     unkno_rie(1:nq,2) = unkno2(1:nq)   
!
!    print*,'xg,ygiface', ifa, xp(2, 1:nptfa)
!
     call riemann_2dhllc(unkno_rie, vnorm, flux)
!
!...adding to the RHS...
    do ideg =1 ,ndegr 
!       rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
!       rhsel(ideg,1:nq,ier)=rhsel(ideg,1:nq,ier) + flux(1:nq)*br(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
    enddo
!
    do ideg =1 ,ndegr 
       rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*wi*djaco
       rhsel(ideg,1:nq,ier)=rhsel(ideg,1:nq,ier) + flux(1:nq)*br(ideg)*wi*djaco
    enddo    
!    
   enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nafac

end subroutine rhsinvifacedg_curv2
!
!...subroutine: domain integral for cubic curved triangle..
!
subroutine rhsinvdomndg_tria_curv2(geoel,unkno,intfac,rhsel,coord,inpoel)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nelem+nbfac),intent(in)::geoel 
 real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
 real*8,dimension(1:ndegr,1:nq,1:nelem),intent(inout)::rhsel
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer*4,dimension(1:nvtri,1:nelem)::inpoel
 integer*4,dimension(1:nifai,1:nafac)::intfac
!
!...local array
 real*8::volel,wi
 real*8:: r, s, t, djaco
 real*8:: dxdr,dxds,dydr,dyds 
 real*8:: eps,c00,c10,c05,c20
 real*8::xc,yc,xg,yg,dx,dy
 real*8::dbdx(6), dbdy(6), b(6)
 real*8::fluxd(1:ndegr,1:nq)
 real*8::unknod(1:nq)
 real*8,dimension(1:2, 1:nptri)::xp 
 real*8,dimension(1:nptri)::shp, dspr, dsps
 real*8:: weigh(ngausd), posi(2, ngausd)
 real*8::rho,rhou,rhov,rhoe,rhom,uadv,vadv,pres
 integer*4::ig,ideg,ielem,ishp
!
  data eps / 1.0d-06 / 
  data c00 / 0.0d0 / 
  data c10 / 1.0d0 / 
  data c05 / 0.5d0 / 
  data c20 / 2.0d0 / 
!
 call rutope(2, ngausd, posi, weigh)
!
!  print*,'domn',posi
    dbdx = 0.d0
    dbdy = 0.d0

   do 3000 ielem=1,nelem !...(1)ie=1,nelem
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
!...
    xc    = geoel(1, ielem)
    yc    = geoel(2, ielem)
    volel = geoel(3, ielem)
    dx    = geoel(4, ielem)*0.5d0
    dy    = geoel(5, ielem)*0.5d0
!
    dbdx(1)= 0.d0
    dbdx(2)= 1.d0/dx
    dbdx(3)= 0.d0

    dbdy(1)= 0.d0
    dbdy(2)= 0.d0
    dbdy(3)= 1.0/dy
!
    do ig=1,ngausd !...(2)ngausd
!
     r  = posi(1,ig)
     s  = posi(2,ig)
!
     t = 1.d0-r-s
!
    wi  = weigh(ig)
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
     b(6) = b(2)*b(3)   - geoel(8,ielem)
!
     dbdx(4)= (xg-xc)/dx/dx
     dbdx(5)= 0.d0
     dbdx(6)= (yg-yc)/dy/dx         

     dbdy(4)= 0.d0
     dbdy(5)= (yg-yc)/dy/dy
     dbdy(6)= (xg-xc)/dx/dy       
!  
     unknod = 0.d0
!
     do ideg =1,mdegr
        unknod(1:nq) = unknod(1:nq)+unkno(ideg,1:nq,ielem)*b(ideg)
     enddo
!
     rho  = unknod(1)
     rhou = unknod(2)
     rhov = unknod(3)
     rhoe = unknod(4)
!
     rhom = 1.d0/rho
     uadv = rhou*rhom
     vadv = rhov*rhom
     pres = (gamma-1.d0)*(rhoe-0.5d0*rho*(uadv*uadv+vadv*vadv))   
!
     fluxd(1,1) = rhou*dbdx(1) + rhov*dbdy(1) 
     fluxd(2,1) = rhou*dbdx(2) + rhov*dbdy(2) 
     fluxd(3,1) = rhou*dbdx(3) + rhov*dbdy(3) 
!
     if(npoly.eq.2)then
      fluxd(4,1) = rhou*dbdx(4) + rhov*dbdy(4) 
      fluxd(5,1) = rhou*dbdx(5) + rhov*dbdy(5) 
      fluxd(6,1) = rhou*dbdx(6) + rhov*dbdy(6) 
     endif

     fluxd(1,2) = (uadv* rhou + pres)*dbdx(1) + vadv* rhou*dbdy(1)
     fluxd(2,2) = (uadv* rhou + pres)*dbdx(2) + vadv* rhou*dbdy(2)
     fluxd(3,2) = (uadv* rhou + pres)*dbdx(3) + vadv* rhou*dbdy(3)
     if(npoly.eq.2)then
      fluxd(4,2) = (uadv* rhou + pres)*dbdx(4) + vadv* rhou*dbdy(4)
      fluxd(5,2) = (uadv* rhou + pres)*dbdx(5) + vadv* rhou*dbdy(5)
      fluxd(6,2) = (uadv* rhou + pres)*dbdx(6) + vadv* rhou*dbdy(6)
     endif

     fluxd(1,3) = uadv* rhov*dbdx(1) + (vadv* rhov + pres)*dbdy(1)
     fluxd(2,3) = uadv* rhov*dbdx(2) + (vadv* rhov + pres)*dbdy(2)
     fluxd(3,3) = uadv* rhov*dbdx(3) + (vadv* rhov + pres)*dbdy(3)
     if(npoly.eq.2)then
      fluxd(4,3) = uadv* rhov*dbdx(4) + (vadv* rhov + pres)*dbdy(4)
      fluxd(5,3) = uadv* rhov*dbdx(5) + (vadv* rhov + pres)*dbdy(5)
      fluxd(6,3) = uadv* rhov*dbdx(6) + (vadv* rhov + pres)*dbdy(6)
     endif

     fluxd(1,4) = uadv*(rhoe+pres)*dbdx(1) + vadv*(rhoe+pres)*dbdy(1)
     fluxd(2,4) = uadv*(rhoe+pres)*dbdx(2) + vadv*(rhoe+pres)*dbdy(2)
     fluxd(3,4) = uadv*(rhoe+pres)*dbdx(3) + vadv*(rhoe+pres)*dbdy(3)
     if(npoly.eq.2)then
      fluxd(4,4) = uadv*(rhoe+pres)*dbdx(4) + vadv*(rhoe+pres)*dbdy(4)
      fluxd(5,4) = uadv*(rhoe+pres)*dbdx(5) + vadv*(rhoe+pres)*dbdy(5)
      fluxd(6,4) = uadv*(rhoe+pres)*dbdx(6) + vadv*(rhoe+pres)*dbdy(6)
     endif
!
!
!finally, scatter the contribution to the RHS
!    
    do ideg = 1,ndegr
     rhsel(ideg,1:nq,ielem)=rhsel(ideg,1:nq,ielem)+fluxd(ideg,1:nq)*djaco
    enddo
!
   enddo !...(2)ngausd

3000 enddo !...(1)ie=1,nelem

    return
end subroutine rhsinvdomndg_tria_curv2

!
!...HLLC solver...
!
subroutine riemann_2dhllc(unint, vnorm, flux)
 use constant
 implicit none
 !real*8,dimension(1:6,1:nelem+nbfac)::unkno
 real*8,dimension(nq, 2),  intent(in)::unint 
 real*8,dimension(ndimn),  intent(in)::vnorm 
 real*8,dimension(nq),     intent(out)::flux
 integer:: iq
 
 real*8::gam1
 real*8 c00, c05, c10, c20, eps
 real*8  rho1,rhou1,rhov1,rhoe1
 real*8  rho2,rhou2,rhov2,rhoe2
 real*8  qadv1, qadv2
 real*8  rhom1,uadv1,vadv1,pres1,enth1,csou1
 real*8  rhom2,uadv2,vadv2,pres2,enth2,csou2
 real*8  rhsca,rmden,rhoav,uaver,vaver,entav
 real*8  qave5,cssa2,cssav,vdon1,vdon2,vnave
 real*8   s1,s2,dsv1,dsv2,sm,omeg1,omeg2,prsta,prst1,prst2
 real*8   rhol,rhoul,rhovl,rhoel
 real*8   rhor,rhour,rhovr,rhoer
 real*8::dwav1,dwav2
!
  data eps   / 1.0d-06 /
  data c00   / 0.0d0 /
  data c10   / 1.0d0 /
  data c05   / 0.5d0 /
  data c20   / 2.0d0 /
!
     gam1 = gamma- c10
!
!...Zero flux...
!
     flux = 0.d0
!
     dwav1 = vnorm(1)
     dwav2 = vnorm(2) 
!
!     print*,'hllc',dwav1, dwav2   
!
!get the fi+ iel
!
     rho1  = unint(1,1)
     rhou1 = unint(2,1)
     rhov1 = unint(3,1)
     rhoe1 = unint(4,1)
!
     rho2  = unint(1,2)
     rhou2 = unint(2,2)
     rhov2 = unint(3,2)
     rhoe2 = unint(4,2)    
!
     rhom1 = c10/rho1
     uadv1 = rhou1*rhom1
     vadv1 = rhov1*rhom1
     qadv1 = uadv1*uadv1 + vadv1*vadv1 
     pres1 = max(eps, gam1*(rhoe1 -                                &
                       c05*rho1*qadv1))
     enth1 = (rhoe1 + pres1)*rhom1
     csou1 = sqrt( max( eps,gamma*pres1*rhom1) )
!
     rhom2 = c10/rho2
     uadv2 = rhou2*rhom2
     vadv2 = rhov2*rhom2
     qadv2 = uadv2*uadv2 + vadv2*vadv2 
     pres2 = max(eps, gam1*( rhoe2 -                               &
                       c05*rho2*qadv2))
     enth2 = (rhoe2 + pres2)*rhom2
     csou2 = sqrt( max( eps,gamma*pres2*rhom2) )
!
!     print*,'hllc', pres1*rhom1*gamma
!     stop
!
!...  get the so-called Roe averaged variables
!
     rhsca = sqrt(rho2*rhom1)
      rmden = c10/(rhsca+c10)
!
      rhoav = rhsca*rho1
      uaver = (rhsca*uadv2 + uadv1)*rmden
      vaver = (rhsca*vadv2 + vadv1)*rmden
      entav = (rhsca*enth2 + enth1)*rmden
!
!...  averaged speed of sound
!
      qave5 = c05*(uaver*uaver + vaver*vaver)
      cssa2 = max(eps,gam1*(entav-qave5))
      cssav = sqrt(cssa2)
!
!...  compute the eigenvalues at the left, right, and Roe states
!
      vdon1 = uadv1*dwav1 + vadv1*dwav2  
      vdon2 = uadv2*dwav1 + vadv2*dwav2  
      vnave = uaver*dwav1 + vaver*dwav2  
!
!...  get the S_L and S_R defined by Eq. (13).
!
      s1    = min( vdon1-csou1, vnave-cssav )
      s2    = max( vdon2+csou2, vnave+cssav )
!
      dsv1  = s1 - vdon1
      dsv2  = s2 - vdon2
!
!...  compute the S_M defined by Eq. (12).
!
      sm    = (rho2*vdon2*dsv2 - rho1*vdon1*dsv1 + pres1-pres2)/  &
              (rho2*dsv2 - rho1*dsv1)
!
!...  compute the Omega_l, Omega_r, and p^*
!
      omeg1 = c10/(s1-sm)
      omeg2 = c10/(s2-sm)
      prsta = rho1*dsv1*(sm-vdon1) + pres1
!
      prst1 = prsta - pres1
      prst2 = prsta - pres2
!
!...  compute the U_l^\star, U_r^\star
!
      rhol  = omeg1* dsv1*rho1
      rhoul = omeg1*(dsv1*rhou1 + prst1*dwav1)   
      rhovl = omeg1*(dsv1*rhov1 + prst1*dwav2)   
      rhoel = omeg1*(dsv1*rhoe1 - pres1*vdon1 + prsta*sm )
!
      rhor  = omeg2* dsv2*rho2
      rhour = omeg2*(dsv2*rhou2 + prst2*dwav1)   
      rhovr = omeg2*(dsv2*rhov2 + prst2*dwav2)   
      rhoer = omeg2*(dsv2*rhoe2 - pres2*vdon2 + prsta*sm )
!
!...  Compute the fluxes according to the wave speed
!
      if(s1 .gt. c00) then
!
      flux(1)  = vdon1* rho1
      flux(2)  = vdon1* rhou1 + pres1*dwav1
      flux(3)  = vdon1* rhov1 + pres1*dwav2
      flux(4)  = vdon1*(rhoe1 + pres1)
!
      else if(s1 .le. c00 .and. sm .gt. c00) then
      flux(1)  = sm* rhol
      flux(2)  = sm* rhoul + prsta*dwav1
      flux(3)  = sm* rhovl + prsta*dwav2
      flux(4)  = sm*(rhoel + prsta)
!
      else if(sm .le. c00 .and. s2 .ge. c00) then
      flux(1)  = sm* rhor
      flux(2)  = sm* rhour + prsta*dwav1
      flux(3)  = sm* rhovr + prsta*dwav2
      flux(4)  = sm*(rhoer + prsta)
!
      else if(s2 .lt. c00) then
!
      flux(1)  = vdon2* rho2
      flux(2)  = vdon2* rhou2 + pres2*dwav1
      flux(3)  = vdon2* rhov2 + pres2*dwav2
      flux(4)  = vdon2*(rhoe2 + pres2)
!
      else
      write(*,*)'------------------------------------------------------'
      write(*,*)'Weird wave speed occured in sub riemann_2dhllc...'
      write(*,*)'Please check before continuing'
      write(*,*) 'rho1  = ', rho1
      write(*,*) 'rhou1 = ', rhou1
      write(*,*) 'rhov1 = ', rhov1
      write(*,*) 'rhoe1 = ', rhoe1
      write(*,*) 'pres1 = ', pres1
      write(*,*) 'enth1 = ', enth1
      write(*,*) 'csou1 = ', csou1
      write(*,*) 'rho2  = ', rho2
      write(*,*) 'rhou2 = ', rhou2
      write(*,*) 'rhov2 = ', rhov2
      write(*,*) 'rhoe2 = ', rhoe2
      write(*,*) 'pres2 = ', pres2
      write(*,*) 'enth2 = ', enth2
      write(*,*) 'csou2 = ', csou2
      write(*,*) 'vdon1 =', vdon1, uadv1, dwav1, vadv1, dwav2
      write(*,*) 'vdon2 =', vdon2
      write(*,*) 'vnave =', vnave
      write(*,*) 'cssav =', cssav
      write(*,*) 's1 = ', s1
      write(*,*) 's2 = ', s2
      write(*,*) 'sm = ', sm
      write(*,*) 'omeg1 =', omeg1
      write(*,*) 'omeg2 =', omeg2
      write(*,*) 'prsta =', prsta
      stop
    endif
end subroutine riemann_2dhllc
!
!...Lax-Friedrich solver...
!
subroutine riemann_2dlf(unint, vnorm, flux)
 use constant
 implicit none
 !real*8,dimension(1:6,1:nelem+nbfac)::unkno
 real*8,dimension(nq, 2),  intent(in)::unint 
 real*8,dimension(ndimn),  intent(in)::vnorm 
 real*8,dimension(nq),     intent(out)::flux
 integer:: iq
 
 real*8::gam1
 real*8 c00, c05, c10, c20, eps
 real*8  rho1,rhou1,rhov1,rhoe1
 real*8  rho2,rhou2,rhov2,rhoe2
 real*8  qadv1, qadv2
 real*8  rhom1,uadv1,vadv1,pres1,enth1,csou1
 real*8  rhom2,uadv2,vadv2,pres2,enth2,csou2
 real*8  vdon1,vdon2
 real*8::dwav1,dwav2
!
 real*8:: sradi, fl1, fl2, fl3, fl4, pre12
!
  data eps   / 1.0d-06 /
  data c00   / 0.0d0 /
  data c10   / 1.0d0 /
  data c05   / 0.5d0 /
  data c20   / 2.0d0 /
!
     gam1 = gamma- c10
!
!...Zero flux...
!
     flux = 0.d0
!
     dwav1 = vnorm(1)
     dwav2 = vnorm(2) 
!
!     print*,'hllc',dwav1, dwav2   
!
!get the fi+ iel
!
     rho1  = unint(1,1)
     rhou1 = unint(2,1)
     rhov1 = unint(3,1)
     rhoe1 = unint(4,1)
!
     rho2  = unint(1,2)
     rhou2 = unint(2,2)
     rhov2 = unint(3,2)
     rhoe2 = unint(4,2)    
!
     rhom1 = c10/rho1
     uadv1 = rhou1*rhom1
     vadv1 = rhov1*rhom1
     qadv1 = uadv1*uadv1 + vadv1*vadv1 
     pres1 = max(eps, gam1*(rhoe1 -                                &
                       c05*rho1*qadv1))
     enth1 = (rhoe1 + pres1)*rhom1
     csou1 = sqrt( max( eps,gamma*pres1*rhom1) )
!
     rhom2 = c10/rho2
     uadv2 = rhou2*rhom2
     vadv2 = rhov2*rhom2
     qadv2 = uadv2*uadv2 + vadv2*vadv2 
     pres2 = max(eps, gam1*( rhoe2 -                               &
                       c05*rho2*qadv2))
     enth2 = (rhoe2 + pres2)*rhom2
     csou2 = sqrt( max( eps,gamma*pres2*rhom2) )
!
     vdon1 = uadv1*dwav1 + vadv1*dwav2 
     vdon2 = uadv2*dwav1 + vadv2*dwav2 
!
!     print*,'Lax-friedrich', pres1*rhom1*gamma
!     stop
     sradi  = max(abs(vdon1)+csou1, abs(vdon2)+csou2)
!
     pre12 = pres1+pres2
!!
     fl1 = vdon1* rho1           + vdon2* rho2
     fl2 = vdon1* rhou1          + vdon2* rhou2 + dwav1*pre12
     fl3 = vdon1* rhov1          + vdon2* rhov2 + dwav2*pre12
     fl4 = vdon1*(rhoe1 + pres1) + vdon2*(rhoe2 + pres2)
!!
    flux(1) = 0.5d0*(fl1 - sradi*(rho2 -rho1 ))
    flux(2) = 0.5d0*(fl2 - sradi*(rhou2-rhou1))
    flux(3) = 0.5d0*(fl3 - sradi*(rhov2-rhov1))
    flux(4) = 0.5d0*(fl4 - sradi*(rhoe2-rhoe1))
!
end subroutine riemann_2dlf
!
!...Roe-Pike
!
subroutine riemann_2droe(unint, vnorm, flux)
  use constant
  implicit none
  !real*8,dimension(1:6,1:nelem+nbfac)::unkno
  real*8,dimension(nq, 2),  intent(in)::unint 
  real*8,dimension(ndimn),  intent(in)::vnorm 
  real*8,dimension(nq),     intent(out)::flux
  integer:: iq
 
  real*8::gam1
  real*8 c00, c05, c10, c20, eps
  real*8  rho1,rhou1,rhov1,rhoe1
  real*8  rho2,rhou2,rhov2,rhoe2
  real*8  qadv1, qadv2
  real*8  rhom1,uadv1,vadv1,pres1,enth1,csou1
  real*8  rhom2,uadv2,vadv2,pres2,enth2,csou2
  real*8  vdon1,vdon2
  real*8::dwav1,dwav2
!
  real*8  dpres, drho, duadv, dvadv, dvdon, dvdt1
  real*8  rhsca,rmden,rhoav,uaver,vaver,entav,csour,qaver,vdonr
!
  real*8,dimension(nq)::flux1,flux2,lamdr,lamdi,lamdj,dw,lamda1,lamda2,epsi
  real*8,dimension(nq)::alfa
  real*8,dimension(nq,nq)::vectr
!
  data eps   / 1.0d-06 /
  data c00   / 0.0d0 /
  data c10   / 1.0d0 /
  data c05   / 0.5d0 /
  data c20   / 2.0d0 /
!
     gam1 = gamma - c10
!
!...Zero flux...
!
     flux = 0.d0
!
     dwav1 = vnorm(1)
     dwav2 = vnorm(2)    
!
!get the fi+ iel
!
     rho1  = unint(1,1)
     rhou1 = unint(2,1)
     rhov1 = unint(3,1)
     rhoe1 = unint(4,1)
!
     rho2  = unint(1,2)
     rhou2 = unint(2,2)
     rhov2 = unint(3,2)
     rhoe2 = unint(4,2)    
!
     rhom1 = c10/rho1
     uadv1 = rhou1*rhom1
     vadv1 = rhov1*rhom1
     qadv1 = uadv1*uadv1 + vadv1*vadv1
     pres1 = max(eps, gam1*(rhoe1 -                                &
                       c05*rho1*qadv1))
     enth1 = (rhoe1 + pres1)*rhom1
     csou1 = sqrt( max( eps,gamma*pres1*rhom1) )
!
     rhom2 = c10/rho2
     uadv2 = rhou2*rhom2
     vadv2 = rhov2*rhom2
     qadv2 = uadv2*uadv2 + vadv2*vadv2 
     pres2 = max(eps, gam1*( rhoe2 -                               &
                       c05*rho2*qadv2))
     enth2 = (rhoe2 + pres2)*rhom2
     csou2 = sqrt( max( eps,gamma*pres2*rhom2) )
!
!...  get the so-called Roe averaged variables
!
      rhsca = sqrt(rho2*rhom1)
      rmden = c10/(rhsca+c10)
!
      rhoav = rhsca*rho1
      uaver = (rhsca*uadv2 + uadv1)*rmden
      vaver = (rhsca*vadv2 + vadv1)*rmden
      entav = (rhsca*enth2 + enth1)*rmden
!
!...  averaged speed of sound
!
      qaver = uaver*uaver + vaver*vaver  
      csour = sqrt(max(eps,gam1*(entav-c05*qaver)))
!
      vdon1 = uadv1*dwav1 + vadv1*dwav2  
      vdon2 = uadv2*dwav1 + vadv2*dwav2  
      vdonr = uaver*dwav1 + vaver*dwav2  
!
      lamdi(1) = vdon1-csou1
      lamdi(2) = vdon1
      lamdi(3) = vdon1
      lamdi(4) = vdon1+csou1

      lamdj(1) = vdon2-csou2
      lamdj(2) = vdon2
      lamdj(3) = vdon2
      lamdj(4) = vdon2+csou2

      lamdr(1) = vdonr-csour
      lamdr(2) = vdonr
      lamdr(3) = vdonr
      lamdr(4) = vdonr+csour
!       
!...Entropy fix...
!  
     do iq = 1, nq
        epsi(iq)=max(0.d0,lamdr(iq)-lamdi(iq),lamdj(iq)-lamdr(iq))
        if(abs(lamdr(iq)).ge.epsi(iq))then
         lamdr(iq)=lamdr(iq)
        else
         lamdr(iq)=epsi(iq)
       endif
     enddo
!
!...Right eigenvector...
!
      vectr(1,1) = 1.d0
      vectr(1,2) = uaver-csour*dwav1
      vectr(1,3) = vaver-csour*dwav2
      vectr(1,4) = entav-csour*vdonr
      
      vectr(2,1) = 1.d0
      vectr(2,2) = uaver
      vectr(2,3) = vaver
      vectr(2,4) = c05*qaver

      vectr(3,1) = 0.d0 
      vectr(3,2) = dwav2
      vectr(3,3) = -dwav1 
      vectr(3,4) = uaver*dwav2-vaver*dwav1
!
      vectr(4,1) = 1.d0
      vectr(4,2) = uaver+csour*dwav1
      vectr(4,3) = vaver+csour*dwav2
      vectr(4,4) = entav+csour*vdonr
!     
      drho   = rho2-rho1
      duadv  = uadv2 - uadv1
      dvadv  = vadv2 - vadv1
      dpres  = pres2 - pres1
      dvdt1  = duadv*dwav2-dvadv*dwav1 
      dvdon  = vdon2-vdon1
!      
      alfa(1) = c05*(dpres-rhoav*csour*dvdon)/csour/csour
      alfa(2) = drho-dpres/csour/csour
      alfa(3) = rhoav*dvdt1 
      alfa(4) = c05*(dpres+rhoav*csour*dvdon)/csour/csour
!
      flux1(1)  = vdon1* rho1
      flux1(2)  = vdon1* rhou1 + pres1*dwav1
      flux1(3)  = vdon1* rhov1 + pres1*dwav2
      flux1(4)  = vdon1*(rhoe1 + pres1)
!
      flux2(1)  = vdon2* rho2
      flux2(2)  = vdon2* rhou2 + pres2*dwav1
      flux2(3)  = vdon2* rhov2 + pres2*dwav2
      flux2(4)  = vdon2*(rhoe2 + pres2)
!
      do iq = 1, nq
        flux(:)=flux(:)+abs(lamdr(iq))*alfa(iq)*vectr(iq,:)
      enddo
  !    enddo
      
     do iq = 1, nq
       flux(iq) = 0.5d0*(flux1(iq)+flux2(iq)-flux(iq)) 
     enddo
!   
!  print*,rhsel(1:4,1:10)
!  stop
end subroutine riemann_2droe
!
!...Barth limiter...
!
subroutine barthlimit(intfac, ltelem, inpoel, unkno ,geoel,coord)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ngeel,1:nelem+nbfac), intent(in)::geoel
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord
integer*4,dimension(1:3,1:nelem),            intent(in)::ltelem
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
!
!...Local
!
integer, dimension(3)::estri

integer:: ie, iv, iest, iq, ideg,ind, ifa,ig
real*8:: xp(nvtri+ngausd), yp(nvtri+ngausd)
real*8:: aflim(1:nq, 1:nelem), alfal(1:nelem)
real*8:: unctr(1:nq)
real*8,  dimension(1:nq, 1:nvtri+ngausd)::alfa
real*8:: xv(nvtri+ngausd), yv(nvtri+ngausd)
real*8:: b(3, nvtri+ngausd)
real*8,  dimension(1:nq, 1:nvtri)::unvtx
real*8:: unmax(1:nq), unmin(1:nq), dunk(1:nq)
real*8,  dimension(1:nq,1:3)     ::unest
real*8,dimension(1:nq,  1:nvtri+ngausd) ::unknv
real*8::uncom(4)
real*8,allocatable::weigh(:), posi(:,:)
!
real*8:: shp1, shp2, shp3
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: rc, sc, dr, ds, fiy,dx,dy,xc,yc
!
eps = 1.e-6
!
!...shape functions
!
!
allocate (weigh(ngausd), posi(2,ngausd))
call rutope(2, ngausd, posi, weigh)
!
posi(1,1) = (-0.577350269189626d0+1.d0)/2.d0; posi(2,1) = 0.d0;
posi(1,2) = (0.577350269189626d0+1.d0)/2.d0; posi(2,2) = 0.0d0;
posi(2,3) = (-0.577350269189626d0+1.d0)/2.d0; posi(1,3) = 0.d0;
posi(2,4) = (0.577350269189626d0+1.d0)/2.d0; posi(1,4) = 0.0d0;!
posi(1,5) = posi(1,1); posi(2,5) = (0.577350269189626d0+1.d0)/2.d0;
posi(1,6) = (0.577350269189626d0+1.d0)/2.d0; posi(2,6) = (-0.577350269189626d0+1.d0)/2.d0;!
!
! if(ie==1) print*,unkno(1, 1:4, ie), unknv(1:4, 1)
!
do ie = 1, nelem
!
estri(1:3) = ltelem(1:3, ie)
!
unctr(1:nq) = unkno(1, 1:nq, ie)
!
! if(ie==1) print*,'update',posi
!
xp(1:nvtri) = coord(1,inpoel(1:nvtri,ie))
yp(1:nvtri) = coord(2,inpoel(1:nvtri,ie))
!
do ig=1,ngausd !...(2)ngausd
!
shp1=posi(1,ig)
shp2=posi(2,ig)
shp3=1.d0-shp1-shp2
!
xp(nvtri + ig) = xp(1)*shp1 + xp(2)*shp2 + xp(3)*shp3
yp(nvtri + ig) = yp(1)*shp1 + yp(2)*shp2 + yp(3)*shp3
!
enddo

!
xc = geoel(1,ie)
yc = geoel(2,ie)
!
dx = geoel(4,ie)*0.5d0
dy = geoel(5,ie)*0.5d0
!
do iv =1 ,nvtri+ngausd
!...Left cell + intfac(3,ifa)
b(1, iv) = 1.d0
b(2, iv) = (xp(iv)-xc)/dx
b(3, iv) = (yp(iv)-yc)/dy
enddo
!
!...zero out unknv
!
unknv = 0.d0
!
do iv   = 1,nvtri+ngausd
do ideg = 1,mdegr
unknv(1:nq, iv) = unknv(1:nq, iv) + unkno(ideg,1:nq,ie)*b(ideg, iv)
enddo
enddo
!
!
! if(ie==1) print*,'update2',unctr, unknv(1:4, 1),estri
!
!...Find the unknown values at the adjacent cells...
!...Pres is stored in the original posiiton for total energy...
!
do iest = 1, 3
!
if(estri(iest).gt.nelem)then
!  if(ie==1) print*,'estri', ie,estri(iest),iest
unest(1:nq, iest) = unctr(1:nq) !-1.d13  !  !Set unknown as a very small number
!  if(ie==1) print*,'estri2', unest(1:nq, iest)
else
unest(1:nq, iest) = unkno(1, 1:nq, estri(iest))
endif
enddo
!
!  if(ie==1) print*,'estri3', unest(1:nq,1), unest(1:nq,2), unest(1:nq,3)
!
!...Find the maximum values...
!
do iq =1, nq
!
uncom(1) = maxval(unest(iq, 1:3))
uncom(2) = maxval(unest(iq, 1:3))
uncom(3) = maxval(unest(iq, 1:3))
uncom(4) = maxval(unest(iq, 1:3))

unmax(iq) = max(uncom(1), unctr(iq))
!
enddo
!
!  if(ie==1) print*,'estri4', unmax
!
do iest = 1, 3
!
if(estri(iest).gt.nelem)then
unest(1:nq, iest) =unctr(1:nq) !1.d13 !  unctr(1:nq) ! !Set unknown as a very small number
endif
enddo
!
!...Find the maximum values...
!
do iq =1, nq
!
uncom(1) = minval(unest(iq, 1:3))
uncom(2) = minval(unest(iq, 1:3))
uncom(3) = minval(unest(iq, 1:3))
uncom(4) = minval(unest(iq, 1:3))
!
unmin(iq) = min(uncom(1),unctr(iq))
!
enddo
!  if(ie==1) print*,'estri4min', unmin
!
!...Impose limiter
!
do iv = 1, nvtri+ngausd
!
do iq = 1,nq
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq) - unctr(iq)),(unmax(iq) - unctr(iq))/dunk(iq),(unmin(iq) - unctr(iq)),&
!                        (unmin(iq) - unctr(iq))/dunk(iq)
!
if(dunk(iq).gt.1.d-8)then

fiy = (unmax(iq) - unctr(iq))/dunk(iq)
alfa(iq, iv) = min(1.d0, 1.0d0*fiy)
!
!  alfa(iq, iv) = max(min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)), 0.d0)
elseif(dunk(iq).lt.-1.d-8)then

fiy = (unmin(iq) - unctr(iq))/dunk(iq)
alfa(iq, iv) = min(1.d0, 1.0d0*fiy)
!
!   alfa(iq, iv) = max(min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)), 0.d0)
else
!
alfa(iq, iv) = 1.d0
endif
!
! if(ie==1) print*,'dunk alfa',iq,iv,alfa(iq,iv)
!
enddo
enddo
!
!...Get the minimum alfa...
!
do iq = 1,nq
aflim(iq, ie) = minval(alfa(iq, 1:nvtri+ngausd))
enddo
!
! alfal(ie) = minval(aflim(1:nq, ie))
!
! if(ie==80) print*,'alfa',aflim(2:3,ie)!,alfa(1, 1:nvtri),alfa(2, 1:nvtri),alfa(3, 1:nvtri),alfa(4, 1:nvtri)
!
!
enddo
!
!
! do ifa =1 ,nbfac
!    aflim(1:nq, intfac(3, ifa)) = 1.d0
! enddo
!
! aflim = 0.d0
!
!
do ie = 1, nelem
!
unkno(1, 1:nq, ie) = unkno(1, 1:nq, ie)
!
unkno(2:3, 1, ie) = unkno(2:3, 1, ie)*aflim(1,ie)
unkno(2:3, 2, ie) = unkno(2:3, 2, ie)*aflim(2,ie)
unkno(2:3, 3, ie) = unkno(2:3, 3, ie)*aflim(3,ie)
unkno(2:3, 4, ie) = unkno(2:3, 4, ie)*aflim(4,ie)
!
enddo
end subroutine barthlimit

!...FVM primitive rhs
!
!
!...subroutine: inviscid flow --face integral for internal face...
!
subroutine rhsinvbface_fvm_prim(geofa,geoel,bface,uchar,unkno,intfac,inpoel,rhsel,coord)
use constant
implicit none
real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel
real*8,intent(in)::uchar(1:nq)
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndegr,1:nq,1:nelem),intent(out)::rhsel
real*8,dimension(1:ndimn,1:npoin)::coord
integer::ifa,iel,ier,ie,ifb,k,m,ig,ideg,jdeg
integer*4,dimension(1:nifai,1:nafac)::intfac
integer*4,dimension(1:nbfai,nbfac)::bface
integer*4,dimension(1:nvtri,1:nelem)::inpoel
real*8,dimension(1:ngefa,1:nafac)::geofa
real*8,dimension(1:nq)::unkno1,unkno2
!
!...local array
real*8::flux(1:nq)
real*8::bl(6)
!...local real number
real*8::xcl,ycl,dxl,dyl,rnx,rny
real*8::x1,x2,y1,y2,xi,shp1,shp2,xg,yg,wi
real*8::rho1,rhou1,rhov1,rhoe1,rho2,rhou2,rhov2,rhoe2
real*8::vn,uadv2,vadv2,dwav1,dwav2
real*8::rhom1,uadv1,vadv1,pres1
real*8::vnorm(ndimn),unkno_rie(1:nq,1:2)

real*8,allocatable::weigh(:), posi(:,:)

allocate (weigh(ngausf), posi(1,ngausf))
call rutope(1, ngausf, posi, weigh)
!
!...zero out rhs...
!
rhsel = 0.d0
!
!...Start loop of boundary faces for RHS....
!
do 250 ifa=1,nbfac !...(1)ifa=1,nafac
!
iel=intfac(1,ifa)
ier=intfac(2,ifa)
!
x1 = coord(1,intfac(3,ifa))
y1 = coord(2,intfac(3,ifa))

x2 = coord(1,intfac(4,ifa))
y2 = coord(2,intfac(4,ifa))
!
xcl = geoel(1,iel)
ycl = geoel(2,iel)
dxl = geoel(4,iel)*0.5d0
dyl = geoel(5,iel)*0.5d0
!
rnx = geofa(1,ifa)
rny = geofa(2,ifa)
!
!...2nd gauss loop...
!
do ig =1,ngausf !...(4)ig=1,ngausf
!
xi =posi(1, ig)

shp1 = 0.5d0*(1.d0-xi)
shp2 = 0.5d0*(1.d0+xi)

xg = x1*shp1 + x2*shp2
yg = y1*shp1 + y2*shp2

bl(1) = 1.d0
bl(2) = (xg-xcl)/dxl
bl(3) = (yg-ycl)/dyl
bl(4) = 0.5d0*bl(2)**2-geoel(6,iel)
bl(5) = 0.5d0*bl(3)**2-geoel(7,iel)
bl(6) = bl(2)*bl(3)   -geoel(8,iel)
!
!...
!...zero out unkno1, unkno2
unkno1 = 0.d0
unkno2 = 0.d0
!
do ideg = 1,ndegr
unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
enddo
!
rho1  = unkno1(1)
uadv1 = unkno1(2)
vadv1 = unkno1(3)
pres1 = unkno1(4)
rhoe1 = pres1/(gamma-1.d0) + 0.5d0*rho1*(uadv1**2 + vadv1**2)
!
unkno1(1) = rho1
unkno1(2) = rho1*uadv1
unkno1(3) = rho1*vadv1
unkno1(4) = rhoe1
!
!...Specify unkno2 for different boundary condition...
if(bface(3,ifa).eq.2)then
!
!...Inviscid boundary wall
!
!
!      rnx = -xg/sqrt(xg**2+yg**2)
!      rny = -yg/sqrt(xg**2+yg**2)
!
!...  unknows for the ghost state
rho2  = rho1
rhoe2 = rhoe1
vn    = uadv1*rnx+ vadv1*rny
uadv2 = uadv1 - 2.d0*vn*rnx
vadv2 = vadv1 - 2.d0*vn*rny
rhou2 = rho2*uadv2
rhov2 = rho2*vadv2
!
unkno2(1) = rho2
unkno2(2) = rhou2
unkno2(3) = rhov2
unkno2(4) = rhoe2
elseif(bface(3,ifa).eq.4)then
!
!...Farfield
unkno2(1:4) = uchar(1:4)
!
endif
!
!...call Riemann solver...
dwav1 = rnx
dwav2 = rny
vnorm(1) = dwav1
vnorm(2) = dwav2
!
unkno_rie(1:nq,1) = unkno1(1:nq)
unkno_rie(1:nq,2) = unkno2(1:nq)
!
call riemann_2dhllc(unkno_rie, vnorm, flux)
!
!...adding to the RHS...
do ideg =1 ,ndegr
rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
enddo
!
enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nbfac

end subroutine rhsinvbface_fvm_prim
!
!...subroutine: inviscid flow --face integral for internal face...
!
subroutine rhsinviface_fvm_prim(geofa,geoel,bface,unkno,intfac,inpoel,rhsel,coord)
use constant
implicit none
real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
real*8,dimension(1:ndegr,1:nq,1:nelem),intent(out)::rhsel
real*8,dimension(1:ndimn,1:npoin)::coord
integer::ifa,iel,ier,ie,ifb,k,m,ig,ideg,jdeg
integer*4,dimension(1:nifai,1:nafac)::intfac
integer*4,dimension(1:nbfai,nbfac)::bface
integer*4,dimension(1:nvtri,1:nelem)::inpoel
real*8,dimension(1:ngefa,1:nafac)::geofa
real*8,dimension(1:nq)::unkno1,unkno2
!...local array
real*8::unkno_rie(1:nq,1:2)
real*8::vnorm(1:2)
real*8::flux(1:nq)
real*8::bl(6),br(6)
!...local real number
real*8::xcl,ycl,xcr,ycr,dxl,dyl,dxr,dyr,rnx,rny
real*8::x1,x2,y1,y2,xi,shp1,shp2,xg,yg,wi
real*8::dwav1,dwav2
real*8::rho1,uadv1,vadv1,pres1,rhoe1
real*8::rho2,uadv2,vadv2,pres2,rhoe2

real*8,allocatable::weigh(:), posi(:,:)

allocate (weigh(ngausf), posi(1,ngausf))
call rutope(1, ngausf, posi, weigh)
!
!
!
do 250 ifa=nbfac+1,nafac !...(1)ifa=1,nafac
!
iel=intfac(1,ifa)
ier=intfac(2,ifa)
!
x1 = coord(1,intfac(3,ifa))
y1 = coord(2,intfac(3,ifa))

x2 = coord(1,intfac(4,ifa))
y2 = coord(2,intfac(4,ifa))
!
xcl = geoel(1,iel)
ycl = geoel(2,iel)
dxl = geoel(4,iel)*0.5d0
dyl = geoel(5,iel)*0.5d0
!
xcr = geoel(1,ier)
ycr = geoel(2,ier)
dxr = geoel(4,ier)*0.5d0
dyr = geoel(5,ier)*0.5d0
!
rnx = geofa(1,ifa)
rny = geofa(2,ifa)
!
!...2nd gauss loop...
!
do ig =1,ngausf !...(4)ig=1,ngausf
!
xi =posi(1, ig)

shp1 = 0.5d0*(1.d0-xi)
shp2 = 0.5d0*(1.d0+xi)

xg = x1*shp1 + x2*shp2
yg = y1*shp1 + y2*shp2

bl(1) = 1.d0
bl(2) = (xg-xcl)/dxl
bl(3) = (yg-ycl)/dyl
bl(4) = 0.5d0*bl(2)**2-geoel(6,iel)
bl(5) = 0.5d0*bl(3)**2-geoel(7,iel)
bl(6) = bl(2)*bl(3) - geoel(8,iel)
!
br(1) = 1.d0
br(2) = (xg-xcr)/dxr
br(3) = (yg-ycr)/dyr
br(4) = 0.5d0*br(2)**2-geoel(6,ier)
br(5) = 0.5d0*br(3)**2-geoel(7,ier)
br(6) = br(2)*br(3) - geoel(8,ier)
!
!...
!...zero out unkno1, unkno2
unkno1 = 0.d0
unkno2 = 0.d0
!
do ideg = 1,mdegr
unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
unkno2(1:nq) = unkno2(1:nq) + unkno(ideg,1:nq,ier)*br(ideg)
enddo
!
rho1  = unkno1(1)
uadv1 = unkno1(2)
vadv1 = unkno1(3)
pres1 = unkno1(4)
rhoe1 = pres1/(gamma-1.d0) + 0.5d0*rho1*(uadv1**2 + vadv1**2)
!
unkno1(1) = rho1
unkno1(2) = rho1*uadv1
unkno1(3) = rho1*vadv1
unkno1(4) = rhoe1
!
rho2  = unkno2(1)
uadv2 = unkno2(2)
vadv2 = unkno2(3)
pres2 = unkno2(4)
rhoe2 = pres2/(gamma-1.d0) + 0.5d0*rho2*(uadv2**2 + vadv2**2)
!
unkno2(1) = rho2
unkno2(2) = rho2*uadv2
unkno2(3) = rho2*vadv2
unkno2(4) = rhoe2
!
!
!...call Riemann solver...
!
dwav1 = rnx
dwav2 = rny
vnorm(1) = dwav1
vnorm(2) = dwav2
!
unkno_rie(1:nq,1) = unkno1(1:nq)
unkno_rie(1:nq,2) = unkno2(1:nq)
!
call riemann_2dhllc(unkno_rie, vnorm, flux)
!
!...adding to the RHS...
do ideg =1 ,ndegr
rhsel(ideg,1:nq,iel)=rhsel(ideg,1:nq,iel) - flux(1:nq)*bl(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
rhsel(ideg,1:nq,ier)=rhsel(ideg,1:nq,ier) + flux(1:nq)*br(ideg)*weigh(ig)*geofa(5,ifa)/2.d0
enddo
!
enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nafac

end subroutine rhsinviface_fvm_prim
