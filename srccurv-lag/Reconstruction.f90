!
!...subroutine: least square for FVM through loop over faces...
!
     subroutine getgraduls(intfac,geoel ,geofa, unkno ,gradu,coord)
      use constant
      implicit none 
      real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel 
      real*8,dimension(1:7,1:nafac),intent(in)::geofa
      real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac)::unkno
      real*8,dimension(1:2,1:3,1:nelem),intent(out)::gradu
      real*8,dimension(1:ndimn,1:npoin)::coord 
      integer*4,dimension(1:8,1:nafac)::intfac
      integer*4::iel,ier,ifa,ncs,ie
      real*8::c00
      real*8::geogh(1:2,1:nbfac)
      real*8::dx,dy,weigh
      real*8::dxwei,dywei
      real*8::du1we,du2we,du3we,du4we
      real*8::rxx,ryy,rxy
      real*8::rco1,rco2,rco3,rco4
      real*8::ru1x,ru2x,ru3x,ru4x
      real*8::ru1y,ru2y,ru3y,ru4y
      real*8::rjac1
      real*8::x1,y1,x2,y2,x0,y0
      real*8::la,lb,lc,abp2,abm2,xg,yg
      integer::ip1,ip2
      real*8, allocatable::dmets(:,:)
!
!...  This sub computes the gradient using the least-square approach
!
      allocate(dmets(1:3,1:nelem+nbfac))

      c00 = 0.0
      dmets(:,:) = 0.0
      ncs=0
!
      gradu = 0.0
!
    do 1400 ifa = 1, nafac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
      ier   = intfac(2,ifa)
     
      if(ifa.le.nbfac) then
        dx  = geofa(3,ifa) - geoel(1,iel)
        dy  = geofa(4,ifa) - geoel(2,iel)
      else 
        dx  = geoel(1,ier) - geoel(1,iel)
        dy  = geoel(2,ier) - geoel(2,iel)
      endif
!
!...  weighting for this edge
!
!     weigh = 1.0
      weigh = 1.0/sqrt(dx**2 + dy**2)
!
      dxwei = weigh*dx
      dywei = weigh*dy
!
      dmets(1,iel) = dmets(1,iel) + dxwei*dxwei 
      dmets(2,iel) = dmets(2,iel) + dywei*dywei 
      dmets(3,iel) = dmets(3,iel) + dxwei*dywei 

!
      dmets(1,ier) = dmets(1,ier) + dxwei*dxwei 
      dmets(2,ier) = dmets(2,ier) + dywei*dywei 
      dmets(3,ier) = dmets(3,ier) + dxwei*dywei
    

     du1we = weigh*(unkno(1,1,ier) - unkno(1,1,iel))
     du2we = weigh*(unkno(1,2,ier) - unkno(1,2,iel))
     du3we = weigh*(unkno(1,3,ier) - unkno(1,3,iel))

!
     gradu(1,1,iel) = gradu(1,1,iel) + dxwei*du1we 
     gradu(1,2,iel) = gradu(1,2,iel) + dxwei*du2we 
     gradu(1,3,iel) = gradu(1,3,iel) + dxwei*du3we 

     gradu(2,1,iel) = gradu(2,1,iel) + dywei*du1we 
     gradu(2,2,iel) = gradu(2,2,iel) + dywei*du2we 
     gradu(2,3,iel) = gradu(2,3,iel) + dywei*du3we 
!  excluding the ghost element gradient
    if(ifa.gt.nbfac) then
     gradu(1,1,ier) = gradu(1,1,ier) + dxwei*du1we 
     gradu(1,2,ier) = gradu(1,2,ier) + dxwei*du2we 
     gradu(1,3,ier) = gradu(1,3,ier) + dxwei*du3we 

     gradu(2,1,ier) = gradu(2,1,ier) + dywei*du1we 
     gradu(2,2,ier) = gradu(2,2,ier) + dywei*du2we 
     gradu(2,3,ier) = gradu(2,3,ier) + dywei*du3we 
    endif
!
 1400 continue
!
!
!
      do 2000 ie  = 1, nelem
!
      rxx         = dmets(1,ie)
      ryy         = dmets(2,ie)
      rxy         = dmets(3,ie)
!
      rco1  = ryy 
      rco2  = -rxy
!nafac
      rco3  = -rxy
      rco4  = rxx

      rjac1       = 1.0/(rxx*ryy - rxy*rxy)
!
!
     ru1x        = gradu(1,ncs+1,ie)
     ru1y        = gradu(2,ncs+1,ie)

     ru2x        = gradu(1,ncs+2,ie)
     ru2y        = gradu(2,ncs+2,ie)
!
     ru3x        = gradu(1,ncs+3,ie)
     ru3y        = gradu(2,ncs+3,ie)
!

     gradu(1,ncs+1,ie) = (ru1x*rco1 + ru1y*rco2 )*rjac1 
     gradu(2,ncs+1,ie) = (ru1x*rco3 + ru1y*rco4 )*rjac1 


     gradu(1,ncs+2,ie) = (ru2x*rco1 + ru2y*rco2 )*rjac1 
     gradu(2,ncs+2,ie) = (ru2x*rco3 + ru2y*rco4 )*rjac1 


     gradu(1,ncs+3,ie) = (ru3x*rco1 + ru3y*rco2 )*rjac1 
     gradu(2,ncs+3,ie) = (ru3x*rco3 + ru3y*rco4 )*rjac1 
!
      if(ie==0)then
        print*,'output'
        print*,gradu(1:2,ncs+1,ie)
        print*,ru1x,ru1y
        print*,rco1,rco2,rco3,rco4  
      endif
 2000 continue
!

      deallocate(dmets)
      return     
      end subroutine getgraduls
!
!...subroutine: least square for FVM through loop over elements...
!
     subroutine getgraduls3(intfac,geoel ,geofa, unkno ,gradu,coord,esuel)
      use constant
      implicit none 
      real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel 
      real*8,dimension(1:7,1:nafac),intent(in)::geofa
      real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac)::unkno
      real*8,dimension(1:2,1:3,1:nelem),intent(out)::gradu
      real*8,dimension(1:ndimn,1:npoin)::coord 
      integer:: esuel(1:3,1:nelem)
      integer*4,dimension(1:8,1:nafac)::intfac
      integer*4::iel,ier,ifa,ncs,ie
      real*8::c00
      real*8::geogh(1:2,1:nbfac)
!
      real*8::dx,dy,weigh
      real*8::dxwei,dywei
      real*8::du1we,du2we,du3we,du4we
      real*8::rxx,ryy,rxy
      real*8::rco1,rco2,rco3,rco4
      real*8::ru1x,ru2x,ru3x,ru4x
      real*8::ru1y,ru2y,ru3y,ru4y
      real*8::rjac1
      real*8::x1,y1,x2,y2,x0,y0
      real*8::la,lb,lc,abp2,abm2,xg,yg
      integer::ip1,ip2,inf
      real*8, allocatable::dmets(:,:)
!
!...  This sub computes the gradient using the least-square approach
!
      allocate(dmets(1:3,1:nelem+nbfac))

      c00 = 0.0
      dmets(:,:) = 0.0
      ncs=0
!
      gradu = 0.0
!
      do 1600 ifa = 1, nbfac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
      ier   = intfac(2,ifa)
     
      dx  = geofa(3,ifa) - geoel(1,iel)
      dy  = geofa(4,ifa) - geoel(2,iel)
!
!...  weighting for this edge
!
      weigh = 1.0
!     weigh = 1.0/sqrt(dx**2 + dy**2)
!
      dxwei = weigh*dx
      dywei = weigh*dy
!
      dmets(1,iel) = dmets(1,iel) + dxwei*dxwei 
      dmets(2,iel) = dmets(2,iel) + dywei*dywei 
      dmets(3,iel) = dmets(3,iel) + dxwei*dywei 
!
     du1we = weigh*(unkno(1,1,ier) - unkno(1,1,iel))
     du2we = weigh*(unkno(1,2,ier) - unkno(1,2,iel))
     du3we = weigh*(unkno(1,3,ier) - unkno(1,3,iel))
!
     gradu(1,1,iel) = gradu(1,1,iel) + dxwei*du1we 
     gradu(1,2,iel) = gradu(1,2,iel) + dxwei*du2we 
     gradu(1,3,iel) = gradu(1,3,iel) + dxwei*du3we 

     gradu(2,1,iel) = gradu(2,1,iel) + dywei*du1we 
     gradu(2,2,iel) = gradu(2,2,iel) + dywei*du2we 
     gradu(2,3,iel) = gradu(2,3,iel) + dywei*du3we 
!
 1600 continue
!
!
!
     do 1500 ie = 1, nelem
 
     do inf = 1,3
      
        iel = ie
        ier = esuel(inf,ie)
        if(ier.le.nelem)then
        dx  = geoel(1,ier) - geoel(1,iel)
        dy  = geoel(2,ier) - geoel(2,iel)

!      weigh = 1.0/sqrt(dx**2 + dy**2)
       weigh =1.0
!
      dxwei = weigh*dx
      dywei = weigh*dy
!
      dmets(1,iel) = dmets(1,iel) + dxwei*dxwei 
      dmets(2,iel) = dmets(2,iel) + dywei*dywei 
      dmets(3,iel) = dmets(3,iel) + dxwei*dywei 
    

      du1we = weigh*(unkno(1,1,ier) - unkno(1,1,iel))
      du2we = weigh*(unkno(1,2,ier) - unkno(1,2,iel))
      du3we = weigh*(unkno(1,3,ier) - unkno(1,3,iel))

!
      gradu(1,1,iel) = gradu(1,1,iel) + dxwei*du1we 
      gradu(1,2,iel) = gradu(1,2,iel) + dxwei*du2we 
      gradu(1,3,iel) = gradu(1,3,iel) + dxwei*du3we 

      gradu(2,1,iel) = gradu(2,1,iel) + dywei*du1we 
      gradu(2,2,iel) = gradu(2,2,iel) + dywei*du2we 
      gradu(2,3,iel) = gradu(2,3,iel) + dywei*du3we 
      endif      
     enddo
1500 enddo
!
! print*,'du1we', dmets(1:3,1)
!
      do 2000 ie  = 1, nelem
!
      rxx         = dmets(1,ie)
      ryy         = dmets(2,ie)
      rxy         = dmets(3,ie)
!
      rco1  = ryy 
      rco2  = -rxy
!nafac
      rco3  = -rxy
      rco4  = rxx

      rjac1       = 1.0/(rxx*ryy - rxy*rxy)
!
!
     ru1x        = gradu(1,ncs+1,ie)
     ru1y        = gradu(2,ncs+1,ie)

     ru2x        = gradu(1,ncs+2,ie)
     ru2y        = gradu(2,ncs+2,ie)
!
     ru3x        = gradu(1,ncs+3,ie)
     ru3y        = gradu(2,ncs+3,ie)
!

     gradu(1,ncs+1,ie) = (ru1x*rco1 + ru1y*rco2 )*rjac1 
     gradu(2,ncs+1,ie) = (ru1x*rco3 + ru1y*rco4 )*rjac1 


     gradu(1,ncs+2,ie) = (ru2x*rco1 + ru2y*rco2 )*rjac1 
     gradu(2,ncs+2,ie) = (ru2x*rco3 + ru2y*rco4 )*rjac1 


     gradu(1,ncs+3,ie) = (ru3x*rco1 + ru3y*rco2 )*rjac1 
     gradu(2,ncs+3,ie) = (ru3x*rco3 + ru3y*rco4 )*rjac1 
 !     print*,'output'
!
!
 2000 continue
!
!...
!
     deallocate(dmets)
      if(nreco.lt.0)then
         do ie =1,nelem
            gradu(1,1,ie) = unkno(1,2,ie)
            gradu(2,1,ie) = unkno(1,3,ie)
         enddo
!
!...nreco==-3, DG(P0P1) + DG(P0)
!
        if(nreco==-3)then
         do ie =1,nelem
            gradu(1:2,2:3,ie) = 0.d0
         enddo
         endif 
      endif
      return  
  end subroutine getgraduls3
!
!...subroutine: least square for rDG(P1P2)...
!
 subroutine getgradulsp1p2(intfac,geoel ,geofa, unkno ,coord,invmat)
  use constant
  implicit none 
 real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel 
 real*8,dimension(1:7,1:nafac),intent(in)::geofa
 real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac)::unkno
 real*8,dimension(6,nelem),intent(in)::invmat
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer*4,dimension(1:8,1:nafac)::intfac
 integer*4::iel,ier,ifa,ncs,ie,ngausl,ig
 real*8::c00
!
  real*8::xcl,ycl,dxl,dyl,rxxl,ryyl,rxyl
  real*8::xcr,ycr,dxr,dyr,rxxr,ryyr,rxyr
  real*8::b2l, b3l, b4l, b5l, b6l 
  real*8::b2r, b3r, b4r, b5r, b6r
  real*8::pl,pxl,pyl,ql,qxl,qyl,xg,yg
  real*8::pr,pxr,pyr,qr,qxr,qyr
  real*8::rl11, rl12, rl21, rl22, rl31, rl32
  real*8::rr11, rr12, rr21, rr22, rr31, rr32
  real*8::rhs11, rhs21, rhs31, rhs12, rhs22, rhs32
  real*8::x1,x2,y1,y2,xi,shp1,shp2
  real*8,allocatable::weigh(:), posi(:,:)
  real*8,allocatable::gradu(:,:,:)
!... 
   ngausl = 2 !...NO. of gauss point used to complete the P1P2 reconstrcution...
!...
  allocate(gradu(1:3,1:2,1:nelem))
  allocate (weigh(ngausl), posi(1,ngausl))
  call rutope(1, ngausl, posi, weigh)

!      print*,'xxxxxxxxxx'
!
!...  This sub computes the gradient using the least-square approach
!
      c00 = 0.0
!
      gradu = 0.0
!
!
!
      do 1200 ifa = 1, nbfac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
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
     rxxl = geoel(11,iel)
     ryyl = geoel(12,iel)
     rxyl = geoel(13,iel)
!
     do ig=1,ngausl

!
     xi =posi(1, ig)
 
     shp1 = 0.5d0*(1.d0-xi)
     shp2 = 0.5d0*(1.d0+xi)
!
     xg = x1*shp1 + x2*shp2
     yg = y1*shp1 + y2*shp2
  
!     xg = geofa(3,ifa)
!     yg = geofa(4,ifa)
!
!...basis function..
!
     b2l = (xg-xcl)/dxl
     b3l = (yg-ycl)/dyl
     b4l =  0.5d0*b2l**2-geoel(11,iel)
     b5l =  0.5d0*b3l**2-geoel(12,iel) 
     b6l =   b2l*b3l - geoel(13,iel)
!
!
     pl  = unkno(1,2,iel) 
     pxl = unkno(2,2,iel)     
     pyl = unkno(3,2,iel)   
   
     ql  = unkno(1,3,iel) 
     qxl = unkno(2,3,iel)     
     qyl = unkno(3,3,iel)           
!
     pr  = (cosh(pi*xg)*pi*sin(pi*yg)+sinh(pi*yg)*pi*cos(pi*xg))/sinh(pi)
     pxr = (sinh(pi*xg)*pi**2*sin(pi*yg)-sinh(pi*yg)*pi**2*sin(pi*xg))/sinh(pi)
     pyr = (cosh(pi*xg)*pi**2*cos(pi*yg)+cosh(pi*yg)*pi**2*cos(pi*xg))/sinh(pi)  
!
     qr  = (sinh(pi*xg)*pi*cos(pi*yg)+pi*cosh(pi*yg)*sin(pi*xg))/sinh(pi)
     qxr = pyr    
     qyr = (-sinh(pi*xg)*pi**2*sin(pi*yg)+sinh(pi*yg)*pi**2*sin(pi*xg))/sinh(pi)   
!
     rl11 = pr-(pl + pxl*b2l + pyl*b3l) 
     rl12 = qr-(ql + qxl*b2l + qyl*b3l)     
!
     rl21 = pxr*dxl-pxl
     rl22 = qxr*dxl-qxl
!
     rl31 = pyr*dyl-pyl
     rl32 = qyr*dyl-qyl
!
!
!     gradu(1,1,iel) =  gradu(1,1,iel) + rl11*b4l          
!     gradu(2,1,iel) =  gradu(2,1,iel) + rl11*b5l          
!     gradu(3,1,iel) =  gradu(3,1,iel) + rl11*b6l  

!     gradu(1,2,iel) =  gradu(1,2,iel) + rl12*b4l  
!     gradu(2,2,iel) =  gradu(2,2,iel) + rl12*b5l    
!     gradu(3,2,iel) =  gradu(3,2,iel) + rl12*b6l  
!
     gradu(1,1,iel) =  gradu(1,1,iel) + rl11*b4l + rl21*b2l         
     gradu(2,1,iel) =  gradu(2,1,iel) + rl11*b5l + rl31*b3l         
     gradu(3,1,iel) =  gradu(3,1,iel) + rl11*b6l + rl21*b3l + rl31*b2l  
!
     gradu(1,2,iel) =  gradu(1,2,iel) + rl12*b4l + rl22*b2l         
     gradu(2,2,iel) =  gradu(2,2,iel) + rl12*b5l + rl32*b3l         
     gradu(3,2,iel) =  gradu(3,2,iel) + rl12*b6l + rl22*b3l + rl32*b2l   

     enddo
!
!...  weighting for this edge
!
 1200 continue

!
      do 1400 ifa = nbfac+1, nafac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
      ier   = intfac(2,ifa)
!
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
     rxxl = geoel(11,iel)
     ryyl = geoel(12,iel)
     rxyl = geoel(13,iel)
!
     xcr = geoel(1,ier)
     ycr = geoel(2,ier)
     dxr = geoel(4,ier)*0.5d0
     dyr = geoel(5,ier)*0.5d0
     rxxr = geoel(11,ier)
     ryyr = geoel(12,ier)
     rxyr = geoel(13,ier)
!...basis function..
!
     b2l = (xcr-xcl)/dxl
     b3l = (ycr-ycl)/dyl
     b4l = (dxr/dxl)**2*rxxr + 0.5d0*b2l**2-rxxl
     b5l = (dyr/dyl)**2*ryyr + 0.5d0*b3l**2-ryyl 
     b6l = (dxr*dyr)/(dxl*dyl)*rxyr + b2l*b3l - rxyl
!
     b2r = (xcl-xcr)/dxr
     b3r = (ycl-ycr)/dyr
     b4r = (dxl/dxr)**2*rxxl + 0.5d0*b2r**2 - rxxr
     b5r = (dyl/dyr)**2*ryyl + 0.5d0*b3r**2 - ryyr 
     b6r = (dxl*dyl)/(dxr*dyr)*rxyl + b2r*b3r - rxyr
!
     pl  = unkno(1,2,iel) 
     pxl = unkno(2,2,iel)     
     pyl = unkno(3,2,iel)      
     ql  = unkno(1,3,iel) 
     qxl = unkno(2,3,iel)     
     qyl = unkno(3,3,iel)           
!
     pr  = unkno(1,2,ier) 
     pxr = unkno(2,2,ier)     
     pyr = unkno(3,2,ier)      
     qr  = unkno(1,3,ier) 
     qxr = unkno(2,3,ier)     
     qyr = unkno(3,3,ier)   
!
     rl11 = pr-(pl + pxl*b2l + pyl*b3l) 
     rl12 = qr-(ql + qxl*b2l + qyl*b3l)     
!
     rr11 = pl-(pr + pxr*b2r + pyr*b3r)
     rr12 = ql-(qr + qxr*b2r + qyr*b3r)
!
     rl21 = pxr*dxl/dxr-pxl
     rl22 = qxr*dxl/dxr-qxl
!
     rr21 = pxl*dxr/dxl-pxr
     rr22 = qxl*dxr/dxl-qxr
!
     rl31 = pyr*dyl/dyr-pyl
     rl32 = qyr*dyl/dyr-qyl
!
     rr31 = pyl*dyr/dyl-pyr
     rr32 = qyl*dyr/dyl-qyr
!
     gradu(1,1,iel) =  gradu(1,1,iel) + rl11*b4l + rl21*b2l         
     gradu(2,1,iel) =  gradu(2,1,iel) + rl11*b5l + rl31*b3l         
     gradu(3,1,iel) =  gradu(3,1,iel) + rl11*b6l + rl21*b3l + rl31*b2l  
!
     gradu(1,2,iel) =  gradu(1,2,iel) + rl12*b4l + rl22*b2l         
     gradu(2,2,iel) =  gradu(2,2,iel) + rl12*b5l + rl32*b3l         
     gradu(3,2,iel) =  gradu(3,2,iel) + rl12*b6l + rl22*b3l + rl32*b2l      
!
!
     gradu(1,1,ier) =  gradu(1,1,ier) + rr11*b4r + rr21*b2r         
     gradu(2,1,ier) =  gradu(2,1,ier) + rr11*b5r + rr31*b3r         
     gradu(3,1,ier) =  gradu(3,1,ier) + rr11*b6r + rr21*b3r + rr31*b2r  
!
     gradu(1,2,ier) =  gradu(1,2,ier) + rr12*b4r + rr22*b2r         
     gradu(2,2,ier) =  gradu(2,2,ier) + rr12*b5r + rr32*b3r         
     gradu(3,2,ier) =  gradu(3,2,ier) + rr12*b6r + rr22*b3r + rr32*b2r         
!
!...  weighting for this edge
!
 1400 continue
!
!print*,'du1we', dmets(1:3,1)
! print*,'gradu', gradu(1,1:3,1), gradu(2,1:3,1)
!
      do 2000 ie  = 1, nelem
!
    rhs11 = gradu(1,1,ie)
    rhs21 = gradu(2,1,ie)
    rhs31 = gradu(3,1,ie)
!
    rhs12 = gradu(1,2,ie)
    rhs22 = gradu(2,2,ie)
    rhs32 = gradu(3,2,ie)
!
    unkno(4,2,ie) = invmat(1,ie)*rhs11 + invmat(2,ie)*rhs21 + invmat(3,ie)*rhs31 
    unkno(5,2,ie) = invmat(2,ie)*rhs11 + invmat(4,ie)*rhs21 + invmat(5,ie)*rhs31 
    unkno(6,2,ie) = invmat(3,ie)*rhs11 + invmat(5,ie)*rhs21 + invmat(6,ie)*rhs31 
!
    unkno(4,3,ie) = invmat(1,ie)*rhs12 + invmat(2,ie)*rhs22 + invmat(3,ie)*rhs32
    unkno(5,3,ie) = invmat(2,ie)*rhs12 + invmat(4,ie)*rhs22 + invmat(5,ie)*rhs32
    unkno(6,3,ie) = invmat(3,ie)*rhs12 + invmat(5,ie)*rhs22 + invmat(6,ie)*rhs32

!
!    unkno(4,3,ie) = unkno(6,2,ie)
!    unkno(5,3,ie) = invmat(2,ie)*rhs12 + invmat(4,ie)*rhs22 + invmat(5,ie)*rhs32
!    unkno(6,3,ie) = unkno(5,2,ie)
!
!
 2000 continue
!
  !   print*,'gradient'
  !    print*,unkno(4:6,2,128),invmat(1:6,128)
  deallocate(gradu)
      return
      
      end subroutine getgradulsp1p2
!
!
!
 subroutine getlsinvmat(intfac,geoel ,geofa,coord,invmat)
 use constant
 implicit none 
  real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel 
  real*8,dimension(1:7,1:nafac),intent(in)::geofa
  real*8,dimension(1:3,1:2,1:nelem)::gradu
  real*8,dimension(6,nelem)::invmat
  real*8,dimension(1:ndimn,1:npoin)::coord 
  integer*4,dimension(1:8,1:nafac)::intfac
  integer*4::iel,ier,ifa,ncs,ie,ig
  real*8::c00
!
  real*8::xcl,ycl,dxl,dyl,rxxl,ryyl,rxyl
  real*8::xcr,ycr,dxr,dyr,rxxr,ryyr,rxyr
  real*8::b2l, b3l, b4l, b5l, b6l 
  real*8::b2r, b3r, b4r, b5r, b6r
  real*8::pl,pxl,pyl,ql,qxl,qyl
  real*8::pr,pxr,pyr,qr,qxr,qyr,det,xg,yg
  real*8::rl11, rl12, rl21, rl22, rl31, rl32
  real*8::rr11, rr12, rr21, rr22, rr31, rr32
  real*8::rhs11, rhs21, rhs31, rhs12, rhs22, rhs32
  real*8::a(3,3), c(3,3),b55(3)
  real*8,allocatable::weigh(:), posi(:,:)
  integer::ip1,ip2,ngausl
  real*8::x1,x2,y1,y2,xi,shp1,shp2
!....
   ngausl = 2 !...NO. of gauss point used to complete the P1P2 reconstrcution...
!...
  allocate (weigh(ngausl), posi(1,ngausl))
 call rutope(1, ngausl, posi, weigh)

!     print*,'xxxxxxxxxx'
!
!...  This sub computes the gradient using the least-square approach
!
      c00 = 0.0
!
      invmat = 0.0
!
     do 1300 ifa = 1, nbfac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
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
     rxxl = geoel(11,iel)
     ryyl = geoel(12,iel)
     rxyl = geoel(13,iel)
!
     do ig=1,ngausl

!
     xi =posi(1, ig)
 
     shp1 = 0.5d0*(1.d0-xi)
     shp2 = 0.5d0*(1.d0+xi)
 
     xg = x1*shp1 + x2*shp2
     yg = y1*shp1 + y2*shp2
!
!...basis function..
!
     b2l = (xg-xcl)/dxl
     b3l = (yg-ycl)/dyl
     b4l =  0.5d0*b2l**2-geoel(11,iel)
     b5l =  0.5d0*b3l**2-geoel(12,iel) 
     b6l =   b2l*b3l - geoel(13,iel)
!
!
!     invmat(1,iel) = invmat(1,iel) + b4l*b4l 
!     invmat(2,iel) = invmat(2,iel) + b4l*b5l 
!     invmat(3,iel) = invmat(3,iel) + b4l*b6l
!     invmat(4,iel) = invmat(4,iel) + b5l*b5l 
!     invmat(5,iel) = invmat(5,iel) + b5l*b6l 
!     invmat(6,iel) = invmat(6,iel) + b6l*b6l 

     invmat(1,iel) = invmat(1,iel) + b4l*b4l + b2l*b2l
     invmat(2,iel) = invmat(2,iel) + b4l*b5l 
     invmat(3,iel) = invmat(3,iel) + b4l*b6l + b2l*b3l
     invmat(4,iel) = invmat(4,iel) + b5l*b5l + b3l*b3l
     invmat(5,iel) = invmat(5,iel) + b5l*b6l + b2l*b3l
     invmat(6,iel) = invmat(6,iel) + b6l*b6l + b2l*b2l + b3l*b3l 
     enddo
!
!...  weighting for this edge
!
 1300 continue

!
      do 1400 ifa = nbfac+1, nafac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
      ier   = intfac(2,ifa)
!
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
     rxxl = geoel(11,iel)
     ryyl = geoel(12,iel)
     rxyl = geoel(13,iel)
!
     xcr = geoel(1,ier)
     ycr = geoel(2,ier)
     dxr = geoel(4,ier)*0.5d0
     dyr = geoel(5,ier)*0.5d0
     rxxr = geoel(11,ier)
     ryyr = geoel(12,ier)
     rxyr = geoel(13,ier)
!...basis function..
!
     b2l = (xcr-xcl)/dxl
     b3l = (ycr-ycl)/dyl
     b4l = (dxr/dxl)**2*rxxr+0.5d0*b2l**2-rxxl
     b5l = (dyr/dyl)**2*ryyr+0.5d0*b3l**2-ryyl 
     b6l = (dxr*dyr)/(dxl*dyl)*rxyr+b2l*b3l-rxyl
!
     b2r = (xcl-xcr)/dxr
     b3r = (ycl-ycr)/dyr
     b4r = (dxl/dxr)**2*rxxl+0.5d0*b2r**2-rxxr
     b5r = (dyl/dyr)**2*ryyl+0.5d0*b3r**2-ryyr 
     b6r = (dxl*dyl)/(dxr*dyr)*rxyl+b2r*b3r-rxyr
!
!
     invmat(1,iel) = invmat(1,iel) + b4l*b4l + b2l*b2l
     invmat(2,iel) = invmat(2,iel) + b4l*b5l 
     invmat(3,iel) = invmat(3,iel) + b4l*b6l + b2l*b3l
     invmat(4,iel) = invmat(4,iel) + b5l*b5l + b3l*b3l
     invmat(5,iel) = invmat(5,iel) + b5l*b6l + b2l*b3l
     invmat(6,iel) = invmat(6,iel) + b6l*b6l + b2l*b2l + b3l*b3l
!
     invmat(1,ier) = invmat(1,ier) + b4r*b4r + b2r*b2r
     invmat(2,ier) = invmat(2,ier) + b4r*b5r 
     invmat(3,ier) = invmat(3,ier) + b4r*b6r + b2r*b3r
     invmat(4,ier) = invmat(4,ier) + b5r*b5r + b3r*b3r
     invmat(5,ier) = invmat(5,ier) + b5r*b6r + b2r*b3r
     invmat(6,ier) = invmat(6,ier) + b6r*b6r + b2r*b2r + b3r*b3r
!
 
!
!...  weighting for this edge
!
 1400 continue
!
! print*,'gradu', gradu(1,1:3,1), gradu(2,1:3,1)
!
    do 2000 ie  = 1, nelem
!
    a= 0.d0
    c= 0.d0
!
     a(1,1) = invmat(1,ie)
     a(1,2) = invmat(2,ie)
     a(1,3) = invmat(3,ie)

     a(2,1) = a(1,2)
     a(2,2) = invmat(4,ie)
     a(2,3) = invmat(5,ie)

     a(3,1) = a(1,3)
     a(3,2) = a(2,3)
     a(3,3) = invmat(6,ie)
!
!     call inverse(a,c,3,det)
!
     call getinvmat(3, a, c, b55)
     invmat(1,ie) = c(1,1)
     invmat(2,ie) = c(1,2) 
     invmat(3,ie) = c(1,3) 
     invmat(4,ie) = c(2,2) 
     invmat(5,ie) = c(2,3)
     invmat(6,ie) = c(3,3) 
!
 2000 continue
!
!     print*,'geoel57'
!      print*,invmat(1:6,92),invmat(1:6,110),invmat(1:6,128)
      return
      
  end subroutine getlsinvmat
!
!...WENOP2 for DG(P0P3) + rDG(P1P2) on skewed mesh...
!
      subroutine wenop2(ltelem,geoel,unkno)
!
      use constant
      implicit real*8 (a-h,o-z)
!
      integer ltelem(3,nelem)
!
      real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel
      real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac)::unkno
!
!...  local arrays
!
      parameter (nsten = 4)
!
      real*8  weigh(2,nsten)
      real*8  gradv(3,2,nsten)
      real*8  os(2,nsten)
      real*8, allocatable :: gradu(:,:,:)
!
      data c00   / 0.0d0    /
      data c05   / 0.5d0    /
      data c10   / 1.0d0    /
      data c20   / 2.0d0    /
      data epsil / 1.0d-6   /   
!     data rpowe / 5.0d0    /   ! works for lilia case 
!     data rpowe / 2.0d0    /   ! works for lilia case 
      data rpowe / 2.0d0    /   ! works for lilia case 
!
      allocate(gradu(3,2,nelem))
!
!...  This sub reconstructs the second derivatives 
!     using a weno reconstruction
!
!...  local integer constant 
!
      iq1 = 2
      iq2 = 3
!
!...  loop over the hexahedral cells
!
      do ie = 1, nelem
!
      ielem = ie
!
!...  get the shape functions and the jacobian of this element
!
      xc    = geoel(1,ielem)
      yc    = geoel(2,ielem)
      volel = geoel(3,ielem)
      dxm   = c20/geoel(4,ielem)  !1/(0.5*dx)
      dym   = c20/geoel(5,ielem)  !1/(0.5*dy)
!
!...  a. curvatures for this cell itself
!
      isten = 1
      gradv(1,1:2,isten) = unkno( 4,iq1:iq2,ie)*dxm*dxm
      gradv(2,1:2,isten) = unkno( 5,iq1:iq2,ie)*dym*dym
      gradv(3,1:2,isten) = unkno( 6,iq1:iq2,ie)*dxm*dym
!
!...  b. curvatures for the face-neighboring cells
!
      do ies = 1, 3
        !
        je = ltelem(ies,ie)
        !
        if(je .le. nelem) then
          !
          isten = isten+1
          !
          dxm = c20/geoel(4,je)
          dym = c20/geoel(5,je)
          !
          gradv(1,1:2,isten) = unkno( 4,iq1:iq2,je)*dxm*dxm
          gradv(2,1:2,isten) = unkno( 5,iq1:iq2,je)*dym*dym
          gradv(3,1:2,isten) = unkno( 6,iq1:iq2,je)*dxm*dym

          !
        endif
        !
      enddo
!
!...  2. compute the oscillation indicators
!    
      weig1 = c00
      weig2 = c00
!
      const = c10
!
      do is = 1, isten
        !
        do iq = 1, 2
          !
          os(iq,is) = gradv(1,iq,is)*gradv(1,iq,is) &
                    + gradv(2,iq,is)*gradv(2,iq,is) &
                    + 2.d0*gradv(3,iq,is)*gradv(3,iq,is) 
          !
        enddo
        !
      enddo
!
      do is = 1, isten
        !
        if(is .eq. 1) then
          !
          wi = 20.d0
          !
        else
          !
          wi = c10
          !
        endif
        !
        weigh(1,is) = wi/(epsil+sqrt(os(1,is)*const))**rpowe
        weigh(2,is) = wi/(epsil+sqrt(os(2,is)*const))**rpowe
        !
        weig1 = weig1 + weigh(1,is)
        weig2 = weig2 + weigh(2,is)
        !
      enddo
!
      do is = 1, isten
        !
        weigh(1,is) = weigh(1,is)/weig1
        weigh(2,is) = weigh(2,is)/weig2
        !
      enddo
!
!...  3. finally, we can reconstruct the curvatures
!
      gradu(1:3,1:2,ie) = c00
!
      do is = 1, isten
!
      gradu(1,1:2,ie) = gradu(1,1:2,ie) + weigh(1:2,is)*gradv(1,1:2,is)
      gradu(2,1:2,ie) = gradu(2,1:2,ie) + weigh(1:2,is)*gradv(2,1:2,is)
      gradu(3,1:2,ie) = gradu(3,1:2,ie) + weigh(1:2,is)*gradv(3,1:2,is)
!
      enddo
!
!...  end of the loop over the hexahedral cells
!
      enddo
!
!...  modify the curvatures
!
      do ie  = 1, nelem
        !
        ielem = ie
        !
        dx = geoel(4,ielem)/c20
        dy = geoel(5,ielem)/c20
        !
        unkno( 4,iq1:iq2,ielem) = gradu(1,1:2,ie)*dx*dx
        unkno( 5,iq1:iq2,ielem) = gradu(2,1:2,ie)*dy*dy
        unkno( 6,iq1:iq2,ielem) = gradu(3,1:2,ie)*dx*dy
        !
      enddo
!
      deallocate(gradu)
!
      return
      end
!
!...Reconstruction P1P2 for BR2 phi
!
 subroutine getgradulsp1p2br2(intfac,geoel ,geofa, unkno ,coord,invmat)
  use constant
  implicit none 
 real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel 
 real*8,dimension(1:7,1:nafac),intent(in)::geofa
 real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac)::unkno
 real*8,dimension(6,nelem),intent(in)::invmat
 real*8,dimension(1:ndimn,1:npoin)::coord 
 integer*4,dimension(1:8,1:nafac)::intfac
 integer*4::iel,ier,ifa,ncs,ie,ngausl,ig
 real*8::c00
!
  real*8::xcl,ycl,dxl,dyl,rxxl,ryyl,rxyl
  real*8::xcr,ycr,dxr,dyr,rxxr,ryyr,rxyr
  real*8::b2l, b3l, b4l, b5l, b6l 
  real*8::b2r, b3r, b4r, b5r, b6r
  real*8::pl,pxl,pyl,ql,qxl,qyl,xg,yg
  real*8::pr,pxr,pyr,qr,qxr,qyr
  real*8::rl11, rl12, rl21, rl22, rl31, rl32
  real*8::rr11, rr12, rr21, rr22, rr31, rr32
  real*8::rhs11, rhs21, rhs31, rhs12, rhs22, rhs32
  real*8::x1,x2,y1,y2,xi,shp1,shp2
  real*8,allocatable::weigh(:), posi(:,:)
  real*8,allocatable::gradu(:,:,:)
!... 
   ngausl = 2 !...NO. of gauss point used to complete the P1P2 reconstrcution...
!...
  allocate(gradu(1:3,1,1:nelem))
  allocate (weigh(ngausl), posi(1,ngausl))
  call rutope(1, ngausl, posi, weigh)

!      print*,'xxxxxxxxxx'
!
!...  This sub computes the gradient using the least-square approach
!
      c00 = 0.0
!
      gradu = 0.0
!
!
!
      do 1200 ifa = 1, nbfac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
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
     rxxl = geoel(11,iel)
     ryyl = geoel(12,iel)
     rxyl = geoel(13,iel)
!
     do ig=1,ngausl

!
     xi =posi(1, ig)
 
     shp1 = 0.5d0*(1.d0-xi)
     shp2 = 0.5d0*(1.d0+xi)
!
     xg = x1*shp1 + x2*shp2
     yg = y1*shp1 + y2*shp2
  
!     xg = geofa(3,ifa)
!     yg = geofa(4,ifa)
!
!...basis function..
!
     b2l = (xg-xcl)/dxl
     b3l = (yg-ycl)/dyl
     b4l =  0.5d0*b2l**2-geoel(11,iel)
     b5l =  0.5d0*b3l**2-geoel(12,iel) 
     b6l =   b2l*b3l - geoel(13,iel)
!
!
     pl  = unkno(1,1,iel) 
     pxl = unkno(2,1,iel)     
     pyl = unkno(3,1,iel)           
!
     pr  = (sinh(pi*xg)*sin(pi*yg)+sinh(pi*yg)*sin(pi*xg))/sinh(pi)
     pxr = (cosh(pi*xg)*pi*sin(pi*yg)+sinh(pi*yg)*pi*cos(pi*xg))/sinh(pi)  
     pyr = (sinh(pi*xg)*pi*cos(pi*yg)+pi*cosh(pi*yg)*sin(pi*xg))/sinh(pi)
!
!
     rl11 = pr-(pl + pxl*b2l + pyl*b3l) 
!
     rl21 = pxr*dxl-pxl
!
     rl31 = pyr*dyl-pyl
!
!
!     gradu(1,1,iel) =  gradu(1,1,iel) + rl11*b4l          
!     gradu(2,1,iel) =  gradu(2,1,iel) + rl11*b5l          
!     gradu(3,1,iel) =  gradu(3,1,iel) + rl11*b6l  

!     gradu(1,2,iel) =  gradu(1,2,iel) + rl12*b4l  
!     gradu(2,2,iel) =  gradu(2,2,iel) + rl12*b5l    
!     gradu(3,2,iel) =  gradu(3,2,iel) + rl12*b6l  
!
     gradu(1,1,iel) =  gradu(1,1,iel) + rl11*b4l + rl21*b2l         
     gradu(2,1,iel) =  gradu(2,1,iel) + rl11*b5l + rl31*b3l         
     gradu(3,1,iel) =  gradu(3,1,iel) + rl11*b6l + rl21*b3l + rl31*b2l  
     enddo
!
!...  weighting for this edge
!
 1200 continue

!
      do 1400 ifa = nbfac+1, nafac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
      ier   = intfac(2,ifa)
!
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
     rxxl = geoel(11,iel)
     ryyl = geoel(12,iel)
     rxyl = geoel(13,iel)
!
     xcr = geoel(1,ier)
     ycr = geoel(2,ier)
     dxr = geoel(4,ier)*0.5d0
     dyr = geoel(5,ier)*0.5d0
     rxxr = geoel(11,ier)
     ryyr = geoel(12,ier)
     rxyr = geoel(13,ier)
!...basis function..
!
     b2l = (xcr-xcl)/dxl
     b3l = (ycr-ycl)/dyl
     b4l = (dxr/dxl)**2*rxxr + 0.5d0*b2l**2-rxxl
     b5l = (dyr/dyl)**2*ryyr + 0.5d0*b3l**2-ryyl 
     b6l = (dxr*dyr)/(dxl*dyl)*rxyr + b2l*b3l - rxyl
!
     b2r = (xcl-xcr)/dxr
     b3r = (ycl-ycr)/dyr
     b4r = (dxl/dxr)**2*rxxl + 0.5d0*b2r**2 - rxxr
     b5r = (dyl/dyr)**2*ryyl + 0.5d0*b3r**2 - ryyr 
     b6r = (dxl*dyl)/(dxr*dyr)*rxyl + b2r*b3r - rxyr
!
     pl  = unkno(1,1,iel) 
     pxl = unkno(2,1,iel)     
     pyl = unkno(3,1,iel)           
!
     pr  = unkno(1,1,ier) 
     pxr = unkno(2,1,ier)     
     pyr = unkno(3,1,ier)      
!
     rl11 = pr-(pl + pxl*b2l + pyl*b3l) 
!
     rr11 = pl-(pr + pxr*b2r + pyr*b3r)
!
     rl21 = pxr*dxl/dxr-pxl
!
     rr21 = pxl*dxr/dxl-pxr
!
     rl31 = pyr*dyl/dyr-pyl
!
     rr31 = pyl*dyr/dyl-pyr
!
     gradu(1,1,iel) =  gradu(1,1,iel) + rl11*b4l + rl21*b2l         
     gradu(2,1,iel) =  gradu(2,1,iel) + rl11*b5l + rl31*b3l         
     gradu(3,1,iel) =  gradu(3,1,iel) + rl11*b6l + rl21*b3l + rl31*b2l  
!
     gradu(1,1,ier) =  gradu(1,1,ier) + rr11*b4r + rr21*b2r         
     gradu(2,1,ier) =  gradu(2,1,ier) + rr11*b5r + rr31*b3r         
     gradu(3,1,ier) =  gradu(3,1,ier) + rr11*b6r + rr21*b3r + rr31*b2r       
!
!...  weighting for this edge
!
 1400 continue
!
!print*,'du1we', dmets(1:3,1)
! print*,'gradu', gradu(1,1:3,1), gradu(2,1:3,1)
!
      do 2000 ie  = 1, nelem
!
    rhs11 = gradu(1,1,ie)
    rhs21 = gradu(2,1,ie)
    rhs31 = gradu(3,1,ie)
!
!
    unkno(4,1,ie) = invmat(1,ie)*rhs11 + invmat(2,ie)*rhs21 + invmat(3,ie)*rhs31 
    unkno(5,1,ie) = invmat(2,ie)*rhs11 + invmat(4,ie)*rhs21 + invmat(5,ie)*rhs31 
    unkno(6,1,ie) = invmat(3,ie)*rhs11 + invmat(5,ie)*rhs21 + invmat(6,ie)*rhs31 
!
!
!    unkno(4,3,ie) = unkno(6,2,ie)
!    unkno(5,3,ie) = invmat(2,ie)*rhs12 + invmat(4,ie)*rhs22 + invmat(5,ie)*rhs32
!    unkno(6,3,ie) = unkno(5,2,ie)
!
!
 2000 continue
!
  !   print*,'gradient'
  !    print*,unkno(4:6,2,128),invmat(1:6,128)
  deallocate(gradu)
      return
      
      end subroutine getgradulsp1p2br2
!
!
!
 subroutine getlsinvmatbr2(intfac,geoel ,geofa,coord,invmat)
 use constant
 implicit none 
  real*8,dimension(1:ngeel,1:nelem+nbfac)::geoel 
  real*8,dimension(1:7,1:nafac),intent(in)::geofa
  real*8,dimension(1:3,1:2,1:nelem)::gradu
  real*8,dimension(6,nelem)::invmat
  real*8,dimension(1:ndimn,1:npoin)::coord 
  integer*4,dimension(1:8,1:nafac)::intfac
  integer*4::iel,ier,ifa,ncs,ie,ig
  real*8::c00
!
  real*8::xcl,ycl,dxl,dyl,rxxl,ryyl,rxyl
  real*8::xcr,ycr,dxr,dyr,rxxr,ryyr,rxyr
  real*8::b2l, b3l, b4l, b5l, b6l 
  real*8::b2r, b3r, b4r, b5r, b6r
  real*8::pl,pxl,pyl,ql,qxl,qyl
  real*8::pr,pxr,pyr,qr,qxr,qyr,det,xg,yg
  real*8::rl11, rl12, rl21, rl22, rl31, rl32
  real*8::rr11, rr12, rr21, rr22, rr31, rr32
  real*8::rhs11, rhs21, rhs31, rhs12, rhs22, rhs32
  real*8::a(3,3), c(3,3),b55(3)
  real*8,allocatable::weigh(:), posi(:,:)
  integer::ip1,ip2,ngausl
  real*8::x1,x2,y1,y2,xi,shp1,shp2
!....
   ngausl = 2 !...NO. of gauss point used to complete the P1P2 reconstrcution...
!...
  allocate (weigh(ngausl), posi(1,ngausl))
 call rutope(1, ngausl, posi, weigh)

!     print*,'xxxxxxxxxx'
!
!...  This sub computes the gradient using the least-square approach
!
      c00 = 0.0
!
      invmat = 0.0
!
     do 1300 ifa = 1, nbfac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
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
     rxxl = geoel(11,iel)
     ryyl = geoel(12,iel)
     rxyl = geoel(13,iel)
!
     do ig=1,ngausl

!
     xi =posi(1, ig)
 
     shp1 = 0.5d0*(1.d0-xi)
     shp2 = 0.5d0*(1.d0+xi)
 
     xg = x1*shp1 + x2*shp2
     yg = y1*shp1 + y2*shp2
!
!...basis function..
!
     b2l = (xg-xcl)/dxl
     b3l = (yg-ycl)/dyl
     b4l =  0.5d0*b2l**2-geoel(11,iel)
     b5l =  0.5d0*b3l**2-geoel(12,iel) 
     b6l =   b2l*b3l - geoel(13,iel)
!
!
!     invmat(1,iel) = invmat(1,iel) + b4l*b4l 
!     invmat(2,iel) = invmat(2,iel) + b4l*b5l 
!     invmat(3,iel) = invmat(3,iel) + b4l*b6l
!     invmat(4,iel) = invmat(4,iel) + b5l*b5l 
!     invmat(5,iel) = invmat(5,iel) + b5l*b6l 
!     invmat(6,iel) = invmat(6,iel) + b6l*b6l 

     invmat(1,iel) = invmat(1,iel) + b4l*b4l + b2l*b2l
     invmat(2,iel) = invmat(2,iel) + b4l*b5l 
     invmat(3,iel) = invmat(3,iel) + b4l*b6l + b2l*b3l
     invmat(4,iel) = invmat(4,iel) + b5l*b5l + b3l*b3l
     invmat(5,iel) = invmat(5,iel) + b5l*b6l + b2l*b3l
     invmat(6,iel) = invmat(6,iel) + b6l*b6l + b2l*b2l + b3l*b3l 
     enddo
!
!...  weighting for this edge
!
 1300 continue

!
      do 1400 ifa = nbfac+1, nafac
!
!...  end-elements of the faces
!
      iel   = intfac(1,ifa)
      ier   = intfac(2,ifa)
!
     xcl = geoel(1,iel)
     ycl = geoel(2,iel)
     dxl = geoel(4,iel)*0.5d0
     dyl = geoel(5,iel)*0.5d0
     rxxl = geoel(11,iel)
     ryyl = geoel(12,iel)
     rxyl = geoel(13,iel)
!
     xcr = geoel(1,ier)
     ycr = geoel(2,ier)
     dxr = geoel(4,ier)*0.5d0
     dyr = geoel(5,ier)*0.5d0
     rxxr = geoel(11,ier)
     ryyr = geoel(12,ier)
     rxyr = geoel(13,ier)
!...basis function..
!
     b2l = (xcr-xcl)/dxl
     b3l = (ycr-ycl)/dyl
     b4l = (dxr/dxl)**2*rxxr+0.5d0*b2l**2-rxxl
     b5l = (dyr/dyl)**2*ryyr+0.5d0*b3l**2-ryyl 
     b6l = (dxr*dyr)/(dxl*dyl)*rxyr+b2l*b3l-rxyl
!
     b2r = (xcl-xcr)/dxr
     b3r = (ycl-ycr)/dyr
     b4r = (dxl/dxr)**2*rxxl+0.5d0*b2r**2-rxxr
     b5r = (dyl/dyr)**2*ryyl+0.5d0*b3r**2-ryyr 
     b6r = (dxl*dyl)/(dxr*dyr)*rxyl+b2r*b3r-rxyr
!
!
     invmat(1,iel) = invmat(1,iel) + b4l*b4l + b2l*b2l
     invmat(2,iel) = invmat(2,iel) + b4l*b5l 
     invmat(3,iel) = invmat(3,iel) + b4l*b6l + b2l*b3l
     invmat(4,iel) = invmat(4,iel) + b5l*b5l + b3l*b3l
     invmat(5,iel) = invmat(5,iel) + b5l*b6l + b2l*b3l
     invmat(6,iel) = invmat(6,iel) + b6l*b6l + b2l*b2l + b3l*b3l
!
     invmat(1,ier) = invmat(1,ier) + b4r*b4r + b2r*b2r
     invmat(2,ier) = invmat(2,ier) + b4r*b5r 
     invmat(3,ier) = invmat(3,ier) + b4r*b6r + b2r*b3r
     invmat(4,ier) = invmat(4,ier) + b5r*b5r + b3r*b3r
     invmat(5,ier) = invmat(5,ier) + b5r*b6r + b2r*b3r
     invmat(6,ier) = invmat(6,ier) + b6r*b6r + b2r*b2r + b3r*b3r
!
 
!
!...  weighting for this edge
!
 1400 continue
!
! print*,'gradu', gradu(1,1:3,1), gradu(2,1:3,1)
!
    do 2000 ie  = 1, nelem
!
    a= 0.d0
    c= 0.d0
!
     a(1,1) = invmat(1,ie)
     a(1,2) = invmat(2,ie)
     a(1,3) = invmat(3,ie)

     a(2,1) = a(1,2)
     a(2,2) = invmat(4,ie)
     a(2,3) = invmat(5,ie)

     a(3,1) = a(1,3)
     a(3,2) = a(2,3)
     a(3,3) = invmat(6,ie)
!
!     call inverse(a,c,3,det)
!
     call getinvmat(3, a, c, b55)
     invmat(1,ie) = c(1,1)
     invmat(2,ie) = c(1,2) 
     invmat(3,ie) = c(1,3) 
     invmat(4,ie) = c(2,2) 
     invmat(5,ie) = c(2,3)
     invmat(6,ie) = c(3,3) 
!
 2000 continue
!
!     print*,'geoel57'
!      print*,invmat(1:6,92),invmat(1:6,110),invmat(1:6,128)
      return
      
  end subroutine getlsinvmatbr2

!
!...ls square...
!
subroutine getgraduls_euler(bface, intfac,geoel ,unkno ,coord ,geofa, uchar)
use constant
implicit none
integer*4,dimension(1:nbfai,nbfac)::bface
real*8,dimension(1:ngeel,1:nelem+nbfac), intent(in)::geoel
real*8,dimension(1:ngefa,1:nafac),       intent(in)::geofa
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),    intent(in)::coord
integer*4,dimension(1:nifai,1:nafac), intent(in)::intfac
real*8::uchar(1:nq)
!
real*8,dimension(1:ndimn,1:nq,1:nelem)::gradu
!
!...  local arrays
!
integer*4::ie,ifa,ncs,iq, iel, ier
real*8::c00
real*8::dmets(1:3,1:nelem+nbfac)
!
real*8::dx,dy,weigh
real*8::dxwei,dywei
real*8::du1we,du2we,du3we,du4we
real*8::rxx,ryy,rxy
real*8::rco1,rco2,rco3,rco4
real*8::ru1x,ru2x,ru3x,ru4x
real*8::ru1y,ru2y,ru3y,ru4y
real*8::rjac1
real*8::x1,y1,x2,y2,x0,y0
real*8::la,lb,lc,abp2,abm2,xg,yg
integer::ip1,ip2
!
real*8::rho1,rhou1,rhov1,rhoe1,rho2,rhou2,rhov2,rhoe2
real*8::vn,uadv2,vadv2,rnx,rny, rhom1,uadv1, vadv1
!
!...  This sub computes the gradient using the least-square approach
!
c00 = 0.0
dmets = 0.0
ncs=0
!
gradu = 0.0
!
unkno(2:3, :, :) = 0.d0
!
do 1400 ifa = 1, nafac
!
!...  end-elements of the faces
!
iel   = intfac(1,ifa)
ier   = intfac(2,ifa)

if(ifa.le.nbfac) then
dx=geofa(3,ifa)-geoel(1,iel)
dy=geofa(4,ifa)-geoel(2,iel)
!
if(bface(3,ifa).eq.2)then
!   print*,'sym ifa',ifa,intfac(3:4,ifa)
!...Inviscid boundary wall
rho1  = unkno(1, 1, iel)
rhou1 = unkno(1, 2, iel)
rhov1 = unkno(1, 3, iel)
rhoe1 = unkno(1, 4, iel)
!
!...  derived variables
rhom1 = 1.d0/rho1
uadv1 = rhou1*rhom1
vadv1 = rhov1*rhom1
!
rnx =  geofa(1, ifa)
rny =  geofa(2, ifa)
!
!...  unknows for the ghost state
rho2  = rho1
rhoe2 = rhoe1
vn    = uadv1*rnx+ vadv1*rny
uadv2 = uadv1 - vn*rnx
vadv2 = vadv1 - vn*rny
rhou2 = rho2*uadv2
rhov2 = rho2*vadv2
!
unkno(1, 1, ier) = rho2
unkno(1, 2, ier) = rhou2
unkno(1, 3, ier) = rhov2
unkno(1, 4, ier) = rhoe2
!
!...Riemann solver for farfield...
!
elseif(bface(3,ifa).eq.4)then
!   print*,'far',ifa,intfac(3:4,ifa)
!...Farfield
unkno(1, 1:nq, ier) = uchar(1:4)
!
endif
!
else
dx= geoel(1,ier) - geoel(1,iel)
dy= geoel(2,ier) - geoel(2,iel)
endif
!
!...  weighting for this edge
!
!      weigh = 1.0
weigh = 1.0/sqrt(dx*dx + dy*dy )!+ dz*dz)
!
dxwei = weigh*dx
dywei = weigh*dy
!
dmets(1,iel) = dmets(1,iel) + dxwei*dxwei
dmets(2,iel) = dmets(2,iel) + dywei*dywei
dmets(3,iel) = dmets(3,iel) + dxwei*dywei

!
dmets(1,ier) = dmets(1,ier) + dxwei*dxwei
dmets(2,ier) = dmets(2,ier) + dywei*dywei
dmets(3,ier) = dmets(3,ier) + dxwei*dywei
!
du1we = weigh*(unkno(1, ncs+1,ier) - unkno(1, ncs+1,iel))
du2we = weigh*(unkno(1, ncs+2,ier) - unkno(1, ncs+2,iel))
du3we = weigh*(unkno(1, ncs+3,ier) - unkno(1, ncs+3,iel))
du4we = weigh*(unkno(1, ncs+4,ier) - unkno(1, ncs+4,iel))
!
gradu(1,ncs+1,iel) = gradu(1,ncs+1,iel) + dxwei*du1we
gradu(1,ncs+2,iel) = gradu(1,ncs+2,iel) + dxwei*du2we
gradu(1,ncs+3,iel) = gradu(1,ncs+3,iel) + dxwei*du3we
gradu(1,ncs+4,iel) = gradu(1,ncs+4,iel) + dxwei*du4we


gradu(2,ncs+1,iel) = gradu(2,ncs+1,iel) + dywei*du1we
gradu(2,ncs+2,iel) = gradu(2,ncs+2,iel) + dywei*du2we
gradu(2,ncs+3,iel) = gradu(2,ncs+3,iel) + dywei*du3we
gradu(2,ncs+4,iel) = gradu(2,ncs+4,iel) + dywei*du4we

!  excluding the ghost element gradient
if(ifa.gt.nbfac) then
gradu(1,ncs+1,ier) = gradu(1,ncs+1,ier) + dxwei*du1we
gradu(1,ncs+2,ier) = gradu(1,ncs+2,ier) + dxwei*du2we
gradu(1,ncs+3,ier) = gradu(1,ncs+3,ier) + dxwei*du3we
gradu(1,ncs+4,ier) = gradu(1,ncs+4,ier) + dxwei*du4we

gradu(2,ncs+1,ier) = gradu(2,ncs+1,ier) + dywei*du1we
gradu(2,ncs+2,ier) = gradu(2,ncs+2,ier) + dywei*du2we
gradu(2,ncs+3,ier) = gradu(2,ncs+3,ier) + dywei*du3we
gradu(2,ncs+4,ier) = gradu(2,ncs+4,ier) + dywei*du4we
endif
!
1400 continue
!
do 2000 ie  = 1, nelem
!
rxx   = dmets(1,ie)
ryy   = dmets(2,ie)
rxy   = dmets(3,ie)
!
rco1  = ryy
rco2  = -rxy
!nafac
rco3  = -rxy
rco4  = rxx

rjac1 = 1.0/(rxx*ryy - rxy*rxy)
!
!
ru1x  = gradu(1,ncs+1,ie)
ru1y  = gradu(2,ncs+1,ie)

ru2x  = gradu(1,ncs+2,ie)
ru2y  = gradu(2,ncs+2,ie)
!
ru3x  = gradu(1,ncs+3,ie)
ru3y  = gradu(2,ncs+3,ie)
!
ru4x  = gradu(1,ncs+4,ie)
ru4y  = gradu(2,ncs+4,ie)
!
gradu(1,ncs+1,ie) = (ru1x*rco1 + ru1y*rco2 )*rjac1
gradu(2,ncs+1,ie) = (ru1x*rco3 + ru1y*rco4 )*rjac1

gradu(1,ncs+2,ie) = (ru2x*rco1 + ru2y*rco2 )*rjac1
gradu(2,ncs+2,ie) = (ru2x*rco3 + ru2y*rco4 )*rjac1

gradu(1,ncs+3,ie) = (ru3x*rco1 + ru3y*rco2 )*rjac1
gradu(2,ncs+3,ie) = (ru3x*rco3 + ru3y*rco4 )*rjac1
!
gradu(1,ncs+4,ie) = (ru4x*rco1 + ru4y*rco2 )*rjac1
gradu(2,ncs+4,ie) = (ru4x*rco3 + ru4y*rco4 )*rjac1
!     print*,'output'
2000 continue
!
do ie = 1, nelem
!
dx = geoel(4,ie)*0.5d0
dy = geoel(5,ie)*0.5d0
!
do iq =1, nq
unkno(2, iq, ie) = gradu(1, iq, ie)*dx
unkno(3, iq, ie) = gradu(2, iq, ie)*dy
enddo

enddo
!   print*,'gradient'
!    print*,gradu(1:2,ncs+1,9),gradu(1:2,ncs+1,10)
return

end subroutine getgraduls_euler
!
!
!
subroutine getgraduls_euler_primitive(bface, intfac,geoel ,unkno ,coord ,geofa, uchar, unpri)
use constant
implicit none
integer*4,dimension(1:nbfai,nbfac)::bface
real*8,dimension(1:ngeel,1:nelem+nbfac), intent(in)::geoel
real*8,dimension(1:ngefa,1:nafac),       intent(in)::geofa
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(inout)::unkno
real*8,dimension(1:ndimn,1:npoin),    intent(in)::coord
integer*4,dimension(1:nifai,1:nafac), intent(in)::intfac
real*8::uchar(1:nq)
!
real*8,dimension(1:ndimn,1:nq,1:nelem)::gradu
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac), intent(out) ::unpri
!
!...  local arrays
!
integer*4::ie,ifa,ncs,iq, iel, ier
real*8::c00
real*8::dmets(1:3,1:nelem+nbfac)
!
real*8::dx,dy,weigh
real*8::dxwei,dywei
real*8::du1we,du2we,du3we,du4we
real*8::rxx,ryy,rxy
real*8::rco1,rco2,rco3,rco4
real*8::ru1x,ru2x,ru3x,ru4x
real*8::ru1y,ru2y,ru3y,ru4y
real*8::rjac1
real*8::x1,y1,x2,y2,x0,y0
real*8::la,lb,lc,abp2,abm2,xg,yg
integer::ip1,ip2
!
real*8::rho1,rhou1,rhov1,rhoe1,rho2,rhou2,rhov2,rhoe2
real*8::vn,uadv2,vadv2,rnx,rny, rhom1,uadv1, vadv1,pres1,pres2
!
!...  This sub computes the gradient using the least-square approach
!
c00 = 0.0
dmets = 0.0
ncs=0
!
gradu = 0.0
!
unkno(2:3, :, :) = 0.d0
!
!...Transform conservative to primitive...
!
do ie =1 ,nelem
!
rho1  = unkno(1, 1, ie)
rhou1 = unkno(1, 2, ie)
rhov1 = unkno(1, 3, ie)
rhoe1 = unkno(1, 4, ie)
!
rhom1 = 1.d0/rho1
uadv1 = rhou1*rhom1
vadv1 = rhov1*rhom1
!
pres1 = (gamma-1.d0)*(rhoe1-0.5d0*rho1*(uadv1**2 + vadv1**2))
!
unpri(1, 1, ie) = rho1
unpri(1, 2, ie) = uadv1
unpri(1, 3, ie) = vadv1
unpri(1, 4, ie) = pres1
!
enddo
!
do 1400 ifa = 1, nafac
!
!...  end-elements of the faces
!
iel   = intfac(1,ifa)
ier   = intfac(2,ifa)

if(ifa.le.nbfac) then
dx=geofa(3,ifa)-geoel(1,iel)
dy=geofa(4,ifa)-geoel(2,iel)
!
if(bface(3,ifa).eq.2)then
!   print*,'sym ifa',ifa,intfac(3:4,ifa)
!...Inviscid boundary wall
rho1  = unkno(1, 1, iel)
rhou1 = unkno(1, 2, iel)
rhov1 = unkno(1, 3, iel)
rhoe1 = unkno(1, 4, iel)
!
!...  derived variables
rhom1 = 1.d0/rho1
uadv1 = rhou1*rhom1
vadv1 = rhov1*rhom1
!
rnx =  geofa(1, ifa)
rny =  geofa(2, ifa)
!
!...  unknows for the ghost state
rho2  = rho1
rhoe2 = rhoe1
vn    = uadv1*rnx+ vadv1*rny
uadv2 = uadv1 - vn*rnx
vadv2 = vadv1 - vn*rny
rhou2 = rho2*uadv2
rhov2 = rho2*vadv2
!
pres2 = (gamma-1.d0)*(rhoe2-0.5d0*rho2*(uadv2**2 + vadv2**2))
!
unpri(1, 1, ier) = rho2
unpri(1, 2, ier) = uadv2
unpri(1, 3, ier) = vadv2
unpri(1, 4, ier) = pres2
!
!...Riemann solver for farfield...
!
elseif(bface(3,ifa).eq.4)then
!   print*,'far',ifa,intfac(3:4,ifa)
!...Farfield
!
rho1  = uchar(1)
rhou1 = uchar(2)
rhov1 = uchar(3)
rhoe1 = uchar(4)
!
rhom1 = 1.d0/rho1
uadv1 = rhou1*rhom1
vadv1 = rhov1*rhom1
!
pres1 = (gamma-1.d0)*(rhoe1-0.5d0*rho1*(uadv1**2 + vadv1**2))
!
unpri(1, 1, ier) = rho1
unpri(1, 2, ier) = uadv1
unpri(1, 3, ier) = vadv1
unpri(1, 4, ier) = pres1
!
endif
!
else
dx= geoel(1,ier) - geoel(1,iel)
dy= geoel(2,ier) - geoel(2,iel)
endif
!
!...  weighting for this edge
!
      weigh = 1.0
!weigh = 1.0/sqrt(dx*dx + dy*dy )!+ dz*dz)
!
dxwei = weigh*dx
dywei = weigh*dy
!
dmets(1,iel) = dmets(1,iel) + dxwei*dxwei
dmets(2,iel) = dmets(2,iel) + dywei*dywei
dmets(3,iel) = dmets(3,iel) + dxwei*dywei

!
dmets(1,ier) = dmets(1,ier) + dxwei*dxwei
dmets(2,ier) = dmets(2,ier) + dywei*dywei
dmets(3,ier) = dmets(3,ier) + dxwei*dywei
!
du1we = weigh*(unpri(1, ncs+1,ier) - unpri(1, ncs+1,iel))
du2we = weigh*(unpri(1, ncs+2,ier) - unpri(1, ncs+2,iel))
du3we = weigh*(unpri(1, ncs+3,ier) - unpri(1, ncs+3,iel))
du4we = weigh*(unpri(1, ncs+4,ier) - unpri(1, ncs+4,iel))
!
gradu(1,ncs+1,iel) = gradu(1,ncs+1,iel) + dxwei*du1we
gradu(1,ncs+2,iel) = gradu(1,ncs+2,iel) + dxwei*du2we
gradu(1,ncs+3,iel) = gradu(1,ncs+3,iel) + dxwei*du3we
gradu(1,ncs+4,iel) = gradu(1,ncs+4,iel) + dxwei*du4we


gradu(2,ncs+1,iel) = gradu(2,ncs+1,iel) + dywei*du1we
gradu(2,ncs+2,iel) = gradu(2,ncs+2,iel) + dywei*du2we
gradu(2,ncs+3,iel) = gradu(2,ncs+3,iel) + dywei*du3we
gradu(2,ncs+4,iel) = gradu(2,ncs+4,iel) + dywei*du4we

!  excluding the ghost element gradient
if(ifa.gt.nbfac) then
gradu(1,ncs+1,ier) = gradu(1,ncs+1,ier) + dxwei*du1we
gradu(1,ncs+2,ier) = gradu(1,ncs+2,ier) + dxwei*du2we
gradu(1,ncs+3,ier) = gradu(1,ncs+3,ier) + dxwei*du3we
gradu(1,ncs+4,ier) = gradu(1,ncs+4,ier) + dxwei*du4we

gradu(2,ncs+1,ier) = gradu(2,ncs+1,ier) + dywei*du1we
gradu(2,ncs+2,ier) = gradu(2,ncs+2,ier) + dywei*du2we
gradu(2,ncs+3,ier) = gradu(2,ncs+3,ier) + dywei*du3we
gradu(2,ncs+4,ier) = gradu(2,ncs+4,ier) + dywei*du4we
endif
!
1400 continue
!
do 2000 ie  = 1, nelem
!
rxx   = dmets(1,ie)
ryy   = dmets(2,ie)
rxy   = dmets(3,ie)
!
rco1  = ryy
rco2  = -rxy
!nafac
rco3  = -rxy
rco4  = rxx

rjac1 = 1.0/(rxx*ryy - rxy*rxy)
!
!
ru1x  = gradu(1,ncs+1,ie)
ru1y  = gradu(2,ncs+1,ie)

ru2x  = gradu(1,ncs+2,ie)
ru2y  = gradu(2,ncs+2,ie)
!
ru3x  = gradu(1,ncs+3,ie)
ru3y  = gradu(2,ncs+3,ie)
!
ru4x  = gradu(1,ncs+4,ie)
ru4y  = gradu(2,ncs+4,ie)
!
gradu(1,ncs+1,ie) = (ru1x*rco1 + ru1y*rco2 )*rjac1
gradu(2,ncs+1,ie) = (ru1x*rco3 + ru1y*rco4 )*rjac1

gradu(1,ncs+2,ie) = (ru2x*rco1 + ru2y*rco2 )*rjac1
gradu(2,ncs+2,ie) = (ru2x*rco3 + ru2y*rco4 )*rjac1

gradu(1,ncs+3,ie) = (ru3x*rco1 + ru3y*rco2 )*rjac1
gradu(2,ncs+3,ie) = (ru3x*rco3 + ru3y*rco4 )*rjac1
!
gradu(1,ncs+4,ie) = (ru4x*rco1 + ru4y*rco2 )*rjac1
gradu(2,ncs+4,ie) = (ru4x*rco3 + ru4y*rco4 )*rjac1
!     print*,'output'
2000 continue
!
do ie = 1, nelem
!
dx = geoel(4,ie)*0.5d0
dy = geoel(5,ie)*0.5d0
!
do iq =1, nq
unpri(2, iq, ie) = gradu(1, iq, ie)*dx
unpri(3, iq, ie) = gradu(2, iq, ie)*dy
enddo

enddo
!   print*,'gradient'
!    print*,gradu(1:2,ncs+1,9),gradu(1:2,ncs+1,10)
return

end subroutine getgraduls_euler_primitive
!
!...Barth limiter
!
subroutine barthlimit_fvm(inpoel, intfac, ltelem, fsuel, geoel, unkno, coord)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(inout)::unkno
real*8,dimension(1:ngeel,1:nelem+nbfac),    intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),          intent(in)::coord
integer*4,dimension(1:3,1:nelem),            intent(in)::ltelem, fsuel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
!
!...Local
!
integer, dimension(3)::estri
integer:: ip(nvtri)
integer:: ie, iv, iest, iq, ideg,ind, ifa
real*8:: aflim(1:nq, 1:nelem), alfal(1:nelem)
real*8:: unctr(1:nq)
real*8,  dimension(1:nq, 1:nvtri)::alfa
real*8:: xv(3), yv(3)
real*8:: b(3, 1:nvtri)
real*8:: unmax(1:nq), unmin(1:nq), dunk(1:nq)
real*8,  dimension(1:nq,1:3)     ::unest
real*8,dimension(1:nq,  1:nvtri) ::unknv
real*8::uncom(4)
!
real*8:: rho, uvtx, vvtx, evtx, pvtx
real*8:: uctr, vctr, ectr, pctr, eps
real*8:: xc, yc, dx, dy, fiy
!
eps = 1.e-6
!
!  print*,'nvtri',nvtri
!
do ie = 1, nelem
!
estri(1:3) = ltelem(1:3, ie)
ip(1:nvtri) = inpoel(1:nvtri, ie)
!
xv(1:nvtri) = coord(1, ip(1:nvtri))
yv(1:nvtri) = coord(2, ip(1:nvtri))
!

xv(1) = 0.5d0*(coord(1, ip(1)) + coord(1, ip(2)))
xv(2) = 0.5d0*(coord(1, ip(2)) + coord(1, ip(3)))
xv(3) = 0.5d0*(coord(1, ip(1)) + coord(1, ip(3)))
!
yv(1) = 0.5d0*(coord(2, ip(1)) + coord(2, ip(2)))
yv(2) = 0.5d0*(coord(2, ip(2)) + coord(2, ip(3)))
yv(3) = 0.5d0*(coord(2, ip(1)) + coord(2, ip(3)))

!
dx = geoel(4,ie)*0.5d0
dy = geoel(5,ie)*0.5d0
!
xc = geoel(1, ie)
yc = geoel(2, ie)
!
b(1, 1:nvtri) = 1.d0
b(2, 1:nvtri) = (xv(1:nvtri)-xc)/dx
b(3, 1:nvtri) = (yv(1:nvtri)-yc)/dy
!
unctr(1:nq) = unkno(1, 1:nq, ie)
!
! if(ie==1) print*,'update',unctr, unknv(1:4, 1)
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

unmax(iq) = max(uncom(1), unctr(iq))
!
enddo
!
!  if(ie==1) print*,'estri4', unmax
!
!
!...Find the maximum values...
!
do iq =1, nq
!
uncom(1) = minval(unest(iq, 1:3))
!
unmin(iq) = min(uncom(1),unctr(iq))
!
enddo
!  if(ie==1) print*,'estri4min', unmin
!
!...Impose limiter
!
do iv = 1, nvtri
!
do iq = 1,nq
!
dunk(iq) = unknv(iq, iv) - unctr(iq)
!
! if(ie==1) print*,'dunk',iq,iv,dunk(iq),(unmax(iq) - unctr(iq)),(unmax(iq) - unctr(iq))/dunk(iq),(unmin(iq) - unctr(iq)),&
!                        (unmin(iq) - unctr(iq))/dunk(iq)
!
if(dunk(iq).gt.1.d-12)then

fiy = (unmax(iq) - unctr(iq))/dunk(iq)
alfa(iq, iv) = max(min(1.d0, 1.0d0*fiy), 0.d0)
!
!  alfa(iq, iv) = max(min(1.d0, (fiy**2+2.d0*fiy)/(fiy**2 + fiy +2.d0)), 0.d0)
elseif(dunk(iq).lt.-1.d-12)then

fiy = (unmin(iq) - unctr(iq))/dunk(iq)
alfa(iq, iv) = max(min(1.d0, 1.0d0*fiy), 0.d0)
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
aflim(iq, ie) = minval(alfa(iq, 1:nvtri))
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
!print*,'before implementation ',unkno(2:3, 1:4, 1)
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
!

!print*,'after implementation ',unkno(2:3, 1:4, 1)
end subroutine barthlimit_fvm
!
!...Main subroutine for calling barth
!
subroutine barthlimit_li(inpoel, intfac, ltelem, fsuel, geoel, unkno, coord)
use constant
implicit none
!...Input arrays
real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(inout)::unkno
real*8,dimension(1:ngeel,1:nelem+nbfac),    intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),          intent(in)::coord
integer*4,dimension(1:3,1:nelem),            intent(in)::ltelem, fsuel
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvtri,1:nelem),        intent(in)::inpoel
!
!...Local
!
integer:: iq
!real*8::unknb(1:ndegr, 1, 1:nelem+nbfac)
!
!print*,'before implementation ',unkno(2:3, 1:4, 1)

! unknb = 0.d0
!
do iq =1, nq
!
!unknb(1:ndegr, 1, 1:nelem+nbfac) = unkno(1:ndegr, iq, 1:nelem+nbfac)
!
call lim_bjas2(ltelem,fsuel,intfac,coord,geoel,unkno(1:ndegr, iq, 1:nelem+nbfac))
!
!unkno(2, iq, 1:nelem+nbfac) = unknb(2, 1, 1:nelem+nbfac)
!unkno(3, iq, 1:nelem+nbfac) = unknb(3, 1, 1:nelem+nbfac)

enddo
!
!print*,'after implementation ',unkno(2:3, 1:4, 1)

!
end subroutine barthlimit_li
!
!...Barth limiter for every variable...
!
subroutine lim_bjas(esuel,fsuel,intfac,coord,geoel,unk)
use constant
implicit none
integer,intent(in):: esuel(3,nelem),intfac(nifai,nafac),fsuel(3,nelem)
real*8, intent(in)::coord(ndimn,npoin),geoel(ngeel,nelem+nbfac)
real*8, dimension(ndegr, 1, nelem+nbfac), intent(inout) :: unk
!
integer, parameter::ngausl=2
real*8 :: shp1,shp2,xi
real*8 :: xc, yc
real*8 :: weight(ngausf),wi,posgp(1,ngausf)
real*8 :: phig(ngausf), phif(3), phie
real*8 :: xg,yg,xp1,xp2,yp1,yp2,ui,uj,uxi,uyi, umin,umax,dumin,dumax,dui
real*8 :: dx, dy
integer:: ie,je,ig,ifc,j
!
!--- Quadrature weights for the faces:
!
call rutope(1, ngausf, posgp, weight)
!
do ie = 1,nelem

!--- centroids of current cell:
xc = geoel(1,ie)
yc = geoel(2,ie)
!
dx = geoel(4, ie)*0.5d0
dy = geoel(5, ie)*0.5d0

!--- Polynomial in the current cell:
ui  = unk(1,1,ie)
uxi = unk(2,1,ie)/dx
uyi = unk(3,1,ie)/dy

phif = 0.d0

!--- min/max values in the neighborhood:
umin = ui
umax = ui
do j = 1,3
je = esuel(j,ie)
if (je.le.nelem) then
uj = unk(1,1,je)
umin = min(umin,uj)
umax = max(umax,uj)
end if
end do !j
dumin = umin - ui
dumax = umax - ui

!--- Loop over the faces of the current cell:
do j = 1, 3

je = esuel(j,ie)

!--- only internal neighboring-cells considered:
if (je.le.nelem) then

ifc = fsuel(j,ie)
xp1 = coord(1,intfac(3,ifc))
xp2 = coord(1,intfac(4,ifc))
yp1 = coord(2,intfac(3,ifc))
yp2 = coord(2,intfac(4,ifc))

phig = 0.d0

do ig = 1,ngausf
!
xi = posgp(1, ig)
!
shp1 = 0.5d0*(1.d0-xi)
shp2 = 0.5d0*(1.d0+xi)
!
! Gauss point coordinates
!
xg = xp1*shp1 + xp2*shp2
yg = yp1*shp1 + yp2*shp2
! delta r_ie: the position vector of ig from cc(ie)
! delta u_ie
dui = (xg - xc)*uxi + (yg - yc)*uyi

! slope limiter for this gauss point

! dui = 0
if (dabs(dui).le.1.d-8) then
phig(ig) = 1.d0

! dui < 0
elseif (sign(1.d0,dui).lt.0.d0) then
phig(ig) = min(1.d0,(dumin/dui))

! dui > 0
elseif (sign(1.d0,dui).gt.0.d0) then
phig(ig) = min(1.d0,(dumax/dui))
end if

!print*,dui,'|',sign(1.d0,dui),dumin/dui,dumax/dui

end do !ig

!--- min value on the face j:
phif(j) = minval(phig)

else
!--- for ghost-cell neighbors:
phif(j) = 1.d0

end if

end do !j

!stop

!--- min value over all faces of ie:
phie = minval(phif)

!--- limited derivatives:
unk(2,1,ie) = phie*unk(2,1,ie)
unk(3,1,ie) = phie*unk(3,1,ie)

!if(ie==1) print*,ie,dxb(ie)

end do !ie

!stop

end subroutine lim_bjas
!
!...Variation of barth limiter
!
!
!...Barth limiter for every variable...
!
subroutine lim_bjas2(esuel,fsuel,intfac,coord,geoel,unk)
use constant
implicit none
integer,intent(in):: esuel(3,nelem),intfac(nifai,nafac),fsuel(3,nelem)
real*8, intent(in)::coord(ndimn,npoin),geoel(ngeel,nelem+nbfac)
real*8, dimension(ndegr, 1, nelem+nbfac), intent(inout) :: unk
!
integer, parameter::ngausl=1
real*8 :: shp1,shp2,xi
real*8 :: xc, yc
real*8 :: weight(ngausl),wi,posgp(1,ngausl)
real*8 :: phig(ngausl), phif(3), phie
real*8 :: xg,yg,xp1,xp2,yp1,yp2,ui,uj,uxi,uyi, umin,umax,dumin,dumax,dui
real*8 :: dx, dy
integer:: ie,je,ig,ifc,j
!
!--- Quadrature weights for the faces:
!
call rutope(1, ngausl, posgp, weight)
!
!posgp(1,1) = -1.d0
!posgp(1,2) = 1.d0
!posgp(1,3) = 0.d0
!
do ie = 1,nelem

!--- centroids of current cell:
xc = geoel(1,ie)
yc = geoel(2,ie)
!
dx = geoel(4, ie)*0.5d0
dy = geoel(5, ie)*0.5d0

!--- Polynomial in the current cell:
ui  = unk(1,1,ie)
uxi = unk(2,1,ie)/dx
uyi = unk(3,1,ie)/dy

phif = 0.d0

!--- min/max values in the neighborhood:
umin = ui
umax = ui
do j = 1,3
je = esuel(j,ie)
if (je.le.nelem) then
uj = unk(1,1,je)
umin = min(umin,uj)
umax = max(umax,uj)
end if
end do !j
dumin = umin - ui
dumax = umax - ui

!--- Loop over the faces of the current cell:
do j = 1, 3

je = esuel(j,ie)

!--- only internal neighboring-cells considered:
if (je.le.nelem) then

ifc = fsuel(j,ie)
xp1 = coord(1,intfac(3,ifc))
xp2 = coord(1,intfac(4,ifc))
yp1 = coord(2,intfac(3,ifc))
yp2 = coord(2,intfac(4,ifc))

phig = 0.d0

do ig = 1,ngausl
!
xi = posgp(1, ig)
!
shp1 = 0.5d0*(1.d0-xi)
shp2 = 0.5d0*(1.d0+xi)
!
! Gauss point coordinates
!
xg = xp1*shp1 + xp2*shp2
yg = yp1*shp1 + yp2*shp2
! delta r_ie: the position vector of ig from cc(ie)
! delta u_ie
dui = (xg - xc)*uxi + (yg - yc)*uyi

! slope limiter for this gauss point

! dui = 0
if (dabs(dui).le.1.d-8) then
phig(ig) = 1.d0

! dui < 0
elseif (sign(1.d0,dui).lt.0.d0) then
phig(ig) = min(1.d0,(dumin/dui))

! dui > 0
elseif (sign(1.d0,dui).gt.0.d0) then
phig(ig) = min(1.d0,(dumax/dui))
end if

!print*,dui,'|',sign(1.d0,dui),dumin/dui,dumax/dui

end do !ig

!--- min value on the face j:
phif(j) = minval(phig)

else
!--- for ghost-cell neighbors:
phif(j) = 1.d0

end if

end do !j

!stop

!--- min value over all faces of ie:
phie = minval(phif)

!--- limited derivatives:
unk(2,1,ie) = phie*unk(2,1,ie)
unk(3,1,ie) = phie*unk(3,1,ie)

!if(ie==1) print*,ie,dxb(ie)

end do !ie

!stop

end subroutine lim_bjas2
