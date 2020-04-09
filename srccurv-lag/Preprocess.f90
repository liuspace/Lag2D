!
!...subroutine: 1st preprocessing to get some information, i.e., nbfac, nelem...
!
subroutine preprocess1
  use constant
  implicit none
  integer*4::nn,i,j,nc
!  integer*4,dimension(1:nepfa)::lhelp
  integer*4::ielem,inode,ipoin,istor,ifael,nnofj,iface,iftri,jftri
  integer*4::jelem,icoun,jnofa,jpoin,ie,je,numfac,ip1,ip2,ip3
  integer:: ifqua, jfqua
! integer*4,dimension(1:2,1:3)::lpofa
  real*8::x1,x2,y1,y2,x3,y3,dx1,dx2,dx3,dy1,dy2,dy3
  real*8:: tinpt
  integer*4,dimension(1:2)::ip
!
integer*4,allocatable::inpoel(:,:)
  integer*4,allocatable::bface(:,:)
  integer*4,allocatable::esup1(:),esup2(:)
  integer*4,allocatable::lpoin(:), lhelp(:)
  integer*4,allocatable::esuel(:,:), estri(:, :), esqua(:, :)
  integer,  allocatable::ipqua(:,:), iptri(:,:)
  integer,  allocatable::lpofa(:,:), lpqua(:,:)
  real*8,allocatable::coord(:,:)
!
!...Step 1: Reading the input file-dgflo2d.ctrl...
!
  open(5,file='dgflo2d.ctrl')
  read(5,*)
  read(5,*)amach,cfl,nflux,breta,dtfix
  read(5,*)
  read(5,*)nmatel
  read(5,*)
  read(5,*)gamlg,cdrho,tend,cimpd
  read(5,*)
  read(5,*)ntime,ndegr,mdegr,npoly, nreco
  read(5,*)
  read(5,*)nlimi,ngaus,ngausf,ngausd,nint
  read(5,*)
  read(5,*)ngausf_geo, ngausd_geo
  read(5,*)
  read(5,*)ncurv
  read(5,*)
  read(5,*)nresi
  read(5,*)
  read(5,*)nmeth, nq, ndens
  read(5,*)
  read(5,*)ndt, nstag
  read(5,*)
  read(5,*)ncase, nriem, nfint
  read(5,*)
  read(5,*)nrz
  read(5,*)
  read(5,*)nout, ndump
  close(5)
!
!...Definition of some parameters...
!...xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!...nepfa: No. of end points of one face. nepfa=2 for 2D face for both linear and curved.
!...nvfac: No. of vertices of one face.  nvfac=2 for linear face, and nvfac=3 for curved face.
!...nvtri: No. of vertices of one triangle.  nvfac=2 for linear tri, and nvfac=3 for curved tri.
!...nptfa: No. of points of one face, linear face(nvfac=2) could be psedo curevd face(nptfa==3)
!
!...nptri: No. of points of one triangle, linear tri(nvtri=3) could be psedo curved tri(nptri==6)
!...nptfa and nptri are parameters to valiate the curved subroutine...

!...nftri: No. of faces of one triangle...nftri==3 for any knid of triangle...

!...ntria: No. of triangle...
!...nquad: No. of quadrilateral...

!...xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
   if(ncurv == 0)then
    nptfa = 3
    nptri = 6;  npqua = 9
    nepfa = 2;  nftri = 3; nfqua = 4
!
    nvfac = 2; nvtri = 3;  nvqua = 4 !...nvfac : No. of vertex for one face; nvtri: No. of vertex for one triangle
!
    nbfai = 5
!
    allocate(lpofa(nepfa, nftri), lpqua(nepfa, nfqua))
!...triangle
    lpofa(1,1) = 2
    lpofa(2,1) = 3
    lpofa(1,2) = 3
    lpofa(2,2) = 1
    lpofa(1,3) = 1
    lpofa(2,3) = 2
!...quad
    lpqua(1,1) = 1
    lpqua(2,1) = 2
    lpqua(1,2) = 2
    lpqua(2,2) = 3
    lpqua(1,3) = 3
    lpqua(2,3) = 4
    lpqua(1,4) = 4
    lpqua(2,4) = 1
!
   elseif(ncurv == 1)then
    nptfa = 3
    nptri = 6; npqua = 9
    nepfa = 2; nftri = 3;  nfqua = 4
    nvfac = 3; nvtri = 6;  nvqua = 9
!
    allocate(lpofa(nepfa, nftri), lpqua(nepfa, nfqua))
!...triangle
    lpofa(1,1) = 2
    lpofa(2,1) = 3
    lpofa(1,2) = 3
    lpofa(2,2) = 1
    lpofa(1,3) = 1
    lpofa(2,3) = 2
!...quad
    lpqua(1,1) = 1
    lpqua(2,1) = 2
    lpqua(1,2) = 2
    lpqua(2,2) = 3
    lpqua(1,3) = 3
    lpqua(2,3) = 4
    lpqua(1,4) = 4
    lpqua(2,4) = 1
!
    nbfai = 5
!
   elseif(ncurv == 2)then
    nptfa = 4
    nptri = 10; npqua = 12
    nepfa = 2; nftri = 3; nfqua = 4
    nvfac = 4; nvtri = 10; nvqua = 12 !...no nodes in cell...
!
    allocate(lpofa(nepfa, nftri), lpqua(nepfa, nfqua))
!...triangle
    lpofa(1,1) = 2
    lpofa(2,1) = 3
    lpofa(1,2) = 3
    lpofa(2,2) = 1
    lpofa(1,3) = 1
    lpofa(2,3) = 2
!...quad
    lpqua(1,1) = 1
    lpqua(2,1) = 2
    lpqua(1,2) = 2
    lpqua(2,2) = 3
    lpqua(1,3) = 3
    lpqua(2,3) = 4
    lpqua(1,4) = 4
    lpqua(2,4) = 1
!
    nbfai = 6
   endif
!
!...Skip the first nn rows
!
  open(10,file='dgflo2d.domn')
  read(10,*)nn
  do i=1,nn
    read(10,*)!Skip the first nn rows
  enddo
!   
  read(10,*)
  read(10,*)ndimn,ntype
  read(10,*)
  read(10,*)ntria,npoin,nbfac,tinpt,nquad
!
!...Consistent with previous...
!
  nelem = ntria !
  ncell = ntria + nquad !...No. of triangles and quads...
  nsize = ncell + nbfac !...No. of all cells including ghost cell..

!
  allocate(lhelp(1:nepfa))
  allocate(inpoel(1:nvtri,1:nelem), iptri(1:nvtri,1:ntria), ipqua(1:nvqua, 1:nquad))
  allocate(bface(1:nbfai,1:nbfac))
!
  allocate(esuel(1:nftri,1:nelem), estri(1:nftri,1:ntria), esqua(1:nfqua, 1:nquad))
  allocate(lpoin(1:npoin))
  allocate(esup2(1:npoin+1))
!
  allocate(coord(1:ndimn,1:npoin))
!  allocate(geoel(1:ngeel,1:nelem+nbfac))
!
!...read element-connectivity...
!
  read(10,*)
!...triangle
  do i=1,ntria
   read(10,*)j,iptri(1:nvtri,i)
  enddo
!...Consistent with previous...
  inpoel = iptri
!...quad
 do i=1,nquad
   read(10,*)j,ipqua(1:nvqua,i)
 enddo
!
!...read coordinates...
!
  read(10,*)
  do i=1,npoin
   read(10,*)j,coord(1:ndimn,i)
  enddo
!
!...useless input...
!
  read(10,*)
  do i=1,npoin
   read(10,*)
  enddo
!
!...read boundary face connectivity...
!
  read(10,*)
  do i=1,nbfac
!   read(10,*)j,bface(1:nbfai,i)
   read(10,*)j,bface(1:nbfai,i)
!   print*,bface(1:5,i)
  enddo
  close(10)
!
!...end of reading input file...(step 1)
!

!
!...Step 2: find the boundary cell host element
!...it is easy to use the 2 end-points to find the host cell...
  nn=0
  do i=1,nbfac
!...triangle
  do ie=1,ntria
     nc=0
    do j=1, 2 !...refers to the two end points of one line...
     if((bface(j,i).eq.iptri(1,ie)).or.(bface(j,i).eq.iptri(2,ie)).or.(bface(j,i).eq.iptri(3,ie))) then
      nc=nc+1
     endif
    enddo
    if(nc==2)then
     bface(4,i)=ie
     nn=nn+1  !...one indicator to check whether all the host elments of boundary faces are found... 
    endif
  enddo
!...quad...
  do ie=1,nquad
    nc=0
   do j=1, 2 !...refers to the two end points of one line...
    if((bface(j,i).eq.ipqua(1,ie)).or.(bface(j,i).eq.ipqua(2,ie)).or.&
       (bface(j,i).eq.ipqua(3,ie)).or.(bface(j,i).eq.ipqua(4,ie))) then
      nc=nc+1
    endif
   enddo
   if(nc==2)then
    bface(4,i)= ie + ntria
    nn=nn+1  !...one indicator to check whether all the host elments of boundary faces are found...
   endif
  enddo
!
  enddo
!
!  print*, 'No. of boundary faces finding the host cell =', nn,&
!          'No. of boundary faces =', nbfac
!
!...end of step 2...
!
!
!...Step 3: Finding elements surrounding elments...
!
!...3.1: find esup1, esup2 (elements surrounding points 2)...
   esup2(1:npoin+1)=0
!...triangle
   do ielem = 1, ntria
   do inode = 1, nvtri
      esup2(iptri(inode,ielem)+1)=esup2(iptri(inode,ielem)+1)+1
   enddo
   enddo
!...quad
  do ielem = 1, nquad
  do inode = 1, nvqua
     esup2(ipqua(inode,ielem)+1)=esup2(ipqua(inode,ielem)+1)+1
  enddo
  enddo
!
   do ipoin = 2, npoin+1
     esup2(ipoin)=esup2(ipoin)+esup2(ipoin-1)
   enddo
!
!...finding esup1...
   allocate (esup1(esup2(npoin+1)))
!
npoin2 = npoin + 1
npoin1 = esup2(npoin+1)

!
!...tria
   do ielem = 1, ntria
   do inode = 1, nvtri
    ipoin=iptri(inode,ielem)
    istor=esup2(ipoin)+1 !...for point j, storage from esup2(j)+1 to esup2(j+1)
    esup2(ipoin)=istor
    esup1(istor)=ielem
   enddo
   enddo
!
!...quad..
  do ielem = 1, nquad
  do inode = 1, nvqua
    ipoin=ipqua(inode,ielem)
    istor=esup2(ipoin)+1 !...for point j, storage from esup2(j)+1 to esup2(j+1)
    esup2(ipoin)=istor
    esup1(istor)=ielem + ntria !...Global numbering...
  enddo
  enddo

!...restore esup2...

   do ipoin=npoin+1,2,-1
    esup2(ipoin)=esup2(ipoin-1)  
   enddo
    esup2(1)=0
!
!...step 3.2:finding esuel...    
    lpoin=0
    esuel=0
    estri = 0
    esqua = 0
!
!...nftri:No. of faces of one element...
!
do 10 ie    = 1, ntria
do 20 iftri = 1, nftri
!
ielem = ie
!
lhelp(1:nepfa)=iptri(lpofa(1:nepfa,iftri),ie)
lpoin(lhelp(1:nepfa))=1
ipoin=lhelp(1)
do 30 istor=esup2(ipoin)+1,esup2(ipoin+1)
jelem=esup1(istor)
if (jelem.ne.ielem)then !...35 if...
!
if(jelem.le.ntria)then!...triangle
!
je = jelem !...local No. for triangle...
!
do 40 jftri=1,nftri
nnofj=2  !...fixed value for 2d, because the points No. of one face if always equal to 2...
if(nnofj.eq.nepfa)then !...36 if...
icoun=0
do 50 jnofa=1,nepfa
jpoin=iptri(lpofa(jnofa,jftri),jelem)
if(lpoin(jpoin).eq.1) icoun=icoun+1
50 enddo
if(icoun.eq.nepfa)then
estri(iftri,ie   )=jelem
estri(jftri,je   )=ielem
if(ielem.eq.33)print*,'badt',estri(iftri,ielem),estri(jftri,jelem),jelem
endif
endif   !...36 if...
40 enddo
!
elseif(jelem.gt.ntria)then!...quad
!
je = jelem - ntria !...local No. for quads
!
do 90 jfqua=1,nfqua
nnofj=2  !...fixed value for 2d, because the points No. of one face if always equal to 2...
if(nnofj.eq.nepfa)then
icoun=0
do 100 jnofa=1,nepfa
jpoin=ipqua(lpqua(jnofa,jfqua),je)
if(lpoin(jpoin).eq.1) icoun=icoun+1
100 enddo
if(icoun.eq.nepfa)then
estri(iftri,ie  )=jelem
esqua(jfqua,je  )=ielem
if(ielem.eq.33)print*,'badt2',estri(iftri,ielem),esqua(jfqua,jelem-ntria),jelem,iftri
endif
endif
90 enddo
endif
!
endif !...35 if...
30 enddo
lpoin(lhelp(1:nepfa))=0
20 enddo
10 enddo
!
!...quad
!
do 110 ie    = 1, nquad
do 120 ifqua = 1, nfqua
!
ielem = ie + ntria
!
lhelp(1:nepfa)=ipqua(lpqua(1:nepfa,ifqua),ie)
lpoin(lhelp(1:nepfa))=1
ipoin=lhelp(1)
do 130 istor=esup2(ipoin)+1,esup2(ipoin+1)
jelem=esup1(istor)
!
if (jelem.ne.ielem)then !...35 if...
!
if(jelem.le.ntria)then!...triangle
!
je = jelem !...local No. for triangle...
!
do 140 jftri=1,nftri
nnofj=2  !...fixed value for 2d, because the points No. of one face if always equal to 2...
if(nnofj.eq.nepfa)then !...36 if...
icoun=0
do 150 jnofa=1,nepfa
jpoin=iptri(lpofa(jnofa,jftri),je)
if(lpoin(jpoin).eq.1) icoun=icoun+1
150 enddo
if(icoun.eq.nepfa)then
esqua(ifqua,ie )=jelem
estri(jftri,je )=ielem
if(ielem.eq.1)print*,'bad',esqua(ifqua,ielem),estri(jftri,jelem),jelem,jftri
endif
endif   !...36 if...
140 enddo
!
elseif(jelem.gt.ntria)then!...quad
!
je = jelem - ntria !...local No. for quads
!
do 190 jfqua=1,nfqua
nnofj=2  !...fixed value for 2d, because the points No. of one face if always equal to 2...
if(nnofj.eq.nepfa)then
icoun=0
do 1100 jnofa=1,nepfa
jpoin=ipqua(lpqua(jnofa,jfqua),je)
if(lpoin(jpoin).eq.1) icoun=icoun+1
1100 enddo
if(icoun.eq.nepfa)then
esqua(ifqua,ie )=jelem
esqua(jfqua,je )=ielem
!if(ielem.eq.1)print*,'badq',esqua(ifqua,ielem),esqua(jfqua,jelem-ntria),jelem,ifqua,jfqua
endif
endif
190 enddo
endif
!
endif !...35 if...
130 enddo
lpoin(lhelp(1:nepfa))=0
120 enddo
110 enddo

!print*,'estri',estri(:,64),esqua(:,1)
!
!...Step 4: finding No. of all faces...
!
     numfac=nbfac
!...triangle
      do ielem = 1,ntria
      do iftri = 1,nftri
          jelem=estri(iftri,ielem)
          if ((jelem .gt. ielem) .and. (jelem .le. (ntria+nquad))) then
            numfac=numfac+1
          endif
      enddo
      enddo
!...quad
   do ielem = 1,nquad
   do ifqua = 1,nfqua
     jelem=esqua(ifqua,ielem)
    if ((jelem .gt. ielem+ntria) .and. (jelem .le. (ntria+nquad))) then
     numfac=numfac+1
    endif
   enddo
   enddo
!
!...No. of all the faces: nafac
   nafac = numfac
!
!  print*, 'No. of all faces =', nafac,&
!          'No. of boundary faces =', nbfac
!
!
  deallocate(inpoel, lhelp)
  deallocate(coord)
  deallocate(bface)
  deallocate(esuel)
  deallocate(lpoin)
  deallocate(esup2)
  deallocate(estri, esqua, iptri, ipqua)

  return
end subroutine preprocess1
!
!...subroutine: 2nd preprocessing to get required arrays, i.e., geoel, bface...
!
subroutine preprocess2(intfac,inpoel,iptri, ipqua, geofa,geoel,bface,coord, esuv1, esuv2)
  use constant
  implicit none
 ! integer*4::nelem,npoin,nbfac,ndimn,ntype
  real*8::omega,l1,l2
  integer*4::nn,i,j,nc,ifa,iel,ier
  integer*4,dimension(1:nvtri,1:nelem)::inpoel
  integer*4,dimension(1:nvtri,1:ntria)::iptri
  integer*4,dimension(1:nvqua,1:nquad)::ipqua
  integer*4,dimension(1:nifai,1:nafac)::intfac
  integer*4,dimension(1:nbfai,1:nbfac)::bface
  integer*4,dimension(1:nbfac)::bflag
  integer*4,dimension(1:nepfa)::lhelp
integer*4, intent(out)::esuv1(npoin1),esuv2(npoin2)
!
  integer*4,allocatable::esup1(:),esup2(:)
  integer*4::ielem,inode,ipoin,istor,iftri,jftri,nnofj,iface
  integer*4::jelem,icoun,jnofa,jpoin,ie,je,numfac,ip1,ip2,ip3
  integer:: ifqua, jfqua
  integer*4,allocatable::lpofac(:, :),lpofa(:, :), lpqua(:, :)
  integer*4,allocatable::lpoin(:)
  integer*4,allocatable::esuel(:,:), estri(:, :), esqua(:, :)
!  integer*4,dimension(1:2,1:3)::lpofa
  real*8::x1,x2,y1,y2,x3,y3,dx1,dx2,dx3,dy1,dy2,dy3
  real*8::x0,y0,la,lb,lc,abp2,abm2,xg,yg
  integer*4,dimension(1:2)::ip
!
  real*8,dimension(1:ndimn,1:npoin)::coord
  real*8,dimension(1:ngefa,1:nafac)::geofa
  real*8,dimension(1:ngeel,1:nsize)::geoel

   allocate(esuel(1:nftri,1:nelem), estri(1:nftri,1:ntria), esqua(1:nfqua, 1:nquad))
   allocate(lpoin(1:npoin))
   allocate(esup2(1:npoin+1))

   if(ncurv==0)then
    allocate(lpofa(nvfac, nftri), lpqua(nvfac, nfqua))  !...need changing...
!
    lpofa(1,1) = 2
    lpofa(2,1) = 3
    lpofa(1,2) = 3
    lpofa(2,2) = 1
    lpofa(1,3) = 1
    lpofa(2,3) = 2
!...quad
    lpqua(1,1) = 1
    lpqua(2,1) = 2
    lpqua(1,2) = 2
    lpqua(2,2) = 3
    lpqua(1,3) = 3
    lpqua(2,3) = 4
    lpqua(1,4) = 4
    lpqua(2,4) = 1
!
   elseif(ncurv==1)then
    allocate(lpofa(nvfac, nftri), lpqua(nvfac, nfqua) )  !...need changing...
!...triangle
    lpofa(1,1) = 2
    lpofa(2,1) = 3
    lpofa(3,1) = 5

    lpofa(1,2) = 3
    lpofa(2,2) = 1
    lpofa(3,2) = 6

    lpofa(1,3) = 1
    lpofa(2,3) = 2
    lpofa(3,3) = 4
!...quad
    lpqua(1,1) = 1
    lpqua(2,1) = 2
    lpqua(3,1) = 5

    lpqua(1,2) = 2
    lpqua(2,2) = 3
    lpqua(3,2) = 6

    lpqua(1,3) = 3
    lpqua(2,3) = 4
    lpqua(3,3) = 7

    lpqua(1,4) = 4
    lpqua(2,4) = 1
    lpqua(3,4) = 8
!
   elseif(ncurv==2)then
    allocate(lpofa(nvfac, nftri), lpqua(nvfac, nfqua) )  !...need changing...
    lpofa(1,1) = 2
    lpofa(2,1) = 3
    lpofa(3,1) = 5
    lpofa(4,1) = 8

    lpofa(1,2) = 3
    lpofa(2,2) = 1
    lpofa(3,2) = 6
    lpofa(4,2) = 9

    lpofa(1,3) = 1
    lpofa(2,3) = 2
    lpofa(3,3) = 4
    lpofa(4,3) = 7
!
!...quad
    lpqua(1,1) = 1
    lpqua(2,1) = 2
    lpqua(3,1) = 5
    lpqua(4,1) = 9

    lpqua(1,2) = 2
    lpqua(2,2) = 3
    lpqua(3,2) = 6
    lpqua(4,2) = 10

    lpqua(1,3) = 3
    lpqua(2,3) = 4
    lpqua(3,3) = 7
    lpqua(4,3) = 11

    lpqua(1,4) = 4
    lpqua(2,4) = 1
    lpqua(3,4) = 8
    lpqua(4,4) = 12

   endif
!
!...Skip the first nn rows
!
    open(10,file='dgflo2d.domn')
    read(10,*)nn
    do i=1,nn
     read(10,*)!Skip the first nn rows
    enddo
!
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
!
!...read element-connectivity...
!
!...triangle
  do i=1,ntria
     read(10,*)j,iptri(1:nvtri,i)
  enddo
!...quad
  do i=1,nquad
    read(10,*)j,ipqua(1:nvqua,i)
  enddo
!...Consistent with previous...
inpoel = iptri
!
!...read coordinates...
!
  read(10,*)
  do i=1,npoin 
    read(10,*)j,coord(1:ndimn,i)
  enddo
!
!coord(2,:) = coord(2,:) + 0.1d0
!
!...useless input...
!
  read(10,*)
  do i=1,npoin
    read(10,*)
  enddo
!
!...read boundary face connectivity...
!
  read(10,*)
  do i=1,nbfac
!    read(10,*)j,bface(1:nbfai,i)
    read(10,*)j,bface(1:nbfai,i)
  enddo
  close(10)
!
if(nmeth.eq.2)then
bflag(1:nbfac) = bface(4,1:nbfac)
endif
!
!...end of reading input file...(step 1)
!
!
!...Step 2: find the boundary cell host element
!...it is easy to use the 2 end-points to find the host cell...
  nn=0
 do i=1,nbfac
!...triangle
 do ie=1,ntria
    nc=0
    do j=1, 2 !...refers to the two end points of one line...
      if((bface(j,i).eq.iptri(1,ie)).or.(bface(j,i).eq.iptri(2,ie)).or.(bface(j,i).eq.iptri(3,ie))) then
        nc=nc+1
      endif
    enddo
  if(nc==2)then
    bface(4,i)=ie
    nn=nn+1  !...one indicator to check whether all the host elments of boundary faces are found...
  endif
 enddo
!...quad...
 do ie=1,nquad
    nc=0
    do j=1, 2 !...refers to the two end points of one line...
      if((bface(j,i).eq.ipqua(1,ie)).or.(bface(j,i).eq.ipqua(2,ie)).or.&
        (bface(j,i).eq.ipqua(3,ie)).or.(bface(j,i).eq.ipqua(4,ie))) then
        nc=nc+1
      endif
    enddo
   if(nc==2)then
     bface(4,i)= ie + ntria
     nn=nn+1  !...one indicator to check whether all the host elments of boundary faces are found...
   endif
 enddo
!
enddo
!
!...end of step 2...
!

!
!...Step 3: Finding elements surrounding elments...
!
!...3.1: find esup1, esup2 (elements surrounding points 2)...
   esup2(1:npoin+1)=0
!...triangle
do ielem = 1, ntria
do inode = 1, nvtri
esup2(iptri(inode,ielem)+1)=esup2(iptri(inode,ielem)+1)+1
enddo
enddo
!...quad
do ielem = 1, nquad
do inode = 1, nvqua
esup2(ipqua(inode,ielem)+1)=esup2(ipqua(inode,ielem)+1)+1
enddo
enddo
!
!   
   do ipoin=2,npoin+1
     esup2(ipoin)=esup2(ipoin)+esup2(ipoin-1)
   enddo
!
!...finding esup1...
   allocate (esup1(esup2(npoin+1)))
!
!...tria
do ielem = 1, ntria
do inode = 1, nvtri
ipoin=iptri(inode,ielem)
istor=esup2(ipoin)+1 !...for point j, storage from esup2(j)+1 to esup2(j+1)
esup2(ipoin)=istor
esup1(istor)=ielem
enddo
enddo
!
!...quad..
do ielem = 1, nquad
do inode = 1, nvqua
ipoin=ipqua(inode,ielem)
istor=esup2(ipoin)+1 !...for point j, storage from esup2(j)+1 to esup2(j+1)
esup2(ipoin)=istor
esup1(istor)=ielem + ntria !...Global numbering...
enddo
enddo
!
!...restore esup2...
   do ipoin=npoin+1,2,-1
    esup2(ipoin)=esup2(ipoin-1)  
   enddo
    esup2(1)=0
!
!...step 3.2:finding esuel... 
    lpoin=0
    esuel=0
    estri = 0
    esqua = 0
!         
do 10 ie    = 1, ntria
do 20 iftri = 1, nftri
!
ielem = ie
!
lhelp(1:nepfa)=iptri(lpofa(1:nepfa,iftri),ie)
lpoin(lhelp(1:nepfa))=1
ipoin=lhelp(1)
do 30 istor=esup2(ipoin)+1,esup2(ipoin+1)
jelem=esup1(istor)
if (jelem.ne.ielem)then !...35 if...
!
if(jelem.le.ntria)then!...triangle
!
je = jelem !...local No. for triangle...
!
do 40 jftri=1,nftri
nnofj=2  !...fixed value for 2d, because the points No. of one face if always equal to 2...
if(nnofj.eq.nepfa)then !...36 if...
icoun=0
do 50 jnofa=1,nepfa
jpoin=iptri(lpofa(jnofa,jftri),jelem)
if(lpoin(jpoin).eq.1) icoun=icoun+1
50 enddo
if(icoun.eq.nepfa)then
estri(iftri,ie   )=jelem
estri(jftri,je   )=ielem
if(ielem.eq.33)print*,'badt',estri(iftri,ielem),estri(jftri,jelem),jelem
endif
endif   !...36 if...
40 enddo
!
elseif(jelem.gt.ntria)then!...quad
!
je = jelem - ntria !...local No. for quads
!
do 90 jfqua=1,nfqua
nnofj=2  !...fixed value for 2d, because the points No. of one face if always equal to 2...
if(nnofj.eq.nepfa)then
icoun=0
do 100 jnofa=1,nepfa
jpoin=ipqua(lpqua(jnofa,jfqua),je)
if(lpoin(jpoin).eq.1) icoun=icoun+1
100 enddo
if(icoun.eq.nepfa)then
estri(iftri,ie  )=jelem
esqua(jfqua,je  )=ielem
if(ielem.eq.33)print*,'badt2',estri(iftri,ielem),esqua(jfqua,jelem-ntria),jelem,iftri
endif
endif
90 enddo
endif
!
endif !...35 if...
30 enddo
lpoin(lhelp(1:nepfa))=0
20 enddo
10 enddo
!
!...quad
!
do 110 ie    = 1, nquad
do 120 ifqua = 1, nfqua
!
ielem = ie + ntria
!
lhelp(1:nepfa)=ipqua(lpqua(1:nepfa,ifqua),ie)
lpoin(lhelp(1:nepfa))=1
ipoin=lhelp(1)
do 130 istor=esup2(ipoin)+1,esup2(ipoin+1)
jelem=esup1(istor)
!
if (jelem.ne.ielem)then !...35 if...
!
if(jelem.le.ntria)then!...triangle
!
je = jelem !...local No. for triangle...
!
do 140 jftri=1,nftri
nnofj=2  !...fixed value for 2d, because the points No. of one face if always equal to 2...
if(nnofj.eq.nepfa)then !...36 if...
icoun=0
do 150 jnofa=1,nepfa
jpoin=iptri(lpofa(jnofa,jftri),je)
if(lpoin(jpoin).eq.1) icoun=icoun+1
150 enddo
if(icoun.eq.nepfa)then
esqua(ifqua,ie )=jelem
estri(jftri,je )=ielem
if(ielem.eq.1)print*,'bad',esqua(ifqua,ielem),estri(jftri,jelem),jelem,jftri
endif
endif   !...36 if...
140 enddo
!
elseif(jelem.gt.ntria)then!...quad
!
je = jelem - ntria !...local No. for quads
!
do 190 jfqua=1,nfqua
nnofj=2  !...fixed value for 2d, because the points No. of one face if always equal to 2...
if(nnofj.eq.nepfa)then
icoun=0
do 1100 jnofa=1,nepfa
jpoin=ipqua(lpqua(jnofa,jfqua),je)
if(lpoin(jpoin).eq.1) icoun=icoun+1
1100 enddo
if(icoun.eq.nepfa)then
esqua(ifqua,ie )=jelem
esqua(jfqua,je )=ielem
!if(ielem.eq.1)print*,'badq',esqua(ifqua,ielem),esqua(jfqua,jelem-ntria),jelem,ifqua,jfqua
endif
endif
190 enddo
endif
!
endif !...35 if...
130 enddo
lpoin(lhelp(1:nepfa))=0
120 enddo
110 enddo

!
!...Step 4: finding intfac...
!
!
!...4.1: finding intfac for boundary face...
!     do iface=1,nbfac
!        ip(1)=bface(1,iface)
!        ip(2)=bface(2,iface)
!        lelem(:)=0
!        do ipoin=esup2(ip(1))+1,esup2(ip(1)+1)
!          lelem(esup1(ipoin))=1
!        enddo
!        do ipoin=esup2(ip(2))+1,esup2(ip(2)+1)
!          if (lelem(esup1(ipoin)) .eq. 1) then
!            intfac(1,iface)=esup1(ipoin)
!            intfac(2,iface)=ntria + nquad + iface
!            intfac(3,iface)=ip(1)
!            intfac(4,iface)=ip(2)
!            if(ncurv==1) then
!             intfac(5,iface) =  bface(5,iface)  !...The third high-order node...
!            elseif(ncurv==2) then
!             intfac(5:6,iface) =  bface(5:6,iface)  !...The third high-order node...
!            endif
!          endif
!        enddo
!      enddo
!
!
do iface=1,nbfac
ip(1)=bface(1,iface)
ip(2)=bface(2,iface)
!
intfac(1,iface)=bface(4,iface)
intfac(2,iface)=ncell + iface
intfac(3,iface)=ip(1)
intfac(4,iface)=ip(2)
if(ncurv==1) then
intfac(5,iface) =  bface(5,iface)  !...The third high-order node...
elseif(ncurv==2) then
intfac(5:6,iface) =  bface(5:6,iface)  !...The third high-order node...
endif

enddo
!
!...4.2: finding intfac for internal face...
     numfac=nbfac
!
do ie   =1,ntria
do iftri=1,nftri
!
ielem = ie
!
jelem=estri(iftri,ie)
if ((jelem .gt. ielem) .and. (jelem .le. ncell)) then !...triangle
numfac=numfac+1
intfac(1,numfac)=ielem
intfac(2,numfac)=jelem
intfac(3,numfac)=iptri(lpofa(1,iftri),ie)
intfac(4,numfac)=iptri(lpofa(2,iftri),ie)
if(ncurv==1)then
intfac(5,numfac)   = iptri(lpofa(3,iftri),ie)    !...The third high-order node...
elseif(ncurv==2)then
intfac(5:6,numfac) = iptri(lpofa(3:4,iftri),ie)    !...The third high-order node...
endif
endif
enddo
enddo
!...quad
do ie   =1,nquad
do ifqua=1,nfqua
!
ielem = ie + ntria
!
jelem=esqua(ifqua, ie)
if ((jelem .gt. ielem) .and. (jelem .le. ncell)) then !...triangle
numfac=numfac+1
intfac(1,numfac)=ielem
intfac(2,numfac)=jelem
intfac(3,numfac)=ipqua(lpqua(1,ifqua),ie)
intfac(4,numfac)=ipqua(lpqua(2,ifqua),ie)
if(ncurv==1)then
intfac(5,numfac) = ipqua(lpqua(3,ifqua),ie)    !...The third high-order node...
elseif(ncurv==2)then
intfac(5:6,numfac) = ipqua(lpqua(3:4,ifqua),ie)    !...The third high-order node...
endif
endif
enddo
enddo
!
!
if(nmeth.eq.2)then
bface(4,1:nbfac) = bflag(1:nbfac)  
endif
!
!...Debugging: compare intfac with bface
 !  do i=1,nbfac
 !  print*,i
 !  print*,intfac(1:4,i)
 !  enddo
!
!...4.3: find the the local number of one point in one element.
!
esuv1 = esup1
esuv2 = esup2
!...Removed...
 !  print*,'!successfully complete the preprocessing'

   deallocate(esuel)
   deallocate(lpoin)
   deallocate(esup2)

end subroutine preprocess2
!
!...subroutine: Calculate the mass matrix for DGMs, npoly=1 for DGP1, and npoly=2 for DGP2...
!
subroutine  getamatr(amatr,geoel,coord,inpoel)
 use constant
 implicit none 
 real*8,dimension(1:ngeel,1:nsize)::geoel
 real*8,dimension(1:ndimn,1:npoin),intent(in)::coord 
 real*8,dimension(1:nmatr,1:nelem),intent(out)::amatr
 integer,allocatable::ip
 integer:: inpoel(1:nvtri,1:nelem)
 integer,parameter::ngausm = 13 !...ngausm is used to calculate the geometry information...
 integer :: ie, ig 
 real*8::xp(1:3), yp(1:3)
 real*8:: xc,yc,xg,yg,dx,dy,f1,f2,f3,sp1,sp2,sp3
 real*8::wi,djac, volel, f4
 real*8::f5,f6,f7,f8
 real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
 real*8::b2,b3,b4,b5,b6
 real*8,allocatable::x5(:,:),b55(:),mmatr(:,:)
 real*8,allocatable::weigh(:), posi(:,:) 

   if(npoly==2) allocate(x5(5,5), mmatr(5,5), b55(5))
!
!...1st part: Get the geometry information(ngausm gauss points)...
!
   allocate(weigh(ngausm), posi(2,ngausm))
   call rutope(ndimn, ngausm, posi, weigh)
!
!...get p2 term for geoel...
!
 do ie = 1,nelem !...(1)ie = 1,nelem
!
    xp(1:3) = coord(1,inpoel(1:3,ie))
    yp(1:3) = coord(2,inpoel(1:3,ie))
  
    xc = geoel(1, ie)
    yc = geoel(2, ie)

    dx = geoel(4, ie)*0.5d0
    dy = geoel(5, ie)*0.5d0

    volel = geoel(3, ie)

    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0  
!
!...
    f4 = 0.d0
    f5 = 0.d0
    f6 = 0.d0
    f7 = 0.d0  

  do ig =1,ngausm
     sp1 = posi(1,ig)
     sp2 = posi(2,ig)
     sp3 = 1.d0-sp1-sp2
  
     wi = weigh(ig)*volel

     xg = sp1*xp(1) + sp2*xp(2) + sp3*xp(3)
     yg = sp1*yp(1) + sp2*yp(2) + sp3*yp(3)

     f1 = f1+(xg-xc)/dx*(xg-xc)/dx/2.d0*wi
     f2 = f2+(xg-xc)/dx*(yg-yc)/dy*wi
     f3 = f3+(yg-yc)/dy*(yg-yc)/dy/2.d0*wi
!
!
     f4 = f4 + ((xg-xc)/dx)**3/6.d0*wi
     f5 = f5 + ((yg-yc)/dy)**3/6.d0*wi
     f6 = f6 + ((xg-xc)/dx)**2*((yg-yc)/dy)/2.d0*wi
     f7 = f7 + ((xg-xc)/dx)*((yg-yc)/dy)**2/2.d0*wi

  enddo
     
     geoel(11,ie) = f1/volel
     geoel(12,ie) = f3/volel
     geoel(13,ie) = f2/volel
!
     geoel(14,ie) = f4/volel
     geoel(15,ie) = f5/volel
     geoel(16,ie) = f6/volel
     geoel(17,ie) = f7/volel
! 
 enddo !...(1)ie = 1,nelem
!
!...2nd part: Get the mass matrix(ngausd gauss points which is same as the domain integral)...
!
  deallocate(weigh, posi)
  allocate(weigh(ngausd), posi(2,ngausd))
  call rutope(ndimn, ngausd, posi, weigh)
!
!...get amatr...
!
 do ie = 1,nelem !...(2)ie = 1,nelem
!
    xp(1:3) = coord(1,inpoel(1:3,ie))
    yp(1:3) = coord(2,inpoel(1:3,ie))
  
    xc = geoel(1, ie)
    yc = geoel(2, ie)

    dx = geoel(4, ie)*0.5d0
    dy = geoel(5, ie)*0.5d0

    volel = geoel(3, ie)

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
     sp1 = posi(1,ig)
     sp2 = posi(2,ig)
     sp3 = 1.d0-sp1-sp2
  
     wi = weigh(ig)*volel

     xg = sp1*xp(1) + sp2*xp(2) + sp3*xp(3)
     yg = sp1*yp(1) + sp2*yp(2) + sp3*yp(3)
       
     b2 = (xg-xc)/dx
     b3 = (yg-yc)/dy
     b4 = 0.5d0*b2**2-geoel(11,ie)
     b5 = 0.5d0*b3**2-geoel(12,ie)
     b6 = b2*b3 - geoel(13,ie)

     f1 = f1+(xg-xc)/dx*(xg-xc)/dx*wi
     f2 = f2+(xg-xc)/dx*(yg-yc)/dy*wi
     f3 = f3+(yg-yc)/dy*(yg-yc)/dy*wi

     if(npoly==2)then
     f22 = f22 + b2*b2*wi
     f23 = f23 + b2*b3*wi
     f24 = f24 + b2*b4*wi
     f25 = f25 + b2*b5*wi
     f26 = f26 + b2*b6*wi

     f33 = f33 + b3*b3*wi
     f34 = f34 + b3*b4*wi
     f35 = f35 + b3*b5*wi
     f36 = f36 + b3*b6*wi

     f44 = f44 + b4*b4*wi
     f45 = f45 + b4*b5*wi
     f46 = f46 + b4*b6*wi

     f55 = f55 + b5*b5*wi
     f56 = f56 + b5*b6*wi

     f66 = f66 + b6*b6*wi
     endif
!   if(ie==1) print*,'gauss',xg, yg
  enddo

  if(npoly==1)then
     djac = f1*f3-f2**2
     
     amatr(1,ie) = f3/djac
     amatr(2,ie) = -f2/djac
     amatr(3,ie) = f1/djac 
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
     
  endif
 enddo !...(2)ie = 1,nelem

  deallocate(weigh, posi)
end subroutine  getamatr 
!
!...subroutine: Calculate the mass matrix of curved elements for DGMs, npoly=1 for DGP1, and npoly=2 for DGP2...
!
 subroutine  getamatr_curv(amatr,geoel,coord,inpoel,iptri, ipqua)
 use constant
 implicit none
!
 real*8,dimension(1:nmatr,1:ncell),intent(out)::amatr
 real*8,dimension(1:ngeel, 1:nsize), intent(inout)::geoel
 real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
!
 integer*4,dimension(1:nvtri, 1:nelem),intent(in)::inpoel
 integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
 integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
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
  real*8:: x5(5,5), mmatr(5,5), b55(5)
!...local real number
  real*8:: dxdr,dxds,dydr,dyds 
  real*8:: eps,c00,c10,c05,c20
  real*8:: r, s, djaco, volel
  real*8:: wi, xc, yc, xgaus, ygaus, dxc, dyc
  real*8:: f1, f2, f3, f4
  real*8::f5,f6,f7,f8
  real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
  real*8::b2,b3,b4,b5,b6, djac
!
  integer::ielem, igaus, ishp
  integer:: ie
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
 do ie = 1, ntria !...(1)ie = 1,nelem
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
    xc    = geoel(1, ie)
    yc    = geoel(2, ie)
    volel = geoel(3, ie)
    dxc   = geoel(4, ie)*0.5d0
    dyc   = geoel(5, ie)*0.5d0
!
!...
    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0
!
  do igaus =1, ngausd_geo
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
     xgaus = 0.d0
     ygaus = 0.d0     
!
     do ishp = 1, nptri  
      xgaus = xgaus + shp(ishp)*xp(1,ishp) 
      ygaus = ygaus + shp(ishp)*xp(2,ishp) 
     enddo  
!
     f1 = f1 + (xgaus-xc)/dxc*(xgaus-xc)/dxc/2.d0*djaco
     f2 = f2 + (xgaus-xc)/dxc*(ygaus-yc)/dyc*djaco
     f3 = f3 + (ygaus-yc)/dyc*(ygaus-yc)/dyc/2.d0*djaco
! 
  enddo
!
    geoel(6, ielem) = f1/volel
    geoel(7, ielem) = f3/volel
    geoel(8, ielem) = f2/volel
!
!    print*,'old geoel',ielem, geoel(11:13, ielem)
!    print*,'new geoel',ielem, f1/volel, f3/volel, f2/volel
! 
 enddo !...(1)ie = 1,nelem
!
!...2nd loop for quads...
!
do ie = 1, nquad !...(1)ie = 1,nelem
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
endif

!
!...
xc    = geoel(1, ielem)
yc    = geoel(2, ielem)
volel = geoel(3, ielem)
dxc   = geoel(4, ielem)*0.5d0
dyc   = geoel(5, ielem)*0.5d0
!
!...
f1 = 0.d0
f2 = 0.d0
f3 = 0.d0
!
do igaus =1, ngausd_geoq
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
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo
!
f1 = f1 + (xgaus-xc)/dxc*(xgaus-xc)/dxc/2.d0*djaco
f2 = f2 + (xgaus-xc)/dxc*(ygaus-yc)/dyc*djaco
f3 = f3 + (ygaus-yc)/dyc*(ygaus-yc)/dyc/2.d0*djaco
!
enddo
!
geoel(6, ielem) = f1/volel
geoel(7, ielem) = f3/volel
geoel(8, ielem) = f2/volel
!
!    print*,'old geoel',ielem, geoel(11:13, ielem)
!    print*,'new geoel',ielem, f1/volel, f3/volel, f2/volel
!
enddo !...(1)ie = 1,nelem

!
!...3rd loop to find amatr...
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
!...Cell center...
    xc    = geoel(1, ielem)
    yc    = geoel(2, ielem)
    volel = geoel(3, ielem)
    dxc   = geoel(4, ielem)*0.5d0
    dyc   = geoel(5, ielem)*0.5d0
!
!...
    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0    
    f4 = 0.d0

    f22 = 0.d0
    f23 = 0.d0
    f24 = 0.d0
    f25 = 0.d0
    f26 = 0.d0


    f33 = 0.d0
    f34 = 0.d0
    f35 = 0.d0
    f36 = 0.d0

    f44 = 0.d0
    f45 = 0.d0
    f46 = 0.d0

    f55 = 0.d0
    f56 = 0.d0

    f66 = 0.d0

   do igaus =1, ngausd_geo
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
     xgaus = 0.d0
     ygaus = 0.d0     
!
     do ishp = 1, nptri  
      xgaus = xgaus + shp(ishp)*xp(1,ishp) 
      ygaus = ygaus + shp(ishp)*xp(2,ishp) 
     enddo  
       
     b2 = (xgaus-xc)/dxc
     b3 = (ygaus-yc)/dyc
     b4 = 0.5d0*b2**2-geoel(6,ielem)
     b5 = 0.5d0*b3**2-geoel(7,ielem)
     b6 = b2*b3 - geoel(8,ielem)

     f1 = f1+(xgaus-xc)/dxc*(xgaus-xc)/dxc*djaco
     f2 = f2+(xgaus-xc)/dxc*(ygaus-yc)/dyc*djaco
     f3 = f3+(ygaus-yc)/dyc*(ygaus-yc)/dyc*djaco

     if(npoly==2)then
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
  enddo

  if(npoly==1)then
!
     djac = f1*f3-f2**2    
     amatr(1,ielem) = f3/djac
     amatr(2,ielem) = -f2/djac
     amatr(3,ielem) = f1/djac 
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
     
  endif
!
!    print*,'old amatr',ielem, amatr(1:15, ielem)
!    print*,'new amatr',ielem, x5
!
 enddo !...(2)ie = 1,nelem
!
!...4th loop to find amatr for quads...
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
endif!
!...Cell center...
xc    = geoel(1, ielem)
yc    = geoel(2, ielem)
volel = geoel(3, ielem)
dxc   = geoel(4, ielem)*0.5d0
dyc   = geoel(5, ielem)*0.5d0
!
!...
f1 = 0.d0
f2 = 0.d0
f3 = 0.d0
f4 = 0.d0

f22 = 0.d0
f23 = 0.d0
f24 = 0.d0
f25 = 0.d0
f26 = 0.d0


f33 = 0.d0
f34 = 0.d0
f35 = 0.d0
f36 = 0.d0

f44 = 0.d0
f45 = 0.d0
f46 = 0.d0

f55 = 0.d0
f56 = 0.d0

f66 = 0.d0

do igaus =1, ngausd_geoq
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
xgaus = 0.d0
ygaus = 0.d0
!
do ishp = 1, npqua
xgaus = xgaus + shpq(ishp)*xpq(1,ishp)
ygaus = ygaus + shpq(ishp)*xpq(2,ishp)
enddo

b2 = (xgaus-xc)/dxc
b3 = (ygaus-yc)/dyc
b4 = 0.5d0*b2**2-geoel(6,ielem)
b5 = 0.5d0*b3**2-geoel(7,ielem)
b6 = b2*b3 - geoel(8,ielem)

f1 = f1+(xgaus-xc)/dxc*(xgaus-xc)/dxc*djaco
f2 = f2+(xgaus-xc)/dxc*(ygaus-yc)/dyc*djaco
f3 = f3+(ygaus-yc)/dyc*(ygaus-yc)/dyc*djaco

if(npoly==2)then
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
enddo

if(npoly==1)then
!
djac = f1*f3-f2**2
amatr(1,ielem) = f3/djac
amatr(2,ielem) = -f2/djac
amatr(3,ielem) = f1/djac
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

endif
!
!    print*,'old amatr',ielem, amatr(1:15, ielem)
!    print*,'new amatr',ielem, x5
!
enddo !...(2)ie = 1,nelem
!
print*,'getamatr',amatr(:,1)
!
 end subroutine  getamatr_curv
!
!...subroutine: Calculate the mass matrix of cubic elements for DGMs, npoly=1 for DGP1, and npoly=2 for DGP2...
!
 subroutine  getamatr_curv2(amatr,geoel,coord,inpoel)
 use constant
 implicit none
!
 real*8,dimension(1:nmatr,1:nelem),intent(out)::amatr
 real*8,dimension(1:ngeel, 1:nelem), intent(inout)::geoel
 real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
!
 integer*4,dimension(1:nvtri, 1:nelem),intent(in)::inpoel
!
!...local array...
!
  real*8,dimension(1:2, 1:3)::xpin
  real*8,dimension(1:2, 1:nptri)::xp 
  real*8,dimension(1:nptri)::shp, dspr, dsps
  real*8:: weigh(ngausd_geo), posi(2, ngausd_geo)
  real*8:: x5(5,5), mmatr(5,5), b55(5)
!...local real number
  real*8:: dxdr,dxds,dydr,dyds 
  real*8:: eps,c00,c10,c05,c20
  real*8:: r, s, t, djaco, volel
  real*8:: wi, xc, yc, xgaus, ygaus, dxc, dyc
  real*8:: f1, f2, f3, f4
  real*8::f5,f6,f7,f8
  real*8::f22,f23,f24,f25,f26,f33,f34,f35,f36,f44,f45,f46,f55,f56,f66
  real*8::b2,b3,b4,b5,b6, djac
!
  integer::ielem, igaus, ishp
!
  data eps / 1.0d-06 / 
  data c00 / 0.0d0 / 
  data c10 / 1.0d0 / 
  data c05 / 0.5d0 / 
  data c20 / 2.0d0 / 
!
!...Find weight and position for gauss points...
  call rutope(2, ngausd_geo, posi, weigh)
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
!...
    xc    = geoel(1, ielem)
    yc    = geoel(2, ielem)
    volel = geoel(3, ielem)
    dxc   = geoel(4, ielem)*0.5d0
    dyc   = geoel(5, ielem)*0.5d0
!
!...
    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0
!
  do igaus =1, ngausd_geo
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
     xgaus = 0.d0
     ygaus = 0.d0     
!
     do ishp = 1, nptri  
      xgaus = xgaus + shp(ishp)*xp(1,ishp) 
      ygaus = ygaus + shp(ishp)*xp(2,ishp) 
     enddo  
!
     f1 = f1 + (xgaus-xc)/dxc*(xgaus-xc)/dxc/2.d0*djaco
     f2 = f2 + (xgaus-xc)/dxc*(ygaus-yc)/dyc*djaco
     f3 = f3 + (ygaus-yc)/dyc*(ygaus-yc)/dyc/2.d0*djaco
! 
  enddo
!
    geoel(6, ielem) = f1/volel
    geoel(7, ielem) = f3/volel
    geoel(8, ielem) = f2/volel
!
!    print*,'old geoel',ielem, geoel(11:13, ielem)
!    print*,'new geoel',ielem, f1/volel, f3/volel, f2/volel
! 
 enddo !...(1)ie = 1,nelem
!
!...2st loop to find amatr...
!
   do ielem = 1,nelem !...(2)ie = 1,nelem
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
!...Cell center...
    xc    = geoel(1, ielem)
    yc    = geoel(2, ielem)
    volel = geoel(3, ielem)
    dxc   = geoel(4, ielem)*0.5d0
    dyc   = geoel(5, ielem)*0.5d0
!
!...
    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0    
    f4 = 0.d0

    f22 = 0.d0
    f23 = 0.d0
    f24 = 0.d0
    f25 = 0.d0
    f26 = 0.d0


    f33 = 0.d0
    f34 = 0.d0
    f35 = 0.d0
    f36 = 0.d0

    f44 = 0.d0
    f45 = 0.d0
    f46 = 0.d0

    f55 = 0.d0
    f56 = 0.d0

    f66 = 0.d0

   do igaus =1, ngausd_geo
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
     xgaus = 0.d0
     ygaus = 0.d0     
!
     do ishp = 1, nptri  
      xgaus = xgaus + shp(ishp)*xp(1,ishp) 
      ygaus = ygaus + shp(ishp)*xp(2,ishp) 
     enddo  
       
     b2 = (xgaus-xc)/dxc
     b3 = (ygaus-yc)/dyc
     b4 = 0.5d0*b2**2-geoel(6,ielem)
     b5 = 0.5d0*b3**2-geoel(7,ielem)
     b6 = b2*b3 - geoel(8,ielem)

     f1 = f1+(xgaus-xc)/dxc*(xgaus-xc)/dxc*djaco
     f2 = f2+(xgaus-xc)/dxc*(ygaus-yc)/dyc*djaco
     f3 = f3+(ygaus-yc)/dyc*(ygaus-yc)/dyc*djaco

     if(npoly==2)then
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
  enddo

  if(npoly==1)then
!
     djac = f1*f3-f2**2    
     amatr(1,ielem) = f3/djac
     amatr(2,ielem) = -f2/djac
     amatr(3,ielem) = f1/djac 
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
     
  endif
!
!    print*,'old amatr',ielem, amatr(1:15, ielem)
!    print*,'new amatr',ielem, x5
!
 enddo !...(2)ie = 1,nelem

 end subroutine  getamatr_curv2
!
!...Subroutine: Check the mesh quality (triangle)...
!
subroutine getmeshquality(inpoel,coord,geoel)
use constant
implicit none
integer*4,dimension(1:nvtri,1:ntria)::inpoel
real*8,dimension(1:ndimn,1:npoin)::coord
real*8,dimension(1:ngeel,1:nsize)::geoel
integer:: ie,ip1,ip2,ip3,ip
real*8::xp(2,3)
real*8::meshrat(ncell)
real*8::l1,l2,l3,arg1,arg2,arg3,arg4,ratio,rincir,routcir,volel
!
!...Initialization...
!
  meshrat = 0.d0

  do ie =1 ,nelem
!
     ip1=inpoel(1,ie)
     ip2=inpoel(2,ie)
     ip3=inpoel(3,ie) 
!
     xp(1:2,1) = coord(1:2,ip1)
     xp(1:2,2) = coord(1:2,ip2) 
     xp(1:2,3) = coord(1:2,ip3)   
!
     volel = geoel(3, ie)
!...Length of 3 sides...
     l1 = sqrt((xp(1,1)-xp(1,2))**2 + (xp(2,1)-xp(2,2))**2 )
     l2 = sqrt((xp(1,3)-xp(1,2))**2 + (xp(2,3)-xp(2,2))**2 )
     l3 = sqrt((xp(1,1)-xp(1,3))**2 + (xp(2,1)-xp(2,3))**2 )    
!
     rincir= volel*2.d0/(l1+l2+l3)
     arg1 = l1+l2+l3
     arg2 = l2+l3-l1
     arg3 = l1+l3-l2
     arg4 = l1+l2-l3
     routcir = l1*l2*l3/sqrt(arg1*arg2*arg3*arg4)
!...
     ratio = routcir/rincir
!
!...find the difference between this triangle and the perfect triagnle..
!     
    meshrat(ie) = abs((ratio-2.d0)/2.d0)
!    if(meshrat(ie).gt.5.0) print*,'ie',ie
  enddo
!
!...Output tecplot format...
!
  open(11,file='meshquality.dat')
  write(11,*)'VARIABLES = "X","Y","R"'
  write(11,*)'ZONE'
  write(11,*)'     T="2D"'
  write(11,*)'     DATAPACKING=BLOCK'
  write(11,*)'     ZONETYPE=FEtriangle'
  write(11,*)'     N=',npoin
  write(11,*)'     E=',nelem
  write(11,*)'     VARLOCATION=([3]=CELLCENTERED)'

  do ip =1,npoin
     write(11,*)coord(1,ip)
  enddo
  do ip=1,npoin
    write(11,*)coord(2,ip)
  enddo
  do ie =1,nelem
    write(11,*)meshrat(ie)
  enddo
    write(11,*)
  do ie=1,nelem
    write(11,*)inpoel(1:3,ie) 
  enddo
 close(11)

end subroutine getmeshquality
!
!...Find curved face center in geofa...
!
 subroutine getgeofa(intfac, geofa, coord)
 use constant
 implicit none
 integer*4,dimension(1:nifai,1:nafac),intent(in)::intfac
!
  real*8,dimension(1:ngefa,1:nafac), intent(inout)::geofa
  real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
!
!...local array...
!
  real*8,dimension(1:2, 1:2)::xpin
  real*8,dimension(1:2, 1:nptfa)::xp 
  real*8,dimension(1:nptfa)::shp, dshpr
  real*8:: weigh(ngausf_geo), posi(1, ngausf_geo)
  real*8:: r, djaco, dxdr, dydr, volel, l1
  real*8:: wi, xcent, ycent
!
  integer::iface, igaus, ishp
!
!...Find weight and position for gauss points...
  call rutope(1, ngausf_geo, posi, weigh)
!  
  do iface = 1, nafac
!
  if(ncurv==0)then
  xpin(1, 1:2) = coord(1, intfac(3:4, iface)) !...One face constitutes of 'nptfa' points...
  xpin(2, 1:2) = coord(2, intfac(3:4, iface))
!
  l1=sqrt((xpin(2, 2)-xpin(2, 1))**2+(xpin(1, 2)-xpin(1, 1))**2)
  geofa(1,iface)= (xpin(2, 2)-xpin(2, 1))/l1 !...The normal vector
  geofa(2,iface)=-(xpin(1, 2)-xpin(1, 1))/l1 
  geofa(5,iface)= l1
!
  xp(1, 1:2) = xpin(1, 1:2)
  xp(2, 1:2) = xpin(2, 1:2)
  xp(1, 3) = 0.5d0*(xpin(1, 1) + xpin(1, 2))
  xp(2, 3) = 0.5d0*(xpin(2, 1) + xpin(2, 2))
 elseif(ncurv==1)then
  xp(1, 1:nptfa) = coord(1, intfac(3:(2+nptfa), iface))  
  xp(2, 1:nptfa) = coord(2, intfac(3:(2+nptfa), iface)) 
!
 endif
!
  volel = 0.d0
  djaco = 0.d0
  xcent = 0.d0
  ycent = 0.d0
! 
   do igaus = 1, ngausf_geo
    
    r  = posi(1, igaus) 
    wi = weigh(igaus)
!
    shp(1) =  -0.5d0*(1.d0-r)*r  
    shp(2) =   0.5d0*(1.d0+r)*r  
    shp(3) =         (1.d0+r)*(1.d0-r) 
!    shp(1) =   0.5d0*(1.d0-r)  
!    shp(2) =   0.5d0*(1.d0+r)  
!
    dshpr(1) = -0.5d0 + r
    dshpr(2) =  0.5d0 + r
    dshpr(3) = -2.d0*r
!    dshpr(1) = -0.5d0  
!    dshpr(2) =  0.5d0 
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
    do ishp = 1, nptfa
     xcent = xcent + shp(ishp)*xp(1, ishp)*djaco*wi
     ycent = ycent + shp(ishp)*xp(2, ishp)*djaco*wi
    enddo
!
    volel = volel + djaco*wi
!
   enddo   
!
   geofa(3, iface) = xcent/volel
   geofa(4, iface) = ycent/volel
!   print*,'old geofa',iface, geofa(3:4, iface), geofa(5,iface)
!   print*,'new geofa',iface, xcent/volel, ycent/volel, volel 
!
  enddo

 end subroutine getgeofa
!
!...Find face center in geoel for curved cell...
!
 subroutine getgeoel(inpoel, iptri, ipqua, geoel, coord)
 use constant
 implicit none
 integer*4,dimension(1:nvtri, 1:ntria),intent(in)::inpoel
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
  real*8:: wi, xcent, ycent
!
  integer::ie
  integer::ielem, igaus, ishp
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
!...1st loop to find center and volulme for triangles...
!
 do ie = 1, ntria !...(1)ie = 1,nelem
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
    xcent = 0.d0
    ycent = 0.d0
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
     do ishp = 1, nptri  
      xcent = xcent + shp(ishp)*xp(1,ishp)*djaco
      ycent = ycent + shp(ishp)*xp(2,ishp)*djaco   
     enddo  
! 
  enddo
!
    geoel(1, ielem) = xcent/volel
    geoel(2, ielem) = ycent/volel
    geoel(3, ielem) = volel
!
    geoel(4, ielem) = maxval(xp(1, 1:3))-minval(xp(1, 1:3))
    geoel(5, ielem) = maxval(xp(2, 1:3))-minval(xp(2, 1:3))
!
!    print*,'old geoel',ielem, geoel(1:5, ielem)
!    print*,'new geoel',ielem, xcent/volel, ycent/volel, volel, maxval(xp(1, 1:3))-minval(xp(1, 1:3)), &
!                         maxval(xp(2, 1:3))- minval(xp(2, 1:3))
! 
 enddo !...(1)ie = 1,nelem
!
!...2nd loop to find center and volulme for quads...
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
xcent = 0.d0
ycent = 0.d0
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
do ishp = 1, npqua
xcent = xcent + shpq(ishp)*xpq(1,ishp)*djaco
ycent = ycent + shpq(ishp)*xpq(2,ishp)*djaco
enddo
!
enddo
!
geoel(1, ielem) = xcent/volel
geoel(2, ielem) = ycent/volel
geoel(3, ielem) = volel
!
geoel(4, ielem) = maxval(xpq(1, 1:4))-minval(xpq(1, 1:4))
geoel(5, ielem) = maxval(xpq(2, 1:4))-minval(xpq(2, 1:4))
!
!    print*,'old geoel',ielem, geoel(1:5, ielem)
!    print*,'new geoel',ielem, xcent/volel, ycent/volel, volel, maxval(xp(1, 1:3))-minval(xp(1, 1:3)), &
!                         maxval(xp(2, 1:3))- minval(xp(2, 1:3))
!
 enddo !...(1)ie = 1,nelem

 end subroutine getgeoel
!
!...Find face center in geoel for cubic cell...
!
 subroutine getgeoel_curv2(inpoel, geoel, coord)
 use constant
 implicit none
 integer*4,dimension(1:nvtri, 1:nelem),intent(in)::inpoel
!
  real*8,dimension(1:ngeel, 1:nelem), intent(inout)::geoel
  real*8,dimension(1:ndimn,1:npoin),  intent(in)::coord
!
!...local array...
!
  real*8,dimension(1:2, 1:3)::xpin
  real*8,dimension(1:2, 1:nptri)::xp 
  real*8,dimension(1:nptri)::shp, dspr, dsps
  real*8:: weigh(ngausd_geo), posi(2, ngausd_geo)
!...local real number
  real*8:: dxdr,dxds,dydr,dyds 
  real*8:: eps,c00,c10,c05,c20
  real*8:: r, s, t, djaco, volel
  real*8:: wi, xcent, ycent
!
  integer::ielem, igaus, ishp
!
  data eps / 1.0d-06 / 
  data c00 / 0.0d0 / 
  data c10 / 1.0d0 / 
  data c05 / 0.5d0 / 
  data c20 / 2.0d0 / 
!
!...Find weight and position for gauss points...
  call rutope(2, ngausd_geo, posi, weigh)
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
!...
    xcent = 0.d0
    ycent = 0.d0
    volel = 0.d0
!
  do igaus =1,ngausd_geo
!
     r  = posi(1,igaus)
     s  = posi(2,igaus)
     t = 1.d0-r-s
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
     volel = volel + djaco
!     
     do ishp = 1, nptri  
      xcent = xcent + shp(ishp)*xp(1,ishp)*djaco
      ycent = ycent + shp(ishp)*xp(2,ishp)*djaco   
     enddo  
! 
  enddo
!
    geoel(1, ielem) = xcent/volel
    geoel(2, ielem) = ycent/volel
    geoel(3, ielem) = volel
!
    geoel(4, ielem) = maxval(xp(1, 1:3))-minval(xp(1, 1:3))
    geoel(5, ielem) = maxval(xp(2, 1:3))-minval(xp(2, 1:3))
!
!    print*,'old geoel',ielem, geoel(1:5, ielem)
!    print*,'new geoel',ielem, xcent/volel, ycent/volel, volel, maxval(xp(1, 1:3))-minval(xp(1, 1:3)), &
!                         maxval(xp(2, 1:3))- minval(xp(2, 1:3))
! 
 enddo !...(1)ie = 1,nelem

 end subroutine getgeoel_curv2
!
!...Read weight matrix from input file
!
subroutine  getWeight_ML
use constant
real*8, dimension(25,8)::weigh0
real*8, dimension(8, 6)::weigh1
real*8, dimension(6, 4)::weigh2
real*8, dimension(4, 1)::weigh3
!
integer::imt, jmt

!...Read weigh0
open(12,file='w0out.txt')
do imt = 1, 25
read(12,*) weigh0(imt,1:8)
enddo
close(12)

!...Read weigh1
open(12,file='w1out.txt')
do imt = 1, 8
read(12,*) weigh1(imt,1:6)
enddo
close(12)

!...Read weigh2
open(12,file='w2out.txt')
do imt = 1, 6
read(12,*) weigh2(imt,1:4)
enddo
close(12)

!...Read weigh3
open(12,file='w3out.txt')
do imt = 1, 4
read(12,*) weigh3(imt,1)
enddo
close(12)

!...Matrix transpose
weigh0t = transpose(weigh0)
weigh1t = transpose(weigh1)
weigh2t = transpose(weigh2)
weigh3t = transpose(weigh3)

end subroutine  getWeight_ML
