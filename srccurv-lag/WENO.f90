!
!...WENO stnecils
!
subroutine weno_stncl_quadp1(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoph,&
coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:9,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)

integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(5)::jelaj
integer,dimension(9)::ielaj
integer,  dimension(1:nvqua):: ipq
integer::ielst1(4,6),ielst2(4,4)
integer::idxiv(3)
!
real*8,dimension(1:2, 1:nvqua)::xpq,xpqj
real*8, allocatable::matst(:, :), mtsti(:, :), binv(:), rhs(:,:), vcrhs(:,:), mtsym(:,:),rhssym(:,:)

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dxc, dyc, xmc, ymc, xcj, ycj
real*8:: dxcrt, dycrt,dxj,dyj, dyrt, dxrt
real*8:: detma
!
integer:: ie,ies, jelem,  ishp, ieaj, istor,iv,nvaj, ivaj, id, ifaj
integer:: idegr, is, iq
integer:: idxod

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Allocate memory
allocate(matst(2, 2), mtsti(2, 2), binv(2), rhs(nq, 2) ,vcrhs(2, 2))
allocate(mtsym(2, 2), rhssym(nq, 2))
!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Scaling parameter
dxc = maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Adjacent cells
ielaj = 0
ielaj(1) = ielem
ielaj(2:5) = esqua(1:4, ielem)

!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
do iv = 1, 4

jelaj = 0
nvaj = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielem)then
nvaj = nvaj +1
jelaj(nvaj) = jelem
endif
enddo !istor
!...Stencil No.
idxod = 1
idxiv = 0
do ivaj=1,3
do ifaj=1,4
if(jelaj(ivaj).eq.ielaj(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
!if(ielem.eq.1)print*,'ielem',iv,nvaj,ivaj,ifaj,jelaj(ivaj),ielaj(ifaj+1),idxod
endif
enddo
enddo

!...Order the adjacent cells
do ivaj=1,3
if(idxiv(ivaj) == 0)then
if(idxod.eq.2)then
ielaj(6) = jelaj(ivaj)
elseif(idxod.eq.6)then
ielaj(7) = jelaj(ivaj)
elseif(idxod.eq.12)then
ielaj(8) = jelaj(ivaj)
elseif(idxod.eq.4)then
ielaj(9) = jelaj(ivaj)
else
!print*,'Wrong stencil!',ielem,idxod
!stop
endif
endif
enddo
!
enddo !do iv = 1, 4

!if(ielem.eq.26)print*,'ielem',ielem,iv,nvaj,ielaj(1:9)

!
if(minval(ielaj(1:9)).eq.0)then
isten = 1

else
!...Part 2: Eliminate cells without 9 adjacient cells
ielst1(1, 1:6) =(/1,2,5,6,9,8/)
ielst1(2, 1:6) =(/1,3,2,7,6,9/)
ielst1(3, 1:6) =(/1,4,3,8,7,6/)
ielst1(4, 1:6) =(/1,5,4,9,8,7/)
!
ielst2(1, 1:4) =(/1,2,6,9/)
ielst2(2, 1:4) =(/1,3,7,6/)
ielst2(3, 1:4) =(/1,4,8,7/)
ielst2(4, 1:4) =(/1,5,9,8/)

!...Zero out other stencils
unkpe(2:ndegr,:,2:nsten) = 0.d0

!...1st stencil
do is =1 ,4
!
isten = isten +1
!
do ies =1, 2
!
jelem = ielaj(ielst1(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dxcrt = (xcj - xmc)/dxc
dycrt = (ycj - ymc)/dyc
!
dxrt = dxj/dxc
dyrt = dyj/dyc

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)

matst(ies, 1)  = (xcj-xmc)/dxc
matst(ies, 2)  = (ycj-ymc)/dyc
!

enddo

!...Symmetrize square matrix
call  matrix_sym(2, 2, nq, matst, rhs, mtsym, rhssym)

!...Invert matrix


detma = mtsym(1, 1)*mtsym(2, 2) - mtsym(1, 2)*mtsym(2, 1)
mtsti(1, 1) = mtsym(2, 2)
mtsti(1, 2) =-mtsym(1, 2)
mtsti(2, 1) =-mtsym(2, 1)
mtsti(2, 2) = mtsym(1, 1)
!
 mtsti =mtsti/detma
!else
!binv = 0.d0
!mtsti = 0.d0
!call getinvmat(ndegr-1, matst, mtsti, binv)
!endif


!...Update the stencil polynomial
do idegr = 2, 3
do ies   = 1, 2
unkpe(idegr, 1:nq, is+1) =  unkpe(idegr, 1:nq, is+1) + mtsti(idegr-1, ies)*rhssym(1:nq, ies)
enddo
enddo
!
!if(ielem.eq.1369)print*,'stencil', is,mtsti,rhs(1, 1:2)
!
enddo
endif !...sencils
!!...Release memory
deallocate(matst, mtsti, binv,mtsym,rhssym)

return
end subroutine weno_stncl_quadp1
!
!...Non least-squares Hermite WENO P2 stnecil, in contrast with subroutine weno_stncl_quadv2,
!...store required quantity in the array geoph to simplify the integral process.
!
subroutine weno_stncl_quad(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoel,geoph,&
coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:ngeel,1:nsize),        intent(in) ::geoel
real*8,dimension(1:9,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)

integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(5)::jelaj
integer,dimension(9)::ielaj
integer,  dimension(1:nvqua):: ipq
integer::ielst1(4,6),ielst2(4,4)
integer::idxiv(3)
!
real*8,dimension(1:2, 1:nvqua)::xpq,xpqj
real*8, allocatable::matst(:, :), mtsti(:, :), binv(:), rhs(:,:), vcrhs(:,:), mtsym(:,:),rhssym(:,:)

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dxc, dyc, xmc, ymc, xcj, ycj
real*8:: dxcrt, dycrt,dxj,dyj, dyrt, dxrt
real*8:: detma
!
integer:: ie,ies, jelem,  ishp, ieaj, istor,iv,nvaj, ivaj, id, ifaj
integer:: idegr, is, iq
integer:: idxod

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Allocate memory
 allocate(matst(ndegr-1, ndegr-1), mtsti(ndegr-1, ndegr-1), binv(ndegr-1), rhs(nq, ndegr-1) ,vcrhs(ndegr-1, ndegr-1))
 allocate(mtsym(ndegr-1, ndegr-1), rhssym(nq, ndegr-1))
!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Scaling parameter
dxc = sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Adjacent cells
ielaj = 0
ielaj(1) = ielem
ielaj(2:5) = esqua(1:4, ielem)

!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
do iv = 1, 4

jelaj = 0
nvaj = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielem)then
nvaj = nvaj +1
jelaj(nvaj) = jelem
endif
enddo !istor
!...Stencil No.
idxod = 1
idxiv = 0
do ivaj=1,3
do ifaj=1,4
if(jelaj(ivaj).eq.ielaj(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
!if(ielem.eq.1)print*,'ielem',iv,nvaj,ivaj,ifaj,jelaj(ivaj),ielaj(ifaj+1),idxod
endif
enddo
enddo

!...Order the adjacent cells
do ivaj=1,3
if(idxiv(ivaj) == 0)then
if(idxod.eq.2)then
ielaj(6) = jelaj(ivaj)
elseif(idxod.eq.6)then
ielaj(7) = jelaj(ivaj)
elseif(idxod.eq.12)then
ielaj(8) = jelaj(ivaj)
elseif(idxod.eq.4)then
ielaj(9) = jelaj(ivaj)
else
!print*,'Wrong stencil!',ielem,idxod
!stop
endif
endif
enddo
!
enddo !do iv = 1, 4

!if(ielem.eq.1153)print*,'ielem',ielem,iv,nvaj,ielaj(1:9)

!
if(minval(ielaj(1:9)).eq.0)then
 isten = 1

else
!...Part 2: Eliminate cells without 9 adjacient cells
ielst1(1, 1:6) =(/1,2,5,6,9,8/)
ielst1(2, 1:6) =(/1,3,2,7,6,9/)
ielst1(3, 1:6) =(/1,4,3,8,7,6/)
ielst1(4, 1:6) =(/1,5,4,9,8,7/)
!
ielst2(1, 1:4) =(/1,2,6,9/)
ielst2(2, 1:4) =(/1,3,7,6/)
ielst2(3, 1:4) =(/1,4,8,7/)
ielst2(4, 1:4) =(/1,5,9,8/)

!...Zero out other stencils
unkpe(2:ndegr,:,2:nsten) = 0.d0

!...1st stencil
do is =1 ,4
!
  isten = isten +1
!
do ies =1, ndegr-1
!
jelem = ielaj(ielst1(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = sqrt(geoel(3, jelem))!maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = dxj!maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dxcrt = (xcj - xmc)/dxc
dycrt = (ycj - ymc)/dyc
!
dxrt = dxj/dxc
dyrt = dyj/dyc

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)

matst(ies, 1)  = (xcj-xmc)/dxc
matst(ies, 2)  = (ycj-ymc)/dyc
!
if(npoly.eq.2)then
matst(ies, 3)  = dxrt**2*geoph(3, jelem)   + 0.5d0*dxcrt**2 - geoph(3, ielem)
matst(ies, 4)  = dyrt**2*geoph(4, jelem)   + 0.5d0*dycrt**2 - geoph(4, ielem)
matst(ies, 5)  = dxrt*dyrt*geoph(5, jelem) +   dxcrt*dycrt  - geoph(5, ielem)
endif

!if(ielem.eq.1153)then
! print*,'weno',ielem,dxc,dyc,dxj,dyj
! print*,'wen2',is,xcj,ycj,xmc,ymc
! print*,'wen3',ies,dxc,dyc,dxj,dyj
! print*,'wen4',ies,geoph(3:5, jelem),geoph(3:5, ielem)
! print*,'wen5',ies,dxrt,dyrt,dxcrt,dycrt
!endif

enddo

!...Symmetrize square matrix
!  call  matrix_sym(ndegr-1, ndegr-1, nq, matst, rhs, mtsym, rhssym)

!binv = 0.d0
!mtsti = 0.d0
!call getinvmat(ndegr-1, mtsym, mtsti, binv)

call matinv225(matst,mtsti,ndegr-1)

!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
 unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhs(1:nq, ies)
enddo
enddo
!
!if(ielem.eq.1248)print*,'matr',ielem,is+1,unkpe(2:3, 1:nq, is+1)
!print*,'mat2',ielem,unkpe(2:ndegr,1,is+1)


enddo

!...2nd stencil
!...The 2st 3 Lagrangian polynomial
if(npoly.eq.-2)then

do is =1 ,4
!...Increasing stencils
 isten = isten +1
!
do ies=1, 3
!
jelem = ielaj(ielst2(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dxcrt = (xcj - xmc)/dxc
dycrt = (ycj - ymc)/dyc

dxrt = dxj/dxc
dyrt = dyj/dyc

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)

matst(ies, 1)  = (xcj-xmc)/dxc
matst(ies, 2)  = (ycj-ymc)/dyc
!
if(npoly.eq.2)then
matst(ies, 3)  = dxrt**2*geoph(3, jelem)  + 0.5d0*dxcrt**2 - geoph(3, ielem)
matst(ies, 4)  = dyrt**2*geoph(4, jelem)  + 0.5d0*dycrt**2 - geoph(4, ielem)
matst(ies, 5)  = dxrt*dyrt*geoph(5, jelem)+   dxcrt*dycrt  - geoph(5, ielem)
endif

enddo

!...The 2 Hermite polynomials from the 2nd stencil
rhs(1:nq, 4:5) = 0.d0

!
jelem = ielaj(ielst2(is, 2))

!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))

!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dxcrt = (xcj - xmc)/dxc
dycrt = (ycj - ymc)/dyc
!
dxrt = dxj/dxc
dyrt = dyj/dyc
!
vcrhs(4, 1) = 2.d0*geoph(3, jelem)
vcrhs(4, 2) =      geoph(5, jelem)
vcrhs(4, 3) = 3.d0*geoph(6, jelem)
vcrhs(4, 4) =      geoph(9, jelem)
vcrhs(4, 5) = 2.d0*geoph(8, jelem)

!
vcrhs(5, 1) =      geoph(5, jelem)
vcrhs(5, 2) = 2.d0*geoph(4, jelem)
vcrhs(5, 3) =      geoph(8, jelem)
vcrhs(5, 4) = 3.d0*geoph(7, jelem)
vcrhs(5, 5) = 2.d0*geoph(9, jelem)
!
matst(4, 1)  = 2.d0*dxrt*geoph(3, jelem)
matst(4, 2)  = dyrt*geoph(5, jelem)
matst(4, 3)  = 3.d0*(dxrt**2)*geoph(6, jelem) + 2.d0*dxrt*dxcrt*geoph(3, jelem)
matst(4, 4)  = (dyrt)**2*geoph(9, jelem) +    dyrt*dycrt*geoph(5, jelem)
matst(4, 5)  = 2.d0*dxrt*dyrt*geoph(8, jelem) + dyrt*dxcrt*geoph(5, jelem) +&
2.d0*dycrt*dxrt*geoph(3, jelem)

!
matst(5, 1)  = dxrt*geoph(5, jelem)
matst(5, 2)  = 2.d0*dyrt*geoph(4, jelem)
matst(5, 3)  = dxrt**2*geoph(8, jelem) + dxrt*dxcrt*geoph(5, jelem)
matst(5, 4)  = 3.d0*(dyrt**2)*geoph(7, jelem) +  2.d0*dyrt*dycrt*geoph(4, jelem)
matst(5, 5)  = 2.d0*dxrt*dyrt*geoph(9, jelem) +  2.d0*dyrt*dxcrt*geoph(4, jelem) +&
dycrt*dxrt*geoph(5, jelem)

do ies = 4, 5
do idegr = 2, ndegr
rhs(1:nq, ies) = rhs(1:nq, ies) + unknp(idegr, 1:nq, jelem)*vcrhs(ies, idegr-1)
enddo
enddo

!...Symmetrize square matrix
call  matrix_sym(ndegr-1, ndegr-1, nq, matst, rhs, mtsym, rhssym)

!...Invert matrix
binv = 0.d0
mtsti = 0.d0
call getinvmat(ndegr-1, mtsym, mtsti, binv)


!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
   unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhssym(1:nq, ies)
enddo
enddo

enddo

endif !if(npoly.eq.2)then

endif !...sencils
!!...Release memory
deallocate(matst, mtsti, binv,mtsym,rhssym)

return
end subroutine weno_stncl_quad
!
!...conventional WENO limiter on characteristic field...
!
subroutine weno_char_quad(ipqua, esqua, unkno, geoel, coord, cooro, esuv1, esuv2)
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
integer, parameter:: nsten=13,nfprj=10
integer:: ipq(nvqua)
integer:: mapfe(1:2,1:nfqua)
!
real*8,dimension(1:2, 1:nvqua)::xpq, xpqi, xpq_j, xpqi_j
real*8,dimension(1:ndegr,1:nq,1:nsize)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten)::unkf_phy,unksc, unksl,unksf_cha, unkf_cha
real*8,dimension(1:ndegr,1:nq)::unk_phy

real*8,dimension(ndegr-1, ndegr-1)::bb,bbq,bbi
real*8,dimension(1:9,1:nsize) ::geoph
real*8,dimension(1:nq,1:ndegr)::rhs
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
real*8::a11, a12, a21, a22, matra,matrb
real*8::ai11, ai12, ai21, ai22
real*8::dnxl,dnyl,rhoxi,rhoet,etxi,etet,uvel,vvel
real*8::ulxi,ulet,vlxi,vlet
real*8::mapmt(2,2,ncell),mapmt2(2,2,ncell)
integer::jelaj(3)
real*8::lamda1,lamda2,lmat1,lmat2,mapt,mapd,matrc,matrd
real*8::unknc(4)
real*8::rhomc,uc,vc,etc,dxc,dyc,dxj,dyj
real*8::os1(nsten)
real*8::ctsft(5), smthid(5)
real*8, allocatable::bbqi(:,:), mbbq(:, :), mbbi(:, :)
real*8, allocatable::binv(:)
!
real*8:: c00, c05, c10, c20, c16,epsil
real*8:: dr, ds, rc, sc, r, s
real*8:: xmc, ymc, xmc_aj, ymc_aj
real*8:: masel, volel
real*8:: xcrho, ycrho
real*8:: xgaus, ygaus, xgausi, ygausi
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
integer:: ie, ies, iq, is, isten, ntrsf, ielem, jelem, ifa, igaus, ishp,ifprj,neaj,ieaj,istor,iv,id,ib
integer:: idegr,jdegr,ijmat,im,jm,jb
integer:: nweno
real*8::  rpowe
!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c16   / 0.16666666666666666d0/
data c20   / 2.0d0    /
data epsil / 1.0d-6   /
data rpowe / .5d0    /   ! works for lilia case
!
if(npoly.eq.1)then
 allocate(bbqi(2,2), binv(2), mbbq(3, nsize), mbbi(4, nsize))
elseif(npoly==2)then
 allocate(bbqi(5,5), binv(5), mbbq(15, nsize), mbbi(25, nsize))

endif
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

!...I.2: Shift terms for DG(P2)
if(npoly.eq.2)then
!
geoph(3:9, :) = 0.d0
!
do ie = 1, nquad
!
ielem = ie + ntria
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
dxc = sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

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
bqp(1) = 1.d0
bqp(2) = (xgaus-xmc)/dxc
bqp(3) = (ygaus-ymc)/dyc

!...Matrix
geoph(3, ielem) = geoph(3, ielem) + 0.5d0*rhogi*djacoi*bqp(2)*bqp(2)
geoph(4, ielem) = geoph(4, ielem) + 0.5d0*rhogi*djacoi*bqp(3)*bqp(3)
geoph(5, ielem) = geoph(5, ielem) +       rhogi*djacoi*bqp(2)*bqp(3)
geoph(6, ielem) = geoph(6, ielem) +   c16*rhogi*djacoi*bqp(2)**3
geoph(7, ielem) = geoph(7, ielem) +   c16*rhogi*djacoi*bqp(3)**3
geoph(8, ielem) = geoph(8, ielem) + 0.5d0*rhogi*djacoi*bqp(2)**2*bqp(3)
geoph(9, ielem) = geoph(9, ielem) + 0.5d0*rhogi*djacoi*bqp(3)**2*bqp(2)
!
enddo

!...Shift terms for P2
geoph(3:9, ielem)= geoph(3:9, ielem)/masel

enddo !do ie = 1, nquad

endif !if(npoly.eq.2)then

!...I.3 Mass matrix
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
dxc = sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

!...mass
masel = geoel(4, ielem)

!...Initial zero
bb = 0.d0
bbq= 0.d0
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
bqp(2) = (xgaus-xmc)/dxc
bqp(3) = (ygaus-ymc)/dyc
!
if(npoly.eq.2)then
bq(4) = 0.5d0*bq(2)*bq(2) - geoel(19, ielem)
bq(5) = 0.5d0*bq(3)*bq(3) - geoel(20, ielem)
bq(6) =       bq(2)*bq(3) - geoel(21, ielem)

bqp(4) = 0.5d0*bqp(2)*bqp(2) - geoph(3, ielem)
bqp(5) = 0.5d0*bqp(3)*bqp(3) - geoph(4, ielem)
bqp(6) =       bqp(2)*bqp(3) - geoph(5, ielem)
endif

!...Matrix(5x5)
do idegr =2, ndegr
do jdegr =2, ndegr
bb(idegr-1, jdegr-1)  =  bb(idegr-1, jdegr-1) + rhogi*djacoi*bq(idegr)*bqp(jdegr)
bbq(idegr-1, jdegr-1) = bbq(idegr-1, jdegr-1) + rhogi*djacoi*bqp(idegr)*bqp(jdegr)
enddo
enddo
!
enddo

!...RHS
rhs = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
rhs(1:nq, idegr) = rhs(1:nq, idegr) + unkno(jdegr, 1:nq, ielem)*bb(jdegr-1, idegr-1)
enddo
enddo

!...Store the array for later use
ijmat = 0
do im = 2, ndegr
do jm = im, ndegr
ijmat = ijmat +1
mbbq(ijmat, ielem) = bbq(im-1, jm-1)
enddo
enddo

!...Inverse of bbq
if(npoly.eq.1)then
detma = bbq(1, 1)*bbq(2, 2)-bbq(1, 2)*bbq(2, 1)
bbqi(1, 1) = bbq(2, 2)
bbqi(1, 2) =-bbq(1, 2)
bbqi(2, 1) =-bbq(2, 1)
bbqi(2, 2) = bbq(1, 1)

bbqi = bbqi/detma
!
detma = bb(1, 1)*bb(2, 2)-bb(1, 2)*bb(2, 1)
bbi(1, 1) = bb(2, 2)
bbi(1, 2) =-bb(1, 2)
bbi(2, 1) =-bb(2, 1)
bbi(2, 2) = bb(1, 1)

bbi = bbi/detma
elseif(npoly.eq.2)then
!binv = 0.d0
!bbqi = 0.d0
!call getinvmat(5, bbq, bbqi, binv)

!binv = 0.d0
!bbi = 0.d0
!call getinvmat(5, bb, bbi, binv)
endif

!
call matinv225(bbq, bbqi, ndegr-1)
call matinv225(bb, bbi, ndegr-1)

!...LHS
unknp(:, 1:nq, ielem) = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
unknp(idegr, 1:nq, ielem) = unknp(idegr, 1:nq, ielem) + rhs(1:nq, jdegr)*bbqi(jdegr-1,idegr-1)
enddo
enddo

!...Cell average on physical domain
unknp(1, 1:nq, ielem) = unkno(1, 1:nq, ielem)

!....
!if(ielem.eq.1153)then
! print*,'bad2',ielem,unknp(1:ndegr,1, ielem)
! print*,'bad3',ielem,unkno(1:ndegr,1,ielem)
!endif

!...Store the array for later use
ijmat = 0
do im = 2, ndegr
do jm = 2, ndegr
ijmat = ijmat +1
mbbi(ijmat, ielem) = bbi(im-1, jm-1)
enddo
enddo

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
!matra = unkno(1, 2, ielem)
!matrb = unkno(1, 3, ielem)
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

!...No local system
mapmt(1, 1, ielem) = 1.d0
mapmt(1, 2, ielem) = 0.d0
mapmt(2, 1, ielem) = 0.d0
mapmt(2, 2, ielem) = 1.d0

enddo !do ie = 1, nquad
!
!...Part II: Recontruction
!
do ie = 1, nquad
!
ielem = ie + ntria
!
if(geoel(10, ielem).lt.10) cycle
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...physical mass center
xmc =geoph(1, ielem)
ymc =geoph(2, ielem)

!...scale
dxc = sqrt(geoel(3,ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))
volel = geoel(3, ielem)

!...High-order shift terms
ctsft(1) = geoph(3, ielem)
ctsft(2) = geoph(4, ielem)
ctsft(3) = geoph(5, ielem)
ctsft(4) = geoph(6, ielem)
ctsft(5) = geoph(7, ielem)

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
xpq_j(1, 1:4) = coord(1, ipqua(1:4,jelem-ntria))
xpq_j(2, 1:4) = coord(2, ipqua(1:4,jelem-ntria))
!
dxj = sqrt(geoel(3, jelem))!maxval(xpq_j(1, 1:4)) - minval(xpq_j(1, 1:4))
dyj = dxj!maxval(xpq_j(2, 1:4)) - minval(xpq_j(2, 1:4))
!
isten = isten + 1
unksc(1, 1:nq, isten) = unknp(1, 1:nq, ielem)
unksc(2, 1:nq, isten) = unknp(2, 1:nq, jelem)*dxc/dxj
unksc(3, 1:nq, isten) = unknp(3, 1:nq, jelem)*dyc/dyj
!
 if(npoly.eq.2)then
  unksc(4, 1:nq, isten) = unknp(4, 1:nq, jelem)*dxc**2/dxj**2
  unksc(5, 1:nq, isten) = unknp(5, 1:nq, jelem)*dyc**2/dyj**2
  unksc(6, 1:nq, isten) = unknp(6, 1:nq, jelem)*dyc*dxc/dyj/dxj
 endif
endif

enddo !do ies = 1, 4

Case (2)

!...Weno
isten = 1
unksc(:, 1:nq, isten) = unknp(:, 1:nq, ielem)

!...WENO
!call weno_stncl_quad(isten, nsten, ielem, ipqua, esqua, unksc, unknp, geoel,geoph,&
!coord, cooro, esuv1, esuv2)

call weno_stncl_quad_Taylor(isten, nsten, ielem, ipqua, esqua, unksc, unknp, geoel,geoph,&
coord, cooro, esuv1, esuv2, mapmt(:,:,ielem), ctsft)

end select


!...Local system
call getunkno_local2(0, 1, mapmt(:,:,ielem), unksc(:,:,1))
!
!if(ielem.eq.265.or.ielem.eq.266)then
!do is=1, isten
!  print*,'unkno',ielem,is,unksc(2:ndegr,1,is)
!enddo
!endif

!...Update local stencile values
unksl = unksc

!...  b. curvatures for the face-neighboring cells
ifprj = 0
unkf_cha = 0.d0
unkf_phy = 0.d0

!...Loop over the faces
do ies = 1, 4

jelem = esqua(ies,ie)
if(jelem .le. ncell) then

!...Find the projection vector
ifprj = ifprj +1

!...normal vector only using linear face even for curved meshes
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
call getmatrix_prj2(qmat, qinvm, dnx, dny, unknl,unknr,ielem)

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
!
!call getvector_local(mapmt(:,:,ielem), dxc,dyc)

!...Smooth indicator for every stencle
if(npoly.eq.1)then
 do is= 1, isten
 do iq= 1, nq
    os(iq, is) = ((unksf_cha(2, iq, is)/dxc)**2 +&
                  (unksf_cha(3, iq, is)/dyc)**2)!/(abs(unksf_cha(1, iq, 1))+epsil)**1
 enddo
 enddo
elseif(npoly.eq.2)then
 do is= 1, isten
 do iq= 1, nq
!...
    smthid(1) = 1.d0/(dxc**2)*(unksf_cha(2, iq, is)**2 +&
                 2.d0*ctsft(1)*unksf_cha(4, iq, is)**2 +&
                 2.d0*ctsft(2)*unksf_cha(6, iq, is)**2 +&
                 2.d0*ctsft(3)*unksf_cha(4, iq, is)*unksf_cha(6, iq, is))

    smthid(2) = 1.d0/(dyc**2)*(unksf_cha(3, iq, is)**2  +&
                 2.d0*ctsft(2)*unksf_cha(5, iq, is)**2 +&
                 2.d0*ctsft(1)*unksf_cha(6, iq, is)**2 +&
                 2.d0*ctsft(3)*unksf_cha(5, iq, is)*unksf_cha(6, iq, is))

    smthid(3) = 1.d0/(dxc**4)*unksf_cha(4, iq, is)**2
    smthid(4) = 1.d0/(dyc**4)*unksf_cha(5, iq, is)**2
    smthid(5) = 1.d0/(dyc*dxc)**2*unksf_cha(6, iq, is)**2
!
    os(iq, is) = (smthid(1) + smthid(2)) +&
                 (smthid(3) + smthid(4) + smthid(5))*volel
 enddo
 enddo

endif

!...Mometum
!os1(:) = 0.25d0*(os(1, :)+os(2, :)+os(3, :)+os(4, :))
os1(:) = max(os(1, :),os(2, :),os(3, :),os(4, :))


do is= 1, isten
os(1:nq, is) = os1(is)
enddo

!...Linear weight for central and biased stencle
weigl(1)=.2d0;
weigl(2:isten)= 0.8d0/(isten-1.d0)


weigt = 0.d0
do is= 1, isten
do iq =1,nq
weigh(iq, is) = weigl(is)/(epsil+os(iq,is))**rpowe
enddo
enddo

!if(ielem.eq.71)then
!print*,'ielem12',ielem,qmat
!print*,'ielem13',ifprj,weigh(:,1)
!endif
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
unk_phy(2:ndegr, iq) =unk_phy(2:ndegr, iq) + weige(iq, ifa)*unkf_phy(2:ndegr,iq,ifa)
enddo
enddo

!...Local system
call getunkno_local2(1, 1, mapmt(:,:,ielem), unk_phy)

!...Transfer back the reference configuration
ijmat = 0
do im = 2, ndegr
do jm = im, ndegr
ijmat = ijmat +1
bbq(im-1, jm-1) = mbbq(ijmat, ielem)
bbq(jm-1, im-1) = mbbq(ijmat, ielem)
enddo
enddo

ijmat = 0
do im = 2, ndegr
do jm = 2, ndegr
ijmat = ijmat +1
bbi(im-1, jm-1) = mbbi(ijmat, ielem)
enddo
enddo

!...RHS
rhs = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
rhs(1:nq, idegr) = rhs(1:nq, idegr) + unk_phy(jdegr, 1:nq)*bbq(jdegr-1, idegr-1)
enddo
enddo
!
!...LHS
unkno(2:ndegr, 1:nq, ielem) = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
unkno(idegr, 1:nq, ielem) = unkno(idegr, 1:nq, ielem) + rhs(1:nq, jdegr)*bbi(jdegr-1, idegr-1)
enddo
enddo

!...  end of the loop over the quad cells
enddo
!...Release memory
deallocate(bbqi, binv, mbbq, mbbi)
return
end subroutine weno_char_quad
!...WENO stnecils
!
subroutine weno_stncl_quad_orthop1(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoel, geoph, geopj,&
coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:ngeel,1:nsize),        intent(in) ::geoel
real*8,dimension(1:16,1:nsize),        intent(in) ::geoph
real*8,dimension(1:3, 1:nsize),        intent(in) ::geopj
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)

integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(5)::jelaj
integer,dimension(9)::ielaj
integer,  dimension(1:nvqua):: ipq
integer::ielst1(4,6),ielst2(4,4)
integer::idxiv(3)
!
real*8,dimension(1:2, 1:nvqua)::xpq,xpqj
real*8, allocatable::matst(:, :), mtsti(:, :), binv(:), rhs(:,:), vcrhs(:,:), mtsym(:,:),rhssym(:,:)

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dxc, dyc, xmc, ymc, xcj, ycj
real*8:: dxcrt, dycrt,dxj,dyj, dyrt, dxrt
real*8:: dai, daj, bqxj, bqyj, dart
real*8:: detma
!
integer:: ie,ies, jelem,  ishp, ieaj, istor,iv,nvaj, ivaj, id, ifaj
integer:: idegr, is, iq
integer:: idxod

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Allocate memory
allocate(matst(ndegr-1, ndegr-1), mtsti(ndegr-1, ndegr-1), binv(ndegr-1), rhs(nq, ndegr-1) ,vcrhs(ndegr-1, ndegr-1))
allocate(mtsym(ndegr-1, ndegr-1), rhssym(nq, ndegr-1))
!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Scaling parameter
dxc = sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Adjacent cells
ielaj = 0
ielaj(1) = ielem
ielaj(2:5) = esqua(1:4, ielem)

!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
do iv = 1, 4

jelaj = 0
nvaj = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielem)then
nvaj = nvaj +1
jelaj(nvaj) = jelem
endif
enddo !istor
!...Stencil No.
idxod = 1
idxiv = 0
do ivaj=1,3
do ifaj=1,4
if(jelaj(ivaj).eq.ielaj(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
!if(ielem.eq.1)print*,'ielem',iv,nvaj,ivaj,ifaj,jelaj(ivaj),ielaj(ifaj+1),idxod
endif
enddo
enddo

!...Order the adjacent cells
do ivaj=1,3
if(idxiv(ivaj) == 0)then
if(idxod.eq.2)then
ielaj(6) = jelaj(ivaj)
elseif(idxod.eq.6)then
ielaj(7) = jelaj(ivaj)
elseif(idxod.eq.12)then
ielaj(8) = jelaj(ivaj)
elseif(idxod.eq.4)then
ielaj(9) = jelaj(ivaj)
else
!print*,'Wrong stencil!',ielem,idxod
!stop
endif
endif
enddo
!
enddo !do iv = 1, 4

!if(ielem.eq.1153)print*,'ielem',ielem,iv,nvaj,ielaj(1:9)

!
if(minval(ielaj(1:9)).eq.0)then
isten = 1

else
!...Part 2: Eliminate cells without 9 adjacient cells
ielst1(1, 1:6) =(/1,2,5,6,9,8/)
ielst1(2, 1:6) =(/1,3,2,7,6,9/)
ielst1(3, 1:6) =(/1,4,3,8,7,6/)
ielst1(4, 1:6) =(/1,5,4,9,8,7/)
!
ielst2(1, 1:4) =(/1,2,6,9/)
ielst2(2, 1:4) =(/1,3,7,6/)
ielst2(3, 1:4) =(/1,4,8,7/)
ielst2(4, 1:4) =(/1,5,9,8/)

!...Zero out other stencils
unkpe(2:ndegr,:,2:nsten) = 0.d0

!...1st stencil
do is =1 ,4
!
isten = isten +1
!
do ies =1, ndegr-1
!
jelem = ielaj(ielst1(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = sqrt(geoel(3, jelem)) !maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = dxj !maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dai  = geoel(3, ielem)
daj  = geoel(3, jelem)
dart = daj/dai

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)
!
bqxj = (xcj-xmc)/dxc
bqyj = (ycj-ymc)/dyc

matst(ies, 1)  = bqxj
matst(ies, 2)  = geoph(3, ielem)*bqxj + bqyj !+ geoph(4, ielem)

!if(ielem.eq.1153)then
! print*,'weno',ielem,dxc,dyc,dxj,dyj
! print*,'wen2',is,xcj,ycj,xmc,ymc
! print*,'wen3',ies,dxc,dyc,dxj,dyj
! print*,'wen4',ies,geoph(3:5, jelem),geoph(3:5, ielem)
! print*,'wen5',ies,dxrt,dyrt,dxcrt,dycrt
!endif

enddo

!...Invert matrix

if(npoly.eq.1)then
detma = matst(1, 1)*matst(2, 2) - matst(1, 2)*matst(2, 1)
mtsti(1, 1) = matst(2, 2)
mtsti(1, 2) =-matst(1, 2)
mtsti(2, 1) =-matst(2, 1)
mtsti(2, 2) = matst(1, 1)
!
mtsti = mtsti/detma
endif


!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhs(1:nq, ies)
enddo
enddo
!
!if(ielem.eq.1248)print*,'matr',ielem,is+1,unkpe(2:3, 1:nq, is+1)
!print*,'mat2',ielem,unkpe(2:ndegr,1,is+1)


enddo
endif
!!...Release memory
deallocate(matst, mtsti, binv,mtsym,rhssym)

return
end subroutine weno_stncl_quad_orthop1
!
!...WENO stnecils
!
subroutine weno_stncl_quad_ortho(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoel, geoph, geopj,&
coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:ngeel,1:nsize),        intent(in) ::geoel
real*8,dimension(1:16,1:nsize),        intent(in) ::geoph
real*8,dimension(1:3, 1:nsize),        intent(in) ::geopj
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)

integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(5)::jelaj
integer,dimension(9)::ielaj
integer,  dimension(1:nvqua):: ipq
integer::ielst1(4,7),ielst2(4,4)
integer::idxiv(3)
!
real*8,dimension(1:2, 1:nvqua)::xpq,xpqj
real*8, allocatable::matst(:, :), mtsti(:, :), binv(:), rhs(:,:), vcrhs(:,:), mtsym(:,:),rhssym(:,:)

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dxc, dyc, xmc, ymc, xcj, ycj
real*8:: dxcrt, dycrt,dxj,dyj, dyrt, dxrt
real*8:: dai, daj, bqxj, bqyj, dart
real*8:: detma
!
integer:: ie,ies, jelem,  ishp, ieaj, istor,iv,nvaj, ivaj, id, ifaj
integer:: idegr, is, iq
integer:: idxod

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Allocate memory
allocate(matst(ndegr-1, ndegr-1), mtsti(ndegr-1, ndegr-1), binv(ndegr-1), rhs(nq, ndegr-1) ,vcrhs(ndegr-1, ndegr-1))
allocate(mtsym(ndegr-1, ndegr-1), rhssym(nq, ndegr-1))
!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Scaling parameter
dxc = sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Adjacent cells
ielaj = 0
ielaj(1) = ielem
ielaj(2:5) = esqua(1:4, ielem)

!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
do iv = 1, 4

jelaj = 0
nvaj = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielem)then
nvaj = nvaj +1
jelaj(nvaj) = jelem
endif
enddo !istor
!...Stencil No.
idxod = 1
idxiv = 0
do ivaj=1,3
do ifaj=1,4
if(jelaj(ivaj).eq.ielaj(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
!if(ielem.eq.1)print*,'ielem',iv,nvaj,ivaj,ifaj,jelaj(ivaj),ielaj(ifaj+1),idxod
endif
enddo
enddo

!...Order the adjacent cells
do ivaj=1,3
if(idxiv(ivaj) == 0)then
if(idxod.eq.2)then
ielaj(6) = jelaj(ivaj)
elseif(idxod.eq.6)then
ielaj(7) = jelaj(ivaj)
elseif(idxod.eq.12)then
ielaj(8) = jelaj(ivaj)
elseif(idxod.eq.4)then
ielaj(9) = jelaj(ivaj)
else
!print*,'Wrong stencil!',ielem,idxod
!stop
endif
endif
enddo
!
enddo !do iv = 1, 4

!if(ielem.eq.1153)print*,'ielem',ielem,iv,nvaj,ielaj(1:9)

!
if(minval(ielaj(1:9)).eq.0)then
isten = 1

else
!...Part 2: Eliminate cells without 9 adjacient cells
ielst1(1, 1:7) =(/1,2,5,6,9,8, 3/)
ielst1(2, 1:7) =(/1,3,2,7,6,9, 5/)
ielst1(3, 1:7) =(/1,4,3,8,7,6,5/)
ielst1(4, 1:7) =(/1,5,4,9,8,7,3/)
!
ielst2(1, 1:4) =(/1,2,6,9/)
ielst2(2, 1:4) =(/1,3,7,6/)
ielst2(3, 1:4) =(/1,4,8,7/)
ielst2(4, 1:4) =(/1,5,9,8/)

!...Zero out other stencils
unkpe(2:ndegr,:,2:nsten) = 0.d0

!...1st stencil
do is =1 ,4
!
isten = isten +1
!
do ies =1, ndegr-1
!
jelem = ielaj(ielst1(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = sqrt(geoel(3, jelem)) !maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = dxj !maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dai  = geoel(3, ielem)
daj  = geoel(3, jelem)
dart = daj/dai

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)
!
bqxj = (xcj-xmc)/dxc
bqyj = (ycj-ymc)/dyc

matst(ies, 1)  = bqxj
matst(ies, 2)  = geoph(3, ielem)*bqxj + bqyj + geoph(4, ielem)
!
if(npoly.eq.2)then
matst(ies, 3)  = dart*geopj(1, jelem) + bqxj**2  + geoph(5, ielem)*bqxj +  geoph(6, ielem)*bqyj + geoph(7, ielem)
matst(ies, 4)  = dart*geopj(1, jelem)*geoph(8, ielem)   + bqxj**2*geoph(8, ielem) +&
                 dart*geopj(3, jelem) + bqxj*bqyj + geoph(9, ielem)*bqxj + geoph(10, ielem)*bqyj + geoph(11, ielem)
matst(ies, 5)  = dart*geopj(1, jelem)*geoph(12, ielem) +   bqxj**2*geoph(12, ielem) +&
                 geoph(13, ielem)*dart*geopj(3, jelem) + geoph(13, ielem)*bqxj*bqyj +&
                 dart*geopj(2, jelem) + bqyj**2 + geoph(14, ielem)*bqxj + geoph(15, ielem)*bqyj + geoph(16, ielem)
endif

!...Weighted
!
!rhs(1:nq, ies) = rhs(1:nq, ies)/sqrt(bqxj**2+bqyj**2)
!matst(ies, 1:5) =  matst(ies, 1:5)/sqrt(bqxj**2+bqyj**2)

!if(ielem.eq.88)then
! print*,'weno',ielem,dxc,dyc,dxj,dyj,sqrt(bqxj**2+bqyj**2)
! print*,'wen2',is,rhs(1:nq, ies)
! print*,'wen3',ies,dxc,dyc,dxj,dyj
! print*,'wen4',ies,geoph(3:5, jelem),geoph(3:5, ielem)
! print*,'wen5',ies,dxrt,dyrt,dxcrt,dycrt
!endif

enddo

!...Symmetrize square matrix
!call  matrix_sym(ndegr-1, ndegr-1, nq, matst, rhs, mtsym, rhssym)

!...Invert small matrix directly
call Matinv225(matst, mtsti, ndegr-1)

!
!do ies =1, ndegr-1
!print*,'matr',ielem,matst(ies,1:ndegr-1)
!enddo

!
!if(ielem.eq.1170)then
!do ies =1, ndegr-1
!print*,'masm',ielem,mtsym(1:ndegr-1 ,ies)
!enddo
!endif
!print*,'mtsm',ielem,mtsym
!...Invert matrix

!if(npoly.eq.1)then
! detma = matst(1, 1)*matst(2, 2) - matst(1, 2)*matst(2, 1)
! mtsti(1, 1) = matst(2, 2)
! mtsti(1, 2) =-matst(1, 2)
! mtsti(2, 1) =-matst(2, 1)
! mtsti(2, 2) = matst(1, 1)
!
! mtsti = mtsti/detma
!else
!binv = 0.d0
!mtsti = 0.d0
!call getinvmat(ndegr-1, matst, mtsti, binv)
!endif


!...Update the stencil polynomial
!
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhs(1:nq, ies)
enddo
enddo
!
!if(ielem.eq.1248)print*,'matr',ielem,is+1,unkpe(2:3, 1:nq, is+1)
!print*,'mat2',ielem,unkpe(2:ndegr,1,is+1)


enddo

!...2nd stencil
!...The 2st 3 Lagrangian polynomial
if(npoly.eq.-2)then

do is =1 ,4
!...Increasing stencils
isten = isten +1
!
do ies=1, 3
!
jelem = ielaj(ielst2(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dxcrt = (xcj - xmc)/dxc
dycrt = (ycj - ymc)/dyc

dxrt = dxj/dxc
dyrt = dyj/dyc

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)

matst(ies, 1)  = (xcj-xmc)/dxc
matst(ies, 2)  = (ycj-ymc)/dyc
!
if(npoly.eq.2)then
matst(ies, 3)  = dxrt**2*geoph(3, jelem)  + 0.5d0*dxcrt**2 - geoph(3, ielem)
matst(ies, 4)  = dyrt**2*geoph(4, jelem)  + 0.5d0*dycrt**2 - geoph(4, ielem)
matst(ies, 5)  = dxrt*dyrt*geoph(5, jelem)+   dxcrt*dycrt  - geoph(5, ielem)
endif

enddo

!...The 2 Hermite polynomials from the 2nd stencil
rhs(1:nq, 4:5) = 0.d0

!
jelem = ielaj(ielst2(is, 2))

!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))

!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dxcrt = (xcj - xmc)/dxc
dycrt = (ycj - ymc)/dyc
!
dxrt = dxj/dxc
dyrt = dyj/dyc
!
vcrhs(4, 1) = 2.d0*geoph(3, jelem)
vcrhs(4, 2) =      geoph(5, jelem)
vcrhs(4, 3) = 3.d0*geoph(6, jelem)
vcrhs(4, 4) =      geoph(9, jelem)
vcrhs(4, 5) = 2.d0*geoph(8, jelem)

!
vcrhs(5, 1) =      geoph(5, jelem)
vcrhs(5, 2) = 2.d0*geoph(4, jelem)
vcrhs(5, 3) =      geoph(8, jelem)
vcrhs(5, 4) = 3.d0*geoph(7, jelem)
vcrhs(5, 5) = 2.d0*geoph(9, jelem)
!
matst(4, 1)  = 2.d0*dxrt*geoph(3, jelem)
matst(4, 2)  = dyrt*geoph(5, jelem)
matst(4, 3)  = 3.d0*(dxrt**2)*geoph(6, jelem) + 2.d0*dxrt*dxcrt*geoph(3, jelem)
matst(4, 4)  = (dyrt)**2*geoph(9, jelem) +    dyrt*dycrt*geoph(5, jelem)
matst(4, 5)  = 2.d0*dxrt*dyrt*geoph(8, jelem) + dyrt*dxcrt*geoph(5, jelem) +&
2.d0*dycrt*dxrt*geoph(3, jelem)

!
matst(5, 1)  = dxrt*geoph(5, jelem)
matst(5, 2)  = 2.d0*dyrt*geoph(4, jelem)
matst(5, 3)  = dxrt**2*geoph(8, jelem) + dxrt*dxcrt*geoph(5, jelem)
matst(5, 4)  = 3.d0*(dyrt**2)*geoph(7, jelem) +  2.d0*dyrt*dycrt*geoph(4, jelem)
matst(5, 5)  = 2.d0*dxrt*dyrt*geoph(9, jelem) +  2.d0*dyrt*dxcrt*geoph(4, jelem) +&
dycrt*dxrt*geoph(5, jelem)

do ies = 4, 5
do idegr = 2, ndegr
rhs(1:nq, ies) = rhs(1:nq, ies) + unknp(idegr, 1:nq, jelem)*vcrhs(ies, idegr-1)
enddo
enddo

!...Symmetrize square matrix
call  matrix_sym(ndegr-1, ndegr-1, nq, matst, rhs, mtsym, rhssym)

!...Invert matrix
binv = 0.d0
mtsti = 0.d0
call getinvmat(ndegr-1, mtsym, mtsti, binv)


!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhssym(1:nq, ies)
enddo
enddo

enddo

endif !if(npoly.eq.2)then

endif !...sencils
!!...Release memory
deallocate(matst, mtsti, binv,mtsym,rhssym)

return
end subroutine weno_stncl_quad_ortho
!
!...WENO limiter on characteristic field using the local orthogonal basis...
!
subroutine weno_char_quad_orthgp1(ipqua, esqua, unkno, geoel, coord, cooro, esuv1, esuv2)
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
real*8,dimension(1:2, 1:nvqua)::xpq, xpqi, xpq_j, xpqi_j
real*8,dimension(1:ndegr,1:nq,1:nsize)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten)::unkf_phy,unksc, unksl,unksf_cha, unkf_cha
real*8,dimension(1:ndegr,1:nq)::unk_phy

real*8,dimension(ndegr-1, ndegr-1)::bb,bbq,bbi
real*8,dimension(1:16,1:nsize) ::geoph
real*8,dimension(1:3,1:nsize) ::geopj
real*8,dimension(1:nq,1:ndegr)::rhs
real*8,dimension(1:4, 1:4)::qmat, qinvm
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:ndegr)::bq,bqp, bqp_aj, bqo
real*8::xpf(1:2, 1:2)
real*8, dimension(1:nq)::unknl,unknr
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: wi
real*8   weigh(1:nq, nsten)
real*8:: weigt(1:nq)
real*8   os(1:nq, nsten)
real*8:: weigl(nsten)
real*8:: weige(nq, nfprj)
real*8::a11, a12, a21, a22, matra,matrb
real*8::ai11, ai12, ai21, ai22
real*8::dnxl,dnyl,rhoxi,rhoet,etxi,etet,uvel,vvel
real*8::ulxi,ulet,vlxi,vlet
real*8::mapmt(2,2,ncell)
integer::jelaj(3)
real*8::lamda1,lamda2,lmat1,lmat2,mapt,mapd,matrc,matrd
real*8::unknc(4)
real*8::rhomc,uc,vc,etc,dxc,dyc,dxj,dyj
real*8::os1(nsten)
real*8::ctsft(5), smthid(5)
real*8, allocatable::bbqi(:,:), mbbq(:, :), mbbi(:, :),vcoef(:),cfmat(:,:), cfmti(:,:)
real*8, allocatable::binv(:), cfb(:)
!
real*8:: c00, c05, c10, c20, c16,c124,epsil
real*8:: dr, ds, rc, sc, r, s
real*8:: xmc, ymc, xmc_aj, ymc_aj
real*8:: masel, volel
real*8:: xcrho, ycrho
real*8:: xgaus, ygaus, xgausi, ygausi
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
integer:: ie, ies, iq, is, isten, ntrsf, ielem, jelem, ifa, igaus, ishp,ifprj,neaj,ieaj,istor,iv,id,ib
integer:: idegr,jdegr,ijmat,im,jm,jb
integer:: nweno
real*8::  rpowe
!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c16   / 0.16666666666666666d0/
data c124  / 0.04166666666666666d0/
data c20   / 2.0d0    /
data epsil / 1.0d-6   /
data rpowe / 2.d0    /   ! works for lilia case
!
if(npoly.eq.1)then
allocate(bbqi(2,2), binv(2), mbbq(3, nsize), mbbi(4, nsize), vcoef(12), cfmat(6, 6), cfb(6), cfmti(6,6))
elseif(npoly==2)then
allocate(bbqi(5,5), binv(5), mbbq(15, nsize), mbbi(25, nsize), vcoef(12), cfmat(6, 6), cfb(6), cfmti(6,6))
endif
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
!

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

!...I.2: Shift terms for DG(P2)
!if(npoly.eq.2)then

geoph(3:16,:) = 0.d0
geopj(1:3,:)  = 0.d0
!
do ie = 1, nquad
!
ielem = ie + ntria
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
!dxc =maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
!dyc =maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

dxc = sqrt(geoel(3, ielem))
dyc = dxc

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

!...mass
masel = geoel(4, ielem)

!...Zero out vcoef
vcoef=0.d0
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
bqp(1) = 1.d0
bqp(2) = (xgaus-xmc)/dxc
bqp(3) = (ygaus-ymc)/dyc

!...Matrix
vcoef(1)  = vcoef( 1) + rhogi*djacoi*bqp(2)*bqp(2)
vcoef(2)  = vcoef( 2) + rhogi*djacoi*bqp(3)*bqp(3)
vcoef(3)  = vcoef( 3) + rhogi*djacoi*bqp(2)*bqp(3)

!
enddo

!...Shift terms for P2
vcoef = vcoef/masel

!...Store some terms
geopj(1:3, ielem) =vcoef(1:3)


!masel = 1.d0
!...Solve the coefficients
!...a31 a32
geoph(3, ielem) =-vcoef(3)/vcoef(1)
geoph(4, ielem) = 0.d0

enddo !do ie = 1, nquad

!print*,'test1',ielem

!endif !if(npoly.eq.2)then

!...I.3 Mass matrix
do ie = 1,nquad
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
dxc =sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc =dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

!...mass
masel = geoel(4, ielem)

!...Initial zero
bb = 0.d0
bbq= 0.d0
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
bqp(2) = (xgaus-xmc)/dxc
bqp(3) = (ygaus-ymc)/dyc

!...Orthogonal basis functions
bqo(1) = 1.d0
bqo(2) = bqp(2)
bqo(3) = geoph(3, ielem)*bqp(2) + bqp(3) !+ geoph(4, ielem)

!...Matrix(5x5)
do idegr =2, ndegr
do jdegr =2, ndegr
bb(idegr-1, jdegr-1)  =  bb(idegr-1, jdegr-1) + rhogi*djacoi*bq(idegr)*bqo(jdegr)
bbq(idegr-1, jdegr-1) = bbq(idegr-1, jdegr-1) + rhogi*djacoi*bqo(idegr)*bqo(jdegr)
enddo
enddo
!
enddo

!...RHS
rhs = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
rhs(1:nq, idegr) = rhs(1:nq, idegr) + unkno(jdegr, 1:nq, ielem)*bb(jdegr-1, idegr-1)
enddo
enddo

!...Store the array for later use
ijmat = 0
do im = 2, ndegr
do jm = im, ndegr
ijmat = ijmat +1
mbbq(ijmat, ielem) = bbq(im-1, jm-1)
enddo
enddo

!...Check
!if(ielem.eq.1)print*,'mat',bbq

!...Inverse of bbq
if(npoly.eq.1)then
detma = bbq(1, 1)*bbq(2, 2)-bbq(1, 2)*bbq(2, 1)
bbqi(1, 1) = bbq(2, 2)
bbqi(1, 2) =-bbq(1, 2)
bbqi(2, 1) =-bbq(2, 1)
bbqi(2, 2) = bbq(1, 1)

bbqi = bbqi/detma
!
detma = bb(1, 1)*bb(2, 2)-bb(1, 2)*bb(2, 1)
bbi(1, 1) = bb(2, 2)
bbi(1, 2) =-bb(1, 2)
bbi(2, 1) =-bb(2, 1)
bbi(2, 2) = bb(1, 1)

bbi = bbi/detma

endif

!...LHS
unknp(:, 1:nq, ielem) = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
unknp(idegr, 1:nq, ielem) = unknp(idegr, 1:nq, ielem) + rhs(1:nq, jdegr)*bbqi(jdegr-1,idegr-1)
enddo
enddo

!...Cell average on physical domain
unknp(1, 1:nq, ielem) = unkno(1, 1:nq, ielem)

!....
!if(ielem.eq.1153)then
! print*,'bad2',ielem,unknp(1:ndegr,1, ielem)
! print*,'bad3',ielem,unkno(1:ndegr,1,ielem)
!endif

!...Store the array for later use
ijmat = 0
do im = 2, ndegr
do jm = 2, ndegr
ijmat = ijmat +1
mbbi(ijmat, ielem) = bbi(im-1, jm-1)
enddo
enddo
enddo !do ie = 1, nquad
!
!...Part II: Recontruction
!
do ie = 1,nquad
!
ielem = ie + ntria
!
!if(geoel(10, ielem).lt.10) cycle
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...physical mass center
xmc =geoph(1, ielem)
ymc =geoph(2, ielem)

!...scale
dxc = sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))
volel = geoel(3, ielem)

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
xpq_j(1, 1:4) = coord(1, ipqua(1:4,jelem-ntria))
xpq_j(2, 1:4) = coord(2, ipqua(1:4,jelem-ntria))
!
dxj = sqrt(geoel(3, jelem))!maxval(xpq_j(1, 1:4)) - minval(xpq_j(1, 1:4))
dyj = dxj!maxval(xpq_j(2, 1:4)) - minval(xpq_j(2, 1:4))
!
isten = isten + 1
unksc(1, 1:nq, isten) = unknp(1, 1:nq, ielem)
unksc(2, 1:nq, isten) = unknp(2, 1:nq, jelem)*dxc/dxj
unksc(3, 1:nq, isten) = unknp(3, 1:nq, jelem)*dyc/dyj
!
if(npoly.eq.2)then
unksc(4, 1:nq, isten) = unknp(4, 1:nq, jelem)*dxc**2/dxj**2
unksc(5, 1:nq, isten) = unknp(5, 1:nq, jelem)*dyc**2/dyj**2
unksc(6, 1:nq, isten) = unknp(6, 1:nq, jelem)*dyc*dxc/dyj/dxj
endif
endif

enddo !do ies = 1, 4

Case (2)

!...Weno
isten = 1
unksc(:, 1:nq, isten) = unknp(:, 1:nq, ielem)

!...WENO
call weno_stncl_quad_orthop1(isten, nsten, ielem, ipqua, esqua, unksc, unknp, geoel,geoph, geopj,&
coord, cooro, esuv1, esuv2)

end select

!...Update local stencile values
unksl = unksc

!...Degenearte to P1
!if(npoly.eq.2)unksl(4:ndegr,:,:) = 0.d0

!...  b. curvatures for the face-neighboring cells
ifprj = 0
unkf_cha = 0.d0
unkf_phy = 0.d0

!...Loop over the faces
do ies = 1, 4

jelem = esqua(ies,ie)
if(jelem .le. ncell) then

!...Find the projection vector
ifprj = ifprj +1

!...normal vector only using linear face even for curved meshes
xpf(1, 1:2) = xpq(1, mapfe(1:2, ies))
xpf(2, 1:2) = xpq(2, mapfe(1:2, ies))

dtx = xpf(1, 2) - xpf(1, 1)
dty = xpf(2, 2) - xpf(2, 1)

dlgt = sqrt(dtx**2 + dty**2)

dnx = dty/dlgt
dny =-dtx/dlgt
!
!call getvector_local(mapmt(:,:,ielem), dnx,dny)

!...left and right unknowns
unknl(1:nq) = unknp(1, 1:nq, ielem)
unknr(1:nq) = unknp(1, 1:nq, jelem)
!
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
call getmatrix_prj2(qmat, qinvm, dnx, dny, unknl,unknr,ielem)

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
!
!call getvector_local(mapmt(:,:,ielem), dxc,dyc)

!...Smooth indicator for every stencle
if(npoly.eq.1)then
do is= 1, isten
do iq= 1, nq
os(iq, is) = (((unksf_cha(2, iq, is)+geoph(3, ielem)*unksf_cha(3, iq, is))/dxc)**2 +&
(unksf_cha(3, iq, is)/dyc)**2)!/(abs(unksf_cha(1, iq, 1))+epsil)**1
enddo
enddo
elseif(npoly.eq.2)then
do is= 1, isten
do iq= 1, nq
!...
smthid(1) = 1.d0/(dxc**2)*(unksf_cha(2, iq, is)**2 +&
2.d0*ctsft(1)*unksf_cha(4, iq, is)**2 +&
2.d0*ctsft(2)*unksf_cha(6, iq, is)**2 +&
2.d0*ctsft(3)*unksf_cha(4, iq, is)*unksf_cha(6, iq, is))
!
smthid(2) = 1.d0/(dyc**2)*(unksf_cha(3, iq, is)**2 +&
2.d0*ctsft(2)*unksf_cha(5, iq, is)**2 +&
2.d0*ctsft(1)*unksf_cha(6, iq, is)**2 +&
2.d0*ctsft(3)*unksf_cha(5, iq, is)*unksf_cha(6, iq, is))
!
smthid(3) = 1.d0/(dxc**4)*unksf_cha(4, iq, is)**2
smthid(4) = 1.d0/(dyc**4)*unksf_cha(5, iq, is)**2
smthid(5) = 1.d0/(dyc*dxc)**2*unksf_cha(6, iq, is)**2
!
os(iq, is) = (smthid(1) + smthid(2)) +&
(smthid(3) + smthid(4) + smthid(5))*volel
enddo
enddo

endif

!...Mometum
!os1(:) = 0.25d0*(os(1, :)+os(2, :)+os(3, :)+os(4, :))
os1(:) = max(os(1, :),os(2, :),os(3, :),os(4, :))


do is= 1, isten
!os(1:nq, is) = os1(is)
enddo

!...Linear weight for central and biased stencle
weigl(1)=.9d0;
weigl(2:isten)= 0.1d0/(isten-1.d0)


weigt = 0.d0
do is= 1, isten
do iq =1,nq
weigh(iq, is) = weigl(is)/(epsil+os(iq,is))**rpowe
enddo
enddo

!if(ielem.eq.1)then
!print*,'ielem12',ielem,qinvm(1,:)
!print*,'ielem13',ifprj,os(1:4,5),epsil+os(1:4,5)
!endif
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
unk_phy(2:ndegr, iq) =unk_phy(2:ndegr, iq) + weige(iq, ifa)*unkf_phy(2:ndegr,iq,ifa)
enddo
enddo

!...Local system
!call getunkno_local(1, 1, mapmt(:,:,ielem), unk_phy)

!...Transfer back the reference configuration
ijmat = 0
do im = 2, ndegr
do jm = im, ndegr
ijmat = ijmat +1
bbq(im-1, jm-1) = mbbq(ijmat, ielem)
bbq(jm-1, im-1) = mbbq(ijmat, ielem)
enddo
enddo

ijmat = 0
do im = 2, ndegr
do jm = 2, ndegr
ijmat = ijmat +1
bbi(im-1, jm-1) = mbbi(ijmat, ielem)
enddo
enddo

!...RHS
rhs = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
rhs(1:nq, idegr) = rhs(1:nq, idegr) + unk_phy(jdegr, 1:nq)*bbq(jdegr-1, idegr-1)
enddo
enddo
!
!...LHS
unkno(2:ndegr, 1:nq, ielem) = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
unkno(idegr, 1:nq, ielem) = unkno(idegr, 1:nq, ielem) + rhs(1:nq, jdegr)*bbi(jdegr-1, idegr-1)
enddo
enddo
!
!if(npoly.eq.2)unkno(4:ndegr, 1:nq, ielem) = 0.d0
!
!if(ielem.eq.1)then
!print*,'ielem1',ielem,unkno(2:ndegr, 1, ielem)
!print*,'ielem2',ielem,bbi(1:ndegr-1,1),rhs(1, 2:ndegr)
!endif
!...  end of the loop over the quad cells
enddo

!stop
!...Release memory
deallocate(bbqi, binv, mbbq, mbbi,cfmat,cfmti,cfb)
return
end subroutine weno_char_quad_orthgp1
!
!...WENO stnecils
!
subroutine weno_stncl_quad_ortho2(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoel, geoph, geopj,&
coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:ngeel,1:nsize),        intent(in) ::geoel
real*8,dimension(1:16,1:nsize),        intent(in) ::geoph
real*8,dimension(1:3, 1:nsize),        intent(in) ::geopj
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)

integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(5)::jelaj
integer,dimension(9)::ielaj
integer,  dimension(1:nvqua):: ipq
integer::ielst1(4,7),ielst2(4,4)
integer::idxiv(3)

real*8::work(ndegr-1)
integer::ipiv(ndegr-1)
integer::info
!
real*8,dimension(1:2, 1:nvqua)::xpq,xpqj
real*8, allocatable::matst(:, :), mtsti(:, :), binv(:), rhs(:,:), vcrhs(:,:), mtsym(:,:),rhssym(:,:)

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dxc, dyc, xmc, ymc, xcj, ycj
real*8:: dxcrt, dycrt,dxj,dyj, dyrt, dxrt
real*8:: dai, daj, bqxj, bqyj, dart
real*8:: detma
!
integer:: ie,ies, jelem,  ishp, ieaj, istor,iv,nvaj, ivaj, id, ifaj
integer:: idegr, is, iq
integer:: idxod

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Allocate memory
allocate(matst(ndegr-1, ndegr-1), mtsti(ndegr-1, ndegr-1), binv(ndegr-1), rhs(nq, ndegr-1) ,vcrhs(ndegr-1, ndegr-1))
allocate(mtsym(ndegr-1, ndegr-1), rhssym(nq, ndegr-1))
!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Scaling parameter
dxc = sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Adjacent cells
ielaj = 0
ielaj(1) = ielem
ielaj(2:5) = esqua(1:4, ielem)

!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
do iv = 1, 4

jelaj = 0
nvaj = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielem)then
nvaj = nvaj +1
jelaj(nvaj) = jelem
endif
enddo !istor
!...Stencil No.
idxod = 1
idxiv = 0
do ivaj=1,3
do ifaj=1,4
if(jelaj(ivaj).eq.ielaj(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
!if(ielem.eq.1)print*,'ielem',iv,nvaj,ivaj,ifaj,jelaj(ivaj),ielaj(ifaj+1),idxod
endif
enddo
enddo

!...Order the adjacent cells
do ivaj=1,3
if(idxiv(ivaj) == 0)then
if(idxod.eq.2)then
ielaj(6) = jelaj(ivaj)
elseif(idxod.eq.6)then
ielaj(7) = jelaj(ivaj)
elseif(idxod.eq.12)then
ielaj(8) = jelaj(ivaj)
elseif(idxod.eq.4)then
ielaj(9) = jelaj(ivaj)
else
!print*,'Wrong stencil!',ielem,idxod
!stop
endif
endif
enddo
!
enddo !do iv = 1, 4

!if(ielem.eq.1153)print*,'ielem',ielem,iv,nvaj,ielaj(1:9)

!
if(minval(ielaj(1:9)).eq.0)then
isten = 1

else
!...Part 2: Eliminate cells without 9 adjacient cells
ielst1(1, 1:7) =(/1,2,5,6,9,8, 3/)
ielst1(2, 1:7) =(/1,3,2,7,6,9, 5/)
ielst1(3, 1:7) =(/1,4,3,8,7,6,5/)
ielst1(4, 1:7) =(/1,5,4,9,8,7,3/)
!
ielst2(1, 1:4) =(/1,2,6,9/)
ielst2(2, 1:4) =(/1,3,7,6/)
ielst2(3, 1:4) =(/1,4,8,7/)
ielst2(4, 1:4) =(/1,5,9,8/)

!...Zero out other stencils
unkpe(2:ndegr,:,2:nsten) = 0.d0

!...1st stencil
do is =1 ,4
!
isten = isten +1
!
do ies =1, ndegr-1
!
jelem = ielaj(ielst1(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = sqrt(geoel(3, jelem)) !maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = dxj !maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dai  = geoel(3, ielem)
daj  = geoel(3, jelem)
dart = daj/dai

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)
!
bqxj = (xcj-xmc)/dxc
bqyj = (ycj-ymc)/dyc

matst(ies, 1)  = bqxj
matst(ies, 2)  = geoph(3, ielem)*bqxj + bqyj + geoph(4, ielem)
!
if(npoly.eq.2)then
matst(ies, 3)  = dart*geopj(1, jelem) + bqxj**2  + geoph(5, ielem)*bqxj +  geoph(6, ielem)*bqyj + geoph(7, ielem)
matst(ies, 4)  = dart*geopj(1, jelem)*geoph(8, ielem)   + bqxj**2*geoph(8, ielem) +&
dart*geopj(3, jelem) + bqxj*bqyj + geoph(9, ielem)*bqxj + geoph(10, ielem)*bqyj + geoph(11, ielem)
matst(ies, 5)  = dart*geopj(1, jelem)*geoph(12, ielem) +   bqxj**2*geoph(12, ielem) +&
geoph(13, ielem)*dart*geopj(3, jelem) + geoph(13, ielem)*bqxj*bqyj +&
dart*geopj(2, jelem) + bqyj**2 + geoph(14, ielem)*bqxj + geoph(15, ielem)*bqyj + geoph(16, ielem)
endif

!if(ielem.eq.16.and.ies.eq.5)then
! print*,'weno',ielem,dxc,dyc,dxj,dyj
! print*,'wen2',is,xcj,ycj,xmc,ymc
! print*,'wen3',ies,geoel(3,jelem)
! print*,'wen4',ies,geoph(3:5, jelem),geoph(3:5, ielem)
! print*,'wen5',jelem,dxrt,dyrt,dxcrt,dycrt
!endif

enddo

!
!if(ielem.eq.7)then
!do ies =1, ndegr-1
!print*,'mat2',matst(ies,1:ndegr-1),ielaj(ielst1(1, 2:6))
!enddo
!endif
mtsym = matst
call matinv225(matst, mtsti, 5)
!call DGETRF(ndegr-1,ndegr-1,mtsym,ndegr-1,ipiv,info)
!call DGETRI(ndegr-1,mtsym,ndegr-1,ipiv,work,ndegr-1,info)
!...Symmetrize square matrix
!call  matrix_sym(ndegr-1, ndegr-1, nq, matst, rhs, mtsym, rhssym,ielem)

!if(ielem.eq.7)then
!do ies =1, ndegr-1
!print*,'mtsy',mtsym(ies,1:ndegr-1)
!enddo
!endif
!
!if(ielem.eq.1170)then
!do ies =1, ndegr-1
!print*,'masm',ielem,mtsym(1:ndegr-1 ,ies)
!enddo
!endif
!print*,'mtsm',ielem,mtsym
!...Invert matrix

!if(npoly.eq.1)then
! detma = matst(1, 1)*matst(2, 2) - matst(1, 2)*matst(2, 1)
! mtsti(1, 1) = matst(2, 2)
! mtsti(1, 2) =-matst(1, 2)
! mtsti(2, 1) =-matst(2, 1)
! mtsti(2, 2) = matst(1, 1)
!
! mtsti = mtsti/detma
!else
!binv = 0.d0
!mtsti = 0.d0
!call getinvmat(ndegr-1, matst, mtsti, binv)
!endif

!if(ielem.eq.7.and.isten.eq.2)then
!do ies =1, ndegr-1
!print*,'mtin',mtsti(ies,1:ndegr-1),rhs(3,ies)
!enddo
!endif

!...Update the stencil polynomial
!
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhs(1:nq, ies)
enddo
enddo
!
!if(ielem.eq.7)print*,'sten',jelem,isten,unkpe(2:ndegr, 2, isten)


enddo

!...2nd stencil
!...The 2st 3 Lagrangian polynomial
if(npoly.eq.-2)then

do is =1 ,4
!...Increasing stencils
isten = isten +1
!
do ies=1, 3
!
jelem = ielaj(ielst2(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dxcrt = (xcj - xmc)/dxc
dycrt = (ycj - ymc)/dyc

dxrt = dxj/dxc
dyrt = dyj/dyc

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)

matst(ies, 1)  = (xcj-xmc)/dxc
matst(ies, 2)  = (ycj-ymc)/dyc
!
if(npoly.eq.2)then
matst(ies, 3)  = dxrt**2*geoph(3, jelem)  + 0.5d0*dxcrt**2 - geoph(3, ielem)
matst(ies, 4)  = dyrt**2*geoph(4, jelem)  + 0.5d0*dycrt**2 - geoph(4, ielem)
matst(ies, 5)  = dxrt*dyrt*geoph(5, jelem)+   dxcrt*dycrt  - geoph(5, ielem)
endif

enddo

!...The 2 Hermite polynomials from the 2nd stencil
rhs(1:nq, 4:5) = 0.d0

!
jelem = ielaj(ielst2(is, 2))

!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))

!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dxcrt = (xcj - xmc)/dxc
dycrt = (ycj - ymc)/dyc
!
dxrt = dxj/dxc
dyrt = dyj/dyc
!
vcrhs(4, 1) = 2.d0*geoph(3, jelem)
vcrhs(4, 2) =      geoph(5, jelem)
vcrhs(4, 3) = 3.d0*geoph(6, jelem)
vcrhs(4, 4) =      geoph(9, jelem)
vcrhs(4, 5) = 2.d0*geoph(8, jelem)

!
vcrhs(5, 1) =      geoph(5, jelem)
vcrhs(5, 2) = 2.d0*geoph(4, jelem)
vcrhs(5, 3) =      geoph(8, jelem)
vcrhs(5, 4) = 3.d0*geoph(7, jelem)
vcrhs(5, 5) = 2.d0*geoph(9, jelem)
!
matst(4, 1)  = 2.d0*dxrt*geoph(3, jelem)
matst(4, 2)  = dyrt*geoph(5, jelem)
matst(4, 3)  = 3.d0*(dxrt**2)*geoph(6, jelem) + 2.d0*dxrt*dxcrt*geoph(3, jelem)
matst(4, 4)  = (dyrt)**2*geoph(9, jelem) +    dyrt*dycrt*geoph(5, jelem)
matst(4, 5)  = 2.d0*dxrt*dyrt*geoph(8, jelem) + dyrt*dxcrt*geoph(5, jelem) +&
2.d0*dycrt*dxrt*geoph(3, jelem)

!
matst(5, 1)  = dxrt*geoph(5, jelem)
matst(5, 2)  = 2.d0*dyrt*geoph(4, jelem)
matst(5, 3)  = dxrt**2*geoph(8, jelem) + dxrt*dxcrt*geoph(5, jelem)
matst(5, 4)  = 3.d0*(dyrt**2)*geoph(7, jelem) +  2.d0*dyrt*dycrt*geoph(4, jelem)
matst(5, 5)  = 2.d0*dxrt*dyrt*geoph(9, jelem) +  2.d0*dyrt*dxcrt*geoph(4, jelem) +&
dycrt*dxrt*geoph(5, jelem)

do ies = 4, 5
do idegr = 2, ndegr
rhs(1:nq, ies) = rhs(1:nq, ies) + unknp(idegr, 1:nq, jelem)*vcrhs(ies, idegr-1)
enddo
enddo

!...Symmetrize square matrix
call  matrix_sym(ndegr-1, ndegr-1, nq, matst, rhs, mtsym, rhssym)

!...Invert matrix
binv = 0.d0
mtsti = 0.d0
call getinvmat(ndegr-1, mtsym, mtsti, binv)


!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhssym(1:nq, ies)
enddo
enddo

enddo

endif !if(npoly.eq.2)then

endif !...sencils
!!...Release memory
deallocate(matst, mtsti, binv,mtsym,rhssym)

return
end subroutine weno_stncl_quad_ortho2
!
!...WENO limiter on characteristic field using the local orthogonal basis...
!
subroutine weno_char_quad_orthg2(ipqua, esqua, unkno, geoel, coord, cooro, esuv1, esuv2)
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
real*8,dimension(1:2, 1:nvqua)::xpq, xpqi, xpq_j, xpqi_j
real*8,dimension(1:ndegr,1:nq,1:nsize)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten)::unkf_phy,unksc, unksl,unksf_cha, unkf_cha
real*8,dimension(1:ndegr,1:nq)::unk_phy

real*8,dimension(ndegr-1, ndegr-1)::bb,bbq,bbi
real*8,dimension(1:16,1:nsize) ::geoph
real*8,dimension(1:3,1:nsize) ::geopj
real*8,dimension(1:nq,1:ndegr)::rhs
real*8,dimension(1:4, 1:4)::qmat, qinvm
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:ndegr)::bq,bqp, bqp_aj, bqo
real*8::xpf(1:2, 1:2)
real*8, dimension(1:nq)::unknl,unknr
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8:: wi
real*8   weigh(1:nq, nsten)
real*8:: weigt(1:nq)
real*8   os(1:nq, nsten)
real*8:: weigl(nsten)
real*8:: weige(nq, nfprj)
real*8::a11, a12, a21, a22, matra,matrb
real*8::ai11, ai12, ai21, ai22
real*8::dnxl,dnyl,rhoxi,rhoet,etxi,etet,uvel,vvel
real*8::ulxi,ulet,vlxi,vlet
real*8::mapmt(2,2,ncell)
integer::jelaj(3)
real*8::lamda1,lamda2,lmat1,lmat2,mapt,mapd,matrc,matrd
real*8::unknc(4)
real*8::rhomc,uc,vc,etc,dxc,dyc,dxj,dyj
real*8::os1(nsten)
real*8::ctsft(5), smthid(5)
real*8::work(ndegr-1)
real*8, allocatable::bbqi(:,:), mbbq(:, :), mbbi(:, :),vcoef(:),cfmat(:,:), cfmti(:,:)
real*8, allocatable::binv(:), cfb(:),cfmat3(:,:), cfmti3(:,:),binv3(:)
real*8, allocatable::cfmat4(:,:), cfmti4(:,:),binv4(:), cfb3(:), cfb4(:)
real*8, allocatable::rhs1(:,:),rhs1sym(:,:),bbsymi(:,:),bbsym(:,:)
real*8::cothg(6,6)
!
real*8:: c00, c05, c10, c20, c16,c124,epsil
real*8:: dr, ds, rc, sc, r, s
real*8:: xmc, ymc, xmc_aj, ymc_aj
real*8:: masel, volel
real*8:: xcrho, ycrho
real*8:: xgaus, ygaus, xgausi, ygausi
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
integer::info,ipiv(ndegr-1)
integer:: ie, ies, iq, is, isten, ntrsf, ielem, jelem, ifa, igaus, ishp,ifprj,neaj,ieaj,istor,iv,id,ib
integer:: idegr,jdegr,ijmat,im,jm,jb
integer:: nweno
real*8::  rpowe
!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c16   / 0.16666666666666666d0/
data c124  / 0.04166666666666666d0/
data c20   / 2.0d0    /
data epsil / 1.0d-6   /
data rpowe / 2.d0    /   ! works for lilia case
!
if(npoly.eq.1)then
allocate(bbqi(2,2), binv(2), mbbq(3, nsize), mbbi(4, nsize), vcoef(12), cfmat(6, 6), cfb(6), cfmti(6,6))
elseif(npoly==2)then
allocate(bbqi(5,5), binv(5), mbbq(15, nsize), mbbi(25, nsize), vcoef(12), cfmat(6, 6), cfb(6), cfmti(6,6))
allocate(cfmat3(3, 3), cfb3(3), cfmti3(3,3), binv3(3))
allocate(cfmat4(4, 4), cfb4(4), cfmti4(4,4), binv4(4))


endif

allocate(rhs1(1:nq,1:ndegr-1),rhs1sym(1:nq,1:ndegr-1),bbsymi(ndegr-1,ndegr-1),bbsym(ndegr-1,ndegr-1))
!
!...Part 0: Basis parameters setup
!
!...Specify the weno
nweno =2

!print*,'coord',coord(:,20)
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
!

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

!
enddo !do ie = 1, nquad

!...I.2: Shift terms for DG(P2)
!if(npoly.eq.2)then

geoph(3:16,:) = 0.d0
geopj(1:3,:)  = 0.d0
!
do ie = 1, nquad
!
ielem = ie + ntria
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
!dxc =maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
!dyc =maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

dxc = sqrt(geoel(3, ielem))
dyc = dxc

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

!...mass
masel = geoel(4, ielem)

!...Zero out vcoef
vcoef=0.d0
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
bqp(1) = 1.d0
bqp(2) = (xgaus-xmc)/dxc
bqp(3) = (ygaus-ymc)/dyc

!...Matrix
vcoef(1)  = vcoef( 1) + rhogi*djacoi*bqp(2)*bqp(2)
vcoef(2)  = vcoef( 2) + rhogi*djacoi*bqp(3)*bqp(3)
vcoef(3)  = vcoef( 3) + rhogi*djacoi*bqp(2)*bqp(3)
vcoef(4)  = vcoef( 4) + rhogi*djacoi*bqp(2)**3
vcoef(5)  = vcoef( 5) + rhogi*djacoi*bqp(3)**3
vcoef(6)  = vcoef( 6) + rhogi*djacoi*bqp(2)**2*bqp(3)
vcoef(7)  = vcoef( 7) + rhogi*djacoi*bqp(3)**2*bqp(2)
vcoef(8)  = vcoef( 8) + rhogi*djacoi*bqp(2)**4
vcoef(9)  = vcoef( 9) + rhogi*djacoi*bqp(3)**4
vcoef(10) = vcoef(10) + rhogi*djacoi*bqp(2)**3*bqp(3)
vcoef(11) = vcoef(11) + rhogi*djacoi*(bqp(2)*bqp(3))**2
vcoef(12) = vcoef(12) + rhogi*djacoi*bqp(3)**3*bqp(2)
!
enddo

!...Shift terms for P2
vcoef = vcoef/masel

!...Store some terms
geopj(1:3, ielem) =vcoef(1:3)


!masel = 1.d0
!...Solve the coefficients
!...a31 a32
geoph(3, ielem) =-vcoef(3)/vcoef(1)
geoph(4, ielem) = 0.d0

if(npoly.eq.2)then
!...a41, a42, a43
cfmat(1, 1) = vcoef(1); cfmat(1, 2) = vcoef(3);
cfmat(2, 1) = vcoef(3); cfmat(2, 2) = vcoef(2);

cfb(1) = -vcoef(4); cfb(2) =-vcoef(6)
!
detma = cfmat(1, 1)*cfmat(2, 2) - cfmat(2, 1)*cfmat(1, 2)

cfmti(1, 1) = vcoef(2); cfmti(1, 2) =-vcoef(3);
cfmti(2, 1) =-vcoef(3); cfmti(2, 2) = vcoef(1);

cfmti = cfmti/detma

geoph(5, ielem) = cfmti(1, 1)*cfb(1) + cfmti(1, 2)*cfb(2)
geoph(6, ielem) = cfmti(2, 1)*cfb(1) + cfmti(2, 2)*cfb(2)
geoph(7, ielem) =-vcoef(1)

!...a51, a52, a53, a54
cfmat3(1, 1) = vcoef(4); cfmat3(1, 2) = vcoef(1); cfmat3(1, 3) = vcoef(3);
cfmat3(2, 1) = vcoef(6); cfmat3(2, 2) = vcoef(3); cfmat3(2, 3) = vcoef(2);

cfmat3(3, 1) = vcoef(8) - vcoef(1)**2;
cfmat3(3, 2) = vcoef(4);
cfmat3(3, 3) = vcoef(6);

cfb3(1) = -vcoef(6); cfb3(2) =-vcoef(7); cfb3(3) =-vcoef(10) + vcoef(3)*vcoef(1)

!
binv3 = 0.d0
cfmti3 = 0.d0
!call getinvmat(3, cfmat3(1:3,1:3), cfmti3(1:3, 1:3), binv3(1:3))

call Matinv225(cfmat3,cfmti3,3)

!
!if(ielem.eq.1)print*,'ielem',
!
do im =1, 3
geoph(8, ielem) = geoph(8, ielem) + cfmti3(1, im)*cfb3(im)
geoph(9, ielem) = geoph(9, ielem) + cfmti3(2, im)*cfb3(im)
geoph(10, ielem)= geoph(10, ielem)+ cfmti3(3, im)*cfb3(im)
enddo

geoph(11, ielem)=-geoph(8, ielem)*vcoef(1) - vcoef(3)

!...a61, a62, a63, a64, a65
cfmat4(1, 1) = vcoef(4); cfmat4(1, 2) = vcoef(6); cfmat4(1, 3) = vcoef(1);  cfmat4(1, 4) = vcoef(3);
cfmat4(2, 1) = vcoef(6); cfmat4(2, 2) = vcoef(7); cfmat4(2, 3) = vcoef(3);  cfmat4(2, 4) = vcoef(2);

cfmat4(3, 1) = vcoef(8) - vcoef(1)**2;
cfmat4(3, 2) = vcoef(10)- vcoef(1)*vcoef(3);
cfmat4(3, 3) = vcoef(4);
cfmat4(3, 4) = vcoef(6);

cfmat4(4, 1) = vcoef(10) - vcoef(1)*vcoef(3);
cfmat4(4, 2) = vcoef(11) - vcoef(3)**2;
cfmat4(4, 3) = vcoef(6);
cfmat4(4, 4) = vcoef(7);

cfb4(1) = -vcoef(7); cfb4(2) =-vcoef(5);
cfb4(3) = -vcoef(11) + vcoef(1)*vcoef(2)
cfb4(4) = -vcoef(12) + vcoef(2)*vcoef(3)
!
cfmti4 = 0.d0
binv4  = 0.d0

!if(ielem.eq.16)then
!do idegr=1,4
! print*,'iem',cfmat4(idegr,:)
!enddo
!endif

cfmti4 = cfmat4
!call DGETRF(4,4,cfmti4,4,ipiv(1:4),info)
!call DGETRI(4,cfmti4,4,ipiv(1:4),work(1:4),4,info)

call Matinv225(cfmat4,cfmti4,4)

!call getinvmat(4, cfmat4(1:4,1:4), cfmti4(1:4, 1:4), binv4(1:4))

!
do im =1, 4
geoph(12, ielem) = geoph(12, ielem) + cfmti4(1, im)*cfb4(im)
geoph(13, ielem) = geoph(13, ielem) + cfmti4(2, im)*cfb4(im)
geoph(14, ielem) = geoph(14, ielem) + cfmti4(3, im)*cfb4(im)
geoph(15, ielem) = geoph(15, ielem) + cfmti4(4, im)*cfb4(im)
enddo
!if(ielem.eq.16) print*,'iem2,cfmti4',cfmti4
geoph(16, ielem) = -(geoph(12, ielem)*vcoef(1) + geoph(13, ielem)*vcoef(3) + vcoef(2))
!
endif

enddo !do ie = 1, nquad

!print*,'test1',ielem

!endif !if(npoly.eq.2)then

!...I.3 Mass matrix
do ie = 1,nquad
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
dxc =sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc =dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)

!...mass
masel = geoel(4, ielem)

!...Initial zero
bb = 0.d0
bbq= 0.d0
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
bqp(2) = (xgaus-xmc)/dxc
bqp(3) = (ygaus-ymc)/dyc

!...Orthogonal basis functions
bqo(1) = 1.d0
bqo(2) = bqp(2)
bqo(3) = geoph(3, ielem)*bqp(2) + bqp(3) + geoph(4, ielem)

!
if(npoly.eq.2)then
bq(4) = 0.5d0*bq(2)*bq(2) - geoel(19, ielem)
bq(5) = 0.5d0*bq(3)*bq(3) - geoel(20, ielem)
bq(6) =       bq(2)*bq(3) - geoel(21, ielem)

bqp(4) = bqp(2)*bqp(2)
bqp(5) = bqp(3)*bqp(3)
bqp(6) = bqp(2)*bqp(3)

bqo(4) = bqp(4) + geoph(5, ielem)*bqp(2) + geoph(6, ielem)*bqp(3) + geoph(7, ielem)
bqo(5) = geoph(8, ielem)*bqp(4) + bqp(6) + geoph(9, ielem)*bqp(2) +&
geoph(10, ielem)*bqp(3) + geoph(11, ielem)
bqo(6) = geoph(12, ielem)*bqp(4) + geoph(13, ielem)*bqp(6) + bqp(5) +&
geoph(14, ielem)*bqp(2) + geoph(15, ielem)*bqp(3) + geoph(16, ielem)
endif

!...Matrix(5x5)
do idegr =2, ndegr
do jdegr =2, ndegr
bb(idegr-1, jdegr-1)  =  bb(idegr-1, jdegr-1) + rhogi*djacoi*bq(idegr)*bqo(jdegr)
bbq(idegr-1, jdegr-1) = bbq(idegr-1, jdegr-1) + rhogi*djacoi*bqo(idegr)*bqo(jdegr)
enddo
enddo
!
enddo

!...RHS
rhs = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
rhs(1:nq, idegr) = rhs(1:nq, idegr) + unkno(jdegr, 1:nq, ielem)*bb(jdegr-1, idegr-1)
enddo
enddo

!...Store the array for later use
ijmat = 0
do im = 2, ndegr
do jm = im, ndegr
ijmat = ijmat +1
mbbq(ijmat, ielem) = bbq(im-1, jm-1)
enddo
enddo

!...Check
!if(ielem.eq.16)then
!do idegr=1,ndegr-1
!print*,'mat7',bbq(idegr,:),geoph(12:16,ielem)
!enddo
!endif

!...Inverse of bbq
if(npoly.eq.1)then
detma = bbq(1, 1)*bbq(2, 2)-bbq(1, 2)*bbq(2, 1)
bbqi(1, 1) = bbq(2, 2)
bbqi(1, 2) =-bbq(1, 2)
bbqi(2, 1) =-bbq(2, 1)
bbqi(2, 2) = bbq(1, 1)

bbqi = bbqi/detma
!
detma = bb(1, 1)*bb(2, 2)-bb(1, 2)*bb(2, 1)
bbi(1, 1) = bb(2, 2)
bbi(1, 2) =-bb(1, 2)
bbi(2, 1) =-bb(2, 1)
bbi(2, 2) = bb(1, 1)

bbi = bbi/detma
elseif(npoly.eq.2)then
binv = 0.d0
bbqi = 0.d0
call matinv225(bbq, bbqi, 5)


!binv = 0.d0
!bbi = 0.d0
!call getinvmat(5, bb, bbi, binv)
endif
!

!...LHS
unknp(:, 1:nq, ielem) = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
unknp(idegr, 1:nq, ielem) = unknp(idegr, 1:nq, ielem) + rhs(1:nq, jdegr)*bbqi(jdegr-1,idegr-1)
enddo
enddo

!...Cell average on physical domain
unknp(1, 1:nq, ielem) = unkno(1, 1:nq, ielem)

!....
!if(ielem.eq.16)then
!do idegr=1,ndegr-1
! print*,'bad3',ielem,bbq(idegr,:)
!enddo
!print*,'unkno',ielem,unknp(2:ndegr, 3, ielem)
!endif

!...Store the array for later use
ijmat = 0
do im = 2, ndegr
do jm = 2, ndegr
ijmat = ijmat +1
mbbi(ijmat, ielem) = bb(im-1, jm-1)
enddo
enddo

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
mapmt(1, 1, ielem) = 1.d0
mapmt(1, 2, ielem) = 0.d0
mapmt(2, 1, ielem) = 0.d0
mapmt(2, 2, ielem) = 1.d0

enddo !do ie = 1, nquad
!
!...Part II: Recontruction
!
do ie = 1,nquad
!
ielem = ie + ntria
!
if(geoel(10, ielem).lt.10) cycle
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
!
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...physical mass center
xmc =geoph(1, ielem)
ymc =geoph(2, ielem)

!...scale
dxc = sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))
volel = geoel(3, ielem)

!...High-order orthogonal coefficients
cothg(3,1) = geoph(3, ielem)
cothg(3,2) = geoph(4, ielem)

cothg(4,1) = geoph(5, ielem)
cothg(4,2) = geoph(6, ielem)
cothg(4,3) = geoph(7, ielem)

cothg(5,1) = geoph(8, ielem)
cothg(5,2) = geoph(9, ielem)
cothg(5,3) = geoph(10, ielem)
cothg(5,4) = geoph(11, ielem)

cothg(6,1) = geoph(12, ielem)
cothg(6,2) = geoph(13, ielem)
cothg(6,3) = geoph(14, ielem)
cothg(6,4) = geoph(15, ielem)
cothg(6,5) = geoph(16, ielem)

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
xpq_j(1, 1:4) = coord(1, ipqua(1:4,jelem-ntria))
xpq_j(2, 1:4) = coord(2, ipqua(1:4,jelem-ntria))
!
dxj = sqrt(geoel(3, jelem))!maxval(xpq_j(1, 1:4)) - minval(xpq_j(1, 1:4))
dyj = dxj!maxval(xpq_j(2, 1:4)) - minval(xpq_j(2, 1:4))
!
isten = isten + 1
unksc(1, 1:nq, isten) = unknp(1, 1:nq, ielem)
unksc(2, 1:nq, isten) = unknp(2, 1:nq, jelem)*dxc/dxj
unksc(3, 1:nq, isten) = unknp(3, 1:nq, jelem)*dyc/dyj
!
if(npoly.eq.2)then
unksc(4, 1:nq, isten) = unknp(4, 1:nq, jelem)*dxc**2/dxj**2
unksc(5, 1:nq, isten) = unknp(5, 1:nq, jelem)*dyc**2/dyj**2
unksc(6, 1:nq, isten) = unknp(6, 1:nq, jelem)*dyc*dxc/dyj/dxj
endif
endif

enddo !do ies = 1, 4

Case (2)

!...Weno
isten = 1
unksc(:, 1:nq, isten) = unknp(:, 1:nq, ielem)

!...WENO
call weno_stncl_quad_ortho4(isten, nsten, ielem, ipqua, esqua, unksc, unknp, geoel,geoph, geopj,&
coord, cooro, esuv1, esuv2)

end select
!
if(ielem.eq.1921)then
!print*,'ielem1921',unksc(2,1,8)
!print*,'stl-v',ielem,isten,unksc(2:ndegr,3,1:isten)
!print*,'stl-e',ielem,isten,unksc(2:ndegr,4,2:5)
endif


!...Local system
!call getunkno_local(0, isten, mapmt(:,:,ielem), unksc(:,:,1:isten))

!...Update local stencile values
unksl = unksc

!...  b. curvatures for the face-neighboring cells
ifprj = 0
unkf_cha = 0.d0
unkf_phy = 0.d0

!...Loop over the faces
do ies = 1, 4

jelem = esqua(ies,ie)
if(jelem .le. ncell) then

!...Find the projection vector
ifprj = ifprj +1

!...normal vector only using linear face even for curved meshes
xpf(1, 1:2) = xpq(1, mapfe(1:2, ies))
xpf(2, 1:2) = xpq(2, mapfe(1:2, ies))

dtx = xpf(1, 2) - xpf(1, 1)
dty = xpf(2, 2) - xpf(2, 1)

dlgt = sqrt(dtx**2 + dty**2)

dnx = dty/dlgt
dny =-dtx/dlgt
!
!call getvector_local(mapmt(:,:,ielem), dnx,dny)

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
!
!call getvector_local(mapmt(:,:,ielem), dxc,dyc)

!...Smooth indicator for every stencle
if(npoly.le.1)then
do is= 1, isten
do iq= 1, nq
os(iq, is) = (((unksf_cha(2, iq, is)+cothg(3, 1)*unksf_cha(3, iq, is))/dxc)**2 +&
(unksf_cha(3, iq, is)/dyc)**2)!/(abs(unksf_cha(1, iq, 1))+epsil)**1
enddo
enddo
elseif(npoly.eq.2)then
do is= 1, isten
do iq= 1, nq
!...
smthid(1) = 1.d0/(dxc**2)*(unksf_cha(2, iq, is) +&
               cothg(3, 1)*unksf_cha(3, iq, is) +&
               cothg(4, 1)*unksf_cha(4, iq, is) +&
               cothg(5, 2)*unksf_cha(5, iq, is) +&
               cothg(6, 3)*unksf_cha(6, iq, is))**2 +&
1.d0/(dxc**2)*geopj(1, ielem)*(2.d0*unksf_cha(4, iq, is)+&
                               2.d0*cothg(5, 1)*unksf_cha(5, iq, is)+&
                               2.d0*cothg(6, 1)*unksf_cha(6, iq, is))**2 +&
1.d0/(dxc**2)*geopj(2, ielem)*(unksf_cha(5, iq, is)+&
                               cothg(6, 2)*unksf_cha(6, iq, is))**2 +&
2.d0/(dxc**2)*geopj(3, ielem)*(unksf_cha(5, iq, is)+&
                               cothg(6, 2)*unksf_cha(6, iq,is))*(2.d0*unksf_cha(4, iq, is)+&
                               2.d0*cothg(5, 1)*unksf_cha(5, iq, is)+&
                               2.d0*cothg(6, 1)*unksf_cha(6, iq, is))
!
smthid(2) = 1.d0/(dyc**2)*(unksf_cha(3, iq, is) +&
                           cothg(4, 2)*unksf_cha(4, iq, is) +&
                           cothg(5, 3)*unksf_cha(5, iq, is) +&
                           cothg(6, 4)*unksf_cha(6, iq, is))**2+&
            1.d0/(dxc**2)*geopj(1, ielem)*(unksf_cha(5, iq, is)+&
                                           cothg(6, 2)*unksf_cha(6, iq, is))**2 +&
            1.d0/(dxc**2)*geopj(2, ielem)*(2.d0*unksf_cha(6, iq, is))**2 +&
            2.d0/(dxc**2)*geopj(3, ielem)*(unksf_cha(5, iq, is)+&
                            cothg(6, 2)*unksf_cha(6, iq,is))*(2.d0*unksf_cha(6, iq, is))

!
smthid(3) = 1.d0/(dxc**4)*(2.d0*unksf_cha(4, iq, is) + 2.d0*cothg(5, 1)*unksf_cha(5, iq, is)+&
                           2.d0*cothg(6, 1)*unksf_cha(6, iq, is) )**2
smthid(4) = 1.d0/(dyc**4)*(2.d0*unksf_cha(6, iq, is))**2
smthid(5) = 1.d0/(dyc*dxc)**2*(unksf_cha(5, iq, is) + cothg(6, 2)*unksf_cha(6, iq, is))**2
!
os(iq, is) = (smthid(1) + smthid(2)) +&
(smthid(3) + smthid(4) + smthid(5))*volel

!os(iq, is) =min((smthid(1) + smthid(2)),(smthid(3) + smthid(4) + smthid(5))*volel)
!if(ielem.eq.7.and.is.eq.1) print*,'derivative',iq,smthid(1:5)
enddo
enddo

endif

!...Mometum
!os1(:) = 0.25d0*(os(1, :)+os(2, :)+os(3, :)+os(4, :))
os1(:) = max(os(1, :),os(2, :),os(3, :),os(4, :))


do is= 1, isten
!os(1:nq, is) = os1(is)
enddo

!...Linear weight for central and biased stencle
weigl(1)=.5d0;
weigl(2:isten)= 0.5d0/(isten-1.d0)


weigt = 0.d0
do is= 1, isten
do iq =1,nq
weigh(iq, is) = weigl(is)/(epsil+os(iq,is))**rpowe
enddo
enddo

!if(ielem.eq.77)then
!print*,'ielem12',ielem,qinvm(1,:)
!print*,'ielem13',ifprj,os(1:4,5),epsil+os(1:4,5)
!endif
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
unk_phy(2:ndegr, iq) =unk_phy(2:ndegr, iq) + weige(iq, ifa)*unkf_phy(2:ndegr,iq,ifa)
enddo
enddo

!...Local system
!call getunkno_local(1, 1, mapmt(:,:,ielem), unk_phy)

!...Transfer back the reference configuration
ijmat = 0
do im = 2, ndegr
do jm = im, ndegr
ijmat = ijmat +1
bbq(im-1, jm-1) = mbbq(ijmat, ielem)
bbq(jm-1, im-1) = mbbq(ijmat, ielem)
enddo
enddo

ijmat = 0
do im = 2, ndegr
do jm = 2, ndegr
ijmat = ijmat +1
bb(im-1, jm-1) = mbbi(ijmat, ielem)
enddo
enddo

!...RHS
rhs = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
rhs(1:nq, idegr) = rhs(1:nq, idegr) + unk_phy(jdegr, 1:nq)*bbq(jdegr-1, idegr-1)
enddo
enddo
!
rhs1(1, 1:ndegr-1) = rhs(1,2:ndegr)
rhs1(2, 1:ndegr-1) = rhs(2,2:ndegr)
rhs1(3, 1:ndegr-1) = rhs(3,2:ndegr)
rhs1(4, 1:ndegr-1) = rhs(4,2:ndegr)
!
!
!if(ielem.eq.16)then
!do idegr=1,ndegr-1
!print*,'ielem7',bbq(idegr, :),unk_phy(idegr+1,3)
!enddo
!endif
!
binv = 0.d0
bbsymi = bb
!call getinvmat(ndegr-1, bb, bbsymi, binv)
!
!call DGETRF(ndegr-1,ndegr-1,bbsymi,ndegr-1,ipiv,info)
!call DGETRI(ndegr-1,bbsymi,ndegr-1,ipiv,work,ndegr-1,info)

call matinv225(bb, bbsymi, ndegr-1)

!
!if(info.ne.0) stop'Matrix inverse fail'
!
!if(ielem.eq.16)then
!do idegr=1,ndegr-1
!print*,'ielem8',bb(idegr,:)
!enddo
!endif

!call RANDOM_NUMBER(binv)
!
!if(ielem.eq.16)then
!do idegr=1,ndegr-1
! print*,'ielem8',bbsymi(idegr,:),rhs1(3, idegr)
!enddo
!endif
!...LHS
unkno(2:ndegr, 1:nq, ielem) = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
unkno(idegr, 1:nq, ielem) = unkno(idegr, 1:nq, ielem) + rhs1(1:nq, jdegr-1)*bbsymi(jdegr-1, idegr-1)
enddo
enddo
!
!
!if(npoly.eq.2)unkno(4:ndegr, 1:nq, ielem) = 0.d0
!
!if(ielem.eq.7)then
!print*,'ielem1',ielem,unkf_phy(2:ndegr, 3,1:4)
!print*,'ielem2',ielem,unkf_phy(2:ndegr,3, 1:4)
!print*,'ielem3',ielem,unk_phy(2:ndegr, 3)
!print*,'ielem4',ielem,unk_phy(2:ndegr, 4)
!print*,'ielem2',ielem,unkno(2:ndegr, 3, ielem)
!endif
!...  end of the loop over the quad cells
enddo
!...Release memory
deallocate(bbqi, binv, mbbq, mbbi,cfmat,cfmti,cfb)
return
end subroutine weno_char_quad_orthg2
!
!
!...WENO stnecils
!
subroutine weno_stncl_quad_ortho4(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoel, geoph, geopj,&
coord, cooro, esuv1, esuv2)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:ngeel,1:nsize),        intent(in) ::geoel
real*8,dimension(1:16,1:nsize),        intent(in) ::geoph
real*8,dimension(1:3, 1:nsize),        intent(in) ::geopj
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)

integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(5)::jelaj
integer,dimension(9)::ielaj
integer,  dimension(1:nvqua):: ipq
integer::ielst1(4,7),ielst2(4,4),ielst3(4,4)
integer::idxiv(3)

real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:2, 1:nvqua)::xpq,xpqj,xpqi,xpqji
real*8,dimension(1:ndegr)::bqp,bqo,bqoj,bqpj
real*8::work(ndegr-1)
integer::ipiv(ndegr-1)
integer::info
!
real*8, allocatable::matst(:, :), mtsti(:, :), binv(:), rhs(:,:), vcrhs(:,:), mtsym(:,:),rhssym(:,:)

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dxc, dyc, xmc, ymc, xcj, ycj
real*8:: dxcrt, dycrt,dxj,dyj, dyrt, dxrt
real*8:: dai, daj, bqxj, bqyj, dart
real*8:: detma
real*8:: dxdri,dxdsi,dydri,dydsi,wi,r,s,xcrho,ycrho,rhogi
real*8:: rcj,scj,djacoi,xgausi,ygausi,xgaus,ygaus
!
integer:: ie,ies, jelem,  ishp, ieaj, istor,iv,nvaj, ivaj, id, ifaj
integer:: idegr, is, iq,igaus, jdegr
integer:: idxod

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Allocate memory
allocate(matst(ndegr-1, ndegr-1), mtsti(ndegr-1, ndegr-1), binv(ndegr-1), rhs(nq, ndegr-1) ,vcrhs(ndegr-1, ndegr-1))
allocate(mtsym(ndegr-1, ndegr-1), rhssym(nq, ndegr-1))

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Scaling parameter
dxc = sqrt(geoel(3, ielem))!maxval(xpq(1, 1:4)) - minval(xpq(1, 1:4))
dyc = dxc!maxval(xpq(2, 1:4)) - minval(xpq(2, 1:4))

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Adjacent cells
ielaj = 0
ielaj(1) = ielem
ielaj(2:5) = esqua(1:4, ielem)

!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
do iv = 1, 4

jelaj = 0
nvaj = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielem)then
nvaj = nvaj +1
jelaj(nvaj) = jelem
endif
enddo !istor
!...Stencil No.
idxod = 1
idxiv = 0
do ivaj=1,3
do ifaj=1,4
if(jelaj(ivaj).eq.ielaj(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
!if(ielem.eq.1)print*,'ielem',iv,nvaj,ivaj,ifaj,jelaj(ivaj),ielaj(ifaj+1),idxod
endif
enddo
enddo

!...Order the adjacent cells
do ivaj=1,3
if(idxiv(ivaj) == 0)then
if(idxod.eq.2)then
ielaj(6) = jelaj(ivaj)
elseif(idxod.eq.6)then
ielaj(7) = jelaj(ivaj)
elseif(idxod.eq.12)then
ielaj(8) = jelaj(ivaj)
elseif(idxod.eq.4)then
ielaj(9) = jelaj(ivaj)
else
!print*,'Wrong stencil!',ielem,idxod
!stop
endif
endif
enddo
!
enddo !do iv = 1, 4

!if(ielem.eq.1153)print*,'ielem',ielem,iv,nvaj,ielaj(1:9)

!
if(minval(ielaj(1:9)).eq.0)then
isten = 1

else
!...Part 2: Eliminate cells without 9 adjacient cells
ielst1(1, 1:7) =(/1,2,5,6,9,8, 3/)
ielst1(2, 1:7) =(/1,3,2,7,6,9, 5/)
ielst1(3, 1:7) =(/1,4,3,8,7,6,5/)
ielst1(4, 1:7) =(/1,5,4,9,8,7,3/)
!
ielst2(1, 1:4) =(/1,2,6,9/)
ielst2(2, 1:4) =(/1,3,7,6/)
ielst2(3, 1:4) =(/1,4,8,7/)
ielst2(4, 1:4) =(/1,5,9,8/)

ielst3(1, 1:4) =(/1,2,5,9/)
ielst3(2, 1:4) =(/1,3,2,6/)
ielst3(3, 1:4) =(/1,4,3,7/)
ielst3(4, 1:4) =(/1,5,4,8/)

!...Zero out other stencils
unkpe(2:ndegr,:,2:nsten) = 0.d0

!...1st stencil
do is =1 ,4
!
isten = isten +1
!
do ies =1, ndegr-1
!
jelem = ielaj(ielst1(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = sqrt(geoel(3, jelem)) !maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = dxj !maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dai  = geoel(3, ielem)
daj  = geoel(3, jelem)
dart = daj/dai

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)
!
bqxj = (xcj-xmc)/dxc
bqyj = (ycj-ymc)/dyc

matst(ies, 1)  = bqxj
matst(ies, 2)  = geoph(3, ielem)*bqxj + bqyj + geoph(4, ielem)
!
if(npoly.eq.2)then
matst(ies, 3)  = dart*geopj(1, jelem) + bqxj**2  + geoph(5, ielem)*bqxj +  geoph(6, ielem)*bqyj + geoph(7, ielem)
matst(ies, 4)  = dart*geopj(1, jelem)*geoph(8, ielem)   + bqxj**2*geoph(8, ielem) +&
dart*geopj(3, jelem) + bqxj*bqyj + geoph(9, ielem)*bqxj + geoph(10, ielem)*bqyj + geoph(11, ielem)
matst(ies, 5)  = dart*geopj(1, jelem)*geoph(12, ielem) +   bqxj**2*geoph(12, ielem) +&
geoph(13, ielem)*dart*geopj(3, jelem) + geoph(13, ielem)*bqxj*bqyj +&
dart*geopj(2, jelem) + bqyj**2 + geoph(14, ielem)*bqxj + geoph(15, ielem)*bqyj + geoph(16, ielem)
endif

enddo

mtsym = matst
call matinv225(matst, mtsti, ndegr-1)

!
!if(ielem.eq.1921.and.isten.eq.2)then
!print*,'local',geoph(1:16, ielem)
!do idegr=1,5
!print*,'mat',ielem,isten,matst(idegr, :), rhs(1,idegr)
!enddo

!do idegr=1,5
!print*,'min',ielem,isten,mtsti(idegr, :)
!enddo
!endif

!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhs(1:nq, ies)
enddo
enddo
!
enddo

!...2nd stencil
!...The 2st 3 Lagrangian polynomial
if(npoly.eq.2)then

do is =1 ,4
!...Increasing stencils
isten = isten +1

!...The 2 Hermite polynomials from the 2nd stencil
rhs   = 0.d0
matst = 0.d0

!...1st part
do ies =1, 3
!
jelem = ielaj(ielst3(is, ies+1))
!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = sqrt(geoel(3, jelem)) !maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = dxj !maxval(xpqj(2, 1:4)) - minval(xpqj(2, 1:4))
!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)
!
dai  = geoel(3, ielem)
daj  = geoel(3, jelem)
dart = daj/dai

rhs(1:nq, ies) = unknp(1, 1:nq, jelem)- unknp(1, 1:nq, ielem)
!
bqxj = (xcj-xmc)/dxc
bqyj = (ycj-ymc)/dyc

matst(ies, 1)  = bqxj
matst(ies, 2)  = geoph(3, ielem)*bqxj + bqyj + geoph(4, ielem)
!
if(npoly.eq.2)then
matst(ies, 3)  = dart*geopj(1, jelem) + bqxj**2  + geoph(5, ielem)*bqxj +  geoph(6, ielem)*bqyj + geoph(7, ielem)
matst(ies, 4)  = dart*geopj(1, jelem)*geoph(8, ielem)   + bqxj**2*geoph(8, ielem) +&
dart*geopj(3, jelem) + bqxj*bqyj + geoph(9, ielem)*bqxj + geoph(10, ielem)*bqyj + geoph(11, ielem)
matst(ies, 5)  = dart*geopj(1, jelem)*geoph(12, ielem) +   bqxj**2*geoph(12, ielem) +&
geoph(13, ielem)*dart*geopj(3, jelem) + geoph(13, ielem)*bqxj*bqyj +&
dart*geopj(2, jelem) + bqyj**2 + geoph(14, ielem)*bqxj + geoph(15, ielem)*bqyj + geoph(16, ielem)
endif

enddo

!...2nd part
jelem = ielaj(ielst3(is, 2))

!...Coordinates
xpqj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))
!
xpqji(1, 1:nvqua) = cooro(1, ipqua(1:nvqua,jelem))
xpqji(2, 1:nvqua) = cooro(2, ipqua(1:nvqua,jelem))
!...Scaling parameter
dxj = sqrt(geoel(3, jelem)) !maxval(xpqj(1, 1:4)) - minval(xpqj(1, 1:4))
dyj = dxj
!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!
xcj = geoph(1, jelem)
ycj = geoph(2, jelem)

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqji, rcj, scj, xcrho, ycrho)
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
dxdri = dxdri + dsprq(ishp)*xpqji(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqji(1,ishp)

dydri = dydri + dsprq(ishp)*xpqji(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqji(2,ishp)
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
xgausi = xgausi + shpq(ishp)*xpqji(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqji(2,ishp)

xgaus  = xgaus + shpq(ishp)*xpqj(1,ishp)
ygaus  = ygaus + shpq(ishp)*xpqj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions

bqp(1) = 1.d0
bqp(2) = (xgaus-xmc)/dxc
bqp(3) = (ygaus-ymc)/dyc

bqpj(1) = 1.d0
bqpj(2) = (xgaus-xcj)/dxj
bqpj(3) = (ygaus-ycj)/dyj


!...Orthogonal basis functions
bqo(1) = 1.d0
bqo(2) = bqp(2)
bqo(3) = geoph(3, ielem)*bqp(2) + bqp(3) + geoph(4, ielem)

bqoj(1) = 1.d0
bqoj(2) = bqpj(2)
bqoj(3) = geoph(3, jelem)*bqpj(2) + bqpj(3) + geoph(4, jelem)

!
if(npoly.eq.2)then
bqp(4) = bqp(2)*bqp(2)
bqp(5) = bqp(3)*bqp(3)
bqp(6) = bqp(2)*bqp(3)

bqo(4) = bqp(4) + geoph(5, ielem)*bqp(2) + geoph(6, ielem)*bqp(3) + geoph(7, ielem)
bqo(5) = geoph(8, ielem)*bqp(4) + bqp(6) + geoph(9, ielem)*bqp(2) +&
geoph(10, ielem)*bqp(3) + geoph(11, ielem)
bqo(6) = geoph(12, ielem)*bqp(4) + geoph(13, ielem)*bqp(6) + bqp(5) +&
geoph(14, ielem)*bqp(2) + geoph(15, ielem)*bqp(3) + geoph(16, ielem)

!bqo(4) = bqp(4) + geoph(7, ielem)
!bqo(5) = bqp(6)
!bqo(6) = bqp(5) + geoph(16, ielem)


!
bqpj(4) = bqpj(2)*bqpj(2)
bqpj(5) = bqpj(3)*bqpj(3)
bqpj(6) = bqpj(2)*bqpj(3)

bqoj(4) = bqpj(4) + geoph(5, jelem)*bqpj(2) + geoph(6, jelem)*bqpj(3) + geoph(7, jelem)
bqoj(5) = geoph(8, jelem)*bqpj(4) + bqpj(6) + geoph(9, jelem)*bqpj(2) +&
geoph(10, jelem)*bqpj(3) + geoph(11, jelem)
bqoj(6) = geoph(12, jelem)*bqpj(4) + geoph(13, jelem)*bqpj(6) + bqpj(5) +&
geoph(14, jelem)*bqpj(2) + geoph(15, jelem)*bqpj(3) + geoph(16, jelem)

!bqoj(4) = bqpj(4) + geoph(7, jelem)
!bqoj(5) = bqpj(6)
!bqoj(6) = bqpj(5) + geoph(16, jelem)

endif
!
!do idegr=2, ndegr
!rhs(1:nq, idegr-1) = rhs(1:nq, idegr-1) + rhogi*djacoi*bqoj(idegr)**2*unknp(idegr,1:nq,jelem)
!enddo

rhs(1:nq, 4) = rhs(1:nq, 4) + rhogi*djacoi*bqoj(2)**2*unknp(2,1:nq,jelem)
rhs(1:nq, 5) = rhs(1:nq, 5) + rhogi*djacoi*bqoj(3)**2*unknp(3,1:nq,jelem)


!...Matrix(5x5)
!do idegr =2, ndegr
!do jdegr =2, ndegr
!matst(idegr-1, jdegr-1)  =  matst(idegr-1, jdegr-1) + rhogi*djacoi*bqo(jdegr)*bqoj(idegr)
!enddo
!enddo

do jdegr =2, ndegr
matst(4, jdegr-1)  =  matst(4, jdegr-1) + rhogi*djacoi*bqo(jdegr)*bqoj(2)
matst(5, jdegr-1)  =  matst(5, jdegr-1) + rhogi*djacoi*bqo(jdegr)*bqoj(3)
enddo
!
enddo
!
call matinv225(matst, mtsti, 5)
!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhs(1:nq, ies)
enddo
enddo

!
!if(ielem.eq.1921.and.isten.eq.8)then
!do idegr=1,5
!print*,'mat',ielem,isten,matst(idegr, :), rhs(1,idegr)
!enddo

!do idegr=1,5
!print*,'min',ielem,isten,mtsti(idegr, :)
!enddo
!endif

enddo


endif !if(npoly.eq.2)then

endif !...sencils

!!...Release memory
deallocate(matst, mtsti, binv,mtsym,rhssym)

return
end subroutine weno_stncl_quad_ortho4
!
!...Rigt vector matrix
!
subroutine getunkno_local2(idtrf, isten, mapmt, unkps)
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
real*8::rhoxi2,rhoet2,rhoxiet,etxi2,etet2,etxiet
real*8::ulxi2,ulet2,ulxiet,vlxi2,vlet2,vlxiet
real*8::uxxi,uyxi,uxet,uyet,vxxi,vyxi,vxet,vyet
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

!...For high-order
if(npoly.eq.2)then
!
rhoxi2 = (unkps(4, 1, is)*a11 + unkps(6, 1, is)*a21)*a11 + &
(unkps(6, 1, is)*a11 + unkps(5, 1, is)*a21)*a21
rhoet2 = (unkps(4, 1, is)*a12 + unkps(6, 1, is)*a22)*a12 + &
(unkps(6, 1, is)*a12 + unkps(5, 1, is)*a22)*a22
rhoxiet= (unkps(4, 1, is)*a12 + unkps(6, 1, is)*a22)*a11 + &
(unkps(6, 1, is)*a12 + unkps(5, 1, is)*a22)*a21
!
uxxi  =   unkps(4, 2, is)*a11 + unkps(6, 2, is)*a21
uyxi  =   unkps(6, 2, is)*a11 + unkps(5, 2, is)*a21

uxet  =   unkps(4, 2, is)*a12 + unkps(6, 2, is)*a22
uyet  =   unkps(6, 2, is)*a12 + unkps(5, 2, is)*a22

vxxi  =   unkps(4, 3, is)*a11 + unkps(6, 3, is)*a21
vyxi  =   unkps(6, 3, is)*a11 + unkps(5, 3, is)*a21

vxet  =   unkps(4, 3, is)*a12 + unkps(6, 3, is)*a22
vyet  =   unkps(6, 3, is)*a12 + unkps(5, 3, is)*a22

!
ulxi2 = (ai11*uxxi + ai12*vxxi)*a11 + &
(ai11*uyxi + ai12*vyxi)*a21
ulet2 = (ai11*uxet + ai12*vxet)*a12 + &
(ai11*uyet + ai12*vyet)*a22
ulxiet =(ai11*uxet + ai12*vxet)*a11+ &
(ai11*uyet + ai12*vyet)*a21
!
vlxi2 = (ai21*uxxi + ai22*vxxi)*a11 + &
(ai21*uyxi + ai22*vyxi)*a21
vlet2 = (ai21*uxet + ai22*vxet)*a12 + &
(ai21*uyet + ai22*vyet)*a22
vlxiet =(ai21*uxet + ai22*vxet)*a11+ &
(ai21*uyet + ai22*vyet)*a21
!
etxi2 = (unkps(4, 4, is)*a11 + unkps(6, 4, is)*a21)*a11 + &
(unkps(6, 4, is)*a11 + unkps(5, 4, is)*a21)*a21
etet2 = (unkps(4, 4, is)*a12 + unkps(6, 4, is)*a22)*a12 + &
(unkps(6, 4, is)*a12 + unkps(5, 4, is)*a22)*a22
etxiet= (unkps(4, 4, is)*a12 + unkps(6, 4, is)*a22)*a11 + &
(unkps(6, 4, is)*a12 + unkps(5, 4, is)*a22)*a21

!
unkps(4, 1, is) =rhoxi2
unkps(5, 1, is) =rhoet2
unkps(6, 1, is) =rhoxiet

unkps(4, 2, is) =ulxi2
unkps(5, 2, is) =ulet2
unkps(6, 2, is) =ulxiet

unkps(4, 3, is) =vlxi2
unkps(5, 3, is) =vlet2
unkps(6, 3, is) =vlxiet

unkps(4, 4, is) =etxi2
unkps(5, 4, is) =etet2
unkps(6, 4, is) =etxiet


endif

enddo

end subroutine getunkno_local2
!
!...Non least squares Hermite WENO P2 stnecils, in contrast with subroutine weno_stncl_quad
!...Here, in a stencile, numerical intergal is used to recover the adjacent cell's average and 1st-order derivative
!
subroutine weno_stncl_quadv2(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoel,geoph,&
coord, cooro, esuv1, esuv2, mapmt, ctsft)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:ngeel,1:nsize),        intent(in) ::geoel
real*8,dimension(1:9,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8,dimension(1:ndimn,1:ndimn),    intent(in) :: mapmt
real*8,dimension(1:3),               intent(out) :: ctsft


integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(5)::jelaj
integer,dimension(9)::ielaj
integer,  dimension(1:nvqua):: ipq
integer::ielst1(4,6),ielst2(4,4)
integer::idxiv(3)
!
real*8:: dxdri,dxdsi,dydri,dydsi
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:2, 1:nvqua)::xpqij,xpqlj
real*8, dimension(1:8,1:9)::geopl
real*8, dimension(1:ndegr,1:nq,1:9)::unkpl
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:ndegr)::bqpl, bqplj
real*8, allocatable::matst(:, :), mtsti(:, :),rhs(:,:), vcrhs(:,:)

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dxcl, dycl, xmc, ymc, xclj, yclj
real*8:: xmcl, ymcl
real*8:: dxclj,dyclj, djacoi
real*8:: xcrho, ycrho
real*8:: xgausl, ygausl, masel, xgausi, ygausi, rhogi
real*8:: wi, r, s, rcj, scj
!
integer:: ie,ies, jelem,  ishp, ieaj, istor,iv,nvaj, ivaj, id, ifaj
integer:: idegr, is, iq, igaus, je
integer:: idxod

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Allocate memory
allocate(matst(ndegr-1, ndegr-1), mtsti(ndegr-1, ndegr-1),rhs(nq, ndegr-1) ,vcrhs(ndegr-1, ndegr-1))

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Scaling parameter for target cell doesn't change without rotation
dxcl = sqrt(geoel(3, ielem))
dycl = dxcl

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Adjacent cells
ielaj = 0
ielaj(1) = ielem
ielaj(2:5) = esqua(1:4, ielem)

!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
do iv = 1, 4

jelaj = 0
nvaj = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielem)then
nvaj = nvaj +1
jelaj(nvaj) = jelem
endif
enddo !istor
!...Stencil No.
idxod = 1
idxiv = 0
do ivaj=1,3
do ifaj=1,4
if(jelaj(ivaj).eq.ielaj(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
!if(ielem.eq.1)print*,'ielem',iv,nvaj,ivaj,ifaj,jelaj(ivaj),ielaj(ifaj+1),idxod
endif
enddo
enddo

!...Order the adjacent cells
do ivaj=1,3
if(idxiv(ivaj) == 0)then
if(idxod.eq.2)then
ielaj(6) = jelaj(ivaj)
elseif(idxod.eq.6)then
ielaj(7) = jelaj(ivaj)
elseif(idxod.eq.12)then
ielaj(8) = jelaj(ivaj)
elseif(idxod.eq.4)then
ielaj(9) = jelaj(ivaj)
else
!print*,'Wrong stencil!',ielem,idxod
!stop
endif
endif
enddo
!
enddo !do iv = 1, 4
!
if(minval(ielaj(1:9)).eq.0)then
isten = 1

else
!...Part 2: Eliminate cells without 9 adjacient cells
ielst1(1, 1:6) =(/1,2,5,6,9,8/)
ielst1(2, 1:6) =(/1,3,2,7,6,9/)
ielst1(3, 1:6) =(/1,4,3,8,7,6/)
ielst1(4, 1:6) =(/1,5,4,9,8,7/)
!
ielst2(1, 1:4) =(/1,2,6,9/)
ielst2(2, 1:4) =(/1,3,7,6/)
ielst2(3, 1:4) =(/1,4,8,7/)
ielst2(4, 1:4) =(/1,5,9,8/)

!
ielst2(1, 1:4) =(/1,9, 2, 5/)
ielst2(2, 1:4) =(/1,6, 2, 3/)
ielst2(3, 1:4) =(/1,7, 4, 3/)
ielst2(4, 1:4) =(/1,8, 4, 5/)

!...Zero out other stencils
unkpe(2:ndegr,:,2:nsten) = 0.d0

!...Localize unknowns
unkpl(:,:,1:9)  = unknp(:, :, ielaj(1:9))
call getunkno_local2(0, 9, mapmt, unkpl)

!...Localize geometry
do ies = 1,9
jelem = ielaj(ies)
!...Local quad No.
je = jelem - ntria
!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,je)

xpqij(1, 1:nvqua) = cooro(1, ipq(1:nvqua))
xpqij(2, 1:nvqua) = cooro(2, ipq(1:nvqua))
!
xpqlj(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpqlj(2, 1:nvqua) = coord(2, ipq(1:nvqua))

xclj = geoph(1, jelem)
yclj = geoph(2, jelem)

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

!...mass
masel = geoel(4, jelem)

!...Local the cell center
call getvector_local(mapmt, xclj,yclj)
!...
geopl(1, ies) = xclj
geopl(2, ies) = yclj

!...Scale doesn't change with rotation
dxclj = sqrt(geoel(3, jelem))
dyclj = dxclj
!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...Zero out
geopl(3:8, ies) = 0.d0
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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqpl(1) = 1.d0
bqpl(2) = (xgausl-xclj)/dxclj
bqpl(3) = (ygausl-yclj)/dyclj

!...Matrix
geopl(3, ies) = geopl(3, ies) + 0.5d0*rhogi*djacoi*bqpl(2)*bqpl(2)
geopl(4, ies) = geopl(4, ies) + 0.5d0*rhogi*djacoi*bqpl(3)*bqpl(3)
geopl(5, ies) = geopl(5, ies) +       rhogi*djacoi*bqpl(2)*bqpl(3)
!
enddo
!
geopl(3:5, ies) =  geopl(3:5, ies)/masel
enddo

!...Record ctfst
ctsft(1:3) = geopl(3:5, 1)

!...1st stencil
do is =1 ,4
!
isten = isten +1

!...Initialize
matst = 0.d0
rhs = 0.d0

do ies =1, ndegr-1
!
jelem = ielaj(ielst1(is, ies+1))

!...Coordinates
xpqij(1, 1:nvqua) = cooro(1, ipqua(1:nvqua,jelem))
xpqij(2, 1:nvqua) = cooro(2, ipqua(1:nvqua,jelem))
!
xpqlj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqlj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))

!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...Local cell center for Target cell
xmcl = geopl(1, 1)
ymcl = geopl(2, 1)

!...RHS
rhs(1:nq, ies) = unkpl(1, 1:nq, ielst1(is, ies+1))- unkpl(1, 1:nq, 1)

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqpl(1) = 1.d0
bqpl(2) = (xgausl-xmcl)/dxcl
bqpl(3) = (ygausl-ymcl)/dycl

!...Matrix
matst(ies, 1)  = matst(ies, 1) + rhogi*djacoi*bqpl(2)/geoel(4, jelem)
matst(ies, 2)  = matst(ies, 2) + rhogi*djacoi*bqpl(3)/geoel(4, jelem)
!
if(npoly.eq.2)then

bqpl(4) = 0.5d0*bqpl(2)**2 - geopl(3, 1)
bqpl(5) = 0.5d0*bqpl(3)**2 - geopl(4, 1)
bqpl(6) =  bqpl(2)*bqpl(3) - geopl(5, 1)

matst(ies, 3)  = matst(ies, 3)  +  rhogi*djacoi*bqpl(4)/geoel(4, jelem)
matst(ies, 4)  = matst(ies, 4)  +  rhogi*djacoi*bqpl(5)/geoel(4, jelem)
matst(ies, 5)  = matst(ies, 5)  +  rhogi*djacoi*bqpl(6)/geoel(4, jelem)
endif
!
enddo

enddo

!
call matinv225(matst,mtsti,ndegr-1)

!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhs(1:nq, ies)
enddo
enddo

enddo

!
!if(ielem.eq.9.or.ielem.eq.16)then
!print*,'ielem',ielem,unkpe(3,1,1:isten)
!endif

!...2nd stencil
!...The 2st 3 Lagrangian polynomial
if(npoly.eq.2)then

do is =1 ,4
!...Increasing stencils
isten = isten +1
!...Initial zero 
matst = 0.d0
rhs =0.d0
do ies=1, 3
!
jelem = ielaj(ielst2(is, ies+1))
!...Coordinates
xpqij(1, 1:nvqua) = cooro(1, ipqua(1:nvqua,jelem))
xpqij(2, 1:nvqua) = cooro(2, ipqua(1:nvqua,jelem))
!
xpqlj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqlj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))

!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...Scaling parameter
dxclj = sqrt(geoel(3, jelem))
dyclj = dxclj
!
xclj = geopl(1, ielst2(is, ies+1))
yclj = geopl(2, ielst2(is, ies+1))
!
xmcl = geopl(1, ielst2(is, 1))
ymcl = geopl(2, ielst2(is, 1))

!...RHS
rhs(1:nq, ies) = unkpl(1, 1:nq, ielst2(is, ies+1))- unkpl(1, 1:nq, 1)

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqpl(1) = 1.d0
bqpl(2) = (xgausl-xmcl)/dxcl
bqpl(3) = (ygausl-ymcl)/dycl

!...Matrix
matst(ies, 1)  = matst(ies, 1) + rhogi*djacoi*bqpl(2)
matst(ies, 2)  = matst(ies, 2) + rhogi*djacoi*bqpl(3)
!
if(npoly.eq.2)then

bqpl(4) = 0.5d0*bqpl(2)**2 - geopl(3, 1)
bqpl(5) = 0.5d0*bqpl(3)**2 - geopl(4, 1)
bqpl(6) =  bqpl(2)*bqpl(3) - geopl(5, 1)

matst(ies, 3)  = matst(ies, 3)  +  rhogi*djacoi*bqpl(4)
matst(ies, 4)  = matst(ies, 4)  +  rhogi*djacoi*bqpl(5)
matst(ies, 5)  = matst(ies, 5)  +  rhogi*djacoi*bqpl(6)
endif

!...Hermite polynomial
if(ies.eq.1)then
bqplj(1) = 1.d0
bqplj(2) = (xgausl-xclj)/dxclj
bqplj(3) = (ygausl-yclj)/dyclj

!...Matrix
matst(4, 1)  = matst(4, 1) + rhogi*djacoi*bqpl(2)*bqplj(2)
matst(4, 2)  = matst(4, 2) + rhogi*djacoi*bqpl(3)*bqplj(2)

matst(5, 1)  = matst(5, 1) + rhogi*djacoi*bqpl(2)*bqplj(3)
matst(5, 2)  = matst(5, 2) + rhogi*djacoi*bqpl(3)*bqplj(3)
!
if(npoly.eq.2)then

bqpl(4) = 0.5d0*bqpl(2)**2 - geopl(3, 1)
bqpl(5) = 0.5d0*bqpl(3)**2 - geopl(4, 1)
bqpl(6) =  bqpl(2)*bqpl(3) - geopl(5, 1)

bqplj(4) = 0.5d0*bqplj(2)**2 - geopl(3, ielst2(is, 2))
bqplj(5) = 0.5d0*bqplj(3)**2 - geopl(4, ielst2(is, 2))
bqplj(6) =  bqplj(2)*bqplj(3)- geopl(5, ielst2(is, 2))

matst(4, 3)  = matst(4, 3)  +  rhogi*djacoi*bqpl(4)*bqplj(2)
matst(4, 4)  = matst(4, 4)  +  rhogi*djacoi*bqpl(5)*bqplj(2)
matst(4, 5)  = matst(4, 5)  +  rhogi*djacoi*bqpl(6)*bqplj(2)

matst(5, 3)  = matst(5, 3)  +  rhogi*djacoi*bqpl(4)*bqplj(3)
matst(5, 4)  = matst(5, 4)  +  rhogi*djacoi*bqpl(5)*bqplj(3)
matst(5, 5)  = matst(5, 5)  +  rhogi*djacoi*bqpl(6)*bqplj(3)
endif
!
do idegr = 2, ndegr
rhs(1:nq, 4) =  rhs(1:nq, 4) + rhogi*djacoi*bqplj(idegr)*bqplj(2)*unkpl(idegr,1:nq,ielst2(is, 2))
rhs(1:nq, 5) =  rhs(1:nq, 5) + rhogi*djacoi*bqplj(idegr)*bqplj(3)*unkpl(idegr,1:nq,ielst2(is, 2))
enddo

endif
!
enddo

enddo

!...Mass average
matst(1, :) = matst(1, :)/geoel(4, ielaj(ielst2(is, 2)))
matst(2, :) = matst(2, :)/geoel(4, ielaj(ielst2(is, 3)))
matst(3, :) = matst(3, :)/geoel(4, ielaj(ielst2(is, 4)))

matst(4, :) = matst(4, :)/geoel(4, ielaj(ielst2(is, 2)))
matst(5, :) = matst(5, :)/geoel(4, ielaj(ielst2(is, 2)))

rhs(:, 4:5) = rhs(:, 4:5)/geoel(4, ielaj(ielst2(is, 2)))

!...Invert matrix
call matinv225(matst, mtsti, ndegr-1)
!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) =  unkpe(idegr, 1:nq, isten) + mtsti(idegr-1, ies)*rhs(1:nq, ies)
enddo
enddo

!if((ielem.eq.1153))then
!do idegr=1,ndegr-1
 !print*,'mat',matst(idegr,:),rhs(1,idegr)
!enddo
!print*,'iele2',ielem,ielaj(ielst2(is, 1:4))
!endif

enddo


endif !if(npoly.eq.2)then

endif !...sencils
!!...Release memory
deallocate(matst, mtsti, rhs, vcrhs)

return
end subroutine weno_stncl_quadv2
!
!...Local geometry infor
!
subroutine getgeophy_local (ipqua,geoel, geoph, ielaj, mapmt, geopl, geopj, coord, cooro)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),    intent(in) ::geoel
real*8,dimension(1:2,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
real*8,dimension(1:ndimn,1:ndimn),    intent(in) :: mapmt
real*8,dimension(16, 9), intent(out)::geopl
real*8,dimension(3, 9),  intent(out)::geopj
integer,dimension(9), intent(in)::ielaj
!
integer, dimension(1:nvqua)::ipq
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:2, 1:nvqua)::xpqij,xpqlj
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:ndegr)::bqpl

!...Local arrays
real*8::vcoef(12)
real*8,dimension(2,2)::cfmat,cfmti
real*8,dimension(3,3)::cfmat3,cfmti3
real*8,dimension(4,4)::cfmat4,cfmti4
real*8::cfb(2),cfb3(3), cfb4(4)

real*8:: c00, c05, c10, c20, epsil
real*8:: xclj, yclj
real*8:: dxclj,dyclj
real*8:: xcrho, ycrho
real*8:: xgausl, ygausl, masel, xgausi, ygausi
real*8:: djacoi,rhogi
real*8:: wi, r, s, rcj, scj
real*8:: dxdri,dxdsi,dydri,dydsi
!
integer:: ies, jelem, je, ishp
integer:: igaus,im

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...Initial zero
geopl(:,1:9) = 0.d0
!
do ies = 1, 9
!
jelem = ielaj(ies)

!...Pure quad No.
je = jelem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,je)

!...Initila coordinates
xpqij(1, 1:nvqua) = cooro(1, ipq(1:nvqua))
xpqij(2, 1:nvqua) = cooro(2, ipq(1:nvqua))
!
xpqlj(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpqlj(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...global physical mass center
xclj =geoph(1, jelem)
yclj =geoph(2, jelem)

!...Local the mass center
call getvector_local(mapmt, xclj, yclj)

!...Store
geopl(1, ies) = xclj
geopl(2, ies) = yclj

!...scale
dxclj = sqrt(geoel(3, jelem))
dyclj = dxclj

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

!...mass
masel = geoel(4, jelem)

!...Zero out vcoef
vcoef=0.d0
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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqpl(1) = 1.d0
bqpl(2) = (xgausl-xclj)/dxclj
bqpl(3) = (ygausl-yclj)/dyclj

!...Matrix
vcoef(1)  = vcoef( 1) + rhogi*djacoi*bqpl(2)*bqpl(2)
vcoef(2)  = vcoef( 2) + rhogi*djacoi*bqpl(3)*bqpl(3)
vcoef(3)  = vcoef( 3) + rhogi*djacoi*bqpl(2)*bqpl(3)
vcoef(4)  = vcoef( 4) + rhogi*djacoi*bqpl(2)**3
vcoef(5)  = vcoef( 5) + rhogi*djacoi*bqpl(3)**3
vcoef(6)  = vcoef( 6) + rhogi*djacoi*bqpl(2)**2*bqpl(3)
vcoef(7)  = vcoef( 7) + rhogi*djacoi*bqpl(3)**2*bqpl(2)
vcoef(8)  = vcoef( 8) + rhogi*djacoi*bqpl(2)**4
vcoef(9)  = vcoef( 9) + rhogi*djacoi*bqpl(3)**4
vcoef(10) = vcoef(10) + rhogi*djacoi*bqpl(2)**3*bqpl(3)
vcoef(11) = vcoef(11) + rhogi*djacoi*(bqpl(2)*bqpl(3))**2
vcoef(12) = vcoef(12) + rhogi*djacoi*bqpl(3)**3*bqpl(2)
!
enddo

!...Shift terms for P2
vcoef = vcoef/masel

!
geopj(1:3, ies) = vcoef(1:3)

!...Solve the coefficients
!...a31 a32
geopl(3, ies) =-vcoef(3)/vcoef(1)
geopl(4, ies) = 0.d0

if(npoly.eq.2)then
!...a41, a42, a43
cfmat(1, 1) = vcoef(1); cfmat(1, 2) = vcoef(3);
cfmat(2, 1) = vcoef(3); cfmat(2, 2) = vcoef(2);

cfb(1) = -vcoef(4); cfb(2) =-vcoef(6)

!...Invert 2x2 matrix
call matinv225(cfmat,cfmti,2)

geopl(5, ies) = cfmti(1, 1)*cfb(1) + cfmti(1, 2)*cfb(2)
geopl(6, ies) = cfmti(2, 1)*cfb(1) + cfmti(2, 2)*cfb(2)
geopl(7, ies) =-vcoef(1)

!...a51, a52, a53, a54
cfmat3(1, 1) = vcoef(4); cfmat3(1, 2) = vcoef(1); cfmat3(1, 3) = vcoef(3);
cfmat3(2, 1) = vcoef(6); cfmat3(2, 2) = vcoef(3); cfmat3(2, 3) = vcoef(2);

cfmat3(3, 1) = vcoef(8) - vcoef(1)**2;
cfmat3(3, 2) = vcoef(4);
cfmat3(3, 3) = vcoef(6);

cfb3(1) = -vcoef(6); cfb3(2) =-vcoef(7); cfb3(3) =-vcoef(10) + vcoef(3)*vcoef(1)

!...Invert 3x3 matrix
call matinv225(cfmat3,cfmti3,3)
!
do im =1, 3
geopl(8, ies) = geopl(8, ies) + cfmti3(1, im)*cfb3(im)
geopl(9, ies) = geopl(9, ies) + cfmti3(2, im)*cfb3(im)
geopl(10, ies)= geopl(10, ies)+ cfmti3(3, im)*cfb3(im)
enddo

geopl(11, ies)=-geopl(8, ies)*vcoef(1) - vcoef(3)

!...a61, a62, a63, a64, a65
cfmat4(1, 1) = vcoef(4); cfmat4(1, 2) = vcoef(6); cfmat4(1, 3) = vcoef(1);  cfmat4(1, 4) = vcoef(3);
cfmat4(2, 1) = vcoef(6); cfmat4(2, 2) = vcoef(7); cfmat4(2, 3) = vcoef(3);  cfmat4(2, 4) = vcoef(2);

cfmat4(3, 1) = vcoef(8) - vcoef(1)**2;
cfmat4(3, 2) = vcoef(10)- vcoef(1)*vcoef(3);
cfmat4(3, 3) = vcoef(4);
cfmat4(3, 4) = vcoef(6);

cfmat4(4, 1) = vcoef(10) - vcoef(1)*vcoef(3);
cfmat4(4, 2) = vcoef(11) - vcoef(3)**2;
cfmat4(4, 3) = vcoef(6);
cfmat4(4, 4) = vcoef(7);

cfb4(1) = -vcoef(7); cfb4(2) =-vcoef(5);
cfb4(3) = -vcoef(11) + vcoef(1)*vcoef(2)
cfb4(4) = -vcoef(12) + vcoef(2)*vcoef(3)
!
cfmti4 = 0.d0

!...Invert 4x4 matrix
call matinv225(cfmat4,cfmti4,4)

!
do im =1, 4
geopl(12, ies) = geopl(12, ies) + cfmti4(1, im)*cfb4(im)
geopl(13, ies) = geopl(13, ies) + cfmti4(2, im)*cfb4(im)
geopl(14, ies) = geopl(14, ies) + cfmti4(3, im)*cfb4(im)
geopl(15, ies) = geopl(15, ies) + cfmti4(4, im)*cfb4(im)
enddo

geopl(16, ies) = -(geopl(12, ies)*vcoef(1) + geopl(13, ies)*vcoef(3) + vcoef(2))
!

endif
!
enddo !do ies = 1, 9

end subroutine getgeophy_local
!
!...Get local physical unkno
!
!
!...Get local physical unkno
!
subroutine getunknp_local(ielem, nsten,ipqua,geoel, geoph, ielaj, mapmt, geopl, coord, cooro, unkpl,unkrl,&
                          mtbbq, mtbbi)
!
use constant
implicit none
!
integer,   intent(in):: ielem, nsten
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),    intent(in) ::geoel
real*8,dimension(1:2,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
real*8,dimension(1:ndimn,1:ndimn),    intent(in) :: mapmt
real*8,dimension(16, 9),              intent(in)::geopl
integer,dimension(9),                 intent(in)::ielaj
real*8,dimension(1:ndegr,1:nq,1:9),  intent(in)::unkrl
real*8, dimension(1:ndegr,1:nq,1:9), intent(inout)::unkpl
real*8, dimension(ndegr-1,ndegr-1), intent(out)::mtbbq, mtbbi
!
integer, dimension(1:nvqua)::ipq
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:2, 1:nvqua)::xpqij,xpqlj
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:ndegr)::bq,bqpl,bqol
real*8,dimension(ndegr-1, ndegr-1)::bb, bbq, bbi,bbqi
real*8,dimension(1:nq, 1:ndegr)::rhs

real*8:: c00, c05, c10, c20, epsil
real*8:: xclj, yclj
real*8:: dxdri,dxdsi,dydri,dydsi, djacoi, rhogi
real*8:: dxclj,dyclj,detma
real*8:: xcrho, ycrho
real*8:: xgausl, ygausl, masel, xgausi, ygausi
real*8:: wi, r, s, rcj, scj,dr,ds
!
integer:: ies, jelem,  ishp
integer:: idegr, jdegr, igaus, je,im

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...I.3 Mass matrix
dr = 1.d0
ds = 1.d0

do ies = 1, 9

!...Element No.
jelem = ielaj(ies)

!...Local quad No.
je = jelem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,je)

!...Initial
xpqij(1, 1:nvqua) = cooro(1, ipq(1:nvqua))
xpqij(2, 1:nvqua) = cooro(2, ipq(1:nvqua))

!...Current
xpqlj(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpqlj(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...physical mass center
xclj =geoph(1, jelem)
yclj =geoph(2, jelem)

!...Local the mass center
call getvector_local(mapmt, xclj, yclj)

!...scale
dxclj =sqrt(geoel(3, jelem))
dyclj =dxclj

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

!...mass
masel = geoel(4, jelem)

!...Initial zero
bb = 0.d0
bbq= 0.d0
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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bq(1) = 1.d0
bq(2) = (r-rcj)/dr
bq(3) = (s-scj)/ds


bqpl(1) = 1.d0
bqpl(2) = (xgausl-xclj)/dxclj
bqpl(3) = (ygausl-yclj)/dyclj

!...Orthogonal basis functions
bqol(1) = 1.d0
bqol(2) = bqpl(2)
bqol(3) = geopl(3, ies)*bqpl(2) + bqpl(3) + geopl(4, ies)

!
if(npoly.eq.2)then
bq(4) = 0.5d0*bq(2)*bq(2) - geoel(19, jelem)
bq(5) = 0.5d0*bq(3)*bq(3) - geoel(20, jelem)
bq(6) =       bq(2)*bq(3) - geoel(21, jelem)

bqpl(4) = bqpl(2)*bqpl(2)
bqpl(5) = bqpl(3)*bqpl(3)
bqpl(6) = bqpl(2)*bqpl(3)

bqol(4) = bqpl(4) + geopl(5, ies)*bqpl(2) + geopl(6, ies)*bqpl(3) + geopl(7, ies)
bqol(5) = geopl(8, ies)*bqpl(4) + bqpl(6) + geopl(9, ies)*bqpl(2) +&
          geopl(10, ies)*bqpl(3) + geopl(11, ies)
bqol(6) = geopl(12, ies)*bqpl(4) + geopl(13, ies)*bqpl(6) + bqpl(5) +&
          geopl(14, ies)*bqpl(2) + geopl(15, ies)*bqpl(3) + geopl(16, ies)
endif

!...Matrix(5x5)
do idegr =2, ndegr
do jdegr =2, ndegr
bb(idegr-1, jdegr-1)  =  bb(idegr-1, jdegr-1) + rhogi*djacoi*bq(idegr)*bqol(jdegr)
bbq(idegr-1, jdegr-1) = bbq(idegr-1, jdegr-1) + rhogi*djacoi*bqol(idegr)*bqol(jdegr)
enddo
enddo
!
enddo

!...RHS
rhs = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
rhs(1:nq, idegr) = rhs(1:nq, idegr) + unkrl(jdegr, 1:nq, ies)*bb(jdegr-1, idegr-1)
enddo
enddo


!...Invert matrix
call matinv225(bbq, bbqi, ndegr-1)
call matinv225(bb,  bbi, ndegr-1)

!...LHS
unkpl(:, 1:nq, ies) = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
unkpl(idegr, 1:nq, ies) = unkpl(idegr, 1:nq, ies) + rhs(1:nq, jdegr)*bbqi(jdegr-1,idegr-1)
enddo
enddo

!
!if((ielem.eq.9.or.ielem.eq.16).and.ies.eq.1)then
! do idegr=1,ndegr-1
!   print*,'matbb',bb(idegr,:)
! enddo

!print*,'unknp',ielem,unkpl(:,1,ies)

!do idegr=1,ndegr-1
!  print*,'matin',bbi(idegr,:)
!enddo
!endif

!...Store the mapping matrix
if(ies.eq.1)then
 mtbbq = bbq
 mtbbi = bbi
endif

!...Cell average on physical domain
unkpl(1, 1:nq, ies) = unkrl(1, 1:nq, ies)

enddo

end subroutine getunknp_local
!
!...Rigt vector matrix
!
subroutine getvector_local2(idtrf, mapmt, dnx, dny)
use constant
implicit none
!...Input arrays
integer,   intent(in) ::idtrf
real*8,dimension(1:2,1:2),intent(in)::mapmt
real*8,intent(inout) ::dnx, dny


!...local array
real*8:: a11, a12, a21, a22
real*8:: dnxl, dnyl

!...Local system
if(idtrf.eq.0)then
a11 = mapmt(1, 1)
a12 = mapmt(1, 2)
a21 = mapmt(2, 1)
a22 = mapmt(2, 2)
elseif(idtrf.eq.1)then
!
a11 = mapmt(2, 2)
a12 =-mapmt(1, 2)
a21 =-mapmt(2, 1)
a22 = mapmt(1, 1)
endif

!....Local system
dnxl = dnx*a11 + dny*a12
dnyl=  dnx*a21 + dny*a22
!
dnx = dnxl
dny = dnyl

end subroutine getvector_local2
!
!...WENO limiter on characteristic field using the local orthogonal basis...
!
subroutine weno_char_quad_orthgl(ipqua, esqua, unkno, geoel, coord, cooro, esuv1, esuv2)
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

!...  local arrays
integer:: ipq(nvqua)
!
real*8,dimension(1:2, 1:nvqua)::xpq, xpqi
real*8,dimension(1:2,1:nsize)         ::geoph
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8::weighq(ngausdq), posiq(2, ngausdq)
real*8::mapmt(2,2,ncell)
!
real*8:: c00, c05, c10, c20, c16,epsil
real*8:: dr, ds, rc, sc, r, s
real*8:: xmc, ymc
real*8:: masel, volel
real*8:: xcrho, ycrho, xgaus, ygaus, xgausi, ygausi
real*8:: wi,djacoi,dxdri,dxdsi,dydri,dydsi
real*8:: dudx,dudy,dvdx,dvdy
real*8:: rhogi
real*8::matra,matrb
real*8::lamda1,lamda2,lmat1,lmat2,mapt,mapd,matrc,matrd
!
integer:: ie, isten, ielem, igaus, ishp
integer:: nweno
!
real*8, allocatable::unknr(:,:,:)
!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c16   / 0.16666666666666666d0/
data c20   / 2.0d0    /
data epsil / 1.0d-6   /
!
allocate(unknr(1:ndegr, 1:nq, 1:nsize))
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
!
!...Part I: Get the L2 projection matrix
!

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
!...mass center...
rc= geoel(1, ielem)  
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

!...Local system
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

!...No local system
!mapmt(1, 1, ielem) = 1.d0
!mapmt(1, 2, ielem) = 0.d0
!mapmt(2, 1, ielem) = 0.d0
!mapmt(2, 2, ielem) = 1.d0
!
enddo !do ie = 1, nquad

!
!...Part II: Recontruction
!
unknr = unkno
!
do ie = 1,nquad
!
ielem = ie + ntria

!...Skip good cells
if(geoel(10, ielem).lt.10) cycle

select Case (nweno)

Case (1)
!...Other WENO
print*,"nweno=1 is not supported."
Case (2)
!...WENO
call weno_stncl_quad_orthol3(ielem, ipqua, esqua, unkno, unknr, geoel, geoph, &
coord, cooro, esuv1, esuv2, mapmt(:,:,ielem))

end select

enddo

!...Release memory
deallocate(unknr)
return
end subroutine weno_char_quad_orthgl
!
!...conventional WENO P2 stnecils 
!
subroutine weno_stncl_quad_Taylor(isten, nsten, ielem, ipqua, esqua, unkpe, unknp, geoel,geoph,&
coord, cooro, esuv1, esuv2, mapmt, ctsft)
!
use constant
implicit none
!
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize), intent(in)::unknp
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(inout)::unkpe
real*8,dimension(1:ngeel,1:nsize),        intent(in) ::geoel
real*8,dimension(1:9,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8,dimension(1:ndimn,1:ndimn),    intent(in) :: mapmt
real*8,dimension(1:3),               intent(out) :: ctsft


integer,    intent(in)::nsten, ielem
integer, intent(inout)::isten
integer,dimension(5)::jelaj
!
integer, parameter::nevaj = 21
integer,dimension(nevaj )::ielaj
integer,  dimension(1:nvqua):: ipq
integer::ielst1(12,8),ielst2(4,4)
integer::idxiv(3)
integer::idxel(nsize)
!
real*8:: dxdri,dxdsi,dydri,dydsi
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:2, 1:nvqua)::xpqij,xpqlj
real*8, dimension(1:8,1:nevaj )::geopl
real*8, dimension(1:ndegr,1:nq,1:nevaj )::unkpl
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:ndegr)::bqpl, bqplj
real*8, allocatable::matlg(:, :),rhslg(:,:)
real*8, allocatable::mathm(:, :),rhshm(:,:)
real*8, allocatable::matls(:, :), matin(:, :),rhsls(:,:)
!
real*8:: c00, c05, c10, c20, epsil
real*8:: dxcl, dycl, xmc, ymc, xclj, yclj
real*8:: xmcl, ymcl
real*8:: dxclj,dyclj, djacoi
real*8:: xcrho, ycrho
real*8:: xgausl, ygausl, masel, xgausi, ygausi, rhogi
real*8:: wi, r, s, rcj, scj, masej
!
integer::neqlg, neqhm, neqhm_l, neqhm_h,ieshm
integer:: ie,ies, jelem,  ishp, ieaj, istor,iv,nvaj, ivaj, id, ifaj, ifai
integer:: idegr, is, iq, igaus, je, jeloc
integer:: idxod

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Allocate memory
if(npoly.eq.1)then
neqlg = 3
elseif(npoly.eq.2)then
!...Lagrangian weno
neqlg  = 7

!...Hermite weno
neqhm   = 9
neqhm_l = 3
neqhm_h = neqhm - neqhm_l

endif

allocate(matlg(neqlg, ndegr-1), rhslg(nq, neqlg))
allocate(mathm(neqhm, ndegr-1), rhshm(nq, neqhm))
allocate(matls(ndegr-1, ndegr-1), matin(ndegr-1, ndegr-1),rhsls(nq, ndegr-1))

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Scaling parameter for target cell doesn't change without rotation
dxcl = sqrt(geoel(3, ielem))
dycl = dxcl

!...physcial center
xmc = geoph(1, ielem)
ymc = geoph(2, ielem)

!...Adjacent cells
ielaj = 0
ielaj(1) = ielem
ielaj(2:5) = esqua(1:4, ielem)

!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
call weno_stencils(nevaj,ielem, ipqua, esqua,esuv1, esuv2, ielaj)
!
!if(minval(ielaj(1:nevaj)).eq.0)print*,'ielem',ielem,ielaj

!
!print*,'ielem',ielem,ielaj
!
if(minval(ielaj(1:nevaj)).eq.0.or.maxval(ielaj(1:nevaj)).gt.ncell)then
isten = 1

else

!...Part 2: Eliminate cells without 9 adjacient cells
ielst1(1, 1:8) =(/1,2,5,6,9,8, 10, 13/)
ielst1(2, 1:8) =(/1,3,2,7,6,9, 10, 11/)
ielst1(3, 1:8) =(/1,4,3,8,7,6, 11, 12/)
ielst1(4, 1:8) =(/1,5,4,9,8,7, 12, 13/)
!
ielst1(5, 1:8) =(/1,2,5,9,13,10, 14,15/)
ielst1(6, 1:8) =(/1,2,3,6,11,10, 17,16/)
ielst1(7, 1:8) =(/1,4,3,7, 11,12,18, 19/)
ielst1(8, 1:8) =(/1,4,5,8, 13, 12, 21 ,20/)
!
ielst1(9, 1:8)  =(/1,5,8,9,20,15,21,14/)
ielst1(10, 1:8) =(/1,2, 6,9, 17, 14, 16, 15/)
ielst1(11, 1:8) =(/1,3, 7, 6, 19, 16, 18, 17/)
ielst1(12, 1:8) =(/1,4, 7 ,8, 18, 21, 19, 20/)
!
ielst2(1, 1:4) =(/1,2,6,9/)
ielst2(2, 1:4) =(/1,3,6,7/)
ielst2(3, 1:4) =(/1,4,8,7/)
ielst2(4, 1:4) =(/1,5,8,9/)

!
!ielst2(1, 1:4) =(/1,9, 5, 2/)
!ielst2(2, 1:4) =(/1,6, 3, 2/)
!ielst2(3, 1:4) =(/1,7, 3, 4/)
!ielst2(4, 1:4) =(/1,8, 5, 4/)

!...Zero out other stencils
unkpe(2:ndegr,:,2:nsten) = 0.d0

!...Localize unknowns
unkpl(:,:,1:nevaj)  = unknp(:, :, ielaj(1:nevaj))
call getunkno_local2(0, nevaj, mapmt, unkpl)

!...Localize geometry
do ies = 1,nevaj
jelem = ielaj(ies)
!...Local quad No.
je = jelem - ntria
!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,je)

xpqij(1, 1:nvqua) = cooro(1, ipq(1:nvqua))
xpqij(2, 1:nvqua) = cooro(2, ipq(1:nvqua))
!
xpqlj(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpqlj(2, 1:nvqua) = coord(2, ipq(1:nvqua))

xclj = geoph(1, jelem)
yclj = geoph(2, jelem)

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

!...mass
masel = geoel(4, jelem)

!...Local the cell center
call getvector_local(mapmt, xclj,yclj)
!...
geopl(1, ies) = xclj
geopl(2, ies) = yclj

!...Scale doesn't change with rotation
dxclj = sqrt(geoel(3, jelem))
dyclj = dxclj
!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...Zero out
geopl(3:8, ies) = 0.d0
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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqpl(1) = 1.d0
bqpl(2) = (xgausl-xclj)/dxclj
bqpl(3) = (ygausl-yclj)/dyclj

!...Matrix
geopl(3, ies) = geopl(3, ies) + 0.5d0*rhogi*djacoi*bqpl(2)*bqpl(2)
geopl(4, ies) = geopl(4, ies) + 0.5d0*rhogi*djacoi*bqpl(3)*bqpl(3)
geopl(5, ies) = geopl(5, ies) +       rhogi*djacoi*bqpl(2)*bqpl(3)
!
enddo
!
geopl(3:5, ies) =  geopl(3:5, ies)/masel
enddo

!...Record ctfst
ctsft(1:3) = geopl(3:5, 1)

!...1st stencil
do is =5, 12
!
isten = isten +1

!...Initialize matlg and rhslg
matlg = 0.d0
rhslg = 0.d0

do ies =1, neqlg
!
jelem = ielaj(ielst1(is, ies+1))

!...Coordinates
xpqij(1, 1:nvqua) = cooro(1, ipqua(1:nvqua,jelem))
xpqij(2, 1:nvqua) = cooro(2, ipqua(1:nvqua,jelem))
!
xpqlj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqlj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))

!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...mass
masej = geoel(4, jelem)

!...Local cell center for Target cell
xmcl = geopl(1, 1)
ymcl = geopl(2, 1)

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

!...Initialize matlg
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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqpl(1) = 1.d0
bqpl(2) = (xgausl-xmcl)/dxcl
bqpl(3) = (ygausl-ymcl)/dycl

!...Matrix
matlg(ies, 1)  = matlg(ies, 1) + rhogi*djacoi*bqpl(2)/masej
matlg(ies, 2)  = matlg(ies, 2) + rhogi*djacoi*bqpl(3)/masej
!
if(npoly.eq.2)then

bqpl(4) = 0.5d0*bqpl(2)**2 - geopl(3, 1)
bqpl(5) = 0.5d0*bqpl(3)**2 - geopl(4, 1)
bqpl(6) =  bqpl(2)*bqpl(3) - geopl(5, 1)

matlg(ies, 3)  = matlg(ies, 3)  +  rhogi*djacoi*bqpl(4)/masej
matlg(ies, 4)  = matlg(ies, 4)  +  rhogi*djacoi*bqpl(5)/masej
matlg(ies, 5)  = matlg(ies, 5)  +  rhogi*djacoi*bqpl(6)/masej
endif

enddo !...igaus

!...RHS
rhslg(1:nq, ies) = unkpl(1, 1:nq, ielst1(is, ies+1))- unkpl(1, 1:nq, 1)
enddo

!...Symmetrize square matrix
call  matrix_sym(neqlg, ndegr-1, nq, matlg, rhslg, matls, rhsls)

!
call matinv225(matls,matin,ndegr-1)

!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) = unkpe(idegr, 1:nq, isten) + matin(idegr-1, ies)*rhsls(1:nq, ies)
enddo
enddo

enddo

!...2nd stencil
!...The 2st 3 Lagrangian polynomial
if(npoly.eq.-2)then

do is =1 ,4
!...Increasing stencils
isten = isten +1

!...Initial zero
mathm = 0.d0
rhshm =0.d0

do ies=1, neqhm_l
!...Local No
jeloc = ielst2(is, ies+1)
!...Global No
jelem = ielaj(jeloc)
!...Coordinates
xpqij(1, 1:nvqua) = cooro(1, ipqua(1:nvqua,jelem))
xpqij(2, 1:nvqua) = cooro(2, ipqua(1:nvqua,jelem))
!
xpqlj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,jelem))
xpqlj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,jelem))

!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...mass
masej = geoel(4, jelem)

!...Scaling parameter
dxclj = sqrt(geoel(3, jelem))
dyclj = dxclj
!
xclj = geopl(1, jeloc)
yclj = geopl(2, jeloc)
!
xmcl = geopl(1, 1)
ymcl = geopl(2, 1)

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqpl(1) = 1.d0
bqpl(2) = (xgausl-xmcl)/dxcl
bqpl(3) = (ygausl-ymcl)/dycl

!...Matrix
mathm(ies, 1)  = mathm(ies, 1) + rhogi*djacoi*bqpl(2)/masej
mathm(ies, 2)  = mathm(ies, 2) + rhogi*djacoi*bqpl(3)/masej
!
if(npoly.eq.2)then
bqpl(4) = 0.5d0*bqpl(2)**2 - geopl(3, 1)
bqpl(5) = 0.5d0*bqpl(3)**2 - geopl(4, 1)
bqpl(6) =  bqpl(2)*bqpl(3) - geopl(5, 1)

mathm(ies, 3)  = mathm(ies, 3)  +  rhogi*djacoi*bqpl(4)/masej
mathm(ies, 4)  = mathm(ies, 4)  +  rhogi*djacoi*bqpl(5)/masej
mathm(ies, 5)  = mathm(ies, 5)  +  rhogi*djacoi*bqpl(6)/masej
endif

!...Hermite polynomial
if(ies.le.(neqhm_h/2))then
!
bqplj(1) = 1.d0
bqplj(2) = (xgausl-xclj)/dxclj
bqplj(3) = (ygausl-yclj)/dyclj

!
bqplj(4) = 0.5d0*bqplj(2)**2 - geopl(3, jeloc)
bqplj(5) = 0.5d0*bqplj(3)**2 - geopl(4, jeloc)
bqplj(6) =  bqplj(2)*bqplj(3)- geopl(5, jeloc)

ieshm = neqhm_l + 2*ies-1

!...Matrix
do idegr =2, ndegr
mathm(ieshm,   idegr-1)  = mathm(ieshm,   idegr-1) + &
rhogi*djacoi*bqpl(idegr)*bqplj(2)/masej
mathm(ieshm+1, idegr-1)  = mathm(ieshm+1, idegr-1) + &
rhogi*djacoi*bqpl(idegr)*bqplj(3)/masej
enddo
!
do idegr = 2, ndegr
rhshm(1:nq,   ieshm) =  rhshm(1:nq,   ieshm) + rhogi*djacoi*bqplj(idegr)*bqplj(2)*unkpl(idegr,1:nq,jeloc)/masej
rhshm(1:nq, ieshm+1) =  rhshm(1:nq, ieshm+1) + rhogi*djacoi*bqplj(idegr)*bqplj(3)*unkpl(idegr,1:nq,jeloc)/masej
enddo

endif
!
enddo

!...RHS for the Lagrangian WENO
rhshm(1:nq, ies) = unkpl(1, 1:nq, jeloc)- unkpl(1, 1:nq, 1)

enddo

!...Symmetrize square matrix
call  matrix_sym(neqhm, ndegr-1, nq, mathm, rhshm, matls, rhsls)

!...Invert matrix
call matinv225(matls, matin, ndegr-1)
!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unkpe(idegr, 1:nq, isten) = unkpe(idegr, 1:nq, isten) + matin(idegr-1, ies)*rhsls(1:nq, ies)
enddo
enddo

!if(ielem.eq.265.and.isten.eq.6)then
!do idegr=1,neqhm
!print*,'mat',isten,mathm(idegr,:),rhshm(1,idegr),ielaj(1:9)
!enddo
!endif

enddo


endif !if(npoly.eq.2)then

endif !...sencils
!!...Release memory
deallocate(matlg, rhslg, mathm, rhshm,matls, matin, rhsls)

return
end subroutine weno_stncl_quad_Taylor
!
!...Call the orthogonal basis
!
!
!...Call the orthogonal basis
!
subroutine basis_orthg(npoly,ndegr, nbsis, b_taylor, cfbs, b_out)
implicit none
integer, intent(in):: ndegr, npoly, nbsis
real*8, dimension(ndegr), intent(in)::b_taylor
real*8, dimension(14), intent(in)::cfbs
real*8, dimension(ndegr),intent(out)::b_out

!...Local
real*8::cfthg(6,6)
real*8, dimension(ndegr)::bph

!...High-order orthogonal coefficients
if(nbsis.eq.1)then
!
bph(1:3) = b_taylor(1:3)

b_out(1:3) = bph(1:3)

if(npoly.eq.2)then
b_out(4) = 0.5d0*bph(2)*bph(2) - cfbs(1)
b_out(5) = 0.5d0*bph(3)*bph(3) - cfbs(2)
b_out(6) =       bph(2)*bph(3) - cfbs(3)
endif

elseif(nbsis.eq.2)then
cfthg(3,1) = cfbs(1)
cfthg(3,2) = cfbs(2)

bph(1:3) = b_taylor(1:3)

b_out(1) = 1.d0
b_out(2) = bph(2)
b_out(3) = cfthg(3, 1)*bph(2) + bph(3) + cfthg(3, 2)

if(npoly.eq.2)then
cfthg(4,1) = cfbs(3)
cfthg(4,2) = cfbs(4)
cfthg(4,3) = cfbs(5)

cfthg(5,1) = cfbs(6)
cfthg(5,2) = cfbs(7)
cfthg(5,3) = cfbs(8)
cfthg(5,4) = cfbs(9)

cfthg(6,1) = cfbs(10)
cfthg(6,2) = cfbs(11)
cfthg(6,3) = cfbs(12)
cfthg(6,4) = cfbs(13)
cfthg(6,5) = cfbs(14)
!
bph(4) = bph(2)*bph(2)
bph(5) = bph(3)*bph(3)
bph(6) = bph(2)*bph(3)
!
b_out(4) =             bph(4) + cfthg(4, 1)*bph(2) + cfthg(4, 2)*bph(3) + cfthg(4, 3)
b_out(5) = cfthg(5, 1)*bph(4) +             bph(6) + cfthg(5, 2)*bph(2) +&
cfthg(5, 3)*bph(3) + cfthg(5, 4)
b_out(6) = cfthg(6, 1)*bph(4) + cfthg(6, 2)*bph(6) +             bph(5) +&
cfthg(6, 3)*bph(2) + cfthg(6, 4)*bph(3) + cfthg(6, 5)

endif

endif
end subroutine basis_orthg
!
!...Local geometry infor
!
subroutine getgeophy_local_Taylor (nbsis,ipqua,geoel, geoph, ielaj, mapmt, geopl, geopj, coord, cooro)
!
use constant
implicit none
!
integer, intent(in):: nbsis
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),    intent(in) ::geoel
real*8,dimension(1:2,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
real*8,dimension(1:ndimn,1:ndimn),    intent(in) :: mapmt
real*8,dimension(16, 9), intent(out)::geopl
real*8,dimension(3, 9),  intent(out)::geopj
integer,dimension(9), intent(in)::ielaj
!
integer, dimension(1:nvqua)::ipq
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:2, 1:nvqua)::xpqij,xpqlj
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:ndegr)::bqpl

!...Local arrays
real*8::vcoef(12)
real*8,dimension(2,2)::cfmat,cfmti
real*8,dimension(3,3)::cfmat3,cfmti3
real*8,dimension(4,4)::cfmat4,cfmti4
real*8::cfb(2),cfb3(3), cfb4(4)

real*8:: c00, c05, c10, c20, epsil
real*8:: xclj, yclj
real*8:: dxclj,dyclj
real*8:: xcrho, ycrho
real*8:: xgausl, ygausl, masel, xgausi, ygausi
real*8:: djacoi,rhogi
real*8:: wi, r, s, rcj, scj
real*8:: dxdri,dxdsi,dydri,dydsi
!
integer:: ies, jelem, je, ishp
integer:: igaus,im

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...Initial zero
geopl(:,1:9) = 0.d0
!
do ies = 1, 9
!
jelem = ielaj(ies)

!...Pure quad No.
je = jelem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,je)

!...Initila coordinates
xpqij(1, 1:nvqua) = cooro(1, ipq(1:nvqua))
xpqij(2, 1:nvqua) = cooro(2, ipq(1:nvqua))
!
xpqlj(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpqlj(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...global physical mass center
xclj =geoph(1, jelem)
yclj =geoph(2, jelem)

!...Local the mass center
call getvector_local(mapmt, xclj, yclj)

!...Store
geopl(1, ies) = xclj
geopl(2, ies) = yclj

!...scale
dxclj = sqrt(geoel(3, jelem))
dyclj = dxclj

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

!...mass
masel = geoel(4, jelem)

!...Zero out vcoef
vcoef=0.d0
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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqpl(1) = 1.d0
bqpl(2) = (xgausl-xclj)/dxclj
bqpl(3) = (ygausl-yclj)/dyclj

!...Matrix
vcoef(1)  = vcoef( 1) + rhogi*djacoi*bqpl(2)*bqpl(2)
vcoef(2)  = vcoef( 2) + rhogi*djacoi*bqpl(3)*bqpl(3)
vcoef(3)  = vcoef( 3) + rhogi*djacoi*bqpl(2)*bqpl(3)
vcoef(4)  = vcoef( 4) + rhogi*djacoi*bqpl(2)**3
vcoef(5)  = vcoef( 5) + rhogi*djacoi*bqpl(3)**3
vcoef(6)  = vcoef( 6) + rhogi*djacoi*bqpl(2)**2*bqpl(3)
vcoef(7)  = vcoef( 7) + rhogi*djacoi*bqpl(3)**2*bqpl(2)
vcoef(8)  = vcoef( 8) + rhogi*djacoi*bqpl(2)**4
vcoef(9)  = vcoef( 9) + rhogi*djacoi*bqpl(3)**4
vcoef(10) = vcoef(10) + rhogi*djacoi*bqpl(2)**3*bqpl(3)
vcoef(11) = vcoef(11) + rhogi*djacoi*(bqpl(2)*bqpl(3))**2
vcoef(12) = vcoef(12) + rhogi*djacoi*bqpl(3)**3*bqpl(2)
!
enddo

!...Shift terms for P2
vcoef = vcoef/masel
!
if(nbsis.eq.1)then
!
geopl(3, ies) = 0.5d0*vcoef(1)
geopl(4, ies) = 0.5d0*vcoef(2)
geopl(5, ies) =       vcoef(3)
!
elseif(nbsis.eq.2)then
geopj(1:3, ies) = vcoef(1:3)

!...Solve the coefficients
!...a31 a32
geopl(3, ies) =-vcoef(3)/vcoef(1)
geopl(4, ies) = 0.d0

if(npoly.eq.2)then
!...a41, a42, a43
cfmat(1, 1) = vcoef(1); cfmat(1, 2) = vcoef(3);
cfmat(2, 1) = vcoef(3); cfmat(2, 2) = vcoef(2);

cfb(1) = -vcoef(4); cfb(2) =-vcoef(6)

!...Invert 2x2 matrix
call matinv225(cfmat,cfmti,2)

geopl(5, ies) = cfmti(1, 1)*cfb(1) + cfmti(1, 2)*cfb(2)
geopl(6, ies) = cfmti(2, 1)*cfb(1) + cfmti(2, 2)*cfb(2)
geopl(7, ies) =-vcoef(1)

!...a51, a52, a53, a54
cfmat3(1, 1) = vcoef(4); cfmat3(1, 2) = vcoef(1); cfmat3(1, 3) = vcoef(3);
cfmat3(2, 1) = vcoef(6); cfmat3(2, 2) = vcoef(3); cfmat3(2, 3) = vcoef(2);

cfmat3(3, 1) = vcoef(8) - vcoef(1)**2;
cfmat3(3, 2) = vcoef(4);
cfmat3(3, 3) = vcoef(6);

cfb3(1) = -vcoef(6); cfb3(2) =-vcoef(7); cfb3(3) =-vcoef(10) + vcoef(3)*vcoef(1)

!...Invert 3x3 matrix
call matinv225(cfmat3,cfmti3,3)
!
do im =1, 3
geopl(8, ies) = geopl(8, ies) + cfmti3(1, im)*cfb3(im)
geopl(9, ies) = geopl(9, ies) + cfmti3(2, im)*cfb3(im)
geopl(10, ies)= geopl(10, ies)+ cfmti3(3, im)*cfb3(im)
enddo

geopl(11, ies)=-geopl(8, ies)*vcoef(1) - vcoef(3)

!...a61, a62, a63, a64, a65
cfmat4(1, 1) = vcoef(4); cfmat4(1, 2) = vcoef(6); cfmat4(1, 3) = vcoef(1);  cfmat4(1, 4) = vcoef(3);
cfmat4(2, 1) = vcoef(6); cfmat4(2, 2) = vcoef(7); cfmat4(2, 3) = vcoef(3);  cfmat4(2, 4) = vcoef(2);

cfmat4(3, 1) = vcoef(8) - vcoef(1)**2;
cfmat4(3, 2) = vcoef(10)- vcoef(1)*vcoef(3);
cfmat4(3, 3) = vcoef(4);
cfmat4(3, 4) = vcoef(6);

cfmat4(4, 1) = vcoef(10) - vcoef(1)*vcoef(3);
cfmat4(4, 2) = vcoef(11) - vcoef(3)**2;
cfmat4(4, 3) = vcoef(6);
cfmat4(4, 4) = vcoef(7);

cfb4(1) = -vcoef(7); cfb4(2) =-vcoef(5);
cfb4(3) = -vcoef(11) + vcoef(1)*vcoef(2)
cfb4(4) = -vcoef(12) + vcoef(2)*vcoef(3)
!
cfmti4 = 0.d0

!...Invert 4x4 matrix
call matinv225(cfmat4,cfmti4,4)

!
do im =1, 4
geopl(12, ies) = geopl(12, ies) + cfmti4(1, im)*cfb4(im)
geopl(13, ies) = geopl(13, ies) + cfmti4(2, im)*cfb4(im)
geopl(14, ies) = geopl(14, ies) + cfmti4(3, im)*cfb4(im)
geopl(15, ies) = geopl(15, ies) + cfmti4(4, im)*cfb4(im)
enddo

geopl(16, ies) = -(geopl(12, ies)*vcoef(1) + geopl(13, ies)*vcoef(3) + vcoef(2))
!

endif
!
endif
enddo !do ies = 1, 9

end subroutine getgeophy_local_Taylor
!
!...Get local physical unkno
!
subroutine getunknp_local_Taylor(nbsis, ielem, nsten,ipqua,geoel, geoph, ielaj, mapmt, geopl, coord, cooro, unkpl,unkrl,&
mtbbq, mtbbi)
!
use constant
implicit none
!
integer,   intent(in):: nbsis,ielem, nsten
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
real*8,dimension(1:ngeel,1:nsize),    intent(in) ::geoel
real*8,dimension(1:2,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
real*8,dimension(1:ndimn,1:ndimn),    intent(in) :: mapmt
real*8,dimension(16, 9),              intent(in)::geopl
integer,dimension(9),                 intent(in)::ielaj
real*8,dimension(1:ndegr,1:nq,1:9),  intent(in)::unkrl
real*8, dimension(1:ndegr,1:nq,1:9), intent(inout)::unkpl
real*8, dimension(ndegr-1,ndegr-1), intent(out)::mtbbq, mtbbi
!
integer, dimension(1:nvqua)::ipq
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8,dimension(1:2, 1:nvqua)::xpqij,xpqlj
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:ndegr)::bq,bqpl,bqol
real*8,dimension(ndegr-1, ndegr-1)::bb, bbq, bbi,bbqi
real*8,dimension(1:nq, 1:ndegr)::rhs

real*8:: c00, c05, c10, c20, epsil
real*8:: xclj, yclj
real*8:: dxdri,dxdsi,dydri,dydsi, djacoi, rhogi
real*8:: dxclj,dyclj,detma
real*8:: xcrho, ycrho
real*8:: xgausl, ygausl, masel, xgausi, ygausi
real*8:: wi, r, s, rcj, scj,dr,ds
!
integer:: ies, jelem,  ishp
integer:: idegr, jdegr, igaus, je,im

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...I.3 Mass matrix
dr = 1.d0
ds = 1.d0

do ies = 1, 9

!...Element No.
jelem = ielaj(ies)

!...Local quad No.
je = jelem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,je)

!...Initial
xpqij(1, 1:nvqua) = cooro(1, ipq(1:nvqua))
xpqij(2, 1:nvqua) = cooro(2, ipq(1:nvqua))

!...Current
xpqlj(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpqlj(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...physical mass center
xclj =geoph(1, jelem)
yclj =geoph(2, jelem)

!...Local the mass center
call getvector_local(mapmt, xclj, yclj)

!...scale
dxclj =sqrt(geoel(3, jelem))
dyclj =dxclj

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

!...mass
masel = geoel(4, jelem)

!...Initial zero
bb = 0.d0
bbq= 0.d0
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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bq(1) = 1.d0
bq(2) = (r-rcj)/dr
bq(3) = (s-scj)/ds
!
if(npoly.eq.2)then
bq(4) = 0.5d0*bq(2)*bq(2) - geoel(19, jelem)
bq(5) = 0.5d0*bq(3)*bq(3) - geoel(20, jelem)
bq(6) =       bq(2)*bq(3) - geoel(21, jelem)
endif


bqpl(1) = 1.d0
bqpl(2) = (xgausl-xclj)/dxclj
bqpl(3) = (ygausl-yclj)/dyclj

!...Orthogonal basis functions
call basis_orthg(npoly, ndegr, nbsis, bqpl, geopl(3:16, ies), bqol)

!...Matrix(5x5)
do idegr =2, ndegr
do jdegr =2, ndegr
bb(idegr-1, jdegr-1)  =  bb(idegr-1, jdegr-1) + rhogi*djacoi*bq(idegr)*bqol(jdegr)
bbq(idegr-1, jdegr-1) = bbq(idegr-1, jdegr-1) + rhogi*djacoi*bqol(idegr)*bqol(jdegr)
enddo
enddo
!
enddo

!...RHS
rhs = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
rhs(1:nq, idegr) = rhs(1:nq, idegr) + unkrl(jdegr, 1:nq, ies)*bb(jdegr-1, idegr-1)
enddo
enddo


!...Invert matrix
call matinv225(bbq, bbqi, ndegr-1)
call matinv225(bb,  bbi, ndegr-1)

!...LHS
unkpl(:, 1:nq, ies) = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
unkpl(idegr, 1:nq, ies) = unkpl(idegr, 1:nq, ies) + rhs(1:nq, jdegr)*bbqi(jdegr-1,idegr-1)
enddo
enddo

!
!if((ielem.eq.9.or.ielem.eq.16).and.ies.eq.1)then
! do idegr=1,ndegr-1
!   print*,'matbb',bb(idegr,:)
! enddo

!print*,'unknp',ielem,unkpl(:,1,ies)

!do idegr=1,ndegr-1
!  print*,'matin',bbi(idegr,:)
!enddo
!endif

!...Store the mapping matrix
if(ies.eq.1)then
mtbbq = bbq
mtbbi = bbi
endif

!...Cell average on physical domain
unkpl(1, 1:nq, ies) = unkrl(1, 1:nq, ies)

enddo

end subroutine getunknp_local_Taylor
!
!...WENO procedure for one cell
!
subroutine weno_cell_Taylor(nbsis, isten, nsten, ielem, xpq, mapmt, esqua, geopl, geoel, geopj, unksl, unknr, unk_phy)
!
use constant
implicit none
!
integer,                               intent(in)::nbsis, ielem, isten, nsten
real*8,dimension(1:2, 1:nvqua),        intent(in)::xpq
real*8,dimension(1:2,1:2),             intent(in)::mapmt
integer, dimension(1:nfqua,1:nquad),   intent(in)::esqua
real*8,dimension(1:16),                intent(in)::geopl
real*8,dimension(1:ngeel),             intent(in)::geoel
real*8,dimension(1:3),                 intent(in)::geopj
real*8,dimension(1:ndegr,1:nq,1:nsten),intent(in)::unksl
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unknr
real*8,dimension(1:ndegr,1:nq),        intent(out)::unk_phy

!...Local
!...Local array
integer, parameter:: nfprj=10
integer:: mapfe(1:2,1:nfqua)
real*8::xpf(1:2, 1:2)
real*8, dimension(1:nq)::unkl,unkr
real*8:: weigh(1:nq, nsten)
real*8:: weigt(1:nq)
real*8:: os(1:nq, nsten)
real*8:: weigl(nsten)
real*8:: weigf(nq, nfprj)
real*8,dimension(1:4, 1:4)::qmat, qinvm
real*8,dimension(1:ndegr,1:nq,1:nsten)::unkf_phy,unksf_cha, unkf_cha
!
real*8::os1(nsten)
real*8::smthid(5)
real*8::cothg(6,6)
real*8::ctsft(3)
!
integer:: ifprj, is, iq,ies, jelem,id ,ifa, ie
real*8:: c00
real*8:: dtx, dty, dlgt, dnx, dny
real*8:: dxc, dyc, volel
real*8:: rpowe, epsil

!...Specify No.
c00 = 0.d0
rpowe =2.d0
epsil = 1.d-6
!...
ie = ielem - ntria

!...
volel = geoel(3)
dxc = sqrt(volel)
dyc = dxc

!...High-order orthogonal coefficients
if(nbsis.eq.1)then
ctsft(1:3) = geopl(3:5)

elseif(nbsis.eq.2)then
cothg(3,1) = geopl(3)
cothg(3,2) = geopl(4)

cothg(4,1) = geopl(5)
cothg(4,2) = geopl(6)
cothg(4,3) = geopl(7)

cothg(5,1) = geopl(8)
cothg(5,2) = geopl(9)
cothg(5,3) = geopl(10)
cothg(5,4) = geopl(11)

cothg(6,1) = geopl(12)
cothg(6,2) = geopl(13)
cothg(6,3) = geopl(14)
cothg(6,4) = geopl(15)
cothg(6,5) = geopl(16)
endif

!...Mapping array
mapfe(1, 1) = 1; mapfe(2, 1) = 2;
mapfe(1, 2) = 2; mapfe(2, 2) = 3;
mapfe(1, 3) = 3; mapfe(2, 3) = 4;
mapfe(1, 4) = 4; mapfe(2, 4) = 1;

!...  b. curvatures for the face-neighboring cells
ifprj = 0
unkf_cha = 0.d0
unkf_phy = 0.d0

!...Loop over the faces
do ies = 1, 4

jelem = esqua(ies,ie)
if(jelem .le. ncell) then

!...Find the projection vector
ifprj = ifprj +1

!...normal vector only using linear face even for curved meshes
xpf(1, 1:2) = xpq(1, mapfe(1:2, ies))
xpf(2, 1:2) = xpq(2, mapfe(1:2, ies))

dtx = xpf(1, 2) - xpf(1, 1)
dty = xpf(2, 2) - xpf(2, 1)

dlgt = sqrt(dtx**2 + dty**2)

dnx = dty/dlgt
dny =-dtx/dlgt
!
call getvector_local(mapmt, dnx,dny)

!...left and right unknowns
unkl(1:nq) = unknr(1, 1:nq, ielem)
unkr(1:nq) = unknr(1, 1:nq, jelem)
!
call getvector_local(mapmt, unkl(2), unkl(3))
call getvector_local(mapmt, unkr(2), unkr(3))

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
call getmatrix_prj2(qmat, qinvm, dnx, dny, unkl,unkr,ielem)

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

!...Smooth indicator for every stencle

!...Taylor basis
if(nbsis.eq.1)then

if(npoly.eq.2)then
do is= 1, isten
do iq= 1, nq
os(iq, is) = ((unksf_cha(2, iq, is)/dxc)**2 +&
(unksf_cha(3, iq, is)/dyc)**2)
enddo
enddo
elseif(npoly.eq.1)then
do is= 1, isten
do iq= 1, nq
!...
smthid(1) = 1.d0/(dxc**2)*(unksf_cha(2, iq, is)**2 +&
2.d0*ctsft(1)*unksf_cha(4, iq, is)**2 +&
2.d0*ctsft(2)*unksf_cha(6, iq, is)**2 +&
2.d0*ctsft(3)*unksf_cha(4, iq, is)*unksf_cha(6, iq, is))

smthid(2) = 1.d0/(dyc**2)*(unksf_cha(3, iq, is)**2  +&
2.d0*ctsft(2)*unksf_cha(5, iq, is)**2 +&
2.d0*ctsft(1)*unksf_cha(6, iq, is)**2 +&
2.d0*ctsft(3)*unksf_cha(5, iq, is)*unksf_cha(6, iq, is))

smthid(3) = 1.d0/(dxc**4)*unksf_cha(4, iq, is)**2
smthid(4) = 1.d0/(dyc**4)*unksf_cha(5, iq, is)**2
smthid(5) = 1.d0/(dyc*dxc)**2*unksf_cha(6, iq, is)**2
!
os(iq, is) = (smthid(1) + smthid(2)) +&
(smthid(3) + smthid(4) + smthid(5))*volel
enddo
enddo
endif

!...Orthonomal basis
elseif(nbsis.eq.2)then

if(npoly.eq.1)then
do is= 1, isten
do iq= 1, nq
os(iq, is) = (((unksf_cha(2, iq, is)+cothg(3, 1)*unksf_cha(3, iq, is))/dxc)**2 +&
(unksf_cha(3, iq, is)/dyc)**2)!/(abs(unksf_cha(1, iq, 1))+epsil)**1
enddo
enddo
elseif(npoly.eq.2)then
do is= 1, isten
do iq= 1, nq
!...
smthid(1) = 1.d0/(dxc**2)*(unksf_cha(2, iq, is) +&
cothg(3, 1)*unksf_cha(3, iq, is) +&
cothg(4, 1)*unksf_cha(4, iq, is) +&
cothg(5, 2)*unksf_cha(5, iq, is) +&
cothg(6, 3)*unksf_cha(6, iq, is))**2 +&
1.d0/(dxc**2)*geopj(1)*(2.d0*unksf_cha(4, iq, is)+&
2.d0*cothg(5, 1)*unksf_cha(5, iq, is)+&
2.d0*cothg(6, 1)*unksf_cha(6, iq, is))**2 +&
1.d0/(dxc**2)*geopj(2)*(unksf_cha(5, iq, is)+&
cothg(6, 2)*unksf_cha(6, iq, is))**2 +&
2.d0/(dxc**2)*geopj(3)*(unksf_cha(5, iq, is)+&
cothg(6, 2)*unksf_cha(6, iq,is))*(2.d0*unksf_cha(4, iq, is)+&
2.d0*cothg(5, 1)*unksf_cha(5, iq, is)+&
2.d0*cothg(6, 1)*unksf_cha(6, iq, is))
!
smthid(2) = 1.d0/(dyc**2)*(unksf_cha(3, iq, is) +&
cothg(4, 2)*unksf_cha(4, iq, is) +&
cothg(5, 3)*unksf_cha(5, iq, is) +&
cothg(6, 4)*unksf_cha(6, iq, is))**2+&
1.d0/(dxc**2)*geopj(1)*(unksf_cha(5, iq, is)+&
cothg(6, 2)*unksf_cha(6, iq, is))**2 +&
1.d0/(dxc**2)*geopj(2)*(2.d0*unksf_cha(6, iq, is))**2 +&
2.d0/(dxc**2)*geopj(3)*(unksf_cha(5, iq, is)+&
cothg(6, 2)*unksf_cha(6, iq,is))*(2.d0*unksf_cha(6, iq, is))

!
smthid(3) = 1.d0/(dxc**4)*(2.d0*unksf_cha(4, iq, is) + 2.d0*cothg(5, 1)*unksf_cha(5, iq, is)+&
2.d0*cothg(6, 1)*unksf_cha(6, iq, is) )**2
smthid(4) = 1.d0/(dyc**4)*(2.d0*unksf_cha(6, iq, is))**2
smthid(5) = 1.d0/(dyc*dxc)**2*(unksf_cha(5, iq, is) + cothg(6, 2)*unksf_cha(6, iq, is))**2
!
os(iq, is) = (smthid(1) + smthid(2)) +&
(smthid(3) + smthid(4) + smthid(5))*volel

!os(iq, is) =min((smthid(1) + smthid(2)),(smthid(3) + smthid(4) + smthid(5))*volel)
!if(ielem.eq.7.and.is.eq.1) print*,'derivative',iq,smthid(1:5)
enddo
enddo

endif

endif

!...Mometum
!os1(:) = 0.25d0*(os(1, :)+os(2, :)+os(3, :)+os(4, :))
os1(:) = max(os(1, :),os(2, :),os(3, :),os(4, :))

do is= 1, isten
os(1:nq, is) = os1(is)
enddo

!...Linear weight for central and biased stencle
weigl(1)=.2d0;
weigl(2:isten)= 0.8d0/(isten-1.d0)


weigt = 0.d0
do is= 1, isten
do iq =1,nq
weigh(iq, is) = weigl(is)/(epsil+os(iq,is))**rpowe
enddo
enddo

!if(ielem.eq.77)then
!print*,'ielem12',ielem,qinvm(1,:)
!print*,'ielem13',ifprj,os(1:4,5),epsil+os(1:4,5)
!endif
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
!
!if(ielem.eq.9.or.ielem.eq.16)then
!print*,'weno-weight',ielem,weigh(1,1:isten),os(1,1:isten)
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

!... Compute the average weight for all the faces
weigf=1.d0/ifprj

!... finally, we can reconstruct the physical polynomails for one cell
unk_phy = c00
do ifa = 1, ifprj
do iq = 1, nq
unk_phy(2:ndegr, iq) =unk_phy(2:ndegr, iq) + weigf(iq, ifa)*unkf_phy(2:ndegr,iq,ifa)
enddo
enddo

!if(ielem.eq.1921)print*,'weno',ielem,unksl(2,1, 1:isten)

end subroutine weno_cell_Taylor
!
!...WENO stnecils
!
subroutine weno_stncl_quad_orthol3(ielem, ipqua, esqua, unkno, unknr, geoel, geoph, &
coord, cooro, esuv1, esuv2, mapmt)
!
use constant
implicit none
!
integer,    intent(in):: ielem
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(inout)::unkno
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unknr
real*8,dimension(1:ngeel,1:nsize),    intent(in) ::geoel
real*8,dimension(1:2,1:nsize),        intent(in) ::geoph
real*8,dimension(1:ndimn,1:npoin),    intent(in) :: coord, cooro
integer*4, intent(in)::esuv1(npoin1),esuv2(npoin2)
real*8,dimension(1:ndimn,1:ndimn),    intent(in) :: mapmt

integer, parameter:: nsten=10, nfprj=10, nevaj=9
integer ::isten, neqlg, neqhm, neqhm_l,neqhm_h,ieshm
integer,dimension(5)::jelaj
integer,dimension(9)::ielaj
integer,dimension(1:nvqua):: ipq
integer::ielst1(4,7),ielst2(4,5),ielst3(4,4)
integer::idxiv(3)

real*8,dimension(1:ndegr,1:nq,1:nsten)::unksl
real*8,dimension(1:16,1:nevaj)::geopl
real*8,dimension(1:3, 1:nevaj)::geopj
real*8,dimension(1:ndegr,1:nq,1:nevaj)::unkpl, unkrl
real*8,dimension(1:nvqua)::shpq, dsprq, dspsq
real*8:: weighq(ngausdq), posiq(2, ngausdq)
real*8,dimension(1:2, 1:nvqua)::xpq,xpqij,xpqlj
real*8,dimension(1:ndegr)::bqpl,bqol,bqolj,bqplj
real*8,dimension(1:ndegr,1:nq)::unk_phy
real*8,dimension(1:ndegr-1,1:ndegr-1)::mtbbq,mtbbi
real*8,dimension(1:nq,1:ndegr-1)::rhsuk

real*8, allocatable::matlg(:, :),rhslg(:,:)
real*8, allocatable::mathm(:, :),rhshm(:,:)
real*8, allocatable::matls(:, :), matin(:, :),rhsls(:,:)

!
real*8:: c00, c05, c10, c20, epsil
real*8:: dxcl, dycl,dxclj,dyclj, xmcl, ymcl, xclj, yclj
real*8:: dxdri,dxdsi,dydri,dydsi,wi,r,s,xcrho,ycrho
real*8:: rcj,scj,xgausi,ygausi,xgausl,ygausl
real*8:: djacoi,rhogi, masej
!
integer:: ie,ies, je, jelem,  ishp,  istor,iv,nvaj, ivaj, ifaj
integer:: idegr, is, iq,igaus, jdegr, jeloc
integer:: idxod, nbsis

!
data c00   / 0.0d0    /
data c05   / 0.5d0    /
data c10   / 1.0d0    /
data c20   / 2.0d0    /
data epsil / 1.0d-6   /

!...Allocate memory
nbsis =1

!...Allocate memory
if(npoly.eq.1)then
neqlg = 3
elseif(npoly.eq.2)then
!...Lagrangian weno
neqlg  = 5
!...Hermite weno
neqhm   = 9
neqhm_l = 3
neqhm_h = neqhm - neqhm_l
endif

allocate(matlg(neqlg, ndegr-1), rhslg(nq, neqlg))
allocate(mathm(neqhm, ndegr-1), rhshm(nq, neqhm))
allocate(matls(ndegr-1, ndegr-1), matin(ndegr-1, ndegr-1),rhsls(nq, ndegr-1))

!...Specification: Stecil selection
ielst1(1, 1:7) =(/1,2,5,6,9,8,3/)
ielst1(2, 1:7) =(/1,3,2,7,6,9,5/)
ielst1(3, 1:7) =(/1,4,3,8,7,6,5/)
ielst1(4, 1:7) =(/1,5,4,9,8,7,3/)
!

ielst2(1, 1:5) =(/1,2,5,9,6/)
ielst2(2, 1:5) =(/1,3,2,6,7/)
ielst2(3, 1:5) =(/1,4,3,7,8/)
ielst2(4, 1:5) =(/1,5,4,8,9/)

ielst2(1, 1:4) =(/1,9,5,2/)
ielst2(2, 1:4) =(/1,8,5,4/)
ielst2(3, 1:4) =(/1,6,3,2/)
ielst2(4, 1:4) =(/1,7,3,4/)

ielst2(1, 1:4) =(/1,2,6,9/)
ielst2(2, 1:4) =(/1,3,6,7/)
ielst2(3, 1:4) =(/1,4,8,7/)
ielst2(4, 1:4) =(/1,5,8,9/)

!...Find weight and position for gauss points...
call ruqope(2, ngausdq, posiq, weighq)

!...Local quad No.
ie = ielem - ntria

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ie)

!...Coordinates
xpq(1, 1:nvqua) = coord(1, ipq(1:nvqua))
xpq(2, 1:nvqua) = coord(2, ipq(1:nvqua))

!...Scaling parameter
dxcl = sqrt(geoel(3, ielem))
dycl = dxcl

!...Adjacent cells
ielaj = 0
ielaj(1) = ielem
ielaj(2:5) = esqua(1:4, ielem)

!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
do iv = 1, 4

jelaj = 0
nvaj = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielem)then
nvaj = nvaj +1
jelaj(nvaj) = jelem
endif
enddo !istor
!...Stencil No.
idxod = 1
idxiv = 0
do ivaj=1,3
do ifaj=1,4
if(jelaj(ivaj).eq.ielaj(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
endif
enddo
enddo

!...Order the adjacent cells
do ivaj=1,3
if(idxiv(ivaj) == 0)then
if(idxod.eq.2)then
ielaj(6) = jelaj(ivaj)
elseif(idxod.eq.6)then
ielaj(7) = jelaj(ivaj)
elseif(idxod.eq.12)then
ielaj(8) = jelaj(ivaj)
elseif(idxod.eq.4)then
ielaj(9) = jelaj(ivaj)
else
stop
endif
endif
enddo
!
enddo !do iv = 1, 4
!
if(minval(ielaj(1:9)).eq.0)then
!...Boundary cell
isten = 1

else
!...Store the local reference solution
unkrl(:, :, 1:nevaj) = unknr(:, :, ielaj(1:nevaj))

!...Project the velocity-related unknown to local basis
do ies =1 ,nevaj
do idegr =1, ndegr
call getvector_local(mapmt, unkrl(idegr, 2, ies),unkrl(idegr, 3, ies))
enddo
enddo

!...Get the local orthonoaml coefficients
call getgeophy_local_Taylor(nbsis, ipqua,geoel, geoph, ielaj, mapmt, geopl, geopj, coord, cooro)

!...Get the local unknown
call getunknp_local_Taylor(nbsis, ielem,nsten,ipqua,geoel, geoph, ielaj, mapmt, geopl, coord, cooro, unkpl,unkrl, &
mtbbq, mtbbi)

!..Store the target cell at the 1st stencil
!...Zero out other stencils
unksl = 0.d0
isten =1
!...Cell average
do is = 1, nsten
unksl(1, :, is) =  unkpl(1, :, 1)
enddo

!...1st stencil
do idegr= 1, ndegr
unksl(idegr, :, 1) =  unkpl(idegr, :, 1)
enddo

!...1st stencil
do is =1 ,4
!
isten = isten +1
!...Initialize matlg and rhslg
matlg = 0.d0
rhslg = 0.d0
!
do ies =1, neqlg
!
jelem = ielaj(ielst1(is, ies+1))

je = jelem - ntria

!...Coordinates
xpqij(1, 1:nvqua) = cooro(1, ipqua(1:nvqua,je))
xpqij(2, 1:nvqua) = cooro(2, ipqua(1:nvqua,je))
!
xpqlj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,je))
xpqlj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,je))

!if(ielem.eq.1921.and.isten.eq.2.and.ies.eq.1)then
!  print*,'igau1',xpqlj(:,1:8)
!endif
!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!if(ielem.eq.1921.and.isten.eq.2.and.ies.eq.1)then
!  print*,'igau2',xpqlj(:,1:8)
!endif

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...mass
masej = geoel(4, jelem)

!...Local cell center for Target cell
xmcl = geopl(1, 1)
ymcl = geopl(2, 1)

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqplj(1) = 1.d0
bqplj(2) = (xgausl-xmcl)/dxcl
bqplj(3) = (ygausl-ymcl)/dycl

!...Orthogonal basis functions
call basis_orthg(npoly, ndegr, nbsis, bqplj, geopl(3:16, 1), bqolj)

!...Matrix
matlg(ies, 1)  = matlg(ies, 1) + rhogi*djacoi*bqolj(2)/masej
matlg(ies, 2)  = matlg(ies, 2) + rhogi*djacoi*bqolj(3)/masej
!
if(npoly.eq.2)then
matlg(ies, 3)  = matlg(ies, 3)  +  rhogi*djacoi*bqolj(4)/masej
matlg(ies, 4)  = matlg(ies, 4)  +  rhogi*djacoi*bqolj(5)/masej
matlg(ies, 5)  = matlg(ies, 5)  +  rhogi*djacoi*bqolj(6)/masej
endif
!
enddo !...igaus

!...RHS
rhslg(1:nq, ies) = unkpl(1, 1:nq, ielst1(is, ies+1))- unkpl(1, 1:nq, 1)

enddo

!...Symmetrize square matrix
call  matrix_sym(neqlg, ndegr-1, nq, matlg, rhslg, matls, rhsls)

!...Invert matrix
call matinv225(matlg,matin,ndegr-1)

!...Check

!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unksl(idegr, 1:nq, isten) =  unksl(idegr, 1:nq, isten) + matin(idegr-1, ies)*rhslg(1:nq, ies)
enddo
enddo
enddo

!...2nd stencil
!...The 2st 3 Lagrangian polynomial
if(npoly.eq.2)then

do is =1 ,4
!...Increasing stencils
isten = isten +1

!...Initial zero
mathm = 0.d0
rhshm = 0.d0

do ies=1, neqhm_l
!...Local No
jeloc = ielst2(is, ies+1)
!
jelem = ielaj(jeloc)
je = jelem - ntria
!...Coordinates
xpqij(1, 1:nvqua) = cooro(1, ipqua(1:nvqua,je))
xpqij(2, 1:nvqua) = cooro(2, ipqua(1:nvqua,je))
!
xpqlj(1, 1:nvqua) = coord(1, ipqua(1:nvqua,je))
xpqlj(2, 1:nvqua) = coord(2, ipqua(1:nvqua,je))

!...Local the cell vertex
do ishp = 1, nvqua
call getvector_local(mapmt, xpqlj(1, ishp),xpqlj(2, ishp))
enddo

!...mass center...
rcj= geoel(1, jelem)
scj= geoel(2, jelem)

!...mass
masej = geoel(4, jelem)

!...Scaling parameter
dxclj = sqrt(geoel(3, jelem))
dyclj = dxclj
!
xclj = geopl(1, ielst2(is, ies+1))
yclj = geopl(2, ielst2(is, ies+1))
!
xmcl = geopl(1, ielst2(is, 1))
ymcl = geopl(2, ielst2(is, 1))

!...Cell cenetr for density
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqij, rcj, scj, xcrho, ycrho)

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
dxdri = dxdri + dsprq(ishp)*xpqij(1,ishp)
dxdsi = dxdsi + dspsq(ishp)*xpqij(1,ishp)

dydri = dydri + dsprq(ishp)*xpqij(2,ishp)
dydsi = dydsi + dspsq(ishp)*xpqij(2,ishp)
enddo
!
djacoi = wi*(dxdri*dydsi - dydri*dxdsi)

!...Physical coord for current...
xgausi = 0.d0
ygausi = 0.d0

xgausl = 0.d0
ygausl = 0.d0
!
do ishp = 1, nvqua
xgausi = xgausi + shpq(ishp)*xpqij(1,ishp)
ygausi = ygausi + shpq(ishp)*xpqij(2,ishp)

xgausl  = xgausl + shpq(ishp)*xpqlj(1,ishp)
ygausl  = ygausl + shpq(ishp)*xpqlj(2,ishp)
enddo

!...The initial density at gauss points
call getrhog_initial(rhogi,  xgausi, ygausi, xcrho, ycrho)

!...Basis functions
bqplj(1) = 1.d0
bqplj(2) = (xgausl-xmcl)/dxcl
bqplj(3) = (ygausl-ymcl)/dycl

!...Orthogonal basis functions
call basis_orthg(npoly, ndegr, nbsis, bqplj, geopl(3:16, 1), bqol)

!...Matrix
mathm(ies, 1)  = mathm(ies, 1) + rhogi*djacoi*bqol(2)/masej
mathm(ies, 2)  = mathm(ies, 2) + rhogi*djacoi*bqol(3)/masej
!
if(npoly.eq.2)then
mathm(ies, 3)  = mathm(ies, 3)  +  rhogi*djacoi*bqol(4)/masej
mathm(ies, 4)  = mathm(ies, 4)  +  rhogi*djacoi*bqol(5)/masej
mathm(ies, 5)  = mathm(ies, 5)  +  rhogi*djacoi*bqol(6)/masej
endif

!...Hermite polynomial
if(ies.le.(neqhm_h/2))then
!
bqplj(1) = 1.d0
bqplj(2) = (xgausl-xclj)/dxclj
bqplj(3) = (ygausl-yclj)/dyclj

!...Orthogonal basis functions
call basis_orthg(npoly, ndegr, nbsis, bqplj, geopl(3:16, jeloc), bqolj)

ieshm = neqhm_l + 2*ies-1
!...Matrix
do idegr =2, ndegr
mathm(ieshm,   idegr-1)  = mathm(ieshm,   idegr-1) + &
rhogi*djacoi*bqol(idegr)*bqolj(2)/masej
mathm(ieshm+1, idegr-1)  = mathm(ieshm+1, idegr-1) + &
rhogi*djacoi*bqol(idegr)*bqolj(3)/masej
enddo

!
do idegr =2, ndegr
rhshm(1:nq, ieshm)   =  rhshm(1:nq, ieshm)   + rhogi*djacoi*bqolj(idegr)*bqolj(2)*unkpl(idegr,1:nq,jeloc)/masej
rhshm(1:nq, ieshm+1) =  rhshm(1:nq, ieshm+1) + rhogi*djacoi*bqolj(idegr)*bqolj(3)*unkpl(idegr,1:nq,jeloc)/masej
enddo

endif
!
enddo

!...RHS for Lag WENO
rhshm(1:nq, ies) = unkpl(1, 1:nq, jeloc)- unkpl(1, 1:nq, 1)

enddo

!...Symmetrize square matrix
call  matrix_sym(neqhm, ndegr-1, nq, mathm, rhshm, matls, rhsls)


!...Invert matrix
call matinv225(matls, matin, ndegr-1)

!...Update the stencil polynomial
do idegr = 2, ndegr
do ies   = 1, ndegr-1
unksl(idegr, 1:nq, isten) =  unksl(idegr, 1:nq, isten) + matin(idegr-1, ies)*rhsls(1:nq, ies)
enddo
enddo
!...Check
!if(ielem.eq.265.and.isten.eq.6)then
!do idegr=1,neqhm
!print*,'mat',isten,mathm(idegr,:),rhshm(1,idegr)
!enddo
!endif

enddo

endif !if(npoly.eq.2)then
!
!...WENO part
!
call weno_cell_Taylor(nbsis, isten, nsten,ielem, xpq, mapmt, esqua, geopl(:, 1), geoel(:, ielem), geopj(:, 1),&
unksl, unknr, unk_phy)

!
!...Inverse L2 projection
!
!...RHS
rhsuk = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
rhsuk(1:nq, idegr-1) = rhsuk(1:nq, idegr-1) + unk_phy(jdegr, 1:nq)*mtbbq(jdegr-1, idegr-1)
enddo
enddo

!...LHS
unkrl(2:ndegr, 1:nq, 1) = 0.d0
do idegr = 2, ndegr
do jdegr = 2, ndegr
unkrl(idegr, 1:nq, 1) = unkrl(idegr, 1:nq, 1) + rhsuk(1:nq, jdegr-1)*mtbbi(jdegr-1, idegr-1)
enddo
enddo

!...Check
!if(ielem.eq.289.or.ielem.eq.290)then
!do is=1,isten
!print*,'inweno2',ielem,is,unksl(2:ndegr,1,is)
!enddo
!print*,'inweno3',ielem,unk_phy(2:ndegr,1)
!print*,'geopl',ielem,geopl(1:16,1)
!print*,'theta',ielem,unk_phy(3,1)+geopl(6,1)*unk_phy(4,1)+geopl(10,1)*unk_phy(5,1)+geopl(15,1)*unk_phy(6,1)
!print*,'rthet',ielem,unk_phy(5,1)+geopl(13,1)*unk_phy(6,1)
!endif


!...Transfer back
do idegr =2, ndegr
call getvector_local2(1, mapmt, unkrl(idegr, 2, 1),unkrl(idegr, 3, 1))
enddo

!...Update the unknown in the reference coordinate
unkno(2:ndegr,1:nq,ielem) = unkrl(2:ndegr, 1:nq, 1)
!
endif !...ies, end of sencils


!
!if(ielem.eq.1921) print*,'ielem1921',unksl(2,1,8)

!...Release memory
deallocate(matlg, rhslg, mathm, rhshm,matls, matin, rhsls)

return
end subroutine weno_stncl_quad_orthol3
!
!...Calculate the 1D Riemann solution at the quadrature points...
!
subroutine get1DRiemann(geoel,bface,intfac,ipqua,coord,coold,unkno, ustar, &
afvec, aflim, rhsel, itime)
use constant
implicit none
!...Input arrays
integer*4,dimension(1:nbfai,nbfac),          intent(in)::bface
integer*4,dimension(1:nifai,1:nafac),        intent(in)::intfac
integer*4,dimension(1:nvqua,1:nquad),        intent(in)::ipqua
real*8,dimension(1:ngeel,1:nsize),           intent(in)::geoel
real*8,dimension(1:ndimn,1:npoin),           intent(in)::coord, coold
real*8,dimension(1:ndegr,1:nq,1:nsize),      intent(in)::unkno
real*8,dimension(1:ndimn,1:npoin),           intent(in)::ustar !...nodal velocity
real*8,dimension(1:nq+1,1:nsize),            intent(in)::aflim !...Limiter coef
real*8,dimension(1:2, 1:2, 1:nsize),         intent(in)::afvec
real*8,dimension(1:ndegr,1:nq,1:ncell),      intent(inout)::rhsel
integer*4,                      intent(in)::itime
!...Local integer
integer::ifa,iel,ier,ie,ig
integer::iv,ishp,ivf,ivq,idegr,iq
!...local integer array
integer,dimension(1:nvfac) :: ipf
integer,dimension(2)::fgausl, fgausr
!
real*8, dimension(1:nq)   ::unkgl, unkgr
real*8, dimension(1:nvqua)::rvq, svq
real*8, dimension(1:2,1:nvfac)::xpf(1:2, 1:nvfac)
real*8, dimension(1:nvfac)::shpf,dshprf
real*8, dimension(1:ndegr)::bql ,bqr
real*8, dimension(1:ndimn)::vewgt,ustrg
real*8, dimension(1:ngausf)::posi, weigh


!...Local real number
real*8::eps
real*8::rhomcr,rhoctr,uctrr,vctrr,ectrr,pctrr,sdctrr
real*8::rhomcl,rhoctl,uctrl,vctrl,ectrl,pctrl,sdctrl
real*8::muctl,muctr,deltul,deltur
real*8::rhoml,rhogl,ugl,vgl,egl,pgl,rhomr,rhogr,ugr,vgr,egr,pgr
real*8::vngl,vngr,vnstr, prstrl, prstrr, unstr,unstrr,unstrl
real*8::dr, ds,rcl,scl, rcr, scr,rgr,rgl,sgr,sgl
real*8::dudr,duds,dvdr,dvds
real*8::r, wi
real*8::djaco,dwav1,dwav2
real*8::dxdr,dydr

!...Initial velocity at gauss point
!..Initial Gauss point velocity from the interpolation of the 3-node velocity...
eps = 1.d-6

!...Quadrature point position and weight
! call ruqope_lobatto(1, ngausf, posi, weigh)
 call ruqope(1, ngausf, posi, weigh)

!...Scaling factors
dr = 1.d0
ds = 1.d0
!
!...Part I: Specify some gauss points
!
rvq(1) = -1.d0; svq(1) = -1.d0
rvq(2) =  1.d0; svq(2) = -1.d0
rvq(3) =  1.d0; svq(3) =  1.d0
rvq(4) = -1.d0; svq(4) =  1.d0
!
if(ncurv.eq.1)then
rvq(5) =  0.d0; svq(5) = -1.d0
rvq(6) =  1.d0; svq(6) =  0.d0
rvq(7) =  0.d0; svq(7) =  1.d0
rvq(8) = -1.d0; svq(8) =  0.d0
rvq(9) =  0.d0; svq(9) =  0.d0
endif
!
!...Part II: Get the corrected velocity at the non-vertex gauss point...
!
fgausl = 0.d0
fgausr = 0.d0

!...Loop over faces
do 450 ifa = 1, nafac !...(1)ifa = 1, nafac

!...let and right cell
iel = intfac(1, ifa)
ier = intfac(2, ifa)

!...mass center
rcl = geoel(1, iel)
scl = geoel(2, iel)
!
if(ifa.gt.nbfac)then
rcr = geoel(1, ier)
scr = geoel(2, ier)
endif
!...face vertices
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
xpf(1, 1:nvfac) = coord(1, ipf(1:nvfac))
xpf(2, 1:nvfac) = coord(2, ipf(1:nvfac))

!...local vertex no. for left and right cells
do ivf = 1, 2
do ivq = 1, 4
!
if(ipf(ivf).eq.ipqua(ivq, iel))then
fgausl(ivf) = ivq
endif
!...interior cells
if(ifa.gt.nbfac)then
if(ipf(ivf).eq.ipqua(ivq, ier))then
fgausr(ivf) = ivq
endif
endif
enddo
enddo

!...Cell average
rhomcl = unkno(1, 1, iel)
uctrl  = unkno(1, 2, iel)
vctrl  = unkno(1, 3, iel)
ectrl  = unkno(1, 4, iel)
!
rhoctl = 1.d0/rhomcl
pctrl  = max(eps,(gamlg-1.d0)*rhoctl*(ectrl-0.5d0*(uctrl**2 + vctrl**2)))
sdctrl = sqrt(gamlg*pctrl/rhoctl)

!...acoustic impedance for left cell
muctl = rhoctl*sdctrl

if(ifa.gt.nbfac)then
!...Cell average
rhomcr = unkno(1, 1, ier)
uctrr  = unkno(1, 2, ier)
vctrr  = unkno(1, 3, ier)
ectrr  = unkno(1, 4, ier)
!
rhoctr = 1.d0/rhomcr
pctrr  = max(eps, (gamlg-1.d0)*rhoctr*(ectrr-0.5d0*(uctrr**2 + vctrr**2)))
sdctrr= sqrt(gamlg*pctrr/rhoctr)

!...acoustic impedance for left cell
muctr = rhoctr*sdctrr
endif

!...loop quadrature points
do ig =1, ngausf
!
r   = posi(ig)
wi  = weigh(ig)

!...Shape function and its derivatives
call getshapfct_edge(ncurv,nvfac,shpf, dshprf, r)

!...Riemann velcoity
ustrg = 0.d0
!
do ishp = 1, nvfac
 ustrg(1) = ustrg(1) + shpf(ishp)*ustar(1, ipf(ishp))
 ustrg(2) = ustrg(2) + shpf(ishp)*ustar(2, ipf(ishp))
enddo

!...Jacobian determinant...
dxdr = 0.d0
dydr = 0.d0
!
do ishp = 1, nvfac
dxdr = dxdr + dshprf(ishp)*xpf(1, ishp)
dydr = dydr + dshprf(ishp)*xpf(2, ishp)
enddo

djaco = sqrt(dxdr**2 + dydr**2)

!...normal vector
dwav1 = dydr/djaco
dwav2 =-dxdr/djaco

!...Quadrature points
rgl = 0.d0
sgl = 0.d0
rgr = 0.d0
sgr = 0.d0
!
do ishp = 1, 2
rgl = rgl + shpf(ishp)*rvq(fgausl(ishp))
sgl = sgl + shpf(ishp)*svq(fgausl(ishp))

if(ifa.gt.nbfac)then
rgr = rgr + shpf(ishp)*rvq(fgausr(ishp))
sgr = sgr + shpf(ishp)*svq(fgausr(ishp))
endif
enddo

!...basis functions
bql(1) = 1.d0
bql(2) = (rgl-rcl)/dr
bql(3) = (sgl-scl)/ds
!DGP2
if(npoly.eq.2)then
bql(4) = 0.5d0*bql(2)*bql(2) - geoel(19 ,iel)
bql(5) = 0.5d0*bql(3)*bql(3) - geoel(20 ,iel)
bql(6) =       bql(2)*bql(3) - geoel(21 ,iel)
endif

if(ifa.gt.nbfac)then
bqr(1) = 1.d0
bqr(2) = (rgr-rcr)/dr
bqr(3) = (sgr-scr)/ds
!DGP2
if(npoly.eq.2)then
bqr(4) = 0.5d0*bqr(2)*bqr(2) - geoel(19 ,ier)
bqr(5) = 0.5d0*bqr(3)*bqr(3) - geoel(20 ,ier)
bqr(6) =       bqr(2)*bqr(3) - geoel(21 ,ier)
endif
endif

!...left and right state
unkgl = 0.d0
unkgr = 0.d0

do iq = 1, nq
do idegr =1, ndegr
unkgl(iq) = unkgl(iq) + unkno(idegr, iq, iel)*bql(idegr)
enddo
enddo

if(ifa.gt.nbfac)then
do iq = 1, nq
do idegr =1, ndegr
unkgr(iq) = unkgr(iq) + unkno(idegr, iq, ier)*bqr(idegr)
enddo
enddo
endif

!...Left state
rhoml = unkgl(1)
ugl   = unkgl(2)
vgl   = unkgl(3)
egl   = unkgl(4)
!
rhogl = 1.d0/rhoml
pgl = max(eps, (gamlg-1.d0)*rhogl*(egl-0.5d0*(ugl**2 + vgl**2)))

!...Limiter
if(nlimi.eq.6.and.npoly.ge.1)then
if(ndens.eq.1)then
rhoml = rhomcl + aflim(1, iel)*(rhoml - rhomcl)
rhogl = 1.d0/rhoml
endif
!
dudr = afvec(1, 1, iel)*unkno(2,2,iel) +  afvec(1, 2, iel)*unkno(2,3,iel)
duds = afvec(1, 1, iel)*unkno(3,2,iel) +  afvec(1, 2, iel)*unkno(3,3,iel)
dvdr = afvec(2, 1, iel)*unkno(2,2,iel) +  afvec(2, 2, iel)*unkno(2,3,iel)
dvds = afvec(2, 1, iel)*unkno(3,2,iel) +  afvec(2, 2, iel)*unkno(3,3,iel)
!
ugl = unkno(1,2,iel)  + dudr*bql(2) + duds*bql(3)
vgl = unkno(1,3,iel)  + dvdr*bql(2) + dvds*bql(3)
!
pgl= pctrl + aflim(4, iel)*(pgl - pctrl)
!
endif

!...Right state
if(ifa.le.nbfac)then

if(bface(3, ifa).eq.21)then
!...acoustic impedance for left cell
muctr = muctl

rhomr = rhoml
ugr   = ugl
vgr   = vgl
!
rhogr = 1.d0/rhomr
pgr = pgl

elseif(bface(3, ifa).eq.25)then
!...acoustic impedance for left cell
muctr = muctl

rhomr = rhoml
ugr   = 0.d0
vgr   = 0.d0
!
rhogr = 1.d0/rhomr
pgr = pgl

elseif(bface(3, ifa).eq.22)then
!...acoustic impedance for left cell
muctr = muctl

if(bface(4, ifa).eq.221)then
rhomr = rhoml
ugr   =  ugl
vgr   = -vgl
elseif(bface(4, ifa).eq.222)then
rhomr = rhoml
ugr   = -ugl
vgr   =  vgl
endif
!
rhogr = 1.d0/rhomr
pgr = pgl

endif

elseif(ifa.gt.nbfac)then

rhomr = unkgr(1)
ugr   = unkgr(2)
vgr   = unkgr(3)
egr   = unkgr(4)
!
rhogr = 1.d0/rhomr
pgr = max(eps, (gamlg-1.d0)*rhogr*(egr-0.5d0*(ugr**2 + vgr**2)))

!...Limiter
if(nlimi.eq.6.and.npoly.ge.1)then
 if(ndens.eq.1)then
  rhomr = rhomcr + aflim(1, ier)*(rhomr - rhomcr)
  rhogr = 1.d0/rhomr
 endif
!
 dudr = afvec(1, 1, ier)*unkno(2,2,ier) +  afvec(1, 2, ier)*unkno(2,3,ier)
 duds = afvec(1, 1, ier)*unkno(3,2,ier) +  afvec(1, 2, ier)*unkno(3,3,ier)
 dvdr = afvec(2, 1, ier)*unkno(2,2,ier) +  afvec(2, 2, ier)*unkno(2,3,ier)
 dvds = afvec(2, 1, ier)*unkno(3,2,ier) +  afvec(2, 2, ier)*unkno(3,3,ier)
!
 ugr = unkno(1,2,ier)  + dudr*bqr(2) + duds*bqr(3)
 vgr = unkno(1,3,ier)  + dvdr*bqr(2) + dvds*bqr(3)
!
 pgr= pctrr + aflim(4, ier)*(pgr - pctrr)
!
endif
!
endif

!...normal velocity at left and right state
vngl     = ugl*dwav1 + vgl*dwav2
vngr     =-ugr*dwav1 - vgr*dwav2

!...normal Riemann velocity from Lag Riemann solver
unstr    = ustrg(1)*dwav1 + ustrg(2)*dwav2

!...acoustic impedance
deltul = abs(unstr-vngl)
!deltul = sqrt((ustrg(1)-ugl)**2 + (ustrg(2)-vgl)**2)
muctl = rhoctl*sdctrl + cimpd*rhoctl*slpdu*deltul

deltur =abs(-unstr-vngr)
!deltur = sqrt((ustrg(1)-ugr)**2 + (ustrg(2)-vgr)**2)
muctr = rhoctr*sdctrr + cimpd*rhoctr*slpdu*deltur

!...1D Riemann solution
vewgt(1) = (muctl*ugl + muctr*ugr)/(muctl + muctr)
vewgt(2) = (muctl*vgl + muctr*vgr)/(muctl + muctr)

!...normal Rieman velocity from 1D Riemann solver
vnstr    = vewgt(1)*dwav1 + vewgt(2)*dwav2 - (pgr-pgl)/(muctl + muctr)

!...BC for smooth noh
if(ifa.le.nbfac)then
if(bface(3, ifa).eq.27)then
 vnstr = unstr
endif
endif

!...normal Riemann velocity from Lag at left and right state
unstrl = vnstr
unstrr =-vnstr

!...Riemann pressure from 1D Riemann solver
prstrl    = pgl - muctl*(unstrl - vngl)
prstrr    = pgr - muctr*(unstrr - vngr)

!...Add the gauss distribution for rhs
djaco = djaco*wi

rhsel(1:ndegr, 1, iel) =  rhsel(1:ndegr, 1, iel) + unstr*djaco*bql(1:ndegr)
rhsel(1:ndegr, 2, iel) =  rhsel(1:ndegr, 2, iel) - prstrl*dwav1*djaco*bql(1:ndegr)
rhsel(1:ndegr, 3, iel) =  rhsel(1:ndegr, 3, iel) - prstrl*dwav2*djaco*bql(1:ndegr)
rhsel(1:ndegr, 4, iel) =  rhsel(1:ndegr, 4, iel) - prstrl*unstr*djaco*bql(1:ndegr)

!...Right cell
if(ifa.gt.nbfac)then
rhsel(1:ndegr, 1, ier) =  rhsel(1:ndegr, 1, ier) - unstr*djaco*bqr(1:ndegr)
rhsel(1:ndegr, 2, ier) =  rhsel(1:ndegr, 2, ier) + prstrr*dwav1*djaco*bqr(1:ndegr)
rhsel(1:ndegr, 3, ier) =  rhsel(1:ndegr, 3, ier) + prstrr*dwav2*djaco*bqr(1:ndegr)
rhsel(1:ndegr, 4, ier) =  rhsel(1:ndegr, 4, ier) + prstrr*unstr*djaco*bqr(1:ndegr)
endif

!...Add the gauss distribution for rhs of velocity
!if(iel.eq.1276.or.ier.eq.1276)print*,'ustar',ifa,iel,ier,djaco,dwav1,dwav2,wi,ig,prstrr,prstrl,unstr,vnstr

!
enddo

450 enddo
!
end subroutine get1DRiemann
!
!...subroutine: Get the nodal velocity U_p^* and pressure for hybrid meshes with RSF...
!
subroutine getRiemvtx_lag_rfmfem(gflag,gelag,gelagq,geoel,bface,intfac,inpoel,iptri,ipqua,&
coord, coold, unkno,ustar, fstar, fstarq, aflim, afvec, rhsel, itime)
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
real*8,dimension(1:ndegr,1:nq,1:ncell),     intent(out)::rhsel
integer,                                      intent(in):: itime
!...Local integer
integer::ie,ig,ideg,jdeg, iv, ifa, ipoin,ielem, iloop
!...local integer array
integer,dimension(1:nvtri) :: ipt
integer,dimension(1:nvqua) :: ipq
integer,dimension(1:nvfac) :: ipf
integer,dimension(npoin)   ::idxpt

!...local real array
real*8::munaci(2, 2)
real*8,dimension(1:3,1:2,1:nvqua)::anmq
real*8,dimension(1:ndimn,1:npoin)::vlave
real*8::radbf(nvfac),radbfx(nvfac), radbfy(nvfac),verad(nvfac)
!...local real number
real*8::eps,c00,c05,c10,c20
real*8::anx,any
real*8::rc, sc, dr, ds
real*8::detma,rhsu1,rhsu2,dtime
!
real*8,allocatable::sgmst(:,:,:),lhsrm(:),lhsrmu(:),rhust(:,:)
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
!...Part I: Specify some preliminary work
!
allocate(sgmst(2,2,npoin),lhsrm(npoin),lhsrmu(npoin),rhust(2,npoin))
allocate (munacn(1:2, 1:2, 1:npoin))
allocate (munacu(1:ndimn, 1:npoin), snsigm(1:ndimn, 1:npoin))
allocate (munaclq(1:2, 1:2, 1:2, 1:nvqua, 1:nquad), munaulq(1:ndimn, 1:2, 1:nvqua,  1:nquad),&
snsigmlq(1:ndimn, 1:2,  1:nvqua,  1:nquad))

!...Zero out vlave (averaged velocity)
vlave = 0.d0
idxpt = 0

!...Mark the boundary nodes...
!...Shockless Noh...
if(ncase.eq.2)then
do ifa = 1, nbfac
ipf(1:nvfac) = intfac(3:(2+nvfac), ifa)
idxpt(ipf(1:nvfac)) = 1
enddo
endif
!
!...Part II: Loop to get the information from Riemann solver
!
do iloop= 1, 1

!...Give vlave
vlave= ustar

!...Tria

!...Quad
!if(nquad.gt.0) call getriem_quad_rsf(ipqua, geoel, gelagq, vlave, unkno, ustar, rhust,&
!sgmst,lhsrm,lhsrmu,coord, coold, aflim, afvec)

!do ipoin = 1, npoin
! sgmst(1:2,1:2, ipoin) = sgmst(1:2,1:2, ipoin)/lhsrm(ipoin)
!enddo
!call getriem_quad_velo_rsf2(ipqua, geoel, gelagq, vlave, unkno, ustar,rhust,&
!sgmst,lhsrmu,coord, coold, aflim, afvec)

!...Update the velocity and stress at the vertex...
!do ipoin = 1, npoin
! if(idxpt(ipoin).eq.0)then
!  ustar(1:2, ipoin)     = rhust(1:2, ipoin)/lhsrmu(ipoin)
! endif
!enddo

!...Initialize munacn, munacu, snsigm
munacn  = 0.d0
munacu  = 0.d0
snsigm  = 0.d0

!...Quad
if(nquad.gt.0) call getriem_quad_matrixsym(ipqua, geoel, gelagq, vlave, unkno, munacn, munacu, snsigm,&
munaclq, munaulq, snsigmlq, coord, coold, aflim, afvec)

!...Update the velocity at the vertex...
do ipoin = 1, npoin
if(idxpt(ipoin).eq.0)then
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

!
!...Part 3: Face integral...
!
!...Initialize the RHS
rhsel = 0.d0

call get1DRiemann(geoel,bface,intfac,ipqua,coord,coold,unkno, ustar, &
afvec, aflim, rhsel, itime)
!
!print*,'rhs',rhsel(1,:,1),ustar(:,8)

!
!...Part 4: Get the Riemann forces for face integral...
!
!...Tria
!...'Traingle will be implemented in future!'

!...Quad
do ie = 1, nquad

ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!...Give the normal vector of every face...
anmq(1:3, 1, 1) = gelagq(1:3, 4, ie); anmq(1:3, 2, 1) = gelagq(1:3, 1, ie)
anmq(1:3, 1, 2) = gelagq(1:3, 1, ie); anmq(1:3, 2, 2) = gelagq(1:3, 2, ie)
anmq(1:3, 1, 3) = gelagq(1:3, 2, ie); anmq(1:3, 2, 3) = gelagq(1:3, 3, ie)
anmq(1:3, 1, 4) = gelagq(1:3, 3, ie); anmq(1:3, 2, 4) = gelagq(1:3, 4, ie)
!

!...Riemann forces
do iv = 1, nvqua
do ifa =1, 2

anx = anmq(1, ifa, iv)*anmq(3, ifa, iv)*0.5d0
any = anmq(2, ifa, iv)*anmq(3, ifa, iv)*0.5d0

fstarq(1, ifa, iv, ie) = sgmst(1, 1, ipq(iv))*anx + sgmst(1, 2, ipq(iv))*any
fstarq(2, ifa, iv, ie) = sgmst(2, 1, ipq(iv))*anx + sgmst(2, 2, ipq(iv))*any
enddo
enddo
!

enddo

!...Free the allocatable arrays
deallocate(sgmst, lhsrm, lhsrmu, rhust)
deallocate(munacn, munacu, snsigm, munaclq, snsigmlq, munaulq)

end subroutine getRiemvtx_lag_rfmfem
!
!...Find the WENO stencils
!
subroutine weno_stencils(nevaj,ielem, ipqua, esqua,esuv1, esuv2, ielse)
!
use constant
implicit none
!
integer,  intent(in)::nevaj, ielem
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nfqua,1:nquad),  intent(in)::esqua
integer,    intent(in)::esuv1(npoin1),esuv2(npoin2)
integer, intent(inout)::ielse(nevaj)
!
integer::ielse2(nevaj)
integer,  dimension(1:nvqua):: ipq
integer,dimension(5)::ielsv
integer::idxiv(3)
integer::idxel(nsize)
!
integer:: jelem, istor,iv,nesv, ivaj, id, ifaj, ifai
integer:: idxod,iej

!...Preliminary
ipq(1:nvqua) = ipqua(1:nvqua, ielem)

!
!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
!
do iv = 1, 4
ielsv = 0
nesv = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielem)then
nesv = nesv +1
ielsv(nesv) = jelem
endif
enddo !istor

!...Excluding the face-adjacent cells
idxod = 1
idxiv = 0
do ivaj=1,nesv
do ifaj=1,4
if(ielsv(ivaj).eq.ielse(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
endif
enddo
enddo

!...Order the vertex-adjacent cells
do ivaj=1,nesv
if(idxiv(ivaj) == 0)then
if(idxod.eq.2)then
ielse(6) = ielsv(ivaj)
elseif(idxod.eq.6)then
ielse(7) = ielsv(ivaj)
elseif(idxod.eq.12)then
ielse(8) = ielsv(ivaj)
elseif(idxod.eq.4)then
ielse(9) = ielsv(ivaj)
else
!print*,'Wrong stencil!',ielem,idxod
!stop
endif
endif
enddo
!
enddo !do iv = 1, 4


!...Coloring the adjacent cells
if(nevaj.gt.9.and.minval(ielse(1:9)).ne.0)then
!
idxel = 0
idxel(ielse(1:9)) = 1
!...Next layer
do ifaj=1,4
do ifai=1,4
if(idxel(esqua(ifai, ielse(ifaj+1))).ne.1)then
if(ifaj.eq.1)then
ielse(10) = esqua(ifai, ielse(ifaj+1))
elseif(ifaj.eq.2)then
ielse(11) = esqua(ifai, ielse(ifaj+1))
elseif(ifaj.eq.3)then
ielse(12) = esqua(ifai, ielse(ifaj+1))
elseif(ifaj.eq.4)then
ielse(13) = esqua(ifai, ielse(ifaj+1))
endif
endif
enddo
enddo
!
if(nevaj.gt.13.and.(minval(ielse(1:9)).gt.0.or.maxval(ielse(1:nevaj)).le.ncell))then

!...2nd layer
do iej = 1, 4
!
ielse2(2:5) = esqua(1:4, ielse(iej+1))

!...Local vertices
ipq(1:nvqua) = ipqua(1:nvqua,ielse(iej+1))
!...Part 1: Loop over the 4 corner vertices for one quad to order the adjacent cells
do iv = 1, 4

ielsv = 0
nesv = 0
do istor=esuv2(ipq(iv))+1,esuv2(ipq(iv)+1)
jelem=esuv1(istor)
if(jelem.ne.ielse(iej+1))then
nesv = nesv+1
ielsv(nesv) = jelem
endif
enddo !istor
!...Stencil No.
idxod = 1
idxiv = 0
do ivaj=1,3
do ifaj=1,4
if(ielsv(ivaj).eq.ielse2(ifaj+1))then
idxod = idxod*ifaj
idxiv(ivaj) = 1
endif
enddo
enddo
!
!print*,'ielem',ielem,iv,iej,ielse(iej+1),idxod,ielse2(2:5),ielsv(:)

!...Order the adjacent cells
do ivaj=1,3
if(idxiv(ivaj) == 0)then
if(iej.eq.1)then !...if(iej.eq.1)then

if(idxod.eq.2)then
ielse(16) = ielsv(ivaj)
elseif(idxod.eq.4)then
ielse(15) = ielsv(ivaj)
endif

elseif(iej.eq.2)then !...if(iej.eq.2)then

if(idxod.eq.2)then
ielse(17) = ielsv(ivaj)
elseif(idxod.eq.6)then
ielse(18) = ielsv(ivaj)
endif

elseif(iej.eq.3)then !...if(iej.eq.2)then

if(idxod.eq.6)then
ielse(19) = ielsv(ivaj)
elseif(idxod.eq.12)then
ielse(20) = ielsv(ivaj)
endif

elseif(iej.eq.4)then !...if(iej.eq.2)then

if(idxod.eq.12)then
ielse(21) = ielsv(ivaj)
elseif(idxod.eq.4)then
ielse(14) = ielsv(ivaj)
endif

endif !...if(iej.eq.1)then
endif
enddo
!
enddo !do iv = 1, 4

enddo
endif

endif

end subroutine weno_stencils
