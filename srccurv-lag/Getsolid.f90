!
!...Get the deviatoric stress for solid dynamics with finite strain
!
subroutine GetStress_deviat_finite(rg,sg,rci,sci,bqhoi,sigma_devt, strnq_devtp, unkgd, ielem, istrn)
use constant
implicit none
real*8, intent(in)::rg,sg, rci, sci
real*8,dimension(3),                           intent(in)::bqhoi
real*8, dimension(1:2, 1:2),                  intent(out)::sigma_devt
real*8,dimension(1:ndegr,1:4),                 intent(in)::unkgd
real*8,dimension(1:2,1:2),                  intent(inout)::strnq_devtp
integer, intent(in)::ielem, istrn

!...Local array
real*8::jacbf(2, 2)
real*8::strn_totl(2, 2), strn_devt(2, 2), strn_devte(2, 2)
real*8::sigmatr_devt(2,2)
real*8::bi(ndegr)
real*8, dimension(2, 2):: gdelas_devt, Cstrnelas, Cstrnelas_inv, Cstrnelas_devt, Glagstrn_devt
real*8, dimension(2, 2)::stres_pk2,stres_pk2_devt

!...Local real
real*8::dr,ds
real*8::gmodu, stres_y, stres_eq, eqstres1
real*8::strn_vlum
real*8::scproj, djacoc ,djaco

!...Local integer
integer:: iv,ideg
!
!...Part I:: Prepare some parameters
!
dr = 1.d0
ds = 1.d0

!...The material shear modulus and yield strength
if(ncase.eq.1)then
gmodu   = 0.286d0
stres_y = 0.0026d0
!
elseif(ncase.eq.2)then
gmodu   = 1.5111d0
stres_y = 1.d0

elseif(ncase.eq.3)then
gmodu   = 0.276d0
stres_y = 0.003d0

elseif(ncase.eq.104)then
gmodu   = 1.5111d0
stres_y = 0.0033d0

endif
!
!...Part II:: Get the inifinitesimal strain
!

!...Jacobian transformation matrix
jacbf = 0.d0
!
bi(1) = 1.d0
bi(2) = (rg - rci)/dr
bi(3) = (sg - sci)/ds
!DGP2
if(npoly.eq.2)then
bi(4) = 0.5d0*bi(2)*bi(2)  - bqhoi(1)
bi(5) = 0.5d0*bi(3)*bi(3)  - bqhoi(2)
bi(6) =       bi(2)*bi(3)  - bqhoi(3)
endif

!
jacbf = 0.d0
do ideg = 1, ndegr
jacbf(1,1) = jacbf(1, 1) + unkgd(ideg, 1)*bi(ideg)
jacbf(1,2) = jacbf(1, 2) + unkgd(ideg, 2)*bi(ideg)
jacbf(2,1) = jacbf(2, 1) + unkgd(ideg, 3)*bi(ideg)
jacbf(2,2) = jacbf(2, 2) + unkgd(ideg, 4)*bi(ideg)
enddo
!
djaco = (jacbf(1,1)*jacbf(2,2)-jacbf(1,2)*jacbf(2,1))

!...Deviatoric elastic deformation gradient
gdelas_devt = 1.d0/(djaco)**(1.d0/3.d0)*jacbf

Cstrnelas_devt(1, 1) = gdelas_devt(1,1)*gdelas_devt(1,1) + gdelas_devt(2,1)*gdelas_devt(2,1)
Cstrnelas_devt(1, 2) = gdelas_devt(1,1)*gdelas_devt(1,2) + gdelas_devt(2,1)*gdelas_devt(2,2)

Cstrnelas_devt(2, 1) = gdelas_devt(1,2)*gdelas_devt(1,1) + gdelas_devt(2,2)*gdelas_devt(2,1)
Cstrnelas_devt(2, 2) = gdelas_devt(1,2)*gdelas_devt(1,2) + gdelas_devt(2,2)*gdelas_devt(2,2)

!...Cauchy strain
Cstrnelas(1, 1) = jacbf(1,1)*jacbf(1,1) + jacbf(2,1)*jacbf(2,1)
Cstrnelas(1, 2) = jacbf(1,1)*jacbf(1,2) + jacbf(2,1)*jacbf(2,2)

Cstrnelas(2, 1) = jacbf(1,2)*jacbf(1,1) + jacbf(2,2)*jacbf(2,1)
Cstrnelas(2, 2) = jacbf(1,2)*jacbf(1,2) + jacbf(2,2)*jacbf(2,2)

!
djacoc = (Cstrnelas(1,1)*Cstrnelas(2,2)-Cstrnelas(1,2)*Cstrnelas(2,1))

!...Inverse Cauchy strain
Cstrnelas_inv(1, 1) = Cstrnelas(2, 2)/djacoc
Cstrnelas_inv(1, 2) =-Cstrnelas(1, 2)/djacoc

Cstrnelas_inv(2, 1) =-Cstrnelas(2, 1)/djacoc
Cstrnelas_inv(2, 2) = Cstrnelas(1, 1)/djacoc
!
Glagstrn_devt(1, 1) = 0.5d0*(Cstrnelas_devt(1, 1) - 1.d0)
Glagstrn_devt(1, 2) = 0.5d0*(Cstrnelas_devt(1, 2))
Glagstrn_devt(2, 1) = 0.5d0*(Cstrnelas_devt(2, 1))
Glagstrn_devt(2, 2) = 0.5d0*(Cstrnelas_devt(2, 2) - 1.d0)

!
stres_pk2 = 2.d0*gmodu*Glagstrn_devt
!
scproj = stres_pk2(1, 1)*Cstrnelas(1, 1) +  stres_pk2(1, 2)*Cstrnelas(1, 2)+&
stres_pk2(2, 1)*Cstrnelas(2, 1) +  stres_pk2(2, 2)*Cstrnelas(2, 2)
!
stres_pk2_devt = stres_pk2 - 1.d0/3.d0*scproj*Cstrnelas_inv


!...Trial stress
sigmatr_devt(1, 1) = jacbf(1,1)*jacbf(1,1)*stres_pk2_devt(1, 1) + jacbf(1,1)*jacbf(1,2)*stres_pk2_devt(2, 1) +&
jacbf(1,1)*jacbf(1,2)*stres_pk2_devt(1, 2) + jacbf(1,2)*jacbf(1,2)*stres_pk2_devt(2, 2) !

sigmatr_devt(1, 2) = jacbf(1,1)*jacbf(2,1)*stres_pk2_devt(1, 1) + jacbf(1,2)*jacbf(2,1)*stres_pk2_devt(2, 1) +&
jacbf(1,1)*jacbf(2,2)*stres_pk2_devt(1, 2) + jacbf(1,2)*jacbf(2,2)*stres_pk2_devt(2, 2) !

sigmatr_devt(2, 1) = sigmatr_devt(1, 2)

sigmatr_devt(2, 2) = jacbf(2,1)*jacbf(2,1)*stres_pk2_devt(1, 1) + jacbf(2,1)*jacbf(2,2)*stres_pk2_devt(2, 1) +&
jacbf(2,1)*jacbf(2,2)*stres_pk2_devt(1, 2) + jacbf(2,2)*jacbf(2,2)*stres_pk2_devt(2, 2) !
!if(ielem.eq.1)then
sigmatr_devt  = sigmatr_devt/djaco

!...Branch
!if(stres_eq.le.stres_y)then

sigma_devt = sigmatr_devt
strnq_devtp = 0.d0

!
!if(ielem.eq.1)then
! print*,'stress1',stres_eq,stres_y,strn_vlum
! print*,'stress2',jacbf,rg,sg,rci,sci
! print*,'stress3',unkgd
! print*,'stress4',strn_totl
! print*,'stress5',strn_devt
! print*,'stress6',strnq_devtp
! print*,'stresse',strn_devte
! print*,'sigmatr',sigmatr_devt
! print*,'sigma',sigma_devt
!endif
!...Recover the solid dynamics
!sigma_devt = 0.d0

end subroutine GetStress_deviat_finite

!
!...Get the deviatoric stress for solid dynamics with infinitesimal strain
!
subroutine GetStress_deviat_infsmal(rg,sg,rci,sci,bqhoi,sigma2d_devt, strnq_devtp, unkgd, ielem, istrn)
use constant
implicit none
real*8, intent(in)::rg,sg, rci, sci
real*8,dimension(3),                           intent(in)::bqhoi
real*8, dimension(1:2, 1:2),                  intent(out)::sigma2d_devt
real*8,dimension(1:ndegr,1:4),                 intent(in)::unkgd
real*8,dimension(1:3,1:3),                  intent(inout)::strnq_devtp
integer, intent(in)::ielem, istrn

!...Local array
real*8::jacbf(3, 3)
real*8::strn_totl(3, 3), strn_devt(3, 3), strn_devte(3, 3)
real*8::sigmatr_devt(3,3),sigma_devt(3,3)
real*8::bi(ndegr)

!...Local real
real*8::dr,ds
real*8::gmodu, stres_y, stres_eq, eqstres1
real*8::strn_vlum

!...Local integer
integer:: iv,ideg
!
!...Part I:: Prepare some parameters
!
dr = 1.d0
ds = 1.d0

!...The material shear modulus and yield strength
if(ncase.eq.1)then
gmodu   = 0.286d0
stres_y = 0.0026d0
!
elseif(ncase.eq.2)then
gmodu   = 1.5111d0
stres_y = 1.d0

elseif(ncase.eq.3)then
gmodu   = 0.276d0
stres_y = 0.003d0

elseif(ncase.eq.104)then
gmodu   = 1.5111d0
stres_y = 0.0033d0

endif
!
!...Part II:: Get the inifinitesimal strain
!

!...Jacobian transformation matrix
jacbf = 0.d0
!
bi(1) = 1.d0
bi(2) = (rg - rci)/dr
bi(3) = (sg - sci)/ds
!DGP2
if(npoly.eq.2)then
bi(4) = 0.5d0*bi(2)*bi(2)  - bqhoi(1)
bi(5) = 0.5d0*bi(3)*bi(3)  - bqhoi(2)
bi(6) =       bi(2)*bi(3)  - bqhoi(3)
endif

!
jacbf = 0.d0
do ideg = 1, ndegr
jacbf(1,1) = jacbf(1, 1) + unkgd(ideg, 1)*bi(ideg)
jacbf(1,2) = jacbf(1, 2) + unkgd(ideg, 2)*bi(ideg)
jacbf(2,1) = jacbf(2, 1) + unkgd(ideg, 3)*bi(ideg)
jacbf(2,2) = jacbf(2, 2) + unkgd(ideg, 4)*bi(ideg)
enddo
!
jacbf(3,3) = 1.d0

!...
strn_totl =0.d0
!...inifinitesimal strain
strn_totl(1,1) = jacbf(1,1)-1.d0
strn_totl(1,2) = 0.5d0*(jacbf(1,2)+jacbf(2,1))
strn_totl(2,1) = strn_totl(1,2)
strn_totl(2,2) = jacbf(2,2)-1.d0
!
strn_totl(3,3) = 0.d0
!

!...Volume strain
strn_vlum = 1.d0/3.d0*(strn_totl(1,1) + strn_totl(2,2) + strn_totl(3,3))

!
strn_devt = 0.d0
!...Deviatoric strain
strn_devt(1,1) = strn_totl(1,1) - strn_vlum
strn_devt(1,2) = strn_totl(1,2)
strn_devt(2,1) = strn_totl(2,1)
strn_devt(2,2) = strn_totl(2,2) - strn_vlum

strn_devt(3,3) = - strn_vlum

!...Plastic Deviatoric strain for the last time step

!...Elatic Deviatoric strain
strn_devte = strn_devt - strnq_devtp

!...Trial stress
sigmatr_devt = 2.d0*gmodu*strn_devte
!
eqstres1 = sigmatr_devt(1,1)**2 + sigmatr_devt(1,2)**2 +  sigmatr_devt(1,3)**2 &
         + sigmatr_devt(2,1)**2 + sigmatr_devt(2,2)**2 +  sigmatr_devt(2,3)**2 &
         + sigmatr_devt(3,1)**2 + sigmatr_devt(3,2)**2 +  sigmatr_devt(3,3)**2

!
stres_eq = sqrt(1.5d0*eqstres1)
!
!if(ielem.eq.1)then
!print*,'strain',strnq_devtp
!endif

!...Branch
if(stres_eq.le.stres_y)then

sigma_devt = sigmatr_devt
!
strnq_devtp = strnq_devtp
else

! print*,'bad',ielem,stres_eq,stres_y
! stop
sigma_devt = sigmatr_devt*stres_y/stres_eq
if(istrn.eq.1)strnq_devtp = strn_devt -0.5d0/gmodu*sigma_devt
endif

!...Store the 2D stress
sigma2d_devt(1,1) = sigma_devt(1,1)
sigma2d_devt(1,2) = sigma_devt(1,2)
sigma2d_devt(2,1) = sigma_devt(2,1)
sigma2d_devt(2,2) = sigma_devt(2,2)

!
!if(ielem.eq.1)then
! print*,'stress1',stres_eq,stres_y,strn_vlum
! print*,'stress2',jacbf,rg,sg,rci,sci
! print*,'stress3',unkgd
! print*,'stress4',strn_totl
! print*,'stress5',strn_devt
! print*,'stress6',strnq_devtp
! print*,'stresse',strn_devte
! print*,'sigmatr',sigmatr_devt
! print*,'sigma',sigma_devt
!endif
!...Recover the solid dynamics
!sigma_devt = 0.d0

end subroutine GetStress_deviat_infsmal
!
!...Run the EOS equations for solid
!
subroutine GetEOS(nmatel,ncase,gamlg,rho, eintl, rhoi, pres, sd ,ielem)
!
implicit none
integer, intent(in):: nmatel, ncase,ielem
real*8,  intent(in):: gamlg,rho, eintl,rhoi
real*8, intent(out):: pres,sd


real*8::presh, enrgh, eta
real*8::dphdeta,dehdeta,detadrho
real*8::dphdrho, sd2
real*8::s0,c0,gmodu
real*8::eps
!
eps =1d-6
!
select case (nmatel)

!...Idea gas
case(1)
!...Pres
pres = max(eps, (gamlg-1.d0)*rho*eintl)
!...Sound speed
sd = sqrt(gamlg*pres/rho)

!...Solid
case(2)

 if(ncase.eq.1)then
  s0 = 1.34d0
  c0 = 0.533d0
  gmodu = 0.286d0

 elseif(ncase.eq.2)then
  s0 = 1.124d0
  c0 = 1.287d0
  gmodu = 1.5111d0

 elseif(ncase.eq.3)then
 s0 = 1.338d0
 c0 = 0.533d0
 gmodu = 0.276d0
 elseif(ncase.eq.104)then
 s0 = 1.124d0
 c0 = 1.287d0
 gmodu = 1.5111d0

 endif

!...
eta = 1.d0-rhoi/rho

!...Hugniout pressure and energy
presh = rhoi*c0**2*eta/(1.d0-s0*eta)**2
enrgh = eta*presh/2.d0/rhoi

!...Pressure
pres = presh + gamlg*rho*(eintl-enrgh)

!...Sound speed
detadrho = rhoi/rho**2
dphdeta  = rhoi*c0**2*(1.d0/(1.d0-s0*eta)**2 + 2.d0*eta*s0/(1.d0-s0*eta)**3)
dehdeta   = .5d0/rhoi*(presh+eta*dphdeta)
sd2      = dphdeta*detadrho + gamlg*(eintl-enrgh) - gamlg*dehdeta*detadrho + gamlg*pres/rho
sd2      = sd2 + 4.d0/3.d0*gmodu/rho

sd = sqrt(sd2)
!
if(ielem.eq.-1)then
 print*,'ieleminside',rhoi,rho,eta,c0,s0
 print*,'ielem2',presh,enrgh,eintl,pres
 print*,'ielem3',detadrho,dphdeta,dehdeta,sd2
 print*,'ielem4',dphdeta*detadrho, gamlg*(eintl-enrgh),gamlg*dehdeta*detadrho ,gamlg*pres/rho,&
4.d0/3.d0*gmodu/rho
endif
!
end select
!
end subroutine GetEOS

!
!...Get the deviatoric stress for solid dynamics with infinitesimal strain
!
subroutine GetStress_deviatfem_infsmal(rg,sg,xpq,xpqi,sigma2d_devt, strnq_devtp, ielem, istrn)
use constant
implicit none
real*8, intent(in)::rg,sg
real*8,dimension(1:ndimn, 1:nvqua),           intent(in) :: xpq, xpqi
real*8, dimension(1:2, 1:2),                  intent(out)::sigma2d_devt
real*8,dimension(1:3,1:3),                  intent(inout)::strnq_devtp
integer, intent(in)::ielem, istrn

!...Local array
real*8::jacbf(3, 3)
real*8::strn_totl(3, 3), strn_devt(3, 3), strn_devte(3, 3)
real*8::sigmatr_devt(3,3),sigma_devt(3,3)
real*8::bi(ndegr)
real*8::shpq(nvqua),dspsq(nvqua),dsprq(nvqua)

!...Local real
real*8::dxdr,dxds,dydr,dyds,dxdri,dxdsi,dydri,dydsi
real*8::dr,ds,jacbi,rm,rp,sm,sp,r,s
real*8::gmodu, stres_y, stres_eq, eqstres1
real*8::strn_vlum
real*8::c10

!...Local integer
integer:: iv,ideg,ishp
!
c10 = 1.d0
!
!...Part I:: Prepare some parameters
!
dr = 1.d0
ds = 1.d0
!
r = rg
s = sg

!...The material shear modulus and yield strength
if(ncase.eq.1)then
gmodu   = 0.286d0
stres_y = 0.0026d0
!
elseif(ncase.eq.2)then
gmodu   = 1.5111d0
stres_y = 1.d0

elseif(ncase.eq.3)then
gmodu   = 0.276d0
stres_y = 0.003d0

endif
!
!...Part II:: Get the inifinitesimal strain
!
if(ncurv.eq.0)then
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

elseif(ncurv.eq.1)then
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
jacbi = dxdri*dydsi-dxdsi*dydri
!
jacbf(1,1)= dxdr*dydsi - dxds*dydri
jacbf(1,2)=-dxdr*dxdsi + dxds*dxdri
jacbf(2,1)= dydr*dydsi - dyds*dydri
jacbf(2,2)=-dydr*dxdsi + dyds*dxdri
!
jacbf = jacbf/jacbi

!...inifinitesimal strain
strn_totl(1,1) = jacbf(1,1)-1.d0
strn_totl(1,2) = 0.5d0*(jacbf(1,2)+jacbf(2,1))
strn_totl(2,1) = strn_totl(1,2)
strn_totl(2,2) = jacbf(2,2)-1.d0

!...Volume strain
strn_vlum = 1.d0/3.d0*(strn_totl(1,1) + strn_totl(2,2))

!...Deviatoric strain
strn_devt = 0.d0

strn_devt(1,1) = strn_totl(1,1) - strn_vlum
strn_devt(1,2) = strn_totl(1,2)
strn_devt(2,1) = strn_totl(2,1)
strn_devt(2,2) = strn_totl(2,2) - strn_vlum

strn_devt(3,3) = - strn_vlum

!...Plastic Deviatoric strain for the last time step

!...Elatic Deviatoric strain
strn_devte = strn_devt - strnq_devtp

!...Trial stress
sigmatr_devt = 2.d0*gmodu*strn_devte
eqstres1 = sigmatr_devt(1,1)**2 + sigmatr_devt(1,2)**2 +  sigmatr_devt(1,3)**2 &
         + sigmatr_devt(2,1)**2 + sigmatr_devt(2,2)**2 +  sigmatr_devt(2,3)**2 &
         + sigmatr_devt(3,1)**2 + sigmatr_devt(3,2)**2 +  sigmatr_devt(3,3)**2
!
stres_eq = sqrt(1.5d0*eqstres1)
!
!if(ielem.eq.1)then
!print*,'strain',strnq_devtp
!endif

!...Branch
if(stres_eq.le.stres_y)then

sigma_devt = sigmatr_devt
strnq_devtp = strnq_devtp
else
sigma_devt = sigmatr_devt*stres_y/stres_eq

!...Guarantee the deviatoric stress is calculated only once
if(istrn.eq.1) strnq_devtp = strn_devt -0.5d0/gmodu*sigma_devt
endif

!...Store the 2D stress
sigma2d_devt(1,1) = sigma_devt(1,1)
sigma2d_devt(1,2) = sigma_devt(1,2)
sigma2d_devt(2,1) = sigma_devt(2,1)
sigma2d_devt(2,2) = sigma_devt(2,2)
!
!if(ielem.eq.-1)then
!print*,'stress1',stres_eq,stres_y,strn_vlum
!print*,'stress4',strn_totl
!print*,'stress5',strn_devt
!print*,'stress6',strnq_devtp
!print*,'stresse',strn_devte
!print*,'sigmatr',sigmatr_devt
!print*,'sigma',sigma_devt
!endif
!...Recover the solid dynamics
!sigma_devt = 0.d0

end subroutine GetStress_deviatfem_infsmal
!
!...Get the threshold  for solid dynamics with infinitesimal strain
!
subroutine GetPlast_threshold_infsmal(rg,sg,rci,sci, strnq_devtp, unkgd, threshld_plast, ielem)
use constant
implicit none
real*8, intent(in)::rg,sg, rci, sci
real*8,dimension(1:ndegr,1:4),    intent(in)::unkgd
real*8,dimension(1:3,1:3),        intent(in)::strnq_devtp
real*8,                          intent(out)::threshld_plast
integer, intent(in)::ielem

!...Local array
real*8::jacbf(3, 3)
real*8::strn_totl(3, 3), strn_devt(3, 3), strn_devte(3, 3)
real*8::sigmatr_devt(3,3)
real*8::bi(ndegr)

!...Local real
real*8::dr,ds
real*8::gmodu, stres_y, stres_eq, eqstres1
real*8::strn_vlum

!...Local integer
integer:: iv,ideg
!
!...Part I:: Prepare some parameters
!
dr = 1.d0
ds = 1.d0

!...The material shear modulus and yield strength
if(ncase.eq.1)then
gmodu   = 0.286d0
stres_y = 0.0026d0
!
elseif(ncase.eq.2)then
gmodu   = 1.5111d0
stres_y = 1.d0

elseif(ncase.eq.3)then
gmodu   = 0.276d0
stres_y = 0.003d0
elseif(ncase.eq.104)then
gmodu   = 1.5111d0
stres_y = 0.0033d0
endif
!
!...Part II:: Get the inifinitesimal strain
!

!...Jacobian transformation matrix
jacbf = 0.d0
!
bi(1) = 1.d0
bi(2) = (rg - rci)/dr
bi(3) = (sg - sci)/ds
!
jacbf = 0.d0
do ideg = 1, 3!ndegr
jacbf(1,1) = jacbf(1, 1) + unkgd(ideg, 1)*bi(ideg)
jacbf(1,2) = jacbf(1, 2) + unkgd(ideg, 2)*bi(ideg)
jacbf(2,1) = jacbf(2, 1) + unkgd(ideg, 3)*bi(ideg)
jacbf(2,2) = jacbf(2, 2) + unkgd(ideg, 4)*bi(ideg)
enddo
!
jacbf(3,3) = 1.d0

!...
strn_totl =0.d0
!...inifinitesimal strain
strn_totl(1,1) = jacbf(1,1)-1.d0
strn_totl(1,2) = 0.5d0*(jacbf(1,2)+jacbf(2,1))
strn_totl(2,1) = strn_totl(1,2)
strn_totl(2,2) = jacbf(2,2)-1.d0
!
strn_totl(3,3) = 0.d0

!...Volume strain
strn_vlum = 1.d0/3.d0*(strn_totl(1,1) + strn_totl(2,2) + strn_totl(3,3))

!
strn_devt = 0.d0
!...Deviatoric strain
strn_devt(1,1) = strn_totl(1,1) - strn_vlum
strn_devt(1,2) = strn_totl(1,2)
strn_devt(2,1) = strn_totl(2,1)
strn_devt(2,2) = strn_totl(2,2) - strn_vlum

strn_devt(3,3) = - strn_vlum

!...Plastic Deviatoric strain for the last time step

!...Elatic Deviatoric strain
strn_devte = strn_devt - strnq_devtp

!...Trial stress
sigmatr_devt = 2.d0*gmodu*strn_devte
!
eqstres1 = sigmatr_devt(1,1)**2 + sigmatr_devt(1,2)**2 +  sigmatr_devt(1,3)**2 &
+ sigmatr_devt(2,1)**2 + sigmatr_devt(2,2)**2 +  sigmatr_devt(2,3)**2 &
+ sigmatr_devt(3,1)**2 + sigmatr_devt(3,2)**2 +  sigmatr_devt(3,3)**2
!
stres_eq = sqrt(1.5d0*eqstres1)
!
!if(ielem.eq.1)then
!print*,'strain',strnq_devtp
!endif

!...Plasticity Threshold

threshld_plast = min(stres_eq/stres_y, 1.d0)

end subroutine GetPlast_threshold_infsmal
!
!...Lagrangian cell vaerage value on hybrid meshes for vtk legacy
!
subroutine getoutputsolid_vtk_cellhybrid(coord,coold,geoel,unkno,unkgd,strnq_devtp,iptri, ipqua)
use constant
implicit none
integer::cntr,cntr2,dump,ie,ip
integer::mtype(2)
integer,intent(in)::iptri(nvtri,ntria), ipqua(nvqua, nquad)
real*8, dimension(ndimn,npoin),intent(in)::coord,coold
real*8, intent(in)::geoel(1:ngeel,1:nsize)
real*8,dimension(1:mdegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ndegr,1:4,1:nsize), intent(in)::unkgd
real*8,dimension(1:3,1:3,ngstrnf+ngausdq, nquad), intent(in)::strnq_devtp
!
integer,dimension(1:nvqua) :: ipq
real*8,dimension(1:ndimn, 1:nvqua) :: xpq, xpqi
real*8,dimension(1:ndimn, 1:nvqua+ngausdq) :: rpq
real*8::xgq(1:2,ngausdq_pst,nquad)

real*8,dimension(1:ndegr):: bi
real*8, dimension(1:nvqua):: shpq, dsprq, dspsq
real*8::jacbf(2, 2)
!
real*8::weighq(ngausdq), posiq(2,ngausdq)
real*8::threshld_plast_g(ngstrnf+ngausdq)
real*8::threshld_plast_g2(ngstrnf+ngausdq_pst, nquad)

!
integer::ielem,ig, ideg, ishp
real*8:: rho, uadv, vadv, eadv
real*8:: xcrho,ycrho,rc,sc
real*8:: rhoct, eintc, rhoi, pctr, sdctr, uctr, vctr, ectr
real*8:: r ,s, rci,sci
real*8:: dxdr,dxds,dydr,dyds
real*8:: a11, a12, a21, a22
real*8:: dr, ds, djaco,wi
real*8:: threshld_plast, volel
real*8:: c10
!
c10 =1.d0
dr = 1.d0
ds = 1.d0
!
!...Give gaussian position and weight...
!
call ruqope(2, ngausdq, posiq, weighq)
!
if(ncurv.eq.0)then
rpq(1, 1) = -1.d0; rpq(2, 1) = -1.d0
rpq(1, 2) =  1.d0; rpq(2, 2) = -1.d0
rpq(1, 3) =  1.d0; rpq(2, 3) =  1.d0
rpq(1, 4) = -1.d0; rpq(2, 4) =  1.d0
!
rpq(:, (nvqua+1):(nvqua+ngausdq)) = posiq(:, 1:ngausdq)
elseif(ncurv.eq.1)then

endif

open(22,file='outputsolid_lag.vtk',status='unknown')
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
do ie=1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!...mass center
rc= geoel(1, ielem)
sc= geoel(2, ielem)
!
xpqi(1, 1:nvqua) = coold(1, ipq(1:nvqua))
xpqi(2, 1:nvqua) = coold(2, ipq(1:nvqua))
!
rhoct = 1.d0/unkno(1, 1, ielem)
uctr  = unkno(1, 2, ielem)
vctr  = unkno(1, 3, ielem)
ectr  = unkno(1, 4, ielem)
!
eintc = ectr-0.5d0*(uctr**2 + vctr**2)
!
call GetCellctr_quad_initial (ncurv,ndimn,nvqua,xpqi, rc, sc, xcrho, ycrho)
call getrhog_initial(rhoi,  xcrho, ycrho, xcrho, ycrho)
call GetEOS(nmatel,ncase,gamlg,rhoct, eintc, rhoi, pctr, sdctr, ielem)
!
write(22,'(3e20.8)') pctr
!
end do
write(22,*)
write(22,*)'SCALARS Plasticitythreshold float'
write(22,*)'LOOKUP_TABLE default'
!-xxx-Plasticity threshold
do ie=1,nquad
!
ipq(1:nvqua) = ipqua(1:nvqua,ie)
ielem = ie + ntria
!...mass center
rci= geoel(7, ielem)
sci= geoel(8, ielem)
!
!
!...Points consitituting one element...
!
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
volel = 0.d0
threshld_plast = 0.d0
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
djaco = wi
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
!
call GetPlast_threshold_infsmal(r,s,rci,sci, strnq_devtp(:,:,ngstrnf+ig,ielem), &
unkgd(:,:,ielem), threshld_plast_g(ngstrnf+ig), ielem)
!
threshld_plast = threshld_plast + threshld_plast_g(ngstrnf+ig)*djaco*(a11*a22-a12*a21)
volel = volel + djaco*(a11*a22-a12*a21)
!
enddo
!if(ielem.eq.5) print*,'test',ielem,unkgd(:,:,ielem),strnq_devtp(:,:,:,ielem)
!
write(22,'(3e20.8)') threshld_plast/volel
!
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

write(*,'(1A17,1A35)') 'Data written for vtk file!'

!
!...Part III: Output the stress ratio at the particle positions and variables
!
xgq = 0.d0
!
do ie = 1, nquad
!
ielem = ie + ntria
!
rc= geoel(1, ielem) !...mass center...
sc= geoel(2, ielem)

!...cell center
rci= geoel(7, ielem)
sci= geoel(8, ielem)
!
xpq(1, 1:nvqua) = coord(1,ipqua(1:nvqua, ie))
xpq(2, 1:nvqua) = coord(2,ipqua(1:nvqua, ie))

xpqi(1, 1:nvqua) = coold(1,ipqua(1:nvqua, ie))
xpqi(2, 1:nvqua) = coold(2,ipqua(1:nvqua, ie))

!...Gauss loop
do ig = 1,ngausdq_pst!...(2)ig = 1,ngausd

!...Reference coordinate
r  = posiq(1,ig)
s  = posiq(2,ig)

!... Shape function & its derivative
call getshapfct_quad(ncurv,nvqua,shpq, dsprq, dspsq, r, s)

!...Get the points
do ishp = 1, nvqua
xgq(1, ig, ie) = xgq(1, ig, ie) + shpq(ishp)*xpq(1,ishp)
xgq(2, ig, ie) = xgq(2, ig, ie) + shpq(ishp)*xpq(2,ishp)
enddo

!
call GetPlast_threshold_infsmal(r,s,rci,sci, strnq_devtp(:,:,ngstrnf+ig,ielem), &
unkgd(:,:,ielem), threshld_plast_g2(ngstrnf+ig, ielem), ielem)!

enddo
!
enddo

!...Part III.2: output for vtk
open(22,file='outputsolid-particle.vtk',status='unknown')
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
write(22,*)'SCALARS Yield ratio float'
write(22,*)'LOOKUP_TABLE default'

!-xxx-Density
if(ndens.eq.1)then
do ie=1,nquad
ielem = ie +ntria
do ig = 1,ngausdq_pst
write(22,*)threshld_plast_g2(ngstrnf+ig, ie)
enddo
enddo
endif

!
write(*,'(1A17,1A35)') 'Particle data written for MEM solid!'
end subroutine getoutputsolid_vtk_cellhybrid
!
!...Get different energies  for solid dynamics
!
subroutine GetEnergy_solid(unkno,geoel,enrgy_solid)
use constant
implicit none
real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
real*8,dimension(1:ngeel,1:nsize),     intent(in)::geoel
real*8,dimension(1:3),                intent(out)::enrgy_solid

!...Local integer
integer:: ie,ielem
real*8::uctr,vctr,ectr,eintc
!
enrgy_solid = 0.d0

do ie =1, nquad
ielem = ie + ntria
!
uctr = unkno(1, 2, ielem)
vctr = unkno(1, 3, ielem)
ectr = unkno(1, 4, ielem)
!
eintc  = ectr - 0.5d0*(uctr**2 + vctr**2)
!
enrgy_solid(1) = enrgy_solid(1) + ectr*geoel(4, ielem)
enrgy_solid(2) = enrgy_solid(2) + 0.5d0*(uctr**2 + vctr**2)*geoel(4, ielem)
enddo
!
enrgy_solid(3) = enrgy_solid(1) - enrgy_solid(2)


end subroutine GetEnergy_solid
