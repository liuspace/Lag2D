!
!...Get the shape functions for an curved edge
!
subroutine getshapfct_edge(ncurv,nvfac,shpf, dshprf, r)
implicit none
integer,                   intent(in)::ncurv,nvfac
real*8, dimension(nvfac), intent(out):: shpf
real*8, dimension(nvfac), intent(out):: dshprf
real*8,                    intent(in):: r
!
if(ncurv.eq.0)then
!...Shape functions
shpf(1) =  (1.d0-r)/2.d0
shpf(2) =  (1.d0+r)/2.d0

dshprf(1) =-.5d0
dshprf(2) = .5d0

!...Quadratic cell
elseif(ncurv.eq.1)then
shpf(1) =  -0.5d0*(1.d0-r)*r
shpf(2) =   0.5d0*(1.d0+r)*r
shpf(3) =         (1.d0+r)*(1.d0-r)
!
dshprf(1) = -0.5d0 + r
dshprf(2) =  0.5d0 + r
dshprf(3) = -2.d0*r

!...Cubic cell
elseif(ncurv.eq.2)then
!...Shape functions
shpf(1) =  9.d0/16.d0*(1.d0-r)*(r+1.d0/3.d0)*(r-1.d0/3.d0)
shpf(2) = -9.d0/16.d0*(1.d0+r)*(1.d0/3.d0-r)*(r+1.d0/3.d0)
shpf(3) =  27.d0/16.d0*(r-1.d0)*(r+1.d0)*(r-1.d0/3.d0)
shpf(4) = -27.d0/16.d0*(r-1.d0)*(r+1.d0)*(r+1.d0/3.d0)
!
dshprf(1) = 9.d0/16.d0*(1.d0/9.d0-3.d0*r**2+2.d0*r)
dshprf(2) =-9.d0/16.d0*(1.d0/9.d0-3.d0*r**2-2.d0*r)
dshprf(3) = 27.d0/16.d0*(-1.d0+3.d0*r**2-2.d0/3.d0*r)
dshprf(4) =-27.d0/16.d0*(-1.d0+3.d0*r**2+2.d0/3.d0*r)
!
endif

end subroutine getshapfct_edge
!
!...Get the shape functions for an curved cell
!
subroutine getshapfct_quad(ncurv,nvqua,shpq, dshprq, dshpsq, r, s)
implicit none
integer,                   intent(in)::ncurv,nvqua
real*8, dimension(nvqua), intent(out):: shpq
real*8, dimension(nvqua), intent(out):: dshprq, dshpsq
real*8,                    intent(in):: r, s
!
!...Local array
real*8,dimension(2, nvqua) ::xrq
integer,dimension(1:nvqua) ::idmpp
!
real*8::rsp,ssp
real*8::rp,rm,sp,sm
integer::ishp

real*8::c10
!
c10 =1.d0

!...Qudratic cell
if(ncurv.eq.0)then
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
!
dshprq(1) = -0.25d0   *sm
dshprq(2) =  0.25d0   *sm
dshprq(3) =  0.25d0   *sp
dshprq(4) = -0.25d0   *sp
!
dshpsq(1) =-0.25d0*rm
dshpsq(2) =-0.25d0*rp
dshpsq(3) = 0.25d0*rp
dshpsq(4) = 0.25d0*rm
elseif(ncurv.eq.1)then
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
dshprq(1) = 0.25d0*(2.d0*r-c10)*s*(s-1.d0)
dshprq(2) = 0.25d0*(2.d0*r+c10)*s*(s-1.d0)
dshprq(3) = 0.25d0*(2.d0*r+c10)*s*(s+1.d0)
dshprq(4) = 0.25d0*(2.d0*r-c10)*s*(s+1.d0)
dshprq(5) =   -r*s*(s-c10)
dshprq(6) =  0.5d0*(2.d0*r+c10)*(c10-s*s)
dshprq(7) =   -r*s*(s+c10)
dshprq(8) = 0.5d0*(2.d0*r-c10)*(1.d0-s*s)
dshprq(9) = -2.0d0*r*(1.d0-s*s)
!
dshpsq(1) = 0.25d0*(2.d0*s-c10)*r*(r-1.d0)
dshpsq(2) = 0.25d0*(2.d0*s-c10)*r*(r+1.d0)
dshpsq(3) = 0.25d0*(2.d0*s+c10)*r*(r+1.d0)
dshpsq(4) = 0.25d0*(2.d0*s+c10)*r*(r-1.d0)
dshpsq(5) =  0.5d0*(1.d0-r*r)*(2.d0*s-c10)
dshpsq(6) =   -r*s*(r+c10)
dshpsq(7) =  0.5d0*(1.d0-r*r)*(2.d0*s+c10)
dshpsq(8) = -r*s*(r-c10)
dshpsq(9) = -2.0d0*s*(1.d0-r*r)

!...Cubic cell
elseif(ncurv.eq.2)then
!...Mapping array
idmpp(1:4) =  (/1,2,3,4/)
idmpp(5:8) =  (/5,9,7,11/)
idmpp(9:12) = (/6,10,8,12/)

!...Reference coordinates
xrq(1, 1) = -1.d0; xrq(2, 1) = -1.d0;
xrq(1, 2) =  1.d0; xrq(2, 2) = -1.d0;
xrq(1, 3) =  1.d0; xrq(2, 3) =  1.d0;
xrq(1, 4) = -1.d0; xrq(2, 4) =  1.d0;
xrq(1, 5) = -1.d0/3.d0; xrq(2, 5) = -1.d0;
xrq(1, 6) =  1.d0;      xrq(2, 6) = -1.d0/3.d0;
xrq(1, 7) =  1.d0/3.d0; xrq(2, 7) =  1.d0;
xrq(1, 8) = -1.d0;      xrq(2, 8) =  1.d0/3.d0;
xrq(1, 9) =  1.d0/3.d0; xrq(2, 9) = -1.d0;
xrq(1,10) =  1.d0;      xrq(2,10) =  1.d0/3.d0;
xrq(1,11) = -1.d0/3.d0; xrq(2,11) =  1.d0;
xrq(1,12) = -1.d0;      xrq(2,12) = -1.d0/3.d0;

!...Shape functions

!...Corner nodes
do ishp =1,4
rsp = xrq(1, idmpp(ishp))
ssp = xrq(2, idmpp(ishp))
shpq(idmpp(ishp)) = 1.d0/32.d0*(1.d0+r*rsp)*(1.d0+s*ssp)*&
(9.d0*(r**2+s**2)-10.d0)
enddo

!...Edge nodes
do ishp =5, 8
rsp = xrq(1, idmpp(ishp))
ssp = xrq(2, idmpp(ishp))
shpq(idmpp(ishp)) = 9.d0/32.d0*(1.d0+s*ssp)*(1.d0-r**2)*&
(9.d0*r*rsp+1.d0)
enddo

!...Edge nodes
do ishp =9 ,12
rsp = xrq(1, idmpp(ishp))
ssp = xrq(2, idmpp(ishp))
shpq(idmpp(ishp)) = 9.d0/32.d0*(1.d0+r*rsp)*(1.d0-s**2)*&
(9.d0*s*ssp+1.d0)
enddo

!...Derivative for shape functions

!...Corner nodes
do ishp =1,4
rsp = xrq(1, idmpp(ishp))
ssp = xrq(2, idmpp(ishp))
dshprq(idmpp(ishp)) = 1.d0/32.d0*(1.d0+s*ssp)*&
(27.d0*rsp*r**2+9.d0*rsp*s**2-10.d0*rsp+18.d0*r)
dshpsq(idmpp(ishp)) = 1.d0/32.d0*(1.d0+r*rsp)*&
(27.d0*ssp*s**2+9.d0*ssp*r**2-10.d0*ssp+18.d0*s)
enddo

!...Edge nodes
do ishp =5,8
rsp = xrq(1, idmpp(ishp))
ssp = xrq(2, idmpp(ishp))
dshprq(idmpp(ishp)) = 9.d0/32.d0*(1.d0+s*ssp)*&
(-27.d0*r**2*rsp-2.d0*r+9.d0*rsp)
dshpsq(idmpp(ishp)) = 9.d0/32.d0*ssp*(1.d0-r**2)*&
(9.d0*r*rsp+1.d0)

enddo

!...Edge nodes
do ishp =9,12
rsp = xrq(1, idmpp(ishp))
ssp = xrq(2, idmpp(ishp))
dshprq(idmpp(ishp)) = 9.d0/32.d0*rsp*(1.d0-s**2)*&
(9.d0*s*ssp+1.d0)
dshpsq(idmpp(ishp)) = 9.d0/32.d0*(1.d0+r*rsp)*&
(-27.d0*s**2*ssp-2.d0*s+9.d0*ssp)
enddo

!
endif

end subroutine getshapfct_quad


