!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!...2D-DG by Xiaodong Liu
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 module constant 
  implicit none
  public
  integer::nsize, ncell
  integer*4::nelem,npoin,nbfac,ndimn,ntype,nafac,ntime,ncase,nout,ndump
  integer*4::nrz !...ncart:indicator of cartesian or RZ coordinate
  integer:: nmatel !...the inicator for materials
  integer:: npoin1, npoin2, nofile
  integer*4::ntria,nquad
  integer*4::nflux,ndegr,npoly,nreco,nlimi,ngaus,nmatr,mdegr 
  integer*4::ngausf,ngausd,nint,nmeth,nq,ndt, ndens
  integer*4::ngausf_geo,ngausd_geo
  integer*4::ncurv, nptfa, nptri
  integer*4::nresi
  integer::nfint, nriem
  integer*4::nvfac, nvtri !...No. of vertex for one face and one triangle.
  integer*4::nbfai
  integer*4,save::ngausd_geoq = 36 !...geometry quadrature for quads...
  integer*4,save::ngausdq = 25!...Domain integral for quads...
  integer*4,save::ngausdqi = 25   !...Domain integral for quads...
  integer*4,save::ngausdq_pst =25  !...Domain integral for quads...
  integer*4,save::ngeel = 24
  integer*4,save::ngefa = 7   
  integer*4,save::nifai = 8
  integer*4,save::nmati = 6
  integer*4,save::nstag
  integer*4,save::nepfa, nftri
  integer, save::nfqua, npqua, nvqua
  integer, save::nsubq,nsubt
  integer::rkstg
  integer::ngsft, ngsfq,ngvtxf
  real*8,save::gamma=1.4d0,toler=-15.d0,pi=4.d0*atan(1.d0)
!...Parametrs for deformation gradient
  integer*4,save::ngedg = 3
!...Lagrangian
  integer, save:: ngstrnf
  integer*4,save::ngelg = 12, ngelgq=16
  integer*4,save::ngesgt = 15, ngesgq=16
  integer*4,save::ngflg = 3
  real*8,dimension(1:3)::rkcoe=[0.d0,0.75d0,0.33333333333333d0]
  real*8,dimension(2,3)::rklag
  real*8,dimension(2)::crklg
  real*8,dimension(3)::crklg2
!
  real*8, dimension(8, 25)::weigh0t
  real*8, dimension(6, 8)::weigh1t
  real*8, dimension(4, 6)::weigh2t
  real*8, dimension(1, 4)::weigh3t
  real*8::amach,cfl,beta,breta,dtfix,gamlg,cdrho,tend,cimpd
  real*8::alfrz, epsaw
  real*8::tfcus, tacru
  real*8::slpdu
  real*8::tstart_cpu, tend_cpu
 end module constant
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!...Main program...
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
program main
use constant
implicit none
 integer*4,allocatable::inpoel(:,:),intfac(:,:)
 integer,  allocatable::ipqua(:,:), iptri(:,:)
 real*8,allocatable::coord(:,:),cooro(:,:),coora(:,:),geofa(:,:),geoel(:,:),resnorm(:,:),resinit(:,:),coold(:,:)
 real*8,allocatable::amatr(:,:),amatr2(:,:),ustar(:,:) !mass matrix
 integer*4,allocatable::bface(:,:)
 real*8,allocatable::unkno(:,:,:),uchar(:),rhsel(:,:,:),unold(:,:,:),reunk(:,:,:),unkaw(:,:,:)
 real*8,allocatable::unknp(:,:)
 real*8,allocatable::strnq_devtp(:,:,:,:)
!...The DoFs and geometry parameters for the deformation gradient
 real*8,allocatable::unkgd(:,:,:), unogd(:,:,:), gesgq0(:, :, :), gesgt0(:, :, :), amagd(:,:)
 real*8,allocatable::ujacb(:)
 real*8,allocatable::gflag(:, :)
 integer*4::itime,ie,iq,istag,i,id,ideg
 real*8::re,alfa,xg,yg,dx,dy,b(6),p,q,x1,x2
 real*8,allocatable::pres(:),vxpoin(:),vypoin(:),velo(:),resmom(:,:)
 real*8,allocatable::enrgy_solid(:)
 integer*4,allocatable::ltelem(:,:), fsuel(:,:)
 integer*4,allocatable::esqua(:,:), estri(:, :)
integer, allocatable:: esuv1(:), esuv2(:)
integer, allocatable:: fasup(:,:)
!
!...Step 0 :Machine learning weight
!
!...Output training data csv fromat.
!call getWeight_ML
!
!xxxxxxx Step 1, preprocessing xxxxxxxxx
!
!...1.1: Preprocess to get some information, i.e., nafac, nelem, nbfac...
  call preprocess1
!
!...The No. of the deviatoric plastic stress in one clel
  if(nmeth.eq.2)then
   if(ncurv.eq.0)then
    ngstrnf = 4
   elseif(ncurv.eq.1)then
    ngesgq = 16
    nsubq =4
    nsubt =4 !...This will corrected in future!
    ngstrnf = 9
    if(nfint.eq.5.or.nfint.eq.10) ngstrnf = 17
   elseif(ncurv.eq.2)then
   ngesgq = 24
    nsubq =8
    nsubt =4 !...This will corrected in future!
   endif
  endif

!...Allocate integer array...
    allocate(intfac(1:nifai,1:nafac))
    allocate(inpoel(1:nvtri,1:ntria), iptri(1:nvtri,1:ntria), ipqua(1:nvqua, 1:nquad))
    allocate(bface(1:nbfai,nbfac))

!...Allocate real array...
    allocate(coord(1:ndimn,1:npoin), coold(1:ndimn,1:npoin), coora(1:ndimn,1:npoin))
    allocate(geoel(1:ngeel,1:nsize))
    allocate(geofa(1:ngefa,1:nafac))
    allocate(unkno(1:ndegr,1:nq,1:nsize))
    allocate(unold(1:ndegr,1:nq,1:nsize))
    allocate(unkaw(1:ndegr,1:nq,1:nsize))
    allocate(unknp(1:npoin,1:nq))
    allocate(strnq_devtp(1:3,1:3,ngstrnf+ngausdq, nquad))
    allocate(rhsel(1:ndegr,1:nq,1:ncell))
    allocate(resnorm(1:ndegr,1:nq))
    allocate(resinit(1:ndegr,1:nq))
    allocate(uchar(1:nq))
    allocate(enrgy_solid(3))
    allocate(ltelem(1:3,1:ncell) ,fsuel(1:3, 1:ncell))
    allocate(estri(1:nftri,1:ntria), esqua(1:nfqua, 1:nquad))
    allocate(esuv1(npoin1),esuv2(npoin2))
!...Paramertes for deformation gradient
    allocate(unkgd(1:ndegr,1:4,1:nsize),unogd(1:ndegr,1:4,1:nsize))
    allocate(gesgq0(1:3,1:ngesgq,1:nquad))
    allocate(gesgt0(1:3,1:ngesgt,1:ntria+nbfac))
    allocate(gflag(1:ngflg,1:nbfac))
    allocate(ujacb(1:ncell))
!
    if(nmeth.eq.2)then
      allocate(cooro(1:ndimn,1:npoin))
      allocate(ustar(1:ndimn,1:npoin))
      allocate(fasup(1:4,1:npoin))

!...Coefficient for acoustic impedance
      slpdu = (gamlg+1.d0)/2.d0
!      print*,'slpdu',gamlg
!...Specify cylindrical coordinates or not
      if(nrz.eq.1.or.nrz.eq.2)then !...RZ...
        alfrz = 1.d0
        tfcus = 6.726842539779d-3
      elseif(nrz.eq.0)then !...Cartesian coordinate...
        alfrz = 0.d0     !...alfrz = 0 should recover the original Cartesian....
        tfcus = 7.26483157d-3
      endif
!...The No. of face Gauss points for one triangle or quad...
      ngsft = ngausf*3
      ngsfq = ngausf*4

    endif
!...Zero out the initial time
    tacru = 0.d0

!...Zero out nofile...
     nofile = -1

!...1.2: Preprocess to get required arrays, i.e., coord, bface, geofa...
 call preprocess2(intfac,inpoel,iptri, ipqua, geofa,geoel,bface,coord, esuv1, esuv2)

 if(ncurv.le.1)then
  call getgeofa(intfac, geofa, coord)
  if(nmeth.eq.1)then
    call getgeoel(inpoel,iptri, ipqua, geoel, coord)
  elseif(nmeth.eq.2)then
 !   call  getfasup(intfac, fasup)
 !...For lagrangian, geoel is used to store the mass cenetr for reference cell...
    call getgeoel_initial_lag(iptri, ipqua, geoel, coord)
 !   call getgeoel_lag_mc(inpoel,iptri,ipqua, geoel, coord)
 !...DGP2
    if(npoly==2)then
     call getgeoelho_initial_lag(iptri, ipqua, geoel, coord)
 !    call getgeoel_lagmc_ho(inpoel, iptri, ipqua, geoel, coord)
 !...High-order terms for the deformation gradient
 !    if(nfint.eq.8.or.nfint.eq.9.or.nfint.eq.10) call getgeoel_lag_hogd(iptri, ipqua, geoel, coord)
    endif
  endif
 elseif(ncurv.ge.2)then
  if(nmeth.eq.1)then
   call getgeofa(intfac, geofa, coord)
   call getgeoel_curv2(inpoel, geoel, coord)
  elseif(nmeth.eq.2)then
   !...For lagrangian, geoel is used to store the mass cenetr for reference cell...
   call getgeoel_initial_lag_curv(iptri, ipqua, geoel, coord)
   if(npoly==2)then
    call getgeoelho_initial_lag_curv(iptri, ipqua, geoel, coord)
   endif
  endif
 endif

!...1.3: Check mesh quality...
!call getmeshquality(inpoel,coord,geoel)
!
!...1.4: Calculate the mass matrix for DGMs, npoly=1 for DGP1, and npoly=2 for DGP2...
 if(nmeth.eq.1)then   !...Eulerian framework...
  if(npoly==1)then
    nmatr = 3            !the dimension of mass matrix,p1
  elseif(npoly==2)then
    nmatr = 15           !the dimension of mass matrix,p2
  endif
  allocate(amatr(nmatr, ncell))
!  call getamatr(amatr,geoel,coord,inpoel)
  if(ncurv.le.1)then
   call getamatr_curv(amatr,geoel,coord,inpoel,iptri,ipqua)
  elseif(ncurv.ge.2)then
   call getamatr_curv2(amatr,geoel,coord,inpoel)
  endif
 endif

! call getlogeo(intfac,geoel,geofa,ltelem,fsuel)
 print*,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'Successfully finished preprocessing'
 print*,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

!...1.5: Set initial field
 if(nmeth.eq.1)then
   call inifield(unkno,uchar,geoel,coord,inpoel,iptri,ipqua)
 elseif(nmeth.eq.2)then
!...This will be overridden by 1.7
   if(ncurv.le.1)then
     call inifield_mcnew(unkno,uchar,geoel,coord,inpoel, iptri, ipqua)
   elseif(ncurv.gt.1)then
     call inifield_curv(unkno,uchar,geoel,coord,inpoel,iptri,ipqua)
   endif
 endif
 
 print*,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'Successfully set the initial flow field'
 print*,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

!...1.6: Calculate the mass matrix in Lagrangian frame, npoly=1 for DGP1, and npoly=2 for DGP2...
if(nmeth.eq.2)then

  !...Set the initial nodal velocity
   ustar = 0.d0
   coold = coord
!
  if(npoly==0)then
    nmatr = 1
  elseif(npoly==1)then
    nmatr = 4            !The dimension of mass matrix,p1
  elseif(npoly==2)then
    nmatr = 16           !The dimension of mass matrix,p2
  endif

  allocate(amatr(nmatr, ncell))
  allocate(amagd(nmatr, ncell))

  if(ncurv.le.1)then
   call getamatr_initial_lag(unkno,amatr,geoel,coord,iptri, ipqua)
   call getamatr_lag_dg(amagd,geoel,coord,iptri, ipqua)
  elseif(ncurv.gt.1)then
   call getamatr_initial_lag_curv(unkno,amatr,geoel,coord,iptri, ipqua)

  !...Implement gradient deformation in future
  endif

!...Set initial nodal velocity...
  if(nmatel.eq.1)then
   if(ncase.eq.2.or.ncase.eq.13)then
    call getustar_slnohmc(ustar, bface, intfac, coord)
   endif
  endif

!...1.7: ReSet initial field, overide 1.5 for initial condition
if(npoly.ge.1)then
 if(ncurv.le.1)then
  call getunkno_initial_lag(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
 elseif(ncurv.gt.1)then
  call getunkno_initial_lag_curv(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
 endif
endif
!!   call getunknot0_lag(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
!!   call getunkini_lag(unkno,amatr,geoel,coord,inpoel, iptri, ipqua)
!
!print*,'amatr',amatr(:,1)

!...Get initial deformation gradient
unkgd = 0.d0
unkgd(1, 1, : ) = 1.d0
unkgd(1, 2, : ) = 0.d0
unkgd(1, 3, : ) = 0.d0
unkgd(1, 4, : ) = 1.d0

!...Zero out the plastic deviatoric strain
strnq_devtp = 0.d0

!...Get the initial geometry parameters for the gardient of deformation
if(ncurv.eq.0)then
 call getfnds_lag_hybridgd(geoel,gflag,gesgq0,gesgt0,intfac,iptri,ipqua,coord,unkgd)
elseif(ncurv.eq.1)then
  call  getfnds_laghybrid_gd(geoel,gflag,gesgq0,gesgt0,intfac,iptri,ipqua,coord,unkgd)
! call getfnds_lagsmsif_hybrid(gflag,gesgt0,gesgq0,intfac,iptri,ipqua,coord)
! call getfnds_lagsms_simpsonhybrid(gflag,gesgq0,gesgt0,intfac,inpoel,iptri,ipqua,coord)
endif

!print*,'ustar',unkno(:,:,1),unkno(:,:,2),unkno(:,:,26),unkno(:,:,27)

!...1.8 Color some boundary cells for kidder shell for L2 error calculation in RZ
!  if(ncase.eq.4.or.ncase.eq.10)then
!    call getcellbd(geoel,coord,iptri,ipqua)
!  endif

!...Get element surrounding element
  call getesuel(intfac,ipqua,iptri,estri,esqua)

 endif !if(nmeth.eq.2)then

!...1.9 Define the RK coefficients
rklag(1, 1:3) = [0.d0, 0.5d0, 0.d0]
rklag(2, 1:3) = [0.d0,0.75d0,0.3333333333333333d0]
!
crklg(1:2) = [0.d0, 1.d0]
crklg2(1:3) = [0.d0, 1.d0, 0.5d0]
!
!xxxxxxx Step 2, solver xxxxxxxxx
!
!...Open residual file 'resi.dat'...
   open(7,file='resi.dat')

!...Record the old unknown values...
   unold = unkno
!...Deformation gradient
   unogd = unkgd

!...Record the old coordinates and deformation gradient
   if(nmeth.eq.2)then
      cooro = coord
      unogd = unkgd
   endif

!...Output the inital files for making movie...
if(ndump.ne.0)then
!...Mark the troubled-cells
if(npoly.gt.1) call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)
!
if((ncurv.eq.0.and.nfint.eq.2).or.(ncurv.eq.1.and.nfint.eq.8))then
call getfile_animation_particlevtk(coord,geoel,unkno,iptri, ipqua)
else
call getfile_animation_hybridvtk(inpoel,coord,geoel,unkno,iptri, ipqua)
endif
endif

!...Record the initial time
 call cpu_time(tstart_cpu)

!...time loop...

  do 1000 itime=1,ntime

    !...Zero out resnorm...
    resnorm=0.d0

!...Call DG solver
!...nmeth==1, Eulerian frame; nmeth==2, Lagrangian;
    if(nmeth.eq.1)then
      call getsolvdg(intfac,inpoel,iptri, ipqua,bface,coord,geofa,geoel,uchar,unkno,unold,ltelem,fsuel,itime,amatr)

    elseif(nmeth.eq.2)then
!...Degenerate to p0
!    unold(2:3,:,:)=0.d0;     unkno(2:3,:,:)=0.d0;
      call getsolvdg_lag(intfac,inpoel,iptri,ipqua,bface,cooro,coord,coora,geofa,geoel,gesgt0,gesgq0,uchar,&
                          unkno,unold,unkaw,unkgd,unogd,strnq_devtp,ltelem,estri,esqua,itime,amatr,amagd,ustar,coold, esuv1, esuv2)
    endif

!...Convergence indicator
    do ie=1,ncell
    do iq = 1, nq
       resnorm(1:ndegr, iq)=resnorm(1:ndegr, iq)+((unkno(1:ndegr,iq,ie)-unold(1:ndegr,iq,ie))**2)*geoel(4,ie)
    enddo
    enddo

!if(itime.gt.1000)print*,'bad',itime,geoel(4,20),unold(1,1,19)

!...Record the initial residual...
    if(itime==1) resinit(1:ndegr,1:nq)=resnorm(1:ndegr,1:nq)
    if (mod(itime, nresi) .eq. 0.or.(abs(tacru-tend)/tend.le.1e-12))then
      write (*,'(i8, 12e30.16)') itime,dtfix,tacru,sqrt(resnorm(1,1:nq))
!
      if(nmatel.eq.1)then
        write(7,'(i8, 24e32.16)')itime,0.5d0*log10(resnorm(:,1:nq))
      elseif(nmatel.eq.2)then
        call GetEnergy_solid(unkno,geoel,enrgy_solid)
        write(7,'(1e32.16, 12e32.16)')(itime-1)*dtfix,0.5d0*log10(resnorm(1,1:nq)),enrgy_solid(1:3),coord(1,1),coord(1,21),&
        ustar(2,368)
      endif
    endif

!...For Eulerian framework
    if(resnorm(1, 1).lt.toler) goto 6000
!
!   unold(5:6,:,:) = 0.d0
!   unold(3,:,:) = 0.d0
!   unold(:,3,:) = 0.d0
!   print*,'unkno',rkstg, unold(:,3,2)

!...Update the array unkno...
     do ie =1,ncell
     do iq =1,ndegr
        unkno(iq,1:nq,ie)=unold(iq,1:nq,ie)
     enddo 
     enddo

!...Update the gradient of deformation
    do ie =1,ncell
    do iq =1,ndegr
       unkgd(iq,1:4,ie)=unogd(iq,1:4,ie)
    enddo
    enddo

!...Update the coordinates
    if(nmeth.eq.2) coord = cooro

!...Correct the negative cells
!    if(mod(itime, nresi) .eq. 0)then!itime.gt.12000)then
!    if(ncurv.eq.1) call getmesh_lagcrect(intfac, inpoel, ipqua, coord)
!    endif

!...Output files for making movie...
    if(ndump.ne.0.and.mod(itime,ndump).eq.0)then

!...Output geometry info
    if(ncurv==2)call getoutput_cellgeo(inpoel,coord,geoel,unkno,iptri, ipqua)

     if(nmatel.eq.1)then
      if(ncurv.le.1)then
       if((ncurv.eq.0.and.nfint.eq.2).or.(ncurv.eq.1.and.nfint.eq.8))then
        call getfile_animation_particlevtk(coord,geoel,unkno,iptri, ipqua)
       else
       ! call getSIHO_general_ml(intfac, ipqua, coord, coold, geoel, unkno)
        call getfile_animation_hybridvtk(inpoel,coord,geoel,unkno,iptri, ipqua)
       endif
      else !...High-order mesh
        call getfile_animation_cellvtk_ho(coord,geoel,unkno,iptri, ipqua)
      endif
     else
       call getfile_animation_hybridvtk(inpoel,coord,geoel,unkno,iptri, ipqua)
     endif
    endif

!...Degenerate to p0
!unold(2:3,2:3,:)=0.d0;     unkno(2:3,2:3,:)=0.d0;
!    if(itime.eq.50000.or.itime.eq.100000)then
!     if(ncurv.eq.0)call geterror_lag_hybrid(intfac, inpoel,iptri, ipqua, coord, coold, geoel, unkno)
!     if(ncurv.eq.1)call geterror_lag_curvhybrid(intfac, inpoel, ipqua, coord, coold, unkno, geoel)
!    endif

!...Exit upon the final time
  if(abs(tacru-tend)/tend.le.1e-11)then
   print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXX'
   print*,'Reach the final time!',tacru
   print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXX'
   goto 6000
  endif

1000 enddo

6000 continue

!print*,'unkno-final',unkno(:,:,1),unkno(:,:,2),unkno(:,:,26),unkno(:,:,27)
!print*,'coord-final',coord(:,ipqua(1:4,1)),coord(:,ipqua(1:4,2))


!...Close file resi.dat...
    close(7)

!...Record the final time
    call cpu_time(tend_cpu)

   print*,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
   print*,'Complete the solver using time',tend_cpu-tstart_cpu
!
!xxxxxxx Step 3, Postprocessingxxxxxxxxx
!
!...3.0: output restart file...
   print*,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
   print*,'Output restart file and tecplot file'!,coord(1:2, 36)
!
  open(12,file='restartdg.dat')
    do ie=1,ncell
       write(12,'(3e32.16)')unkno(1:ndegr,1,ie)
       write(12,'(3e32.16)')unkno(1:ndegr,2,ie)
       write(12,'(3e32.16)')unkno(1:ndegr,3,ie)
       write(12,'(3e32.16)')unkno(1:ndegr,4,ie)
    enddo
  close(12)

!...3.1: output tecplot format...
  print*,'Complete the solver'
!  call getpostp1(unkno,inpoel,coord,geoel)
!...Eulerian method
   if(nmeth.eq.1)then
    if(nout.eq.1) call getoutput_cell(inpoel,coord,geoel,unkno)
    if(nout.eq.2) call getoutput_cellpara_euler(inpoel,coord,unkno,iptri,ipqua)
   endif

!....Lagrangian method
   if(nmeth.eq.2)then
    select case (nmatel)

!...Gas
     case (1)
!...Update the cell volume
     call getgeoel_lag_post(iptri, ipqua, geoel, coord)
!     if(ndens.eq.2) call getdensity_lag(unkno)
     if(nout.eq.1) call getoutput_celllag(inpoel,ipqua,coord,geoel,unkno)
     if(nout.eq.2)then
!    call getpostnd_lag(unkno,inpoel,coord,geoel,unknp)
     call getpostnd_lag_hybrid(unkno,inpoel,iptri, ipqua, coord,geoel,unknp)

!...Idenity the troubled-cells from Persson's
    if(npoly.gt.1)then
      call getSIHO_Persson(intfac, ipqua, coord, coold, geoel, unkno)
      !call output_csv(intfac, ipqua, coord, coold, geoel, unkno)
    endif

!...Get the minimum Jacobian
     call getJacobian_lagcurv(iptri, ipqua, coord, ujacb)
!...Output final file
    if(ncurv.le.1)then
    !...vtu file
!      call getoutput_cellpara_hybrid(inpoel,coord,geoel,unkno,ujacb, iptri,ipqua)

!...vtk file
      call getoutput_cellvtk(inpoel,coord,geoel,unkno,iptri, ipqua)
    else
    !...vtk file
      call getoutput_cellvtk_ho(inpoel,coord,geoel,unkno,iptri, ipqua)
    endif

!...Output geometry info
  if(ncurv==2)call getoutput_cellgeo(inpoel,coord,geoel,unkno,iptri, ipqua)

!...vtk file for particle method.
!   if(nfint.eq.2)then
     call getoutput_particlevtk(coord,geoel,unkno,iptri, ipqua)
!   endif
    endif

!...Solid
    case (2)
     call getoutputsolid_vtk_cellhybrid(coord,coold,geoel,unkno,unkgd,strnq_devtp,iptri, ipqua)
    end select
   endif

!...3.2: Calculate the error...
   if(nmeth.eq.1)then
    if(ncurv.le.1)call geterror(inpoel, iptri, ipqua, geoel, coord, unkno)
    if(ncurv.eq.2)call geterror_curv2(inpoel, geoel, coord, unkno)
   elseif(nmeth.eq.2)then
    if(nfint.eq.2.and.ncurv.eq.0)then
      call geterror_lag_hybrid_gd(intfac,iptri, ipqua, coord, coold, geoel, unkno, unkgd)
    else
     if(ncurv.eq.0)call geterror_lag_hybrid(intfac, inpoel,iptri, ipqua, coord, coold, geoel, unkno)
     if(ncurv.eq.1)call geterror_lag_curvhybrid(intfac, inpoel, ipqua, coord, coold, unkno, geoel)
     if(ncurv.eq.2)call geterror_lag_curvhybrid_general(intfac, inpoel, ipqua, coord, coold, unkno, geoel)

   ! if(ncurv.eq.1)call geterror_lagtgvfvm(intfac, iptri, ipqua, coord, coold, geoel, unkno)
   ! call geterror_lag_curv(intfac, inpoel, coord, unkno)
    endif
   endif
end program main
!
!...Subroutine: Set the initial flow field...
!
 subroutine inifield(unkno,uchar,geoel,coord,inpoel, iptri, ipqua)
 use constant
 implicit none
 real*8,dimension(1:mdegr,1:nq,1:nsize),intent(out)::unkno
 real*8,dimension(1:nq)::uchar
 real*8,dimension(1:ngeel,1:nsize)::geoel
 real*8,dimension(1:ndimn, 1:npoin),intent(in)::coord
 integer,dimension(1:nvtri, 1:nelem)::inpoel
 integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
 integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
!...Local
!...local real array
 real*8:: coorp(1:ndimn, 1:nvtri)
 real*8:: shp(nptri), dspr(nptri), dsps(nptri)
!...local real
 real*8:: dp0dx, dp0dy, du0dx,du0dy,dv0dx,dv0dy
 real*8:: du0dr, du0ds, dv0dr, dv0ds,de0dx,de0dy
 real*8:: dr, ds
 real*8:: dxdr,dxds,dydr,dyds
 real*8:: dx21, dy21, dx31, dy31
 real*8:: volel,wi,xc,yc,xg,yg
 real*8:: uvel, vvel, evel, ieve
 real*8:: pini, uini, vini
 real*8:: uctr,vctr

 real*8:: c00,c10,c05,c20,eps
 real*8::r, s
real*8::radie, radii,radie2,radii2,radic2, radig2,paras,rhoi,rhoe
real*8::rho, rho0, rhom, drho0dx, drho0dy, drhomdx, drhomdy, drhomdr, drhomds, pctr
!...local integer
 integer*4::ie,id, ig, ishp
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
!...give the far-field value--(nondimensionalized value)...
!
  uchar(1) = 1.d0
  uchar(2) = amach
  uchar(3) = 0.d0
  uchar(4) = 1.d0/gamma/(gamma-1.d0)+0.5d0*amach**2 !..initial rho*e

!
!...Zero out everything...
!
   unkno = 0.d0
!
!...initilization of flow field...
!
 if(nmeth.eq.1)then !...nmeth=1 for Euleriam frame...
!
  if(nint==1)then
     do ie=1,ncell
        unkno(1,1:nq,ie)=uchar(1:nq)
     enddo
  elseif(nint==2)then
    open(7,file='restartdg.dat')
     do ie=1,ncell
      read(7,*)unkno(1:ndegr,1,ie)
      read(7,*)unkno(1:ndegr,2,ie)
      read(7,*)unkno(1:ndegr,3,ie)
      read(7,*)unkno(1:ndegr,4,ie)
     ! read(7,*)unkno(1:3,1,ie)
     ! read(7,*)unkno(1:3,2,ie)
     ! read(7,*)unkno(1:3,3,ie)
     ! read(7,*)unkno(1:3,4,ie)
     enddo
    close(7)
  endif
!
 elseif(nmeth.eq.2)then !...nmeth=2 for Euleriam frame...
!
  if(nint==1)then
  if(ncase.eq.1)then   !...ncase=1 for Taylor-Green Vortex(Only linear meshes act as initial mesh)
   do ie = 1, nelem
!
!...Physical coordinates
!
   coorp(1, 1:nvtri) = coord(1, inpoel(1:nvtri, ie))
   coorp(2, 1:nvtri) = coord(2, inpoel(1:nvtri, ie))
!
!...Reference coordinates
!
   xc = geoel(1,ie)
   yc = geoel(2,ie)
   volel = geoel(3, ie)
!
   dr = 0.5d0  !...Refrence cell (Delt x)
   ds = 0.5d0
!
   uvel = 0.d0
   vvel = 0.d0
   evel = 0.d0
!
   do ig =1,ngausd
!
   r  = posi(1,ig)
   s  = posi(2,ig)
!  
   shp(1) = 1.d0-r-s
   shp(2) = r
   shp(3) = s
!
   wi = weigh(ig)*volel
!...Linear mesh
   xg = shp(1)*coorp(1, 1) + shp(2)*coorp(1, 2) + shp(3)*coorp(1, 3)
   yg = shp(1)*coorp(2, 1) + shp(2)*coorp(2, 2) + shp(3)*coorp(2, 3)
!
   uini = sin(pi*xg)*cos(pi*yg)
   vini =-cos(pi*xg)*sin(pi*yg)
   pini = 0.25d0*(cos(2.d0*pi*xg) + cos(2.d0*pi*yg)) + 1.d0
!
   uvel = uvel + uini*wi
   vvel = vvel + vini*wi
   evel = evel + (pini/(gamlg-1.d0) + 0.5d0*(uini**2 + vini**2))*wi
!
! if(ie==21.or.ie==5.or.ie==4)print*,'initial ig',ie,ig,uini,-cos(pi*xg)*sin(pi*yg),xg,yg,xc,yc,nvtri,ig,ngausd
!
   enddo
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
   do ishp = 1, 3  !....Only 3 points for linear mesh...
   dxdr = dxdr + dspr(ishp)*coorp(1,ishp)
   dxds = dxds + dsps(ishp)*coorp(1,ishp)

   dydr = dydr + dspr(ishp)*coorp(2,ishp)
   dyds = dyds + dsps(ishp)*coorp(2,ishp)
   enddo
!
!
   dx21 = coorp(1, 2) - coorp(1, 1)
   dx31 = coorp(1, 3) - coorp(1, 1)
   dy21 = coorp(2, 2) - coorp(2, 1)
   dy31 = coorp(2, 3) - coorp(2, 1)
!
   unkno(1, 1, ie) = 1.d0; unkno(2:ndegr, 1, ie) = 0.d0;
!
   unkno(1, 2, ie) = uvel/volel;
!
   du0dx = pi*cos(pi*xc)*cos(pi*yc)
   du0dy =-pi*sin(pi*xc)*sin(pi*yc)
!
 !
   unkno(1, 3, ie) = vvel/volel;
!
   dv0dx = pi*sin(pi*xc)*sin(pi*yc)
   dv0dy =-pi*cos(pi*xc)*cos(pi*yc)
!
   du0dr = du0dx*dxdr + du0dy*dydr
   du0ds = du0dx*dxds + du0dy*dyds
!
   dv0dr = dv0dx*dxdr + dv0dy*dydr
   dv0ds = dv0dx*dxds + dv0dy*dyds
!
   unkno(2, 2, ie) = du0dr
   unkno(3, 2, ie) = du0ds
!
   unkno(2, 3, ie) = dv0dr
   unkno(3, 3, ie) = dv0ds
!
!
   dp0dx = 0.25d0*(-2.d0*pi)*(sin(2.d0*pi*xc))
   dp0dy = 0.25d0*(-2.d0*pi)*(sin(2.d0*pi*yc))
!
   uctr = sin(pi*xc)*cos(pi*yc)
   vctr =-cos(pi*xc)*sin(pi*yc)
!
!  if(ie==1)print*,'initial',du0dx, du0dy,uctr,vctr,xc,yc
!   if(ie==21.or.ie==5.or.ie==4)print*,'initial',ie,unkno(1, 3, ie),-cos(pi*xc)*sin(pi*yc),vvel,volel!+dv0dx*(0.8d0-xc) + dv0dy*(0.d0-yc)
!
   unkno(1, 4, ie) = evel/volel
   de0dx = dp0dx/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
   de0dy = dp0dy/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)
!
   unkno(2, 4, ie) = de0dx*dxdr + de0dy*dydr
   unkno(3, 4, ie) = de0dx*dxds + de0dy*dyds
!
!...Update the dierivative with drr and dsr...
!
   unkno(2, 1:nq, ie) = unkno(2, 1:nq, ie)*dr
   unkno(3, 1:nq, ie) = unkno(3, 1:nq, ie)*ds
!
   enddo
!
elseif(ncase.eq.2)then !...ncase==2 : shockless Noh test problems...
!
do ie = 1, nelem
!
!...Physical coordinates
!
coorp(1, 1:nvtri) = coord(1, inpoel(1:nvtri, ie))
coorp(2, 1:nvtri) = coord(2, inpoel(1:nvtri, ie))
!
!...Reference coordinates
!
xc = geoel(1,ie)
yc = geoel(2,ie)
volel = geoel(3, ie)
!
dr = 0.5d0  !...Refrence cell (Delt x)
ds = 0.5d0
!
uvel = 0.d0
vvel = 0.d0
evel = 0.d0
!
do ig =1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
!
wi = weigh(ig)*volel
!
xg = shp(1)*coorp(1, 1) + shp(2)*coorp(1, 2) + shp(3)*coorp(1, 3)
yg = shp(1)*coorp(2, 1) + shp(2)*coorp(2, 2) + shp(3)*coorp(2, 3)
!
uini = -xg
vini = -yg
pini = (gamlg-1.d0)
!
uvel = uvel + uini*wi
vvel = vvel + vini*wi
evel = evel + (pini/(gamlg-1.d0) + 0.5d0*(uini**2 + vini**2))*wi
!
! if(ie==21.or.ie==5.or.ie==4)print*,'initial ig',ie,ig,uini,-cos(pi*xg)*sin(pi*yg),xg,yg,xc,yc,nvtri,ig,ngausd
!
enddo
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
do ishp = 1, 3  !...Only linear meshes...
dxdr = dxdr + dspr(ishp)*coorp(1,ishp)
dxds = dxds + dsps(ishp)*coorp(1,ishp)

dydr = dydr + dspr(ishp)*coorp(2,ishp)
dyds = dyds + dsps(ishp)*coorp(2,ishp)
enddo
!
!
!
unkno(1, 1, ie) = 1.d0; unkno(2:ndegr, 1, ie) = 0.d0;
!
unkno(1, 2, ie) = uvel/volel;
!
 du0dx =-1.d0
 du0dy = 0.d0
!
!
unkno(1, 3, ie) = vvel/volel;
!
dv0dx =  0.d0
dv0dy = -1.d0
!
 du0dr = du0dx*dxdr + du0dy*dydr
 du0ds = du0dx*dxds + du0dy*dyds
!
 dv0dr = dv0dx*dxdr + dv0dy*dydr
 dv0ds = dv0dx*dxds + dv0dy*dyds
!
 unkno(2, 2, ie) = du0dr
 unkno(3, 2, ie) = du0ds
!
 unkno(2, 3, ie) = dv0dr
 unkno(3, 3, ie) = dv0ds
!
!
 dp0dx = 0.d0
 dp0dy = 0.d0
!
 uctr  = -xc
 vctr  = -yc
!
!  if(ie==1)print*,'initial',du0dx, du0dy,uctr,vctr,xc,yc
!   if(ie==21.or.ie==5.or.ie==4)print*,'initial',ie,unkno(1, 3, ie),-cos(pi*xc)*sin(pi*yc),vvel,volel!+dv0dx*(0.8d0-xc) + dv0dy*(0.d0-yc)
!
 unkno(1, 4, ie) = evel/volel
 de0dx = dp0dx/(gamlg-1.d0) + (uctr*du0dx + vctr*dv0dx)
 de0dy = dp0dy/(gamlg-1.d0) + (uctr*du0dy + vctr*dv0dy)
!
 unkno(2, 4, ie) = de0dx*dxdr + de0dy*dydr
 unkno(3, 4, ie) = de0dx*dxds + de0dy*dyds
!
!...Update the dierivative with drr and dsr...
!
 unkno(2, 1:nq, ie) = unkno(2, 1:nq, ie)*dr
 unkno(3, 1:nq, ie) = unkno(3, 1:nq, ie)*ds
!
enddo
!
elseif(ncase.eq.4)then !...ncase==4 : Kidder shell test problems...
!
!...some parameters
!
radie = 1.d0
radii = 0.9d0
rhoi  = 1.d-3
rhoe  = 1.d-2
paras = 1.d5
!
!print*,'initial ig',ie
!
radie2 = radie**2
radii2 = radii**2
!
do ie = 1, nelem
!
!...Physical coordinates
!
coorp(1, 1:nvtri) = coord(1, inpoel(1:nvtri, ie))
coorp(2, 1:nvtri) = coord(2, inpoel(1:nvtri, ie))
!
!...Reference coordinates
!
xc = geoel(1,ie)
yc = geoel(2,ie)
volel = geoel(3, ie)
!
dr = 0.5d0  !...Refrence cell (Delt x)
ds = 0.5d0
!
uvel = 0.d0
vvel = 0.d0
evel = 0.d0
rhom = 0.d0
!
do ig =1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
!
wi = weigh(ig)*volel
!
xg = shp(1)*coorp(1, 1) + shp(2)*coorp(1, 2) + shp(3)*coorp(1, 3)
yg = shp(1)*coorp(2, 1) + shp(2)*coorp(2, 2) + shp(3)*coorp(2, 3)
!
radig2 = xg**2 + yg**2
!
rho =(radie2-radig2)/(radie2-radii2)*rhoi + (radig2-radii2)/(radie2-radii2)*rhoe !gamma=2 means gamma-1 = 1
!
uini = 0.d0
vini = 0.d0
pini = paras*(rho)**gamlg
!
rhom = rhom + 1.d0/rho*wi
uvel = uvel + uini*wi
vvel = vvel + vini*wi
evel = evel + (pini/rho/(gamlg-1.d0))*wi
!
! if(ie==200)print*,'initial ig',ie,ig,uini,pini,xg,yg,xc,yc,nvtri,gamlg, volel,radie2, radii2, paras
!
enddo
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
!...Density inverse...
unkno(1, 1, ie) = rhom/volel;
!
radic2 = xc**2 + yc**2
!
rho0 = (radie2-radic2)/(radie2-radii2)*rhoi + (radic2-radii2)/(radie2-radii2)*rhoe !...density at cell center...
drho0dx = 2.d0*xc*(rhoe-rhoi)/(radie2-radii2)
drho0dy = 2.d0*yc*(rhoe-rhoi)/(radie2-radii2)
!
drhomdx = -drho0dx/rho0/rho0
drhomdy = -drho0dy/rho0/rho0
!
drhomdr = drhomdx*dxdr + drhomdy*dydr
drhomds = drhomdx*dxds + drhomdy*dyds
!
unkno(2, 1, ie) = drhomdr
unkno(3, 1, ie) = drhomds
!
unkno(1, 2, ie) = 0.d0;
unkno(1, 3, ie) = 0.d0;
!
!
unkno(2, 2, ie) = 0.d0
unkno(3, 2, ie) = 0.d0
!
unkno(2, 3, ie) = 0.d0
unkno(3, 3, ie) = 0.d0
!
!   dp0dx = 0.25d0*(-2.d0*pi)*(sin(2.d0*pi*xc)*(dx21) +  &
!                             sin(2.d0*pi*yc)*(dy21))
!   dp0dy = 0.25d0*(-2.d0*pi)*(sin(2.d0*pi*xc)*(dx31) +  &
!                             sin(2.d0*pi*yc)*(dy31))
!
!dp0dx = paras*gamlg*rho0**(gamlg-1.d0)*drho0dx
!dp0dy = paras*gamlg*rho0**(gamlg-1.d0)*drho0dy
!
dp0dx = -paras*gamlg*rho0**(gamlg+1.d0)*drhomdx
dp0dy = -paras*gamlg*rho0**(gamlg+1.d0)*drhomdy
!
uctr  = 0.d0
vctr  = 0.d0
pctr  = paras*(rho0)**gamlg
!
  if(ie==200)print*,'initial',xc, yc,radic2,rho0,drho0dx,drhomdx,drhomdy,dxdr,dydr,drhomdr,drhomds
  if(ie==200)print*,coorp(1, 1)-xc, coorp(2, 1)-yc, drho0dx, drho0dy, rho0,rho0+drho0dx*(coorp(1, 1)-xc)+drho0dy*(coorp(2, 1)-yc)
  if(ie==200)print*,coorp(1, 1)-xc, coorp(2, 1)-yc, drhomdx, drhomdy, 1.d0/rho0,&
                    unkno(1, 1, ie),1.d0/rho0+drhomdx*(coorp(1, 1)-xc)+drhomdy*(coorp(2, 1)-yc),&
                    1.d0/rho0+drhomdr*(-1.d0/3.d0)+drhomds*(-1.d0/3.d0)
!   if(ie==21.or.ie==5.or.ie==4)print*,'initial',ie,unkno(1, 3, ie),-cos(pi*xc)*sin(pi*yc),vvel,volel!+dv0dx*(0.8d0-xc) + dv0dy*(0.d0-yc)
!
unkno(1, 4, ie) = evel/volel
!de0dx = (dp0dx*rho0 - pctr*drho0dx)/(rho0**2)/(gamlg-1.d0)
!de0dy = (dp0dy*rho0 - pctr*drho0dy)/(rho0**2)/(gamlg-1.d0)
!
de0dx = (dp0dx/rho0 + pctr*drhomdx)/(gamlg-1.d0)
de0dy = (dp0dy/rho0 + pctr*drhomdy)/(gamlg-1.d0)
!
if(ie==200)print*,'initial',pctr, rho0, dp0dx, dp0dy,radic2,de0dx, de0dy
!
unkno(2, 4, ie) = de0dx*dxdr + de0dy*dydr
unkno(3, 4, ie) = de0dx*dxds + de0dy*dyds
!
if(ie==200)print*,coorp(1, 1)-xc, coorp(2, 1)-yc, de0dx, de0dy,&
unkno(1, 4, ie),de0dx*(coorp(1, 1)-xc)+de0dy*(coorp(2, 1)-yc),&
unkno(2, 4, ie)*(-1.d0/3.d0)+unkno(3, 4, ie)*(-1.d0/3.d0)
!
!...Update the dierivative with drr and dsr...
!
unkno(2, 1:nq, ie) = unkno(2, 1:nq, ie)*dr
unkno(3, 1:nq, ie) = unkno(3, 1:nq, ie)*ds
!
enddo
!
elseif(ncase.eq.5)then   !...ncase=5 for Kidder ball
!
do ie = 1, nelem
!
!...Physical coordinates
!
coorp(1, 1:nvtri) = coord(1, inpoel(1:nvtri, ie))
coorp(2, 1:nvtri) = coord(2, inpoel(1:nvtri, ie))
!
!...Reference coordinates
!
xc = geoel(1,ie)
yc = geoel(2,ie)
volel = geoel(3, ie)
!
dr = 0.5d0  !...Refrence cell (Delt x)
ds = 0.5d0
!
rhom = 0.d0
uvel = 0.d0
vvel = 0.d0
evel = 0.d0
!
do ig =1,ngausd
!
r  = posi(1,ig)
s  = posi(2,ig)
!
shp(1) = 1.d0-r-s
shp(2) = r
shp(3) = s
!
wi = weigh(ig)*volel
!
xg = shp(1)*coorp(1, 1) + shp(2)*coorp(1, 2) + shp(3)*coorp(1, 3)
yg = shp(1)*coorp(2, 1) + shp(2)*coorp(2, 2) + shp(3)*coorp(2, 3)
!
rho  = 1.d0/sqrt(2.d0)*exp(-0.5d0*(xg**2 + yg**2))
uini = -0.5d0*xg
vini = -0.5d0*yg
ieve = 3.d0/8.d0
!
rhom = rhom + 1.d0/rho*wi
uvel = uvel + uini*wi
vvel = vvel + vini*wi
evel = evel + (ieve + 0.5d0*(uini**2 + vini**2))*wi
!
! if(ie==21.or.ie==5.or.ie==4)print*,'initial ig',ie,ig,uini,-cos(pi*xg)*sin(pi*yg),xg,yg,xc,yc,nvtri,ig,ngausd
!
enddo
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
do ishp = 1, 3
dxdr = dxdr + dspr(ishp)*coorp(1,ishp)
dxds = dxds + dsps(ishp)*coorp(1,ishp)

dydr = dydr + dspr(ishp)*coorp(2,ishp)
dyds = dyds + dsps(ishp)*coorp(2,ishp)
enddo
!
!
dx21 = coorp(1, 2) - coorp(1, 1)
dx31 = coorp(1, 3) - coorp(1, 1)
dy21 = coorp(2, 2) - coorp(2, 1)
dy31 = coorp(2, 3) - coorp(2, 1)
!
unkno(1, 1, ie) = rhom/volel;
!
rho0 = 1.d0/sqrt(2.d0)*exp(-0.5d0*(xc**2 + yc**2))
!
drho0dx = -xc*rho0
drho0dy = -yc*rho0
!
drhomdx = -drho0dx/rho0**2
drhomdy = -drho0dy/rho0**2
!
drhomdr = drhomdx*dxdr + drhomdy*dydr
drhomds = drhomdx*dxds + drhomdy*dyds
!
unkno(2, 1, ie) = drhomdr
unkno(3, 1, ie) = drhomds
!
unkno(1, 2, ie) = uvel/volel;
du0dx = -0.5d0
du0dy =  0.d0
!
!
unkno(1, 3, ie) = vvel/volel;
dv0dx = 0.d0
dv0dy =-0.5d0
!
du0dr = du0dx*dxdr + du0dy*dydr
du0ds = du0dx*dxds + du0dy*dyds
!
dv0dr = dv0dx*dxdr + dv0dy*dydr
dv0ds = dv0dx*dxds + dv0dy*dyds
!
unkno(2, 2, ie) = du0dr
unkno(3, 2, ie) = du0ds
!
unkno(2, 3, ie) = dv0dr
unkno(3, 3, ie) = dv0ds
!
!   dp0dx = 0.25d0*(-2.d0*pi)*(sin(2.d0*pi*xc)*(dx21) +  &
!                             sin(2.d0*pi*yc)*(dy21))
!   dp0dy = 0.25d0*(-2.d0*pi)*(sin(2.d0*pi*xc)*(dx31) +  &
!                             sin(2.d0*pi*yc)*(dy31))
!
uctr = -0.5d0*xc
vctr = -0.5d0*yc
!
!  if(ie==1)print*,'initial',du0dx, du0dy,uctr,vctr,xc,yc
!   if(ie==21.or.ie==5.or.ie==4)print*,'initial',ie,unkno(1, 3, ie),-cos(pi*xc)*sin(pi*yc),vvel,volel!+dv0dx*(0.8d0-xc) + dv0dy*(0.d0-yc)
!
unkno(1, 4, ie) = evel/volel
de0dx = (uctr*du0dx + vctr*dv0dx)
de0dy = (uctr*du0dy + vctr*dv0dy)
!
unkno(2, 4, ie) = de0dx*dxdr + de0dy*dydr
unkno(3, 4, ie) = de0dx*dxds + de0dy*dyds
!
!...Update the dierivative with drr and dsr...
!
unkno(2, 1:nq, ie) = unkno(2, 1:nq, ie)*dr
unkno(3, 1:nq, ie) = unkno(3, 1:nq, ie)*ds
!
enddo


endif
!
elseif(nint==2)then
!
!Input restarted file...
!
  print*,'...................Notice................... !'
  print*,'Restart function is not allowed by Lagrangian framework right now !'
  stop
 endif
!
 endif
 end subroutine inifield
!
!...Finding element surrounding element...
!
subroutine getesuel(intfac,ipqua,iptri,estri,esqua)
use constant
implicit none
integer*4,dimension(1:nifai,1:nafac), intent(in)::intfac
integer,  dimension(1:nvtri,1:ntria), intent(in):: iptri
integer,  dimension(1:nvqua,1:nquad), intent(in):: ipqua
integer, dimension(1:nftri,1:ntria), intent(out)::estri
integer, dimension(1:nfqua,1:nquad), intent(out)::esqua
integer*4::ie,ifa,numb, ielem
!
!...For tria
!
do ie=1,ntria
!
ielem = ie
!
numb=0
do ifa=1,nafac
if(intfac(1,ifa).eq.ielem)then
numb=numb+1
if(intfac(3, ifa).eq.iptri(1, ie))then
estri(1,ie)=intfac(2,ifa)
elseif(intfac(3, ifa).eq.iptri(2, ie))then
estri(2,ie)=intfac(2,ifa)
elseif(intfac(3, ifa).eq.iptri(3, ie))then
estri(3,ie)=intfac(2,ifa)
else
print*,'We can not find element surrounding the cell in getesuel!'
stop
endif
elseif(intfac(2,ifa).eq.ielem)then
numb=numb+1
if(intfac(4, ifa).eq.iptri(1, ie))then
estri(1,ie)=intfac(1,ifa)
elseif(intfac(4, ifa).eq.iptri(2, ie))then
estri(2,ie)=intfac(1,ifa)
elseif(intfac(4, ifa).eq.iptri(3, ie))then
estri(3,ie)=intfac(1,ifa)
else
print*,'We can not find element surrounding the cell in getesuel!'
stop
endif
endif
enddo
!
enddo
!
!...For quad
!
do ie=1,nquad
!
ielem = ie +ntria
!
numb=0
do ifa=1,nafac
if(intfac(1,ifa).eq.ielem)then
numb=numb+1
 if(intfac(3, ifa).eq.ipqua(1, ie))then
   esqua(1,ie)=intfac(2,ifa)
 elseif(intfac(3, ifa).eq.ipqua(2, ie))then
   esqua(2,ie)=intfac(2,ifa)
 elseif(intfac(3, ifa).eq.ipqua(3, ie))then
   esqua(3,ie)=intfac(2,ifa)
 elseif(intfac(3, ifa).eq.ipqua(4, ie))then
   esqua(4,ie)=intfac(2,ifa)
 else
 print*,'We can not find element surrounding the cell in getesuel!'
 stop
 endif
elseif(intfac(2,ifa).eq.ielem)then
numb=numb+1
 if(intfac(4, ifa).eq.ipqua(1, ie))then
   esqua(1,ie)=intfac(1,ifa)
 elseif(intfac(4, ifa).eq.ipqua(2, ie))then
   esqua(2,ie)=intfac(1,ifa)
 elseif(intfac(4, ifa).eq.ipqua(3, ie))then
   esqua(3,ie)=intfac(1,ifa)
 elseif(intfac(4, ifa).eq.ipqua(4, ie))then
   esqua(4,ie)=intfac(1,ifa)
 else
  print*,'We can not find element surrounding the cell in getesuel!'
  stop
 endif

endif
enddo
enddo
!
!print*,'esqua', esqua(:,23)
end subroutine getesuel

!
!...Subrotuine:get the faces and elements around a element...
!
subroutine getlogeo(intfac,geoel,geofa,ltelem,fsuel,fsuel_q)
use constant
implicit none
integer*4,dimension(1:nifai,1:nafac)::intfac
real*8,dimension(1:ngeel,1:nsize)::geoel
real*8,dimension(1:ngefa,1:nafac)::geofa
integer*4,dimension(1:3,1:ntria), intent(out)::ltelem
integer*4,dimension(1:3,1:ntria), intent(out)::fsuel
integer*4,dimension(1:4,1:nquad), intent(out)::fsuel_q
integer*4::ie,ifa,numb

!...Triangle
do ie=1,ntria
    numb=0
 do ifa=1,nafac
     if(intfac(1,ifa).eq.ie)then
     numb=numb+1
     ltelem(numb,ie)=intfac(2,ifa)
     fsuel(numb, ie)=ifa
     elseif(intfac(2,ifa).eq.ie)then
     numb=numb+1
     ltelem(numb,ie)=intfac(1,ifa)
     fsuel(numb, ie)=ifa
     endif

 enddo
enddo

!...Quad
do ie=1,nquad
numb=0
do ifa=1,nafac
if(intfac(1,ifa).eq.ie)then
numb=numb+1
fsuel_q(numb, ie)=ifa
elseif(intfac(2,ifa).eq.ie)then
numb=numb+1
fsuel_q(numb, ie)=ifa
endif
enddo
enddo

end subroutine getlogeo
!
!...Subrotuine: get local time step...
!
  subroutine getdelt(deltt,geoel,bface,uchar,unkno,intfac,coord)
  use constant
  implicit none
  real*8,dimension(1:ngeel,1:nsize), intent(in)::geoel
  real*8,dimension(1:nq), intent(in)::uchar 
  real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
  real*8,dimension(1:ndimn,1:npoin), intent(in)::coord 
  real*8,dimension(1:ncell), intent(out)::deltt
  integer*4,dimension(1:nifai,1:nafac), intent(in)::intfac
  integer*4,dimension(1:nbfai,1:nbfac), intent(in)::bface
 !
 !...local array
  real*8,dimension(1:ncell)::veloc
 !
  integer:: ielem
 !
 !...1.find the characteristic velocity...
 if(ncurv.le.1) call getveloc_curv(geoel,bface,uchar,unkno,intfac,coord,veloc)
 if(ncurv.eq.2) call getveloc_curv2(geoel,bface,uchar,unkno,intfac,coord,veloc)
 !
 !...2.find the time step...
   do ielem = 1, ncell
      deltt(ielem) = cfl*geoel(3, ielem)/veloc(ielem)
   enddo
 !
  if(ndt.eq.2) deltt = min(dtfix, deltt(1:ncell))
 !
  end subroutine getdelt
!
!...Subroutine to calculate the characteristic velicity...
!
subroutine getveloc_curv(geoel,bface,uchar,unkno,intfac,coord,veloc)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nsize),intent(in)::geoel
 real*8,intent(in)::uchar(1:nq)
 real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
 real*8,dimension(1:ndimn,1:npoin)::coord 
 real*8,dimension(1:ncell), intent(out)::veloc
 integer*4,dimension(1:nifai,1:nafac)::intfac
 integer*4,dimension(1:nbfai,1:nbfac)::bface
!
!...local array
 real*8::xpin(2, 2)
 real*8::xp(2, nptfa)
 real*8,dimension(1:nptfa)::shp, dshpr
 real*8::weigh(ngausf), posi(1,ngausf)
 real*8,dimension(1:nq)::unkno1,unkno2
 real*8::flux(1:nq)
 real*8::bl(6), br(6)
!...local real number
 integer::ifa,iel,ier,ie,ifb,k,m,ig,ishp,ideg,jdeg
 real*8::xcl,ycl,dxl,dyl,rnx,rny
 real*8::xcr,ycr,dxr,dyr 
 real*8:: r, djaco, dxdr, dydr 
 real*8::xg,yg,wi
 real*8::rho1,rhou1,rhov1,rhoe1,rho2,rhou2,rhov2,rhoe2
 real*8::vn,dwav1,dwav2
 real*8::rhom1,uadv1,vadv1,pres1, csou1, vdon1
 real*8::rhom2,uadv2,vadv2,pres2, csou2, vdon2
 real*8::eps
!
 data   eps / 1.0d-06 /
!
 call rutope(1, ngausf, posi, weigh)
!
! print*,'gauss',posi(1,ngausf)
!
!...zero out rhs...
!
    veloc = 0.d0
!
!...Characteristic velocity from boundary faces...
!
    do 250 ifa = 1,nbfac !...(1)ifa=1,nafac
!
     iel=intfac(1,ifa)
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
     xg = 0.d0
     yg = 0.d0
!
   do ishp = 1, nptfa
     xg = xg + shp(ishp)*xp(1, ishp) 
     yg = yg + shp(ishp)*xp(2, ishp) 
    enddo
!
!...Basis function...
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
!
     do ideg = 1,3 !...Calculating time, use only p1....
       unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
     enddo
!
!... Conservative variables for left state...
!
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!... Primitive variable for left state...
      rhom1 = 1.d0/rho1
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
      pres1 = max(eps, (gamma-1.d0)*(rhoe1-0.5d0*rho1*(uadv1*uadv1+vadv1*vadv1)))
      csou1 = sqrt(max(eps,gamma*pres1/rho1))
      vdon1 = uadv1*dwav1+vadv1*dwav2 
!
      veloc(iel) = veloc(iel) + wi*djaco*(abs(vdon1) + csou1) 

!    
   enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nbfac
!
!...Characteristic velocity from internal faces...
!
    do 260 ifa=nbfac+1,nafac !...(1)ifa=1,nafac
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
     do ideg = 1,3  !...Calculating time, use only p1....
       unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
       unkno2(1:nq) = unkno2(1:nq) + unkno(ideg,1:nq,ier)*br(ideg)
     enddo
!
!... Conservative variables for left state...
!
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!... Primitive variable for left state...
      rhom1 = 1.d0/rho1 
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
      pres1 = max(eps, (gamma-1.d0)*(rhoe1-0.5d0*rho1*(uadv1*uadv1+vadv1*vadv1)))
      csou1 = sqrt(max(eps,gamma*pres1/rho1))
      vdon1 = uadv1*dwav1 + vadv1*dwav2 
!
!... Conservative variables for right state...
!
      rho2  = unkno2(1)
      rhou2 = unkno2(2)
      rhov2 = unkno2(3)
      rhoe2 = unkno2(4)
!...Primitive variable for right state...
      rhom2 = 1.d0/rho2 
      uadv2 = rhou2*rhom2
      vadv2 = rhov2*rhom2
      pres2 = max(eps, (gamma-1.d0)*(rhoe2-0.5d0*rho2*(uadv2*uadv2+vadv2*vadv2)))
      csou2 = sqrt(max(eps,gamma*pres2/rho2))
      vdon2 = uadv2*dwav1 + vadv2*dwav2 
!    
      veloc(iel) = veloc(iel) + wi*djaco*max(abs(vdon1) + csou1, abs(vdon2) + csou2) 
      veloc(ier) = veloc(ier) + wi*djaco*max(abs(vdon1) + csou1, abs(vdon2) + csou2) 
!    
   enddo !...(4)ig=1,ngausf
!
260 enddo  !...(1)ifa=1,nafac

end subroutine getveloc_curv
!
!...Subroutine to calculate the characteristic velicity...
!
subroutine getveloc_curv2(geoel,bface,uchar,unkno,intfac,coord,veloc)
 use constant
 implicit none
 real*8,dimension(1:ngeel,1:nelem+nbfac),intent(in)::geoel 
 real*8,intent(in)::uchar(1:nq)
 real*8,dimension(1:ndegr,1:nq,1:nelem+nbfac),intent(in)::unkno
 real*8,dimension(1:ndimn,1:npoin)::coord 
 real*8,dimension(1:nelem), intent(out)::veloc
 integer*4,dimension(1:nifai,1:nafac)::intfac
 integer*4,dimension(1:nbfai,nbfac)::bface 
!
!...local array
 real*8::xpin(2, 2)
 real*8::xp(2, nptfa)
 real*8,dimension(1:nptfa)::shp, dshpr
 real*8::weigh(ngausf), posi(1,ngausf)
 real*8,dimension(1:nq)::unkno1,unkno2
 real*8::flux(1:nq)
 real*8::bl(6), br(6)
!...local real number
 integer::ifa,iel,ier,ie,ifb,k,m,ig,ishp,ideg,jdeg
 real*8::xcl,ycl,dxl,dyl,rnx,rny
 real*8::xcr,ycr,dxr,dyr 
 real*8:: r, djaco, dxdr, dydr 
 real*8::xg,yg,wi
 real*8::rho1,rhou1,rhov1,rhoe1,rho2,rhou2,rhov2,rhoe2
 real*8::vn,dwav1,dwav2
 real*8::rhom1,uadv1,vadv1,pres1, csou1, vdon1
 real*8::rhom2,uadv2,vadv2,pres2, csou2, vdon2
 real*8::eps
!
 data   eps / 1.0d-06 /
!
 call rutope(1, ngausf, posi, weigh)
!
! print*,'gauss',posi(1,ngausf)
!
!...zero out rhs...
!
    veloc = 0.d0
!
!...Characteristic velocity from boundary faces...
!
    do 250 ifa = 1,nbfac !...(1)ifa=1,nafac
!
     iel=intfac(1,ifa)
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
     xg = 0.d0
     yg = 0.d0
!
   do ishp = 1, nptfa
     xg = xg + shp(ishp)*xp(1, ishp) 
     yg = yg + shp(ishp)*xp(2, ishp) 
    enddo
!
!...Basis function...
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
!
     do ideg = 1,ndegr
       unkno1(1:nq) = unkno1(1:nq) + unkno(ideg,1:nq,iel)*bl(ideg)
     enddo
!
!... Conservative variables for left state...
!
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!... Primitive variable for left state...
      rhom1 = 1.d0/rho1
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
      pres1 = (gamma-1.d0)*(rhoe1-0.5d0*rho1*(uadv1*uadv1+vadv1*vadv1))   
      csou1 = sqrt(max(eps,gamma*pres1/rho1))
      vdon1 = uadv1*dwav1+vadv1*dwav2 
!
      veloc(iel) = veloc(iel) + wi*djaco*(abs(vdon1) + csou1) 

!    
   enddo !...(4)ig=1,ngausf
!
250 enddo  !...(1)ifa=1,nbfac
!
!...Characteristic velocity from internal faces...
!
    do 260 ifa=nbfac+1,nafac !...(1)ifa=1,nafac
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
!... Conservative variables for left state...
!
      rho1  = unkno1(1)
      rhou1 = unkno1(2)
      rhov1 = unkno1(3)
      rhoe1 = unkno1(4)
!... Primitive variable for left state...
      rhom1 = 1.d0/rho1 
      uadv1 = rhou1*rhom1
      vadv1 = rhov1*rhom1
      pres1 = (gamma-1.d0)*(rhoe1-0.5d0*rho1*(uadv1*uadv1+vadv1*vadv1))   
      csou1 = sqrt(max(eps,gamma*pres1/rho1))
      vdon1 = uadv1*dwav1 + vadv1*dwav2 
!
!... Conservative variables for right state...
!
      rho2  = unkno2(1)
      rhou2 = unkno2(2)
      rhov2 = unkno2(3)
      rhoe2 = unkno2(4)
!...Primitive variable for right state...
      rhom2 = 1.d0/rho2 
      uadv2 = rhou2*rhom2
      vadv2 = rhov2*rhom2
      pres2 = (gamma-1.d0)*(rhoe2-0.5d0*rho2*(uadv2*uadv2+vadv2*vadv2))   
      csou2 = sqrt(max(eps,gamma*pres2/rho2))
      vdon2 = uadv2*dwav1 + vadv2*dwav2 
!    
      veloc(iel) = veloc(iel) + wi*djaco*max(abs(vdon1) + csou1, abs(vdon2) + csou2) 
      veloc(ier) = veloc(ier) + wi*djaco*max(abs(vdon1) + csou1, abs(vdon2) + csou2) 
!    
   enddo !...(4)ig=1,ngausf
!
260 enddo  !...(1)ifa=1,nafac

end subroutine getveloc_curv2
!
!...Subroutine: Solver for DG...
!
subroutine getsolvdg(intfac,inpoel,iptri, ipqua, bface,coord,geofa,geoel,uchar,unkno,unold,ltelem,fsuel,itime,amatr)
   use constant
   implicit none
   integer*4::istag,itime,ie,iq,ig,id,iunk
   integer*4,dimension(1:nvtri,1:ntria)::inpoel
   integer*4,dimension(1:nvtri,1:ntria)::iptri
   integer*4,dimension(1:nvqua,1:nquad)::ipqua
   integer*4,dimension(1:nifai,1:nafac)::intfac
   real*8,dimension(1:ndimn,1:npoin)::coord
   real*8,dimension(1:ngefa,1:nafac)::geofa
   real*8,dimension(1:ngeel,1:nsize)::geoel
   real*8,dimension(1:nq)::uchar
   real*8,dimension(1:ncell)::deltt
   real*8,dimension(1:ndegr,1:nq,1:nsize)::unold
   real*8,dimension(1:ndegr,1:nq,1:nsize),intent(in)::unkno
   real*8,dimension(1:ndegr,1:nq,1:ncell)::rhsel
   real*8,intent(in)::amatr(nmatr,ncell)
   real*8::m(ndegr, ndegr)
   real*8:: unint(1:nq)
   integer*4,dimension(1:nbfai,1:nbfac)::bface
   integer*4,dimension(1:3,1:ncell)::ltelem, fsuel
   real*8::b(6)
   real*8::alfa,dx,dy,xg,yg ,xc,yc,p,q
   integer::ideg
!
!...find the advancing time step...
!
    call getdelt(deltt,geoel,bface,uchar,unkno,intfac,coord)   !local time step dt,unkno,geoel,cfl
!
    do 2500 istag=1,nstag !...(1)istag=1,nstag

      alfa=rkcoe(istag)
     !
     !...find the advancing time step...
     ! call getdelt(deltt,geoel,bface,uchar,unkno,intfac,coord)   !local time step dt,unkno,geoel,cfl
     ! print*,'mmmmm'
     !print*,'xxxxx1',deltt(1:3)
      call getrhsdg(uchar,bface,unold,rhsel,intfac,inpoel,iptri,ipqua,geofa,geoel,coord,ltelem,fsuel,amatr,itime)
     !
     !
     do ie =1, ncell
     do iq =1, nq  !...nq variables for hyperbolic DGM.     
        rhsel(1:ndegr,iq,ie)=deltt(ie)*rhsel(1:ndegr,iq,ie)
     enddo
     enddo
     ! print*,'xxxxxrhs',rhsel(1:ndegr,2,1)
!
!...for high-order dg...
!
     if(npoly.ge.1)then !...(2)npoly.ge.1
!
      do ie=1,ncell     !...(3)ie=1,nelem
!
        dx = geoel(4, ie)*0.5d0
        dy = geoel(5, ie)*0.5d0
!
        if(npoly==1)then
          m(1,1) = 1.d0/geoel(3,ie)
          m(1,2) = 0.d0
          m(1,3) = 0.d0

          m(2,1) = m(1,2)
          m(2,2) = amatr(1,ie)
          m(2,3) = amatr(2,ie)

          m(3,1) = m(1,3)
          m(3,2) = amatr(2,ie)
          m(3,3) = amatr(3,ie)

        elseif(npoly==2)then

          m(1,1) = 1.d0/geoel(3,ie)
          m(1,2) = 0.d0
          m(1,3) = 0.d0
          m(1,4) = 0.d0
          m(1,5) = 0.d0
          m(1,6) = 0.d0

          m(2,1) = m(1,2)
          m(2,2) = amatr(1,ie)
          m(2,3) = amatr(2,ie)
          m(2,4) = amatr(3,ie)
          m(2,5) = amatr(4,ie)
          m(2,6) = amatr(5,ie)
 
          m(3,1) = m(1,3)
          m(3,2) = m(2,3)
          m(3,3) = amatr(6,ie)
          m(3,4) = amatr(7,ie)
          m(3,5) = amatr(8,ie)
          m(3,6) = amatr(9,ie)

          m(4,1) = m(1,4)
          m(4,2) = m(2,4)
          m(4,3) = m(3,4)
          m(4,4) = amatr(10,ie)
          m(4,5) = amatr(11,ie)
          m(4,6) = amatr(12,ie)

          m(5,1) = m(1,5)
          m(5,2) = m(2,5)
          m(5,3) = m(3,5)
          m(5,4) = m(4,5)
          m(5,5) = amatr(13,ie)
          m(5,6) = amatr(14,ie)

          m(6,1) = m(1,6)
          m(6,2) = m(2,6)
          m(6,3) = m(3,6)
          m(6,4) = m(4,6)
          m(6,5) = m(5,6)
          m(6,6) = amatr(15,ie)
        endif
!
!...solve the mdegr independant varaible...
!
    do id =1,ndegr !...(4)id =1,ndegr
!
!...step 1
!
        unint = 0.d0
     do iunk = 1,ndegr
        unint(1:nq) = unint(1:nq) + m(id, iunk)*rhsel(iunk,1:nq,ie)
     enddo
!
!...step 2
!
        unold(id,1:nq,ie)= alfa*unkno(id,1:nq,ie) + (1.d0-alfa)*(unold(id,1:nq,ie) + unint(1:nq))
!
     enddo !...(4)id =1,ndegr
!  
    enddo !...(3)ie=1,nelem 
!
    endif !...(2)npoly.ge.1
!
!    print*,'xxxxx',unkno(1:3, 2, 1)
!
2500 enddo  !...(1)istag=1,nstag
 !

!
!...Output the solution field after calculation...
!
!
    end subroutine getsolvdg
!
!...Subroutine: Get inverse matrix....
!
subroutine inverse(a,c,n,det)

   implicit none
   integer n
   real*8 a(n,n), c(n,n)
   real*8 det

   ! step 0: initialization for matrices L and U and b
   ! Fortran 90/95 aloows such operations on matrices
   det = a(3, 1)*a(1, 2)*a(2, 3)-a(3, 1)*a(1, 3)*a(2, 2)-a(2, 1)*a(1, 2)*a(3, 3) &
         -a(2, 3)*a(3, 2)*a(1, 1)+a(2, 1)*a(1, 3)*a(3, 2)+a(2, 2)*a(3, 3)*a(1, 1)

   c(1,1) = a(2, 2)*a(3, 3)-a(2, 3)*a(3, 2)
   c(2,1) = -a(2, 1)*a(3, 3)+a(2, 3)*a(3, 1)
   c(3,1) = a(2, 1)*a(3, 2)-a(2, 2)*a(3, 1)

   c(1,2) = -a(1, 2)*a(3, 3)+a(1, 3)*a(3, 2)
   c(2,2) = a(1, 1)*a(3, 3)-a(1, 3)*a(3, 1)
   c(3,2) = -a(1, 1)*a(3, 2)+a(1, 2)*a(3, 1)

   c(1,3) = a(1, 2)*a(2, 3)-a(1, 3)*a(2, 2)
   c(2,3) = -a(1, 1)*a(2, 3)+a(1, 3)*a(2, 1)
   c(3,3) = a(1, 1)*a(2, 2)-a(1, 2)*a(2, 1)

   c = c/det

   if(n.ne.3) print*,'inverse should not be used for  n ne 3'
end
!
!...Subroutine: Gauss quadrature for triangle...
!
subroutine rutope(ndimn, ngaus, posgp, weigp)
!
!     This routine sets up the integration constants of open rules for
!     triangles and tetrahedra
!
!     NDIMN = 1             NDIME = 2             NDIME = 3
!
!     NGAUS  EXACT POL.    NGAUS  EXACT POL.     NGAUS  EXACT POL. 
!     -----  ---------     -----  ----------     -----  ----------
!         1         p1         1       p1            1       p1
!         2         p3         3       p2            4       p2
!         3         p5         4       p3            5       p3
!                              6       p4           11       p4
!                              7       p5           14       p5
!                             13       p9
!
      implicit real*8 (a-h,o-z)
      real*8 posgp(ndimn,*), weigp(*)
!
      istop = 0
!
      goto (100, 200, 300) ndimn
!
!...  Line integral (the same as for brick elements)
!
 100  continue
        if(ngaus .eq. 1) then
          posgp(1,1) = 0.d0
          weigp(  1) = 2.d0
        else if(ngaus .eq. 2) then
          posgp(1,1) =-0.577350269189626d0
          posgp(1,2) = 0.577350269189626d0
          weigp(  1) = 1.d0
          weigp(  2) = 1.d0
        else if(ngaus .eq. 3) then
          posgp(1,1) =-sqrt(.6d0)!-0.774596669241483d0 
          posgp(1,2) = 0.d0
          posgp(1,3) = sqrt(.6d0)!0.774596669241483d0
          weigp(  1) = 5.d0/9.d0!0.555555555555556d0
          weigp(  2) = 8.d0/9.d0!0.888888888888889d0
          weigp(  3) = 5.d0/9.d0!0.555555555555556d0
        else if(ngaus .eq. 4) then
          posgp(1,1) =-0.8611363115940526d0 
          posgp(1,2) =-0.3399810435848563d0
          posgp(1,3) = 0.3399810435848563d0
          posgp(1,4) = 0.8611363115940526d0
          weigp(  1) = 0.3478548451374538d0 
          weigp(  2) = 0.6521451548625461d0
          weigp(  3) = 0.6521451548625461d0
          weigp(  4) = 0.3478548451374538d0 
        else
          istop=1
        end if
      goto 999
!
!...  Area integral (triangles)
!
 200  continue
        if(ngaus .eq. 1) then
          posgp(1,1) = 1.d0/3.d0
          posgp(2,1) = 1.d0/3.d0
          weigp(  1) = 1.d0
        else if(ngaus .eq. 3) then
          posgp(1,1) = 2.d0/3.d0
          posgp(2,1) = 1.d0/6.d0
          posgp(1,2) = 1.d0/6.d0
          posgp(2,2) = 2.d0/3.d0
          posgp(1,3) = 1.d0/6.d0
          posgp(2,3) = 1.d0/6.d0
          weigp(  1) = 1.d0/3.d0
          weigp(  2) = 1.d0/3.d0
          weigp(  3) = 1.d0/3.d0
        else if(ngaus .eq. 4) then
          posgp(1,1) = 1.0/3.0
          posgp(2,1) = 1.0/3.0
          posgp(1,2) = 1.0/5.0
          posgp(2,2) = 1.0/5.0
          posgp(1,3) = 3.0/5.0
          posgp(2,3) = 1.0/5.0
          posgp(1,4) = 1.0/5.0
          posgp(2,4) = 3.0/5.0
          weigp(  1) =-27.0/48.0
          weigp(  2) = 25.0/48.0
          weigp(  3) = 25.0/48.0
          weigp(  4) = 25.0/48.0 
        else if(ngaus .eq. 6) then
          ex1 = 0.816847572980459d0 
          et1 = 0.091576213509771d0
          ez1 = 0.091576213509771d0
          ex2 = 0.108103018168070d0
          et2 = 0.445948490915965d0
          ez2 = 0.445948490915965d0
          posgp(1,1) = ex1
          posgp(2,1) = et1
          posgp(1,2) = et1
          posgp(2,2) = ez1
          posgp(1,3) = ez1
          posgp(2,3) = ex1
          posgp(1,4) = ex2
          posgp(2,4) = et2
          posgp(1,5) = et2
          posgp(2,5) = ez2
          posgp(1,6) = ez2
          posgp(2,6) = ex2
!          a = 0.054975870996713638d0*2.d0
!          b = 0.1116907969117165d0*2.d0   
!
          a = 0.109951743655322d0 
          b = 0.223381589678011d0  
  
          weigp(1)   = a
          weigp(2)   = a
          weigp(3)   = a
          weigp(4)   = b
          weigp(5)   = b
          weigp(6)   = b
        else if(ngaus .eq. 7) then
          posgp(1,1) = 0.333333333333333d0 
          posgp(2,1) = 0.333333333333333d0
          posgp(1,2) = 0.059715871789770d0
          posgp(2,2) = 0.470142064105115d0
          posgp(1,3) = 0.470142064105115d0
          posgp(2,3) = 0.059715871789770d0
          posgp(1,4) = 0.470142064105115d0
          posgp(2,4) = 0.470142064105115d0
          posgp(1,5) = 0.797426985353087d0
          posgp(2,5) = 0.101286507323456d0
          posgp(1,6) = 0.101286507323456d0
          posgp(2,6) = 0.797426985353087d0
          posgp(1,7) = 0.101286507323456d0
          posgp(2,7) = 0.101286507323456d0
          weigp(  1) = 0.225000000000000d0
          weigp(  2) = 0.132394152788506d0
          weigp(  3) = 0.132394152788506d0
          weigp(  4) = 0.132394152788506d0
          weigp(  5) = 0.125939180544827d0
          weigp(  6) = 0.125939180544827d0
          weigp(  7) = 0.125939180544827d0
        else if(ngaus .eq. 13) then
          a = 0.333333333333333d0 
          b = 0.479308067841920d0
          c = 0.869739794195568d0
          d = 0.638444188569810d0
          e = 0.260345966079040d0
          f = 0.065130102902216d0
          g = 0.312865496004874d0
          h = 0.048690315425316d0
          w1=-0.149570044467682d0
          w2= 0.175615257433208d0
          w3= 0.053347235608838d0
          w4= 0.077113760890257d0
          posgp(1, 1)= a
          posgp(2, 1)= a         
          posgp(1, 2)= e
          posgp(2, 2)= e
          posgp(1, 3)= b
          posgp(2, 3)= e        
          posgp(1, 4)= e
          posgp(2, 4)= b        
          posgp(1, 5)= f
          posgp(2, 5)= f        
          posgp(1, 6)= c
          posgp(2, 6)= f        
          posgp(1, 7)= f
          posgp(2, 7)= c        
          posgp(1, 8)= d
          posgp(2, 8)= g        
          posgp(1, 9)= d
          posgp(2, 9)= h        
          posgp(1,10)= g
          posgp(2,10)= d        
          posgp(1,11)= g
          posgp(2,11)= h        
          posgp(1,12)= h
          posgp(2,12)= d        
          posgp(1,13)= h
          posgp(2,13)= g
          weigp( 1) = w1
          weigp( 2) = w2
          weigp( 3) = w2
          weigp( 4) = w2
          weigp( 5) = w3
          weigp( 6) = w3
          weigp( 7) = w3
          weigp( 8) = w4
          weigp( 9) = w4
          weigp(10) = w4
          weigp(11) = w4
          weigp(12) = w4
          weigp(13) = w4
        else
          istop=1
        end if
      goto 999
!
!...  Volume integral ( tetrahedra )
!
 300  continue
        if(ngaus .eq. 1) then
          posgp(1,1)= 1.0/4.0
          posgp(2,1)= 1.0/4.0
          posgp(3,1)= 1.0/4.0
          weigp(1)  = 1.0
        else if(ngaus .eq. 4) then
          a = 0.5854101966249685
          b = 0.1381966011250105
          posgp(1,1) = b
          posgp(2,1) = b
          posgp(3,1) = b
          posgp(1,2) = a
          posgp(2,2) = b
          posgp(3,2) = b
          posgp(1,3) = b
          posgp(2,3) = a
          posgp(3,3) = b
          posgp(1,4) = b
          posgp(2,4) = b
          posgp(3,4) = a
          weigp(  1) = 0.25   
          weigp(  2) = 0.25   
          weigp(  3) = 0.25   
          weigp(  4) = 0.25   
        else if(ngaus .eq. 5) then
          posgp(1,1)= 1.0/4.0
          posgp(2,1)= 1.0/4.0
          posgp(3,1)= 1.0/4.0
          posgp(1,2)= 1.0/6.0
          posgp(2,2)= 1.0/6.0
          posgp(3,2)= 1.0/6.0
          posgp(1,3)= 1.0/2.0
          posgp(2,3)= 1.0/6.0
          posgp(3,3)= 1.0/6.0
          posgp(1,4)= 1.0/6.0
          posgp(2,4)= 1.0/2.0
          posgp(3,4)= 1.0/6.0
          posgp(1,5)= 1.0/6.0
          posgp(2,5)= 1.0/6.0
          posgp(3,5)= 1.0/2.0
          weigp(  1)=-12.0/15.0
          weigp(  2)= 9.0/20.0
          weigp(  3)= 9.0/20.0
          weigp(  4)= 9.0/20.0
          weigp(  5)= 9.0/20.0
        else if(ngaus .eq. 11) then
          a = 0.3994035761667992
          b = 0.1005964238332008
          c = 343.0/7500.0/6.0
          d = 56.0/375.0/6.0
          posgp(1,1) = 1.0/4.0
          posgp(2,1) = 1.0/4.0
          posgp(3,1) = 1.0/4.0
          posgp(1,2) = 11.0/14.0
          posgp(2,2) = 1.0/14.0
          posgp(3,2) = 1.0/14.0
          posgp(1,3) = 1.0/14.0
          posgp(2,3) = 11.0/14.0
          posgp(3,3) = 1.0/14.0
          posgp(1,4) = 1.0/14.0
          posgp(2,4) = 1.0/14.0
          posgp(3,4) = 11.0/14.0
          posgp(1,5) = 1.0/14.0
          posgp(2,5) = 1.0/14.0
          posgp(3,5) = 1.0/14.0
          posgp(1,6) = a
          posgp(2,6) = a
          posgp(3,6) = b
          posgp(1,7) = a
          posgp(2,7) = b
          posgp(3,7) = a
          posgp(1,8) = a
          posgp(2,8) = b
          posgp(3,8) = b
          posgp(1,9) = b
          posgp(2,9) = a
          posgp(3,9) = a
          posgp(1,10)= b
          posgp(2,10)= a
          posgp(3,10)= b
          posgp(1,11)= b
          posgp(2,11)= b
          posgp(3,11)= a
          weigp(1)   =-148.0/1875.0
          weigp(2)   = c*6.0
          weigp(3)   = c*6.0
          weigp(4)   = c*6.0
          weigp(5)   = c*6.0
          weigp(6)   = d*6.0
          weigp(7)   = d*6.0
          weigp(8)   = d*6.0
          weigp(9)   = d*6.0
          weigp(10)  = d*6.0
          weigp(11)  = d*6.0
        else if(ngaus .eq. 14) then
          a = 0.0673422422100983
          b = 0.3108859192633005
          c = 0.7217942490673264
          d = 0.0927352503108912
          e = 0.4544962958743506
          f = 0.0455037041256494
          p = 0.1126879257180162
          q = 0.0734930431163619
          r = 0.0425460207770812
          posgp(1,1) = a
          posgp(2,1) = b
          posgp(3,1) = b
          posgp(1,2) = b
          posgp(2,2) = a
          posgp(3,2) = b
          posgp(1,3) = b
          posgp(2,3) = b
          posgp(3,3) = a
          posgp(1,4) = b
          posgp(2,4) = b
          posgp(3,4) = b
          posgp(1,5) = c
          posgp(2,5) = d
          posgp(3,5) = d
          posgp(1,6) = d
          posgp(2,6) = c
          posgp(3,6) = d
          posgp(1,7) = d
          posgp(2,7) = d
          posgp(3,7) = c
          posgp(1,8) = d
          posgp(2,8) = d
          posgp(3,8) = d
          posgp(1,9) = e
          posgp(2,9) = e
          posgp(3,9) = f
          posgp(1,10)= e
          posgp(2,10)= f
          posgp(3,10)= e
          posgp(1,11)= e
          posgp(2,11)= f
          posgp(3,11)= f
          posgp(1,12)= f
          posgp(2,12)= e
          posgp(3,12)= e
          posgp(1,13)= f
          posgp(2,13)= e
          posgp(3,13)= f
          posgp(1,14)= f
          posgp(2,14)= f
          posgp(3,14)= e
          weigp(1)   = p
          weigp(2)   = p
          weigp(3)   = p
          weigp(4)   = p
          weigp(5)   = q
          weigp(6)   = q
          weigp(7)   = q
          weigp(8)   = q
          weigp(9)   = r
          weigp(10)  = r
          weigp(11)  = r
          weigp(12)  = r
          weigp(13)  = r
          weigp(14)  = r
        else
          istop=1
        end if
      goto 999
!
 999  continue
!      
!...  Errors
!
      if(istop .eq. 1) then
      write(*,*)'RUTOPE: NO  AVAILABLE QUADRATURE'
      write(*,*)'ndimn = ', ndimn
      write(*,*)'ngaus = ', ngaus
      stop
      endif
!
      return
end subroutine rutope
!
!...Subroutine for gauss quadrature for quad and hexa(Gauss-Lobatto)...
!
subroutine ruqope_lobatto(ndime,ngaus,posgp,weigp)
!
! This routine sets up the integration constants of open
! integration rules for brick elements:
!
!         NDIME = 1             NDIME = 2             NDIME = 3
!
!     NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL.
!     -----  ----------     -----  ----------     -----  ----------
!       1      q1           1 x 1     q1          1x1x1     q1
!	2      q3           2 x 2     q3          2x2x2     q3
!	3      q5           3 x 3     q5          3x3x3     q5
!	4      q7           4 x 4     q7          4x4x4     q7
!
!
implicit real*8(a-h,o-z)
real*8 posgl(6), posgp(ndime,ngaus), weigl(6), weigp(ngaus)
!
nlocs=0       !This is needed if NGAUS=0

if(ndime.eq.1) then
nlocs=ngaus
else if(ndime.eq.2) then
if (ngaus .eq. 1) then
nlocs=1
else if (ngaus .eq. 4) then
nlocs=2
else if (ngaus .eq. 9) then
nlocs=3
else if (ngaus .eq. 16) then
nlocs=4
else if (ngaus .eq. 25) then
nlocs=5
else if (ngaus .eq. 36) then
nlocs=6
end if
else
if (ngaus .eq. 1) then
nlocs=1
else if (ngaus .eq. 8) then
nlocs=2
else if (ngaus .eq. 27) then
nlocs=3
else if (ngaus .eq. 64) then
nlocs=4
else if (ngaus .eq. 125) then
nlocs=5
else if (ngaus .eq. 216) then
nlocs=6
end if
end if
!
if(nlocs.eq.1) then
!
print*,'nlocs.eq.1 does not exist for Lobatto'
stop
!
else if(nlocs.eq.2) then
posgl(1)=-1.d0
posgl(2)= 1.d0
weigl(1)= 1.0d0
weigl(2)= 1.0d0
else if(nlocs.eq.3) then
!
posgl(1) =-1.d0
posgl(2) = 0.d0
posgl(3) = 1.d0
weigl(  1) = 1.d0/3.d0
weigl(  2) = 4.d0/3.d0
weigl(  3) = 1.d0/3.d0
!print*,'nlocs.eq.3 does not exist for Lobatto'
!stop
!
else if(nlocs.eq.4)  then
posgl(1) =-1.d0
posgl(2) =-0.4472135954999579392818d0
posgl(3) = 0.4472135954999579392818d0
posgl(4) = 1.d0
weigl(  1) = 0.166666666666666666666d0
weigl(  2) = 0.833333333333333333333d0
weigl(  3) = 0.833333333333333333333d0
weigl(  4) = 0.166666666666666666666d0
else if(nlocs.eq.5)  then
!
posgl(1) =-1.d0
posgl(2) =-0.6546536707079771437983d0
posgl(3) = 0.d0
posgl(4) = 0.6546536707079771437983d0
posgl(5) = 1.d0

weigl(  1) = 0.1d0
weigl(  2) = 0.544444444444444444444d0
weigl(  3) = 0.7111111111111111111111d0
weigl(  4) = 0.544444444444444444444d0
weigl(  5) = 0.1d0
!
else if(nlocs.eq.6)  then
posgl(1) =-1.d0
posgl(2) =-0.765055323929464692851d0
posgl(3) =-0.2852315164806450963142d0
posgl(4) = 0.2852315164806450963142d0
posgl(5) = 0.765055323929464692851d0
posgl(6) = 1.d0
weigl(  1) = 0.06666666666666666666667d0
weigl(  2) = 0.378474956297846980317d0
weigl(  3) = 0.55485837703548635302d0
weigl(  4) = 0.55485837703548635302d0
weigl(  5) = 0.378474956297846980317d0
weigl(  6) = 0.06666666666666666666667d0
!
end if
!
if(ndime.eq.1) then
igaus=0
do ilocs=1,nlocs
igaus=igaus+1
weigp(  igaus)=weigl(ilocs)
posgp(1,igaus)=posgl(ilocs)
end do
else if(ndime.eq.2) then
igaus=0
do ilocs=1,nlocs
do jlocs=1,nlocs
igaus=igaus+1
weigp(  igaus)=weigl(ilocs)*weigl(jlocs)
posgp(1,igaus)=posgl(ilocs)
posgp(2,igaus)=posgl(jlocs)
end do
end do
else if(ndime.eq.3) then
igaus=0
do ilocs=1,nlocs
do jlocs=1,nlocs
do klocs=1,nlocs
igaus=igaus+1
weigp(  igaus)=weigl(ilocs)*weigl(jlocs)*weigl(klocs)
posgp(1,igaus)=posgl(ilocs)
posgp(2,igaus)=posgl(jlocs)
posgp(3,igaus)=posgl(klocs)
end do
end do
end do
end if
return
!
end subroutine ruqope_lobatto
!
!...Subroutine for gauss quadrature for quad and hexa...
!
subroutine ruqope(ndime,ngaus,posgp,weigp)
!
! This routine sets up the integration constants of open
! integration rules for brick elements:
!
!         NDIME = 1             NDIME = 2             NDIME = 3
!
!     NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL.
!     -----  ----------     -----  ----------     -----  ----------
!       1      q1           1 x 1     q1          1x1x1     q1
!	2      q3           2 x 2     q3          2x2x2     q3
!	3      q5           3 x 3     q5          3x3x3     q5
!	4      q7           4 x 4     q7          4x4x4     q7
!
!
implicit real*8(a-h,o-z)
real*8 posgl(6), posgp(ndime,ngaus), weigl(6), weigp(ngaus)
!
nlocs=0       !This is needed if NGAUS=0

if(ndime.eq.1) then
nlocs=ngaus
else if(ndime.eq.2) then
if (ngaus .eq. 1) then
nlocs=1
else if (ngaus .eq. 4) then
nlocs=2
else if (ngaus .eq. 9) then
nlocs=3
else if (ngaus .eq. 16) then
nlocs=4
else if (ngaus .eq. 25) then
nlocs=5
else if (ngaus .eq. 36) then
nlocs=6
end if
else
if (ngaus .eq. 1) then
nlocs=1
else if (ngaus .eq. 8) then
nlocs=2
else if (ngaus .eq. 27) then
nlocs=3
else if (ngaus .eq. 64) then
nlocs=4
else if (ngaus .eq. 125) then
nlocs=5
else if (ngaus .eq. 216) then
nlocs=6
end if
end if
!
if(nlocs.eq.1) then
posgl(1)=0.0d0
weigl(1)=2.0d0
else if(nlocs.eq.2) then
posgl(1)=-0.5773502691896257645091488d0
posgl(2)= 0.5773502691896257645091488d0
weigl(1)= 1.0d0
weigl(2)= 1.0d0
else if(nlocs.eq.3) then
posgl(1)=-0.7745966692414833770358531d0
posgl(2)= 0.0d0
posgl(3)= 0.7745966692414833770358531d0
weigl(1)= 0.5555555555555555555555556d0
weigl(2)= 0.8888888888888888888888889d0
weigl(3)= 0.5555555555555555555555556d0
else if(nlocs.eq.4)  then
posgl(1)=-0.8611363115940525752239465d0
posgl(2)=-0.3399810435848562648026658d0
posgl(3)= 0.3399810435848562648026658d0
posgl(4)= 0.8611363115940525752239465d0
weigl(1)= 0.3478548451374538573730639d0
weigl(2)= 0.6521451548625461426269361d0
weigl(3)= 0.6521451548625461426269361d0
weigl(4)= 0.3478548451374538573730639d0
else if(nlocs.eq.5)  then
posgl(3)= 0.d0
posgl(4)= 0.5384693101056830910363144d0
posgl(2)=-0.5384693101056830910363144d0
posgl(5)= 0.9061798459386639927976269d0
posgl(1)=-0.9061798459386639927976269d0
weigl(3)=  0.5688888888888888888888889d0
weigl(4)=  0.4786286704993664680412915d0
weigl(2)=  0.4786286704993664680412915d0
weigl(5)=  0.2369268850561890875142640d0
weigl(1)=  0.2369268850561890875142640d0

else if(nlocs.eq.6)  then
posgl(1)= 0.2386191860831969086305017d0
posgl(2)= -0.2386191860831969086305017d0
posgl(3)= 0.6612093864662645136613996d0
posgl(4)= -0.6612093864662645136613996d0
posgl(5)=  -0.9324695142031520278123016d0
posgl(6)=  0.9324695142031520278123016d0
weigl(1)=  0.4679139345726910d0!473898703d0
weigl(2)=  0.4679139345726910d0!473898703d0
weigl(3)=  0.3607615730481386d0!075698335d0
weigl(4)=  0.3607615730481386d0!075698335d0
weigl(5)=  0.1713244923791704d0!3450402961d0
weigl(6)=  0.1713244923791704d0!3450402961d0
!
end if
!
if(ndime.eq.1) then
igaus=0
do ilocs=1,nlocs
igaus=igaus+1
weigp(  igaus)=weigl(ilocs)
posgp(1,igaus)=posgl(ilocs)
end do
else if(ndime.eq.2) then
igaus=0
do ilocs=1,nlocs
do jlocs=1,nlocs
igaus=igaus+1
weigp(  igaus)=weigl(ilocs)*weigl(jlocs)
posgp(1,igaus)=posgl(ilocs)
posgp(2,igaus)=posgl(jlocs)
end do
end do
else if(ndime.eq.3) then
igaus=0
do ilocs=1,nlocs
do jlocs=1,nlocs
do klocs=1,nlocs
igaus=igaus+1
weigp(  igaus)=weigl(ilocs)*weigl(jlocs)*weigl(klocs)
posgp(1,igaus)=posgl(ilocs)
posgp(2,igaus)=posgl(jlocs)
posgp(3,igaus)=posgl(klocs)
end do
end do
end do
end if
return
!
end subroutine ruqope
!
!=====================================================================
!
!  this sub. computes the inverse of a matrix of arbitrary size
!  by using LU decomposition
!
!=====================================================================
!
subroutine getinvmat(n, a, x, b )
 !
 implicit   real*8 (a-h,o-z)
 !
 real*8  a(n,n), x(n,n), b(n)
 !
 ! inverse of a matrix
 !
 call ludeco(n, a)
 !
 do i = 1, n
   !
   do j = 1, n
     !
     b(j) = 0.0d0
     !
   enddo
   !
   b(i) = 1.0d0
   !
   call lureso(n, a, b, x(1,i) )
   !
 enddo
 !
 return
end
!
!=====================================================================
!
!  algorithm of Lower-Upper decomposition for a square matrix
!
!=====================================================================
!
subroutine ludeco(n, a)
 !
 implicit   real*8 (a-h,o-z)
 !
 real*8  a(n,n)
 !
 ! LU decomposition by point
 !
 do k = 1, n
   !
   ! decompose the diagonal terms
   !
   do m = 1, k-1
     !
     a(k,k) = a(k,k) - a(k,m)*a(m,k)
     !
   enddo
   !
   adicv = 1.0d0/a(k,k)
   !
   ! decompose the non diagonal terms
   !
   do i = k+1, n
     !
     do m = 1, k-1
       !
       a(k,i) = a(k,i) - a(k,m)*a(m,i)
       a(i,k) = a(i,k) - a(i,m)*a(m,k)
       !
     enddo
     !
     a(i,k) = a(i,k)*adicv
     !
   enddo
   !
 enddo
 !
 return
end
!
!=====================================================================
!
!  algorithm of Lower-Upper resolution for a square matrix
!
!=====================================================================
!
subroutine lureso(n, a, b, x )
 !
 implicit   real*8 (a-h,o-z)
 !
 real*8  a(n,n), b(n), x(n)
 !
 ! LU resolution
 !
 do i = 1, n
   !
   c = 0.0d0
   !
   do j = 1, i-1
     !
     c = c + a(i,j)*x(j)
     !
   enddo
   !
   x(i) = b(i) - c
   !
 enddo
 !
 do i = n, 1, -1
   !
   c = 0.0d0
   !
   do j = i+1, n
     !
     c = c + a(i,j)*x(j)
     !
   enddo
   !
   x(i) = (x(i) - c )/a(i,i)
   !
 enddo
 !
 return
end
