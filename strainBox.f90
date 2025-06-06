program strain_box
use hdf5
!use, intrinsic :: ISO_FORTRAN_ENV
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'
type ragged_array
  !2D array containing vectors of different lengths
  real,allocatable::v(:)
  character(len=12) :: indexName
end type ragged_array
character(len=200) :: filename
real :: modnor
integer,parameter :: maxMom=2
integer,parameter :: nxt=256, nyt=256, nzt=256, boxWid=nint(nzt/10.)
integer, dimension(3), parameter :: nt = (/nxt,nyt,nzt/)
real, parameter :: pi=3.14159265358979
real,parameter :: lx=2*pi, ly=2*pi, lz=2*pi
real, dimension(3), parameter :: l = (/lx,ly,lz/)
real,parameter :: dx=lx/nxt, dy=ly/nxt, dz=lz/nxt
real, dimension(nxt,nyt,nzt) :: kur,u,v,w,phase,dxxPhase
real, dimension(nxt,nyt,nzt) :: chemPot,dxChemPot,dyChemPot,dzChemPot,surfFX,surfFY,surfFZ
real, dimension(3,nxt,nyt,nzt) :: nor, vel
!use int8 for non shared arrays to save memory
integer, dimension(nxt,nyt,nzt) :: s_drop
integer, dimension(0:1,0:1,0:1) :: paint, neigh
integer, dimension(3)  :: last0,last1,pos
integer :: i,j,k,ip,jp,kp,iq,jq,kq,iShifted,mom,ii,jj,ky,kz,im,jm,km
integer :: cols,paintIt,faceOnCorner,genus,onInt,error,intR,nVels,ios
integer :: statU,posiU,veloU,MoInU,topoU,straU,dVelU,listU,specU,forcU,vortU
real, dimension(maxMom,3) :: dropPos, dropVel
real, dimension(3) :: MoIEiVals, dVeldx, strainEiVals, vort
real, dimension(3,3) :: MoI, dVeldxBox, strainBox
real :: diag, deformation, dropArea, dA, Cn, r, time, work(8), We, res
real :: maxNor, kurMean, kurStdDev, kurInv, kurInvSq, invSize, QInva, RInva
type(ragged_array) :: hist(3) !histogram of drop mass in x, y and z directions
!integer(kind=int64) :: dropSize, faces, edges, vertices
integer :: dropSize, faces, edges, vertices
character(len=200) :: fileEnd, fmtstr
integer(hid_t) :: file_id, dset_id
integer(hsize_t) :: dims(3)=(/nxt,nyt,nzt/),  dims1d(1)=(/1/) 
complex, dimension(nxt/2+1,nyt,nzt) :: cHat,chemPotHat,dxChemPotHat,dyChemPotHat,dzChemPotHat
complex, dimension(nxt/2+1,nyt,nzt) :: surfFXHat,surfFYHat,surfFZHat,uHat,vHat,wHat
real, dimension(0:nzt/2) :: FSpec,ESpec
type(C_PTR)  :: plan, plan_inverse
write(*,'(1x,a)') 'starting number of drops calculation                       '
do ii=1,3
  allocate(hist(ii)%v(nt(ii)))
enddo
hist(1)%indexName='1'
hist(2)%indexName='2'
hist(3)%indexName='3'
fileEnd='.txt'
open(newunit=statU,file='./output/statDrops'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
open(newunit=posiU,file='./output/posiDrops'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
open(newunit=veloU,file='./output/veloDrops'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
open(newunit=MoInU,file='./output/MoInDrops'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
open(newunit=dVelU,file='./output/dVeldxBox'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
open(newunit=straU,file='./output/strainBox'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
open(newunit=topoU,file='./output/topoDrops'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
open(newunit=specU,file='./output/ESpec'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
open(newunit=forcU,file='./output/FSpec'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
open(newunit=vortU,file='./output/vortBox'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
!call system("ls /home/alberto.velamartin/drop_time/we_05/run_break_197/field*.h5 > file_list.txt")
!call system("ls /home/alberto.velamartin/drop_time/we_05/run_break_197/field.160.h5 > file_list.txt")
call system("ls /home/alberto.velamartin/drop_time/we_05/run_break_197/field.*90.h5 > file_list.txt")
!Open the generated file list
open(newunit=listU, file="file_list.txt", status="old", action="read")
!Loop through file names
do
  read(listU, '(A)', iostat=ios) filename
  if (ios /= 0) exit
  print *, trim(filename)
  ! Open the file (read-only)
  call h5open_f(error)
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
    call h5dopen_f(file_id, "Cn", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, "We", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, We, dims1d, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, "c", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, phase, dims, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, "res", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, res, dims1d, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, "time", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, time, dims1d, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, "u", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, u, dims, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, "v", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, v, dims, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, "w", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, w, dims, error)
    call h5dclose_f(dset_id, error)
  call h5fclose_f(file_id, error)
  !change box size from 256 to 2 pi
  u=u*lx/nxt
  v=v*ly/nyt
  w=w*lz/nxt
  vel(1,:,:,:)=u
  vel(2,:,:,:)=v
  vel(3,:,:,:)=w
  !write(*,*)'E',0.5*sum(vel**2)/nzt**3
  write(*,*)'rmsVel',sqrt(sum(vel**2)/nzt**3/3)
  !write(*,*)'maxVelE',1.5*maxval(vel)**2
  write(*,*)'maxVel',maxval(vel)
  plan        =fftwf_plan_dft_r2c_3d(nzt, nyt, nxt, phase, cHat, FFTW_ESTIMATE)
  plan_inverse=fftwf_plan_dft_c2r_3d(nzt, nyt, nxt, cHat, phase, FFTW_ESTIMATE)
  call fftwf_execute_dft_r2c(plan, phase, cHat)
  do k=1,nzt
    kz=k-1
    if (kz.gt.nzt/2) kz=kz-nzt
    do j=1,nyt
      ky=j-1
      if (ky.gt.nyt/2) ky=ky-nyt
      do i=1,nxt/2+1
        cHat(i,j,k) = (kz**2 + ky**2 + (i-1)**2) * cHat(i,j,k)
      enddo
    enddo
  enddo
  call fftwf_execute_dft_c2r(plan_inverse, cHat, dxxPhase)
  chemPot = 1/Cn * (phase*phase - 1)*phase - Cn*dxxPhase
  call fftwf_execute_dft_r2c(plan, chemPot, chemPotHat)
  do k=1,nzt
    kz=k-1
    if (kz.gt.nzt/2) kz=kz-nzt
    do j=1,nyt
      ky=j-1
      if (ky.gt.nyt/2) ky=ky-nyt
      do i=1,nxt/2+1
        dxChemPotHat(i,j,k) =(i-1)* chemPotHat(i,j,k)
        dyChemPotHat(i,j,k) = ky * chemPotHat(i,j,k)
        dzChemPotHat(i,j,k) = kz * chemPotHat(i,j,k)
      enddo
    enddo
  enddo
  call fftwf_execute_dft_c2r(plan_inverse, dxChemPotHat, dxChemPot)
  call fftwf_execute_dft_c2r(plan_inverse, dyChemPotHat, dyChemPot)
  call fftwf_execute_dft_c2r(plan_inverse, dzChemPotHat, dzChemPot)
  do k=1,nzt
    do j=1,nyt
      do i=1,nxt
        surfFX(i,j,k) = phase(i,j,k) * dxChemPot(i,j,k)
        surfFY(i,j,k) = phase(i,j,k) * dyChemPot(i,j,k)
        surfFZ(i,j,k) = phase(i,j,k) * dzChemPot(i,j,k)
      enddo
    enddo
  enddo
  call fftwf_execute_dft_r2c(plan, surfFX, surfFXHat)
  call fftwf_execute_dft_r2c(plan, surfFY, surfFYHat)
  call fftwf_execute_dft_r2c(plan, surfFZ, surfFZHat)
  call fftwf_execute_dft_r2c(plan, u, uHat)
  call fftwf_execute_dft_r2c(plan, v, vHat)
  call fftwf_execute_dft_r2c(plan, w, wHat)
  ESpec=0.0
  FSpec=0.0
  do k=1,nzt
    kz=k-1
    if (kz.gt.nzt/2) kz=kz-nzt
    do j=1,nyt
      ky=j-1
      if (ky.gt.nyt/2) ky=ky-nyt
      do i=1,nxt/2+1
        r = sqrt( kz**2 + ky**2 + (i-1)**2 + 0.0)
        !Find the bin number for this radius. Bins have width 1.0/pointsPerWvNum.
        intR = int(r+0.5) 
        if(intR.le.nzt/2) then
          !Factor of half is cancelled when taking the real part of velocity*force
          ESpec(intR) = ESpec(intR) + 1.0/nxt/nyt/nzt* & 
                                    ( abs(uHat(i,j,k))**2 + &
                                      abs(vHat(i,j,k))**2 + &
                                      abs(wHat(i,j,k))**2 )
          FSpec(intR) = FSpec(intR) + 1.0/nxt/nyt/nzt* & 
                                    ( real( conjg(uHat(i,j,k)) * surfFXHat(i,j,k) )+ &
                                      real( conjg(vHat(i,j,k)) * surfFYHat(i,j,k) )+ &
                                      real( conjg(wHat(i,j,k)) * surfFZHat(i,j,k) ))
        endif
      enddo
    enddo
  enddo
  call fftw_destroy_plan(plan)
  call fftw_destroy_plan(plan_inverse)
  call fftw_cleanup()
  dropSize=0
  do k=1,nzt
    do j=1,nyt
      do i=1,nxt
        if(phase(i,j,k).ge.0.9)then
          s_drop(i,j,k)=1
          dropSize=dropSize+1
        else
          s_drop(i,j,k)=0
        endif
        !compute cell-centred normal
        ip=i+1
        jp=j+1
        kp=k+1
        if(ip.gt.nxt) ip=ip-nxt
        if(jp.gt.nyt) jp=jp-nyt
        if(kp.gt.nzt) kp=kp-nzt
        im=i-1
        jm=j-1
        km=k-1
        if(im.lt.1) im=im+nxt
        if(jm.lt.1) jm=jm+nyt
        if(km.lt.1) km=km+nzt
        nor(1,i,j,k)=(phase(ip,j,k)-phase(im,j,k))*(0.5/dx)
        nor(2,i,j,k)=(phase(i,jp,k)-phase(i,jm,k))*(0.5/dy)
        nor(3,i,j,k)=(phase(i,j,kp)-phase(i,j,km))*(0.5/dz)
        modnor=sqrt(nor(1,i,j,k)**2+nor(2,i,j,k)**2+nor(3,i,j,k)**2)
        !outward pointing normal
        nor(:,i,j,k)=-nor(:,i,j,k)/modnor
      enddo
    enddo
  enddo
  do k=1,nzt
    do j=1,nyt
      do i=1,nxt
        im=i-1
        jm=j-1
        km=k-1
        if(im.lt.1) im=im+nxt
        if(jm.lt.1) jm=jm+nyt
        if(km.lt.1) km=km+nzt
        ip=i+1
        jp=j+1
        kp=k+1
        if(ip.gt.nxt) ip=ip-nxt
        if(jp.gt.nyt) jp=jp-nyt
        if(kp.gt.nzt) kp=kp-nzt
        !compute curvature in backward direction
        kur(i,j,k)=(nor(1,ip,j,k)-nor(1,im,j,k))*(0.5/dx)+ &
                   (nor(2,i,jp,k)-nor(2,i,jm,k))*(0.5/dy)+ &
                   (nor(3,i,j,kp)-nor(3,i,j,km))*(0.5/dz)
      enddo
    enddo
  enddo
  !calculate drop statistics if minimum i j k in the drop lies in the processor's domain
  invSize = 1.0/dropSize
  dropArea=0.0
  kurMean=0.0
  kurStdDev=0.0
  kurInv=0.0
  kurInvSq=0.0
  faces=0
  edges=0
  vertices=0
  do ii=1,3
    hist(ii)%v=0.0
  enddo
  dropVel=0.0
  !make histograms of drop mass in each direction so we can avoid overlap with periodic 
  !boundaries in moment of inertia calculations
  do k=1,nzt
    do j=1,nyt
      do i=1,nxt
        if(s_drop(i,j,k).ne.0) then
          hist(1)%v(i)=hist(1)%v(i)+invSize
          hist(2)%v(j)=hist(2)%v(j)+invSize
          hist(3)%v(k)=hist(3)%v(k)+invSize
          do ii=1,3
            do mom = 1,maxMom
              dropVel(mom,ii) = dropVel(mom,ii) + (vel(ii,i,j,k)**mom)*invSize
            enddo
          enddo
        endif
      enddo
    enddo
  enddo 
  dropPos=0.0
  last0=0
  last1=0
  do ii=1,3
    !find beginning and end of drop in x,y,z directions
    do i=1,nt(ii)
      ip=i+1
      if(ip.gt.nt(ii)) ip=ip-nt(ii)
      if(hist(ii)%v(i).gt.0.5*invSize.and.hist(ii)%v(ip).lt.0.5*invSize) last1(ii)=i
      if(hist(ii)%v(i).lt.0.5*invSize.and.hist(ii)%v(ip).gt.0.5*invSize) last0(ii)=i
    enddo
    !if part of the drop touches i=1, shift that part to other side of domain so drop is contiguous
    do i=1,nt(ii)
      iShifted = i
      if(last1(ii).lt.last0(ii).and.i.le.last1(ii)) iShifted = i+nt(ii)
      do mom = 1,maxMom
        dropPos(mom,ii) = dropPos(mom,ii) + hist(ii)%v(i)*(iShifted**mom)
      enddo
    enddo 
  enddo
  MoI=0.0
  do k=1,nzt
    pos(3)=k
    if(last1(3).lt.last0(3).and.k.le.last1(3)) pos(3) = k+nzt
    do j=1,nyt
      pos(2)=j
      if(last1(2).lt.last0(2).and.j.le.last1(2)) pos(2) = j+nyt
      do i=1,nxt
        pos(1)=i
        if(last1(1).lt.last0(1).and.i.le.last1(1)) pos(1) = i+nxt
        !use https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
        !Bunner & Tryggvason JFM 2003 misses out on diag part of MoI definition
        !don't bother moving dropPos inside domain until after the moment of inertia calculations
        if(s_drop(i,j,k).ne.0) then
          !contribution from each cell
          diag=dx**5*(1.0/6.0)
          do ii=1,3
            diag = diag + ( pos(ii)-dropPos(1,ii) )**2 *dx**5
          enddo
          do ii=1,3
            do jj=1,3
              MoI(ii,jj) = MoI(ii,jj) - ( pos(ii)-dropPos(1,ii) )*( pos(jj)-dropPos(1,jj) )*dx**5
              if(ii.eq.jj) MoI(ii,jj) = MoI(ii,jj) + diag
            enddo
          enddo
        endif
        !the interface is here if s_drop changes in any of the 6 directions
        onInt=0
        do ii=1,3
          ip = i
          jp = j
          kp = k
          iq = i
          jq = j
          kq = k
          if (ii.eq.1) then
            ip=i+1
            if(ip.gt.nxt)ip=ip-nxt
            iq=i-1
            if(iq.lt.1)  iq=iq+nxt
          endif
          if (ii.eq.2) then
            jp=j+1
            if(jp.gt.nyt)jp=jp-nyt
            jq=j-1
            if(jq.lt.1)  jq=jq+nyt
          endif
          if (ii.eq.3) then
            kp=k+1
            if(kp.gt.nzt)kp=kp-nzt
            kq=k-1
            if(kq.lt.1)  kq=kq+nzt
          endif
          if (s_drop(ip,jp,kp).ne.s_drop(i,j,k)) faces = faces + 1
          if (s_drop(ip,jp,kp).ne.s_drop(i,j,k).or.s_drop(iq,jq,kq).ne.s_drop(i,j,k)) onInt=1
        enddo
        if(onInt.eq.1)then
          !find direction which is most aligned with interface normal
          maxNor=0.0
          dA=0.0
          do ii=1,3
            if(abs(nor(ii,i,j,k)).gt.maxNor)then
              maxNor = abs(nor(ii,i,j,k))
              if(ii.eq.1)dA=dy*dz
              if(ii.eq.2)dA=dz*dx
              if(ii.eq.3)dA=dx*dy
            endif
          enddo
          !area of interface is 
          dA=dA/(2.0*maxNor)
          dropArea = dropArea + dA
          kurMean = kurMean + dA*kur(i,j,k)
          kurStdDev = kurStdDev + dA*kur(i,j,k)**2
          kurInv = kurInv + dA/kur(i,j,k)
          kurInvSq = kurInvSq + dA/kur(i,j,k)**2
        endif
        !I Cannon 2023 May
        !calculate the genus of the drop using Mendoza et al. 2006 Acta Mater 10.1016/j.actamat.2005.10.010.
        !genus = tunnels - voids 
        !for a shape made of simple polygons:
        !genus = 1 - (faces - edges + corners)/2
        !use a painter algorithm to count the number of contiguous drop corners at this vertex
        !paint the vof with up to 8 colours surrounding a vertex between cells
        cols=0
        neigh=0
        paint=0
        do kq=0,1
          kp = k + kq
          if(kp.gt.nzt) kp=kp-nzt
          do jq=0,1
            jp = j + jq
            if(jp.gt.nyt) jp=jp-nyt
            do iq=0,1
              ip = i + iq
              if(ip.gt.nxt) ip=ip-nxt
              neigh(iq,jq,kq)=s_drop(ip,jp,kp)
              if(s_drop(ip,jp,kp).ne.0) then
                cols=cols+1
                paint(iq,jq,kq)=cols
              endif
            enddo
          enddo
        enddo
        faceOnCorner=0
        !no corners if vertex entirely inside drop
        if(cols.eq.8) cols=0
        !trivial case with all drop except for one corner
        if(cols.eq.7) cols=1
        !if drop only at the corner, there are three faces touching it
        if(cols.eq.1) faceOnCorner=3
        !if just two voids on opposite sides of the vertex, we connect the voids using a tunnel with
        !Euler characteristic = v-e+f = 3-9+6 = 0 = 0-6+6. I choose v=0 to skip the painting loop below
        if(cols.eq.6) then
          cols=0
          faceOnCorner=6
          do jq=0,1
            do iq=0,1
              if(neigh(iq,jq,1).ne.neigh(1-iq,1-jq,0)) then
                cols=6
                faceOnCorner=0
              endif
            enddo
          enddo
        endif
        if(cols.gt.1) then
          !in the worst case, we need 4 iterations to paint the 2x2x2 cube by neighbours
          !loop 4 times instead of recursion
          do paintIt=1,4
            do kq=0,1
              do jq=0,1
                do iq=0,1
                  if(neigh(iq,jq,kq).gt.0) then
                    do ii=1,3
                      ip=iq
                      jp=jq
                      kp=kq
                      if(ii.eq.1) ip=1-iq
                      if(ii.eq.2) jp=1-jq
                      if(ii.eq.3) kp=1-kq
                      !paint any touching cells the same colour
                      if(neigh(iq,jq,kq).eq.neigh(ip,jp,kp)) then
                        if(paint(ip,jp,kp).ne.paint(iq,jq,kq)) then
                          paint(iq,jq,kq) = min(paint(iq,jq,kq), paint(ip,jp,kp))
                          paint(ip,jp,kp) = min(paint(iq,jq,kq), paint(ip,jp,kp))
                        endif
                      !if s_drop has a sign change, there is a face and an edge
                      else if(paintIt.eq.1) then
                        faceOnCorner=faceOnCorner+1
                      endif
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
          !count how many colours remain, this is equal to the number of drop corners at the vertex
          cols=0
          do paintIt=1,8
            colourLoop: do kq=0,1
              do jq=0,1
                do iq=0,1
                  if(paintIt.eq.paint(iq,jq,kq)) then
                    cols=cols+1
                    exit colourLoop
                  endif
                enddo
              enddo
            enddo colourLoop
          enddo
        endif
        vertices = vertices + cols
        !the number of edges touching the vertex is equal to the number of faces
        edges=edges+faceOnCorner
      enddo
    enddo
  enddo
  !every edge is between two corners so we double counted the edges
  edges = edges/2
  if(mod(vertices-edges+faces,2).ne.0)then
    write(*,*)'warning: Euler characteristic is odd, v,e,f=',vertices,edges,faces
  endif
  genus = 1-(vertices-edges+faces)/2
  kurMean = kurMean / dropArea
  !kurStdDev = dsqrt( kurStdDev/dropArea - kurMean**2 )
  kurStdDev = sqrt( kurStdDev/dropArea - kurMean**2 )
  kurInv = kurInv / dropArea
  kurInvSq = kurInvSq / dropArea
  !Calculate eigenvalues of moment of inertia
  if(last0(1).eq.0.or.last0(2).eq.0.or.last0(3).eq.0) then
    write(*,*) 'drop spans all of domain, deformation is undefined'
    deformation=-1.0
  else
    call SSYEV("V","U",3,MoI,3,MoIEiVals,work,8,error)
    if (error.ne.0) write(*,*) 'dropSize', dropSize, 'MoIEiVals', MoIEiVals
    deformation=sqrt(MoIEiVals(3)/MoIEiVals(1))
  endif
  do ii=1,3
    !move drop back inside domain
    if(dropPos(1,ii).gt.nt(ii)) then
      do mom=1,maxMom
        dropPos(mom,ii)=dropPos(mom,ii)-nt(ii)**mom
      enddo
    endif
    pos(ii) = int(dropPos(1,ii))
    !put in units of simulation domain size
    do mom=1,maxMom
      dropPos(mom,ii)=dropPos(mom,ii) * ( l(ii)/nt(ii) )**mom
    enddo
  enddo
  dVeldxBox=0.0
  jj=1
  ip = pos(1) + boxWid
  ip = mod(ip-1, nt(jj)) + 1
  im = pos(1) - boxWid
  im = mod(im-1, nt(jj)) + 1
  nVels = 0
  do k = pos(3)-boxWid, pos(3)+boxWid
    kp = mod(k-1, nzt) + 1
    do j = pos(2)-boxWid, pos(2)+boxWid
      jp = mod(j-1, nyt) + 1
      nVels = nVels + 1
      !dVeldx(:) = (vel(:,ip,jp,kp) - vel(:,im,jp,kp)) / dx / boxWid
      dVeldx(:) = vel(:,ip,jp,kp)
      dVeldxBox(:,jj) = dVeldxBox(:,jj) + ( dVeldx(:) - dVeldxBox(:,jj) ) / nVels
    enddo
  enddo
  write(*,*) 'nVels',nVels
  write(*,*) 'dVeldx',dVeldx
  write(*,*) 'dVeldxBox',dVeldxBox
  jj=2
  jp = pos(jj) + boxWid
  jp = mod(jp-1, nt(jj)) + 1
  jm = pos(jj) - boxWid
  jm = mod(jm-1, nt(jj)) + 1
  nVels = 0
  do k = pos(3)-boxWid, pos(3)+boxWid
    kp = mod(k-1, nzt) + 1
    do i = pos(1)-boxWid, pos(1)+boxWid
      ip = mod(i-1, nxt) + 1
      nVels = nVels + 1
      !dVeldx(:) = (vel(:,ip,jp,kp) - vel(:,ip,jm,kp)) / dy / boxWid
      dVeldx(:) = vel(:,ip,jp,kp)
      dVeldxBox(:,jj) = dVeldxBox(:,jj) + ( dVeldx(:) - dVeldxBox(:,jj) ) / nVels
    enddo
  enddo
  jj=3
  kp = pos(jj) + boxWid
  kp = mod(kp-1, nt(jj)) + 1
  km = pos(jj) - boxWid
  km = mod(km-1, nt(jj)) + 1
  nVels = 0
  do j = pos(2)-boxWid, pos(2)+boxWid
    jp = mod(j-1, nyt) + 1
    do i = pos(1)-boxWid, pos(1)+boxWid
      ip = mod(i-1, nxt) + 1
      nVels = nVels + 1
      !dVeldx(:) = (vel(:,ip,jp,kp) - vel(:,ip,jp,km)) / dz / boxWid
      dVeldx(:) = vel(:,ip,jp,kp)
      dVeldxBox(:,jj) = dVeldxBox(:,jj) + ( dVeldx(:) - dVeldxBox(:,jj) ) / nVels
    enddo
  enddo
  do ii=1,3
    do jj=1,3
      strainBox(ii,jj) = 0.5*(dVeldxBox(ii,jj) + dVeldxBox(jj,ii))
    enddo
  enddo
  !vorticity = curl(u) = del x u
  vort(1) = dVeldxBox(3,2) - dVeldxBox(2,3)
  vort(2) = dVeldxBox(1,3) - dVeldxBox(3,1)
  vort(3) = dVeldxBox(2,1) - dVeldxBox(1,2)
  do j = 1,3
    do i = 1,3
      !Q citerion as equation (4) from Paul22roleOfBreakup
      QInva = QInva - 0.5*dVeldxBox(i,j)*dVeldxBox(j,i)
      do k = 1,3
        !Third invariant of velocty grad tensor as equation (5) from Paul22roleOfBreakup
        !Assumes incompressibility so that 3*det(dVeldxBox)=tr(dVeldxBox^3)
        RInva = RInva - (1.0/3.0)*dVeldxBox(i,j)*dVeldxBox(j,k)*dVeldxBox(k,i) 
      enddo
    enddo
  enddo
  call SSYEV("V","U",3,strainBox,3,strainEiVals,work,8,error)
  write(statU,'(ES16.7E3,i16,4ES16.7E3,6i16,2ES16.7E3)') time,dropSize,dropArea,deformation,kurMean,kurStdDev,&
        last0(:),last1(:),kurInv,kurInvSq
  flush(statU)
  ! generate format string for writing 
  write(fmtstr,'(a,i0,a)') '(',1+maxMom*3,'(ES16.7E3))'
  write(posiU,fmtstr) time, (dropPos(mom,:),mom=1,maxMom)
  flush(posiU)
  write(veloU,fmtstr) time, (dropVel(mom,:),mom=1,maxMom)
  flush(veloU)
  write(MoInU,'(13ES16.7E3)') time, MoIEiVals, MoI
  flush(MoInU)
  write(dVelU,'(13ES16.7E3)') time, QInva, RInva, dVeldxBox
  flush(dVelU)
  write(straU,'(13ES16.7E3)') time, strainEiVals, strainBox
  flush(straU)
  write(topoU,'(ES16.7E3, 4i16)') time,vertices,edges,faces,genus
  flush(topoU)
  write(fmtstr,'(a,i0,a)') '(',nzt/2+2,'(ES16.7E3))'
  write(specU,fmtstr) time,ESpec 
  flush(specU)
  write(forcU,fmtstr) time,FSpec 
  write(vortU,'(4ES16.7E3)') time, vort
  flush(vortU)
enddo
do ii=1,3
  deallocate(hist(ii)%v)
enddo
close(posiU,status='keep')
close(veloU,status='keep')
close(MoInU,status='keep')
close(dVelU,status='keep')
close(straU,status='keep')
close(statU,status='keep')
close(topoU,status='keep')
close(specU,status='keep')
close(forcU,status='keep')
close(listU)
close(vortU)
call system("rm file_list.txt")
write(6,*) 'This is the end'
end program strain_box
