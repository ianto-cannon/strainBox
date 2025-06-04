program strain_box
use hdf5
!use, intrinsic :: ISO_FORTRAN_ENV
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'
type ragged_array
  !2D array containing vectors of different lengths
  real,allocatable::v(:)
  character(len=12) :: iStepChar
end type ragged_array
character(len=200) :: filename
real :: modnor
integer,parameter :: maxMom=2,statU=41,posiU=42,veloU=43,MoInU=44,counU=45,topoU=46
integer,parameter :: nxt=256, nyt=256, nzt=256
real, parameter :: pi=3.14159265358979
real,parameter :: lx=2*pi, ly=2*pi, lz=2*pi
real, dimension(3), parameter :: l = (/lx,ly,lz/)
real,parameter :: dx=lx/nxt, dy=ly/nxt, dz=lz/nxt
real, dimension(nxt,nyt,nzt) :: kur,u,v,w,phase,dxxPhase
real, dimension(nxt,nyt,nzt) :: chemPot,dxChemPot,dyChemPot,dzChemPot,surfFX,surfFY,surfFZ
real, dimension(nxt,nyt,nzt,3) :: nor, vel
!use int8 for non shared arrays to save memory
integer, dimension(nxt,nyt,nzt) :: s_drop
integer, dimension(0:1,0:1,0:1) :: paint, neigh
integer :: i,j,k,ip,jp,kp,iq,jq,kq,direc,iShifted,mom,last0(3),last1(3),pos(3),ii,jj,error,intR
integer :: cols,paintIt,faceOnCorner,genus,onInt,ky,kz,im,jm,km
real, dimension(maxMom,3) :: dropPos, dropVel
real :: MoI(3,3), work(8), eiVals(3)
real :: diag, deformation, dropArea, dA, Cn, r
real :: maxNor, kurMean, kurStdDev, kurInv, kurInvSq, invSize
type(ragged_array) :: hist(3) !histogram of drop mass in x, y and z directions
!integer(kind=int64) :: dropSize, faces, edges, vertices
integer :: dropSize, faces, edges, vertices
character(len=200) :: fileEnd, fmtstr
integer, dimension(3), parameter :: nt = (/nxt,nyt,nzt/)
integer(hid_t) :: file_id, dset_id
integer(hsize_t) :: dims(3)=(/nxt,nyt,nzt/),  dims1d(1)=(/1/) 
complex, dimension(nxt/2+1,nyt,nzt) :: cHat,chemPotHat,dxChemPotHat,dyChemPotHat,dzChemPotHat
complex, dimension(nxt/2+1,nyt,nzt) :: surfFXHat,surfFYHat,surfFZHat,uHat,vHat,wHat
real, dimension(0:nzt/2) :: FSpec,ESpec
type(C_PTR)  :: plan, plan_inverse
write(*,'(1x,a)') '                 starting number of drops calculation                       '
! Open the file (read-only)
call h5open_f(error)
call h5fopen_f("/home/alberto.velamartin/drop_time/we_10/run_break_009/field.015.h5", &
                 H5F_ACC_RDONLY_F, file_id, error)
  call h5dopen_f(file_id, "Cn", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
  call h5dclose_f(dset_id, error)
  write(*,*) 'Cn', Cn
  call h5dopen_f(file_id, "We", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
  call h5dclose_f(dset_id, error)
  write(*,*) 'We', Cn
  call h5dopen_f(file_id, "c", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_REAL, phase, dims, error)
  call h5dclose_f(dset_id, error)
  call h5dopen_f(file_id, "res", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
  call h5dclose_f(dset_id, error)
  write(*,*) 'res', Cn
  call h5dopen_f(file_id, "time", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
  call h5dclose_f(dset_id, error)
  write(*,*) 'time', Cn
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
write(*,*)'E',0.5*sum(u**2+v**2+w**2)
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
write(*,*)'ESpec',sum(ESpec)
filename = './output/spec.txt'
write(6,*) 'saving to ',trim(filename)
open(30,file=trim(filename),form='formatted',action='write')
  do i=0,nzt/2
    write(30,'(3E23.15)') 1.0*i,ESpec(i),FSpec(i)
  enddo
close(30,status='keep')
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
      nor(i,j,k,1)=(phase(ip,j,k)-phase(im,j,k))*(0.5/dx)
      nor(i,j,k,2)=(phase(i,jp,k)-phase(i,jm,k))*(0.5/dy)
      nor(i,j,k,3)=(phase(i,j,kp)-phase(i,j,km))*(0.5/dz)
      !modnor=dsqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2)
      modnor=sqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2)
      !outward pointing normal
      nor(i,j,k,:)=-nor(i,j,k,:)/modnor
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
      kur(i,j,k)=(nor(ip,j,k,1)-nor(im,j,k,1))*(0.5/dx)+ &
                 (nor(i,jp,k,2)-nor(i,jm,k,2))*(0.5/dy)+ &
                 (nor(i,j,kp,3)-nor(i,j,km,3))*(0.5/dz)
    enddo
  enddo
enddo
fileEnd='.dat'
open(statU,file='./output/statDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(statU,status='keep')
open(posiU,file='./output/posiDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(posiU,status='keep')
open(veloU,file='./output/veloDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(veloU,status='keep')
open(MoInU,file='./output/MoInDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(MoInU,status='keep')
open(topoU,file='./output/topoDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(topoU,status='keep')
dropPos=0.0
dropVel=0.0
do direc=1,3
  allocate(hist(direc)%v(nt(direc)))
enddo
hist(1)%iStepChar='1'
hist(2)%iStepChar='2'
hist(3)%iStepChar='3'
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
do direc=1,3
  hist(direc)%v=0.0
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
        do direc=1,3
          do mom = 1,maxMom
            dropVel(mom,direc) = dropVel(mom,direc) + (vel(i,j,k,direc)**mom)*invSize
          enddo
        enddo
      endif
    enddo
  enddo
enddo 
write(*,*) "invSize", invSize
write(*,*) "hist1", hist(1)%v
dropPos=0.0
last0=0
last1=0
do direc=1,3
  !find beginning and end of drop in x,y,z directions
  do i=1,nt(direc)
    ip=i+1
    if(ip.gt.nt(direc)) ip=ip-nt(direc)
    if(hist(direc)%v(i).gt.0.5*invSize.and.hist(direc)%v(ip).lt.0.5*invSize) last1(direc)=i
    if(hist(direc)%v(i).lt.0.5*invSize.and.hist(direc)%v(ip).gt.0.5*invSize) write(*,*) 'last0',direc,i,ip
    if(hist(direc)%v(i).lt.0.5*invSize.and.hist(direc)%v(ip).gt.0.5*invSize) last0(direc)=i
  enddo
  !if part of the drop touches i=1, shift that part to other side of domain so drop is contiguous
  do i=1,nt(direc)
    iShifted = i
    if(last1(direc).lt.last0(direc).and.i.le.last1(direc)) iShifted = i+nt(direc)
    do mom = 1,maxMom
      dropPos(mom,direc) = dropPos(mom,direc) + hist(direc)%v(i)*(iShifted**mom)
    enddo
  enddo 
enddo
write(*,*) '296last0',last0
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
      do direc=1,3
        ip = i
        jp = j
        kp = k
        iq = i
        jq = j
        kq = k
        if (direc.eq.1) then
          ip=i+1
          if(ip.gt.nxt)ip=ip-nxt
          iq=i-1
          if(iq.lt.1)  iq=iq+nxt
        endif
        if (direc.eq.2) then
          jp=j+1
          if(jp.gt.nyt)jp=jp-nyt
          jq=j-1
          if(jq.lt.1)  jq=jq+nyt
        endif
        if (direc.eq.3) then
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
        do direc=1,3
          if(abs(nor(i,j,k,direc)).gt.maxNor)then
            maxNor = abs(nor(i,j,k,direc))
            if(direc.eq.1)dA=dy*dz
            if(direc.eq.2)dA=dz*dx
            if(direc.eq.3)dA=dx*dy
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
                  do direc=1,3
                    ip=iq
                    jp=jq
                    kp=kq
                    if(direc.eq.1) ip=1-iq
                    if(direc.eq.2) jp=1-jq
                    if(direc.eq.3) kp=1-kq
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
do direc=1,3
  !move drop back inside domain
  if(dropPos(1,direc).gt.nt(direc)) then
    do mom=1,maxMom
      dropPos(mom,direc)=dropPos(mom,direc)-nt(direc)**mom
    enddo
  endif
  !put in units of simulation domain size
  do mom=1,maxMom
    dropPos(mom,direc)=dropPos(mom,direc) * ( l(direc)/nt(direc) )**mom
  enddo
enddo
!Calculate eigenvalues of moment of inertia
if(last0(1).eq.0.or.last0(2).eq.0.or.last0(3).eq.0) then
  write(*,*) 'drop spans all of domain, deformation is undefined'
  deformation=-1.0
else
  call DSYEV("N","U",3,MoI,3,eiVals,work,8,error)
  if (error.ne.0) write(*,*) 'dropSize', dropSize, 'eiVals', eiVals
  deformation=sqrt(eiVals(3)/eiVals(1))
endif
open(statU,file='./output/statDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
  write(statU,'(i16,4ES16.7E3)') dropSize,dropArea,deformation,kurMean,kurStdDev
close(statU,status='keep')
! generate format string for writing 
write(fmtstr,'(a,i0,a)') '(i16,',maxMom*3,'(ES16.7E3))'
open(posiU,file='./output/posiDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
  write(posiU,fmtstr) dropSize,(dropPos(mom,:),mom=1,maxMom)
close(posiU,status='keep')
open(veloU,file='./output/veloDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
  write(veloU,fmtstr) dropSize,(dropVel(mom,:),mom=1,maxMom)
close(veloU,status='keep')
open(MoInU,file='./output/MoInDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
  !write(MoInU,'(i16,9ES16.7E3,6i16,2ES16.7E3)') dropSize,(eiVals(i),i=1,3),&
  !  MoI(1,1),MoI(1,2),MoI(1,3),MoI(2,2),MoI(2,3),MoI(3,3),last0(:),last1(:),kurInv,kurInvSq
  write(*,*) '519last0',last0(:)
  write(MoInU,*) last0(:)
close(MoInU,status='keep')
open(topoU,file='./output/topoDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
  write(topoU,'(5i16)') dropSize,vertices,edges,faces,genus
close(topoU,status='keep')
write(6,*) 'This is the end'
end program strain_box
